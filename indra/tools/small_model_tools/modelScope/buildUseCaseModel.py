from copy import deepcopy
from pick import pick
from indra.statements import * 
from indra.mechlinker import MechLinker
from indra.preassembler import Preassembler
from indra.assemblers import PysbAssembler
from indra.assemblers import SifAssembler
import paths_graph as pg
from indra.explanation import model_checker
from indra.tools import assemble_corpus as ac
#from modelContext import enforceCascadeContext as cs
#from modelContext import combinePhosphorylationSites as ptm
from indra.tools.small_model_tools.modelScope import indraDB_query as idb
from indra.tools.small_model_tools.modelContext import buildSmallModel as bsm
import pysb

import logging
logging.getLogger("assemble_corpus").setLevel(logging.WARNING)
logging.getLogger("pre_cfpg").setLevel(logging.WARNING)
logging.getLogger("model_checker").setLevel(logging.WARNING)
logging.getLogger("paths_graph").setLevel(logging.WARNING)
logging.getLogger("grounding_mapper").setLevel(logging.WARNING)
logging.getLogger("preassembler").setLevel(logging.WARNING)
logging.getLogger("pysb_assembler").setLevel(logging.WARNING)

#TODO:
#Remove unused mods that exist on stmts, ex: Complex(HIF1A(mods: (sumoylation, phos_act)), MYC(mods: (ubiquitination, phos_act))),
#Look for bugs in phosphorylation sites being generated.



#Build sif file
def buildDirectedSif(stmts,save=True,fn='./directedSifFile.sif'):

    #Apply small model simplification to list of statements 
    modelStmts = bsm.buildSmallModel(stmts)
#    finalSorafStmts = ac.load_statements('finalSorafStmts_reduced.pkl')
#    finalStmts = modelStmts + finalSorafStmts

    #build sif graph 
    sa = SifAssembler(modelStmts)
    sifModel = sa.make_model(use_name_as_key=True, include_mods=True,include_complexes=True)
    directedSifModel = sa.print_model(include_unsigned_edges=True)
    sa.save_model(fname='./rawSif_tmp.sif',include_unsigned_edges=True)

    newSifModelText = []
    with open('./rawSif_tmp.sif','r') as f:
        sifModelText = f.readlines()
    #For pg: 0 = positive, 1 = negative
    #In raw sif: 1 = positive, 0 = neutral, -1 = negative
    #below, keeping 0 the same turns neutral to positive. Is this best rule of thumb?
    for line in sifModelText:       #check here
        if ' 1 ' in line:
            newline = line.replace(' 1 ', ' 0 ')
        elif ' -1 ' in line:
            newline = line.replace(' -1 ', ' 1 ')
        else:
            newline = line    
        newSifModelText.append(newline)

    if save:
        with open(fn,'w') as f:
            for line in newSifModelText:
                f.write(line)

    return fn, modelStmts


#######################

#Build PG and find paths

def findPaths(sifFile,source,target,length):

    #Make num_samples variable? Best practice?
    g = pg.api.load_signed_sif(sifFile)
    #This checks if there is any path whatsoever, positive or negative, and reports all nodes 
    #is it better to be more mindful of polarities?
    #depending on how/when it's called I can have an 3-part option to enforece either polarity or include both 
    paths1 = pg.api.sample_paths(g,source,target,max_depth=length,num_samples=10,cycle_free=True,signed=False,target_polarity=0) 
    paths2 = pg.api.sample_paths(g,source,target,max_depth=length,num_samples=10,cycle_free=True,signed=False,target_polarity=1) 
    paths = paths1+paths2

    return paths

#########################
#Build model from paths

#def buildModel(paths,stmts,drugStmts,save=False,fn='./modelStmts.pkl'):
def buildModel(paths,stmts,nodes=None,save=False,fn='./modelStmts.pkl'):
    uniqueNodes = []
    for path in paths:
        uniqueNodes = uniqueNodes + [el for el in list(path) if el not in uniqueNodes]

    #Major bug in model building - fails if no ligand 
#    #Need to handle this better, but for now just add relevant ligands
#    #Also need drug and newly added phosphatase in the final model, not sure best place to add that, here for now 
#    ligands = ['PDGF','FLT3LG','PDGFA','SORAFENIB','GenericPhosphatase']
#    uniqueNodes = uniqueNodes + ligands

    #print('Nodes for final model are %s' % uniqueNodes)
    modelStmts = ac.filter_gene_list(stmts,uniqueNodes,'all')


    if save:
        ac.dump_statements(modelStmts,fn)

    modelStmts = Preassembler.combine_duplicate_stmts(modelStmts)   #check here


    pa = PysbAssembler()
    pa.add_statements(modelStmts)
    PySB_Model = pa.make_model()

    return PySB_Model, modelStmts

##########################
#Test experimental findings against candidate model
def testExpStmt(model,expStmt): #take in statements? Sentences?

    testStmts = [expStmt]
    mc = model_checker.ModelChecker(model=model,statements=testStmts)
    results = mc.check_model(max_path_length=8)

    return results, mc


#drug targets and cogante ligands are build into funtions. These should be pulled out and made input variables. 
def runModelCheckingRoutine(stmts,drug,nodeToExplain,expStmts):
    fn, contextStmts = buildDirectedSif(stmts)
    paths = findPaths(fn,drug,nodeToExplain,8)
    if paths:
        PySB_Model,modelStmts = buildModel(paths,contextStmts)
        results, mc = testExpStmt(PySB_Model,expStmts)
        passResult = results[0][1].path_found
    else:
        passResult = None
        modelStmts = None
    return passResult, modelStmts

###########################








def findActiveForm(node,currentNodes):
    newAFStmt = None
    #find all non-mutation active forms for newly added node 
    stmtsDB_AF = idb.queryDB_nonmod(node,'activeform')
    stmtsDB_AF = remove_mutations(stmtsDB_AF)
    stmtsDB_AF = [st for st in stmtsDB_AF if st.is_active]

    #provide list of af's to follow up on?
    #Ideally, check for a bound condition, and see if binder is in model already. If not, check for phos. If neither...
    #Will rework this whole section in the futre

    for st in stmtsDB_AF:

        if st.agent.bound_conditions:
            boundName = newAF[-1].agent.bound_conditions[0].agent.name
            if boundName in currentNodes:
                newAFStmt = st
        elif st.agent.mods:
            if st.agent.mods[0].mod_type == 'phosphorylation':  

                newMod = 'phosphorylation'    
                #Later, filter and feedback instead of generalizing   
                st_general = deepcopy(st)
                st_general.agent.mods[0].residue = None    
                st_general.agent.mods[0].position = None       

                newAFStmt = st_general       #Make sure to take a stmt where is_active=True
                break

    return newAFStmt

##############################################
#Move somewhere else, need now for filtering DB results because I only care about wild type cells currently
def remove_mutations(stmts_in):
    stmts_in = Preassembler.combine_duplicate_stmts(stmts_in)
    stmts_new = stmts_in[:]
    for st in stmts_in:
        for ag in st.agent_list():
                if ag.mutations:
                    stmts_new.remove(st)
    return stmts_new

##############################################



#NEW ISSUE:
#Not all phos are equal. 
#Should filter AF stmts on ones that are true, then filter phos stmts on ones that match these phos sites 
#potential work around: when adding statement, remove res and pos.
#This loses accuracy and generally suck, but would probably lead to functional model in short term. 

#NEED to add failure after certain number of iterations that will present model for examination, can restart and choose different paths. 


def expandModel(expObservations,drug,drugTargets,initialStmts=None,initialNodes=None):
    finalStmts = []
    if initialStmts:
        currentStmts = initialStmts
        currentNodes = initialNodes
    else:
        currentStmts = []
        currentNodes = []

    for key in expObservations:
        finalNode = key
        nodeToExplain = key
        modToExplain = expObservations[key][0]
        expStmt = expObservations[key][1]
        passResult, modelStmts = runModelCheckingRoutine(currentStmts,drug,nodeToExplain,expStmt)

        #Model doesn't satisfy condition, go searching for statements to add to model 
        if passResult:
            finalStmts = finalStmts + modelStmts
        else:
            found = 0
            while found == 0:
                #search for potential nodes and stmts to add to the model 
                stmtsDB = idb.queryDB(nodeToExplain,modToExplain)
                #identify all new nodes that are candidates for the model 
                candidateNodes = idb.getCandidateNodes(stmtsDB)
                candidateNodes = list(set(candidateNodes))

                targetNodes = [nd for nd in candidateNodes if nd in drugTargets]    #Want to ensure any drug targets are in the potential node list, they are included before the cut off
                if len(candidateNodes) >= 20:
                    candidateNodes = candidateNodes[:21] + targetNodes
                #Add exit point
                candidateNodes.append('exit')

                #Test adding each new node, and all of it's
                #See if model passes model checker test with any added node 
                #Could get around slowness of adding lots of nodes by adding one node at a time, then testing.  But then could potentially run the test many, many times 

                    

                #If no single node satisfies condition, pick one node to add to model. 

                    #allow user to pick node for next iteration. Following all nodes will result in a combinatorial explosion
                proteinOptions = candidateNodes
                title = 'Pick the protein to follow up on, for %s on %s' % (modToExplain,nodeToExplain)
                #allMechs = allMechs + '\n' + title
                if proteinOptions:
                    option,index = pick(proteinOptions,title,indicator='=>',default_index=0)    #This really sucks
                else:
                    option = currentNodes[-2]   
                if option == 'exit':
                    found = 1
                else:
                    currentNodes.append(option)
                    print('Testing new node %s' % option)
                    testStmts = ac.filter_gene_list(currentStmts,currentNodes,'all') + ac.filter_gene_list(stmtsDB,option,'one')  
                    #print(ac.filter_gene_list(stmtsDB,node,'one') )
    #                    passResult, modelStmts = runModelCheckingRoutine(testStmts,drug,nodeToExplain,sentence)      
                    passResult, modelStmts = runModelCheckingRoutine(testStmts,drug,finalNode,expStmt)      
                    if passResult:
                        found = 1
                        finalStmts = finalStmts + modelStmts
                        print('Worked!')
                        break

                    if found == 0:
                        #specify stmts, nodes, mods for next loop iteration
                        #This finds all relevant statements returned by DB for the new node. Will often be one, but should likely limit to one carefully selected statement
                        stmtsForNewNode = ac.filter_gene_list(stmtsDB,option,'one')  
                        newAFStmt = findActiveForm(option,currentNodes)                   
                        currentStmts = currentStmts + stmtsForNewNode + [newAFStmt]
                        currentNodes.append(option)

                        nodeToExplain = option  
                        modToExplain = 'phosphorylation' #This need to be fixed/generalized



    return finalStmts              



############





finalStmts = expandModel(expObservations,'SORAFENIB',['FLT3','KDR','PDGFRA'],initialStmts=filteredStmts,initialNodes=nodes)







