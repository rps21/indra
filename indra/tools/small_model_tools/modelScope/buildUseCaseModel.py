from copy import deepcopy
from pick import pick
import tempfile
from indra.statements import * 
from indra.mechlinker import MechLinker
from indra.preassembler import Preassembler
from indra.assemblers import PysbAssembler
from indra.assemblers import SifAssembler
import paths_graph as pg
from indra.explanation import model_checker
from indra.tools import assemble_corpus as ac
from indra.tools.small_model_tools.modelScope import indraDB_query as idb
from indra.tools.small_model_tools.modelContext import buildSmallModel as bsm
from indra.tools.small_model_tools.modelContext import extraModelReductionTools as ex
from indra.tools.small_model_tools.modelContext import enforceCascadeContext as cs
from indra.tools.small_model_tools.modelContext import addImplicitMechs as aim
import pysb
import csv 

import logging
#logging.getLogger("assemble_corpus").setLevel(logging.WARNING)
logging.getLogger("pre_cfpg").setLevel(logging.WARNING)
logging.getLogger("model_checker").setLevel(logging.WARNING)
logging.getLogger("paths_graph").setLevel(logging.WARNING)
logging.getLogger("grounding_mapper").setLevel(logging.WARNING)
logging.getLogger("preassembler").setLevel(logging.WARNING)
logging.getLogger("pysb_assembler").setLevel(logging.WARNING)

#TODO:

#Build sif file
def buildDirectedSif(stmts,drugStmts):

    #Apply small model simplification to list of statements 
    #TEMP HARD CODING
    mutations = {'BRAF': [('V', '600', 'E')]}
    modelStmts = ac.filter_mutation_status(stmts,mutations,deletions=[])
    modelStmts = bsm.buildSmallModel(modelStmts)
    modelStmts = modelStmts + drugStmts
    modelStmts = cs.run_mechlinker_step(modelStmts)

    #build sif graph 
    sa = SifAssembler(modelStmts)
    sifModel = sa.make_model(use_name_as_key=True, include_mods=True,include_complexes=True)
    directedSifModel = sa.print_model(include_unsigned_edges=True)
    tempRawSif = tempfile.NamedTemporaryFile(mode='w+t',delete=True)
    sa.save_model(fname=tempRawSif.name,include_unsigned_edges=True)

    newSifModelText = []
    sifModelText = tempRawSif.readlines()

    #For pg: 0 = positive, 1 = negative
    #In raw sif: 1 = positive, 0 = neutral, -1 = negative
    #below, keeping 0 the same turns neutral to positive. Is this best rule of thumb?
    for line in sifModelText:       
        if ' 1 ' in line:
            newline = line.replace(' 1 ', ' 0 ')
        elif ' -1 ' in line:
            newline = line.replace(' -1 ', ' 1 ')
        else:
            newline = line    
        newSifModelText.append(newline)
    tempRawSif.close()

    directedSif = tempfile.NamedTemporaryFile(mode='w+t',delete=False)
    directedSifFN = directedSif.name
    directedSif.writelines(newSifModelText)
    directedSif.close()

    return directedSifFN, modelStmts


#Build PG and find paths
def findPaths(sifFile,source,target,length):

    #Make num_samples variable? Best practice?
    g = pg.api.load_signed_sif(sifFile)
    #This checks if there is any path whatsoever, positive or negative, and reports all nodes 
    paths1 = pg.api.sample_paths(g,source,target,max_depth=length,num_samples=10,cycle_free=True,signed=False,target_polarity=0) 
    paths2 = pg.api.sample_paths(g,source,target,max_depth=length,num_samples=10,cycle_free=True,signed=False,target_polarity=1) 
    paths = paths1+paths2
    os.remove(sifFile)

    return paths

#########################
#Build model from paths
def buildModel(paths,stmts,drug,nodes=[]):
    uniqueNodes = []
    for path in paths:
        uniqueNodes = uniqueNodes + [el for el in list(path) if el not in uniqueNodes]
    
    necessaryNodes = ['GenericAgent'] + nodes  #Allow for specification of nodes that must be in a model, even if they aren't found on a path. ir
    uniqueNodes = uniqueNodes + necessaryNodes
    modelStmts = ac.filter_gene_list(stmts,uniqueNodes,'all')

    recLigNamePairs = cs.find_rec_lig(modelStmts)[0]
    recLigNames = list(itertools.chain.from_iterable(recLigNamePairs))
    uniqueNodes = uniqueNodes + recLigNames    #Additionally add back any missing ligands, so a receptor always has a pa
    modelStmts = ac.filter_gene_list(stmts,uniqueNodes,'all')

    modelStmts = Preassembler.combine_duplicate_stmts(modelStmts)  
    pa = PysbAssembler()
    pa.add_statements(modelStmts)
    PySB_Model = pa.make_model()

    return PySB_Model, modelStmts

##########################
#Test experimental findings against candidate model
def testExpStmt(model,expStmt): #take in statements? Sentences?

    testStmts = [expStmt]
    mc = model_checker.ModelChecker(model=model,statements=testStmts)
    results = mc.check_model(max_path_length=15)

    return results, mc


def runModelCheckingRoutine(stmts,drug,drugStmts,nodeToExplain,expStmts,extraNodes):
    fn, contextStmts = buildDirectedSif(stmts,drugStmts)
    paths = findPaths(fn,drug,nodeToExplain,15)

    if paths:
        extraNodes.append(drug)
        PySB_Model,modelStmts = buildModel(paths,contextStmts,drug,extraNodes)
        ac.dump_statements(modelStmts,'testingStmts.pkl')
        results, mc = testExpStmt(PySB_Model,expStmts)
        print(results)
        print(paths)

        passResult = results[0][1].path_found
        if passResult:
            pathNodes = list(paths[0])
            pathStmts = ac.filter_gene_list(modelStmts,pathNodes,'all')

            allEvList = []
            for st in pathStmts:
                if st.evidence:
                    evList = [st,st.evidence[0].source_api,st.evidence[0].pmid,st.evidence[0].text]
                    allEvList.append(evList)
            with open("evidence.csv", "a") as f:
                writer = csv.writer(f)
                writer.writerows(allEvList)
#        ac.dump_statements(modelStmts,'testingStmts.pkl')
    else:
        passResult = None
        modelStmts = None
    return passResult, modelStmts

###########################

def findActiveForm(node,currentNodes):
    newAFStmt = None
    #find all non-mutation active forms for newly added node 
    stmtsDB_AF = idb.queryDB_nonmod(node,'activeform')
    stmtsDB_AF = ex.removeMutations(stmtsDB_AF)
    stmtsDB_AF = [st for st in stmtsDB_AF if st.is_active]
    stmtsDB = idb.cleanStatements(stmtsDB_AF)

    #provide list of af's to follow up on?
    #Ideally, check for a bound condition, and see if binder is in model already. If not, check for phos. If neither...
    #Will rework this whole section in the futre

    for st in stmtsDB_AF:
        if st.agent.bound_conditions:
            boundName = st.agent.bound_conditions[0].agent.name
            print(boundName)
            
            print(currentNodes)
            if boundName in currentNodes:
                print('YES')
                newComplexStmt = Complex([st.agent,st.agent.bound_conditions[0].agent])
                newAFStmt = [st,newComplexStmt]
                print("NEW AF STMT")
                print(newAFStmt)
                break
    if not newAFStmt:
        for st in stmtsDB_AF:
            if st.agent.mods:
                if st.agent.mods[0].mod_type == 'phosphorylation':  

                    newMod = 'phosphorylation'    
                    #Later, filter and feedback instead of generalizing   
                    st_general = deepcopy(st)
                    st_general.agent.mods[0].residue = None    
                    st_general.agent.mods[0].position = None       

                    newAFStmt = [st_general]       #Make sure to take a stmt where is_active=True
                    print('NEWAFSTMT 2')
                    print(newAFStmt)
                    break

    if not stmtsDB_AF:
        stmtsDB_AF = idb.queryDB_nonmod(node,'increaseamount')    
        #pick tf
        candidateNodes = idb.getCandidateNodes(stmtsDB_AF)
        #Add exit point
        candidateNodes.append('exit')
        #allow user to pick node to test
        title = 'Pick the TF to follow up on, for %s of %s' % ('transcription',node)
        if candidateNodes:
            option,index = pick(candidateNodes,title,indicator='=>',default_index=0)    #This really sucks
        else:
            print('No options founds')
            found = 1
        #Handle selected option 
        if option == 'exit':
            found = 1
        newAFStmt = ac.filter_gene_list(stmtsDB_AF,[option],'one')
        print('NEWAFSTMT 3')
        print(newAFStmt)
    return newAFStmt


##############################################

def expandModel(expObservations,drug,drugTargets,drugStmts,initialStmts=None,initialNodes=None,extraNodes=[]):
    if initialStmts:
        currentStmts = initialStmts
        currentNodes = initialNodes
    else:
        currentStmts = []
        currentNodes = []

    for key in expObservations:
        #Find better way to track stmts 
        try:
            currentStmts = finalStmts
            currentStmts = [st for st in currentStmts if st not in drugStmts]
        except NameError:
            finalStmts = [] 

        finalNode = key
        nodeToExplain = key
        modToExplain = expObservations[key][0]
        expStmt = expObservations[key][1]
        passResult, modelStmts = runModelCheckingRoutine(currentStmts,drug,drugStmts,nodeToExplain,expStmt,extraNodes)

        #Model doesn't satisfy condition, go searching for statements to add to model 
        if passResult:
            finalStmts = modelStmts + finalStmts
        else:
            #modelStmts = modelStmts.remove(drugStmts)
            found = 0
            while found == 0:

                #search for potential nodes and stmts to add to the model 
                stmtsDB = idb.queryDB(nodeToExplain,modToExplain)
                stmtsDB = idb.cleanStatements(stmtsDB)
                #identify all new nodes that are candidates for the model 
                candidateNodes = idb.getCandidateNodes(stmtsDB)
                candidateNodes = list(set(candidateNodes))
                currentNodes = list(set(currentNodes))
                #check if any candidate nodes are already in the model
                for node in currentNodes:
                    if node in candidateNodes:
                        #If yes, check if satisfy model checker 
                        currentNodes.append(node)
                        print('Testing new node %s' % node)
                        print(ac.filter_gene_list(stmtsDB,node,'one'))
                        testStmts = currentStmts + ac.filter_gene_list(stmtsDB,node,'one')  #+ prior filtered by new node?

                        passResult, modelStmts = runModelCheckingRoutine(testStmts,drug,drugStmts,finalNode,expStmt,extraNodes)      
                        if passResult:
                            found = 1
                            finalStmts = aim.addAll(modelStmts) + finalStmts
                            finalStmts = Preassembler.combine_duplicate_stmts(finalStmts)
                            print('Worked!')
                            break
                        else:
                            currentNodes.pop()  
#                            modelStmts = modelStmts.remove(drugStmts)

                #If we haven't found a node to complete the model yet
                if found == 0:  
                    #Add exit point
                    candidateNodes.append('exit')
                    #allow user to pick node to test
                    title = 'Pick the protein to follow up on, for %s on %s' % (modToExplain,nodeToExplain)
                    if candidateNodes:
                        option,index = pick(candidateNodes,title,indicator='=>',default_index=0)    #This really sucks
                    else:
                        print('No options founds')
                        found = 1
                    #Handle selected option 
                    if option == 'exit':
                        found = 1

                    else:
                        currentNodes.append(option)
                        print('Testing new node %s' % option)
                        testStmts = currentStmts + ac.filter_gene_list(stmtsDB,node,'one')  
                        print(ac.filter_gene_list(stmtsDB,node,'one'))
                        passResult, modelStmts = runModelCheckingRoutine(testStmts,drug,drugStmts,finalNode,expStmt,extraNodes)      
                        if passResult:
                            found = 1
                            finalStmts = modelStmts
                            print('Worked!')
                            break
                        else:
                            #specify stmts, nodes, mods for next loop iteration
                            #This finds all relevant statements returned by DB for the new node. Will often be one, but should likely limit to one carefully selected statement
                            stmtsForNewNode = ac.filter_gene_list(stmtsDB,option,'one')  
                            if not stmtsForNewNode:
                                stmtsForNewNode = []
                            newAFStmt = findActiveForm(option,currentNodes)                   
                            currentStmts = currentStmts + stmtsForNewNode + newAFStmt
                            if len(newAFStmt) == 1:
                                if isinstance(newAFStmt[0],ActiveForm):
                                    currentNodes.append(option) #Is this getting all nodes? With new bound one?
                                    nodeToExplain = option  
                                    modToExplain = 'phosphorylation' #WRONG
                                elif isinstance(newAFStmt[0],IncreaseAmount):
                                    currentNodes.append(option) #
                                    nodeToExplain = newAFStmt[0].subj.name
                                    currentNodes.append(nodeToExplain) #
                                    newAFStmt = findActiveForm(nodeToExplain,currentNodes)                   
                                    currentStmts = currentStmts + stmtsForNewNode + newAFStmt
                            else:           
                                while len(newAFStmt) > 1:                     
                                    currentNodes.append(option)
                                    nodeToExplain = ac.filter_by_type(newAFStmt,ActiveForm)[0].agent.bound_conditions[0].agent.name
                                    currentNodes.append(nodeToExplain) #Is this getting all nodes? With new bound one?
                                    newAFStmt = findActiveForm(nodeToExplain,currentNodes)                   
                                    currentStmts = currentStmts + stmtsForNewNode + newAFStmt
                                modToExplain = 'phosphorylation' #WRONG
                                    #recheck model, get new af

    return finalStmts              





