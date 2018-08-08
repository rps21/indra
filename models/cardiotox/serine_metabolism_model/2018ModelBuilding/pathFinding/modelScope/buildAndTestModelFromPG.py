
from indra.assemblers import PysbAssembler
from indra.assemblers import SifAssembler
import paths_graph as pg
from indra.tools import assemble_corpus as ac
from indra.tools.small_model_tools import enforceCascadeContext as cs
from indra.tools.small_model_tools import combinePhosphorylationSites as ptm
from indra.sources import trips
from indra.explanation import model_checker
from indra.preassembler import Preassembler
import pysb
from copy import deepcopy
from indra.statements import * 
from pick import pick
from modelContext import addRecContext as rc
import pysb
from modelContext import addDephosphorylationAndDegradationRules as add 
from indra.mechlinker import MechLinker
from modelScope import indraDB_query as idb

import logging

logging.getLogger("ac").setLevel(logging.WARNING)
logging.getLogger("indra").setLevel(logging.WARNING)
logging.getLogger("assemble_corpus").setLevel(logging.WARNING)
logging.getLogger("explanation").setLevel(logging.WARNING)
logging.getLogger("paths_graph").setLevel(logging.WARNING)
logging.getLogger("indra/pre_cfpg").setLevel(logging.WARNING)


#TODO:
#check sif direction, if lines are going in both direction for complex, may make things difficult
#make file writing mandatory in sif generation. Decide if want to rm tmp intermediate file

#Build sif file
def buildDirectedSif(stmts,save=True,fn='./directedSifFile.sif'):

    #add rule context changes
    newstmts, uplist, downlist = cs.add_all_af(stmts)     
    newstmts = cs.run_mechlinker_step_reduced(newstmts, uplist, downlist)
    prelimSmallStmts_contextChanges = ptm.coarse_grain_phos(newstmts)
    prelimSmallStmts_contextChanges = Preassembler.combine_duplicate_stmts(prelimSmallStmts_contextChanges)

    #add prespecified drug stmts
    finalSorafStmts = ac.load_statements('finalSorafStmts_reduced.pkl')
    finalStmts = prelimSmallStmts_contextChanges + finalSorafStmts

    #Need to handle receptor active forms and ligand binding/phosphorylations specially, want to better incorporate this into cs/ptm code eventually
    recList = ['PDGFRA','FLT3','KDR']
    newRecPhosStmts = rc.fixRecPhosContext(recList,finalStmts)
    finalStmts = finalStmts + newRecPhosStmts    

    #need corresponding dephosphorylation and degredation for any phosphorylation and transcription if not present in stmt set
    finalStmts = add.add_dephosphorylations(finalStmts)
    finalStmts = add.add_degradations(finalStmts)

    #replace activations if possible, removes unnecessary site and improve context logic
    ml = MechLinker(finalStmts)
    ml.replace_activations()
    finalStmts = ml.statements  
    finalStmts = ac.map_grounding(finalStmts)

    #build sif graph 
    sa = SifAssembler(finalStmts)
    sifModel = sa.make_model(use_name_as_key=True, include_mods=True,include_complexes=True)
    directedSifModel = sa.print_model(include_unsigned_edges=True)
    sa.save_model(fname='./rawSif_tmp.sif',include_unsigned_edges=True)

    newSifModelText = []
    with open('./rawSif_tmp.sif','r') as f:
        sifModelText = f.readlines()
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

    if save:
        with open(fn,'w') as f:
            for line in newSifModelText:
                f.write(line)

    return fn, finalStmts


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
    #Need to handle this better, but for now just add relevant ligands
    #Also need drug and newly added phosphatase in the final model, not sure best place to add that, here for now 
    ligands = ['PDGF','FLT3LG','PDGFA','SORAFENIB','GenericPhosphatase']
    uniqueNodes = uniqueNodes + ligands
    if nodes:
        uniqueNodes = uniqueNodes + nodes
    print('Nodes for final model are %s' % uniqueNodes)
    modelStmts = ac.filter_gene_list(stmts,uniqueNodes,'all')
    modelStmts = ac.map_grounding(modelStmts)

    if save:
        ac.dump_statements(modelStmts,fn)

    modelStmts = Preassembler.combine_duplicate_stmts(modelStmts)
    #FIX THIS BETTER SOMEHOW> FIGURE OUT WHERE ITS COMING FROM
    for st in modelStmts:
        if len(st.agent_list()) > 2:
             st.agent_list()[2:] = []

    pa = PysbAssembler()
    pa.add_statements(modelStmts)
    PySB_Model = pa.make_model()

    return PySB_Model, modelStmts

##########################
#Test experimental findings against candidate model

def testExpStmt(model,sentences): #take in statements? Sentences?

#    #may try to allow stmts or sentences, either through two input variables, or detecting type.
    tp = trips.process_text(sentences)
    testStmts = tp.statements
    testStmts = ac.map_grounding(testStmts)

    mc = model_checker.ModelChecker(model=model,statements=testStmts)
    results = mc.check_model(max_path_length=8)

    return results, mc
#    results[0][1].path_found #True or False [0] index means first tuple, one tuple per sentence checking. [1] is the part of a tuple with failure/success result. 
###########################


#drug targets and cogante ligands are build into funtions. These should be pulled out and made input variables. 
def runModelCheckingRoutine(stmts,drug,nodeToExplain,sentences):
    fn, contextStmts = buildDirectedSif(stmts)
    paths = findPaths(fn,drug,nodeToExplain,8)
    if paths:
        PySB_Model,modelStmts = buildModel(paths,contextStmts)
        results, mc = testExpStmt(PySB_Model,sentences)
        passResult = results[0][1].path_found
    else:
        passResult = None
        modelStmts = None
    return passResult, modelStmts

def findActiveForm(node):
    #find all non-mutation active forms for newly added node 
    stmtsDB_AF = idb.queryDB_nonmod(node,'activeform')
    stmtsDB_AF = remove_mutations(stmtsDB_AF)

    #provide list of af's to follow up on?
    #Ideally, check for a bound condition, and see if binder is in model already. If not, check for phos. If neither...
    #Will rework this whole section in the futre
#    newAFStmts = []
    for st in stmtsDB_AF:
#                    bound_names = []
#                    if st.agent.bound_conditions:
#                        bound_names.append(st.stmts[-1].agent_list_with_bound_condition_agents())
#                    elif st.agent.mods:
        if st.agent.mods:
            if st.agent.mods[0].mod_type == 'phosphorylation':  
                if st.is_active == True:
                    newMod = 'phosphorylation'    
                    #Later, filter and feedback instead of generalizing   
                    st_general = deepcopy(st)
                    st_general.agent.mods[0].res = None    
                    st_general.agent.mods[0].pos = None     
                    st_general.agent.mods[0].residue = None    
                    st_general.agent.mods[0].position = None       
#                            st_general.agent.mods[0].is_active = True
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
        sentence = expObservations[key][1]
        passResult, modelStmts = runModelCheckingRoutine(currentStmts,drug,nodeToExplain,sentence)

        #Model doesn't satisfy condition, go searching for statements to add to model 
        if passResult:
            finalStmts = finalStmts + modelStmts
        else:
            found = 0
            while found == 0:
                #search for potential nodes and stmts to add to the model 
                stmtsDB = idb.queryDB(nodeToExplain,modToExplain)
                stmtsDB = remove_mutations(stmtsDB)
                #identify all new nodes that are candidates for the model 
                candidateNodes = idb.getCandidateNodes(stmtsDB)
                #candidateNodes = [el for el in candidateNodes if el not in currentNodes]
                candidateNodes = list(set(candidateNodes))
                targetNodes = []
                for nd in candidateNodes:
                    if nd in drugTargets:
                        targetNodes.append(nd)

                if len(candidateNodes) >= 20:
                    candidateNodes = candidateNodes[:21] + targetNodes
                #Add exit point
                candidateNodes.append('exit')

                #Test adding each new node, and all of it's associated stmts to the model 
                #See if model passes model checker test with any added node 
                #Could get around slowness of adding lots of nodes by adding one node at a time, then testing.  But then could potentially run the test many, many times 
                for node in candidateNodes:
                    currentNodes.append(node)
                    testStmts = ac.filter_gene_list(currentStmts,currentNodes,'all') + ac.filter_gene_list(stmtsDB,node,'one')  
#                    passResult, modelStmts = runModelCheckingRoutine(testStmts,drug,nodeToExplain,sentence)      
                    passResult, modelStmts = runModelCheckingRoutine(testStmts,drug,finalNode,sentence)      
                    if passResult:
                        found = 1
                        finalStmts = finalStmts + modelStmts
                        print('Worked!')
                        break
                    else:
                        currentNodes.pop()    #Is this necessary?

                #If no single node satisfies condition, pick one node to add to model. 
                if found == 0:
                    #allow user to pick node for next iteration. Following all nodes will result in a combinatorial explosion
                    proteinOptions = candidateNodes
                    title = 'Pick the protein to follow up on, for %s on %s' % (modToExplain,nodeToExplain)
                    #allMechs = allMechs + '\n' + title
                    if proteinOptions:
                        option,index = pick(proteinOptions,title,indicator='=>',default_index=0)    #This really sucks
                    else:
                        option = nodes[1]   #This doesn't seem right 

                    #specify stmts, nodes, mods for next loop iteration
                    #This finds all relevant statements returned by DB for the new node. Will often be one, but should likely limit to one carefully selected statement
                    stmtsForNewNode = ac.filter_gene_list(stmtsDB,option,'one')  
                    newAFStmt = findActiveForm(option)                   
                    currentStmts = currentStmts + stmtsForNewNode + [newAFStmt]
                    currentNodes.append(option)

                    nodeToExplain = option  
                    modToExplain = 'phosphorylation' #This need to be fixed/generalized

                    if option == 'exit':
                        found = 1

    return finalStmts              



############
#TESTING
#load large set of statements, do some messy filtering. 
#Gives a basis to start model building, may not be entirely necessary
#should extend to have the option of starting with some set of statements, or starting from scratch
#Eventually move this to a cardiotox-specific file and just save the filtered statements as the starting point. 
prelimStmts = ac.load_statements('largeModelStmts.pkl') 
stmts = ac.filter_by_type(prelimStmts,Activation)+ac.filter_by_type(prelimStmts,ActiveForm)+ac.filter_by_type(prelimStmts,Phosphorylation) + ac.filter_by_type(prelimStmts,Dephosphorylation) + ac.filter_by_type(prelimStmts,IncreaseAmount)+ac.filter_by_type(prelimStmts,DecreaseAmount) + ac.filter_by_type(prelimStmts,Complex)

large_model_stmts = []
for st in stmts:
    app=1
    for ag in st.agent_list():
        for mod in ag.mods:
            if (mod.mod_type=='sumoylation' or mod.mod_type=='acetylation' or mod.mod_type=='ubiquitination'):
                app = 0 
    if app == 1:
        large_model_stmts.append(st)
        app=0

nodes = ['JUN','STAT1','PKM','RPS6','AURKA','HIF1A','MYC','PDGFRA','KDR','PDGF','FLT3LG','PDGFA','SORAFENIB','GenericPhosphatase']
large_model_stmts = ac.filter_gene_list(large_model_stmts,nodes,'all')


expObservations = {'RPS6':['phosphorylation','SORAFENIB dephosphorylates RPS6'],'PKM':['phosphorylation','SORAFENIB phosphorylates PKM']} #,'AURKA':['phosphorylation','SORAFENIB phosphorylates AURKA'],'HIF1A':['increaseamount','SORAFENIB transcribes HIF1A'],'MYC':['increaseamount','SORAFENIB transcribes MYC'],'JUN':['phosphorylation','SORAFENIB phosphorylates JUN']}#,'STAT1':['phosphorylation','SORAFENIB phosphorylates STAT1']}#,'PDGFRA':'increased'}


finalStmts = expandModel(expObservations,'SORAFENIB',['FLT3','KDR','PDGFRA'],initialStmts=large_model_stmts,initialNodes=nodes)



#key='RPS6'
#paths1 = findPaths(sifFile,'SORAFENIB',key,8)  #checking against original node in question instead
##model, stmts_to_test = buildModel(paths,contextStmts,nodes,save=True,fn='./modelStmts_paths.pkl')
##sentence = effectType[key][1]
##results, mc = testExpStmt(model,sentence)   #PROBLEM: checking partial model, based on paths to current, not final node, against sentence referencing final node. 

#key = 'PKM'
#paths2 = findPaths(sifFile,'SORAFENIB',key,8)  #checking against original node in question instead
##model, stmts_to_test = buildModel(paths,contextStmts,nodes,save=True,fn='./modelStmts_paths.pkl')
##sentence = effectType[key][1]
##results, mc = testExpStmt(model,sentence)   #PROBLEM: checking partial model, based on paths to current, not final node, against sentence referencing final node. 

#key = 'HIF1A'
#paths3 = findPaths(sifFile,'SORAFENIB',key,8)  #checking against original node in question instead
##model, stmts_to_test = buildModel(paths,contextStmts,nodes,save=True,fn='./modelStmts_paths.pkl')
##sentence = effectType[key][1]
##results, mc = testExpStmt(model,sentence)   #PROBLEM: checking partial model, based on paths to current, not final node, against sentence referencing final node. 


#paths = paths1+paths2+paths3
#model, reducedFinalStmts = buildModel(paths,newContextStmts,nodes,save=True,fn='./modelStmts_paths.pkl')


#ac.dump_statements(reducedFinalStmts,'reducedPpotentialFinalGroupMeetingStmts.pkl')

#bngl_model = pysb.export.export(PySB_Model,'bngl')
#with open('final_bngl_model_NEW_NEW_NEW.bngl','w') as f:
#    f.write(bngl_model)





#######################################33
#OLD CODE





#finalModelStmts = []
#finalStmts = []
#for key in expObservations:
#    paths = None
#    results = None
#    passResult = None
#    finalModelStmts = finalModelStmts + finalStmts
#    currentNode = key 
#    if finalModelStmts:
#        firstModel = finalModelStmts
#    else:
#        firstModel = ac.filter_gene_list(large_model_stmts,nodes,'all') 

#    #First check: is there a path between nodes 
#    sifFile,contextStmts = buildDirectedSif(firstModel,save=True,fn='./directedSifFile.sif')
#    paths = findPaths(sifFile,'SORAFENIB',key,8)  
#    
#    #second check, can this path satisfy the mechanism of interest 
#    #10 paths of length 8 for now. This should result in a fairly big/complete model if we're dealing with a small set of base statements
#    if paths:
#        model, stmts_to_test = buildModel(paths,contextStmts,save=True,fn='./modelStmts.pkl')
#        sentence = expObservations[key][1]
#        results, mc = testExpStmt(model,sentence)
#        passResult = results[0][1].path_found
#    else:
#        passResult = False

    



#    effect_tmp = {key:effectType[key] for key in [key]} #this is such a stupid way of doing things, i shouldn't be making/passing entire dict 
#    effect_tmp[key] = effect_tmp[key][0]

#    #Model doesn't satisfy condition, go searching for statements to add to model 
#    if not passResult:
#        found = 0
#        newNodes = []
#        newStmts = []
#        pastNodes = []
#        allMechs = ''

#        nodeToExplain = key
#        modToExplain = effectType[key]

#        while found == 0:



#            stmtsDB = idb.queryDB(nodeToExplain,modToExplain)
#            stmtsDB = remove_mutations(stmtsDB)
#            candidateNodes = idb.getCandidateNodes(stmtsDB)
#            #Here, will have a new set of candidate nodes 

#            #identify all new nodes that are candidates for the model 
#            candidateNodes = [el for el in candidateNodes if el not in pastNodes]
#            candidateNodes = list(set(candidateNodes))


#            #Could get around slowness of adding lots of nodes by adding one node at a time, then testing.  But then could potentially run the test many, many times 
#            modelNodes = nodes+newNodes
#            if len(candidateNodes) >= 20:
##                candidateNodes = candidateNodes[:21]

#            #Can I filter by nodes in model first, and only test with those?
#            #honestly, may not be entirely necessary, but helps automation and enforcing smallness of model by identifying first possible solution
#            for node in candidateNodes:
#                modelNodes.append(node)
#                modelStmts = ac.filter_gene_list(large_model_stmts,modelNodes,'all') + ac.filter_gene_list(stmtsDB,node,'one') + newStmts + finalModelStmts                                                                             
#                sifFile,contextStmts = buildDirectedSif(modelStmts,save=True,fn='./directedSifFile.sif') 
##                ac.dump_statements(contextStmts,'./modelStmts_paths.pkl')

#                paths = findPaths(sifFile,'SORAFENIB',finalNode,8)  #checking against original node in question instead
#                if paths:
#                    model, stmts_to_test = buildModel(paths,contextStmts,nodes,save=True,fn='./modelStmts_paths.pkl')
#                    results, mc = testExpStmt(model,sentence)   #PROBLEM: checking partial model, based on paths to current, not final node, against sentence referencing final node. 
#                    passResult = results[0][1].path_found
#                    if passResult:
#                        found = 1
#                        newNodes.append(node)   #unsure if I need this, can't hurt but may look sloppy to have repeated nodes in final list 
#                        finalStmts = modelStmts
#                        print('Worked!')
#                        break
#                    else:
#                        modelNodes.pop()    #Is this necessary?


#            if found == 0:
#                finalStmts = []
#                #allow user to pick node for next iteration. Following all nodes will result in a combinatorial explosion
#                #maybe try to add more info for the user here, sentence of interest, previous node[s]
#                proteinOptions = candidateNodes
#                title = 'Pick the protein to follow up on, for %s on %s' % (modToExplain,nodeToExplain)
#                allMechs = allMechs + '\n' + title
#                if proteinOptions:
#                    option,index = pick(proteinOptions,allMechs,indicator='=>',default_index=0)    #This really sucks
#                else:
#                    option = nodes[1]   #This doesn't seem right 

#                newStmts = newStmts + ac.filter_gene_list(stmtsDB,option,'one')  
#                stmtsDB_AF = idb.queryDB_nonmod(option,'activeform')
#                stmtsDB_AF = remove_mutations(stmtsDB_AF)
#                #provide list of af's to follow up on?
#                #Ideally, check for a bound condition, and see if binder is in model already. If not, check for phos. If neither...

#                newAFStmts = []
#                for st in stmtsDB_AF:
##                    bound_names = []
##                    if st.agent.bound_conditions:
##                        bound_names.append(st.stmts[-1].agent_list_with_bound_condition_agents())
##                    elif st.agent.mods:
#                    if st.agent.mods:
#                        if st.agent.mods[0].mod_type == 'phosphorylation':  
#                            if st.is_active == True:
#                                newMod = 'phosphorylation'    
#                                #Later, filter and feedback instead of generalizing   
#                                st_general = deepcopy(st)
#                                st_general.agent.mods[0].res = None    
#                                st_general.agent.mods[0].pos = None     
#                                st_general.agent.mods[0].residue = None    
#                                st_general.agent.mods[0].position = None       
#    #                            st_general.agent.mods[0].is_active = True
#                                newAFStmts.append(st_general)       #Make sure to take a stmt where is_active=True
#                                break


#                newStmts = newStmts + newAFStmts
#                newNodes.append(option)
#                nodeToExplain = option 
#                modToExplain = 'phosphorylation' #This need to be fixed/generalized

#    nodes = nodes + newNodes









