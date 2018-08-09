from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from copy import deepcopy
import pickle
from indra.preassembler.hierarchy_manager import hierarchies
from indra.tools import assemble_corpus as ac
from indra.statements import *
from indra.mechlinker import MechLinker
from indra.preassembler import Preassembler
import logging
logging.getLogger("assemble_corpus").setLevel(logging.WARNING)



#TODO: 
#Outstanding question: does this handle multiple receptors and or ligands gracefully?
#Test for infinite loops in while loop 

#def add_receptor_ligand_activeform(stmts):
#    """Identify receptor-ligand pairs and add appropriate active form for receptor.

#    Parameters
#    ----------
#    stmts_in : list[indra.statements.Statement]
#        A list of statements to process.

#    Returns
#    -------
#    stmts_out : list[indra.statements.Statement]
#        A list of statements with added ActiveForm statements for all receptor-ligand pairs.
#    rec_list: list[indra.statements.Agent]
#        A list of receptor agents
#    lig_list: list[indra.statements.Agent]
#        A list of ligand agents
#    """


def get_agent_name(agent):
    name = agent.name
    return name

#Functionalize searching for agent in statements by name
def find_agent(stmts, inputName):
    stmts = ac.strip_agent_context(stmts)
    for st in stmts:
        if None not in st.agent_list():
            for ag in st.agent_list():
                if ag.name == inputName: 
                    output_agent = ag
                    break
    return output_agent


def findAllAgentNames(stmts):
    names = []
    for st in stmts:
        names = names + list(map(lambda obj: obj.name, st.agent_list()))
        names = list(set(names))
    return names 

#    return lig_names, rec_names, lig_list, rec_list

def find_rec_lig(stmts):
    with open('/home/bobby/Dropbox/Sorger_Lab/indra/indra/tools/small_model_tools/lig_rec_dict.pkl','rb') as f:
        receptor_dict = pickle.load(f)
    allRecNames = list(itertools.chain.from_iterable(list(receptor_dict.values())))
    allAgentNames = findAllAgentNames(stmts)
    recInCorpus = [name for name in allAgentNames if name in allRecNames]
    ligRecPairNames = []
    ligRecPairAgents = []
    for receptor in recInCorpus:
        for ligRef, recRef in receptor_dict.items():    
            if receptor in recRef:
                if ligRef in allAgentNames:
                    #Handle names
                    ligRecPair = (ligRef,receptor)  #strings, not agents
                    ligRecPairNames.append(ligRecPair)
                
                    #Handle agents 
                    try:
                        ligAg = find_agent(stmts,ligRef)
                    except UnboundLocalError:
                        print('No agent found with name %s, creating agent' % ligRef)
                        ligAg = Agent(ligRef)
                    try:
                        recAg = find_agent(stmts,receptor)
                    except UnboundLocalError:
                        print('No agent found with name %s' % receptor)
                        recAg = Agent(receptor) #This will be missing db_refs
                    ligRecPairAgent = (ligAg,recAg)  
                    ligRecPairAgents.append(ligRecPairAgent)

                else:
                    backupLigRecPair = (ligRef,receptor) #strings, not agents

        if not ligRecPairNames:
            ligRecPairNames.append(backupLigRecPair)
    ligRecPairNames = list(set(ligRecPairNames))
    return ligRecPairNames, ligRecPairAgents, recInCorpus

#Main idea here: typically, ligand binding is the active form for a receptor
#Exception: if there is phosphorylation of receptor without outside kinase i.e. cis or trans phos following ligand binding 
    #In this case phosphorylation is the active form, but the phosphorylation reaction should require ligand binding. 

#list(map(lambda obj: obj.name, list(itertools.chain.from_iterable(ligRecPairAgentList)))) #Get names form list of tuples of agents

def add_receptor_ligand_activeform(stmts):
    new_af_stmts = []
    ligRecPairs, ligRecPairAgents, recInCorpus = find_rec_lig(stmts)
    rec_lig_stmts = ac.filter_gene_list(stmts,list(itertools.chain.from_iterable(ligRecPairs)),'all')
    
    #Remove old phosphorylation statements involving only ligand and receptor, going to replace with new statements. 
    oldPhosStmts = []
    for st in rec_lig_stmts:
        if isinstance(st,Phosphorylation):
            oldPhosStmts.append(st)
            outputStmts = [st for st in stmts if st not in oldPhosStmts]

    #Need to sort receptors into two groups: those phosphorylated and those not
    allPhosNames = findAllAgentNames(oldPhosStmts)
    phos_rec = []
    non_phos_rec = []
    for rec in recInCorpus:
        if rec in allPhosNames:
            phos_rec.append(rec)
        else:
            non_phos_rec.append(rec)

    #now use rec-lig pairs to build af stmts directly, instead of looping through. handle phos and non phos rec diff
    for pair in ligRecPairAgents:
        if any([el for el in pair if el.name in phos_rec]):

            ligAg = pair[0]     #Fix so that we have agents in pairs 
            recAg = pair[1] 

            kinAg = deepcopy(recAg)
            kinAg.bound_conditions = [BoundCondition(ligAg)]
            newPhosStmt = Phosphorylation(kinAg,recAg)#Is this dangerous/wrong? Should I be checking and modifying phos stmts directly? Minimaly need to remove old ones
            new_af_stmts.append(newPhosStmt)

            afAg = deepcopy(recAg)
            afAg.mods.append(ModCondition(mod_type='phosphorylation'))
            af_stmt = ActiveForm(afAg,activity='activity',is_active=True)
            new_af_stmts.append(af_stmt)

        else: #no phos 
            ligAg = pair[0]     #Fix so that we have agents in pairs 
            recAg = pair[1] 
            recAg.bound_conditions = [BoundCondition(ligAg)]
            af_stmt = ActiveForm(recAg,activity='activity',is_active=True)
            new_af_stmts.append(af_stmt)


    new_af_stmts = Preassembler.combine_duplicate_stmts(new_af_stmts)
    outputStmts = outputStmts + new_af_stmts
    return outputStmts, ligRecPairs, recInCorpus



#identify statements that have receptor but no ligand 
#this allows us to focus on 'next level' 
#i.e. proteins that bind receptor following ligation 
#then build af stmts for these 
#should be able to generalize this beyond rec-lig and cut down a function later

#Work with tuples
#tuple[0] is upstream 
#tuple[1] is downstream 
#find statements with tuple[1], not tuple[0], not *only* tuple[1]
#new tuple 
#tuple[1] becomes tuple[0]
#new tuple[1] comes from found statements 

def findNextPPISet(stmts,oldPairsNames):
    newPairsNames = []
    for pair in oldPairsNames:

        relevantStmts = ac.filter_gene_list(stmts,[pair[1]],'one')
        for st in relevantStmts:
            if pair[0] not in list(map(lambda obj: obj.name, st.agent_list())):
                for ag in st.agent_list():
                    if ag.name != pair[1]:  #this will effectively skip homodimers, which should already be taken care of
                        newPair = (pair[1],ag.name)
                        newPairsNames.append(newPair)
    newPairsNames = list(set(newPairsNames))
    return newPairsNames
#names or agents best?
#Is there a better way for searching and indexing agents than the double loop?
#Maybe just build a separate function that will be cleaner, but same effeciency

def add_active_forms(stmts,newPairsNames):
    new_af_stmts = []
    af_agents = []
    for pair in newPairsNames:
        relevantStmts = ac.filter_gene_list(stmts,list(pair),'all')
        for st in relevantStmts:    
            if isinstance(st, Modification):    
                newStmt = add_modification_active_form(st)
                new_af_stmts.append(newStmt)
                af_agents.append(newStmt.agent.name)
            elif isinstance(st, Complex):
                newStmt = add_complex_active_form(st,  pair[1])
                new_af_stmts.append(newStmt)
                af_agents.append(newStmt.agent.name)
#        elif isinstance(st, IncreaseAmount):
#            new_af_stmts_init = add_transcription_active_form(st, upstream_list)
#            new_af_stmts = new_af_stmts + new_af_stmts_init
     
    new_af_stmts = Preassembler.combine_duplicate_stmts(new_af_stmts) 
    return new_af_stmts, af_agents


def add_modification_active_form(stmt):
    #In this case we know enz is upstream and sub downstream. Shouldn't have to do any matching 
    new_af_stmts = []
    af_agent = deepcopy(stmt.sub)
    af_mods = [stmt._get_mod_condition()] #need to enclose in list 
    af_agent.mods = af_mods
    af_stmt = ActiveForm(af_agent,activity='activity',is_active=True)  #Kinase here is too specific. May be able to take any string
    #new_af_stmts.append(af_stmt) #This leads to a duplicate if there is no new af_stmts (if stmt is not Modification or Compl

    return af_stmt 


def add_complex_active_form(stmt, upstreamElement):
    if all(list(map(lambda obj: obj.name==stmt.agent_list()[0].name, stmt.agent_list()))): #homodimers 
        afAg = deepcopy(stmt.agent_list()[0])
        bindingAg = deepcopy(stmt.agent_list()[0])
    else:
        for ag in stmt.agent_list(): 
            if ag.name == upstreamElement:
                afAg = deepcopy(ag)
            else:
                bindingAg = deepcopy(ag)
    af_boundconditions = [BoundCondition(bindingAg)]
    afAg.bound_conditions = af_boundconditions
    af_stmt = ActiveForm(afAg,activity='activity',is_active=True)

    return af_stmt 





def add_all_af(stmts):
    outputStmts,ligRecPairs,recInCorpus = add_receptor_ligand_activeform(stmts)
    i=0
    afAgents = recInCorpus
    oldPairs = ligRecPairs
    while oldPairs:
        newPairs = findNextPPISet(outputStmts,oldPairs)
        newPairs = [pair for pair in newPairs if pair[1] not in afAgents]
        newAFStmts, newAFAgents = add_active_forms(outputStmts,newPairs)
        outputStmts = outputStmts + newAFStmts
        afAgents = afAgents + newAFAgents
        oldPairs = newPairs

    return outputStmts



def run_mechlinker_step_reduced(stmts):
    ml = MechLinker(stmts)
    ml.gather_explicit_activities() #Why was this commented out?
    ml.gather_modifications()
    ml.require_active_forms() #THIS IS KEY
    ml.require_active_forms_complex() #d.

    updated_stmts = ml.statements
    output_stmts = Preassembler.combine_duplicate_stmts(updated_stmts)   
    return output_stmts






###################
#If an agent has active forms for both bound and modification conditions, only keep mods
def reduce_complex_activeforms(stmts):
    new_af_stmts = []
    af_stmts = ac.filter_by_type(stmts,ActiveForm)
    stmts_to_keep = ac.filter_by_type(stmts,ActiveForm,invert=True)
    af_agents = []
    for st in af_stmts:
        if not any(list(map(lambda obj: obj.entity_matches(st.agent), af_agents))): 
            af_agents.append(st.agent)
    for ag in af_agents:
        ag_stmts = ac.filter_gene_list(stmts,[ag.name],'one')
        ag_af_stmts = ac.filter_by_type(ag_stmts,ActiveForm)
        #all af statements for this agent 
        #Can I replace this whole block with 1-2 one liners?
        if(any([st.agent.mods for st in ag_af_stmts])): #if there are any mod AFs, going to ignore bound_conditions
            for st in ag_af_stmts:
                if st.agent.mods:
                    st.agent.bound_conditions = [] #This may be a terrible idea. Also should be copying
                    stmts_to_keep.append(st)
                elif st.is_active == False:
                    stmts_to_keep.append(st)
        else:
            stmts_to_keep = stmts_to_keep + ag_af_stmts
    output_stmts =  stmts_to_keep
    return output_stmts


#Making multiple phos af's 'and' gated
#Can loop through all proteins (going to be very slow) and pull out all af statements for a given protein

def combine_multiple_phos_activeforms(stmts):
    new_af_stmts = []
    af_stmts = ac.filter_by_type(stmts,ActiveForm)
    output_stmts = ac.filter_by_type(stmts,ActiveForm,invert=True)
    af_agents = []
    for st in af_stmts:
        if not any(list(map(lambda obj: obj.entity_matches(st.agent), af_agents))): 
            af_agents.append(st.agent)
    for ag in af_agents:
        ag_stmts = ac.filter_gene_list(stmts,[ag.name],'one')
        ag_af_stmts = ac.filter_by_type(ag_stmts,ActiveForm)
    
        all_mods = []
        for st in ag_af_stmts:
            if st.agent.bound_conditions:
                output_stmts.append(st)
            elif st.is_active == False:
                output_stmts.append(st)
            else:
                all_mods = all_mods + st.agent.mods
        new_mods = []
        for mod in all_mods:
            if not any(list(map(lambda obj: obj.matches(mod), new_mods))):
                new_mods.append(mod)
        if any([mod.residue for mod in all_mods]):
            #remove any mods with no residue
            new_mods = [mod for mod in all_mods if mod.residue] 

        af_agent = deepcopy(ag)
        af_mods = new_mods
        af_agent.mods = af_mods
        new_afstmt = ActiveForm(af_agent,activity='kinase',is_active=True)
        new_af_stmts.append(new_afstmt)
        output_stmts = new_af_stmts + output_stmts
    return output_stmts

