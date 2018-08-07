from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from copy import deepcopy
import pickle
from indra.preassembler.hierarchy_manager import hierarchies
from indra.tools import assemble_corpus as ac
from indra.statements import *
from indra.mechlinker import MechLinker
from indra.preassembler import Preassembler



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
def find_agent(stmts, name_list):
    output_agents = []
    output_names = []
    stmts = ac.strip_agent_context(stmts)
    for st in stmts:
        if None not in st.agent_list():
            ag_list = [ag for ag in st.agent_list() if ag.name in name_list and not any(list(map(lambda obj: obj.matches(ag), output_agents)))]
            output_agents = output_agents + ag_list #receptor agent 
    return output_agents


#HERE################
#Find pairs not two separate lists
def find_rec_lig(stmts):
    #Requires dict of receptor-ligand pairs, should store this somewhere better
    with open('/home/bobby/Dropbox/Sorger_Lab/indra/indra/tools/small_model_tools/lig_rec_dict.pkl','rb') as f:
        receptor_dict = pickle.load(f)

    ligRecPairList = []
    for sublist in list(receptor_dict.values()):
        rec_list_init = find_agent(stmts,sublist)   #have a list of receptors to loop through. For each receptor, want either all pairs in stmts or first addition from dict if none 
        if rec_list_init:
            for rec in rec_list_init:
                #this pulls first ligand from dict. Should probably allow searching a corpus for ligands already present in st list 
                #break this down rec by rec (may only be one) and combine at end?
                for lig, recs in receptor_dict.items():    
                    if rec.name in recs:
                        #Testing
                        lig_ag = find_agent(stmts,[lig])
                        if lig_ag:
                            ligRecPair = (lig_ag[0].name,rec.name)  #strings, not agents
                            ligRecPairList.append(ligRecPair)
                            break
                        else:
                            backupLigRecPair = (lig,rec.name)  #strings, not agents

                if not ligRecPairList:
                    ligRecPairList.append(backupLigRecPair)
    ligRecPairList = list(set(ligRecPairList))
    return ligRecPairList

#    lig_list, lig_names = find_agent(stmts,list(receptor_dict.keys()))

#    return lig_names, rec_names, lig_list, rec_list







def add_receptor_ligand_activeform(stmts):
    lig_ag = []
    new_af_stmts = []
    lig_names, rec_names, lig_list, rec_list = find_rec_lig(stmts)
    rec_lig_stmts = ac.filter_gene_list(stmts,rec_names+lig_names,'all')
    
#    for st in rec_lig_stmts:
#        if isinstance(st, Modification):   

    for st in rec_lig_stmts:
        for ag in st.agent_list():
            #This step could probably be combined with find_rec_lig(). Should find a way to find pairs of lig-rec instead of individual lists
            if ag.name in lig_names:
                lig_ag = deepcopy(ag) 
            if ag.name in rec_names:
                rec_ag = deepcopy(ag) 
        if lig_ag:
            rec_ag.bound_conditions = [BoundCondition(lig_ag)]
            af_stmt = ActiveForm(rec_ag,activity='kinase',is_active=True)
            new_af_stmts.append(af_stmt)
            lig_ag = None
    new_af_stmts = Preassembler.combine_duplicate_stmts(new_af_stmts)
    output_stmts = stmts + new_af_stmts
    return output_stmts, lig_list, rec_list

#def fixRecPhosContext(rec_names, stmts):
#    phosStmts = ac.filter_by_type(stmts,Phosphorylation)
#    newRecPhosStmts = []
#    for rec in rec_names:
#        recPhosStmts = ac.filter_gene_list(phosStmts,[rec],'one')
#        recPhosStmts_self = []
#        for st in recPhosStmts:
#            if st.enz.name == rec and st.sub.name == rec:
#                recPhosStmts_self.append(st)

#        for st in recPhosStmts_self:
#            newSt = deepcopy(st)
#            newSt.enz.mods = []
#            lig_ag = getLigAgent(rec,stmts)
#            if lig_ag:
#                newSt.enz.bound_conditions = [BoundCondition(lig_ag)]
#                newRecPhosStmts.append(newSt)
#    return newRecPhosStmts



#identify statements that have receptor but no ligand 
#this allows us to focus on 'next level' 
#i.e. proteins that bind receptor following ligation 
#then build af stmts for these 
#should be able to generalize this beyond rec-lig and cut down a function later

def find_next_downstream_list(next_level_stmts,old_upstream_list,old_downstream_list):
    new_downstream_list = []
    new_upstream_list = []
    for st in next_level_stmts:
        for ag1 in st.agent_list():
            for ag2 in old_downstream_list:
                if not ag1.entity_matches(ag2):
                    new_downstream_list.append(ag1)
                else:
                    new_upstream_list.append(ag1)
    new_downstream_list = list(set(new_downstream_list))
    new_upstream_list = list(set(new_upstream_list))
    return new_downstream_list, new_upstream_list

def find_next_cascade_step(stmts,upstream_list,downstream_list): 
    new_downstream_list = []
    downstream_names = []
    upstream_names = []
    for rec in downstream_list:
        downstream_names.append(rec.name)
    for lig in upstream_list:
        upstream_names.append(lig.name)
    downstream_names = list(set(downstream_names))
    upstream_names = list(set(upstream_names))
    downstream_stmts = ac.filter_gene_list(stmts,downstream_names,'one')
    upstream_stmts = ac.filter_gene_list(stmts,upstream_names,'one')
    next_level_stmts = [x for x in downstream_stmts if x not in upstream_stmts] #identifies stmts involing elements in the downstream list but not upstream, i.e. receptors binding adaptors but not ligands
    new_downstream_list, new_upstream_list = find_next_downstream_list(next_level_stmts,upstream_list,downstream_list)
    return next_level_stmts, new_downstream_list, new_upstream_list



def add_modification_active_form(stmt, upstream_list):
    new_af_stmts = []
    for ag in stmt.agent_list(): 
        if any([ag.entity_matches(previous_level) for previous_level in upstream_list]):   
            pass            
        else:
            af_agent = deepcopy(ag)
            af_mods = [stmt._get_mod_condition()] #need to enclose in list 
            af_agent.mods = af_mods
            af_stmt = ActiveForm(af_agent,activity='kinase',is_active=True)  #Kinase here is too specific. May be able to take any string
            new_af_stmts.append(af_stmt) #This leads to a duplicate if there is no new af_stmts (if stmt is not Modification or Complex)
    return new_af_stmts #, new_downstream_list

def add_complex_active_form(stmt, upstream_list):
    af_agent = None
    new_af_stmts = []
    for ag in stmt.agent_list(): 
        if all([ag.entity_matches(previous_level) for previous_level in upstream_list]): #check for homodimers
            previous_agent = deepcopy(ag)
            af_agent = deepcopy(ag)
        elif any([ag.entity_matches(previous_level) for previous_level in upstream_list]):  
            previous_agent = deepcopy(ag)
        else:
            af_agent = deepcopy(ag) 
    af_boundconditions = [BoundCondition(previous_agent)]
    af_agent.bound_conditions = af_boundconditions
    af_stmt = ActiveForm(af_agent,activity='kinase',is_active=True)
    new_af_stmts.append(af_stmt) #This leads to a duplicate if there is no new af_stmts (if stmt is not Modification or Complex)
    return new_af_stmts 

#def add_transcription_active_form(stmt, upstream_list):
#    new_af_stmts = []
#    ag = stmt.subj
#    if any([ag.entity_matches(previous_level) for previous_level in upstream_list]):   
#        pass            
#    else:
#        af_agent = deepcopy(ag)
#        af_stmt = ActiveForm(af_agent,activity='transcription',is_active=True)  #Kinase here is too specific. May be able to take any string
#        new_af_stmts.append(af_stmt) #This leads to a duplicate if there is no new af_stmts (if stmt is not Modification or Complex)
#    return new_af_stmts
#Not sure this is working correctly. Requires knowledge tf is modified, rather than finding modification of tf. 

#Translocation
#filter ones without a from or to location 
#In [9]: trans1.from_location
#Then look for upstream partner
#def add_complex_active_form(stmt, upstream_list):
#    af_agent = None
#    new_af_stmts = []
#    for ag in stmt.agent_list(): 
#        if all([ag.entity_matches(previous_level) for previous_level in upstream_list]): #check for homodimers
#            previous_agent = deepcopy(ag)
#            af_agent = deepcopy(ag)
#        elif any([ag.entity_matches(previous_level) for previous_level in upstream_list]):  
#            previous_agent = deepcopy(ag)
#        else:
#            af_agent = deepcopy(ag) 
#    af_boundconditions = [BoundCondition(previous_agent)]
#    af_agent.bound_conditions = af_boundconditions
#    af_stmt = ActiveForm(af_agent,activity='kinase',is_active=True)
#    new_af_stmts.append(af_stmt) #This leads to a duplicate if there is no new af_stmts (if stmt is not Modification or Complex)
#    return new_af_stmts 




def add_active_forms(next_level_stmts,upstream_list):
    #split to handle modification statements and complex statements separately, notably missing transcription
    #break these into separate functions for clarity
    new_af_stmts = []
    for st in next_level_stmts:
        #is this method going to lead to a lot of duplicates?
        if isinstance(st, Modification):    
            new_af_stmts_init = add_modification_active_form(st, upstream_list)
            new_af_stmts = new_af_stmts + new_af_stmts_init
        elif isinstance(st, Complex):
            new_af_stmts_init = add_complex_active_form(st, upstream_list)
            new_af_stmts = new_af_stmts + new_af_stmts_init
        elif isinstance(st, IncreaseAmount):
            new_af_stmts_init = add_transcription_active_form(st, upstream_list)
            new_af_stmts = new_af_stmts + new_af_stmts_init
     
    new_af_stmts = Preassembler.combine_duplicate_stmts(new_af_stmts) 
    return new_af_stmts

####################################



def add_all_af(stmts):
    upstream_list_total = []
    downstream_list_total = []

    #First, handle receptors
    updated_stmts, lig_list, rec_list = add_receptor_ligand_activeform(stmts)

    upstream_list = lig_list
    downstream_list = rec_list
    
    #CHECK FOR INFINITE LOOPS
    i=0
    outputStmts = []
    while downstream_list:
        next_stmts, new_downstream_list, new_upstream_list = find_next_cascade_step(updated_stmts,upstream_list, downstream_list)
        new_af_stmts = add_active_forms(next_stmts,new_upstream_list)
        updated_stmts = updated_stmts + new_af_stmts

        #fixing loop issue
        newStmts = [st for st in updated_stmts if st not in outputStmts]
        outputStmts = newStmts + outputStmts
        updated_stmts = [st for st in updated_stmts if st not in next_stmts]
        #end fix 

        #Save lists 
        upstream_list_total.append(upstream_list)
        downstream_list_total.append(downstream_list)

        #Reset for next iteration
        downstream_list = new_downstream_list
        upstream_list = new_upstream_list   

    outputStmts = reduce_complex_activeforms(outputStmts) 
    outputStmts = combine_multiple_phos_activeforms(outputStmts)

    return outputStmts, downstream_list_total, upstream_list_total


def run_mechlinker_step_reduced(stmts,downstream_list,upstream_list):
    ml = MechLinker(stmts)
    ml.gather_explicit_activities() #Why was this commented out?
    ml.gather_modifications()
    ml.require_active_forms() #THIS IS KEY

    print(len(downstream_list))
    for i in range(len(downstream_list)):
        ml.require_active_forms_complex(downstream_list[i], upstream_list[i]) #reexamine why input list of agents is needed.
        updated_stmts = ml.statements
        updated_stmts = Preassembler.combine_duplicate_stmts(updated_stmts)      
        ml = MechLinker(updated_stmts)
        ml.gather_explicit_activities() #Why was this commented out?
        ml.gather_modifications()
        ml.require_active_forms() #THIS IS KEY

    output_stmts = ml.statements
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

