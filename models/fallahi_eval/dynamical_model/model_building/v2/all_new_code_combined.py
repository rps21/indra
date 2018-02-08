from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
from indra.preassembler.hierarchy_manager import hierarchies
from indra.tools import assemble_corpus as ac
from indra.statements import *
from indra.mechlinker import MechLinker
from copy import deepcopy
from indra.preassembler import Preassembler
from indra.assemblers import PysbAssembler
import pysb
from indra.sources import trips

#TODO: 
#Outstanding question: does this handle multiple receptors and or ligands gracefully?
#Test for infinite loops in while loop 

def add_receptor_ligand_activeform(stmts):
    with open('/home/bobby/Dropbox/Sorger_Lab/BigMech/reporting_and_eval/phase3_eval_pt2/dynamical_model/model_building/lig_rec_dict.pkl','rb') as f:
        receptor_dict = pickle.load(f)

    new_af_stmts = []
    rec_list = []
    lig_list = []
    for st in stmts:
#    if isinstance(st, Complex):
        for ag in st.agent_list(): #checking every ag in every st against dict
            for lig,rec in receptor_dict.items():
                if ag.name in rec:
#                    if not any(list(map(lambda obj: obj.entity_matches(ag), rec_list))):
                    if not any(list(map(lambda obj: obj.matches(ag), rec_list))):
                        rec_ag = ag
                        rec_list.append(ag)#receptor agent 
                        #This is ineffecient but should work. Improve efficiency later. Literally looping through a list within a loop of that list. 
                    for ag in st.agent_list():
                        #This pulls some incorrect entries. e.g. 'EGF' in 'EGFR'. Either direcly match name or add '('
                        #if lig in str(ag):
                        if lig==ag.name:
                            if not any(list(map(lambda obj: obj.entity_matches(ag), lig_list))):
                                lig_list.append(ag)
#                                lig_ag = ag 
                                lig_ag = deepcopy(ag)
                                #Write af stmt here
#                                rec_ag.bound_conditions = rec_ag.bound_conditions + [BoundCondition(lig_ag)]
#                                new_ag = rec_ag
                                new_ag = deepcopy(rec_ag)
                                new_ag.bound_conditions = [BoundCondition(lig_ag)] #Should I make bound_conditions additive? in rec-lig context maybe only need the one 
                                af_stmt = ActiveForm(new_ag,activity='kinase',is_active=True)
                                new_af_stmts.append(af_stmt)
    output_stmts = stmts+new_af_stmts
    return output_stmts, rec_list, lig_list


#######################
###### NEXT STEP ######
### PRIMARY BINDERS ###
#######################
def add_primary_binders_af(stmts,downstream_list,upstream_list): 
    new_downstream_list = []
    downstream_names = []
    upstream_names = []
    for rec in downstream_list:
        downstream_names.append(rec.name)
    for lig in upstream_list:
        upstream_names.append(lig.name)
    rec_stmts = ac.filter_gene_list(stmts,downstream_names,'one')
    lig_stmts = ac.filter_gene_list(rec_stmts,upstream_names,'one')
    rec_no_lig_stmts = [x for x in rec_stmts if x not in lig_stmts]
    new_af_stmts = []
    for st in rec_no_lig_stmts:
        if isinstance(st, Modification):    
            for ag in st.agent_list(): 
                if any([ag.entity_matches(rec) for rec in downstream_list]):   #not sure how this will scale with mulitple receptors
                    #do I need rec_ag for anything?
                    pass
                else:
                    new_downstream_list.append(ag) 
                    af_agent = deepcopy(ag)
                    af_mods = [st._get_mod_condition()] #need to enclose in list 
                    af_agent.mods = af_mods
                    af_stmt = ActiveForm(af_agent,activity='kinase',is_active=True)  #Kinase here is too specific. May be able to take any string
                    new_af_stmts.append(af_stmt) #This leads to a duplicate if there is no new af_stmts (if stmt is not Modification or Complex)
        elif isinstance(st, Complex):
            #TODO: Dimerization works, partly due to ML change. Double check below is necessary and improve efficiencey
            #QUESTION: Did this really do anything? Isn't problem only in applying AF to stmts with two EGFR? Not in creating AF
            #Somehow popped an error on prim_ag not being assigned
            #This shouldn't happen due to complex statement agent list, but should and dimer check, but should be enforced better. 
            prim_ag = None #HORRIBLE TEMP FIX
            for ag in st.agent_list(): 
#                print(upstream_list) #upstream_list is irrelevant, was filtered out. 
                #PROBLEM: Have issue where passing first if check, not second. Some kind of heterodimer? Results in no prim_ag assignment. 
                if any([ag.entity_matches(rec) for rec in downstream_list]):  
                    rec_ag = deepcopy(ag)
                    #extra check for all instead of any to account for dimerization. maybe reorder so this check is first, using if/elif/else instead of nested if
                    if all([ag.entity_matches(rec) for rec in downstream_list]):
                        rec_ag = deepcopy(ag)
                        prim_ag = deepcopy(ag)
                else:
                    prim_ag = deepcopy(ag) 
                    new_downstream_list.append(deepcopy(ag))
                if not prim_ag:
                    prim_ag = deepcopy(ag)

            af_agent = prim_ag
            af_boundconditions = [BoundCondition(rec_ag)]
            af_agent.bound_conditions = af_boundconditions
            af_stmt = ActiveForm(af_agent,activity='kinase',is_active=True)
            new_af_stmts.append(af_stmt) #This leads to a duplicate if there is no new af_stmts (if stmt is not Modification or Complex)
    new_af_stmts = Preassembler.combine_duplicate_stmts(new_af_stmts) #this is a lazy way to avoid adding duplicates, should fix. Now causing a recursion error? Why?

    output_stmts = stmts + new_af_stmts
    new_upstream_list = downstream_list
    return output_stmts, new_downstream_list, new_upstream_list

def add_all_af(stmts):
    upstream_list_total = []
    downstream_list_total = []

    #First, handle receptors
    updated_stmts, rec_list, lig_list = add_receptor_ligand_activeform(stmts)
    #updated_stmts = add_receptor_ligand_context(updated_stmts,rec_list)
    upstream_list = lig_list
    downstream_list = rec_list
    
    #CHECK FOR INFINITE LOOPS
    i=0
    while downstream_list:
        updated_stmts, new_downstream_list, new_upstream_list = add_primary_binders_af(updated_stmts, downstream_list, upstream_list)
        
        upstream_list_total.append(upstream_list)
        downstream_list_total.append(downstream_list)

##        #Mech linker step
#        ml = MechLinker(updated_stmts)
#        ml.gather_explicit_activities() #Why was this commented out?
#        ml.gather_modifications()
#        ml.require_active_forms() #THIS IS KEY
#        ml.require_active_forms_complex(downstream_list, upstream_list) #reexamine why input list of agents is needed.
#        updated_stmts = ml.statements
#        updated_stmts = Preassembler.combine_duplicate_stmts(updated_stmts)      

        #Reset for next iteration
        downstream_list = new_downstream_list
        upstream_list = new_upstream_list   

        print(i)
        if i >= 30:
            break
        i=i+1

    updated_stmts = Preassembler.combine_duplicate_stmts(updated_stmts)

#    #NEW
#    updated_stmts = reduce_complex_activeforms(updated_stmts)
#    updated_stmts = combine_multiple_phos_activeforms(updated_stmts)
  
    return updated_stmts, downstream_list_total, upstream_list_total


def run_mechlinker_step(stmts,downstream_list,upstream_list):
    ml = MechLinker(stmts)
    for i in range(len(downstream_list)):
#        ml = MechLinker(stmts)
        ml.gather_explicit_activities() #Why was this commented out?
        ml.gather_modifications()
        ml.require_active_forms() #THIS IS KEY
        ml.require_active_forms_complex(downstream_list[i], upstream_list[i]) #reexamine why input list of agents is needed.
        print(downstream_list[i])
        print(upstream_list[i])
        updated_stmts = ml.statements
        updated_stmts = Preassembler.combine_duplicate_stmts(updated_stmts)      
        ml = MechLinker(updated_stmts)
    output_stmts = ml.statements
    return output_stmts


#####################
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
    stmts_to_keep = ac.filter_by_type(stmts,ActiveForm,invert=True)
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
                stmts_to_keep.append(st)
            elif st.is_active == False:
                stmts_to_keep.append(st)
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
        output_stmts = new_af_stmts + stmts_to_keep
    return output_stmts



#dimerization testing


##More testing
#test_stmts = ac.load_statements('../fallahi_eval_pysb_stmts_updated.pkl')
#erk_model_st = ac.filter_gene_list(test_stmts,['EGF','EGFR','SOS1','GRB2','KRAS','BRAF','MAPK1','MAPK3','MAP2K1'],'all')
#erk_model_st_reduced = ac.filter_by_type(erk_model_st,Translocation,invert=True)

#old_af = ac.filter_by_type(erk_model_st_reduced,ActiveForm)
#new_af = ac.filter_by_type(test_output,ActiveForm)

#pa = PysbAssembler('two_step')
#pa.add_statements(test_output)
#model2 = pa.make_model()

#bngl_model_filter = pysb.export.export(model2,'bngl')
#bngl_file_filter = open('rbm/new_version.bngl','w')
#bngl_file_filter.write(bngl_model_filter)
#bngl_file_filter.close()



#Last up in the air case (before more testing):
#Case 1:
# ActiveForm(MAP2K1(mods: (phosphorylation)), kinase, True),
# ActiveForm(MAP2K1(bound: [MTOR, True]), kinase, True),

#Case 1: One phosphorylation active form. One bound active form
#Seems likely the both are not correct
#Either a Branching (likely with error) or just an error 
#If it's a branching, 'or' is correct - works that way now
#If it's an error should revert to phos only - maybe check test cases









##############################################################


#Function added to mechlinker to handle complex in ActiveForm
#    def require_active_forms_receptor_complex(self,rec_list):

#        new_stmts = []
#        for stmt in self.statements:
#            if isinstance(stmt, Complex):
#                #only this if block needs to change
#                #this may be the place to check if either agent is in 'input list'. here, receptors
##                if ag in rec_list for ag in stmt.agent_list(): #set.intersection or list comprehension are reasonable, but probably need to loop since should use ag.entity_matches
#                match=0
#                for ag in stmt.agent_list():
#                    for rec in rec_list:
#                        if ag.entity_matches(rec):
#                            match = 1
#                            rec_ag = rec
#                            rec_ag_index = stmt.agent_list().index(rec_ag)

#                if match == 0: 
#                    new_stmts.append(stmt)
#                    continue


#                rec_base = self._get_base(rec_ag) 
#                active_forms = rec_base.get_active_forms()
#                if not active_forms:
#                    new_stmts.append(stmts)
#                else:
#                    for af in active_forms:
#                        new_stmt = deepcopy(stmt)
#                        new_stmt.uuid = str(uuid.uuid4())
##                        af.apply_to(new_stmt.enz)  #Need clean way to refer to appropriate agent in stmt. 
#                        af.apply_to(new_stmt.agent_list()[rec_ag_index])  #Need clean way to refer to appropriate agent in stmt. 
#                        new_stmts.append(new_stmt)    

#            #rest of this is fine
#            else:
#                new_stmts.append(stmt)
#        self.statements = new_stmts
#        return new_stmts


######################################################################################

from indra.tools import assemble_corpus as ac
from indra.statements import *      

#Things to think about.
#Right now, combining phos sites on substrate, if same enz and substrate
#But combine sites on kinase no matter what
#Ideally should make a stricter version on sub, looser version on kin, and provide option


#Temporary hack to ignore mutations for now
#af_stmts.pop(1)

#rewrites active form statements to coarse grain pY and pS/T sites 
def coarse_grain_phos_on_af(stmts):
    af_stmts = ac.filter_by_type(stmts, ActiveForm)
    non_af_stmts = ac.filter_by_type(stmts, ActiveForm, invert=True)
    new_af_stmts = []
    for st in af_stmts:   
        act_agent = st.agent
        #Preserve inhibitory statements
        if st.is_active == True: 
            #in each AF statement, identify Tyrosine and Serine/Threonine phosphorylations
            phos_mods = [mod for mod in st.agent.mods if mod.mod_type == 'phosphorylation']
            pY_mods = [mod for mod in phos_mods if mod.residue == 'Y']
            pST_mods = [mod for mod in phos_mods if mod.residue in ('S','T')] 
            modlist = []
            try:
                generic_pY = ModCondition(pY_mods[0].mod_type,pY_mods[0].residue)
                modlist.append(generic_pY)
            except IndexError:
                generic_pY = []
            try:
                generic_pST = ModCondition(pST_mods[0].mod_type,'S') #This could appear as S or T but I'm assuming equivaqlency. scratch that, make s, to remove redundancies
                modlist.append(generic_pST)
            except IndexError:
                generic_pST = []
            #Write replacement ActiveForm statement with coarse-grained modifications
            newagent = Agent(act_agent.name,modlist)
            newst = ActiveForm(newagent,st.activity,st.is_active)
            new_af_stmts.append(newst)
        else:
            new_af_stmts.append(st)
    new_af_stmts = ac.Preassembler.combine_duplicate_stmts(new_af_stmts)
    filtered_stmts = non_af_stmts + new_af_stmts
    return(filtered_stmts)



##This may be unnecessary if I'm now further coarse graining AF. Could be useful as a middle road
#        #loop through list of mods on activeform statement
#        for mod in act_agent.mods:
#            #loop through all phosphorylation statements
#            for pst in phos_stmts:
#                #check if mod in af statement matches a phosphorylation statement
#                if act_agent.entity_matches(pst.agent_list()[1]): #Second agent is one being modified
#                    if mod.position == pst.position:
#                        used_stmts.append(pst)

        #Generic, coarse-grained
#        try:
#            newmod = ModCondition(mod_type=act_agent.mods[0].mod_type,residue=act_agent.mods[0].residue) #needs improvement
#        except IndexError:
#            #this is mut instead of mod, maybe other exceptions I haven't come across yet
#            newmod = None
#            #pass
#        try:
#            newmut = ModCondition(mod_type=act_agent.mods[0].mod_type,residue=act_agent.mods[0].residue) #needs improvement
#        except IndexError:
#            #this is mut instead of mod, maybe other exceptions I haven't come across yet
#            newmod = None#[]
#need to handle 'activating' mutations differently
#ActiveForm(HRAS(muts: (G, 12, V)), gtpbound, True) #Mutation condition may break things
#In [82]: af_stmts[1].agent.mutations
#Out[82]: [MutCondition(G, 12, V)]



#Sort through phosphorylation stmts and remove ones with the same enzyme and substrate   
def remove_redundant_phosphorylations(stmts):
    """Filter to statements containing genes only.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    specific_only : Optional[bool]
        If True, only elementary genes/proteins will be kept and families
        will be filtered out. If False, families are also included in the
        output. Default: False
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of coarse grained statements.
    """
    phos_stmts = ac.filter_by_type(stmts, Phosphorylation) + ac.filter_by_type(stmts, Dephosphorylation)
    non_phos_stmts = ac.filter_by_type(stmts, Phosphorylation, invert = True)
    non_phos_stmts = ac.filter_by_type(non_phos_stmts, Dephosphorylation, invert = True)
    #Start new list with coarse grained first statement
    unique_phos_stmts = [phos_stmts[0]]
    unique_phos_stmts[0].position=None
    for st1 in phos_stmts:
        unique = 1
        for st2 in unique_phos_stmts:
            #Comparae enz/sub pair for list of phos stmts with stmts in new list. If enz/sub match already added statement, discard. Otherwise add and coarse-grain (remove phos position)
            if st1.enz.entity_matches(st2.enz) and st1.sub.entity_matches(st2.sub):    
                unique = 0
#            else:
        if unique == 1:
            unique_phos_stmts.append(st1)
            unique_phos_stmts[-1].position = None
    filtered_stmts = non_phos_stmts + unique_phos_stmts
    return filtered_stmts

#Coarse graining kinase context
##rewrites active form statements to coarse grain pY and pS/T sites 
def coarse_grain_kinase_context(stmts):
    phos_stmts = ac.filter_by_type(stmts, Phosphorylation) + ac.filter_by_type(stmts, Dephosphorylation)
    non_phos_stmts = ac.filter_by_type(stmts, Phosphorylation, invert = True)
    non_phos_stmts = ac.filter_by_type(non_phos_stmts, Dephosphorylation, invert = True)
    kinase_with_mods = list(filter(lambda x: x.enz.mods != [], phos_stmts)) #had to add list() for pyth3 combatability
    kinase_wo_mods = list(filter(lambda x: x.enz.mods == [], phos_stmts))
    
    coarse_grained_phos_stmts = []
    for st in kinase_with_mods:
        pY_mods = [mod for mod in st.enz.mods if mod.residue == 'Y']
        pST_mods = [mod for mod in st.enz.mods if mod.residue in ('S','T')]  
        modlist = []
        try:
            generic_pY = ModCondition(pY_mods[0].mod_type,pY_mods[0].residue)
            modlist.append(generic_pY)
        except IndexError:
            generic_pY = []
        try:
            generic_pST = ModCondition(pST_mods[0].mod_type,'S') #This could appear as S or T but I'm assuming equivaqlency. Again, scratch this for S only
            modlist.append(generic_pST)
        except IndexError:
            generic_pST = []
        #Write replacement ActiveForm statement with coarse-grained modifications
        newagent = Agent(st.enz.name,modlist)
        newst = Phosphorylation(newagent,st.sub,st.residue)
        coarse_grained_phos_stmts.append(newst)

    new_phos_stmts = coarse_grained_phos_stmts+kinase_wo_mods
    new_phos_stmts = ac.Preassembler.combine_duplicate_stmts(new_phos_stmts)
    filtered_stmts = non_phos_stmts + new_phos_stmts
    return filtered_stmts

def coarse_grain_substrate_context(stmts):
    phos_stmts = ac.filter_by_type(stmts, Phosphorylation) + ac.filter_by_type(stmts, Dephosphorylation)
    non_phos_stmts = ac.filter_by_type(stmts, Phosphorylation, invert = True)
    non_phos_stmts = ac.filter_by_type(non_phos_stmts, Dephosphorylation, invert = True)
    substrate_with_mods = list(filter(lambda x: x.sub.mods != [], phos_stmts)) #had to add list() for pyth3 combatability
    substrate_wo_mods = list(filter(lambda x: x.sub.mods == [], phos_stmts))
    
    coarse_grained_phos_stmts = []
    for st in substrate_with_mods:
        pY_mods = [mod for mod in st.sub.mods if mod.residue == 'Y']
        pST_mods = [mod for mod in st.sub.mods if mod.residue in ('S','T')]  
        modlist = []
        try:
            generic_pY = ModCondition(pY_mods[0].mod_type,pY_mods[0].residue)
            modlist.append(generic_pY)
        except IndexError:
            generic_pY = []
        try:
            generic_pST = ModCondition(pST_mods[0].mod_type,'S') #This could appear as S or T but I'm assuming equivaqlency. Again, scratch this for S only
            modlist.append(generic_pST)
        except IndexError:
            generic_pST = []
        #Write replacement ActiveForm statement with coarse-grained modifications
        newagent = Agent(st.sub.name,modlist)
        newst = Phosphorylation(st.enz,newagent,st.residue)
        coarse_grained_phos_stmts.append(newst)

    new_phos_stmts = coarse_grained_phos_stmts+substrate_wo_mods
    new_phos_stmts = ac.Preassembler.combine_duplicate_stmts(new_phos_stmts)
    filtered_stmts = non_phos_stmts + new_phos_stmts
    return filtered_stmts

#stmts = my_model_stmts_larger[:]
#stmts = coarse_grain_phos_on_af(stmts)
#stmts = remove_redundant_phosphorylations(stmts)
#stmts = coarse_grain_kinase_context(stmts)

#pa = PysbAssembler('two_step')
#pa.add_statements(my_model_stmts_larger)
#model = pa.make_model()

#bngl_model_filter = pysb.export.export(model,'bngl')
#bngl_file_filter = open('rbm/egfr_small_no_cg.bngl','w')
#bngl_file_filter.write(bngl_model_filter)
#bngl_file_filter.close()


#pa = PysbAssembler('two_step')
#pa.add_statements(stmts)
#model = pa.make_model()

#bngl_model_filter = pysb.export.export(model,'bngl')
#bngl_file_filter = open('rbm/group_meeting_6_5_17/small.bngl','w')
#bngl_file_filter.write(bngl_model_filter)
#bngl_file_filter.close()





#ac.dump_statements('egf_egfr_sos_grb_olig_stmts.pkl',stmts)

######################################################################################################################

def remove_dimers(stmts_in):
    stmts_new = stmts_in[:]
    for st in stmts_in:
        if len(st.agent_list()) > 1:
            if st.agent_list()[0]: #Weird circumstance with a None type in transcription
                if st.agent_list()[0].entity_matches(st.agent_list()[1]):
                    stmts_new.remove(st)
    return stmts_new

def remove_mutations(stmts_in):
    stmts_new = stmts_in[:]
    for st in stmts_in:
        for ag in st.agent_list():
                if ag.mutations:
                    stmts_new.remove(st)
    return stmts_new

