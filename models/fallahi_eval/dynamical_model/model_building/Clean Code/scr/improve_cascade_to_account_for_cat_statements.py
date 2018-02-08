from indra import mechlinker as ml


#Rec-Ligand
#EGF
#EGFR

#Primary binders
#{AKT1(bound: [EGFR, True]): EGFR(bound: [EGF, True]),   #This is from a phosphorylation statement 
# GRB2(bound: [EGFR, True]): EGFR(bound: [EGF, True]),
# EGFR(bound: [EGF, True], bound: [EGFR, True]): EGFR(bound: [EGF, True])}


#Right now find ligand-receptor
#Then find all statements with receptor but now ligand 
#Then add context to theses statements such that the receptor must be bound to the ligand for additional statments to fire.
#This should not change. This allows phosphorylation of akt by egfr to require egf binding. with two-step assembly this will result in free p-akt

#Next, take all statements with receptor and no ligand. These have added context that the receptor is bound to ligand
#Look for all statements that have primary binders as definied based on previous statemnt list 
#To each primary binder, add context that they be bound to receptor to be functional.
#The problem here comes from catalytic statements. If a receptor catalyzes a phosphorylation or other modifications on a protein, they need that modification to be active, not to be bound to receptor. 

#Best idea might be two parallel tracks of similar but distinct functions, which take in different sets of statements.
#Complex goes to set1. 'Else' goes to set2. filter out ActiveForm. Double check there aren't other types to worry about. Ans: Translocation should also be filtered.  
#Split has to happen at first step, before receptor-nonlig binding. End result at this stage should be the same, but will need to inherit info from statement type.
#Functions inheriting Complex statements should behave exactly the same
#Functions inheriting catalytic statements:
#    Receptor-primary statement should function the same. Primary catalysis only happens when ligand is bound to receptor.
#    Next set of primary statements are different. Functionally similar process, but the context we add at the end is different.
#    Instead of inheriting a bound condition, inherit mods list.



from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
from indra.preassembler.hierarchy_manager import hierarchies
from indra.tools import assemble_corpus as ac
from indra.statements import *

#Example function calls:
##Theses functions all need to work together in this order. In future should combine in a better way
##New Version:
stmts_in = ac.load_statements('test_stmts_for_context.pkl')
rec_lig_list = find_rec_lig(stmts_in)
stmts = gather_rec_stmts(stmts_in,rec_lig_list)
stmts = add_ligand_context(stmts,rec_lig_list) #stmts here are list of three lists

#problem is supposed to call this with rec_nolig_st
#now need three sets of statements: complex and phos, both of which are fed separtely through the next two sets, and others which are just maintained
#Which base do I work off of? stmts variable is messy - three lists. 
#I guess stmts[1]
complex_stmts = ac.filter_by_type(stmts[1],Complex)
phos_stmts = ac.filter_by_type(stmts[1],Phosphorylation) + ac.filter_by_type(stmts[1],Autophosphorylation) 
other_stmts = [x for x in stmts[1] if x not in phos_stmts and x not in complex_stmts]

#Complex is failing to include egfr with grb2 key
primary_binders_dict_comp = find_primary_binders_complex(complex_stmts,rec_lig_list) #This is getting sloppy, but middle list in stmts is rec_nolig_st 
primary_binders_dict_phos = find_primary_binders_phos(phos_stmts,rec_lig_list) #This is getting sloppy, but middle list in stmts is rec_nolig_st 


##Still working on these
#my_primary_stmts, my_non_primary_stmts = find_primary_stmts(primary_binders_dict,rec_lig_list,stmts_in)
#final_primary_st = add_context_to_stmts(my_primary_stmts)
#final_stmts = final_primary_st + my_non_primary_stmts + stmts[0] + stmts[2]



#Old Version:
#stmts_in = ac.load_statements('test_stmts_for_context.pkl')
#rec_lig_list = find_rec_lig(stmts_in)
#stmts = gather_rec_stmts(stmts_in,rec_lig_list)
#stmts = add_ligand_context(stmts,rec_lig_list) #stmts here are list of three lists

##problem is supposed to call this with rec_nolig_st
#primary_binders_dict = find_primary_binders(stmts[1],rec_lig_list) #This is getting sloppy, but middle list in stmts is rec_nolig_st 
#my_primary_stmts, my_non_primary_stmts = find_primary_stmts(primary_binders_dict,rec_lig_list,stmts_in)
#final_primary_st = add_context_to_stmts(my_primary_stmts)
#final_stmts = final_primary_st + my_non_primary_stmts + stmts[0] + stmts[2]



def split_stmts_for_context(stmts1):
    stmts = ac.filter_by_type(stmts1,ActiveForm,invert=True)
    stmts = ac.filter_by_type(stmts,Translocation,invert=True)
    comp_stmts = ac.filter_by_type(stmts,Complex)
    cat_stmts = ac.filter_by_type(stmts,Complex,invert=True)
    return comp_stmts, cat_stmts




#Function takes a list of statments as input and outputs a list of two lists, of receptor agents and ligand agents respectively
#This is horribly inefficient and needs to be cleaned up, but at least it works again 
def find_rec_lig(stmts):
    with open('/home/bobby/Dropbox/Sorger_Lab/BigMech/reporting_and_eval/phase3_eval_pt2/dynamical_model/model_building/lig_rec_dict.pkl','rb') as f:
        receptor_dict = pickle.load(f)

    rec_list = []
    lig_list = []
    for st in stmts:
#    if isinstance(st, Complex):
        for ag in st.agent_list(): #checking every ag in every st against dict
            for lig,rec in receptor_dict.items():
                if ag.name in rec:
                    if not any (list(map(lambda obj: obj.entity_matches(ag), rec_list))):
                        rec_list.append(ag)#receptor agent 
                        #This is ineffecient but should work. Improve efficiency later
                    for ag in st.agent_list():
                        #This pulls some incorrect entries. e.g. 'EGF' in 'EGFR'. Either direcly match name or add '('
                        #if lig in str(ag):
                        if lig==ag.name:
                            if not any (list(map(lambda obj: obj.entity_matches(ag), lig_list))):
                                lig_list.append(ag)
    #For now, remove any mods on receptors, just take base agent. may revisit this.
    for rec in rec_list:
        rec.mods=[]
    rec_lig_list =[]
    rec_lig_list.append(rec_list)
    rec_lig_list.append(lig_list)
    #rec_lig_list is now a list of two lists, first is receptor agents, second is ligand agents
    return rec_lig_list


#New version
def gather_rec_stmts(model_stmts,rec_lig_list):
    rec_stmts = []
    non_rec_stmts = []
    rec_list = rec_lig_list[0]
    lig_list = rec_lig_list[1]
    for rec in rec_list:
        rec_stmts = [] #this works since there's one rec for now, need to use different object or list of lists when dealing with multiple recs, but I'm on a deadline
        for st in model_stmts:
            for ag in st.agent_list():    #rec is string with no () now
                if ag:
                    if rec.entity_matches(ag):
                        rec_stmts.append(st)
            else:
                non_rec_stmts.append(st)

    #This is all receptor statements, both with and without ligand, next step will split those
    rec_stmts = ac.Preassembler.combine_duplicate_stmts(rec_stmts)
    non_rec_stmts = ac.Preassembler.combine_duplicate_stmts(non_rec_stmts)

    #This can almost certainly be folded into above to save a step
    rec_lig_stmts = []
    rec_nolig_stmts = []
    for rec in rec_list: #just egfr 
        for lig in lig_list:
            for st in rec_stmts:
    #            if rec in str(st.agent_list()):
#                if rec in st.agent_list(): #this is guaranteed
                if str(lig).split(')')[0] in str(st.agent_list()):
                    
                    rec_lig_stmts.append(st)
                else:
                    rec_nolig_stmts.append(st)

    rec_lig_stmts = ac.Preassembler.combine_duplicate_stmts(rec_lig_stmts)
    rec_nolig_stmts = ac.Preassembler.combine_duplicate_stmts(rec_nolig_stmts)
    #Now have three sets of statements, ones without receptor, ones with receptor and ligand, ones with receptor and no ligand
    output_stmts = []
    output_stmts.append(rec_lig_stmts)
    output_stmts.append(rec_nolig_stmts)
    output_stmts.append(non_rec_stmts)
    return output_stmts #list of three lists


#New function for generarting context in statments
def add_ligand_context(stmts,rec_lig_list):
    rec_lig_stmts = stmts[0]
    rec_nolig_stmts = stmts[1]
    non_rec_stmts = stmts[2]
    rec_list = rec_lig_list[0]
    lig_list = rec_lig_list[1]

    for st1 in rec_lig_stmts:
        for ag in st1.agent_list():
            if ag in lig_list: #find the ligand for a given receptor so we can add it as context for all receptor stmts. There is probably a smarter way to do this by dropping a function and folding this into above
                lig_ag = ag #need ligand agent so we can add it as context to additional receptor binding statements
                for st2 in rec_nolig_stmts:
                    for ag in st2.agent_list():
                        if ag in rec_list: 
                            rec_ag = ag
                            rec_ag_index = st2.agent_list().index(ag)

                    st2.agent_list()[rec_ag_index].bound_conditions = st2.agent_list()[rec_ag_index].bound_conditions + [BoundCondition(lig_ag)] #have to make additive to get 'or' behavior
    #Make sure each conext element is present only once. Definitely better way to do this. 
    new_bound_conditions = []
    for st in rec_nolig_stmts:
        for ag in st.agent_list():
            if len(ag.bound_conditions) > 1:
                new_bound_conditions.append(ag.bound_conditions[0])
                for bc in ag.bound_conditions[1:]:
                    if bc.agent.entity_matches(new_bound_conditions[-1].agent):
                        pass
                    else:
                        new_bound_conditions.append(bc)
                ag.bound_conditions = new_bound_conditions

    #Now want a single list of usable statements. rec_nolig_stmts have been rewritten to add context. other two sets should remain the same.
    output_stmts = []
    output_stmts.append(rec_lig_stmts)
    output_stmts.append(rec_nolig_stmts)
    output_stmts.append(non_rec_stmts)
    return output_stmts #list of three lists



#This works for complex statments
def find_primary_binders_complex(stmts,rec_lig_list): 
    rec_list = rec_lig_list[0]
    key = None #are these necessary?
    val = None
    primary_binders_dict = {}
    for st in stmts:
        key = None
        for ag in st.agent_list():
            for rec in rec_list:
                if ag.entity_matches(rec):
                    val = ag
                else:
                    key = ag
        if key: #why was this necessary?
            primary_binders_dict[key] = val
#            key = None
    return primary_binders_dict

#This works now. For Phosphorylation and Autophosphorylation statements
def find_primary_binders_phos(stmts,rec_lig_list): 
    rec_list = rec_lig_list[0]
    key = None #are these necessary?
    val = None
    primary_binders_dict = {}
    for st in stmts:
        key = None
        for ag in st.agent_list():
            if ag not in rec_list:
                key = ag 
            res = st.residue    #Note: res and pos can be None
            pos = st.position 
            val = ModCondition(mod_type='phosphorylation',residue=res,position=pos)
#            print(val)
        if key: #why was this necessary?
            primary_binders_dict[key] = val
    return primary_binders_dict

#######################
#### SPLIT HERE #######
#######################

def find_primary_stmts(primary_binders_dict,rec_lig_list,model_stmts):
    rec_list = rec_lig_list[0]
    #finding relevant stmts
    primary_st=[]
    non_primary_st = []
    #loop through stmts for proteins at this level (in this case direct binders to receptor
    for prot in primary_binders_dict.keys():
        flag = 0
        #loop through all stmts to find stmts that contain this protein. is there a faster way to do this?
        for st in model_stmts:                 ##########VARIABLE STATEMENT LIST HERE. Contains all statements in model 
            flag = max([flag+1 if prot.entity_matches(ag) else flag+0 for ag in st.agent_list()])
            if flag == 0:
                non_primary_st.append(st)
            else:
                primary_st.append(st)   #should be checking for receptor, not sure how
                flag = 0
                #may need to check for edge cases where flag > 1, but I think they still fit in this category

    #being replaced, ignore for now
    #removing stmts with rec
    #Should this be an additional condition on the flag above rather than a separate step?
    primary_st_test = primary_st[:]
    for st in primary_st:
        for ag in st.agent_list():
            if ag in rec_list:
                primary_st_test.remove(st) #removing stmt where primary binds ligand, so i can add context to the others
                non_primary_st.append(st)
#    output = primary_st_test + 
    return primary_st_test, non_primary_st


