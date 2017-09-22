#def find_rec_lig_pair(stmts_in):
import pickle

with open('/home/bobby/Dropbox/Sorger_Lab/BigMech/reporting_and_eval/phase3_eval_pt2/dynamical_model/model_building/lig_rec_dict.pkl','rb') as f:
    receptor_dict = pickle.load(f)

#send list of stmts, loop through stmts inside function, check for complex stmts, loop through all agents, check in if in ligand list
#stmts = ac.load_statements('my_reading/stmts_preassembled_fixed_be_filtering.pkl')
#stmts_old = ac.load_statements('/home/bobby/Dropbox/Sorger_Lab/BigMech/reporting_and_eval/phase3_eval_pt2/dynamical_model/model_building/secondpass_model_stmts.pkl')

#stmts = ac.load_statements('')


#AUTOPHOS IS IMPORTANT TO HANDLE NOW

stmts_in = stmts_final 




rec_list = []
lig_list_pre = []
for st in stmts_in:
#    if isinstance(st, Complex):
        for ag in st.agent_list():
            for lig,rec in receptor_dict.items():
                if str(ag.name).strip('()') in rec:
#                    print(str(ag))
                    lig_list_pre.append(lig+'()')
                    rec_list.append(ag)

lig_list = []
for lig in lig_list_pre:
    for st in stmts_in:
        if lig in str(st.agent_list()):
            lig_list.append(lig)



rec_list = list(set(rec_list)) #list off all receptors in stmt list 
lig_list = list(set(lig_list)) #list of all ligands associated with these receptors, not necessarily in stmts.

rec_stmts = []
non_rec_stmts = []
for rec in rec_list:
    rec_stmts = [] #this works since there's one rec for now, need to use different object or list of lists when dealing with multiple recs, but I'm on a deadline
    for st in stmts_in:
        for ag in st.agent_list():    #rec is string with no () now
            if rec.entity_matches(ag):
                rec_stmts.append(st)
        else:
            non_rec_stmts.append(st)

rec_stmts = ac.Preassembler.combine_duplicate_stmts(rec_stmts)
non_rec_stmts = ac.Preassembler.combine_duplicate_stmts(non_rec_stmts)


#rec_list has all rec (in this case one)
#rec_stmts has all stmts for each rec (in this case it's a single list, may become list of list or a different data structure


rec_lig_st = []
rec_nolig_st = []
for rec in rec_list: #just egfr 
    for lig in lig_list:
        for st in rec_stmts:
#            if rec in str(st.agent_list()):
            if rec in st.agent_list():
                if lig in str(st.agent_list()):
                    rec_lig_st.append(st)

for st in rec_stmts:
    if st not in rec_lig_st:
        rec_nolig_st.append(st) #still getting hbegf? this comes from looping through ligs

rec_lig_st = list(set(rec_lig_st))
rec_nolig_st = list(set(rec_nolig_st))
 #ligand bound, effectively making an or gate on all ligands


new_rec_nolig_st = []
rec_nolig_st_copy = rec_nolig_st[:]
rec_ag_index=[]
rec_ag=[]
for st1 in rec_lig_st:
    for ag in st1.agent_list():
        if str(ag) in lig_list: #find the ligand for a given receptor so we can add it as context for all receptor stmts
            lig_ag = ag 
    for st2 in rec_nolig_st_copy:
        for ag in st2.agent_list():
            if ag in rec_list: 
                rec_ag = ag
                rec_ag_index = st2.agent_list().index(ag)
            #for now this assums only one receptor agent 
            #this fails on receptor dimerizations, in which both should have ligand-binding context. 
            #not a terrible failure for now because cascade will be reliant on ligand binding, but ahould be fixed
#        if rec_ag_index:      
#PROBLEM HERE  
        st2.agent_list()[rec_ag_index].bound_conditions = st.agent_list()[rec_ag_index].bound_conditions + [BoundCondition(lig_ag)] #have to make additive to get 'or' behavior
        new_rec_nolig_st.append(st2)


new_bound_conditions = []
for st in new_rec_nolig_st:
    for ag in st.agent_list():
        if len(ag.bound_conditions) > 1:
            new_bound_conditions.append(ag.bound_conditions[0])
            for bc in ag.bound_conditions[1:]:
                if bc.agent.entity_matches(new_bound_conditions[-1].agent):
                    pass
                else:
                    new_bound_conditions.append(bc)
            ag.bound_conditions = new_bound_conditions
            
            

new_rec_nolig_st = ac.Preassembler.combine_duplicate_stmts(new_rec_nolig_st)

final_st = non_rec_stmts + new_rec_nolig_st
final_st = ac.Preassembler.combine_duplicate_stmts(final_st)


