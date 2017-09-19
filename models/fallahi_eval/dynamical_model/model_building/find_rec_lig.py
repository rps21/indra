#def find_rec_lig_pair(stmts_in):
import pickle

with open('/home/bobby/Dropbox/Sorger_Lab/BigMech/reporting_and_eval/phase3_eval_pt2/model_building/lig_rec_dict.pkl','rb') as f:
    receptor_dict = pickle.load(f)

#send list of stmts, loop through stmts inside function, check for complex stmts, loop through all agents, check in if in ligand list
#stmts = ac.load_statements('my_reading/stmts_preassembled_fixed_be_filtering.pkl')
stmts = ac.load_statements('/home/bobby/Dropbox/Sorger_Lab/BigMech/reporting_and_eval/phase3_eval_pt2/model_building/secondpass_model_stmts.pkl')

stmts = ac.load_statements('')




stmts_in = stmts 

#rec_list = []
#for st in stmts_in:
#    if isinstance(st, Complex):
#        for ag in st.agent_list():
#            if ag.name in receptor_dict.keys():
#                lig = str(lig).strip('()')
#                rec = receptor_dict[lig] #list of receptors
#                print(rec)
#                rec_list = rec_list + rec
#rec_list = list(set(rec_list))        

#might be better to loop through values, each being a list, may or may not need the whole list to later recover ligand. How to handle multiple ligands per receptor? Base it on which ligs in list of stmts

rec_list = []
lig_list_pre = []
for st in stmts_in:
#    if isinstance(st, Complex):
        for ag in st.agent_list():
            for lig,rec in receptor_dict.items():
                if str(ag.name).strip('()') in rec:
#                    print(str(ag))
                    lig_list_pre.append(lig+'()')
                    rec_list.append(str(ag))

lig_list = []
for lig in lig_list_pre:
    for st in stmts_in:
        if lig in str(st.agent_list()):
            lig_list.append(lig)



rec_list = list(set(rec_list)) #list off all receptors in stmt list 
lig_list = list(set(lig_list)) #list of all ligands associated with these receptors, not necessarily in stmts.

rec_stmts = []
for rec in rec_list:
    rec_stmts = [] #this works since there's one rec for now, need to use different object or list of lists when dealing with multiple recs, but I'm on a deadline
    for st in stmts_in:
        for ag in st.agent_list():    #rec is string with no () now
            if rec in str(ag):
                rec_stmts.append(st)

rec_stmts = list(set(rec_stmts))
#this is limited for now, didn't included egf in first pass and grb2 is missing for some reason and i had manually added.
#Out[60]: [Complex(EGFR(), EGFR()), Phosphorylation(EGFR(), AKT1())]

#rec_stmts = ac.load_statements('egfr_context_testing.pkl')
#rec_list has all rec (in this case one)
#rec_stmts has all stmts for each rec (in this case it's a single list, may become list of list or a different data structure


rec_lig_st = []
rec_nolig_st = []
for rec in rec_list: #just egfr 
    for lig in lig_list:
        for st in rec_stmts:
            if rec in str(st.agent_list()):
                if lig in str(st.agent_list()):
                    rec_lig_st.append(st)

for st in rec_stmts:
    if st not in rec_lig_st:
        rec_nolig_st.append(st) #still getting hbegf? this comes from looping through ligs

rec_lig_st = list(set(rec_lig_st))
rec_nolig_st = list(set(rec_nolig_st))

#In [81]: rec_lig_st
#Out[81]: [Complex(HBEGF(), EGFR()), Ubiquitination(EGF(), EGFR())]
#where is hbegf coming from? regardless makes for good test example
#why am i missing egf? also need to filter by complex, thought i was already

#pull context(binding site) from rec_lig_st
#acutally, maybe add bound conditions to rec_nolig_st
#rec_lig_st[0].agent_list()[0].bound_conditions

#loop through ligand binding statements.
#take ligand
#loop through no_ligand
#make new list 
#add new st with bound_conditions
#should have one st per ligand bound, effectively making an or gate on all ligands

new_rec_nolig_st = []
rec_nolig_st_copy = rec_nolig_st[:]
for st1 in rec_lig_st:
    for ag in st1.agent_list():
        if str(ag) in lig_list: #find the ligand for a given receptor so we can add it as context for all receptor stmts
            lig_ag = ag 
    for st2 in rec_nolig_st_copy:
        for ag in st2.agent_list():
            if str(ag) in rec_list: 
                rec_ag = ag
                rec_ag_index = rec_list.index(str(ag))
            #for now this assums only one receptor agent 
            #this fails on receptor dimerizations, in which both should have ligand-binding context. 
            #not a terrible failure for now because cascade will be reliant on ligand binding, but ahould be fixed
        st2.agent_list()[rec_ag_index].bound_conditions = st.agent_list()[rec_ag_index].bound_conditions + [BoundCondition(lig_ag)] #have to make additive to get 'or' behavior
        new_rec_nolig_st.append(st2)

new_rec_nolig_st = list(set(new_rec_nolig_st))








#def filter_enzyme_kinase(stmts_in, **kwargs):
#    """Filter Phosphorylations to ones where the enzyme is a known kinase.

#    Parameters
#    ----------
#    stmts_in : list[indra.statements.Statement]
#        A list of statements to filter.
#    save : Optional[str]
#        The name of a pickle file to save the results (stmts_out) into.

#    Returns
#    -------
#    stmts_out : list[indra.statements.Statement]
#        A list of filtered statements.
#    """
#    logger.info('Filtering %d statements to remove ' % len(stmts_in) +
#                'phosphorylation by non-kinases...')
#    path = os.path.dirname(os.path.abspath(__file__))
#    kinase_table = read_unicode_csv(path + '/../resources/kinases.tsv',
#                                    delimiter='\t')
#    gene_names = [lin[1] for lin in list(kinase_table)[1:]]
#    stmts_out = []
#    for st in stmts_in:
#        if isinstance(st, Phosphorylation):
#            if st.enz is not None:
#                if st.enz.name in gene_names:
#                    stmts_out.append(st)
#        else:
#            stmts_out.append(st)
#    logger.info('%d statements after filter...' % len(stmts_out))
#    dump_pkl = kwargs.get('save')
#    if dump_pkl:
#        dump_statements(stmts_out, dump_pkl)
#    return stmts_out
