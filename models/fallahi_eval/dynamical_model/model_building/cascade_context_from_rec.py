#this has to follow rec_lig finding
#everything binding to rec will have ligand context
#continue to carry this down cascade
#every complex reaction must have first ag bound to rec/other adaptor to allow 2nd ag binding
#This will build complexes, phosphorylations should be break points
#have to watch for redundant complex/phosphorylation statments which will allow binding that isn't breaking the complex incorrectly
#this should be taken care of by the mech linker/preassebler but a couple seemed to slip through, not sure why.

#going to have to iteratively build a list
#already can first build list of rec-ligand 
#then rec-nonligand
#need to use that list of nonligand direct binders to get next level of binders
#iterate stopping on an enzymatic (non Complex?) statement 


stmts = ac.load_statements('/home/bobby/Dropbox/Sorger_Lab/BigMech/reporting_and_eval/phase3_eval_pt2/model_building/secondpass_model_stmts.pkl')

#save rec stmts?
#call other function?
#having rec stmts be an input variable is a good idea, will just save and load them for now
#can save rec_list or rec_stmts.

rec_lig_st = list(set(rec_lig_st))
rec_nolig_st = list(set(rec_nolig_st))
#so gather list of 'primary binders' from rec_nolig_st.
#loop through and find all stmts with a primary binder and no receptor, in new list
#change context for these stmts.
#if a Complex stmts, populate new list of 'secondary binders.



rec_nolig_st_tosave = rec_nolig_st[:]
rec_nolig_st_touse = rec_nolig_st[:]


primary_binders = []
primary_binders_dict = {}
key = None
val = None
#for st in rec_nolig_st:
for st in rec_nolig_st_touse:
#    print(st)
    for ag in st.agent_list():
        if ag not in rec_list:
#            primary_binders.append(ag) #maybe here we can save the stmts that go with the primary binder, maybe in a dict? list of lists?
            key = ag 
        else:
            val = ag
    if key:
        primary_binders_dict[key] = val
#    print(key)


primary_st=[]
for prot in primary_binders_dict.keys():
    for st in stmts:
        for ag in st.agent_list():
            if prot.entity_matches(ag):
                primary_st.append(st)   #should be checking for receptor, not sure how
    
#
primary_st_test = primary_st[:]
for st in primary_st:
    for ag in st.agent_list():
        if ag in rec_list:
            primary_st_test.remove(st)



for st in primary_st_test:
    for ag in st.agent_list():
#        print(ag)
#        if ag in primary_binders_dict.keys():
        for ag2 in primary_binders_dict.keys():
            if ag.entity_matches(ag2):
                print(ag)
#context here. should be primary binder is bound to rec 
                ag.bound_conditions = ag.bound_conditions + [BoundCondition(primary_binders_dict[ag2])] #FIX REC REFERENCE actually just need coresponding receptor, not whole stmt 



#should continue this on with mtor - Complex(AKT1(bound: [EGFR, True]), MTOR()),
#all other stmts are phosphorylation
#need to encapsulate this into some type of while loop that continues while there are options for expansion
#While loop that calls a function as many times as necessary, after encoding this in a funciton








#rec_stmts = []
#for rec in rec_list:
#    rec_stmts = [] #this works since there's one rec for now, need to use different object or list of lists when dealing with multiple recs, but I'm on a deadline
#    for st in stmts_in:
#        for ag in st.agent_list():    #rec is string with no () now
#            if rec in str(ag):
#                rec_stmts.append(st)

#rec_stmts = list(set(rec_stmts))
##this is limited for now, didn't included egf in first pass and grb2 is missing for some reason and i had manually added.
##Out[60]: [Complex(EGFR(), EGFR()), Phosphorylation(EGFR(), AKT1())]

##rec_stmts = ac.load_statements('egfr_context_testing.pkl')
##rec_list has all rec (in this case one)
##rec_stmts has all stmts for each rec (in this case it's a single list, may become list of list or a different data structure


#rec_lig_st = []
#rec_nolig_st = []
#for rec in rec_list: #just egfr 
#    for lig in lig_list:
#        for st in rec_stmts:
#            if rec in str(st.agent_list()):
#                if lig in str(st.agent_list()):
#                    rec_lig_st.append(st)

#for st in rec_stmts:
#    if st not in rec_lig_st:
#        rec_nolig_st.append(st) #still getting hbegf? this comes from looping through ligs

#rec_lig_st = list(set(rec_lig_st))
#rec_nolig_st = list(set(rec_nolig_st))

##In [81]: rec_lig_st
##Out[81]: [Complex(HBEGF(), EGFR()), Ubiquitination(EGF(), EGFR())]
##where is hbegf coming from? regardless makes for good test example
##why am i missing egf? also need to filter by complex, thought i was already

##pull context(binding site) from rec_lig_st
##acutally, maybe add bound conditions to rec_nolig_st
##rec_lig_st[0].agent_list()[0].bound_conditions

##loop through ligand binding statements.
##take ligand
##loop through no_ligand
##make new list 
##add new st with bound_conditions
##should have one st per ligand bound, effectively making an or gate on all ligands

#new_rec_nolig_st = []
#rec_nolig_st_copy = rec_nolig_st[:]
#for st1 in rec_lig_st:
#    for ag in st1.agent_list():
#        if str(ag) in lig_list: #find the ligand for a given receptor so we can add it as context for all receptor stmts
#            lig_ag = ag 
#    for st2 in rec_nolig_st_copy:
#        for ag in st2.agent_list():
#            if str(ag) in rec_list: 
#                rec_ag = ag
#                rec_ag_index = rec_list.index(str(ag))
#            #for now this assums only one receptor agent 
#            #this fails on receptor dimerizations, in which both should have ligand-binding context. 
#            #not a terrible failure for now because cascade will be reliant on ligand binding, but ahould be fixed
#        st2.agent_list()[rec_ag_index].bound_conditions = st.agent_list()[rec_ag_index].bound_conditions + [BoundCondition(lig_ag)] #have to make additive to get 'or' behavior
#        new_rec_nolig_st.append(st2)

#new_rec_nolig_st = list(set(new_rec_nolig_st))



