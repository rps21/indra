#Enforce dephos and degredation
#Check unbinding context 

#TODO:
#Break up multiple phosp sites
#Ex:  Dephosphorylation(None, BRAF(mods: (phosphorylation, S, 446), (phosphorylation, T, 599), (phosphorylation, S, 602))),

#stmts = erk_model_st_reduced

##This will error on non catalytic statements, add try
##For now only calling with list of phos/dephos statements
#def collect_substrates(stmts):
#    substrate_list = []
#    for st in stmts:
#        new_sub = st.sub
##        print(new_sub)
#        if any(list(map(lambda x: x.matches(new_sub),substrate_list))):
#            pass
#        else:
#            substrate_list.append(new_sub)
#    return substrate_list


#def add_dephosphorylations(stmts):
#    output_stmts = deepcopy(stmts)
#    new_dephos = []
#    phos_stmts = ac.filter_by_type(stmts,Phosphorylation)
#    dephos_stmts = ac.filter_by_type(stmts,Dephosphorylation)
#    phos_subs = collect_substrates(phos_stmts)
#    dephos_subs = collect_substrates(dephos_stmts)

#    for sub in phos_subs:
#        if not any(list(map(lambda x: x.matches(sub),dephos_subs))):
#            new_dephos.append(Dephosphorylation(None, deepcopy(sub)))
#    
#    output_stmts = output_stmts + new_dephos
#    return output_stmts

#stmts = erk_model_st_reduced
#stmts = add_dephosphorylations(stmts)




#############################################################3
#Phosphorylation(BRAF(mods: (phosphorylation, S)), MAP2K1(), S),
#Want
#Dephosphorylation(None, MAP2K1(), S)

#Alternative to above, for every phos stmt immediately make a matching, generic dephos stmts. 
#Less elegant, may have repeated stmts of overlap with explicit phosphotase, but both workable.

#>>> dusp6 = Agent('DUSP6')
#>>> erk = Agent('MAPK1')
#>>> dephos = Dephosphorylation(dusp6, erk, 'T', '185')

def add_dephosphorylations(stmts):
    implicit_phosphotase = Agent('ImplicitPhos')
    output_stmts = deepcopy(stmts)
    new_dephos = []
    phos_stmts = ac.filter_by_type(stmts,Phosphorylation)

    for st in phos_stmts:
        dephos_ag = deepcopy(st.sub)
        dephos_res = deepcopy(st.residue)
        dephos_pos = deepcopy(st.position)

        dephos_mod = ModCondition(mod_type='phosphorylation',residue=dephos_res,position=dephos_pos)
        dephos_ag.mods = [dephos_mod]

        dephos_st = Dephosphorylation(implicit_phosphotase, dephos_ag, residue=dephos_res, position=dephos_pos)
        new_dephos.append(dephos_st)
        #print(st)
    new_dephos = Preassembler.combine_duplicate_stmts(new_dephos)
    output_stmts = output_stmts + new_dephos
    return output_stmts

#new_st = add_dephosphorylations(erk_model_st_reduced)


#    dephos_stmts = ac.filter_by_type(stmts,Dephosphorylation)
#    phos_subs = collect_substrates(phos_stmts)
#    dephos_subs = collect_substrates(dephos_stmts)

#    for sub in phos_subs:
#        if not any(list(map(lambda x: x.matches(sub),dephos_subs))):
#            new_dephos.append(Dephosphorylation(None, deepcopy(sub)))
#    
#    output_stmts = output_stmts + new_dephos
#    return output_stmts




















