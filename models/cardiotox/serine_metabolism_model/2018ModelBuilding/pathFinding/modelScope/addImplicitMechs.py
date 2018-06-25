

effectType = {'PKM':['phosphorylation','SORAFENIB phosphorylates PKM'],'RPS6':['dephosphorylation','SORAFENIB phosphorylates RPS6'],'AURKA':['phosphorylation','SORAFENIB phosphorylates AURKA'],'HIF1A':['increaseamount','SORAFENIB transcribes HIF1A'],'MYC':['increaseamount','SORAFENIB transcribes MYC'],'JUN':['phosphorylation','SORAFENIB phosphorylates JUN']}#,'STAT1':['phosphorylation','SORAFENIB phosphorylates STAT1']}#,'PDGFRA':'increased'}


#eventually replace this with something from indra db potentially
#Notes: double check what stmts to run through here. Have to be careful about ptm stuff
implicit_agent = Agent('GenericAgent')
for key in effectType:
    agentName = key
    mechType = effectType[key][0]

    if effectType == 'phosphorylation':
        dephosStmts = ac.filter_by_type(stmts,Dephosphorylation)
        originalDephos = ac.filter_gene_list(dephosStmts,[agentName],'one')
        correctDephos = []
        for st in originalDephos:
            if st.enz.name == agentName:
                correctDephos.append(st)
        if not correctDephos:        
            filtStmts = ac.filter_gene_list(stmts,[agentName],'all') #HERE
            dephos_ag = deepcopy(filtStmts[0].agent_list()[0])
            dephos_mod = ModCondition(mod_type='phosphorylation',residue='phos_act')
            dephos_ag.mods = [dephos_mod]

            dephos_st = Dephosphorylation(implicit_agent, dephos_ag, residue=dephos_res)
            new_dephos.append(dephos_st)


    elif effectType == 'dephosphorylation':
        phosStmts = ac.filter_by_type(stmts,Phosphorylation)
        originalPhos = ac.filter_gene_list(PhosStmts,[agentName],'one')
        correctPhos = []
        for st in originalPhos:
            if st.enz.name == agentName:
                correctPhos.append(st)
        if not correctPhos:        
            filtStmts = ac.filter_gene_list(stmts,[agentName],'all') #HERE
            phos_ag = deepcopy(filtStmts[0].agent_list()[0])
            phos_mod = ModCondition(mod_type='dephosphorylation',residue='phos_act')
            phos_ag.mods = [phos_mod]

            phos_st = Phosphorylation(implicit_agent, phos_ag, residue=phos_res)
            new_phos.append(phos_st)


    elif effectType == 'increaseamount':
        decreaseStmts = ac.filter_by_type(stmts,DecreaseAmount)
        originalDecrease = ac.filter_gene_list(decreaseStmts,[agentName],'one')
        correctDecrease = []
        for st in correctDecrease:
            if st.subj.name == agentName:
                correctDecrease.append(st)
        if not correctDecrease:        
            filtStmts = ac.filter_gene_list(stmts,[agentName],'all')    #HERE
            decrease_ag = deepcopy(filtStmts[0].agent_list()[0])
#            dephos_mod = ModCondition(mod_type='phosphorylation',residue='phos_act')
#            dephos_ag.mods = [dephos_mod]

            decrease_st = DecreaseAmount(implicit_agent, decrease_ag)
            new_decrease.append(decrease_st)

    elif effectType == 'decreaseamount':
        increaseStmts = ac.filter_by_type(stmts,IncreaseAmount)
        originalIncrease = ac.filter_gene_list(increaseStmts,[agentName],'one')
        correctIncrease = []
        for st in correctIncrease:
            if st.subj.name == agentName:
                correctIncrease.append(st)
        if not correctIncrease:        
            filtStmts = ac.filter_gene_list(stmts,[agentName],'all')    #HERE
            increase_ag = deepcopy(filtStmts[0].agent_list()[0])
#            dephos_mod = ModCondition(mod_type='phosphorylation',residue='phos_act')
#            dephos_ag.mods = [dephos_mod]

            increase_st = IncreaseAmount(implicit_agent, increase_ag)
            new_increase.append(increase_st)



