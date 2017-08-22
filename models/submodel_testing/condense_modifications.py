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
                generic_pST = ModCondition(pST_mods[0].mod_type,pST_mods[0].residue) #This could appear as S or T but I'm assuming equivaqlency
                modlist.append(generic_pST)
            except IndexError:
                generic_pST = []
            #Write replacement ActiveForm statement with coarse-grained modifications
            newagent = Agent(act_agent.name,modlist)
            newst = ActiveForm(newagent,st.activity,st.is_active)
            new_af_stmts.append(newst)
        else:
            new_af_stmts.append(st)
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
    phos_stmts = ac.filter_by_type(stmts, Phosphorylation)
    non_phos_stmts = ac.filter_by_type(stmts, Phosphorylation, invert = True)
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
#
def coarse_grain_kinase_context(stmts):
    phos_stmts = ac.filter_by_type(stmts, Phosphorylation)
    non_phos_stmts = ac.filter_by_type(stmts, Phosphorylation, invert = True)
    kinase_with_mods = filter(lambda x: x.enz.mods != [], phos_stmts)
    kinase_wo_mods = filter(lambda x: x.enz.mods == [], phos_stmts)
    
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
            generic_pST = ModCondition(pST_mods[0].mod_type,pST_mods[0].residue) #This could appear as S or T but I'm assuming equivaqlency
            modlist.append(generic_pST)
        except IndexError:
            generic_pST = []
        #Write replacement ActiveForm statement with coarse-grained modifications
        newagent = Agent(st.enz.name,modlist)
        newst = Phosphorylation(newagent,st.sub,st.residue)
        coarse_grained_phos_stmts.append(newst)
    new_phos_stmts = coarse_grained_phos_stmts+kinase_wo_mods
    filtered_stmts = non_phos_stmts + new_phos_stmts
    return filtered_stmts

stmts = my_model_stmts_larger[:]
stmts = coarse_grain_phos_on_af(stmts)
stmts = remove_redundant_phosphorylations(stmts)
stmts = coarse_grain_kinase_context(stmts)

pa = PysbAssembler('two_step')
pa.add_statements(my_model_stmts_larger)
model = pa.make_model()

bngl_model_filter = pysb.export.export(model,'bngl')
bngl_file_filter = open('rbm/egfr_small_no_cg.bngl','w')
bngl_file_filter.write(bngl_model_filter)
bngl_file_filter.close()


pa = PysbAssembler('two_step')
pa.add_statements(stmts)
model = pa.make_model()

bngl_model_filter = pysb.export.export(model,'bngl')
bngl_file_filter = open('rbm/group_meeting_6_5_17/small.bngl','w')
bngl_file_filter.write(bngl_model_filter)
bngl_file_filter.close()





#ac.dump_statements('egf_egfr_sos_grb_olig_stmts.pkl',stmts)

