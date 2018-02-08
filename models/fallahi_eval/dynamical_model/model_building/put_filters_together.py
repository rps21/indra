stmts_start = all_stmts[:]

stmts_new = remove_dimers(stmts_start)
stmts_new = remove_mutations(stmts_new)
stmts_new = coarse_grain_phos_on_af(stmts_new)
stmts_new = remove_redundant_phosphorylations(stmts_new)
stmts_new = coarse_grain_kinase_context(stmts_new)
stmts_new = coarse_grain_substrate_context(stmts_new)
#stmts_new = ac.filter_inconsequential_mods(stmts_new)
#stmts_new = ac.filter_inconsequential_acts(stmts_new)

stmts_final = stmts_new[:]

##egf-egfr complext
##egfr_grb2 complex
##kras-braf activate
##egfr_akt phosphorylate


##ADD MISSING STMTS NECESSARY FOR MODEL
#for st in stmts_old:
#    for ag in st.agent_list():
#        if 'EGF(' in str(ag):
#            egf_agent = ag
#        elif 'EGFR(' in str(ag):
#            egfr_agent = ag
#            egfr_agent.mods = []
#        elif 'GRB2(' in str(ag):
#            grb2_agent = ag
#            grb2_agent.mods = []
#        elif 'KRAS(' in str(ag):
#            kras_agent = ag 
#            kras_agent.mods = []
#        elif 'AKT1(' in str(ag):
#            akt_agent = ag
#            akt_agent.mods = []
#        elif 'BRAF(' in str(ag):
#            braf_agent = ag
#            braf_agent.mods = []



#missing_stmts=[]
#missing_stmts.append(Complex([egf_agent,egfr_agenet]))
#missing_stmts.append(Complex([egfr_agent,grb2_agent]))
#missing_stmts.append(Activation(kras_agent,braf_agent))
#missing_stmts.append(Phosphorylation(egfr_agent,akt_agent))
#stmts = stmts + missing_stmts
#ac.dump_statements(stmts,'fallahi_eval_pysb_stmts_updated.pkl')
#stmts_check = ac.load_statements('fallahi_eval_pysb_stmts_updated.pkl')

erk_model_st = ac.filter_gene_list(final_st,['EGF','EGFR','SOS1','GRB2','KRAS','BRAF','MAPK1','MAPK3','MAP2K1'],'all')
erk_model_st_reduced = ac.filter_by_type(erk_model_st,Translocation,invert=True)




pa = PysbAssembler('two_step')
pa.add_statements(newstmts)
model1 = pa.make_model()

bngl_model_filter = pysb.export.export(model1,'bngl')
bngl_file_filter = open('rbm/testing_nodes.bngl','w')
bngl_file_filter.write(bngl_model_filter)
bngl_file_filter.close()

with open('rbm/final_version/erk_model_final2.txt','w') as f:
    for item in erk_model_st_reduced:
        f.write("%s\n" % item)

ac.dump_statements(erk_model_st_reduced,'rbm/final_version/erk_model_final_stmts.pkl')

#remove hspb1, akt3, 
#rptor, grb2, kras not connectected to anything
