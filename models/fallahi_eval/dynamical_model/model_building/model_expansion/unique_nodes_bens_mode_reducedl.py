['RASGRF2', 'RASA2', 'MAP3K8',  'PLK1', 'BRAF', 'JUN', 'MAPK3', 'MAP2K2', 'RHEB', 'SOS1', 'PRKAA2', 'HRAS', 'KSR1', 'MYC', 'RASGRP4', 'ATF3', 'PAK2', 'ERBB2', 'RALGDS', 'RASAL1', 'RASGRP1', 'TBK1', 'ARHGAP35', 'PLCE1', 'RAC1',  'SHC1', 'MAPK14', 'VAV1', 'DUSP6', 'PAK3', 'DUSP4', 'MAPK10', 'MAP2K1', 'PDPK1', 'MAPK9', 'PRKAA1', 'RASA1', 'MYB', 'KRAS', 'E2F7', 'JUND', 'RASGRF1',  'DUSP1', 'RELA', 'MAPK1', 'PRKAB1', 'NRAS', 'SREBF1', 'ELK1', 'CDKN1B', 'AURKA', 'RAF1', 'RASGRP3',  'MAPK8', 'PIN1',  'PAK1',  'PTPN11', 'EGFR','GRB2','EGF']


erk_model_st = ac.filter_gene_list(final_st,['EGF','EGFR','SOS1','GRB2','KRAS','BRAF','MAPK1','MAPK3','MAP2K1'],'all')


newnodes=['RASA2', 'MAP3K8', 'BRAF', 'JUN', 'MAPK3', 'MAP2K2', 'SOS1', 'PRKAA2', 'HRAS', 'PAK2', 'RASAL1', 'RASGRP1', 'TBK1', 'RAC1',  'SHC1', 'MAPK14', 'DUSP6', 'PAK3', 'DUSP4', 'MAPK10', 'MAP2K1', 'MAPK9', 'PRKAA1', 'RASA1', 'KRAS',  'DUSP1', 'RELA', 'MAPK1', 'PRKAB1', 'NRAS', 'CDKN1B', 'AURKA', 'RAF1', 'RASGRP3',  'MAPK8', 'PIN1',  'PAK1',  'PTPN11', 'EGFR','GRB2','EGF']


erk_expanded_st = ac.filter_gene_list(stmts,newnodes,'all')


new_nodes = filter_dead_end(erk_expanded_st,newnodes)


erk_expanded_st_tweaked = ac.filter_gene_list(stmts,new_nodes,'all')


stmts_start = erk_expanded_st_tweaked[:]

stmts_new = remove_dimers(stmts_start)
stmts_new = remove_mutations(stmts_new)
stmts_new = coarse_grain_phos_on_af(stmts_new)
stmts_new = remove_redundant_phosphorylations(stmts_new)
stmts_new = coarse_grain_kinase_context(stmts_new)
stmts_new = coarse_grain_substrate_context(stmts_new)
#stmts_new = ac.filter_inconsequential_mods(stmts_new)
#stmts_new = ac.filter_inconsequential_acts(stmts_new)

stmts_final = stmts_new[:]


pa = PysbAssembler('two_step')
pa.add_statements(stmts_final)
model1 = pa.make_model()

bngl_model_filter = pysb.export.export(model1,'bngl')
bngl_file_filter = open('rbm/testing_nodes.bngl','w')
bngl_file_filter.write(bngl_model_filter)
bngl_file_filter.close()


#things for further simplification:
#removing redundant family members (criteria?)
#expanding context extrapolation beyond first two layers
#be critical of transcription factors (are they transcribing? if not are they an output? if not remove)
#need EGFR, SOS, GRB stmts from last model. Also need to eventually figure out why they're missing. 

#../../fallahi_eval_pysb_stmts_updated.pkl
old_stmts = ac.load_statements('../../fallahi_eval_pysb_stmts_updated.pkl')

#In [74]: old_stmts[828:]
#Out[74]: 
#[Ubiquitination(EGF(), EGFR()),
# Complex(GRB2(), SOS1()),
# Complex(EGFR(), GRB2()),
# Activation(KRAS(), BRAF()),
# Phosphorylation(EGFR(), AKT1())]

stmts_final = stmts_final + old_stmts[828:]
#after rec_lig now using final_st
#s

pa = PysbAssembler('two_step')
pa.add_statements(final_st)
model1 = pa.make_model()

bngl_model_filter = pysb.export.export(model1,'bngl')
bngl_file_filter = open('rbm/testing_nodes2.bngl','w')
bngl_file_filter.write(bngl_model_filter)
bngl_file_filter.close()







