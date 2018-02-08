import indra

stmts = ac.load_statements('../../fallahi_eval_pysb_stmts_updated.pkl')

#orig_stmts = ac.load_statements('../../fallahi_eval_pysb_stmts.pkl')

erk_model_st = ac.filter_gene_list(stmts,['EGF','EGFR','SOS1','GRB2','KRAS','BRAF','MAPK1','MAPK3','MAP2K1'],'all')
erk_model_st_reduced = ac.filter_by_type(erk_model_st,Translocation,invert=True)

#things to add:
#'DUSP1',
#expand1 = ac.filter_gene_list(stmts,['BRAF','MAPK1','MAPK3','MAP2K1'],'one')
#nodes = filter_unique_nodes(expand1)
#nodes_to_use = ['DUSP1','JUN','FOS','PTPN11','PTPN6','CDK1','STAT3','ELK1']
#expand2 = ac.filter_gene_list(stmts,nodes_to_use,'one')
#nodes2 = filter_unique_nodes(expand2)
#nodes_to_use2 = ['ATF3','CDK1','CDK4','CCNA2','TP53','MYC']

#Remove ATF2 for now. 
#Also p38 (mapk14), SRC, SPRY4.
#SPRY2 included but doesn't do anything
master_nodes = ['EGF','EGFR','SOS1','GRB2','KRAS','BRAF','MAPK1','MAPK3','MAP2K1'] + ['JUN','DUSP1','MYB']#['DUSP1','JUN','FOS','PTPN6','ELK1','ATF2'] + ['CDK1','MYC']
erk_model_st = ac.filter_gene_list(stmts,master_nodes,'all')
erk_model_st_reduced = ac.filter_by_type(erk_model_st,Translocation,invert=True)


sentences = 'MAPK1 phosphorylates MYB. SPRY2 inhibits KRAS. Phosphorylated MYB transcribes SPRY2.'
trips_processor = trips.process_text(sentences)
manual_stmts = trips_processor.statements




def assembly_stuff(stmts1):
    #more stuff
#    stmts = ac.filter_direct(stmts1)
#    stmts = ac.filter_genes_only(stmts)
    stmts = ac.map_grounding(stmts1)
#    stmts = ac.filter_grounded_only(stmts)
#    stmts = ac.filter_human_only(stmts)
#    stmts = ac.expand_families(stmts)    
    stmts = ac.map_sequence(stmts)
    stmts = ac.filter_enzyme_kinase(stmts)
    stmts = ac.filter_mod_nokinase(stmts)
#    stmts = ac.filter_transcription_factor(stmts)
    # Run preassembly and save result
    stmts = ac.run_preassembly(stmts, return_toplevel=False)

    #ML Stuff
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.reduce_activities()
    ml.replace_activations()
    af_stmts = ac.filter_by_type(ml.statements, ActiveForm)
    non_af_stmts = ac.filter_by_type(ml.statements, indra.statements.ActiveForm, invert=True)
    af_stmts = ac.run_preassembly(af_stmts)
    stmts = af_stmts + non_af_stmts
    # Replace activations when possible


#    stmts = ac.filter_belief(stmts, 0.7)
    #stmts = ac.filter_top_level(stmts)
    return stmts 


#manual_stmts = assembly_stuff(manual_stmts)
erk_model_st_reduced = erk_model_st_reduced + manual_stmts




#New Steps
def assemble_new(stmts1):
    stmts = remove_dimers(stmts1)
#    stmts = remove_mutations(stmts)
    stmts = coarse_grain_phos_on_af(stmts)
    stmts = remove_redundant_phosphorylations(stmts)
    stmts = coarse_grain_kinase_context(stmts)
    stmts = coarse_grain_substrate_context(stmts)


    #new version
    stmts, down_list, up_list = add_all_af(stmts)
    stmts = reduce_complex_activeforms(stmts)
    stmts = combine_multiple_phos_activeforms(stmts)

    stmts = run_mechlinker_step(stmts,down_list,up_list)

    stmts = coarse_grain_phos_on_af(stmts)
    stmts = remove_redundant_phosphorylations(stmts)
    stmts = coarse_grain_kinase_context(stmts)
    stmts = coarse_grain_substrate_context(stmts)
    stmts = Preassembler.combine_duplicate_stmts(stmts)

    #Extra new
#    print(len(stmts))
    stmts = add_dephosphorylations(stmts)
#    print(len(stmts))
    return stmts 


#model_stmts_one = assemble_new(stmts_first_assembly_one)
#model_stmts_all = assemble_new(erk_model_st_reduced)

#def assembly_stuff(stmts1):
#    #more stuff
#    stmts = ac.filter_direct(stmts1)
#    stmts = ac.filter_genes_only(stmts)
#    stmts = ac.map_grounding(stmts)
#    stmts = ac.filter_grounded_only(stmts)
#    stmts = ac.filter_human_only(stmts)
#    stmts = ac.expand_families(stmts)    
#    stmts = ac.map_sequence(stmts)
#    stmts = ac.filter_enzyme_kinase(stmts)
#    stmts = ac.filter_mod_nokinase(stmts)
#    stmts = ac.filter_inconsequential_mods(stmts)
#    stmts = ac.filter_inconsequential_acts(stmts)
#    stmts = ac.filter_transcription_factor(stmts)
#     #Run preassembly and save result
#    stmts = ac.run_preassembly(stmts, return_toplevel=False)

#    #ML Stuff
#    ml = MechLinker(stmts)
#    ml.gather_explicit_activities()
#    ml.reduce_activities()
#    af_stmts = ac.filter_by_type(ml.statements, ActiveForm)
#    non_af_stmts = ac.filter_by_type(ml.statements, indra.statements.ActiveForm, invert=True)
#    af_stmts = ac.run_preassembly(af_stmts)
#    stmts = af_stmts + non_af_stmts
#    # Replace activations when possible
#    ml = MechLinker(stmts)
#    ml.gather_explicit_activities()
#    ml.replace_activations()
#    stmts = ml.statements
#    stmts = ac.filter_belief(stmts, 0.7)
#    #stmts = ac.filter_top_level(stmts)
#    return stmts 




gtp_stmts = ac.filter_by_type(stmts,GtpActivation)
kras_braf_gtp_stmts = ac.filter_gene_list(gtp_stmts,['BRAF','KRAS'],'all')
kras_braf_to_keep = deepcopy(kras_braf_gtp_stmts[0])
kras_braf_to_keep.agent_list()[0].mutations = []
#Out[154]: GtpActivation(KRAS(gtpbound: True), BRAF(), kinase)

#Double check Gef stmt is retained, and context on sos1 
#Gef(SOS1, KRAS)
#kras_braf_to_keep

#together give sos1->kras->braf 
#Need to incorporate sos AF - add gef to cascade context
#Need to carry down braf AF 

#May need manual intervention on BRAF. No obvious way to incorporate kinase~active without digging into pysb assembler.
#add a rule to all braf phosphorylations that allows it to happen when kinase~active
#Manual sos change as well, add its grb2 binding context onto gef

#Handle braf inhibition
#At some point go back and check what activating phosphorylations there were 





#Can drop:
#Generic kras activates braf
#egfr phos 
#egfr ubiq



#cycle changes to be made: 
#braf-mapk1-mapk3-map2k1
#ImplicitPhos

#Maybe add spry2 related stmts 
#MYB transcribes SPRY2
#mapk1 activates myb

#sentences = 'MAPK1 phosphorylates and activates MYB. SPRY2 inhibits KRAS. MYB transcribes SPRY2.'
#trips_processor = trips.process_text(sentences)
#trips_processor.statements

#erk_model_st_reduced = erk_model_st_reduced + trips_processor.statements




model_stmts_all = assemble_new(erk_model_st_reduced)
#model_stmts_all = assembly_stuff(model_stmts_all1)

model_stmts_all = model_stmts_all[2:]
#Removing artificial kras-braf activation 
#add v600e back in 
raw_braf_stmts = ac.filter_gene_list(erk_model_st_reduced,['BRAF'],'one')
braf_v600e_stmts = [deepcopy(raw_braf_stmts[18]),deepcopy(raw_braf_stmts[34]),deepcopy(raw_braf_stmts[42])]
for st in braf_v600e_stmts:
    st.residue=None
    st.position=None

model_stmts_all = model_stmts_all + braf_v600e_stmts



pa = PysbAssembler('two_step')
pa.add_statements(model_stmts_all)
model2 = pa.make_model()
 
bngl_model_filter = pysb.export.export(model2,'bngl')
bngl_file_filter = open('rbm/dusp_spry_feedback_v600e.bngl','w')
bngl_file_filter.write(bngl_model_filter)
bngl_file_filter.close()




#Runs quick, check dynamics 
#get gml and cycle check 
#import networkx as nx
#g = nx.read_gml('improved_reading_model_reduced_contactmap.gml')
#nx.cycle_basis(g)     





#ac.dump_statements(model_stmts_all,'model_stmts/model_stmts_all_new.pkl')

#Really weird bug, when TF has activity and bound condition it seems to cause PysbAssembler to error
#Need to investigate further

#edited_model_stmts_all = deepcopy(model_stmts_all)
#for st in edited_model_stmts_all:
#    if st.agent_list()[0].activity:
#        if st.agent_list()[0].activity.activity_type=='transcription':
#            st.agent_list()[0].bound_conditions = []





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






#pa = PysbAssembler('two_step')
#pa.add_statements(model_stmts_all)
#model1 = pa.make_model()

#bngl_model_filter = pysb.export.export(model1,'bngl')
#bngl_file_filter = open('rbm/secondpass_v1.bngl','w')
#bngl_file_filter.write(bngl_model_filter)
#bngl_file_filter.close()

#with open('rbm/final_version/erk_model_final2.txt','w') as f:
#    for item in erk_model_st_reduced:
#        f.write("%s\n" % item)

#ac.dump_statements(erk_model_st_reduced,'rbm/final_version/erk_model_final_stmts.pkl')

#remove hspb1, akt3, 
#rptor, grb2, kras not connectected to anything
