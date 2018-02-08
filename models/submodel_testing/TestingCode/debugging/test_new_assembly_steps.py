import indra
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
from indra.preassembler.hierarchy_manager import hierarchies
from indra.tools import assemble_corpus as ac
from indra.statements import *
from indra.mechlinker import MechLinker
from copy import deepcopy
from indra.preassembler import Preassembler
from indra.assemblers import PysbAssembler
import pysb
from indra.sources import trips

###################################
############ IMPORTANT ############
########### NEW IMPORTS ###########
###################################





stmts = ac.load_statements('/home/bobby/Dropbox/Sorger_Lab/indra/models/fallahi_eval/dynamical_model/fallahi_eval_pysb_stmts_updated.pkl')

master_nodes = ['EGF','EGFR','SOS1','GRB2','KRAS','BRAF','MAPK1','MAPK3','MAP2K1'] + ['JUN','DUSP1','MYB']#['DUSP1','JUN','FOS','PTPN6','ELK1','ATF2'] + ['CDK1','MYC']
erk_model_st = ac.filter_gene_list(stmts,master_nodes,'all')
erk_model_st_reduced = ac.filter_by_type(erk_model_st,Translocation,invert=True)

sentences = 'MAPK1 phosphorylates MYB. SPRY2 inhibits KRAS. Phosphorylated MYB transcribes SPRY2.'
trips_processor = trips.process_text(sentences)
manual_stmts = trips_processor.statements

sentences = ['Phosphorylated AKT1 phosphorylates and activates RAF1. Phosphorylated RAF1 phosphorylates and activates JNK. Phosphorylated JNK phosphorylates and activates JUN on Serine.']
trips_processor = trips.process_text(sentences)
manual_stmts2 = trips_processor.statements


base_statements = erk_model_st_reduced + manual_stmts + manual_stmts2


def phos_simplification(stmts1):
    stmts = coarse_grain_phos_on_af(stmts)
    stmts = remove_redundant_phosphorylations(stmts)
    stmts = coarse_grain_kinase_context(stmts)
    stmts = coarse_grain_substrate_context(stmts)
    stmts = Preassembler.combine_duplicate_stmts(stmts)


def context_simplification(stmts1):
    #new version
    stmts, down_list, up_list = add_all_af(stmts)
    stmts = reduce_complex_activeforms(stmts)
    stmts = combine_multiple_phos_activeforms(stmts)
    stmts = run_mechlinker_step(stmts,down_list,up_list)


phos_simp_stmts = phos_simplification(base_statements)
context_simp_stmts = context_simplification(base_statements)



#rerun phos on new statements?
    #does it make a difference if we run before/after/both?
#dimer/mutations
#phosphatase
#gml - add something around it 


#New Steps
def assemble_new(stmts1):
    stmts = remove_dimers(stmts1)

    stmts = add_dephosphorylations(stmts)
#    print(len(stmts))
    return stmts 






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


old_assembly = assembly_stuff(erk_model_st_reduced)

gtp_stmts = ac.filter_by_type(stmts,GtpActivation)
kras_braf_gtp_stmts = ac.filter_gene_list(gtp_stmts,['BRAF','KRAS'],'all')
kras_braf_to_keep = deepcopy(kras_braf_gtp_stmts[0])
kras_braf_to_keep.agent_list()[0].mutations = []
#Out[154]: GtpActivation(KRAS(gtpbound: True), BRAF(), kinase)
erk_model_st_reduced.append(kras_braf_to_keep)


#erk_model_st_reduced = erk_model_st_reduced + trips_processor.statements




model_stmts_all = assemble_new(erk_model_st_reduced)
model_stmts_all = model_stmts_all[1:]


model_stmts_all = model_stmts_all[2:]
#Removing artificial kras-braf activation 
#add v600e back in 
raw_braf_stmts = ac.filter_gene_list(erk_model_st_reduced,['BRAF'],'one')
braf_v600e_stmts = [deepcopy(raw_braf_stmts[18]),deepcopy(raw_braf_stmts[34]),deepcopy(raw_braf_stmts[42])]
for st in braf_v600e_stmts:
    st.residue=None
    st.position=None

model_stmts_all = model_stmts_all + braf_v600e_stmts

ac.dump_statements(model_stmts_all,'../../ResultsandWriteup/v2/FinalResults/Final_Statements.pkl')

pa = PysbAssembler('two_step')
pa.add_statements(model_stmts_all)
model2 = pa.make_model()
 
bngl_model_filter = pysb.export.export(model2,'bngl')
bngl_file_filter = open('rbm/new_final_model.bngl','w')
bngl_file_filter.write(bngl_model_filter)
bngl_file_filter.close()







