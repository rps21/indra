from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from copy import deepcopy
import itertools
import indra
from indra.statements import *      
from indra.mechlinker import MechLinker
from indra.sources import trips
from indra.tools.gene_network import GeneNetwork
from indra.sources import biopax
from indra.assemblers.pysb_assembler import PysbAssembler
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.tools import assemble_corpus as ac
from indra.statements import Phosphorylation, Agent, Evidence
import pysb
from pysb.export import export



def grouped_biopax_query(gene_names, database_filter, block_size=60):
    gene_blocks = [gene_names[i:i+block_size] for i in
                   range(0, len(gene_names), block_size)]
    stmts = []
    # Run pathsfromto between pairs of blocks and pathsbetween
    # within each block. This breaks up a single call with N genes into
    # (N/block_size)*(N/blocksize) calls with block_size genes
    for genes1, genes2 in itertools.product(gene_blocks, repeat=2):
        if genes1 == genes2:
            bp = biopax.process_pc_pathsbetween(genes1,
                                            database_filter=database_filter)
        else:
            bp = biopax.process_pc_pathsfromto(genes1, genes2,
                                           database_filter=database_filter)
        if bp:
            stmts += bp.statements
    return stmts


def build_prior(gene_names):
    """Build a corpus of prior Statements from PC and BEL."""
    gn = GeneNetwork(gene_names)
    #bel_stmts = gn.get_bel_stmts(filter=False)
    #ac.dump_statements(bel_stmts, prefixed_pkl('bel'))
    # This call results in 504 error currently
    #biopax_stmts = gn.get_biopax_stmts(filter=False)
    database_filter = ['reactome', 'psp', 'kegg', 'pid']
    biopax_stmts = grouped_biopax_query(gene_names, database_filter) #This requires other function
#    ac.dump_statements(biopax_stmts, prefixed_pkl('biopax'))
#    ac.dump_statements(biopax_stmts, prefixed_pkl('biopax'))
    return biopax_stmts

#################
#Building model 
#################


#model_genes = ['FLT3', 'FLT3LG', 'STAT3', 'JAK2', 'PKM', 'HIF1A', 'MYC', 'CDKN1A', 'EPAS1']
#soraf_targets = ['VEGFR','KDR','VEGFR1','PDGFR','BRAF','FLT3','VEGF','PDGF','FLT3LG','PDGFA','PDGFRA','VEGFC']
model_genes = ['FLT3', 'FLT3LG', 'STAT3', 'JAK2', 'PKM', 'HIF1A', 'MYC', 'CDKN1A']
soraf_targets = ['KDR','VEGFR1','FLT3','FLT1','PDGF','FLT3LG','PDGFA','PDGFRA','VEGFC']
pkm_genes = ['JAK2','FLT3','ABL1','AGL1','PTPN6','SHP1']
#lapat_targets = ['EGFR','ERBB2','EGF']
extra_genes = ['RPS6','RPS6KB1','PIK3CA','AKT1','AURKA']


model_genes = model_genes+soraf_targets+pkm_genes
model_genes = model_genes+extra_genes

#biopax_stmts_new = build_prior(model_genes)
#ac.dump_statements(biopax_stmts,'model_stmts/biopax_output_new.pkl')


#In [25]: len(statementsAll)
#Out[25]: 135

#In [26]: len(statementsOne)
#Out[26]: 2600






#Identify and separate reach and trips pkls
reach_stmts = ac.load_statements('model_stmts/reach_output.pkl')
trips_stmts = ac.load_statements('model_stmts/trips_output.pkl')
biopax_stmts = ac.load_statements('model_stmts/biopax_output_new.pkl')

raw_stmts = reach_stmts + trips_stmts + biopax_stmts


filtered_stmts_one = ac.filter_gene_list(raw_stmts,model_genes,'one')
filtered_stmts_all = ac.filter_gene_list(raw_stmts,model_genes,'all')


##ac.dump_statements(filtered_stmts_one,'model_stmts/filtered_stmts_one.pkl')
#ac.dump_statements(filtered_stmts_all,'model_stmts/filtered_stmts_all_new.pkl')



#replace agl1 with PTPN6
for st in filtered_stmts_one:
    for ag in st.agent_list():
        if ag:
            if ag.name == 'PTPN6':
                ptpn_ag = deepcopy(ag)
                break

for st in filtered_stmts_all:
    for ag in st.agent_list():
        if ag:
            if ag.name == 'PTPN6':
                ptpn_ag = deepcopy(ag)
                break


#for st in ac.filter_by_type(filtered_stmts_one,Dephosphorylation):
to_remove = []
for st in ac.filter_gene_list(filtered_stmts_one,['AGL1'],'one'):
    for ag in st.agent_list():
        if ag:
            if ag.name == 'AGL1':
                ag_index = st.agent_list().index(ag)
                newst = deepcopy(st)
                newst_ag = deepcopy(newst.agent_list())
                newst_ag[ag_index] = deepcopy(ptpn_ag)
                newst.set_agent_list(newst_ag)
                filtered_stmts_one.append(newst)
                to_remove.append(st)
filtered_stmts_one = [st for st in filtered_stmts_one if st not in to_remove]


#for st in ac.filter_by_type(filtered_stmts_all,Dephosphorylation):
to_remove = []
for st in ac.filter_gene_list(filtered_stmts_all,['AGL1'],'one'):
    for ag in st.agent_list():
        if ag:
            if ag.name == 'AGL1':
                ag_index = st.agent_list().index(ag)
                newst = deepcopy(st)
                newst_ag = deepcopy(newst.agent_list())
                newst_ag[ag_index] = deepcopy(ptpn_ag)
                newst.set_agent_list(newst_ag)
                filtered_stmts_one.append(newst)
                to_remove.append(st)
filtered_stmts_one = [st for st in filtered_stmts_one if st not in to_remove]

def assembly_stuff(stmts1):
    #more stuff
    stmts = ac.filter_direct(stmts1)
    stmts = ac.filter_genes_only(stmts)
    stmts = ac.map_grounding(stmts)
    stmts = ac.filter_grounded_only(stmts)
    stmts = ac.filter_human_only(stmts)
    stmts = ac.expand_families(stmts)    
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
    af_stmts = ac.filter_by_type(ml.statements, ActiveForm)
    non_af_stmts = ac.filter_by_type(ml.statements, indra.statements.ActiveForm, invert=True)
    af_stmts = ac.run_preassembly(af_stmts)
    stmts = af_stmts + non_af_stmts
    # Replace activations when possible
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.replace_activations()
    stmts = ml.statements
    stmts = ac.filter_belief(stmts, 0.7)
    #stmts = ac.filter_top_level(stmts)
    return stmts 


stmts_first_assembly_one = assembly_stuff(filtered_stmts_one)
stmts_first_assembly_all = assembly_stuff(filtered_stmts_all)

#Issue with transcription with None as TF
#Probably need ac.filter_transcription_factor(stmts)
to_remove = []
for st in filtered_stmts_one:
     if st.agent_list()[0] == None:
         to_remove.append(st)    

for st in to_remove:
     filtered_stmts_one.remove(st)
     
statementsOne = [] 
threeag = []
for st in filtered_stmts_one:
    if len(st.agent_list()) < 3:
        statementsOne.append(st)
    else:
        threeag.append(st)

to_remove = []
for st in filtered_stmts_all:
     if st.agent_list()[0] == None:
         to_remove.append(st)    

for st in to_remove:
     filtered_stmts_all.remove(st)
     
statementsAll = [] 
threeag = []
for st in filtered_stmts_all:
    if len(st.agent_list()) < 3:
        statementsAll.append(st)
    else:
        threeag.append(st)

##Moving filtering below first assembly, experiment with running trips on all stmts above. maybe need an intermediate filtering first
##filtered_stmts_one = ac.filter_gene_list(raw_stmts,model_genes,'one')
#filtered_stmts_all = ac.filter_gene_list(raw_stmts,model_genes,'all')
##ac.dump_statements(filtered_stmts_one,'model_stmts/filtered_stmts_one.pkl')
#ac.dump_statements(filtered_stmts_all,'model_stmts/filtered_stmts_all_new.pkl')

#ac.dump_statements(stmts_first_assembly_one,'model_stmts/stmts_first_assembly_one.pkl')
#ac.dump_statements(stmts_first_assembly_all,'model_stmts/stmts_first_assembly_all_new.pkl')


from indra.sources import trips
manualSentences = 'FLT3 binds JAK2. PDGFRA binds JAK2'
trips_processor = trips.process_text(manualSentences)
statementsAll = statementsAll + trips_processor.statements 
statementsOne = statementsOne + trips_processor.statements 



soraf_targets = ['KDR','VEGFR1','FLT3','PDGFRA']
sorafSentences = 'Sorafenib inhibits KDR. Sorafenib inhibits VEGFR1. Sorafenib inhibits FLT3. Sorafenib inhibits PDGFRA.'
trips_processor = trips.process_text(sorafSentences)
soraf_stmts = trips_processor.statements 


experimentalSentences = 'Sorafenib phosphorylates cJun. Sorafenib phosphorylates STAT1. Sorafenib phosphorylates PKM2. Sorafenib dephosphorylates RPS6. Sorafenib phosphorylates Aurora kinase A. Sorafenib transcribes HIF1A. Sorafenib transcribes cMyc. Sorafenib transcribes PDGFRA.'  
trips_processor = trips.process_text(experimentalSentences)
exp_stmts = trips_processor.statements 


statementsAll = statementsAll + soraf_stmts + exp_stmts
statementsOne = statementsOne + soraf_stmts + exp_stmts


##Assemble model without new assembly steps
pa = PysbAssembler('one_step')
pa.add_statements(statementsOne)
largeModel = pa.make_model()

bngl_model_filter = pysb.export.export(largeModel,'bngl')
bngl_file_filter = open('largeModel.bngl','w')
bngl_file_filter.write(bngl_model_filter)
bngl_file_filter.close()

pa = PysbAssembler('one_step')
pa.add_statements(statementsAll)
smallModel = pa.make_model()

bngl_model_filter = pysb.export.export(smallModel,'bngl')
bngl_file_filter = open('smallModel.bngl','w')
bngl_file_filter.write(bngl_model_filter)
bngl_file_filter.close()

#Pickle stmts, pysb model, bngl model for both 
ac.dump_statements(statementsAll,'smallStatementSet.pkl')
ac.dump_statements(statementsOne,'largeStatementSet.pkl')

import pickle
with open('largePYSBModel.pkl','wb') as f:
    pickle.dump(largeModel,f)

with open('smallPYSBModel.pkl','wb') as f:
    pickle.dump(smallModel,f)



##New Steps
#def assemble_new(stmts1):
#    stmts = remove_dimers(stmts1)
##    stmts = remove_mutations(stmts)
#    stmts = coarse_grain_phos_on_af(stmts)
#    stmts = remove_redundant_phosphorylations(stmts)
#    stmts = coarse_grain_kinase_context(stmts)
#    stmts = coarse_grain_substrate_context(stmts)


#    #new version
#    stmts, down_list, up_list = add_all_af(stmts)
#    stmts = reduce_complex_activeforms(stmts)
#    stmts = combine_multiple_phos_activeforms(stmts)

#    stmts = run_mechlinker_step(stmts,down_list,up_list)

#    stmts = coarse_grain_phos_on_af(stmts)
#    stmts = remove_redundant_phosphorylations(stmts)
#    stmts = coarse_grain_kinase_context(stmts)
#    stmts = coarse_grain_substrate_context(stmts)
#    stmts = Preassembler.combine_duplicate_stmts(stmts)

#    return stmts 


#model_stmts_one = assemble_new(stmts_first_assembly_one)
#model_stmts_all = assemble_new(stmts_first_assembly_all)

##ac.dump_statements(model_stmts_one,'model_stmts/model_stmts_one.pkl')
#ac.dump_statements(model_stmts_all,'model_stmts/model_stmts_all_new.pkl')



##Really weird bug, when TF has activity and bound condition it seems to cause PysbAssembler to error
##Need to investigate further

#edited_model_stmts_all = deepcopy(model_stmts_all)
#for st in edited_model_stmts_all:
#    if st.agent_list()[0].activity:
#        if st.agent_list()[0].activity.activity_type=='transcription':
#            st.agent_list()[0].bound_conditions = []



#model_stmts_all = ac.load_statements('model_stmts/model_stmts_all.pkl')
#modelStmts = ac.filter_gene_list(raw_stmts,modelGenes,'all')

#pa = PysbAssembler('two_step')
#pa.add_statements(edited_model_stmts_all)
#model2 = pa.make_model()

#bngl_model_filter = pysb.export.export(model2,'bngl')
#bngl_file_filter = open('multi_target_model_v2.bngl','w')
#bngl_file_filter.write(bngl_model_filter)
#bngl_file_filter.close()


