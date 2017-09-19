from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
from indra.assemblers import PysbAssembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.tools import assemble_corpus as ac
from indra.statements import *
import pysb
from pysb.export import export
import indra


othergenes = ['EGFR','GRB2','SOS1','KRAS','SOS']
observablegenes = ['MAP2K1','MAPK1','MAPK3','RPS6KA1','AKT1','AKT3','MTOR','RPTOR','RICTOR','RPS6KB1','RPS6','PRKAA1','PRKAA2','PRKAB1','PRKAB2','PRKAG1','PRKAG2','PRKAG3','MAPK8','MAPK9','MAPK10','JUN','MAPK14','MAPK11','HSPB1','RELA','NFKB1','NFKB2','BCL2L11','PARP1','PARP3','HIST3H3','CDKN1B']
targetgenes = ['RAF1','BRAF','KDR','DDR2','LYN','FLT1','FLT3','CSF1R','MAP2K1']
reduced_targets = ['RAF1','BRAF']
allgenes = othergenes+targetgenes+observablegenes
allgenes_reduced = othergenes+reduced_targets+observablegenes

#to remove after first pass analysis
#for now, based on model alone and not biology or exp, mapk8 (JNK), mapk9 (JNK2), mapk14 (p38), NFKB1, RELA
allgenes_reduced.remove('MAPK8')
allgenes_reduced.remove('MAPK9')
allgenes_reduced.remove('MAPK14')
allgenes_reduced.remove('NFKB1')
allgenes_reduced.remove('RELA')


stmts = ac.load_statements('stmts_preassembled_fixed_be_filtering.pkl')

#all_stmts = ac.filter_gene_list(stmts, allgenes, 'all')
all_stmts = ac.filter_gene_list(stmts, allgenes_reduced, 'all')
#drug_stmts = ac.filter_gene_list(stmts, targetgenes, 'one')
drug_stmts = ac.filter_gene_list(stmts, reduced_targets, 'one')
exp_stmts = ac.filter_gene_list(stmts, observablegenes, 'one')



pa = PysbAssembler('two_step')
pa.add_statements(all_stmts)
model1 = pa.make_model()

bngl_model_filter = pysb.export.export(model1,'bngl')
bngl_file_filter = open('rbm/second_pass.bngl','w')
bngl_file_filter.write(bngl_model_filter)
bngl_file_filter.close()


##to remove
#grb2, sos1, nfkb1, rela, akt3, mapk9, prkaa1, parp3

##to investigate
#cdkn1b phos
#parp1/3/cleavage
