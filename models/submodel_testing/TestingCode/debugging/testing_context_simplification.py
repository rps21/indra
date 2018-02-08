#from __future__ import absolute_import, print_function, unicode_literals
#from builtins import dict, str
#import pickle
#from indra.preassembler.hierarchy_manager import hierarchies
from indra.statements import *
#from indra.mechlinker import MechLinker
from copy import deepcopy
#from indra.preassembler import Preassembler
#from indra.assemblers import PysbAssembler
#import pysb
from indra.sources import trips
from indra.tools import assemble_corpus as ac
#import pickle

import context_simplification as cs
import ptm_simplification as ptm

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


#KEY
#newstmts, uplist, downlist = cs.add_all_af(base_statements)     
#final_stmts = cs.run_mechlinker_step_reduced(newstmts, uplist, downlist)

##smallTest = ac.filter_gene_list(final_stmts,['BRAF','MAP2K1','MAPK1'],'all')
#stmts_to_use = deepcopy(final_stmts)
#newStmts = ptm.coarse_grain_phos(stmts_to_use)

#bigger testing
#KEY
newstmts, uplist, downlist = cs.add_all_af(stmts)     
final_stmts = cs.run_mechlinker_step_reduced(newstmts, uplist, downlist)

#smallTest = ac.filter_gene_list(final_stmts,['BRAF','MAP2K1','MAPK1'],'all')
stmts_to_use = deepcopy(final_stmts)
newStmts = ptm.coarse_grain_phos(stmts_to_use)




###########
#loop testing
#smallset = ac.filter_gene_list(stmts,['BRAF','MAP2K1','MAPK1'],'all')
#newstmts, uplist, downlist = cs.add_all_af(smallset)     
#final_stmts = cs.run_mechlinker_step_reduced(newstmts, uplist, downlist)
