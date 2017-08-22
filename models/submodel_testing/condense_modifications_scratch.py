from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
import re
import os
import pandas as pd
from indra.databases import relevance_client, hgnc_client
from indra.assemblers import PysbAssembler
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.tools import assemble_corpus as ac
from indra.statements import Phosphorylation, Agent, Evidence
import pysb
from pysb.export import export
import indra
from indra.statements import *      
from indra.mechlinker import MechLinker
from indra.tools import assemble_corpus as ac
from indra.statements import *   

stmts = my_model_stmts_larger
stmts = coarse_grain_phos_on_af(stmts)
stmts = remove_redundant_phosphorylations(stmts)
stmts = coarse_grain_kinase_context(stmts)

#pa = PysbAssembler('two_step')
#pa.add_statements(stmts)
#model = pa.make_model()

#bngl_model_filter = pysb.export.export(model,'bngl')
#bngl_file_filter = open('rbm/egf_egfr_sos_grb_olig.bngl','w')
#bngl_file_filter.write(bngl_model_filter)
#bngl_file_filter.close()

#ac.dump_statements('egf_egfr_sos_grb_olig_stmts.pkl',stmts)

