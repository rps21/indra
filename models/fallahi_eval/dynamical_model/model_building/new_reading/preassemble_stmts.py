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

stmts1 = ac.load_statements('stmts.pkl')

stmts = ac.filter_direct(stmts1)
#stmts = ac.filter_belief(stmts, 0.95)
stmts = ac.filter_top_level(stmts)
stmts = ac.filter_enzyme_kinase(stmts)
stmts = ac.filter_mod_nokinase(stmts)
stmts = ac.filter_transcription_factor(stmts)
stmts = ac.filter_grounded_only(stmts)
stmts = ac.filter_genes_only(stmts)
stmts =ac.filter_human_only(stmts)
stmts = ac.run_preassembly(stmts)  #can split into deduplication and related

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

ac.dump_statements(stmts,'stmts_preassembled_fixed_be_filtering.pkl')



