import pickle
from indra.mechlinker import MechLinker
from indra.statements import *
from indra.sources import signor
from indra.sources import biopax
from indra.sources import trips
from indra.tools.gene_network import GeneNetwork
from indra.tools import assemble_corpus as ac
from indra.preassembler import Preassembler
from indra.util import _require_python3



def add_database_stmts(model_genes):
    """Build a corpus of prior Statements from PC and BEL."""
    gn = GeneNetwork(model_genes)

    # Read BEL Statements
    bel_stmts = gn.get_bel_stmts(filter=False)
    # Read Pathway Commons Statements
    database_filter = ['reactome', 'kegg', 'pid']
    biopax_stmts = gn.get_biopax_stmts(database_filter=database_filter)
    # Eliminate blacklisted interactions
    tmp_stmts = []
    biopax_blacklist = \
        ['http://pathwaycommons.org/pc2/Catalysis_13953a072d6f992f2388d85f9059a475']
    for stmt in biopax_stmts:
        source_ids = [ev.source_id for ev in stmt.evidence]
        if set(source_ids) & set(biopax_blacklist):
            continue
        tmp_stmts.append(stmt)
    biopax_stmts = tmp_stmts
    # Read Signor Statements
    sp = signor.api.process_from_web()
    signorStmtsAll = sp.statements
    signor_stmts = ac.filter_gene_list(signorStmtsAll, model_genes, policy='one')#,save='signor_stmts_list.pkl')

    priorStmts = bel_stmts + biopax_stmts + signor_stmts
    return priorStmts


def build_raw_prior(model_genes,dbs=True,additional_stmts_files=[]):
    if dbs:
        db_stmts = add_database_stmts(model_genes)
    else:
        db_stmts = []
    additional_stmts = []
    for fn in additional_stmts_files:
        additional_stmts = ac.load_statements(fn) + additional_stmts  
    raw_prior_stmts = db_stmts + additional_stmts
    return raw_prior_stmts


def assemble_prior(stmts):
    stmts = ac.map_grounding(stmts)
    stmts = ac.filter_grounded_only(stmts)#, save='intermediateStmts/grounded.pkl')
    stmts = ac.filter_genes_only(stmts)
    stmts = ac.filter_human_only(stmts)
    stmts = ac.filter_direct(stmts)
    stmts = ac.run_preassembly(stmts)#, save='intermediateStmts/preassembled.pkl')
    stmts = ac.filter_belief(stmts, 0.50)
    stmts = ac.filter_top_level(stmts)
    stmts = ac.filter_enzyme_kinase(stmts)
    stmts = ac.filter_mod_nokinase(stmts)
    stmts = ac.filter_transcription_factor(stmts)   
    # Simplify activity types
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.reduce_activities()
    ml.gather_modifications()
    ml.reduce_modifications()
    # Replace activations when possible
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.replace_activations()
    # Require active forms
    ml.require_active_forms()

    stmts = ml.statements
    return stmts


def filter_prior(stmts,model_genes,model_types):
    filtered_by_genes = ac.filter_gene_list(stmts,model_genes,'all')
    filtered_stmts = []
    for ty in model_types:
        filtered_stmts = filtered_stmts + ac.filter_by_type(filtered_by_genes,ty)
    return filtered_stmts
    
def build_prior(model_genes,model_types,dbs=True,additional_stmts_files=[]):
    raw_model_stmts =  build_raw_prior(model_genes,dbs,additional_stmts_files)
    filtered_stmts = filter_prior(raw_model_stmts,model_genes,model_types)
    assembled_stmts = assemble_prior(filtered_stmts)
    return assembled_stmts









