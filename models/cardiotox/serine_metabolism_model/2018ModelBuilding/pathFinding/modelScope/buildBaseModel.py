import pickle
import itertools
from indra.tools import assemble_corpus as ac
from indra.assemblers import PysbAssembler
from indra.mechlinker import MechLinker
from indra.statements import *
from indra.sources import signor
from indra.sources import biopax
from indra.sources import trips
from indra.tools.gene_network import GeneNetwork
from indra.preassembler import Preassembler
from indra.util import _require_python3
#from util import prefixed_pkl, based, basen

#process_databases.py 
phosphosite_owl_file = 'sources/Kinase_substrates.owl'


#def read_phosphosite_owl(fname=phosphosite_owl_file):
#    bp = biopax.process_owl(fname)
#    for stmt in bp.statements:
#        for ev in stmt.evidence:
#            ev.source_api = 'phosphosite'
#            ev.epistemics = {'direct': True}
#    return bp.statements


def build_prior(gene_names):
    """Build a corpus of prior Statements from PC and BEL."""
    gn = GeneNetwork(gene_names)
    # Read BEL Statements
    bel_stmts = gn.get_bel_stmts(filter=False)
    #ac.dump_statements(bel_stmts, prefixed_pkl('bel'))
    # Read Pathway Commons Statements
    database_filter = ['reactome', 'kegg', 'pid']
    biopax_stmts = gn.get_biopax_stmts(database_filter=database_filter)
    # Eliminate blacklisted interactions
    tmp_stmts = []
    for stmt in biopax_stmts:
        source_ids = [ev.source_id for ev in stmt.evidence]
        if set(source_ids) & set(biopax_blacklist):
            continue
        tmp_stmts.append(stmt)
    biopax_stmts = tmp_stmts
    #ac.dump_statements(biopax_stmts, prefixed_pkl('biopax'))
    # Read Phosphosite Statements
   # phosphosite_stmts = read_phosphosite_owl(phosphosite_owl_file)
    #ac.dump_statements(phosphosite_stmts, prefixed_pkl('phosphosite'))

    priorStmts = bel_stmts + biopax_stmts 
#    priorStmts = biopax_stmts #+ phosphosite_stmts
    return priorStmts

biopax_blacklist = \
    ['http://pathwaycommons.org/pc2/Catalysis_13953a072d6f992f2388d85f9059a475']


def cleanStatements(stmts):
    stmts = ac.map_grounding(stmts)
    stmts = ac.filter_grounded_only(stmts, save='intermediateStmts/grounded.pkl')
    stmts = ac.filter_genes_only(stmts)
    stmts = ac.filter_human_only(stmts)
    stmts = ac.run_preassembly(stmts, save='intermediateStmts/preassembled.pkl')
    stmts = ac.filter_belief(stmts, 0.50)
    stmts = ac.filter_top_level(stmts)
    stmts = ac.filter_enzyme_kinase(stmts)
    stmts = ac.filter_mod_nokinase(stmts)
    stmts = ac.filter_transcription_factor(stmts)   #Any downside?
    # Simplify activity types
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.reduce_activities()
    ml.gather_modifications()
    ml.reduce_modifications()
#    stmts = normalize_active_forms(ml.statements)
    # Replace activations when possible
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.replace_activations()
    # Require active forms
    ml.require_active_forms()
    num_stmts = len(ml.statements)
    """
    while True:
        # Remove inconsequential PTMs
        ml.statements = ac.filter_inconsequential_mods(ml.statements,
                                                       get_mod_whitelist())
        ml.statements = ac.filter_inconsequential_acts(ml.statements,
                                                       get_mod_whitelist())
        if num_stmts <= len(ml.statements):
            break
        num_stmts = len(ml.statements)
    """
    stmts = ml.statements
    return stmts


model_genes = ['EGF','EGFR','GRB2','SOS1','KRAS','BRAF']

dbStmts = build_prior(model_genes)

#sp = signor.processor.SignorProcessor()
#signorStmtsAll = sp.statements
#signorStmts = ac.filter_gene_list(signorStmtsAll, model_genes, policy='one',save='signor_stmts_list.pkl')

#reach_stmts = ac.load_statements('model_stmts/filter_gene_list.pkl')
#trips_stmts = ac.load_statements('model_stmts/trips_output.pkl')


stmts = dbStmts #+ signorStmts  
filteredStmts = ac.filter_gene_list(stmts,model_genes,'all')
filteredStmts = Preassembler.combine_duplicate_stmts(filteredStmts) 


#AT THIS POINT:
#ac.dump_statements(filteredStmts,'testingModelBuilding.pkl')



##New filter, because of errors that arrise from none agents
#removeNone = [st for st in filteredStmts if None not in st.agent_list()]

# Gef(SOS1, KRAS),
# GtpActivation(KRAS(gtpbound), BRAF(), kinase),

stmtTypesToKeep = [Phosphorylation,Dephosphorylation,Activation,Inhibition,ActiveForm,IncreaseAmount,DecreaseAmount,Gef,GtpActivation]#Check sos and raf stmts 
finalStmts = []
for ty in stmtTypesToKeep:
    finalStmts = finalStmts + ac.filter_by_type(filteredStmts,ty)

ac.dump_statements(finalStmts,'testingModelBuilding.pkl')


#remove mutants, dimers 
#modelStmts = cleanStatements(removeNone)























