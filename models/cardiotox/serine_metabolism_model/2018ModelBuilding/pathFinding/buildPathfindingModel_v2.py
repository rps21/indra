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


def normalize_active_forms(stmts):
    af_stmts = ac.filter_by_type(stmts, ActiveForm)
    relevant_af_stmts = []
    for stmt in af_stmts:
        if (not stmt.agent.mods) and (not stmt.agent.mutations):
            continue
        relevant_af_stmts.append(stmt)
    print('%d relevant ActiveForms' % len(relevant_af_stmts))
    non_af_stmts = ac.filter_by_type(stmts, ActiveForm, invert=True)
    af_stmts = ac.run_preassembly(relevant_af_stmts)
    stmts = af_stmts + non_af_stmts
    return stmts

model_genes = [
    'FLT3','STAT3', 'JAK2', 'PKM', 'HIF1A', 'MYC',
    'CDKN1A','KDR','VEGFR1','FLT3','FLT1','PDGF','FLT3LG','PDGFA','PDGFRA',
    'JAK2','FLT3','ABL1','PTPN6','RPS6','RPS6KB1','PIK3CA','AKT1','AURKA'
    ]

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
    return biopax_stmts

biopax_stmts = build_prior(model_genes)
ac.dump_statements(biopax_stmts,'model_stmts/biopax_output.pkl')
trips_stmts = ac.load_statements('model_stmts/trips_output.pkl')

sp = signor.SignorProcessor()
signorStmts = sp.statements
#filteredSignorStmts = ac.filter_gene_list(signorStmts, model_genes, policy='one',save='signor_stmts_list.pkl')
filteredSignorStmts = ac.load_statements('model_stmts/signor_stmts_list.pkl')


#raw_stmts = ac.load_statements('heart_reading_stmts.pkl')
#stmts = ac.filter_no_hypothesis(raw_stmts)
#stmts = ac.filter_gene_list(stmts, model_genes, policy='one',
#                            save='filter_gene_list.pkl')

reach_stmts = ac.load_statements('model_stmts/filter_gene_list.pkl')
stmts = reach_stmts + biopax_stmts + trips_stmts + filteredSignorStmts

def cleanStatements(stmts):
    stmts = ac.map_grounding(stmts)
    stmts = ac.filter_grounded_only(stmts, save='intermediateStmts/grounded.pkl')
    stmts = ac.filter_genes_only(stmts)
    stmts = ac.filter_human_only(stmts)
    stmts = ac.run_preassembly(stmts, save='intermediateStmts/preassembled.pkl')
    stmts = ac.filter_belief(stmts, 0.90)
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
    stmts = normalize_active_forms(ml.statements)
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


largeModelStmts = cleanStatements(stmts)
smallModelRawStmts = ac.filter_gene_list(stmts,model_genes,'all')
smallModelStmts = cleanStatements(smallModelRawStmts)




#Coarse grain phos and add binding context, interpretting as activationsfrom indra.tools.small_model_tools import enforceCascadeContext as cs
from indra.tools.small_model_tools import enforceCascadeContext as cs
from indra.tools.small_model_tools import combinePhosphorylationSites as ptm

#newstmts, uplist, downlist = cs.add_all_af(largeModelStmts)     
#newstmts = cs.run_mechlinker_step_reduced(newstmts, uplist, downlist)
#largeModelStmts_contextChanges = ptm.coarse_grain_phos(newstmts)

newstmts, uplist, downlist = cs.add_all_af(smallModelStmts)     
smallModelStmts_contextChanges = cs.run_mechlinker_step_reduced(newstmts, uplist, downlist)
smallModelStmts_contextChanges = ptm.coarse_grain_phos(newstmts)


#Add statements for Sorafenib Targets
drugTargetStmts = ac.load_statements('sorafenibTargetStmts.pkl')
largeModelStmts = largeModelStmts + drugTargetStmts
#largeModelStmts_contextChanges = largeModelStmts_contextChanges + drugTargetStmts
smallModelStmts = smallModelStmts + drugTargetStmts
smallModelStmts_contextChanges = smallModelStmts_contextChanges + drugTargetStmts


largeModelStmts = Preassembler.combine_duplicate_stmts(largeModelStmts)
#largeModelStmts_contextChanges = Preassembler.combine_duplicate_stmts(largeModelStmts_contextChanges)
#largeModelStmts_contextChanges = ac.filter_inconsequential_mods(largeModelStmts_contextChanges)

smallModelStmts = Preassembler.combine_duplicate_stmts(smallModelStmts)
smallModelStmts_contextChanges = Preassembler.combine_duplicate_stmts(smallModelStmts_contextChanges)
smallModelStmts_contextChanges = ac.filter_inconsequential_mods(smallModelStmts_contextChanges)



#Save four final statment lists
ac.dump_statements(largeModelStmts,'finalLargeModelStmts.pkl')
ac.dump_statements(smallModelStmts,'finalSmallModelStmts.pkl')
#ac.dump_statements(largeModelStmts,'largeModelStmts_contextChanges.pkl')
ac.dump_statements(smallModelStmts_contextChanges,'smallModelStmts_contextChanges.pkl')

#Assemble and save four PySB models
pa = PysbAssembler()
pa.add_statements(largeModelStmts)
largeModel = pa.make_model()
with open('largePYSBModel.pkl','wb') as f:
    pickle.dump(largeModel,f)

pa = PysbAssembler()
pa.add_statements(smallModelStmts)
smallModel = pa.make_model()
with open('smallPYSBModel.pkl','wb') as f:
    pickle.dump(smallModel,f)

#pa = PysbAssembler()
#pa.add_statements(largeModelStmts_contextChanges)
#largeModel = pa.make_model()
#with open('largePYSBModel_contextChanges.pkl','wb') as f:
#    pickle.dump(largeModel,f)

pa = PysbAssembler()
pa.add_statements(smallModelStmts_contextChanges)
smallModelContext = pa.make_model()
with open('smallPYSBModel_contextChanges.pkl','wb') as f:
    pickle.dump(smallModelContext,f)

