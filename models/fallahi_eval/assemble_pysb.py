from indra.statements import *
from indra.mechlinker import MechLinker
import indra.tools.assemble_corpus as ac
from indra.assemblers import PysbAssembler, IndexCardAssembler

def assemble_pysb(stmts, data_genes):
    # Filter the INDRA Statements to be put into the model
    stmts = ac.filter_by_type(stmts, Complex, invert=True)
    stmts = ac.filter_direct(stmts)
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)
    stmts = ac.filter_gene_list(stmts, data_genes, 'all')
    stmts = ac.filter_enzyme_kinase(stmts)
    stmts = ac.filter_mod_nokinase(stmts)
    stmts = ac.filter_transcription_factor(stmts)
    # Simplify activity types
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.reduce_activities()
    ml.gather_modifications()
    ml.reduce_modifications()
    af_stmts = ac.filter_by_type(ml.statements, ActiveForm)
    non_af_stmts = ac.filter_by_type(ml.statements, ActiveForm, invert=True)
    af_stmts = ac.run_preassembly(af_stmts)
    stmts = af_stmts + non_af_stmts
    # Replace activations when possible
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.replace_activations()
    # Require active forms
    ml.require_active_forms()
    num_stmts = len(ml.statements)
    while True:
        # Remove inconsequential PTMs
        ml.statements = ac.filter_inconsequential_mods(ml.statements,
                                                       get_mod_whitelist())
        ml.statements = ac.filter_inconsequential_acts(ml.statements,
                                                       get_mod_whitelist())
        if num_stmts <= len(ml.statements):
            break
        num_stmts = len(ml.statements)
    stmts = ml.statements

    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    return model

def contextualize_model(model, cell_line):
    # Here we just make a PysbAssembler to be able
    # to apply set_context on the model being passed in
    pa = PysbAssembler()
    pa.model = model
    pa.set_context(cell_line)

def get_mod_whitelist():
    # TODO: populate this with the actual antibody sites
    return {'MAP2K1': [('phosphorylation', 'S', '222')]}
