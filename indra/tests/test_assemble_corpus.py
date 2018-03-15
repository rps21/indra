from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
from indra.tools import assemble_corpus as ac
from indra.statements import *

a = Agent('a', db_refs={'HGNC': '1234', 'TEXT': 'a'})
b = Agent('b', db_refs={'UP': 'P15056', 'TEXT': 'b'})
c = Agent('c', db_refs={'FPLX': 'XXX', 'TEXT': 'c'})
d = Agent('d', db_refs={'TEXT': 'd'})
e = Agent('e', db_refs={'CHEBI': 'CHEBI:1234', 'TEXT': 'e'})
f = Agent('b', db_refs={'UP': 'P28028', 'TEXT': 'b'})
g = Agent('g', db_refs={'FPLX': 'ERK'})
h = Agent('g', mods=['x', 'y'], mutations=['x', 'y'], activity='x',
               location='nucleus', bound_conditions=['x', 'y', 'z'])
mapk1 = Agent('MAPK1', db_refs={'HGNC':'6871', 'UP':'P28482'})
erk = Agent('ERK', db_refs={'FPLX': 'ERK'})
st1 = Phosphorylation(a, b)
st2 = Phosphorylation(a, d)
st3 = Phosphorylation(c, d)
st4 = Phosphorylation(b, e)
st5 = Phosphorylation(None, b)
st6 = Phosphorylation(None, d)
st7 = Phosphorylation(None, e)
st8 = Phosphorylation(b, f)
st9 = Phosphorylation(None, f)
st10 = Phosphorylation(None, g)
st11 = Phosphorylation(None, h)
st12 = Phosphorylation(a, b, evidence=[Evidence(epistemics={'direct': True})])
st13 = Phosphorylation(a, b, evidence=[Evidence(epistemics={'direct': False})])
st14 = Activation(a, b, 'activity')
st15 = Activation(a, b, 'kinase')
st14.supports = [st15]
st15.supported_by = [st14]
st16 = Phosphorylation(a, mapk1)
st17 = Phosphorylation(a, erk)
st1.belief = 0.9
st2.belief = 0.8
st3.belief = 0.7

def test_load_stmts():
    with open('_test.pkl', 'wb') as fh:
        pickle.dump([st1], fh, protocol=2)
    st_loaded = ac.load_statements('_test.pkl')
    assert(len(st_loaded) == 1)
    assert(st_loaded[0].equals(st1))

def test_dump_stmts():
    ac.dump_statements([st1], '_test.pkl')
    st_loaded = ac.load_statements('_test.pkl')
    assert(len(st_loaded) == 1)
    assert(st_loaded[0].equals(st1))

def test_filter_grounded_only():
    st_out = ac.filter_grounded_only([st1, st4])
    assert len(st_out) == 2
    st_out = ac.filter_grounded_only([st3])
    assert len(st_out) == 0

def test_filter_uuid_list():
    st_out = ac.filter_uuid_list([st1, st4], [st1.uuid])
    assert len(st_out) == 1

def test_filter_genes_only():
    st_out = ac.filter_genes_only([st1, st5])
    assert len(st_out) == 2
    st_out = ac.filter_genes_only([st6, st7])
    assert len(st_out) == 0
    st_out = ac.filter_genes_only([st4])
    assert len(st_out) == 0
    st_out = ac.filter_genes_only([st3], specific_only=True)
    assert len(st_out) == 0

def test_filter_human_only():
    st_out = ac.filter_human_only([st1, st5])
    assert len(st_out) == 2
    st_out = ac.filter_human_only([st8, st9])
    assert len(st_out) == 0

def test_filter_gene_list_one():
    st_out = ac.filter_gene_list([st1, st2], ['a'], 'one')
    assert(len(st_out) == 2)
    st_out = ac.filter_gene_list([st1, st2], ['a'], 'all')
    assert(len(st_out) == 0)
    st_out = ac.filter_gene_list([st1, st2], ['a', 'b'], 'all')
    assert(len(st_out) == 1)
    st_out = ac.filter_gene_list([st1, st2], ['a', 'b'], 'invalid')
    assert(len(st_out) == 2)

def test_filter_gene_list_families():
    stmts_out = ac.filter_gene_list([st16, st17], ['MAPK1'], 'one',
                                    allow_families=False)
    assert len(stmts_out) == 1
    assert stmts_out[0] == st16
    stmts_out = ac.filter_gene_list([st16, st17], ['MAPK1'], 'one',
                                    allow_families=True)
    assert len(stmts_out) == 2
    assert st16 in stmts_out
    assert st17 in stmts_out

def test_run_preassembly():
    st_out = ac.run_preassembly([st1, st3, st5, st6])
    assert(len(st_out) == 2)

def test_run_preassembly_all_stmts():
    st_out = ac.run_preassembly([st1, st3, st5, st6], return_toplevel=False)
    assert(len(st_out) == 4)

def test_expand_families():
    st_out = ac.expand_families([st10])
    assert(len(st_out) == 2)

def test_strip_agent_context():
    st_out = ac.strip_agent_context([st11])
    assert(len(st_out) == 1)
    assert(not st_out[0].sub.mods)
    assert(not st_out[0].sub.mutations)
    assert(not st_out[0].sub.bound_conditions)
    assert(not st_out[0].sub.activity)
    assert(not st_out[0].sub.location)

def test_filter_direct():
    st_out = ac.filter_direct([st12])
    assert(len(st_out) == 1)
    st_out = ac.filter_direct([st13])
    assert(len(st_out) == 0)

def test_filter_belief():
    st_out = ac.filter_belief([st1, st2, st3], 0.75)
    assert(len(st_out) == 2)

def test_reduce_activities():
    st_out = ac.reduce_activities([st14, st15])
    assert(st_out[0].obj_activity == 'kinase')
    assert(st_out[1].obj_activity == 'kinase')

def test_filter_source():
    ev1 = Evidence(source_api='bel')
    ev2 = Evidence(source_api='biopax')
    ev3 = Evidence(source_api='reach')
    st1 = Activation(Agent('a'), Agent('b'), evidence=[ev3])
    st2  = Activation(Agent('a'), Agent('b'), evidence=[ev1, ev2])
    st3 = Activation(Agent('a'), Agent('b'), evidence=[ev1, ev3])
    st_out = ac.filter_evidence_source([st1, st2], ['reach'], 'one')
    assert(len(st_out) == 1)
    st_out = ac.filter_evidence_source([st1, st2, st3], ['reach'], 'all')
    assert (len(st_out) == 2)
    st_out = ac.filter_evidence_source([st1, st2, st3], ['bel', 'biopax'],
                                       'one')
    assert (len(st_out) == 2)
    st_out = ac.filter_evidence_source([st1, st2, st3], ['bel', 'biopax'],
                                       'all')
    assert (len(st_out) == 1)
    st_out = ac.filter_evidence_source([st1, st2, st3], ['bel', 'biopax'],
                                       'none')
    assert (len(st_out) == 1)

def test_map_grounding():
    a = Agent('MEK', db_refs={'TEXT': 'MEK'})
    b = Agent('X', db_refs={'TEXT': 'ERK'})
    st = Activation(a, b)
    st_out = ac.map_grounding([st], do_rename=False)
    assert(len(st_out) == 1)
    assert(st_out[0].subj.db_refs.get('FPLX'))
    assert(st_out[0].obj.db_refs.get('FPLX'))
    assert(st_out[0].obj.name == 'X')
    st_out = ac.map_grounding([st], do_rename=True)
    assert(len(st_out) == 1)
    assert(st_out[0].subj.db_refs.get('FPLX'))
    assert(st_out[0].obj.db_refs.get('FPLX'))
    assert(st_out[0].obj.name == 'ERK')

def test_map_sequence():
    a = Agent('MAPK1', db_refs={'UP': 'P28482', 'HGNC': '6871'})
    st1 = Phosphorylation(None, a, 'T', '182')
    st2 = Phosphorylation(None, a, 'T', '185')
    st3 = Phosphorylation(None, a, 'Y', '999')
    st_out = ac.map_sequence([st1])
    assert(len(st_out) == 1)
    assert(st_out[0].position == '185')
    st_out = ac.map_sequence([st2])
    assert(len(st_out) == 1)
    assert(st_out[0].position == '185')
    st_out = ac.map_sequence([st3])
    assert(len(st_out) == 0)

def test_filter_by_type():
    st_out = ac.filter_by_type([st1, st14], Phosphorylation)
    assert(len(st_out) == 1)

def test_filter_top_level():
    st_out = ac.filter_top_level([st14, st15])
    assert(len(st_out) == 1)

def test_filter_no_hypothesis():
    a = Agent('MAPK1')
    ev1 = Evidence(epistemics={'hypothesis': True})
    ev2 = Evidence(epistemics={'hypothesis': False})
    st1 = Phosphorylation(None, a, evidence=[ev1, ev2])
    st2 = Phosphorylation(None, a, evidence=[ev1, ev1])
    st_out = ac.filter_no_hypothesis([st1, st2])

def test_belief_cut_plus_filter_top():
    st1 = Phosphorylation(None, Agent('a'))
    st2 = Phosphorylation(Agent('b'), Agent('a'))
    st1.supports = [st2]
    st2.supported_by = [st1]
    st1.belief = 0.9
    st2.belief = 0.1
    st_high_belief = ac.filter_belief([st1, st2], 0.5)
    st_top_level = ac.filter_top_level(st_high_belief)
    assert(len(st_top_level) == 1)

def test_filter_inconsequential_mods():
    mc = ModCondition('phosphorylation', None, None, True)
    st1 = Phosphorylation(None, Agent('a'))
    st2 = Phosphorylation(Agent('a', mods=[mc]), Agent('b'))
    st_out = ac.filter_inconsequential_mods([st1, st2])
    assert(len(st_out) == 1)
    whitelist = {'b': [('phosphorylation', None, None)]}
    st_out = ac.filter_inconsequential_mods([st1, st2], whitelist=whitelist)
    assert(len(st_out) == 2)

def test_filter_inconsequential_mods2():
    st1 = Phosphorylation(Agent('a'), Agent('b'), 'S', '315')
    whitelist = {'b': [('phosphorylation', 'S', '315')]}
    st_out = ac.filter_inconsequential_mods([st1, st2], whitelist=whitelist)
    assert(len(st_out) == 1)

def test_filter_inconsequential_activities():
    st1 = Activation(Agent('a', activity=ActivityCondition('kinase', True)),
                     Agent('b'), 'activity')
    st2 = Activation(Agent('c'), Agent('a'), 'kinase')
    st_out = ac.filter_inconsequential_acts([st1, st2])
    assert(len(st_out) == 1)
    st_out = ac.filter_inconsequential_acts(st_out)
    assert(len(st_out) == 0)

def test_filter_mutation_status():
    braf_mut = Agent('BRAF', mutations=MutCondition('600', 'V', 'E'))
    braf_other_mut = Agent('BRAF', mutations=MutCondition('555', 'K', 'G'))
    st1 = Phosphorylation(braf_mut, Agent('a'))
    st2 = Phosphorylation(braf_other_mut, Agent('a'))
    mutations = {'BRAF': [('V', '600', 'E')]}
    deletions = []
    st_out = ac.filter_mutation_status([st1, st2], mutations, deletions)
    assert(len(st_out) == 1)
    mutations = {}
    deletions = ['a']
    st_out = ac.filter_mutation_status([st1, st2], mutations, deletions)
    assert(len(st_out) == 0)

def test_get_unreachable_mods():
    st1 = Phosphorylation(Agent('X'), Agent('Y'), 'S', '222')
    mcs = [ModCondition('phosphorylation', 'S', '218', True),
           ModCondition('phosphorylation', 'S', '222', True)]
    st2 = ActiveForm(Agent('Y', mods=mcs), 'activity', True)
    res = ac.get_unreachable_mods([st1, st2])
    assert 'Y' in res, res
    assert res['Y'] == set([('phosphorylation', 'S', '218')])


def test_rename_db_ref():
    x = Agent('X', db_refs={'BE': 'X'})
    y = Agent('Y', db_refs={'FPLX': 'Y'})
    st1 = Phosphorylation(x, y)
    stmts = ac.rename_db_ref([st1], 'BE', 'FPLX')
    assert len(stmts) == 1
    assert stmts[0].enz.db_refs.get('FPLX') == 'X'
    assert 'BE' not in stmts[0].enz.db_refs
    assert stmts[0].sub.db_refs.get('FPLX') == 'Y'
