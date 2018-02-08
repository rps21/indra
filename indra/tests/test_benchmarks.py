from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import dirname, abspath, join
from indra.benchmarks import bioprocesses as bp
# from indra.benchmarks import complexes as cp
from indra.benchmarks import phosphorylations as phos
from indra.util import unicode_strs
from nose.plugins.attrib import attr
import unittest

eval_file = join(dirname(abspath(__file__)),
                 '../benchmarks/assembly_eval/batch4/reach/' +
                 'reach_stmts_batch_4_eval.pkl')

# Removing due to QuickGO web service reliability issues
#def test_bioprocesses():
#    """Smoke test to see if bioprocesses analysis works."""
#    bp.analyze(eval_file)

#def test_bioprocesses_get_genes():
#    gene_set = bp.get_genes_for_go_id('GO:0006915')
#    assert gene_set
#    assert unicode_strs(gene_set)

@attr('nonpublic', 'webservice')
@unittest.skip('Complex analysis has been removed, test should be too, later.')
def test_complexes():
    """Smoke test to see if complexes analysis works."""
    cp.analyze(eval_file)

