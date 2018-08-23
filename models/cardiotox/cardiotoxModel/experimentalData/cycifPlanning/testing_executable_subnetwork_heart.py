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


def get_subnetwork(statements, nodes, relevance_network=None,
                   relevance_node_lim=10):

    if relevance_network is not None:
        relevant_nodes = _find_relevant_nodes(nodes, relevance_network,
                                              relevance_node_lim)
        all_nodes = nodes + relevant_nodes
    else:
        all_nodes = nodes

#    stmts = ac.filter_by_type(stmts1, Complex, invert=True)
    stmts = ac.filter_direct(stmts1)
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)
#    stmts = ac.filter_gene_list(stmts, data_genes, 'all')
    #Double check what these do, may be some overlap with functions below.
    stmts = ac.filter_enzyme_kinase(stmts)
    stmts = ac.filter_mod_nokinase(stmts)
    stmts = ac.filter_transcription_factor(stmts)
    stmts = ac.run_preassembly(stmts)  #definitely check this one. does it include others?

#    my_direct = ac.map_sequence(my_direct)


    stmts_filtered = ac.filter_gene_list(stmts, nodes, 'all')


    pa = PysbAssembler('two_step')
    pa.add_statements(stmts_filtered)
    model = pa.make_model()
    return model, stmts_filtered
#    return stmts_filtered

def _find_relevant_nodes(query_nodes, relevance_network, relevance_node_lim):
    """Return a list of nodes that are relevant for the query.

    ParametersSix months ago our washer started leaking into the apartment below ours so we had to stop using it. For two months I had to do laundry at my brothers while the other two spent $50 a week getting th
    ----------
    query_nodes : list[str]
        A list of node names to query for.
    relevance_network : str
        The UUID of the NDEx network to query relevance in.
    relevance_node_lim : int
        The number of top relevant nodes to return.

    Returns
    -------
    nodes : list[str]
        A list of node names that are relevant for the query.
    """
    all_nodes = relevance_client.get_relevant_nodes(relevance_network,
                                                    query_nodes)
    nodes = [n[0] for n in all_nodes[:relevance_node_lim]]
    return nodes

#from indra.benchmarks import phosphorylations
#def get_kinase_activities():
##    kinase_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../indra/resources/kinases.tsv') doesn't work interactively
#    kinase_file = '../../indra/resources/kinases.tsv' 
#    kinases = []
#    with open(kinase_file, 'rt') as fh:
#        lines = [l.strip() for l in fh.readlines()]
#        for lin in lines[1:]:
#            up_id, hgnc_name, _, _ = lin.split('\t')
#            hgnc_id = hgnc_client.get_hgnc_id(hgnc_name)
#            agent = Agent(hgnc_name, db_refs={'UP': up_id, 'HGNC': hgnc_id})
#            kinases.append(agent)
#    kin_activities = []
#    from indra.statements import HasActivity
#    for kin in kinases:
#        stmt = HasActivity(kin, 'kinase', True)
#        kin_activities.append(stmt)
#    return kin_activities


#def filter_enzkinase(stmts):
#    kinase_activities = get_kinase_activities()
#    stmts_enzkinase = []
#    for stmt in stmts:
#        is_kinase = False
#        for kin in kinase_activities:
#            if stmt.enz.entity_matches(kin.agent):
#                is_kinase = True
#                break
#            if kin.agent.refinement_of(stmt.enz, hierarchies):
#                is_kinase = True
#                break
#        if is_kinase:
#            stmts_enzkinase.append(stmt)
#    return stmts_enzkinase

if __name__ == '__main__':
#    genes = ['EGF', 'EGFR', 'ERBB2', 'GRB2', 'SOS1', 'HRAS', 'RAF1',
#            'MAP2K1', 'MAPK1']
#    genes = ['MAP2K1', 'MAPK1']
#    ac.filter_gene_list(stmts, data_genes, 'all')

    with open('processed_node_names.pkl','rb') as f:
        genes = pickle.load(f)
    with open('putative_targets.txt','r') as f2:
        targets = f2.readlines()
    targets_pre = targets[0].strip().split(',')
#    with open('metabolism_genes.txt','r') as f3:
#        metgenes = f3.readlines()
#    metenes = metgenes[0].

    metgenes = ['CMDN','TRPN','CAMK','MYO','CSR','JSR','BSR','BSL','CSQN','ALR','HK','PFK','ALD','GAPDH','PGK','ENOL','PK','GP','PGM','G6PDH','6PGDH','R5PI','RuPE','Tkt','TAL','XyDH','GR','GSH','CAT','SOD']
    targets = []
    for thing in targets_pre:
        thing = thing.strip('"')
        targets.append(thing)        


    gene_list_pre = targets + genes
    gene_list = []
    for gene in gene_list_pre:
        gene = gene.strip('"')
        gene = gene.strip("'")
        gene_list.append(gene)
        
    stmts1 = ac.load_statements('updated_preassembled_stmts.pkl')

    rasmachine_network = '50e3dff7-133e-11e6-a039-06603eb7f303'
    my_model, my_stmts = get_subnetwork(stmts1, gene_list)
#    my_stmts = get_subnetwork(stmts1, gene_list)

bngl_model_filter = pysb.export.export(my_model,'bngl')
bngl_file_filter = open('kdr_new_heart_model_small_all_superfiltered.bngl','w')
bngl_file_filter.write(bngl_model_filter)
bngl_file_filter.close()




