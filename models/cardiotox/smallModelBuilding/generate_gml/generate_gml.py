from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
import os
import sys
import subprocess
import argparse
import shlex
import pysb
from pysb.export import export
import indra
from indra.databases import relevance_client, hgnc_client
from indra.assemblers import PysbAssembler
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.statements import Phosphorylation, Agent, Evidence

#python generate_gml.py --list list_of_gene_names --stmts indra_created_pkl_file"

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--list', help='delimited list input', type=str)
parser.add_argument('-s', '--stmts', help='pickle file of indra statements')
args = parser.parse_args()

def get_subnetwork(statements, nodes, relevance_network=None,
                   relevance_node_lim=10):

    if relevance_network is not None:
        relevant_nodes = _find_relevant_nodes(nodes, relevance_network,
                                              relevance_node_lim)
        all_nodes = nodes + relevant_nodes
    else:
        all_nodes = nodes

    filtered_statements = _filter_statements(statements, all_nodes)
    direct_stmts = filter_direct(filtered_statements)

	#adding preassembler to resolve hierarchies
    pra = Preassembler(hierarchies, direct_stmts)
    stmts_unique = pra.combine_related()

    pa = PysbAssembler('two_step')
    pa.add_statements(stmts_unique)
    model = pa.make_model()
    return model


def _filter_statements(statements, agents):
    """Return INDRA Statements which have Agents in the given list.

    Only statements are returned in which all appearing Agents as in the
    agents list.

    Parameters
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to filter.
    agents : list[str]
        A list of agent names that need to appear in filtered statements.

    Returns
    -------
    filtered_statements : list[indra.statements.Statement]
        The list of filtered INDRA Statements.
    """
    filtered_statements = []
    for s in statements:
        if all([a is not None for a in s.agent_list()]) and \
            all([a.name in agents for a in s.agent_list()]):
            filtered_statements.append(s)
    return filtered_statements

def _find_relevant_nodes(query_nodes, relevance_network, relevance_node_lim):
    """Return a list of nodes that are relevant for the query.

    Parameters
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

def get_is_direct(stmt):
    '''Returns true if there is evidence that the statement is a direct
    interaction. If any of the evidences associated with the statement
    indicates a direct interatcion then we assume the interaction
    is direct. If there is no evidence for the interaction being indirect
    then we default to direct.'''
    any_indirect = False
    if isinstance(stmt,indra.statements.Complex):
        return True
    for ev in stmt.evidence:
        if ev.epistemics.get('direct') is True:
            return True
        elif ev.epistemics.get('direct') is False:
            # This guarantees that we have seen at least
            # some evidence that the statement is indirect
            any_indirect = True
    if any_indirect:
        return False
    if isinstance(stmt,indra.statements.Complex):
        return True
    return True

def filter_direct(stmts):
    direct_stmts = []
    indirect_stmts = []
    for stmt in stmts:
        if get_is_direct(stmt):
            direct_stmts.append(stmt)
        else:
            indirect_stmts.append(stmt)
    indirect_phos_stmts = []
    for st in indirect_stmts:
        if 'Phos' in str(st):
            indirect_phos_stmts.append(st)
    direct_phos = filter_enzkinase(indirect_phos_stmts)
    direct_stmts.extend(direct_phos)
    return direct_stmts

def get_kinase_activities():
#    kinase_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../indra/resources/kinases.tsv') doesn't work interactively
    kinase_file = 'kinases.tsv' 
    kinases = []
    with open(kinase_file, 'rt') as fh:
        lines = [l.strip() for l in fh.readlines()]
        for lin in lines[1:]:
            up_id, hgnc_name, _, _ = lin.split('\t')
            hgnc_id = hgnc_client.get_hgnc_id(hgnc_name)
            agent = Agent(hgnc_name, db_refs={'UP': up_id, 'HGNC': hgnc_id})
            kinases.append(agent)
    kin_activities = []
    from indra.statements import HasActivity
    for kin in kinases:
        stmt = HasActivity(kin, 'kinase', True)
        kin_activities.append(stmt)
    return kin_activities


def filter_enzkinase(stmts):
    kinase_activities = get_kinase_activities()
    stmts_enzkinase = []
    for stmt in stmts:
        is_kinase = False
        for kin in kinase_activities:
            if stmt.enz.entity_matches(kin.agent):
                is_kinase = True
                break
            if kin.agent.refinement_of(stmt.enz, hierarchies):
                is_kinase = True
                break
        if is_kinase:
            stmts_enzkinase.append(stmt)
    return stmts_enzkinase

if __name__ == '__main__':
    #genes = ['EGF', 'EGFR', 'ERBB2', 'GRB2', 'SOS1', 'HRAS', 'RAF1', 'MAP2K1', 'MAPK1']
    genes=args.list

    #with open('model-2016-11-30-10-18-57.pkl', 'rb') as f:
    with open(args.stmts, 'rb') as f:
        model = pickle.load(f)
    stmts = []
    for k, v in model.items():
        stmts += v

    rasmachine_network = '50e3dff7-133e-11e6-a039-06603eb7f303'
    model = get_subnetwork(stmts, genes)


bngl_model = pysb.export.export(model,'bngl')
bngl_file = open('model.bngl','w')
bngl_file.write(bngl_model)
bngl_file.close()

args_str = "perl Perl2/Visualization/visualize.pl --type contactmap --bngl model.bngl"
args = shlex.split(args_str)
subprocess.Popen(args)



