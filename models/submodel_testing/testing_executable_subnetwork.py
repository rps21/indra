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


#model, my_direct, my_indirect, my_stmts = get_subnetwork(stmts1, genes)
def get_subnetwork(statements, nodes, relevance_network=None,
                   relevance_node_lim=10):

    if relevance_network is not None:
        relevant_nodes = _find_relevant_nodes(nodes, relevance_network,
                                              relevance_node_lim)
        all_nodes = nodes + relevant_nodes
    else:
        all_nodes = nodes

    filtered_statements = _filter_statements(statements, all_nodes)
    my_direct, my_indirect = filter_direct(filtered_statements)

	#adding preassembler to resolve hierarchies
    my_direct = ac.run_preassembly(my_direct)
    my_direct = ac.map_sequence(my_direct)
    pra = Preassembler(hierarchies, my_direct)
    stmts_unique = pra.combine_related()

    pa = PysbAssembler('two_step')
    pa.add_statements(stmts_unique)
    model = pa.make_model()
    return model, my_direct, my_indirect, stmts_unique

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
    #THIS SHOULD BE DONE BETTER
    indirect_phos_stmts = []
    for st in indirect_stmts:
        if 'Phos' in str(st):
            indirect_phos_stmts.append(st)
    direct_phos = filter_enzkinase(indirect_phos_stmts)
    direct_stmts.extend(direct_phos)
    return direct_stmts, indirect_stmts

#from indra.benchmarks import phosphorylations
def get_kinase_activities():
#    kinase_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../indra/resources/kinases.tsv') doesn't work interactively
    kinase_file = '../../indra/resources/kinases.tsv' 
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
#    genes = ['EGF', 'EGFR', 'ERBB2', 'GRB2', 'SOS1', 'HRAS', 'RAF1',
#            'MAP2K1', 'MAPK1']
    genes = ['MAP2K1', 'MAPK1']

    with open('reading/model-2016-11-30-10-18-57.pkl', 'rb') as f:
#    with open('phase3_preassembled.pkl', 'rb') as f:
        model = pickle.load(f)
    stmts1 = []
    for k, v in model.items():
        stmts1 += v

    #preassembly
    #stmts1 = ac.run_preassembly(stmts1)
    #stmts1 = ac.map_sequence(stmts1)

    rasmachine_network = '50e3dff7-133e-11e6-a039-06603eb7f303'
    #model_nofilter = get_subnetwork(stmts, genes)
    model, my_direct, my_indirect, my_stmts = get_subnetwork(stmts1, genes)


bngl_model_filter = pysb.export.export(model,'bngl')
bngl_file_filter = open('rbm/mapk_sitemap.bngl','w')
bngl_file_filter.write(bngl_model_filter)
bngl_file_filter.close()


##Loop through list of indra statements marked as direct
##Create new lists based on mechanism type (phosphorylation, binding, activation)
#direct_phos = []
#direct_complex = []
#direct_activation = []
#for st in my_direct:
#	if 'Phos' in str(st):
#		direct_phos.append(st)
#	elif 'Complex' in str(st):
#		direct_complex.append(st)
#	elif 'Act' in str(st):
#		direct_activation.append(st)


##Recreate total list of direct mechanisms, now grouped by mechanism type.
##There's probably a more efficient way to reorder statments
#direct_total = []
#for st in direct_complex:
#	direct_total.append(st)
#for st in direct_phos:
#	direct_total.append(st)
#for st in direct_activation:
#	direct_total.append(st)


##Loop through list of indra statements marked as indirect
##Create new lists based on mechanism type (phosphorylation, binding, activation)
#indirect_phos = []
#indirect_complex = []
#indirect_activation = []
#for st in my_indirect:
#	if 'Phos' in str(st):
#		indirect_phos.append(st)
#	elif 'Complex' in str(st):
#		indirect_complex.append(st)
#	elif 'Act' in str(st):
#		indirect_activation.append(st)


##Recreate total list of indirect mechanisms, now grouped by mechanism type.
##There's probably a more efficient way to reorder statments
#indirect_total = []
#for st in indirect_complex:
#	indirect_total.append(st)
#for st in indirect_phos:
#	indirect_total.append(st)
#for st in indirect_activation:
#	indirect_total.append(st)


##Create separate litsts of statement, evidence, and epistemics for each direct mechanism
##Group into a list of lists with common ordering
##Think this format may be beneficial for pandas
#dir_statement = []
#dir_evidence = []
#dir_epistemics = []
#for st in direct_total:
#	dir_statement.append(st)
#	dir_evidence.append(st.evidence)
#	dir_epistemics.append(st.evidence[0].epistemics)

##Clean up evidence strings
#dir_rule = []
#dir_sentence = []
#for ev in dir_evidence:
#	#for evidence info, first split on opening {
#	list1 = re.split(r'{',str(ev))
#	
#	#Discard first of two elements in list, contains method(reach) and PMID
#	#Split second element on closing }
#	list2 = re.split(r'}',list1[1])
#	
#	#First element of new list contains reach rule, but needs additional parsing. Second element contains supporting 	sentence
#	#Split first element on found by tag that indicates reach rule
#	list3=re.split(r"found_by': u'",list2[0])
#	
#	#Split again on species tag to remove trailing information
#	rule_list = re.split(r"', 'sp",list3[1])
#	#Rule name is first element of resulting list
#	dir_rule.append(rule_list[0])
#	
#	#Strip special characters from sentence in list2[1]
#	dir_sentence.append(re.sub('[^A-Za-z0-9 ]+', '', list2[1]))


##Compile lists into dict, make dataframe with pandas
#d = {'Sentence':dir_sentence ,'Rule':dir_rule, 'Epistemics':dir_epistemics, 'Statement':dir_statement}
#df = pd.DataFrame(data=d)

##reorder columns in df
#cols = df.columns.tolist()
#cols = [cols[3],cols[0],cols[1],cols[2]]
#df = df[cols]
##df.to_csv('all_direct_statements.csv')


##Create separate litsts of statement, evidence, and epistemics for each indirect mechanism
##Group into a list of lists with common ordering
#indir_statement = []
#indir_evidence = []
#indir_epistemics = []
#for st in indirect_total:
#	indir_statement.append(st)
#	indir_evidence.append(st.evidence)
#	indir_epistemics.append(st.evidence[0].epistemics)

##Clean up evidence strings
#indir_rule = []
#indir_sentence = []
#for ev in indir_evidence:
#	#for evidence info, first split on opening {
#	list1 = re.split(r'{',str(ev))
#	
#	#Discard first of two elements in list, contains method(reach) and PMID
#	#Split second element on closing }
#	list2 = re.split(r'}',list1[1])
#	
#	#First element of new list contains reach rule, but needs additional parsing. Second element contains supporting 	sentence
#	#Split first element on found by tag that indicates reach rule
#	list3=re.split(r"found_by': u'",list2[0])
#	
#	#Split again on species tag to remove trailing information
#	rule_list = re.split(r"', 'sp",list3[1])
#	#Rule name is first element of resulting list
#	indir_rule.append(rule_list[0])
#	
#	#Strip special characters from sentence in list2[1]
#	indir_sentence.append(re.sub('[^A-Za-z0-9 ]+', '', list2[1]))


##Compile lists into dict, make dataframe with pandas
#d = {'Sentence':indir_sentence ,'Rule':indir_rule, 'Epistemics':indir_epistemics, 'Statement':indir_statement}
#df = pd.DataFrame(data=d)
##reorder columns in df
#cols = df.columns.tolist()
#cols = [cols[3],cols[0],cols[1],cols[2]]
#df = df[cols]
##df.to_csv('all_indirect_statements.csv')






