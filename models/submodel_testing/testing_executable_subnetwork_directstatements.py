from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
from indra.databases import relevance_client
from indra.assemblers import PysbAssembler
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
import pysb
from pysb.export import export
import pandas as pd
import re

def get_subnetwork(statements, nodes, relevance_network=None,
                   relevance_node_lim=10):

    if relevance_network is not None:
        relevant_nodes = _find_relevant_nodes(nodes, relevance_network,
                                              relevance_node_lim)
        all_nodes = nodes + relevant_nodes
    else:
        all_nodes = nodes

    filtered_statements = _filter_statements(statements, all_nodes)
    #filtered_phos_statements = filter_phos(filtered_statements)
    filtered_egfr_statements = filter_egf_egfr(filtered_statements)
    my_direct, my_indirect = filter_direct(filtered_egfr_statements)


    pa = PysbAssembler()
    pa.add_statements(filtered_egfr_statements)
    model = pa.make_model()
    return model, my_direct, my_indirect

def filter_egf_egfr(statements):
    
    filtered_first = []
    for st in statements:
        if 'EGF(' in str(st):
            filtered_first.append(st)
        elif 'egf(' in str(st):
            filtered_first.append(st)

    egf_egfr_stmts = []
    for st in filtered_first:
        if 'EGFR(' in str(st):
            egf_egfr_stmts.append(st)
        elif 'egfr(' in str(st):
            egf_egfr_stmts.append(st)

    return(egf_egfr_stmts)

def filter_phos(statements):
    filtered_phos = []
    for st in statements:
        if 'Phos' in str(st):
            filtered_phos.append(st)
    return(filtered_phos)


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
    for ev in stmt.evidence:
        if ev.epistemics.get('direct') is True:
            return True
        elif ev.epistemics.get('direct') is False:
            # This guarantees that we have seen at least
            # some evidence that the statement is indirect
            any_indirect = True
    if any_indirect:
        return False
    return True

def filter_direct(stmts):
    direct_stmts = []
    indirect_stmts = []
    for stmt in stmts:
        if get_is_direct(stmt):
            direct_stmts.append(stmt)
        else:
            indirect_stmts.append(stmt)
    return direct_stmts, indirect_stmts

if __name__ == '__main__':
    genes = ['EGF', 'EGFR', 'ERBB2', 'GRB2', 'SOS1', 'HRAS', 'RAF1',
            'MAP2K1', 'MAPK1']

    with open('reading/model.pkl', 'rb') as f:
        model = pickle.load(f)
    stmts1 = []
    for k, v in model.items():
        stmts1 += v
   # stmts2 = filter_direct(stmts1)
	

    rasmachine_network = '50e3dff7-133e-11e6-a039-06603eb7f303'
    #model_nofilter = get_subnetwork(stmts, genes)
    #model, my_filtered_egfr_statements, my_filtered_direct_egfr_statments = get_subnetwork(stmts1, genes)
    model, my_direct, my_indirect = get_subnetwork(stmts1, genes)


#bngl_model = pysb.export.export(model_nofilter,'bngl')
#bngl_file = open('rbm/ras_nofilter.bngl','w')
#bngl_file.write(bngl_model)
#bngl_file.close()

#bngl_model_filter = pysb.export.export(model,'bngl')
#bngl_file_filter = open('rbm/ras_newsmall.bngl','w')
#bngl_file_filter.write(bngl_model_filter)
#bngl_file_filter.close()


#Loop through list of indra statements marked as direct
#Create new lists based on mechanism type (phosphorylation, binding, activation)
direct_phos = []
direct_complex = []
direct_activation = []
for st in my_direct:
	if 'Phos' in str(st):
		direct_phos.append(st)
	elif 'Complex' in str(st):
		direct_complex.append(st)
	elif 'Act' in str(st):
		direct_activation.append(st)


#Recreate total list of direct mechanisms, now grouped by mechanism type.
#There's probably a more efficient way to reorder statments
direct_total = []
for st in direct_complex:
	direct_total.append(st)
for st in direct_phos:
	direct_total.append(st)
for st in direct_activation:
	direct_total.append(st)


#Loop through list of indra statements marked as indirect
#Create new lists based on mechanism type (phosphorylation, binding, activation)
indirect_phos = []
indirect_complex = []
indirect_activation = []
for st in my_indirect:
	if 'Phos' in str(st):
		indirect_phos.append(st)
	elif 'Complex' in str(st):
		indirect_complex.append(st)
	elif 'Act' in str(st):
		indirect_activation.append(st)


#Recreate total list of indirect mechanisms, now grouped by mechanism type.
#There's probably a more efficient way to reorder statments
indirect_total = []
for st in indirect_complex:
	indirect_total.append(st)
for st in indirect_phos:
	indirect_total.append(st)
for st in indirect_activation:
	indirect_total.append(st)


#Create separate litsts of statement, evidence, and epistemics for each direct mechanism
#Group into a list of lists with common ordering
#Think this format may be beneficial for pandas
dir_statement = []
dir_evidence = []
dir_epistemics = []
for st in direct_total:
	dir_statement.append(st)
	dir_evidence.append(st.evidence)
	dir_epistemics.append(st.evidence[0].epistemics)

#Clean up evidence strings
dir_rule = []
dir_sentence = []
for ev in dir_evidence:
	#for evidence info, first split on opening {
	list1 = re.split(r'{',str(ev))
	
	#Discard first of two elements in list, contains method(reach) and PMID
	#Split second element on closing }
	list2 = re.split(r'}',list1[1])
	
	#First element of new list contains reach rule, but needs additional parsing. Second element contains supporting 	sentence
	#Split first element on found by tag that indicates reach rule
	list3=re.split(r"found_by': u'",list2[0])
	
	#Split again on species tag to remove trailing information
	rule_list = re.split(r"', 'sp",list3[1])
	#Rule name is first element of resulting list
	dir_rule.append(rule_list[0])
	
	#Strip special characters from sentence in list2[1]
	dir_sentence.append(re.sub('[^A-Za-z0-9 ]+', '', list2[1]))


#Compile lists into dict, make dataframe with pandas
d = {'Sentence':dir_sentence ,'Rule':dir_rule, 'Epistemics':dir_epistemics, 'Statement':dir_statement}
df = pd.DataFrame(data=d)

#reorder columns in df
cols = df.columns.tolist()
cols = [cols[3],cols[0],cols[1],cols[2]]
df = df[cols]
#df.to_csv('direct_statements_egf_egfr.csv')


#Create separate litsts of statement, evidence, and epistemics for each indirect mechanism
#Group into a list of lists with common ordering
indir_statement = []
indir_evidence = []
indir_epistemics = []
for st in indirect_total:
	indir_statement.append(st)
	indir_evidence.append(st.evidence)
	indir_epistemics.append(st.evidence[0].epistemics)

#Clean up evidence strings
indir_rule = []
indir_sentence = []
for ev in indir_evidence:
	#for evidence info, first split on opening {
	list1 = re.split(r'{',str(ev))
	
	#Discard first of two elements in list, contains method(reach) and PMID
	#Split second element on closing }
	list2 = re.split(r'}',list1[1])
	
	#First element of new list contains reach rule, but needs additional parsing. Second element contains supporting 	sentence
	#Split first element on found by tag that indicates reach rule
	list3=re.split(r"found_by': u'",list2[0])
	
	#Split again on species tag to remove trailing information
	rule_list = re.split(r"', 'sp",list3[1])
	#Rule name is first element of resulting list
	indir_rule.append(rule_list[0])
	
	#Strip special characters from sentence in list2[1]
	indir_sentence.append(re.sub('[^A-Za-z0-9 ]+', '', list2[1]))


#Compile lists into dict, make dataframe with pandas
d = {'Sentence':indir_sentence ,'Rule':indir_rule, 'Epistemics':indir_epistemics, 'Statement':indir_statement}
df = pd.DataFrame(data=d)
#reorder columns in df
cols = df.columns.tolist()
cols = [cols[3],cols[0],cols[1],cols[2]]
df = df[cols]
#df.to_csv('indirect_statements_egf_egfr.csv')






