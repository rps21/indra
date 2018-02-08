import itertools
from indra.assemblers import EnglishAssembler
from indra.explanation.reporting import stmts_from_path
from indra.explanation.model_checker import ModelChecker
from util import pklload, pkldump
from process_data import *


def get_task_5(data, inverse=False):
    """Return the test cases to be explained for Task 5."""
    # TASK 5: We observe a dose-dependent increase in total c-Jun
    # in the cell line RVH421 but a dose-dependent decrease in total
    # c-Jun in the cell line C32.
    antibody_agents = get_antibody_agents()
    obs_agents = antibody_agents['Total c-Jun']
    # We fix the time point to 24 hours
    time = 24
    # We look at doses at least 0.1 since the effect is only observed there
    dose_lower_bound = 0.1
    # Structure: cell line / drug / dose / time
    stmts_to_check = {}
    for cell_line  in ('C32', 'RVH421'):
        stmts_to_check[cell_line] = {}
        for drug in drug_names.keys():
            stmts_to_check[cell_line][drug] = {}
            drug_agent = Agent(drug, db_refs=drug_grounding[drug])
            for dose in drug_doses:
                if dose < dose_lower_bound:
                    continue
                values = get_agent_values_for_condition(data, cell_line,
                                                        drug, time, dose)
                stmts_to_check[cell_line][drug][dose] = [[], values]
                for obs in obs_agents:
                    if (cell_line == 'C32' and not inverse) or \
                        (cell_line == 'RVH421' and inverse):
                        st = DecreaseAmount(drug_agent, obs)
                    else:
                        st = IncreaseAmount(drug_agent, obs)
                    stmts_to_check[cell_line][drug][dose][0].append(st)
    return stmts_to_check


##################
# Refactor these in reusable chunks into explanations/reporting
def export_paths(scored_paths, model, stmts):
    """Export paths for pathway map in JSON-like format."""
    conc = 1.0
    time = 10
    paths = {}
    for cell_line in scored_paths.keys():
        for drug in scored_paths[cell_line].keys():
            scpaths = scored_paths[cell_line][drug]
            path, score = scpaths[0]
            label = '%s_%s_%s_%s' % (drug, time, conc, cell_line)
            paths[label] = {'meta': [], 'path': []}
            path_stmts = stmts_from_path(path, model, stmts)
            uuids = [stmt.uuid for stmt in path_stmts]
            paths[label]['path'] = uuids
    return paths


def report_paths(scored_paths, model, stmts, cell_line):
    """Report paths for a specific cell line."""
    citations = {}
    citation_count = 1
    ab_name = 'Total c-Jun'
    if cell_line == 'C32':
        pol = 'decreased'
    else:
        pol = 'increased'
    for drug in scored_paths.keys():
        paths = scored_paths[drug]
        for path, score in paths[:1]:
            title = 'How does %s treatment result in %s %s' % \
                (drug, pol, ab_name)
            title += ' in %s cells?' % cell_line
            print(title)
            print('=' * len(title))
            path_stmts = stmts_from_path(path, model, stmts)
            sentences = []
            for i, stmt in enumerate(path_stmts):
                if i == 0:
                    target = stmt.agent_list()[0].name
                    sentences.append('%s is a target of %s.' %
                                     (target, drug))
                # Make citations
                pmids = [ev.pmid for ev in stmt.evidence if ev.pmid]
                cit_nums = []
                for pmid in pmids:
                    cit_num = citations.get(pmid)
                    if cit_num is None:
                        citations[pmid] = citation_count
                        cit_num = citation_count
                        citation_count += 1
                    cit_nums.append(cit_num)
                if cit_nums:
                    cit_nums = sorted(list(set(cit_nums)))
                    cit_str = ' [%s]' % (','.join([str(c) for c
                                                  in cit_nums]))
                else:
                    cit_str = ''
                ea = EnglishAssembler([stmt])
                sentence = ea.make_model()
                sentence = sentence[:-1] + cit_str + '.'
                sentences.append(sentence)
            sentences[-1] = sentences[-1][:-1] + \
                ', which is measured by %s.' % ab_name
            print(' '.join(sentences))
            print()
    references = 'References\n==========\n'
    for k, v in sorted(citations.items(), key=lambda x: x[1]):
        references += '[%d] https://www.ncbi.nlm.nih.gov/pubmed/%s\n' % (v, k)
    print(references)
#################

def get_agent_values(antibody_agents, values):
    agent_values = {}
    for ab, value in values.items():
        for agent in antibody_agents[ab]:
            agent_values[agent] = value
    return agent_values


def get_global_mc(model, stmts_to_check, agents_to_observe):
    all_stmts_condition = []
    for cell_line in stmts_to_check.keys():
        for drug in stmts_to_check[cell_line].keys():
            stmts_condition, _ = stmts_to_check[cell_line][drug][1.0]
            all_stmts_condition += stmts_condition
    mc = ModelChecker(model, all_stmts_condition, agents_to_observe, do_sampling=True)
    mc.prune_influence_map()
    return mc


def flatten(x): return list(itertools.chain.from_iterable(x))


if __name__ == '__main__':
    INVERSE = False
    if INVERSE:
        print('Running Task 5 in INVERSE mode as control')
        print('=========================================')
    data = read_rppa_data()
    antibody_agents = get_antibody_agents()
    agents_to_observe = []
    for agents in antibody_agents.values():
        agents_to_observe += agents

    # Get all the data for Task 5
    stmts_to_check = get_task_5(data, INVERSE)
    dose = 1.0
    scored_paths = {}
    models = {}
    global_mc = None
    # Run model checking per cell line
    for cell_line in stmts_to_check.keys():
        print('Cell line: %s\n=============' % cell_line)
        model = pklload('pysb_model_%s' % cell_line)
        scored_paths[cell_line] = {}
        # Make a Model Checker with the
        # - model contextualized to the cell line
        # - the statements for the given condition
        # - agents for which observables need to be made
        global_mc = get_global_mc(model, stmts_to_check, agents_to_observe)
        for drug in stmts_to_check[cell_line].keys():
            print('Drug: %s\n=============' % drug)
            # Get all the statements and values for the condition
            stmts_condition, values_condition = \
                stmts_to_check[cell_line][drug][dose]
            path_results = []
            for stmt in stmts_condition:
                pr = global_mc.check_statement(stmt, 100, 7)
                path_results.append(pr)
            # Get a dict of measured values by INDRA Agents for this condition
            agent_values = get_agent_values(antibody_agents, values_condition)
            # Get a single list of paths for this condition
            # The length of this will be the max path number times the number
            # of Statements being checked
            paths = flatten([pr.paths for pr in path_results])
            # Run score paths with
            # - the list of paths the model checker produced
            # - the values of observed agents in the current condition
            # - loss of function mode since these are inhibitory drugs
            # - sigma corresponding to the log2 normalized measurements
            scored_paths_condition = global_mc.score_paths(paths, agent_values,
                                                    loss_of_function=True,
                                                    sigma=0.5)
            scored_paths[cell_line][drug] = scored_paths_condition
        models[cell_line] = global_mc.model

    # Dump results in standard folder
    fname = 'task5_scored_paths'
    if INVERSE:
        fname += '_inverse'
    pkldump(fname, (scored_paths, models))
