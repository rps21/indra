import requests
from indra.statements import stmts_from_json
from indra.tools import assemble_corpus as ac 
from indra.mechlinker import MechLinker
from indra.statements import *


def queryDB(gene,mod):
    timeout = 1
    timeout_counter = 0
    while timeout == 1:

        resp = requests.get('https://lsm6zea7gg.execute-api.us-east-1.amazonaws.com/production/statements/',
                            headers={'x-api-key': 'XH36SVAGBN9L3SA8XNuHu5hvJ9v3j9mq8PTkPYjG'},
                            params={'object': gene,
                                    'type': mod})
        if 'message' in resp.json():
            print(resp.json().get('message'))
            timeout_counter = timeout_counter + 1
            print('Number of timeouts: %s' % timeout_counter)
        else:
            timeout=0

    stmts_json = resp.json()
    stmts = stmts_from_json(stmts_json['statements'])
    stmts = cleanStatements(stmts)
    return stmts

def queryDB_nonmod(gene,mod):
    resp = requests.get('https://lsm6zea7gg.execute-api.us-east-1.amazonaws.com/production/statements/',
                        headers={'x-api-key': 'XH36SVAGBN9L3SA8XNuHu5hvJ9v3j9mq8PTkPYjG'},
                        params={'agent': gene,
                                'type': mod})
    stmts_json = resp.json()
    stmts = stmts_from_json(stmts_json['statements'])
    stmts = cleanStatements(stmts)
    return stmts

def getCandidateNodes(stmts):
    candidateNodes = []
    for st in stmts:
        if isinstance(st,Modification):
            candidateNodes.append(st.enz.name)
        if isinstance(st,RegulateAmount):
            candidateNodes.append(st.subj.name)
    return candidateNodes

#isinstance(st, stmt_type)
#- :py:class:`Complex`
#- :py:class:`Modification`
#- :py:class:`SelfModification`
#- :py:class:`RegulateActivity`
#- :py:class:`RegulateAmount`
#- :py:class:`ActiveForm`

def rmNoneStmts(stmts):
    outputStmts = []
    for st in stmts:
        if any([ag == None for ag in st.agent_list()]):
            pass 
        else:
            outputStmts.append(st)
    return outputStmts

def cleanStatements(stmts):
    stmts = ac.map_grounding(stmts)
#    stmts = ac.filter_grounded_only(stmts, save='intermediateStmts/grounded.pkl')
    stmts = ac.filter_grounded_only(stmts)
    stmts = ac.filter_genes_only(stmts)
    stmts = ac.filter_human_only(stmts)
#    stmts = ac.run_preassembly(stmts, save='intermediateStmts/preassembled.pkl')
    stmts = ac.run_preassembly(stmts)
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

    stmts = ml.statements
#    stmts = rmNoneStmts(stmts)
    return stmts



##Testing
##gene='PKM'
##mod = ['phosphorylation']#,'dephosphorylation','increaseamount','decreaseamount'] #,'complex']  
##for mo in mod:
##    stmts = idb.queryDB(gene,mo)
##    nodes = idb.getCandidateNodes(stmts)
##    print(nodes)
