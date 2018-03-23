from copy import deepcopy
from indra.statements import *
from indra.tools import assemble_corpus as ac

def remove_dimers(stmts_in):
    stmts_to_remove = []
    for st in stmts_in:
        if len(st.agent_list()) > 1:
            if st.agent_list()[0]: #Can have a None type in agent list, better workaround probably exists 
                if st.agent_list()[0].entity_matches(st.agent_list()[1]):
                    stmts_to_remove.append(st)
    stmts_new = [st for st in stmts_in if st not in stmts_to_remove]
    return stmts_new

def remove_mutations(stmts_in):
    stmts_to_remove = []
    for st in stmts_in:
        for ag in st.agent_list():
                if ag.mutations:
                    stmts_to_remove.append(st)
                    break
    stmts_new = [st for st in stmts_in if st not in stmts_to_remove]
    return stmts_new

def remove_bad_translocations(stmts_in):
    tran_stmts = ac.filter_by_type(stmts_in,Translocation)
    stmts_to_remove = []
    for st in tran_stmts:
        if st.to_location is None or st.from_location is None:
            stmts_to_remove.append(st)
    stmts_new = [st for st in stmts_in if st not in stmts_to_remove]
    return stmts_new


