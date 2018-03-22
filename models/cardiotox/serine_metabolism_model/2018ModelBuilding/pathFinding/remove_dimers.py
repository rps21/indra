def remove_dimers(stmts_in):
    stmts_new = stmts_in[:]
    for st in stmts_in:
        if len(st.agent_list()) > 1:
            if st.agent_list()[0]: #Weird circumstance with a None type in transcription
                if st.agent_list()[0].entity_matches(st.agent_list()[1]):
                    stmts_new.remove(st)
    return stmts_new

