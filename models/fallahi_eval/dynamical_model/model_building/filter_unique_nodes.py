#my_model_stmts_larger_all

def filter_unique_species(stmts):
    nodes = []
    for st in stmts:
        for ag in st.agent_list():
            nodes.append(str(ag))

    unique_species = list(set(nodes))
    return unique_species

def filter_unique_nodes(stmts):
    nodes = []
    for st in stmts:
        for ag in st.agent_list():
            nodes.append(str(ag.name))

    unique_nodes = list(set(nodes))
    return unique_nodes



