from indra.statements import *
from indra.tools.small_model_tools import enforceCascadeContext as cs
from indra.tools.small_model_tools import combinePhosphorylationSites as ptm
from indra.tools import assemble_corpus as ac 

#recList = ['PDGFRA','FLT3','KDR']
#ligList = ['PDGFA','FLT3LG','VEGFA']

#def getLigAgent(lig_name,stmts):
#    st in stmts:
#        for ag in st.agent_list():
#            if ag.name == lig_name:
#                lig_agent = ag
#    return lig_agent

def getLigAgent(rec,stmts):
    lig_names, rec_names, lig_list, rec_list = cs.find_rec_lig(stmts)
    lig_agent = None
    for st in stmts:
        for ag in st.agent_list():
            if ag.name == rec:
                for ag in st.agent_list():
                    for lig in lig_names:
                        if ag.name == lig:
                            lig_agent = ag  
    return lig_agent




def fixRecPhosContext(rec_names, stmts):
    phosStmts = ac.filter_by_type(stmts,Phosphorylation)
    newRecPhosStmts = []
    for rec in rec_names:
        recPhosStmts = ac.filter_gene_list(phosStmts,[rec],'one')
        recPhosStmts_self = []
        for st in recPhosStmts:
            if st.enz.name == rec and st.sub.name == rec:
                recPhosStmts_self.append(st)

        for st in recPhosStmts_self:
            newSt = deepcopy(st)
            newSt.enz.mods = []
            lig_ag = getLigAgent(rec,stmts)
            if lig_ag:
                newSt.enz.bound_conditions = [BoundCondition(lig_ag)]
                newRecPhosStmts.append(newSt)
    return newRecPhosStmts


#Testing
#recList = ['PDGFRA','FLT3','KDR']
#newRecPhosStmts = fixRecPhosContext(recList,finalStmts)

#lig_names, rec_names, lig_list, rec_list = cs.find_rec_lig(finalStmts)



#def find_rec_lig(stmts):
#    #Requires dict of receptor-ligand pairs, should store this somewhere better
##    dir_path = os.path.dirname(os.path.realpath(__file__))
#    dir_path = '/home/bobby/Dropbox/Sorger_Lab/indra/indra/tools/small_model_tools/'
#    with open(dir_path+'/lig_rec_dict.pkl','rb') as f:
#        receptor_dict = pickle.load(f)

#    new_af_stmts = []
#    lig_list, lig_names = find_agent(stmts,list(receptor_dict.keys()))
#    lig_stmts = ac.filter_gene_list(stmts,lig_names,'one')

#    rec_list = []
#    rec_names = []
#    for sublist in list(receptor_dict.values()):
#        rec_list_init, rec_names_init = find_agent(lig_stmts,sublist)
#        rec_list = rec_list + rec_list_init #change this to avoid repeats
#        rec_names = rec_names + rec_names_init
#    lig_list = list(set(lig_list))
#    rec_list = list(set(rec_list))
#    return lig_names, rec_names, lig_list, rec_list
