from indra.statements import *
from indra.tools.small_model_tools import enforceCascadeContext as cs
from indra.tools.small_model_tools import combinePhosphorylationSites as ptm
from indra.tools import assemble_corpus as ac 

#recList = ['PDGFRA','FLT3','KDR']
#ligList = ['PDGFA','FLT3LG','VEGFA']



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


