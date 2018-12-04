#from indra.statements import *
#from indra.databases import hgnc_client

#def buildDrugStmt(drug,target,stmtType):
#    #ground drug agent 
#    drugAg = Agent(drug)
#    drugAg


#    #ground target agent 
#    targetAg = Agent(target)
#    hgnc_id = hgnc_client.get_hgnc_id('RPS6')
#    if hgnc_id:
#        ag2.db_refs['HGNC'] = hgnc_id
#    standard_up_id = hgnc_client.get_uniprot_id(hgnc_id)
#    if standard_up_id:
#        ag2.db_refs['UP'] = standard_up_id

#    expStmt = stType(drugAg,targetAg)
#    return expStmt

##'SORAFENIB dephosphorylates RPS6'
#drug='Sorafenib'
#target = 'RPS6'
#stmtType = Dephosphorylation
#ag1.db_refs = {'CHEBI': 'CHEBI:50924', 'PUBCHEM': '216239', 'TEXT': 'sorafenib'}
#expStmt = buildDrugStmt(drug,target,stmtType)

#Try to reformulate these in a way to not rely on trips and take in lists of genes instead of nl sentences 


import pickle
from indra.sources import trips
from indra.tools.small_model_tools.modelContext import enforceCascadeContext as cs
from indra.tools.small_model_tools.modelContext import combinePhosphorylationSites as ptm
from indra.tools.small_model_tools.modelContext import addImplicitMechs as aim

def buildExperimentalStatements(sentences):
    tp = trips.process_text(sentences)
    expStmts = tp.statements
    return expStmts

def buildDrugTargetStmts(sentences):
    tp = trips.process_text(sentences)
    drugStmts = tp.statements
    drugStmts = cs.add_all_af(drugStmts)
    drugStmts = ptm.coarse_grain_phos(drugStmts)
    drugStmts = cs.combine_multiple_activeforms(drugStmts)
    return drugStmts







