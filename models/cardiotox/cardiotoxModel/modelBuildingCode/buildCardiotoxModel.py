from indra.tools.small_model_tools.modelScope.buildPriorModel import build_prior
from indra.statements import *
from indra.tools import assemble_corpus as ac
import pickle
from indra.tools.small_model_tools.modelScope.buildDrugTreatmentStatements import buildExperimentalStatements
from indra.tools.small_model_tools.modelScope.buildDrugTreatmentStatements import buildDrugTargetStmts



#model_types = [Phosphorylation,Dephosphorylation,ActiveForm,IncreaseAmount,DecreaseAmount,Complex] #Gef,GtpActivation]#Check sos and raf stmts 
#model_genes = ['JUN','STAT1','PKM','RPS6','AURKA','HIF1A','MYC','PDGFRA','KDR','PDGF','PDGFA','FLT3LG','SRC','SORAFENIB']
##model_genes = ['JUN','STAT1','PKM','RPS6','AURKA','HIF1A','MYC','PDGFRA','KDR','PDGF','FLT3LG','PDGFA','SORAFENIB','GenericPhosphatase','SRC']
##reach_stmts = ac.load_statements('../intermediateStmts/rawStmts.pkl')
#reach_stmts = './intermediateStmts/rawStmts.pkl'


#prior_stmts = build_prior(model_genes,model_types,dbs=True,additional_stmts_files=['./intermediateStmts/rawStmts.pkl'])
#drugSentences = 'SORAFENIB dephosphorylates PDGFRA. SORAFENIB dephosphorylates KDR. SORAFENIB dephosphorylates FLT3'
#drug_stmts = buildDrugTargetStmts(drugSentences)

#prior_model_stmts = prior_stmts + drug_stmts


#ac.dump_statements(prior_model_stmts,'cardiotoxPrior.pkl')

###########################################################################

prior_model_stmts = ac.load_statements('cardiotoxPrior.pkl')
expSentences = 'SORAFENIB dephosphorylates RPS6. SORAFENIB phosphorylates PKM. SORAFENIB transcribes HIF1A.'
exp_stmts = buildExperimentalStatements(expSentences)

initialStmts = prior_model_stmts
initialNodes = []
for st in prior_model_stmts:
    for ag in st.agent_list():
        if ag.name not in initialNodes:
            initialNodes.append(ag.name)


expObservations = {'RPS6':['phosphorylation',exp_stmts[0]]}
drug = 'SORAFENIB'
drugTargets = ['FLT3','PDGFRA','KDR']
stmts = prior_model_stmts


modelStmts = expandModel(expObservations,drug,drugTargets,initialStmts,initialNodes)

#Outstanding issues:
#Ligands in statements.
#   tmp fix: include as input
#   better fix: identify missing ligands during context generation, partially implemented
#phosphorylation states in final model - phos_act and kinase, False?
#Filter ubiquintation
