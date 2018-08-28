import pickle
import tempfile
from indra.tools.small_model_tools.modelScope.buildPriorModel import build_prior
from indra.statements import *
from indra.tools import assemble_corpus as ac
from indra.tools.small_model_tools.modelScope.buildDrugTreatmentStatements import buildExperimentalStatements
from indra.tools.small_model_tools.modelScope.buildDrugTreatmentStatements import buildDrugTargetStmts
from indra.tools.small_model_tools.modelScope.buildUseCaseModel import expandModel


model_types = [Phosphorylation,Dephosphorylation,ActiveForm,IncreaseAmount,DecreaseAmount,Complex] #Gef,GtpActivation]#Check sos and raf stmts 
drugTargets = ['PDGFRA','KDR','FLT3']
ligands = ['PDGFA','PDGF','VEGF','VEGFA','VEGFC']
modifiedNodes = ['RPS6','PKM','HIF1A']# ['JUN','STAT1','PKM','RPS6','AURKA','HIF1A','MYC']
otherNodes = []
model_genes = drugTargets + ligands + modifiedNodes + otherNodes


#reach_stmts = './intermediateStmts/rawStmts.pkl'
reach_stmts = './heart_reading/raw_stmts/reach_output.pkl'
trips_stmts = './heart_reading/raw_stmts/trips_output.pkl'
prior_stmts = build_prior(model_genes,model_types,dbs=True,additional_stmts_files=[reach_stmts,trips_stmts])
drugSentences = 'SORAFENIB dephosphorylates PDGFRA. SORAFENIB dephosphorylates KDR. SORAFENIB dephosphorylates FLT3'
drug_stmts = buildDrugTargetStmts(drugSentences)

prior_model_stmts = prior_stmts + drug_stmts

ac.dump_statements(prior_model_stmts,'cardiotoxPrior.pkl')
#ac.dump_statements(prior_model_stmts,'cardiotoxPrior_testing.pkl')

###########################################################################

#prior_model_stmts = ac.load_statements('cardiotoxPrior.pkl')
prior_model_stmts = ac.load_statements('cardiotoxPrior_testing.pkl')
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
#drugTargets = ['FLT3','PDGFRA','KDR']
drugTargets = ['PDGFRA']
stmts = prior_model_stmts


modelStmts = expandModel(expObservations,drug,drugTargets,ligands,initialStmts,initialNodes)

#Outstanding issues:
#Ligands in statements.
#   tmp fix: include as input
#   better fix: identify missing ligands during context generation, partially implemente





