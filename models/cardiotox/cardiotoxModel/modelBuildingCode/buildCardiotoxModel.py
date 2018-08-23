from indra.tools.small_model_tools.modelScope.buildPriorModel import build_prior
from indra.statements import *


##TESTING
model_types = [Phosphorylation,Dephosphorylation,ActiveForm,IncreaseAmount,DecreaseAmount,Complex] #Gef,GtpActivation]#Check sos and raf stmts 
model_genes = ['JUN','STAT1','PKM','RPS6','AURKA','HIF1A','MYC','PDGFRA','KDR','PDGF','PDGFA','FLT3LG','SRC','SORAFENIB']
#model_genes = ['JUN','STAT1','PKM','RPS6','AURKA','HIF1A','MYC','PDGFRA','KDR','PDGF','FLT3LG','PDGFA','SORAFENIB','GenericPhosphatase','SRC']
#reach_stmts = ac.load_statements('../intermediateStmts/rawStmts.pkl')
reach_stmts = '../intermediateStmts/rawStmts.pkl'


prior_stmts = build_prior(model_genes,model_types,dbs=False,additional_stmts_files=['../intermediateStmts/rawStmts.pkl'])




