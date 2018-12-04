import pickle
import tempfile
from indra.tools.small_model_tools.modelScope.buildPriorModel import build_prior
from indra.statements import *
from indra.tools import assemble_corpus as ac
from indra.tools.small_model_tools.modelScope.buildDrugTreatmentStatements import buildExperimentalStatements
from indra.tools.small_model_tools.modelScope.buildDrugTreatmentStatements import buildDrugTargetStmts
from indra.tools.small_model_tools.modelScope.buildUseCaseModel import expandModel
from indra.tools.small_model_tools.modelScope.addBNGLParameters import addObservables
from indra.tools.small_model_tools.modelScope.addBNGLParameters import addSimParamters
from indra.tools.small_model_tools.modelContext import extraModelReductionTools as ex
from indra.assemblers import PysbAssembler
import pysb
from indra.sources import trips
from indra.tools.small_model_tools.modelContext import enforceCascadeContext as cs
from indra.tools.small_model_tools.modelContext import combinePhosphorylationSites as ptm
from indra.tools.small_model_tools.modelContext import addImplicitMechs as aim




#testingStmts = ac.load_statements('braf_mut_context_testing.pkl')
testingStmts = ac.load_statements('af_false_test.pkl')
contextStmts = cs.add_all_af(testingStmts)
ptmStmts = ptm.coarse_grain_phos(contextStmts)
context2 = cs.combine_multiple_phos_activeforms(ptmStmts)
finalStmts = aim.addAll(context2)

drug_stmts = ac.load_statements('testing_drug_mut.pkl')
modelStmts = finalStmts + drug_stmts

pa = PysbAssembler()
pa.add_statements(modelStmts)
pysbModel = pa.make_model()
#WARNING: [2018-08-29 16:12:59] indra/pysb_assembler - HIF1A transcribes itself, skipping???

model_obs = addObservables(pysbModel,bound=True)
drug='VEMURAFENIB'    

bngl_model_filter = pysb.export.export(model_obs,'bngl')
bngl_file_filter = open('vemurafenib_bngl_model.bngl','w')
bngl_file_filter.write(bngl_model_filter+'\n')
bngl_file_filter.close()
model_actions = addSimParamters(method='ode',equil=True,equilSpecies=['%s()' % drug],viz=True)
bngl_file_filter = open('vemurafenib_bngl_model.bngl','a')
bngl_file_filter.write(model_actions)
bngl_file_filter.close()


