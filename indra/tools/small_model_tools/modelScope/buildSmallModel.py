import pickle
from indra.tools import assemble_corpus as ac
from indra.sources import trips
from modelContext import enforceCascadeContext as cs
from modelContext import combinePhosphorylationSites as ptm
from modelContext import addImplicitMechs as aim
from indra.assemblers import PysbAssembler
import pysb
from addBNGLParameters import addSimParamters
from addBNGLParameters import addObservables
from modelContext import extraModelReductionTools as ex

def buildSmallModel(stmts):

    finalStmts = cs.add_all_af(stmts)
    finalStmts2 = ptm.coarse_grain_phos(finalStmts)
    finalStmts3 = aim.addAll(finalStmts2)
    finalStmts4 = ex.removeDimers(finalStmts3)
    finalStmts5 = ex.removeMutations(finalStmts4)

    return finalStmts5

#def writeBNGLFile(stmts):
#    pa = PysbAssembler()
#    pa.add_statements(finalStmts5)
#    modifiedModel = pa.make_model()
#    modifiedModel = addObservables(modifiedModel,bound=True)

#    if writeBNGL:
#        bngl_model = pysb.export.export(modifiedModel,'bngl')
#        with open('bnglModelTesting/modifiedModel.bngl','w') as bnglFile:
#            bnglFile.write(bngl_model)
#            actionsBlock = addSimParamters('ode',True,['EGF(erbb)'])
#            bnglFile.write(actionsBlock)
