import pickle
from indra.tools import assemble_corpus as ac
from indra.sources import trips
from indra.tools.small_model_tools.modelContext import enforceCascadeContext as cs
from indra.tools.small_model_tools.modelContext import combinePhosphorylationSites as ptm
from indra.tools.small_model_tools.modelContext import addImplicitMechs as aim
from indra.assemblers import PysbAssembler
import pysb
#from indra.tools.small_model_toolsaddBNGLParameters import addSimParamters
#from indra.tools.small_model_toolsaddBNGLParameters import addObservables
from indra.tools.small_model_tools.modelContext import extraModelReductionTools as ex

def buildSmallModel(stmts):

    finalStmts = cs.add_all_af(stmts)
    finalStmts2 = ptm.coarse_grain_phos(finalStmts)
#    print(finalStmts2)
    finalStmts3 = cs.combine_multiple_phos_activeforms(finalStmts2)
#    print(finalStmts3)
    finalStmts4 = aim.addAll(finalStmts3)
    finalStmts5 = ex.removeDimers(finalStmts4)
    finalStmts6 = ex.removeMutations(finalStmts5)

    return finalStmts6

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
