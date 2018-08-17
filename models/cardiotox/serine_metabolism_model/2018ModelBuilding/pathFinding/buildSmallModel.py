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


#sentences = 'EGF binds EGFR. EGFR phosphorylates EGFR. EGFR binds EGFR.EGFR binds GRB2. GRB2 bounds SOS1, EGFR binds SOS1. SOS1 phosphorylates KRAS. KRAS phosphorylates BRAF. KDR binds VEGFC. KDR binds SRC. KDR phosphorylates SRC. SRC phosphorylates AKT1.'
##sentences = 'EGF binds EGFR. EGFR phosphorylates EGFR. EGFR binds EGFR. EGFR binds GRB2. KDR binds VEGFC. KDR binds SRC'

#tp = trips.process_text(sentences)
#testStmts = tp.statements


testStmts = ac.load_statements('contextTestStmts.pkl')
finalStmts = cs.add_all_af(testStmts)
finalStmts2 = ptm.coarse_grain_phos(finalStmts)
finalStmts3 = aim.addAll(finalStmts2)

pa = PysbAssembler()
pa.add_statements(testStmts)
originalModel = pa.make_model()
orignalModel = addObservables(originalModel,bound=True)

bngl_model = pysb.export.export(originalModel,'bngl')
with open('bnglModelTesting/originalModel.bngl','w') as bnglFile:
    bnglFile.write(bngl_model)
    actionsBlock = addSimParamters('ode',True,['EGF(erbb)'])
    bnglFile.write(actionsBlock)


pa = PysbAssembler()
pa.add_statements(finalStmts3)
modifiedModel = pa.make_model()
modifiedModel = addObservables(modifiedModel,bound=True)

bngl_model = pysb.export.export(modifiedModel,'bngl')
with open('bnglModelTesting/modifiedModel.bngl','w') as bnglFile:
    bnglFile.write(bngl_model)
    actionsBlock = addSimParamters('ode',True,['EGF(erbb)'])
    bnglFile.write(actionsBlock)


