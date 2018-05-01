from indra.sources import trips
from indra.explanation import model_checker
import indra.tools.assemble_corpus as ac
import pickle
from indra.assemblers import PysbAssembler


###############


#expStmts = ac.load_statements('sorafenibExperimentalStmts.pkl')
#mc = model_checker.ModelChecker(model=smallModelContext,statements=expStmts)
#results = mc.check_model()

#mc2 = model_checker.ModelChecker(model=largeModel,statements=expStmts)
#results2 = mc2.check_model()



##############################
##More Sorafenib testing
#In [34]: fltStmts
#Out[34]: 
#[Activation(FLT3LG(), FLT3()),
# Activation(FLT3(mods: (phosphorylation, Y, 589), (phosphorylation, Y, 591), (phosphorylation, Y, 597), (phosphorylation, Y, 599)), AKT1()),
# ActiveForm(FLT3(mods: (phosphorylation, Y, 589), (phosphorylation, Y, 591), (phosphorylation, Y, 597), (phosphorylation, Y, 599)), kinase, True),
# Complex(FLT3LG(), FLT3()),
# Inhibition(SORAFENIB(), FLT3()),
# Phosphorylation(FLT3(mods: (phosphorylation, Y, 589), (phosphorylation, Y, 591), (phosphorylation, Y, 597), (phosphorylation, Y, 599)), FLT3(), Y, 589),
# Phosphorylation(FLT3(mods: (phosphorylation, Y, 589), (phosphorylation, Y, 591), (phosphorylation, Y, 597), (phosphorylation, Y, 599)), FLT3(), Y, 591),
# Phosphorylation(FLT3(mods: (phosphorylation, Y, 589), (phosphorylation, Y, 591), (phosphorylation, Y, 597), (phosphorylation, Y, 599)), FLT3(), Y, 597),
# Phosphorylation(FLT3(mods: (phosphorylation, Y, 589), (phosphorylation, Y, 591), (phosphorylation, Y, 597), (phosphorylation, Y, 599)), FLT3(), Y, 599)]

#In [38]: sorafTestStmts
#Out[38]: 
#[Inhibition(SORAFENIB(), FLT3()),
# Activation(FLT3(mods: (phosphorylation, Y, 589), (phosphorylation, Y, 591), (phosphorylation, Y, 597), (phosphorylation, Y, 599)), AKT1())]

#sorafTestStmts = []
#sorafTestStmts.append(fltStmts[4])
#sorafTestStmts.append(fltStmts[1])

#pa = PysbAssembler()
#pa.add_statements(sorafTestStmts)
#testModel = pa.make_model()
#with open('testModel.pkl','wb') as f:
#    pickle.dump(testModel,f)

with open('smallPYSBModel_contextChanges.pkl','rb') as f:
    model = pickle.load(f)

sentence = 'SORAFENIB inhibits FLT3. SORAFENIB inhibits AKT1'
tp = trips.process_text(sentence)
testStmt = tp.statements

mc = model_checker.ModelChecker(model=model,statements=testStmt)
results = mc.check_model()

#Still Failing
# Dephosphorylation(SORAFENIB(), FLT3(), phos_act)
# ActiveForm(FLT3(mods: (phosphorylation, phos_act)), kinase, True),
# Activation(FLT3(mods: (phosphorylation, phos_act)), AKT1()),

#testStmts
#[Inhibition(SORAFENIB(), FLT3()), Inhibition(SORAFENIB(), AKT1())]

#########
#Testing Phos Site
fltStmts2 = ac.filter_gene_list(smallModelStmts,['FLT3'],'one')
smallModelStmts_contextChanges.append(fltStmts2[1])
smallModelStmts_contextChanges.append(fltStmts2[5])

pa = PysbAssembler()
pa.add_statements(smallModelStmts_contextChanges)
smallModelContext = pa.make_model()
with open('smallPYSBModel_contextChanges_phosTest.pkl','wb') as f:
    pickle.dump(smallModelContext,f)

#Model Monomer for flt3 
# Monomer('FLT3', ['activity', 'phos_act', 'flt3lg', 'phospho', 'Y589'], {'activity': ['inactive', 'active'], 'Y589': ['u', 'p'], 'phospho': ['u', 'p'], 'phos_act': ['u', 'p']}),









#with open('smallPYSBModel_contextChanges.pkl','rb') as f:
with open('smallPYSBModel_contextChanges.pkl','rb') as f:
    model = pickle.load(f)

sentence = 'SORAFENIB inhibits FLT3. SORAFENIB inhibits AKT1. Sorafenib dephosphorylates FLT3. Sorafenib dephosphorylates AKT1'
tp = trips.process_text(sentence)
testStmt = tp.statements

mc = model_checker.ModelChecker(model=model,statements=testStmt)
results = mc.check_model()


# Dephosphorylation(SORAFENIB(), FLT3(), phos_act),
# ActiveForm(FLT3(mods: (phosphorylation, phos_act)), kinase, True),
# Activation(FLT3(mods: (phosphorylation, phos_act)), AKT1()),




########################

#modelSentences = 'FLT3 activates AKT1. Sorafenib inhibits FLT3.'   #First works, second fails
#modelSentences = 'Active FLT3 activates AKT1. Sorafenib inhibits FLT3.' #Works
modelSentences = 'Phosphorylated FLT3 is active. Active FLT3 activates AKT1. Sorafenib inhibits FLT3.' #Both Fail
#modelSentences = 'Phosphorylated FLT3 is active. Active FLT3 activates AKT1. Sorafenib dephosphorylates FLT3.'  #Works
#modelSentences = 'Phosphorylated FLT3 is active. Phosphorylated FLT3 activates AKT1. Sorafenib dephosphorylates FLT3.'  #Works
#modelSentences = 'Phosphorylated FLT3 activates AKT1. Sorafenib dephosphorylates FLT3.'  #Error - No 'activity' in FLT3
#modelSentences = 'Phosphorylated FLT3 activates AKT1. Sorafenib dephosphorylates FLT3. FLT3LG activates FLT3'  #First fails, second works
tp = trips.process_text(modelSentences)
pa = PysbAssembler()
pa.add_statements(tp.statements)
testModel = pa.make_model()


testSentences = 'SORAFENIB inhibits FLT3. SORAFENIB inhibits AKT1'
tp = trips.process_text(testSentences)
testStmts = tp.statements

mc = model_checker.ModelChecker(model=testModel,statements=testStmts)
results = mc.check_model()











