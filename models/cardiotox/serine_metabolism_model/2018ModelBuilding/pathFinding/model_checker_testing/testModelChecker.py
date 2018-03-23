from indra.sources import trips
from indra.explanation import model_checker
import indra.tools.assemble_corpus as ac
import pickle

with open('smallPYSBModel_contextChanges.pkl','rb') as f:
    model = pickle.load(f)

modelStmts = ac.load_statements('smallModelStmts_contextChanges.pkl')

#sentence = 'JAK2 phosphorylates STAT3. JAK2 activates STAT3'
sentence = 'SORAFENIB inhibits FLT3. SORAFENIB inhibits AKT1'
tp = trips.process_text(sentence)
testStmt = tp.statements

mc = model_checker.ModelChecker(model=model,statements=testStmt)
results = mc.check_model()

# fltStmts = ac.filter_gene_list(modelStmts,['FLT3'],'one')

# Inhibition(SORAFENIB(), FLT3()),
# Activation(FLT3(mods: (phosphorylation, Y, 589), (phosphorylation, Y, 591), (phosphorylation, Y, 597), (phosphorylation, Y, 599)), AKT1()),
# Rule('SORAFENIB_deactivates_FLT3_activity', SORAFENIB() + FLT3(activity='active') >> SORAFENIB() + FLT3(activity='inactive'), kf_sf_act_2),



###############


expStmts = ac.load_statements('sorafenibExperimentalStmts.pkl')
mc = model_checker.ModelChecker(model=smallModelContext,statements=expStmts)
results = mc.check_model()

mc2 = model_checker.ModelChecker(model=largeModel,statements=expStmts)
results2 = mc2.check_model()



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

sorafTestStmts = []
sorafTestStmts.append(fltStmts[4])
sorafTestStmts.append(fltStmts[1])



pa = PysbAssembler()
pa.add_statements(sorafTestStmts)
testModel = pa.make_model()
with open('sorafTestPYSBModel.pkl','wb') as f:
    pickle.dump(testModel,f)

sentence = 'SORAFENIB inhibits FLT3. SORAFENIB inhibits AKT1'
tp = trips.process_text(sentence)
testStmt = tp.statements

mc = model_checker.ModelChecker(model=testModel,statements=testStmt)
results = mc.check_model()

###

sentence = 'FLT3 activates AKT1'
tp = trips.process_text(sentence)
sorafTestStmts = sorafTestStmts + tp.statements

pa = PysbAssembler()
pa.add_statements(sorafTestStmts)
testModel = pa.make_model()

mc = model_checker.ModelChecker(model=testModel,statements=testStmt)
results = mc.check_model()







