from indra.sources import trips
from indra.explanation import model_checker


#Test Stmts

#sentence = 'JAK2 phosphorylates STAT3. JAK2 activates STAT3'
sentence = 'SORAFENIB inhibits FLT3. SORAFENIB inhibits AKT1'
tp = trips.process_text(sentence)
testStmt = tp.statements

mc = model_checker.ModelChecker(model=smallModelContext,statements=testStmt)
mc.get_im()
results = mc.check_model()

# Inhibition(SORAFENIB(), FLT3()),
# Activation(FLT3(mods: (phosphorylation, Y, 589), (phosphorylation, Y, 591), (phosphorylation, Y, 597), (phosphorylation, Y, 599)), AKT1()),
# Rule('SORAFENIB_deactivates_FLT3_activity', SORAFENIB() + FLT3(activity='active') >> SORAFENIB() + FLT3(activity='inactive'), kf_sf_act_2),



###############


expStmts = ac.load_statements('sorafenibExperimentalStmts.pkl')
mc = model_checker.ModelChecker(model=smallModelContext,statements=expStmts)
results = mc.check_model()

mc2 = model_checker.ModelChecker(model=largeModel,statements=expStmts)
results2 = mc2.check_model()
