from indra.sources import trips
from indra.explanation import model_checker


#Test Stmts

sentence = 'JAK2 phosphorylates STAT3'
tp = trips.process_text(sentence)
testStmt = tp.statements

mc = model_checker.ModelChecker(model=smallModel,statements=testStmt)
mc.get_im()


results = mc.check_model()
#Returns list of tuples, likely 1 tupler per statement, current list has len 1 for single statement


expStmts = ac.load_statements('sorafenibExperimentalStmts.pkl')
mc = model_checker.ModelChecker(model=smallModel,statements=expStmts)
results = mc.check_model()

mc2 = model_checker.ModelChecker(model=largeModel,statements=expStmts)
results2 = mc2.check_model()
