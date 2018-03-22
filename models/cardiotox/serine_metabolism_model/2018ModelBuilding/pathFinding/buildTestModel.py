from indra.assemblers import PysbAssembler
from indra.sources import trips

#testModelSentences = 'EGF activates EGFR. Active EGFR activates GRB2. GRB2 binds SOS1. SOS1 is active when bound to GRB2'
testModelSentences = 'EGF binds EGFR. EGFR binds GRB2. GRB2 binds SOS1.'
trips_processor = trips.process_text(testModelSentences)
testModelStmts = trips_processor.statements 


pa = PysbAssembler()
pa.add_statements(testModelStmts)
testModel = pa.make_model()
#with open('testPYSBModel.pkl','wb') as f:
#    pickle.dump(testModel,f)


from indra.explanation import model_checker
sent_to_check = 'EGF activates SOS1. EGF activates GRB2.'
trips_processor = trips.process_text(sent_to_check)
stmt_to_check = trips_processor.statements
stmt_to_check



mc = model_checker.ModelChecker(model=testModel,statements=stmt_to_check)
results = mc.check_model()



#/home/bobby/Dropbox/Sorger_Lab/indra/indra/explanation/model_checker.py in _check_regulate_activity(self, stmt, max_paths, max_path_length)
#    355         # appropriate type
#    356         obs_names = self.stmt_to_obs[stmt]
#--> 357         for obs_name in obs_names:
#    358             return self._find_im_paths(subj_mp, obs_name, target_polarity,
#    359                                        max_paths, max_path_length)

#TypeError: 'NoneType' object is not iterable

#ipdb> self.stmt_to_obs
#{Activation(EGF(), GRB2()): ['GRB2_activity_active_obs'], Activation(EGF(), SOS1()): None}

#In Pysb model SOS1 has 'grb2' domain but no 'activity' domain  



##################################################



from indra.tools.small_model_tools import enforceCascadeContext as cs
from indra.tools.small_model_tools import combinePhosphorylationSites as ptm

newstmts, uplist, downlist = cs.add_all_af(testModelStmts)     
final_stmts = cs.run_mechlinker_step_reduced(newstmts, uplist, downlist)

pa = PysbAssembler()
pa.add_statements(final_stmts)
testModel2 = pa.make_model()

mc2 = model_checker.ModelChecker(model=testModel2,statements=stmt_to_check)
results = mc2.check_model()

#stmts_to_use = deepcopy(final_stmts)
#newStmts = ptm.coarse_grain_phos(stmts_to_use)


#ISSUES WITH OLD CODE
#Repeatedly adding AF statements. Should simplify these within the base code
#GRB2 gets EGFR in it's active form, doesn't translate to GRB2-SOS1 binding. 
# ActiveForm(GRB2(bound: [EGFR, True]), kinase, True),
# Complex(GRB2(), SOS1()),

