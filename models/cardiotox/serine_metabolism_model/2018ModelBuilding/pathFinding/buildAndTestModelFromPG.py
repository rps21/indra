from indra.assemblers import PysbAssembler
from indra.assemblers import SifAssembler
import paths_graph as pg
from indra.tools import assemble_corpus as ac
from indra.tools.small_model_tools import enforceCascadeContext as cs
from indra.tools.small_model_tools import combinePhosphorylationSites as ptm
from indra.sources import trips
from indra.explanation import model_checker


g = pg.api.load_signed_sif('newSifModel.sif')
source = 'FLT3LG'
target = 'AKT1'
target_polarity=0
length=1

paths = pg.api.sample_paths(g,source,target,max_depth=2,num_samples=2,cycle_free=True,signed=False,target_polarity=0) #Works well


testPath = list(paths[0])
largeModelStmts = ac.load_statements('finalLargeModelStmts.pkl')
testModel = ac.filter_gene_list(largeModelStmts,testPath,'all')


newstmts, uplist, downlist = cs.add_all_af(testModel)     
newstmts = cs.run_mechlinker_step_reduced(newstmts, uplist, downlist)
testModel_contextChanges = ptm.coarse_grain_phos(newstmts)
#ac.dump_statements(testModel_contextChanges,'testModel_contextChanges.pkl')

pa = PysbAssembler()
pa.add_statements(testModel_contextChanges)
testModelPySB = pa.make_model()


sentence = 'FLT3LG activates AKT1.'
tp = trips.process_text(sentence)
testStmt = tp.statements

mc = model_checker.ModelChecker(model=testModelPySB,statements=testStmt)
results = mc.check_model()

results[0][1].path_found #True or False
