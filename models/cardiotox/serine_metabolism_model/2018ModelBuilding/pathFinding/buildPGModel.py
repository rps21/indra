from indra.assemblers import PysbAssembler
from indra.assemblers import SifAssembler
import paths_graph as pg
from indra.tools import assemble_corpus as ac

#with open('smallModelStmts_contextChanges.pkl','rb') as f:
#    model = pickle.load(f)

model_stmts = ac.load_statements('smallModelStmts_contextChanges.pkl')


sa = SifAssembler(model_stmts)
sifModel = sa.make_model(use_name_as_key=True, include_mods=True,include_complexes=True)
testSifModel = sa.print_model(include_unsigned_edges=True)
sa.save_model(fname='testSifModel.sif',include_unsigned_edges=True)

with open('testSifModel.sif','r') as f:
    sifModelText = f.readlines()

newSifModelText = []
#For pg - 0 = positive, 1 = negative
#In raw sif - 1 = positive, 0 = neutral, -1 = negative
for line in sifModelText:
    if ' 1 ' in line:
        newline = line.replace(' 1 ', ' 0 ')
    elif ' -1 ' in line:
        newline = line.replace(' -1 ', ' 1 ')
    else:
        newline = line    
    newSifModelText.append(newline)


with open('newSifModel.sif','w') as f:
    for line in newSifModelText:
        f.write(line)


g = pg.api.load_signed_sif('newSifModel.sif')
source = 'FLT3LG'
target = 'AKT1'
target_polarity=0
length=1


paths = pg.api.sample_paths(g,source,target,max_depth=2,num_samples=2,cycle_free=True,signed=False,target_polarity=0) #Works well
#paths is a list of tuples of strings.
#Could feed strings back into a filtering by genes routine
#Could start from bigger, less biases stmt set. Sample paths, build model with strict filtering from paths, use model checker on model. If fails, new set of sampled paths. 

#to investigate:
#class CombinedCFPG(object):
#    """Combine a set of CFPGs for different lengths into a single super-CFPG.



#cfpg = pg.CFPG.from_graph(g,source,target,target_polarity,length)   #doesn't seem to work 


#large stmt set 
#may need to figure out how to run context tweaks on large set of stmts 
large_model_stmts = ac.load_statements('finalLargeModelStmts.pkl')


sa = SifAssembler(large_model_stmts)
sifModelLarge = sa.make_model(use_name_as_key=True, include_mods=True,include_complexes=True)
sa.save_model(fname='newSifModelLarge.sif',include_unsigned_edges=True)

with open('newSifModelLarge.sif','r') as f:
    sifModelText = f.readlines()

newSifModelText = []
#For pg - 0 = positive, 1 = negative
#In raw sif - 1 = positive, 0 = neutral, -1 = negative
for line in sifModelText:
    if ' 1 ' in line:
        newline = line.replace(' 1 ', ' 0 ')
    elif ' -1 ' in line:
        newline = line.replace(' -1 ', ' 1 ')
    else:
        newline = line    
    newSifModelText.append(newline)


with open('newSifModelLarge.sif','w') as f:
    for line in newSifModelText:
        f.write(line)


g = pg.api.load_signed_sif('newSifModelLarge.sif')
source = 'FLT3'
target = 'AKT1'
target_polarity=0
length=1

#import timeit
#timeit.timeit(paths = pg.api.sample_paths(g,source,target,max_depth=10,num_samples=10,cycle_free=True,signed=False,target_polarity=0))
paths = pg.api.sample_paths(g,source,target,max_depth=10,num_samples=10,cycle_free=True,signed=False,target_polarity=0)
#check sif directionality. May be adding lines in both directions for complexes. Finding a way around this may be tricky, but would be a big  help. 
#timeit failed - probably can't set a variable. But took roughly 25 min by indra output


#This may crash
newstmts, uplist, downlist = cs.add_all_af(largeModelStmts)     
newstmts = cs.run_mechlinker_step_reduced(newstmts, uplist, downlist)
largeModelStmts_contextChanges = ptm.coarse_grain_phos(newstmts)
ac.dump_statements(largeModelStmts_contextChanges,'newLargeModelStmts_contextChanges.pkl')


#sentence = 'SORAFENIB inhibits FLT3. SORAFENIB inhibits AKT1'
#tp = trips.process_text(sentence)
#testStmt = tp.statements

#mc = model_checker.ModelChecker(model=model,statements=testStmt)
#results = mc.check_model()







