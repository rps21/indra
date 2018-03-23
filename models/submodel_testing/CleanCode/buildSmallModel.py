import removeDimersMutants as rdm
import combinePhosphorylationSites as cps
import enforceCascadeContext as ecc
import addDephosphorylationAndDegradationRules as add
import findOligomers as fo
from indra.tools import assemble_corpus as ac 
from indra.assemblers import PysbAssembler
from indra.statements import *
import pysb

def buildSmallModel(stmts):
    stmts = rdm.remove_mutations(stmts)
    stmts = rdm.remove_dimers(stmts)
    stmts = rdm.remove_bad_translocations(stmts)
    stmts, uplist, downlist = ecc.add_all_af(stmts)     
    stmts = ecc.run_mechlinker_step_reduced(stmts, uplist, downlist)
    stmts = cps.coarse_grain_phos(stmts)    #Can have extra input here, generic = T/F
    stmts = add.add_dephosphorylations(stmts)
    stmts = add.add_deegradations(stmts)
    cycles = fo.findCycles(stmts)

    pa = PysbAssembler('two_step')
    pa.add_statements(stmts)
    model = pa.make_model()
     
    bngl_model = pysb.export.export(model,'bngl')
    bngl_file = open('smallModel.bngl','w')  #Add optional input for naming output file?
    bngl_file.write(bngl_model)
    bngl_file.close()

    return stmts, cycles

stmts = ac.load_statements('/home/bobby/Dropbox/Sorger_Lab/indra/models/submodel_testing/TestingCode/debugging/fallahi_eval_pysb_stmts_updated.pkl')
egfrGenes = ['EGF','EGFR','GRB2','SOS1','KRAS','BRAF','MAP2K1','MAPK1']
egfrStmts = ac.filter_gene_list(stmts,egfrGenes,'all')
#egfrStmts = ac.filter_by_type(egfrStmts,Autophosphorylation,invert=True)
newStmts, cycles = buildSmallModel(egfrStmts)



#autophosphoyrlation - not being caught in phosphorylation adjustment. Keep? Throw away?
#Additionally, activeforms with phos get precedent. If there's no accompanying phosphorylation statement, problem
#EGFR hits both these 


#DONE
#generic phosphorylations with no residue sneaking through
#   Comes from dephosphorylation. Probably doesn't like changed residue names. 
#   Had to edit statements.py to not check if a residue is valid. This seems super dangerous, need a better workaround.
#Remove translocations, minimially ones with a None component
#check if all cycles are being returned, or only one 
