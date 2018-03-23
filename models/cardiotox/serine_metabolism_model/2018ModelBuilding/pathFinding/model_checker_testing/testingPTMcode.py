import pickle
import itertools
from indra.assemblers import PysbAssembler
from indra.mechlinker import MechLinker
from indra.statements import *
from indra.tools.gene_network import GeneNetwork
from indra.preassembler import Preassembler
from indra.explanation import model_checker
import indra.tools.assemble_corpus as ac
from indra.tools.small_model_tools import enforceCascadeContext as cs
from indra.tools.small_model_tools import combinePhosphorylationSites as ptm


#From buildPathfindingModel_v2.py
#ac.dump_statements(smallModelStmts,'model_checker_testing/smallModel_noContext.pkl')

#modelStmts = ac.load_statements('model_checker_testing/smallModel_noContext.pkl')

#newstmts, uplist, downlist = cs.add_all_af(modelStmts)     
#modelStmts_cascade = cs.run_mechlinker_step_reduced(newstmts, uplist, downlist)
#modelStmts_cascade_ptm = ptm.coarse_grain_phos(modelStmts_cascade)

#testStmts = ac.filter_gene_list(modelStmts_cascade_ptm,['FLT3LG','FLT3','JAK2','STAT3','ABL1'],'all')
#ac.dump_statements(testStmts,'model_checker_testing/ptm_testing.pkl')
testStmts = ac.load_statements('model_checker_testing/ptm_testing.pkl')

#Three Problems
#1. Homodimers Fail
#2. Activity conflicts:
#     ActiveForm(ABL1(mods: (phosphorylation, phos_act)), activity, False),
#     ActiveForm(ABL1(mods: (phosphorylation, phos_act)), kinase, True),
#3. Other rare failures, unsure of reason:
#     Phosphorylation(ABL1(mods: (phosphorylation, phos_act)), JAK2(), Y, 1007),
#4. Check that activities are correct. JAK2 looks like it could be wrong. Also:
#     ActiveForm(JAK2(mods: (phosphorylation, phos_inhib)), kinase, True),
#     ActiveForm(JAK2(mods: (phosphorylation, phos_inhib)), activity, False),
#     ActiveForm(JAK2(mods: (ubiquitination, phos_inhib)), activity, False),


#PROBLEM from source
# ActiveForm(JAK2(mods: (phosphorylation, Y, 1007)), activity, False),
# ActiveForm(JAK2(mods: (phosphorylation, Y, 1007)), activity, True),
# ActiveForm(JAK2(mods: (phosphorylation, Y, 1008)), activity, False),
# ActiveForm(JAK2(mods: (phosphorylation, Y, 1008)), activity, True),
# All from Signor 

#ActiveForm(JAK2(mods: (phosphorylation, Y, 1007)), activity, False),
#5 pieces of evidence, all from Signor
#4 look like mis readings 
#But 1 says 1007 phos allows for binding of SOCS1, ubiquitination, degradtion of JAK2 
#So no catalyical inhibition, but does decrease JAK2 activity 
#We report that SHP-2 dephosphorylates tyrosine (Tyr-1007) of Jak2 kinase, a critical recruitment site for the ubiquitin ligase-associated inhibitory protein suppressor of cytokine signaling-1 (SOCS-1), thereby contributing to Jak2 stability. Inactivation of SHP-2 function by blocking receptor/SHP-2 association or by using a catalytically inactive mutant of SHP-2 led to a marked increase in Jak2 ubiquitination/degradation, Jak2 phosphorylation on Tyr-1007, and Jak2/SOCS-1 association),


#jakStmts = ac.filter_gene_list(smallModelStmts,['JAK2'],'all')
#testStmts1 = jakStmts[1]
#testStmt2 = jakStmts[14]

#In [61]: relevantStatements
#Out[61]: 
#[Phosphorylation(JAK2(mods: (phosphorylation, Y, 1007)), JAK2(), S, 523),
# ActiveForm(JAK2(mods: (phosphorylation, Y, 1007)), activity, False)]

