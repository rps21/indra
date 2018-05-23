import pickle
import itertools
from indra.tools import assemble_corpus as ac
from indra.assemblers import PysbAssembler
from indra.mechlinker import MechLinker
from indra.statements import *
from indra.sources import signor
from indra.sources import biopax
from indra.sources import trips
from indra.tools.gene_network import GeneNetwork
from indra.preassembler import Preassembler


def normalize_active_forms(stmts):
    af_stmts = ac.filter_by_type(stmts, ActiveForm)
    relevant_af_stmts = []
    for stmt in af_stmts:
        if (not stmt.agent.mods) and (not stmt.agent.mutations):
            continue
        relevant_af_stmts.append(stmt)
    print('%d relevant ActiveForms' % len(relevant_af_stmts))
    non_af_stmts = ac.filter_by_type(stmts, ActiveForm, invert=True)
    af_stmts = ac.run_preassembly(relevant_af_stmts)
    stmts = af_stmts + non_af_stmts
    return stmts

model_genes = [
    'FLT3','STAT3', 'JAK2', 'PKM', 'HIF1A', 'MYC',
    'CDKN1A','KDR','VEGFR1','FLT3','FLT1','PDGF','FLT3LG','PDGFA','PDGFRA',
    'JAK2','FLT3','ABL1','PTPN6','RPS6','RPS6KB1','PIK3CA','AKT1','AURKA'
    ]


#biopax_stmts = build_prior(model_genes)
#ac.dump_statements(biopax_stmts,'model_stmts/biopax_output.pkl')


#sp = signor.SignorProcessor()
#signorStmts = sp.statements
#filteredSignorStmts = ac.filter_gene_list(signorStmts, model_genes, policy='one',save='signor_stmts_list.pkl')



#raw_stmts = ac.load_statements('heart_reading_stmts.pkl')
#stmts = ac.filter_no_hypothesis(raw_stmts)
#stmts = ac.filter_gene_list(stmts, model_genes, policy='one',
#                            save='filter_gene_list.pkl')


reach_stmts = ac.load_statements('model_stmts/filter_gene_list.pkl')
trips_stmts = ac.load_statements('model_stmts/trips_output.pkl')
biopax_stmts = ac.load_statements('model_stmts/biopax_output.pkl')
filteredSignorStmts = ac.load_statements('model_stmts/signor_stmts_list.pkl')
stmts = reach_stmts + biopax_stmts + trips_stmts + filteredSignorStmts  #Probably need to up date signor stmts 




def cleanStatements(stmts):
    stmts = ac.map_grounding(stmts)
    stmts = ac.filter_grounded_only(stmts, save='intermediateStmts/grounded.pkl')
    stmts = ac.filter_genes_only(stmts)
    stmts = ac.filter_human_only(stmts)
    stmts = ac.run_preassembly(stmts, save='intermediateStmts/preassembled.pkl')
    stmts = ac.filter_belief(stmts, 0.90)
    stmts = ac.filter_top_level(stmts)
    stmts = ac.filter_enzyme_kinase(stmts)
    stmts = ac.filter_mod_nokinase(stmts)
    stmts = ac.filter_transcription_factor(stmts)   #Any downside?
    # Simplify activity types
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.reduce_activities()
    ml.gather_modifications()
    ml.reduce_modifications()
    stmts = normalize_active_forms(ml.statements)
    # Replace activations when possible
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.replace_activations()
    # Require active forms
    ml.require_active_forms()
    num_stmts = len(ml.statements)
    """
    while True:
        # Remove inconsequential PTMs
        ml.statements = ac.filter_inconsequential_mods(ml.statements,
                                                       get_mod_whitelist())
        ml.statements = ac.filter_inconsequential_acts(ml.statements,
                                                       get_mod_whitelist())
        if num_stmts <= len(ml.statements):
            break
        num_stmts = len(ml.statements)
    """
    stmts = ml.statements
    return stmts


#New filter, because of errors that arrise from none agents
stmts = [st for st in stmts if None not in st.agent_list()]

largeModelStmts = cleanStatements(stmts)
smallModelRawStmts = ac.filter_gene_list(stmts,model_genes,'all')
smallModelStmts = cleanStatements(smallModelRawStmts)

###############################################
#Add context modifications to stmt lists
#Need more work to assess necessity/affect on sif file generation 
#################################################

#Coarse grain phos and add binding context, interpretting as activationsfrom indra.tools.small_model_tools import enforceCascadeContext as cs
from indra.tools.small_model_tools import enforceCascadeContext as cs
from indra.tools.small_model_tools import combinePhosphorylationSites as ptm

#newstmts, uplist, downlist = cs.add_all_af(largeModelStmts)     
#newstmts = cs.run_mechlinker_step_reduced(newstmts, uplist, downlist)
#largeModelStmts_contextChanges = ptm.coarse_grain_phos(newstmts)

newstmts, uplist, downlist = cs.add_all_af(smallModelStmts)     
smallModelStmts_contextChanges = cs.run_mechlinker_step_reduced(newstmts, uplist, downlist)
smallModelStmts_contextChanges = ptm.coarse_grain_phos(newstmts)






#######################################################
#Add manual stmts for Sorafenib interactions 
#May want to check these for complex vs inhibition
########################################################


#Add statements for Sorafenib Targets
drugTargetStmts = ac.load_statements('sorafenibTargetStmts.pkl')
largeModelStmts = largeModelStmts + drugTargetStmts
#largeModelStmts_contextChanges = largeModelStmts_contextChanges + drugTargetStmts
smallModelStmts_contextChanges = smallModelStmts_contextChanges + drugTargetStmts



#####################################################
#Final processing/cleanup
#Deduplication, removing inconsequential mods
######################################################

largeModelStmts = Preassembler.combine_duplicate_stmts(largeModelStmts)
#largeModelStmts_contextChanges = Preassembler.combine_duplicate_stmts(largeModelStmts_contextChanges)
#largeModelStmts_contextChanges = ac.filter_inconsequential_mods(largeModelStmts_contextChanges)

smallModelStmts_contextChanges = Preassembler.combine_duplicate_stmts(smallModelStmts_contextChanges)
smallModelStmts_contextChanges = ac.filter_inconsequential_mods(smallModelStmts_contextChanges)



#######################
# Save stmts 
#######################

##Save four final statment lists
#ac.dump_statements(largeModelStmts,'finalLargeModelStmts.pkl')
##ac.dump_statements(largeModelStmts,'largeModelStmts_contextChanges.pkl')
#ac.dump_statements(smallModelStmts_contextChanges,'smallModelStmts_contextChanges.pkl')


#############################
#Build and save PySB models 
##############################


#pa = PysbAssembler()
#pa.add_statements(largeModelStmts)
#largeModel = pa.make_model()
#with open('largePYSBModel.pkl','wb') as f:
#    pickle.dump(largeModel,f)


##pa = PysbAssembler()
##pa.add_statements(largeModelStmts_contextChanges)
##largeModel = pa.make_model()
##with open('largePYSBModel_contextChanges.pkl','wb') as f:
##    pickle.dump(largeModel,f)

#pa = PysbAssembler()
#pa.add_statements(smallModelStmts_contextChanges)
#smallModelContext = pa.make_model()
#with open('smallPYSBModel_contextChanges.pkl','wb') as f:
#    pickle.dump(smallModelContext,f)

