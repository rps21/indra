import pickle
import tempfile
from indra.tools.small_model_tools.modelScope.buildPriorModel import build_prior
from indra.statements import *
from indra.tools import assemble_corpus as ac
from indra.tools.small_model_tools.modelScope.buildDrugTreatmentStatements import buildExperimentalStatements
from indra.tools.small_model_tools.modelScope.buildDrugTreatmentStatements import buildDrugTargetStmts
from indra.tools.small_model_tools.modelScope.buildUseCaseModel import expandModel
from indra.tools.small_model_tools.modelScope.addBNGLParameters import addObservables
from indra.tools.small_model_tools.modelScope.addBNGLParameters import addSimParamters
from indra.tools.small_model_tools.modelContext import extraModelReductionTools as ex
from indra.assemblers import PysbAssembler
import pysb


#Receptor Names = [ERBB2, ERBB3, ERBB4, EGFR, MET, AXL, Tyro3]

#Changes in expression levels of receptors at DAY4 (fold change relative to untreated cells): 
#ERBB2 - 3.5
#ERBB3 - 1.5
#ERBB4 - 3.5
#EGFR - 0.4
#MET - 0.3
#AXL - 0.3
#Tyro3 - 0.5
#I only know that inhibiting EGFR, ERBB2, ERBB3, and ERBB4 is important for the survival of these cells. 
#My hypothesis is that there is signaling through the ERBB2/3 dimer or ERBB3/4 dimer, but it could of course be wrong. 




#User Input
drugTargets = ['BRAF']
modifiedNodes = ['ERBB2','EGFR','ERBB3','MET']
#ERBB4 failed. 
otherNodes =['KRAS','MAP2K1','SOS1','GRB2']
model_genes = drugTargets + modifiedNodes + otherNodes

#Add drug-target interactions 
drug = 'VEMURAFENIB'
#Inhibitor target = BRAFV600E mutant; does not inhibit formation and signaling through BRAF/CRAF dimers.
#drugSentences = 'VEMURAFENIB dephosphorylates BRAF V600E.' 
expSentences = 'VEMURAFENIB degrades EGFR. VEMURAFENIB degrades MET. VEMURAFENIB transcribes ERBB2. VEMURAFENIB transcribes ERBB3.'# VEMURAFENIB transcribes ERBB4.'

mutations = {'BRAF': [('Val', '600', 'Glu'),('V', '600', 'E')]}
#drugSentences = 'VEMURAFENIB dephosphorylates BRAF V600E.' 

drugSentences = 'BRAF V600E bound to VEMURAFENIB is inactive. VEMURAFENIB BINDS BRAF V600E' 
drug_stmts = buildDrugTargetStmts(drugSentences)




def unify_mutations(stmts):
    mutation_stmts = []
    output_stmts = []
    for st in stmts:
        if any(list(map(lambda obj: obj.mutations, st.agent_list()))): 
            mutation_stmts.append(st)
        else:
            output_stmts.append(st)

    for st in mutation_stmts:
        for ag in st.agent_list():
            for mut in ag.mutations:
                if mut.residue_from.lower() in amino_acids_reverse:
                    mut.residue_from = amino_acids_reverse.get(mut.residue_from.lower())
                    mut.residue_to = amino_acids_reverse.get(mut.residue_to.lower())

    output_stmts = output_stmts + mutation_stmts
    return output_stmts






def buildExpObs(expSentences):
    exp_stmts = buildExperimentalStatements(expSentences)
    expObservations = {}
    for st in exp_stmts:
        if isinstance(st,Phosphorylation):
            expObservations[st.sub.name] = ['dephosphorylation',st] #works for mods 
        elif isinstance(st,Dephosphorylation):
            expObservations[st.sub.name] = ['phosphorylation',st] #works for mods 
        elif isinstance(st,IncreaseAmount):
            expObservations[st.obj.name] = ['decreaseamount',st] #works for mods 
        elif isinstance(st,DecreaseAmount):
            expObservations[st.obj.name] = ['increaseamount',st] #works for mods 
    return expObservations


expObservations = buildExpObs(expSentences)     #Need to save this or separate it out

#Expand model to match experimental observations
def buildModel(prior_model_stmts,expObservations,drug,drugTargets,drug_stmts,otherNodes):
    #prior_model_stmts = ac.load_statements('modelPrior.pkl')
    initialStmts = prior_model_stmts
    initialNodes = []
    for st in prior_model_stmts:
        for ag in st.agent_list():
            if ag.name not in initialNodes:
                initialNodes.append(ag.name)


    modelStmts = expandModel(expObservations,drug,drugTargets,drug_stmts,initialStmts,initialNodes,otherNodes)
    return modelStmts




###############################################################





#Run 
prior_model_stmts = ac.load_statements('modelPrior_large.pkl')
#prior_model_stmts,expObservations = buildAndCleanPrior(modifiedNodes,model_genes,mutations,expSentences)
prior_model_stmts = unify_mutations(prior_model_stmts)
modelStmts  = buildModel(prior_model_stmts,expObservations,drug,drugTargets,drug_stmts,otherNodes)


#Temp inhibitory hack
finalModelStmts = []
for st in modelStmts:
    if 'phos_inhib' not in str(st):
        finalModelStmts.append(st)


ac.dump_statements(fname='vem_model_stmts.pkl',stmts=finalModelStmts)


#def buildAndCleanPrior(drugTargets,modifiedNodes,model_genes,drug,drugSentences,expSentences):
#def expandModel(prior_model_stmts,expObservations,drug,drugTargets,otherNodes):


###############################################################
#Build Prior 

#expObservations = {}
#for st in exp_stmts:
#    if isinstance(st,Modification):
#        expObservations[st.sub.name] = [st._get_mod_condition().mod_type,st] #works for mods 
#    elif isinstance(st,RegulateAmount):
#        expObservations[st.obj.name] = [type(st).__name__.lower(),st] #works for mods 
##Build function for inverse of statements

def buildAndCleanPrior(modifiedNodes,model_genes,mutations,expSentences):

    exp_stmts = buildExperimentalStatements(expSentences)

    expObservations = {}
    for st in exp_stmts:
        if isinstance(st,Phosphorylation):
            expObservations[st.sub.name] = ['dephosphorylation',st] #works for mods 
        elif isinstance(st,Dephosphorylation):
            expObservations[st.sub.name] = ['phosphorylation',st] #works for mods 
        elif isinstance(st,IncreaseAmount):
            expObservations[st.obj.name] = ['decreaseamount',st] #works for mods 
        elif isinstance(st,DecreaseAmount):
            expObservations[st.obj.name] = ['increaseamount',st] #works for mods 


    #reach_stmts = './indraReading/raw_stmts/reach_output.pkl'
    #trips_stmts = './indraReading/raw_stmts/trips_output.pkl'
    #additional_stmts_files = [reach_stmts,trips_stmts]
    model_types = [Phosphorylation,Dephosphorylation,ActiveForm,IncreaseAmount,DecreaseAmount,Complex] #Gef,GtpActivation]#Check sos and raf stmts 
    prior_stmts = build_prior(model_genes,model_types,dbs=True,additional_stmts_files=[])

    #Prior reduction

    #Map to model_types
    mods_to_keep = list(map(lambda obj: obj.__name__.lower(), model_types))
    stmts_to_remove = []
    for st in prior_stmts:
        full_mods = list(map(lambda obj: obj.mods, st.agent_list())) #gives mods for a stmt, in list of lists. Actually need mod types. [(phosphorylation, Y, 589)]
        flattened_mods = list(itertools.chain.from_iterable(full_mods))
        #mod_types = list(map(lambda obj: obj.mod_type, mods) #mod types in lol 
        mod_types = list(map(lambda obj: obj.mod_type, flattened_mods))
        if any([el for el in mod_types if el not in mods_to_keep]):
            stmts_to_remove.append(st)
    prior_model_stmts = [st for st in prior_stmts if st not in stmts_to_remove]


    stmts_to_remove = []
    for st in prior_model_stmts:
        for ag in st.agent_list():
            for mod in ag.mods:
                if mod.mod_type not in mods_to_keep:
                    stmts_to_remove.append(st)

    prior_model_stmts = [st for st in prior_model_stmts if st not in stmts_to_remove]


#    prior_model_stmts = ex.removeMutations(prior_model_stmts)
#    prior_model_stmts = ex.removeDimers(prior_model_stmts)
#    prior_model_stmts = prior_stmts + drug_stmts
    prior_model_stmts = ac.filter_mutation_status(prior_model_stmts,mutations,deletions=[])
    prior_model_stmts = unify_mutations(prior_model_stmts)


    #save prior model 
    ac.dump_statements(prior_model_stmts,'modelPrior_large.pkl')


    return prior_model_stmts, expObservations

###########################################################################




###########################################################################

#Build and write bngl model 

pa = PysbAssembler()
pa.add_statements(finalModelStmts)
pysbModel = pa.make_model()
#WARNING: [2018-08-29 16:12:59] indra/pysb_assembler - HIF1A transcribes itself, skipping???

model_obs = addObservables(pysbModel,bound=True)


bngl_model_filter = pysb.export.export(model_obs,'bngl')
bngl_file_filter = open('vemurafenib_bngl_model.bngl','w')
bngl_file_filter.write(bngl_model_filter+'\n')
bngl_file_filter.close()
model_actions = addSimParamters(method='ode',equil=True,equilSpecies=['%s()' % drug],viz=True)
bngl_file_filter = open('vemurafenib_bngl_model.bngl','a')
bngl_file_filter.write(model_actions)
bngl_file_filter.close()



