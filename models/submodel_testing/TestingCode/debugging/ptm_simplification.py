from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from copy import deepcopy
import pickle
from indra.tools import assemble_corpus as ac
from indra.preassembler import Preassembler
from indra.statements import *      

#TODO: ADD OTHER MODS
#COARSE GRAIN AF PROBLEMS:
#NEW PLAN BELOW
#First get context working then try to collect all modification sites. see what it looks like.
#Check to see about specific sites that don't belong to an activation or inhibitory group. Keep? Assume activating? Drop?

##############################
##NEW THOUGHTS
##GATHER ALL SITES FIRST
########################## ###
#gather 'list' of all activating and inhibitory sites (maybe use different data structure)
#so two 'lists' for every active form with a modification
#Compare lists, see if there's overlap - potentially throw out overlapping sites, or through a warning and make a decision then
#Then check for all Y, S/T, generic phos on a single protein, join into grouped activating and inhibitory sites
#check all other statements for these sites, and rewrite, with new coarse grained site 

#potentially should do this after other context building
#make sure there is an active form for every protein first, rather than relying on ones already in stmt corpus


#Things to think about.
#Right now, combining phos sites on substrate, if same enz and substrate
#But combine sites on kinase no matter what
#Ideally should make a stricter version on sub, looser version on kin, and provide option


#Temporary hack to ignore mutations for now
#af_stmts.pop(1)

#if stmt.agent.mutations:


#loop through af stmts_in
#check is_active 
#pull agent.modifications
   #Thought, already combined af for a given protein into one statement, unless is_active is different 


#Maybe different approach
#check if prtoein has two activeforms
#if yes, handle separately, make sure coarse grained sites are well albeled
#if no, should'nt matter much. could still add a act/inhib label for clarity

#maybe do this anyway, could simplify code my removing the check for repeated active forms and model may read better
    #If I do this, I have to be sure to apply it to every time this site appears in the model. 




#collect list of proteins, list of all sites on that protein. Compare to list of sites on af stmts with that protein 
def collect_agent_name_list(stmts):
    agent_names = []
    for st in stmts:
        for ag in st.agent_list():
            if ag.name not in agent_names:        
                agent_names.append(ag.name)
    return agent_names


def find_all_ptm_sites(stmts,name):
    relevantStatements = ac.filter_gene_list(stmts,[name],'one')
    uniqueMods = []
    for st in relevantStatements:
        for ag in st.agent_list():
            if ag.name == name:
                for mod in ag.mods:
                    if not any([mod.matches(umod) for umod in uniqueMods]):
                        uniqueMods.append(mod)
    return uniqueMods



#Two sites in full set, not in af set:  (phosphorylation, T, 599), (phosphorylation, S, 602)]

def find_act_inhib_sites(stmts,name):
    relevantStatements = ac.filter_gene_list(stmts,[name],'one')
    af_stmts = ac.filter_by_type(relevantStatements,ActiveForm)
    actMods = []
    inhibMods = []
    for st in af_stmts:
        for ag in st.agent_list():
            if ag.name == name:
                for mod in ag.mods:
                    if st.is_active:
                        if not any([mod.matches(umod) for umod in actMods]):
                            actMods.append(mod)
                    else:
                        if not any([mod.matches(umod) for umod in inhibMods]):
                            inhibMods.append(mod)
    allMods = find_all_ptm_sites(stmts, name)
    neutralMods = [mod for mod in allMods if mod not in actMods and mod not in inhibMods]
    return actMods,inhibMods,neutralMods




#Now, first check all of these against sites in active form statements 

#compareing firstModSet vs mods produced from same name, but only af_stmts statement set



#allNames = collect_agent_name_list(base_statements)
#brafMods = find_all_ptm_sites(base_statements,'BRAF')   #This gives a list of all unique modifications for a single protein in the statement list


#actMods, inhibMods, neutralMods = find_act_inhib_sites(final_stmts,'BRAF')

#loop through act mods 
#S/T -> S_T_act
#Y -> Y_act 
#generic -> check if matches any above, and fold in or keep 
    #CAN'T make unique if going with this approach, these may need to be handled separately 
#Repeat with inhib and neut 

#Will want a second method that makes all act and inhib residues a single generic 


#can write this to take stmts and mod lists for a single protein 
#or can search and build those lists if given a name 
#Start with 2nd for now - architecture idea - replace_ptms as master call, calls build dict, calls find sites
#for catching unique:  not any([mod.matches(umod) for umod in uniqueMods]):
def build_ptm_dict(stmts,name):
    actMods, inhibMods, neutralMods = find_act_inhib_sites(stmts,name)
    actModsDict = {}
    inhibModsDict = {}
    neutralModsDict = {}
    for mod in actMods:
        #might need a catch for if there is no residue 
        if mod.residue == 'S' or mod.residue == 'T':
            actModsDict[str(mod)] = 'S_T_act'
        elif mod.residue == 'Y':
            actModsDict[str(mod)] = 'Y_act'
    for mod in inhibMods:
        #might need a catch for if there is no residue 
        if mod.residue == 'S' or mod.residue == 'T':
            inhibModsDict[str(mod)] = 'S_T_inhib'
        elif mod.residue == 'Y':
            inhibModsDict[str(mod)] = 'Y_inhib'
    for mod in neutralMods:
        #might need a catch for if there is no residue 
        if mod.residue == 'S' or mod.residue == 'T':
            neutralModsDict[str(mod)] = 'S_T_neutral'
        elif mod.residue == 'Y':
            neutralModsDict[str(mod)] = 'Y_neutral'
    return [actModsDict, inhibModsDict, neutralModsDict]

#Going to need a check for sites that overlap multiple dicts. 

#maybe call replace_ptms with name 
#build dicts 

def replace_ptms(stmts,name,dictionary):
    relevantStatements = ac.filter_gene_list(stmts,[name],'one')
    otherStatements = [st for st in stmts if st not in relevantStatements] 
    for st in relevantStatements:
        for ag in st.agent_list():
            #change context on any agents 
            if ag.name == name:
                newmodlist = []
                for mod in ag.mods:
                    if str(mod) in list(dictionary.keys()):
#                    if any([checkmod for checkmod in list(dictionary.keys()) if checkmod.matches(mod)]):
                        newmod = deepcopy(mod)
                        newmod.residue = dictionary[str(mod)]
                        newmod.position = None
                        if not any([mod for mod in newmodlist if newmod.matches(mod)]):
                            newmodlist.append(newmod)
                    else:
                        if not any([mod1 for mod1 in newmodlist if mod.matches(mod1)]):
                            newmodlist.append(mod)
                ag.mods = newmodlist 
        #change phosphorylated residues on a substrate 
        if isinstance(st,Modification):
            for entry in list(dictionary.keys()):
#                if st.sub.name == name and st.residue == entry.residue and st.position == entry.position:   #Don't love all the string matching here
                if st.sub.name == name and st.residue == entry.split(',')[1].strip() and st.position == entry.split(',')[2].strip().strip(')'):   #Don't love all the string matching here
                    st.residue = dictionary[str(entry)]
                    st.position = None
    outputStmts = relevantStatements + otherStatements
    return outputStmts 

def remove_dup_phos(stmts):
    outputStmts = deepcopy(stmts)
    for st in outputStmts:
        for ag in st.agent_list():
            new_mod_list = []
            for mod in ag.mods:
                if not any([newmod for newmod in new_mod_list if newmod.matches(mod)]):
                    new_mod_list.append(mod)
            ag.mods = new_mod_list
    return outputStmts




def coarse_grain_phos(stmts):
    stmts = remove_dup_phos(stmts)
    allNames = collect_agent_name_list(stmts)
    for name in allNames:
        dictList = build_ptm_dict(stmts,name)
        for dictionary in dictList:
            outputStmts = replace_ptms(stmts,name,dictionary)
    outputStmts = Preassembler.combine_duplicate_stmts(outputStmts)
    return outputStmts

#Going to need to rerun combining statements
#Add a deepcopy somewhere to not overwrite original statements 




#########################################################
#OLD


##rewrites active form statements to coarse grain pY and pS/T sites 
#def coarse_grain_phos_on_af(stmts):
#    af_stmts = ac.filter_by_type(stmts, ActiveForm)
#    non_af_stmts = ac.filter_by_type(stmts, ActiveForm, invert=True)
#    new_af_stmts = []
#    for st in af_stmts:   
#        act_agent = st.agent
#        #Preserve inhibitory statements
#        if st.is_active == True: 
#            #skip mutants
#            #Make this a separate function and call it?
#            if st.agent.mutations:
#                new_af_stmts.append(st)
#            else:
#                #in each AF statement, identify Tyrosine and Serine/Threonine phosphorylations
#                phos_mods = [mod for mod in st.agent.mods if mod.mod_type == 'phosphorylation']
#                pY_mods = [mod for mod in phos_mods if mod.residue == 'Y']
#                pST_mods = [mod for mod in phos_mods if mod.residue in ('S','T')] 
#                modlist = []
#                try:
#                    generic_pY = ModCondition(pY_mods[0].mod_type,pY_mods[0].residue)
#                    modlist.append(generic_pY)
#                except IndexError:
#                    generic_pY = []
#                try:
#                    generic_pST = ModCondition(pST_mods[0].mod_type,'S') #This could appear as S or T but I'm assuming equivaqlency. scratch that, make s, to remove redundancies
#                    modlist.append(generic_pST)
#                except IndexError:
#                    generic_pST = []
#                #Write replacement ActiveForm statement with coarse-grained modifications
#                newagent = Agent(act_agent.name,modlist)
#                newst = ActiveForm(newagent,st.activity,st.is_active)
#                new_af_stmts.append(newst)
#        else:
#            new_af_stmts.append(st)
#    new_af_stmts = ac.Preassembler.combine_duplicate_stmts(new_af_stmts)
#    filtered_stmts = non_af_stmts + new_af_stmts
#    return(filtered_stmts)



###This may be unnecessary if I'm now further coarse graining AF. Could be useful as a middle road
##        #loop through list of mods on activeform statement
##        for mod in act_agent.mods:
##            #loop through all phosphorylation statements
##            for pst in phos_stmts:
##                #check if mod in af statement matches a phosphorylation statement
##                if act_agent.entity_matches(pst.agent_list()[1]): #Second agent is one being modified
##                    if mod.position == pst.position:
##                        used_stmts.append(pst)

#        #Generic, coarse-grained
##        try:
##            newmod = ModCondition(mod_type=act_agent.mods[0].mod_type,residue=act_agent.mods[0].residue) #needs improvement
##        except IndexError:
##            #this is mut instead of mod, maybe other exceptions I haven't come across yet
##            newmod = None
##            #pass
##        try:
##            newmut = ModCondition(mod_type=act_agent.mods[0].mod_type,residue=act_agent.mods[0].residue) #needs improvement
##        except IndexError:
##            #this is mut instead of mod, maybe other exceptions I haven't come across yet
##            newmod = None#[]
##need to handle 'activating' mutations differently
##ActiveForm(HRAS(muts: (G, 12, V)), gtpbound, True) #Mutation condition may break things
##In [82]: af_stmts[1].agent.mutations
##Out[82]: [MutCondition(G, 12, V)]



##Sort through phosphorylation stmts and remove ones with the same enzyme and substrate   
#def remove_redundant_phosphorylations(stmts):

#    
#    phos_stmts = ac.filter_by_type(stmts, Phosphorylation) + ac.filter_by_type(stmts, Dephosphorylation)
#    non_phos_stmts = ac.filter_by_type(stmts, Phosphorylation, invert = True)
#    non_phos_stmts = ac.filter_by_type(non_phos_stmts, Dephosphorylation, invert = True)
#    #Start new list with coarse grained first statement
#    unique_phos_stmts = [phos_stmts[0]]
#    unique_phos_stmts[0].position=None
#    for st1 in phos_stmts:
#        unique = 1
#        for st2 in unique_phos_stmts:
#            #Comparae enz/sub pair for list of phos stmts with stmts in new list. If enz/sub match already added statement, discard. Otherwise add and coarse-grain (remove phos position)
#            if st1.enz.entity_matches(st2.enz) and st1.sub.entity_matches(st2.sub):    
#                unique = 0
##            else:
#        if unique == 1:
#            unique_phos_stmts.append(st1)
#            unique_phos_stmts[-1].position = None
#    filtered_stmts = non_phos_stmts + unique_phos_stmts
#    return filtered_stmts

##Coarse graining kinase context
###rewrites active form statements to coarse grain pY and pS/T sites 
#def coarse_grain_kinase_context(stmts):
#    phos_stmts = ac.filter_by_type(stmts, Phosphorylation) + ac.filter_by_type(stmts, Dephosphorylation)
#    non_phos_stmts = ac.filter_by_type(stmts, Phosphorylation, invert = True)
#    non_phos_stmts = ac.filter_by_type(non_phos_stmts, Dephosphorylation, invert = True)
#    kinase_with_mods = list(filter(lambda x: x.enz.mods != [], phos_stmts)) #had to add list() for pyth3 combatability
#    kinase_wo_mods = list(filter(lambda x: x.enz.mods == [], phos_stmts))
#    
#    coarse_grained_phos_stmts = []
#    for st in kinase_with_mods:
#        pY_mods = [mod for mod in st.enz.mods if mod.residue == 'Y']
#        pST_mods = [mod for mod in st.enz.mods if mod.residue in ('S','T')]  
#        modlist = []
#        try:
#            generic_pY = ModCondition(pY_mods[0].mod_type,pY_mods[0].residue)
#            modlist.append(generic_pY)
#        except IndexError:
#            generic_pY = []
#        try:
#            generic_pST = ModCondition(pST_mods[0].mod_type,'S') #This could appear as S or T but I'm assuming equivaqlency. Again, scratch this for S only
#            modlist.append(generic_pST)
#        except IndexError:
#            generic_pST = []
#        #Write replacement ActiveForm statement with coarse-grained modifications
#        newagent = Agent(st.enz.name,modlist)
#        newst = Phosphorylation(newagent,st.sub,st.residue)
#        coarse_grained_phos_stmts.append(newst)

#    new_phos_stmts = coarse_grained_phos_stmts+kinase_wo_mods
#    new_phos_stmts = ac.Preassembler.combine_duplicate_stmts(new_phos_stmts)
#    filtered_stmts = non_phos_stmts + new_phos_stmts
#    return filtered_stmts

#def coarse_grain_substrate_context(stmts):
#    phos_stmts = ac.filter_by_type(stmts, Phosphorylation) + ac.filter_by_type(stmts, Dephosphorylation)
#    non_phos_stmts = ac.filter_by_type(stmts, Phosphorylation, invert = True)
#    non_phos_stmts = ac.filter_by_type(non_phos_stmts, Dephosphorylation, invert = True)
#    substrate_with_mods = list(filter(lambda x: x.sub.mods != [], phos_stmts)) #had to add list() for pyth3 combatability
#    substrate_wo_mods = list(filter(lambda x: x.sub.mods == [], phos_stmts))
#    
#    coarse_grained_phos_stmts = []
#    for st in substrate_with_mods:
#        pY_mods = [mod for mod in st.sub.mods if mod.residue == 'Y']
#        pST_mods = [mod for mod in st.sub.mods if mod.residue in ('S','T')]  
#        modlist = []
#        try:
#            generic_pY = ModCondition(pY_mods[0].mod_type,pY_mods[0].residue)
#            modlist.append(generic_pY)
#        except IndexError:
#            generic_pY = []
#        try:
#            generic_pST = ModCondition(pST_mods[0].mod_type,'S') #This could appear as S or T but I'm assuming equivaqlency. Again, scratch this for S only
#            modlist.append(generic_pST)
#        except IndexError:
#            generic_pST = []
#        #Write replacement ActiveForm statement with coarse-grained modifications
#        newagent = Agent(st.sub.name,modlist)
#        newst = Phosphorylation(st.enz,newagent,st.residue)
#        coarse_grained_phos_stmts.append(newst)

#    new_phos_stmts = coarse_grained_phos_stmts+substrate_wo_mods
#    new_phos_stmts = ac.Preassembler.combine_duplicate_stmts(new_phos_stmts)
#    filtered_stmts = non_phos_stmts + new_phos_stmts
#    return filtered_stmts

