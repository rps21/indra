from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from copy import deepcopy
import pickle
from indra.tools import assemble_corpus as ac
from indra.preassembler import Preassembler
from indra.statements import *      

#TODO: ADD OTHER MODS



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



def build_ptm_dict_genericres(stmts,name):
    actMods, inhibMods, neutralMods = find_act_inhib_sites(stmts,name)
    actModsDict = {}
    inhibModsDict = {}
    neutralModsDict = {}
    for mod in actMods:
        actModsDict[str(mod)] = 'phos_act'
    for mod in inhibMods:
        inhibModsDict[str(mod)] = 'phos_inhib'
    for mod in neutralMods:
        neutralModsDict[str(mod)] = 'phos_act'
    return [actModsDict, inhibModsDict, neutralModsDict]

def build_ptm_dict_keepres(stmts,name):
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
        #Assume 'neutral' are activating for now. May want to revisit this. 
        if mod.residue == 'S' or mod.residue == 'T':
            neutralModsDict[str(mod)] = 'S_T_act'
        elif mod.residue == 'Y':
            neutralModsDict[str(mod)] = 'Y_act'
    return [actModsDict, inhibModsDict, neutralModsDict]




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
                try:
                    pos = entry.split(',')[2].strip().strip(')')
                except IndexError:
                    pos = None
                try: 
                    res = entry.split(',')[1].strip()
                except IndexError:
                    res = None

                if st.sub.name == name and st.residue == res and st.position == pos:   
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




def coarse_grain_phos(stmts,generic=True):
    stmts = remove_dup_phos(stmts)
    allNames = collect_agent_name_list(stmts)
    for name in allNames:
        if generic == True:
            dictList = build_ptm_dict_genericres(stmts,name)
        elif generic == False:
            dictList = build_ptm_dict_keepres(stmts,name)
        for dictionary in dictList:
            outputStmts = replace_ptms(stmts,name,dictionary)
    outputStmts = Preassembler.combine_duplicate_stmts(outputStmts)
    return outputStmts


