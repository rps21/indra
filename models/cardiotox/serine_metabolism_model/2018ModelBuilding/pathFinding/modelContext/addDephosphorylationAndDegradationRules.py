from indra.statements import *
from copy import deepcopy
from indra.preassembler import Preassembler
from indra.tools import assemble_corpus as ac

#Enforce dephos and degredation
#Check unbinding context 

#TODO:
#Break up multiple phosp sites
#Ex:  Dephosphorylation(None, BRAF(mods: (phosphorylation, S, 446), (phosphorylation, T, 599), (phosphorylation, S, 602))),



def add_dephosphorylations(stmts):
    implicit_phosphotase = Agent('GenericPhosphatase')
    output_stmts = deepcopy(stmts)
    new_dephos = []
    phos_stmts = ac.filter_by_type(stmts,Phosphorylation)

    for st in phos_stmts:
        dephos_ag = deepcopy(st.sub)
        dephos_res = deepcopy(st.residue)
        dephos_pos = deepcopy(st.position)

        #check for existing dephos 
        dephosStmts = ac.filter_by_type(stmts,Dephosphorylation)
        originalDephos = ac.filter_gene_list(dephosStmts,[dephos_ag.name],'one')
        correctDephos = []
        for st in originalDephos:
            if st.enz == dephos_ag.name:
                correctDephos.append(st)
        if not correctDephos:        

            dephos_mod = ModCondition(mod_type='phosphorylation',residue=dephos_res,position=dephos_pos)
            dephos_ag.mods = [dephos_mod]

            dephos_st = Dephosphorylation(implicit_phosphotase, dephos_ag, residue=dephos_res, position=dephos_pos)
            new_dephos.append(dephos_st)
            print(dephos_st)
    new_dephos = Preassembler.combine_duplicate_stmts(new_dephos)
    output_stmts = output_stmts + new_dephos
    return output_stmts

def add_degradations(stmts):
    output_stmts = deepcopy(stmts)
    new_deg = []
    transcription_stmts = ac.filter_by_type(stmts,IncreaseAmount)

    for st in transcription_stmts:
        deg_ag = deepcopy(st.obj)
        deg_st = DecreaseAmount(None, deg_ag)
        new_deg.append(deg_st)
    new_deg = Preassembler.combine_duplicate_stmts(new_deg)
    output_stmts = output_stmts + new_deg
    return output_stmts





















