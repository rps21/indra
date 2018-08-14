from indra.statements import *
from copy import deepcopy
from indra.preassembler import Preassembler
from indra.tools import assemble_corpus as ac

#TODO:


def addPhosDephos(phosStmts,dephosStmts):
    phosSubstrates = list(map(lambda stmt: stmt.sub.name, phosStmts))
    dephosSubstrates = list(map(lambda stmt: stmt.sub.name, dephosStmts))
    implicitAgent = Agent('GenericAgent')
    outputStmts = []

    for st in phosStmts:
        dephos_ag = deepcopy(st.sub)
        dephos_res = deepcopy(st.residue)
        dephos_pos = deepcopy(st.position)

        if dephos_ag.name not in dephosSubstrates: #This means there is no statement dephosphorylating an agent that is phosphorylated
            dephos_mod = ModCondition(mod_type='phosphorylation',residue=dephos_res,position=dephos_pos)
            dephos_ag.mods = [dephos_mod]
            dephos_st = Dephosphorylation(implicitAgent, dephos_ag, residue=dephos_res, position=dephos_pos)
            outputStmts.append(dephos_st)

    for st in dephosStmts:
        phos_ag = deepcopy(st.sub)
        phos_res = deepcopy(st.residue)
        phos_pos = deepcopy(st.position)

        if phos_ag.name not in phosSubstrates: #This means there is no statement dephosphorylating an agent that is phosphorylated
            phos_mod = ModCondition(mod_type='phosphorylation',residue=phos_res,position=phos_pos)
            phos_ag.mods = [phos_mod]
            phos_st = Phosphorylation(implicitAgent, phos_ag, residue=phos_res, position=phos_pos)
            outputStmts.append(phos_st)

    outputStmts = Preassembler.combine_duplicate_stmts(outputStmts)
    return outputStmts



def addIncreaseDecrease(increaseStmts,decreaseStmts):
    increaseObjects = list(map(lambda stmt: stmt.obj.name, increaseStmts))
    decreaseObjects = list(map(lambda stmt: stmt.obj.name, decreaseStmts))
    implicitAgent = Agent('GenericAgent')
    outputStmts = []

    for st in increaseStmts:
        deg_ag = deepcopy(st.obj)
        if deg_ag.name not in decreaseObjects:
            deg_st = DecreaseAmount(implicitAgent, deg_ag)
            outputStmts.append(deg_st)

    for st in decreaseStmts:
        inc_ag = deepcopy(st.obj)
        if inc_ag.name not in increaseObjects:
            inc_st = IncreaseAmount(implicitAgent, inc_ag)
            outputStmts.append(inc_st)

    outputStmts = Preassembler.combine_duplicate_stmts(outputStmts)
    return outputStmts



def addActivateDeactivate(activateStmts,deactivateStmts):
    activateObjects = list(map(lambda stmt: stmt.obj.name, activateStmts))
    deactivateObjects = list(map(lambda stmt: stmt.obj.name, deactivateStmts))
    implicitAgent = Agent('GenericAgent')
    outputStmts = []

    for st in activateStmts:
        deact_ag = deepcopy(st.obj)
        if deact_ag.name not in deactivateObjects:
            deact_st = DecreaseAmount(implicitAgent, deact_ag)
            outputStmts.append(deact_st)

    for st in deactivateStmts:
        act_ag = deepcopy(st.obj)
        if act_ag.name not in activateObjects:
            act_st = IncreaseAmount(implicitAgent, act_ag)
            output_stmts.append(act_st)

    outputStmts = Preassembler.combine_duplicate_stmts(outputStmts)
    return outputStmts



def addAll(stmts):
    newStmts = []
    phosStmts = ac.filter_by_type(stmts,Phosphorylation) 
    dephosStmts = ac.filter_by_type(stmts,Dephosphorylation)
    newStmts = addPhosDephos(phosStmts,dephosStmts) + newStmts

    increaseStmts = ac.filter_by_type(stmts,IncreaseAmount)
    decreaseStmts = ac.filter_by_type(stmts,DecreaseAmount)
    newStmts = addIncreaseDecrease(increaseStmts,decreaseStmts) + newStmts

    activateStmts = ac.filter_by_type(stmts,Activation)
    deactivateStmts = ac.filter_by_type(stmts,Inhibition)
    newStmts = addActivateDeactivate(activateStmts,deactivateStmts) + newStmts

    outputStmts = stmts + newStmts
    
    return outputStmts








