In [3]: old_af
Out[3]: 
[ActiveForm(BRAF(mods: (phosphorylation, T, 753), (phosphorylation, S, 750), (phosphorylation, S, 151), (phosphorylation, T, 401)), kinase, False),
 ActiveForm(BRAF(mods: (phosphorylation, S, 446), (phosphorylation, S, 729)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, S, 579)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 373)), kinase, True),
 ActiveForm(BRAF(muts: (V, 600, E)), kinase, True),

 ActiveForm(EGFR(mods: (phosphorylation, Y)), kinase, True),

 ActiveForm(KRAS(muts: (G, 12, C), (G, 12, V), (G, 13, D), (Q, 61, H), (Q, 61, L)), gtpbound, True),
 ActiveForm(KRAS(muts: (G, 12, V), (G, 12, C), (G, 13, D), (Q, 61, H), (Q, 61, L)), gtpbound, True),
 ActiveForm(KRAS(muts: (G, 13, D), (G, 12, C), (G, 12, V), (Q, 61, H), (Q, 61, L)), gtpbound, True),
 ActiveForm(KRAS(muts: (Q, 61, H), (G, 12, C), (G, 12, V), (G, 13, D), (Q, 61, L)), gtpbound, True),
 ActiveForm(KRAS(muts: (Q, 61, L), (G, 12, C), (G, 12, V), (G, 13, D), (Q, 61, H)), gtpbound, True),

 ActiveForm(MAP2K1(mods: (phosphorylation, S, 218), (phosphorylation, S, 222), (phosphorylation, T, 286), (phosphorylation, T, 292)), kinase, True),
 ActiveForm(MAPK1(mods: (phosphorylation, Y, 187), (phosphorylation, T, 185)), kinase, True),
 ActiveForm(MAPK3(mods: (phosphorylation, T, 202), (phosphorylation, Y, 204)), kinase, True)]


#without last step
[ActiveForm(BRAF(mods: (phosphorylation, T, 753), (phosphorylation, S, 750), (phosphorylation, S, 151), (phosphorylation, T, 401)), kinase, False),
 ActiveForm(BRAF(mods: (phosphorylation, S, 446), (phosphorylation, S, 729)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, S, 579)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 373)), kinase, True),
 ActiveForm(BRAF(muts: (V, 600, E)), kinase, True),

 ActiveForm(EGFR(mods: (phosphorylation, Y)), kinase, True),
 ActiveForm(EGFR(mods: (phosphorylation, Y), bound: [EGF, True]), kinase, True),

 ActiveForm(GRB2(bound: [EGFR, True]), kinase, True),

 ActiveForm(KRAS(muts: (Q, 61, H), (G, 12, C), (G, 12, V), (G, 13, D), (Q, 61, L)), gtpbound, True),

 ActiveForm(MAP2K1(mods: (phosphorylation, S, 218), (phosphorylation, S, 222), (phosphorylation, T, 286), (phosphorylation, T, 292)), kinase, True),
 ActiveForm(MAPK1(mods: (phosphorylation, Y, 187), (phosphorylation, T, 185)), kinase, True),
 ActiveForm(MAPK3(mods: (phosphorylation, T, 202), (phosphorylation, Y, 204)), kinase, True),

 ActiveForm(SOS1(bound: [GRB2, True]), kinase, True)]

#Added:
 ActiveForm(EGFR(mods: (phosphorylation, Y), bound: [EGF, True]), kinase, True),
 ActiveForm(GRB2(bound: [EGFR, True]), kinase, True),
 ActiveForm(SOS1(bound: [GRB2, True]), kinase, True)]



#Problems: 
#Loss of KRAS mutant statments.  - This comes from Preassembler.combine_duplicate_stmts, combines all KRAS mutant statements for some reason. Not my bug.
#2nd EGFR statment having both phos and bound, only want bound 
# ActiveForm(EGFR(mods: (phosphorylation, Y), bound: [EGF, True]), kinase, True),


