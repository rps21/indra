

drugSentences = 'SORAFENIB dephosphorylates PDGFRA. SORAFENIB dephosphorylates KDR. SORAFENIB dephosphorylates FLT3'
drug_stmts = buildDrugTargetStmts(drugSentences)

#testingStmts=drug_stmts
#contextStmts = cs.add_all_af(testingStmts)
#ptmStmts = ptm.coarse_grain_phos(contextStmts)
#context2 = cs.combine_multiple_phos_activeforms(ptmStmts)



drugSentences = 'BRAF V600E bound to VEMURAFENIB is inactive. VEMURAFENIB BINDS BRAF V600E' 
drug_stmts = buildDrugTargetStmts(drugSentences)


 
[Dephosphorylation(SORAFENIB(), PDGFRA()),
 Dephosphorylation(SORAFENIB(), KDR()),
 Dephosphorylation(SORAFENIB(), FLT3())]

In [8]: contextStmts = cs.add_all_af(testingStmts)
   ...: ptmStmts = ptm.coarse_grain_phos(contextStmts)
   ...: 

In [9]: ptmStmts
Out[9]: 
[ActiveForm(FLT3(mods: (phosphorylation, phos_act)), activity, True),
 ActiveForm(KDR(mods: (phosphorylation, phos_act)), activity, True),
 ActiveForm(PDGFRA(mods: (phosphorylation, phos_act)), activity, True),
 Complex(FLT3LG(), FLT3()),
 Complex(VEGF(), KDR()),
 Complex(PDGFA(), PDGFRA()),
 Dephosphorylation(SORAFENIB(), FLT3(), phos_act),
 Dephosphorylation(SORAFENIB(), KDR(), phos_act),
 Dephosphorylation(SORAFENIB(), PDGFRA(), phos_act),
 Phosphorylation(FLT3(bound: [FLT3LG, True]), FLT3(), phos_act),
 Phosphorylation(KDR(bound: [VEGF, True]), KDR(), phos_act),
 Phosphorylation(PDGFRA(bound: [PDGFA, True]), PDGFRA(), phos_act)]




['MAPK3', 'BRAF', 'JUN', 'KRAS', 'EGFR', 'ERBB2', 'MAP2K1', 'EGFR']
expSentences = 'VEMURAFENIB degrades EGFR.'# VEMURAFENIB transcribes ERBB2'

#inhibit a tf
Jun transcribes EGFR
Jun activated by MAPK3
MAPK3 activated by MAP2K1
BRAFV600E activates map2k1


mutations = {'BRAF': [('Val', '600', 'Glu'),('V', '600', 'E')]}
#drugSentences = 'VEMURAFENIB dephosphorylates BRAF V600E.' 

drugSentences = 'BRAF V600E bound to VEMURAFENIB is inactive. VEMURAFENIB BINDS BRAF V600E' 
