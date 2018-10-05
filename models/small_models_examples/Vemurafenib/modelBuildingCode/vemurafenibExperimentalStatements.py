
from indra.sources import trips

soraf_targets = ['']

sorafStmts = []

for target in soraf_targets:

    targetStmts = ac.filter_gene_list(smallModelStmts_contextChanges,[target],'one')    #Need model stmts here
    targetMods = []
    for st in targetStmts:
        for ag in st.agent_list():
            if ag.name == target:
                for mod in ag.mods:
                    if not any([m for m in targetMods if m.matches(mod)]):
                        targetMods.append(mod)

    baseSentence = 'Sorafenib dephosphorylates %s' % target
    trips_processor = trips.process_text(baseSentence)
    baseStmt = trips_processor.statements[0]
    for mod in targetMods:
        newStmt = deepcopy(baseStmt)
        newStmt.residue = mod.residue
        newStmt.position = mod.position
        sorafStmts.append(newStmt)
ac.dump_statements(sorafStmts,'sorafenibTargetStmts.pkl')



experimentalSentences = 'Sorafenib phosphorylates cJun. Sorafenib phosphorylates STAT1. Sorafenib phosphorylates PKM2. Sorafenib dephosphorylates RPS6. Sorafenib phosphorylates Aurora kinase. Sorafenib transcribes HIF1A. Sorafenib transcribes cMyc. Sorafenib transcribes PDGFRA.'  

exp_stmts = []
trips_processor = trips.process_text(experimentalSentences)
exp_stmts = exp_stmts + trips_processor.statements 


#############


