from indra.sources import trips

targetSentences = 'Sorafenib dephosphorylates FLT3. Sorafenib inhibits FLT3. Sorafenib dephosphorylates KDR. Sorafenib inhibits KDR. Sorafenib dephosphorylates PDGFRA. Sorafenib inhibits PDGFRA.'
target_stmts = []
trips_processor = trips.process_text(targetSentences)
target_stmts = target_stmts + trips_processor.statements 
ac.dump_statements(target_stmts, 'sorafenibTargetStmts2.pkl')




