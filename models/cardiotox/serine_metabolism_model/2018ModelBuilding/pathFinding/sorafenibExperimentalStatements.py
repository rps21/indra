#Sorafenib increases cJun phosphorylation
#Sorafenib increases STAT1 phosphorylation
#Sorafenib increases PKM2 phosphorylation
#Sorafenib decreases phosphorylation of RPS6
#Sorafenib increases Aurora kinase phosphorylation - impossible in model 

#Sorafenib increases the amount of HIF1A
#Sorafenib increases cMyc abundance.
#Sorafenib increases PDGFRA abundance - impossible in model 


from indra.sources import trips

soraf_targets = ['KDR','VEGFR1','FLT3','PDGFRA']

sorafSentences = 'Sorafenib inhibits KDR. Sorafenib inhibits VEGFR1. Sorafenib inhibits FLT3. Sorafenib inhibits PDGFRA.'


soraf_stmts = []
trips_processor = trips.process_text(sorafSentences)
soraf_stmts = soraf_stmts + trips_processor.statements 


experimentalSentences = 'Sorafenib phosphorylates cJun. Sorafenib phosphorylates STAT1. Sorafenib phosphorylates PKM2. Sorafenib dephosphorylates RPS6. Sorafenib phosphorylates Aurora kinase. Sorafenib transcribes HIF1A. Sorafenib transcribes cMyc. Sorafenib transcribes PDGFRA.'  

exp_stmts = []
trips_processor = trips.process_text(experimentalSentences)
exp_stmts = exp_stmts + trips_processor.statements 


#############


