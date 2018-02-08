#preparing for read.
#st.evidence is a list.
In [83]: stmts_in[-1].evidence[0].source_api
Out[83]: 'reach'
In [84]: stmts_in[-1].evidence[0].text

In [106]: stmts_in[-1].evidence[0].source_api=='reach'
Out[106]: True


In [114]: sentences = []
     ...: for st in testst:
     ...:     for ev in st.evidence:
     ...:         if ev.source_api == 'reach':
     ...:             #print(ev.source_api)
     ...:             sentences.append(ev.text)



from indra.sources import trips
sentence = 'MAP2K1 phosphorylates MAPK3 at Thr-202 and Tyr-204'
trips_processor = trips.process_text(sentence)
 trips_processor.statements
#gives sentence, which we can accumulate and save
#probably easy to build list of sentences for all stmts, feed to trips
#will have to be cautious of computation time, probably only do this after filtering initial stmts.
#BUT before and coarsegraining or other editing of stmt context 
