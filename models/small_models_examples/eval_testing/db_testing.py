import requests
from indra.statements import stmts_from_json
from indra.tools import assemble_corpus as ac 
from indra.mechlinker import MechLinker
from indra.statements import *
  
resp = requests.get('https://lsm6zea7gg.execute-api.us-east-1.amazonaws.com/production/statements/',
                headers={'x-api-key': 'XH36SVAGBN9L3SA8XNuHu5hvJ9v3j9mq8PTkPYjG'},
                params={'agent': 'KRAS','type':'activeform'})
stmts_json = resp.json()
stmts = stmts_from_json(stmts_json['statements'])
stmts = cleanStatements(stmts)


resp = requests.get('https://lsm6zea7gg.execute-api.us-east-1.amazonaws.com/production/statements/',
                headers={'x-api-key': 'XH36SVAGBN9L3SA8XNuHu5hvJ9v3j9mq8PTkPYjG'},
                params={'subject':'PTEN','object': 'MAPK3',
                        'type': 'dephosphorylation'})
stmts_json = resp.json()
stmts = stmts_from_json(stmts_json['statements'])


resp = requests.get('https://lsm6zea7gg.execute-api.us-east-1.amazonaws.com/production/statements/',
                headers={'x-api-key': 'XH36SVAGBN9L3SA8XNuHu5hvJ9v3j9mq8PTkPYjG'},
                params={'object': 'PIP-3'})
stmts_json = resp.json()
stmts = stmts_from_json(stmts_json['statements'])
stmts = cleanStatements(stmts)



from indra.sources import trips
sentences = 'PI3K phosphorylates PIP2.'
tp = trips.process_text(sentences)
stmts = tp.statements

sentences = 'PTEN dephosphorylates PIP3.'
tp = trips.process_text(sentences)
stmts = tp.statements


db_rest_api readme
  `subject=6871@HGNC`.

CHEBI:16618 (pip3)
agent=ERK@TEXT


resp = requests.get('https://lsm6zea7gg.execute-api.us-east-1.amazonaws.com/production/statements/',
                headers={'x-api-key': 'XH36SVAGBN9L3SA8XNuHu5hvJ9v3j9mq8PTkPYjG'},
                params={'agent': '16618@CHEBI'})
stmts_json = resp.json()
stmts = stmts_from_json(stmts_json['statements'])
stmts = cleanStatements(stmts)


resp = requests.get('https://lsm6zea7gg.execute-api.us-east-1.amazonaws.com/production/statements/',
                headers={'x-api-key': 'XH36SVAGBN9L3SA8XNuHu5hvJ9v3j9mq8PTkPYjG'},
                params={'agent': 'PIP-3@TEXT'})
stmts_json = resp.json()
stmts = stmts_from_json(stmts_json['statements'])
stmts = cleanStatements(stmts)
