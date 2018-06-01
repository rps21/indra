import requests
from indra.statements import stmts_from_json

resp = requests.get('https://lsm6zea7gg.execute-api.us-east-1.amazonaws.com/production/statements/',
                    headers={'x-api-key': 'XH36SVAGBN9L3SA8XNuHu5hvJ9v3j9mq8PTkPYjG'},
                    params={'subject': 'MAP2K1',
                            'object': 'MAPK1',
                            'type': 'phosphorylation'})
stmts_json = resp.json()
stmts = stmts_from_json(stmts_json)


#params=u'agent=SMAD2&agent=SMURF2')
resp = requests.get('https://lsm6zea7gg.execute-api.us-east-1.amazonaws.com/production/statements/',
                    headers={'x-api-key': 'XH36SVAGBN9L3SA8XNuHu5hvJ9v3j9mq8PTkPYjG'},
                    params=u'agent=MAP2K1&agent=MAPK1')
stmts_json = resp.json()
stmts = stmts_from_json(stmts_json)





