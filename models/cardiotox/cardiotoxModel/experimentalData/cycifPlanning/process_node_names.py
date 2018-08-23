import pickle

with open('node_names.csv','r') as f:
   raw_nodes = f.readlines()

raw_nodes = raw_nodes[1:len(raw_nodes)-1]
processed_nodes = []
for node in raw_nodes:
    processed_nodes.append(node.split(',')[1].strip().upper())

with open('processed_node_names.pkl','wb') as f:
    pickle.dump(processed_nodes,f)

