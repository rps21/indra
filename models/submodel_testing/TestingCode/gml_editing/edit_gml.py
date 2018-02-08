import re
from copy import deepcopy
import networkx as nx

with open('secondpass_v2_contactmap.gml','r') as f:
    fullGraph = f.readlines()

nodes = []
edges = []

for line in fullGraph:
    if 'edge' in line:
        edges.append(line)
    else:
        nodes.append(line)

#Fix edges
edges_copy = []
for edge in edges:
    oldEdgeIndexes = re.findall('\d+',edge)
    oldEdgeIndex1 = oldEdgeIndexes[0]
    oldEdgeIndex2 = oldEdgeIndexes[1]

    if oldEdgeIndexes:
        for num, line in enumerate(nodes, 1):
            if oldEdgeIndex1 in line:
                linenum = num
                for i in range(num,0,-1):
                    if len(nodes[i]) > 7: #This is so sloppy
                        if 'gid' not in nodes[i]:
                            newEdgeIndex1 = re.search('\d+',nodes[i]).group(0)
                            break
            elif oldEdgeIndex2 in line:
                linenum = num
                for i in range(num,0,-1):
                    if len(nodes[i]) > 7: #This is so sloppy
                        if 'gid' not in nodes[i]:
                            newEdgeIndex2 = re.search('\d+',nodes[i]).group(0)
                            break

    edge = edge.replace(oldEdgeIndex1,newEdgeIndex1)
    edge = edge.replace(oldEdgeIndex2,newEdgeIndex2)
    edges_copy.append(edge)


#Remove subnodes

nodes_to_remove = []
nodes_copy = deepcopy(nodes)
for num, node in enumerate(nodes,1):
    if len(node) > 7:
        if 'gid' in node:
            nodes_to_remove = nodes_to_remove + [num-2,num-1,num]

for node in reversed(nodes_to_remove):
    del nodes_copy[node]

#Remove 'directed' and remove new line before open bracket in first line. 
del nodes_copy[2]
nodes_copy[0] = nodes[0].replace('\n',' ')

edges_copy = list(set(edges_copy)) #remove duplicate edges

del nodes_copy[-1] #don't close graph without edges. Should be smoother way to do this upstream
newGraph = nodes_copy + edges_copy
newGraph.append(']') #Close graph definition



with open('newGraph.gml','w') as f:
    for line in newGraph:
        f.write("%s\n" % line)

##Find cycles
#g = nx.read_gml('newGraph.gml')
#nx.cycle_basis(g)

#[['ImplicitPhos', 'MAPK3', 'MAPK1'],
# ['BRAF', 'MAPK3', 'MAPK1'],
# ['MAP2K1', 'MAPK3', 'MAPK1'],
# ['ImplicitPhos', 'MAP2K1', 'MAPK1'],
# ['BRAF', 'MAP2K1', 'MAPK1'],
# ['ImplicitPhos', 'BRAF', 'MAPK1']]


