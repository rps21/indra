


#########################

node [
id 14 label "GRB2" isGroup 1  graphics [ type "roundrectangle" fill "#D2D2D2" outline "#000000"  ] LabelGraphics [ text "GRB2" anchor "t" fontStyle "bold"  ] 
 ]
node [
id 15 label "catalytic" isGroup 1 gid 14  graphics [ type "roundrectangle" fill "#FFFFFF" outline "#000000"  ] LabelGraphics [ text "catalytic" anchor "t"  ] 
 ]
node [
id 16 label "inactive" gid 15  graphics [ type "roundrectangle" fill "#FFCC00" outline "#000000"  ] LabelGraphics [ text "inactive" anchor "c"  ] 
 ]
node [
id 17 label "active" gid 15  graphics [ type "roundrectangle" fill "#FFCC00" outline "#000000"  ] LabelGraphics [ text "active" anchor "c"  ] 
 ]

edge [ source 8 target 8  graphics [ fill "#000000"  ] ]
edge [ source 19 target 9  graphics [ fill "#000000"  ] ]
edge [ source 10 target 28  graphics [ fill "#000000"  ] ]
edge [ source 1 target 7  graphics [ fill "#000000"  ] ]
edge [ source 18 target 2  graphics [ fill "#000000"  ] ]
edge [ source 23 target 29  graphics [ fill "#000000"  ] ]


#don't care about graphics
#isGroup means there will be subnodes. gid not yet specified
#id 4 label "Y" isGroup 1 gid 3
#sub-nodes with further subnodes have isGroup, but also gid.
#anything with gid is subnode - may be able to simply remove all of these

#readline()
#if line contains node
#check next line for gid
#if yes, remove line, next line, 3rd line
#don't remove, just don't write to new list for new file

from bisect import bisect_left

with open('egf_egfr_sos_grb_olig_contactmap.gml','r') as f:
    original_graph = f.readlines()

#edges
major_nodes = []

i=0
new_graph = ['graph [']
for line in original_graph:
    if 'node' in line:
        if 'gid' not in original_graph[i+1]:
            major_nodes.append(int(original_graph[i+1].split(' ')[1]))
#            new_graph.append(original_graph[i:i+2])
            new_graph = new_graph+list(map(str.strip,original_graph[i:i+3]))
#    elif 'edge' in line:
#        new_graph.append(line)
    i=i+1

for line in original_graph:
    if 'edge' in line:
        src_subnode = line.split(' ')[3]
        tar_subnode = line.split(' ')[5]
        src_node = major_nodes[bisect_left(major_nodes,int(src_subnode))-1]
        tar_node = major_nodes[bisect_left(major_nodes,int(tar_subnode))-1]
    
        edge_line = 'edge [ source %s target %s  graphics [ fill "#000000"  ] ]' % (src_node,tar_node)
        new_graph.append(edge_line)


#edge [ source 23 target 29  graphics [ fill "#000000"  ] ]


new_graph.append(']')

with open('test_reduced_graph.gml','w') as newfile:
    for line in new_graph:
      newfile.write("%s\n" % line)

import networkx as nx
test_graph = nx.read_gml('test_reduced_graph.gml')
cyc = nx.cycle_basis(test_graph)

#Can't copy the edges directly
#Were between subnodes.
#Need to rewrite them

#go to edge first
#fine node id
#find first node upstream with no gid

