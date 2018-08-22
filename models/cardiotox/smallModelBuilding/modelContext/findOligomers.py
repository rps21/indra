import re
import shutil
from copy import deepcopy
import networkx as nx
import pysb
from pysb import bng
from indra.assemblers import PysbAssembler

def generateOriginalGML(stmts):
    console = bng.BngConsole()
    outputdir = console.base_directory

    pa = PysbAssembler('two_step')
    pa.add_statements(stmts)
    model = pa.make_model()

    bngl_file_path = console.base_directory+'/model.bngl'
    bngl_model = pysb.export.export(model,'bngl')
    bngl_file = open(bngl_file_path,'w')
    bngl_file.write(bngl_model)
    bngl_file.close()

    console.load_bngl(bngl_file_path)
    console.action('visualize',type='contactmap')

    gml_file_path = outputdir+'/model_contactmap.gml'
    shutil.move(gml_file_path, './model_contactmap.gml')
    final_gml_file_path = './model_contactmap.gml'
    return final_gml_file_path



def generateModifiedGML(gml_file_path):
    with open(gml_file_path,'r') as f:
        fullGraph = f.readlines()

    nodes = []
    edges = []

    for line in fullGraph:
        if 'edge' in line:
            edges.append(line)
        else:
            nodes.append(line)

    #Fix edges
    edgesNew = []
    edgeIndexList = []
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
        newEdgePair = set([newEdgeIndex1,newEdgeIndex2])
        if newEdgePair not in edgeIndexList:
            edgeNew = edge.replace(oldEdgeIndex1,newEdgeIndex1).replace(oldEdgeIndex2,newEdgeIndex2)
            edgesNew.append(edgeNew)
            edgeIndexList.append(newEdgePair)
#            print(newEdgePair)
#            print(edgePairList)

    #Remove subnodes
    nodes_to_remove = []
    nodesNew = deepcopy(nodes)
    for num, node in enumerate(nodes,1):
        if len(node) > 7:
            if 'gid' in node:
                nodes_to_remove = nodes_to_remove + [num-2,num-1,num]

    for node in reversed(nodes_to_remove):
        del nodesNew[node]

    #Remove 'directed' and remove new line before open bracket in first line. 
    del nodesNew[2]
    nodesNew[0] = nodes[0].replace('\n',' ')

    edgesNew = list(set(edgesNew)) #remove duplicate edges

    del nodesNew[-1] #don't close graph without edges. Should be smoother way to do this upstream
    newGraph = nodesNew + edgesNew
    newGraph.append(']') #Close graph definition


    newGMLFile = 'fixed_model_contactmap.gml'
    with open(newGMLFile,'w') as f:
        for line in newGraph:
            f.write("%s\n" % line)

    return newGMLFile


def processNetworkx(gml_file_path):
    g = nx.read_gml(gml_file_path)
    cycles = nx.cycle_basis(g)
    return cycles

def findCycles(stmts):
    origGMLPath = generateOriginalGML(stmts)
    newGMLPath = generateModifiedGML(origGMLPath)
    cycles = processNetworkx(newGMLPath)
    return cycles

#cycles = findCycles(braf_rawStmts)



