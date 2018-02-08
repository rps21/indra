from indra.tools import assemble_corpus as ac
from indra.tools.assemble_corpus import filter_by_type
from indra import statements
from indra.statements import *
from filter_unique_nodes import filter_unique_species

metgenes = ['HK','GPI','PFK','ALDO','GAPDH','PGK','PGM','ENO','PKM2','PKM','G6PC','FBP','PEPCK','PC','PGHDH','PSAT1','PSPH']
targetgenes = ['VEGFR','KDR','PDGFR','CSF1R','FLT2','KIT','BRAF','EGFR','ERBB2']
expgenes = ['MLC2V','MYL2','MYLK3','p53','TP53','STAT5','STAT5A','STAT5B','SRC','cJUN','JUN','MEK','MAP2K1','STAT1','cMyc','MYC','CAMKII','CAMK2N2','MLC2a','MYL7','STAT3','CTNNB1','HIF1A','14-3-3','YWHAZ','LC3','MAP1LC3A','MAP1LC3A','mTOR','ACTA1','HCN4','FLT','FLT3','MCM6','p44','MAPK3','p42','PSMC6','S6','RPS6KB1','RPS6KA1','RPS6','PCNA','EGFR','p21','CDKN1A','ACTB','COX4','COX4I1','p38','MAPK11','MAPK12','MAPK13','MAPK14','SAPK','JNK1','NFkB','NFKB1','NFKB2','RELA','REL','RELB','PKM2','PKM','cMYC','MYC','CASP3','TUBB','Aurora','AURKA','AURKB','PDGFR','PDGFRA','PDGFRB','PIK3R1']

#genes_to_keep = metgenes + targetgenes + expgenes

genes_to_keep = newnodes

##basic implementation. includes nodes that are in more than one statement, but may be dead-ends in that they have only one binding partner
##i.e. they are in 2 statements - one binding to partner, two phosphorylation by partner
#def filter_dead_end(stmts,nodes):
#    final_stmts = stmts[:]
#    for nd in nodes:
#        counter = 0
#        for st in stmts:
#            if nd in str(st):
#                temp_st = st
#                counter = counter + 1
#        if counter < 2:
##            print(temp_st)
#            if temp_st in final_stmts:
#                final_stmts.remove(temp_st)
#        temp_st = None
#    return final_stmts
erknodes = ['EGF','EGFR','SOS1','GRB2','KRAS','BRAF','MAPK1','MAPK3','MAP2K1']

def filter_dead_end(stmts,nodes):
    nodes_to_remove = []
    final_stmts = stmts[:]
    new_nodes = nodes[:]
#    binding_stmts = ac.filter_by_type(stmts,Complex,invert=False)
    for nd in nodes:
        counter = 0
        for st in stmts:
            if nd in str(st):
                temp_st = st
                counter = counter + 1
        if counter < 2:
            #This is where change needs to happen
            if str(nd).split('(')[0] not in erknodes:
                nodes_to_remove.append(nd)
                new_nodes.remove(nd)
                print(nodes_to_remove)
#        temp_st = None
#    for nd in nodes_to_remove:
#        for st in stmts:
#            if str(nd) in str(st):
#                if st in final_stmts:
#                    final_stmts.remove(st)
#    return final_stmts
    return new_nodes


