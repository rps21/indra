from indra.literature import *
from indra.databases import hgnc_client
import random

# Get a list of PMIDs for the gene
#with open('heartmachine_new_searchterms.txt','r') as f:
with open('drugtargets.txt','r') as f:
   targetgenes = f.readlines()
with open('experimental_observables.txt','r') as f:
   experimentgenes = f.readlines()
with open('other.txt','r') as f:
   othergenes = f.readlines()
pre_genelist = targetgenes+experimentgenes+othergenes

with open('drugnames.txt','r') as f:
   drugnames = f.readlines()

genelist = []
problem = []
for gene in pre_genelist:
    hgnc_id = hgnc_client.get_hgnc_id(gene.strip())
    if hgnc_id is None:
        problem.append(gene.strip())
    else:
        genelist.append(gene.strip())

with open('combined_gene_list.txt', 'w') as f:
   for gene in genelist:
    f.write('%s\n' % gene)

#with open('final_lists/prelim_not_gene_list.txt', 'w') as f:
#   for notgene in not_gene_list:
#    f.write('%s\n' % notgene)


with open('combined_gene_list.txt') as f:
    gene_search_terms = f.readlines()

#with open('final_not_gene_list.txt') as f:
#   non_gene_search_terms = f.readlines()   

pmids = []
for term in gene_search_terms:
    print(term)
    pmids = pmids + pubmed_client.get_ids_for_gene(term.strip(),retmax=2500)

for term in drugnames:
    print(term)
    pmids = pmids + pubmed_client.get_ids(term.strip(),retmax=2500)

# Get the PMIDs that have XML in PMC
#pmids_oa_xml = pmc_client.filter_pmids(pmids, 'oa_xml')

pmids = sorted(list(set(pmids)))
random.shuffle(pmids)

# Write the results to a file
with open('phase3eval_pt2_pmids.txt', 'w') as f:
   for pmid in pmids:
    f.write('%s\n' % pmid)
