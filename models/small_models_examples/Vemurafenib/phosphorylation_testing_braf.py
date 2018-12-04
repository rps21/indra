import pickle
from indra.sources import trips
from indra.tools.small_model_tools.modelContext import enforceCascadeContext as cs
from indra.tools.small_model_tools.modelContext import combinePhosphorylationSites as ptm
from indra.tools.small_model_tools.modelContext import addImplicitMechs as aim
from indra.tools import assemble_corpus as ac 

prior_model_stmts = ac.load_statements('modelPrior.pkl')
brafStmts = ac.filter_gene_list(prior_model_stmts,['BRAF'],'all')
contextStmts = cs.add_all_af(brafStmts)
ptmStmts = ptm.coarse_grain_phos(contextStmts)
context2 = cs.combine_multiple_activeforms(ptmStmts)


