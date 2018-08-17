from pysb import *
from pysb.core import MonomerPattern
from copy import deepcopy

def addSimParamters(method='ode',equil=True,equilSpecies=[],viz=True):
    if method == "ode":
        newText = ("generate_network({overwrite=>1})\n"
                  "writeMexfile()\n")               
        if equil:
            for species in equilSpecies:
                newText+= 'setConcentration("%s",0)\n' % species
            newText+='simulate({method=>"ode",t_end=>1e7,n_steps=>500,atol=>1e-8,rtol=>1e-8})\n'
            for species in equilSpecies:
                newText+= 'setConcentration("%s",1e4)\n' % species
            newText+='simulate({suffix=>"equil",method=>"ode",t_start=>0,t_end=>36000,n_steps=>1000,atol=>1e-8,rtol=>1e-8});\n'
        else:
            newText = 'simulate({method=>"ode",t_start=>0,t_end=>36000,n_steps=>1000,atol=>1e-8,rtol=>1e-8});\n'

    elif method == "ssa":
       newText = "generate_network({overwrite=>1})\n"
       if equil:
            for species in equilSpecies:
                newText+= 'setConcentration("%s",0)\n' % species
            newText+='simulate({suffix=>"equil",method=>"ssa",t_end=>1e7,n_steps=>500})\n'
            for species in equilSpecies:
                newText+= 'setConcentration("%s",1e4)\n' % species
            newText+='simulate({method=>"ssa",t_start=>0,t_end=>36000,n_steps=>1000);\n'
       else:
            newText+='simulate({method=>"ssa",t_start=>0,t_end=>36000,n_steps=>1000});\n'

    elif method == "nf":
        if equil:
            for species in equilSpecies:
                newText+= 'setConcentration("%s",0)\n' % species
            newText+='simulate({suffix=>"equil",method=>"nf",t_end=>1e7,n_steps=>500})'
            for species in equilSpecies:
                newText+='setConcentration("%s",1e4)\n' % species
            newText+='simulate({method=>"nf",t_start=>0,t_end=>36000,n_steps=>1000);\n'
        else:
            newText = 'simulate{method=>"nf",t_start=>0,t_end=>36000,n_steps=>1000});\n'

    else:
        print("Invalid Simulation Method Called")
        
    if viz:
        newText+=('visualize({type=>"contactmap"})\n'
                 'visualize({type=>"regulatory",groups=>1,collapse=>1})\n')

    return newText
       
def addObservables(pysbModel,bound=False):
    i=1
    newModel = deepcopy(pysbModel)
    for monomer in newModel.monomers:
        for site in monomer.sites: 
            try:
                i+=1
                monomer.site_states[site]  #means a site has a flippable state 
                pattern = MonomerPattern(compartment=None,monomer=monomer,site_conditions={site:monomer.site_states[site][1]})   #Should see if [1] is always the perturbed/activated state
                newModel.observables.add(Observable(name='%s_%s_phos'%(monomer.name,site),reaction_pattern=pattern))

            #Bound observables are optionally available
            except KeyError:
                if bound:
                    pattern = MonomerPattern(compartment=None,monomer=monomer,site_conditions={site:WILD})
                    newModel.observables.add(Observable(name='%s_%s_bound'%(monomer.name,site),reaction_pattern=pattern))
                else:
                    pass

    return newModel










