def addSimParamters(method='ode',equil=True,equilSpecies=[],viz=True):

    if method == "ode":
        newText = """generate_network({overwrite=>1})
        writeMexfile()\n"""
        if equil:
            for species in equilSpecies:
                newText+= 'setConcentration("%s",0)\n' % species
            newText+='simulate({method=>"ode",t_end=>1e7,n_steps=>500,atol=>1e-8,rtol=>1e-8})\n'
            for species in equilSpecies:
                newText+= 'setConcentration("%s",1e4)\n' % species
            newText+='simulate({method=>"ode",t_start=>0,t_end=>36000,n_steps=>1000,atol=>1e-8,rtol=>1e-8});\n'
        else:
            newText = 'simulate({method=>"ode",t_start=>0,t_end=>36000,n_steps=>1000,atol=>1e-8,rtol=>1e-8});\n'

    elif method == "ssa":
       newText = "generate_network({overwrite=>1})\n"
       if equil:
            for species in equilSpecies:
                newText+= 'setConcentration("%s",0)\n' % species
            newText+='simulate({method=>"ssa",t_end=>1e7,n_steps=>500})\n'
            for species in equilSpecies:
                newText+= 'setConcentration("%s",1e4)\n' % species
            newText+='simulate({method=>"ssa",t_start=>0,t_end=>36000,n_steps=>1000);\n'
       else:
            newText+='simulate({method=>"ssa",t_start=>0,t_end=>36000,n_steps=>1000});\n'

    elif method == "nf":
        if equil:
            for species in equilSpecies:
                newText+= 'setConcentration("%s",0)\n' % species
            newText+='simulate({method=>"nf",t_end=>1e7,n_steps=>500})'
            for species in equilSpecies:
                newText+='setConcentration("%s",1e4)\n' % species
            newText+='simulate({method=>"nf",t_start=>0,t_end=>36000,n_steps=>1000);\n'
        else:
            newText = 'simulate{method=>"nf",t_start=>0,t_end=>36000,n_steps=>1000});\n'

    else:
        print("Invalid Simulation Method Called")
        
    if viz:
        newText+="""visualize({type=>"contactmap"})
        visualize({type=>"regulatory",groups=>1,collapse=>1})\n"""

    return newText


#Issue will be defining species for equilibration, observables.
#Can pass as a string, but not ideal for 
#Could write a function to parse seed speices based on species name 

#In [24]: species_codes = [pysb.generator.bng.format_complexpattern(cp) for cp, param in originalModel.initial_conditions]

#In [25]: species_codes
#Out[25]: 
#['GRB2(erbb)',
# 'SOS1(erbb)',

#species codes - should be useful for grabbing species for setConcentration



#Want as observables:
#SRC(vegfr?~p)   SRC(phospho='p')
#Maybe 
#GRB2(erbb!+)    GRB2(erbb=1)
#SOS1(erbb!+)    SOS1(erbb=1)

from pysb import *
from pysb.core import MonomerPattern
#This is acting on direct issue, which we don't want
#Tweak naming        
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


#srcMon = originalModel.monomers[2]
#srcPhosSite = srcMon.sites[1]   #[0] is binding 
#srcPhosState = srcMon.site_states[srcPhosSite][1] #[0] is unphos
#srcPhosPattern = pysb.core.MonomerPattern(compartment=None,monomer=srcMon,site_conditions={srcPhosSite:srcPhosState})

#srcBindSite = srcMon.sites[0]
##srcBindState = srcMon.site_states[srcBindSite]  #doesn't exist
#srcBindPattern = pysb.core.MonomerPattern(compartment=None,monomer=srcMon,site_conditions={srcBindSite:1})

















