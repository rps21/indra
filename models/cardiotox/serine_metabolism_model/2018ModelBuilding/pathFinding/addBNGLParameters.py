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



