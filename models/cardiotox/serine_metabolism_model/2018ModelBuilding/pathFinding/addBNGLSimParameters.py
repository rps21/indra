def addSimParamters(method='ode',equil=True,equilSpecies=[],viz=True):

    if method == "ode":
        newText = """generate_network({overwrite=>1})
        writeMexfile()\n"""
        if equil:
            for species in equilSpecies:
                newText+= "setConcentration(%s,0)\n" % species
            newText+="""simulate({method=>'ode',t_end=>1e7,n_steps=>500,atol=>1e-8,rtol=>1e-8})
            setConcentration(,1e4)
            simulate({method=>'ode',t_start=>0,t_end=>36000,n_steps=>1000,atol=>1e-8,rtol=>1e-8});\n"""
        else:
            newText = "simulate({method=>'ode',t_start=>0,t_end=>36000,n_steps=>1000,atol=>1e-8,rtol=>1e-8});\n"

    elif method == "ssa":
       newText = "generate_network({overwrite=>1})\n"
       if equil:
            for species in equilSpecies:
                newText+= "setConcentration(%s,0)\n" % species
            newText+="""simulate({method=>'ssa',t_end=>1e7,n_steps=>500})
            setConcentration(,1e4)
            simulate({method=>'ssa',t_start=>0,t_end=>36000,n_steps=>1000);\n"""
       else:
            newText+="simulate({method=>'ssa',t_start=>0,t_end=>36000,n_steps=>1000});\n"

    elif method == "nf":
        if equil:
            for species in equilSpecies:
                newText+= "setConcentration(%s,0)\n" % species
            newText+="""setConcentration(,0)
            "simulate({method=>'nf',t_end=>1e7,n_steps=>500})
            "setConcentration(,1e4)
            "simulate_({method=>'nf',t_start=>0,t_end=>36000,n_steps=>1000);\n"""
        else:
            newText = "simulate{method=>'nf',_start=>0,t_end=>36000,n_steps=>1000});\n"

    else:
        print("Invalid Simulation Method Called")
        
    if viz:
        newText+="visualize({type=>'contactmap'})\n"
        newText+="visualize({type=>'regulatory',groups=>1,collapse=>1})\n"

    return newText


#Issue will be defining species for equilibration, observables.
#Can pass as a string, but not ideal for 
#Could write a function to parse seed speices based on species name 
