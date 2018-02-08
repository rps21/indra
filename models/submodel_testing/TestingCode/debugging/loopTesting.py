   
downstream_list = stmts[-4].sub
downstream_list = [downstream_list]
upstream_list = stmts[-2].agent
upstream_list = [upstream_list]
updated_stmts = smallestset

#upstream_list
#Out[55]: BRAF()

#downstream_list
#MAP2K1(mods: (phosphorylation, S, 222), (phosphorylation, S, 218))


i=0
outputStmts = []
while downstream_list:
    next_stmts, new_downstream_list, new_upstream_list = cs.find_next_cascade_step(updated_stmts,upstream_list, downstream_list)
    print(upstream_list)
    print(downstream_list)
    new_af_stmts = cs.add_active_forms(next_stmts,new_upstream_list)
    updated_stmts = updated_stmts + new_af_stmts

#    #fixing loop issue
    newStmts = [st for st in updated_stmts if st not in outputStmts]
    outputStmts = newStmts + outputStmts
    updated_stmts = [st for st in updated_stmts if st not in next_stmts]
#    #end fix 

    upstream_list_total.append(upstream_list)
    downstream_list_total.append(downstream_list)

    #Reset for next iteration
    downstream_list = new_downstream_list
    upstream_list = new_upstream_list   

    print(i)
    if i >= 30: #Need better option to handle feedbacks leading to infinite loops
        break
    i=i+1


[Phosphorylation(BRAF(mods: (phosphorylation, T, 373)), MAP2K1(), S, 222),
 Phosphorylation(MAP2K1(mods: (phosphorylation, S, 218), (phosphorylation, S, 222), (phosphorylation, T, 286), (phosphorylation, T, 292)), MAPK1(), T, 185),
 Phosphorylation(MAPK1(mods: (phosphorylation, Y, 187), (phosphorylation, T, 185)), BRAF(mods: (phosphorylation, S, 446), (phosphorylation, T, 599), (phosphorylation, S, 602)), T, 401),
 ActiveForm(MAPK1(mods: (phosphorylation, T, 185)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 401)), kinase, True),
 ActiveForm(MAP2K1(mods: (phosphorylation, S, 222)), kinase, True),
 ActiveForm(MAPK1(mods: (phosphorylation, T, 185)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 401)), kinase, True),
 ActiveForm(MAP2K1(mods: (phosphorylation, S, 222)), kinase, True),
 ActiveForm(MAPK1(mods: (phosphorylation, T, 185)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 401)), kinase, True),
 ActiveForm(MAP2K1(mods: (phosphorylation, S, 222)), kinase, True),
 ActiveForm(MAPK1(mods: (phosphorylation, T, 185)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 401)), kinase, True),
 ActiveForm(MAP2K1(mods: (phosphorylation, S, 222)), kinase, True),
 ActiveForm(MAPK1(mods: (phosphorylation, T, 185)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 401)), kinase, True),
 ActiveForm(MAP2K1(mods: (phosphorylation, S, 222)), kinase, True),
 ActiveForm(MAPK1(mods: (phosphorylation, T, 185)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 401)), kinase, True),
 ActiveForm(MAP2K1(mods: (phosphorylation, S, 222)), kinase, True),
 ActiveForm(MAPK1(mods: (phosphorylation, T, 185)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 401)), kinase, True),
 ActiveForm(MAP2K1(mods: (phosphorylation, S, 222)), kinase, True),
 ActiveForm(MAPK1(mods: (phosphorylation, T, 185)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 401)), kinase, True),
 ActiveForm(MAP2K1(mods: (phosphorylation, S, 222)), kinase, True),
 ActiveForm(MAPK1(mods: (phosphorylation, T, 185)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 401)), kinase, True),
 ActiveForm(MAP2K1(mods: (phosphorylation, S, 222)), kinase, True),
 ActiveForm(MAPK1(mods: (phosphorylation, T, 185)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 401)), kinase, True),
 ActiveForm(MAP2K1(mods: (phosphorylation, S, 222)), kinase, True),
 ActiveForm(MAPK1(mods: (phosphorylation, T, 185)), kinase, True)]


########################
#NEW WIN
#AFTER FIX:
 [73]: outputStmts
Out[73]: 
[ActiveForm(MAP2K1(mods: (phosphorylation, S, 222)), kinase, True),
 ActiveForm(BRAF(mods: (phosphorylation, T, 401)), kinase, True),
 Phosphorylation(BRAF(mods: (phosphorylation, T, 373)), MAP2K1(), S, 222),
 Phosphorylation(MAP2K1(mods: (phosphorylation, S, 218), (phosphorylation, S, 222), (phosphorylation, T, 286), (phosphorylation, T, 292)), MAPK1(), T, 185),
 Phosphorylation(MAPK1(mods: (phosphorylation, Y, 187), (phosphorylation, T, 185)), BRAF(mods: (phosphorylation, S, 446), (phosphorylation, T, 599), (phosphorylation, S, 602)), T, 401),
 ActiveForm(MAPK1(mods: (phosphorylation, T, 185)), kinase, True)]

