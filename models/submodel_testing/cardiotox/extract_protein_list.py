#loop through stmts
#pull both agents and add to list

#remove obvious things
#redundancies, mistakes, disconnected mini networks

#see how big list is

#try to remove things not in proteomics

new_lists = []
for st in fresh_stmts:
    new_lists = new_lists + st.agent_list()

str_list = []
for element in new_lists:
    str_list.append(str(element))

unique_protein_list = list(set(str_list))

with open('protein_list.txt','w') as newfile:
    newfile.write('\n'.join(unique_protein_list))


