import pandas as pd
import pickle

with open('dlrp_symbolonly.csv','r') as f:
    wholelist=f.readlines()

receptor_dict = {}

firstlig = 1
lig = 0
key1 = wholelist[0].strip()
val1 = []
newval = []
for line in wholelist[1:]:
    if firstlig == 1:
        if line != '\n':
            val1.append(line.strip())
        else:
            receptor_dict[key1] = val1
            firstlig = 0
    #Should never trip else until finished with the first ligand, and past the blank line    
    else:
        if lig == 0:
            newkey = line.strip()
            lig = 1
        else:
            if line != '\n':
                newval.append(line.strip())
            else:
                receptor_dict[newkey] = newval
                lig = 0
                newval = []


to_remove = []
for key in receptor_dict.keys():
    if receptor_dict[key] == []:
        to_remove.append(key)
for key in to_remove:
    receptor_dict.pop(key)

df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in receptor_dict.items() ]))
df = df.transpose()
df.to_csv('lig_rec.csv')

with open('lig_rec_df.pkl','wb') as f:
    testdict = pickle.dump(df,f)

with open('lig_rec_df.pkl','rb') as f:
    testdict = pickle.load(f)


with open('lig_rec_dict.pkl','wb') as f:
    testdict = pickle.dump(receptor_dict,f)

with open('lig_rec_dict.pkl','rb') as f:
    testdict = pickle.load(f)


