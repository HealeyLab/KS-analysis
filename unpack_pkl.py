import pickle
import sys
import numpy as np
path =  sys.argv[1]

with open(path, 'rb') as f:
    data = pickle.load(f)

keys, items = [], []
for key, val in data.items():
    keys.append(key[0])
    items.append(val)
    
print(np.array([keys, items]))
