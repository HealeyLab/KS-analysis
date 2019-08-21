import pickle
import sys
import numpy as np
# sys.argv[1] #
path =  'D:\DJP thesis sorting\mda 8 14 19\mda 8 14 19 habituate no fem_190814_165845_Kilosort1850\.phy\memcache\phycontrib.template.gui.get_best_channel.pkl' # sys.argv[1]

with open(path, 'rb') as f:
    data = pickle.load(f)

keys, items = [], []
for key, val in data.items():
    keys.append(key[0])
    items.append(val)
    
print(np.array([keys, items]))
