import pickle
import pandas as pd
import scipy as sp
import os


objfolder ='.'

def dump_object(obj, filename, path=objfolder):
    filename = os.path.join(path, filename)
    print(filename)
    with open(filename, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_object(filename, path=objfolder):
    filename = os.path.join(path, filename)
    with open(filename, 'rb') as f:
        obj = pickle.load(f)
    return obj


res = load_object('res_2.pkl')


print(res.describe())

print(res['phi0'])
print(res['phi1'])
print(res['ncycle'])