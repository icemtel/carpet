import  scipy as sp
import  pickle
import scipy.linalg as lin
import matplotlib.pyplot as plt

with open('phi_2.pkl', 'rb') as f:
    phis = pickle.load(f)


print(sp.array(phis).shape)


plt.plot(lin.norm(phis, axis=1))
plt.show()