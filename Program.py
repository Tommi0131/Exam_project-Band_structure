import numpy as np
import pandas as pd
from numpy import linalg
#from numpy import linalg
import timeit
import matplotlib.pyplot as plt
%matplotlib inline
plt.style.use('seaborn-whitegrid')


Es= -8.0 #-8.0
Ep= 0.0 #0.0
Vss = -2.0 #-2.0
Vpp_sigma = +2.2
Vpp_pi = -1.8
Vsp = -2.1

#if you need to add orbitals, just add the (spx,spy...) and the orbitals in "orbitals"

def K_pathes(*path): #put a way to have just the G,X,M letters inserted!
    """function that takes 'G', 'X' OR 'M' letters (with up-commas)
    and return a tuple with k*a values along x and y directions"""

    path_dict = {'G' : (0,0), 'X' : (np.pi,0), 'M' : (np.pi,np.pi)}
    k_path_x=[]
    k_path_y=[]
    precision=150 #how many steps for each couple of letters
    #path_dict['G'][0]
    for index in range(len(path)-1): #we have 2 arrays, with the total k path in x and y dim. There is surely a better way!
        k_path_x = np.concatenate((k_path_x, np.linspace(path_dict[path[index]][0],path_dict[path[index+1]][0], precision)))
        k_path_y = np.concatenate((k_path_y, np.linspace(path_dict[path[index]][1],path_dict[path[index+1]][1], precision)))
    return k_path_x, k_path_y


def Ham_matrix(kx_a, ky_a): #qua potrei creare un oggetto matrice, cioè fare una classe.. invece che inizializzare ogni volta N matrici
    """it takes the comined orbitals as arguments (ss,spx,spy..)
    then creates a square matrix with the length of the orbitals
    then it becomes a DataFrame, with indexes and columns equal to orbitals names
    it inizializes with a double cycle each matrix element with the corresponding value (uses the name of the indexes)
    return the DataFrame inizialized
    return a matrix with values inizialized"""

    ss = Es + 2*Vss*(np.cos(kx_a)+np.cos(ky_a))
    spx = 2*Vsp*np.sin(kx_a)
    spy = 2*Vsp*np.sin(ky_a)
    pxpx = Ep + 2*(Vpp_sigma*np.cos(kx_a)+Vpp_pi*np.cos(ky_a))
    pypy = Ep + 2*(Vpp_sigma*np.cos(ky_a)+Vpp_pi*np.cos(kx_a))
    pzpz = Ep + 2*Vpp_pi*(np.cos(kx_a)+np.cos(ky_a))
    #all above could be put out

    orbitals = np.array(['s','px','py','pz'])
    basic_matrix = np.zeros(shape = (len(orbitals), len(orbitals)), dtype = float)
    final_matrix = pd.DataFrame(basic_matrix, index=orbitals, columns=orbitals)

    #now inizialize every matrix elements with the names "spx,spy.."
    """top funziona"""
    for i in orbitals:
        for j in orbitals:
            try:
                final_matrix.loc[i,j] = eval(i+j)
                final_matrix.loc[j,i] = eval(i+j)
            except NameError: #i think that this shouldn't be done.. it's nasty
                pass
    return final_matrix

Ham_matrix(0,0)

#k_a=np.linspace(0,np.pi,100) this is old!

def eigenvalues_path():
    """calculate the k pathes given the X,G,M labels (just those!)"""
    (kx_a, ky_a)= K_pathes('X','G','M','X','G') #forse è meglio tenerlo fuori dalla funzione
    eigenvals = []

    for i in range(len(kx_a)):
        """calculate eigenvalues for each matrix with a certain kx_a, ky_a values
        append the values to 'eigenvals'
        return eigenvalues"""
        eigenvals_temporary = linalg.eigvalsh(Ham_matrix(kx_a[i], ky_a[i]))
        eigenvals = np.concatenate((eigenvals, eigenvals_temporary))

    return eigenvals, len(kx_a)

bands, k_length = eigenvalues_path()
bands = np.reshape(bands, newshape=(k_length,4))
bands[34:39, :]
np.shape(bands)
k_length

#%% #to run use ctrl+alt+enter
"""graphical part"""
fig = plt.figure()
ax = plt.axes()

plt.title("Band structure")
plt.xlabel("k path")
plt.ylabel("Energy bands (eV)");

x=[0,100,200,300,400]
my_xticks = ['X','G','M','X','M']
plt.xticks(x, my_xticks)

ax.plot(range(k_length), bands[:,0], color='blue');
ax.plot(range(k_length), bands[:,1], color='blue');
ax.plot(range(k_length), bands[:,2], color='blue');
ax.plot(range(k_length), bands[:,3], color='blue');
plt.legend()
ax.set(ylim=(-24, 12))
#%%
