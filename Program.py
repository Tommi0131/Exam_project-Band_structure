
import numpy as np
import pandas as pd
from numpy import linalg
import matplotlib.pyplot as plt
%matplotlib inline
plt.style.use('seaborn-whitegrid')

"""these need to be put in a different script"""
Es= -8.0
Ep= 0.0
Vss = -2.0
Vpp_sigma = +2.2
Vpp_pi = -1.8
Vsp = -2.1


def K_pathes(*path): #put a way to have just the G,X,M letters inserted!
    """takes 'G', 'X' OR 'M' as chars (with up-commas)
    create through 'path_dict' the k pathes.
    return a tuple with k*a values along x and y directions"""

    path_dict = {'G' : (0,0), 'X' : (np.pi,0), 'M' : (np.pi,np.pi)}

    k_path_x=[]
    k_path_y=[]
    precision=100 #how many steps for each couple of letters
    for index in range(len(path)-1):
        k_path_x = np.concatenate((k_path_x, np.linspace(path_dict[path[index]][0],path_dict[path[index+1]][0], precision)))
        k_path_y = np.concatenate((k_path_y, np.linspace(path_dict[path[index]][1],path_dict[path[index+1]][1], precision)))
    return k_path_x, k_path_y


def Ham_matrix(kx_a, ky_a): #qua potrei creare un oggetto matrice, cioè fare una classe.. invece che inizializzare ogni volta N matrici
    """takes the k_pathes along x and y (multiplied by Ham vector of reciprocal lattice, 'a')
    inizitalizes the parameters (Hamiltonian matrix elements) with kx_a and ky_a values
    creates a DataFrame with indices and columns as 'orbitals' (s, px, py, pz)
    inizializes the Dataframe by orbitals' names with the parameters
    return the DataFrame"""

    #parameters to be set in Hamiltonian matrix
    orbitals_dict = {
    'ss' : Es + 2*Vss*(np.cos(kx_a)+np.cos(ky_a)),
    'spx' : 2*Vsp*np.sin(kx_a),
    'spy' : 2*Vsp*np.sin(ky_a),
    'pxpx' : Ep + 2*(Vpp_sigma*np.cos(kx_a)+Vpp_pi*np.cos(ky_a)),
    'pypy' : Ep + 2*(Vpp_sigma*np.cos(ky_a)+Vpp_pi*np.cos(kx_a)),
    'pzpz' : Ep + 2*Vpp_pi*(np.cos(kx_a)+np.cos(ky_a))
    }

    #create the Dataframe that have the parameters
    orbitals = np.array(['s','px','py','pz'])
    matrix_order = len(orbitals)
    matrix = np.zeros(shape = (matrix_order, matrix_order))
    matrix = pd.DataFrame(matrix, index=orbitals, columns=orbitals)

    #inizialize Hamiltonian matrix thorugh parameters' name
    for i in orbitals:
        for j in orbitals:
            mixed_orbitals = i+j
            matrix_elements = orbitals_dict.get(mixed_orbitals)
            if (matrix_elements == None):
                print(f"{mixed_orbitals} has null values")
            else:
                matrix.loc[i,j] = orbitals_dict[mixed_orbitals]
                matrix.loc[j,i] = orbitals_dict[mixed_orbitals]
                """
            try:
                matrix.loc[i,j] = orbitals_dict[mixed_orbitals]
                matrix.loc[j,i] = orbitals_dict[mixed_orbitals]
            except KeyError: #i think that this shouldn't be done.. it's nasty
                pass
                """
    return matrix


def eigenvalues_path(*critical_points): #put a variable that can be inserted by keyboard or script
    """take char values (only 'X','G' or 'M')
    calculate the k pathes
    calculate the eigenvalues for each k point and concatenate them
    reshape such that each row has eingevalues relative to just a k point
    return eigenvalues and length of k path"""
    #I DO NOT LIKE THAT INPUT IN K_pathes IS THE SAME OF THE MAIN FUNCTION! BAD!
    #create k pathes
    (kx_a, ky_a)= K_pathes(*critical_points)

    #calculate the eigenvalues for each k point
    eigenvalues = []
    for i in range(len(kx_a)):
        eigenvalues_temporary = linalg.eigvalsh(Ham_matrix(kx_a[i], ky_a[i]))
        eigenvalues = np.concatenate((eigenvalues, eigenvalues_temporary))

    n_of_k_points = len(kx_a)
    n_orbitals = 4 #maybe i can put 'orbitals' global
    eigenvalues = np.reshape(eigenvalues, newshape=(n_of_k_points, n_orbitals)) #each row has eigenvalues corresponding to a k point

    return eigenvalues, n_of_k_points

eigenvalues, n_of_k_points = eigenvalues_path('X','G','M','X','G')

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

ax.plot(range(n_of_k_points), eigenvalues[:,0], color='blue');
ax.plot(range(n_of_k_points), eigenvalues[:,1], color='blue');
ax.plot(range(n_of_k_points), eigenvalues[:,2], color='blue');
ax.plot(range(n_of_k_points), eigenvalues[:,3], color='blue');
plt.legend()
#ax.set(ylim=(-24, 12))
#%%
