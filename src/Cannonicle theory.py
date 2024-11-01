# Runtime on Ryzen 9 5900x: 4.6 seconds

import numpy as np
import matplotlib.pyplot as plt
import random
import scipy.odr as odr
from mpl_toolkits.mplot3d import Axes3D
import scipy.special as sp
import math

#Dit script berekend de kansen op bepaalde macrostates en bepaalt daaruit de cannoncical S

# Aantal quanta en bins
N = np.array(range(0, 100, 1))  # waarden N
M = np.array(range(1, 100, 1))  # waarden M
K_B = 1.380649e-23

# Initialiseer de volledige matrix
Fullmatrix = []
N_grid, M_grid = np.meshgrid(N, M)
# Itereert over gekozen N en M
for Mchosen in M:
    Row = []
    for Nchosen in N:
        Pn = [] # a list of chances 
        
        # Itereer over alle mogelijke n waarden (0 tot Nchosen)
        for n in range(Nchosen + 1):
            # Bereken de multipliciteit voor deze specifieke n
            omega = sp.comb(Mchosen + Nchosen - 1, Mchosen - 1, exact=False)
            if omega > 0:
                Pn_value = sp.comb(Mchosen + Nchosen - 2 - n, Mchosen - 2, exact=False) / omega #This is calculates the likilyhood of a macrostate
            else:
                Pn_value = 0
            Pn.append(Pn_value)
        Row.append(sum(-K_B*np.array(Pn)*np.log(np.array(Pn))))  #This adds a value of S to the row 
    Fullmatrix.append(Row)  # Voeg de rij toe aan de volledige matrix


# Convert Fullmatrix naar een numpy array als je dat nodig hebt
Fullmatrix = np.array(Fullmatrix)

#%% The plotting part
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')




# Plot het fitte vlak
ax.plot_surface(M_grid, N_grid, Fullmatrix,  cmap='magma')
ax.view_init(elev=30, azim= -65) 
# Labels voor de assen
ax.set_xlabel('M', size =12)
ax.set_ylabel('N', size =12)
ax.set_zlabel(r'$S$(J/K)   .',size =12, labelpad= 6) #labelpad to prevent clipping with the axes
#increases size of the axes numbers
ax.tick_params(labelsize = 11)
#zooms out to prevent clipping
ax.set_box_aspect(None, zoom =0.80)
# Laat de plot zien
plt.show()
# Print de grootte van de matrix
print("Shape of Fullmatrix:", Fullmatrix.shape)
print("Fullmatrix:")
print(Fullmatrix)

def sliceplotNwise(n): #functie om een slice van de 3d plot weer te geven
    plt.plot(M,Fullmatrix[n,:])  #plots the slice of the matrix S    
    plt.figtext(0.12, -0.05,('Entropy plotted against bins for ' + str(n) + ' quanta'))
    #labels
    plt.ylabel('Entropy ' r'$S$(J/K)')
    plt.xlabel('Amount of bins')
    plt.show()
    
def sliceplotMwise(m): #functie om een slice van de 3d plot weer te geven
    plt.plot(N,Fullmatrix[:,m])  #plots the slice of the matrix S
    #labels
    plt.ylabel('Entropy ' r'$S$(J/K)')
    plt.xlabel('Amount of quanta')
    plt.figtext(0.12, -0.05,('Entropy plotted against quanta for ' + str(m) + ' bins'))
    plt.show()
    
   