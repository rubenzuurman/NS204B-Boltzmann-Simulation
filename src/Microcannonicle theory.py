# Runtime on Ryzen 9 5900x: <1 second

# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:30:04 2024

@author: Quint
"""


#Dit script berekend Stotaal zoals in het microcanonicle
import numpy as np
import matplotlib.pyplot as plt
import random
import scipy.odr as odr
from mpl_toolkits.mplot3d import Axes3D
import math
import scipy.special as sp
random.seed(1234)

WaardenN = range(0,100, 1)
WaardenM = range(1,100, 1)

Ns = 100 #aantal simulatie sweeps
Nth = 50 # aantal thermalizatie sweeps
K_B = 1.380649e-23
# K_B = 1
S = []

#loop setup to itarate over N and M                   
for index, N in enumerate(WaardenN):
    Srow = []# row in S
    print(index/len(WaardenN)*100, '%')
    for M in WaardenM:
        Srow.append(K_B*np.log(sp.comb(N + M -1 , M - 1)))  #calculates S for given N and add to row in S        
    S.append(Srow)#voegt Srow toe aan entropie matrix

    
S = np.array(S)

#%% The plotting part

M, N = np.meshgrid(WaardenM, WaardenN) #creates meshgrid to for 3d plot based on WaardenM and WaardenN




# Stap 5: Plot het oppervlak
# for i in range(40, 360, 20):
#     print(i)
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     # ax.plot_surface(X, Y, correlatieS, cmap='magma')  # gebruik de 'viridis' kleurenkaart
#     # ax.plot_surface(X, Y, z_data, cmap='viridis') 
#     ax.plot_surface(X, Y, S, cmap='inferno') 
#     ax.view_init(elev=30, azim=i) 
#     # Labels voor de assen
#     ax.set_xlabel('M as')
#     ax.set_ylabel('N as')
#     ax.set_zlabel('S as')
#     ax.dist = 12
#     # Laat de plot zien
#     plt.show()
        
fig = plt.figure(constrained_layout=True)  
ax = fig.add_subplot(111, projection='3d')
 
ax.plot_surface(M, N, S, cmap='magma') 
ax.view_init(elev=30, azim= -65) 
# Labels voor de assen
ax.set_xlabel('M', size =12)
ax.set_ylabel('N', size =12)
ax.set_zlabel(r'$S$(J/K)',size =12, labelpad= 10) #labelpad to prevent clipping with the axes
#increases size of the axes numbers
ax.tick_params(labelsize = 11)
#zooms out to prevent clipping
ax.set_box_aspect(None, zoom =0.8)
# Laat de plot zien
# plt.tight_layout()
plt.show()



def sliceplotNwise(n): #function to slice 3d plot
    plt.plot(WaardenM,S[n,:]) #plots the slice of the matrix S
    plt.figtext(0.12, -0.05,'Entropy plotted against bins for ' + str(n) + ' quanta') #description of plot
    #labels for the axes
    plt.ylabel('Entropy ' r'$S$(J/K)')
    plt.xlabel('Amount of bins')
    plt.show() #shows plot
    
def sliceplotMwise(m): #function to slice 3d plot
    plt.plot(WaardenN,S[:,m]) #plots the slice of the matrix S
    #labels fot the axes
    plt.ylabel('Entropy ' r'$S$(J/K)')
    plt.xlabel('Amount of quanta')
    plt.figtext(0.12, -0.05,f'Entropy plotted against quanta for {str(m)} bins') #description of plot
    plt.show() #shows plot
