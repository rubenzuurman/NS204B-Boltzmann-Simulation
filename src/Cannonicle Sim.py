

# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 11:20:47 2024

@author: Quint
"""
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

S = [] #matrix that gets filled with entropy values
SD_matrix = [] #matrix that gets filled with standart deviation entropy values

def thermalization(M, N, Nth): #thermalizing function
    for i in range(0,Nth):
        for i in range(0,M): 
            selector_1 = random.randint(0,M-1) #chooses donating bin
            selector_2 = random.randint(0,M-1) #chooses accepting bin
            if Main_lijst[selector_1] != 0: #checks if donating bin is not empty
                Main_lijst[selector_1] = Main_lijst[selector_1] - 1 #donating
                Main_lijst[selector_2] = Main_lijst[selector_2] + 1 #reciving


def randommovement(M,N,Ns): #movement algorithm
    for i in range(0,Ns):
        bincounts.append(np.bincount(Main_lijst)) #takes bincounts(macrostates) per sweep
        for i in range(0,M):
            selector_1 = random.randint(0,M-1) #chooses donating bin
            selector_2 = random.randint(0,M-1) #chooses accepting bin
            if Main_lijst[selector_1] != 0:  #checks if donating bin is not empty
                Main_lijst[selector_1] = Main_lijst[selector_1] - 1 #donating
                Main_lijst[selector_2] = Main_lijst[selector_2] + 1 #reciving
                

   
#loop setup                    
for index, N in enumerate(WaardenN):
    Srow = [] #rows in matrix S
    SDrow= [] #rows in matrix SD_matrix
    print(index/len(WaardenN)*100, '%')
    for M in WaardenM:
        Main_lijst = [] #lattice is created
        bincounts = [] # list where bincounts are collected
        
        Main_lijst = [0] * M #lattice is filled with lattice spots
        # print(Main_lijst)
        for i in range(0,N):
            Main_lijst[i % M] += 1
        
        thermalization(M, N, Nth) #thermelizes
        randommovement(M, N, Ns) #algorithm random movement

               

        

        
        max_bincount_length = max([len(nparr) for nparr in bincounts]) #takes bincount length 
        bincounts_summed = {i: 0 for i in range(max_bincount_length)} #dictonary maken
        
        for bincount in bincounts: #summs the bincounts
             for index, value in enumerate(bincount):
                bincounts_summed[index] += value
        
            
        max_bincount_length = max([len(nparr) for nparr in bincounts])
        bincounts2_summed = {i: 0 for i in range(max_bincount_length)} #dictonary maken
        for bincount in bincounts: #summs the bincounts quadraticly
            for index, value in enumerate(bincount):
                bincounts2_summed[index] += (value**2)

        
         
        
        N_av = np.array(list(bincounts_summed.values()))/Ns  #deelt de som van de histogrammen door het aantal histogrammen(sweeps)
        N2_av = np.array(list(bincounts2_summed.values()))/Ns #deelt de kwadraten som van de histogrammen door het aantal histogrammen(sweeps)
        SD = np.sqrt(N2_av - N_av**2) #berekend de standaart de
        Pn = N_av/M
        
        x = np.array(list(bincounts_summed.keys())) #aantal quanta bin 
        y = np.array(list(bincounts_summed.values())) #aantal bins die per hoeveelheid quanta
        Pn = [p for p in Pn if p != 0] #removes 0 values from consideration 0ln0 = 0 
        Srow.append(-K_B*sum(Pn*np.log(Pn))) #berekning entropie en voegt element toe aan row
        SDrow.append(SD)
    S.append(Srow)#voegt Srow toe aan entropie matrix
    SD_matrix.append(SDrow)#voegt SDrow toe aan SD matrix
    
S = np.array(S)

#%%

M, N = np.meshgrid(WaardenM, WaardenN)




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

# Makes 3D plot figure
fig = plt.figure(constrained_layout=True) 
ax = fig.add_subplot(111, projection='3d') 
 # Plot a 3D surface plot with M, N, and S data arrays, and use a color map
ax.plot_surface(M, N, S, cmap='inferno') 
# Set the view angle for better visualization of the 3D plot
ax.view_init(elev=30, azim= -65) 
# Labels voor de assen
ax.set_xlabel('M', size =12)
ax.set_ylabel('N', size =12)
ax.set_zlabel(r'$S$(J/K)',size =12, labelpad= 10)
# Adjust tick label size
ax.tick_params(labelsize = 11) 
ax.set_box_aspect(None, zoom =0.8) # to prevent cuttof of S axes
# Laat de plot zien
# plt.tight_layout()
plt.show()



def sliceplotNwise(n): #function to slice 3d plot
    plt.plot(WaardenM,S[n,:])
    plt.figtext(0.12, -0.05,'Entropy plotted against bins for ' + str(n) + ' quanta')
    plt.ylabel('Entropy ' r'$S$(J/K)')
    plt.xlabel('Amount of bins')
    plt.show()
    
def sliceplotMwise(m):  #function to slice 3d plot
    plt.plot(WaardenN,S[:,m])
    plt.ylabel('Entropy ' r'$S$(J/K)')
    plt.xlabel('Amount of quanta')
    plt.figtext(0.12, -0.05,f'Entropy plotted against quanta for {str(m)} bins')
    plt.show()
    
    
