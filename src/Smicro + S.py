# Names: Quinten Schuttevaer, Ruben Zuurman
# Runtime on Ryzen 9 5900x: 52.5 seconds

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 22:53:28 2024

@author: Quint
"""
#Dit script berekend Smicro het verschil tussen Stotaal,zoals in het microcanonicle, en S zoals in het canonicle
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
Smicromatrix = []

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
    Srow = []
    Smicrow= []
    print(index/len(WaardenN)*100, '%')
    for M in WaardenM:
        Main_lijst = []   
        bincounts = []
        
        Main_lijst = [0] * M
        # print(Main_lijst)
        for i in range(0,N):
            Main_lijst[i % M] += 1
        
        thermalization(M, N, Nth) #thermalizes the 
        randommovement(M, N, Ns) #fills bincounts with the bincounts of the sim between sweeps

               

        

        
        max_bincount_length = max([len(nparr) for nparr in bincounts])
        bincounts_summed = {i: 0 for i in range(max_bincount_length)} #dictonary maken
        
        for bincount in bincounts:
             for index, value in enumerate(bincount):
                bincounts_summed[index] += value
        # print(bincounts_summed)
        # print(1)        
        
        
            
        max_bincount_length = max([len(nparr) for nparr in bincounts])
        bincounts2_summed = {i: 0 for i in range(max_bincount_length)}
        for bincount in bincounts:
            for index, value in enumerate(bincount):
                bincounts2_summed[index] += (value**2)
        # print(bincounts2_summed)
        
         
        
        N_av = np.array(list(bincounts_summed.values()))/Ns  #deelt de som van de histogrammen door het aantal histogrammen(sweeps)
        N2_av = np.array(list(bincounts2_summed.values()))/Ns #deelt de kwadraten som van de histogrammen door het aantal histogrammen(sweeps)
        Pn = N_av/M
        Smicro = []
        
        for n in range(len(Pn)):
            Smicro.append(Pn[n]*K_B*np.log(sp.comb(N + M -1 - n, M - 2)))            
        Smicro = sum(Smicro)
        x = np.array(list(bincounts_summed.keys())) #aantal quanta bin 
        y = np.array(list(bincounts_summed.values())) #aantal bins die per hoeveelheid quanta
        Pn = [p for p in Pn if p != 0] #removes 0 values from consideration 0ln0 = 0 
        Srow.append(-K_B*sum(Pn*np.log(Pn))) #berekning entropie en voegt element toe aan row
        Smicrow.append(Smicro)
    S.append(Srow)#voegt Srow toe aan entropie matrix
    Smicromatrix.append(Smicrow)#voegt SDrow toe aan SD matrix
    
S = np.array(S)
Smicromatrix = np.array(Smicromatrix)
#%% The plotting part

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
      
#Creates 3D plot with tight layout to prevent clipping
fig = plt.figure(constrained_layout=True) 
ax = fig.add_subplot(111, projection='3d')
 
#Actualy plots the surface based on Stot=S+Smicromatrix
ax.plot_surface(M, N, S + Smicromatrix , cmap='magma') 
#Sets the viewing angle
ax.view_init(elev=30, azim= -65) 
# Labels voor de assen
ax.set_xlabel('M', size =12)
ax.set_ylabel('N', size =12)
ax.set_zlabel(r'$S$(J/K)',size =12, labelpad= 10)#labelpad to prevent clipping with the axes
ax.tick_params(labelsize = 11)
ax.set_box_aspect(None, zoom =0.8)
# Laat de plot zien
# plt.tight_layout()
plt.show()



def sliceplotNwise(n): #function to slice 3d plot
    plt.plot(WaardenM,S[n,:])
    plt.figtext(0.12, -0.05,'Entropy plotted against bins for ' + str(n) + ' quanta') #description of plot
    #labels
    plt.ylabel('Entropy ' r'$S$(J/K)')
    plt.xlabel('Amount of bins')
    plt.show() #shows plot
    
def sliceplotMwise(m): #function to slice 3d plot
    plt.plot(WaardenN,S[:,m])
    #labels
    plt.ylabel('Entropy ' r'$S$(J/K)')
    plt.xlabel('Amount of quanta')
    plt.figtext(0.12, -0.05,f'Entropy plotted against quanta for {str(m)} bins') #description of plot
    plt.show() #shows plot
