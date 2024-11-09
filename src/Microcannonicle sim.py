# Names: Quinten Schuttevaer, Ruben Zuurman
# Runtime on Ryzen 9 5900x: 2.8 seconds

# -*- coding: utf-8 -*-


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
from collections import Counter

# random.seed(1234)

#Dit script simuleerd en systeem om te proberen alle microstates te identificeren en zo de microcannonicle S te berekenen


WaardenN = range(0,5, 1) #genration of a list of quanta
WaardenM = range(1,5, 1) #genration of a list of lattice sites

Ns = 50000 #aantal simulatie sweeps
Nth = 50 # aantal thermalizatie sweeps
K_B = 1.380649e-23

S = [] #matrix that gets filled with entropy values

def thermalization(M, N, Nth): #thermalizing function
    for i in range(0,Nth):
        for i in range(0,M): 
            selector_1 = random.randint(0,M-1) #chooses donating bin
            selector_2 = random.randint(0,M-1) #chooses accepting bin
            if Main_lijst[selector_1] != 0: #checks if donating bin is not empty
                Main_lijst[selector_1] = Main_lijst[selector_1] - 1 #donating
                Main_lijst[selector_2] = Main_lijst[selector_2] + 1 #reciving


def randommovement(M,N,Ns):
    bincounts = [] #list that gets filled with current states between sweeps
    for i in range(0,Ns):
        current_state = Main_lijst[:] #takes the mainlist as current state between sweeps
        bincounts.append(current_state)#takes current state and adds it to bincounts
        for i in range(0,M):
            selector_1 = random.randint(0,M-1) #chooses donating bin
            selector_2 = random.randint(0,M-1) #chooses accepting bin
            if Main_lijst[selector_1] != 0: #checks if donating bin is not empty
                Main_lijst[selector_1] = Main_lijst[selector_1] - 1 #the donating
                Main_lijst[selector_2] = Main_lijst[selector_2] + 1 #the reciving
    return bincounts           


#loop setup                    
for index, N in enumerate(WaardenN):
    Srow = [] #rows in S matrix
    SDrow= []  #rows in SD_matrix
    print(index/len(WaardenN)*100, '%')
    for M in WaardenM:
        Main_lijst = [] #lattice is created        
        Main_lijst = [0] * M #latice sites are created
        for i in range(0,N): #fills sites with quanta
            Main_lijst[i % M] += 1
        thermalization(M, N, Nth) #thermilezes distrabution
        
        seen = set() #makes a set in which already seen microstates will be recorded
        microstates = [] #list of unique microstates
        alle_states = randommovement(M, N, Ns) #all states given by the random movement function
        for sublist in alle_states:
            # Voeg de tuple van de sublijst toe aan de set en de sublijst zelf aan de unieke lijsten
            tuple_sublist = tuple(sublist)
            if tuple_sublist not in seen: #if the microstate is not in seen and thus new it is added to seen and to the list of unique microstates
                seen.add(tuple_sublist)
                microstates.append(sublist)
               

        

        
      
        omega = len(microstates) #omega is the amount of microstates
        Srow.append(K_B*np.log(omega)) #berekning entropie en voegt element toe aan row
    S.append(Srow)#voegt Srow toe aan entropie matrix
    
S = np.array(S) #makes S a np.rray



#%% The plotting part


# for i in range(0,360, 20):
# Stap 1: x- en y-waarden definiëren
x = WaardenM
y = WaardenN

# Stap 2: Maak een grid voor x en y
X, Y = np.meshgrid(x, y)


#Maak de 3D-plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#plots the surface
ax.plot_surface(X, Y, S, cmap='inferno') 
#Geeft de weergave hoek van de 3D plot
ax.view_init(elev=30, azim= -65) 
# Labels voor de assen
ax.set_xlabel('M', size =12)
ax.set_ylabel('N', size =12)
ax.set_zlabel(r'$S$(J/K)  .',size =12, labelpad= 10) #labelpad to prevent clipping with the axes
#Makes the axes numbers larger
ax.tick_params(labelsize = 11)
#zooms out to prevent clipping
ax.set_box_aspect(None, zoom =0.80)
# Laat de plot zien
# plt.tight_layout()
plt.show()


def sliceplotNwise(n): #function to slice 3d plot
    plt.plot(WaardenM,S[n,:])  #plots the slice of the matrix S
    plt.figtext(0.12, -0.05,'Entropy plotted against bins for ' + str(n) + ' quanta') #description of plot
    #labels
    plt.ylabel('Entropy ' r'$S$(J/K)')
    plt.xlabel('Amount of bins')
    plt.show() #shows plot
    
def sliceplotMwise(m): #function to slice 3d plot
    plt.plot(WaardenN,S[:,m])  #plots the slice of the matrix S
    #labels
    plt.ylabel('Entropy ' r'$S$(J/K)')
    plt.xlabel('Amount of quanta')
    plt.figtext(0.12, -0.05,f'Entropy plotted against quanta for {str(m)} bins') #description of plot
    plt.show() #shows plot
    