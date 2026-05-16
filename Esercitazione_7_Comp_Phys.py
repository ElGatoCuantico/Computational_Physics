#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 21:44:08 2026

@author: alberto
"""

import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science', 'notebook', 'grid'])   #migliora la visibilità



#Esercizio 1

N = 100_000

x = np.random.rand(N)

N_valori = np.arange(1, N + 1)

media_cumulativa = np.cumsum(x) / N_valori

errore_teorico = 1/np.sqrt(12*N_valori)

media_teorica = 0.5

plt.plot(N_valori, media_cumulativa, label='Media cumulativa numerica')
plt.axhline(y=media_teorica, label='Media teorica = 1/2')
plt.plot(N_valori, media_teorica + errore_teorico, label=r'Media teorica + 1/$\sqrt{12N}$ ')
plt.plot(N_valori, media_teorica - errore_teorico, label=r'Media teorica - 1/$\sqrt{12N}$ ')
plt.xlabel('N')
plt.xscale('log')
plt.ylabel('Media')
plt.legend()
plt.show()


#Esercizio 2


def f(x):
    
    return 4*np.sqrt(1-x**2)

x = np.random.rand(N)

N_valori = np.arange(1, N + 1)

stima_integrale = (1/N_valori) * np.cumsum(f(x))

dev_stnd_f = np.std(f(x))
sigma_I = dev_stnd_f/np.sqrt(N_valori)

pi_th = np.pi
    
plt.plot(N_valori, stima_integrale, label=r'Stima MC $\int_{0}^{1}4\sqrt{1-x^2}dx$ ')
plt.axhline(y=pi_th, label=r'$\pi$ teorico')
plt.plot(N_valori, pi_th + sigma_I, label=r'$\pi$ teorico + $\sigma_f/\sqrt{N}$ ')
plt.plot(N_valori, pi_th - sigma_I, label=r'$\pi$ teorico - $\sigma_f/\sqrt{N}$ ')
plt.xlabel('N')
plt.xscale('log')
plt.ylabel(r'$\pi$')
plt.legend()
plt.show()


#Esercizio 3


def inizializza_reticolo(L, tipo="random"):
    
    reticolo = 0
    
    if tipo == "random":
        
        reticolo = np.random.choice([1, -1], size=(L, L)) 
        #permette di estrarre elementi casuali da un 
        #array fornito in input, e accetta anche un parametro size 
        #per definire la forma della matrice di uscita
        
    elif tipo == "ordinato":
        
        reticolo = np.ones((L,L))
    
    return reticolo


def calcola_magnetizzazione(reticolo):
    
    return np.sum(reticolo)
    




#np.roll()   
#permette di far "scorrere" gli elementi di una matrice lungo un asse, 
#riportando quelli che escono da un lato all'inizio dell'altro



def calcola_energia(reticolo, J):
    
    #sposta tutto a sinistra di 1, così il vicino a destra finisce sopra lo spin centrale
    vicini_destra = np.roll(reticolo, shift=-1, axis=1) 
    
    #sposta tutto in alto di 1
    vicini_sotto = np.roll(reticolo, shift=-1, axis=0)
    
    H = -J*np.sum(reticolo*vicini_sotto + reticolo*vicini_destra)
    
    return H 



def delta_energia(reticolo, i, j, J, L):
    
    S_vicini = reticolo[(i+1)%L, j] + reticolo[(i-1)%L, j] + reticolo[i, (j+1)%L] + reticolo[i, (j-1)%L]

    return 2*J*reticolo[i,j] * S_vicini



def passo_metropolis(reticolo, J, T, L):
    
    for _ in range(L**2):
        i, j = np.random.randint(0, L, size=2)  #prendo randomicamente gli indici dello spin
        dE = delta_energia(reticolo, i, j, J, L) #calcolo l'energia
        
        if dE <= 0:
            reticolo[i, j] *= -1 #accetto: flip dello spin
        else:
            r = np.random.rand()
            if r < np.exp(-dE/T):
                reticolo[i, j] *= -1
                
    return reticolo
    
N_passi = 2000

E = np.zeros(N_passi)
M = np.zeros(N_passi)  #energia e magnetizzazione



L = 20
J = 1.0
T = 2.0



reticolo = inizializza_reticolo(L, tipo='random')

    
for step in range(N_passi):
    #passo di Metropolis per aggiornare il reticolo
    reticolo = passo_metropolis(reticolo, J, T, L)
    
    #Calcolo i valori correnti
    energia_corrente = calcola_energia(reticolo, J)
    mag_corrente = calcola_magnetizzazione(reticolo)
    
    E[step] = energia_corrente
    M[step] = mag_corrente


t = np.arange(0, N_passi)

plt.plot(t, E, label='E(t)')

plt.plot(t, M, label='M(t)')

plt.legend()
plt.xlabel('t')


plt.show()


#Esercizio 4



#reticolo = inizializza_reticolo(L, tipo='random')

T_values = [1.0, 1.5, 2.0, 2.2, 2.3, 2.4, 2.5, 3.0, 4.0]

N_passi = 2000

E = np.zeros(N_passi)
M = np.zeros(N_passi)


E_medie = []
M_medie = []
N_spin = L**2


for T in T_values:
    # A ogni nuova temperatura, ripartiamo da un reticolo nuovo
    reticolo = inizializza_reticolo(L, tipo="random")
    
    E_validi = []
    M_validi = []
    
    for step in range(N_passi):
        reticolo = passo_metropolis(reticolo, J, T, L)
        
        #Dopo 500 passi di termalizzazione, misura
        if step >= 500:
            E_validi.append(calcola_energia(reticolo, J))
            M_validi.append(calcola_magnetizzazione(reticolo))
            
    #Finiti i passi per questa temperatura, calcolo i valori medi per spin
    E_medie.append(np.mean(E_validi) / N_spin)
    M_medie.append(np.mean(np.abs(M_validi)) / N_spin)
    
    plt.imshow(reticolo.T, aspect='auto', origin='lower', cmap='inferno')
    plt.colorbar(label='Spin')
    plt.title(f'Reticolo a temperatura T = {T}')
    plt.show()




T_c = 2.269

temp_array = np.array(T_values)


fig, axs = plt.subplots(1, 2, figsize=(18, 8))


axs[0].plot(temp_array, E_medie, label='<E>/N')    
axs[0].set_xlabel('T')
axs[0].set_ylabel('<E>/N')
axs[0].axvline(x=T_c, color='r', linestyle='--')
axs[0].legend()

axs[1].plot(temp_array, M_medie, label='<|M|>/N')    
axs[1].set_xlabel('T')
axs[1].set_ylabel('<|M|>/N')
axs[1].axvline(x=T_c, color='r', linestyle='--')
axs[1].legend()


plt.show()

























