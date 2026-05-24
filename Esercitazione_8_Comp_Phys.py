#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 21:18:56 2026

@author: alberto
"""

import numpy as np
import matplotlib.pyplot as plt
#import scienceplots
from numba import njit
#plt.style.use(['science', 'notebook', 'grid', 'ieee'])   #migliora la visibilità#

def main():
    
    
    
    #Esercizio 1
    
    
    L = 40
    T_0 = 5.0
    T_f = 0.01
    N_passi = 2000
    N_spin = L**2
    
    J_v, J_h = genera_accoppiamenti(L)
    
    
    visualizza_accoppiamenti(J_v, r'Accoppiamenti per $J_v$')
    visualizza_accoppiamenti(J_h, r'Accoppiamenti per $J_h$')
    
    n_fer_J_v, n_antifer_J_v = conta_spin(J_v)
    n_fer_J_h, n_antifer_J_h = conta_spin(J_h)
    
    print(f"Numero di +1 per J_v: {n_fer_J_v}")
    print(f"Numero di -1 per J_v: {n_antifer_J_v}")
    print(f"Numero di +1 per J_h: {n_fer_J_h}")
    print(f"Numero di -1 per J_h: {n_antifer_J_h}")
    

    #Esercizio 2
    
    J_v, J_h = genera_accoppiamenti(L)
    E_finale_N_greedy = []
    M_finale_N_greedy = []
    
    for _ in range(10):
        
        reticolo = inizializza_reticolo(L, 'random')
        
    
        for step in range(N_passi):
            
            reticolo = passo_greedy(reticolo, J_h, J_v, L)
            
            energia_corrente = calcola_energia(reticolo, J_h, J_v)
            mag_corrente = calcola_magnetizzazione(reticolo)
            
        E_finale_N_greedy.append(energia_corrente/N_spin)
        M_finale_N_greedy.append(np.abs(mag_corrente)/N_spin)
            
    
    E_finale_N_greedy = np.array(E_finale_N_greedy)
    M_finale_N_greedy = np.array(M_finale_N_greedy)
    E_massimo = np.max(E_finale_N_greedy)
    E_minimo = np.min(E_finale_N_greedy)

    print(f"Valore medio dell'energia per reticoli random: {np.mean(E_finale_N_greedy)} +/- {np.std(E_finale_N_greedy)}")
    print(f"Valore medio della magnetizzazione per reticoli random: {np.mean(M_finale_N_greedy)} +/- {np.std(M_finale_N_greedy)}")
    print(f"Valore massimo dell'energia per reticoli random: {E_massimo}")  
    print(f"Valore minimo dell'energia per reticoli random: {E_minimo}")

    #Esercizio 3


    J_v, J_h = genera_accoppiamenti(L)
    E_finale_N_metro = []
    M_finale_N_metro = []
    
    for _ in range(10):
        
        reticolo = inizializza_reticolo(L, 'random')
        
    
        for step in range(N_passi):
            
            reticolo = passo_metropolis(reticolo, J_h, J_v, L, T_0)
            
            energia_corrente = calcola_energia(reticolo, J_h, J_v)
            mag_corrente = calcola_magnetizzazione(reticolo)
            
        E_finale_N_metro.append(energia_corrente/N_spin)
        M_finale_N_metro.append(np.abs(mag_corrente)/N_spin)


    E_finale_N_metro = np.array(E_finale_N_metro)
    M_finale_N_metro = np.array(M_finale_N_metro)
    E_massimo = np.max(E_finale_N_metro)
    E_minimo = np.min(E_finale_N_metro)

    print(f"Valore medio dell'energia per reticoli random: {np.mean(E_finale_N_metro)} +/- {np.std(E_finale_N_metro)}")
    print(f"Valore medio della magnetizzazione per reticoli random: {np.mean(M_finale_N_metro)} +/- {np.std(M_finale_N_metro)}")
    print(f"Valore massimo dell'energia per reticoli random: {E_massimo}")  
    print(f"Valore minimo dell'energia per reticoli random: {E_minimo}")

      
    #Esercizio 4

    t_max = N_passi
    T_storia_esponenziale = np.zeros(N_passi)
    E_storia_esponenziale = np.zeros(N_passi)
    reticolo = inizializza_reticolo(L, 'random')
    
    #schedule esponenziale

    for step in range(N_passi):
        
        T_t_esponenziale = schedule_esponenziale(step, T_0, T_f, t_max)
        reticolo = passo_metropolis(reticolo, J_h, J_v, L, T_t_esponenziale)
        energia_corrente = calcola_energia(reticolo, J_h, J_v)
        E_storia_esponenziale[step] = energia_corrente/N_spin
        T_storia_esponenziale[step] = T_t_esponenziale


    T_storia_lineare = np.zeros(N_passi)
    E_storia_lineare = np.zeros(N_passi)
    reticolo = inizializza_reticolo(L, 'random')
    
    #schedule lineare

    for step in range(N_passi):
        
        T_t_lineare = schedule_lineare(step, T_0, T_f, t_max)
        reticolo = passo_metropolis(reticolo, J_h, J_v, L, T_t_lineare)
        energia_corrente = calcola_energia(reticolo, J_h, J_v)
        E_storia_lineare[step] = energia_corrente/N_spin
        T_storia_lineare[step] = T_t_lineare

    passi = np.arange(N_passi)


    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].plot(passi, E_storia_esponenziale, label='Schedule Esponenziale')
    axes[0].plot(passi, E_storia_lineare, linestyle='--', label='Schedule Lineare')
    axes[0].set_title('Evoluzione dell\'energia')
    axes[0].set_xlabel('Passi')
    axes[0].set_ylabel('Energia per spin')  
    axes[0].grid()
    axes[0].legend()

    axes[1].plot(passi, T_storia_esponenziale, linestyle='--', label='Schedule Esponenziale')
    axes[1].plot(passi, T_storia_lineare, label='Schedule Lineare')
    axes[1].set_title('Evoluzione della temperatura')
    axes[1].set_xlabel('Passi')
    axes[1].set_ylabel('Temperatura')
    axes[1].set_yscale('log')
    axes[1].grid()
    axes[1].legend()
    plt.tight_layout()
    plt.show()


    E_finale_N_exp = []

    for _ in range(10):
        
        reticolo = inizializza_reticolo(L, 'random')
        
    
        for step in range(N_passi):
            
            T_t_esponenziale = schedule_esponenziale(step, T_0, T_f, t_max)
            reticolo = passo_metropolis(reticolo, J_h, J_v, L, T_t_esponenziale)
            
            energia_corrente = calcola_energia(reticolo, J_h, J_v)
            mag_corrente = calcola_magnetizzazione(reticolo)
            
        E_finale_N_exp.append(energia_corrente/N_spin)
        
    E_finale_N_exp = np.array(E_finale_N_exp)
    E_massimo_exp = np.max(E_finale_N_exp)
    E_minimo_exp = np.min(E_finale_N_exp)

    print(f"Valore medio dell'energia per schedule esponenziale: {np.mean(E_finale_N_exp)} +/- {np.std(E_finale_N_exp)}")
    print(f"Valore massimo dell'energia per schedule esponenziale: {E_massimo_exp}")  
    print(f"Valore minimo dell'energia per schedule esponenziale: {E_minimo_exp}")


    E_finale_N_lineare = []

    for _ in range(10):
        
        reticolo = inizializza_reticolo(L, 'random')
        
    
        for step in range(N_passi):
            
            T_t_lineare = schedule_lineare(step, T_0, T_f, t_max)
            reticolo = passo_metropolis(reticolo, J_h, J_v, L, T_t_lineare)
            
            energia_corrente = calcola_energia(reticolo, J_h, J_v)
            
        E_finale_N_lineare.append(energia_corrente/N_spin)


    E_finale_N_lineare = np.array(E_finale_N_lineare)
    E_massimo_lineare = np.max(E_finale_N_lineare)
    E_minimo_lineare = np.min(E_finale_N_lineare)

    print(f"Valore medio dell'energia per schedule lineare: {np.mean(E_finale_N_lineare)} +/- {np.std(E_finale_N_lineare)}")
    print(f"Valore massimo dell'energia per schedule lineare: {E_massimo_lineare}")  
    print(f"Valore minimo dell'energia per schedule lineare: {E_minimo_lineare}")

    if np.mean(E_finale_N_exp) < np.mean(E_finale_N_lineare):
        print("La schedule esponenziale ha un'energia media più bassa rispetto alla schedule lineare.")
    elif np.mean(E_finale_N_exp) > np.mean(E_finale_N_lineare):
        print("La schedule lineare ha un'energia media più bassa rispetto alla schedule esponenziale.")
    else:
        print("Entrambe le schedule hanno  la stessa energia media.")


    #Esercizio 5

    media_greedy = np.mean(E_finale_N_greedy)
    std_greedy = np.std(E_finale_N_greedy)
    media_metro = np.mean(E_finale_N_metro)
    std_metro = np.std(E_finale_N_metro)
    media_lin = np.mean(E_finale_N_lineare)
    std_lin = np.std(E_finale_N_lineare)
    media_exp = np.mean(E_finale_N_exp) 
    std_exp = np.std(E_finale_N_exp)


    medie = [media_greedy, media_metro, media_lin, media_exp]
    errori = [std_greedy, std_metro, std_lin, std_exp]
    metodi = ['Greedy', 'Metropolis', 'Lineare', 'Esponenziale']

    plt.bar(metodi, medie, width=0.5, yerr=errori, capsize=5)
    plt.ylabel('Energia per spin')
    plt.title('Confronto tra i metodi di ottimizzazione')
    plt.show()

    print("Confronto tra i metodi di ottimizzazione:")
    print(f"Greedy: {media_greedy:.4f} +/- {std_greedy:.4f}")
    print(f"Metropolis: {media_metro:.4f} +/- {std_metro:.4f}")
    print(f"Schedule Lineare: {media_lin:.4f} +/- {std_lin:.4f}")
    print(f"Schedule Esponenziale: {media_exp:.4f} +/- {std_exp:.4f}")


    passi_list = [200, 500, 1000, 2000]

    medie_passi = []
    std_passi = []

    for N_passi in passi_list:
        
        E_annealing = []
        
        for _ in range(4):
            
            reticolo = inizializza_reticolo(L, 'random')
            
        
            for step in range(N_passi):
                
                T_t_esponenziale = schedule_esponenziale(step, T_0, T_f, N_passi)
                reticolo = passo_metropolis(reticolo, J_h, J_v, L, T_t_esponenziale)
                energia_corrente = calcola_energia(reticolo, J_h, J_v)
                
                    
            E_annealing.append(energia_corrente/N_spin)
        
        E_annealing = np.array(E_annealing)
        medie_passi.append(np.mean(E_annealing))
        std_passi.append(np.std(E_annealing))



    plt.errorbar(passi_list, medie_passi, yerr=std_passi, fmt='o-', capsize=5)
    plt.xscale('log')
    plt.xlabel('Numero di passi (log scale)')
    plt.ylabel('Energia per spin')
    plt.grid()
    plt.title('Effetto del numero di passi sull\'energia finale (Schedule Esponenziale)')
    plt.show()

    reticolo_partenza = inizializza_reticolo(L, 'random')

    #faccio delle copie del reticolo di partenza per ogni metodo, in modo da confrontarli a parità di condizioni iniziali
    
    reticolo_greedy = np.copy(reticolo_partenza)
    for _ in range(N_passi):
        reticolo_greedy = passo_greedy(reticolo_greedy, J_h, J_v, L)

    
    reticolo_metro = np.copy(reticolo_partenza)
    for _ in range(N_passi):
        reticolo_metro = passo_metropolis(reticolo_metro, J_h, J_v, L, 0.5)

    
    reticolo_annealing = np.copy(reticolo_partenza)
    for step in range(N_passi):
        T_t = schedule_esponenziale(step, T_0, T_f, N_passi)
        reticolo_annealing = passo_metropolis(reticolo_annealing, J_h, J_v, L, T_t)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    #Reticolo Greedy
    im = axes[0].imshow(reticolo_greedy, cmap='coolwarm')
    axes[0].set_title('Greedy')
    axes[0].axis('off') 

    #Reticolo Metropolis
    axes[1].imshow(reticolo_metro, cmap='coolwarm')
    axes[1].set_title('Metropolis (T=0.5)')
    axes[1].axis('off')

    #Reticolo Simulated Annealing
    axes[2].imshow(reticolo_annealing, cmap='coolwarm')
    axes[2].set_title('Annealing Esponenziale')
    axes[2].axis('off')

    #Colorbar collegata a 'im' 
    fig.colorbar(im, ax=axes, label='Spin (+1 / -1)', fraction=0.02, pad=0.04)

    plt.show()

    #Calcolo delle magnetizzazioni finali per spin |M|/N
    mag_greedy = np.abs(calcola_magnetizzazione(reticolo_greedy)) / N_spin
    mag_metro = np.abs(calcola_magnetizzazione(reticolo_metro)) / N_spin
    mag_annealing = np.abs(calcola_magnetizzazione(reticolo_annealing)) / N_spin

    print("Magnetizzazione finale per spin |M|/N:")
    print(f"- Greedy: {mag_greedy:.4f}")
    print(f"- Metropolis a T=0.5: {mag_metro:.4f}")
    print(f"- Annealing Esponenziale: {mag_annealing:.4f}")




















def genera_accoppiamenti(L):

    J_v = np.random.choice([1, -1], size=(L, L))
    J_h = np.random.choice([1, -1], size=(L, L))

    return J_v, J_h


def conta_spin(J):

    n_fer = np.sum(J==1)  #True == 1, False == 0 in python
    n_antifer = np.sum(J != 1)
    
    return n_fer, n_antifer

def visualizza_accoppiamenti(J, titolo):

    plt.imshow(J.T, aspect='auto', origin='lower', cmap='coolwarm')
    plt.colorbar(label='Interazione')
    plt.title(titolo)
    plt.show()

@njit   #non uso np.roll perché numba mi velocizza, ma vuole i cicli for
def calcola_energia(reticolo, J_h, J_v):
    L = reticolo.shape[0]
    H = 0.0
    
    for i in range(L):
        for j in range(L):
            # Moltiplico lo spin (i,j) per il vicino sotto e il vicino a destra
            H += -(reticolo[i,j] * reticolo[(i+1)%L, j] * J_v[i,j] + 
                   reticolo[i,j] * reticolo[i, (j+1)%L] * J_h[i,j])
            
    return H


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


@njit
def delta_energia(reticolo, i, j, J_h, J_v, L):

    delta_E = 2*reticolo[i,j] * (J_h[i,j]*reticolo[i, (j+1)%L] + 
                                 J_h[i,(j-1)%L]*reticolo[i, (j-1)%L] + 
                                 J_v[i,j]*reticolo[(i+1)%L,j] + 
                                 J_v[(i-1)%L, j]*reticolo[(i-1)%L, j])
    
    return delta_E
    
@njit
def passo_greedy(reticolo, J_h, J_v, L):
    
    for _ in range(L**2):
        i, j = np.random.randint(0, L, size=2)  #prendo randomicamente gli indici dello spin
        dE = delta_energia(reticolo, i, j, J_h, J_v, L) #calcolo l'energia
        
        if dE <= 0:
            reticolo[i, j] *= -1 #accetto: flip dello spin
         
    return reticolo

def calcola_magnetizzazione(reticolo):
    
    return np.sum(reticolo)

@njit
def passo_metropolis(reticolo, J_h, J_v, L, T):
    
    for _ in range(L**2):
        i, j = np.random.randint(0, L, size=2)
        dE = delta_energia(reticolo, i, j, J_h, J_v, L)
        
        if dE <= 0 or np.random.rand() < np.exp(-dE/T):
            reticolo[i, j] *= -1
            
    return reticolo

def schedule_lineare(t, T_0, T_f, t_max):

    T_t = T_0 - (T_0 - T_f)*(t/t_max)

    return T_t

def schedule_esponenziale(t, T_0, T_f, t_max):

    T_t = T_0*(T_f/T_0)**(t/t_max)

    return T_t

main()






















































































































