#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 15:53:19 2026

@author: alberto
"""

import numpy as np
import matplotlib.pyplot as plt


def main():
    
    #Esercizio 1
    
    
    L = 1  #m
    alpha = 0.01  #m^2/s
    N_x = 50
    dx = L/N_x
    pi = np.pi
    dt_max = (dx**2)/(2*alpha)
    
    dt = 0.4 * dt_max

    
    T_0 = 1.0 #°C
    
    t_max = 5.0 #s
    
    x = np.linspace(0, L, N_x + 1)  #len 51
    t = np.arange(0, t_max + dt/2, dt)
    
    tt, xx = np.meshgrid(t,x) #meshgrid crea len(t) colonne per le t e len(x) righe per le x
    
    
    #condizione iniziale
    
    T_x_0 = T_0 * np.sin(pi * x / L)
    
    T_left = 0.0
    T_right = 0.0
    
    T_analitica = T_0 * np.exp(-alpha * tt * (pi / L)**2) * np.sin(pi * xx / L)

    T_x_t = FTCS(T_x_0, x, t, dx, dt, alpha, T_left, T_right)
        
    times = [0.0, 0.5, 1.0, 2.0, 5.0]
    
    fig, axs = plt.subplots(1, 2, figsize=(18, 8))
    
    for t_target in times:
        
        n = int(t_target / dt) #recupero l'indice
        
        axs[0].plot(x, T_x_t[:, n], label=f't = {t_target} s')
        
       
        
    T_analitica_finale = T_analitica[:, -1]
    epsilon = np.max(np.abs(T_analitica_finale - T_x_t[:,-1]))
    
    
    
    
    axs[0].plot(x, T_analitica_finale, 'r.', linewidth = 0.6, label='Sol analitica a t = 5.0s')    
    axs[0].set_xlabel('x [m]')
    axs[0].set_ylabel('T(x) [°C]')
    axs[0].legend()
    axs[0].grid(True)
        
    axs[1].imshow(T_x_t.T, aspect='auto', origin='lower', cmap='inferno')
    axs[1].set_xlabel('x [m]')
    axs[1].set_ylabel('T(x) [°C]')
    
    fig.suptitle(r'Profilo iniziale T(x,0) = $T_0$ sin($\pi$x/L) ', fontsize=18, fontweight='bold')
    
    
    
    
    print(f'Errore massimo a t_max = {t_max}s: {epsilon}')
    #Errore massimo a t_max = 5.0s: 1.9819753150351893e-05

    

    #Esercizio 2
    
    #cond iniz
    
    sigma = 0.1 * L
    T_x_0 = T_0 * np.exp(-(x - L/2)**2 / (2*sigma**2)) #gaussiana

    
    T_left = 0.0
    T_right = 0.0

    t_max = 10.0 #s

    t = np.arange(0, t_max + dt/2, dt)
    
    t_values = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]
    
    T_x_t = FTCS(T_x_0, x, t, dx, dt, alpha, T_left, T_right)
    
    
    fig, axs = plt.subplots(1, 3, figsize=(20, 8))
    
    for time in t_values:
        
        n = int(time/dt) #indice
        
        
        
        axs[0].plot(x, T_x_t[:, n], label=f't = {time} s')
        axs[0].set_xlabel('x [m]')
        axs[0].set_ylabel('T(x) [°C]')
        axs[0].legend()
        axs[0].grid(True)
        
        
        
    

    T_max_t = np.max(T_x_t, axis = 0) #cerca i massimi sull'asse 0, ovvero su quello dello spazio
    
    
    axs[1].plot(t, T_max_t, color = 'purple')
    axs[1].set_xlabel('t [s]')
    axs[1].set_ylabel(r'$T_{max}$(t) [°C]')
    axs[1].grid()
        
    
    
    axs[2].imshow(T_x_t.T, aspect='auto', origin='lower', cmap='inferno')
    axs[2].set_xlabel('x [m]')
    axs[2].set_ylabel('T(x) [°C]')
    
    fig.suptitle('Distribuzione iniziale gaussiana', fontsize=18, fontweight='bold')
    
    
    plt.tight_layout() #evita che le etichette si sovrappongano
    plt.show()
    
    #Esercizio 3
    
    #cond iniz
    
    sigma = 0.03 * L
    T_x_0 = T_0 * np.exp(-(x - L/2)**2 / (2*sigma**2)) #gaussiana stretta

    
    T_left = 0.0
    T_right = 0.0

    t_max = 0.5 #s

    r_values = [0.4, 0.5, 0.6, 0.8, 1.0, 1.2]
    
    
    fig, axs = plt.subplots(2, 3, figsize=(10, 6))
    axs_flat = axs.flatten()
    
    for i, r_value in enumerate(r_values):
        
        dt = (r_value * dx**2) / alpha
        
        t = np.arange(0, t_max + dt/2, dt)
        
        T_x_t = FTCS(T_x_0, x, t, dx, dt, alpha, T_left, T_right)
        
        T_x_t_finale = T_x_t[:,-1]
        
        
        axs_flat[i].plot(x, T_x_t_finale, label=f'r = {r_value}', color = 'purple')
        axs_flat[i].set_xlabel('x [m]')
        axs_flat[i].set_ylabel('T(x) [°C]')
        axs_flat[i].legend()
        axs_flat[i].grid(True)
   
        #da r=0.6 iniziano le oscillazioni; da r= 0.8 il profilo gaussiano si perde completamente
        
    fig.suptitle(r'Distribuzione iniziale gaussiana $\sigma = 0.03$ al variare di r', fontsize=18, fontweight='bold')


    #Esercizio 4
    
    
    
    T_x_0 = 50 #°C
    
    T_left = 100
    T_right = 20
    
    t_max = 50.0
    
    dt = 0.4 * dt_max
    
    t = np.arange(0, t_max + dt/2, dt)
    
    T_staz = 100.0 - 80.0*(x/L)
    
    T_x_t = FTCS(T_x_0, x, t, dx, dt, alpha, T_left, T_right)
    
    times = [0.0, 1.0, 5.0, 10.0, 20.0, 50.0]
    
    
    fig, axs = plt.subplots(1, 3, figsize=(20, 8))
    
    for time in times:
        
        n = int(time / dt) #recupero l'indice
        
        axs[0].plot(x, T_x_t[:, n], label=f't = {time} s')
        axs[0].set_xlabel('x [m]')
        axs[0].set_ylabel('T(x) [°C]')
        axs[0].legend()
        axs[0].grid(True)
    
    
    
    d_t = np.max(np.abs(T_staz[:, None] - T_x_t), axis=0)
    
    axs[1].plot(t, d_t, color = 'purple')
    axs[1].set_xlabel('t [s]')
    axs[1].set_ylabel(r'd(t) [°C]')
    axs[1].set_yscale('log')
    axs[1].grid()
    
    axs[2].imshow(T_x_t.T, aspect='auto', origin='lower', cmap='inferno')
    axs[2].set_xlabel('x [m]')
    axs[2].set_ylabel('T(x) [°C]')
    
    fig.suptitle('Profilo iniziale costante', fontsize=16, fontweight='bold')
    
    plt.tight_layout() #evita che le etichette si sovrappongano
    plt.show()
    
    
    
    
    


def FTCS(f_x_0, x, t, dx, dt, alpha, f_left, f_right):
    
    r = alpha*(dt / dx**2)
    
    f_x_t = np.zeros((len(x), len(t))) #inizializzo a zero (ho len(x) righe e len(t) colonne)
    
    f_x_t[:, 0] = f_x_0  #condizione a t_0 = 0; per ogni valore di riga, la colonna 0 è inizializzata a f_x_0
    
    #condizioni al contorno
    
    f_x_t[0, :] = f_left  #per ogni valore di colonna, le righe 0 e -1 sono inizializzate dai valori al contorno
    f_x_t[-1, :] = f_right
    
    for n in range(0, len(t) - 1):
        
        for i in range(1, len(x) - 1):  #lascio i limiti (bordi dx e sx)
            
            f_x_t[i, n+1] = f_x_t[i, n] + r * (f_x_t[i+1, n] - 2*f_x_t[i, n] + f_x_t[i-1, n])
            
            
    return f_x_t




main()



