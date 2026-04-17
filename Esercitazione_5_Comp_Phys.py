#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 12:51:52 2026

@author: alberto
"""

import numpy as np
import matplotlib.pyplot as plt




def main():

    


#Esercizio 1


    N = 50
    L = 1.0
    h = L/N
    rho_0 = 1.0
    pi = np.pi
    #epsilon_0 = 8.85418782 * 10**(-12) #m^-3 kg^-1 s^4 A^2
    epsilon_0 = 1.0 #in unità naturali!
    
    x = np.linspace(0, L, N+1)  # 51 punti, includo x=L
    
    rho = rho_0 * np.sin((pi*x)/L)
    
    V_a, V_b = 0.0, 0.0

    mat_A = A(N, h) #matrice A
    vec_b = b(rho, epsilon_0, V_a, V_b, h) #termine noto


    

    
    #Risoluzione
    
    phi_analitica = (rho_0/epsilon_0)*((L/pi)**2)*np.sin((pi*x)/L)

    phi_internal = np.linalg.solve(mat_A, vec_b)

    phi_completa = np.concatenate([[V_a], phi_internal, [V_b]])  #51 elementi
    
    errore_max = np.max(np.abs(phi_completa - phi_analitica))
    
    errore_rel = errore_max/np.max(np.abs(phi_analitica)) * 100  
    
    E_analitico = -(rho_0/epsilon_0)*(L/pi)*np.cos((pi*x)/L)

    E_num = -derivata_prima_f(phi_completa, x, h)
    
    
    
    #figura con 1 riga e 2 colonne
    fig, axs = plt.subplots(1, 2, figsize=(10, 6))
    
    #Potenziale
    axs[0].plot(x, phi_completa, '-', color = 'purple', label=r'$\phi_{num}$(x)')
    axs[0].plot(x, phi_analitica, '--', color = 'orange', label=r'$\phi_{analitica}$(x)')
    axs[0].set_title('Potenziale Elettrostatico')
    axs[0].set_xlabel('x [m]')
    axs[0].set_ylabel(r'$\phi$(x) [V]')
    axs[0].legend()
    axs[0].grid()

    #Campo Elettrico
    axs[1].plot(x, E_num, '-', color='purple', label=r'$E_{num}$(x)')
    axs[1].plot(x, E_analitico, '--', color='orange', label=r'$E_{analitico}$(x)')
    axs[1].set_title('Campo Elettrico')
    axs[1].set_xlabel('x [m]')
    axs[1].set_ylabel('E(x) [V/m]')
    axs[1].legend()
    axs[1].grid()

    fig.suptitle(r'Distribuzione di carica sinusoidale: $\rho (x) = \rho_0 sin(\pi x/L)$', fontsize=16, fontweight='bold')

    plt.tight_layout() #evita che le etichette si sovrappongano
    plt.show()
    
    print(f"Errore massimo: {errore_max}") #3.3339914100094514e-05
    print(f'epsilon_rel = {errore_rel} %') #0.03290517629342339 %



#Esercizio 2


    N = [10, 20, 40, 80, 160]
    errori = []
    h_values = []
    
    
    for n in N:
        
        h = L/n
        h_values.append(h)
        
        x = np.linspace(0, L, n+1)
        
        rho = rho_0 * np.sin((pi*x)/L)
        
        mat_A = A(n, h) #matrice A
        vec_b = b(rho, epsilon_0, V_a, V_b, h) #termine noto
        
        phi_analitica = (rho_0/epsilon_0)*((L/pi)**2)*np.sin((pi*x)/L)

        phi_internal = np.linalg.solve(mat_A, vec_b)

        phi_completa = np.concatenate([[V_a], phi_internal, [V_b]])  #51 elementi
        
        errore_max = np.max(np.abs(phi_completa - phi_analitica))
        
        errori.append(errore_max)
        
    
    
    
    
    #conversione in array
    
    h_array = np.array(h_values)
    
    #retta di riferimento
    linea_riferimento =  (h_array**2)  #in scala log-log, mi darà una retta  con pendenza 2
    
    plt.plot(h_values, errori, 'o-', color = 'purple', label='Errore Massimo')
    plt.plot(h_values, linea_riferimento, '--', color = 'orange',label=r'Riferimento $h^2$')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    
    plt.legend()
    plt.title('Errore al variare di h', fontweight='bold')
    plt.show()
        
   
    
#Esercizio 3
    
    N = 50
    L = 1.0
    h = L/N    
    
    x = np.linspace(0, L, N+1)  # 51 punti, includo x=L
    
    rho = rho_0*np.ones(len(x))
    
    V_a, V_b = 0.0, 0.0
    
    mat_A = A(N, h) #matrice A
    vec_b = b(rho, epsilon_0, V_a, V_b, h) #termine noto
    
    
    
    
    
    #Risoluzione
    
    phi_analitica = (rho_0/(2*epsilon_0))*x*(L-x)
    
    phi_internal = np.linalg.solve(mat_A, vec_b)
    
    phi_completa = np.concatenate([[V_a], phi_internal, [V_b]])  #51 elementi
    
    errore_max = np.max(np.abs(phi_completa - phi_analitica))
    
    errore_rel = errore_max/np.max(np.abs(phi_analitica)) * 100 

    E_analitico = -(rho_0 / (2 * epsilon_0)) * (L - 2*x)
    
    E_num = -derivata_prima_f(phi_completa, x, h)


    #figura con 1 riga e 2 colonne
    fig, axs = plt.subplots(1, 2, figsize=(10, 6))
    
    #Potenziale
    axs[0].plot(x, phi_completa, '-', color = 'purple', label=r'$\phi_{num}$(x)')
    axs[0].plot(x, phi_analitica, '--', color = 'orange', label=r'$\phi_{analitica}$(x)')
    axs[0].set_title('Potenziale Elettrostatico')
    axs[0].set_xlabel('x [m]')
    axs[0].set_ylabel(r'$\phi$(x) [V]')
    axs[0].legend()
    axs[0].grid()

    #Campo Elettrico
    axs[1].plot(x, E_num, '-', color='purple', label=r'$E_{num}$(x)')
    axs[1].plot(x, E_analitico, '--', color='orange', label=r'$E_{analitico}$(x)')
    axs[1].set_title('Campo Elettrico')
    axs[1].set_xlabel('x [m]')
    axs[1].set_ylabel('E(x) [V/m]')
    axs[1].legend()
    axs[1].grid()

    fig.suptitle(r'Distribuzione di carica costante: $\rho(x) = \rho_0$ ', fontsize=16, fontweight='bold')

    plt.tight_layout() #evita che le etichette si sovrappongano
    plt.show()
    

    print(f'epsilon_rel = {errore_rel} %') #5.10702591327572e-13 %

    #max potenziale per x = L/2 => phi_max = (rho_0*L**2)/(8*epsilon_0)
    
    phi_max_th = (rho_0*L**2)/(8*epsilon_0)
    
    phi_max_num = np.max(phi_completa)
    
    err_rel_phi_max = (np.abs(phi_max_num - phi_max_th)/phi_max_th)*100
    
    print(f'epsilon_phi_max = {err_rel_phi_max} %') #4.218847493575595e-13 %
    
    
    
#Esercizio 4



    N = 50
    L = 1.0
    h = L/N    
    
    x = np.linspace(0, L, N+1)  # 51 punti, includo x=L
    
    rho = np.zeros(len(x))
    
    V_a, V_b = 0.0, 10.0
    
    mat_A = A(N, h) #matrice A
    vec_b = b(rho, epsilon_0, V_a, V_b, h) #termine noto
    
    
    
    
    
    #Risoluzione
    
    phi_analitica = V_a + ((V_b - V_a)*x)/L
    
    phi_internal = np.linalg.solve(mat_A, vec_b)
    
    phi_completa = np.concatenate([[V_a], phi_internal, [V_b]])  #51 elementi
    
    errore_max = np.max(np.abs(phi_completa - phi_analitica))
    
    errore_rel = errore_max/np.max(np.abs(phi_analitica)) * 100 

    E_analitico = (-(V_b - V_a)/L)*np.ones(len(x))
    
    E_num = -derivata_prima_f(phi_completa, x, h)

    
    #figura con 1 riga e 2 colonne
    fig, axs = plt.subplots(1, 2, figsize=(10, 6))
    
    #Potenziale
    axs[0].plot(x, phi_completa, '-', color = 'purple', label=r'$\phi_{num}$(x)')
    axs[0].plot(x, phi_analitica, '--', color = 'orange', label=r'$\phi_{analitica}$(x)')
    axs[0].set_title('Potenziale Elettrostatico')
    axs[0].set_xlabel('x [m]')
    axs[0].set_ylabel(r'$\phi$(x) [V]')
    axs[0].legend()
    axs[0].grid()

    #Campo Elettrico
    axs[1].plot(x, E_num, '-', color='purple', label=r'$E_{num}$(x)')
    axs[1].plot(x, E_analitico, '--', color='orange', label=r'$E_{analitico}$(x)')
    axs[1].set_title('Campo Elettrico')
    axs[1].set_xlabel('x [m]')
    axs[1].set_ylim(-15.0, -5.0)
    axs[1].set_ylabel('E(x) [V/m]')
    axs[1].legend()
    axs[1].grid()

    fig.suptitle('Condensatore a facce piane parallele', fontsize=16, fontweight='bold')

    plt.tight_layout() #evita che le etichette si sovrappongano
    plt.show()
    
    
    
    print(f'epsilon_rel = {errore_rel} %') #4.352074256530614e-13

#Esercizio 5



    N = 100
    L = 1.0
    h = L/N    
    
    x = np.linspace(0, L, N+1)  # 51 punti, includo x=L
    
    sigma = 0.1 * L
    
    rho = rho_0*np.exp((-(x-L/2)**2)/(2*sigma**2))
    
    V_a, V_b = 0.0, 0.0
    
    mat_A = A(N, h) #matrice A
    vec_b = b(rho, epsilon_0, V_a, V_b, h) #termine noto
    
    
    
    
    
    #Risoluzione
    
    
    phi_internal = np.linalg.solve(mat_A, vec_b)
    
    phi_completa = np.concatenate([[V_a], phi_internal, [V_b]])  #51 elementi
    
    E_num = -derivata_prima_f(phi_completa, x, h)

    phi_max = np.max(phi_completa)    
    
    E_abs = np.abs(E_num)          
      
    E_max_abs = np.max(E_abs)     
        
    indice_E_max = np.argmax(E_abs)
    
    x_E_max = x[indice_E_max]
    
    indice_phi_max = np.argmax(phi_completa)
    
    x_phi_max = x[indice_phi_max]
    
    
    print(f"E_max = {E_max_abs:.4f} V/m in x = {x_E_max:.4f} m") #E_max = 0.1253 V/m in x = 1.0000 m
    print(f"phi_max = {phi_max:.4f} V in x = {x_phi_max:.4f} m") #phi_max = 0.0527 V in x = 0.5000 m
    
    
    
    
    #E nullo in x = L/2

    #figura con 1 riga e 3 colonne
    fig, axs = plt.subplots(1, 3, figsize=(15, 7))

    #Densità di carica
    axs[0].plot(x, rho, '-')
    axs[0].set_title('Densità di carica')
    axs[0].set_xlabel('x [m]')
    axs[0].set_ylabel(r'$\rho$(x) [C/m]')
    axs[0].grid()

    #Potenziale
    axs[1].plot(x, phi_completa, '-', color='purple')
    axs[1].set_title('Potenziale Elettrostatico')
    axs[1].set_xlabel('x [m]')
    axs[1].set_ylabel(r'$\phi$(x) [V]')
    axs[1].grid()

    #Campo Elettrico
    axs[2].plot(x, E_num, '-', color='orange')
    axs[2].set_title('Campo Elettrico')
    axs[2].set_xlabel('x [m]')
    axs[2].set_ylabel('E(x) [V/m]')
    axs[2].grid()

    fig.suptitle('Distribuzione di carica gaussiana', fontsize=16, fontweight='bold')

    plt.tight_layout() #evita che le etichette si sovrappongano
    plt.show()







# matrice A

def A(N, h):

    # Diagonale principale: 2/h^2
    
    diag_main = np.ones(N-1) * 2.0 / h**2
    
    # Diagonali laterali: -1/h^2
    
    diag_off = np.ones(N-2) * (-1.0) / h**2
    
    # Matrice tridiagonale
    
    A = np.diag(diag_main) + np.diag(diag_off, k=1) + np.diag(diag_off, k=-1)

    return A

# Termine noto

def b(rho, epsilon_0, V_a, V_b, h):

    b = rho[1:-1] / epsilon_0 #taglio i valori al bordo di rho, lunghezza N-1 = 49
    
    b[0] += V_a / h**2 # condizione al contorno sinistra
    
    b[-1] += V_b / h**2 # condizione al contorno destra
    
    return b



def derivata_prima_f(f, x, h):
    df_dx = np.zeros(len(x))
    
    #parto da 1 e mi fermo al penultimo elemento
    for i in range(1, len(x)-1):
        
        df_dx[i] = (f[i+1] - f[i-1]) / (2*h)
        
    df_dx[0] = (f[1] - f[0]) / h  #bordo sx
        
    df_dx[-1] = (f[-1] - f[-2]) / h   #bordo dx
        
    return df_dx
    
    



main()








