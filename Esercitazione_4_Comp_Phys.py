#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 17:10:05 2026

@author: alberto
"""

import numpy as np
from matplotlib import pyplot as plt



def main():
    
    
    #Esercizio 1
    
    
    h = 0.001
    t_max = 20.0
    
    t = np.arange(0, t_max, h)
    
    
    Y_sol = np.zeros((len(t), 4))  #matrice di zeri per theta_1, theta_2, omega_1, omega_2
    
    Y_sol[0] = [0.1, 0.1, 0.0, 0.0]  # valori iniziali per theta_i, omega_i
    
    for i in range(len(t) - 1):
        
        Y_sol[i+1] = RK4(dY_dt, t[i], Y_sol[i], h) #la funzione calcola il passo basandosi 
        #esclusivamente sulle condizioni dell'istante corrente
        
        
    
    theta_1 = Y_sol[:, 0] #slicing: prendo tutti gli elementi (:) della riga 0
    theta_2 = Y_sol[:, 1] #slicing: prendo tutti gli elementi (:) della riga 1

    

    plt.plot(t, theta_1, 'r-')
    plt.xlabel('t [s]')
    plt.ylabel(r'$\theta_1$(t) [rad]')
    plt.title(r'$\theta_1$ vs t')
    plt.grid()
    
    plt.show()


    plt.plot(t, theta_2, 'b-')
    plt.xlabel('t [s]')
    plt.ylabel(r'$\theta_2$(t) [rad]')
    plt.title(r'$\theta_2$ vs t')
    plt.grid()
    
    plt.show()

    #analisi di fourier per le ampiezze dei theta(t)  
    #rivela le frequenze di fluttuazione numerica


    
    freq_positive, ampiezze_theta_1, nu_max, A_max = analisi_fourier(theta_1, h, label=r'$\theta_1$')
    
    freq_positive, ampiezze_theta_2, nu_max, A_max = analisi_fourier(theta_2, h, label=r'$\theta_2$')

    
    plt.plot(freq_positive, ampiezze_theta_1, color='red')
    plt.xlabel(r'$\nu$ $[s^{-1}]$')
    plt.ylabel(r'Ampiezza $\theta_1$')
    plt.xscale('log')
    plt.title(r'Ampiezza $\theta_1$ vs $\nu$')
    plt.grid()
    plt.show()
    
    plt.plot(freq_positive, ampiezze_theta_2, color='blue')
    plt.xlabel(r'$\nu$ $[s^{-1}]$')
    plt.ylabel(r'Ampiezza $\theta_2$')
    plt.xscale('log')
    plt.title(r'Ampiezza $\theta_2$ vs $\nu$')
    plt.grid()
    plt.show()
    
    #Massimo di ampiezza per theta_1: 679.8080014775281 rad a frequenza nu_1 = 0.4 Hz
    #Massimo di ampiezza per theta_2: 967.0443980801607 rad a frequenza nu_2 = 0.4 Hz
    
    
    
    
    
    
    
    


    #Esercizio 2
    
    
    Y_sol = np.zeros((len(t), 4))  #matrice di zeri per theta_1, theta_2, omega_1, omega_2
    
    Y_sol[0] = [np.pi/2, np.pi/2, 0.0, 0.0]  # valori iniziali per theta_i, omega_i
    
    for i in range(len(t) - 1):
        
        Y_sol[i+1] = RK4(dY_dt, t[i], Y_sol[i], h)



    theta_1 = Y_sol[:, 0] #slicing: prendo tutti gli elementi (:) della riga 0
    theta_2 = Y_sol[:, 1] #slicing: prendo tutti gli elementi (:) della riga 1

    

    plt.plot(t, theta_1, 'r-')
    plt.xlabel('t [s]')
    plt.ylabel(r'$\theta_1$(t) [rad]')
    plt.title(r'$\theta_1$ vs t')
    plt.grid()
    
    plt.show()


    plt.plot(t, theta_2, 'b-')
    plt.xlabel('t [s]')
    plt.ylabel(r'$\theta_2$(t) [rad]')
    plt.title(r'$\theta_2$ vs t')
    plt.grid()
    
    plt.show()
    
    #coordinate cartesiane
    
    L_1, L_2 = 1.0, 1.0
    
    x_1 =  L_1*np.sin(theta_1)
    y_1 = -L_1*np.cos(theta_1)
    
    x_2 = x_1 + L_2*np.sin(theta_2)
    y_2 = y_1 - L_2*np.cos(theta_2)
    
    
    plt.scatter(0,0, s = 50, color='black', label='Pivot')
    plt.plot(x_1, y_1, linewidth = 0.5, alpha = 0.5, label='massa 1')
    plt.plot(x_2, y_2, color='purple', linewidth = 0.5, alpha = 0.5, label='massa 2')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Piano cartesiano')
    plt.grid()
    
    plt.show()

    

    #analisi di fourier per le ampiezze dei theta(t)    
    
    
    freq_positive, ampiezze_theta_1, nu_max, A_max = analisi_fourier(theta_1, h, label=r'$\theta_1$')
    
    freq_positive, ampiezze_theta_2, nu_max, A_max = analisi_fourier(theta_2, h, label=r'$\theta_2$')

    
    plt.plot(freq_positive, ampiezze_theta_1, color='red')
    plt.xlabel(r'$\nu$ $[s^{-1}]$')
    plt.ylabel(r'Ampiezza $\theta_1$')
    plt.xscale('log')
    plt.title(r'Ampiezza $\theta_1$ vs $\nu$')
    plt.grid()
    plt.show()
    
    plt.plot(freq_positive, ampiezze_theta_2, color='blue')
    plt.xlabel(r'$\nu$ $[s^{-1}]$')
    plt.ylabel(r'Ampiezza $\theta_2$')
    plt.xscale('log')
    plt.title(r'Ampiezza $\theta_2$ vs $\nu$')
    plt.grid()
    plt.show()
    
     
 #Massimo di ampiezza per $\theta_1$: 7799.4437 rad  a frequenza ν = 0.4000 Hz
 #Massimo di ampiezza per $\theta_2$: 92607.6687 rad  a frequenza ν = 0.0500 Hz

    #Esercizio 3
    
    
    h = 0.001
    t_max = 30.0
    
    L_1, L_2 = 1.0, 1.0
    m_1, m_2 = 1.0, 1.0
    
    t = np.arange(0, t_max, h)
    
    
    Y_sol = np.zeros((len(t), 4))  #matrice di zeri per theta_1, theta_2, omega_1, omega_2
    
    Y_sol[0] = [np.pi/3, np.pi/4, 0.0, 0.0]  # valori iniziali per theta_i, omega_i
    
    for i in range(len(t) - 1):
        
        Y_sol[i+1] = RK4(dY_dt, t[i], Y_sol[i], h)
        
    #def gli array dei dati    
    
    theta_1 = Y_sol[:, 0] 
    theta_2 = Y_sol[:, 1] 
    omega_1 = Y_sol[:, 2] 
    omega_2 = Y_sol[:, 3] 

    E = energia(m_1, m_2, L_1, L_2, theta_1, theta_2, omega_1, omega_2)
    
    E_0 = E[0]  #energia al tempo t = 0
    
    epsilon = ((E - E_0)/np.abs(E_0)) #discrepanza 
    
    epsilon_max = np.max(epsilon) #1.6772095370500467e-12
    
    print(f"Variazione massima osservata: {epsilon_max}")
            
    
    plt.plot(t, epsilon, linewidth = 0.5,color='purple')
    plt.xlabel('t [s]')
    plt.ylabel(r'$\Delta E / E_0$')
    plt.title(r'$\Delta E / E_0$ vs t')
    plt.grid()
    
    plt.show()
    
    
    
    #fft
    
    
    freq, amp, nu_max, A_max = analisi_fourier(epsilon, h, label=r'$\epsilon$')
    
    plt.plot(freq, amp, linewidth = 0.5, color='purple')
    plt.xlabel(r'$\nu$ $[s^{-1}]$')
    plt.ylabel(r'Ampiezza $\Delta E / E_0$')
    plt.xscale('log')
    plt.title(r'$\Delta E / E_0$ vs $\nu$')
    plt.grid()
    plt.show()
    
    
            
    

    #Esercizio 4
    
    
    h = 0.001
    t_max = 30.0
    
    L_1, L_2 = 1.0, 1.0
    m_1, m_2 = 1.0, 1.0
    
    t = np.arange(0, t_max, h)
    
    
    Y_sol_A = np.zeros((len(t), 4))  
    Y_sol_A[0] = [np.pi/2, np.pi/2, 0.0, 0.0]
    
    Y_sol_B = np.zeros((len(t), 4))  
    Y_sol_B[0] = [np.pi/2 + 10**(-6), np.pi/2, 0.0, 0.0]
    
    for i in range(len(t) - 1):
        
        Y_sol_A[i+1] = RK4(dY_dt, t[i], Y_sol_A[i], h)
        Y_sol_B[i+1] = RK4(dY_dt, t[i], Y_sol_B[i], h)
        
    #def gli array dei dati    
    
    theta_1_A = Y_sol_A[:, 0] 
    
    theta_2_A = Y_sol_A[:, 1] 
     
    theta_1_B = Y_sol_B[:, 0]
    
    theta_2_B = Y_sol_B[:, 1]
    
    abs_delta_theta_1 = np.abs(theta_1_A - theta_1_B)
    
    abs_delta_theta_2 = np.abs(theta_2_A - theta_2_B)

    indice_lyapunov = np.argmax(abs_delta_theta_1 > 1) #cerca il primo momento in cui 
    #la condizione è vera e ti restituisce l'indice
    
    print(f"Tempo di Lyapunov: {t[indice_lyapunov]} s") #18.666 s
    
    plt.plot(t, theta_1_A, linewidth = 0.5, label=r'$\theta_{1, A}$ (t)')
    plt.plot(t, theta_1_B, color='purple', linewidth = 0.5, label=r'$\theta_{1, B}$ (t)')
    plt.legend()
    plt.xlabel('t [s]')
    plt.ylabel(r'$\theta$ (t) [rad]')
    plt.title('Plot')
    plt.grid()
    
    plt.show()
    
    plt.plot(t, theta_2_A, linewidth = 0.5, label=r'$\theta_{1, A}$ (t)')
    plt.plot(t, theta_2_B, color='purple', linewidth = 0.5, label=r'$\theta_{1, B}$ (t)')
    plt.legend()
    plt.xlabel('t [s]')
    plt.ylabel(r'$\theta$ (t) [rad]')
    plt.title('Plot')
    plt.grid()
    
    plt.show()
    
    plt.plot(t, abs_delta_theta_1, linewidth = 0.5, label = r'|$\Delta \theta_1$|')
    plt.plot(t, abs_delta_theta_2, linewidth = 0.5, color='purple', label = r'|$\Delta \theta_2$|')
    plt.legend()
    plt.xlabel('t [s]')
    plt.ylabel(r'|$\Delta \theta$|')
    plt.title(r'|$\Delta \theta$| vs t')
    #plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    
    plt.show()
    
    #Esercizio 5
    
    
    h = 0.001
    t_max = 50.0
    
    t = np.arange(0, t_max, h)
    
    Y_sol = np.zeros((len(t), 4))  #matrice di zeri per theta_1, theta_2, omega_1, omega_2
    
    Y_sol[0] = [np.pi/2, np.pi/2, 0.0, 0.0]  # valori iniziali per theta_i, omega_i
    
    for i in range(len(t) - 1):
        
        Y_sol[i+1] = RK4(dY_dt, t[i], Y_sol[i], h)



    theta_1 = Y_sol[:, 0] 
    theta_2 = Y_sol[:, 1] 
    omega_1 = Y_sol[:, 2] 
    omega_2 = Y_sol[:, 3]
    
    plt.plot(theta_1, omega_1, 'r-', linewidth = 0.5)
    plt.ylabel(r'$\omega_1$ [rad/s]')
    plt.xlabel(r'$\theta_1$ [rad]')
    plt.title(r'Spazio delle fasi ($\theta_1, \omega_1$)')
    plt.grid()
    
    plt.show()


    plt.plot(theta_2, omega_2, 'b-', linewidth = 0.5)
    plt.ylabel(r'$\omega_2$ [rad/s]')
    plt.xlabel(r'$\theta_2$ [rad]')
    plt.title(r'Spazio delle fasi ($\theta_2, \omega_2$)')
    plt.grid()
    
    plt.show()




#mi serve un vettore Y = [theta_1, theta_2, omega_1, omega_2].
#ne  calcolo la derivata prima



def dY_dt(t, Y):
    
    theta_1, theta_2, omega_1, omega_2 = Y
    
    g = 9.81
    L1, L2 = 1.0, 1.0
    m1, m2 = 1.0, 1.0
    
    delta = theta_1 - theta_2
    
    den1 = L1**2 * ((m1 + m2) - m2 * np.cos(delta)**2)
    
    
    alpha_1 = (- m2 * L1 * omega_1**2 * np.sin(delta) * np.cos(delta) 
              + m2 * g * np.sin(theta_2) * np.cos(delta) 
              - m2 * L2 * omega_2**2 * np.sin(delta) 
              - (m1 + m2) * g * np.sin(theta_1)) / den1
              
    alpha_2 = (L1 * omega_1**2 * np.sin(delta) - g * np.sin(theta_2) 
               - L1 * np.cos(delta) * alpha_1) / L2
    
    return np.array([omega_1, omega_2, alpha_1, alpha_2]) #array della derivata prima


def RK4(dY_dt, t, Y, h):  #funzione generale per RK4
    
    
    k_1 = dY_dt(t, Y)
    k_2 = dY_dt(t + (h/2), Y+k_1*(h/2))
    k_3 = dY_dt(t + (h/2), Y+k_2*(h/2))
    k_4 = dY_dt(t + h, Y+k_3*h) 
    
    return Y + (h/6)*(k_1 + 2*k_2 + 2*k_3 + k_4) #mi ritorna la media pesata




def energia(m_1, m_2, L_1, L_2, theta_1, theta_2, omega_1, omega_2, g=9.81):
    
    T = 0.5*m_1*(L_1*omega_1)**2 + 0.5*m_2*((L_1*omega_1)**2 + (L_2*omega_2)**2 
                + 2*L_1*L_2*omega_1*omega_2*np.cos(theta_1 - theta_2))

    V = -(m_1 + m_2)*g*L_1*np.cos(theta_1) - m_2*g*L_2*np.cos(theta_2)

    E = T + V
    
    return E


def analisi_fourier(vector, h, label='θ'):
    """
    Esegue l'analisi di Fourier su un segnale angolare.
    
    parametri da inserire
    
    vector : array numpy; Array dei valori dipendenti dall tempo.
    h : float Passo temporale di campionamento [s]
        
    label : str Etichetta
        
    """
    
    N = len(vector)
    
    trasf_fourier = np.fft.fft(vector)
    freq = np.fft.fftfreq(N, h) #campionamento frequenze: numero n dei campionamenti len(t), passo d del campionamento h
    
    # Solo frequenze positive, escluso il picco DC a freq=0
    freq_positive  = freq[1 : N//2]
    ampiezze       = np.abs(trasf_fourier)[1 : N//2]
    
    indice_max  = np.argmax(ampiezze) #calcola l'indice della freq in corrispondenza del max
    freq_max    = freq_positive[indice_max]
    ampiezza_max = ampiezze[indice_max]
    
    print(rf'Massimo di ampiezza per {label}: {ampiezza_max:.4f} rad  '
          rf'a frequenza $\nu$ = {freq_max:.4f} Hz')
    
    return freq_positive, ampiezze, freq_max, ampiezza_max




main()








