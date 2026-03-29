#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 21:36:47 2026

@author: alberto
"""
import numpy as np
from matplotlib import pyplot as plt




def main():
    
    #Esercizio 1
    
    A = 0.0
    x_0 = 1.0
    v_0 = 0.0
    omega_0 = 1.0
    omega = 0.0
    gamma = 0.1
    m = 1.0
    h = 0.01
    t_max = 50
    
    
    t = np.arange(0, t_max, h)
    
    x, v = RK4(t, h, x_0, v_0, omega_0, gamma, A, omega)

    decad = np.exp(-gamma*t/2)  #decadimento esponenziale

    E = energia(x, v, m, omega_0)

    plt.plot(t, x, 'r-' ,label='x(t)')
    plt.plot(t, decad, 'b--', label=r'$e^{\frac{-\gamma t}{2}}$')
    plt.xlabel("t[s]")
    plt.ylabel("x(t)[m]")
    plt.legend()
    plt.grid(True)
    plt.title("Posizione in funzione del tempo in assenza di forzante")
    
    plt.show()

    plt.plot(t, E)
    plt.xlabel("t[s]")
    plt.ylabel("E(t)[J]")
    plt.grid(True)
    plt.title("Energia in funzione del tempo") #decadimento exp

    plt.show()

    #Esercizio 2
    
    
    A = 1.0
    x_0 = 0.0
    v_0 = 0.0
    omega_0 = 1.0
    omega = 0.8
    gamma = 0.2
    m = 1.0
    h = 0.01
    t_max = 100
    
    A_staz = A_s(A, omega_0, omega, gamma) #valore teorico 2.5384
    
    t = np.arange(0, t_max, h)
    
    x, v = RK4(t, h, x_0, v_0, omega_0, gamma, A, omega)


    E = energia(x, v, m, omega_0)
    
    
    plt.plot(t, x, 'r-')
    plt.xlabel("t [s]")
    plt.ylabel("x(t) [m]")
    plt.grid(True)
    plt.title("Posizione in funzione del tempo con forzante")
    
    plt.show()

    plt.plot(t, E)
    plt.xlabel("t [s]")
    plt.ylabel("E(t) [J]")
    plt.grid(True)
    plt.title("Energia in funzione del tempo")

    plt.show()
    
    #dal plot dell'energia (più facile da vedere) sembra che il transitorio sia di circa 40s
    
    #indice corrispondente a t = 40s
    #taglia il transitorio a 5 tau secondi, ma se 5 tau è troppo 
    #grande rispetto alla simulazione, taglia ai 2/3 del tempo totale
    tau = (2 / gamma)
    i_trans = min(int(5 * tau / h), int(len(t) * 2 / 3)) #divido per h per avere l'indice

    #ampiezza numerica = max di |x| dopo il transitorio
    A_num = np.max(np.abs(x[i_trans:]))

    print("A_s th  = %.4f" % A_staz) #val th 2.5384
    print("A_s num = %.4f" % A_num)  #val num 2.5503
    print("Errore relativo = %.2f %%" % (np.abs(A_staz - A_num)/A_staz * 100)) # 0.47 %
    
    
    #Esercizio 3
    
    
    A = 1.0
    x_0 = 0.0
    v_0 = 0.0
    omega_0 = 1.0
    omega_values = np.arange(0.5, 1.5, 0.01) 
    gamma = 0.1
    m = 1.0
    h = 0.01
    t_max = 200
    
    t = np.arange(0, t_max, h)
    
    A_s_values_th  = np.zeros(len(omega_values))
    A_s_values_num = np.zeros(len(omega_values))
    epsilon = np.zeros(len(omega_values))

    for i in range(len(omega_values)):
    
        omega = omega_values[i]  #prendo l'elemento i-esimo e lo sostituisco a 0
        
        A_staz_th = A_s(A, omega_0, omega, gamma)
        A_s_values_th[i] = A_staz_th  # lo stesso qui
        
        x = RK4(t, h, x_0, v_0, omega_0, gamma, A, omega)[0]
        tau = (2 / gamma)
        i_trans = min(int(5 * tau / h), int(len(t) * 2 / 3))
        A_num = np.max(np.abs(x[i_trans:]))
        A_s_values_num[i] = A_num  # anche qui
        
        errore_rel = (np.abs(A_staz_th - A_num) / A_staz_th * 100)
        epsilon[i] = errore_rel  # e qui
        
    
    err_rel_medio = np.mean(epsilon)
    
    print("Errore relativo medio: %.2f %%" % err_rel_medio) #0.37%
    
    plt.plot(omega_values, A_s_values_num, 'b-', label='numerico')
    plt.plot(omega_values, A_s_values_th, 'r--', label='teorico')
    plt.xlabel(r"$\omega$")
    plt.ylabel(r"$A_s (\omega)$")
    plt.grid(True)
    plt.legend()
    plt.title(r"$A_s$ vs $\omega$ ($\gamma = 0.1$)")
    
    plt.show()
    
    
    
    A = 1.0
    x_0 = 0.0
    v_0 = 0.0
    omega_0 = 1.0
    omega_values = np.arange(0.5, 1.5, 0.01) 
    gamma = 0.3
    m = 1.0
    h = 0.01
    t_max = 200
    
    t = np.arange(0, t_max, h)
    
    A_s_values_th  = np.zeros(len(omega_values))
    A_s_values_num = np.zeros(len(omega_values))
    epsilon = np.zeros(len(omega_values))

    for i in range(len(omega_values)):
    
        omega = omega_values[i]  #prendo l'elemento i-esimo e lo sostituisco a 0
        
        A_staz_th = A_s(A, omega_0, omega, gamma)
        A_s_values_th[i] = A_staz_th  # lo stesso qui
        
        x = RK4(t, h, x_0, v_0, omega_0, gamma, A, omega)[0]
        tau = (2 / gamma)
        i_trans = min(int(5 * tau / h), int(len(t) * 2 / 3))
        A_num = np.max(np.abs(x[i_trans:]))
        A_s_values_num[i] = A_num  # anche qui
        
        errore_rel = (np.abs(A_staz_th - A_num) / A_staz_th * 100)
        epsilon[i] = errore_rel  # e qui
        
    
    
    err_rel_medio = np.mean(epsilon)
    
    print("Errore relativo medio: %.2f %%" % err_rel_medio) #0.24%
    
    plt.plot(omega_values, A_s_values_num, 'b-', label='numerico')
    plt.plot(omega_values, A_s_values_th, 'r--', label='teorico')
    plt.xlabel(r"$\omega$")
    plt.ylabel(r"$A_s (\omega)$")
    plt.grid(True)
    plt.legend()
    plt.title(r"$A_s$ vs $\omega$ ($\gamma = 0.3$)")
    
    plt.show()
    
    #Esercizio 4
    
    h = 0.01
    t_max = 150
    
    t = np.arange(0, t_max, h)
    
    x, v = RK4(t, h, x_0=1.0, v_0=0.0, omega_0=1.0, gamma=0.1, A=0.0, omega=0.0)
    x_forz, v_forz = RK4(t, h, x_0, v_0, omega_0, gamma=0.2, A=1, omega=1.0)

    tau = (2 / 0.2) #tau forzante
    i_trans = min(int(5 * tau / h), int(len(t) * 2 / 3)) #transitorio forzante

    plt.plot(x, v, 'b--')
    
    plt.xlabel("x(t) [m]")
    plt.ylabel("v(t) [m/s]")
    plt.grid(True)
    plt.title("Spazio delle fasi in assenza di forzante")
    
    plt.show()
    
    plt.plot(x_forz[i_trans:], v_forz[i_trans:], 'r--', label='Forzante')
    
    plt.xlabel("x(t) [m]")
    plt.ylabel("v(t) [m/s]")
    plt.grid(True)
    plt.title("Spazio delle fasi con forzante dopo il transitorio")
    
    plt.show()
    
    
    
    
    
    
    
def RK4(t, h, x_0, v_0, omega_0, gamma, A, omega):
    
    x = np.zeros(len(t))  #def gli array pieni di zeri
    v = np.zeros(len(t))
    
    x[0] = x_0   # def i valori iniziali per pos e vel
    v[0] = v_0
    
    for i in range(len(t) - 1):
    
        k_1_x = v[i]
        k_1_v = -(omega_0**2)*x[i] - gamma*v[i] + A*np.cos(omega*t[i])
        
        k_2_x = v[i] + (h/2)*k_1_v
        k_2_v = -(omega_0**2)*(x[i] + (h/2)*k_1_x) - gamma*(v[i] + (h/2)*k_1_v) + A*np.cos(omega*(t[i] + (h/2)))
        
        k_3_x = v[i] + (h/2)*k_2_v
        k_3_v = -(omega_0**2)*(x[i] + (h/2)*k_2_x) - gamma*(v[i] + (h/2)*k_2_v) + A*np.cos(omega*(t[i] + (h/2)))

        k_4_x = v[i] + h * k_3_v
        k_4_v = -(omega_0**2)*(x[i] + h*k_3_x) - gamma*(v[i] + h*k_3_v) + A*np.cos(omega*(t[i] + h))
        
        x[i+1] = x[i] + (h/6)*(k_1_x + 2*k_2_x + 2*k_3_x + k_4_x)
        v[i+1] = v[i] + (h/6)*(k_1_v + 2*k_2_v + 2*k_3_v + k_4_v)
        
        
    return x, v



def energia(x, v, m, omega_0):
    
    E = 0.5*m*(v**2 + (omega_0*x)**2)
    
    return E


def A_s(A, omega_0, omega, gamma):
    
    A_staz = A/np.sqrt((omega_0**2 - omega**2)**2 + (gamma*omega)**2)
    
    return A_staz

    
main()








                    
                    