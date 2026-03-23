#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 09:20:25 2026

@author: alberto
"""


import numpy as np
from matplotlib import pyplot as plt


#Esercizio 1


def velocity_verlet(t, h, theta_0, omega_in, a_0, L, g):
    
    theta = np.zeros(len(t))  #def gli array pieni di zeri
    omega = np.zeros(len(t))
    a     = np.zeros(len(t))
    
    theta[0] = theta_0   # def i valori iniziali per theta e omega
    omega[0] = omega_in
    a[0]     = a_0
    
    for i in range(len(t) - 1):
        
        theta[i+1] = theta[i] + h*omega[i] + 0.5*(h**2)*a[i]
        a[i+1]     = -(g/L)*np.sin(theta[i+1])
        omega[i+1] = omega[i] + 0.5*h*(a[i]+a[i+1])
        
    
    return theta, omega, a #restituisce una tupla
    



def RK4(t, h, theta_0, omega_in, L, g):
    
    theta = np.zeros(len(t))  #def gli array pieni di zeri
    omega = np.zeros(len(t))
    
    theta[0] = theta_0   # def i valori iniziali per theta e omega
    omega[0] = omega_in
    
    for i in range(len(t) - 1):
    
        k_1_theta = omega[i]
        k_1_omega = -(g/L)*np.sin(theta[i])
        
        k_2_theta = omega[i] + (h/2)*k_1_omega
        k_2_omega = -(g/L)*np.sin(theta[i] + (h/2)*k_1_theta)
        
        k_3_theta = omega[i] + (h/2)*k_2_omega
        k_3_omega = -(g/L)*np.sin(theta[i] + (h/2)*k_2_theta)

        k_4_theta = omega[i] + h*k_3_omega
        k_4_omega = -(g/L)*np.sin(theta[i] + h*k_3_theta)

        theta[i+1] = theta[i] + (h/6)*(k_1_theta + 2*k_2_theta + 2*k_3_theta + k_4_theta)
        omega[i+1] = omega[i] + (h/6)*(k_1_omega + 2*k_2_omega + 2*k_3_omega + k_4_omega)
        
        
    return theta, omega
        
h = 0.01 #s
t_max = 20.0 #s
L = 1 #m
theta_0 = 0.5 #rad
g = 9.81 #m/s^2
t_0 = 0.0 #s
omega_in = 0.0 #rad/s
a_0 = -(g/L)*np.sin(theta_0)

t = np.arange(t_0, t_max, h)


theta_verlet, omega_verlet, a_verlet = velocity_verlet(t, h, theta_0, omega_in, a_0, L, g)

theta_RK4, omega_RK4 = RK4(t, h, theta_0, omega_in, L, g)

#blocco di plot

plt.plot(t, theta_verlet, label="Velocity Verlet")
plt.plot(t, theta_RK4, "r--", label="RK4")

plt.xlabel("t (s)")
plt.ylabel(r"$\theta (t)$ (rad)")
plt.legend()
plt.grid(True)

plt.show()


##############################################################################

#Esercizio 2



h = 0.01 #s
t_max = 500.0 #s
L = 1 #m
theta_0 = 0.5 #rad
g = 9.81 #m/s^2
t_0 = 0.0 #s
omega_in = 0.0 #rad/s
a_0 = -(g/L)*np.sin(theta_0)
m = 1 #massa unitaria
t = np.arange(t_0, t_max, h)


E_0 = 0.5*m*(L*omega_in)**2 + m*g*L*(1 - np.cos(theta_0))

def energia(theta, omega, L, m, g):
    
    return 0.5*m*(L*omega)**2 + m*L*g*(1 - np.cos(theta))


theta_verlet, omega_verlet, a_verlet = velocity_verlet(t, h, theta_0, omega_in, a_0, L, g)
theta_RK4, omega_RK4 = RK4(t, h, theta_0, omega_in, L, g)


energia_verlet = energia(theta_verlet, omega_verlet, L, m, g)
energia_RK4    = energia(theta_RK4, omega_RK4, L, m, g)

delta_E_verlet = E_0 - energia_verlet
delta_E_RK4    = E_0 - energia_RK4 

U_V = (delta_E_verlet/E_0)*100
U_RK4 = (delta_E_RK4/E_0)*100

plt.plot(t, U_V)

plt.xlabel("t (s)")
plt.ylabel(r"$\Delta E/E_0$")
plt.grid(True)
plt.title('Variazione Verlet')

plt.show()

plt.plot(t, U_RK4)
plt.xlabel("t (s)")
plt.ylabel(r"$\Delta E/E_0$")
plt.title('Variazione RK4')
plt.grid(True)

plt.show()

##############################################################################

#Esercizio 3

#diagrammi delle fasi: riutilizzo gli stessi di prima



plt.plot(theta_verlet, omega_verlet, 'r-')

plt.xlabel(r"$\theta (t)$ (rad)")
plt.ylabel(r"$\omega (t)$ (rad/s)")
plt.grid(True)
plt.title(r'$\omega$ vs $\theta$ Verlet')

plt.show()

plt.plot(theta_RK4, omega_RK4, 'b-')

plt.xlabel(r"$\theta (t)$ (rad)")
plt.ylabel(r"$\omega (t)$ (rad/s)")
plt.grid(True)
plt.title(r'$\omega$ vs $\theta$ RK4')

plt.show()



#ordine di convergenza


h_values = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]

errore_verlet = []
errore_RK4    = []

h_c = 0.0001
t_finale = 10.0
N_c   = round(t_finale / h_c)          # numero di passi esatto
t_c   = np.linspace(0, t_finale, N_c + 1)
theta_0 = 0.5
omega_in = 0
L = 1
g = 9.81
a_0 = -(g/L)*np.sin(theta_0)

theta_RK4_c = RK4(t_c, h_c, theta_0, omega_in, L, g)[0] 

theta_RK4_c_last = theta_RK4_c[-1]  #valore campione


for h in h_values:
    
    N = round(t_finale / h) # numero di passi esatto, serve a evitare che 
    #divisioni float tipo 10.0/0.1 diano 99.9999... invece di 100
    t = np.linspace(0, t_finale, N + 1) # t[-1] = 10.0 garantito
    
    theta_verlet = velocity_verlet(t, h, theta_0, omega_in, a_0, L, g)[0]
    theta_RK4 = RK4(t, h, theta_0, omega_in, L, g)[0]
    
    diff_verlet = theta_verlet[-1] - theta_RK4_c_last
    diff_RK4    = theta_RK4[-1] - theta_RK4_c_last
    
    errore_verlet.append(np.abs(diff_verlet))
    errore_RK4.append(np.abs(diff_RK4))


plt.plot(h_values, errore_verlet, 'r-o', label='Errore Verlet')

plt.plot(h_values, errore_RK4, 'b-o', label='Errore RK4')


plt.legend()
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$h$")
plt.ylabel(r"$\varepsilon$")
plt.show()

# Stima ordine di convergenza 

log_h = np.log(h_values)
ordine_verlet = np.polyfit(log_h, np.log(errore_verlet), 1)[0]
ordine_RK4    = np.polyfit(log_h, np.log(errore_RK4),    1)[0]
print(f"Ordine Verlet: {ordine_verlet:.2f}")  # atteso 2, esce 1.99
print(f"Ordine RK4:    {ordine_RK4:.2f}") # atteso 4, esce 4.04



