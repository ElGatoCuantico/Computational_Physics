#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 09:57:22 2026

@author: alberto
"""

import numpy as np
from matplotlib import pyplot as plt

#Esercizio 1

def eulero(t, h, theta_0, omega_in, L, g):
    
    theta = np.zeros(len(t))  #def gli array pieni di zeri
    omega = np.zeros(len(t))
    
    theta[0] = theta_0   # def i valori iniziali per theta e omega
    omega[0] = omega_in
    
    for i in range(len(t) - 1):
        
        theta[i+1] = theta[i] + h*omega[i]
        omega[i+1] = omega[i] - h*(g/L)*np.sin(theta[i])
        
    
    return theta, omega #restituisce una tupla
    
    
def eulero_cromer(t, h, theta_0, omega_in, L, g):
    
    theta = np.zeros(len(t))  #def gli array pieni di zeri
    omega = np.zeros(len(t))
    
    theta[0] = theta_0   # def i valori iniziali per theta e omega
    omega[0] = omega_in
    
    for i in range(len(t) - 1):
        
        omega[i+1] = omega[i] - h*(g/L)*np.sin(theta[i])
        theta[i+1] = theta[i] + h*omega[i+1]
        
        
    
    return theta, omega

#def le costanti da usare

h = 0.01 #s
t_max = 20.0 #s
L = 1 #m
theta_0 = 0.1 #rad
omega_in = 0.0 #rad/s
g = 9.81 #m/s^2
t_0 = 0.0 #s
w_0 = np.sqrt(g/L) #s^-1


t = np.arange(t_0, t_max, h)


#def array per i due metodi + confronto sol analitica

theta_eulero, omega_eulero = eulero(t, h, theta_0, omega_in, L, g)

theta_eulero_cromer, omega_eulero_cromer = eulero_cromer(t, h, theta_0, omega_in, L, g)

theta_analitica = theta_0*np.cos(w_0*t) + (omega_in/w_0)*np.sin(w_0*t)

#blocco di plot

plt.plot(t, theta_eulero, label="Eulero")
plt.plot(t, theta_eulero_cromer, label="Eulero-Cromer")
plt.plot(t, theta_analitica, label="Analitica")

plt.xlabel("t (s)")
plt.ylabel(r"$\theta (t)$ (rad)")
plt.legend()
plt.grid(True)

plt.show()


#errori 


e_glob_theta_eulero = np.abs(theta_analitica - theta_eulero) #vettore degli errori di discretizzazione globale al passo i-esimo
e_medio_theta_eulero = np.mean(e_glob_theta_eulero) #errore medio
e_max_theta_eulero = np.max(e_glob_theta_eulero)


print("Errore globale medio (Eulero): %.4f" % e_medio_theta_eulero)
print("Errore globale massimo (Eulero) : %.4f" % e_max_theta_eulero)

e_glob_theta_eulero_cromer = np.abs(theta_analitica - theta_eulero_cromer) #vettore degli errori di discretizzazione globale al passo i-esimo
e_medio_theta_eulero_cromer = np.mean(e_glob_theta_eulero_cromer) #errore medio
e_max_theta_eulero_cromer = np.max(e_glob_theta_eulero_cromer)

print("Errore globale medio (Eulero-Cromer): %.4f" % e_medio_theta_eulero_cromer)
print("Errore globale massimo (Eulero-Cromer) : %.4f" % e_max_theta_eulero_cromer)


#variazione di h: plot di e(h) 


h_values = np.array([0.1, 0.05, 0.025, 0.01, 0.005,0.0025, 0.001])

errori_max_values_eulero = []
errori_max_values_eulero_cromer = []

for h in h_values:   #ciclo sui valori di h
    
    t = np.arange(t_0, t_max, h) #ricalcolo il vettore dei tempi
    
    #ricalcolo le varie theta
    
    theta_eulero = eulero(t, h, theta_0, omega_in, L, g)[0]
    theta_eulero_cromer = eulero_cromer(t, h, theta_0, omega_in, L, g)[0]
    
    theta_analitica = theta_0*np.cos(w_0*t) + (omega_in/w_0)*np.sin(w_0*t)
    
    #calcolo vettore dell'errore globale al passo n
    
    errore_eulero = np.abs(theta_analitica - theta_eulero)
    errore_eulero_cromer = np.abs(theta_analitica - theta_eulero_cromer)
    
    e_max_eulero = np.max(errore_eulero)
    e_max_eulero_cromer = np.max(errore_eulero_cromer)

    #append del valore massimo dell'errore all'array degli errori
    
    errori_max_values_eulero.append(e_max_eulero)
    errori_max_values_eulero_cromer.append(e_max_eulero_cromer)

#plot 

plt.title("$e_{n, max}$ max in funzione di h (Eulero)")
plt.scatter(h_values, errori_max_values_eulero, marker="o")
plt.xlabel("h")
plt.ylabel(r"$e_{n, max}$")
plt.grid(True)
plt.show()


plt.title("$e_{n, max}$ max in funzione di h (Eulero-Cromer)")
plt.scatter(h_values, errori_max_values_eulero_cromer, marker="o")
plt.xlabel("h")
plt.ylabel(r"$e_{n, max}$")
plt.grid(True)

plt.show()


#####################################################################


#Esercizio 2

#def le costanti da usare

theta_0 = 0.5

h = 0.01 #s
L = 1 #m
omega_in = 0.0 #rad/s
g = 9.81 #m/s^2
t_0 = 0.0 #s
w_0 = np.sqrt(g/L) #s^-1
T = 2*np.pi/w_0

t_max = 100*T

t = np.arange(t_0, t_max, h)

m = 1 #suppongo la massa uguale a 1kg

#energia meccanica 

def energia(theta, omega, L, m, g):
    
    return 0.5*m*(L*omega)**2 + m*L*g*(1 - np.cos(theta))

#def gli array di theta e omega


theta_eulero, omega_eulero = eulero(t, h, theta_0, omega_in, L, g)
theta_eulero_cromer, omega_eulero_cromer = eulero_cromer(t, h, theta_0, omega_in, L, g)


#def gli array dell'energia in entrambi i metodi


energia_eulero = energia(theta_eulero, omega_eulero, L, m, g)
energia_eulero_cromer = energia(theta_eulero_cromer, omega_eulero_cromer, L, m, g)


#plot

plt.plot(t, energia_eulero, label="Eulero")
plt.plot(t, energia_eulero_cromer, label="Eulero-Cromer")

plt.xlabel("t (s)")
plt.ylabel("E (J)")
plt.title("Energia in funzione del tempo")

plt.ylabel("E (J)")
plt.grid(True)
plt.legend()

plt.show()

#errore percentuale


err_eulero = ((energia_eulero[-1] - energia_eulero[0])/energia_eulero[0])*100
err_eulero_cromer = ((energia_eulero_cromer[-1] - energia_eulero_cromer[0])/energia_eulero_cromer[0])*100


print("Errore percentuale Eulero: %.5f" % err_eulero)
print("Errore percentuale Eulero-Cromer: %.5f" % err_eulero_cromer)



#####################################################################

#Esercizio 3

#def le costanti da usare


theta_0 = [np.pi/2, 3]  # due valori per theta_0. conviene salvarli in una lista

h = 0.01 #s
L = 1 #m
omega_in = 0.0 #rad/s
g = 9.81 #m/s^2
t_0 = 0.0 #s
w_0 = np.sqrt(g/L) #s^-1
T = 2*np.pi/w_0

t_max = 5*T

t = np.arange(t_0, t_max, h)


for theta in theta_0:


    theta_eulero, omega_eulero = eulero(t, h, theta, omega_in, L, g)

    theta_eulero_cromer, omega_eulero_cromer = eulero_cromer(t, h, theta, omega_in, L, g)

    theta_analitica = theta*np.cos(w_0*t) + (omega_in/w_0)*np.sin(w_0*t)


    plt.title(fr"$\theta_0$ = {theta: .2f}")
    plt.plot(t, theta_eulero, label="Eulero")
    plt.plot(t, theta_eulero_cromer, label="Eulero-Cromer")
    plt.plot(t, theta_analitica, label="Analitica")


    plt.xlabel("t (s)")
    plt.ylabel(r"$\theta (t)$ (rad)")
    plt.legend()
    plt.grid(True)

    plt.show()










