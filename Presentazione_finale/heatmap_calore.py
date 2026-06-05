#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 11:23:53 2026

@author: alberto
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter




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


   
    
T_analitica_finale = T_analitica[:, -1]
epsilon = np.max(np.abs(T_analitica_finale - T_x_t[:,-1]))



# --- animazione heatmap (la sbarra che si raffredda) ---
salto = 4
frames = range(0, len(t), salto)
altezza = 8                       # righe impilate per dare spessore alla striscia

def striscia(n):
    return np.tile(T_x_t[:, n], (altezza, 1))

fig_h, ax_h = plt.subplots(figsize=(9, 9))
img = ax_h.imshow(striscia(0), aspect='auto', origin='lower',
                      cmap='inferno', vmin=0.0, vmax=T_0,
                      extent=[0, L, 0, 1])
titolo_h = ax_h.set_title('t = 0.000 s')
ax_h.set_xlabel('x [m]'); ax_h.set_yticks([])
fig_h.colorbar(img, ax=ax_h, label='T [°C]')

def update_h(n):
    img.set_array(striscia(n))    # il pupazzo è l'immagine: riscrivo i pixel
    titolo_h.set_text(f't = {t[n]:.3f} s')
    return img, titolo_h

ani_h = FuncAnimation(fig_h, update_h, frames=frames, interval=30, blit=False)
ani_h.save("calore_heatmap.gif", writer=PillowWriter(fps=30), dpi=90)







