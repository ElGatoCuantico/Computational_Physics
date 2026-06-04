#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:01:42 2026

@author: alberto
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import PowerNorm
import scienceplots


plt.style.use(['science', 'notebook', 'grid'])   #migliora la visibilità



def main():
    

    a  = 20.0
    dr = 0.05
    r  = np.arange(0, a, dr)
    V  = np.zeros(len(r))   #dentro la buca il potenziale e' zero
     
    #gaussiana + impulso k 
    r0, sigma, k = 6.0, 1.0, 5.0
    psi0 = np.exp(-(r - r0)**2 / (2*sigma**2)) * np.exp(1j*k*r)   #gaussiana * fase
    u0   = r * psi0                                              #u = r * psi
    u0[0] = 0.0; u0[-1] = 0.0                                    # i due muri
    u0   = u0 / np.sqrt(np.sum(np.abs(u0)**2) * dr)             # normalizzo: probabilita' = 1
     
    R0, I0 = np.real(u0), np.imag(u0)
    
    E_max   = 2/dr**2 + np.max(np.abs(V))     #autovalore max di H
    dt      = 0.5 * (2/E_max)                 #meta' del limite di stabilita'
    T       = 25.0
    n_passi = int(T/dt)
     
    #inizializzazione
    
    R, I = leapfrog(R0, I0, dr, V, dt, n_passi)
    prob = calcola_probabilita(R, I, u0)
    
    #stesso setup, stesso dt; T piccolo perché FTCS esplode in fretta
    
    T_cmp = 1.0
    n_cmp = int(T_cmp / dt)
    
    R, I    = leapfrog(R0, I0, dr, V, dt, n_cmp)
    norm_lf = np.sum(calcola_probabilita(R, I, u0), axis=0) * dr      # Visscher
    
    psi       = FTCS_schrodinger(u0, dr, V, dt, n_cmp)
    norm_ftcs = np.sum(np.abs(psi)**2, axis=0) * dr                   # |psi|^2 diretto
    
    t = np.arange(n_cmp) * dt
    plt.semilogy(t, norm_lf,   label="Leapfrog (Visscher)")
    plt.semilogy(t, norm_ftcs, label="FTCS (Eulero esplicito)")
    plt.xlabel("tempo"); plt.ylabel("probabilità totale (log)"); plt.legend()
    plt.show()
    
    
    #plot 1D
    
    fig, ax = plt.subplots(figsize=(8, 5))
    linea, = ax.plot(r, prob[:, 0], color='purple')            #disegno la curva UNA volta
    ax.set_xlim(0, a)
    ax.set_ylim(0, prob.max())
    ax.set_xlabel("r"); ax.set_ylabel(r"$|u|^2$")
    ax.set_title("Particella in buca infinita")
     
    def aggiorna(n):
        linea.set_ydata(prob[:, n])            #cambio solo i valori verticali
        ax.set_title("Particella in buca infinita   t = %.2f" % (n*dt))
        return (linea,)
     
    ani = FuncAnimation(fig, aggiorna,
                        frames=range(0, n_passi, 20),
                        interval=20, blit=False)
    ani.save("film.gif", writer=PillowWriter(fps=30), dpi=80)
        
        
    #plot 2d

    N2d = 300
    asse = np.linspace(-a, a, N2d)      #valori da -a a +a
    X, Y = np.meshgrid(asse, asse)      #due matrici: X = coord. x, Y = coord. y di ogni pixel    
    
    #ho creato una griglia 300 x 300 da -a ad a
    
    R = np.sqrt(X**2 + Y**2)
    
    r2 = r[:, None]**2    #r^2 come colonna che si allinea a ogni colonna-tempo
    prob_vera = np.divide(prob, r2,
                      out=np.zeros_like(prob), #dove non divido, lascio zero 
                      where=r2 > 1e-12)    #divido solo dove r != 0
    
    vmax = np.percentile(prob_vera[r > 0.3], 99.5)    # scala colore robusta al rumore

    
    
    def frame_2d(n):
        Z = np.interp(R, r, prob_vera[:, n])          # leggo la curva radiale a ogni distanza R
        Z[R > r.max()] = 0                            # fuori dalla sfera: niente
        return Z

    
    fig2, ax2 = plt.subplots(figsize=(6, 6))
    mappa = ax2.imshow(frame_2d(0), extent=[-a, a, -a, a], origin="lower",
                       cmap="inferno", norm=PowerNorm(gamma=0.4, vmin=0, vmax=vmax))
    ax2.set_xlabel("x"); ax2.set_ylabel("y")
    fig2.colorbar(mappa, ax=ax2, label=r"$|\psi|^2$", fraction=0.046)
    
    
    def aggiorna_2d(n):
        mappa.set_array(frame_2d(n))
        ax2.set_title("t = %.2f" % (n*dt))
        return (mappa,)
    
    ani2 = FuncAnimation(fig2, aggiorna_2d, frames=range(0, n_passi, 20),
                         interval=20, blit=False)
    
    ani2.save("heatmap.gif", writer=PillowWriter(fps=30), dpi=80)
    
    #plot 3D
    
    gamma  = 0.4                                   #alza i valori deboli 
    n_surf = 180
    asse_s = np.linspace(-a, a, n_surf)            #griglia piu' rada perché il 3D e' costoso
    Xs, Ys = np.meshgrid(asse_s, asse_s)
    Rs     = np.sqrt(Xs**2 + Ys**2)
    z_cap  = vmax**gamma
    
    fig3 = plt.figure(figsize=(8, 6))
    ax3  = fig3.add_subplot(111, projection='3d')  # assi 3D
    
    def aggiorna_surf(n):
        ax3.clear()                                 #la superficie NON si aggiorna ma si ridisegna
        Z = np.interp(Rs, r, prob_vera[:, n])
        Z[Rs > r.max()] = 0
        Z = np.clip(Z, 0, vmax)**gamma              #sollevo i valori deboli
        ax3.plot_surface(Xs, Ys, Z, cmap="inferno", vmin=0, vmax=z_cap,
                         linewidth=0, antialiased=False)
        ax3.set_zlim(0, z_cap * 1.2)
        ax3.set_xlabel("x"); ax3.set_ylabel("y"); ax3.set_zlabel(r"$|\psi|^2$")
        ax3.set_title("t = %.2f" % (n*dt))
        return ()
    
    ani3 = FuncAnimation(fig3, aggiorna_surf, frames=range(0, n_passi, 30),
                         interval=60, blit=False)
    
    ani3.save("heatmap_3D.gif", writer=PillowWriter(fps=30), dpi=80)
            
    


def derivata_seconda(u, dr):
    
    d2u_dr2 = np.zeros_like(u)
    
    d2u_dr2[1 : -1] = (1/dr**2)*(u[2 : ] - 2 * u[1 : -1] + u[ : -2]) 
    
    return d2u_dr2


def applica_H(u, dr, V, h_bar=1, m=1):
    
    d2u_dr2 = derivata_seconda(u, dr)
    
    H_u = -(h_bar**2 / (2*m)) * d2u_dr2 + V * u
    
    return H_u


def leapfrog(R0, I0, dr, V, dt, n_passi, h_bar=1):
    R = np.zeros((len(R0), n_passi))   #  R
    I = np.zeros((len(I0), n_passi))   #  I

    R[:, 0] = R0
    I[:, 0] = I0 - (dt/(2*h_bar)) * applica_H(R0, dr, V)   #avviamento con Eulero

    for n in range(n_passi - 1):
        R[:, n+1] = R[:, n] + (dt/h_bar)*applica_H(I[:, n],   dr, V)
        I[:, n+1] = I[:, n] - (dt/h_bar)*applica_H(R[:, n+1], dr, V)

    return R, I


def calcola_probabilita(R, I, u0):
    
    prob = np.zeros_like(R)
    prob[:, 0]  = np.abs(u0)**2                      
    prob[:, 1:] = R[:, 1:]**2 + I[:, 1:]*I[:, :-1]
    
    return prob


def FTCS_schrodinger(u0, dr, V, dt, n_passi, h_bar=1):
    
    u_r_t = np.zeros((len(u0), n_passi), dtype=complex)  
    
    u_r_t[:, 0] = u0
    
    for n in range(n_passi - 1):
        
        u_r_t[:, n+1] = u_r_t[:, n] - (1j*dt/h_bar)*applica_H(u_r_t[:, n], dr, V)

    return u_r_t



main()





