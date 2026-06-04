#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 24 17:09:57 2026

@author: alberto
"""

import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science', 'notebook', 'grid'])   



def figura_spirale(E=1.0, dt=0.45, N=9):
    
    
    t = np.linspace(0, 2*np.pi/E, 400)
    psi_esatta = np.exp(-1j*E*t)
    psi = 1 + 0j
    traj = [psi]
    
    
    for _ in range(N):
        
        psi = psi * (1 - 1j*E*dt)
        traj.append(psi)
        
        
    traj = np.array(traj)

    fig, ax = plt.subplots(figsize=(6, 6))
    
    ax.plot(psi_esatta.real, psi_esatta.imag, color="#1D9E75", lw=2.5,
            label="Esatta: $|\\psi|$ costante (cerchio)")
    
    ax.plot(traj.real, traj.imag, color="#C0392B", lw=2, marker='o', ms=4,
            label="FTCS/Eulero: $|\\psi|$ cresce (spirale)")
    
    ax.set_aspect('equal')
    
    ax.axhline(0, color='gray', lw=0.4)
    ax.axvline(0, color='gray', lw=0.4)
    ax.set_xlabel("Re$(\\psi)$")
    ax.set_ylabel("Im$(\\psi)$")
    ax.legend(loc='upper left')
    
    fig.savefig("spirale_eulero.png", dpi=110)  
    plt.show()

figura_spirale()











