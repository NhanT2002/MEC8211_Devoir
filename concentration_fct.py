# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 18:27:11 2025

@author: dorid
"""

import numpy as np

def concentration_anal(prm):
    r_range = np.linspace(0, prm.R, prm.N)
    # C_anal = 0.25*prm.S/prm.D_eff*prm.R*prm.R*(r_range*r_range/(prm.R*prm.R) - 1) + prm.C_e
    C_anal = 0.25*prm.S/prm.D_eff*prm.R**2*(r_range**2/(prm.R**2) - 1) + prm.C_e

    return C_anal

def concentration_mdf_o1(prm):
    N = prm.N
    r = np.linspace(0, prm.R, N)
    dr = prm.R/(N-1)
    
    #Initilisation de la matrice LHS
    A = np.zeros((N, N))
    # Conditions limites
    # Utilisation de l'approximation avant de la dérivée 1ère
    A[0,0] = -1
    A[0,1] = 1
    A[-1,-1] = 1
    
    #Initialisation du vecteur RHS
    B = np.zeros(N)
    
    #Remplir les coefficient
    for i in range(1, N-1): #Exclure les premiers et derniers noeuds déjà définis
        #Remplir les coefficients de la matrice A
        A[i,i-1] = 1/(dr*dr)
        A[i,i] = -(2/(dr*dr) + 1/(r[i]*dr))
        A[i,i+1] = 1/(dr*dr) + 1/(r[i]*dr)
        #Remplir les coefficients du vecteur B
        B[i] = prm.S/prm.D_eff
    
    B[-1] = prm.C_e
    
    #Résolution dy système matriciel
    C = np.linalg.solve(A, B)
    
    return r, dr, C

def concentration_mdf_o2(prm):
    N = prm.N
    r = np.linspace(0, prm.R, N)
    dr = prm.R/(N-1)
    
    
    #Initilisation de la matrice LHS
    A = np.zeros((N, N))
    # Conditions limites
    # Utilisation de l'approximation avant "Gear" de la dérivée 1ère
    A[0,0] = -3/(2*dr)
    A[0,1] = 4/(2*dr)
    A[0,2] = -1/(2*dr)
    A[-1,-1] = 1
    #Initialisation du vecteur RHS
    B = np.zeros(N)
    
    #Remplir la matrice LHS
    for i in range(1, N-1): #Exclure les premier et derniers noeuds
        A[i,i-1] = 1/(dr*dr) - 1/(r[i]*2*dr)
        A[i,i] = -2/(dr*dr)
        A[i,i+1] = 1/(dr*dr) + 1/(r[i]*2*dr)
        
        B[i] = prm.S/prm.D_eff
        
    B[-1] = prm.C_e
    #Resolution dy systeme matriciel
    C = np.linalg.solve(A, B)
    
    
    return r, dr, C