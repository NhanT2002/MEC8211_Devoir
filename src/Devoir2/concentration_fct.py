#Importation des modules
import numpy as np
import matplotlib.pyplot as plt

#Définition des fonctions de concentration
def concentration_anal_MMS(prm):
    '''
    Fonction qui calcule, selon la solution analytique, 
    la concentration de sel dans un pilier

    Parameters
    ----------
    prm : Objet class parametres()
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]
        - k : Constante de réaction [1/s]
        - Ce : Concentration en sel à la surface du poteau [mol/m^3]
        - K : Nombre de pas de temps [-]
        - T : Temps total de simulation [s]
    r : float
        Position radiale
    t : float
        Temps

    Returns
    -------
    C_anal : Vecteur (array) composée de la concentration en sel selon la 
             position [mol/L]

    '''

    # Définition de la discrétisation spatiale
    N = prm.N
    r = np.linspace(0, prm.R, N)

    # Définition de la discrétisation temporelle
    K = prm.K
    t = np.linspace(0, prm.T, K)

    r, t = np.meshgrid(r, t)

    return prm.C_e*(np.exp(1.0e-10*t) - 1)*np.cos(np.pi*r/(2*prm.R)) + prm.C_e

def source_MMS(prm, r, t):
    '''
    Fonction qui calcule le terme source selon la solution MMS

    Parameters
    ----------
    prm : Objet class parametres()
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]
        - k : Constante de réaction [1/s]
        - Ce : Concentration en sel à la surface du poteau [mol/m^3]
        - K : Nombre de pas de temps [-]
        - T : Temps total de simulation [s]
    r : float
        Position radiale à laquelle le terme source est calculé [m]
    t : float
        Temps à laquelle le terme source est calculé [s]

    Returns
    -------
    f : float
        Valeur du terme source à la position r et au temps t [mol/m^3/s]

    '''
    C_e = prm.C_e
    R = prm.R
    k = prm.k
    Deff = prm.D_eff
    exp = np.exp
    cos = np.cos
    sin = np.sin
    pi = np.pi
    f = 1.0e-10*C_e*exp(1.0e-10*t)*cos(pi*r/(2*R)) \
    - Deff*(-pi*C_e*(exp(1.0e-10*t) - 1)*sin(pi*r/(2*R))/(2*R) \
    - pi**2*C_e*r*(exp(1.0e-10*t) - 1)*cos(pi*r/(2*R))/(4*R**2))/r \
    + k*(C_e*(exp(1.0e-10*t) - 1)*cos(pi*r/(2*R)) + C_e)
    return f

def source_zero(prm, r, t):
    '''
    Fonction qui calcule le terme source nul

    Parameters
    ----------
    prm : Objet class parametres()
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]
        - k : Constante de réaction [1/s]
        - Ce : Concentration en sel à la surface du poteau [mol/m^3]
        - K : Nombre de pas de temps [-]
        - T : Temps total de simulation [s]
    r : float
        Position radiale à laquelle le terme source est calculé [m]
    t : float
        Temps à laquelle le terme source est calculé [s]

    Returns
    -------
    f : float
        Valeur du terme source à la position r et au temps t [mol/m^3/s]

    '''
    return 0
    

def concentration_mdf_o2(prm, source):
    '''
    Fonction qui simule, selon la méthode des différence finies avec un 
    schéma d'approximation d'ordre 2, la concentration de sel dans un pilier dans le temps
    

    Parameters
    ----------
    prm : Objet class parametres()
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]
        - k : Constante de réaction [1/s]
        - C_e : Concentration en sel à la surface du poteau [mol/m^3]
        - K : Nombre de pas de temps [-]
        - T : Temps total de simulation [s]
        - C_i : Concentration initiale [mol/m^3]
    
    source : Fonction
        Fonction qui calcule le terme source selon la solution

    Returns
    -------
    r : Vecteur (array) de dimension N composé de la position radiale à 
        laquelle les concentrations sont calculées, où N le nombre de noeuds.
    dr : float représentant de la pas de discrétisation spatial selon le rayon 
        du pilier
    C : Vecteur (array) composée de la concentration en sel selon la 
        position [mol/L]

    '''
    # Définition de la discrétisation spatiale
    N = prm.N
    r = np.linspace(0, prm.R, N)
    dr = prm.R/(N-1)

    # Définition de la discrétisation temporelle
    K = prm.K
    t = np.linspace(0, prm.T, K)
    dt = prm.T/(K-1)

    # Liste pour stocker les concentrations
    C_all = []
        
    #Initilisation de la matrice LHS
    A = np.zeros((N, N))
    # Conditions limites
    # Utilisation de l'approximation avant "Gear" de la dérivée 1ère
    A[0,0] = -3
    A[0,1] = 4
    A[0,2] = -1
    A[-1,-1] = 1
    #Initialisation du vecteur RHS
    B = np.zeros(N)
    B[-1] = prm.C_e #Condition limite en r=R

    C_t = prm.C_i*np.ones(N) # Condition initiale
    C_all.append(C_t)
    for k in range(1, K) :
        # print(f"itération {k}")
        #Remplir les coefficient de la matrice LHS et le vecteur RHS
        f_array = np.zeros(N-2)
        for i in range(1, N-1): #Exclure les premiers et derniers noeuds
            A[i,i-1] = prm.D_eff*dt*(1/dr**2 - 1/(r[i]*2*dr))
            A[i,i] = -2*prm.D_eff*dt/dr**2 - 1 - prm.k*dt
            A[i,i+1] = prm.D_eff*dt*(1/dr**2 + 1/(r[i]*2*dr))
            #Remplir les coefficients du vecteur B
            f = source(prm, r[i], t[k])
            f_array[i-1] = f
            B[i] = -C_t[i] - dt*f
        #Resolution dy systeme matriciel
        C = np.linalg.solve(A, B)

        #Sauvegarder la concentration
        C_all.append(C)

        #Mise à jour de la concentration
        C_t = C

    
    return r, dr, C_all