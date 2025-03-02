#Importation des modules
import numpy as np

#Définition des fonctions de concentration
def concentration_anal(prm):
    '''
    Fonction qui calcule, selon la solution analytique, 
    la concentration de sel dans un pilier

    Parameters
    ----------
    prm : Objet class parametres()
        - S : Terme source constant [mol/m^3/s]/
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]

    Returns
    -------
    C_anal : Vecteur (array) composée de la concentration en sel selon la 
             position [mol/L]

    '''
    r_range = np.linspace(0, prm.R, prm.N)
    # C_anal = 0.25*prm.S/prm.D_eff*prm.R*prm.R*(r_range*r_range/(prm.R*prm.R) - 1) + prm.C_e
    C_anal = 0.25*prm.S/prm.D_eff*prm.R**2*(r_range**2/(prm.R**2) - 1) + prm.C_e

    return C_anal

def concentration_mdf_o1(prm):
    '''
    Fonction qui simule, selon la méthode des différence finies avec un 
    schéma d'approximation d'ordre 1, la concentration de sel dans un pilier
    

    Parameters
    ----------
    prm : Objet class parametres()
        - S : Terme source constant [mol/m^3/s]
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]

    Returns
    -------
    r : Vecteur (array) de dimension N composé de la position radiale à 
        laquelle les concentrations sont calculées, où N le nombre de noeuds.
    dr : float représentant de la pas de discrétisation spatial selon le rayon 
        du pilier
    C : Vecteur (array) composée de la concentration en sel selon la 
        position [mol/L]

    '''
    # Définition de la discrétisation temporelle
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
    
    #Remplir les coefficient de la matrice LHS et le vecteur RHS
    for i in range(1, N-1): #Exclure les premiers et derniers noeuds
        #Remplir les coefficients de la matrice A
        A[i,i-1] = 1/(dr*dr)
        A[i,i] = -(2/(dr*dr) + 1/(r[i]*dr))
        A[i,i+1] = 1/(dr*dr) + 1/(r[i]*dr)
        #Remplir les coefficients du vecteur B
        B[i] = prm.S/prm.D_eff

    B[-1] = prm.C_e #Condition limite en r=R
    
    #Résolution dy système matriciel
    C = np.linalg.solve(A, B)
    
    return r, dr, C

def concentration_mdf_o2(prm):
    '''
    Fonction qui simule, selon la méthode des différence finies avec un 
    schéma d'approximation d'ordre 2, la concentration de sel dans un pilier
    

    Parameters
    ----------
    prm : Objet class parametres()
        - S : Terme source constant [mol/m^3/s]
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]

    Returns
    -------
    r : Vecteur (array) de dimension N composé de la position radiale à 
        laquelle les concentrations sont calculées, où N le nombre de noeuds.
    dr : float représentant de la pas de discrétisation spatial selon le rayon 
        du pilier
    C : Vecteur (array) composée de la concentration en sel selon la 
        position [mol/L]

    '''
    # Définition de la discrétisation temporelle
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
    
    #Remplir les coefficient de la matrice LHS et le vecteur RHS
    for i in range(1, N-1): #Exclure les premiers et derniers noeuds
        A[i,i-1] = 1/(dr*dr) - 1/(r[i]*2*dr)
        A[i,i] = -2/(dr*dr)
        A[i,i+1] = 1/(dr*dr) + 1/(r[i]*2*dr)
        #Remplir les coefficients du vecteur B
        B[i] = prm.S/prm.D_eff
        
    B[-1] = prm.C_e #Condition limite en r=R
    #Resolution dy systeme matriciel
    C = np.linalg.solve(A, B)
    
    return r, dr, C