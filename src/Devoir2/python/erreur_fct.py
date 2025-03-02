#Importation des modules
import numpy as np

#Définition des fonctions d'erreurs
def erreur_L1(prm, u, u_ref):
    '''
    Fonction qui calcule la norme L_1 discrète de l'erreur associée à une 
    résolution numérique par la méthode de différence finies (MDF).

    Parameters
    ----------
    prm : Objet class parametres()
        - S : Terme source constant [mol/m^3/s]/
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]
    u : Vecteur (array) de la solution numérique par la méthode MDF
    u_ref : Vecteur (array) de la solution de référence, soit celle analytique

    Returns
    -------
    L1 : float64 de la norme L_1 discrète de l'erreur associée à une  
    résolution numérique par la MDF.

    
    '''
    L1 = 1/(prm.N*prm.K)*np.sum(np.abs(u-u_ref))
    return L1

def erreur_L2(prm, u, u_ref):
    '''
    Fonction qui calcule la norme L_2 discrète de l'erreur associée à une 
    résolution numérique par la méthode de différence finies (MDF).

    Parameters
    ----------
    prm : Objet class parametres()
        - S : Terme source constant [mol/m^3/s]/
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]
    u : Vecteur (array) de la solution numérique par la méthode MDF
    u_ref : Vecteur (array) de la solution de référence, soit celle analytique

    Returns
    -------
    L2 : float64 de la norme L_2 discrète de l'erreur associée à une  
    résolution numérique par la MDF.

    
    '''
    L2 = np.sqrt(1/(prm.N*prm.K)*np.sum(np.abs(u-u_ref)**2))
    return L2

def erreur_inf(u, u_ref):
    '''
    Fonction qui calcule la norme L_infinie discrète de l'erreur associée à une 
    résolution numérique par la méthode de différence finies (MDF).

    Parameters
    ----------
    u : Vecteur (array) de la solution numérique par la méthode MDF
    u_ref : Vecteur (array) de la solution de référence, soit celle analytique

    Returns
    -------
    L_inf: float64 de la norme L_infinie discrète de l'erreur associée à une  
    résolution numérique par la MDF.

    
    '''
    L_inf = np.max(np.abs(u-u_ref))
    return L_inf