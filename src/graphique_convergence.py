#Importation des modules
import numpy as np
import matplotlib.pyplot as plt
from erreur_fct import *
from concentration_fct import *

#Définition de la fonction de graphique de convergence
def graphique_convergence_ordre(prm, sol_anal, mdf_oi, ordre_schema):
    '''
    Fonction qui trace le graphique de l'ordre de covergence de l'erreur 
    observée selon la norme L_1, L_2 et L_infinie pour la résolution d'une 
    solution par la  méthode des différences finies (MDF), 
    pour une taille d'élément donnée.

    Parameters
    ----------
    prm : Objet class parametres()
        - S : Terme source constant [mol/m^3/s]
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]
    sol_anal : fonction qui calcule la solution analytique 
    mdf_oi : fonction qui simule la solution numérique par MDF d'ordre "i"
            spécifiée.
    ordre_schema : float64 qui définit l'ordre de convergence du schéma formel

    Returns
    -------
    None.

    '''
    # Mailles differentes
    N_range = np.logspace(0.5, 3, 10, dtype=int)
    #Initialisation pour le schema d'ordre spécifié
    L1_range = []
    L2_range = []
    L_inf_range = []
    dr_range = []
    for n in N_range:
        prm.N = n
        #Calcul de la solution selon le maillage et l'ordre de la MDF
        C_anal = sol_anal(prm)
        r, dr, C_mdf_oi = mdf_oi(prm)
        dr_range.append(dr)
        
        #Calcul des normes des erreurs selon le schema d'ordre spécifié
        L1_oi = erreur_L1(prm, C_mdf_oi, C_anal)
        L2_oi = erreur_L2(prm, C_mdf_oi, C_anal)
        L_inf_oi = erreur_inf(C_mdf_oi, C_anal)
        
        L1_range.append(L1_oi)
        L2_range.append(L2_oi)
        L_inf_range.append(L_inf_oi)
    
    # Loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
    order_L1, intercept_L1 = np.polyfit(np.log(dr_range), np.log(L1_range), 1)
    order_L2, intercept_L2 = np.polyfit(np.log(dr_range), np.log(L2_range), 1)
    order_Linf, intercept_Linf = np.polyfit(np.log(dr_range), np.log(L_inf_range), 1)
    
    C_L1 = np.exp(intercept_L1)
    C_L2 = np.exp(intercept_L2)
    C_Linf = np.exp(intercept_Linf)

    # Plot computed reference lines
    ref_x = np.array([dr_range[0], dr_range[-1]])
    ref_y1 = C_L1 * ref_x**order_L1  # Computed reference line for L1
    ref_y2 = C_L2 * ref_x**order_L2  # Computed reference line for L2
    ref_yinf = C_Linf * ref_x**order_Linf  # Computed reference line for Linf
    
    # Tracer le graphique en échelle log-log avec des points et la courbe de régression
    plt.figure()
    plt.plot(dr_range, L1_range, 'o', label=f'$L_1$ numérique', color='tab:blue')
    plt.plot(dr_range, L2_range, 'o', label=f'$L_2$ numérique', color='tab:orange')
    plt.plot(dr_range, L_inf_range, 'o', label=f'$L_\infty$ numérique', color='tab:green')

    plt.plot(ref_x, ref_y1, "--", color='tab:blue')
    plt.plot(ref_x, ref_y2, "--", color='tab:orange')
    plt.plot(ref_x, ref_yinf, "--", color='tab:green')
    plt.plot([], [], "--", color='black', label='Régression en loi de puissance')

    # Afficher l'équation de la régression en loi de puissance
    equation_text = f'$L_1 = {C_L1:.4f} \\times Δr^{{{order_L1:.4f}}}$\n$L_2 = {C_L2:.4f} \\times Δr^{{{order_L2:.4f}}}$\n$L_\infty = {C_Linf:.4f} \\times Δr^{{{order_Linf:.4f}}}$'
    equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')
    equation_text_obj.set_position((0.5, 0.2))
    
    # Ajouter des étiquettes et un titre au graphique
    plt.legend()
    plt.xlabel(f'$\Delta r$ [m]')
    plt.ylabel('Erreur [mol/$m^3$]')
    plt.title(f"Convergence de l'erreur en fonction de $\Delta r$ pour le schema d'ordre {ordre_schema}")
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    
    # Rendre les axes plus gras
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)
    plt.tick_params(width=2, which='both', direction='in', top=True, right=True, length=6)
    
    return
