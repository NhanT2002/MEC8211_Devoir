# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 00:46:25 2025

@author: dorid
"""
import numpy as np
import matplotlib.pyplot as plt

def erreur_L1(prm, u, u_ref):
    L1 = 1/prm.N*sum(abs(u-u_ref))
    return L1

def erreur_L2(prm, u, u_ref):
    L2 = np.sqrt(1/prm.N*sum(abs(u-u_ref)**2))
    return L2

def erreur_inf(u, u_ref):
    L_inf = max(abs(u-u_ref))
    return L_inf

def erreur_ordre_observe(h_range, Li_range):
    # Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
    coefficients = np.polyfit(np.log(h_range[:]), np.log(Li_range[:]), 1)
    exponent = coefficients[0]
    
    # Fonction de régression en termes de logarithmes
    fit_function_log = lambda x: exponent * x + coefficients[1]

    # Fonction de régression en termes originaux
    fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

    # Extrapoler la valeur prédite pour la dernière valeur de h_range
    extrapolated_value = fit_function(h_range[-1])
    
    # Tracer le graphique en échelle log-log avec des points et la courbe de régression extrapolée
    plt.figure(figsize=(8, 6))
    plt.scatter(h_range, Li_range, marker='o', color='b', label='Données numériques obtenues')
    plt.plot(h_range, fit_function(h_range), linestyle='--', color='r', label='Régression en loi de puissance')

    # Marquer la valeur extrapolée
    #plt.scatter(h_range[-1], extrapolated_value, marker='x', color='g', label='Extrapolation')

    # Ajouter des étiquettes et un titre au graphique
    plt.title('Convergence d\'ordre 2\n de l\'erreur $L_2$ en fonction de $Δx$',
              fontsize=14, fontweight='bold', y=1.02)  # Le paramètre y règle la position verticale du titre

    plt.xlabel('Taille de maille $h_{max}$ ou $Δx$ (cm)', fontsize=12, fontweight='bold')  # Remplacer "h" par "Δx"
    plt.ylabel('Erreur $L_2$ (m/s)', fontsize=12, fontweight='bold')

    # Rendre les axes plus gras
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)

    # Placer les marques de coche à l'intérieur et les rendre un peu plus longues
    plt.tick_params(width=2, which='both', direction='in', top=True, right=True, length=6)

    # Afficher l'équation de la régression en loi de puissance
    equation_text = f'$L_2 = {np.exp(coefficients[1]):.4f} \\times Δx^{{{exponent:.4f}}}$'
    equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

    # Déplacer la zone de texte
    equation_text_obj.set_position((0.5, 0.4))


    # Afficher le graphique
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.show()
    
    return
