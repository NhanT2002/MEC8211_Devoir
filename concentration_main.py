# Importation de modules
import numpy as np
import matplotlib.pyplot as plt
from concentration_fct import *
#%% Assisgnation des paramètres
class parametres():
    C_e = 20         # Concentration en sel à la surface du poteau [mol/m^3]
    S = 2e-8         # Terme source constant [mol/m^3/s]
    D_eff = 1e-10   # Coefficient de diffusion effectif du sel dans le béton [1/s]
    R = 0.5          # Rayon du pilier [m]
    N = 20           # Nombre de noeuds [-]

prm = parametres()

#%% Résolution du problème
C_anal = concentration_anal(prm)
r, dr, C_mdf_o1 = concentration_mdf_o1(prm)
r, dr, C_mdf_o2 = concentration_mdf_o2(prm)

# Plot des concentrations
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(r, C_mdf_o1, '-o', label='MDF: Schéma d\'ordre 1')
ax.plot(r, C_mdf_o2, '-o', label='MDF: Schéma d\'ordre 2')
ax.plot(r, C_anal, 'k--', lw = 1, label='Solution analytique')
ax.set_xlabel(r'Rayon du pilier $r$ [m]')
ax.set_ylabel(r'Concentration $C$ [mol/m$^3$]')
ax.set_title('Profil de concentration de sel dans le pilier, \n à l\'état stationnaire')
plt.legend()
plt.grid(zorder=0)
plt.show()

#%% Analyse des erreurs
# Mailles differentes
N_range = np.arange(5,105,10)

#Initialisation pour le schema d'ordre 1
L1_range_o1 = []
L2_range_o1 = []
L_inf_range_o1 = []

#Initialisation pour le schema d'ordre 2
L1_range_o2 = []
L2_range_o2 = []
L_inf_range_o2 = []
dr_range = []

for i in range(len(N_range)):
    prm.N = N_range[i]
    
    #Calcul de la solution selon le maillage
    C_anal = concentration_anal(prm)
    r, dr, C_mdf_o1 = concentration_mdf_o1(prm)
    r, dr, C_mdf_o2 = concentration_mdf_o2(prm)
    dr_range.append(dr)
    
    #Calcul des normes des erreurs - schema d'ordre 1
    L1_o1 = erreur_L1(prm, C_mdf_o1, C_anal)
    L2_o1 = erreur_L2(prm, C_mdf_o1, C_anal)
    L_inf_o1 = erreur_inf(C_mdf_o1, C_anal)
    
    L1_range_o1.append(L1_o1)
    L2_range_o1.append(L2_o1)
    L_inf_range_o1.append(L_inf_o1)
    
    #Calcul des normes des erreurs - schema d'ordre 2
    L1_o2 = erreur_L1(prm, C_mdf_o2, C_anal)
    L2_o2 = erreur_L2(prm, C_mdf_o2, C_anal)
    L_inf_o2 = erreur_inf(C_mdf_o2, C_anal)
        
    L1_range_o2.append(L1_o2)
    L2_range_o2.append(L2_o2)
    L_inf_range_o2.append(L_inf_o2)
    
#%% Plot des graphiques

#Plot de l'ordre de convergence
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(dr_range, L1_range_o1, '-o', label='$L_1$, données numériques')
ax.plot(dr_range, L2_range_o1, '-o', label='$L_2$, données numériques')
ax.plot(dr_range, L_inf_range_o1, 'k--', lw = 1, label=r'$L_{\inf}$, données numériques')
# Reste à rajouter les polyfit.

ax.set_xlabel(r'Taille de maille $\Delta r$ [m]')
ax.set_xscale('log')
ax.set_ylabel(r'Erreur ||$u-u_{ref}$|| [mol/m$^3$]')
ax.set_yscale('log')
ax.set_title('Convergence d\'ordre XXX \n des erreurs $L_1$, $L_2$ et $L_{inf}$ en fonction de $\Delta r$')
plt.legend()
plt.grid(zorder=0)
plt.show()

#Plot de l'ordre de convergence
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(dr_range, L1_range_o2, '-o', label='$L1$')
ax.plot(dr_range, L2_range_o2, '-o', label='$L2$')
ax.plot(dr_range, L_inf_range_o2, 'k--', lw = 1, label=r'$L_{\inf}$')
# Reste à rajouter les polyfit.

ax.set_xlabel(r'Taille de maille $\Delta r$ [m]')
ax.set_xscale('log')
ax.set_ylabel(r'Erreur ||$u-u_{ref}$|| [mol/m$^3$]')
ax.set_yscale('log')
ax.set_title('Convergence d\'ordre XXX \n des erreurs $L_1$, $L_2$ et $L_{inf}$ en fonction de $\Delta r$')
plt.legend()
plt.grid(zorder=0)
plt.show()

#Plot des erreurs