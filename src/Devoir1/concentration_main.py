# Importation de modules
import numpy as np
import matplotlib.pyplot as plt
from concentration_fct import *
from erreur_fct import *
from graphique_convergence import *
#%% Assisgnation des paramètres
class parametres():
    C_e = 20         # Concentration en sel à la surface du poteau [mol/m^3]
    S = 2e-8         # Terme source constant [mol/m^3/s]
    D_eff = 1e-10    # Coefficient de diffusion effectif du sel dans le béton [1/s]
    R = 0.5          # Rayon du pilier [m]
    N = 20           # Nombre de noeuds [-]

prm = parametres()

#%% Résolution du problème
C_anal = concentration_anal(prm) # solution analytique
r, dr, C_mdf_o1 = concentration_mdf_o1(prm) # solution numérique, schémas d'ordre 1
r, dr, C_mdf_o2 = concentration_mdf_o2(prm) # solution numérique, schémas d'ordre 2

# Plot des concentrations : solutions analytique et numérique
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(r, C_mdf_o1, 'o', label='MDF: Schéma d\'ordre 1')
ax.plot(r, C_mdf_o2, 'o', label='MDF: Schéma d\'ordre 2')
ax.plot(r, C_anal, 'k--', lw = 1, label='Solution analytique')
ax.set_xlabel(r'Rayon du pilier $r$ [m]')
ax.set_ylabel(r'Concentration $C$ [mol/m$^3$]')
ax.set_title(f'Profil de concentration de sel dans le pilier, \n à l\'état stationnaire (N = {prm.N})')
plt.legend()
plt.grid(zorder=0)
plt.show()


#%%Vérification de code
# Plot des ordres de convergence observés
graphique_convergence_ordre(prm, concentration_anal, concentration_mdf_o1, 1) #ordre 1
graphique_convergence_ordre(prm, concentration_anal, concentration_mdf_o2, 2) #ordre 2

# Test de symétrie
prm.N = 20
prm.R = -prm.R
C_anal_sym = concentration_anal(prm) # solution analytique
r_sym, dr_sym, C_mdf_o1_sym = concentration_mdf_o1(prm) # solution numérique, schémas d'ordre 1
r_sym, dr_sym, C_mdf_o2_sym, = concentration_mdf_o2(prm) # solution numérique, schémas d'ordre 2

# Plot des concentrations : solutions analytique et numérique
plt.figure()
plt.plot(r, C_mdf_o1, "o", label=f'Solution numérique')
plt.plot(r_sym, C_mdf_o1_sym, "o", label=f'Solution numérique (test symétrie)')
plt.plot(r, C_anal, color="tab:green", label='Solution analytique')
plt.plot(r_sym, C_anal_sym, color="tab:green")
plt.xlabel('r [m]')
plt.ylabel('C [mol/$m^3$]')
plt.title(f"Profil de concentration à l'état stationnaire  (N={prm.N})\n Ordre 1")
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.plot(r, C_mdf_o2, "o", label=f'Solution numérique')
plt.plot(r_sym, C_mdf_o2_sym, "o", label=f'Solution numérique (test symétrie)')
plt.plot(r, C_anal, color="tab:green", label='Solution analytique')
plt.plot(r_sym, C_anal_sym, color="tab:green")
plt.xlabel('r [m]')
plt.ylabel('C [mol/$m^3$]')
plt.title(f"Profil de concentration à l'état stationnaire  (N={prm.N})\n Ordre 2")
plt.legend()
plt.grid()
plt.show()