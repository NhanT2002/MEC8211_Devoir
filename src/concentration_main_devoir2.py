# Importation de modules
import numpy as np
import matplotlib.pyplot as plt
from concentration_fct import *
from erreur_fct import *
from graphique_convergence import *
#%% Assisgnation des paramètres
class parametres():
    C_e = 20         # Concentration en sel à la surface du poteau [mol/m^3]
    k = 4e-9         # Constante de réaction [1/s]
    D_eff = 1e-10    # Coefficient de diffusion effectif du sel dans le béton [1/s]
    R = 0.5          # Rayon du pilier [m]
    N = 11           # Nombre de noeuds [-]
    ti = 0           # Temps initial [s]
    tf = 4*10**9     # Temps final [s]
    dt = 1*365*24*3600 # Pas de temps (1 an) [s]

prm = parametres()

#%% Résolution du problème
r, dr, dt, temps, C = concentration_mdf_o2_implicite(prm)
temps = temps/prm.dt #en année

#%% Plot de la solution
# Plot des concentrations : solutions numérique en temps
fig = plt.figure()
ax = fig.add_subplot()
for i in range(int(temps[0]),int(temps[-1]),1): 
    ax.plot(r, C[:,i], '-o', label=f't = {temps[i]}')
ax.set_xlabel(r'Rayon du pilier $r$ [m]')
ax.set_ylabel(r'Concentration $C$ [mol/m$^3$]')
ax.set_title(f'Profil de concentration de sel dans le pilier, \n à l\'état transitoire (N = {prm.N})')
#plt.legend()
plt.grid(zorder=0)
plt.show()

# Représentation 2D: maillage temporel et spatial
ti, ri = np.meshgrid(temps, r, indexing='ij')

# Tracé des résultats
plt.figure()
plt.contourf(ri, ti, C.transpose(), levels=200)
plt.colorbar()
plt.title(f'Profil de concentration de sel dans le pilier, \n à l\'état transitoire (N = {prm.N})')
plt.xlabel('R [m]')
plt.ylabel('t [s]')
plt.show()