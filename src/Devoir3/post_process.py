import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

seeds = [101, 102, 103, 104, 105]

dx = []
poro = []
Re = []
k = []

# Load the data
for seed in seeds:
    data = pd.read_csv(f"results_seed_{seed}.csv")
    dx.append(data["dx_i"])
    poro.append(data["poro_eff"])
    Re.append(data["Re"])
    k.append(data["k_in_micron2"])

dx_mean = np.mean(dx, axis=0)
poro_mean = np.mean(poro, axis=0)
Re_mean = np.mean(Re, axis=0)
k_mean = np.mean(k, axis=0)

p_hat = np.log((k_mean[1]-k_mean[3])/(k_mean[3]-k_mean[7]))/np.log(2)

error = np.abs(k_mean[-1]-k_mean)/k_mean[-1]
order, intercept = np.polyfit(np.log(dx_mean[4:-1]), np.log(error[4:-1]), 1)

ref_y = np.exp(intercept) * dx_mean**order

# Tracer le graphique en échelle log-log avec des points et la courbe de régression
plt.figure()
plt.plot(dx_mean, error, 'o', label=r'$|\frac{k_{finest}-k_{\Delta x}}{k_{finest}}|$ numérique', color='tab:blue')
plt.plot(dx_mean, ref_y, "--", color='tab:blue')
plt.plot([], [], "--", color='black', label='Régression en loi de puissance')

# Afficher l'équation de la régression en loi de puissance
equation_text = r'$|\frac{k_{finest}-k_{\Delta x}}{k_{finest}}|$' + f'$= {intercept:.4f} \\times Δx^{{{order:.4f}}} $'
equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')
equation_text_obj.set_position((0.4, 0.2))

# Ajouter des étiquettes et un titre au graphique
plt.legend()
plt.xlabel(f'$\Delta x$ [s]')
plt.ylabel(r'$|\frac{k_{finest}-k_{\Delta x}}{k_{finest}}|$ [$\mu m^2$]')
plt.title('Convergence of the permeability $k$')
plt.xscale('log')
plt.yscale('log')
plt.grid()

# Rendre les axes plus gras
plt.gca().spines['bottom'].set_linewidth(2)
plt.gca().spines['left'].set_linewidth(2)
plt.gca().spines['right'].set_linewidth(2)
plt.gca().spines['top'].set_linewidth(2)
plt.tick_params(width=2, which='both', direction='in', top=True, right=True, length=6)







