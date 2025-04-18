#%%Importation des modules
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import os

#%% Importation des données
#Seeds generated from LBM from MATLAB
seeds = [101, 102, 103, 104, 105]

#Initialisation of vector ranges for data
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
#%% Mean values of the parameters taken from the different seeds
dx_mean = np.mean(dx, axis=0)
n_raf = np.arange(len(dx_mean)) + 1
nx = 100*n_raf
poro_mean = np.mean(poro, axis=0)
Re_mean = np.mean(Re, axis=0)
k_mean = np.mean(k, axis=0)

#%% Vérification de la convergence
# Vérifier si l'on est dans une région asymptotique
error = np.abs(k_mean[-1]-k_mean)/k_mean[-1]
order, intercept = np.polyfit(np.log(dx_mean[:-1]), np.log(error[:-1]), 1)

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
plt.xlabel(f'$\Delta x$ [m]')
plt.ylabel(r'$|\frac{k_{finest}-k_{\Delta x}}{k_{finest}}|$ [-]')
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
plt.savefig("convergence_k.svg", dpi=300, bbox_inches='tight')

#%% Vérification de solution, après avoir vu que l'on converge vers une solution
# Ordre de convergence obeservé
# pour NX = 100;200;400 pour avoir le ratio de r = 2
p_hat_a = np.log((k_mean[0]-k_mean[1])/(k_mean[1]-k_mean[3]))/np.log(2)

# pour NX = 200;400;800 pour avoir le ratio de r = 2
p_hat_b = np.log((k_mean[1]-k_mean[3])/(k_mean[3]-k_mean[7]))/np.log(2)

# Calcul du GCI selon la valeur obtenue de p_hat_b, soit celui avec le maillage plus raffiné
FS = 3.0 #Facteur de sécurité
r = 2.0 # ratio de raffinement de maillage
# on voudra faire l'analyse de propagation avec NX = 200, donc on reste consitent
f2 = k_mean[0] # NX = 100, maillage le moins fin
f1 = k_mean[1] # NX = 200, maillage le plus fin
GCI = FS/(r**p_hat_b-1)*np.abs(f2-f1)
#%% --------------------------------------------------- u_input ---------------------------------------------------
# Load data from the random seed generated with LBM in MATLAB
data = pd.read_csv("results_seed_0.csv")

# Plot the different statistical graphs (PFD, CFD) with their respective parameters
# plot histogram of the diameter
plt.figure()
plt.hist(data["d_equivalent"], bins=30, color='tab:blue', alpha=0.7, density=True, edgecolor='black')
plt.xlabel(r'$d_{eq}$ [$\mu m$]')
plt.ylabel('PDF')
plt.title('Histogram of Equivalent Diameter')
mu, std = stats.norm.fit(data["d_equivalent"])
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = stats.norm.pdf(x, mu, std)
plt.plot(x, p, color="tab:orange", linewidth=2, label='Normal fit, $\mu$={:.3f}, $\sigma$={:.3f}'.format(mu, std))
plt.legend()
plt.grid()
plt.savefig("histogram_d_equivalent.svg", dpi=300, bbox_inches='tight')

# plot histogram of the porosity
plt.figure()
plt.hist(data["poro_eff"], bins=30, color='tab:blue', alpha=0.7, density=True, edgecolor='black')
plt.xlabel('Porosity')
plt.ylabel('PDF')
plt.title('Histogram of Porosity')
mu_poro, std_poro = stats.norm.fit(data["poro_eff"]) 
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = stats.norm.pdf(x, mu_poro, std_poro)
plt.plot(x, p, color="tab:orange", linewidth=2, label='Normal fit, $\mu$={:.3f}, $\sigma$={:.3e}'.format(mu, std))
plt.legend()
plt.grid()
plt.savefig("histogram_poro_eff.svg", dpi=300, bbox_inches='tight')

# plot cummulative distribution function of porosity
plt.figure()
plt.hist(data["poro_eff"], bins=30, color='tab:blue', alpha=0.7, density=True, cumulative=True, edgecolor='black')
plt.xlabel('Porosity')
plt.ylabel('CDF')
plt.title('Cumulative Distribution Function of Porosity')
cdf = stats.norm.cdf(x, mu, std)
plt.plot(x, cdf, color="tab:orange", linewidth=2)
plt.grid()
plt.savefig("cdf_poro_eff.svg", dpi=300, bbox_inches='tight')

# plot histogram of the permeability
plt.figure()
plt.hist(data["k_in_micron2"], bins=30, color='tab:blue', alpha=0.7, density=True, edgecolor='black')
plt.xlabel(r'$k$ [$\mu m^2$]')
plt.ylabel('PDF')
plt.title('Histogram of Permeability')
mu, std = stats.norm.fit(np.log(data["k_in_micron2"])) # Fit a log-normal distribution to the data
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
lognorm_pdf = stats.lognorm.pdf(x, s=std, scale=np.exp(mu))  # Log-normal PDF
plt.plot(x, lognorm_pdf, color="tab:orange", linewidth=2, label='Log-normal fit, $\mu$={:.3f}, $\sigma$={:.3f}'.format(mu, std))
plt.legend()
plt.grid()
log_mean = np.exp(mu + std**2 / 2)
log_stddev = log_mean * np.sqrt(np.exp(std**2) - 1)
plt.text(0.45, 0.7, f'Median = $e^\mu = $ {np.exp(mu):.3f} $\mu m^2$', fontsize=12, transform=plt.gca().transAxes, color='k')
plt.text(0.45, 0.65, f'FVG = $e^\sigma = $ {np.exp(std):.3f} $\mu m^2$', fontsize=12, transform=plt.gca().transAxes, color='k')
plt.text(0.45, 0.55, r'$E[X] = e^{\mu+\frac{1}{2}\sigma^2}=$' + f'{log_mean:.4f}', fontsize=12, transform=plt.gca().transAxes, color='k')
plt.text(0.45, 0.45, r'$SD[X] = E[X]\sqrt{e^{\sigma^2}-1}=$' + f'{log_stddev:.4f}', fontsize=12, transform=plt.gca().transAxes, color='k')
plt.savefig("histogram_k.svg", dpi=300, bbox_inches='tight')

# plot cummulative distribution function of permeability
plt.figure()
plt.hist(data["k_in_micron2"], bins=30, color='tab:blue', alpha=0.7, density=True, cumulative=True, edgecolor='black')
plt.xlabel(r'$k$ [$\mu m^2$]')
plt.ylabel('CDF')
plt.title('Cumulative Distribution Function of Permeability')
lognorm_cdf = stats.lognorm.cdf(x, s=std, scale=np.exp(mu))
plt.plot(x, lognorm_cdf, color="tab:orange", linewidth=2)
plt.grid()
plt.savefig("cdf_k.svg", dpi=300, bbox_inches='tight')

# Calcul des intervalles autour de la médiane avec le FVG
Median = np.exp(mu)
FVG = np.exp(std)
FVG_lower = Median/FVG
FVG_upper = Median*FVG
#%% --------------------------------------------- Erreur du modèle ---------------------------------------------
# Incertittudes numériques
u_num = GCI/2
# Incertitude des données d'entrées
u_input_lower = Median - FVG_lower
u_input_upper = FVG_upper - Median
# Incertitude des données expérimentales
u_d = np.sqrt(14.7**2 + 10.0**2)
# Incertitude du modèle
u_val_upper = np.sqrt(u_num**2 + u_input_upper**2 + u_d**2)
u_val_lower = np.sqrt(u_num**2 + u_input_lower**2 + u_d**2)
# Erreur du modèle
S = Median
D = 80.6 # [10e-3 m^2]
E = S - D 

# Intervalles de confiance
k_confiance = 2. #Car on veut à 95.4% de confiance
U_95_4_lower = k_confiance * u_val_lower
U_95_4_upper = k_confiance * u_val_upper

# Vraie erreur du modèle bornée
delta_model_lower = E - U_95_4_lower
delta_model_upper = E + U_95_4_upper

# Vraie erreur du modèle sans u_input
u_val_no_input = np.sqrt(u_num**2 + u_d**2)
delta_model_lower_no_input = E - k_confiance * u_val_no_input
delta_model_upper_no_input = E + k_confiance * u_val_no_input

# plot distance de l'erreur (avec résultats de 200x200)
fig=plt.figure()
ax=fig.add_subplot(111)
plt.xlabel(r'Porosity [-]')
plt.ylabel(r'$E = S-D$ [$\mu m^2$]')
plt.title('Relative model error for mesh of 200x200')
ax.set_xlim(xmin = 0.8, xmax = 1)
ax.set_ylim(ymin = -100, ymax = 10)
ax.plot(mu_poro, E, 'o')
ax.hlines(0, 0.8, 1.0, 'k')
ax.errorbar(mu_poro, E, yerr = ([U_95_4_lower], [U_95_4_upper]), xerr = std_poro, ecolor = "tab:orange", capsize = 5, lw = 3, label = r'$u_{input} \neq 0$', zorder = 0)
ax.errorbar(mu_poro, E, yerr = k_confiance*u_val_no_input, xerr = std_poro, ecolor = "tab:green", capsize = 5, label = r'$u_{input} = 0$', zorder = 0)
ax.legend(loc='lower right')
plt.savefig("relative_model_error.svg", dpi=300, bbox_inches='tight')
