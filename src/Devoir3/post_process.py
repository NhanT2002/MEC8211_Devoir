import numpy as np
from scipy import stats
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
n_raf = np.arange(len(dx_mean)) + 1
nx = 100*n_raf
poro_mean = np.mean(poro, axis=0)
Re_mean = np.mean(Re, axis=0)
k_mean = np.mean(k, axis=0)

p_hat = np.log((k_mean[0]-k_mean[1])/(k_mean[1]-k_mean[3]))/np.log(2)
p_hat = np.log((k_mean[1]-k_mean[3])/(k_mean[3]-k_mean[7]))/np.log(2)
GCI = 3.0/(2**p_hat-1)*np.abs(k_mean[1]-k_mean[0])

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
plt.xlabel(f'$\Delta x$ [m]')
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
plt.savefig("convergence_k.svg", dpi=300, bbox_inches='tight')


#%% --------------------------------------------------- u_input ---------------------------------------------------

data = pd.read_csv("results_seed_0.csv")


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


#%% --------------------------------------------- Error of model ---------------------------------------------
u_num = GCI/2
u_input = log_stddev
u_d = np.sqrt(14.7**2 + 10.0**2)
u_val = np.sqrt(u_num**2 + u_input**2 + u_d**2)
E = np.exp(mu) - 80.6
delta_model_lower = E - 2. * u_val
delta_model_upper = E + 2. * u_val

# plot distance of error (with 200x200 results)
fig=plt.figure()
ax=fig.add_subplot(111)
plt.xlabel(r'Porosity [-]')
plt.ylabel(r'$E = S-D$ [$\mu m^2$]')
plt.title('Realtive model error')
ax.set_xlim(xmin = 0.8, xmax = 1)
ax.set_ylim(ymin = -100, ymax = 10)
ax.plot(mu_poro, E, 'o')
ax.hlines(0, 0.8, 1.0, 'k')
ax.errorbar(mu_poro, E, yerr = 2*u_val, xerr = std_poro, ecolor = "tab:orange", capsize = 5, zorder = 0)

