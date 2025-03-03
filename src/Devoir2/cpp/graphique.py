import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def plot(filename) :
    data = pd.read_csv(filename)
    r = np.float64(np.array(data.columns)[1:])
    plt.figure()
    for i in np.linspace(0, len(data)-1, 10, dtype=int):
        plt.plot(r, np.array(data.iloc[i, 1:]), label=f't = {data.iloc[i, 0]:.2e}')
    plt.xlabel(r'Rayon du pilier $r$ [m]')
    plt.ylabel(r'Concentration $C$ [mol/m$^3$]')
    plt.title(f'Profil de concentration de sel dans le pilier (N = {len(r)})')
    plt.legend()
    plt.grid(zorder=0)
    plt.show()

plot("output_MMS.csv")
plot("output.csv")


data_space = pd.read_csv("convergence_results_space.csv")
data_time = pd.read_csv("convergence_results_time.csv")

def plot_convergence(data, title, space=True):
    L1_range = np.array(data['L1_error'])
    L2_range = np.array(data['L2_error'])
    L_inf_range = np.array(data['L_inf_error'])
    dr_range = 0.5/(np.array(data.iloc[:,0])-1)
    # Loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
    if space:
        order_L1, intercept_L1 = np.polyfit(np.log(dr_range[1:5]), np.log(L1_range[1:5]), 1)
        order_L2, intercept_L2 = np.polyfit(np.log(dr_range[1:5]), np.log(L2_range[1:5]), 1)
        order_Linf, intercept_Linf = np.polyfit(np.log(dr_range[1:5]), np.log(L_inf_range[1:5]), 1)
    else:
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
    if space:
        plt.xlabel(f'$\Delta r$ [m]')
    else:
        plt.xlabel(f'$\Delta t$ [s]')
    plt.ylabel('Erreur [mol/$m^3$]')
    plt.title(title)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()

    # Rendre les axes plus gras
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)
    plt.tick_params(width=2, which='both', direction='in', top=True, right=True, length=6)

delta_t = 4e9/(10000-1)
title = f"Convergence de l'erreur en fonction de $\Delta r$ pour $\Delta t$ = {delta_t:.2e} s"
plot_convergence(data_space, title)

delta_r = 0.5/(1000-1)
title = f"Convergence de l'erreur en fonction de $\Delta t$ pour $\Delta r$ = {delta_r:.2e} m"
plot_convergence(data_time, title, space=False)