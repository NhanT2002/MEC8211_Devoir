import numpy as np
import matplotlib.pyplot as plt

Ce = 20 # [mol/m3]
D = 1 # [m]
R = D/2 # [m]
D_eff = 1E-10 # [m2/s]
k = 4E-9 # [s^-1]
S = 2E-8 # [mol/m3/s]

def C_order_1(N) :
    r = np.linspace(0, R, N) # domain
    dr = r[1] - r[0] # grid spacing

    A = np.zeros((N, N))
    B = np.zeros(N)

    for i in range(N) :
        if i == 0 :
            A[i, i] = -1
            A[i, i+1] = 1
            B[i] = 0
        elif i == N-1 :
            A[i, i] = 1
            B[i] = Ce
        else :
            A[i, i-1] = 1/dr**2
            A[i, i] = -2/dr**2 - 1/(dr*r[i])
            A[i, i+1] = 1/dr**2 + 1/(dr*r[i])
            B[i] = S/D_eff

    C = np.linalg.solve(A, B)
    C_ana = S/(4*D_eff)*R**2*(r**2/R**2 -1) + Ce
    return C, C_ana, r

def C_order_2(N) :
    r = np.linspace(0, R, N) # domain
    dr = r[1] - r[0] # grid spacing

    A = np.zeros((N, N))
    B = np.zeros(N)

    for i in range(N) :
        if i == 0 :
            A[i, i] = -3
            A[i, i+1] = 4
            A[i, i+2] = -1
            B[i] = 0
        elif i == N-1 :
            A[i, i] = 1
            B[i] = Ce
        else :
            A[i, i-1] = 1/dr**2 - 1/(2*dr*r[i])
            A[i, i] = -2/dr**2
            A[i, i+1] = 1/dr**2 + 1/(2*dr*r[i])
            B[i] = S/D_eff

    C = np.linalg.solve(A, B)
    C_ana = S/(4*D_eff)*R**2*(r**2/R**2 -1) + Ce
    return C, C_ana, r

N = 20 # number of grid points
C_o1, C_ana, r = C_order_1(N)
C_o2, C_ana, r = C_order_2(N)

plt.figure()
plt.plot(r, C_o1, "o", label=f'Solution numérique, schéma d\'ordre 1')
plt.plot(r, C_o2, "o", label=f'Solution numérique, schéma d\'ordre 2')
plt.plot(r, C_ana, label='Solution analytique')
plt.xlabel('r [m]')
plt.ylabel('C [mol/$m^3$]')
plt.title(f"Profil de concentration à l'état stationnaire  (N={N})")
plt.legend()
plt.grid()
plt.show()

def convergence_plot(scheme_order) :
    L1_array = []
    L2_array = []
    L_inf_array = []
    dr_array = []
    N_range = np.logspace(0.5, 3, 10, dtype=int)
    for n in N_range :
        C, C_ana, r = scheme_order(n)
        L1 = np.sum(np.abs(C - C_ana))/n
        L2 = np.sqrt(np.sum((C - C_ana)**2)/n)
        L_inf = np.max(np.abs(C - C_ana))

        L1_array.append(L1)
        L2_array.append(L2)
        L_inf_array.append(L_inf)
        dr_array.append(R/n)

    L1_array = np.array(L1_array)
    L2_array = np.array(L2_array)
    L_inf_array = np.array(L_inf_array)
    dr_array = np.array(dr_array)

    order_L1, intercept_L1 = np.polyfit(np.log(dr_array), np.log(L1_array), 1)
    order_L2, intercept_L2 = np.polyfit(np.log(dr_array), np.log(L2_array), 1)
    order_Linf, intercept_Linf = np.polyfit(np.log(dr_array), np.log(L_inf_array), 1)

    C_L1 = np.exp(intercept_L1)
    C_L2 = np.exp(intercept_L2)
    C_Linf = np.exp(intercept_Linf)

    # Plot computed reference lines
    ref_x = np.array([dr_array[0], dr_array[-1]])
    ref_y1 = C_L1 * ref_x**order_L1  # Computed reference line for L1
    ref_y2 = C_L2 * ref_x**order_L2  # Computed reference line for L2
    ref_yinf = C_Linf * ref_x**order_Linf  # Computed reference line for Linf

    plt.figure()
    plt.plot(dr_array, L1_array, 'o', label=f'$L_1$ numérique', color='tab:blue')
    plt.plot(dr_array, L2_array, 'o', label=f'$L_2$ numérique', color='tab:orange')
    plt.plot(dr_array, L_inf_array, 'o', label=f'$L_\infty$ numérique', color='tab:green')

    plt.plot(ref_x, ref_y1, "--", color='tab:blue')
    plt.plot(ref_x, ref_y2, "--", color='tab:orange')
    plt.plot(ref_x, ref_yinf, "--", color='tab:green')
    plt.plot([], [], "--", color='black', label='Régression en loi de puissance')

    equation_text = f'$L_1 = {C_L1:.4f} \\times Δr^{{{order_L1:.4f}}}$\n$L_2 = {C_L2:.4f} \\times Δr^{{{order_L2:.4f}}}$\n$L_\infty = {C_Linf:.4f} \\times Δr^{{{order_Linf:.4f}}}$'
    equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')
    equation_text_obj.set_position((0.5, 0.2))

    plt.legend()
    plt.xlabel(f'$\Delta r$ [m]')
    plt.ylabel('Erreur [mol/$m^3$]')
    plt.title("Convergence de l'erreur en fonction de $\Delta r$")
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()

    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)
    plt.tick_params(width=2, which='both', direction='in', top=True, right=True, length=6)

convergence_plot(C_order_1)
convergence_plot(C_order_2)