import sympy as sp
import matplotlib.pyplot as plt
import numpy as np

# Indepentant variables
r, t = sp.symbols('r t')

# Parameters
R, C_e, k, Deff = sp.symbols('R C_e k Deff')

# Analytical solution
C_ana = C_e*(sp.exp(1e-10*t) - 1) * sp.cos((sp.pi / 2) * (r / R)) + C_e

dCdt = sp.diff(C_ana, t)
dCdr = sp.diff(C_ana, r)
d2Cdr2 = sp.diff(dCdr, r)

# Boundary and initial conditions
dCdr_0 = sp.diff(C_ana, r).subs(r, 0)
C_R = C_ana.subs(r, R)
C_t0 = C_ana.subs(t, 0)

# Source term
f = dCdt - Deff*(1/r*sp.diff(r*dCdr, r)) + k*C_ana

# print
sp.pprint(f"Analytical solution: {C_ana}")
sp.pprint(f"Boundary condition left : {dCdr_0}")
sp.pprint(f"Boundary condition right : {C_R}")
sp.pprint(f"Initial condition : {C_t0}")
sp.pprint(f"Source term : {f}")

# sp.pprint(Deff*(1/r*sp.diff(r*dCdr, r)) - k*C_ana + f)

# Plot
r_range = np.linspace(0, 0.5, 5)
t_range = np.linspace(0, 4e9, 10)

C = sp.lambdify((r, t, R, C_e), C_ana, 'numpy')

plt.figure()
for t_val in t_range:
    plt.plot(r_range, C(r_range, t_val, 0.5, 20), label=f"t = {t_val}")
plt.legend()
plt.xlabel('r [m]')
plt.ylabel('C [mol/$m^3$]')

# plot source
f_source = sp.lambdify((r, t, R, C_e, Deff, k), f, 'numpy')
plt.figure()
for t_val in t_range:
    plt.plot(r_range, f_source(r_range, t_val, 0.5, 20, 1e-10, 4e-9), label=f"t = {t_val}")
plt.legend()
plt.xlabel('r [m]')
plt.ylabel('f [mol/$m^3$/s]')