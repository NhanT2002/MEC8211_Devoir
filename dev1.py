import numpy as np
import matplotlib.pyplot as plt

Ce = 20 # [mol/m3]
D = 1 # [m]
R = D/2 # [m]
D_eff = 1E-10 # [m2/s]
k = 4E-9 # [s^-1]
S = 2E-8 # [mol/m3/s]

N = 100 # number of grid points

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

plt.figure()
plt.plot(r, C, label='Solution numérique')
plt.plot(r, C_ana, label='Solution analytique')
plt.xlabel('r [m]')
plt.ylabel(f'C [mol/$m^3$]')
plt.legend()
plt.grid()
plt.show()