import numpy as np
import matplotlib.pyplot as plt

# Définition des paramètres de l'équation à résoudre
D_eff = 10**(-10)
S = 2*10**(-8)
Ce = 20

# Définition du domaine
N_tot = 5
R = 0.5


# Système matriciel d'ordre 1
def mat1(n):
    r = np.zeros(n)
    dr = R/(n - 1)
    C_ana = np.zeros(n)
    for i in range(n):
        r[i] = i*dr
        C_ana[i] = S*R**2/(4*D_eff)*(r[i]**2/R**2 - 1) + Ce
        
    coef = np.zeros([n,n])
    coef[0,0] = -1
    coef[0,1] = 1
    coef[-1,-1] = 1
    RHS = np.zeros(n)
    RHS[-1] = Ce
    for i in range(1,n-1):
        coef[i,i-1] = 1/dr**2
        coef[i,i] = -1/(r[i]*dr) - 2/dr**2
        coef[i,i+1] = 1/(r[i]*dr) + 1/dr**2
        RHS[i] = S/D_eff
    C = np.linalg.solve(coef,RHS)
    
    return r, dr, C, C_ana


# Système matriciel d'ordre 2
def mat2(n):
    r = np.zeros(n)
    dr = R/(n - 1)
    C_ana = np.zeros(n)
    for i in range(n):
        r[i] = i*dr
        C_ana[i] = S*R**2/(4*D_eff)*(r[i]**2/R**2 - 1) + Ce
        
    coef = np.zeros([n,n])
    coef[0,0] = -3
    coef[0,1] = 4
    coef[0,2] = -1
    coef[-1,-1] = 1
    RHS = np.zeros(n)
    RHS[-1] = Ce
    for i in range(1,n-1):
        coef[i,i-1] = 1/dr**2 - 1/(2*r[i]*dr)
        coef[i,i] = -2/dr**2
        coef[i,i+1] = 1/dr**2 + 1/(2*r[i]*dr)
        RHS[i] = S/D_eff
    C = np.linalg.solve(coef,RHS)
    
    return r,dr,C,C_ana


# Définition de l'erreur
def erreurs(C,C_ana):
    N = len(C)
    L1 = 0
    L2 = 0
    for i in range(N):
        L1 += abs(C[i]-C_ana[i])/N
        L2 += abs(C[i]-C_ana[i])**2/N
    L2 = L2**0.5
    L_inf = max(abs(C-C_ana))
    
    return L1,L2,L_inf



### QUESTION D)

# Affichage des solutions analytique et numériques
r1,dr1,C1,C_ana1 = mat1(N_tot)
plt.plot(r1,C_ana1,"k--",label='Solution analytique')
plt.plot(r1,C1,label='Solution numérique')
plt.xlabel("Position [m]")
plt.xlim([0,R])
plt.ylabel(r"Concentration [$mol/m^3$]")
plt.title("Concentration en sel en fonction de la position radiale dans le pilier \n Schéma avant de la dérivée première")
plt.legend()
plt.grid()
plt.show()

# Calcul de l'erreur pour différents maillages, ordre 1
it = 6
dr = np.zeros(it)
L1 = np.zeros(it)
L2 = np.zeros(it)
L_inf = np.zeros(it)
N = 4
for i in range(it):
    _,dr[i],C,C_ana = mat1(N+1)
    L1[i],L2[i],L_inf[i] = erreurs(C,C_ana)
    N *= 2
p1 = np.zeros(it-1)
p2 = np.zeros(it-1)
p_inf = np.zeros(it-1)
for i in range(it-1):
    p1[i] = np.log(L1[i]/L1[i+1])/np.log(dr[i]/dr[i+1])
    p2[i] = np.log(L2[i]/L2[i+1])/np.log(dr[i]/dr[i+1])
    p_inf[i] = np.log(L_inf[i]/L_inf[i+1])/np.log(dr[i]/dr[i+1])
print("Ordre de l'erreur L1 : ", np.mean(p1))
print("Ordre de l'erreur L2 : ", np.mean(p2))
print("Ordre de l'erreur L_inf : ", np.mean(p_inf))
plt.loglog(dr,L1,label=r'Erreur $L_1$')
plt.loglog(dr,L2,label=r'Erreur $L_2$')
plt.loglog(dr,L_inf,label=r'Erreur $L_{\infty}$')
plt.xlabel("Pas de distance [m]")
plt.ylabel("Erreur")
plt.title("Erreurs \n Schéma avant de la dérivée première")
plt.legend()
plt.grid()
plt.show()



### QUESTION E)

# Calcul de l'erreur pour différents maillages, ordre 2
r2,dr2,C2,C_ana2 = mat2(N_tot)
plt.plot(r2,C_ana2,"k--",label='Solution analytique')
plt.plot(r1,C1,label='Solution numérique, schéma avant')
plt.plot(r2,C2,label='Solution numérique, schéma centré')
plt.xlabel("Position [m]")
plt.xlim([0,R])
plt.ylabel(r"Concentration [$mol/m^3$]")
plt.title("Concentration en sel en fonction de la position radiale dans le pilier \n Différents schémas de la dérivée première")
plt.legend()
plt.grid()
plt.show()

it = 6
dr = np.zeros(it)
L1 = np.zeros(it)
L2 = np.zeros(it)
L_inf = np.zeros(it)
N = 4
for i in range(it):
    _,dr[i],C,C_ana = mat2(N+1)
    L1[i],L2[i],L_inf[i] = erreurs(C,C_ana)
    N *= 2
p1 = np.zeros(it-1)
p2 = np.zeros(it-1)
p_inf = np.zeros(it-1)
for i in range(it-1):
    p1[i] = np.log(L1[i]/L1[i+1])/np.log(dr[i]/dr[i+1])
    p2[i] = np.log(L2[i]/L2[i+1])/np.log(dr[i]/dr[i+1])
    p_inf[i] = np.log(L_inf[i]/L_inf[i+1])/np.log(dr[i]/dr[i+1])
print("Ordre de l'erreur L1 : ", np.mean(p1))
print("Ordre de l'erreur L2 : ", np.mean(p2))
print("Ordre de l'erreur L_inf : ", np.mean(p_inf))
plt.loglog(dr,L1,label=r'Erreur $L_1$')
plt.loglog(dr,L2,label=r'Erreur $L_2$')
plt.loglog(dr,L_inf,label=r'Erreur $L_{\infty}$')
plt.xlabel("Pas de distance [m]")
plt.ylabel("Erreur")
plt.title("Erreurs \n Schéma centré de la dérivée première")
plt.legend()
plt.grid()
plt.show()






