import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


D_eff = 10**(-10)
k = 4*10**(-9)
Ce = 20

N_tot = 5
R = 0.5

t_fin = 4*10**9


# Définition du terme source
def source(r,t):
    c1 = np.pi*r**2/2/R**2
    c2 = 10**(-10)
    f = 1/R**4*(Ce*(-np.pi*D_eff*(1 - np.exp(c2*t))*(2*R**2*np.sin(c1) + np.pi*r**2*np.cos(c1))
                    + R**4*(k*(-1*(1 - np.exp(c2*t))*np.cos(c1) + 1) + c2*np.exp(c2*t)*np.cos(c1))))
    return f



def mat(n,dt,C_ini):
    ri = np.zeros(n)
    dr = R/(n - 1)
    nt = t_fin//dt
    if t_fin%dt != 0:
        nt += 1
    for i in range(n):
        ri[i] = i*dr
    
    k = 4*10**(-9)
    alpha = D_eff*dt/dr
    
    coef = np.zeros([n,n])
    coef[0,0] = -3
    coef[0,1] = 4
    coef[0,2] = -1
    coef[-1,-1] = 1
    RHS = np.zeros(n)
    RHS[-1] = Ce
    for i in range(1,n-1):
        coef[i,i-1] = alpha*(1/(2*ri[i]) - 1/dr)
        coef[i,i] = 1 + 2*alpha/dr + k*dt
        coef[i,i+1] = -alpha*(1/(2*ri[i]) + 1/dr)
    
    Ct0 = RHS
    Ct1 = np.zeros(n)
    t_actu = 0
    for k in range(int(nt)):
        t_actu += dt
        RHS[0] = 0
        for i in range(1,n-1):
            RHS[i] = Ct0[i]
        RHS[-1] = Ce
        Ct1 = np.linalg.solve(coef,RHS)
        Ct0 = Ct1
    
    for i in range(len(ri)):
        C_ana[i] = CMMS(ri[i],t_fin)
    
    return ri,dr,Ct0,C_ana



def matsource(n,dt,C_ini):
    ri = np.zeros(n)
    dr = R/(n - 1)
    nt = t_fin//dt
    if t_fin%dt != 0:
        nt += 1
    for i in range(n):
        ri[i] = i*dr
    
    k = 4*10**(-9)
    alpha = D_eff*dt/dr
    
    coef = np.zeros([n,n])
    coef[0,0] = -3
    coef[0,1] = 4
    coef[0,2] = -1
    coef[-1,-1] = 1
    RHS = np.zeros(n)
    RHS[-1] = Ce
    for i in range(1,n-1):
        coef[i,i-1] = alpha*(1/(2*ri[i]) - 1/dr)
        coef[i,i] = 1 + 2*alpha/dr + k*dt
        coef[i,i+1] = -alpha*(1/(2*ri[i]) + 1/dr)
    
    Ct0 = RHS
    Ct1 = np.zeros(n)
    t_actu = 0
    for k in range(int(nt)):
        t_actu += dt
        RHS[0] = 0
        for i in range(1,n-1):
            RHS[i] = Ct0[i] + dt*source(ri[i],t_actu)
        RHS[-1] = Ce
        Ct1 = np.linalg.solve(coef,RHS)
        Ct0 = Ct1
    
    C_ana = np.zeros(len(ri))
    for i in range(len(ri)):
        C_ana[i] = CMMS(ri[i],t_fin)
    
    return ri,dr,Ct0,C_ana



def CMMS(r,t):
    sol = Ce*((np.exp(10**(-10)*t) - 1)*np.cos(np.pi*r**2/2/R**2) + 1)
    return sol



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




r1,dr1,C1,C_ana = matsource(N_tot,60*60*24*365,np.zeros(N_tot))
plt.plot(r1,C_ana,"k--",label='Solution analytique')
plt.scatter(r1,C1,label='Solution numérique')
plt.xlabel("Position [m]")
plt.xlim([0,R])
plt.ylabel(r"Concentration [$mol/m^3$]")
plt.title("Concentration en sel en fonction de la position radiale dans le pilier \n Schéma avant de la dérivée première")
plt.legend()
plt.grid()
plt.show()



# Analyse de convergence en espace
it = 6
dr = np.zeros(it)
L1 = np.zeros(it)
L2 = np.zeros(it)
L_inf = np.zeros(it)
N = 4
for i in range(it):
    _,dr[i],C,C_ana = matsource(N,400000,np.zeros(N_tot))
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
plt.title("Analyse de convergence, espace")
plt.legend()
plt.grid()
plt.show()



# Analyse de convergence en temps
it = 4
dt = np.zeros(it)
L1 = np.zeros(it)
L2 = np.zeros(it)
L_inf = np.zeros(it)
N = 1
for i in range(it):
    _,_,C,C_ana = matsource(61,60*60*365/N,np.zeros(11))
    dt[i] = 60*60*24*365/N
    L1[i],L2[i],L_inf[i] = erreurs(C,C_ana)
    N *= 2
p1 = np.zeros(it-1)
p2 = np.zeros(it-1)
p_inf = np.zeros(it-1)
for i in range(it-1):
    p1[i] = np.log(L1[i]/L1[i+1])/np.log(dt[i]/dt[i+1])
    p2[i] = np.log(L2[i]/L2[i+1])/np.log(dt[i]/dt[i+1])
    p_inf[i] = np.log(L_inf[i]/L_inf[i+1])/np.log(dt[i]/dt[i+1])
print("Ordre de l'erreur L1 : ", np.mean(p1))
print("Ordre de l'erreur L2 : ", np.mean(p2))
print("Ordre de l'erreur L_inf : ", np.mean(p_inf))
plt.loglog(dt,L1,label=r'Erreur $L_1$')
plt.loglog(dt,L2,label=r'Erreur $L_2$')
plt.loglog(dt,L_inf,label=r'Erreur $L_{\infty}$')
plt.xlabel("Pas de temps [s]")
plt.ylabel("Erreur")
plt.title("Analyse de convergence, temps")
plt.legend()
plt.grid()
plt.show()








