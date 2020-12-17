import numpy as np
from math import exp
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
#----------------Constants----------------
T=298 #K
n=2
R=8.314 #kg m2/molK
F=96485 #C/mol
V_an=15e-3 #L
M=0.24541 #kg struvite/mol
y=1e-2 #m
density=1700 #kg/m3 struvite
eps= 7.1e-10 # permittivity of water, F/m
kb= 1.38065e-23 # Boltzmann Constant m2kg s-2 -1
e=1.602e-19 # C
#---------Variables-------------
C_an=np.array([.12, .12,0, 10**-7]) #NH4, PO4, Mg, H+. Molar
z_k=np.array([1,-3,2,1])
C_ca=C_an
D_k=np.array([1.98e-9, 0.612e-9,0.705e-9, 9.310e-9])#m2/s from aquion.de
phi_eq_an=-2.357 #V SHE
phi_eq_ca=0

# C_NH4=C_an[0]
# C_PO4=C_an[1]
# C_H=C_an[3]
# -----Water Chemistry-----------------
K_eq=10**-13 # Solubility product
C=K_eq**(1/3) # concentration at saturation


# ----Initialize some stuff-------
t=np.arange(1.0,5000.0,1)
Q=np.zeros_like(t)
N_Mg=np.zeros_like(t)
C_Mg=np.zeros_like(t)
C_NH4=np.zeros_like(t)
C_PO4=np.zeros_like(t)
N_struvite=np.zeros_like(t)
L=np.zeros_like(t)
iplot=np.zeros_like(t)
Cbulkplot=np.zeros_like(t)
A_zero=np.ones_like(t)
A=4e-3*A_zero

phi_an = -1.0  # V
beta = 0.5
A[0] = 4e-4  # m2
# -----------Anode--------------------
# Calculate i_an and Coulombs for each time step


C_Mg[0] = 1e-50  # Cannot start at 0
C_Mgbefore = C_Mg[0]  # Variable for [Mg] from previous time step in loop

C_NH4[0] = 0.12
C_PO4[0] = 0.12
N_struviteprev = 0
eta_an = phi_an - phi_eq_an

Cbulk=0.12
count = 0
Aold=4e-3
for time in t:
    A[count]=Aold
    k0_an = 1e-6  # rate constant of Mg2+ from A.Chadwick, et al. J. Electrochem Soc. 163 A1813 2016
    # i = 100*A[count] * n * F * (exp(-beta * n * F * eta_an / (R * T)))
    # i = A*i_o*(exp((1-beta)*F*eta_an/(R*T))-exp(-beta*F*eta_an/(R*T)))
    i_o = n * F * k0_an * C_Mgbefore

    i = Aold*i_o*(exp(-beta * n * F * eta_an / (R * T)))
    iplot[count]=i
    #print(iplot[count])
    Q[count] = i * time
    N_Mg[count] = Q[count] / (n * F)  # Faraday's Law to calculate moles of Mg2+
    C_Mg[count] = N_Mg[count] / V_an
    C_Mgbefore = C_Mg[count]  # make new [Mg] the old value for next iteration of loop

    Cbulkplot[count]=Cbulk
    kappa1=((2*F*e*Cbulk)/(kb*T*eps))**0.5 #Debye length
    #phi_bulk=phi_an/exp(-1*kappa1*2) #Poisson-Boltzmann equation to calculate Potential at anode surface, r=kappa
    C_NH4[count]=Cbulk*exp(-z_k[0]*e*phi_an/(kb*T)) # Concentration at surface of anode, Boltzmann

    kappa2 = ((2 * (9) * F* e * Cbulk) / (kb * T * eps)) ** 0.5  # K, where 1/K is the Debye length
    #phi_bulk = phi_an / exp(-1*kappa2 * 2)  # Poisson-Boltzmann equation to calculate Potential after dbl, r=kappa
    C_PO4[count] = Cbulk * exp(-z_k[1] *e * phi_an / (kb * T))  # Concentration at surface of anode, Boltzmann


    Ksp = C_Mg[count] * C_NH4[count] * C_PO4[count]

    if C_Mg[count] >= C:  # Condition for precipitation since [Mg] is limiting
        N_struvite[count] = ((C_Mg[count] - C) * V_an) + N_struviteprev #moles of struvite

        C_Mg[count] = C
        C_NH4[count] = C
        C_PO4[count] = C
        L[count] = M * Q[count] / (n * F * density * Aold)  # length of precipitate covering anode area


        Aold = Aold - (L[count] * y)  # New, uncovered area assuming uniform precipitation along y axis
        if Aold<0: #Cannot have negative area
           Aold=0

        N_struviteprev = N_struvite[count]
        Cbulk=Cbulk-C #Cannot have negative concentration
        if Cbulk<0:
            Cbulk=0
    else:
        # C_NH4[count] = 0.12
        # C_PO4[count] = 0.12
        N_struvite[count] =  N_struviteprev
        N_struviteprev = N_struvite[count]
        Cbulk=Cbulk
    # print(N_struviteprev)
    count += 1



plt.plot(t, N_struvite)
plt.xlabel('Time (s)')
plt.ylabel('Struvite moles')
plt.show()
#
# plt.plot(t, A)
# plt.xlabel('Time (s)')
# plt.ylabel('A')
# plt.show()
