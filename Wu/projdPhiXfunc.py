

import numpy as np
from math import exp
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt

def fxn(phi_an, timeMax):
    # ----------------Constants----------------
    T = 298  # K
    n = 2
    R = 8.314  # kg m2/molK
    F = 96485  # C/mol
    V_an = 15e-3/1000  # m3
    M = 0.24541  # kg struvite/mol
    y = 1e-2  # m
    density = 1700  # kg/m3 struvite
    eps = 7.1e-10  # permittivity of water, F/m
    kb = 1.38065e-23  # Boltzmann Constant m2kg s-2 -1
    e = 1.602e-19  # C
    # ---------Variables-------------
    C_an = np.array([.12, .12, 0, 10 ** -7])  # NH4, PO4, Mg, H+. Molar
    z_k = np.array([1, -3, 2, 1])
    C_ca = C_an
    D_k = np.array([1.98e-9, 0.612e-9, 0.705e-9, 9.310e-9])  # m2/s from aquion.de
    phi_eq_an = -2.357  # V SHE
    phi_eq_ca = 0

    # C_NH4=C_an[0]
    # C_PO4=C_an[1]
    # C_H=C_an[3]
    # -----Water Chemistry-----------------
    K_eq = 10 ** -13  # Solubility product
    C = K_eq ** (1 / 3)  # concentration at saturation

    # ----Initialize some stuff-------
    # t = np.arange(1.0, 5000.0, 1)
    tMAX = timeMax
    t = np.linspace(1.0, tMAX, tMAX)
    t_xplot = [1]
    Q = np.zeros_like(t)
    N_Mg = np.zeros_like(t)
    C_Mg = np.zeros_like(t)
    C_NH4 = np.zeros_like(t)
    C_PO4 = np.zeros_like(t)
    N_struvite = np.zeros_like(t)
    L = np.zeros_like(t)
    iplot = np.zeros_like(t)
    Cbulkplot = np.zeros_like(t)
    A_zero = np.ones_like(t)
    A = 4e-3 * A_zero
    x = np.linspace(0, 5e-8, 10)
    phi_bulk = np.zeros((len(t_xplot), len(x)))
    phi_bulk2 = np.zeros((len(t_xplot), len(x)))
    phi_bulk1 = np.zeros((len(t_xplot), len(x)))
    C_Hbulk = np.zeros_like(t)
    C_H_surf = np.zeros_like(t)
    C_NH4plot=np.zeros(((len(t_xplot), len(x))))

    beta = 0.5
    A[0] = 4e-4  # m2
    # -----------Anode--------------------
    # Calculate i_an and Coulombs for each time step

    # def fxn(phi_an,t):
    C_Mg[0] = 1e-50  # Cannot start at 0
    C_Mgbefore = C_Mg[0]  # Variable for [Mg] from previous time step in loop

    C_NH4[0] = 0.12/1000 #m3
    C_PO4[0] = 0.12/1000
    N_struviteprev = 0
    eta_an = phi_an - phi_eq_an

    Cbulk = 0.12/1000 #m3
    C_Hbulk[0] = (10 ** -7)/1000 #pH=7
    count = 0
    t_xcount = 0
    Aold = 4e-3
    for time in t:
        xcount = 0
        A[count] = Aold
        k0_an = 1e-6  # rate constant of Mg2+ from A.Chadwick, et al. J. Electrochem Soc. 163 A1813 2016
        # i = 100*A[count] * n * F * (exp(-beta * n * F * eta_an / (R * T)))
        # i = A*i_o*(exp((1-beta)*F*eta_an/(R*T))-exp(-beta*F*eta_an/(R*T)))
        i_o = n * F * k0_an * C_Mgbefore

        i = Aold * 1e25 * (exp(-beta * n * F * eta_an / (R * T)))
        iplot[count] = i
        # print(iplot[count])
        Q[count] = i * time
        N_Mg[count] = Q[count] / (n * F)  # Faraday's Law to calculate moles of Mg2+
        C_Mg[count] = N_Mg[count] / V_an
        C_Mgbefore = C_Mg[count]  # make new [Mg] the old value for next iteration of loop

        Cbulkplot[count] = Cbulk
        kappa1 = (2 * F * e * Cbulk/ (kb * T * eps)) ** 0.5  # 1/Debye length
        C_NH4[count] = Cbulk * exp(1 * e * phi_an / (kb * T))  # Concentration at surface of anode, Boltzmann

        kappa2 = ((2 * (9) * F * e * Cbulk) / (kb * T * eps)) ** 0.5  # K, where 1/K is the Debye length
        C_PO4[count] = Cbulk * exp(-z_k[1] * e * phi_an / (kb * T))  # Concentration at surface of anode, Boltzmann

        Ksp = C_Mg[count] * C_NH4[count] * C_PO4[count]

        if C_Mg[count] >= C:  # Condition for precipitation
            N_struvite[count] = ((C_Mg[count] - C) * V_an) + N_struviteprev  # moles of struvite

            C_Mg[count] = C
            C_NH4[count] = C
            C_PO4[count] = C
            L[count] = M * Q[count] / (n * F * density * Aold)  # length of precipitate covering anode area

            Aold = Aold - (L[count] * y)  # New, uncovered area assuming uniform precipitation along y axis
            if Aold < 0:  # Cannot have negative area
                Aold = 0

            N_struviteprev = N_struvite[count]
            Cbulk = C_NH4[count]/exp( e * phi_an / (kb * T))  #
            if Cbulk < 0: # Cannot have negative concentration
                Cbulk = 0
        else:
            # C_NH4[count] = 0.12
            # C_PO4[count] = 0.12
            N_struvite[count] = N_struviteprev
            N_struviteprev = N_struvite[count]


        # Plot phi change with distance from electrode
        for element in t_xplot:
            if time == element:
                # Calculate phi change with distance
                for distance in x:
                    phi_bulk1[t_xcount][xcount] = phi_an / exp(
                        1 * kappa1 * distance)  # Poisson-Boltzmann equation to calculate Potential as a fxn of x
                    phi_bulk2[t_xcount][xcount] = phi_an / exp(
                        1 * kappa2 * distance)  # Poisson-Boltzmann equation to calculate Potential

                    phi_bulk[t_xcount][xcount] = phi_bulk2[t_xcount][xcount]
                    C_NH4plot[t_xcount][xcount] = Cbulk*exp(1 * e * phi_bulk[t_xcount][xcount] / (kb * T))
                    xcount += 1
                t_xcount += 1

        count += 1
    # -----------------Plot things vs X-------------------
    for things in range(len(t_xplot)):
# -------- Comment this out to plot phi- -----------
        plt.plot(x, C_NH4plot[things])
        plt.xlabel('Distance from anode, m')
        plt.ylabel('C_NH4, mol/m3')
# ------------- Comment this out to plot C_NH4----------
        # plt.plot(x, phi_bulk[things])
        # plt.xlabel('Distance from anode,m')
        # plt.ylabel('phi')

    return N_struvite, C_NH4plot, phi_bulk,A


# ---------------------------------------------------------------------------------------------------
phi_an2 = [-1,-0.9,-0.8,-0.7,-0.6]
N_struvite_plot = {}
C_NH4plot2 = {}
phi_bulkplot = {}
Aplot={}
x = np.linspace(0, 60e-9, 10)
j = 0   # Loop counter
tmax = 10000
t = np.linspace(1.0, tmax, tmax)

#Plot by calling function
for phi in phi_an2:
    N_struvite_plot[j], C_NH4plot2[j], phi_bulkplot[j] ,Aplot[j]= fxn(phi, tmax)
# ------------Comment this out to plot CH4 or phi (x)------------------
    plt.plot(t,N_struvite_plot[j])
    plt.plot(t, Aplot[j])
# ------------------------------------------------------
plt.legend(phi_an2)
# plt.xlabel('Time,s')
# plt.ylabel('Moles Struvite')
plt.show()









