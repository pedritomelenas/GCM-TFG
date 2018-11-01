import random
import arff
import numpy as np
import scipy.optimize as op
import data as dt
import commonFunctions as cf
from calLimits import calLimits
from intervalMinim import intervalMin
from varphiMinim import varphiMin
from WeighProd import WeighProd
from collections import deque
import matplotlib.pyplot as plt


#####################________ REVISAR PRECISIÓN __________######################
# print(1/((n-nu)*errs[0]**2))
# print(1/((n-nu)*errs[1]**2))
# print(1/((n-nu)*errs[2]**2))
# print(weights)
#####################_____________________________________#######################
galaxdata = {
    "radii": np.array([]),
    "vrot": np.array([]),
    "vbary": np.array([]),
    "weights": np.array([]),
    "CteDim": dt.CteDim,
    "totalnullvbary": False,
    "somenullvbary": False,
    "vones": np.array([]),
    "vv": np.array([]),
    "vvbary": np.array([]),
    "profile": ''
}
profiles = ['ISO', 'BUR', 'NFW']

for i in dt.galaxlist:
    print(" ****************** GALAXY: ", i, " **********************")
    radii = dt.galaxies[i]["R"]
    galaxdata["radii"] = radii
    #print("  R = ", radii)
    vrot = dt.galaxies[i]["vrot"]
    galaxdata["vrot"] = vrot
    #print("  Vrotacional = ", vrot)
    vbary = dt.galaxies[i]["vbary"]
    galaxdata["vbary"] = vbary
    #print("  Vbariónica = ", vbary)
    n = len(radii)
    vones = np.ones(n)
    galaxdata["vones"] = vones
    weights = 1 / ((n - dt.nu) * dt.galaxies[i]["errs"] ** 2)
    galaxdata["weights"] = weights
    totalnullvbary = np.sum(vbary) == 0
    galaxdata["totalnullvbary"] = totalnullvbary
    somenullvbary = round(np.prod(vbary)) == 0
    galaxdata["somenullvbary"] = somenullvbary
    vv = cf.vv(galaxdata)  # vrot * np.diag(weights) * vrot
    galaxdata["vv"] = vv
    vvbary = cf.vvbary(galaxdata)  # vbary * np.diag(weights) * vbary
    galaxdata["vvbary"] = vvbary

    for p in profiles:
        print(" ********** PROFILE: ", p, " ************")
        galaxdata["profile"] = p

        if n - dt.nu > 0:

            ### CALCULATION OF THE LIMITS ###
            limits = calLimits(galaxdata)
            varphiLim0 = limits[0]
            varphiLimInf = limits[1]
            #print("Xinf = ", Xinf)
            #print("X0 = ", X0)
            print("varphiLimInf", varphiLimInf)
            print("varphiLim0", varphiLim0)

            ### INTERVAL MINIMIZATION ###
            interval = intervalMin(varphiLim0, varphiLimInf, galaxdata)
            intervalinf = interval[0][0]
            intervalsup = interval[0][1]
            print("[", intervalinf, ", ", intervalsup, "]")
            Xi = interval[1]
            Yi = interval[2]
            X = np.logspace(np.log10(intervalinf), np.log10(intervalsup), 8)

            ### VARPHI MINIMIZATION ###
            varphimin = varphiMin(intervalinf, intervalsup, galaxdata)
            minvarphi = varphimin[0]
            minvarphiX = varphimin[1]
            Xj = varphimin[2]
            print("len(Xj) = ", len(Xj))
            Yj = varphimin[3]
            print("len(Yj) = ", len(Yj))
            #print("minphi = ", interval[3])    ## mínimo valor encontrado en la minimización del intervalo
            print("minvarphi = ", minvarphi)

            plt.semilogx()
            plt.scatter(X, np.zeros(len(X)), marker=3)
            plt.scatter(Xi, Yi, marker='.')
            plt.scatter(Xj, Yj, c='r', marker='.')
            plt.scatter(minvarphiX, minvarphi, c='y', marker='s')
            plt.hlines(varphiLimInf, 10 ** -2, intervalsup)
            plt.hlines(varphiLim0, intervalinf, 10)
            plt.show()
'''
# print(WeighProd(np.array([1,2,3]), np.array([1,2,3]), np.array([1,1,1])))
radii = dt.galaxies["DDO101"]["R"]
vrot = dt.galaxies["DDO101"]["vrot"]
vbary = dt.galaxies["DDO101"]["vbary"]
errs = dt.galaxies["DDO101"]["errs"]
n = len(radii)
weights = 1 / ((n - dt.nu) * dt.galaxies["DDO101"]["errs"] ** 2)
totalnullvbary = np.sum(vbary) == 0
somenullvbary = round(np.prod(vbary)) == 0

vv = WeighProd(vrot, vrot, weights)  # vrot * np.diag(weights) * vrot
vvbary = WeighProd(vbary, vbary, weights)  # vbary * np.diag(weights) * vbary
vones = np.ones(n)

if n - dt.nu > 0:

    #########################################################
    ##### Computation of the limit values for varphi(s) #####
    #########################################################

    def ginf(x):
        return np.ones(len(x))


    # Function defining equation #34 --> #35
    def eqVLimInf(t):
        return WeighProd(ginf(radii), vones, weights) - \
               WeighProd(vrot, (ginf(radii) / np.sqrt(t * ginf(radii) + vbary ** 2)), weights)
        # (ginf(radii) * np.diag(weights) * vones) / \
        # (vrot * np.diag(weights) * (ginf(radii) / np.sqrt(t * ginf(radii) + vbary ** 2)))
        # WeighProd(ginf(radii),vones,weights) - WeighProd(vrot,ginf(radii)./sqrt(t*ginf(radii)+(vbary.^2)),weights)


    def g0(x):
        return x ** 2


    # Function defining equation #38 --> #40
    def eqVLim0(t):
        # WeighProd(g0(radii),vones,weights) - WeighProd(vrot,g0(radii)./sqrt(t.*g0(radii)+(vbary.^2)),weights)
        return WeighProd(g0(radii), vones, weights) - WeighProd(vrot, g0(radii) / np.sqrt(t * g0(radii) + vbary ** 2),
                                                                weights)

    # print(WeighProd(vrot, np.sqrt(ginf(radii)), weights))
    # print(WeighProd(ginf(radii), vones, weights))
    # print(WeighProd(vrot, np.sqrt(ginf(radii)), weights)/WeighProd(ginf(radii), vones, weights))
    # print(ginf(radii) * np.diag(weights) * vones)
    # print(np.shape(ginf(radii) * np.diag(weights) * vones))

    # The solution to #34 --> #35 with vbary=0, which is an upper bound for the solutions to #34 --> #35
    XVbaryNullinf = (WeighProd(vrot, np.sqrt(ginf(radii)), weights) / WeighProd(ginf(radii), vones, weights)) ** 2
    print("XVbaryNullinf = ", XVbaryNullinf)
    # The solution to #40 is an upper bound for the solutions to #40
    XVbaryNull0 = (WeighProd(vrot, np.sqrt(g0(radii)), weights) / WeighProd(g0(radii), vones, weights)) ** 2
    print("XVbaryNull0 = ", XVbaryNull0)
    # print(XVbaryNullinf)
    # print(XVbaryNull0)

    if totalnullvbary:
        # print("if")
        Xinf = XVbaryNullinf
        X0 = XVbaryNull0
    elif somenullvbary: ####################################################################### NO ENTIENDO ESTE CASO ##########
        # print(round(np.prod(vbary)) == 0)
        jinf = -3
        # print(eqVLimInf(10 ** (j) * XVbaryNullinf))
        # print(eqVLim0(10 ** jinf * XVbaryNullinf))
        while eqVLimInf(10 ** jinf * XVbaryNullinf) > 0:
            jinf -= 1
        j0 = -3
        while eqVLim0(10 ** j0 * XVbaryNull0) > 0:
            j0 -= 1
        Xinf = op.brentq(eqVLimInf, 10 ** jinf * XVbaryNullinf, XVbaryNullinf)  # Brent's Method (escalar)
        # Xinf = op.fsolve(eqVLimInf, (10 ** jinf * XVbaryNullinf + XVbaryNullinf) / 2)
        X0 = op.brentq(eqVLim0, 10 ** j0 * XVbaryNull0, XVbaryNull0)  # Brent's Method (escalar)
        # X0 = op.fsolve(eqVLim0, (10 ** j0 * XVbaryNull0 + XVbaryNull0) / 2)  ## USAR ESCALAR
        # X = fzero(equationVLimInf, [10 ^ (j) * XVbaryNull, XVbaryNull]); # Solves equation  # 34 --> #35

    else:
        if eqVLimInf(0) >= 0:   # >= ? #########################################################################################
            Xinf = 0
        else:
            Xinf = op.fsolve(eqVLimInf, XVbaryNullinf / 2)  # Solves equation --> #35
        if eqVLim0(0) >= 0:
            X0 = 0
        else:
            X0 = op.fsolve(eqVLim0, XVbaryNull0 / 2)  # Solves equation #38 --> #40

    ##### Xinf y X0 calculado con brentq y con f solve. ################################################################
    ##### En el caso de somenullvbary, ¿por qué usar brentq y no fsolve?
    ##### (brentq es ESCALAR, fsolve es MULTIDIMENSIONAL)
    ##### brentq devuelve UNA RAÍZ de f, fsolve devuelve LAS RAÍCES de f

    ######### PRUEBA PARA LA DUDA ENTRE fsolve Y brentq ###########
    jinf = -3
    while eqVLimInf(10 ** jinf * XVbaryNullinf) > 0:
        jinf -= 1
    #print("Xinf with brentq = ", op.brentq(eqVLimInf, 10 ** jinf * XVbaryNullinf, XVbaryNullinf))
    #print("Xinf with fsolve = ", op.fsolve(eqVLimInf, XVbaryNullinf / 2))
    j0 = -3
    while eqVLim0(10 ** j0 * XVbaryNull0) > 0:
        j0 -= 1
    #print("X0 with brentq = ", op.brentq(eqVLim0, 10 ** j0 * XVbaryNull0, XVbaryNull0))
    #print("X0 with fsolve = ", op.fsolve(eqVLim0, XVbaryNull0 / 2))

    ####################################################################################################################

    ## Calculation of the limit value by using Lemma 2.1 and development #16 ##
    varphiLimInf = vv + vvbary - 2 * WeighProd(vrot, np.sqrt(Xinf * ginf(radii) + (vbary ** 2)), weights) + \
                   Xinf * WeighProd(ginf(radii), vones, weights)
    # varphiLimInf = vv + vvbary - 2. * WeighProd(vrot, sqrt(X * ginf(radii) + (vbary. ^ 2)), weights) +
    # X. * WeighProd(ginf(radii), vones, weights)

    ## Calculation of the limit value by using Lemma 2.2 (2.1?) and development #16 ##
    varphiLim0 = vv + vvbary + X0 * WeighProd(g0(radii), vones, weights) \
                 - 2 * WeighProd(vrot, np.sqrt(X0 * g0(radii) + (vbary ** 2)), weights)

    #print("varphiLimInf = ", varphiLimInf)
    #print("varphiLim0 = ", varphiLim0)



    print("phi(", (Xinf + X0) / 2, ")=", phi(np.array([(Xinf + X0) / 2])))
    print("Xinf = ", Xinf)
    print("phi(Xinf)=", phi(Xinf))
    print("phi(1000)=", phi(np.array([1000.0])))
    print("phi(X0)=", phi(X0))
    print("phi(0.00000000001)=", phi(np.array([0.00000000001])))  # 10^-11
    print("phi(0.00331078)=", phi(np.array([0.00331078])))  #
    # print("phi(0.000000000001)=", phi(np.array([0.000000000001]))) #10^-12


    # print("****", alphaMVISO(10**(-3+np.array([-0.2000,-0.1000,-0.0000,0.1000,0.2000]))))
    # print("----", phi(10**(-3+np.array([-0.2000,-0.1000,-0.0000,0.1000,0.2000]))))


    ###############################################
    ##### Minimization interval determination #####
    ###############################################

    ##### Extremo inferior #####
    #einf = -3
    #dir = -1
    #iter = 0
    #print("s antes de alpha = ", (10**(-3 + np.array([-0.2,-0.1,0.0,0.1,0.2]))))
    #print("alphaMVISO", alphaMVISO(10 ** (-3 + np.array([-0.2, -0.1, 0.0, 0.1, 0.2]))))
    #print("vbary", vbary)
    #print("eval = ", eval)

    #while iter < 100 and dir != 0:
     #   iter += 1
      #  eval = abs(
       #             vv + vvbary + alphaMVISO(10 ** (-3 + np.array([-0.2, -0.1, 0.0, 0.1, 0.2]))) - varphiLim0) / varphiLim0
        #izq = sum(eval[0:3])  # mejor hasta 2 ([0:2]) y luego eval[3]+eval[4] ??


    ###############################################
    ############### Tabu Search ###################
    ###############################################

    s0 = random.uniform(0, 10**3)   # solución inicial aleatoria
    Nv = 5  # mínimo número de soluciones en la vecindad actual
    Nmin = 10   # mínimo número de soluciones en la lista tabú
    Nmax = 100  # número máximo de soluciones en la lista tabú
    K = 15000   # máximo número de iteraciones
    e = 10**(-3)    # umbral de precisión
    M = 10  # máximo número de soluciones encontradas
    N = 30  # máximo número de iteraciones sin que la solución óptima cambie

    tabulist = deque()
    #print(tabulist)
'''