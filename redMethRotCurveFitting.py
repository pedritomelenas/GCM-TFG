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
import time


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
    "profile": '',
    "graphic": False
}
profiles = ['ISO', 'BUR', 'NFW']
start_time = time.time()
for i in dt.galaxlist:
    istart = time.time()
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
    galaxdata["graphic"] = True

    for p in profiles:
        print(" ********** PROFILE: ", p, " ************")
        galaxdata["profile"] = p
        pstart = time.time()
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
            if galaxdata["graphic"]:
                Xi = interval[1]
                Yi = interval[2]
                intminvarphi = interval[3]
                intminvarphiX = interval[4]
            else:
                intminvarphi = interval[1]
                intminvarphiX = interval[2]
            X = np.logspace(np.log10(intervalinf), np.log10(intervalsup), 8)

            ### VARPHI MINIMIZATION ###
            pmin = varphiMin(intervalinf, intervalsup, galaxdata)
            minvarphi = pmin[0]
            minvarphiX = pmin[1]
            if galaxdata["graphic"]:
                Xj = pmin[2]
                Yj = pmin[3]
                forkpoints = pmin[4]
            else:
                forkpoints = pmin[2]
            #print("minphi = ", interval[3])    ## mínimo valor encontrado en la minimización del intervalo
            print("minvarphi = ", minvarphi)
            if galaxdata["graphic"]:
                plt.semilogx()
                plt.title("Galaxia "+i+" con perfil "+p)
                plt.scatter(X, np.zeros(len(X)), marker=3)
                plt.scatter(forkpoints, np.zeros(len(forkpoints)), c='r', marker=3)
                plt.scatter(Xi, Yi, marker='.')
                plt.scatter(Xj, Yj, c='r', marker='.')
                plt.scatter(minvarphiX, minvarphi, c='y', marker='s')
                plt.hlines(varphiLimInf, 10 ** -2, intervalsup)
                plt.hlines(varphiLim0, intervalinf, 10)
                #plt.show()
                plt.savefig("C:/Users/marin/PycharmProjects/TFG/galaxies/graphics/" + p + '-' + i + '.png')
                plt.gcf().clear()
        pend = time.time()
        print("Tiempo para el perfil ", p, " para la galaxia ", i, " = ", pend - pstart, " segundos")
    iend = time.time()
    print("Tiempo para la galaxia ", i, " = ", iend - istart, " segundos")
endtime = time.time()
print("Tiempo total = ", endtime - start_time, " segundos")
'''

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