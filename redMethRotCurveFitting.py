import numpy as np
import data as dt
import commonFunctions as cf
from calLimits import calLimits
from intervalMinim import intervalMin
from varphiMinim import varphiMin
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from galaxygraphic3D import generate3Dgraphic


galaxdata = {                          # Colocar en data.py?
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

start_time = time.time()
for i in dt.galaxlist:
    istart = time.time()
    print(" ****************** GALAXY: ", i, " **********************")
    radii = dt.galaxies[i]["R"]
    galaxdata["radii"] = radii
    vrot = dt.galaxies[i]["vrot"]
    galaxdata["vrot"] = vrot
    vbary = dt.galaxies[i]["vbary"]
    galaxdata["vbary"] = vbary
    n = len(radii)
    vones = np.ones(n)
    galaxdata["vones"] = vones
    weights = 1 / ((n - dt.nu) * dt.galaxies[i]["errs"] ** 2)
    galaxdata["weights"] = weights
    totalnullvbary = np.sum(vbary) == 0
    galaxdata["totalnullvbary"] = totalnullvbary
    somenullvbary = round(np.prod(vbary)) == 0
    galaxdata["somenullvbary"] = somenullvbary
    vv = cf.vv(galaxdata)
    galaxdata["vv"] = vv
    vvbary = cf.vvbary(galaxdata)
    galaxdata["vvbary"] = vvbary
    galaxdata["graphic"] = True

    for p in dt.profiles:
        print(" ********** PROFILE: ", p, " ************")
        galaxdata["profile"] = p
        pstart = time.time()
        if n - dt.nu > 0:
            ### CALCULATION OF THE LIMITS ###
            limits = calLimits(galaxdata)
            varphiLim0 = limits[0]
            varphiLimInf = limits[1]
            #print("varphiLimInf", varphiLimInf)
            #print("varphiLim0", varphiLim0)

            ### INTERVAL MINIMIZATION ###
            intervalmin_start = time.time()
            interval = intervalMin(varphiLim0, varphiLimInf, galaxdata)
            intervalinf = interval[0][0]
            intervalsup = interval[0][1]
            print("[", intervalinf, ", ", intervalsup, "]")
            if galaxdata["graphic"]:
                Xi = interval[1]
                Yi = interval[2]
                intinfmin = interval[3]
                intsupmin = interval[4]
            else:
                intinfmin = interval[1]
                intsupmin = interval[2]
            intervalmin_end = time.time()
            print("Tiempo de minimización del intervalo = ", intervalmin_end - intervalmin_start, " segundos")
            ### VARPHI MINIMIZATION ###
            pmin = varphiMin(varphiLim0, varphiLimInf, intinfmin, intsupmin, intervalinf, intervalsup, galaxdata)
            minvarphi = pmin[0]
            minrho = pmin[1]
            minvarphiX = pmin[2]
            if galaxdata["graphic"]:
                Xj = pmin[3]
                Yj = pmin[4]
                forkpoints = pmin[5]
                X = pmin[6]
            else:
                forkpoints = pmin[3]
                X = pmin[4]
            print("minvarphi = ", minvarphi)
            print("para s = ", minvarphiX)
            print("con rho(s) = ", minrho)

            if galaxdata["graphic"]:
                plt.semilogx()
                #X = np.logspace(-6, 20, 5000)               # Gráfica de ejemplo
                plt.title("Galaxia "+i+" con perfil "+p)
                #phiX, rho = cf.phi(X, galaxdata)           # Gráfica de ejemplo
                #plt.plot(X, phiX, 'k')                     # Gráfica de ejemplo
                plt.xlabel("s (parámetro de escala)")      # Gráfica de ejemplo
                plt.ylabel(r"$\varphi(s)$")                # Gráfica de ejemplo
                #plt.vlines(intervalsup, minvarphi, varphiLimInf,
                #           colors='b', linestyles='dashed')            # Ejemplo: diferencia entre condiciones de convergencia
                #plt.gca().set_xticks([intervalsup])                    # Ejemplo: diferencia entre condiciones de convergencia
                #plt.gca().set_xticklabels([""+str(intervalsup)+""])    # Ejemplo: diferencia entre condiciones de convergencia
                plt.scatter(intervalinf, 0, c='black', marker=3)
                plt.scatter(intervalsup, 0, c='black', marker=3)
                #plt.vlines(intervalinf, -0.1, 0.1)
                #plt.vlines(intervalsup, -0.1, 0.1)
                plt.hlines(0, intervalinf, intervalsup)
                plt.scatter(X, np.zeros(len(X)), color='black', marker=3)
                #plt.scatter(forkpoints, np.zeros(len(forkpoints)), c='r', marker=3)
                plt.scatter(Xi, Yi, c='r', marker='.')                      # Exploración de la minimización del intervalo
                plt.scatter(Xj, Yj, c='b', marker='.', linewidths=0.01)     # Exploración de la minimización de varphi
                #plt.scatter(minvarphiX, minvarphi, c='c', marker='*', linewidths=2)     # Mínimo de varphi
                #plt.hlines(varphiLimInf, 10 ** -2, intervalsup)        # Límite en infinito
                #plt.hlines(varphiLim0, intervalinf, 10)                # Límite en 0
                plt.show()
                #plt.savefig("galaxies/graphics/" + p + '-' + i + '-interval-improvement.png')
                #plt.gcf().clear()
                #generate3Dgraphic(i, minvarphiX, minvarphi, minrho, galaxdata)     # Generador de gráficas 3D para ejemplos
        pend = time.time()
        print("Tiempo para el perfil ", p, " para la galaxia ", i, " = ", pend - pstart, " segundos")
    iend = time.time()
    print("Tiempo para la galaxia ", i, " = ", iend - istart, " segundos")
endtime = time.time()
print("Tiempo total = ", endtime - start_time, " segundos")
