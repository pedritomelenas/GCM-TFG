import numpy as np
import data as dt
import commonFunctions as cf
from calLimits import calLimits
from intervalMinim import intervalMin
from varphiMinim import varphiMin
import matplotlib.pyplot as plt
import time

'''
Proceso de ajuste de curvas de rotación:
    1) Cálculo de límites
    2) Minimización del intervalo de búsqueda
    3) Minimización de la función varphi
'''

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

start_time = time.time()
for i in dt.galaxlist:
    istart = time.time()
    print("\n")
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
    # galaxdata["graphic"] = True

    for p in dt.profiles:
        print(" ********** PROFILE: ", p, " ************")
        galaxdata["profile"] = p
        pstart = time.time()
        if n - dt.nu > 0:
            ### Cálculo de límites ###
            limits = calLimits(galaxdata)
            varphiLim0 = limits[0]
            varphiLimInf = limits[1]

            ### Minimización del intervalo de búsqueda ###
            intervalmin_start = time.time()
            interval = intervalMin(varphiLim0, varphiLimInf, galaxdata)
            intervalmin_end = time.time()
            print("Tiempo de minimización del intervalo = ", intervalmin_end - intervalmin_start, " segundos")
            intervalinf = interval[0][0]
            intervalsup = interval[0][1]
            print("Intervalo de búsqueda inicial = [", intervalinf, ", ", intervalsup, "]")
            if galaxdata["graphic"]:
                Xi = interval[1]
                Yi = interval[2]
                intinfmin = interval[3]
                intsupmin = interval[4]
            else:
                intinfmin = interval[1]
                intsupmin = interval[2]

            ### Minimización de la función varphi ###
            varphimin_start = time.time()
            pmin = varphiMin(varphiLim0, varphiLimInf, intinfmin, intsupmin, intervalinf, intervalsup, galaxdata)
            varphimin_end = time.time()
            print("Tiempo de minimización de varphi = ", varphimin_end - varphimin_start, " segundos")
            minvarphi = pmin[0]
            minrho = pmin[1]
            minvarphiX = pmin[2]
            old_intervalsup = intervalsup
            if galaxdata["graphic"]:
                Xj = pmin[3]
                Yj = pmin[4]
                forkpoints = pmin[5]
                X = pmin[6]
                intervalinf = pmin[7]
                intervalsup = pmin[8]
            else:
                forkpoints = pmin[3]
                X = pmin[4]
                intervalinf = pmin[5]
                intervalsup = pmin[6]
            print("Nuevo intervalo de búsqueda = [", intervalinf, ", ", intervalsup, "]")
            print("Mínimo de varphi = ", minvarphi)
            print("para s = ", minvarphiX)
            print("con rho(s) = ", minrho)
            if galaxdata["graphic"]:
                plt.semilogx()
                plt.title("Galaxia " + i + " con perfil " + p)
                plt.xlabel("s (parámetro de escala)")
                plt.ylabel(r"$\varphi(s)$")

                ## Gráfica de ejemplo ##
                # Xe = np.logspace(-5, 3, 5000)
                # phiX, rho = cf.phi(Xe, galaxdata)
                # plt.plot(Xe, phiX, 'k')

                ## Ejemplo: diferencia entre condiciones de convergencia ##
                # plt.vlines(intervalsup, minvarphi, varphiLimInf, colors='b', linestyles='dashed')
                # plt.gca().set_xticks([intervalsup])
                # plt.gca().set_xticklabels([""+str(intervalsup)+""])

                plt.scatter(intervalinf, 0, c='black', marker=3)                        # Extremo inferior del intervalo
                plt.scatter(intervalsup, 0, c='black', marker=3)                        # Extremo superior del intervalo
                plt.hlines(0, intervalinf, intervalsup)                                 # Intervalo
                plt.scatter(X, np.zeros(len(X)), color='black', marker=3)               # División del intervalo
                # plt.scatter(forkpoints, np.zeros(len(forkpoints)), c='r', marker=3)   # Puntos en situación de fork

                ## Exploración de la minimización del intervalo ##
                plt.scatter(Xi, Yi, c='r', marker='.')
                ## Exploración de la minimización de varphi ##
                plt.scatter(Xj, Yj, c='b', marker='.', linewidths=0.01)
                plt.scatter(minvarphiX, minvarphi, c='yellow', marker='.', linewidths=2)            # Mínimo de varphi
                plt.hlines(varphiLimInf, 0, old_intervalsup, colors='g', linestyles='dotted')       # Límite en infinito
                plt.hlines(varphiLim0, 0, intervalsup, colors='g', linestyles='dotted')             # Límite en 0
                plt.show()
                #plt.savefig("galaxies/graphics/" + p + '-' + i + '-all-improvements.png')
                #plt.gcf().clear()
                # generate3Dgraphic(i, minvarphiX, minvarphi, minrho, galaxdata)     # Generador de gráficas 3D para ejemplos
        pend = time.time()
        print("Tiempo para el perfil ", p, " para la galaxia ", i, " = ", pend - pstart, " segundos")

        ## Para comprobar si podemos asegurar la existencia del mínimo ##
        #if (abs(minvarphiX - intervalinf) < 10 ** (-12) and (abs(minvarphi - varphiLim0)/varphiLim0 < 10**(-2))) \
        #        or (abs(minvarphiX - intervalsup) < 10 ** (-12) and (abs(minvarphi - varphiLimInf)/varphiLimInf < 10**(-2))):
        #    print("En este caso NO PODEMOS ASEGURAR LA EXISTENCIA DEL MÍNIMO")

    iend = time.time()
    print("Tiempo para la galaxia ", i, " = ", iend - istart, " segundos")
endtime = time.time()
print("Tiempo total = ", endtime - start_time, " segundos")
