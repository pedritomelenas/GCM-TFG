import numpy as np
import data as dt
import commonFunctions as cf
from calLimits import calLimits
from intervalMinim import intervalMin
from varphiMinim import varphiMin
import matplotlib.pyplot as plt
import time

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
profiles = ['ISO', 'BUR', 'NFW']       # Colocar en data.py?
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
    #galaxdata["graphic"] = True

    for p in profiles:
        print(" ********** PROFILE: ", p, " ************")
        galaxdata["profile"] = p
        pstart = time.time()
        if n - dt.nu > 0:
            ### CALCULATION OF THE LIMITS ###
            limits = calLimits(galaxdata)
            varphiLim0 = limits[0]
            varphiLimInf = limits[1]
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
            print("minvarphi = ", minvarphi)
            if galaxdata["graphic"]:
                plt.semilogx()
                plt.title("Galaxia "+i+" con perfil "+p)
                plt.vlines(intervalinf, -0.1, 0.1)
                plt.vlines(intervalsup, -0.1, 0.1)
                plt.hlines(-0.05, intervalinf, intervalsup)
                plt.scatter(X, np.zeros(len(X)), color='black', marker=3)
                #plt.scatter(forkpoints, np.zeros(len(forkpoints)), c='r', marker=3)
                plt.scatter(Xi, Yi, c='r', marker='.')
                plt.scatter(Xj, Yj, c='b', marker='.', linewidths=0.01)
                #plt.scatter(minvarphiX, minvarphi, c='c', marker='*', linewidths=2)
                #plt.hlines(varphiLimInf, 10 ** -2, intervalsup)
                #plt.hlines(varphiLim0, intervalinf, 10)
                plt.savefig("galaxies/graphics/" + p + '-' + i + '-varphimin.png')
                plt.gcf().clear()
        pend = time.time()
        print("Tiempo para el perfil ", p, " para la galaxia ", i, " = ", pend - pstart, " segundos")
    iend = time.time()
    print("Tiempo para la galaxia ", i, " = ", iend - istart, " segundos")
endtime = time.time()
print("Tiempo total = ", endtime - start_time, " segundos")
