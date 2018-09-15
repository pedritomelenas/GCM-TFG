import random
import numpy as np
from commonFunctions import phi
from redMethRotCurveFitting import intervalinf, intervalsup, galaxdata

#def varphiMin(intervalinf, intervalsup, galaxdata):
random.seed(1)
tol = 10**-2
subint = np.logspace(np.log10(intervalinf), np.log10(intervalsup), 8)
bestphi = 10**4
for i in np.arange(len(subint) - 1):
    intizq = subint[i]
    intder = subint[i + 1]
    M = 1
    lastM = 0
    minphi = 10**4
    k = 0
    while abs(M - lastM) > tol:
        print(k)
        lastM = M
        print("intizq = ", intizq, "; intder = ", intder)
        m = (intder + intizq) / 2
        print("m = ", m)
        M = phi(np.array([m]), 'ISO', galaxdata)
        print("M = ", M)
        i = random.uniform(intizq, m)
        print("i = ", i)
        I = phi(np.array([i]), 'ISO', galaxdata)
        print("I = ", I)
        d = random.uniform(m, intder)
        print("d = ", d)
        D = phi(np.array([d]), 'ISO', galaxdata)
        print("D = ", D)
        if I < M < D:
            intder = m
        elif I > M > D:
            intizq = m
        else:
            print("FORK")
        if M < minphi:
            minphi = M
        k += 1
        #print("intervalo = [", intizq, ", ", intder, "]")
    if minphi < bestphi:
        bestphi = minphi

    print("BESTPHI = ", bestphi)





#return 0