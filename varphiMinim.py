import random
import numpy as np
from commonFunctions import phi
from redMethRotCurveFitting import intervalinf, intervalsup, galaxdata
'''
def getIMD(intizq, intder, galaxdata):
    m = (intder + intizq) / 2
    M = phi(np.array([m]), 'ISO', galaxdata)
    i = random.uniform(intizq, m)
    I = phi(np.array([i]), 'ISO', galaxdata)
    d = random.uniform(m, intder)
    D = phi(np.array([d]), 'ISO', galaxdata)

    return {'I': I, 'M': M, 'D': D}

def forkRoutine(subint1, subint2, bestphi, tol, galaxdata):
    bestphi1 = subintervalRoutine(subint1, 0, bestphi, tol, False, galaxdata)
    bestphi2 = subintervalRoutine(subint2, 0, bestphi, tol, False, galaxdata)
    bestphi = min(bestphi1, bestphi2)

    return bestphi

def subintervalRoutine(subint, i, bestphi, tol, fork, galaxdata):
    print("SUBINT = ", i)
    intizq = subint[i]
    intder = subint[i + 1]
    M = 1
    lastM = 0
    minphi = 10 ** 4
    k = 0
    while abs(M - lastM) > tol:
        print(k)
        lastM = M
        print("intizq = ", intizq, "; intder = ", intder)
        IMD = getIMD(intizq, intder, galaxdata)
        M = IMD['M']
        print("M = ", M)
        I = IMD['I']
        print("I = ", I)
        D = IMD['D']
        print("D = ", D)
        if I < M < D:
            intder = m
        elif I > M > D:
            intizq = m
        else:
            if fork:
                intizqf = intizq
                intderf = intder
                forkRoutine([intizqf, m], [m, intderf])
            else:
                if I < D:
                    intder = m
                else:
                    intizq = m
        if min(IMD) < minphi:
            minphi = min(IMD)
        k += 1
        # print("intervalo = [", intizq, ", ", intder, "]")
    if minphi < bestphi:
        bestphi = minphi

    return bestphi
'''
#def varphiMin(intervalinf, intervalsup, galaxdata):
random.seed(1)
tol = 10**-2
subint = np.asarray(np.logspace(np.log10(intervalinf), np.log10(intervalsup), 8))
bestphi = 10**4

for s in np.arange(len(subint) - 1):
    print("SUBINT = ", s)
    print("len(SUBINT) = ", len(subint))
    print("subint = ", type(subint))
    intizq = subint[s]
    intder = subint[s + 1]
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
            intder = m
            print("SUBINT LAST = ", subint)
            print("s+1 = ", s+1)
            subint = np.insert(subint, s+1, m)
            print("SUBINT NEW = ", subint)
        if M < minphi:
            minphi = M
        k += 1
        #print("intervalo = [", intizq, ", ", intder, "]")
    if minphi < bestphi:
        bestphi = minphi

    print("BESTPHI = ", bestphi)





#return 0