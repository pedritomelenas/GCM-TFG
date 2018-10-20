import random
import numpy as np
from commonFunctions import phi
random.seed(1)

def getIMD(intizq, intder, galaxdata):
    m = (intder + intizq) / 2
    M = phi(np.array([m]), galaxdata)
    i = random.uniform(intizq, m)
    I = phi(np.array([i]), galaxdata)
    d = random.uniform(m, intder)
    D = phi(np.array([d]), galaxdata)

    return {'I': I, 'i': i, 'M': M, 'm': m, 'D': D, 'd': d}

def varphiMin(intervalinf, intervalsup, galaxdata):
    tol = 10**-8
    subint = np.asarray(np.logspace(np.log10(intervalinf), np.log10(intervalsup), 8))
    bestphi = 10**4

    for s in np.arange(len(subint) - 1):
        #print("SUBINT = ", s)
        #print("len(SUBINT) = ", len(subint))
        #print("subint = ", type(subint))
        intizq = subint[s]
        intder = subint[s + 1]
        M = 1
        lastM = 0
        minphi = 10**4
        k = 0
        fork = 0
        while abs(M - lastM) > tol:
            #print(k)
            lastM = M
            #print("intizq = ", intizq, "; intder = ", intder)
            IMD = getIMD(intizq, intder, galaxdata)
            m = IMD['m']
            M = IMD['M']
            #print("M = ", M)
            i = IMD['i']
            I = IMD['I']
            #print("I = ", I)
            d = IMD['d']
            D = IMD['D']
            #print("D = ", D)
            if I < M < D:
                intder = m
            elif I > M > D:
                intizq = m
            else:
                fork += 1
                intder = m
                #print("SUBINT LAST = ", subint)
                #print("s+1 = ", s+1)
                subint = np.insert(subint, s+1, m)
                subint = np.inser(subint, s+2, d)       ## sin esto no hace el fork bien! Falta probar
                #print("SUBINT NEW = ", subint)
            #print("min = ", min([I, M, D]))
            if min([I, M, D]) < minphi:
                minphi = min([I, M, D])
            #print("minphi = ", minphi)
            k += 1
            #print("intervalo = [", intizq, ", ", intder, "]")
        if minphi < bestphi:
            bestphi = minphi
        #print("FORK = ", fork)
        #print("BESTPHI = ", bestphi)
    return bestphi