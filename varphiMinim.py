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

    return [[i, I], [m, M], [d, D]]

def varphiMin(intervalinf, intervalsup, galaxdata):
    tol = 10**-8
    subint = np.asarray(np.logspace(np.log10(intervalinf), np.log10(intervalsup), 8))
    bestphi = 10**4
    if galaxdata["graphic"]:
        X = []
        Y = []
    s = 0
    nfork = 0
    forkpoints = []
    while s < len(subint) - 1:
        intizq = subint[s]
        intder = subint[s+1]
        if s == nfork + 1:
            nfork = 0
        M = 1
        lastM = 0
        izqphi = phi(np.array([intizq]), galaxdata)
        derphi = phi(np.array([intder]), galaxdata)
        if izqphi < derphi:
            minphi = izqphi
            minphiX = intizq
        else:
            minphi = derphi
            minphiX = intder
        k = 0
        while abs(M - lastM) > tol:
            lastM = M
            IMD = getIMD(intizq, intder, galaxdata)
            m = IMD[1][0]
            M = IMD[1][1]
            I = IMD[0][1]
            D = IMD[2][1]
            if galaxdata["graphic"]:
                Y.append(M)
                Y.append(I)
                Y.append(D)
                X.append(m)
                X.append(IMD[0][0])
                X.append(IMD[2][0])
            if I < M < D:
                intder = m
            elif I > M > D:
                intizq = m
            else:
                if nfork + 1 > 3:
                    if random.randint(0, 1):
                        intder = m
                    else:
                        intizq = m

                else:
                    intder = m
                    subint = np.insert(subint, s+1, m)
                    forkpoints.append(m)
                    nfork += 1
            sorted_IMD = sorted(IMD, key=lambda tup: tup[1])
            if sorted_IMD[0][1] < minphi:
                minphi = sorted_IMD[0][1]
                minphiX = sorted_IMD[0][0]
            k += 1
        if minphi < bestphi:
            bestphi = minphi
            bestphiX = minphiX
        s += 1
        if galaxdata["graphic"]:
            sol = [bestphi, bestphiX, X, Y, forkpoints]
        else:
            sol = [bestphi, bestphiX, forkpoints]
    return sol
