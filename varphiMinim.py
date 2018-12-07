import random
import numpy as np
from commonFunctions import phi
random.seed(1)

## VARPHI MINIMIZATION ##

def getIMD(intizq, intder, galaxdata):
    m = (intder + intizq) / 2
    M, rhoM = phi(np.array([m]), galaxdata)
    i = random.uniform(intizq, m)
    I, rhoI = phi(np.array([i]), galaxdata)
    d = random.uniform(m, intder)
    D, rhoD = phi(np.array([d]), galaxdata)

    return [[i, I, rhoI], [m, M, rhoM], [d, D, rhoD]]

def reductionInterval(varphiLim0, varphiLimInf, intinfmin, intsupmin, intervalinf, intervalsup):
    intinf = intervalinf
    if varphiLim0 > varphiLimInf:
        if abs(varphiLim0 - intinfmin[1]) < abs(varphiLimInf - intinfmin[1]):
            intinf = intinfmin[0]
    else:
        if abs(varphiLim0 - intinfmin[1]) < 0.1:
            intinf = intinfmin[0]
    intsup = intervalsup
    if varphiLimInf > varphiLim0:
        if abs(varphiLimInf - intsupmin[1]) < abs(varphiLim0 - intsupmin[1]):
            intsup = intsupmin[0]
    else:
        if abs(varphiLimInf - intsupmin[1]) < 0.1:
            intsup = intsupmin[0]
    return intinf, intsup

def varphiMin(varphiLim0, varphiLimInf, intinfmin, intsupmin, intervalinf, intervalsup, galaxdata):
    tol = 10**-8
    intervalinf, intervalsup = reductionInterval(varphiLim0, varphiLimInf,
                                                 intinfmin, intsupmin, intervalinf, intervalsup)    # Mejora
    subint = np.asarray(np.logspace(np.log10(intervalinf), np.log10(intervalsup), 8))
    Xs = subint
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
        izqphi, rhoi = phi(np.array([intizq]), galaxdata)
        derphi, rhod = phi(np.array([intder]), galaxdata)
        if izqphi < derphi:
            minphi = izqphi
            rho = rhoi
            minphiX = intizq
        else:
            minphi = derphi
            rho = rhod
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
                rho = sorted_IMD[0][2]
                minphiX = sorted_IMD[0][0]
            k += 1
        if minphi < bestphi:
            bestphi = minphi
            bestrho = rho
            bestphiX = minphiX
        s += 1
        if galaxdata["graphic"]:
            sol = [bestphi, bestrho, bestphiX, X, Y, forkpoints, Xs, intervalinf, intervalsup]
        else:
            sol = [bestphi, bestrho, bestphiX, forkpoints, Xs, intervalinf, intervalsup]
    return sol
