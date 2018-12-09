import random
import numpy as np
from commonFunctions import phi
random.seed(1)

## VARPHI MINIMIZATION ##


##  Dados los extremos del intervalo, calcula el valor de varphi en el punto medio y el valor de varphi de dos puntos
##  escogidos aleatoriamente, uno a la izquierda del punto medio y otro a la derecha.
#   param:
#       intizq: extremo inferior del intervalo
#       intder: extremo superior del intervalo
#       galaxdata: datos de la galaxia
#   return:
#       un vector de tuplas:
#           m: punto medio
#           M: varphi(m)
#           rhoM: el valor de rho_0 de la Proposición 1 para s = m
#           i: punto aleatorio a la izquierda de m
#           I: varphi(i)
#           rhoI: el valor de rho_0 de la Proposición 1 para s = i
#           d: punto aleatorio a la derecha de m
#           D: varphi(d)
#           rhoD: el valor de rho_0 de la Proposición 1 para s = d
def getIMD(intizq, intder, galaxdata):
    m = (intder + intizq) / 2
    M, rhoM = phi(np.array([m]), galaxdata)
    i = random.uniform(intizq, m)
    I, rhoI = phi(np.array([i]), galaxdata)
    d = random.uniform(m, intder)
    D, rhoD = phi(np.array([d]), galaxdata)

    return [[i, I, rhoI], [m, M, rhoM], [d, D, rhoD]]


##  Realiza la mejora propuesta en la memoria.
#   param:
#       varphiLim0: límite de varphi cuando s tiende a 0
#       varphiLimInf: límite de varphi cuando s tiende a infinito
#       intinfmin: punto mínimo encontrado en la exploración del intervalo inferior
#       intsupmin: punto mínimo encontrado en la exploración del intervalo superior
#       intervalinf: extremo inferior del intervalo calculado en intervalMinim.py
#       intervalsup: extremo superior del intervalo calculado en intervalMinim.py
#   return:
#       intinf: nuevo extremo inferior del intervalo
#       intsup: nuevo extremo superior del intervalo
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


##  Realiza la minimización de varphi.
#   param:
#       varphiLim0: límite de varphi cuando s tiende a 0
#       varphiLimInf: límite de varphi cuando s tiende a infinito
#       intinfmin: punto mínimo encontrado en la exploración del intervalo inferior
#       intsupmin: punto mínimo encontrado en la exploración del intervalo superior
#       intervalinf: extremo inferior del intervalo calculado en intervalMinim.py
#       intervalsup: extremo superior del intervalo calculado en intervalMinim.py
#       galaxdata: datos de la galaxia
#   return:
#       vector solución
#           bestphi: devuelve el valor mínimo encontrado para varphi(s)
#           bestrho: devuelve el valor de rho_0(s) tal que s es bestphiX
#           bestphiX: devuelve el valor de s tal que varphi(s) es bestphi
#           forkpoints: devuelve los puntos m en situación de fork
#           Xs: división del intervalo inicial, sin los subintervalos añadidos de situaciones fork
#           intervalinf: extremo inferior del intervalo calculado en intervalMinim.py
#           intervalsup: extremo superior del intervalo calculado en intervalMinim.py
def varphiMin(varphiLim0, varphiLimInf, intinfmin, intsupmin, intervalinf, intervalsup, galaxdata):
    tol = 10**-8
    intervalinf, intervalsup = reductionInterval(varphiLim0, varphiLimInf,
                                                 intinfmin, intsupmin, intervalinf, intervalsup)    # Mejora propuesta
    subint = np.asarray(np.logspace(np.log10(intervalinf), np.log10(intervalsup), 8))
    Xs = subint
    bestphi = 10**4
    if galaxdata["graphic"]:
        X = []
        Y = []
    s = 0
    nfork = 0
    cfork = 0
    forkpoints = []
    while s < len(subint) - 1:
        intizq = subint[s]
        intder = subint[s+1]
        if cfork == nfork + 1:
            nfork = 0
            cfork = 0
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
        if minphi < bestphi:
            bestphi = minphi
            bestrho = rho
            bestphiX = minphiX
        s += 1
        cfork += 1
        if galaxdata["graphic"]:
            sol = [bestphi, bestrho, bestphiX, X, Y, forkpoints, Xs, intervalinf, intervalsup]
        else:
            sol = [bestphi, bestrho, bestphiX, forkpoints, Xs, intervalinf, intervalsup]
    return sol
