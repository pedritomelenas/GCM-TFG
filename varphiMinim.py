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
    #return {'I': I, 'i': i, 'M': M, 'm': m, 'D': D, 'd': d}

def varphiMin(intervalinf, intervalsup, galaxdata):
    tol = 10**-8
    subint = np.asarray(np.logspace(np.log10(intervalinf), np.log10(intervalsup), 8))
    initialsubint = subint
    bestphi = 10**4
    X = []
    Y = []
    subintsize = np.arange(len(subint) - 1)
    s = 0
    #step = 1
    nfork = 0
    forkpoints = []
    while s < len(subint) - 1:
        #print("subintervalo ", s)
        intizq = subint[s]
        intder = subint[s+1]
        if s == nfork + 1:
            nfork = 0
        M = 1
        lastM = 0
        minphi = 10**4
        k = 0
        #fork = False
        while abs(M - lastM) > tol:
            #print(k)
            lastM = M
            #print("intizq = ", intizq, "; intder = ", intder)
            IMD = getIMD(intizq, intder, galaxdata)
            m = IMD[1][0]
            M = IMD[1][1]
            #print("M = ", M)
            i = IMD[0][0]
            I = IMD[0][1]
            #print("I = ", I)
            d = IMD[2][0]
            D = IMD[2][1]
            #print("D = ", D)
            Y.append(M)
            Y.append(I)
            Y.append(D)
            X.append(m)
            X.append(i)
            X.append(d)
            '''if s_i >= 2:
                print("[m, M] = [", m, ", ", M, "]")
                print("[i, I] = [", i, ", ", I, "]")
                print("[d, D] = [", d, ", ", D, "]")'''
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
                    #fork = True
                    intder = m
                    #print("fork = ", fork)
                    #print("SUBINT LAST = ", subint)
                    #print("s+1 = ", s+1)
                    subint = np.insert(subint, s+1, m)
                    forkpoints.append(m)
                    nfork += 1
                    #step += 1
                    #print("nfork = ", nfork)
                    #print("SUBINT NEW = ", subint)
                    #print("len(subint) = ", len(subint))
            #print("min = ", min([I, M, D]))
            sorted_IMD = sorted(IMD, key=lambda tup: tup[1])
            #print("SORTED_IMD = ", sorted_IMD)
            if sorted_IMD[0][1] < minphi:
                minphi = sorted_IMD[0][1]
                minphiX = sorted_IMD[0][0]
            k += 1
            #print("intervalo = [", intizq, ", ", intder, "]")
        if minphi < bestphi:
            bestphi = minphi
            bestphiX = minphiX
        #print("FORK = ", fork)
        #print("BESTPHI = ", bestphi)
        s += 1
    return [bestphi, bestphiX, X, Y, forkpoints]
