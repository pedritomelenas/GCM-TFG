import random
import numpy as np
from redMethRotCurveFitting import phi, varphiLimInf
import time

random.seed(10)
s0 = [random.uniform(0, 10000)]
minphi = phi(np.asarray(s0))
mins = s0
bestphi = minphi
bests = mins
bestsols = []
bestsols.append(s0[0])
close = True
M = 100000.0
m = 0
bests0 = 0
S = 10
dist = 100
k = 0
n = 0
N = 5000
d = 0.0
farest = 0.0
tol = 10**(-2)
print("s0 = ", s0)
print("phi(np.asarray([M])) = ", phi(np.asarray([M])))
print("varphiLimInf = ", varphiLimInf)
print("abs(phi(np.asarray([M])) - varphiLimInf) / varphiLimInf = ", abs(phi(np.asarray([M])) - varphiLimInf) / varphiLimInf)
start_time = time.time()
while abs(phi(np.asarray([M])) - varphiLimInf) / varphiLimInf >= tol:
    print("k = ", k)
    print("n = ", n)
    n = 0
    second_start = time.time()
    while len(bestsols) < S:    # and n < N:
        while close:
            s0[0] = random.uniform(0, M)
            i = 0
            for s in bestsols:
                if abs(s-s0[0]) <= dist:
                    i += 1
            if i <= int(len(bestsols)/2)+1 or len(bestsols) == 1:
                #sols.append(s0[0])
                close = False
        close = True
        phi0 = phi(np.asarray(s0))
        if phi0 < minphi:
            print("mejor *******")
            minphi = phi0
            print("bestphi = ", bestphi, " para ", s0[0])
            bestsols.append(s0[0])
            mins = s0[0]
            #print(bestsols)
            #print(len(bestsols))
            n = 0
        else:
            n += 1
            if n > N:
                if dist/10 > 0.0:
                    dist /= 10
                if S - 1 > 1:
                    S -= 1
                n = 0
                print("S_ = ", S)
    print(time.time() - second_start, "SEGUNDOS EN LLENAR bestsols CON ", len(bestsols), " VALORES EN LA ITERACIÓN k = ", k)
    if (minphi < bestphi):
        bestphi = minphi
        bests = mins
    media = sum(bestsols)/len(bestsols)
    print("media = ", media)
    farest = max(bestsols)
    for s in bestsols:
        if abs(s - media) >= d:
            farest = s
            d = abs(farest - media)
    print("farest = ", farest)
    M = abs(farest - media) + media
    #if abs(farest - media) == 0.0:
    #        m = 0
    #elif media - abs(farest - media) >= 0:
    #    m = media - abs(farest - media)
    print("bestsols = ", bestsols)
    print("M = ", M)
    print("abs(m - M) = ", abs(M - m))
   # if S - 1 > 1:
    #    S -= 1
    #N *= 2
    rsol = random.uniform(media, M)     #bestsols[0] + 1   #(bestsols[0] + bestsols[1]) / 2  #random.choice(bestsols)
    bestsols.clear()
    print("rsol = ", rsol)
    bestsols.append(rsol)
    bestphi = phi(np.asarray([rsol]))   ## AQUÍ ESTÁ PERDIENDO LA MEJOR SOLUCIÓN
    n = 0
    #if dist/10 > 0.0:
    #    dist /= 10
    #print("dist = ", dist)
    k += 1
    print("S = ", S)
    #print("dist = ", dist)

print(time.time() - start_time, " SEGUNDOS")
print("interval = [", m, ", ", M, "]")
print("interval distance = ", abs(m-M))
print("bestphi = ", bestphi)
print("bests0 = ", bests0)