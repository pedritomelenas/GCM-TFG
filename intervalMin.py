import random
import numpy as np
from redMethRotCurveFitting import phi
import time
random.seed(10)
s0 = [random.uniform(0, 10000)]
bestphi = phi(np.asarray(s0))
bestsols = []
bestsols.append(s0[0])
close = True
M = 10000
m = 0
bests0 = 0
S = 10
dist = 100
k = 0
n = 0
N = 5000

start_time = time.time()
while abs(M - m) >= 10.0:
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
        if phi0 < bestphi:
            #print("mejor")
            bestphi = phi0
            bestsols.append(s0[0])
            bests0 = s0[0]
            #print(bestsols)
            #print(len(bestsols))
            n = 0
        else:
            n += 1
            if n > N and dist/10 > 0.0 and S - 1 > 1:
                dist /= 10
                S -= 1
                print("S_ = ", S)
    print(time.time() - second_start, "SEGUNDOS EN LLENAR bestsols CON ", len(bestsols), " VALORES EN LA ITERACIÓN k = ", k)
    media = sum(bestsols)/len(bestsols)
    print("media = ", media)
    d = 0.0
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
    N *= 2
    rsol = (bestsols[0] + bestsols[1]) / 2  #random.choice(bestsols)
    bestsols.clear()
    print("rsol = ", rsol)
    bestsols.append(rsol)
    bestphi = phi(np.asarray([rsol]))
    if dist/10 > 0.0:
        dist /= 10
    #print("dist = ", dist)
    k += 1
    print("S = ", S)
    #print("dist = ", dist)

print(time.time() - start_time, " SEGUNDOS")
print("interval = [", m, ", ", M, "]")
print("interval distance = ", abs(m-M))
print("bestphi = ", bestphi)
print("bests0 = ", bests0)