import random
import numpy as np
from redMethRotCurveFitting import phi

random.seed(10)
s0 = [random.uniform(0, 10000)]
bestphi = phi(np.asarray(s0))
bestsols = []
bestsols.append(s0[0])
close = True
M = 5000
S = 10
dist = 1000
k = 0
n = 0

while k < 5 and S > 0:
    print("k = ", k)
    print("n = ", n)
    n = 0
    while len(bestsols) < S and n < 5000:
        while close:
            s0[0] = random.uniform(0, 2*M)
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
            print("mejor")
            bestphi = phi0
            bestsols.append(s0[0])
            bests0 = s0[0]
            print(bestsols)
            print(len(bestsols))
            n = 0
        else:
            n += 1
            if n > 500 and dist/10 > 0.0:
                dist /= 10

    M = sum(bestsols)/len(bestsols)
    if S - 2 > 0:
        S -= 2
    #sols.clear()
    bestsols.clear()
    #sols.append(bests0)
    bestsols.append(bests0)
    if dist/10 > 0.0:
        dist /= 10
    print("dist = ", dist)
    k += 1
    print("M =", M)
    print("S = ", S)
    #print("dist = ", dist)
print("interval = ", 2*M)
print("bestphi = ", bestphi)
print("bests0 = ", bests0)