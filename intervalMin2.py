import random
import numpy as np
from redMethRotCurveFitting import phi, varphiLimInf, varphiLim0
import time
import matplotlib.pyplot as plt

tol = 10**-2
maxiter = 0
direction = -1
intervalinf = -3
intervalsup = 3

k = 0
lastint = 0.0
dir = []
stop = False
i = 0.0
X = []
Y = []
random.seed(1)

#eval = (phi(10**(intervalinf + np.array([-0.2, -0.1, 0.0, 0.1, 0.2]))) - varphiLim0) / varphiLim0
#print(eval)
#print(eval[2:5])

while maxiter < 100 and direction != 0 and k < 50:
    maxiter += 1
    varphi = phi(10**(intervalinf + np.array([-0.2, -0.1, 0.0, 0.1, 0.2])))
    eval = abs(varphi - varphiLim0) / varphiLim0
    X.append(10**(intervalinf + np.array([-0.2, -0.1, 0.0, 0.1, 0.2])))
    Y.append(varphi)
    # [2:5] [0,2]
    # [5:10] [0:5]
    #test1 = sum(eval[2:5]) < tol
    #test2 = sum(eval[0:3]) < tol
    test1 = True                    # En vez de considerar la suma de los elementos de eval, se considera
    test2 = True                    # cada uno de forma independiente
    for e in eval[2:5]:
        if e > tol:
            test1 = False
    for e in eval[0:2]:
        if e > tol:
            test2 = False

    if test1 == 0 and test2 == 1:
        opcion = 1
        stop = True                 # Podemos convertir stop en un contador, en lugar de un booleano, para
        direction = -1              # establecer un patrón más específico
        i = intervalinf
        intervalinf = intervalinf + random.uniform(0.2, 0.6) * direction        # En cada iteración, intervalinf
        #direction = 0                                                          # varía con una dirección determinada
    elif test1 == 0 and test2 == 0:                                             # pero con módulo aleatorio
        opcion = 2
        if stop:
            stop = False
        direction = -1
        intervalinf = intervalinf + random.uniform(0.2, 0.3) * direction
        #intervalinf = intervalinf + 0.3 * direction
    elif test1 == 1 and test2 == 1:
        opcion = 3
        if stop:
            direction = 0
            intervalinf = i
            diferencia = abs(phi(np.asarray([10 ** intervalinf])) - varphiLim0)
        else:
            direction = 1
            intervalinf = intervalinf + random.uniform(0.2, 0.3) * direction
        #intervalinf = intervalinf + 0.3 * direction
    else:
        opcion = 4
        if stop:
            stop = False
        direction = -1
        intervalinf = intervalinf + random.uniform(0.2, 0.3) * direction
        #intervalinf = intervalinf + 0.3 * direction
    if abs(10**intervalinf - lastint) < tol:
        k += 1
        if k >= 25 and abs(phi(np.asarray([10**intervalsup])) - varphiLimInf) >= 1:
            k = 0
        #print("opcion = ", opcion)
        #print("intervalinf = ", intervalinf)
        #dir.append(direction)
    else:
        k = 0
    lastint = 10**intervalinf
    dir.append(direction)

print("maxiter = ", maxiter)
print("k = ", k)
print("dir = ", dir)
print("direction = ", direction)
print("diferencia = ", diferencia)
direction = 1
maxiter = 0
k = 0
lastint = 0.0
dir.clear()
stop = False
i = 0.0
dir.clear()

## IDEA: APROVECHAR ESTE ALGORITMO PARA DETECTAR, SI ES QUE SE DA, LA EXISTENCIA DEL MÍNIMO ##
    ## Si varphiLim0 < varphiLimInf, basta con que exista un s tal que varphi(s) < varphiLim0 para asegurar la existencia del límite
    ## Si varphiLimInf < varphiLim0, basta con que exista un s tal que varphi(s) < varphiLimInf para asegurar la existencia del límite

while maxiter < 100 and direction != 0 and k < 50:  # Nueva condición de parada --> relacionar con el límite
    maxiter += 1
    # [-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.10, 0.15, 0.2]
    varphi = phi(10**(intervalsup + np.array([-0.2, -0.1, 0.0, 0.1, 0.2])))
    eval = abs(varphi - varphiLimInf) / varphiLimInf    ##  CORREGIDO
    X.append(10**(intervalsup + np.array([-0.2, -0.1, 0.0, 0.1, 0.2])))
    Y.append(varphi)
    # [3:5] [0,3]
    # [5:10] [0:5]
    #test1 = sum(eval[0:3]) < tol
    #test2 = sum(eval[2:5]) < tol
    test1 = True        # En vez de considerar la suma de los elementos de eval, se considera
    test2 = True        # cada uno de forma independiente
    for e in eval[3:5]:
        if e > tol:
            test2 = False
    for e in eval[0:3]:
        if e > tol:
            test1 = False
    if test1 == 0 and test2 == 1:
        stop = True         # Podemos convertir stop en un contador, en lugar de un booleano, pera
        direction = 1       # establecer un patrón más específico --> agrandar el boleano: 0.3 - 0.6
        i = intervalsup
        intervalsup = intervalsup + random.uniform(0.2, 0.6) * direction            # En cada iteración, intervalsup
        # direction = 0                                                             # varía con una dirección determinada
    elif test1 == 0 and test2 == 0:
        if stop:
            stop = False
        direction = 1
        intervalsup = intervalsup + random.uniform(0.2, 0.3) * direction
    elif test1 == 1 and test2 == 1:
        if stop:
            direction = 0
            intervalsup = i
            diferencia = abs(phi(np.asarray([10**intervalsup])) - varphiLimInf)
        else:
            direction = -1
            intervalsup = intervalsup + random.uniform(0.2, 0.3) * direction
    else:
        if stop:
            stop = False
        direction = 1
        intervalsup = intervalsup + random.uniform(0.2, 0.3) * direction
    if abs(10**intervalsup - lastint) < tol:
        k += 1
        if k >= 25 and abs(phi(np.asarray([10**intervalsup])) - varphiLimInf) >= 1:
            k = 0
        #print("intervalsup = ", intervalsup)
        #dir.append(direction)
    else:
        k = 0
    lastint = 10**intervalinf
    dir.append(direction)

print("k = ", k)
print("diferencia = ", diferencia)
#print("dir = ", dir)
print("maxiter = ", maxiter)
print("direction = ", direction)
print("dir = ", dir)
Xa = np.asarray(X)
Ya = np.asarray(Y)
plt.scatter(X, Y)
plt.hlines(varphiLimInf, 0, 1000)
plt.hlines(varphiLim0, 0, 100)
plt.vlines(10**intervalsup, 0, 50)
plt.vlines(10**intervalinf, 50, 70)
plt.show()

if intervalinf > intervalsup:
    print("The interval selection module has not worked properly")
    intervalinf = -3
    intervalsup = 5

print("[", 10**intervalinf, ", ", 10**intervalsup, "]")
