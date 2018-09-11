import random
import numpy as np
from redMethRotCurveFitting import phi, varphiLimInf, varphiLim0
import time

tol = 10**-6
maxiter = 0
direction = -1
intervalinf = -3
intervalsup = 3

k = 0
lastint = 0.0
dir = []
stop = False
i = 0.0

random.seed(1)

#eval = (phi(10**(intervalinf + np.array([-0.2, -0.1, 0.0, 0.1, 0.2]))) - varphiLim0) / varphiLim0
#print(eval)
#print(eval[2:5])

while maxiter < 100 and direction != 0 and k < 50:
    maxiter += 1
    eval = abs(phi(10**(intervalinf + np.array([-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.10, 0.15, 0.2]))) - varphiLim0) / varphiLim0
    #print("eval[2] = ", eval[2])
    #print("eval[3] = ", eval[3])
    #print("eval[4] = ", eval[4])
    #print("sum(eval[2:5]) = ", sum(eval[2:5]))
    #print("sum(eval[0:3]) = ", sum(eval[0:3]))
    # [2:5] [0,3]
    # [5:10] [0:5]
    #test1 = sum(eval[2:5]) < tol
    #test2 = sum(eval[0:3]) < tol
    test1 = True                    # En vez de considerar la suma de los elementos de eval, se considera
    test2 = True                    # cada uno de forma independiente
    for e in eval[5:10]:
        if e > tol:
            test1 = False
    for e in eval[0:5]:
        if e > tol:
            test2 = False

    if test1 == 0 and test2 == 1:
        opcion = 1
        stop = True                 # Podemos convertir stop en un contador, en lugar de un booleano, pera
        direction = -1              # establecer un patrón más específico
        i = intervalinf
        intervalinf = intervalinf + random.uniform(0.2, 0.3) * direction        # En cada iteración, intervalinf
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
        print("opcion = ", opcion)
        print("intervalinf = ", intervalinf)
        #dir.append(direction)
    else:
        k = 0
    lastint = 10**intervalinf

print("maxiter = ", maxiter)
print("k = ", k)
#print("dir = ", dir)
print("direction = ", direction)
direction = 1
maxiter = 0
k = 0
lastint = 0.0
dir.clear()
stop = False
i = 0.0

while maxiter < 100 and direction != 0 and k < 50:  # Nueva condición de parada --> relacionar con el límite
    maxiter += 1
    # [-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.10, 0.15, 0.2]
    eval = (phi(10**(intervalsup + np.array([-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.10, 0.15, 0.2]))) - varphiLimInf) / varphiLimInf
    # [2:5] [0,3]  corregir
    # [5:10] [0:5]
    #test1 = sum(eval[0:3]) < tol
    #test2 = sum(eval[2:5]) < tol
    test1 = True        # En vez de considerar la suma de los elementos de eval, se considera
    test2 = True        # cada uno de forma independiente
    for e in eval[5:10]:
        if e > tol:
            test2 = False
    for e in eval[0:5]:
        if e > tol:
            test1 = False
    if test1 == 0 and test2 == 1:
        stop = True         # Podemos convertir stop en un contador, en lugar de un booleano, pera
        direction = 1       # establecer un patrón más específico --> agrandar el boleano: 0.3 - 0.6
        i = intervalsup
        intervalsup = intervalsup + random.uniform(0.2, 0.3) * direction            # En cada iteración, intervalsup
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
        print("intervalsup = ", intervalsup)
        #dir.append(direction)
    else:
        k = 0
    lastint = 10**intervalinf

print("k = ", k)
#print("dir = ", dir)
print("maxiter = ", maxiter)
print("direction = ", direction)

if intervalinf > intervalsup:
    print("The interval selection module has not worked properly")
    intervalinf = -3
    intervalsup = 5

print("[", 10**intervalinf, ", ", 10**intervalsup, "]")
