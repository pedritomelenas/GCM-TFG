import random
import numpy as np
from math import trunc
from commonFunctions import phi

tol = 10 ** -2

def inftestElementwise(eval):
    test1 = True
    test2 = True
    l = len(eval)
    for e in eval[trunc(l/2):l]:
        if e > tol:
            test1 = False
    for e in eval[0:trunc(l/2)]:
        if e > tol:
            test2 = False
    return test1, test2

def suptestElementwise(eval):
    test1 = True
    test2 = True
    l = len(eval)
    for e in eval[round(l/2):l]:
        if e > tol:
            test2 = False
    for e in eval[0:round(l/2)]:
        if e > tol:
            test1 = False
    return test1, test2

def inftestElementsum(eval):
    l = len(eval)
    test1 = sum(eval[trunc(l/2):l]) < tol
    test2 = sum(eval[0:trunc(l/2)]) < tol
    return test1, test2

def suptestElementsum(eval):
    l = len(eval)
    test1 = sum(eval[0:round(l/2)]) < tol
    test2 = sum(eval[round(l/2):l]) < tol
    return test1, test2

def infConditions(test1, test2, intervalinf, stop, i):
    if (not test1) and test2:
        stop = True  # Podemos convertir stop en un contador, en lugar de un booleano, para
        direction = -1  # establecer un patrón más específico
        i = intervalinf
        intervalinf = intervalinf + random.uniform(0.2, 0.6) * direction  # En cada iteración, intervalinf
        # direction = 0                                                          # varía con una dirección determinada
    elif (not test1) and (not test2):  # pero con módulo aleatorio
        if stop:
            stop = False
        direction = -1
        intervalinf = intervalinf + random.uniform(0.2, 0.3) * direction
        # intervalinf = intervalinf + 0.3 * direction
    elif test1 and test2:
        if stop:
            direction = 0
            intervalinf = i
        else:
            direction = 1
            intervalinf = intervalinf + random.uniform(0.2, 0.3) * direction
    else:
        if stop:
            stop = False
        direction = -1
        intervalinf = intervalinf + random.uniform(0.2, 0.3) * direction

    return intervalinf, direction, stop, i

def supConditions(test1, test2, intervalsup, stop,i):
    if (not test1) and test2:
        stop = True  # Podemos convertir stop en un contador, en lugar de un booleano, pera
        direction = 1  # establecer un patrón más específico --> agrandar el boleano: 0.3 - 0.6
        i = intervalsup
        intervalsup = intervalsup + random.uniform(0.2, 0.6) * direction  # En cada iteración, intervalsup
        #  varía con una dirección determinada
    elif (not test1) and (not test2):
        if stop:
            stop = False
        direction = 1
        intervalsup = intervalsup + random.uniform(0.2, 0.3) * direction
    elif test1 and test2:
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

    return intervalsup, direction, stop, i

def jumpCondition(twoclosevar, varLimdistance, interval, direction, k):
    jump = False
    if twoclosevar:
        k += 1
        if k > 5 and varLimdistance < 1:
            jump = True
            interval = interval + random.uniform(0.7, 0.8) * direction
            k = 0
        if k >= 25 and varLimdistance >= 1:
            k = 0
    else:
        k = 0
    return jump, interval, k

## IDEA:
## 1) APROVECHAR ESTE ALGORITMO PARA DETECTAR, SI ES QUE SE DA, LA EXISTENCIA DEL MÍNIMO ##
## Si varphiLim0 < varphiLimInf, basta con que exista un s tal que varphi(s) < varphiLim0 para asegurar la existencia del límite
## Si varphiLimInf < varphiLim0, basta con que exista un s tal que varphi(s) < varphiLimInf para asegurar la existencia del límite
## 2) APROVECHAR ESTE ALGORITMO PARA ALMACENAR EL MÍNIMO EVALUADO HASTA EL MOMENTO <-- Hecho (quizás no haya merecido la pena => Revisar)
def intervalMin(varphiLim0, varphiLimInf, galaxdata):
    if galaxdata["graphic"]:
        X = []
        Y = []
    minphi = 10**4  # para devolver el mínimo encontrado en la exploración
    minx = 0        # para devolver el mínimo encontrado en la exploración
    random.seed(1)

    # INTERVALO INFERIOR #
    maxiter = 0
    direction = -1
    intervalinf = -3
    k = 0
    lastint, rho = phi(np.array([float(10 ** intervalinf)]), galaxdata)
    dir = []
    stop = False
    i = 0.0

    while maxiter < 100 and direction != 0 and k < 50:
        maxiter += 1
        s = 10**(intervalinf + np.array([-0.2, -0.1, 0.0, 0.1, 0.2]))
        varphi, rho = phi(s, galaxdata)
        if min(varphi) < minphi:
            minphi = min(varphi)
            pos = (varphi.tolist()).index(minphi)
            minx = s[pos]
        if galaxdata["graphic"]:
            X.append(10**(intervalinf + np.array([-0.2, -0.1, 0.0, 0.1, 0.2])))
            Y.append(varphi)

        eval = abs(varphi - varphiLim0) / varphiLim0
        test1, test2 = inftestElementwise(eval)
        intervalinf, direction, stop, i = infConditions(test1, test2, intervalinf, stop, i)
        var, rho = phi(np.asarray([10 ** intervalinf]), galaxdata)
        twoclosevar = abs(var - lastint) < tol
        varLimdistance = abs(var - varphiLim0)
        jump, intervalinf, k = jumpCondition(twoclosevar, varLimdistance, intervalinf, direction, k)
        if jump:
            lastint, rho = phi(np.asarray([float(10 ** intervalinf)]), galaxdata)
        else:
            lastint = var


    # INTERVALO SUPERIOR #
    maxiter = 0
    direction = 1
    intervalsup = 3
    k = 0
    lastint, rho = phi(np.array([float(10**intervalsup)]), galaxdata)
    dir.clear()
    stop = False
    i = 0.0

    while maxiter < 100 and direction != 0 and k < 50:  # Nueva condición de parada --> relacionar con el límite
        maxiter += 1
        # [-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.10, 0.15, 0.2]
        s = 10**(intervalsup + np.array([-0.2, -0.1, 0.0, 0.1, 0.2]))
        varphi, rho = phi(s, galaxdata)
        if min(varphi) < minphi:
            minphi = min(varphi)
            pos = (varphi.tolist()).index(minphi)
            minx = s[pos]
        if galaxdata["graphic"]:
            X.append(10**(intervalsup + np.array([-0.2, -0.1, 0.0, 0.1, 0.2])))
            Y.append(varphi)

        #eval = abs(varphi - varphiLimInf) / varphiLimInf
        if galaxdata["profile"] == 'ISO' or galaxdata["profile"] == 'BUR':
            eval = abs(varphi - varphiLimInf) / (2*varphiLimInf)
        else:
            eval = abs(varphi - varphiLimInf) / varphiLimInf
        test1, test2 = suptestElementwise(eval)
        intervalsup, direction, stop, i = supConditions(test1, test2, intervalsup, stop, i)
        var, rho = phi(np.asarray([10 ** intervalsup]), galaxdata)
        twoclosevar = abs(var - lastint) < tol
        varLimdistance = abs(var - varphiLimInf)
        jump, intervalsup, k = jumpCondition(twoclosevar, varLimdistance, intervalsup, direction, k)
        if jump:
            lastint, rho = phi(np.asarray([float(10 ** intervalsup)]), galaxdata)
        else:
            lastint = var

    if intervalinf > intervalsup:
        print("WARNING: The interval selection module has not worked properly")
        intervalinf = -3
        intervalsup = 5

    interval = [10**intervalinf, 10**intervalsup]
    if galaxdata["graphic"]:
        sol = [interval, X, Y, minphi, minx]
    else:
        sol = [interval, minphi, minx]

    return sol
