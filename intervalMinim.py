import random
import numpy as np
from math import trunc
from commonFunctions import phi

tol = 10 ** -2
random.seed(1)


'''
Comprueba si cada uno de los valores del vector pasado por parámetro cumple el
criterio analítico de convergencia para el límite en cero, indicado en (40).
   param:
       eval: vector de puntos vecinos al candidato a extremo inferior del intervalo.
   return:
       test1: booleano, indica si los puntos de la derecha cumplen (40)
       test2: booleano, indica si los puntos de la izquierda cumplen (40)
'''
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


'''
Comprueba si cada uno de los valores del vector pasado por parámetro cumple el
criterio analítico de convergencia para el límite en infinito, indicado en (39).
   param:
       eval: vector de puntos vecinos al candidato a extremo superior del intervalo.
   return:
       test1: booleano, indica si los puntos de la izquierda cumplen (39)
       test2: booleano, indica si los puntos de la derecha cumplen (39)
'''
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


'''
Comprueba si la suma de los valores del vector pasado por parámetro cumple el
criterio analítico de convergencia para el límite en cero, indicado en (40).
   param:
       eval: vector de puntos vecinos al candidato a extremo inferior del intervalo.
   return:
       test1: booleano, indica si la suma de los puntos de la derecha cumplen (40)
       test2: booleano, indica si la suma de los puntos de la izquierda cumplen (40)
'''
def inftestElementsum(eval):
    l = len(eval)
    test1 = sum(eval[trunc(l/2):l]) < tol
    test2 = sum(eval[0:trunc(l/2)]) < tol
    return test1, test2


'''
Comprueba si la suma de los valores del vector pasado por parámetro cumple el
criterio analítico de convergencia para el límite en infinito, indicado en (39).
   param:
       eval: vector de puntos vecinos al candidato a extremo superior del intervalo.
   return:
       test1: booleano, indica si la suma de los puntos de la izquierda cumplen (39)
       test2: booleano, indica si la suma de los puntos de la derecha cumplen (39)
'''
def suptestElementsum(eval):
    l = len(eval)
    test1 = sum(eval[0:round(l/2)]) < tol
    test2 = sum(eval[round(l/2):l]) < tol
    return test1, test2


'''
Comprueba la situación del punto candidato a extremo inferior del intervalo y la
de sus puntos vecinos. Decide si está en condición óptima y hacia qué dirección
moverse.
   param:
       test1: booleano, indica si los puntos de la derecha (o su suma) cumplen (40)
       test2: booleano, indica si los puntos de la izquierda (o su suma) cumplen (40)
       intervalinf: candidato a extremo inferior del intervalo
       stop: controla la condición de parada
       i: almacena el anterior candidato a extremo inferior
   return:
       intervalinf: candidato a extremo inferior del intervalo
       direction: dirección de movimiento (-1, 0 o 1)
       stop: controla la condición de parada
       i: devuelve el anterior candidato a extremo inferior
'''
def infConditions(test1, test2, intervalinf, stop, i):
    if (not test1) and test2:
        stop = True
        direction = -1
        i = intervalinf
        intervalinf = intervalinf + random.uniform(0.2, 0.6) * direction
    elif (not test1) and (not test2):
        if stop:
            stop = False
        direction = -1
        intervalinf = intervalinf + random.uniform(0.2, 0.3) * direction
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


'''
Comprueba la situación del punto candidato a extremo superior del intervalo y la
de sus puntos vecinos. Decide si está en condición óptima y hacia qué dirección
moverse.
   param:
       test1: booleano, indica si los puntos de la izquierda (o su suma) cumplen (39)
       test2: booleano, indica si los puntos de la derecha (o su suma) cumplen (39)
       intervalsup: candidato a extremo superior del intervalo
       stop: controla la condición de parada
       i: almacena el anterior candidato a extremo superior
   return:
       intervalinf: candidato a extremo superior del intervalo
       direction: dirección de movimiento (-1, 0 o 1)
       stop: controla la condición de parada
       i: devuelve el anterior candidato a extremo superior
'''
def supConditions(test1, test2, intervalsup, stop, i):
    if (not test1) and test2:
        stop = True
        direction = 1
        i = intervalsup
        intervalsup = intervalsup + random.uniform(0.2, 0.6) * direction
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


'''
Comprueba si se da la condición de salto.
   param:
       twoclosevar: booleano, indica si los dos últimos candidatos están "cerca"
       varLimdistance: indica a qué distancia está el candidato del valor del límite
       interval: candidato a intervalo inferior o superior
       direction: dirección de movimiento (-1, 0 o 1)
       k: contador de la condición de salto
   return:
       jump: booleano, indica si ha habido salto
       interval: candidato a intervalo inferior o superior
       k: contador de la condición de salto
'''
def jumpCondition(twoclosevar, varLimdistance, interval, direction, k):
    jump = False
    if twoclosevar:
        k += 1
        if k > 5 and varLimdistance < 1:
            jump = True
            interval = interval + random.uniform(0.7, 0.8) * direction
            k = 0
        if k >= 7 and varLimdistance >= 1+tol:
            k = 0
    else:
        k = 0
    return jump, interval, k


'''
Realiza la minimización del intervalo de búsqueda: primero busca el extremo inferior
que cumple alguna condición satisfactoria y luego el extremo superior, análogamente.
   param:
       varphiLim: Límite de varphi en 0
       varphiLimInf: Límite de varphi en infinito
       galaxdata: diccionario de datos de la galaxia
   return:
       vector solución
           interval: intervalo deducido
           infmin: mínimo valor encontrado en la exploración de la parte izquierda
           supmin: mínimo valor encontrado en la exploración de la parte derecha
           X, Y: datos para la elaboración de gráficas (opcional)
'''
def intervalMin(varphiLim0, varphiLimInf, galaxdata):
    if galaxdata["graphic"]:
        X = []
        Y = []
    infminphi = 10**4  # para devolver el mínimo encontrado en la exploración
    infminx = 0        # para devolver el mínimo encontrado en la exploración

    # INTERVALO INFERIOR #
    maxiter = 0
    direction = -1
    intervalinf = -3
    k = 0
    lastint, rho = phi(np.array([float(10 ** intervalinf)]), galaxdata)
    stop = False
    i = 0.0

    while maxiter < 100 and direction != 0 and k < 10:
        maxiter += 1
        s = 10**(intervalinf + np.array([-0.2, -0.1, 0.0, 0.1, 0.2]))
        varphi, rho = phi(s, galaxdata)
        if min(varphi) < infminphi:
            infminphi = min(varphi)
            pos = (varphi.tolist()).index(infminphi)
            infminx = s[pos]
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
    stop = False
    i = 0.0
    supminphi = 10**4   # para devolver el mínimo encontrado en la exploración
    supminx = 0         # para devolver el mínimo encontrado en la exploración

    while maxiter < 100 and direction != 0 and k < 50:
        maxiter += 1
        s = 10**(intervalsup + np.array([-0.2, -0.1, 0.0, 0.1, 0.2]))
        varphi, rho = phi(s, galaxdata)
        if min(varphi) < supminphi:
            supminphi = min(varphi)
            pos = (varphi.tolist()).index(supminphi)
            supminx = s[pos]
        if galaxdata["graphic"]:
            X.append(10**(intervalsup + np.array([-0.2, -0.1, 0.0, 0.1, 0.2])))
            Y.append(varphi)

        #eval = abs(varphi - varphiLimInf) / varphiLimInf
        if galaxdata["profile"] == 'BUR' or galaxdata["profile"] == 'NFW':      # Mejora propuesta en la memoria
            eval = abs(varphi - varphiLimInf) / (2*varphiLimInf)                # Mejora propuesta en la memoria
        else:                                                                   # Mejora propuesta en la memoria
            eval = abs(varphi - varphiLimInf) / varphiLimInf                    # Mejora propuesta en la memoria

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
    infmin = [infminx, infminphi]
    supmin = [supminx, supminphi]
    if galaxdata["graphic"]:
        sol = [interval, X, Y, infmin, supmin]
    else:
        sol = [interval, infmin, supmin]

    return sol
