import numpy as np
from data import radii, weights, CteDim, totalnullvbary, somenullvbary, vbary, vrot
import scipy.optimize as op

def WeighProd(x, y, sigmas):
    m = x.transpose()*sigmas    #np.multiply(x, sigmas)
    n = m.transpose()*y     #np.multiply(m, y)
    return n.sum(axis=0)

def ginf(x, model):
    if model == 'ISO':
        ginf = np.ones(len(x))
    elif model == 'BUR' or model == 'NFW':
        ginf = x**(-1)
    return ginf

def g0(x, model):
    if model == 'NFW':
        g0 = x
    elif model == 'BUR' or model == 'ISO':
        g0 = x**2
    return g0

def v(x, s, model):
    if model == 'ISO':
        v = np.sqrt(4 * np.pi * (np.outer(x, s) - np.arctan(np.outer(x, s))) / np.outer(x, np.ones(len(s))))
        # aux=sqrt(4*pi*(x*s-atan(x*s))./(x*ones(size(s))));
    elif model == 'BUR':
        v = np.sqrt(np.pi * (np.log((1+(np.outer(x, s))**2) * ((1 + np.outer(x, s))**2)) - 2*np.arctan(np.outer(x, s)))
                    / np.outer(x, np.ones(len(s))))
        # aux=sqrt(pi*(log((1+(x*s).^2).*((1+x*s).^2))-2*atan(x*s))./(x*ones(size(s))));
    elif model == 'NFW':
        v = np.sqrt(4 * np.pi * (np.log(1 + np.outer(x, s) / np.outer(x, np.ones(len(s))) - np.outer(np.ones(len(s)), s)
                                        / (1 + np.outer(x, s)))))
        # aux = sqrt(4 * pi * (log(1 + x * s). / (x * ones(size(s))) - (ones(size(x)) * s). / (1 + x * s)));
    return v

def rho(s, model):
    aux = 0 * s
    if model == 'ISO':
        vHalos = v(radii, s, 'ISO')
    elif model == 'BUR':
        vHalos = v(radii, s, 'BUR')
    elif model == 'NFW':
        vHalos = v(radii, s, 'NFW')
    # print("vHalos = ", vHalos) # Para s = 10^-12, vHalos = 0
    rhs = WeighProd(vHalos, vHalos, weights)  # rhs Eq  # 19 and #20 up to multiplicative constant 1/(s^3 CteDim) .
    # print("RHS = ", rhs)
    rhoVbaryNull = CteDim * (s ** 3) * (WeighProd(np.dot(np.atleast_2d(vrot).T, np.atleast_2d(np.ones(len(s)))),
                                                  vHalos, weights) / rhs) ** 2  # Eq #21
    if totalnullvbary:
        aux = rhoVbaryNull
    elif somenullvbary:
        rango = np.nonzero(rhs)
        for i in rango:
            def rhoequation(t):
                return rhs[i] - WeighProd(vrot, (vHalos[:][i] ** 2) / np.sqrt(t *
                                            (vHalos[:][i] ** 2) / ((s[i] ** 3) * CteDim) +
                                                                                   (vbary ** 2)), weights)
            j = -3
            while rhoequation((10 ** j) * rhoVbaryNull[i]) > 0:
                j -= 1
            aux[i] = op.fsolve(rhoequation, ((10 ** j) * rhoVbaryNull[i] + rhoVbaryNull[i]) / 2) ## MIRAR MULTIDIMENSIONAL

    else:
        lhs = WeighProd(np.dot(np.atleast_2d(vrot).T, np.atleast_2d(np.ones(len(s)))), (vHalos ** 2) /
                        (np.dot(np.atleast_2d(vbary).T, np.atleast_2d(np.ones(len(s))))), weights)
        #print("LHS = ", lhs)
        rango = np.where(rhs < lhs)
        #print("RANGO = ", rango)
        for i in rango[0]:      ## Intentar VECTORIZAR
            #print(i)
            #print(vHalos)
            #print(vHalos[:, i])
            #print(rhs[i])
            #print(s[i])
            def rhoequation(t):
                return rhs[i] - WeighProd(vrot, (np.square(vHalos[:, i])) / np.sqrt(t * np.square((vHalos[:, i])) / ((s[i] ** 3) * CteDim) +
                                                                                   np.square(vbary)), weights)
            aux[i] = op.fsolve(rhoequation, rhoVbaryNull[i] / 2)
    return aux

def alphaMV(s, model):
    if model == 'ISO':
        rhoaux = rho(s, 'ISO')
        vaux = v(radii, s, 'ISO')
    elif model == 'BUR':
        rhoaux = rho(s, 'BUR')
        vaux = v(radii, s, 'BUR')
    elif model == 'NFW':
        rhoaux = rho(s, 'NFW')
        vaux = v(radii, s, 'NFW')

    eval = rhoaux * WeighProd(vaux, vaux, weights) / (CteDim * s ** 3)
    eval -= 2 * (WeighProd(np.dot(np.atleast_2d(vrot).T, np.atleast_2d(np.ones(len(s)))),
                           np.sqrt(np.square(vaux) * (np.ones((len(radii), 1)) *
                                                      (rhoaux / (CteDim * s ** 3))) +
                                   np.dot(np.atleast_2d(np.square(vbary)).T, np.atleast_2d(np.ones(len(s))))), weights))
    return eval

def phi(s, vv, vvbary, model):
    if model == 'ISO':
        phi = vv + vvbary + alphaMV(s, 'ISO')
    elif model == 'BUR':
        phi = vv + vvbary + alphaMV(s, 'BUR')
    elif model == 'NFW':
        phi = vv + vvbary + alphaMV(s, 'NFW')
    return phi