import numpy as np
import scipy.optimize as op

def WeighProd(x, y, sigmas):
    m = x.transpose()*sigmas
    n = m.transpose()*y
    return n.sum(axis=0)

def ginf(x, model):
    if model == 'ISO':
        ginf = np.ones(len(x))
    elif model == 'BUR' or model == 'NFW':
        ginf = x**(-1)
    return ginf

def eqVLimInf(t, ginf, galaxdata):
    return WeighProd(ginf, galaxdata["vones"], galaxdata["weights"]) - \
           WeighProd(galaxdata["vrot"], (ginf / np.sqrt(t * ginf + galaxdata["vbary"] ** 2)), galaxdata["weights"])

def g0(x, model):
    if model == 'NFW':
        g0 = x
    elif model == 'BUR' or model == 'ISO':
        g0 = x**2
    return g0


def eqVLim0(t, g0, galaxdata):
    return WeighProd(g0, galaxdata["vones"], galaxdata["weights"]) - \
           WeighProd(galaxdata["vrot"], g0 / np.sqrt(t * g0 + galaxdata["vbary"] ** 2), galaxdata["weights"])

def v(x, s, model):
    if model == 'ISO':
        v = np.sqrt(4 * np.pi * (np.outer(x, s) - np.arctan(np.outer(x, s))) / np.outer(x, np.ones(len(s))))
    elif model == 'BUR':
        v = np.sqrt(np.pi * (np.log((1+(np.outer(x, s))**2) * ((1 + np.outer(x, s))**2)) - 2*np.arctan(np.outer(x, s)))
                    / np.outer(x, np.ones(len(s))))
    elif model == 'NFW':
        v = np.sqrt(4 * np.pi * (np.log(1 + np.outer(x, s)) / np.outer(x, np.ones(len(s))) - (np.outer(np.ones(len(x)), s)
                                        / (1 + np.outer(x, s)))))
    return v

def rho(s, galaxdata):
    aux = 0 * s
    vHalos = v(galaxdata["radii"], s, galaxdata["profile"])     # Para s = 10^-12, vHalos = 0
    rhs = WeighProd(vHalos, vHalos, galaxdata["weights"])  # rhs Eq  # 19 and #20 up to multiplicative constant 1/(s^3 CteDim)
    rhoVbaryNull = galaxdata["CteDim"] * (s ** 3) * (WeighProd(np.dot(np.atleast_2d(galaxdata["vrot"]).T,
                                                                      np.atleast_2d(np.ones(len(s)))),
                                                               vHalos, galaxdata["weights"]) / rhs) ** 2  # Eq #21
    if galaxdata["totalnullvbary"]:
        aux = rhoVbaryNull
    elif galaxdata["somenullvbary"]:
        rango = np.nonzero(rhs)
        for i in rango:
            def rhoequation(t):
                return rhs[i] - WeighProd(galaxdata["vrot"], (vHalos[:][i] ** 2) / np.sqrt(t *
                                            (vHalos[:][i] ** 2) / ((s[i] ** 3) * galaxdata["CteDim"]) +
                                                                                   (galaxdata["vbary"] ** 2)), galaxdata["weights"])
            j = -3
            while rhoequation((10 ** j) * rhoVbaryNull[i]) > 0:
                j -= 1
            aux[i] = op.brentq(rhoequation, (10 ** j) * rhoVbaryNull[i], rhoVbaryNull[i])

    else:
        lhs = WeighProd(np.dot(np.atleast_2d(galaxdata["vrot"]).T, np.atleast_2d(np.ones(len(s)))), (vHalos ** 2) /
                        (np.dot(np.atleast_2d(galaxdata["vbary"]).T, np.atleast_2d(np.ones(len(s))))), galaxdata["weights"])
        rango = np.where(rhs < lhs)
        for i in rango[0]:
            def rhoequation(t):
                return rhs[i] - WeighProd(galaxdata["vrot"], (np.square(vHalos[:, i])) /
                                          np.sqrt(t * np.square((vHalos[:, i])) / ((s[i] ** 3) * galaxdata["CteDim"]) +
                                                                np.square(galaxdata["vbary"])), galaxdata["weights"])
            aux[i] = op.brentq(rhoequation, 0, rhoVbaryNull[i])
    return aux

def alphaMV(s, galaxdata):
    rhoaux = rho(s, galaxdata)
    vaux = v(galaxdata["radii"], s, galaxdata["profile"])
    eval = rhoaux * WeighProd(vaux, vaux, galaxdata["weights"]) / (galaxdata["CteDim"] * s ** 3)
    eval -= 2 * (WeighProd(np.dot(np.atleast_2d(galaxdata["vrot"]).T, np.atleast_2d(np.ones(len(s)))),
                           np.sqrt(np.square(vaux) * (np.ones((len(galaxdata["radii"]), 1)) *
                                                      (rhoaux / (galaxdata["CteDim"] * s ** 3))) +
                                   np.dot(np.atleast_2d(np.square(galaxdata["vbary"])).T, np.atleast_2d(np.ones(len(s))))),
                           galaxdata["weights"]))
    return eval

def vv(galaxdata):
    v = galaxdata["vrot"]
    return WeighProd(v, v, galaxdata["weights"])

def vvbary(galaxdata):
    vbary = galaxdata["vbary"]
    return WeighProd(vbary, vbary, galaxdata["weights"])

def phi(s, galaxdata):
    phi = vv(galaxdata) + vvbary(galaxdata) + alphaMV(s, galaxdata)
    return phi
