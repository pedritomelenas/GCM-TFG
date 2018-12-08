import numpy as np
import scipy.optimize as op

## SOME COMMON FUNCTIONS ##


##  Cálculo del producto escalar pesado definido en (15)
#   param:
#       x: array
#       y: array
#       sigmas: array de pesos
#   return:
#       producto escalar pesado de x e y, con pesos sigmas
def WeighProd(x, y, sigmas):
    m = x.transpose()*sigmas
    n = m.transpose()*y
    return n.sum(axis=0)


##  Cálculo del valor de g cuando s tiende a infinito, definida en la Tabla 2
#   param:
#       x: array
#       model: perfil d e densidad de materia oscura ('ISO', 'BUR', o 'NFW')
#   return:
#       el valor de g cuando s tiende a infinito
def ginf(x, model):
    if model == 'ISO':
        ginf = np.ones(len(x))
    elif model == 'BUR' or model == 'NFW':
        ginf = x**(-1)
    return ginf


##  Cálculo de la ecuación definida en (33)
#   param:
#       t: parámetro de la ecuación
#       ginf: valor de g cuando s tiende a infinito
#       galaxdata: datos de la galaxia
#   return:
#       ecuación definida en (33)
def eqVLimInf(t, ginf, galaxdata):
    return WeighProd(ginf, galaxdata["vones"], galaxdata["weights"]) - \
           WeighProd(galaxdata["vrot"], (ginf / np.sqrt(t * ginf + galaxdata["vbary"] ** 2)), galaxdata["weights"])


##  Cálculo del valor de g cuando s tiende a 0, definida en la Tabla 2
#   param:
#       x: array
#       model: perfil d e denisdad de materia oscura ('ISO', 'BUR', o 'NFW')
#   return:
#       el valor de g cuando s tiende a 0
def g0(x, model):
    if model == 'NFW':
        g0 = x
    elif model == 'BUR' or model == 'ISO':
        g0 = x**2
    return g0


##  Cálculo de la ecuación definida en (35)
#   param:
#       t: parámetro de la ecuación
#       ginf: valor de g cuando s tiende a 0
#       galaxdata: datos de la galaxia
#   return:
#       ecuación definida en (35)
def eqVLim0(t, g0, galaxdata):
    return WeighProd(g0, galaxdata["vones"], galaxdata["weights"]) - \
           WeighProd(galaxdata["vrot"], g0 / np.sqrt(t * g0 + galaxdata["vbary"] ** 2), galaxdata["weights"])


##  Cálculo de W^s(r), ecuación (18)
#   param:
#       r: array de radios
#       s: array de parámetros de inverso de escala
#       model: perfil d e densidad de materia oscura ('ISO', 'BUR', o 'NFW')
#   return:
#       el valor de la ecuación (18) para los parámetros especificados
def v(r, s, model):
    if model == 'ISO':
        v = np.sqrt(4 * np.pi * (np.outer(r, s) - np.arctan(np.outer(r, s))) / np.outer(r, np.ones(len(s))))
    elif model == 'BUR':
        v = np.sqrt(np.pi * (np.log((1+(np.outer(r, s))**2) * ((1 + np.outer(r, s))**2)) - 2*np.arctan(np.outer(r, s)))
                    / np.outer(r, np.ones(len(s))))
    elif model == 'NFW':
        v = np.sqrt(4 * np.pi * (np.log(1 + np.outer(r, s)) / np.outer(r, np.ones(len(s))) - (np.outer(np.ones(len(r)), s)
                                        / (1 + np.outer(r, s)))))
    return v


##  Cálculo de chi^2(rho, s), ecuación (16)
#   param:
#       rho: array de parámetros de densidad central
#       s: array de parámetros de inverso de escala
#       galaxdata: datos de la galaxia
#   return:
#       el valor de la ecuación (16) para los parámetros especificados
def chiquad(rho, s, galaxdata):
    vaux = v(galaxdata["radii"], s, galaxdata["profile"])
    eval = rho * WeighProd(vaux, vaux, galaxdata["weights"]) / (galaxdata["CteDim"] * s ** 3)
    eval -= 2 * (WeighProd(np.dot(np.atleast_2d(galaxdata["vrot"]).T, np.atleast_2d(np.ones(len(s)))),
                           np.sqrt(np.square(vaux) * (np.ones((len(galaxdata["radii"]), 1)) *
                                                      (rho / (galaxdata["CteDim"] * s ** 3))) +
                                   np.dot(np.atleast_2d(np.square(galaxdata["vbary"])).T,
                                          np.atleast_2d(np.ones(len(s))))),
                           galaxdata["weights"]))
    chi = vv(galaxdata) + vvbary(galaxdata) + eval
    return chi


##  Cálculo de rho_0(s), el mínimo de la función chi^2(·, s), Proposición 1
#   param:
#       s: array de parámetros de escala
#       galaxdata: datos de la galaxia
#   return:
#       el valor de rho_0(s) para los parámetros especificados
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


##  Cálculo de g(rho_0,s), ecuación (24)
#   param:
#       s: array de parámetros de escala
#       galaxdata: datos de la galaxia
#   return:
#       el valor de g(rho_0,s) para los parámetros especificados
def alphaMV(s, galaxdata):
    rhoaux = rho(s, galaxdata)
    vaux = v(galaxdata["radii"], s, galaxdata["profile"])
    eval = rhoaux * WeighProd(vaux, vaux, galaxdata["weights"]) / (galaxdata["CteDim"] * s ** 3)
    eval -= 2 * (WeighProd(np.dot(np.atleast_2d(galaxdata["vrot"]).T, np.atleast_2d(np.ones(len(s)))),
                           np.sqrt(np.square(vaux) * (np.ones((len(galaxdata["radii"]), 1)) *
                                                      (rhoaux / (galaxdata["CteDim"] * s ** 3))) +
                                   np.dot(np.atleast_2d(np.square(galaxdata["vbary"])).T, np.atleast_2d(np.ones(len(s))))),
                           galaxdata["weights"]))
    return eval, rhoaux


##  Cálculo del producto escalar pesado de la velocidad rotacional para una determinada galaxia
#   param:
#       galaxdata: datos de la galaxia
#   return:
#       el valor del producto escalar pesado de la velocidad rotacional
def vv(galaxdata):
    v = galaxdata["vrot"]
    return WeighProd(v, v, galaxdata["weights"])


##  Cálculo del producto escalar pesado de la velocidad debida a la materia bariónica para una determinada galaxia
#   param:
#       galaxdata: datos de la galaxia
#   return:
#       el valor del producto escalar pesado de la velocidad debida a la materia bariónica
def vvbary(galaxdata):
    vbary = galaxdata["vbary"]
    return WeighProd(vbary, vbary, galaxdata["weights"])


##  Cálculo de varphi(s) para una determinada galaxia
#   param:
#       s: parámetro de inverso de escala
#       galaxdata: datos de la galaxia
#   return:
#       phi: el valor de varphi para la s especificada
#       rho: el valor de rho_0(s) definido en la Proposición 1
def phi(s, galaxdata):
    alpha, rho = alphaMV(s, galaxdata)
    phi = vv(galaxdata) + vvbary(galaxdata) + alpha
    return phi, rho
