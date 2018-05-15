from vISO import vISO
from data import radii, weights, CteDim, totalnullvbary, somenullvbary, vbary, vrot
from WeighProd import WeighProd
import numpy as np
import scipy.optimize as op


# Numerical approximation of the function rho_0(s) defined in Proposition 2.1 considering the Isothermal model
#
# s is a row vector containing scale factors 
#
# rhoISO(s) returns a row vector whose ith component is rho_0(s(i))
#
# This function requires global variables defined in the script data.m

def rhoISO(s):
    aux = 0 * s
    vHalos = vISO(radii, s)
    rhs = WeighProd(vHalos, vHalos, weights)  # rhs Eq  # 19 and #20 up to multiplicative constant 1/(s^3 CteDim) .
    #print("RHS = ", rhs)
    rhoVbaryNull = CteDim * (s ** 3) * (WeighProd(np.dot(np.atleast_2d(vrot).T, np.atleast_2d(np.ones(len(s)))),
                                                       vHalos, weights) / rhs) ** 2     # Eq #21

   # print(np.dot(np.atleast_2d(vrot).T, np.atleast_2d(np.ones(len(s)))))       BIEN: (9x1)*(1x5)=(9x5)

    ####################        D U D A         ############################
    # This rootfinding can not be vectorized because of the scalar character of arguments of routine fzero
    # Those cases with vbary null at some radii have to be treated particularly
    ########################################################################

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
        for i in rango[0]:
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
