from vISO import vISO
from rhoISO import rhoISO
from data import radii, CteDim, weights, vbary, vrot
from WeighProd import WeighProd
import numpy as np


# Evaluation of the function g(rho_0(s),s), where g is defined in equation #24
# considering the Isothermal model. By expression #16 this function provides 
# varphi(s) adding appropriate terms. 
#
# s is a row vector containing scale factors 
#
# alphaMVISO(s) returns a column vector whose ith component is g(rho_0(s(i)),s(i))
#
# This function requires global variables defined in the script redMethRotCurveFitting.py
def alphaMVISO(s):
    rhoaux = rhoISO(s)
    #print("rhoaux", rhoaux)
    vaux = vISO(radii, s)
    eval = rhoaux * WeighProd(vaux, vaux, weights) / (CteDim * s ** 3)
    eval -= 2 * (WeighProd(np.dot(np.atleast_2d(vrot).T, np.atleast_2d(np.ones(len(s)))),
                           np.sqrt(np.square(vaux) * (np.ones((len(radii), 1)) *
                                                    (rhoaux / (CteDim * s ** 3))) +
                                   np.dot(np.atleast_2d(np.square(vbary)).T, np.atleast_2d(np.ones(len(s))))), weights))
    return eval
