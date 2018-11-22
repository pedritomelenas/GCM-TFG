import numpy as np
import commonFunctions as cf
import scipy.optimize as op

## CALCULATION OF VARPHI LIMITS ##

def calLimits(galaxdata):
    ginf = cf.ginf(galaxdata["radii"], galaxdata["profile"])
    g0 = cf.g0(galaxdata["radii"], galaxdata["profile"])
    # The solution to #34 --> #35 with vbary=0, which is an upper bound for the solutions to #34 --> #35
    XVbaryNullinf = (cf.WeighProd(galaxdata["vrot"], np.sqrt(ginf), galaxdata["weights"]) /
                     cf.WeighProd(ginf, galaxdata["vones"], galaxdata["weights"])) ** 2
    # The solution to #40 is an upper bound for the solutions to #40
    XVbaryNull0 = (cf.WeighProd(galaxdata["vrot"], np.sqrt(g0), galaxdata["weights"]) /
                   cf.WeighProd(g0, galaxdata["vones"], galaxdata["weights"])) ** 2

    if galaxdata["totalnullvbary"]:
        Xinf = XVbaryNullinf
        X0 = XVbaryNull0
    elif galaxdata["somenullvbary"]:
        jinf = -3
        while cf.eqVLimInf(10 ** jinf * XVbaryNullinf, ginf, galaxdata) > 0:
            jinf -= 1
        j0 = -3
        while cf.eqVLim0(10 ** j0 * XVbaryNull0, g0, galaxdata) > 0:
            j0 -= 1
        Xinf = op.brentq(cf.eqVLimInf, 10 ** jinf * XVbaryNullinf, XVbaryNullinf)  # Brent's Method (escalar)
        X0 = op.brentq(cf.eqVLim0, 10 ** j0 * XVbaryNull0, XVbaryNull0)  # Brent's Method (escalar)

    else:
        if cf.eqVLimInf(0, ginf, galaxdata) >= 0:
            Xinf = 0
        else:
            Xinf = op.brentq(cf.eqVLimInf, 0, XVbaryNullinf, (ginf, galaxdata))  # Solves equation --> #35
        if cf.eqVLim0(0, g0, galaxdata) >= 0:
            X0 = 0
        else:
            X0 = op.brentq(cf.eqVLim0, 0, XVbaryNull0, (g0, galaxdata))  # Solves equation #38 --> #40

    ## Calculation of the limit value by using Lemma and development #16 ##
    varphiLimInf = galaxdata["vv"] + galaxdata["vvbary"] + Xinf * cf.WeighProd(ginf, galaxdata["vones"], galaxdata["weights"])\
                   - 2 * cf.WeighProd(galaxdata["vrot"], np.sqrt(Xinf * ginf + (galaxdata["vbary"] ** 2)), galaxdata["weights"])

    ## Calculation of the limit value by using Lemma and development #16 ##
    varphiLim0 = galaxdata["vv"] + galaxdata["vvbary"] + X0 * cf.WeighProd(g0, galaxdata["vones"], galaxdata["weights"]) \
                 - 2 * cf.WeighProd(galaxdata["vrot"], np.sqrt(X0 * g0 + (galaxdata["vbary"] ** 2)), galaxdata["weights"])

    limits = [varphiLim0, varphiLimInf]
    return limits
