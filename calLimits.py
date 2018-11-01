import numpy as np
import commonFunctions as cf
import scipy.optimize as op

def calLimits(galaxdata):
    ginf = cf.ginf(galaxdata["radii"], galaxdata["profile"])
    g0 = cf.g0(galaxdata["radii"], galaxdata["profile"])
    # The solution to #34 --> #35 with vbary=0, which is an upper bound for the solutions to #34 --> #35
    XVbaryNullinf = (cf.WeighProd(galaxdata["vrot"], np.sqrt(ginf), galaxdata["weights"]) /
                     cf.WeighProd(ginf, galaxdata["vones"], galaxdata["weights"])) ** 2
    # print("XVbaryNullinf = ", XVbaryNullinf)
    # The solution to #40 is an upper bound for the solutions to #40
    XVbaryNull0 = (cf.WeighProd(galaxdata["vrot"], np.sqrt(g0), galaxdata["weights"]) /
                   cf.WeighProd(g0, galaxdata["vones"], galaxdata["weights"])) ** 2
    # print("XVbaryNull0 = ", XVbaryNull0)

    if galaxdata["totalnullvbary"]:
        # print("if")
        Xinf = XVbaryNullinf
        X0 = XVbaryNull0
    elif galaxdata["somenullvbary"]:
        # print(round(np.prod(vbary)) == 0)
        jinf = -3
        # print(eqVLimInf(10 ** (j) * XVbaryNullinf))
        # print(eqVLim0(10 ** jinf * XVbaryNullinf))
        # if galaxdata["profile"] == 'BUR':
        #    print("ERROR 1: = ", 10 ** jinf * XVbaryNullinf)
        while cf.eqVLimInf(10 ** jinf * XVbaryNullinf, ginf, galaxdata) > 0:
            jinf -= 1
        j0 = -3
        while cf.eqVLim0(10 ** j0 * XVbaryNull0, g0, galaxdata) > 0:
            j0 -= 1
        Xinf = op.brentq(cf.eqVLimInf, 10 ** jinf * XVbaryNullinf, XVbaryNullinf)  # Brent's Method (escalar)
        X0 = op.brentq(cf.eqVLim0, 10 ** j0 * XVbaryNull0, XVbaryNull0)  # Brent's Method (escalar)

    else:
        if cf.eqVLimInf(0, ginf,
                        galaxdata) >= 0:  # <= ? #########################################################################################
            Xinf = 0
        else:
            # Xinf, info, ier, msg = op.fsolve(cf.eqVLimInf, XVbaryNullinf / 2, (ginf, galaxdata), full_output=True)  # Solves equation --> #35
            Xinf = op.brentq(cf.eqVLimInf, 0, XVbaryNullinf, (ginf, galaxdata))
            '''
            if galaxdata["profile"] == 'BUR':
                print("ERROR 2")
                print("ier = ", ier)
                print("msg = ", msg)        ######## fsolve no encuentra solución !!!!!!!
                '''
        if cf.eqVLim0(0, g0, galaxdata) >= 0:
            X0 = 0
        else:
            X0 = op.brentq(cf.eqVLim0, 0, XVbaryNull0, (g0, galaxdata))
            # X0 = op.fsolve(cf.eqVLim0, XVbaryNull0 / 2, (g0, galaxdata))  # Solves equation #38 --> #40

    ## Calculation of the limit value by using Lemma 2.1 and development #16 ##
    varphiLimInf = galaxdata["vv"] + galaxdata["vvbary"] + Xinf * cf.WeighProd(ginf, galaxdata["vones"], galaxdata["weights"])\
                   - 2 * cf.WeighProd(galaxdata["vrot"], np.sqrt(Xinf * ginf + (galaxdata["vbary"] ** 2)), galaxdata["weights"])
    # varphiLimInf = vv + vvbary - 2. * cf.WeighProd(vrot, sqrt(X * ginf(radii) + (vbary. ^ 2)), weights) +
    # X. * cf.WeighProd(ginf(radii), vones, weights)

    ## Calculation of the limit value by using Lemma 2.2 (2.1?) and development #16 ##
    varphiLim0 = galaxdata["vv"] + galaxdata["vvbary"] + X0 * cf.WeighProd(g0, galaxdata["vones"], galaxdata["weights"]) \
                 - 2 * cf.WeighProd(galaxdata["vrot"], np.sqrt(X0 * g0 + (galaxdata["vbary"] ** 2)), galaxdata["weights"])

    limits = [varphiLim0, varphiLimInf]
    return limits
