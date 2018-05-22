import arff
import numpy as np
import scipy.optimize as op
from data import radii, vrot, errs, vbary, n, nu, weights, CteDim, totalnullvbary, somenullvbary
from WeighProd import WeighProd
from alphaMVISO import alphaMVISO


#####################________ REVISAR PRECISIÓN __________######################
# print(1/((n-nu)*errs[0]**2))
# print(1/((n-nu)*errs[1]**2))
# print(1/((n-nu)*errs[2]**2))
# print(weights)
#####################_____________________________________#######################


# print(WeighProd(np.array([1,2,3]), np.array([1,2,3]), np.array([1,1,1])))

vv = WeighProd(vrot, vrot, weights)  # vrot * np.diag(weights) * vrot
vvbary = WeighProd(vbary, vbary, weights)  # vbary * np.diag(weights) * vbary
vones = np.ones(n)

if n - nu > 0:

    #########################################################
    ##### Computation of the limit values for varphi(s) #####
    #########################################################

    def ginf(x):
        return np.ones(len(x))


    # Function defining equation #34 --> #35
    def eqVLimInf(t):
        return WeighProd(ginf(radii), vones, weights) - \
               WeighProd(vrot, (ginf(radii) / np.sqrt(t * ginf(radii) + vbary ** 2)), weights)
        # (ginf(radii) * np.diag(weights) * vones) / \
        # (vrot * np.diag(weights) * (ginf(radii) / np.sqrt(t * ginf(radii) + vbary ** 2)))
        # WeighProd(ginf(radii),vones,weights) - WeighProd(vrot,ginf(radii)./sqrt(t*ginf(radii)+(vbary.^2)),weights)


    def g0(x):
        return x ** 2


    # Function defining equation #38 --> #40
    def eqVLim0(t):
        # WeighProd(g0(radii),vones,weights) - WeighProd(vrot,g0(radii)./sqrt(t.*g0(radii)+(vbary.^2)),weights)
        return WeighProd(g0(radii), vones, weights) - WeighProd(vrot, g0(radii) / np.sqrt(t * g0(radii) + vbary ** 2),
                                                                weights)

    # print(WeighProd(vrot, np.sqrt(ginf(radii)), weights))
    # print(WeighProd(ginf(radii), vones, weights))
    # print(WeighProd(vrot, np.sqrt(ginf(radii)), weights)/WeighProd(ginf(radii), vones, weights))
    # print(ginf(radii) * np.diag(weights) * vones)
    # print(np.shape(ginf(radii) * np.diag(weights) * vones))

    # The solution to #34 --> #35 with vbary=0, which is an upper bound for the solutions to #34 --> #35
    XVbaryNullinf = (WeighProd(vrot, np.sqrt(ginf(radii)), weights) / WeighProd(ginf(radii), vones, weights)) ** 2

    # The solution to #40 is an upper bound for the solutions to #40
    XVbaryNull0 = (WeighProd(vrot, np.sqrt(g0(radii)), weights) / WeighProd(g0(radii), vones, weights)) ** 2

    # print(XVbaryNullinf)
    # print(XVbaryNull0)

    if totalnullvbary:
        # print("if")
        Xinf = XVbaryNullinf
        X0 = XVbaryNull0
    elif somenullvbary:
        # print(round(np.prod(vbary)) == 0)
        jinf = -3
        # print(eqVLimInf(10 ** (j) * XVbaryNullinf))
        # print(eqVLim0(10 ** jinf * XVbaryNullinf))
        while eqVLimInf(10 ** jinf * XVbaryNullinf) > 0:
            jinf -= 1
        j0 = -3
        while eqVLim0(10 ** j0 * XVbaryNull0) > 0:
            j0 -= 1
        Xinf = op.brentq(eqVLimInf, 10 ** jinf * XVbaryNullinf, XVbaryNullinf)  # Brent's Method (escalar)
        # Xinf = op.fsolve(eqVLimInf, (10 ** jinf * XVbaryNullinf + XVbaryNullinf) / 2)
        X0 = op.brentq(eqVLim0, 10 ** j0 * XVbaryNull0, XVbaryNull0)  # Brent's Method (escalar)
        # X0 = op.fsolve(eqVLim0, (10 ** j0 * XVbaryNull0 + XVbaryNull0) / 2)  ## USAR ESCALAR
        # X = fzero(equationVLimInf, [10 ^ (j) * XVbaryNull, XVbaryNull]); # Solves equation  # 34 --> #35

    else:
        if eqVLimInf(0) >= 0:
            Xinf = 0
        else:
            Xinf = op.fsolve(eqVLimInf, XVbaryNullinf / 2)  # Solves equation #38 --> #40
        if eqVLim0(0) >= 0:
            X0 = 0
        else:
            X0 = op.fsolve(eqVLim0, XVbaryNull0 / 2)  # Solves equation #38 --> #40

    print("Xinf = ", Xinf)
    print("X0 = ", X0)

    ## Calculation of the limit value by using Lemma 2.1 and development #16 ##
    varphiLimInf = vv + vvbary - 2 * WeighProd(vrot, np.sqrt(Xinf * ginf(radii) + (vbary ** 2)), weights) + \
                   Xinf * WeighProd(ginf(radii), vones, weights)
    # varphiLimInf = vv + vvbary - 2. * WeighProd(vrot, sqrt(X * ginf(radii) + (vbary. ^ 2)), weights) +
    # X. * WeighProd(ginf(radii), vones, weights)

    ## Calculation of the limit value by using Lemma 2.2 (2.1?) and development #16 ##
    varphiLim0 = vv + vvbary + X0 * WeighProd(g0(radii), vones, weights) \
                 - 2 * WeighProd(vrot, np.sqrt(X0 * g0(radii) + (vbary ** 2)), weights)

    print("varphiLimInf = ", varphiLimInf)
    print("varphiLim0 = ", varphiLim0)

    ###############################################
    ##### Minimization interval determination #####
    ###############################################
    #print(radii)
    #print(10**(-3+np.array([-3.2000,-3.1000,-3.0000,-2.9000,-2.8000])))

    def phi(s):
        return vv+vvbary+alphaMVISO(s)

    print("phi(",(Xinf+X0)/2,")=",phi(np.array([(Xinf+X0)/2])))

    print("phi(Xinf)=", phi(Xinf))
    print("phi(1000)=", phi(np.array([1000.0])))
    print("phi(X0)=",phi(X0))
    print("phi(0.00000000001)=",phi(np.array([0.00000000001]))) # 10 ceros
    print("phi(0.0000000001)=", phi(np.array([0.0000000001]))) # 9 ceros
    #print("****", alphaMVISO(10**(-3+np.array([-0.2000,-0.1000,-0.0000,0.1000,0.2000]))))
    #print("----", phi(10**(-3+np.array([-0.2000,-0.1000,-0.0000,0.1000,0.2000]))))
