import numpy as np

# Isothermal rotation velocity (up to constants, see Eq.#45 --> A1 ?).
#
# x is a column vector containing radii.
# s is a row vector containing scale factors
#
# vISO(x,s) returns a matrix where the (i,j) element
# is the rotation velocity due to Isothemal density profile
# at radious x(i) with scalar factor s(j) and central density equal to 1
# up to multiplicative constants.
# See J. Gunn, J.R. Gott, Astrophys. J. 176 (1972) for more details.


def vISO(x, s):
    #print("arctan(outer(x, s)) = ", np.arctan(np.outer(x, s)))
    #print("outer(x, s) = ", np.outer(x, s))
    # Para s = 10^-12, arctan(outer(x,s)) = outer(x,s)
    return np.sqrt(4 * np.pi * (np.outer(x, s) - np.arctan(np.outer(x, s))) / np.outer(x, np.ones(len(s))))
