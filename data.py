import arff
import numpy as np

fp = open('DDO101.arff')
data = arff.load(fp)
data = np.array(data['data'])

radii = data[:, 0] * 1000
vrot = abs(data[:, 1])

errs = data[:, 3]
vbary = np.sqrt(data[:, 4] ** 2 + data[:, 5] ** 2)

n = len(radii)
nu = 2
weights = 1 / ((n - nu) * errs ** 2)
CteDim = 10000/((3.0856776 ** 2) * 4.51697)
totalnullvbary = np.sum(vbary) == 0
somenullvbary = round(np.prod(vbary)) == 0  ###  sustituir por  "0.0 in vbary"