#import arff                            CON PYTHON 3.5
from scipy.io import arff
import numpy as np

## GALAXIES DATA READING AND CONSTANTS ##

galaxlist = ["DDO43", "DDO46", "DDO47", "DDO52", "DDO53", "DDO70", "DDO87", "DDO101", "DDO126",
             "DDO133", "DDO154", "DDO168", "DDO210", "DDO216", "F564_v3", "haro29", "haro36",
             "ic10", "ic1613", "NGC1569", "NGC2366", "NGC3738", "UGC8508"]
profiles = ['ISO', 'BUR', 'NFW']
galaxies = {}
for i in galaxlist:
    fp = open("galaxies/"+i+".arff")
    # data = arff.load(fp)              # CON PYTHON 3.5
    # data = np.array(data['data'])     # CON PYTHON 3.5
    dt, metadt = arff.loadarff(fp)      # CON PYTHON 3.7
    data = []                           # CON PYTHON 3.7
    for d in dt.tolist():               # CON PYTHON 3.7
        data.append(np.asarray(d))
    data = np.asarray(data)             # CON PYTHON 3.7
    galaxies[i] = {
        "R": data[:, 0] * 1000,
        "vrot": abs(data[:, 1]),
        "errs": data[:, 3],
        "vbary": np.sqrt(data[:, 4] ** 2 + data[:, 5] ** 2)
    }
    fp.close()
nu = 2  # SÓLO PARA ISO, BUR Y NFW (PARA EIN nu = 3)
CteDim = 10000 / ((3.0856776 ** 2) * 4.51697)
# 1 kpc = 3.0856776 * 10^{16} km; 1pc = 0.001 kpc = 3.0856776 * 10^{13} km
# M_{solar} = 1.98847 * 10^{30} kg