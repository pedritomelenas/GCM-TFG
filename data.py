import arff
import numpy as np

galaxlist = ["DDO101", "DDO126", "DDO133", "DDO154", "DDO168"]
galaxies = {}
for i in galaxlist:
    fp = open("galaxies/"+i+".arff")
    data = arff.load(fp)
    data = np.array(data['data'])
    galaxies[i] = {
        "R": data[:, 0] * 1000,
        "vrot": abs(data[:, 1]),
        "errs": data[:, 3],
        "vbary": np.sqrt(data[:, 4] ** 2 + data[:, 5] ** 2)
    }
    fp.close()
nu = 2
CteDim = 10000 / ((3.0856776 ** 2) * 4.51697)
