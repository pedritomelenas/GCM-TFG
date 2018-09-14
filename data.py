import arff
import numpy as np

galaxlist = ["DDO101", "DDO126" , "DDO133"]#, "DDO154", "DDO168"]
galaxies = {}
for i in galaxlist:
    fp = open(i+".arff")
    data = arff.load(fp)
    data = np.array(data['data'])
    galaxies[i] = {
        "R": data[:, 0] * 1000,
        "vrot": abs(data[:, 1]),
        "errs": data[:, 3],
        "vbary": np.sqrt(data[:, 4] ** 2 + data[:, 5] ** 2)
    }
    fp.close()
#fp = open('DDO101.arff')
#data = arff.load(fp)
#data = np.array(data['data'])
#print(data)
#galaxies = {
 #   "DDO101": {
  #      "R": data[:, 0] * 1000,
   #     "vrot": abs(data[:, 1]),
    #    "errs": data[:, 3],
     #   "vbary": np.sqrt(data[:, 4] ** 2 + data[:, 5] ** 2)
    #}
#}
#DDO101 = dict(R=data[:, 0] * 1000, vrot=abs(data[:, 1]), errs=data[:, 3], vbary=np.sqrt(data[:, 4] ** 2 + data[:, 5] ** 2))
#print(DDO101)

#radii = data[:, 0] * 1000
#print(radii)
#print(galaxies["DDO101"]["R"])
#3vrot = abs(data[:, 1])
#p3333333rint(vrot)
#print(galaxies["DDO101"]["vrot"])
#errs = data[:, 3]
#print(errs)
#print(galaxies["DDO101"]["errs"])
#vbary = np.sqrt(data[:, 4] ** 2 + data[:, 5] ** 2)
#print(vbary)
#print(galaxies["DDO101"]["vbary"])

#n = len(radii)
nu = 2
#weights = 1 / ((n - nu) * galaxies["DDO101"]["errs"] ** 2)
CteDim = 10000 / ((3.0856776 ** 2) * 4.51697)
#totalnullvbary = np.sum(vbary) == 0
#somenullvbary = round(np.prod(vbary)) == 0  ###  sustituir por  "0.0 in vbary"
