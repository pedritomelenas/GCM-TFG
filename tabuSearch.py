from collections import deque
import random
import matplotlib.pyplot as plt
import numpy as np

from commonFunctions import phi

###############################################
############### Tabu Search ###################
###############################################

Nv = 5  # mínimo número de soluciones en la vecindad actual
Nmin = 3  # mínimo número de soluciones en la lista tabú
Nmax = 5  # número máximo de soluciones en la lista tabú
K = 15000   # máximo número de iteraciones
e = 10**(-2)    # umbral de precisión
M = 20  # máximo número de soluciones encontradas
N = 200  # máximo número de iteraciones sin que la solución óptima cambie

s0 = []
random.seed(10)
s0.append([np.array(np.array(random.uniform(500, 600)))]) # solución inicial aleatoria
tabulist = deque()  # lista tabú
s_opt = s0   # solución óptima
#print(np.asarray(s0))
phi_opt = phi(np.asarray(s0))
#print(phi_opt)
tabulist.append(500)
m = 0   #  contador de soluciones encontradas
n = 0   # contador de iteraciones sin que sopt cambie
k = 0   # número de iteraciones
vicinity = []
previcinity = []
dic_spm = {}
X = []
Y = []

while m < M and n < N and k < 10000: #cambiar a K
    while len(vicinity) < Nv:
        s_pm = random.choice(np.arange(0.0, 1.0, 0.1))
        if not (s_pm in tabulist) and not (s_pm in previcinity):
            #print("s_pm = ", s_pm)
            previcinity.append(s_pm)
            s_pmf = random.uniform(s_pm, s_pm+0.09)
            dic_spm[s_pmf] = s_pm
            vicinity.append([s_pmf])
    #print("dic = ", dic_spm)
    previcinity.clear()
    #print("vicinity = ", vicinity)
    #print("*******************************")
    if len(vicinity) != 0:
        avicinity = np.asarray(vicinity)
        test =[]
        #print("av = ", avicinity)
        for v in avicinity:
            #print("phi(",v,") = ", phi(v))
            test.append([phi(v)])
        #print("test = ", test)
        #print("min(test) = ", min(test))
        phi_s0 = np.asarray(min(test))
        s0 = avicinity[test.index(phi_s0)]
        #print("s0 = ", s0[0])
        X.append(s0[0])
        Y.append(phi_s0[0][0])
        tabulist.append(dic_spm[s0[0]])
        dic_spm.clear()
        #print("phis0 = ", phi_s0[0][0])
        #print("tabulist = ", tabulist)
        #print(phi_s0, " < ", phi_opt)
        if phi_s0 < phi_opt:
            s_opt = s0
            n = 0
            #print(phi_s0, " - ", phi_opt, " = ", abs(phi_s0 - phi_opt))
            #print(abs(phi_s0 - phi_opt) < e)
            if abs(phi_s0 - phi_opt) < e:
                m += 1
                #print("suma m = ", m)
            phi_opt = phi_s0
            if len(tabulist) > Nmin:
                tabulist.popleft()
            #print(tabulist)
        else:
            n += 1
            #print("suma n = ", n)
            #print("len(tabulist) = ", len(tabulist))
            while len(tabulist) > Nmax:
                tabulist.popleft()
            m = 0
        k += 1
    vicinity.clear()
print("s_opt = ", s_opt," ; phi_opt = ", phi_opt)
print("m = ", m, "; n = ", n, "; k = ", k)
Xa = np.asarray(X)
Ya = np.asarray(Y)
plt.scatter(X, Y)
plt.show()




