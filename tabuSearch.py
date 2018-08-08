from collections import deque
import random

import numpy as np

from redMethRotCurveFitting import phi

###############################################
############### Tabu Search ###################
###############################################

Nv = 5  # mínimo número de soluciones en la vecindad actual
Nmin = 2  # mínimo número de soluciones en la lista tabú
Nmax = 40  # número máximo de soluciones en la lista tabú
K = 15000   # máximo número de iteraciones
e = 1    # umbral de precisión
M = 200  # máximo número de soluciones encontradas
N = 1000  # máximo número de iteraciones sin que la solución óptima cambie


s0 = []
random.seed(4)
s0.append([np.array(np.array(random.uniform(0, 10**3)))]) # solución inicial aleatoria
tabulist = deque()  # lista tabú
s_opt = s0   # solución óptima
print(s0)
#print(np.asarray(s0))
phi_opt = phi(np.asarray(s0))
#print(phi_opt)
tabulist.append(s_opt)
m = 0   #  contador de soluciones encontradas
n = 0   # contador de iteraciones sin que sopt cambie
k = 0   # número de iteraciones
vicinity = []

while m < M and n < N and k < 10: #cambiar a K
    while len(vicinity) < Nv:
        s_pm = random.choice(range(0, 100, 10))
        list_spm = np.linspace(s_pm, s_pm+10)
        print(list_spm)
        if not (s_pm in tabulist):
            vicinity.append([s_pm])
    #print(vicinity)
    #print("*******************************")
    if len(vicinity) != 0:
        avicinity = np.asarray(vicinity)
        test =[]
        print("av = ", avicinity)
        '''
        for v in avicinity:
            #print("phi(",v,") = ", phi(v))
            test.append([phi(v)])
        #print("test = ", test)
        #print("min(test) = ", min(test))
        phi_s0 = np.asarray(min(test))
        s0 = avicinity[test.index(phi_s0)]
        #print(s0)
        #print(tabulist)
        tabulist.extend(range(int(s0 - 10), int(s0 + 10)))
        #print(tabulist)
        #print(phi_s0, " < ", phi_opt)
        if phi_s0 < phi_opt:
            print("s0 = ", s0)
            s_opt = s0
            n = 0
            #print(phi_s0, " - ", phi_opt)
            if abs(phi_s0 - phi_opt) < e:
                m += 1
                print("suma m = ", m)
            phi_opt = phi_s0
            if len(tabulist) > Nmin:
                tabulist.popleft()
            #print(tabulist)
        else:
            n += 1
            while len(tabulist) > Nmax:
                tabulist.popleft()
            m = 0'''
        k += 1
    vicinity.clear()
#print("s_opt = ", s_opt," ; phi_opt = ", phi_opt)
#print("m = ", m, "; n = ", n, "; k = ", k)'''




