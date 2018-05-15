from WeighProd import WeighProd
import numpy as np

v = np.array([0.1,0.1])
m = np.array([[1,3],[2,4]])

x1 = np.arange(3.0)
x2 = np.arange(3.0)

a = np.array([5,4])[np.newaxis]


x = np.ones(5)
print(x)
print(x.T)
y = np.atleast_2d(x)
print(y)
print(y.T)



