import numpy as np
from matplotlib import pyplot as plt 

def g(x, dx, t, l):

    x1 = x + dx
    x2 = x - dx
    f_1 = f(x1, t, l)
    f_2 = f(x2, t, l)

    der_f = (f_1 - f_2)/(2.0 * dx)
    return der_f

def f(x, t, l):

    sum = 0.0
    for n in range(1,l+1):
        sum += np.exp(-n * t**2) * np.cos(n * x)
    sum *= 2.0
    sum += 1.0

    return np.log(sum)

x = np.power(1e-10, 1/3)
print(x)

'''
m = 100
dx=1e-10
t=1e-10
l=10000
x_lst = np.linspace(0.0, 2.0 * np.pi, m)
g_lst = np.zeros(m)
for i in range(len(x_lst)):
    g_lst[i] = t**3 * g(x_lst[i], dx, t, l)

plt.figure()
plt.plot(x_lst, g_lst)
plt.show()
'''




