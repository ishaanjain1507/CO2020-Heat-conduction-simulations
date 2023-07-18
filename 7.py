import numpy as np
import matplotlib.pyplot as plt

def f(X, t):
    a = np.array([X[1], X[2], -0.5*X[0]*X[2]])
    return a


def RK2(X0, t0, tn, h):
    Xi = np.zeros((int((tn-t0)/h) + 1, 3))
    ti = np.zeros(int((tn-t0)/h)+1)
    for i in range(0, int((tn-t0)/h)+1):
        X = X0 + h * f(X0, t0 + i*h)
        X0 = X
        Xi[i] = X0
        ti[i] = t0 + i*h
    return Xi, ti

Xi, ti = RK2(np.array([0, 0, 1]), 0, 10, 0.001)
plt.plot(ti, Xi)

plt.show()
print(RK2(np.array([0, 0, 1]), 0, 20, 0.01))
