import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('seaborn-white')

Lmin = 3
Lmax = 4.5
H1 = np.array([1, 0, 0, 0])
H2 = np.array([0, 0, 0, 1])


# 1
def f(x, y):
    A = np.array([[3.0, 0, 0, 0],
                  [0, 3.0, 0, 0],
                  [0, 0, 4.0, 0],
                  [0, 0, 0, 4.5]])
    b = np.array([5, 5, 5, 5])
    c = 6

    z = np.ones((len(x), len(y)), dtype=float)
    for i in range(len(x)):
        for j in range(len(y)):
            stolb = H1 * x[i] + y[j] * H2
            z[i][j] = (stolb.transpose().dot(A.dot(stolb)) - b.dot(stolb) + c)
    return z


x = np.linspace(-5, 5, 500)
y = np.linspace(-5, 5, 500)
Z = f(x, y)
plt.imshow(Z, extent=[-5, 5, -5, 5], origin='lower', cmap='RdGy', alpha=0.5)
# 2
# От 1 до 2 Задаётся задний фон поля
with open(r'C:\Users\user\CLionProjects\SLAE2.0\Data_FPI.txt') as f:
    x1 = []
    y1 = []
    for line in f:
        string = line.split()
        x1.append(float(string[0]))
        y1.append(float(string[1]))
with open(r'C:\Users\user\CLionProjects\SLAE2.0\Data_CG.txt') as f:
    x2 = []
    y2 = []
    for line in f:
        string = line.split()
        x2.append(float(string[0]))
        y2.append(float(string[1]))
with open(r'C:\Users\user\CLionProjects\SLAE2.0\Data_GD.txt') as f:
    x3 = []
    y3 = []
    for line in f:
        string = line.split()
        x3.append(float(string[0]))
        y3.append(float(string[1]))
with open(r'C:\Users\user\CLionProjects\SLAE2.0\Data_FPIo.txt') as f:
    x4 = []
    y4 = []
    for line in f:
        string = line.split()
        x4.append(float(string[0]))
        y4.append(float(string[1]))
with open(r'C:\Users\user\CLionProjects\SLAE2.0\Data_FPIc.txt') as f:
    x5 = []
    y5 = []
    for line in f:
        string = line.split()
        x5.append(float(string[0]))
        y5.append(float(string[1]))


plt.plot(x1, y1, 's-c', alpha=0.7, label="first", lw=0.5, mec='c', mew=2, ms=2)
plt.plot(x2, y2, 'o-m', alpha=0.7, label="first", lw=0.5, mec='m', mew=2, ms=2)
plt.plot(x3, y3, 'X-r', alpha=0.7, label="first", lw=0.5, mec='r', mew=2, ms=2)
plt.plot(x4, y4, '^-g', alpha=0.7, label="first", lw=0.5, mec='g', mew=2, ms=2)
plt.plot(x5, y5, 'o-k', alpha=0.7, label="first", lw=0.5, mec='k', mew=2, ms=2)

plt.show()
