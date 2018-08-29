
import numpy as np


def fit(x, y, deg):
    X = np.zeros([deg+1, deg+1])
    for i in range(0, len(x)):
        for j in range(0, deg+1):
            X[i, j] = x[i]**(deg-j)
    P = np.dot(np.linalg.pinv(X), y)

    return P


def gen(t, p):
    val = np.zeros(len(t))
    deg = len(p)-1
    t_mat = np.zeros(deg+1)
    for i in range(0, len(t)):
        for j in range(0, deg+1):
            t_mat[j] = t[i]**(deg-j)
        val[i] = np.dot(t_mat, p)

    return val


def dev(t, p):
    deg = len(p)-1
    t_mat = np.zeros(deg)
    for i in range(1, deg):
        t_mat[i] = (deg-i+1)*t**(deg-i)
    val = np.dot(t_mat, p[0:-1])

    return val
