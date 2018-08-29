
import numpy as np
from scipy.linalg import pascal
from scipy.special import comb
from math import hypot


#  Generate Bernstein Matrix
#  deg: degree of Bezeir curve
def mxbern(t, deg):
    elements = np.size(t)

    if type(t) == np.ndarray:
        n = len(t)
        m = elements/n
        if m != 1:
            raise ValueError('Input t must be a column vector')
        elif min(t) < 0 or max(t) > 1.0:
            raise ValueError('Input nodes t must be within [0, 1]')

        ct = 1.0 - t
        B = np.zeros([n, deg+1])

        for i in range(0, deg+1):
            B[:, i] = (t**i)*(ct**(deg-i))

        if deg < 23:
            lower = (np.cumprod(np.arange(deg, 0.9, -1)) /
                     np.cumprod(np.arange(1, deg+0.1, 1)))
            full = lower.tolist()
            full.insert(0, 1.0)
            B = np.dot(B, np.diag(full))
        else:
            B = np.dot(B, np.diag(np.diag(np.fliplr(pascal(deg+1)))))

    elif type(t) == float:
        ct = 1.0 - t
        B = np.zeros(deg+1)

        for i in range(0, deg+1):
            B[i] = (t**i)*(ct**(deg-i))

        if deg < 23:
            lower = (np.cumprod(np.arange(deg, 0.9, -1)) /
                     np.cumprod(np.arange(1, deg+0.1, 1)))
            full = lower.tolist()
            full.insert(0, 1.0)
            B = np.dot(B, np.diag(full))
        else:
            B = np.dot(B, np.diag(np.diag(np.fliplr(pascal(deg+1)))))

    else:
        raise ValueError('Wrong data type: Only numpy.ndarray or float is \
accepted')

    return B


def aff_angle(X):
    n = len(X)
    meanX = np.ones([n, 2])
    meanX[:, 0] = meanX[:, 0]*np.mean(X, axis=0)[0]
    meanX[:, 1] = meanX[:, 1]*np.mean(X, axis=0)[1]

    Xbar = X - np.dot(np.ones([n, 2]), np.diag(np.mean(X, axis=0)))
    Xcov = np.dot(np.transpose(Xbar), Xbar/n)
    A = np.linalg.inv(Xcov)

    V = X[1:n, :] - X[0:n-1, :]
    t = np.diag(np.dot(np.dot(V, A), np.transpose(V)))**0.5
    V2 = X[2:n, :] - X[0:n-2, :]
    t2 = np.diag(np.dot(np.dot(V2, A), np.transpose(V2)))

    theta = np.zeros(n-1)
    for j in range(1, n-1):
        theta[j] = min(np.pi - np.arccos((t[j-1]**2 + t[j]**2 - t2[j-1]) /
                                         (2*t[j]*t[j-1])), np.pi/2)

    h = np.zeros(n-1)

    h[0] = t[0]*(1.0 + (1.5*theta[1]*t[1])/(t[0] + t[1]))
    for j in range(1, n-2):
        h[j] = t[j]*(1.0 + (1.5*theta[j]*t[j-1])/(t[j-1] + t[j]) +
                     (1.5*theta[j+1]*t[j+1])/(t[j] + t[j+1]))

    h[n-2] = t[n-2]*(1.0 + (1.5*theta[n-2]*t[n-3])/(t[n-3] + t[n-2]))

    h = h.tolist()
    h.insert(0, 0.0)
    h = np.array(h)
    h = np.cumsum(h)
    h = h/h[n-1]

    return h


#  epsilon: the stopping criteria
def fit(data, deg, epsilon):
    #  Step 1

    i = len(data)
    j = deg + 1
    t = aff_angle(data)

    bez_mat = mxbern(t, deg)
    p = np.dot(np.linalg.pinv(bez_mat), data)  # initial guess

    #  Step2
    counter = 0
    resid_old = 0
    resid_new = np.dot(bez_mat, p) - data

    error = np.linalg.norm(resid_new - resid_old)/max(1,
                                                      np.linalg.norm(resid_new))
    while error > epsilon:
        deriv = deg*np.dot(mxbern(t, deg-1), p[1:j, :] - p[0:j-1, :])
        t = t - ((np.dot(deriv[:, 0], resid_new[:, 0]) +
                  np.dot(deriv[:, 1], resid_new[:, 1])) /
                 (deriv[:, 0]**2 + deriv[:, 1]**2))

        t = -min(t)*np.ones(i) + t
        t = t/max(t)

        bez_mat = mxbern(t, deg)
        p = np.dot(np.linalg.pinv(bez_mat), data)
        resid_old = resid_new
        resid_new = np.dot(bez_mat, p) - data

        counter = counter + 1
        error = (np.linalg.norm(resid_new - resid_old) /
                 max(1, np.linalg.norm(resid_new)))

    return [p, t]


def gen(t, deg, p):
    if type(t) == np.ndarray or type(t) == list:
        bernstein = np.zeros(deg+1)
        val = np.zeros([len(t), 2])
        for i in range(0, len(t)):
            for j in range(0, deg+1):
                bernstein[j] = comb(deg, j)*t[i]**j*(1.0-t[i])**(deg-j)
            val[i, :] = np.dot(bernstein, p)
    elif type(t) == float:
        bernstein = np.zeros(deg+1)
        for j in range(0, deg+1):
            bernstein[j] = comb(deg, j)*t**j*(1.0-t)**(deg-j)
        val = np.dot(bernstein, p)

    return val


def get_t(deg, P, R):
    t = 0.3  # Initial guess of t
    error = 1.0
    while error > 10**-4:
        f = gen(t, deg, P)
        error = hypot(f[0] - R[0], f[1] - R[1])
        resid = [f[0]-R[0], f[1]-R[1]]

        if t <= 0.0 or t >= 1.0:
            t = 0.1
        df = deg*np.dot(mxbern(t, deg-1), P[1:deg+1, :] - P[0:deg, :])

        t = t - np.dot(f, resid)/hypot(df[0], df[1])
        t = float(t)

    return t
