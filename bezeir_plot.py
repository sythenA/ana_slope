
import bezeir_fit as bf
import matplotlib.pyplot as plt
import numpy as np


def run():
    t = np.array([0.0, 0.2, 1.0])
    B = bf.mxbern(t, 2)
    P = [[0.0, 5.0], [5.0, 8.0], [15.0, 5.0]]

    P = np.array(P)
    b = np.dot(B, P)

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(6, 6))
    axes[0, 0].plot(P[:, 0], P[:, 1], 'r')
    axes[0, 0].plot(b[:, 0], b[:, 1], 'b')
    axes[0, 0].set_title('$t_1 = 0.2$')

    t = np.array([0.0, 0.4, 1.0])
    B = bf.mxbern(t, 2)
    b = np.dot(B, P)
    axes[0, 1].plot(P[:, 0], P[:, 1], 'r')
    axes[0, 1].plot(b[:, 0], b[:, 1])
    axes[0, 1].set_title('$t_1 = 0.4$')

    t = np.array([0.0, 0.5, 1.0])
    B = bf.mxbern(t, 2)
    b = np.dot(B, P)
    axes[0, 2].plot(P[:, 0], P[:, 1], 'r')
    axes[0, 2].plot(b[:, 0], b[:, 1])
    axes[0, 2].set_title('$t_1 = 0.5$')

    t = np.array([0.0, 0.6, 1.0])
    B = bf.mxbern(t, 2)
    b = np.dot(B, P)
    axes[1, 0].plot(P[:, 0], P[:, 1], 'r')
    axes[1, 0].plot(b[:, 0], b[:, 1])
    axes[1, 0].set_title('$t_1 = 0.6$')

    t = np.array([0.0, 0.8, 1.0])
    B = bf.mxbern(t, 2)
    b = np.dot(B, P)
    axes[1, 1].plot(P[:, 0], P[:, 1], 'r')
    axes[1, 1].plot(b[:, 0], b[:, 1])
    axes[1, 1].set_title('$t_1 = 0.8$')

    t = np.array([0.0, 0.9, 1.0])
    B = bf.mxbern(t, 2)
    b = np.dot(B, P)
    axes[1, 2].plot(P[:, 0], P[:, 1], 'r')
    axes[1, 2].plot(b[:, 0], b[:, 1])
    axes[1, 2].set_title('$t_1 = 0.9$')

    plt.show()
