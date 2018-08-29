
import pickle
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
import bezeir_fit as bf


def writeTxt(file_path, data):

    try:
        f = open(file_path, 'w')
        for i in range(0, len(data)):
            line = ""
            if type(data[i]) == str:
                line = data[i] + "\n"
            else:
                for j in range(0, len(data[i])):
                    line = line + str(data[i][j]) + " "
                line = line + "\n"

            f.write(line)
    except(IOError):
        print ("Can't create text file %s, please check folder or \
writing authorizations." % file_path)


def reverse_t(t):
    new_t = np.zeros(len(t))
    for i in range(0, len(t)):
        new_t[i] = 1.0 - t[i]
    return new_t


def mat_flip(mat):
    rev_mat = np.fliplr([mat])[0]

    return rev_mat


def max_outflow(ie):
    ie = float(ie)/3600/1000
    return 100.0*ie


def max_sto(ie):
    ie = float(ie)/3600/1000
    alpha = 1.0/0.1*sqrt(0.01)
    length = 100.0
    m = 5.0/3
    sto = m/(m+1)*(ie*length**(m+1)/alpha)**(1/m)

    return sto


def run():
    _file = "/home/room801/anaoverland/Steady_fit/ks1_Rise_CurveRec.flow"
    RiseRec = pickle.load(open(_file, 'r'))
    _file = "/home/room801/anaoverland/Steady_fit/ks1_Fall_CurveRec.flow"
    FallRec = pickle.load(open(_file, 'r'))

    RiseParam = pickle.load(
        open("/home/room801/anaoverland/Steady_fit/fit_1_RiseLimb.pick", 'r'))
    FallParam = pickle.load(
        open("/home/room801/anaoverland/Steady_fit/fit_1_FallLimb.pick", 'r'))
    plt.figure

    ks = [1.01, 1.05, 1.1, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0,
          7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0,
          90.0, 100.0]
    #  ks = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0]
    k_and_NS = list()
    for k in ks:
        keyStr = 'k=' + str(float(k))
        out_s = RiseRec['S'][keyStr]
        out_q = RiseRec['Q'][keyStr]

        #  Normalized curve to [0, 1]
        out_s[:, 1] = (out_s[:, 1] - min(out_s[:, 1]))/(max(out_s[:, 1]) -
                                                        min(out_s[:, 1]))
        out_q[:, 1] = (out_q[:, 1] - min(out_q[:, 1]))/(max(out_q[:, 1]) -
                                                        min(out_q[:, 1]))

        #  Interpolation on flowrate record
        standard_s = np.arange(0, 1.001, 0.01)
        standard_q = np.interp(standard_s, out_s[0:len(out_q), 1], out_q[:, 1])
        rec_sq = np.zeros([len(standard_s), 2])
        rec_sq[:, 0] = standard_s
        rec_sq[:, 1] = standard_q

        t = bf.aff_angle(rec_sq)
        P = RiseParam.ref_points(k)
        sq = bf.gen(t, 6, P)

        R2 = 0.0
        for j in range(0, len(sq)):
            R2 = R2 + ((sq[j, 0] - rec_sq[j, 0])**2 +
                       (sq[j, 1] - rec_sq[j, 1])**2)
        R2 = sqrt(R2)

        plt.plot(out_q[:, 1], out_s[0:len(out_q), 1], color='aqua',
                 linewidth='1.2')
        plt.plot(sq[:, 1], sq[:, 0], '-b')
        print k, R2
        k_and_NS.append([k, R2])

    for k in ks:
        keyStr = 'k=' + str(float(k))
        out_s = FallRec['S'][keyStr]
        out_q = FallRec['Q'][keyStr]

        #  Normalization to [0, 1]
        out_s[:, 1] = (out_s[:, 1] - min(out_s[:, 1]))/(max(out_s[:, 1]) -
                                                        min(out_s[:, 1]))
        out_q[:, 1] = (out_q[:, 1] - min(out_q[:, 1]))/(max(out_q[:, 1]) -
                                                        min(out_q[:, 1]))

        standard_s = np.arange(0, 1.001, 0.01)
        standard_q = np.interp(standard_s, mat_flip(out_s[0:len(out_q), 1]),
                               mat_flip(out_q[:, 1]))

        rec_sq = np.zeros([len(standard_s), 2])
        rec_sq[:, 0] = mat_flip(standard_s)
        rec_sq[:, 1] = mat_flip(standard_q)

        t = bf.aff_angle(rec_sq)
        P = FallParam.ref_points(1./k)
        P[0] = [1.0, 1.0]
        P[-1] = [0.0, 0.0]
        sq = bf.gen(t, 6, P)

        R2 = 0.0
        for j in range(0, len(sq)):
            R2 = R2 + ((sq[j, 0] - rec_sq[j, 0])**2 +
                       (sq[j, 1] - rec_sq[j, 1])**2)
        R2 = sqrt(R2)

        plt.plot(out_q[:, 1], out_s[0:len(out_q), 1], color='darksalmon',
                 lw=1.2)
        plt.plot(sq[:, 1], sq[:, 0], '-r')

        k_and_NS.append([1./k, R2])

    writeTxt('k_and_NS.txt', k_and_NS)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('$q^{{}*{}}$')
    plt.ylabel('$s^{{}*{}}$')
    plt.show()
