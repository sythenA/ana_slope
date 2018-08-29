
import matplotlib.pyplot as plt
import numpy as np
import ana_new as an
import pickle
from math import sqrt


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
    plt.figure

    _file = "/home/room801/anaoverland/Steady_fit/ks1_Rise_CurveRec.flow"
    RiseRec = pickle.load(open(_file, 'r'))
    _file = "/home/room801/anaoverland/Steady_fit/ks1_Fall_CurveRec.flow"
    FallRec = pickle.load(open(_file, 'r'))

    ks = [1.01, 1.6, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 40.0,
          50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
    for k in ks:
        keystring = 'k=' + str(k)
        out_q = RiseRec['Q'][keystring]
        out_s = RiseRec['S'][keystring]

        #  Normalized curve to [0, 1]
        out_s[:, 1] = (out_s[:, 1] - min(out_s[:, 1]))/(max(out_s[:, 1]) -
                                                        min(out_s[:, 1]))
        out_q[:, 1] = (out_q[:, 1] - min(out_q[:, 1]))/(max(out_q[:, 1]) -
                                                        min(out_q[:, 1]))

        plt.plot(out_q[:, 1], out_s[0:len(out_q), 1], color='dodgerblue',
                 linestyle='-')

    for k in ks:
        keystring = 'k=' + str(k)
        out_q = FallRec['Q'][keystring]
        out_s = FallRec['S'][keystring]

        #  Normalization to [0, 1]
        out_s[:, 1] = (out_s[:, 1] - min(out_s[:, 1]))/(max(out_s[:, 1]) -
                                                        min(out_s[:, 1]))
        out_q[:, 1] = (out_q[:, 1] - min(out_q[:, 1]))/(max(out_q[:, 1]) -
                                                        min(out_q[:, 1]))

        plt.plot(out_q[:, 1], out_s[0:len(out_q), 1], color='tomato',
                 linestyle='-')

    aa = an.trans_ana(10.0, 10.0*200, 100.0, 201, 0.5, 0.1, 0.01)
    aa.run()
    out_s = aa.out_s
    out_q = aa.out_q

    out_s[:, 1] = (out_s[:, 1] - min(out_s[:, 1]))/(max(out_s[:, 1]) -
                                                    min(out_s[:, 1]))
    out_q[:, 1] = (out_q[:, 1] - min(out_q[:, 1]))/(max(out_q[:, 1]) -
                                                    min(out_q[:, 1]))
    plt.plot(out_q[:, 1]/max(out_q[:, 1]),
             out_s[0:len(out_q), 1]/max(out_s[:, 1]), color='dodgerblue',
             linestyle='-', label='Rising Limb')

    ab = an.trans_ana(10.0*200, 10.0, 100.0, 201, 0.5, 0.1, 0.01)
    ab.run()
    out_s = ab.out_s
    out_q = ab.out_q

    out_s[:, 1] = (out_s[:, 1] - min(out_s[:, 1]))/(max(out_s[:, 1]) -
                                                    min(out_s[:, 1]))
    out_q[:, 1] = (out_q[:, 1] - min(out_q[:, 1]))/(max(out_q[:, 1]) -
                                                    min(out_q[:, 1]))
    plt.plot(out_q[:, 1]/max(out_q[:, 1]),
             out_s[0:len(out_q), 1]/max(out_s[:, 1]), color='tomato',
             linestyle='-', label='Falling Limb')

    #  Rising Limb limit
    ac = an.progression(200.0, 100.0, 201, 0.5, 0.1, 0.01)
    ac.run()
    out_q = ac.out_q
    out_s = ac.out_s

    out_s[:, 1] = (out_s[:, 1] - min(out_s[:, 1]))/(max(out_s[:, 1]) -
                                                    min(out_s[:, 1]))
    out_q[:, 1] = (out_q[:, 1] - min(out_q[:, 1]))/(max(out_q[:, 1]) -
                                                    min(out_q[:, 1]))

    plt.plot(out_q[:, 1], out_s[0:len(out_q), 1], 'k-', lw=1.5)

    #  Falling Limb limit
    ad = an.Recession(200.0, 100.0, 201, 0.5, 0.1, 0.01)
    ad.run()
    out_q = ad.out_q
    out_s = ad.out_s

    out_s[:, 1] = (out_s[:, 1] - min(out_s[:, 1]))/(max(out_s[:, 1]) -
                                                    min(out_s[:, 1]))
    out_q[:, 1] = (out_q[:, 1] - min(out_q[:, 1]))/(max(out_q[:, 1]) -
                                                    min(out_q[:, 1]))
    plt.plot(out_q[:, 1], out_s[0:len(out_q), 1], linestyle='--',
             color='dimgray', lw=1.5)

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.legend(loc=2)
    plt.xlabel('$q^{{}*{}}$')
    plt.ylabel('$s^{{}*{}}$')
    plt.show()
