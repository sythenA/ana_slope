
import two_step
import ana_new
import numpy as np
import pickle
import matplotlib.pyplot as plt
import bezeir_fit
import matplotlib as mpl
from math import sqrt

std_SQ = pickle.load(open('Steady_fit/SQ_Steady.pick', 'r'))
SQ = np.zeros([len(std_SQ['steady_s']), 2])
SQ[:, 0] = std_SQ['steady_s']
SQ[:, 1] = std_SQ['steady_q']


def farPoint(st_curve, qt_curve):
    SQ = np.zeros([len(std_SQ['steady_s']), 2])
    SQ[:, 0] = std_SQ['steady_s']
    SQ[:, 1] = std_SQ['steady_q']

    Max_diff_Q = 0.0
    Max_index = 0
    for i in range(1, len(st_curve)-1):
        Q_std = abs(np.interp(st_curve[i, 1], SQ[:, 0], SQ[:, 1])
                    - qt_curve[i-1, 1])
        if Q_std > Max_diff_Q:
            Max_diff_Q = Q_std
            Max_index = i

    print('Storage:', st_curve[Max_index, 0], st_curve[Max_index, 1])
    print('Outflow:', qt_curve[Max_index-1, 0], qt_curve[Max_index-1, 1])
    print(np.interp(st_curve[Max_index, 1], std_SQ['steady_s'],
                    std_SQ['ie'])*1000.0*3600)


def run_farPoint(L, dt, n, S, x_grds, ie0, ie1):
    aa = ana_new.trans_ana(ie0, ie1, L, x_grds, dt, n, S)
    aa.run()
    farPoint(aa.out_s, aa.out_q)


def steadyPoint(st_curve, qt_curve, B_idx, D_idx):
    Min_diff_Q = 100.0
    Min_index = 0
    for i in range(B_idx, D_idx):
        Q_std = abs(np.interp(st_curve[i, 1], SQ[:, 0], SQ[:, 1])
                    - qt_curve[i-1, 1])
        if Q_std < Min_diff_Q:
            Min_diff_Q = Q_std
            Min_index = i
    return Min_index


def clipSecondCurve(ab, bc):
    b_idx, c_idx, d_idx = bcdPoints(ab, bc)

    print(b_idx, c_idx, d_idx)
    SC = np.zeros([d_idx-b_idx, 3])
    SC[:, 0] = ab.q_curve[b_idx:d_idx, 0]
    SC[:, 1] = ab.q_curve[b_idx:d_idx, 1]
    SC[:, 2] = ab.s_curve[b_idx:d_idx, 1]

    return SC


def SWRFitTest():
    c0 = 100.0
    c1 = 5.0
    c2 = [5.5, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
          110.0, 130.0, 150.0, 170.0, 200.0]

    a_curves = list()
    b_curves = list()
    SWR_curves = list()
    SWR_fit = list()
    for i in range(0, len(c2)):
        ca = two_step.two_step_overland(c0, c1, c2[i], 100.0, 201, 0.5, 0.1,
                                        0.01)
        ca.ie1_duration(405)
        cb = ana_new.trans_ana(c1, c2[i], 100.0, 201, 0.5, 0.1, 0.01)
        ca.run()
        cb.run()
        SWR = clipSecondCurve(ca, cb)

        a_curves.append(ca)
        b_curves.append(cb)
        SWR_curves.append(SWR)

    plt.figure
    for j in range(0, len(c2)):
        plt.plot(a_curves[j].q_curve[:, 1], a_curves[j].s_curve[:, 1],
                 color='forestgreen', lw=1.0, ls='-')
        plt.plot(b_curves[j].out_q[:, 1], b_curves[j].out_s[:, 1],
                 color='chocolate', lw=1.0, ls='-')
        plt.plot(SWR_curves[j][:, 1], SWR_curves[j][:, 2],
                 label='{:4.3f}'.format(SWR_curves[j][-1, 1]/SWR_curves[j][0, 1]))
    plt.legend(loc=1)
    plt.xlabel('Storage ($m^2$)')
    plt.ylabel('Outflow ($m^2/s$)')
    plt.savefig('multipleSWR.png')
    plt.close()

    #                                      #
    #  Fit Second Wetting Curve From Above #
    #                                      #
    SWR_org = list()
    SWR_fit = list()
    for j in range(0, len(c2)):
        SWR = SWR_curves[j]

        # Normallization by QS(0, 0)
        SWR[:, 0] = (SWR[:, 0] - SWR[0, 0])/(SWR[-1, 0] - SWR[0, 0])
        SWR[:, 1] = SWR[:, 1]/SWR[0, 1]
        SWR[:, 2] = SWR[:, 2]/SWR[0, 2]

        t = np.arange(0, 1.01, 1./20)

        _SWR = np.zeros([21, 2])
        _SWR[:, 0] = np.interp(t, SWR[:, 0], SWR[:, 1])
        _SWR[:, 1] = np.interp(t, SWR[:, 0], SWR[:, 2])
        [p, z] = bezeir_fit.fit(_SWR, 2, 1.0E-8)
        p[0] = [1.0, 1.0]
        n_SWR = bezeir_fit.gen(t, 2, p)

        SWR_fit.append(n_SWR)
        SWR_org.append(_SWR)

    plt.figure(figsize=(10, 8))
    ax = plt.subplot(111)
    for k in range(0, len(c2)):
        ax.plot(SWR_fit[k][:, 0], SWR_fit[k][:, 1], lw=2.0, ls='-',
                label='{:4.3f}'.format(SWR_fit[k][-1, 0]))
        ax.plot(SWR_org[k][:, 0], SWR_org[k][:, 1], lw=1.0, ls='--',
                label='{:4.3f}'.format(SWR_org[k][-1, 0]))

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), ncol=2)
    ax.set_xlabel('q')
    ax.set_ylabel('s')
    plt.title('Secondary Wetting Curve - Normalized and Fitting')
    plt.savefig('SWR_fitting.png')
    plt.close()


def SDRFit(c0, c1, c2, ToC, L, x_space, dt, n, slope, deg):
    a_curves = list()
    b_curves = list()
    SDR_curves = list()
    SDR_fit = list()
    for i in range(0, len(c2)):
        ca = two_step.two_step_overland(c0, c1, c2[i], L, x_space, dt, n,
                                        slope)
        ca.ie1_duration(ToC)
        cb = ana_new.trans_ana(c1, c2[i], L, x_space, dt, n, slope)
        ca.run()
        cb.run()
        SDR = clipSecondCurve(ca, cb)

        a_curves.append(ca)
        b_curves.append(cb)
        SDR_curves.append(SDR)

    plt.figure
    for j in range(0, len(c2)):
        plt.plot(a_curves[j].q_curve[:, 1], a_curves[j].s_curve[:, 1],
                 color='forestgreen', lw=1.0, ls='-')
        plt.plot(b_curves[j].out_q[:, 1], b_curves[j].out_s[:, 1],
                 color='chocolate', lw=1.0, ls='-')
        plt.plot(SDR_curves[j][:, 1], SDR_curves[j][:, 2],
                 label='{:4.3f}'.format(SDR_curves[j][-1, 1]/SDR_curves[j][0, 1]))
    plt.legend(loc=1)
    plt.xlabel('Storage ($m^2$)')
    plt.ylabel('Outflow ($m^2/s$)')
    plt.savefig('multipleSDR.png')
    plt.close()

    #                                      #
    #  Fit Second Wetting Curve From Above #
    #                                      #
    SDR_org = list()
    SDR_fit = list()

    header = list()
    ptList = list()
    for j in range(0, len(c2)):
        SDR = SDR_curves[j]

        # Normallization by QS(0, 0)
        SDR[:, 0] = (SDR[:, 0] - SDR[0, 0])/(SDR[-1, 0] - SDR[0, 0])
        SDR[:, 1] = SDR[:, 1]/SDR[0, 1]
        SDR[:, 2] = SDR[:, 2]/SDR[0, 2]

        t = np.arange(0, 1.01, 1./20)

        _SDR = np.zeros([21, 2])
        _SDR[:, 0] = np.interp(t, SDR[:, 0], SDR[:, 1])
        _SDR[:, 1] = np.interp(t, SDR[:, 0], SDR[:, 2])
        [p, z] = bezeir_fit.fit(_SDR, deg, 1.0E-8)
        p[0] = [1.0, 1.0]
        n_SDR = bezeir_fit.gen(t, deg, p)

        SDR_fit.append(n_SDR)
        SDR_org.append(_SDR)
        header.append(n_SDR[-1, 0]/n_SDR[0, 0])
        ptList.append([p, z])

    plt.figure(figsize=(10, 8))
    ax = plt.subplot(111)
    for k in range(0, len(c2)):
        ax.plot(SDR_fit[k][:, 0], SDR_fit[k][:, 1], lw=2.0, ls='-',
                label='{:4.3f}'.format(SDR_fit[k][-1, 0]))
        ax.plot(SDR_org[k][:, 0], SDR_org[k][:, 1], lw=1.0, ls='--',
                label='{:4.3f}'.format(SDR_org[k][-1, 0]))

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), ncol=2)
    ax.set_xlabel('q')
    ax.set_ylabel('s')
    plt.title('Secondary Drying Curve - Normalized and Fitting')
    plt.savefig('SDR_fitting.png')
    plt.close()


def bcdPoints(ab, ac):
    #  ab: Two-step transition overland flow object
    #  ac: Two rainfall transition overland flow object
    point_B_idx = np.where(ab.s_curve[:, 0] == ab.t1)[0][0]

    point_D_idx = 10
    for i in range(point_B_idx, len(ab.q_curve)):
        if ab.ie1 > ab.ie2:
            q_1 = np.interp(ab.s_curve[i, 1],
                            np.fliplr([ac.out_s[:, 1]])[0],
                            np.fliplr([ac.out_q[:, 1]])[0])
        else:
            q_1 = np.interp(ab.s_curve[i, 1],
                            ac.out_s[:, 1],
                            ac.out_q[:, 1])
        if abs(q_1 - ab.q_curve[i, 1]) <= 1.0e-6:
            point_D_idx = i
            break

    point_C_idx = steadyPoint(ab.s_curve, ab.q_curve, point_B_idx, point_D_idx)
    #  Index of the point on steady-state SQ

    return point_B_idx, point_C_idx, point_D_idx


def Normalize(q_clip, s_clip):
    q_clip[:, 0] = (q_clip[:, 0] - q_clip[0, 0])/(q_clip[-1, 0] - q_clip[0, 0])
    s_clip[:, 0] = (s_clip[:, 0] - s_clip[0, 0])/(s_clip[-1, 0] - s_clip[0, 0])

    q_clip[:, 1] = (q_clip[:, 1] - min(q_clip[:, 1]))/(max(q_clip[:, 1]) -
                                                       min(q_clip[:, 1]))

    s_clip[:, 1] = (s_clip[:, 1] - min(s_clip[:, 1]))/(max(s_clip[:, 1]) -
                                                       min(s_clip[:, 1]))

    return [q_clip, s_clip]


def Normalize2(q_clip, s_clip):
    q_clip[:, 0] = (q_clip[:, 0] - q_clip[0, 0])/(q_clip[-1, 0] - q_clip[0, 0])
    s_clip[:, 0] = (s_clip[:, 0] - s_clip[0, 0])/(s_clip[-1, 0] - s_clip[0, 0])

    q_clip[:, 1] = q_clip[:, 1]/q_clip[0, 1]

    s_clip[:, 1] = s_clip[:, 1]/s_clip[0, 1]

    return [q_clip, s_clip]


class SDRCurves:
    def __init__(self, L, x_space, dt, n, slope):
        self.alpha = 1./n*sqrt(slope)
        self.L = L
        self.n = n
        self.slope = slope
        self.dt = dt

    def runSDRFit(self):
        ie0 = 5.0
        ie1 = 300.0
        ie2 = [280.0, 260.0, 240.0, 220.0, 200.0, 180.0, 160.0, 140.0, 120.0,
               100.0, 80, 60, 40, 20.0, 10.0]
        TimeOfChange = 400.0
        deg = 2
        header, PTs = SDRFit(ie0, ie1, ie2, TimeOfChange, self.L, self.x_space,
                             self.n, self.slope, deg)
        self.header = header
        self.PTs = PTs

    def fitPT(self):
        PTs = self.PTs
        header = self.header
