
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


def SWRFit():
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

    header = list()
    ptList = list()
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
        header.append(n_SWR[-1, 0]/n_SWR[0, 0])
        ptList.append([p, z])

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

    f = open('ks1_SWRrec.pick', 'w')
    pickle.dump({'header': header, 'PT': ptList}, f)

    return header, ptList


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

    f = open('ks1_SDRrec.pick', 'w')
    pickle.dump({'header': header, 'PT': ptList}, f)

    return header, ptList


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


def testFitSDR():
    a = pickle.load(open('SDRrec.pick', 'r'))
    header = a['header']
    PT = a['PT']

    Plist = list()
    p2List = list()
    p3List = list()
    for j in range(0, len(PT)):
        Plist.append(PT[j][0])
    for k in range(0, len(Plist)):
        p2List.append(Plist[k][1])
        p3List.append(Plist[k][2])
    p2List = np.array(p2List)
    p3List = np.array(p3List)

    P2, T2 = bezeir_fit.fit(p2List, 4, 1.0E-8)
    P3, T3 = bezeir_fit.fit(p3List, 4, 1.0E-8)

    T = np.arange(0., 1.001, 1./50)
    P2_fit = bezeir_fit.gen(T, 4, P2)
    P3_fit = bezeir_fit.gen(T, 4, P3)

    plt.figure
    plt.plot(P2_fit[:, 0], P2_fit[:, 1], label='P2', lw=3.0)
    plt.plot(p2List[:, 0], p2List[:, 1], ls='--', lw=3.0)
    plt.plot(P3_fit[:, 0], P3_fit[:, 1], label='P3')
    plt.plot(p3List[:, 0], p3List[:, 1], ls='--')
    plt.legend()
    plt.savefig('p2p3fit.png')
    plt.close()

    t2 = np.interp([1.277], np.fliplr([header])[0], np.fliplr([T2])[0])
    t3 = np.interp([1.277], np.fliplr([header])[0], np.fliplr([T3])[0])
    _p2 = bezeir_fit.gen(t2, 4, P2)[0]
    _p3 = bezeir_fit.gen(t3, 4, P3)[0]

    _dat1 = bezeir_fit.gen(T, 2, [[1, 1], _p2, _p3])
    _dat2 = bezeir_fit.gen(T, 2, PT[9][0])
    plt.figure
    plt.plot(_dat1[:, 0], _dat1[:, 1], label='gen')
    plt.plot(_dat2[:, 0], _dat2[:, 1], label='rec')
    plt.legend()
    plt.show()


class SDRCurves:
    def __init__(self, L, x_space, dt, n, slope, degC, degP, **kwargs):
        # degC: Degree of bezeir curve in fitting SDR
        # degP: Degree of bezeir curve in fitting points generated in fitting
        # SDR.
        self.alpha = 1./n*sqrt(slope)
        self.L = L
        self.n = n
        self.slope = slope
        self.x_space = x_space
        self.degC = degC
        self.degP = degP
        self.dt = dt

        if 'useFile' in kwargs.keys():
            recfile = kwargs['useFile']
            self.readRec(recfile)
        else:
            self.runSDRFit()

    def runSDRFit(self):
        ie0 = 5.0
        ie1 = 300.0
        ie2 = [280.0, 260.0, 240.0, 220.0, 200.0, 180.0, 160.0, 140.0, 120.0,
               100.0, 80, 60, 40, 20.0]
        TimeOfChange = 400.0
        header, PTs = SDRFit(ie0, ie1, ie2, TimeOfChange, self.L, self.x_space,
                             self.dt, self.n, self.slope, self.degC)
        self.header = header
        self.PTs = PTs

    def readRec(self, recFile):
        a = pickle.load(open(recFile, 'r'))
        self.header = a['header']
        self.PTs = a['PT']

    def fitPT(self):
        PTs = self.PTs

        degC = len(PTs[0][0])-1
        self.degC = degC
        pointsContainer = list()
        for i in range(0, degC):
            pointsList = list()
            for j in range(0, len(PTs)):
                pointsList.append(PTs[j][0][i+1])
            pointsContainer.append(pointsList)

        Tlist = list()
        PPList = list()  # Container of points for generating interpolated curve
        for i in range(0, degC):
            p, t = bezeir_fit.fit(np.array(pointsContainer[i]), self.degP,
                                  1.0E-8)
            PPList.append(p)
            Tlist.append(t)

        self.PPlist = PPList
        self.Tlist = Tlist

    def interpCurve(self, gamma):
        header = self.header
        PPlist = self.PPlist
        Tlist = self.Tlist

        curvePoints = [[1, 1]]
        for i in range(0, len(PPlist)):
            t = np.interp([gamma], np.fliplr([header])[0],
                          np.fliplr([Tlist[i]])[0])[0]
            p = bezeir_fit.gen([t], self.degP, PPlist[i])[0]
            curvePoints.append(p)
        T = np.arange(0., 1.0001, 1./100)
        curve = bezeir_fit.gen(T, self.degC, np.array(curvePoints))
        return curve, T



class SWRCurves:
    def __init__(self, L, x_space, dt, n, slope, degP, **kwargs):
        # degC: Degree of bezeir curve in fitting SWR
        # degP: Degree of bezeir curve in fitting points generated in fitting
        # SWR.
        self.alpha = 1./n*sqrt(slope)
        self.L = L
        self.n = n
        self.slope = slope
        self.x_space = x_space
        self.degP = degP
        self.dt = dt

        if 'useFile' in kwargs.keys():
            recfile = kwargs['useFile']
            self.readRec(recfile)
        else:
            self.runSWRFit()

    def runSWRFit(self):
        ie0 = 5.0
        ie1 = 300.0
        ie2 = [280.0, 260.0, 240.0, 220.0, 200.0, 180.0, 160.0, 140.0, 120.0,
               100.0, 80, 60, 40, 20.0]
        TimeOfChange = 400.0
        header, PTs = SWRFit(ie0, ie1, ie2, TimeOfChange, self.L, self.x_space,
                             self.dt, self.n, self.slope, self.degC)
        self.header = header
        self.PTs = PTs

    def readRec(self, recFile):
        a = pickle.load(open(recFile, 'r'))
        self.header = a['header']
        self.PTs = a['PT']

    def fitPT(self):
        PTs = self.PTs

        degC = len(PTs[0])-1
        self.degC = degC
        for i in range(0, degC):
            pointsList = list()
            for j in range(0, len(PTs)):
                pointsList.append(PTs[j][0][i+1])

        Tlist = list()
        PPList = list()  # Container of points for generating interpolated curve
        for i in range(0, degC):
            p, t = bezeir_fit.fit(pointsList[i], self.degP, 1.0E-8)
            PPList.append(p)
            Tlist.append(t)

        self.PPlist = PPList
        self.Tlist = Tlist

    def interpCurve(self, gamma):
        header = self.header
        PPlist = self.PPlist
        Tlist = self.Tlist

        curvePoints = [[1, 1]]
        for i in range(0, len(PPlist)):
            t = np.interp(gamma, header, Tlist[i])[0]
            p = bezeir_fit.gen(t, self.degP, PPlist[i])[0][0]
            curvePoints.append(p)
        T = np.arange(0., 1.0001, 1./100)
        curve = bezeir_fit.gen(T, self.degC, curvePoints)
        return curve, T


def interpS(S, S_curve, T):
    # Return all possible value in t corresponds to storage in a SDR curve.
    sol = list()
    for i in range(1, len(T)):
        if (S_curve[i-1] - S)*(S_curve[i] - S) < 0:
            t0 = T[i-1]
            t1 = T[i]
            S1 = S_curve[i]
            S0 = S_curve[i-1]
            # linear interpolation
            t = t0 + (S-S0)/(S1-S0)*(t1-t0)
            sol.append(t)

    return sol


def interpQtoT(Q, t0, Q_curve, T):
    # Return the closiest value of q in bezeir-curve-fit.
    sol = list()
    for i in range(1, len(T)):
        if (Q_curve[i-1]-Q)*(Q_curve[i]-Q) < 0:
            t1 = T[i-1]
            t2 = T[i]
            Q1 = Q_curve[i-1]
            Q2 = Q_curve[i]

            t = t1 + (Q-Q1)/(Q2-Q1)*(t2-t1)
            sol.append(t)
    sol = np.array(sol)
    print sol
    sol = sol - t0
    sol.sort()
    print sol
    idx = np.where(sol > 0)

    print idx

    return sol[idx[0]]+t0  # returns T(reference position on bezeir-curve)


def SCQuerry(curve, T, Sc, Q0, t0):
    """
    Sc: the storage to find the outflow
    Q0: last outflow
    T: reference knots
    curve: standard generated bezeir curve of outflow - storage relation.
    t0: the last t-value(indicator of position on a bezeir-curve)
    """
    q = curve[:, 0]
    s = curve[:, 1]

    t0 = interpQtoT(Q0, t0, q, T)
    s_to_t = interpS(Sc, s, T)

    if len(s_to_t) > 1:
        r = np.array(s_to_t) - t0
        r = r.sort()
        idx = np.where(r > 0)
        idx = idx[0]

        t = r[idx] + t0
        Q1 = np.interp(t, T, q)
    elif len(s_to_t) == 1:
        t = s_to_t[0]
        Q1 = np.interp(t, T, q)

    return Q1, t
