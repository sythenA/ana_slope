
import two_step
import ana_new
import numpy as np
import pickle
from steady_fit import fit_steady
from scipy.interpolate import interp1d, interp2d, Rbf
import matplotlib.pyplot as plt
from operator import itemgetter
import bezeir_fit
from random import randint, uniform
from math import sqrt

with open('Steady_fit/SQ_Steady.pick', 'r') as f:
    std_SQ = pickle.load(f)

SQ = np.zeros([len(std_SQ[b'steady_s']), 2])
SQ[:, 0] = std_SQ[b'steady_s']
SQ[:, 1] = std_SQ[b'steady_q']


def farPoint(st_curve, qt_curve):
    SQ = np.zeros([len(std_SQ[b'steady_s']), 2])
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


def randColor():
    color = (uniform(0., 1.0), uniform(0., 1.0), uniform(0., 1.0))
    return color


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


def SWRFit(i0, i1, i2, L, x_space, dt, n, slope, deg):
    c1 = list()
    qBDRatio = list()
    sBDRatio = list()
    ptList = list()

    plt.figure
    plt.figure(figsize=(10, 8))
    ax = plt.subplot(111)
    dat_And_Indicator = list()
    # Index 0: qB/qD
    # Index 1: sB/sD
    # Index 2: C1
    # Index 3: [p, t]  Secondary Curve Fit
    for i in range(0, len(i1)):
        ie1 = i1[i]
        N = 5./3  # Manning Eq
        alpha = 1./n*sqrt(slope)  # Manning n
        steady_time = ((ie1/3600./1000.)**(1-N)*L/alpha)**(1/N)
        ToC = [steady_time*0., steady_time*0.001, steady_time*0.01,
               steady_time*0.05, steady_time*0.1, steady_time*0.3,
               steady_time*0.5, steady_time*0.55, steady_time*0.6,
               steady_time*0.65, steady_time*0.7, steady_time*0.75,
               steady_time*0.8, steady_time*0.9, steady_time*0.97,
               steady_time*0.99]
        for j in range(0, len(i2[i])):
            ie2 = i2[i][j]
            q2 = ie2/3600./1000.*L
            s2 = 1./(N+1)*alpha/(ie2/3600./1000.)*(
                (ie2/3600./1000.)*L/alpha)**((N+1)/N)
            for k in range(0, len(ToC)):
                d = ToC[k]
                ca = two_step.two_step_overland(i0, ie1, ie2, L, x_space, dt,
                                                n, slope)
                ca.ie1_duration(d)
                ca.run()

                SWR = bcdPoints(ca)  # Clip Secondary Curve
                SWR = np.array(SWR)
                qBD = SWR[-1, 1]/SWR[0, 1]
                sBD = SWR[-1, 2]/SWR[0, 2]

                #                                      #
                #  Fit Second Drying Curve From Above  #
                #                                      #
                SWR[:, 0] = (SWR[:, 0] - SWR[0, 0])/(SWR[-1, 0] - SWR[0, 0])
                SWR[:, 1] = SWR[:, 1]/SWR[0, 1]
                SWR[:, 2] = SWR[:, 2]/SWR[0, 2]

                t = np.arange(0, 1.01, 1./20)

                _SWR = np.zeros([21, 2])
                _SWR[:, 0] = np.interp(t, SWR[:, 0], SWR[:, 1])
                _SWR[:, 1] = np.interp(t, SWR[:, 0], SWR[:, 2])

                # Fitting
                [p, z] = bezeir_fit.fit(_SWR, deg, 1.0E-8)
                p[0] = [1.0, 1.0]
                p[-1] = [SWR[-1, 1], SWR[-1, 2]]
                n_SWR = bezeir_fit.gen(t, deg, p)

                dat_And_Indicator.append([qBD,
                                          sBD,
                                          ie1/ie2,
                                          [p, z]])

                color = randColor()
                labelString = ('$q_B/q_D$=' +
                               '{:4.3f}'.format(qBD) + ', ' +
                               '$c_1$=' +
                               '{:4.3f}'.format(ie1/ie2))
                print(labelString)

                ax.plot(SWR[:, 1], SWR[:, 2], lw=2.0, ls='-',
                        label=labelString, color=color)
                ax.plot(n_SWR[:, 0], n_SWR[:, 1], lw=1.0, ls='--', color=color)

    # dat_And_Indicator = sorted(dat_And_Indicator, key=itemgetter(2))
    for i in range(0, len(dat_And_Indicator)):
        qBDRatio.append(dat_And_Indicator[i][0])
        sBDRatio.append(dat_And_Indicator[i][1])
        c1.append(dat_And_Indicator[i][2])
        ptList.append(dat_And_Indicator[i][3])

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), ncol=2)
    ax.set_xlabel('q')
    ax.set_ylabel('s')
    plt.title('Secondary Wetting Curve - Normalized and Fitting')
    plt.savefig('SWR_fitting.png')
    plt.close()

    f = open('ks1_SWRrec.pick', 'wb')
    pickle.dump({'qBD': qBDRatio, 'sBD': sBDRatio, 'PT': ptList, 'C1': c1}, f)
    f.close()


def SDRFit(i0, i1, i2, L, x_space, dt, n, slope, deg):
    endSlope = list()
    qBi1Ratio = list()
    SBi1Ratio = list()
    qDi1Ratio = list()
    SDi1Ratio = list()
    ptList = list()

    plt.figure
    plt.figure(figsize=(10, 8))
    ax = plt.subplot(111)
    dat_And_Indicator = list()
    # Index 0: qB/qD
    # Index 1: sB/sD
    # Index 2: EndSlope
    # Index 3: [p, t]  Secondary Curve Fit
    for i in range(0, len(i1)):
        ie1 = i1[i]
        N = 5./3  # Manning Eq
        alpha = 1./n*sqrt(slope)  # Manning n
        steady_time = ((ie1/3600./1000.)**(1-N)*L/alpha)**(1/N)
        ToC = [steady_time*0., steady_time*0.001, steady_time*0.01,
               steady_time*0.05, steady_time*0.1, steady_time*0.3,
               steady_time*0.5, steady_time*0.55, steady_time*0.6,
               steady_time*0.65, steady_time*0.7, steady_time*0.75,
               steady_time*0.8, steady_time*0.9, steady_time*0.97,
               steady_time*0.99]
        q1 = ie1/3600./1000.*L
        s1 = N/(N+1)*alpha/(ie1/3600./1000.)*(
            (ie1/3600./1000.)*L/alpha)**((N+1)/N)
        for j in range(0, len(i2[i])):
            ie2 = i2[i][j]
            for k in range(0, len(ToC)):
                d = ToC[k]
                ca = two_step.two_step_overland(i0, ie1, ie2, L, x_space, dt,
                                                n, slope)
                ca.ie1_duration(d)
                ca.run()

                SDR, _endSlope = bcdPoints(ca)  # Clip Secondary Curve
                SDR = np.array(SDR)
                qBi1 = SDR[0, 1]/q1
                SBi1 = SDR[0, 2]/s1
                qDi1 = SDR[-1, 1]/q1
                SDi1 = SDR[-1, 2]/s1
                qBD = SDR[0, 1]/SDR[-1, 1]

                #                                      #
                #  Fit Second Drying Curve From Above  #
                #                                      #
                SDR[:, 0] = (SDR[:, 0] - SDR[0, 0])/(SDR[-1, 0] - SDR[0, 0])
                SDR[:, 1] = SDR[:, 1]/SDR[0, 1]
                SDR[:, 2] = SDR[:, 2]/SDR[0, 2]

                t = np.arange(0, 1.01, 1./20)

                _SDR = np.zeros([21, 2])
                _SDR[:, 0] = np.interp(t, SDR[:, 0], SDR[:, 1])
                _SDR[:, 1] = np.interp(t, SDR[:, 0], SDR[:, 2])

                # Fitting
                [p, z] = bezeir_fit.fit(_SDR, deg, 1.0E-8)
                p[0] = [1.0, 1.0]
                p[-1] = [SDR[-1, 1], SDR[-1, 2]]
                n_SDR = bezeir_fit.gen(t, deg, p)

                dat_And_Indicator.append([qBi1,
                                          SBi1,
                                          qDi1,
                                          SDi1,
                                          _endSlope,
                                          [p, z]])

                color = randColor()
                labelString = ('$q_B/q_D$=' +
                               '{:4.3f}'.format(qBD) + ', ' +
                               '$c_1$=' +
                               '{:4.3f}'.format(ie1/ie2))
                print(labelString)

                ax.plot(SDR[:, 1], SDR[:, 2], lw=2.0, ls='-',
                        label=labelString, color=color)
                ax.plot(n_SDR[:, 0], n_SDR[:, 1], lw=1.0, ls='--', color=color)

    # dat_And_Indicator = sorted(dat_And_Indicator, key=itemgetter(2))
    for i in range(0, len(dat_And_Indicator)):
        qBi1Ratio.append(dat_And_Indicator[i][0])
        SBi1Ratio.append(dat_And_Indicator[i][1])
        qDi1Ratio.append(dat_And_Indicator[i][2])
        SDi1Ratio.append(dat_And_Indicator[i][3])
        endSlope.append(dat_And_Indicator[i][4])
        ptList.append(dat_And_Indicator[i][5])

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), ncol=2)
    ax.set_xlabel('q')
    ax.set_ylabel('s')
    plt.title('Secondary Drying Curve - Normalized and Fitting')
    plt.savefig('SDR_fitting.png')
    plt.close()

    f = open('ks1_SDRrec.pick', 'wb')
    pickle.dump({'qBi1': qBi1Ratio, 'SBi1': SBi1Ratio, 'qDi1': qDi1Ratio,
                 'SDi1': SDi1Ratio, 'PT': ptList, 'endSlope': endSlope}, f)
    f.close()


def bcdPoints(ab):
    #  ab: Two-step transition overland flow object
    #  ac: Two rainfall transition overland flow object
    def getBIdx(ab, d):
        minIdx = 0
        minErr = 1000
        for j in range(0, len(ab.q_curve)):
            err = abs(ab.q_curve[j, 0] - d)
            if err < minErr:
                minIdx = j
                minErr = err

        return minIdx

    i1 = ab.ie1
    i2 = ab.ie2

    d = ab.d
    m = 5./3
    alpha = ab.alpha
    L = ab.L
    c1 = i1/i2

    tr = (i2**(1-m)/alpha*(L+(i1*d)**m*(alpha/i2 - alpha/i1)))**(1/m)-c1*d

    y = i1*d + i2*tr
    qD = alpha*y**m
    sD = np.interp(d+tr, ab.s_curve[:, 0], ab.s_curve[:, 1])

    qB = np.interp(d, ab.q_curve[:, 0], ab.q_curve[:, 1])
    sB = np.interp(d, ab.s_curve[:, 0], ab.s_curve[:, 1])

    qDp1 = np.interp(d+tr+0.5, ab.q_curve[:, 0], ab.q_curve[:, 1])
    sDp1 = np.interp(d+tr+0.5, ab.s_curve[:, 0], ab.s_curve[:, 1])
    endSlope = (sDp1 - sD)/(qDp1 - qD)

    B_idx = getBIdx(ab, d)
    D_idx = getBIdx(ab, d+tr)

    t_clip = ab.q_curve[B_idx:D_idx+1, 0]
    s_curve = ab.s_curve[B_idx:D_idx+1, 1]
    q_curve = ab.q_curve[B_idx:D_idx+1, 1]
    s_curve[0] = sB
    s_curve[-1] = sD
    q_curve[0] = qB
    q_curve[-1] = qD

    sq_curve = [zip(t_clip, q_curve, s_curve)][0]

    return list(sq_curve), endSlope


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

    def plotFit(self):
        PPlist = self.PPlist
        qBDRatio = self.qBDRatio
        sBDRatio = self.sBDRatio
        PTs = self.PTs
        C1 = self.C1

        for j in range(0, len(PPlist)):
            p = list()
            org_P = list()
            fx = PPlist[j][0]
            fy = PPlist[j][1]
            for i in range(0, len(qBDRatio)):
                qBD = qBDRatio[i]
                sBD = sBDRatio[i]
                val_C1 = C1[i]

                p.append([fx(qBD, val_C1), fy(qBD, sBD)])
                org_P.append(PTs[i][0][j+1])
            p = np.array(p)
            org_P = np.array(org_P)

            plt.figure
            plt.plot(p[:, 0], p[:, 1], 'r+')
            plt.plot(org_P[:, 0], org_P[:, 1], 'bo', markerfacecolor='None')
            plt.savefig('point_'+str(j)+'.png', dpi=150)
            plt.close()

    def runSDRFit(self):
        ie0 = 10.0
        ie1 = [15., 20., 50., 100.0, 200.0, 300.0]
        ie2 = [[10., 5.],
               [12., 10., 5., 3.],
               [45., 40., 35., 30., 25., 20., 15., 10., 5.],
               [80, 60.0, 40., 20., 10., 5.],
               [180., 160., 140., 120., 100., 80., 60., 40., 20., 10., 5.],
               [280, 260, 240, 220, 200, 180, 160, 140.0, 120.0, 100.0, 80, 60,
               40, 20., 10., 5.]]
        SDRFit(ie0, ie1, ie2, self.L, self.x_space, self.dt,
               self.n, self.slope, self.degC)

        dat = pickle.load(open('ks1_SDRrec.pick', 'rb'))
        self.qBi1Ratio = dat['qBi1']
        self.SBi1Ratio = dat['SBi1']
        self.qBi1Ratio = dat['qDi1']
        self.DBi1Ratio = dat['SDi1']
        self.endSlope = dat['endSlope']
        self.PTs = dat['PT']

    def readRec(self, recFile):
        a = pickle.load(open(recFile, 'rb'))
        self.qBi1Ratio = a[b'qBi1']
        self.SBi1Ratio = a[b'SBi1']
        self.qDi1Ratio = a[b'qDi1']
        self.SDi1Ratio = a[b'SDi1']
        self.endSlope = a[b'endSlope']
        self.PTs = a[b'PT']

    def fitPT(self):
        PTs = self.PTs

        degC = len(PTs[0][0])-1
        pointsContainer = list()
        for i in range(0, degC):
            pointsList_X = list()
            pointsList_Y = list()
            for j in range(0, len(PTs)):
                pointsList_X.append(PTs[j][0][i+1][0])
                pointsList_Y.append(PTs[j][0][i+1][1])
            pointsContainer.append([pointsList_X, pointsList_Y])

        # Container of points for generating interpolated curve
        PPList = list()
        for i in range(0, degC):
            X_list = pointsContainer[i][0]
            Y_list = pointsContainer[i][1]
            # Two-dimensional interpolation function
            fX = interp2d(self.qBDRatio, self.sBDRatio, X_list, kind='cubic')
            fY = interp2d(self.qBDRatio, self.sBDRatio, Y_list, kind='cubic')

            PPList.append([fX, fY])

        self.PPlist = PPList

    def nearestNeighbor(self, p):
        distance = np.zeros(len(self.qBi1Ratio))
        idx = np.arange(0, len(self.qBi1Ratio))
        Pts = zip(self.qBi1Ratio, self.SBi1Ratio, self.qDi1Ratio,
                  self.SDi1Ratio, distance, idx)
        Pts = np.array(Pts)

        for i in range(0, len(Pts)):
            dist = sqrt((p[0] - Pts[i][0])**2 + (p[1] - Pts[i][1])**2 +
                        (p[2] - Pts[i][2])**2 + (p[3] - Pts[i][3])**2)
            Pts[i][4] = dist

        midPoint = list()
        Pts = sorted(Pts, key=itemgetter(4))
        for j in range(0, 10):
            idx = Pts[j][5]
            midPoint.append([self.PTs[int(idx)][0][1][0],
                             self.PTs[int(idx)][0][1][1]])

        return Pts[0:10], midPoint

    def buildInterp(self):
        Px = list()
        Py = list()
        for i in range(0, len(self.PTs)):
            Px.append(self.PTs[i][0][1][0])
            Py.append(self.PTs[i][0][1][1])
        Px = np.array(Px)
        Py = np.array(Py)

        iX = Rbf(self.qBi1Ratio, self.SBi1Ratio, self.qDi1Ratio,
                 self.SDi1Ratio, Px)
        iY = Rbf(self.qBi1Ratio, self.SBi1Ratio, self.qDi1Ratio,
                 self.SDi1Ratio, Py)

        self.iX = iX
        self.iY = iY

    def inverseDist(self, p, Pts, midPoints):
        distList = list()
        for i in range(0, len(Pts)):
            dist = sqrt((p[0]-Pts[i][0])**2 + (p[1]-Pts[i][1])**2 +
                        (p[2]-Pts[i][2])**2 + (p[3]-Pts[i][3])**2)
            distList.append(dist)

        distList = np.array(distList)
        weights = 1./distList*0.5
        weights /= weights.sum(axis=0)

        rX = sum(weights*midPoints[:, 0])
        rY = sum(weights*midPoints[:, 1])

        return [rX, rY]

    def interpCurve(self, p):
        Pts, midPoint = self.nearestNeighbor(p)
        Pts = np.array(Pts)

        midPoint = np.array(midPoint)
        """
        xrbf = Rbf(Pts[:, 0], Pts[:, 1], Pts[:, 2], Pts[:, 3], midPoint[:, 0])
        yrbf = Rbf(Pts[:, 0], Pts[:, 1], Pts[:, 2], Pts[:, 3], midPoint[:, 1])

        P1X = xrbf(p[0], p[1], p[2], p[3])*p[0]
        P1Y = yrbf(p[0], p[1], p[2], p[3])*p[1]"""

        P1X, P1Y = self.inverseDist(p, Pts, midPoint)

        P1X = P1X*p[0]
        P1Y = P1Y*p[1]

        print(midPoint)
        print(P1X, P1Y)

        curvePoints = [[p[0], p[1]], [P1X, P1Y], [p[2], p[3]]]

        T = np.arange(0., 1.0001, 1./100)
        curve = bezeir_fit.gen(T, self.degC, np.array(curvePoints))
        return curve, T


class SWRCurves:
    def __init__(self, L, x_space, dt, n, slope, degC, degP, **kwargs):
        # degC: Degree of bezeir curve in fitting SWR
        # degP: Degree of bezeir curve in fitting points generated in fitting
        # SWR.
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
            self.runSWRFit()

    def runSWRFit(self):
        ie0 = 100.0
        ie1 = [5.0, 20.0, 40.0, 60.0, 80.0]
        ie2 = [[6.0, 10.0, 30.0, 50.0, 70.0, 90.0, 100.0, 120.0, 140.0, 160.0,
                180.0, 200.0],
               [22.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0,
                190.0, 210.0],
               [45.0, 60.0, 80.0, 100.0, 120.0, 140.0, 150.0, 170.0, 190.0,
                210.0],
               [65.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0],
               [85.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0]]
        SWRFit(ie0, ie1, ie2, self.L, self.x_space, self.dt, self.n,
               self.slope, self.degC)

    def readRec(self, recFile):
        a = pickle.load(open(recFile, 'rb'))
        self.qBDRatio = a[b'qBD']
        self.sBDRatio = a[b'sBD']
        self.C1 = a[b'C1']
        self.PTs = a[b'PT']

    def fitPT(self):
        PTs = self.PTs

        degC = len(PTs[0][0])-1
        pointsContainer = list()
        for i in range(0, degC):
            pointsList_X = list()
            pointsList_Y = list()
            for j in range(0, len(PTs)):
                pointsList_X.append(PTs[j][0][i+1][0])
                pointsList_Y.append(PTs[j][0][i+1][1])
            pointsContainer.append([pointsList_X, pointsList_Y])

        # Container of points for generating interpolated curve
        PPList = list()
        for i in range(0, degC):
            X_list = pointsContainer[i][0]
            Y_list = pointsContainer[i][1]
            # Two-dimensional interpolation function
            fX = interp2d(self.qBDRatio, self.sBDRatio, X_list, kind='cubic')
            fY = interp2d(self.qBDRatio, self.sBDRatio, Y_list, kind='cubic')

            PPList.append([fX, fY])

        self.PPlist = PPList

    def interpCurve(self, qBD, sBD):
        PPlist = self.PPlist

        curvePoints = [[1, 1]]
        for i in range(0, len(PPlist)):
            fX = PPlist[i][0]
            fY = PPlist[i][1]

            p = [fX(qBD, sBD), fY(qBD, sBD)]
            curvePoints.append(p)

        # p[-1] = [qBD, sBD]
        T = np.arange(0., 1.0001, 1./100)

        curve = bezeir_fit.gen(T, self.degC, np.array(curvePoints))
        return curve, T


def randColor():
    r = randint(0, 255)
    g = randint(0, 255)
    b = randint(0, 255)

    return '#%02x%02x%02x' % (r, g, b)


def compareSDRCurve():
    L = 100.0
    slope = 0.01
    n = 0.1
    dt = 0.5
    x_space = 201
    a = SDRCurves(L, x_space, dt, n, slope, 2, 4, useFile='ks1_SDRrec.pick')
    a.fitPT()

    # Run curves
    ie0 = 5.0
    ie1 = 300.0
    ie2 = [280.0, 260.0, 240.0, 220.0, 200.0, 180.0, 160.0, 140.0, 120.0,
           100.0, 80, 60, 40, 20.0]
    TimeOfChange = 400.0

    SDR_curves = list()
    gen_Curves = list()
    Qratios = list()
    for i in range(0, len(ie2)):
        ca = two_step.two_step_overland(ie0, ie1, ie2[i], L, x_space, dt, n,
                                        slope)
        ca.ie1_duration(TimeOfChange)
        cb = ana_new.trans_ana(ie1, ie2[i], L, x_space, dt, n, slope)
        ca.run()
        cb.run()
        SDR = bcdPoints(ca)
        SDR[:, 0] = (SDR[:, 0] - SDR[0, 0])/(SDR[-1, 0] - SDR[0, 0])
        SDR[:, 1] = SDR[:, 1]/SDR[0, 1]
        SDR[:, 2] = SDR[:, 2]/SDR[0, 2]

        SDR_curves.append(SDR)

        Qratio = SDR[-1, 1]/SDR[0, 1]
        curve, T = a.interpCurve(Qratio)
        gen_Curves.append(curve)
        Qratios.append(Qratio)

    plt.figure(figsize=(10, 8))
    ax = plt.subplot(111)
    for k in range(0, len(ie2)):
        color = randColor()
        ax.plot(SDR_curves[k][:, 1], SDR_curves[k][:, 2], lw=1.0, ls='--',
                label='{:4.3f}'.format(Qratios[k]), color=color)
        ax.plot(gen_Curves[k][:, 0], gen_Curves[k][:, 1], lw=2.0, ls='-',
                label='{:4.3f}'.format(gen_Curves[k][-1, 0]),
                color=color)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), ncol=2)
    ax.set_xlabel('q')
    ax.set_ylabel('s')
    plt.title('Secondary Drying Curve - Original and Generated')
    plt.savefig('SDR_gen_compare.png')
    plt.close()


def compareSWRCurve():
    ie0 = 100.0
    ie1 = 5.0
    ie2 = [5.5, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0,
           100.0, 110.0, 130.0, 150.0, 170.0, 200.0]
    L = 100.0
    slope = 0.01
    n = 0.1
    dt = 0.5
    x_space = 201
    a = SWRCurves(L, x_space, dt, n, slope, 2, 4, useFile='ks1_SWRrec.pick')
    a.fitPT()

    # Run curves
    TimeOfChange = 405.0

    SWR_curves = list()
    gen_Curves = list()
    Qratios = list()
    for i in range(0, len(ie2)):
        ca = two_step.two_step_overland(ie0, ie1, ie2[i], L, x_space, dt, n,
                                        slope)
        ca.ie1_duration(TimeOfChange)
        cb = ana_new.trans_ana(ie1, ie2[i], L, x_space, dt, n, slope)
        ca.run()
        cb.run()
        SWR = bcdPoints(ca)
        SWR[:, 0] = (SWR[:, 0] - SWR[0, 0])/(SWR[-1, 0] - SWR[0, 0])
        SWR[:, 1] = SWR[:, 1]/SWR[0, 1]
        SWR[:, 2] = SWR[:, 2]/SWR[0, 2]

        SWR_curves.append(SWR)

        Qratio = SWR[-1, 1]/SWR[0, 1]
        curve, T = a.interpCurve(Qratio)
        gen_Curves.append(curve)
        Qratios.append(Qratio)

    plt.figure(figsize=(10, 8))
    ax = plt.subplot(111)
    for k in range(0, len(ie2)):
        color = randColor()
        ax.plot(SWR_curves[k][:, 1], SWR_curves[k][:, 2], lw=1.0, ls='--',
                label='{:4.3f}'.format(Qratios[k]), color=color)
        ax.plot(gen_Curves[k][:, 0], gen_Curves[k][:, 1], lw=2.0, ls='-',
                label='{:4.3f}'.format(gen_Curves[k][-1, 0]),
                color=color)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), ncol=2)
    ax.set_xlabel('q')
    ax.set_ylabel('s')
    plt.title('Secondary Wetting Curve - Original and Generated')
    plt.savefig('SWR_gen_compare.png')
    plt.close()
