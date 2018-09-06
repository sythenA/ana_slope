
import numpy as np
import matplotlib.pyplot as plt
from SQsearch import RisingLimb, FallingLimb
import pickle
import ana_new as an
import steady_fit as st
from farPoint import SWRCurves, SDRCurves


def rainParse(rain_rec, current_t):
    return np.interp(current_t, rain_rec[:, 0], rain_rec[:, 1])/3600.0*10**-3


def rainParse2(rain_rec, current_t, dt):
    t = np.arange(rain_rec[0, 0], rain_rec[-1, 0]+1, dt)
    new_rec = np.interp(t, rain_rec[:, 0], rain_rec[:, 1])
    if current_t == 0.0:
        return rain_rec[0, 1]/3600.*10**-3
    elif current_t >= rain_rec[-1, 0]:
        return rain_rec[-1, 1]/3600.*10**-3

    for i in range(0, len(new_rec)-1):
        if current_t >= t[i] and current_t < t[i+1]:
            return new_rec[i]/3600.*10**-3


def ie_select(t1, t2, rain_rec):
    ie0 = rainParse(rain_rec, t1)
    ie1 = rainParse(rain_rec, t2)

    if ie0 == ie1:
        return ['nan', ie1]
    else:
        return [ie0, ie1]


def h1(length, ie, k, normS, dt, std, Fall, Rise):
    max_Q = std.ie_to_outflow(ie)
    dS = length*dt*ie
    if k >= 1.0:
        H1 = dS - Rise.s_to_q(k, normS)*max_Q*dt
    elif k < 1.0:
        H1 = dS - Fall.s_to_q(k, normS)*max_Q*dt

    return H1


def h2(length, ie, k, normS, h1, dt, std, Fall, Rise):
    max_S = std.ie_to_storage(ie)
    max_Q = std.ie_to_outflow(ie)
    dI = length*(1 + 1.0/3)*dt*ie
    if k >= 1.0:
        H2 = dI - Rise.s_to_q(k, normS + (1.0/3)*h1/max_S)*max_Q*dt
    elif k < 1.0:
        H2 = dI - Fall.s_to_q(k, normS + (1.0/3)*h1/max_S)*max_Q*dt

    return H2


def h3(length, ie, k, normS, h2, dt, std, Fall, Rise):
    max_S = std.ie_to_storage(ie)
    max_Q = std.ie_to_outflow(ie)
    dI = length*(1 + 2.0/3)*dt*ie
    if k >= 1.0:
        H3 = dI - Rise.s_to_q(k, normS + (2.0/3)*h2/max_S)*max_Q*dt
    elif k < 1.0:
        H3 = dI - Fall.s_to_q(k, normS + (2.0/3)*h2/max_S)*max_Q*dt

    return H3


def r_h1(length, ie, S, dt, std):
    dS = length*ie*dt
    _ie = std.s_to_ie(S)
    q = std.ie_to_outflow(_ie)

    h1 = dS - q*dt

    return h1


def r_h2(length, ie, S, dt, h1, std):
    dS = length*ie*dt
    _ie = std.s_to_ie(S+1.0/3*h1)
    q = std.ie_to_outflow(_ie)

    h2 = dS*(1 + 1.0/3.0) - q*dt

    return h2


def r_h3(length, ie, S, dt, h2, std):
    dS = length*ie*dt*(1.0 + 2.0/3)
    _ie = std.s_to_ie(S+2.0/3*h2)
    q = std.ie_to_outflow(_ie)

    h3 = dS - q*dt

    return h3


def initialk(rain_rec, dt, std, Rise, Fall):
    # Set initial setting before simulation start
    pre_ie = rainParse(rain_rec, 0.0)
    c_ie = rainParse(rain_rec, dt)
    if pre_ie > 0:
        k = rainParse(rain_rec, dt)/pre_ie  # Initial k setting
        max_S = std.ie_to_storage(max([c_ie, pre_ie]))
        max_Q = std.ie_to_outflow(max([c_ie, pre_ie]))
        if k > 1.:
            currentCurve, t = Rise.getCurve(k, max_S, max_Q)
        elif k < 1.:
            currentCurve, t = Fall.getCurve(k, max_S, max_Q)
        elif k == 0:
            k = 1./100
            currentCurve, t = Fall.getCurve(k, max_S, max_Q)
        elif k == 1:
            currentCurve, t = std.bzCurve()
    elif pre_ie == 0 and c_ie > 0:
        max_S = std.ie_to_storage(max([c_ie, pre_ie]))
        max_Q = std.ie_to_outflow(max([c_ie, pre_ie]))
        k = 100.0
        currentCurve, t = Rise.getCurve(k, max_S, max_Q)
    elif pre_ie == 0 and c_ie == 0:
        currentCurve, t = std.bzCurve()

    return k, currentCurve, t


def curveFromSteadyK(k, maxS, maxQ, Rise, Fall):


def getK(S, Q, pre_ie, c_ie, std, Rise, Fall, SDR, SWR):
    relative_ie = std.s_to_ie(S)
    relative_Q = std.ie_to_outflow(relative_ie)
    maxS = std.ie_to_storage(max([c_ie, pre_ie]))
    maxQ = std.ie_to_outflow(max([c_ie, pre_ie]))
    k1 = c_ie/pre_ie
    if Q < relative_Q*0.999:  # S, Q on the rising side
        if k1 > 1:
            k, t_star = Rise.SQtok(
                S/max_S, Q/max_Q, (std.ie_to_storage(c_ie)/max_S)**(5.0/3))
            curve, T = Rise.getCurve(k, maxS, maxQ)
        elif k1 < 1:  # Call SDC

    elif Q > relative_Q*1.001:  # S, Q on the falling side
        if k1 < 1:
            k, t_star = Fall.SQtok(
                S/max_S, Q/max_Q,
                (std.ie_to_storage(c_ie)/max_S)**(5.0/3) - 0.01)
            curve, T = Fall.getCurve(k, maxS, maxQ)
        elif k1 > 1:  # Call SWC
    print('k=%f' % k)


def calc(rain_rec, length):
    dt = 20.0
    Fall = FallingLimb(rec_name='ks1', folder='Steady_fit')
    Rise = RisingLimb(rec_name='ks1', folder='Steady_fit')

    # Secondary curves
    SDR = SDRCurves(100.0, 201, 0.5, 0.1, 0.01, 2, 4,
                    useFile='ks1_SDRrec.pick')
    SDR.fitPT()
    SWR = SWRCurves(100.0, 201, 0.5, 0.1, 0.01, 2, 4,
                    useFile='ks1_SWRrec.pick')
    SWR.fitPT()

    std = st.fit_steady(save_name='ks1', folder='Steady_fit')
    std.load()

    rain_rec = np.array(rain_rec)
    max_t = max(rain_rec[:, 0])

    q_curve = list()
    s_curve = list()

    bkw_q_curve = list()
    bkw_s_curve = list()

    S0 = std.ie_to_storage(rainParse2(rain_rec, 0.0, dt))
    Q0 = std.ie_to_outflow(rainParse2(rain_rec, 0.0, dt))
    max_S = S0
    max_Q = Q0

    q_curve.append([0.0, Q0])
    s_curve.append([0.0, S0])
    bkw_q_curve.append([0.0, Q0])
    bkw_s_curve.append([0.0, S0])

    t = 0.0
    S = S0
    Q = Q0
    bkw_S = S0
    bkw_Q = Q0

    pre_ie = rainParse(rain_rec, 0.0)
    k, c_curve, t = initialk(rain_rec, dt, std, Rise, Fall)
    steay = True

    while t < max_t:
        t = t + dt

        # Current rainfall
        c_ie = rainParse2(rain_rec, t, dt)

        # Compare current rainfall and last rainfall
        if c_ie != pre_ie:
            print('rainfall change to %f' % c_ie*3600*1000)
            max_S = std.ie_to_storage(c_ie)
            max_Q = std.ie_to_outflow(c_ie)

            if steady:
                k = c_ie/pre_ie
            else:



        dS1 = h1(length, c_ie, k, S/max_S, dt, std, Fall, Rise)
        dS2 = h2(length, c_ie, k, S/max_S, dS1, dt, std, Fall, Rise)
        dS3 = h3(length, c_ie, k, S/max_S, dS2, dt, std, Fall, Rise)

        dS = 0.25*dS1 + 0.75*dS3

        if k > 1.0:
            Q = Rise.s_to_q(k, (S + dS)/max_S)*max_Q
        elif k < 1.0:
            Q = Fall.s_to_q(k, (S + dS)/max_S)*max_Q

        #  print k, (S + dS)/max_S
        S = S + length*c_ie*dt - Q*dt
        s_curve.append([t, S])
        q_curve.append([t, Q])

        bkw_h1 = r_h1(length, c_ie, bkw_S, dt, std)
        bkw_h2 = r_h2(length, c_ie, bkw_S, dt, bkw_h1, std)
        bkw_h3 = r_h3(length, c_ie, bkw_S, dt, bkw_h2, std)

        dS = 0.25*bkw_h1 + 0.75*bkw_h3
        bkw_Q = std.ie_to_outflow(std.s_to_ie(bkw_S+dS))
        bkw_S = bkw_S + length*c_ie*dt - bkw_Q*dt

        bkw_s_curve.append([t, bkw_S])
        bkw_q_curve.append([t, bkw_Q])

        if c_ie != pre_ie:
            pre_ie = c_ie

    q_curve = np.array(q_curve)
    s_curve = np.array(s_curve)
    bkw_s_curve = np.array(bkw_s_curve)
    bkw_q_curve = np.array(bkw_q_curve)

    # Add reference line of steady-state S-Qs.
    ie_list = np.arange(0.1, 160.1, 1.0)
    std_s = std.ie_to_storage(ie_list/3600*10**-3)
    std_o = std.ie_to_outflow(ie_list/3600*10**-3)

    # Plot rainfall and outflow
    t = np.arange(0.0, max(rain_rec[:, 0]) + 0.1, 60.0)
    rain_bar = np.interp(t, rain_rec[:, 0], rain_rec[:, 1])

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 10))
    axes[0].bar(t/60, rain_bar, facecolor='k', width=1.3, edgecolor='b',
                linewidth=0.1)
    axes[0].set_ylim([0.0, max(rain_rec[:, 1])])
    axes[0].set_ylabel('Rain (mm/hr)')
    axes[0].set_xlabel('Time (min)')
    axes[0].set_ylim([220.0, 0.0])
    axes[0].set_title('Rain - Test Case')

    axes[1].plot(q_curve[:, 0]/60, q_curve[:, 1], color='dodgerblue',
                 label='mBKW', lw=2.0)
    axes[1].plot(bkw_q_curve[2:-1, 0]/60, bkw_q_curve[2:-1, 1], color='r',
                 label='BKW', lw=1.5)
    axes[1].set_title('Outflow - Test Case')
    axes[1].legend(loc='upper right')
    axes[1].set_xlabel('time (s)')
    axes[1].set_ylabel('Outflow ($m^2/s$)')
    plt.savefig('ks2_ShowCase_Outflow_Compare.png')
    plt.close()

    plt.figure(figsize=(10, 5))
    plt.plot(s_curve[:, 0], s_curve[:, 1], color='dodgerblue', label='mBKW',
             lw=2.0)
    plt.plot(bkw_s_curve[2:-1, 0], bkw_s_curve[2:-1, 1], color='r',
             label='BKW', lw=1.5)
    plt.legend(loc='upper right')
    plt.title('Storage - Test Case')
    plt.xlabel('time (s)')
    plt.ylabel('Storage ($m^3$)')
    plt.savefig('ks2_ShowCase_Storage_Compare.png')
    plt.close()

    plt.figure
    plt.plot(std_o, std_s, '-r', linewidth=2.0)
    plt.plot(q_curve[:, 1], s_curve[:, 1], '--b', label='mBKW', lw=2.5)
    plt.plot(bkw_q_curve[2:-1, 1], bkw_s_curve[2:-1, 1], '-k', label='BKW')
    plt.legend(loc=2)
    plt.title('S-Q Curve')
    plt.xlabel('Outflow ($m^3/s$)')
    plt.ylabel('Storage ($m^3$)')
    plt.savefig('ks2_storage-outflow_ShowCase.png')
    plt.close()


def calc_testk(rain_rec, length):
    dt = 5.0
    Fall = sh.FallingLimb()
    Rise = sh.RisingLimb()
    std = pickle.load(open('Steady_fit/SQ_Steady.pick', 'r'))

    rain_rec = np.array(rain_rec)
    max_t = max(rain_rec[:, 0])

    q_curve = list()
    s_curve = list()

    q2_curve = list()
    s2_curve = list()

    S0 = std.ie_to_storage(rainParse2(rain_rec, 0.0))
    Q0 = std.ie_to_outflow(rainParse2(rain_rec, 0.0))
    max_S = S0
    max_Q = Q0

    q_curve.append([0.0, Q0])
    s_curve.append([0.0, S0])
    q2_curve.append([0.0, Q0])
    s2_curve.append([0.0, S0])

    t = 0.0
    S = S0
    Q = Q0
    S2 = S0
    Q2 = Q0
    k_rec = list()

    k = rainParse(rain_rec, dt)/rainParse(rain_rec, 0.0)  # Initial k setting
    fixed_k = k
    while t < max_t:
        t = t + dt

        c_ie = rainParse2(rain_rec, t)
        pre_ie = rainParse2(rain_rec, t-dt)
        if c_ie != pre_ie:
            print(c_ie)
            max_S = std.ie_to_storage(c_ie)
            max_Q = std.ie_to_outflow(c_ie)

        print(t, Q)

        if t % 60.0 == 0 or c_ie != pre_ie:
            if S < std.ie_to_storage(c_ie):
                k, t_star = Rise.SQtok(
                    S/max_S, Q/max_Q, (std.ie_to_storage(c_ie)/max_S)**(5.0/3))
                #  k, t_star = Rise.SQtok(S/max_S, Q/max_Q, 5.0)
            elif S > std.ie_to_storage(c_ie):
                k, t_star = Fall.SQtok(S/max_S, Q/max_Q,
                                       (std.ie_to_storage(c_ie)/max_S)**(5.0/3)
                                       - 0.01)
            print(k)
            k_rec.append(k)

        dS1 = h1(length, c_ie, k, S/max_S, dt, std, Fall, Rise)
        dS2 = h2(length, c_ie, k, S/max_S, dS1, dt, std, Fall, Rise)
        dS3 = h3(length, c_ie, k, S/max_S, dS2, dt, std, Fall, Rise)

        dS = 0.25*dS1 + 0.75*dS3

        dS1_2 = h1(length, c_ie, fixed_k, S2/max_S, dt, std, Fall, Rise)
        dS2_2 = h2(length, c_ie, fixed_k, S2/max_S, dS1_2, dt, std, Fall, Rise)
        dS3_2 = h3(length, c_ie, fixed_k, S2/max_S, dS2_2, dt, std, Fall, Rise)

        dS_2 = 0.25*dS1_2 + 0.75*dS3_2

        if k > 1.0:
            Q = Rise.s_to_q(k, (S + dS)/max_S)*max_Q
        elif k < 1.0:
            Q = Fall.s_to_q(k, (S + dS)/max_S)*max_Q

        if fixed_k > 1.0:
            Q2 = Rise.s_to_q(fixed_k, (S2 + dS_2)/max_S)*max_Q
        elif fixed_k < 1.0:
            Q2 = Fall.s_to_q(fixed_k, (S2 + dS_2)/max_S)*max_Q

        #  print k, (S + dS)/max_S
        S = S + length*c_ie*dt - Q*dt
        S2 = S2 + length*c_ie*dt - Q2*dt
        s_curve.append([t, S])
        q_curve.append([t, Q])
        s2_curve.append([t, S2])
        q2_curve.append([t, Q2])

    q_curve = np.array(q_curve)
    s_curve = np.array(s_curve)
    q2_curve = np.array(q2_curve)
    s2_curve = np.array(s2_curve)

    q_curve[:, 0] = q_curve[:, 0] - 60.0
    s_curve[:, 0] = s_curve[:, 0] - 60.0
    q2_curve[:, 0] = q2_curve[:, 0] - 60.0
    s2_curve[:, 0] = s2_curve[:, 0] - 60.0

    ab = an.trans_ana(54.3, 3.11, 100.0, 201, 1.0, 0.1, 0.01)
    ab.run()

    standard_k = np.ones(len(k_rec))
    standard_k = standard_k*18.31

    ie_list = np.arange(0.1, 160.1, 1.0)
    std_s = std.ie_to_storage(ie_list/3600*10**-3)
    std_o = std.ie_to_outflow(ie_list/3600*10**-3)

    plt.figure
    plt.plot(q_curve[:, 0], q_curve[:, 1], '-b', label='k-AutoSearch')
    plt.plot(q2_curve[:, 0], q2_curve[:, 1], '-k', label='k-fixed')
    ab.plot_QT(LT='--r', plt_label='Analytic')
    plt.title('Outflow - Test k 54.3 mm/hr to 3.11 mm/hr')
    plt.legend(loc=4)
    plt.xlabel('time (s)')
    plt.ylabel('Outflow ($m^3/s$)')
    plt.xlim([0, t])
    plt.savefig('modified-BKW_testk3.png')
    plt.close()

    plt.figure
    plt.plot(s_curve[:, 0], s_curve[:, 1], '-b', label='k-AutoSearch')
    plt.plot(s2_curve[:, 0], s2_curve[:, 1], '-k', label='k-fixed')
    ab.plot_ST(LT='--r', plt_label='Analytic')
    plt.legend(loc=4)
    plt.title('Storage - Test k 54.3 mm/hr to 3.11 mm/hr')
    plt.xlabel('time (s)')
    plt.ylabel('Storage ($m^3$)')
    plt.xlim([0, t])
    plt.savefig('modified-BKW_storage_testk3.png')
    plt.close()

    plt.figure
    plt.plot(std_o, std_s, '-r', linewidth=2.0)
    plt.plot(q_curve[:, 1], s_curve[:, 1], '-b', label='k-AutoSearch')
    plt.plot(q2_curve[:, 1], s2_curve[:, 1], '-m', label='k-fixed')
    plt.legend(loc=2)
    plt.title('S-Q - Test Case k 54.3 mm/hr to 3.11 m/hr')
    plt.xlabel('Outflow ($m^3/s$)')
    plt.ylabel('Storage ($m^3$)')
    plt.savefig('storage-outflow_testk3.png')
    plt.close()


def run():
    rain_rec = [[0.0, 5.0], [876.0, 5.0], [876.01, 20.0], [1379.0, 20.0],
                [1379.01, 100.0], [1644.0, 100.0], [1644.01, 200.0],
                [1845.0, 200.0], [1845.01, 160], [2065.0, 160], [2065.01, 80],
                [2355.0, 80], [2355.01, 5], [3231.0, 5.0], [3231.01, 15.0],
                [3797.0, 15.0], [3797.01, 2.0], [5062.0, 2.0],
                [5062.01, 0.0], [7200.0, 0.0]]
    rain_rec = np.array(rain_rec)

    t = np.arange(0.0, max(rain_rec[:, 0]) + 0.1, 2.0)
    rain_bar = np.interp(t, rain_rec[:, 0], rain_rec[:, 1])

    plt.figure
    plt.bar(t, rain_bar, facecolor='k', width=1.3, edgecolor='b',
            linewidth=0.1)
    plt.ylim([0.0, max(rain_rec[:, 1])])
    plt.ylabel('Rain (mm/hr)')
    plt.xlabel('Time (min)')
    plt.title('Rain')

    plt.savefig('TestCase_Rain 5.png')
    plt.close()

    # rain_rec[:, 0] = rain_rec[:, 0]*60
    calc(rain_rec, 10.0)

    """
rain_rec1 = [[0.0, 10.0], [0.01, 73.62], [30.0, 73.62]]
rain_rec2 = [[0.0, 50.0], [0.01, 10.0], [100.0, 10.0]]

    rain_rec = [[0.0, 10.0], [0.01, 20.0], [10.0, 10.0], [15.0, 60.0],
                [30.0, 60.0], [35.0, 60.0], [35.01, 25.0], [50.0, 25.0]]

    rain_rec = [[0.0, 132.0], [0.01, 85.12], [10.0, 65.12], [10.01, 65.12],
                [30.0, 53.41], [30.01, 53.41], [50.0, 86.17], [60.0, 195.21],
                [70.0, 195.21], [90.0, 35.10]]
    """


def run_testk():
    rain_rec = [[0.0, 54.3], [0.01, 3.11], [90.0, 3.11]]
    rain_rec = np.array(rain_rec)

    t = np.arange(0.0, max(rain_rec[:, 0]) + 0.1, 2.0)
    rain_bar = np.interp(t, rain_rec[:, 0], rain_rec[:, 1])

    plt.figure
    plt.bar(t, rain_bar, facecolor='k', width=1.3, edgecolor='b', linewidth=0.1)
    plt.ylim([0.0, 100.0])
    plt.ylabel('Rain (mm/hr)')
    plt.xlabel('Time (min)')
    plt.title('Rain')

    plt.savefig('TestCase_Rain_test_K.png')
    plt.close()

    rain_rec[:, 0] = rain_rec[:, 0]*60
    calc_testk(rain_rec, 100.0)
