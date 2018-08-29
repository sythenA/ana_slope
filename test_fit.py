
import pickle
import bezeir_fit as bf
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from ana_new import progression, Recession


def run():
    REC = pickle.load(open('Steady_fit/ks1_Rise_CurveRec.flow'))
    param = pickle.load(open('Steady_fit/fit_1_RiseLimb.pick'))
    key_str = 'k=10.0'
    sto = REC['S'][key_str]
    out = REC['Q'][key_str]

    sto = sto[0:len(out), :]

    out[:, 1] = (out[:, 1] - min(out[:, 1]))/(1.0 - min(out[:, 1]))
    sto[:, 1] = (sto[:, 1] - min(sto[:, 1]))/(1.0 - min(sto[:, 1]))
    lst = np.arange(1, len(out), len(out)/10)
    lst = lst.tolist()
    lst.pop(-1)

    plt.figure
    plt.plot(out[:, 1], sto[0:len(out), 1], color='dodgerblue', lw=3.0,
             label='Analytic')

    p = param['param']
    p = p['k=10.0']

    new_t = np.arange(0, 1.001, 0.01)
    new_t = new_t.tolist()
    new_t.append(1.0)
    new_t = np.array(new_t)
    y2 = bf.gen(new_t, 6, p)

    plt.plot(y2[:, 1], y2[:, 0], '--r', label='Bezeir-fit')
    counter = 0
    for point in p:
        plt.plot(point[0], point[1], marker='s', color='k')
        plt.text(point[0]+0.02, point[1]+0.01, 'P'+str(counter))
        counter += 1
    plt.legend(loc=4)
    plt.xlabel('$Q_{{}*{}}$')
    plt.ylabel('$S_{{}*{}}$')
    plt.title('k=10 Analytical S-Q fit (n=6)')
    plt.savefig('s-q fit test.png')
    plt.close()


# To exibit the multiple correspndence of storage to outflow
def run2():
    a1 = progression(10.0, 10.0, 21, 1.0, 0.1, 0.01)
    a2 = Recession(10.0, 10.0, 21, 1.0, 0.1, 0.01)

    b1 = progression(30.0, 10.0, 21, 1.0, 0.1, 0.01)
    b2 = Recession(30.0, 10.0, 21, 1.0, 0.1, 0.01)

    c1 = progression(50.0, 10.0, 21, 1.0, 0.1, 0.01)
    c2 = Recession(50.0, 10.0, 21, 1.0, 0.1, 0.01)

    d1 = progression(75.0, 10.0, 21, 1.0, 0.1, 0.01)
    d2 = Recession(75.0, 10.0, 21, 1.0, 0.1, 0.01)

    e1 = progression(100.0, 10.0, 21, 1.0, 0.1, 0.01)
    e2 = Recession(100.0, 10.0, 21, 1.0, 0.1, 0.01)

    a1.run()
    a2.run()
    b1.run()
    b2.run()
    c1.run()
    c2.run()
    d1.run()
    d2.run()
    e1.run()
    e2.run()

    steady_SQ = list()
    m = 5.0/3
    for i in range(1, 121):
        ie = float(i)/3600.0*1.0E-3
        alpha = 1/0.1*sqrt(0.01)
        q = ie*10.0
        s = m/(m+1)*(ie/alpha)**(1/m)*10.0**(1 + 1./m)
        steady_SQ.append([q, s])

    steady_SQ = np.array(steady_SQ)

    plt.figure
    plt.plot(a1.out_q[:, 1], a1.out_s[:, 1], '-b', lw=1.5)
    plt.plot(a2.out_q[:, 1], a2.out_s[:, 1], '-b', lw=1.5)
    plt.plot(b1.out_q[:, 1], b1.out_s[:, 1], '-b', lw=1.5)
    plt.plot(b2.out_q[:, 1], b2.out_s[:, 1], '-b', lw=1.5)
    plt.plot(c1.out_q[:, 1], c1.out_s[:, 1], '-b', lw=1.5)
    plt.plot(c2.out_q[:, 1], c2.out_s[:, 1], '-b', lw=1.5)
    plt.plot(d1.out_q[:, 1], d1.out_s[:, 1], '-b', lw=1.5)
    plt.plot(d2.out_q[:, 1], d2.out_s[:, 1], '-b', lw=1.5)
    plt.plot(e1.out_q[:, 1], e1.out_s[:, 1], '-b', lw=1.5)
    plt.plot(e2.out_q[:, 1], e2.out_s[:, 1], '-b', lw=1.5)

    plt.plot(steady_SQ[:, 0], steady_SQ[:, 1], '--r', lw=2.0)

    plt.text(steady_SQ[10, 0], steady_SQ[10, 1]+0.001, '10.0')
    plt.text(steady_SQ[30, 0], steady_SQ[30, 1]+0.001, '30.0')
    plt.text(steady_SQ[50, 0], steady_SQ[50, 1]+0.001, '50.0')
    plt.text(steady_SQ[75, 0], steady_SQ[75, 1]+0.001, '75.0')
    plt.text(steady_SQ[100, 0], steady_SQ[100, 1]+0.001, '100.0')
    plt.xlabel('q ($m^2$/s)')
    plt.ylabel('s ($m^2$)')
    plt.savefig('Hydrolodical Loop.png')
    plt.close()


# To show the "looping" in hydrological phenomena
def run3():
    a1 = progression(100.0, 10.0, 21, 1.0, 0.1, 0.01)
    a2 = Recession(100.0, 10.0, 21, 1.0, 0.1, 0.01)

    a1.run()
    a2.run()

    steady_SQ = list()
    m = 5.0/3
    for i in range(1, 121):
        ie = float(i)/3600.0*1.0E-3
        alpha = 1/0.1*sqrt(0.01)
        q = ie*10.0
        s = m/(m+1)*(ie/alpha)**(1/m)*10.0**(1 + 1./m)
        steady_SQ.append([q, s])

    steady_SQ = np.array(steady_SQ)

    plt.figure
    plt.plot(a1.out_q[:, 1], a1.out_s[:, 1], '-', color='dodgerblue', lw=2.5)
    plt.plot(a2.out_q[:, 1], a2.out_s[:, 1], '-', color='dodgerblue', lw=2.5)

    plt.plot(steady_SQ[:, 0], steady_SQ[:, 1], '--k', lw=2.0)

    plt.xlabel('q ($m^2$/s)')
    plt.ylabel('s ($m^2$)')
    plt.savefig('one Loop.png')
    plt.close()


#  To draw the characteristic curve movement on slope
def run4():
    m = 5.0/3
    alpha = 1./0.1*sqrt(0.01)  # 1% slope, n=0.1
    ie = 100.0/3600.0*1.0E-3  # effective rainfall 100.0 mm/hr

    trace = np.zeros([10, 267])
    depth = np.zeros([267, 101])
    x = np.arange(0, 10.1, 0.1)

    ref = np.zeros(267)
    trace[:, 0] = np.arange(0, 10, 1.)
    for i in range(0, 10):
        for j in range(1, 267):
            trace[i, j] = trace[i, 0] + alpha*ie**(m-1)*j**m
    for k in range(0, 267):
        ref[k] = alpha*ie**(m-1)*k**m

    for i in range(0, 267):
        for j in range(0, 101):
            c_front = alpha*ie**(m-1)*i**m
            if x[j] < c_front:
                depth[i, j] = (ie*x[j]/alpha)**(1./m)
            else:
                depth[i, j] = ie*i

    for i in range(0, 267):
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))
        axes[0].plot(trace[0, 0:i], ref[0:i], '-r', lw=2.0)
        for j in range(1, 10):
            axes[0].plot(trace[j, 0:i], ref[0:i], '-b')
        axes[0].set_xlim([0, 10.0])
        axes[0].set_ylim([0, 12.0])
        axes[0].set_xlabel('Wave Position (m)')
        axes[0].set_ylabel('Distance Traveled (m)')
        axes[0].set_title('Characteristic Curve')

        axes[1].plot(x, depth[i, :], '-r', lw=1.5)
        axes[1].set_xlabel('Length (m)')
        axes[1].set_ylabel('Depth (m)')
        axes[1].set_ylim([0.0, 0.01])
        axes[1].set_xlim([0.0, 10.0])
        axes[1].set_title('Flow Depth')

        fig.suptitle('Time ' + str(i) + ' s')
        plt.savefig('car_curve/' + str(i).zfill(3) + '.png')
        plt.close()
