
import ana_new as ana
import pickle
import matplotlib.pyplot as plt


def run():
    aa = pickle.load(open('Steady_fit/fit_1_DownToUp.pick', 'r'))
    t, out = aa.out_curve(4.3)

    ab = ana.trans_ana(10.0, 43.0, 100.0, 1001, 1.0, 0.1, 0.01)
    ab.run2()
    out = out*ab.max_q
    t = t*ab.max_t

    plt.figure
    plt.plot(t, out, '-b', label='Interpolation')
    ab.plot_QT(plt_label='Analytical', LT='--r')
    plt.xlabel('time (s)')
    plt.ylabel('flowrate ($m^3/s$)')
    plt.title('10 mm/hr to 43 mm/hr, k=4.3')
    plt.legend(loc=4)
    plt.savefig('10to43.png')
    plt.close()

    ac = pickle.load(open('Steady_fit/fit_1_UpToDown.pick', 'r'))
    t, out = ac.out_curve(0.34)

    ad = ana.trans_ana(100.0, 34.0, 100.0, 1001, 1.0, 0.1, 0.01)
    ad.run2()
    out = out*ad.min_q
    t = t*ad.max_t

    plt.figure
    plt.plot(t, out, '-b', label='Interpolation')
    ad.plot_QT(plt_label='Analytical', LT='--r')
    plt.xlabel('time (s)')
    plt.ylabel('flowrate ($m^3/s$)')
    plt.title('100 mm/hr to 34 mm/hr, k=0.34')
    plt.legend(loc=4)
    plt.savefig('100to34.png')
    plt.close()

    plt.figure
    aa.plot_C(0)
    plt.title('$C_1$')
    plt.legend(loc=3)
    plt.savefig('C1.png')
    plt.close()

    plt.figure
    aa.plot_C(1)
    plt.title('$C_2$')
    plt.legend(loc=4)
    plt.savefig('C2.png')
    plt.close()

    plt.figure
    aa.plot_C(2)
    plt.title('$C_3$')
    plt.legend(loc=1)
    plt.savefig('C3.png')
    plt.close()

    plt.figure
    aa.plot_C(3)
    plt.title('$C_4$')
    plt.legend(loc=1)
    plt.savefig('C4.png')
    plt.close()

    plt.figure
    aa.plot_C(4)
    plt.title('$C_5$')
    plt.legend(loc=1)
    plt.savefig('C5.png')
    plt.close()
