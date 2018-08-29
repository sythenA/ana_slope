
import sys
import steady_fit as st
import bezeir_fit as bf
import pickle
import numpy as np
from math import fabs, hypot
import matplotlib.pyplot as plt


class RisingLimb:
    """
    To generate curve by ie1/ie0 or search by normalized storage, outflow,
    and ie1/ie0, with normalized storage and k, normalized outflow can be
    specified.

    To build RisingLimb class:
        folder: folder of preserved fitting results.
    """
    def __init__(self, folder='Steady_fit', rec_name='ks1_Rise'):

        aa = st.RisingLimb_fit(rec_name=rec_name, save_name='fit_1',
                               folder=folder)
        aa.load()
        """
        file_name = folder + '/' + 'fit_1_RiseLimb.pick'
        aa = pickle.load(open(file_name, 'r'))"""

        self.ref = aa

    def load_sq(self, k):
        ref = self.ref
        min_k = min(ref.ks)
        max_k = max(ref.ks)

        if k >= min_k and k <= max_k:
            val = ref.SQ_curve(k)
        elif k > max_k:
            val = ref.SQ_curve(max_k)
        else:
            val = ref.SQ_curve(min_k)

        return val

    def mat_flip(self, mat):
        rev_mat = np.fliplr([mat])[0]

        return rev_mat

    def s_to_q(self, k, normS):
        nQ = self.ref.NS_to_NQ(normS, k)
        return nQ

    def test_k(self, k0, k1, normS, normQ):
        Q_star1 = self.s_to_q(1.0/k0, normS)
        Q_star2 = self.s_to_q(1.0/k1, normS)
        f1 = Q_star1 - normQ
        f2 = Q_star2 - normQ

        return f1*f2

    def get_k(self, normS, normQ):
        err = 10*6
        k0 = 1.0
        k1 = 200.0
        while err > 10**-5:
            t1 = self.test_k(k0, k1, normS, normQ)
            if t1 < 0:
                t2 = self.test_k(k0, 0.5*(k0+k1), normS, normQ)
                t3 = self.test_k(0.5*(k0+k1), k1, normS, normQ)
                if t2 < 0 and t3 > 0:
                    k0 = k0
                    k1 = 0.5*(k0+k1)
                elif t3 < 0 and t2 > 0:
                    k0 = 0.5*(k0+k1)
                    k1 = k1
                elif t2 > 0 and t3 > 0:
                    raise ValueError('No solution in k')
                elif t2 == 0:
                    return 0.5*(k0 + k1)
                elif t3 == 0:
                    return k1
            else:
                raise ValueError('No solution in k')
            err = fabs(k0 - k1)

        return k1

    def derivate(self, inv_k, t_star):
        seq = self.seq
        max_k = self.max_k

        k1 = 1.0/(1.0/inv_k+0.05)
        if k1 > 1.0 or k1 < 0.0:
            k1 = 0.6
        k2 = 1.0/(1.0/inv_k-0.05)
        if k2 > 1.0 or k2 < 0.0:
            k2 = 0.5
        for i in range(1, len(max_k)):
            if 1.0/k1 > max_k[i-1] and 1.0/k1 <= max_k[i]:
                Dos1 = seq[i-1].outflow_D(k1)
                print k1
            if 1.0/k2 > max_k[i-1] and 1.0/k2 <= max_k[i]:
                Dos2 = seq[i-1].outflow_D(k2)
                print k2
        print max_k

        dev_Ds = (Dos2 - Dos1)/(k2 - k1)
        fun_dev_D = np.poly1d(dev_Ds)

        return fun_dev_D(t_star)

    def SQtok(self, normS, normQ, init_k):
        ref = self.ref
        err = 1.0
        k0 = init_k
        counter = 0

        while err > 10**-3:
            kp = k0 + 0.001
            km = k0 - 0.001
            Q_star = ref.NS_to_NQ(normS, k0)
            Q_star1 = ref.NS_to_NQ(normS, kp)
            Q_star2 = ref.NS_to_NQ(normS, km)

            #  deg = len(ref.ref_points(k0)) - 1
            t_star = ref.NS_to_t(normS, k0)
            print t_star
            f = fabs(normQ - Q_star)
            """
            print ref.derivate(t_star)
            df = bf.gen(t_star, deg, ref.derivate(t_star))
            """
            df = (fabs(Q_star1 - normQ) - fabs(Q_star2 - normQ))/0.002

            """
            df = hypot(df[0], df[1])
            """
            print df
            if df == 0:
                break

            #  k1 = k0 - f/hypot(df[0], df[1])
            k1 = k0 - f/df
            if k1 < min(ref.ks) or k1 > max(ref.ks):
                k1 = 5.0
            """
            elif normQ < 1.0/k1:
                k1 = k1 + 0.5
                """
            err = fabs(k1-k0)
            print err
            k0 = k1
            counter += 1
            if counter > 50:
                break

        return k0, self.s_to_q(k0, normS)


class FallingLimb:
    def __init__(self, folder='Steady_fit', rec_name='ks1_Fall'):

        aa = st.FallingLimb_fit(rec_name=rec_name, save_name='fit_1',
                                folder=folder)
        aa.load()
        """
        file_name = folder + '/' + 'fit_1' '_FallLimb.pick'
        print file_name
        aa = pickle.load(open(file_name, 'r'))"""

        self.ref = aa
        self.max_k = max(aa.ks)
        self.min_k = min(aa.ks)

    def load_sq(self, k):
        max_k = self.max_k
        min_k = self.min_k

        if k <= max_k and k >= self.min_k:
            val = self.ref.SQ_curve(1./k)
        elif k > max_k:
            val = self.ref.SQ_curve(1./max_k)
        elif k < min_k:
            val = self.ref.SQ_curve(1./min_k)

        return val

    def mat_flip(self, mat):
        rev_mat = np.fliplr([mat])[0]

        return rev_mat

    def s_to_q(self, inv_k, normS):

        nQ = self.ref.NS_to_NQ(normS, inv_k)
        return nQ

    def test_k(self, k0, k1, normS, normQ):
        Q_star1 = self.s_to_q(k0, normS)
        Q_star2 = self.s_to_q(k1, normS)
        f1 = Q_star1 - normQ
        f2 = Q_star2 - normQ

        return f1*f2

    def get_k(self, min_k, max_k, normS, normQ):
        err = 10*6
        k0 = min_k
        k1 = 200.0
        while err > 10**-3:
            t1 = self.test_k(k0, k1, normS, normQ)
            if t1 < 0:
                t2 = self.test_k(k0, 0.5*(k0+k1), normS, normQ)
                t3 = self.test_k(0.5*(k0+k1), k1, normS, normQ)
                if t2 < 0 and t3 > 0:
                    k0 = k0
                    k1 = 0.5*(k0+k1)
                elif t3 < 0 and t2 > 0:
                    k0 = 0.5*(k0+k1)
                    k1 = k1
                elif t2 > 0 and t3 > 0:
                    raise ValueError('No solution in k')
                elif t2 == 0:
                    return 0.5*(k0 + k1)
                elif t3 == 0:
                    return k1
            else:
                raise ValueError('No solution in k')
            err = fabs(self.s_to_q(k0, normS) - normQ)

        return k0

    def SQtok(self, normS, normQ, init_k):
        ref = self.ref
        err = 1.0
        k0 = 1./init_k
        counter = 0

        while err > 10**-3:
            kp = k0 + 0.01
            km = k0 - 0.01
            Q_star = ref.NS_to_NQ(normS, 1./k0)
            Q_star1 = ref.NS_to_NQ(normS, 1./kp)
            Q_star2 = ref.NS_to_NQ(normS, 1./km)

            """
            deg = len(ref.ref_points(1./k0)) - 1
            t_star = ref.NS_to_t(normS, 1./k0)
            """
            f = fabs(normQ - Q_star)
            #  df = bf.gen(t_star, deg, ref.derivate(1./k0))
            df = (fabs(normQ - Q_star2) - fabs(normQ - Q_star1))/(1./km - 1./kp)

            """
            if hypot(df[0], df[1]) == 0:
                break
            """

            #  k1 = k0 - f/hypot(df[0], df[1])
            k1 = k0 - f/-df
            if k1 < min(ref.ks) or k1 > max(ref.ks):
                k1 = 5.0
            """
            elif normQ < 1.0/k1:
                k1 = k1 + 0.5
                """
            err = fabs(k1 - k0)
            k0 = k1
            counter += 1
            if counter > 50:
                break
            print err

        return 1./k0, self.s_to_q(k0, normS)


def run_Steady_gen(length, dt, n, slope):
    aa = st.fit_steady(save_name='SQ')
    aa.param(length, n, dt, slope)
    aa.gen()
    aa.SaveResult()


def run_Rising_gen(length, dt, n, slope, gen_dat=False):
    ks1 = [1.01, 1.05, 1.1, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0,
           7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0,
           90.0, 100.0]

    aa = st.RisingLimb_gen(length, dt, n, slope, ks1, rec_name='ks1_Rise',
                           gen=gen_dat)
    aa.run()
    plt.figure
    aa.plot_norm_SQ()
    plt.show()


def run_RisingLimb_fit():
    ks1 = [1.01, 1.05, 1.1, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0,
           7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0,
           90.0, 100.0]
    ab = st.RisingLimb_fit(ks1, rec_name='ks1_Rise', save_name='fit_1')
    ab.run()
    ab.SaveResult()

    val1 = ab.SQ_curve(42.0)
    val2 = ab.SQ_curve(7.4)
    plt.figure
    plt.plot(val1[:, 1], val1[:, 0], '-b')
    plt.plot(val2[:, 1], val2[:, 0], '-r')
    plt.show()


def run_Falling_gen(length, dt, n, slope, gen_dat=False):
    ks1 = [1.01, 1.05, 1.1, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0,
           7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0,
           90.0, 100.0]

    aa = st.FallingLimb_gen(length, dt, n, slope, ks1, rec_name='ks1_Fall',
                            gen=gen_dat)
    aa.run()
    plt.figure
    aa.plot_norm_SQ()
    plt.show()


def run_FallingLimb_fit():
    ks1 = [1.01, 1.05, 1.1, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0,
           7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0,
           90.0, 100.0]

    ab = st.FallingLimb_fit(ks1, rec_name='ks1_Fall', save_name='fit_1')
    ab.run()
    ab.SaveResult()
    val1 = ab.SQ_curve(0.32)
    val2 = ab.SQ_curve(0.015)
    plt.figure
    plt.plot(val1[:, 1], val1[:, 0], '-b')
    plt.plot(val2[:, 1], val2[:, 0], '-r')
    plt.show()
