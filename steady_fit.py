
import os
import os.path
import numpy as np
import ana_new as aw
import matplotlib.pyplot as plt
from math import sqrt
import pickle
import bezeir_fit as bf
import polyfit as pf


def mat_flip(mat):
    rev_mat = np.fliplr([mat])[0]

    return rev_mat


class FallingLimb_gen:
    def __init__(self, length, dt, n, slope, ks, rec_name='',
                 folder='ana_slope/Steady_fit', gen=False):
        self.length = length
        self.dt = dt
        self.n = n
        self.slope = slope
        self.alpha = 1.0/n*sqrt(slope)
        self.ks = ks
        self.rec_name = rec_name
        self.folder = folder
        self.gen = gen

    def load(self):
        rec_file = self.folder + '/' + self.rec_name + '_CurveRec.flow'
        REC = pickle.load(open(rec_file, 'r'))
        self.REC = REC

    def fit(self):
        ks = self.ks
        param = dict()
        for k in ks:
            key_str = 'k=' + str(float(k))
            out_s = self.REC['S'][key_str]
            out_q = self.REC['Q'][key_str]

            """
            lst = np.arange(len(out_q)*0.01, len(out_q)*0.95, len(out_q)/15.)
            for j in range(0, len(lst)):
                lst[j] = int(lst[j])
            """
            sq = np.zeros([16, 2])
            sq[0:16, 0] = np.arange(0, 1.01, 1./15)
            sq[-1, 0] = 1.0
            sq[:, 1] = np.interp(sq[:, 0], mat_flip(out_s[0:len(out_q), 1]),
                                 mat_flip(out_q)[:, 1])
            sq[:, 0] = mat_flip(sq[:, 0])
            sq[:, 1] = mat_flip(sq[:, 1])

            """
            sq = list()
            sq.append([1.0, 1.0])
            for j in range(0, len(lst)):
                sq.append([out_s[lst[j], 1], out_q[lst[j], 1]])
            sq = np.array(sq)"""

            p, t = bf.fit(sq, 6, 10**-6)
            p[0, :] = [1.0, 1.0]
            p[-1, :] = [0.0, 0.0]

            param.update({'k='+str(k): p})

        self.param = param
        self.param.update({'k': ks})

        if self.rec_name:
            pickle.dump(param, open(self.folder + '/' +
                                    self.rec_name + '_Fall_param.flow', 'w'))

    def run(self):
        if self.gen:
            self.gen_rec()
            self.fit()
        else:
            self.load()
            self.fit()

    def gen_rec(self):
        dt = self.dt
        n = self.n
        slope = self.slope
        length = self.length
        ks = self.ks

        out_rec = dict()
        sto_rec = dict()

        self.Steady_SQ()
        for k in range(0, len(self.ks)):
            ie1 = 10.0
            aa = aw.trans_ana(ie1*self.ks[k], ie1, length, 201, dt, n, slope)
            aa.run()

            out_s = aa.out_s
            out_q = aa.out_q
            out_s[:, 0] = out_s[:, 0]/aa.max_t
            out_s[:, 1] = (out_s[:, 1] - aa.min_s)/(aa.max_s - aa.min_s)
            out_q[:, 0] = out_q[:, 0]/aa.max_t
            out_q[:, 1] = (out_q[:, 1] - aa.min_q)/(aa.max_q - aa.min_q)

            out_rec.update({'k='+str(ks[k]): out_q})
            sto_rec.update({'k='+str(ks[k]): out_s})

        REC = {'S': sto_rec, 'Q': out_rec}
        self.REC = REC

        if self.rec_name:
            pickle.dump(REC, open(self.folder + '/' +
                                  self.rec_name + '_CurveRec.flow', 'w'))

    def Steady_SQ(self):
        length = self.length
        alpha = self.alpha
        m = 5.0/3
        max_k = max(self.ks)

        sq = list()
        for ie in range(1, int(5.0*max_k)):
            ie = float(ie)/3600*10**-3
            _q = alpha*(ie*length/alpha)**(1/m)
            _s = m/(m+1)*(ie*length**(m+1)/alpha)**(1/m)

            sq.append([_s, _q])

        self.SQ_std = np.array(sq)

    def plot_SQ(self):
        ks = self.ks
        param = self.param
        length = self.length
        alpha = self.alpha
        m = 5.0/3
        ie0 = 5.0

        dot_s = np.arange(0, 1.0, 0.01)
        for k in range(0, len(ks)):
            ie1 = ie0/ks[k]
            min_q = ie1*length
            max_q = ie0*length
            min_s = m/(m+1)*(ie1*length**(m+1)/alpha)*(1/m)
            max_s = m/(m+1)*(ie0*length**(m+1)/alpha)*(1/m)

            _str = 'k=' + str(ks[k])
            p = param[_str]
            sq = bf.gen(dot_s, p)
            sq[:, 1] = sq[:, 1]*(max_q - min_q) + min_q
            sq[:, 0] = sq[:, 0]*(max_s - min_s) + min_s

            plt.plot(sq[:, 1], sq[:, 0])

    def plot_norm_SQ(self):
        ks = self.ks
        param = self.param
        REC = self.REC

        dot_s = np.arange(0, 0.901, 0.01)
        dot_s = dot_s.tolist()
        for k in range(0, len(ks)):
            _str = 'k=' + str(ks[k])
            out_s = REC['S'][_str][:, 1]
            out_q = REC['Q'][_str][:, 1]
            p = param[_str]
            sq = bf.gen(dot_s, 6, p)

            sq = sq.tolist()
            sq.append([0.0, 0.0])
            sq = np.array(sq)
            """
            q_fun = np.poly1d(p)
            dot_q = q_fun(dot_s)
            dot_q[0] = 0.0
            dot_q[-1] = 1.0
            """
            plt.plot(sq[:, 1], sq[:, 0], 'r')
            #  plt.plot(dot_s, dot_q, '-r')
            plt.plot(out_q, out_s[0: len(out_q)], '--b')


class FallingLimb_fit:
    """Load parameters of falling limb s-q curve, then build interpolation
    function of curve."""
    def __init__(self, rec_name='', save_name='',
                 folder='ana_slope/Steady_fit', **kwargs):

        param = pickle.load(
            open(folder + '/' + rec_name + '_Fall_param.flow', 'r'))

        self.param = param
        self.folder = folder
        self.ks = param['k']
        self.save_name = save_name

    def run(self):
        self.p_fit()
        self.SaveResult()

    def load(self):
        save_name = self.save_name
        folder = self.folder

        tmp_dict = pickle.load(open(folder + '/' + save_name +
                               '_FallLimb.pick', 'rb'))
        self.__dict__.update(tmp_dict)

    def mat_flip(self, mat):
        rev_mat = np.fliplr([mat])[0]

        return rev_mat

    def derivate(self, inv_k):
        #  t_star must be float, bring t_star with k.
        px = self.zx
        py = self.zy

        k = 1./inv_k

        dx = np.zeros(len(px)+2)
        dy = np.zeros(len(px)+2)
        for i in range(0, len(px)):
            dx[i+1] = pf.dev(k, px[i])
            dy[i+1] = pf.dev(k, py[i])

        dPdk = np.zeros([len(dx), 2])
        dPdk[:, 0] = dx
        dPdk[:, 1] = dy
        dPdk[0, :] = [0.0, 0.0]

        return dPdk

    def plot_p(self, r):
        param = self.param

        _px = list()
        _py = list()
        for k in range(0, len(self.ks)):
            key = 'k=' + str(self.ks[k])
            p = param[key]

            _px.append([1.0/self.ks[k], p[r, 0]])
            _py.append([1.0/self.ks[k], p[r, 1]])

        _px = np.array(_px)
        _py = np.array(_py)
        plt.plot(_px[:, 0], _px[:, 1], '-r', label='x')
        plt.plot(_py[:, 0], _py[:, 1], '-k', label='y')

    def p_fit(self):
        param = self.param
        deg = len(param['k=5.0'])-1
        ks = self.ks

        _px = np.zeros([len(ks), deg])
        _py = np.zeros([len(ks), deg])
        for k in range(0, len(ks)):
            key = 'k=' + str(ks[k])
            p = param[key]

            _px[k, 0] = ks[k]
            _py[k, 0] = ks[k]
            _px[k, 1:deg] = p[1:deg, 0]
            _py[k, 1:deg] = p[1:deg, 1]

        self._px = _px
        self._py = _py

        self.zx = self.fit_px()
        self.zy = self.fit_py()

        zx_fun = list()
        zy_fun = list()
        for i in range(0, len(self.zx)):
            zx_fun.append(np.poly1d(self.zx[i]))
            zy_fun.append(np.poly1d(self.zy[i]))
        self.zx_fun = zx_fun
        self.zy_fun = zy_fun

    def fit_px(self):
        _px = self._px
        deg = len(_px[0])

        zx = list()
        for i in range(1, deg):
            zx.append(np.polyfit(_px[:, 0], _px[:, i], 12))

        return zx

    def fit_py(self):
        _py = self._py
        deg = len(_py[0])

        zy = list()
        for i in range(1, deg):
            zy.append(np.polyfit(_py[:, 0], _py[:, i], 12))

        return zy

    def plot_fit_p(self):
        zx_fun = self.zx_fun
        zy_fun = self.zy_fun
        _px = self._px
        _py = self._py

        k = np.arange(min(self.ks), max(self.ks)+0.1, 0.1)

        plt.figure
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(6, 6))

        axes[0, 0].plot(k, zx_fun[0](k), '-r')
        axes[0, 0].plot(k, zy_fun[0](k), '-k')
        axes[0, 0].plot(_px[:, 0], _px[:, 1], '--r')
        axes[0, 0].plot(_py[:, 0], _py[:, 1], '--k')
        axes[0, 0].set_title('$P_1$')

        axes[0, 1].plot(k, zx_fun[1](k), '-r')
        axes[0, 1].plot(k, zy_fun[1](k), '-k')
        axes[0, 1].plot(_px[:, 0], _px[:, 2], '--r')
        axes[0, 1].plot(_py[:, 0], _py[:, 2], '--k')
        axes[0, 1].set_title('$P_2$')

        axes[0, 2].plot(k, zx_fun[2](k), '-r')
        axes[0, 2].plot(k, zy_fun[2](k), '-k')
        axes[0, 2].plot(_px[:, 0], _px[:, 3], '--r')
        axes[0, 2].plot(_py[:, 0], _py[:, 3], '--k')
        axes[0, 2].set_title('$P_3$')

        axes[1, 0].plot(k, zx_fun[3](k), '-r')
        axes[1, 0].plot(k, zy_fun[3](k), '-k')
        axes[1, 0].plot(_px[:, 0], _px[:, 4], '--r')
        axes[1, 0].plot(_py[:, 0], _py[:, 4], '--k')
        axes[1, 0].set_title('$P_4$')

        axes[1, 1].plot(k, zx_fun[4](k), '-r')
        axes[1, 1].plot(k, zy_fun[4](k), '-k')
        axes[1, 1].plot(_px[:, 0], _px[:, 5], '--r')
        axes[1, 1].plot(_py[:, 0], _py[:, 5], '--k')
        axes[1, 1].set_title('$P_5$')

        plt.show()

    def px(self, inv_k):
        k = 1./inv_k
        zx_fun = self.zx_fun
        deg = len(zx_fun) + 2
        _px = np.zeros(deg)
        _px[-1] = 1.0
        for j in range(1, deg-1):
            _px[j] = zx_fun[j-1](k)

        return _px

    def py(self, inv_k):
        k = 1./inv_k
        zy_fun = self.zy_fun
        deg = len(zy_fun) + 2
        _py = np.zeros(deg)
        _py[-1] = 1.0
        for j in range(1, deg-1):
            _py[j] = zy_fun[j-1](k)

        return _py

    def ref_points(self, inv_k):
        x = self.px(inv_k)
        y = self.py(inv_k)

        P = np.zeros([len(x), 2])
        P[:, 0] = x
        P[:, 1] = y

        return P

    def SQ_curve(self, inv_k):
        k = 1./inv_k
        zx_fun = self.zx_fun
        zy_fun = self.zy_fun

        deg = len(zy_fun) + 1
        P = np.zeros([deg+1, 2])
        P[0, :] = [1.0, 1.0]
        P[-1, :] = [0.0, 0.0]
        for i in range(0, deg-1):
            P[i+1, 0] = zx_fun[i](k)
            P[i+1, 1] = zy_fun[i](k)

        t = np.arange(0.0, 1.01, 0.01)
        val = bf.gen(t, deg, P)
        return val

    def NS_to_NQ(self, nS, inv_k):
        #  Normalized storage to find normalized outflow with assigned
        #  k-inverse
        val = self.SQ_curve(inv_k)
        val[:, 0] = val[:, 0]*(1.0 - inv_k**(3./5))/inv_k**(3./5) + 1.  # S
        val[:, 1] = val[:, 1]*(1.0 - inv_k)/inv_k + 1.  # Q

        nQ = np.interp(nS, self.mat_flip(val[:, 0]), self.mat_flip(val[:, 1]))

        return nQ

    def NS_to_t(self, nS, inv_k):
        val = self.SQ_curve(inv_k)
        val[:, 0] = val[:, 0]*(1.0 - (inv_k)**(3./5))/(inv_k)**(3./5) + 1.
        t = np.arange(0.0, 1.01, 0.01)
        nt = np.interp(nS, self.mat_flip(val[:, 0]), self.mat_flip(t))

        return nt

    def SaveResult(self):
        pickle.dump(self.__dict__,
                    open(self.folder + '/' + self.save_name +
                         '_FallLimb.pick', 'w'))


class RisingLimb_gen:
    def __init__(self, length, dt, n, slope, ks, rec_name='',
                 folder='ana_slope/Steady_fit', gen=False):
        self.length = length
        self.dt = dt
        self.n = n
        self.slope = slope
        self.ks = ks
        self.rec_name = rec_name
        self.alpha = 1.0/n*sqrt(slope)
        self.folder = folder
        self.gen = gen

    def load(self):
        _file = self.folder + '/' + self.rec_name + '_CurveRec.flow'
        REC = pickle.load(open(_file, 'r'))
        self.REC = REC

    def fit(self):
        REC = self.REC

        bez_deg = 6
        self.bez_deg = bez_deg

        ks = self.ks
        param = dict()
        for k in ks:
            key_str = 'k=' + str(float(k))
            out_s = REC['S'][key_str]
            out_q = REC['Q'][key_str]
            """
            lst = np.arange(len(out_q)*0.1, len(out_q)*0.96, len(out_q)/20.)
            for j in range(0, len(lst)):
                lst[j] = int(lst[j])
            """
            sq = np.zeros([15, 2])
            sq[:, 0] = np.arange(0, 1.001, 1./14)
            sq[:, 1] = np.interp(sq[:, 0], out_s[0:len(out_q), 1], out_q[:, 1])

            """
            sq = list()
            sq.append([0.0, 0.0])
            for j in range(0, len(lst)):
                sq.append([out_s[lst[j], 1], out_q[lst[j], 1]])
            sq = np.array(sq)
            """

            p, t = bf.fit(sq, bez_deg, 10**-6)
            p[0, :] = [0.0, 0.0]
            p[-1, :] = [1.0, 1.0]

            param.update({'k='+str(k): p})

        self.param = param
        self.param.update({'k': ks})

        if self.rec_name:
            pickle.dump(param, open(self.folder + '/' +
                                    self.rec_name + '_Rise_param.flow', 'w'))

    def gen_rec(self):
        dt = self.dt
        n = self.n
        slope = self.slope
        length = self.length
        ks = self.ks

        out_rec = dict()
        sto_rec = dict()

        for k in range(0, len(ks)):
            ie0 = 10.0
            aa = aw.trans_ana(ie0, ie0*ks[k], length, 201, dt, n, slope)
            aa.run()

            out_s = aa.out_s
            out_q = aa.out_q
            out_s[:, 0] = out_s[:, 0]/aa.max_t
            out_s[:, 1] = (out_s[:, 1] - min(out_s[:, 1]))/(max(out_s[:, 1]) -
                                                            min(out_s[:, 1]))
            out_q[:, 0] = out_q[:, 0]/aa.max_t
            out_q[:, 1] = (out_q[:, 1] - min(out_q[:, 1]))/(max(out_q[:, 1]) -
                                                            min(out_q[:, 1]))

            out_rec.update({'k='+str(ks[k]): out_q})
            sto_rec.update({'k='+str(ks[k]): out_s})

        REC = {'S': sto_rec, 'Q': out_rec}

        self.REC = REC

        if self.rec_name:
            pickle.dump(REC, open(self.folder + '/' +
                                  self.rec_name + '_CurveRec.flow', 'w'))

    def run(self):
        if self.gen:
            self.gen_rec()
            self.fit()
        else:
            self.load()
            self.fit()

    def plot_SQ(self):
        #  max_k = self.max_k
        param = self.param
        length = self.length
        alpha = self.alpha
        ks = self.ks
        m = 5.0/3
        ie0 = 5.0

        dot_s = np.arange(0, 1.0, 0.01)
        for k in range(0, len(ks)):
            ie1 = ie0*ks[k]
            max_q = ie1*length
            max_s = m/(m+1)*(ie1*length**(m+1)/alpha)*(1/m)
            min_s = m/(m+1)*(ie0*length**(m+1)/alpha)*(1/m)

            _str = 'k=' + str(ks[k])
            p = param[_str]
            dot_s = np.arange(0, 1.001, 0.01)
            dot_q = bf.gen(dot_s, 6, p)

            dot_s = (dot_s*(max_s - min_s)/max_s + min_s/max_s)*max_s
            q = (dot_q*(1.0 - 1.0/k) + 1.0/k)*max_q

            plt.plot(q, dot_s)

    def plot_norm_SQ(self):
        #  max_k = self.max_k
        param = self.param
        ks = self.ks
        REC = self.REC

        dot_s = np.arange(0, 0.981, 0.01)
        dot_s = dot_s.tolist()
        dot_s.append(1.0)
        dot_s - np.array(dot_s)
        for k in range(0, len(ks)):
            _str = 'k=' + str(ks[k])
            p = param[_str]
            out_s = REC['S'][_str]
            out_q = REC['Q'][_str]
            dot_q = bf.gen(dot_s, self.bez_deg, p)

            plt.plot(dot_q[:, 1], dot_q[:, 0], '-b')
            plt.plot(out_q[:, 1], out_s[0:len(out_q[:, 1]), 1], '--r')

    def Steady_SQ(self):
        n = self.n
        slope = self.slope
        length = self.length
        m = 5.0/3
        ks = self.ks

        sq = list()
        alpha = 1.0/n*sqrt(slope)
        self.alpha = alpha
        for ie in range(1, int(5.0*ks[-1])):
            ie = float(ie)/3600*10**-3
            _q = alpha*(ie*length/alpha)**(1/m)
            _s = m/(m+1)*(ie*length**(m+1)/alpha)**(1/m)

            sq.append([_s, _q])

        self.SQ_std = np.array(sq)


class RisingLimb_fit:
    """Load parameters of rising limb s-q curve, then build interpolation
    function of curve."""
    def __init__(self, rec_name='', save_name='',
                 folder='ana_slope/Steady_fit', **kwargs):

        param = pickle.load(
            open(folder + '/' + rec_name + '_Rise_param.flow', 'r'))

        self.param = param
        self.folder = folder
        self.ks = param['k']
        self.save_name = save_name

    def run(self):
        self.p_fit()
        self.SaveResult()

    def load(self):
        save_name = self.save_name
        folder = self.folder

        tmp_dict = pickle.load(open(folder + '/' + save_name +
                                    '_RiseLimb.pick', 'r'))
        self.__dict__.update(tmp_dict)

    def mat_flip(self, mat):
        rev_mat = np.fliplr([mat])[0]

        return rev_mat

    def plot_p(self, r):
        param = self.param  # Bezeir curve reference points
        ks = self.ks

        px = list()
        py = list()

        for k in range(0, len(ks)):
            key = 'k=' + str(ks[k])
            p = param[key]

            px.append([ks[k], p[r, 0]])
            py.append([ks[k], p[r, 1]])

        px = np.array(px)
        py = np.array(py)
        plt.plot(px[:, 0], px[:, 1], '-r', label='x')
        plt.plot(py[:, 0], py[:, 1], '-k', label='y')

    def p_fit(self):
        param = self.param
        ks = self.ks

        deg = len(param['k=5.0']) - 1

        _px = np.zeros([len(ks), deg])
        _py = np.zeros([len(ks), deg])
        for k in range(0, len(ks)):
            key = 'k=' + str(ks[k])
            p = param[key]

            _px[k, 0] = ks[k]
            _px[k, 1:deg] = p[1:deg, 0]
            _py[k, 0] = ks[k]
            _py[k, 1:deg] = p[1:deg, 1]

        self._px = _px
        self._py = _py

        self.zx = self.fit_px()
        self.zy = self.fit_py()

        zx_fun = list()
        zy_fun = list()
        for i in range(0, len(self.zx)):
            zx_fun.append(np.poly1d(self.zx[i]))
            zy_fun.append(np.poly1d(self.zy[i]))
        self.zx_fun = zx_fun
        self.zy_fun = zy_fun

    def fit_px(self):
        _px = self._px

        z = list()
        for i in range(1, len(_px[0])):
            z.append(np.polyfit(_px[:, 0], _px[:, i], 6))

        return z

    def fit_py(self):
        _py = self._py

        z = list()
        for i in range(1, len(_py[0])):
            z.append(np.polyfit(_py[:, 0], _py[:, i], 6))

        return z

    def plot_fit_p(self):

        zx_fun = self.zx_fun
        zy_fun = self.zy_fun
        _px = self._px
        _py = self._py

        k = np.arange(min(self.ks), max(self.ks)+0.1, 0.1)

        plt.figure
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(6, 6))

        axes[0, 0].plot(k, zx_fun[0](k), '-r')
        axes[0, 0].plot(k, zy_fun[0](k), '-k')
        axes[0, 0].plot(_px[:, 0], _px[:, 1], '--r')
        axes[0, 0].plot(_py[:, 0], _py[:, 1], '--k')
        axes[0, 0].set_title('$P_1$')

        axes[0, 1].plot(k, zx_fun[1](k), '-r')
        axes[0, 1].plot(k, zy_fun[1](k), '-k')
        axes[0, 1].plot(_px[:, 0], _px[:, 2], '--r')
        axes[0, 1].plot(_py[:, 0], _py[:, 2], '--k')
        axes[0, 1].set_title('$P_2$')

        axes[0, 2].plot(k, zx_fun[2](k), '-r')
        axes[0, 2].plot(k, zy_fun[2](k), '-k')
        axes[0, 2].plot(_px[:, 0], _px[:, 3], '--r')
        axes[0, 2].plot(_py[:, 0], _py[:, 3], '--k')
        axes[0, 2].set_title('$P_3$')

        axes[1, 0].plot(k, zx_fun[3](k), '-r')
        axes[1, 0].plot(k, zy_fun[3](k), '-k')
        axes[1, 0].plot(_px[:, 0], _px[:, 4], '--r')
        axes[1, 0].plot(_py[:, 0], _py[:, 4], '--k')
        axes[1, 0].set_title('$P_4$')

        axes[1, 1].plot(k, zx_fun[4](k), '-r')
        axes[1, 1].plot(k, zy_fun[4](k), '-k')
        axes[1, 1].plot(_px[:, 0], _px[:, 5], '--r')
        axes[1, 1].plot(_py[:, 0], _py[:, 5], '--k')
        axes[1, 1].set_title('$P_5$')

        plt.show()

    def SQ_curve(self, k):
        zx_fun = self.zx_fun
        zy_fun = self.zy_fun

        deg = len(zy_fun) + 1
        P = np.zeros([deg+1, 2])
        P[0, :] = [0.0, 0.0]
        P[-1, :] = [1.0, 1.0]
        for i in range(0, deg-1):
            P[i+1, 0] = zx_fun[i](k)
            P[i+1, 1] = zy_fun[i](k)

        t = np.arange(0.0, 1.01, 0.01)
        val = bf.gen(t, deg, P)
        return val

    def NS_to_t(self, nS, k):
        val = self.SQ_curve(k)
        val[:, 0] = val[:, 0]*(1.0 - (1./k)**(3./5)) + (1./k)**(3./5)
        t = np.arange(0.0, 1.01, 0.01)
        nt = np.interp(nS, val[:, 0], t)

        return nt

    def NS_to_NQ(self, nS, k):
        #  Normalized storage to find normalized outflow with assigned k
        val = self.SQ_curve(k)
        val[:, 0] = val[:, 0]*(1.0 - (1./k)**(3./5)) + (1./k)**(3./5)
        val[:, 1] = val[:, 1]*(1.0 - 1./k) + 1./k

        nQ = np.interp(nS, val[:, 0], val[:, 1])

        return nQ

    def NSNQ_to_norm(self, nS, nQ, k):
        norm_nS = (nS - (1./k)**(3./5))/(1.0 - (1./k)**(3./5))
        norm_nQ = (nQ - 1./k)/(1.0 - 1./k)

        return norm_nS, norm_nQ

    def px(self, k):
        zx_fun = self.zx_fun
        deg = len(zx_fun) + 2
        _px = np.zeros(deg)
        _px[-1] = 1.0
        for j in range(1, deg-1):
            _px[j] = zx_fun[j-1](k)

        return _px

    def py(self, k):
        zy_fun = self.zy_fun
        deg = len(zy_fun) + 2
        _py = np.zeros(deg)
        _py[-1] = 1.0
        for j in range(1, deg-1):
            _py[j] = zy_fun[j-1](k)

        return _py

    def ref_points(self, k):
        x = self.px(k)
        y = self.py(k)

        P = np.zeros([len(x), 2])
        P[:, 0] = x
        P[:, 1] = y

        return P

    def derivate(self, t_star):
        px = self.zx
        py = self.zy

        dx = np.zeros(len(px)+2)
        dy = np.zeros(len(px)+2)
        for i in range(0, len(px)):
            dx[i+1] = pf.dev(t_star, px[i])
            dy[i+1] = pf.dev(t_star, py[i])

        dPdk = np.zeros([len(dx), 2])
        dPdk[:, 0] = dx[:]
        dPdk[:, 1] = dy[:]
        dPdk[-1, :] = [0.0, 0.0]
        return dPdk

    def SaveResult(self):
        pickle.dump(self.__dict__, open(self.folder + '/' + self.save_name +
                                        '_RiseLimb.pick', 'w'))


class fit_steady:
    def __init__(self, save_name,
                 folder='ana_slope/Steady_fit'):

        self.save_name = save_name
        self.folder = folder

    def param(self, length, n, dt, slope):
        self.alpha = 1.0/n*sqrt(slope)
        self.L = length
        self.dt = dt
        self.m = 5.0/3
        self.steady_SQ()
        self.steady_time()

    def load(self):
        tmp_dict = pickle.load(open(self.folder + '/' +
                                    self.save_name + '_Steady.pick', 'r'))
        self.__dict__.update(tmp_dict)

    def bzCurve(self):
        # Start from zero
        self.steady_q.tolist().insert(0, 0.)
        self.steady_s.tolist().insert(0, 0.)

        t = np.arange(0., 1.0001, 1./10)
        Q = np.interp(t, self.ie/max(self.ie), self.steady_q)
        S = np.interp(t, self.ie/max(self.ie), self.steady_s)

        std_SQ = np.array([zip(Q, S)])[0]
        p, t = bf.fit(std_SQ, 2, 1.0E-8)
        print p

        p[0] = [0, 0]
        p[-1] = std_SQ[-1]
        t = np.arange(0., 1.0001, 1./100)
        f_std_SQ = bf.gen(t, 2, p)

        return f_std_SQ, t

    def steady_SQ(self):
        length = self.L
        m = self.m
        alpha = self.alpha
        ie = np.arange(1.0, 2000.1, 1.0)
        ie = ie/3600*10**-3
        steady_sto = np.zeros(2000)
        steady_out = np.zeros(2000)

        for i in range(0, len(ie)):
            steady_sto[i] = m/(m+1)*(ie[i]*length**(m+1)/alpha)**(1/m)
            steady_out[i] = ie[i]*length

        self.steady_q = steady_out
        self.steady_s = steady_sto
        self.ie = ie
        self.SaveResult()

    def steady_time(self):
        L = self.L
        alpha = self.alpha
        ie = self.ie
        m = self.m

        steady_time = np.zeros(len(ie))
        for i in range(0, len(ie)):
            steady_time[i] = (L*ie[i]**(1-m)/alpha)**(1/m)

        self.steady_t = steady_time

    def ie_to_storage(self, ie_input):
        ie = self.ie
        steady_s = self.steady_s

        return np.interp(ie_input, ie, steady_s)

    def ie_to_outflow(self, ie_input):
        ie = self.ie
        steady_q = self.steady_q

        return np.interp(ie_input, ie, steady_q)

    def ie_to_t(self, ie_input):
        ie = self.ie
        steady_t = self.steady_t

        return np.interp(ie_input, ie, steady_t)

    def s_to_ie(self, s):
        ie = self.ie
        steady_s = self.steady_s

        return np.interp(s, steady_s, ie)

    def gen(self):
        L = self.L
        m = self.m
        dt = 1.0
        ie = 100.0/3600*10**-3
        alpha = self.alpha

        t = 0.0
        time = list()
        outflow = list()
        storage = list()
        while t < self.ie_to_t(ie):
            h = ie*t
            q = alpha*h**m
            xs = alpha*ie**(m-1)*t**m
            s = m/(m+1)*(ie/alpha)**(1/m)*xs**((m+1)/m) + (L - xs)*ie*t
            time.append(t)
            outflow.append(q)
            storage.append(s)

            t = t + dt

        time = np.array(time)
        outflow = np.array(outflow)
        storage = np.array(storage)

        time = time/self.ie_to_t(ie)
        outflow = outflow/self.ie_to_outflow(ie)
        storage = storage/self.ie_to_storage(ie)

        Co = np.polyfit(time, outflow, 4)
        Cs = np.polyfit(time, storage, 4)

        self.Co = Co
        self.Cs = Cs

    def SaveResult(self):
        save_name = self.save_name
        folder = self.folder

        pickle.dump(self.__dict__, open(folder + '/' + save_name +
                                        '_Steady.pick', 'w'))
