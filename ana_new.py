
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, fabs, floor
from parameters import main_folder
import time
import os
import copy

"""
Solving the analytical solution of kinematic wave overland flow from a
steady-state under steady rainfall ie0 to another steady-state under rainfall
ie1.
                                                                            """


class trans_ana:

    def __init__(self, ie0, ie1, t_dist, x_space, dt, n, slope, **kwargs):
        """
        ie0: effective rainfall before rainfall intensity shifts.
        ie1: effective rainfall after rainfall intensity shifts.
        t_dist: total distance of slope(from 0 to L).
        x_space: number of grid points(boundaries included).
        dt: time space.
        n: manning coefficient.
        slope: slope for manning equation.
        """
        self.t_dist = t_dist  # Total distance of slope
        self.x_space = x_space  # Number of grid points(include boundary)
        self.dt = dt  # Time space
        alpha = 1.0/n*sqrt(slope)
        self.alpha = alpha  # paramter of manning
        self.m = 5.0/3  # power of Manning equation

        if 'folder' in kwargs:
            self.folder = kwargs['folder']
            if not os.path.isdir(kwargs['folder']):
                os.mkdir(kwargs['folder'])
                os.mkdir(kwargs['folder']+'/depth_change')
        else:
            folder = 'calc_result/' + str(ie0) + 'to' + str(ie1)
            if not os.path.isdir(folder):
                os.mkdir(folder)
                os.mkdir(folder + '/depth_change')
            self.folder = folder

        ie0 = ie0*10**-3/3600
        ie1 = ie1*10**-3/3600
        #  Change rainfall intensity unit from mm/hr to m/s.

        self.ie0 = ie0  # Initial rainfall
        self.ie1 = ie1  # Secondary rainfall

        if 'plot' in kwargs:
            if kwargs['plot']:
                self.plot = True
            else:
                self.plot = False
        else:
            self.plot = False

        self.get_max_q()
        self.get_max_s()

        self.method = 'simpson'
        if (x_space-1) % 2 != 0:
            print('Number total column-1 is not even, using trapezoid method.')
            self.method = 'trapezoid'
        if 'method' in kwargs:
            method = kwargs['method']
            self.method = method

    def get_max_q(self):
        ie0 = self.ie0
        ie1 = self.ie1

        if ie0 > ie1:
            self.max_q = self.t_dist*ie0
            self.min_q = self.t_dist*ie1
        else:
            self.max_q = self.t_dist*ie1
            self.min_q = self.t_dist*ie0

    def get_max_s(self):
        ie0 = self.ie0
        ie1 = self.ie1
        m = self.m
        alpha = self.alpha

        if ie0 > ie1:
            self.max_s = m/(m+1)*(ie0*self.t_dist**(m+1)/alpha)**(1/m)
            self.min_s = m/(m+1)*(ie1*self.t_dist**(m+1)/alpha)**(1/m)
        else:
            self.max_s = m/(m+1)*(ie1*self.t_dist**(m+1)/alpha)**(1/m)
            self.min_s = m/(m+1)*(ie0*self.t_dist**(m+1)/alpha)**(1/m)

    def grids(self):
        t_dist = self.t_dist
        x_space = self.x_space
        dx = t_dist/(x_space-1)
        x = np.zeros(x_space)
        for i in range(0, x_space):
            x[i] = i*dx
        self.grd = x

    def initial_state(self):
        ie0 = self.ie0
        grd = self.grd
        alpha = self.alpha
        d0 = np.zeros(len(grd))

        for i in range(0, len(d0)):
            d0[i] = alpha**-1.5*(ie0*(grd[i] - grd[0]))**1.5
        self.d0 = d0

    def test_bi(self, xa, t, rh, ih):
        #  Function of bisection method, test if solution between rh and ih.
        alpha = self.alpha
        m = self.m
        ie0 = self.ie0
        ie1 = self.ie1

        f1 = (xa - alpha/ie0*rh**m - alpha/ie1*(rh + ie1*t)**m
              + alpha/ie1*rh**m)
        f2 = (xa - alpha/ie0*ih**m - alpha/ie1*(ih + ie1*t)**m
              + alpha/ie1*ih**m)
        return f1*f2

    def bisec(self, xa, t, ih, rh):
        #  Bisection method
        f1f2 = self.test_bi(xa, t, rh, ih)

        if f1f2 < 0:
            n_rh = (rh + ih)*0.5
            f1f2_2 = self.test_bi(xa, t, n_rh, ih)
            if f1f2_2 < 0:
                return [n_rh, ih]
            elif f1f2_2 > 0:
                return [n_rh, rh]
        elif f1f2 > 0:
            print('Error, no solution available.')
            time.sleep(10)
            raise(ValueError)

    def iteration_counter(self, ih, f, df):
        counter = self.counter
        if counter < 30:
            rh = ih - f/df
        elif counter >= 30:
            rh = ih - f/df*0.5
        elif counter >= 50:
            rh = ih - f/df*0.1
        elif counter >= 200:
            rh = ih - f/df*0.01
        elif counter >= 500:
            rh = ih - f/df*0.001

        return rh

    def bisec_iter(self, xa, t, allow_err):
        #  Bisection method iteration.
        ie1 = self.ie1
        ie0 = self.ie0
        alpha = self.alpha
        m = self.m
        counter = self.counter

        ih = 0
        if ie0 > ie1:
            rh = (ie0*xa/alpha)**(1/m)
        else:
            rh = (ie1*xa/alpha)**(1/m)
        err = 100.0
        while err > allow_err:
            try:
                [rh, ih] = self.bisec(xa, t, rh, ih)
            except(ValueError):
                return None
            err = fabs(rh - ih)
            counter = counter + 1
            self.counter = counter

            if self.counter > 500:
                print(rh, ih)

            if self.counter > 1000:
                print(rh, ih)
                time.sleep(300)

        return rh

    def iteration(self, xa, t, ih):
        """
        Input assigned position and time , backtrace to its starting position
        and depth at rainfall intensity shift."""
        #  ih is the initial guess of starting position depth.
        ie0 = self.ie0
        ie1 = self.ie1
        alpha = self.alpha
        m = self.m
        err = 100.0
        last_err = err + 1.0
        self.counter = 0

        while err > 10**-8:
            #  Newton method iteration.
            try:
                f = (xa - alpha/ie0*ih**m - alpha/ie1*(ih + ie1*t)**m
                     + alpha/ie1*ih**m)
                df = (-(alpha/ie0 + alpha/ie1)*m*ih**(m-1)
                      - alpha/ie1*m*(ih + ie1*t)**m-1)

                rh = self.iteration_counter(ih, f, df)

                err = fabs(rh - ih)
            except(ValueError):
                rh = self.bisec_iter(xa, t, 10**-8)
                return rh
            """
            Newton method easily goes divergent if df is too small, and not
            getting accurate value, thus if error between current and last
            iteration expands, use bisection method instead.
            """
            if last_err > err:
                ih = rh
                last_err = err
                self.counter = self.counter + 1
            else:
                rh = self.bisec_iter(xa, t, 10**-8)
                return rh

        return rh

    def gen_depth_graph(self, h_on_x, t):
        x_grd = self.grd
        ie0 = self.ie0
        ie1 = self.ie1
        alpha = self.alpha
        m = self.m
        folder = self.folder
        out_q = self.out_q

        balance1 = (ie0*x_grd/alpha)**(1/m)
        balance2 = (ie1*x_grd/alpha)**(1/m)

        ie0 = ie0*3600*10**3
        ie1 = ie1*3600*10**3

        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))
        axes[0].plot(x_grd, h_on_x)
        axes[0].plot(x_grd, balance1, 'r--', linewidth=1.0,
                     label='$i_e=$'+str(ie0)+'  Steady-state')
        axes[0].plot(x_grd, balance2, 'g--', linewidth=1.0,
                     label='$i_e=$'+str(ie1)+'  Steady-state')
        axes[0].legend(loc=2)
        axes[0].set_xlabel('m')
        axes[0].set_ylabel('m')
        axes[0].set_title('Depth Along Slope')

        c_q = np.interp(t, out_q[:, 0], out_q[:, 1])
        axes[1].plot(out_q[:, 0], out_q[:, 1], '-b')
        axes[1].plot(t, c_q, 'bo')
        axes[1].set_xlabel('Time (s)')
        axes[1].set_ylabel('Flowrate ($m^3/s$)')
        axes[1].set_xlim([0, max(out_q[:, 0])])
        axes[1].set_title('Flowrate - Time')
        fig.suptitle('t = ' + str(t))
        filename = str(t)+'.png'
        fig.savefig(folder + '/depth_change/' + filename.zfill(10))
        plt.close()

    def in_process(self):
        x_grd = self.grd
        L = self.t_dist

        dt = self.dt
        ie0 = self.ie0
        ie1 = self.ie1
        alpha = self.alpha
        m = self.m
        folder = self.folder

        h_on_x = np.zeros(len(x_grd))

        balance_t = (L*ie1**(1-m)/alpha)**(1/m)
        self.max_t = balance_t

        xr_t = list()
        xf_t = list()
        t_rec = list()
        depth_rec = list()

        out_h = list()
        out_q = list()
        out_s = list()
        t = 0.0
        diff = np.zeros(len(x_grd))
        balance2 = (ie1*x_grd/alpha)**(1/m)
        hl = (ie0*L/alpha)**(1/m)
        xf = 0
        s0 = (ie0/alpha)**(1/m)*m/(m+1)*L**((m+1)/m)
        q0 = ie0*L
        out_h.append([0.0, (ie0*L/alpha)**(1/m)])
        out_s.append([0.0, s0])
        out_q.append([0.0, q0])
        c_s = s0

        while t < balance_t or xf <= L:
            t = t + dt
            hl = self.iteration(L, t, hl)

            xf = alpha/ie1*(ie1*t)**m
            if float(floor(t)) == t:
                print(t, xf)

            for i in range(1, len(h_on_x)-1):
                if x_grd[i] > xf:
                    h_on_x[i] = self.iteration(x_grd[i], t, ie1*t) + ie1*t
                else:
                    h_on_x[i] = (ie1*x_grd[i]/alpha)**(1/m)

            try:
                h_on_x[-1] = hl + ie1*t
            except(TypeError):
                break

        #  Record on xr(balance front) position and xf(start from x=0)
        #  position.
            xr_list = list()
            for i in range(0, len(x_grd)):
                diff[i] = fabs(h_on_x[i] - balance2[i])

            for i in range(0, len(x_grd)):
                if diff[i] <= 10**-9:
                    xr_list.append(x_grd[i])
            xr = max(xr_list)
            xr_t.append(xr)
            xf_t.append(xf)
            t_rec.append(t)
            self.xr = xr_t
            self.xf = xf
            self.t = t_rec

            depth_rec.append([t, copy.copy(h_on_x)])

            out_h.append([t, hl])
            out_q.append([t, alpha*h_on_x[-1]**m])

            c_s = self.simpson(h_on_x)
            out_s.append([t, c_s])

        out_h = np.array(out_h)
        out_q = np.array(out_q)
        out_s = np.array(out_s)

        self.out_q = out_q
        self.out_s = out_s

        #  Plot xr-xf comparison and catchment outflow rate.
        if self.plot:
            self.plot_depth(depth_rec)
            plt.figure
            plt.plot(t_rec, xr_t, label='$x_r$')
            plt.plot(t_rec, xf_t, label='$x_f$')
            plt.legend(loc=2)
            plt.xlabel('t')
            plt.ylabel('x (m)')
            plt.title(str(ie0*3600*1000) + ' to ' + str(ie1*3600*1000) +
                      '$x_f$ and $x_r$')
            plt.savefig(folder + '/xr_and_xf.png')
            plt.close()

            plt.figure
            plt.plot(out_q[:, 0], out_q[:, 1])
            plt.title('Outflow rate')
            plt.xlabel('t')
            plt.ylabel('q $(m^3)/s$')
            plt.title(str(ie0*3600*1000) + ' to ' + str(ie1*3600*1000) +
                      ' catchment outflow rate')
            plt.savefig(folder + '/outflow.png')
            plt.close()

            plt.figure
            plt.plot(out_s[:, 0], out_s[:, 1])
            plt.title('Storage')
            plt.xlabel('t')
            plt.ylabel('S $(m^3)$')
            plt.title(str(ie0*3600*1000) + ' to ' + str(ie1*3600*1000) +
                      ' catchment storage')
            plt.savefig(folder + '/storage.png')
            plt.close()

    def plot_depth(self, depth_rec):
        for i in range(0, len(depth_rec)):
            if str(depth_rec[i][0])[-1] == '0':
                self.gen_depth_graph(depth_rec[i][1], depth_rec[i][0])

    def simple_iter(self, h_on_x, t):
        alpha = self.alpha
        ie0 = self.ie0
        ie1 = self.ie1
        m = self.m
        x = self.grd
        grds = self.x_space
        xz_pos = np.zeros(grds)
        init_h = (x*ie0/alpha)**(1/m)

        xf = alpha*ie1**(m-1)*t**m
        for i in range(0, grds):
            xz_pos[i] = (x[i] + alpha/ie1*(init_h[i] + ie1*t)**m -
                         alpha/ie1*init_h[i]**m)
        for j in range(1, grds):
            c_x = np.interp(x[j], xz_pos, x)
            if x[j] > xf:
                h_on_x[j] = (c_x*ie0/alpha)**(1/m) + ie1*t
            else:
                h_on_x[j] = (x[j]*ie1/alpha)**(1/m)

        return h_on_x

    def in_process2(self):
        x_grd = self.grd
        L = self.t_dist

        dt = self.dt
        ie0 = self.ie0
        ie1 = self.ie1
        alpha = self.alpha
        m = self.m
        folder = self.folder

        h_on_x = np.zeros(len(x_grd))

        balance_t = (L*ie1**(1-m)/alpha)**(1/m)
        self.max_t = balance_t

        xr_t = list()
        xf_t = list()
        t_rec = list()

        out_h = list()
        out_q = list()
        out_s = list()
        t = 0.0
        diff = np.zeros(len(x_grd))
        balance2 = (ie1*x_grd/alpha)**(1/m)
        xf = 0
        s0 = (ie0/alpha)**(1/m)*m/(m+1)*L**((m+1)/m)
        out_s.append([0.0, s0])
        c_s = s0

        while t < balance_t or xf < L:
            t = t + dt

            xf = alpha/ie1*(ie1*t)**m
            print(t, xf)

            h_on_x = self.simple_iter(h_on_x, t)

            if str(t)[-1] == '0' and self.plot:
                self.gen_depth_graph(h_on_x, t, xf)

        #  Record on xr(balance front) position and xf(start from x=0)
        #  position.
            xr_list = list()
            for i in range(0, len(x_grd)):
                diff[i] = fabs(h_on_x[i] - balance2[i])

            for i in range(0, len(x_grd)):
                if diff[i] <= 10**-9:
                    xr_list.append(x_grd[i])
            xr = max(xr_list)
            xr_t.append(xr)
            xf_t.append(xf)
            t_rec.append(t)
            self.xr = xr_t
            self.xf = xf
            self.t = t_rec

            out_h.append([t, h_on_x[-1]])
            out_q.append([t, alpha*h_on_x[-1]**m])

            c_s = self.storage(h_on_x)
            out_s.append([t, c_s])

        out_h = np.array(out_h)
        out_q = np.array(out_q)
        out_s = np.array(out_s)

        self.out_q = out_q
        self.out_s = out_s

        #  Plot xr-xf comparison and catchment outflow rate.
        plt.figure
        plt.plot(t_rec, xr_t, label='$x_r$')
        plt.plot(t_rec, xf_t, label='$x_f$')
        plt.legend(loc=2)
        plt.xlabel('t')
        plt.ylabel('x (m)')
        plt.title(str(ie0*3600*1000) + ' to ' + str(ie1*3600*1000) + '$x_f$ and \
$x_r$')
        plt.savefig(folder + '/xr_and_xf.png')
        plt.close()

        plt.figure
        plt.plot(out_q[:, 0], out_q[:, 1])
        plt.title('Outflow rate')
        plt.xlabel('t')
        plt.ylabel('q $(m^3)/s$')
        plt.title(str(ie0*3600*1000) + ' to ' + str(ie1*3600*1000) + ' catchment\
 outflow rate')
        plt.savefig(folder + '/outflow.png')
        plt.close()

        plt.figure
        plt.plot(out_s[:, 0], out_s[:, 1])
        plt.title('Storage')
        plt.xlabel('t')
        plt.ylabel('S $(m^3)$')
        plt.title(str(ie0*3600*1000) + ' to ' + str(ie1*3600*1000) + ' catchment\
 storage')
        plt.savefig(folder + '/storage.png')
        plt.close()

    def run(self):
        """Operation function"""
        self.grids()
        self.initial_state()
        self.in_process()

    def run2(self):
        """Operation function"""
        self.grids()
        self.initial_state()
        self.in_process2()

    def storage(self, h_on_x):
        L = self.t_dist
        x_space = self.x_space
        dx = L/(x_space-1)

        if self.method == 'trapezoid' or self.method == 'trape':
            Sto = 0.0
            for i in range(1, len(h_on_x)):
                Sto = Sto + (h_on_x[i] + h_on_x[i-1])*0.5*dx
        elif self.method == 'simpson':
            Sto = self.simpson(h_on_x)

        return Sto

    def simpson(self, h_on_x):
        x_space = self.x_space
        n = x_space - 1
        L = self.t_dist
        dx = L/n

        mult = np.zeros(x_space)
        mult[0] = 1.0
        mult[-1] = 1.0
        #  Multiplier of simpson's 1/2 method
        for i in range(1, x_space-1):
            if i % 2 == 0:
                mult[i] = 2.0
            elif (i-1) % 2 == 0:
                mult[i] = 4.0

        return dx/3.*np.dot(mult, h_on_x)

    def plot_ST(self, LT='-r', plt_label=''):
        S = self.out_s
        plt.plot(S[:, 0], S[:, 1], LT, label=plt_label)

    def plot_QT(self, LT='-r', plt_label=''):
        t = self.t
        Q = self.out_q
        plt.plot(t, Q[:, 1], LT, label=plt_label)


class unsteady_trans:
    def __init__(self, ie0, ie1, t_dist, x_space, x_s, dt, n, slope, **kwargs):

        self.L = t_dist
        self.grds = x_space
        self.x_s = x_s
        self.dt = dt
        alpha = 1.0/n*sqrt(slope)
        self.alpha = alpha
        self.m = 5.0/3

        if 'folder' in kwargs:
            self.folder = kwargs['folder']
            if not os.path.isdir(kwargs['folder']):
                os.mkdir(kwargs['folder'])
                os.mkdir(kwargs['folder']+'/depth_change')
        else:
            folder = (main_folder + '/' + str(ie0) + '_unsteady_to_' +
                      str(ie1) + '_xs_' + str(x_s))
            if not os.path.isdir(folder):
                os.mkdir(folder)
                os.mkdir(folder + '/depth_change')
            self.folder = folder

        ie0 = ie0/3600.0/1000
        ie1 = ie1/3600.0/1000
        #  x_s: balance position at rainfall intensity shifts.
        self.ie0 = ie0
        self.ie1 = ie1

        if 'plot' in kwargs:
            if kwargs['plot']:
                self.plot = True
            else:
                self.plot = False
        else:
            self.plot = False

        self.get_max_q()
        self.get_max_s()
        self.xs_to_L()

    def storage(self, h_on_x):
        L = self.L
        x_space = self.grds
        dx = L/(x_space-1)

        Sto = 0.0
        for i in range(1, len(h_on_x)):
            Sto = Sto + (h_on_x[i] + h_on_x[i-1])*0.5*dx

        return Sto

    def xs_to_L(self):
        alpha = self.alpha
        xs = self.x_s
        ie1 = self.ie1
        m = self.m
        L = self.L

        t = np.arange(1, (L*ie1**(1-m)/alpha)**(1/m), 1.0)
        hs = (xs*ie1/alpha)**(1/m)
        r_err = xs - L + alpha/ie1*(hs + ie1*t)**m - alpha/ie1*hs**m
        t_L = np.interp(0.0, r_err, t)

        ie0 = self.ie0
        h_L = (ie0*xs/alpha)**(1/m) + ie1*t_L
        self.q_xs = alpha*h_L**m

        self.t_xs = t_L

    def x_grds(self):
        L = self.L
        grds = self.grds
        grd = np.zeros(grds)
        for i in range(0, len(grd)):
            grd[i] = L/(grds-1)*i
        self.grd = grd

    def initial_cond(self):
        ie0 = self.ie0
        x_s = self.x_s
        grd = self.grd
        m = self.m
        alpha = self.alpha

        d0 = np.zeros(len(grd))
        for i in range(0, len(grd)):
            if grd[i] <= x_s:
                d0[i] = (ie0*grd[i]/alpha)**(1/m)
            else:
                d0[i] = (ie0*x_s/alpha)**(1/m)

        self.d0 = d0

    def bisec(self, xa, t, ih, rh):
        #  Bisection method
        f1f2 = self.test_bi(xa, t, rh, ih)

        if f1f2 < 0:
            n_rh = (rh + ih)*0.5
            f1f2_2 = self.test_bi(xa, t, n_rh, ih)
            if f1f2_2 < 0:
                return [n_rh, ih]
            elif f1f2_2 > 0:
                return [n_rh, rh]
        elif f1f2 > 0:
            print('Error, no solution available.')
            time.sleep(10)
            raise ValueError

    def iteration_counter(self, ih, f, df):
        counter = self.counter
        if counter < 30:
            rh = ih - f/df
        elif counter >= 30:
            rh = ih - f/df*0.5
        elif counter >= 50:
            rh = ih - f/df*0.1
        elif counter >= 200:
            rh = ih - f/df*0.01
        elif counter >= 500:
            rh = ih - f/df*0.001

        return rh

    def test_bi(self, xa, t, rh, ih):
        #  Function of bisection method, test if solution between rh and ih.
        alpha = self.alpha
        m = self.m
        ie0 = self.ie0
        ie1 = self.ie1

        f1 = (xa - alpha/ie0*rh**m - alpha/ie1*(rh + ie1*t)**m
              + alpha/ie1*rh**m)
        f2 = (xa - alpha/ie0*ih**m - alpha/ie1*(ih + ie1*t)**m
              + alpha/ie1*ih**m)
        return f1*f2

    def bisec_iter(self, xa, t, allow_err):
        #  Bisection method iteration.
        ie1 = self.ie1
        ie0 = self.ie0
        alpha = self.alpha
        m = self.m
        counter = self.counter

        ih = 0
        if ie0 > ie1:
            rh = (ie0*xa/alpha)**(1/m)
        else:
            rh = (ie1*xa/alpha)**(1/m)
        err = 100.0
        while err > allow_err:
            try:
                [rh, ih] = self.bisec(xa, t, rh, ih)
            except(ValueError):
                return None
            err = fabs(rh - ih)
            counter = counter + 1
            self.counter = counter

            if self.counter > 500:
                print(rh, ih)

            if self.counter > 1000:
                print(rh, ih)
                time.sleep(300)

        return rh

    def iteration(self, xa, t, ih):
        ie0 = self.ie0
        ie1 = self.ie1
        alpha = self.alpha
        m = self.m
        err = 100.0
        last_err = err + 1.0
        self.counter = 0

        while err > 10**-8:
            #  Newton method iteration.
            try:
                f = (xa - alpha/ie0*ih**m - alpha/ie1*(ih + ie1*t)**m
                     + alpha/ie1*ih**m)
                df = (-(alpha/ie0 + alpha/ie1)*m*ih**(m-1)
                      - alpha/ie1*m*(ih + ie1*t)**m-1)
            except(ValueError):
                rh = self.bisec_iter(xa, t, 10**-8)
                return rh

            rh = self.iteration_counter(ih, f, df)

            err = fabs(rh - ih)
            #  iteration
            if last_err > err:
                try:
                    ih = rh
                    last_err = err
                    self.counter = self.counter + 1
                except(ValueError):
                    rh = self.bisec_iter(xa, t, 10**-8)
                    return rh
            else:
                rh = self.bisec_iter(xa, t, 10**-8)
                return rh
        return rh

    def h_on_x(self, x_grd, xf, x_sta, t):
        ie0 = self.ie0
        ie1 = self.ie1
        alpha = self.alpha
        x_s = self.x_s
        m = self.m

        if x_grd > xf:
            if x_grd >= x_sta:
                depth = (ie0*x_s/alpha)**(1/m) + ie1*t
            else:
                depth = self.iteration(x_grd, t, ie1*t) + ie1*t
        else:
            depth = (ie1*x_grd/alpha)**(1/m)

        return depth

    def get_max_q(self):
        ie0 = self.ie0
        ie1 = self.ie1

        if ie0 > ie1:
            self.max_q = self.L*ie0
        else:
            self.max_q = self.L*ie1

    def get_max_s(self):
        ie0 = self.ie0
        ie1 = self.ie1
        m = self.m
        alpha = self.alpha

        if ie0 > ie1:
            self.max_s = m/(m+1)*(ie0*self.L**(m+1)/alpha)**(1/m)
        else:
            self.max_s = m/(m+1)*(ie1*self.L**(m+1)/alpha)**(1/m)

    def calc_process(self):
        x_grd = self.grd
        L = self.L

        dt = self.dt
        ie0 = self.ie0
        ie1 = self.ie1
        alpha = self.alpha
        m = self.m
        x_s = self.x_s
        folder = self.folder

        h_on_x = np.zeros(len(x_grd))

        balance_t = (L*ie1**(1-m)/alpha)**(1/m)
        self.max_t = balance_t

        out_h = list()
        out_q = list()
        out_s = list()
        t = 0.0
        hl = (ie0*L/alpha)**(1/m)

        xr_t = list()
        xf_t = list()
        t_rec = list()
        xf = 0.0

        while t < balance_t or xf < L:
            t = t + dt
            hl = self.iteration(L, t, hl)
            x_sta = self.x_staMove(t)

            diff = np.zeros(len(x_grd))
            balance2 = (ie1*x_grd/alpha)**(1/m)
            xf = alpha/ie1*(ie1*t)**m
            print(t, xf)
            for i in range(1, len(h_on_x)-1):
                h_on_x[i] = self.h_on_x(x_grd[i], xf, x_sta, t)

            try:
                if x_sta >= L:
                    h_on_x[-1] = hl + ie1*t
                else:
                    h_on_x[-1] = (ie0*x_s/alpha)**(1/m) + ie1*t
            except(TypeError):
                break

            sto = self.storage(h_on_x)

            if str(t)[-1] == '0' and self.plot:
                self.gen_depth_graph(h_on_x, t, xf)

            xr_list = list()
            for i in range(0, len(x_grd)):
                diff[i] = h_on_x[i] - balance2[i]

            for i in range(0, len(x_grd)):
                if diff[i] == 0:
                    xr_list.append(x_grd[i])
            xr = max(xr_list)
            xr_t.append(xr)
            xf_t.append(xf)
            t_rec.append(t)

            out_h.append([t, h_on_x[-1]])
            out_q.append([t, alpha*h_on_x[-1]**m])
            out_s.append([t, sto])

        out_h = np.array(out_h)
        out_q = np.array(out_q)
        out_s = np.array(out_s)

        self.out_q = out_q
        self.out_s = out_s
        self.t = t_rec

        plt.figure
        plt.plot(t_rec, xr_t, label='$x_r$')
        plt.plot(t_rec, xf_t, label='$x_f$')
        plt.legend(loc=2)
        plt.xlabel('t')
        plt.ylabel('x (m)')
        plt.title('Unsteady ' + str(ie0*3600*1000) + ' to ' + str(ie1*3600*1000)
                  + '$x_f$ and  $x_r$')
        plt.savefig(folder + '/xr_and_xf.png')
        plt.close()

        plt.figure
        plt.plot(out_q[:, 0], out_q[:, 1])
        plt.title('Outflow rate')
        plt.xlabel('t')
        plt.ylabel('q $(m^3)/s$')
        plt.title('Unsteady ' + str(ie0*3600*1000) + ' to ' + str(ie1*3600*1000)
                  + 'catchment outflow rate')
        plt.savefig(folder + '/outflow.png')
        plt.close()

    def run(self):
        self.x_grds()
        self.initial_cond()
        self.calc_process()

    def gen_depth_graph(self, h_on_x, t, xf):
        x_grd = self.grd
        ie0 = self.ie0
        ie1 = self.ie1
        L = self.L
        alpha = self.alpha
        m = self.m
        x_s = self.x_s
        folder = self.folder

        balance1 = np.zeros(len(x_grd))
        for i in range(0, len(x_grd)):
            if x_grd[i] <= x_s:
                balance1[i] = (ie0*x_grd[i]/alpha)**(1/m)
            else:
                balance1[i] = (ie0*x_s/alpha)**(1/m)
        balance2 = (ie1*x_grd/alpha)**(1/m)

        hl_lim = max((L*ie0/alpha)**(1/m), (L*ie1/alpha)**(1/m))

        ie0 = ie0*3600*10**3
        ie1 = ie1*3600*10**3

        plt.figure
        plt.plot(x_grd, h_on_x)
        plt.plot(x_grd, balance1, 'r--', linewidth=1.0,
                 label='$i_e=$'+str(ie0)+' $x_s=$'+str(x_s))
        plt.plot(x_grd, balance2, 'g--', linewidth=1.0,
                 label='$i_e=$'+str(ie1)+'  Steady-state')
        plt.axvline(x=xf, linestyle='dashed', color='k')
        plt.ylim(0.0, hl_lim)
        plt.title(str(t))
        plt.legend(loc=2)
        filename = str(t)+'.png'
        plt.savefig(folder + '/depth_change/' + filename.zfill(10))
        plt.close()

    def x_staMove(self, t):
        ie0 = self.ie0
        ie1 = self.ie1
        x_s = self.x_s
        alpha = self.alpha
        m = self.m

        h_s = (ie0*x_s/alpha)**(1/m)
        x_sta = x_s + alpha/ie1*(h_s + ie1*t)**m - ie0/ie1*x_s

        return x_sta

    def plot_ST(self, LT='-r', LW=1.0):
        ie0 = self.ie0
        ie1 = self.ie1
        x_s = self.x_s
        t = self.t
        S = self.out_s
        plt.plot(t, S[:, 1], LT, label=str(ie0) + 'unsteady to' + str(ie1) +
                 ' $x_s=$'+str(x_s))

    def plot_QT(self, LT='-r', LW=1.0):
        ie0 = self.ie0
        ie1 = self.ie1
        x_s = self.x_s
        t = self.t
        Q = self.out_q
        plt.plot(t, Q[:, 1], LT, label=str(ie0) + 'unsteady to' + str(ie1) +
                 ' $x_s=$'+str(x_s))


class progression:
    def __init__(self, ie, L, x_space, dt, n, slope, **kwargs):
        self.alpha = 1.0/n*sqrt(slope)
        self.m = 5.0/3
        self.L = L  # Length of slope
        self.x_space = x_space  # Number of grids
        self.dt = dt

        if 'folder' in kwargs:
            self.folder = kwargs['folder']
            if not os.path.isdir(kwargs['folder']):
                os.mkdir(kwargs['folder'])
                os.mkdir(kwargs['folder']+'/depth_change')
        else:
            folder = main_folder + '/calc_result/' + str(ie) + '_pro'
            if not os.path.isdir(folder):
                os.mkdir(folder)
                os.mkdir(folder + '/depth_change')
            self.folder = folder

        if 'plot' in kwargs:
            self.plot = kwargs['plot']
        else:
            self.plot = False

        self.ie = float(ie)/3600*10**-3

    def run(self):
        self.calc_process()

    def calc_process(self):
        ie = self.ie
        m = self.m
        L = self.L
        grds = self.x_space
        dt = self.dt
        alpha = self.alpha

        dx = L/(grds - 1)
        balance_t = (L*ie**(1-m)/alpha)**(1/m)
        self.max_t = balance_t
        self.max_q = L*ie
        x = np.arange(0, L+dx, dx)

        h_on_x = np.zeros(grds)

        out_q = list()
        out_s = list()
        depth_rec = list()

        t = 0.0
        while t < balance_t:
            xs = alpha*ie**(m-1)*t**m  # position of balance front
            for i in range(0, grds):
                if x[i] <= xs:
                    h_on_x[i] = (x[i]*ie/alpha)**(1/m)
                else:
                    h_on_x[i] = ie*t

            sto = m/(m+1)*(ie*xs**(m+1)/alpha)**(1/m) + (L - xs)*ie*t
            q = alpha*h_on_x[-1]**m

            out_q.append([t, q])
            out_s.append([t, sto])
            depth_rec.append([t, copy.copy(h_on_x)])
            t = t + dt

        out_q = np.array(out_q)
        out_s = np.array(out_s)
        self.out_q = out_q
        self.out_s = out_s

        if self.plot:
            self.plot_depth(depth_rec)

    def plot_QT(self, LT='-b'):
        ie = self.ie

        out_q = self.out_q
        plt.plot(out_q[:, 0], out_q[:, 1], LT,
                 label='$i_e$ = ' + str(ie) + ' Rising')

    def plot_ST(self, LT='-b'):
        ie = self.ie
        out_s = self.out_s

        plt.plot(out_s[:, 0], out_s[:, 1], LT,
                 label='$i_e$ = ' + str(ie) + ' Rising')

    def plot_depth(self, depth_rec):
        for i in range(0, len(depth_rec)):
            if str(depth_rec[i][0])[-1] == '0':
                self.gen_depth_graph(depth_rec[i][1], depth_rec[i][0])

    def gen_depth_graph(self, h_on_x, t):
        L = self.L
        grds = self.x_space
        out_q = self.out_q
        max_h = (self.ie*L/self.alpha)**(1./self.m)
        dx = L/(grds-1)

        x_grd = np.zeros(grds)
        for i in range(1, grds):
            x_grd[i] = dx*i

        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))
        axes[0].plot(x_grd, h_on_x, 'r')
        axes[0].set_xlabel('m')
        axes[0].set_ylabel('m')
        axes[0].set_ylim([0, 2*max_h])
        axes[0].set_title('Depth Along Slope')

        c_q = np.interp(t, out_q[:, 0], out_q[:, 1])
        axes[1].plot(out_q[:, 0], out_q[:, 1], '-b')
        axes[1].plot(t, c_q, 'bo')
        axes[1].set_title('Flowrate - Time')
        axes[1].set_xlabel('Time (s)')
        axes[1].set_ylabel('Flowrate ($m^3/s$)')
        axes[1].set_xlim([0, max(out_q[:, 0])])

        fig.suptitle('$i_e=100.0$ Rising - t = ' + str(t))
        filename = str(t)+'.png'
        fig.savefig(self.folder + '/depth_change/' + filename.zfill(10))
        plt.close()


class Recession:
    def __init__(self, ie, L, x_space, dt, n, slope, **kwargs):
        self.L = L
        self.grds = x_space
        self.dt = dt
        self.alpha = 1./n*sqrt(slope)
        self.m = 5.0/3

        if 'folder' in kwargs:
            self.folder = kwargs['folder']
            if not os.path.isdir(kwargs['folder']):
                os.mkdir(kwargs['folder'])
                os.mkdir(kwargs['folder']+'/depth_change')
        else:
            folder = main_folder + '/calc_result/' + str(ie) + '_res'
            if not os.path.isdir(folder):
                os.mkdir(folder)
                os.mkdir(folder + '/depth_change')
            self.folder = folder

        if 'plot' in kwargs:
            self.plot = kwargs['plot']
        else:
            self.plot = False

        self.ie = float(ie)/3600*10**-3
        self.initial_state()

    def initial_state(self):
        alpha = self.alpha
        L = self.L
        dx = L/(self.grds-1)
        x = np.zeros(self.grds)
        for i in range(1, self.grds):
            x[i] = i*dx

        h_on_x = np.zeros(self.grds)
        for i in range(0, self.grds):
            h_on_x[i] = (self.ie*x[i]/alpha)**1/self.m

        self.init_h = h_on_x
        self.x = x

    def simpson(self, h_on_x):
        x_space = self.grds
        n = x_space - 1
        L = self.L
        dx = L/n

        mult = np.zeros(x_space)
        mult[0] = 1.0
        mult[-1] = 1.0
        #  Multiplier of simpson's 1/2 method
        for i in range(1, x_space-1):
            if i % 2 == 0:
                mult[i] = 2.0
            elif (i-1) % 2 == 0:
                mult[i] = 4.0

        return dx/3.*np.dot(mult, h_on_x)

    def run(self):
        self.calc_process()

    def calc_process(self):
        h_on_x = self.init_h
        ie = self.ie
        m = self.m
        alpha = self.alpha
        L = self.L
        grds = self.grds
        x = self.x
        dt = self.dt

        q0 = ie*L
        q = q0
        S0 = (ie/alpha)**(1/m)*m/(m+1)*L**((m+1)/m)
        t = 0.

        out_q = list()
        out_s = list()
        depth_rec = list()

        out_q.append([t, q0])
        out_s.append([t, S0])
        depth_rec.append([t, copy.copy(h_on_x)])

        while q > q0*10**-3:
            t = t + dt
            for i in range(1, grds):
                h_on_x[i] = (ie*self.iteration(x[i], t, x[i])/alpha)**(1/m)
            q = alpha*h_on_x[-1]**m
            S = self.simpson(h_on_x)
            out_q.append([t, q])
            out_s.append([t, S])
            depth_rec.append([t, copy.copy(h_on_x)])
            if str(t)[-1] == '0':
                print(q, t, q/q0)

        self.out_q = np.array(out_q)
        self.out_s = np.array(out_s)
        if self.plot:
            for i in range(0, len(depth_rec)):
                if str(depth_rec[i][0])[-1] == '0':
                    self.gen_depth_graph(depth_rec[i][1], depth_rec[i][0])

    def iteration(self, xa, t, ix):
        """
        Input assigned position and time , backtrace to its starting position
        and depth at rainfall intensity shift."""
        """ ih is the initial guess of starting position depth."""
        ie = self.ie
        alpha = self.alpha
        m = self.m
        err = 100.0
        last_err = err + 1.0
        self.counter = 0

        while err > 10**-6:
            #  Newton method iteration.
            try:
                f = xa - ix - alpha*m*(ie*ix/alpha)**((m-1)/m)*t
                df = -1 + (1-m)*alpha**(1/m)*ie**((m-1)/m)*ix**(-1/m)

                rx = self.iteration_counter(ix, f, df)
            except(ValueError):
                rx = self.bisec_iter(xa, t, 10**-6)
                return rx

            err = fabs(rx - ix)
            """
            Newton method easily goes divergent if df is too small, and not
            getting accurate value, thus if error between current and last
            iteration expands, use bisection method instead.
            """
            if last_err > err:
                ix = rx
                self.counter = self.counter + 1
            else:
                rx = self.bisec_iter(xa, t, 10**-6)
                return rx

        return rx

    def iteration_counter(self, ix, f, df):
        counter = self.counter
        if counter < 30:
            rx = ix - f/df
        elif counter >= 30:
            rx = ix - f/df*0.5
        elif counter >= 50:
            rx = ix - f/df*0.1
        elif counter >= 200:
            rx = ix - f/df*0.01
        elif counter >= 500:
            rx = ix - f/df*0.001

        return rx

    def test_bi(self, xa, t, rx, ix):
        #  Function of bisection method, test if solution between rh and ih.
        alpha = self.alpha
        m = self.m
        ie = self.ie

        f1 = (xa - rx - alpha*m*(ie*rx/alpha)**((m-1)/m)*t)
        f2 = (xa - ix - alpha*m*(ie*ix/alpha)**((m-1)/m)*t)
        return f1*f2

    def bisec(self, xa, t, ix, rx):
        #  Bisection method
        f1f2 = self.test_bi(xa, t, rx, ix)

        if f1f2 < 0:
            n_rx = (rx + ix)*0.5
            f1f2_2 = self.test_bi(xa, t, n_rx, ix)
            if f1f2_2 < 0:
                return [n_rx, ix]
            elif f1f2_2 > 0:
                return [n_rx, rx]
        elif f1f2 > 0:
            print('Error, no solution available.')
            time.sleep(10)
            raise(ValueError)

    def bisec_iter(self, xa, t, allow_err):
        #  Bisection method iteration.
        counter = self.counter

        ix = 0.
        rx = xa
        err = 100.0
        while err > allow_err:
            try:
                [rx, ix] = self.bisec(xa, t, rx, ix)
            except(ValueError):
                return None
            err = fabs(rx - ix)
            counter = counter + 1
            self.counter = counter

            if self.counter > 500:
                print(rx, ix)

            if self.counter > 1000:
                print(rx, ix)
                time.sleep(300)

        return rx

    def plot_QT(self, LT='-b'):
        ie = self.ie

        out_q = self.out_q
        plt.plot(out_q[:, 0], out_q[:, 1], LT,
                 label='$i_e$ = ' + str(ie) + ' Falling')

    def plot_ST(self, LT='-b'):
        ie = self.ie
        out_s = self.out_s

        plt.plot(out_s[:, 0], out_s[:, 1], LT,
                 label='$i_e$ = ' + str(ie) + ' Falling')

    def gen_depth_graph(self, h_on_x, t):
        L = self.L
        grds = self.grds
        out_q = self.out_q
        max_h = (self.ie*L/self.alpha)**(1./self.m)
        dx = L/(grds-1)

        x_grd = np.zeros(grds)
        for i in range(1, grds):
            x_grd[i] = dx*i

        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))
        axes[0].plot(x_grd, h_on_x, 'r')
        axes[0].set_xlabel('m')
        axes[0].set_ylabel('m')
        axes[0].set_ylim([0, 2*max_h])
        axes[0].set_title('Depth Along Slope')

        c_q = np.interp(t, out_q[:, 0], out_q[:, 1])
        axes[1].plot(out_q[:, 0], out_q[:, 1], '-b')
        axes[1].plot(t, c_q, 'bo')
        axes[1].set_title('Flowrate - Time')
        axes[1].set_xlabel('Time (s)')
        axes[1].set_ylabel('Flowrate ($m^3/s$)')
        axes[1].set_xlim([0, max(out_q[:, 0])])

        fig.suptitle('$i_e=100.0$ Falling - t = ' + str(t))
        filename = str(t)+'.png'
        fig.savefig(self.folder + '/depth_change/' + filename.zfill(10))
        plt.close()
