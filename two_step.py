
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from parameters import main_folder
import os
import os.path


class two_step_overland:

    def __init__(self, ie0, ie1, ie2, t_dist, x_space, dt, n, slope, **kwargs):
        self.ie0 = float(ie0)/3600*10**-3
        self.ie1 = float(ie1)/3600*10**-3
        self.ie2 = float(ie2)/3600*10**-3
        self.L = float(t_dist)
        self.grds = x_space
        self.dt = float(dt)
        self.alpha = 1.0/n*sqrt(slope)
        self.n = n
        self.slope = slope
        self.m = 5.0/3

        self.grids()

        if 'folder' in kwargs:
            self.folder = kwargs['folder']
            if not os.path.isdir(kwargs['folder']):
                folder_str = 'mkdir -p ' + kwargs['folder']
                os.system(folder_str)
                os.system(folder_str + '/depth_change')
        else:
            folder = ('calc_result' + '/'
                      + str(ie0) + 'to' + str(ie1) + 'to' + str(ie2))
            if not os.path.isdir(folder):
                os.mkdir(os.getcwd() + '/' + folder)
                os.mkdir(os.getcwd() + '/' + folder + '/depth_change')
            self.folder = folder

        self.initial_state()

        if 'plot' in kwargs:
            if kwargs['plot']:
                self.plot = True
            else:
                self.plot = False
        else:
            self.plot = False

        self.steady_sq()

    def ie1_duration(self, t1):
        self.t1 = t1
        self.state_at_t1()
        self.ie1_xs()
        self.ie1_term_time()
        self.xs_to_L()

    def ie1_xs(self):
        ie1 = self.ie1
        m = self.m
        alpha = self.alpha
        t1 = self.t1

        xs = alpha/ie1*(ie1*t1)**m
        self.xs1 = xs

    def storage(self, h_on_x):
        L = self.L
        grds = self.grds
        dx = L/(grds-1)

        Sto = 0.0
        for i in range(1, len(h_on_x)):
            Sto = Sto + (h_on_x[i] + h_on_x[i-1])*0.5*dx

        return Sto

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

    def grids(self):
        grds = self.grds  # Number of grids
        L = self.L  # Total length
        dx = L/(grds-1)

        x = np.arange(0.0, L + dx, dx)
        self.x = x

    def initial_state(self):
        ie0 = self.ie0
        grds = self.grds
        alpha = self.alpha
        init_h = np.zeros(grds)
        x = self.x
        m = self.m

        for i in range(0, grds):
            init_h[i] = (x[i]*ie0/alpha)**(1/m)

        self.init_h = init_h

    def ie1_term_time(self):
        #  Steady time required in ie1 rainfall
        L = self.L
        ie1 = self.ie1
        alpha = self.alpha
        m = self.m
        print (L*ie1**(1-m)/alpha)**(1/m)
        self.t_ie1 = (L*ie1**(1-m)/alpha)**(1/m)

    def plot_depth(self, depth_rec):
        x = self.x
        folder = self.folder
        max_h = (max(self.q_curve[:, 1])/self.alpha)**(1/self.m)

        for recs in depth_rec:
            t = recs[0]
            if int(t) % 5 == 0:
                print t
                fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))
                axes[0].plot(x, recs[1], '-r')
                axes[0].set_xlabel('m')
                axes[0].set_ylabel('m')
                axes[0].set_title('Depth Along Slope')
                axes[0].set_ylim([0, max_h])

                c_q = np.interp(t, self.q_curve[:, 0], self.q_curve[:, 1])
                axes[1].plot(self.q_curve[:, 0], self.q_curve[:, 1], '-b')
                axes[1].plot(t, c_q, 'bo')
                axes[1].set_xlabel('Time (s)')
                axes[1].set_ylabel('Flowrate ($m^3/s$)')
                axes[1].set_xlim([0, max(self.q_curve[:, 0])])
                axes[1].set_title('Flowrate - Time')
                fig.suptitle('t = ' + str(t))
                filename = str(t)+'.png'
                fig.savefig(folder + '/depth_change/' + filename.zfill(10))
                plt.close()

    def state_at_t1(self):
        t1 = self.t1
        x = self.x
        grds = self.grds
        ie0 = self.ie0
        ie1 = self.ie1
        alpha = self.alpha
        init_h = self.init_h
        m = self.m

        xz_pos = np.zeros(grds)
        h_on_x = np.zeros(grds)

        xf = alpha*ie1**(m-1)*t1**m
        for i in range(0, grds):
            xz_pos[i] = (x[i] + alpha/ie1*(init_h[i] + ie1*t1)**m -
                         alpha/ie1*init_h[i]**m)
        for j in range(1, grds):
            c_x = np.interp(x[j], xz_pos, x)
            if x[j] > xf:
                h_on_x[j] = (c_x*ie0/alpha)**(1/m) + ie1*t1
            else:
                h_on_x[j] = (x[j]*ie1/alpha)**(1/m)
        self.h_on_t1 = h_on_x

    def step1(self, h_on_x, t):
        alpha = self.alpha
        ie0 = self.ie0
        ie1 = self.ie1
        m = self.m
        init_h = self.init_h
        x = self.x
        grds = self.grds
        xz_pos = np.zeros(grds)

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

    def step2(self, h_on_x, t):
        grds = self.grds
        x = self.x
        ie2 = self.ie2
        alpha = self.alpha
        m = self.m
        t1 = self.t1
        h_on_t1 = self.h_on_t1

        xz2_pos = np.zeros(grds)

        xf = alpha*ie2**(m-1)*(t - t1)**m
        for i in range(0, grds):
            xz2_pos[i] = (x[i] + alpha/ie2*(h_on_t1[i] + ie2*(t - t1))**m
                          - alpha/ie2*h_on_t1[i]**m)
        for j in range(1, grds):
            c_x = np.interp(x[j], xz2_pos, x)

            if x[j] > xf:
                h_on_x[j] = np.interp(c_x, x, h_on_t1) + ie2*(t - t1)
            else:
                h_on_x[j] = (x[j]*ie2/alpha)**(1/m)

        return h_on_x

    def run(self):
        alpha = self.alpha
        ie2 = self.ie2
        t1 = self.t1
        L = self.L
        m = self.m
        grds = self.grds
        dt = self.dt
        depth_rec = list()

        term_time = t1 + (L*ie2**(1-m)/alpha)**(1/m)
        t = 0.0

        h_on_x = np.zeros(grds)
        q = list()  # Outflow curve
        s = list()  # Storage curve

        while t < term_time:
            if t <= t1:
                h_on_x = self.step1(h_on_x, t)
            elif t > t1 and t < term_time:
                h_on_x = self.step2(h_on_x, t)

            q.append([t, alpha*h_on_x[-1]**m])
            s.append([t, self.simpson(h_on_x)])

            depth_rec.append([t, h_on_x.copy()])

            t = t + dt

        q = np.array(q)
        s = np.array(s)

        self.q_curve = q
        self.s_curve = s

        if self.plot:
            self.plot_depth(depth_rec)

    def plot_q(self, LT='-r'):
        plt.plot(self.q_curve[:, 0], self.q_curve[:, 1], LT)
        plt.xlabel('t (s)')
        plt.ylabel('flowrate ($m^3/s$)')
        title = (str(self.ie0*1000*3600) + ' to ' + str(self.ie1*1000*3600)
                 + 'to ' + str(self.ie0*1000*3600) + ' $t_1$ = ' + str(self.t1)
                 + '\n' + 'outflow rate')
        plt.title(title)

    def plot_s(self, LT='-k'):
        plt.plot(self.s_curve[:, 0], self.s_curve[:, 1], LT)
        plt.xlabel('t (s)')
        plt.ylabel('storage ($m^3$)')
        title = (str(self.ie0*1000*3600) + ' to ' + str(self.ie1*1000*3600)
                 + 'to ' + str(self.ie0*1000*3600) + ' $t_1$ = ' + str(self.t1)
                 + '\n' + 'storage')
        plt.title(title)

    def ie1_to_ie2(self):
        ie1 = self.ie1
        ie2 = self.ie2
        alpha = self.alpha
        m = self.m
        L = self.L
        grds = self.grds
        x = self.x
        dt = self.dt

        term_time = (L*ie2**(1-m)/alpha)**(1/m)
        #  initial state
        init_h = np.zeros(grds)
        for i in range(0, grds):
            init_h[i] = (x[i]*ie1/alpha)**(1/m)

        t = 0.0
        h_on_x = np.zeros(grds)
        xz_pos = np.zeros(grds)

        q_curve = list()
        s_curve = list()
        while t < term_time:
            for j in range(1, grds):
                xz_pos[j] = (x[j] + alpha/ie2*(init_h[j] + ie2*t)**m
                             - alpha/ie2*init_h[j]**m)

            xf = alpha/ie2*(ie2*t)**m
            for z in range(1, grds):
                c_x = np.interp(x[z], xz_pos, x)
                if x[z] > xf:
                    h_on_x[z] = (c_x*ie1/alpha)**(1/m) + ie2*t
                else:
                    h_on_x[z] = (x[z]*ie2/alpha)**(1/m)
            q_curve.append([t, alpha*h_on_x[-1]**m])
            s_curve.append([t, self.simpson(h_on_x)])

            t = t + dt

        self.ie1_to_ie2_q = np.array(q_curve)
        self.ie1_to_ie2_s = np.array(s_curve)

    def plot_ie1_to_ie2_q(self, LT='-b'):
        ie1_to_ie2_q = self.ie1_to_ie2_q
        ie1_to_ie2_q[:, 0] = ie1_to_ie2_q[:, 0] + self.t1

        plt.plot(ie1_to_ie2_q[:, 0], ie1_to_ie2_q[:, 1], LT)

    def plot_ie1_to_ie2_Normq(self, LT='-b'):
        ie1_to_ie2_q = self.ie1_to_ie2_q
        ie1_to_ie2_q[:, 0] = ie1_to_ie2_q[:, 0]

        plt.plot(ie1_to_ie2_q[:, 0], ie1_to_ie2_q[:, 1], LT)

    def norm(self):
        L = self.L
        ie2 = self.ie2
        alpha = self.alpha
        m = self.m
        t1 = self.t1

        ie2_q = L*ie2
        ie2_s = (m/(m+1))*L**((m+1)/m)*(ie2/alpha)**(1/m)
        ie2_time = (L*ie2**(1-m)/alpha)**(1/m)

        self.q_curve[:, 1] = self.q_curve[:, 1]/ie2_q
        self.q_curve[:, 0] = (self.q_curve[:, 0] - t1)/ie2_time
        self.s_curve[:, 1] = self.s_curve[:, 1]/ie2_s
        self.s_curve[:, 0] = (self.s_curve[:, 0] - t1)/ie2_time

        self.ie1_to_ie2_q[:, 1] = self.ie1_to_ie2_q[:, 1]/ie2_q
        self.ie1_to_ie2_q[:, 0] = self.ie1_to_ie2_q[:, 0]/ie2_time

        self.ie1_to_ie2_s[:, 1] = self.ie1_to_ie2_s[:, 1]/ie2_s
        self.ie1_to_ie2_s[:, 0] = self.ie1_to_ie2_s[:, 0]/ie2_time

    def xs_to_L(self):
        ie1 = self.ie1
        ie2 = self.ie2
        t1 = self.t1
        alpha = self.alpha
        m = self.m
        L = self.L

        xs1 = alpha/ie1*(ie1*t1)**m

        term_time = (L*ie2**(1-m)/alpha)**(1/m)

        t = np.arange(0.0, term_time, 1.0)
        xz_pos = np.zeros(len(t))

        for i in range(0, len(t)):
            xz_pos[i] = (xs1 + alpha/ie2*(ie1*t1 + ie2*t[i])**m
                         - alpha/ie2*(ie1*t1)**m)

        t_xs = np.interp(L, xz_pos, t)  # Time required for xs1 to reach L.
        print t_xs
        q1 = alpha*(ie1*t1 + ie2*(t_xs))**m
        q2 = alpha*(ie1*t1)**m

        print ('xs1 to reach L: ' + str(t_xs))
        print ('flowrate: ' + str(alpha*(ie1*t1 + ie2*(t_xs))**m))
        print ('flowate at t1: ' + str(alpha*(ie1*t1)**m))
        print ('z = ' + str(q2/q1))

        self.t_xs = t_xs
        self.t_term = term_time

    def transition(self):
        ie1 = self.ie1
        ie2 = self.ie2
        h_on_t1 = self.h_on_t1
        alpha = self.alpha
        m = self.m
        t_xs = self.t_xs
        grds = self.grds
        x = self.x
        dt = self.dt

        q_curve = list()
        s_curve = list()
        t = 0.0
        h_on_x = np.zeros(grds)
        xz_pos = np.zeros(grds)
        while t < t_xs:
            for j in range(1, grds):
                xz_pos[j] = (x[j] + alpha/ie2*(h_on_t1[j] + ie2*t)**m
                             - alpha/ie2*h_on_t1[j]**m)

            xf = alpha/ie2*(ie2*t)**m
            for z in range(1, grds):
                c_x = np.interp(x[z], xz_pos, x)
                if x[z] > xf:
                    h_on_x[z] = (c_x*ie1/alpha)**(1/m) + ie2*t
                else:
                    h_on_x[z] = (x[z]*ie2/alpha)**(1/m)
            q_curve.append([t, alpha*h_on_x[-1]**m])
            s_curve.append([t, self.simpson(h_on_x)])

            t = t + dt

        q_curve = np.array(q_curve)
        s_curve = np.array(s_curve)

        self.transition_q = q_curve
        self.transition_s = s_curve

        self.steady_intersection()

    def plot_transition_q(self, LT='-r'):
        plt.plot(self.transition_q[:, 0], self.transition_q[:, 1], LT)

    def plot_transition_s(self, LT='-r'):
        plt.plot(self.transition_s[:, 0], self.transition_s[:, 1], LT)

    def steady_sq(self):
        #  Generate storage and outflow under different rainfall intensity
        #  at steady-state.
        m = self.m
        alpha = self.alpha
        L = self.L

        ie = np.arange(0.0, 500.0, 0.5)
        ie = ie/1000.0/3600.0
        storage = np.zeros(len(ie))
        outflow = np.zeros(len(ie))

        for i in range(0, len(ie)):
            storage[i] = m/(m+1)*(ie[i]/alpha)**(1/m)*L**(1+1/m)
            outflow[i] = ie[i]*L

        steady_ref = np.zeros([len(ie), 3])
        steady_ref[:, 0] = ie
        steady_ref[:, 1] = storage
        steady_ref[:, 2] = outflow

        #  The output
        #  1st column of steady_ref : rainfall intensity in mm/hr
        #  2nd column of steady-ref : storage at steady-state
        #  3rd column of steady-ref : outflow at steady-state
        self.steady_ref = steady_ref

    def steady_intersection(self):
        steady_sq = self.steady_ref
        q_rec = self.q_curve
        s_rec = self.s_curve

        min_error = 1.0
        seq = 0
        counter = 0
        #  Find the intersection of hydrograph S-Q and steady-state S-Q
        #  Using q=f(s), the closest estimated q and hydrograph q, is the
        #  intersection point.
        for i in range(1, len(q_rec)-1000):
            q_est = np.interp(s_rec[i, 1], steady_sq[:, 1], steady_sq[:, 2])
            error = (q_rec[i, 1] - q_est)**2
            counter = counter + 1
            if error < min_error and q_rec[i, 0] > self.t1:
                seq = i
                min_error = error
                counter = 0
            if q_rec[i, 0] > self.t1 and counter > 1000:
                break

        #  Find the closest point to t1 (transition point of ie1 -> ie2)
        seq_t1 = 0
        min_t1_err = 100.0
        for i in range(0, len(q_rec)):
            err = (q_rec[i, 0] - self.t1)**2
            if err < min_t1_err:
                seq_t1 = i
                min_t1_err = err

        print seq, min_error
        self.q_AC = q_rec[seq_t1:seq]
        self.s_AC = s_rec[seq_t1:seq]
        self.trans_ie = np.interp(s_rec[seq, 1], steady_sq[:, 1],
                                  steady_sq[:, 0])*1000.0*3600.0

        gamma = q_rec[seq, 1]/np.interp(self.t1, q_rec[:, 0], q_rec[:, 1])
        print gamma
