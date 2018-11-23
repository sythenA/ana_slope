
# coding: utf-8

# In[3]:



from SQsearch import RisingLimb, FallingLimb
from farPoint import SWRCurves
from ana_new import trans_ana
from two_step import two_step_overland
from steady_fit import fit_steady
import numpy as np
import matplotlib.pyplot as plt

Rise = RisingLimb(rec_name='ks1', folder='Steady_fit')
Fall = FallingLimb(rec_name='ks1', folder='Steady_fit')

SWR = SWRCurves(100.0, 201, 0.5, 0.1, 0.01, 2, 4, useFile='ks1_SWRrec.pick')
SWR.fitPT()

std = fit_steady(save_name='ks1', folder='Steady_fit')
std.load()

ie0 = 100.0
ie1 = 25.0
ie2 = 70.0

a = trans_ana(100.0, 25.0, 100.0, 201, 0.5, 0.1, 0.01)
a.run()

b = trans_ana(25.0, 70.0, 100.0, 201, 0.5, 0.1, 0.01)
b.run()

timeOfChange = 650.0

c = two_step_overland(100.0, 25.0, 70.0, 100.0, 201, 0.5, 0.1, 0.01)
c.ie1_duration(timeOfChange)
c.run()

S = np.interp(timeOfChange, a.out_s[:, 0], a.out_s[:, 1])
Q = np.interp(timeOfChange, a.out_q[:, 0], a.out_q[:, 1])
"""
maxS = std.ie_to_storage(70.0/1000/3600)
maxQ = std.ie_to_outflow(70.0/1000/3600)
curve, T = Rise.getCurve(70.0/25.0, maxS, maxQ)"""

curve = np.array([zip(b.out_q[:, 1], b.out_s[:, 1])])[0]

gQ, gS = SWR.endFit2(S, Q, curve)
print gQ/Q

swc, t = SWR.interpCurve(gQ/Q)
swc[:, 0] = swc[:, 0]*Q
swc[:, 1] = swc[:, 1]*S

plt.figure
plt.plot(a.out_q[:, 1], a.out_s[:, 1], '--k')
plt.plot(b.out_q[:, 1], b.out_s[:, 1], '--k')
plt.plot(c.q_curve[:, 1], c.s_curve[:, 1], 'r')
plt.plot(np.array(SWR.header)*Q, np.array(SWR.Sheader)*S)
plt.plot(swc[:, 0], swc[:, 1])
plt.xlabel('q($m^2/s$)')
plt.ylabel('s(m^2)')
plt.title('100 - 25 - 70')
plt.savefig('SWCtest1.png')

