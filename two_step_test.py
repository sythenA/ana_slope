
import two_step as ts
import matplotlib.pyplot as plt


aa = ts.two_step_overland(90.0, 100.0, 40.0, 100.0, 201, 0.5, 0.1, 0.01)
aa.ie1_duration(502.5)
ab = ts.two_step_overland(80.0, 200.0, 40.0, 100.0, 201, 0.5, 0.1, 0.01)
ab.ie1_duration(207.0)
ac = ts.two_step_overland(10.0, 90.0, 5.0, 100.0, 201, 0.5, 0.1, 0.01)
ac.ie1_duration(117.3)

aa.run()
ab.run()
ac.run()

aa.transition()
ab.transition()
ac.transition()

aq = aa.transition_q
bq = ab.transition_q
cq = ac.transition_q

aq[:, 1] = (aq[:, 1] - min(aq[:, 1]))/(1.0 - min(aq[:, 1]))
bq[:, 1] = (bq[:, 1] - min(bq[:, 1]))/(1.0 - min(bq[:, 1]))
cq[:, 1] = (cq[:, 1] - min(cq[:, 1]))/(1.0 - min(cq[:, 1]))

plt.figure
plt.plot(aq[:, 0], aq[:, 1], '-b')
plt.plot(bq[:, 0], bq[:, 1], '--c', lw=3.0)
plt.plot(cq[:, 0], cq[:, 1], '--m', lw=2.5)
plt.show()
