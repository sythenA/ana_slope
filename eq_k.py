
import ana_new as an
import matplotlib.pyplot as plt


def run(k):

    if k > 1.0:
        aa = an.trans_ana(10.0, 10.0*k, 100.0, 201, 0.5, 0.1, 0.01)
        aa.run()
        ab = an.trans_ana(30.0, 30.0*k, 100.0, 201, 0.5, 0.1, 0.01)
        ab.run()
        ac = an.trans_ana(60.0, 60.0*k, 100.0, 201, 0.5, 0.1, 0.01)
        ac.run()
    elif k < 1.0:
        aa = an.trans_ana(50.0, 50.0*k, 100.0, 201, 0.5, 0.1, 0.01)
        aa.run()
        ab = an.trans_ana(100.0, 100.0*k, 100.0, 201, 0.5, 0.1, 0.01)
        ab.run()
        ac = an.trans_ana(200.0, 200.0*k, 100.0, 201, 0.5, 0.1, 0.01)
        ac.run()

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))
    axes[0].plot(ac.out_q[:, 0]/max(ac.out_q[:, 0]),
                 ac.out_q[:, 1]/max(ac.out_q[:, 1]), linestyle='-',
                 color='orange', lw=8.0,
                 label=str(ac.ie0*3600*1000)+' to '+str(ac.ie1*3600*1000))
    axes[0].plot(ab.out_q[:, 0]/max(ab.out_q[:, 0]),
                 ab.out_q[:, 1]/max(ab.out_q[:, 1]), linestyle='-',
                 color='beige', lw=3.5,
                 label=str(ab.ie0*3600*1000)+' to '+str(ab.ie1*3600*1000))
    axes[0].plot(aa.out_q[:, 0]/max(aa.out_q[:, 0]),
                 aa.out_q[:, 1]/max(aa.out_q[:, 1]), linestyle='--',
                 color='steelblue',
                 label=str(aa.ie0*3600*1000)+' to '+str(aa.ie1*3600*1000),
                 lw=2.0)
    if k > 1:
        axes[0].legend(loc=2)
    elif k < 1:
        axes[0].legend(loc=1)
    axes[0].set_xlabel('$t_{{}*{}}$')
    axes[0].set_ylabel('$q_{{}*{}}$')
    axes[0].set_title('Outflow rate - Time')

    axes[1].plot(ac.out_s[:, 0]/max(ac.out_s[:, 0]),
                 ac.out_s[:, 1]/max(ac.out_s[:, 1]), linestyle='-',
                 color='orange', lw=8.0,
                 label=str(ac.ie0*3600*1000)+' to '+str(ac.ie1*3600*1000))
    axes[1].plot(ab.out_s[:, 0]/max(ab.out_s[:, 0]),
                 ab.out_s[:, 1]/max(ab.out_s[:, 1]), linestyle='-',
                 color='beige', lw=3.5,
                 label=str(ab.ie0*3600*1000)+' to '+str(ab.ie1*3600*1000))
    axes[1].plot(aa.out_s[:, 0]/max(aa.out_s[:, 0]),
                 aa.out_s[:, 1]/max(aa.out_s[:, 1]), linestyle='--',
                 color='steelblue',
                 label=str(aa.ie0*3600*1000)+' to '+str(aa.ie1*3600*1000),
                 lw=2.0)
    if k > 1:
        axes[1].legend(loc=2)
    elif k < 1:
        axes[1].legend(loc=1)
    axes[1].set_xlabel('$t_{{}*{}}$')
    axes[1].set_ylabel('$S_{{}*{}}$')
    axes[1].set_title('Storage - Time')
    fig.suptitle('k = '+str(k))
    plt.savefig('graph/eq_k_'+str(k)+'.png')
    plt.close()
