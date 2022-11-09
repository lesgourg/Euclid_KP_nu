import numpy as np
import matplotlib.pyplot as plt
from math import sqrt,log

# plot for non-diagonal coefficient of Fisher matrix:

if True:

    xmax=40

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(6, 6))

    ######################

    ymin=-0.1
    ymax=0.4

    ax1.set_ylim(ymin,ymax)

    ######################

    ax1.set_ylabel(r'error on $F_{w_0\,w_a}$ (%)')
    ax1.set_xlim(0.,0.4)

    DP = np.loadtxt('photometric/pessimistic_DP/d2L_dw0_dwa.dat')
    HP = np.loadtxt('photometric/pessimistic_HP/d2L_dw0_dwa.dat')
    fit = np.polynomial.polynomial.Polynomial.fit(HP[:,0], HP[:,1], 2, domain=[0.1,1.])
    asymptot=fit(0)

    ax1.plot(HP[:xmax,0],100.*(fit(HP[:xmax,0])/asymptot-1.),linestyle='-',color='k')
    ax1.plot(DP[:xmax,0],100.*(DP[:xmax,1]/asymptot-1.),label="photo/pess DP",linestyle='--',marker='+')
    ax1.plot(HP[:xmax,0],100.*(HP[:xmax,1]/asymptot-1.),label="photo/pess HP",linestyle='--',marker='+')
    ax1.legend()

    ######################

    DP = np.loadtxt('photometric/optimistic_DP/d2L_dw0_dwa.dat')
    HP = np.loadtxt('photometric/optimistic_HP/d2L_dw0_dwa.dat')
    fit = np.polynomial.polynomial.Polynomial.fit(HP[:,0], HP[:,1], 2, domain=[0.1,1])
    asymptot=fit(0)

    ax2.plot(HP[:xmax,0],100.*(fit(HP[:xmax,0])/asymptot-1.),linestyle='-',color='k')
    ax2.plot(DP[:xmax,0],100.*(DP[:xmax,1]/asymptot-1.),label="photo/opt DP",linestyle='--', marker='+')
    ax2.plot(HP[:xmax,0],100.*(HP[:xmax,1]/asymptot-1.),label="photo/opt HP",linestyle='--', marker='+')
    ax2.legend()

    ######################

    ymin=-0.1
    ymax=0.4

    ax3.set_ylim(ymin,ymax)

    ######################

    ax3.set_xlabel(r'stepsize/error')
    ax3.set_ylabel(r'error on $F_{w_0\,w_a}$ (%)')

    DP = np.loadtxt('spectroscopic/pessimistic_DP/d2L_dw0_dwa.dat')
    HP = np.loadtxt('spectroscopic/optimistic_HP/d2L_dw0_dwa.dat')
    fit = np.polynomial.polynomial.Polynomial.fit(HP[:,0], HP[:,1], 2, domain=[0.1,1])
    asymptot=fit(0)

    ax3.plot(HP[:xmax,0],100.*(fit(HP[:xmax,0])/asymptot-1.),linestyle='-',color='k')
    ax3.plot(DP[:xmax,0],100.*(DP[:xmax,1]/asymptot-1.),label="spectro/pess DP",linestyle='--', marker='+')
    ax3.plot(HP[:xmax,0],100.*(HP[:xmax,1]/asymptot-1.),label="spectro/pess HP",linestyle='--', marker='+')
    ax3.legend()

    ######################

    ax4.set_xlabel(r'stepsize/error')

    DP = np.loadtxt('spectroscopic/optimistic_DP/d2L_dw0_dwa.dat')
    HP = np.loadtxt('spectroscopic/optimistic_HP/d2L_dw0_dwa.dat')
    fit = np.polynomial.polynomial.Polynomial.fit(HP[:,0], HP[:,1], 2, domain=[0.1,1])
    asymptot=fit(0)

    ax4.plot(HP[:xmax,0],100.*(fit(HP[:xmax,0])/asymptot-1.),linestyle='-',color='k')
    ax4.plot(DP[:xmax,0],100.*(DP[:xmax,1]/asymptot-1.),label="spectro/opt DP",linestyle='--', marker='+')
    ax4.plot(HP[:xmax,0],100.*(HP[:xmax,1]/asymptot-1.),label="spectro/opt HP",linestyle='--', marker='+')
    ax4.legend()

    ######################

    plt.savefig('test_stepsizes.pdf')

    ######################

# plot for diagonal coefficient of Fisher matrix:

if False:

    xmax=5

    fig2, ((ax5, ax6)) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(6, 3))

    #######################

    DP = np.loadtxt('spectroscopic/pessimistic_DP/d2L_d2w0.dat')
    HP = np.loadtxt('spectroscopic/pessimistic_HP/d2L_d2w0.dat')
    fit = np.polynomial.polynomial.Polynomial.fit(HP[:,0], HP[:,1], 2, domain=[0.1,5])
    asymptot=fit(0)

    ax5.plot(HP[:xmax,0],100.*(fit(HP[:xmax,0])/asymptot-1.),linestyle='-',color='k')
    ax5.plot(DP[:xmax,0],100.*(DP[:xmax,1]/asymptot-1.),label="spectro/pess $w_0$ DP",linestyle='--', marker='+')
    ax5.plot(HP[:xmax,0],100.*(HP[:xmax,1]/asymptot-1.),label="spectro/pess $w_0$ HP",linestyle='--', marker='+')
    ax5.legend()

    #######################

    DP = np.loadtxt('spectroscopic/pessimistic_DP/d2L_d2wa.dat')
    HP = np.loadtxt('spectroscopic/pessimistic_HP/d2L_d2wa.dat')
    fit = np.polynomial.polynomial.Polynomial.fit(HP[:,0], HP[:,1], 2, domain=[0.1,5])
    asymptot=fit(0)

    ax6.plot(HP[:xmax,0],100.*(fit(HP[:xmax,0])/asymptot-1.),linestyle='-',color='k')
    ax6.plot(DP[:xmax,0],100.*(DP[:xmax,1]/asymptot-1.),label="spectro/pess $w_0$ DP",linestyle='--', marker='+')
    ax6.plot(HP[:xmax,0],100.*(HP[:xmax,1]/asymptot-1.),label="spectro/pess $w_0$ HP",linestyle='--', marker='+')
    ax6.legend()

    ######################

    ######################

    plt.savefig('test_stepsizes_diag.pdf')
