import numpy as np
import matplotlib.pyplot as plt

def derivative(case,param,step):
    kk_plus = np.loadtxt('../../input_4_cast/output/%s/%s_pl_eps_1p0E-02/k_values_list.txt'%(case,param))
    pl_plus = np.loadtxt('../../input_4_cast/output/%s/%s_pl_eps_1p0E-02/Plin-zk.txt'%(case,param))
    pn_plus = np.loadtxt('../../input_4_cast/output/%s/%s_pl_eps_1p0E-02/Pnonlin-zk.txt'%(case,param))
    kk_minus = np.loadtxt('../../input_4_cast/output/%s/%s_mn_eps_1p0E-02/k_values_list.txt'%(case,param))
    pl_minus = np.loadtxt('../../input_4_cast/output/%s/%s_mn_eps_1p0E-02/Plin-zk.txt'%(case,param))
    pn_minus = np.loadtxt('../../input_4_cast/output/%s/%s_mn_eps_1p0E-02/Pnonlin-zk.txt'%(case,param))
    kk_fid = np.loadtxt('../../input_4_cast/output/%s/fiducial_eps_0/k_values_list.txt'%(case))
    pl_fid = np.loadtxt('../../input_4_cast/output/%s/fiducial_eps_0/Plin-zk.txt'%(case))
    pn_fid = np.loadtxt('../../input_4_cast/output/%s/fiducial_eps_0/Pnonlin-zk.txt'%(case))
    der_pl = (pl_plus-pl_minus)/2./step/pl_fid
    der_pn = (pn_plus-pn_minus)/2./step/pl_fid
    if kk_plus[0] != kk_minus[0] or kk_plus[0] != kk_fid[0]:
        print ("problem with k sampling! Don't trust results! Interpolation needed!")
    if kk_plus[-1] != kk_minus[-1] or kk_plus[-1] != kk_fid[-1]:
        print ("problem with k sampling! Don't trust results! Interpolation needed!")
    if len(kk_plus) != len(kk_minus) or len(kk_plus) != len(kk_fid):
        print ("problem with k sampling! Don't trust results! Interpolation needed!")
    return kk_fid,der_pl,der_pn

#########

params = ['Ob','Om','h','ns','s8','w0','wa']
latex_names = ['$\Omega_b$','$\Omega_m$','$h$','$n_s$','$\sigma_8$','$w_0$','$w_a$']

param_num = len(params)

fig, axs = plt.subplots(param_num,2,figsize=(9,10),sharex=True)

logscale = False

for i, param in enumerate(params):

    if i == 0:
        axs[i,0].set_title(r'$\partial_\alpha \ln P_L(k,z=0)$')
        axs[i,1].set_title(r'$\partial_\alpha \ln P_{NL}(k,z=0)$')
    if i == param_num-1:
        axs[i,0].set_xlabel(r'$k$   (1/Mpc)')
        axs[i,1].set_xlabel(r'$k$   (1/Mpc)')

    axs[i,0].set_ylabel(latex_names[i])

    axs[i,0].set_xlim(1.e-4,50.)
    axs[i,1].set_xlim(1.e-4,50.)

    #if params == 'w0':
    #    axs[i,0].set_ylim(-600,600)
    #    axs[i,1].set_ylim(-600,600)
    #elif params == 'wa':
    #    axs[i,0].set_ylim(-100,100)
    #    axs[i,1].set_ylim(-100,100)
    #else:
    #    axs[i,0].set_ylim(-60000,60000)
    #    axs[i,1].set_ylim(-60000,60000)

    #factor=1
    #if

    #axs[i,0].set_title(str(latex_names[i]),pad=-0.1)
    #axs[i,1].set_title(str(latex_names[i]),pad=-0.1)

    kk,dl,dn = derivative('class_w0wa_DP',param,0.01)
    if logscale == False:
        axs[i,0].semilogx(kk[:],dl[0,:],'k-',label='CLASS DP',linewidth=2)
        axs[i,1].semilogx(kk[:],dn[0,:],'k-',linewidth=2)
    else:
        axs[i,0].loglog(kk[:],dl[0,:],'k-',label='CLASS DP',linewidth=2)
        axs[i,1].loglog(kk[:],dn[0,:],'k-',linewidth=2)
        axs[i,0].loglog(kk[:],-dl[0,:],'k--',linewidth=2)
        axs[i,1].loglog(kk[:],-dn[0,:],'k--',linewidth=2)

    kk,dl,dn = derivative('class_w0wa_HP',param,0.01)
    if logscale == False:
        axs[i,0].semilogx(kk[:],dl[0,:],'g-',label='CLASS HP',linewidth=0.9)
        axs[i,1].semilogx(kk[:],dn[0,:],'g-',linewidth=0.9)
    else:
        axs[i,0].loglog(kk[:],dl[0,:],'g-',label='CLASS HP',linewidth=0.9)
        axs[i,1].loglog(kk[:],dn[0,:],'g-',linewidth=0.9)
        axs[i,0].loglog(kk[:],-dl[0,:],'g--',linewidth=0.9)
        axs[i,1].loglog(kk[:],-dn[0,:],'g--',linewidth=0.9)

    kk,dl,dn = derivative('camb_w0wa_P1',param,0.01)
    if logscale == False:
        axs[i,0].semilogx(kk[:],dl[0,:],'c-',label='CAMB P1',linewidth=0.9)
        axs[i,1].semilogx(kk[:],dn[0,:],'c-',linewidth=0.9)
    else:
        axs[i,0].loglog(kk[:],dl[0,:],'c-',label='CAMB P1',linewidth=0.9)
        axs[i,1].loglog(kk[:],dn[0,:],'c-',linewidth=0.9)
        axs[i,0].loglog(kk[:],-dl[0,:],'c--',linewidth=0.9)
        axs[i,1].loglog(kk[:],-dn[0,:],'c--',linewidth=0.9)

    kk,dl,dn = derivative('camb_w0wa_P2',param,0.01)
    if logscale == False:
        axs[i,0].semilogx(kk[:],dl[0,:],'b-',label='CAMB P2',linewidth=0.9)
        axs[i,1].semilogx(kk[:],dn[0,:],'b-',linewidth=0.9)
    else:
        axs[i,0].loglog(kk[:],dl[0,:],'b-',label='CAMB P2',linewidth=0.9)
        axs[i,1].loglog(kk[:],dn[0,:],'b-',linewidth=0.9)
        axs[i,0].loglog(kk[:],-dl[0,:],'b--',linewidth=0.9)
        axs[i,1].loglog(kk[:],-dn[0,:],'b--',linewidth=0.9)

    kk,dl,dn = derivative('camb_w0wa_P3',param,0.01)
    if logscale == False:
        axs[i,0].semilogx(kk[:],dl[0,:],'r-',label='CAMB P3',linewidth=0.9)
        axs[i,1].semilogx(kk[:],dn[0,:],'r-',linewidth=0.9)
    else:
        axs[i,0].loglog(kk[:],dl[0,:],'r-',label='CAMB P3',linewidth=0.9)
        axs[i,1].loglog(kk[:],dn[0,:],'r-',linewidth=0.9)
        axs[i,0].loglog(kk[:],-dl[0,:],'r--',linewidth=0.9)
        axs[i,1].loglog(kk[:],-dn[0,:],'r--',linewidth=0.9)

    if i==0:
        axs[i,0].legend()

#plt.show()
plt.savefig('derivatives.pdf')
