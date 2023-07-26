# verification of Stupakov's formulas implementation (PRST-AB 2007)
# comparing with M. Zobov - O. Frasciello wakes and impedances
# (ipython commands)

ipython
from Impedance import *
# plot of G function (see paper)
x=np.arange(0.,4.01,0.01)
G=np.zeros((len(x),4))
for i in range(4):
    G[:,i]=G_Stupakov(x,i);

pylab.plot(x,G)
pylab.show()

# check w.r.t. M. Zobov curve of BB impedance (current coll.)
garray=np.arange(1.,12.,0.5)*1e-3;
w=70e-3/2;l=97e-3;delta=17.6e-3;
fr=5e9;Q=1;
R=np.zeros((len(garray),5));
Rexact=np.zeros((len(garray),5));
Rtround=np.zeros(len(garray));
Rlround=np.zeros(len(garray));
for ig,g in enumerate(garray):
    R[ig,:]=2*broadband_imp_taper_flat_Stupakov(g,g+delta,delta/l,w,fr);
    Rexact[ig,:]=2*broadband_imp_taper_flat_Stupakov(g,g+delta,delta/l,w,fr,approx=False);
    Rtround[ig]=2*transverse_imp_taper_round_Yokoya(g,g+delta,delta/l);
    Rlround[ig]=2*longitudinal_imp_taper_round_Yokoya(g,g+delta,delta/l,fr);

fig,ax=init_figure();
ax.semilogy(garray*1e3,R[:,2]/1e3,'-xb',label='Approximation',lw=3.,ms=10.);
ax.semilogy(garray*1e3,Rexact[:,2]/1e3,'-r',label='Exact (with G1 function)',lw=3.,ms=10.);
ax.set_xlabel('Half-gap [mm]');
ax.set_ylabel(r"Constant inductive impedance [k$ \Omega/ $m]")
end_figure(fig,ax)
ax.grid(which='minor')
pylab.show()

# check all components
listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad'];
for icomp,comp in enumerate(listcomp):
    fig,ax=init_figure();
    ax.semilogy(garray*1e3,np.abs(R[:,icomp])/1e3,'-xb',label=comp+', flat, approx.',lw=3.,ms=10.);
    ax.semilogy(garray*1e3,np.abs(Rexact[:,icomp])/1e3,'-r',label=comp+', flat, exact',lw=3.,ms=10.);
    if icomp==0: ax.semilogy(garray*1e3,Rlround/1e3,'-g',label=comp+', round',lw=3.,ms=10.);
    elif icomp<3: ax.semilogy(garray*1e3,Rtround/1e3,'-g',label=comp+', round',lw=3.,ms=10.);
    ax.set_xlabel('Half-gap [mm]');
    ax.set_ylabel(r"Constant inductive impedance [k$ \Omega/ $m]")
    end_figure(fig,ax)
    ax.grid(which='minor')

pylab.show()

# check kick factors (current coll.) vs M. Zobov/O. Fraciello results
garray=np.array([1.,3.,5.])*1e-3;
w=70e-3/2;l=97e-3;delta=17.6e-3;
fr=5e9;Q=1;
kick=np.zeros(len(garray));kickex=np.zeros(len(garray));
for ig,g in enumerate(garray):
    R=2*broadband_imp_taper_flat_Stupakov(g,g+delta,delta/l,w,fr);
    Rexact=2*broadband_imp_taper_flat_Stupakov(g,g+delta,delta/l,w,fr,approx=False);
    imp_mod,wake_mode=imp_model_resonator(R,fr,Q,listcomp=['Zlong', 'Zxdip', 'Zydip', 'Zxquad', 'Zyquad']);
    imp_mod_ex,wake_mode=imp_model_resonator(Rexact,fr,Q,listcomp=['Zlong', 'Zxdip', 'Zydip', 'Zxquad', 'Zyquad']);
    kick[ig]=transverse_kick_factor(imp_mod,7.5e-2,compname='Zydip')
    kickex[ig]=transverse_kick_factor(imp_mod_ex,7.5e-2,compname='Zydip')

print(kick,kickex)

# check kick factors (new coll., BPM cavity) vs M. Zobov/O. Fraciello results
garray=np.array([1.,3.,5.])*1e-3;
l1=(95-37.32-31.9)*1e-3;delta1=8e-3;w1=(33e-3+23.3e-3)/2;
l2=37.32e-3;delta2=10.7e-3;w2=70e-3/2;
fr=5e9;Q=1;
kick=np.zeros(len(garray));
for ig,g in enumerate(garray):
    R1=2*broadband_imp_taper_flat_Stupakov(g,g+delta1,delta1/l1,w1,fr,approx=False);
    R2=2*broadband_imp_taper_flat_Stupakov(g+delta1,g+delta1+delta2,delta2/l2,w2,fr,approx=False);
    imp_mod,wake_mode=imp_model_resonator(R1+R2,fr,Q,listcomp=['Zlong', 'Zxdip', 'Zydip', 'Zxquad', 'Zyquad']);
    kick[ig]=transverse_kick_factor(imp_mod,7.5e-2,compname='Zydip')

print(kick)

# check wake potential formula (after previous one with current coll., g=5mm)
R=22753.30122837;
sigmaz=7.5e-2;zscan=np.arange(-4*sigmaz,4.0005*sigmaz,sigmaz/1000.);
w=resonator_wake_potential(R,fr,Q,sigmaz,zscan,beta=1)
lambdaz=np.exp(-zscan**2/(2.*sigmaz**2)) / (np.sqrt(2*np.pi)*sigmaz); # line density
np.trapz(lambdaz*w,zscan)/kick[-1]
# -1 (due to sign convention in kick factor) -> fine !

# compare Stupakov wake potential to real wake potential from M. Zobov and O. Frasciello
filename='../../LHC_impedance_and_scripts/Coll_settings/TCS_TCT_geometric/GdfidLWakepotential_dip_halfgap5mm_bunchlength2mm_mesh0p2mm.txt'
s=read_ncol_file(filename)
pylab.plot(s[:,0],-s[:,1])
w=resonator_wake_potential(R,5e9,1,2e-3,s[:,0],beta=1)
pylab.plot(s[:,0],w)
wmeas=-s[:,1];z=s[:,0]
pylab.show()

# attempt to fit the real wake potential with resonator
# 3D
from scipy import optimize
from tables_lib import *
def f3D(param,z,wmeas,sigmaz):
    R=param[0];fr=param[1];Q=param[2];
    w=resonator_wake_potential(R,fr,Q,sigmaz,z,beta=1);
    return np.sum((w-wmeas)**2);

# 2D
def f2D(param,Q,z,wmeas,sigmaz):
    R=param[0];fr=param[1];
    w=resonator_wake_potential(R,fr,Q,sigmaz,z,beta=1);
    return np.sum((w-wmeas)**2);

# n D
def f(param,z,wmeas,sigmaz):
    R=param[::3]*1e4;fr=param[1::3]*1e9;Q=param[2::3];
    w=resonator_wake_potential(R,fr,Q,sigmaz,z,beta=1);
    return w-wmeas;

def fbis(param,z,wmeas,sigmaz):
    R=param[::3];fr=param[1::3];Q=param[2::3];
    w=resonator_wake_potential(R,fr,Q,sigmaz,z,beta=1);
    return np.sum((w-wmeas)**2);

def fter(parfit,z,wmeas,sigmaz,parfix,parfixind):
    param=np.zeros(len(parfit)+len(parfix));
    param[parfixind]=parfix;
    param[complementary(parfixind,list(range(len(param))))]=parfit;
    R=param[::3];fr=param[1::3];Q=param[2::3];
    w=resonator_wake_potential(R,fr,Q,sigmaz,z,beta=1);
    return np.sum((w-wmeas)**2);

ind=np.where(z*(z-3e-3)<0)[0];

par=optimize.fmin(f2D,[R,5e9],args=(1,z[ind],wmeas[ind],2e-3),maxfun=20000,maxiter=200)
w=resonator_wake_potential(par[0],par[1],1,2e-3,z,beta=1)
pylab.plot(z,wmeas,z,w)

ind=np.where(z-3e-3>0)[0];

par2=optimize.fmin(f3D,[par[0]/2,par[1]/2,10],args=(z[ind],wmeas[ind]-w[ind],2e-3),maxfun=20000,maxiter=200)
w2=resonator_wake_potential(par2[0],par2[1],par2[2],2e-3,z,beta=1)
pylab.plot(z,w+w2)

par,inf=optimize.leastsq(f,np.ones(6),args=(z,wmeas,2e-3),maxfev=20000);
w=resonator_wake_potential(par[::3]*1e4,par[1::3]*1e9,par[2::3],1,2e-3,z,beta=1)
pylab.plot(z,wmeas,z,w)

par=optimize.brute(fbis,(slice(0,2e4,1e3),slice(1e10,1e11,1e10),slice(1,2,1),slice(0,1e4,1e3),slice(5e8,5e9,5e8),slice(1,100,5)),args=(z,wmeas,2e-3),Ns=20);
w=resonator_wake_potential(par[::3],par[1::3],par[2::3],2e-3,z,beta=1)
pylab.plot(z,wmeas,z,w)

# fit very first peak
ind=np.where(z*(z-3e-3)<0)[0];
parfixind=[2];parfix=[1];
parfit=optimize.fmin(fter,[R,5e10],args=(z[ind],wmeas[ind],2e-3,parfix,parfixind),maxfun=20000)
par=np.zeros(len(parfit)+len(parfix));
par[parfixind]=parfix;
par[complementary(parfixind,list(range(len(par))))]=parfit;
w=resonator_wake_potential(par[::3],par[1::3],par[2::3],2e-3,z,beta=1)
pylab.plot(z,wmeas,z,w)
pylab.plot(z,wmeas-w)
wfunc=resonator_wake(par[0],par[1],par[2],z/c*1e9)
pylab.plot(z,wfunc)
pylab.plot(z,wmeas-wfunc.reshape((1,-1))[0])


ind=np.where((z-3e-3)*(z-1e-2)<0)[0];
parfixind=[2];parfix=np.hstack(([1]));
#parfit=optimize.fmin(fter,[par[0],R/10,5e9],args=(z[ind],wmeas[ind]-w[ind],2e-3,parfix,parfixind),maxfun=20000,maxiter=200)
parfit,som,info=optimize.fmin_l_bfgs_b(fter,[par[0],par[1],par[0],5e9,1],args=(z[ind],wmeas[ind],2e-3,parfix,parfixind),bounds=((100,R),(1e9,1e11),(10,par[0]),(1e8,5e10),(1,100)),approx_grad=True,maxfun=20000)
#parfit=optimize.brute(fter,((0,par[0]),(1e9,par[1]),(1,20)),args=(z[ind],wmeas[ind]-w[ind],2e-3,parfix,parfixind),Ns=20,finish=optimize.fmin)
par2=np.zeros(len(parfit)+len(parfix));
par2[parfixind]=parfix;
par2[complementary(parfixind,list(range(len(par2))))]=parfit;
w2=resonator_wake_potential(par2[::3],par2[1::3],par2[2::3],2e-3,z,beta=1)
pylab.plot(z,wmeas,z,w2)
pylab.plot(z,wmeas-w2)

# take first peak larger
ind=np.where(z*(z-2e-2)<0)[0];
parfixind=[2];parfix=[1];
parfit=optimize.fmin(fter,[R,5e10],args=(z[ind],wmeas[ind],2e-3,parfix,parfixind),maxfun=20000)
par=np.zeros(len(parfit)+len(parfix));
par[parfixind]=parfix;
par[complementary(parfixind,list(range(len(par))))]=parfit;
w=resonator_wake_potential(par[::3],par[1::3],par[2::3],2e-3,z,beta=1)
pylab.plot(z,wmeas,z,w)
pylab.plot(z,wmeas-w)
wfunc=resonator_wake(par[0],par[1],par[2],z/c*1e9)
pylab.plot(z,wfunc)
pylab.plot(z,wmeas-wfunc.reshape((1,-1))[0])

ind=np.where((z-5e-2)<0)[0];
parfixind=[0,1,2];parfix=np.hstack((par,[]));
parfit=optimize.fmin(fter,[par[0],5e9,1],args=(z[ind],wmeas[ind],2e-3,parfix,parfixind),maxfun=20000,maxiter=200)
#parfit,som,info=optimize.fmin_l_bfgs_b(fter,[par[0],5e9,1],args=(z[ind],wmeas[ind],2e-3,parfix,parfixind),bounds=((100,par[0]),(1e8,5e10),(1,100)),approx_grad=True,maxfun=20000)
#parfit=optimize.brute(fter,((0,par[0]),(1e9,par[1]),(1,20)),args=(z[ind],wmeas[ind]-w[ind],2e-3,parfix,parfixind),Ns=20,finish=optimize.fmin)
par2=np.zeros(len(parfit)+len(parfix));
par2[parfixind]=parfix;
par2[complementary(parfixind,list(range(len(par2))))]=parfit;
w2=resonator_wake_potential(par2[::3],par2[1::3],par2[2::3],2e-3,z,beta=1)
pylab.plot(z,wmeas,z,w2)
pylab.plot(z,wmeas-w2)
# all of this is a mess... no way to fit correctly apparently...
# use orthogonal polynomials ? (Hermite ? see eqs. (30)-(31) p. 195 in Higher trans. func. (Erdelyi) part 2)

# some tests with fft
c=299792458.;sigmaz=0.002;
fileReZ='../../LHC_impedance_and_scripts/Coll_settings/TCS_TCT_geometric/ReZX_halfgap5mm_bunchlength2mm_mesh0p2mm_offset1mm_sent22082013.txt';
fileImZ='../../LHC_impedance_and_scripts/Coll_settings/TCS_TCT_geometric/ImZX_halfgap5mm_bunchlength2mm_mesh0p2mm_offset1mm_sent22082013.txt';
ReZ=read_ncol_file(fileReZ,ignored_rows=2)
ImZ=read_ncol_file(fileImZ,ignored_rows=2)
pylab.semilogx(ReZ[:,0],ReZ[:,1]*1e3,ImZ[:,0],ImZ[:,1]*1e3); # offset of 1mm for impedance to be taken into account
f=np.arange(len(z))*c/(len(z)*np.diff(z)[0]);
Zt=1j*np.fft.fft(wmeas)/c/np.exp(-(2.*np.pi*f)**2*sigmaz**2/(2.*c**2))
pylab.plot(f[~np.isinf(Zt)],np.real(Zt[~np.isinf(Zt)]),f[~np.isinf(Zt)],np.imag(Zt[~np.isinf(Zt)]));pylab.axis([0,1e11,0,1e5])
# it's a big mess

# plot all dipolar wakes
col=['b','r','g','m','k','c','y'];
fig,ax=init_figure()
halfgapscan=[1,3,5,11.5,20];
for ihg,hg in enumerate(halfgapscan):
    filename='../../LHC_impedance_and_scripts/Coll_settings/TCS_TCT_geometric/GdfidLWakepotential_dip_halfgap'+str(hg)+'mm_bunchlength2mm_mesh0p2mm.txt'
    s=read_ncol_file(filename)
    ax.plot(s[:,0],-s[:,1],'-'+col[ihg])

pylab.show()

# compare Stupakov's dip. wake potentials & functions to these
width=70e-3/2;l=97e-3;delta=17.6e-3;
fr=5e9;Q=1;sigmaz=2e-3;
for ig,g in enumerate(np.array(halfgapscan)*1e-3):
    R=2*broadband_imp_taper_flat_Stupakov(g,g+delta,delta/l,width,fr,approx=False);
    w=resonator_wake_potential(R[2],fr,Q,sigmaz,s[:,0],beta=1);
    ax.plot(s[:,0],w,'--'+col[ig])
    wf=resonator_wake(R[2],fr,Q,s[:,0]*1e9/299792458.);
    ax.plot(s[:,0],wf,':'+col[ig])

pylab.show()

# plot all quadrupolar wakes
col=['b','r','g','m','k','c','y'];
fig,ax=init_figure()
halfgapscan=[1,3,5,11.5,20];
for ihg,hg in enumerate(halfgapscan):
    filename='../../LHC_impedance_and_scripts/Coll_settings/TCS_TCT_geometric/GdfidLWakepotential_quad_halfgap'+str(hg)+'mm_bunchlength2mm_mesh0p2mm.txt'
    s=read_ncol_file(filename)
    ax.plot(s[:,0],-s[:,1],'-'+col[ihg])

pylab.show()

# compare Stupakov's quad. wake potentials & functions to these
width=70e-3/2;l=97e-3;delta=17.6e-3;
fr=5e9;Q=1;sigmaz=2e-3;
for ig,g in enumerate(np.array(halfgapscan)*1e-3):
    R=2*broadband_imp_taper_flat_Stupakov(g,g+delta,delta/l,width,fr,approx=False);
    w=resonator_wake_potential(R[4],fr,Q,sigmaz,s[:,0],beta=1);
    ax.plot(s[:,0],w,'--'+col[ig])
    wf=resonator_wake(R[4],fr,Q,s[:,0]*1e9/299792458.);
    ax.plot(s[:,0],wf,':'+col[ig])

pylab.show()


# plot all longitudinal wakes
col=['b','r','g','m','k','c','y'];
fig,ax=init_figure()
halfgapscan=[1,3,5,11.5];
for ihg,hg in enumerate(halfgapscan):
    filename='../../LHC_impedance_and_scripts/Coll_settings/TCS_TCT_geometric/GdfidLWakepotential_long_halfgap'+str(hg)+'mm_bunchlength2mm_mesh0p2mm.txt'
    s=read_ncol_file(filename)
    ax.plot(s[:,0],-s[:,1])

pylab.show()

# compare Stupakov's long. wake FUNCTIONS to these
width=70e-3/2;l=97e-3;delta=17.6e-3;
fr=5e10;Q=1;sigmaz=2e-3;
for ig,g in enumerate(np.array(halfgapscan)*1e-3):
    R=2*broadband_imp_taper_flat_Stupakov(g,g+delta,delta/l,width,fr,approx=False);
    w=resonator_long_wake_potential(R[0],fr,Q,sigmaz,s[:,0],beta=1)
    ax.plot(s[:,0],w,'--'+col[ig])
    wf=resonator_long_wake(R[0],fr,Q,s[:,0]*1e9/299792458.);
    ax.plot(s[:,0],wf,':'+col[ig])

pylab.show()
