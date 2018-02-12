from pylab import *
from scipy.signal import savgol_filter
from tools import *
import batch as b


truth=mrdfits("model/tophat_dens0_ngrid1024.fits",1)
t=mrdfits("outputs/bench0/clmean.fits",1)
cols=t.dtype.names[1:]
l=t.ell

newfig=True


Nz=5000
Nsamp=41253*Nz*0.1
print("Nz={} mean Ngal={}".format(Nz,Nsamp))
nbar=Nsamp/(4*np.pi)

zval=(0,0.1,0.2,0.3,0.4,0.5)

if newfig:
    fig=figure()
for i in range(1,5):
    if i==1 :
        ax=subplot(2,2,i)
    else:
        subplot(2,2,i,sharex=ax,sharey=ax)
    ishell=i+1
    clt=truth[cols[i-1]]
    clt[0]=0
    cl=t[cols[i-1]]
    res=cl-clt
    var=2*(clt+1/nbar)**2/(2*l+1)
    plot(l,clt/max(clt),'k--')
    plot(res/sqrt(var),'k')
    xlim(1,200)
    ylim(-1,1)
    zmax=zval[ishell]
    zmin=zval[ishell-1]
    title(r"$z\in [{},{}]$".format(zmin,zmax),fontsize=9)
    ax0()
    xlabel(r"$\ell$",fontsize=12)
    ylabel(r"$\frac{\Delta C_\ell}{\sigma_{CV}}$",fontsize=12)
tight_layout()
print("ngrid={}".format(ngrid))
