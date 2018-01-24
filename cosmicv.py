from pylab import *
from scipy.signal import savgol_filter
from tools import *
import batch as b

dens=0
rsd=True
shell=(2,3,4,5)
#ngrid=256
#ngrid=1024
ngrid=512
newfig=True


zval=(0,0.1,0.2,0.3,0.4,0.5)

if newfig:
    fig=figure()
    fig.suptitle("Ngrid={}".format(ngrid), fontsize=14)
for i in range(1,5):
    if i==1 :
        ax=subplot(2,2,i)
    else:
        subplot(2,2,i,sharex=ax,sharey=ax)
    ishell=i+1
    l,clt=b.model(dens_type=dens,ngrid=ngrid,rsd=rsd,ishell=ishell)
    cl=b.get(dens_type=dens,ngrid=ngrid,rsd=rsd,ishell=ishell)
    res=cl-clt
    var=2*clt**2/(2*l+1)
    plot(res/sqrt(var))
    xlim(1,200)
    ylim(-1,1)
    zmax=zval[ishell]
    zmin=zval[ishell-1]
    title(r"$z\in [{},{}]$".format(zmin,zmax),fontsize=9)
    ax0()
    xlabel(r"$\ell$",fontsize=12)
    ylabel(r"$\Delta C_\ell/\sigma_{CV}$",fontsize=12)
tight_layout()
print("ngrid={}".format(ngrid))
