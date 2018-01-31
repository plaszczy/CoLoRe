import batch as b
from pylab import *
from tools import *

dens=0
rsd=True
shell=(2,3,4,5)
ngrid=256
#ngrid=512
ngrid=1024
do_legend=True
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
    plot(l,cl-clt,'k',alpha=.7)
    fill(l,clt/max(clt)*1e-6,'r',alpha=0.7)
    #diff_fac(l,clt,cl,do_legend=do_legend)
    zmax=zval[ishell]
    zmin=zval[ishell-1]
    title(r"$z\in [{},{}]$".format(zmin,zmax),fontsize=9)
    #text(0.4, 0.9,r"$z\in [{},{}]$".format(zmin,zmax), transform=gca().transAxes,fontsize=12)

    ax0()
    xlabel(r"$\ell$")
    xlim(0,500)
    #ylabel(r"$\Delta C_\ell$")
tight_layout()
ylim(-1e-6,1e-6)
print("ngrid={}".format(ngrid))
