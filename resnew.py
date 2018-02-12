import batch as b
from pylab import *
from tools import *

truth=mrdfits("model/tophat_dens0_ngrid1024.fits",1)
t=mrdfits("outputs/bench0/clmean.fits",1)
cols=t.dtype.names[1:]
l=t.ell

do_legend=True
newfig=True

zval=(0,0.1,0.2,0.3,0.4,0.5)

Nz=5000
Nsamp=41253*Nz*0.1
print("Nz={} mean Ngal={}".format(Nz,Nsamp))
nbar=Nsamp/(4*np.pi)


if newfig:
    fig=figure(figsize=(12,10))
    #fig.suptitle("Ngrid={}".format(ngrid), fontsize=14)
for i in range(1,5):
    if i==1 :
        ax=subplot(2,2,i)
    else:
        subplot(2,2,i,sharex=ax,sharey=ax)
    ishell=i+1
    clt=truth[cols[i-1]]
    clt[0]=0
    clrec=t[cols[i-1]]
    cvar=2*(clt+1/nbar)**2/(2*l+1)
    plot(clrec,'k',label=r"$<C_\ell^i>-SN$",lw=3)
    plot(clt,'r',label=r"$C_\ell^{th}$",lw=1)
    plot(clrec-clt,'b',label='residue')
    fill_between(l,clt+sqrt(cvar),clt-sqrt(cvar),color='r',alpha=0.1)
    plot(l,clrec-clt)
    zmax=zval[ishell]
    zmin=zval[ishell-1]
    legend()
    title(r"$z\in [{},{}]$".format(zmin,zmax),fontsize=14)
    #text(0.4, 0.9,r"$z\in [{},{}]$".format(zmin,zmax), transform=gca().transAxes,fontsize=12)

    ax0()
    xlabel(r"$\ell$")
    xlim(0,250)
    ylabel(r"$C_\ell$")
ylim(-5e-5,2e-4)
tight_layout()
