from pylab import *
from tools import *

def ratio(t,truth,lmax=200):
    names=t.names
    [plot(t.ell[1:lmax],t[n][1:lmax]/truth[n][1:lmax],label=n) for n in names[1:]]
    xlabel(r"$\ell$")
    ylabel(r"$C_\ell(rec)/C_\ell(true)$")
    legend()
    ylim(0,3)
    ax1()

#
#data="outputs/bench3/clmean.fits"
#model="model/tophat_clip120.fits"
#gauss
#data="outputs/gauss4/clmean.fits"
#model="model/clgauss_logn.fits"

data="outputs/tophat4/clmean.fits"
model="model/tophat_logn120.fits"

truth=mrdfits(model,1)
t=mrdfits(data,1)

figure()
ratio(t,truth)
title(model)
tight_layout()

#lognormal
data="outputs/bench0/clmean.fits"
model="model/tophat_logn120.fits"

truth=mrdfits(model,1)
t=mrdfits(data,1)

figure()
ratio(t,truth)
title(model)
tight_layout()


