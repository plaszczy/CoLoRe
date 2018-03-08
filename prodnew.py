from pylab import *
from tools import *

def ratio(t,truth,lmax=200):
    names=t.names
    [plot(t.ell[1:lmax],t[n][1:lmax]/truth[n][1:lmax],label=n) for n in names[1:]]
    xlabel(r"$\ell$")
    ylabel(r"$C_\ell(rec)/C_\ell(true)$")
    legend()
    ax1()

#clip
#data="outputs/bench3/clmean.fits"
#model="model/tophat_clip120.fits"
#logn
data="batch/gauss_bench0/clmean.fits"
model="model/clgauss_logn.fits"


truth=mrdfits(model,1)
t=mrdfits(data,1)

figure()
ratio(t,truth)
title(model)
ylim(0.85,1.15)
tight_layout()

#lognormal
data="batch/bench0/clmean.fits"
model="model/tophat_logn120.fits"

truth=mrdfits(model,1)
t=mrdfits(data,1)

figure()
ratio(t,truth)
title(model)
ylim(0.85,1.15)
tight_layout()


