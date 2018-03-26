from pylab import *
from tools import *

def ratio(t,truth,h,lmax=200):
    names=t.names
    for n in names[1:] :
        cl=t[n][1:lmax]
        if h['REMOVESN']==False:
            SN=mean(t[n][500:])
            print("estimating SN for {}={}".format(n,SN))
            #cl-=SN
        plot(t.ell[1:lmax],cl,label=n)
        #plot(t.ell[1:lmax],cl/truth[n][1:lmax],label=n)
        xlabel(r"$\ell$")
        ylabel(r"$C_\ell(rec)/C_\ell(true)$")
        legend()
        #ylim(0,3)
        #ax1()
        ax0()

#
#data="outputs/bench3/clmean.fits"
#model="model/tophat_clip120.fits"

data="outputs/tophatsmear_logn/clmean.fits"
data="outputs/gauss4_dens0/clmean.fits"
model="model/clgauss_logn.fits"

#data="outputs/tophatsmear_clip/clmean.fits"
#model="model/clgauss_clip.fits"


#data="outputs/tophat4/clmean.fits"
#model="model/tophat_logn120.fits"

truth=mrdfits(model,1)
(t,h)=mrdfits(data,1,header=True)


#figure()
ratio(t,truth,h,200)
title(model)
tight_layout()
