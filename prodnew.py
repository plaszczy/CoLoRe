from pylab import *
from tools import *

def ratio(t,truth,h,lmax=200):
    names=t.names
    for n in names[1:] :
        cl=t[n][1:lmax]
        if h['REMOVESN']==False:
            SN=mean(t[n][600:])
            print("estimating SN for {}={}".format(n,SN))
            cl-=SN
        #plot(t.ell[1:lmax],cl,label=n)
        plot(t.ell[1:lmax],cl/truth[n][1:lmax],label=n)
        xlabel(r"$\ell$")
        ylabel(r"$<C_\ell(rec)>/C_\ell(true)$")
        legend()
        ylim(0.8,1.2)
        ax1()
        #ax0()

#LOGN

#tophat
#data="outputs/bench0/clmean.fits"
#model="model/cltophat_logn.fits"
#
#data="outputs/bench0/clmean_norsd.fits"
#model="model/cltophat_logn_norsd.fits"
#gauss
#data="outputs/gauss4_dens0/clmean.fits"
#model="model/clgauss_logn.fits"

###CLIP
###tophat
data="outputs/bench3/clmean.fits"
model="model/cltophat_clip.fits"
#gauss
#data="outputs/gauss4_dens3/clmean.fits"
#model="model/clgauss_clip.fits"

#SMEARED
#data="outputs/tophat4/clmean.fits"
#model="model/tophat_logn120.fits"

tmod=mrdfits(model,1)
(t,h)=mrdfits(data,1,header=True)


figure()
ratio(t,tmod,h,200)
title(model)
tight_layout()
