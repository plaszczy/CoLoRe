from pylab import *

from tools import *

import sys


def dndz_spectra(cl1,cl2):
        #figure()
        subplot(2,1,1)
        
        title(n)
        plot(cl1,label='const')
        plot(cl2,'r',label='red')
        xlim(1,200)
        legend()
        ax0()
        ylabel("Cl")
        rcl=cl2-cl1
        subplot(2,1,2)
        plot(rcl[1:200],'r')
        xlim(1,200)
        ax0()
        ylabel("Cl(red)-Cl(const)")
        
        tight_layout()

def angpow(t,ta,h,lmax=500):
    figure()
    it=gca()._get_lines.prop_cycler
    names=t.names
    for n in names[1:] :
        cl=t[n][1:lmax]
        cla=ta[n][1:lmax]
        if h['REMOVESN']==False:
            Ntot=h["ngal"+n[-1]]
            sn=4*pi/Ntot
            cl-=sn
        col=next(it)['color']
        plot(cl,label=n,color=col)
        plot(cla,'--',color=col)
        plot(cl-cla,color=col)
        ax0()
    legend()

#data
t,h=mrdfits("outputs/bench0/clmean_red.fits",1,header=True)
ta=mrdfits("model/cltophat_logn_red.fits",1)

#t,h=mrdfits("outputs/bench0/clmean.fits",1,header=True)
#ta=mrdfits("model/cltophat_logn.fits",1)

#

angpow(t,ta,h)


#print("||{:04.2e} || {:04.2e} || {:04.2e} ||".format(4*pi/Ntot,S2,S2/(4*pi/Ntot)))
