from pylab import *
from tools import *
from healpy import *

import healpy as hp
import glob,os


#cls
clt=read_cl("cl_R4.fits")
clt[0]=0.

nbar=2066239./(4*pi)
clsn=ones_like(clt)*1/nbar

#p2=getPixWin2(256)

files=glob.glob("batch/cl_map256_cat*.fits")

#no rsd
#clt=read_cl("cl_norsd_R4.fits")
#files=glob.glob("batch/cl_map_norsd_gal*.fits")


cls=[]
for file in files:
    print(file)
    cls.append(hp.read_cl(file))

#[plot(cl) for cl in cls]

clm=mean(cls,0,dtype=float64)

lmax=min([len(clm),len(clt),len(clsn)])-1


clsn=clsn[0:lmax]
clt=clt[0:lmax]
clm=clm[0:lmax]

#c=cov(transpose(cls))
#cc=cov2cor(c)
#sig2=diag(c)
#l=arange(shape(cls)[1])
#errorbar(l,clm,yerr=sqrt(sig2))

#cut clt and put monopole to 0
#clt=clt[0:len(clm)]

res=clm-clt
snlevel=mean(clm[400:]-clt[400:])
print("sn level={}".format(snlevel))

#p2=p2[0:len(clm)]

    
figure()
plot(clt+clsn,'r',label=r"$C_\ell^{th}+SN$")
plot(clm,'k',label=r"$<C_\ell^i>$")
plot(clm-(clt+clsn),label='residue')
#plot(clsn,'r',label='SN')
#plot(clsn*p2,label='SN*pixwin')

axhline(0,color='k',lw=0.5)
legend()
xlabel(r"$\ell$")
ylabel(r"$C_\ell$")
xlim(0,500)
tight_layout()

show()
