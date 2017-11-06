from pylab import *
from healpy import *
import os,sys


dens_type=int(sys.argv[1])
ishell=5



dens=("LogN","1LPT","2LPT","Gaussclip")

dir="batch_{}".format(dens_type)



#truth="clR4_shell{:0d}_flat.txt".format(ishell)
truth="clR4_shell{:0d}_bordersout.txt".format(ishell)
#truth="clR4_shell{:0d}_bordersin.txt".format(ishell)

print("Analyzing batch={} vs {}".format(dir,truth))

f1=os.path.join("batch","dens_type{}".format(dens_type),"shell{}".format(ishell),"clmean.fits")
clrec=read_cl(f1)
l=arange(len(clrec))

lt,clt=loadtxt(truth,unpack=True)
#clt=read_cl("cl_R4.fits")
clt=clt[l]
clt[0]=0


figure()
plot(clt,'r',label=r"$C_\ell^{th}$")
plot(clrec,'k',label=r"$<C_\ell^i>-SN$")
plot(clrec-clt,label='residue')
axhline(0,color='k',lw=0.5)
legend()
xlabel(r"$\ell$")
ylabel(r"$C_\ell$")
xlim(0,500)
ylim(-2e-5,8e-5)
title(dens[dens_type])
tight_layout()

show()
