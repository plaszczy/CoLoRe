from pylab import *
from healpy import *
import os

ishell=5
truth="clR4_shell{:0d}_bordersout.txt".format(ishell)

lt,clt=loadtxt(truth,unpack=True)
clt[0]=0



dens=("LogN","1LPT","2LPT","Gaussclip")

i=0
for dens_type in dens :
    f1=os.path.join("batch","dens_type{}".format(i),"shell{}".format(ishell),"clmean.fits")
    clrec=read_cl(f1)
    l=arange(len(clrec))
    res=clrec-clt[l]
    plot(res,label=dens_type,alpha=0.7)
    i+=1

legend()
axhline(0,color='k',lw=0.5)
xlabel(r"$\ell$")
ylabel(r"residue")
xlim(0,500)
tight_layout()
