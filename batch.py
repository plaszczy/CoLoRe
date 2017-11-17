import numpy as np
from astropy.io import fits
import healpy as hp
import sys,os,glob
from tools import *


from pylab import *



def read_catalog(fname,zcut,rsd=True) :
    data=(fits.open(fname)[1]).data
    print("reading {}".format(fname))
    
    ra=data['RA']
    dec=data['DEC']
    zrec=data['Z_COSMO']
    if rsd:
        zrec+=data['DZ_RSD']

    #cut
    w=np.logical_and(zrec>zcut[0],zrec<zcut[1])
    
    return zrec[w],ra[w],dec[w]

####
def get_path(dens_type=0,ishell=5,ngrid=512):
    return os.path.join("batch","ngrid{}".format(ngrid),"dens_type{}".format(dens_type),"shell{}".format(ishell))

####

def proj(dens_type=0,ishell=5,ngrid=512,nside=256,lmax=750,rsd=True,write=True):

    dir=get_path(dens_type,ishell,ngrid)
    zval=(0,0.1,0.2,0.3,0.4,0.5)

    zmax=zval[ishell]
    zmin=zval[ishell-1]
    
    
    files=glob.glob(os.path.join(dir,"..","cat*.fits"))
    print("Analyzing: {}".format(dir))
    print("shell #{}: z=[{},{}] ".format(ishell,zmin,zmax))
    
    print("there are {} files".format(len(files)))
    
    zcut=[zmin,zmax]
    
    print(" --> slice z=[{},{}] onto nside={} map".format(zcut[0],zcut[1],nside)) 
    npix=hp.nside2npix(nside)
    
    cls=[]
    cpt=1
    for file in files:
        zrec,ra,dec=read_catalog(file,zcut,rsd)
        Nsamp=len(ra)
        nbar=Nsamp/(4*np.pi)
        print(" {} -> Nsamp={}, SN={}".format(cpt,Nsamp,1/nbar))
        mp=np.bincount(hp.ang2pix(nside,np.radians(90-dec),np.radians(ra)),minlength=npix)
        Nmean=mp.mean()
        map=mp.astype(float)/Nmean-1.
        #anafast
        cl=hp.anafast(map,lmax=lmax,iter=0,pol=False,use_weights=True,datapath=os.environ['HEALPIXDATA'])
        #remove SN
        cl-=1./nbar
        cls.append(cl)
        cpt+=1
        
    clm=np.mean(cls,0)
    covmat=np.cov(np.transpose(cls))
        
    if write:
        dirout=os.path.join(dir,"..","shell{:d}".format(ishell))
        os.makedirs(dirout,exist_ok=True)
        f1=os.path.join(dirout,"clmean.fits")
        hp.write_cl(f1,clm,overwrite=True)
        f2=os.path.join(dirout,"covmat.fits")
        hdu=fits.ImageHDU(covmat)
        hdu.writeto(f2,overwrite=True)
        
        print("writing {} and {}".format(f1,f2))
        
    return clm,covmat


def ana(dens_type=0,ishell=5,ngrid=512):

    dens=("LogN","1LPT","2LPT","Gaussclip")

    dir=get_path(dens_type,ishell,ngrid)

    f1=os.path.join(dir,"shell{}".format(ishell),"clmean.fits")
    clrec=hp.read_cl(f1)
    l=arange(len(clrec))

    #truth="clR4_shell{:0d}_bordersout.txt".format(ishell)
    #lt,clt=loadtxt(truth,unpack=True)
    t=mrdfits("colore.fits",1)
    key='cl{}{}'.format(ishell-2,ishell-2)
    lt=t.l
    clt=t[key]
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
    ylim(-2e-5,8e-5)
    title(dens[dens_type])
    
    zval=(0,0.1,0.2,0.3,0.4,0.5)
    zmax=zval[ishell]
    zmin=zval[ishell-1]
    text(0.4, 0.9,r"$z\in [{},{}]$".format(zmin,zmax), transform=gca().transAxes,fontsize=12)

    tight_layout()

    show()


    return clrec-clt


def residues(ishell=5,ngrid=256,dens_types=(3,0,1,2),lmax=500):

    densnames=("LogN","1LPT","2LPT","Gaussclip")
    zval=(0,0.1,0.2,0.3,0.4,0.5)
    zmax=zval[ishell]
    zmin=zval[ishell-1]

    truth="clR4_shell{:0d}_bordersout.txt".format(ishell)
    
    lt,clt=loadtxt(truth,unpack=True)
    clt[0]=0

    figure()
    for i in dens_types :
        f1=os.path.join("batch","ngrid{}".format(ngrid),"dens_type{}".format(i),"shell{}".format(ishell),"clmean.fits")
        clrec=hp.read_cl(f1)
        l=arange(len(clrec))
        res=clrec-clt[l]
        plot(res,label=densnames[i],alpha=0.8)
        
    plot(clt[l]/100,'k--',label=r"$C_\ell^{true}/100$")
 
    
    legend(loc='lower right')
    title(r"$z\in [{},{}]$".format(zmin,zmax))
    axhline(0,color='k',lw=0.5)
    xlabel(r"$\ell$")
    ylabel(r"residue")
    xlim(0,lmax)
    tight_layout()
    show()

