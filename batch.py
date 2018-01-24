import numpy as np
from astropy.io import fits
import healpy as hp
import sys,os,glob
from tools import *


from pylab import *

class Catalog:
    def __init__(self,fname,rsd=True):
        data=(fits.open(fname)[1]).data
        print("reading {}".format(fname))
        self.ra=data['RA']
        self.dec=data['DEC']
        self.zrec=data['Z_COSMO']
        if rsd:
            self.zrec+=data['DZ_RSD']

    def get(self,zcut):
        w=np.logical_and(self.zrec>zcut[0],self.zrec<zcut[1])
        return self.zrec[w],self.ra[w],self.dec[w]

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

def get(dens_type=0,ishell=5,ngrid=512,rsd=True) :
    dirin=get_path(dens_type,ishell,ngrid)
    if rsd:
        f1=os.path.join(dirin,"clmean.fits")
    else:
        f1=os.path.join(dirin,"clmean_norsd.fits")
    assert os.path.exists(f1),"file does not exist:"+f1
    return hp.read_cl(f1)


def model(dens_type=0,ishell=5,ngrid=512,rsd=True):
    #camgal
    name="model/tophat_dens{:d}_ngrid{:d}".format(dens_type,ngrid)
    if not rsd:
        name+="_norsd"
    name+=".fits"
    assert os.path.exists(name),"file does not exist: "+name
    t=mrdfits(name,1)
    key='cl{}{}'.format(ishell-2,ishell-2)
    print("model={} :key={}".format(name,key))

    lt=t['ell']
    clt=t[key]
    clt[0]=0
    return lt,clt
    

def proj(dens_type=0,ishell=5,ngrid=512,nside=256,lmax=750,rsd=True,write=True):

    dirin=get_path(dens_type,ishell,ngrid)
    os.makedirs(dirin,exist_ok=True)
    zval=(0,0.1,0.2,0.3,0.4,0.5)

    zmax=zval[ishell]
    zmin=zval[ishell-1]
    
    files=glob.glob(os.path.join(dirin,"..","cat*.fits"))
    print("Analyzing: {}".format(dirin))
    print("shell #{}: z=[{},{}] ".format(ishell,zmin,zmax))
    
    print("there are {} files".format(len(files)))
    
    zcut=[zmin,zmax]
    
    print(" --> slice z=[{},{}] onto nside={} map".format(zcut[0],zcut[1],nside)) 
    npix=hp.nside2npix(nside)
    
    cls=[]
    cpt=1
    for file in files:
        #zrec,ra,dec=read_catalog(file,zcut,rsd)
        catalog=Catalog(file,rsd)
        zrec,ra,dec=catalog.get(zcut)
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
        dirout=os.path.join(dirin,"..","shell{:d}".format(ishell))
        os.makedirs(dirout,exist_ok=True)
        clname="clmean.fits"
        if not rsd :
            clname="clmean_norsd.fits"
        f1=os.path.join(dirout,clname)
        hp.write_cl(f1,clm,overwrite=True)
        covname="covmat.fits"
        if not rsd :
            clname="covmat_norsd.fits"
        f2=os.path.join(dirout,"covmat.fits")
        hdu=fits.ImageHDU(covmat)
        hdu.writeto(f2,overwrite=True)
        
        print("writing {} and {}".format(f1,f2))
        
    return clm,covmat


def ana(dens_type=0,ishell=5,ngrid=512,rsd=True):

    dens=("LogN","1LPT","2LPT","clipped")

    clrec=get(dens_type,ishell,ngrid,rsd)
    l=arange(len(clrec))
    lt,clt=model(dens_type,ishell,ngrid,rsd)

    #resize
    lmin=min(len(l),len(lt))
    clt=clt[0:lmin]
    clrec=clrec[0:lmin]

    #figure()
    plot(clt,'r',label=r"$C_\ell^{th}$")
    plot(clrec,'k',label=r"$<C_\ell^i>-SN$")
    plot(clrec-clt,label='residue')
    axhline(0,color='k',lw=0.5)
    legend()
    xlabel(r"$\ell$")
    ylabel(r"$C_\ell$")
    #ylim(-2e-5,8e-5)
    title(dens[dens_type]+"  (ngrid={})".format(ngrid))
    
    zval=(0,0.1,0.2,0.3,0.4,0.5)
    zmax=zval[ishell]
    zmin=zval[ishell-1]
    text(0.4, 0.9,r"$z\in [{},{}]$".format(zmin,zmax), transform=gca().transAxes,fontsize=12)

    tight_layout()

    show()
    return clrec-clt


def residues(ishell=5,ngrid=512,dens_types=(3,0,1,2),lmax=500):

    densnames=("LogN","1LPT","2LPT","Gaussclip")
    zval=(0,0.1,0.2,0.3,0.4,0.5)
    zmax=zval[ishell]
    zmin=zval[ishell-1]

    lt,clt=model(dens_type,ishell,ngrid)

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

