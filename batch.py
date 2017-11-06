import numpy as np
from astropy.io import fits
import healpy as hp
import sys,os,glob

def read_fits(fname,zcut) :
    data=(fits.open(fname)[1]).data
    print("reading {}".format(fname))
    
    ra_arr=data['RA']
    dec_arr=data['DEC']
    z0_arr=data['Z_COSMO']
    rsd_arr=data['DZ_RSD']
    
    zrec=z0_arr+rsd_arr
    w=np.logical_and(zrec>zcut[0],zrec<zcut[1])
    
    z0=z0_arr[w]
    zrec=zrec[w]
    ra=ra_arr[w]
    dec=dec_arr[w]
    
    return zrec,z0,ra,dec

####

def proj(dens_type=0,ishell=5,nside=256,lmax=500):

    dir=os.path.join("batch","dens_type{}".format(dens_type))

    zval=(0,0.1,0.2,0.3,0.4,0.5)

    zmax=zval[ishell]
    zmin=zval[ishell-1]
    
    
    files=glob.glob(os.path.join(dir,"cat*.fits"))
    print("Analyzing: {}".format(dir))
    print("shell #{}: z=[{},{}] ".format(ishell,zmin,zmax))
    
    print("there are {} files".format(len(files)))
    
    zcut=[zmin,zmax]
    
    print(" --> slice z=[{},{}] onto nside={} map".format(zcut[0],zcut[1],nside)) 
    npix=hp.nside2npix(nside)
    
    cls=[]
    cpt=1
    for file in files:
        zrec,z0,ra,dec=read_fits(file,zcut)
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
    
    dirout=os.path.join(dir,"shell{:d}".format(ishell))
    os.makedirs(dirout,exist_ok=True)
    
    f1=os.path.join(dirout,"clmean.fits")
    hp.write_cl(f1,clm,overwrite=True)
    
    f2=os.path.join(dirout,"covmat.fits")
    covmat=np.cov(np.transpose(cls))
    hdu=fits.ImageHDU(covmat)
    hdu.writeto(f2,overwrite=True)
    
    print("writing {} and {}".format(f1,f2))
    
    return clm,covmat
