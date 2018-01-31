import numpy as np
from astropy.io import fits
import healpy as hp
import sys,os,glob

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


def proj_all(dens_type=0,ngrid=512,nside=256,lmax=750,rsd=True,write=True):


    npix=hp.nside2npix(nside)
    dirin=os.path.join("batch","ngrid{}".format(ngrid),"dens_type{}".format(dens_type))
    files=glob.glob(os.path.join(dirin,"cat*.fits"))
    ncat=len(files)
    print("{} files in {}".format(len(files),dirin))
    if ncat==0:
        return
    zval=(0,0.1,0.2,0.3,0.4,0.5)

    cpt=1
    cls=np.zeros((4,lmax+1))
    for file in files:
        cat=Catalog(file,rsd)
        print("OK {}".format(cpt))
        for i in range(0,4):
            ishell=i+2
            zmax=zval[ishell]
            zmin=zval[ishell-1]
            zrec,ra,dec=cat.get([zmin,zmax])
            Nsamp=len(ra)
            nbar=Nsamp/(4*np.pi)
            mp=np.bincount(hp.ang2pix(nside,np.radians(90-dec),np.radians(ra)),minlength=npix)
            Nmean=mp.mean()
            map=mp.astype(float)/Nmean-1.
            #anafast
            cl=hp.anafast(map,lmax=lmax,iter=0,pol=False,use_weights=True,datapath=os.environ['HEALPIXDATA'])
            l=np.arange(len(cl))
            s2=np.sum((2*l+1)*cl)/(4*np.pi)
            print("\t -> shell={}: z in [{},{}]  Nsamp={:10.2f}M shot noise={:5.2}".format(ishell,zmin,zmax,Nsamp/1e6,np.sqrt(s2),1/nbar))
            #remove SN
            cl-=1./nbar
            cls[i,:]+=cl
        cpt+=1
    #normalize
    for i in range(0,4):
        cls[i,:]/=ncat

    if write:
        for i in range(0,4):
            ishell=i+2
            dirout=os.path.join(dirin,"shell{:d}".format(ishell))
            os.makedirs(dirout,exist_ok=True)
            clname="clmean.fits"
            if not rsd :
                clname="clmean_norsd.fits"
            f1=os.path.join(dirout,clname)
            print("writing {}".format(f1))
            hp.write_cl(f1,cls[i],overwrite=True)
                

    return cls


