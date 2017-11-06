import numpy as np
from astropy.io import fits
import healpy as hp
import sys,os,glob
import matplotlib.pyplot as plt
import tools

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
assert(len(sys.argv)==3,"Usage: proj_cl dens_type numshell")

dens=sys.argv[1]
ishell=int(sys.argv[2])

dir=os.path.join("batch","dens_type"+dens)

zval=(0,0.1,0.2,0.3,0.4,0.5)

zmax=zval[ishell]
zmin=zval[ishell-1]


files=glob.glob(os.path.join(dir,"cat*.fits"))
print("Analyzing: {}".format(dir))
print("shell #{}: z=[{},{}] ".format(ishell,zmin,zmax))

print("there are {} files".format(len(files)))

zcut=[zmin,zmax]

nside=256
print(" --> slice z=[{},{}] onto nside={} map".format(zcut[0],zcut[1],nside)) 
npix=hp.nside2npix(nside)

cls=[]
for file in files:
    zrec,z0,ra,dec=read_fits(file,zcut)
    Nsamp=len(ra)
    nbar=Nsamp/(4*np.pi)
    print(" -> Nsamp={}, SN={}".format(Nsamp,1/nbar))
    mp=np.bincount(hp.ang2pix(nside,np.radians(90-dec),np.radians(ra)),minlength=npix)
    Nmean=mp.mean()
    map=mp.astype(float)/Nmean-1.
    #anafast
    cl=hp.anafast(map,lmax=500,iter=0,pol=False,use_weights=True,datapath=os.environ['HEALPIXDATA'])
    #remove SN
    cl-=1./nbar
    cls.append(cl)
    
clm=np.mean(cls,0)

dirout=os.path.join(dir,"shell{:d}".format(ishell))
os.makedirs(dirout,exist_ok=True)
    
f1=os.path.join(dirout,"clmean.fits")
hp.write_cl(f1,clm,overwrite=True)

f2=os.path.join(dirout,"cormat.fits")
covmat=np.cov(np.transpose(cls))
cormat=tools.cov2cor(covmat)

print("writing {} and {}".format(f1,f2))

tools.write_fits_image(cormat,f2)

tools.plotcor(cormat,vmin=-0.1,vmax=0.1)
