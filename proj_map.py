import numpy as np
from astropy.io import fits
import healpy as hp
import sys,os,glob
from tools import *

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


file=sys.argv[1]

nside=256
zcut=[0.4,0.5]
print("projecting slice z=[{},{}] onto nside={} map".format(zcut[0],zcut[1],nside)) 
npix=hp.nside2npix(nside)
zrec,z0,ra,dec=read_fits(file,zcut)
mp=np.histogram(hp.ang2pix(nside,np.radians(90-dec),np.radians(ra)),
                bins=npix,range=[0,npix])[0]
plt.figure()
hist_plot(mp,bins=50)
plt.show(block=True)
Nmean=mp.mean()
map=mp.astype(float)/Nmean-1.
print("#samples={}M,sigma={}".format(len(ra)/1e6,map.std()))
d=os.path.dirname(file)
f=os.path.basename(file)
fout=os.path.join(d,"map{}_".format(nside)+f)
print("writing "+fout)
hp.write_map(fout,map,overwrite=True)
