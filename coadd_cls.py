import os,glob,sys
from astropy.io import fits
from astropy.table import Table
import argparse

parser = argparse.ArgumentParser(description='coadd the cls from a directory')
parser.add_argument('clpattern', help='dir/cls*.fits like pattern')

args= parser.parse_args()

files=glob.glob(args.clpattern)
Ntot=len(files)
print(Ntot,files)
hdulist = fits.open(files[0])
h0=hdulist[0].header

h=hdulist[1].header
d=hdulist[1].data
names=d.dtype.names

cls=[]
[cls.append(d[n]) for n in names]
hdulist.close()

#now loop on remainder add
for file in files[1:]:
    hdulist = fits.open(file)
    d=hdulist[1].data
    for cl,n in zip(cls[1:],names[1:]):
        cl+=d[n]

#normalize
for cl in cls[1:]:
    cl/=Ntot

# write bintable with same header
t=Table(cls,names=names)
primary_hdu = fits.PrimaryHDU(header=h0)
hh=fits.BinTableHDU(data=t, header=h)
hdul = fits.HDUList([primary_hdu, hh])

d=os.path.dirname(files[0])
fout=os.path.join(d,"clmean.fits")
print("-> output={}".format(fout))
hdul.writeto(fout,overwrite=True)

