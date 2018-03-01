import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import py_cosmo_mad as csm
from scipy.integrate import quad

zmax=3.0
numz=1024

zarr=np.linspace(0,zmax,numz)
pcs=csm.PcsPar()
pcs.set_verbosity(1)
pcs.background_set(0.315,0.685,0.045,-1.,0.,0.7,2.7255)
gfarr=np.array([pcs.growth_factor(1./(1+z))/pcs.growth_factor(1.) for z in zarr])


zr,nzr=np.loadtxt("nz_red.txt" ,unpack=True);
nzrf=interp1d(zr,nzr,bounds_error=False,fill_value=0);
nzrarr=60*60*nzrf(zarr)
np.savetxt("NzRed.txt",np.transpose([zarr,nzrarr]))
bzrarr=1.5/gfarr
np.savetxt("BzRed.txt",np.transpose([zarr,bzrarr]))

zb,nzb=np.loadtxt("nz_blue.txt" ,unpack=True);
nzbf=interp1d(zb,nzb,bounds_error=False,fill_value=0);
nzbarr=60*60*nzbf(zarr)
np.savetxt("NzBlue.txt",np.transpose([zarr,nzbarr]))
bzbarr=1.0/gfarr
np.savetxt("BzBlue.txt",np.transpose([zarr,bzbarr]))
print " %d %d"%(quad(nzrf,0,2.5)[0]*4*np.pi*(180*60/np.pi)**2,quad(nzbf,0,2.5)[0]*4*np.pi*(180*60/np.pi)**2)

plt.figure()
plt.xlabel('$z$',fontsize=16)
plt.ylabel('$N(z)$',fontsize=16)
plt.plot(zarr,nzrarr,'r-')
plt.plot(zarr,nzbarr,'b-')

plt.figure()
plt.xlabel('$z$',fontsize=16)
plt.ylabel('$b(z)$',fontsize=16)
plt.plot(zarr,bzrarr,'r-')
plt.plot(zarr,bzbarr,'b-')
plt.show()
