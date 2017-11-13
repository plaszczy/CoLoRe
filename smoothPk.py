import numpy as np
import sys
from tools import *

fin=sys.argv[1]
R=float(sys.argv[2])
h=0.679
R/=h

kh,pkh=np.loadtxt(fin,unpack=True)
k,pk=remove_h(kh,pkh,h=h)
pkf=pk*np.exp(-(k*R)**2)
kh,pkfh=add_h(k,pkf,h=h)

np.savetxt(sys.stdout.buffer,np.transpose([kh,pkfh]))
