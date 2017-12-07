import sys
from pylab import *


ngrid=int(sys.argv[1])

h=0.679
Rsm=4/h

Lbox=4234.

l=Lbox/ngrid
print("Ngrid={} cell size l={}".format(ngrid,l))

R1=sqrt(Rsm**2+(l/sqrt(12))**2)

print("Rin={:5.2f} -> R={:5.2f}".format(Rsm,R1))
