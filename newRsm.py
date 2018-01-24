import sys
from pylab import *


ngrid=int(sys.argv[1])

h=0.679
Rsm=4/h

Lbox=4234.

l=Lbox/ngrid


ladd=(l/sqrt(12))

print("Ngrid={}: cell size l={:5.2f} -> adding {:5.2f} Mpc".format(ngrid,l,ladd))

R1=sqrt(Rsm**2+ladd**2)

print("Rin={:5.2f} -> R={:5.2f}".format(Rsm,R1))
