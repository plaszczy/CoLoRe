from pylab import *
import os
from healpy import *

l,cl=loadtxt("colore_cl.txt",unpack=True)

t,ct=loadtxt("colore_ctheta.txt",unpack=True)

figure()
plot(l,cl)
ylim(ymin=0)

print("sigma={}".format(sqrt(sum((2*l+1)*cl)/(4*pi))))
