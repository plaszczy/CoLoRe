
from tools import *


i=1

#dndz=5000
t=mrdfits("outputs/bench0/clmean.fits",1)

#red
tr,hr=mrdfits("outputs/bench0/clmean_red.fits",1,header=True)

colnames=tr.names
n=colnames[i+1]

cl1=t[n]

cl2=tr[n]
#remove SN
S2=mean(cl2[600:])
cl2-=S2

figure()
subplot(2,1,1)

title(n)
plot(cl1,label='const')
plot(cl2,'r',label='red')
xlim(1,200)
legend()
ax0()
ylabel("Cl")

rcl=cl2/cl1

subplot(2,1,2)
plot(rcl[1:200],'r')
xlim(1,200)
ax1()
ylabel("Cl(red)/Cl(const)")


#SHOTNOISE
Ntot=hr["ngal{}".format(i)]

print("||{:04.2e} || {:04.2e} || {:04.2e} ||".format(4*pi/Ntot,S2,S2/(4*pi/Ntot)))
