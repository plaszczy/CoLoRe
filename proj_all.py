
import batch

ngrid=1024
dens=0

rsd=True
for shell in range(2,6):
    batch.proj(dens_type=dens,ishell=shell,ngrid=ngrid,rsd=rsd)

rsd=False
for shell in range(2,6):
    batch.proj(dens_type=dens,ishell=shell,ngrid=ngrid,rsd=rsd)
