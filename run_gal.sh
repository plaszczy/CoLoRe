
NSIM=100

dens_type=3
ngrid=1024

TEMPLATE="template_gal.cfg"
OUTDIR="batch/ngrid$ngrid/dens_type${dens_type}"
mkdir -p $OUTDIR

#loop

for ((i=0 ; i<$NSIM ; i++))
do
seed=$RANDOM
 
awk -v seed=$RANDOM -v dens=$dens_type -v ngrid=$ngrid '{if (/_SEED_/) {print "seed = "seed}  else if (/_DENS_TYPE_/) {print "dens_type="dens} else if (/_N_GRID_/) {print "n_grid="ngrid} else {print}}' template_gal.cfg > tt.cfg

#mpiexec -f ~/machinefile16 -n 8 ./CoLoRe tt.cfg
./CoLoRe tt.cfg

#renommer
cp tt.cfg $OUTDIR/cat$seed.cfg
mv -f out_srcs_s1_0.fits $OUTDIR/cat$seed.fits
rm -f out*
done
