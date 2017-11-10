
NSIM=1

dens_type=0

TEMPLATE="template_gal.cfg"
OUTDIR="batch/dens_type${dens_type}"
mkdir -p $OUTDIR

#loop

for ((i=0 ; i<$NSIM ; i++))
do
seed=$RANDOM
 
awk -v seed=$RANDOM -v dens=$dens_type '{if (/_SEED_/) {print "seed = "seed}  else if (/_DENS_TYPE_/) {print "dens_type="dens} else {print}}' template_gal.cfg > tt.cfg

./CoLoRe tt.cfg

#renommer
mv -f tt.cfg $OUTDIR/cat$seed.cfg
mv -f out_srcs_0.fits $OUTDIR/cat$seed.fits
#rm -f out*
done