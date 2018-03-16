#!/bin/bash

#arguments
nargs=$#
if ! [ $nargs -eq 4 ]; then
echo "usage: run_batch colore.cfg proj.cfg nsims dirout/ "
exit
fi

COLORECONF=$1
SHELLSCONF=$2
NSIMS=$3
OUTDIR=$(readlink -e $4)

#
cp -f $COLORECONF $OUTDIR/$COLORECONF
grep -v file $SHELLSCONF > $OUTDIR/$SHELLSCONF
echo "filein=out_srcs_s1_0.fits" >> $OUTDIR/$SHELLSCONF
echo "fileout=cls.fits" >> $OUTDIR/$SHELLSCONF


range="1-10"
NSIMS_PER_BATCH=$(($NSIMS/10))

#where teh commad was run
HERE=$PWD

#default is running on 8 cores.
#if you wish to change it, set the environment variable NCORES to another value.
if [ -z "$NCORES" ] ; then 
NCORES=8
fi

##some checks
if ! [ -f $HERE/HEAD/amd64_sl6/CoLoRe ] ; then
echo "miss CoLoRe in $HERE/$HEAD/amd64_sl6"
exit
fi
if ! [ -f $HERE/HEAD/amd64_sl7/CoLoRe ] ; then
echo "miss CoLoRe in $HERE/$HEAD/amd64_sl7"
exit
fi
if ! [ -f $HERE/HEAD/amd64_sl6/proj ] ; then
echo "miss proj in $HERE/$HEAD/amd64_sl6"
exit
fi
if ! [ -f $HERE/HEAD/amd64_sl7/proj ] ; then
echo "miss proj in $HERE/$HEAD/amd64_sl6"
exit
fi

# go into OUTDIR
cd $OUTDIR

###################################################
cat > colore_run.sh  <<EOBATCH
#!/bin/bash

cd \$TMPDIR
df -h .


#for icc
COMPILERVARS_ARCHITECTURE=intel64
COMPILERVARS_PLATFORM=linux
source gcc_env.sh 5.2.0
source /usr/local/intel/icc/bin/iccvars.sh

#needs cmt for CMTCONFIG
source $LSSTLIB/CMT/v1r26/mgr/setup.sh

# copies locale execs
cp $HERE/HEAD/$CMTCONFIG/CoLoRe .
cp $HERE/HEAD/$CMTCONFIG/proj .

cp $OUTDIR/$COLORECONF .
cp $OUTDIR/$SHELLSCONF .

cp -r $HERE/data .

ls
#
export OMP_NUM_THREADS=$NCORES

#init ze random from task array
RANDOM=\${SGE_TASK_ID}

for isim in {1..${NSIMS_PER_BATCH}}
do
seed=\$RANDOM
 
awk -v seed=\$RANDOM  '{if (/_SEED_/) {print "seed = "seed}  else {print}}' $COLORECONF > tt.cfg

#this will create the catalog
./CoLoRe tt.cfg

#then projetc it into shells
#change seed 
grep -v "seed" $SHELLSCONF > proj.par
echo "seed=\$seed" >> proj.par

./proj proj.par

#copy& clean
cp cls.fits $OUTDIR/cls_\$seed.fits
#cp tt.cfg $OUTDIR/colore\$seed.cfg

rm -f *.fits
done

qstat -j \${JOB_ID} -nenv
EOBATCH

##############################################################################################
if [ -z "${QSUB_CMD}" ] ; then
QSUB_CMD="qsub -P P_$GROUP -pe multicores $NCORES -q mc_long -R y -j y -l sps=1"
fi
jobsub="${QSUB_CMD} -o $PWD -N $(basename $OUTDIR) -t $range colore_run.sh"

echo "about to run : $jobsub"
echo "OK? [y/n]"
read answer
if [ $answer != 'y' ] ; then
echo "exiting"
exit 1
fi

eval $jobsub

echo ">>>> Ouptuts in: $OUTDIR"

