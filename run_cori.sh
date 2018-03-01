
#Output paths
predir="/global/cscratch1/sd/plaszczy"
rundir=${predir}"/colore"

ingal=galLSST10Y.cfg

#Number of nodes
nnod=32
nthreads=64
#Time limit in minutes
timelim=10
#Queue to use
#which_partition="regular"
which_partition="debug"

#ecriture cfg
parfile=${rundir}/param_colore.cfg
awk -v rundir=$rundir '{ if (/prefix_out/){print "\"prefix_out="rundir"/out\""} else {print}}' $ingal > $parfile

#run script
runfile=${rundir}/run.sh

cat > ${runfile} << EOF
#!/bin/bash -l
#SBATCH --partition ${which_partition}
#SBATCH --nodes ${nnod}
#SBATCH --time=00:${timelim}:00
#SBATCH --job-name=CoLoRe_SP
#SBATCH -C haswell

export OMP_NUM_THREADS=${nthreads}
srun -n ${nnod} -c ${nthreads} $PWD/CoLoRe $parfile

EOF

cat $runfile
ask "launch?"
if [ $? -eq 0 ] ; then
sbatch ${runfile}
echo "output in $rundir"
fi


