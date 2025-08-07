#!/usr/bin/bash -l
#SBATCH -p epyc -c 32 -N 1 -n 1 --mem 300gb  --out logs/gtkdb.log --time 48:00:00
CPUS=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPUS=$SLURM_CPUS_ON_NODE
fi


module load workspace/scratch
module load gtdbtk
VERSION=226
if [ ! -f $VERSION ]; then
	ln -s /srv/projects/db/gtdbtk/$VERSION
fi

RUNDIR=$SCRATCH/$(date +%Y_%m_%d_%H_%M)
mkdir -p $RUNDIR
IN=$(realpath input)
# gotta just do this one at a time
for sample in $(ls $IN)
do
	for bin in $(ls $IN/$sample/bins/*.fa)
	do
		ln -s $bin $RUNDIR
	done
done

COUNT=$(ls -c $RUNDIR | wc -l)
OUTDIR=results/gtkdb
OUT=$OUTDIR/bincount_${COUNT}
mkdir -p $OUT
GTDBTK_DATA_PATH=$VERSION gtdbtk classify_wf --genome_dir $RUNDIR --out_dir $OUT -x fa --prefix gtdbtk --cpus $CPUS \
	--pplacer_cpus 8  --tmpdir $SCRATCH --mash_db $SCRATCH/mashdb.msh
