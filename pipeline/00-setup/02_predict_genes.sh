#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 8 --mem 16gb --out logs/prodigal.%a.log

module load prodigal

INPUT=input
OUTPUT=results/genes

mkdir -p $OUTPUT
CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "cannot run without a number provided either cmdline or --array in sbatch"
  exit
fi
DIR=$(ls $INPUT | sed -n ${N}p)
NAME=$(basename $DIR)
mkdir -p $OUTPUT/$NAME
parallel -j $CPU prodigal -q -i {} -f gff -o $OUTPUT/$NAME/{/.}.gff -a $OUTPUT/$NAME/{/.}.aa.fa -d $OUTPUT/$NAME/{/.}.cds.fa -p single ::: $(ls $INPUT/$DIR/bins/*.fa)

pigz -f $OUTPUT/$NAME/*.gff $OUTPUT/$NAME/*.cds.fa
