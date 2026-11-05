#!/bin/bash
# usage: showcase.sh <design-prefix> <name> ; simulate one tube per pool, analyze, copy to igv site
set -e
S=/tmp/claude-0/-home-user-my-amplicon-project/f58d2bc1-8662-564d-9a2b-eb2192724c3d/scratchpad
export PATH=$S/mamba3/envs/amp3/bin:$PATH
cd /home/user/my_amplicon_project
P=$1; N=$2; G=common_refs/ecoli_genome.fna
mkdir -p $S/sim/$N
for bed in ${P}_pool_*.bed; do
  k=$(basename $bed .bed | sed 's/.*_pool_/pool/')
  python analysis/generate_simulated_fastq.py --genome $G --bed $bed --output-dir $S/sim/$N --sample-name ${N}_$k --coverage 300 --read-length 150 --seed 1 >/dev/null 2>&1
  python analysis/analyze.py --reads "$S/sim/$N/${N}_${k}_R{1,2}.fastq.gz" --bed $bed --primers ${P}_${k#pool}_trimming.fasta --genome $G --outdir $S/an/$N >/dev/null 2>&1 || python analysis/analyze.py --reads "$S/sim/$N/${N}_${k}_R{1,2}.fastq.gz" --bed $bed --primers ${P}_pool_${k#pool}_trimming.fasta --genome $G --outdir $S/an/$N 2>&1 | tail -5
done
ls $S/an/$N
