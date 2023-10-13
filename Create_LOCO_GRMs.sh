#!/bin/sh
#$ -cwd
#$ -l h_rt=2:00:00
#$ -pe mpi 16
#$ -l h_vmem=15G
#$ -o _o_files/
#$ -e _e_files/
#$ -t 1-26

. /etc/profile.d/modules.sh

module load igmm/apps/plink/1.90b4  
module load igmm/apps/dissect/1.15.2c  
module load python/3.4.3  


i=$SGE_TASK_ID

plink --bfile sheep_geno_imputed_Plates1to97_20220627 --sheep --not-chr ${i},X --make-bed --out LOCO_${i}  
  
mpirun -np 8 dissect.mpich  --make-grm --bfile LOCO_${i} --out LOCO_${i}

plink --bfile LOCO_${i} --sheep --exclude LOCO_${i}.badsnps  --make-bed --out LOCO_${i}_2  
  
mpirun -np 8 dissect.mpich  --make-grm --bfile LOCO_${i}_2 --out LOCO_${i}
