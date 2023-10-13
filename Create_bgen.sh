#!/bin/sh

#$ -cwd  
#$ -o _o_files/
#$ -e _e_files/
#$ -l h_vmem=20G
#$ -l h_rt=10:00:00
#$ -t 1-26 
#$ -pe sharedmem 16 

. /etc/profile.d/modules.sh

module load igmm/apps/plink/2.00
module load igmm/apps/shapeit4/4.2.2
module load igmm/apps/bcftools/1.10.2


BFILE=sheep_geno_imputed_Plates1to97_20220627
OUT=Phasing/Chr${SGE_TASK_ID}


plink2 --bfile ${BFILE} \
	--export vcf bgz id-paste=iid \
  --chr ${SGE_TASK_ID} \
  --sheep \
	--threads ${NSLOTS} \
	--out ${BFILE}_${SGE_TASK_ID} 


bcftools index ${BFILE}_${SGE_TASK_ID}.vcf.gz

shapeit4.2 --input ${BFILE}_${SGE_TASK_ID}.vcf.gz \
	--region ${SGE_TASK_ID} \
	--out ${OUT}.vcf.gz \
	--seed 1234 \
	--thread ${NSLOTS} \
	--sequencing


 plink2 --vcf ${OUT}.vcf.gz \
 --export bgen-1.2 'id-paste=iid'\
  --sheep \
  --threads ${NSLOTS} \
  --out ${OUT}
  
rm ${BFILE}_${SGE_TASK_ID}.vcf.gz*
