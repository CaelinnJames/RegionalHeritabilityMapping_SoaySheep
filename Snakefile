def rhm_for_phenotype(phenotype):
  filenames = []
  for chrom in range(1, 27):
    trtfile = "Haps/Haps_{chrom}.txt".format(chrom=chrom)
    with open(trtfile, "r") as tf:
      num_haplotypes = len(tf.readlines())

    for hap in range(1, num_haplotypes + 1):
      filenames.append("{trait}/Output/Chr{chrom}/{trait}_{chrom}_{hap}.reml.tests".format(trait=phenotype, chrom=chrom, hap=hap))

  return filenames

def run_all_phenotypes(file):
  filenames = []
  with open (file, "r") as trtfile:
    trimmedfile = [next(trtfile) for _ in range(11)]
    for line in trimmedfile:
      currentLine = line.split(",")
      pheno = currentLine[0]
      filenames.extend(rhm_for_phenotype(pheno))
  return filenames

rule all:
  input: run_all_phenotypes("/exports/igmm/eddie/haley-soay/traits.txt")
   

rule precorrect:
  input:
    bfile=expand("/exports/igmm/eddie/haley-soay/GRMs/LOCO_P1_97/LOCO_{{chrom}}_2.{ext}", ext=["bim", "bed", "fam"]),
    pheno="{trait}/Precorrecting/{trait}Phenotypes.txt",
    covar="{trait}/Precorrecting/{trait}Assoc.txt",
    qcovar="{trait}/Precorrecting/{trait}QAssoc.txt",
    random="{trait}/Precorrecting/{trait}Random.txt"
  output:
    blups="{trait}/Precorrecting/Chr{chrom}/{trait}_chr{chrom}.E.blup.indiv"
  params:
    bfile="/exports/igmm/eddie/haley-soay/GRMs/LOCO_P1_97/LOCO_{chrom}_2", 
    blups="{trait}/Precorrecting/Chr{chrom}/{trait}_chr{chrom}",
    rt="-l h_rt=04:00:00",
    mem="-l h_vmem=4G",
    pe="-pe mpi 16 -R y"
  priority: 40
  benchmark:
    "_benchmark/precorrect/precorrect{trait}_{chrom}.txt"
  shell:
    """
    mkdir -p {wildcards.trait}/Precorrecting/Chr{wildcards.chrom}
    declare -A dict

    while IFS=',' read -r key value ignore; do
      dict[$key]=$value
    done < /exports/igmm/eddie/haley-soay/traits.txt
    
    module load igmm/mpi/gcc/mpich/3.1.4
    module load igmm/apps/dissect/1.15.2c

  mpirun -np 8 dissect.mpich  --reml \
    --bfile {params.bfile} \
    --pheno {input.pheno} \
    --pheno-col 1 \
    --covar {input.covar} \
    --qcovar {input.qcovar} \
    --random-effects {input.random} \
    --random-effects-cols ${{dict[{wildcards.trait}]}} \
    --indiv-blup \
    --out {params.blups}
    """

rule prep_precorrecting:
  input:
    blups="{trait}/Precorrecting/Chr{chrom}/{trait}_chr{chrom}.E.blup.indiv"
  output:
    tempblup=temp("{trait}/Precorrecting/Chr{chrom}/{trait}_chr{chrom}_temp.txt"),
    precorrected="{trait}/Precorrecting/Chr{chrom}/{trait}_chr{chrom}_Precorrected.txt"
  params:
    pe="",
    rt="-l h_rt=00:10:00",
    mem="-l h_vmem=5G"  
  priority: 30
  benchmark:
    "_benchmark/prep_precorrecting/prep_precorrecting_{trait}_{chrom}.txt"
  shell:
    """
    awk '{{print $2,$2 , $3}}' {input.blups} > {output.tempblup} 
    awk NR\>1 {output.tempblup} > {output.precorrected}
    
    """
    
rule get_hap_snps:
  output:
    "SNPs/Chr{chrom}/regional_snps_{chrom}_{hap}.txt"
  params:
    rt="-l h_rt=00:10:00",
    pe="",
    mem="-l h_vmem=5G"  
  priority: 20
  benchmark:
    "_benchmark/get_hap_snps/get_hap_snps_{chrom}_{hap}.txt"
  shell:
    """
    mkdir -p SNPs
    mkdir -p SNPs/Chr{wildcards.chrom}
    module load igmm/apps/R/4.1.3
    Rscript /exports/igmm/eddie/haley-soay/3_RegionalHaplotypeMapping/SNP_Extractor.R {wildcards.hap} {wildcards.chrom}    
    """
    
rule run_hap_analysis:
  input:
    bgen="/exports/igmm/eddie/haley-soay/Phasing/SHAPEIT4/Chr{chrom}.bgen",
    pheno="{trait}/Precorrecting/Chr{chrom}/{trait}_chr{chrom}_Precorrected.txt",
    region="SNPs/Chr{chrom}/regional_snps_{chrom}_{hap}.txt"
  output:
    analysis="{trait}/Output/Chr{chrom}/{trait}_{chrom}_{hap}.reml.tests"
  params:
    output="{trait}/Output/Chr{chrom}/{trait}_{chrom}_{hap}",
    bgen="/exports/igmm/eddie/haley-soay/Phasing/SHAPEIT4/Chr{chrom}",
    rt="-l h_rt=14:00:00",
    mem="-l h_vmem=8G",
    pe="-pe sharedmem 2 -R y"      
  priority: 5
  benchmark:
    "_benchmark/run_hap_analysis/run_hap_analysis_{trait}_{chrom}_{hap}.txt"
  shell:
    """
    mkdir -p {wildcards.trait}/Output
    mkdir -p {wildcards.trait}/Output/Chr{wildcards.chrom}
    export MKLROOT=/exports/applications/apps/SL7/intel/parallel_studio_xe_2015_update5/composer_xe_2015.5.223/mkl 
    export LD_LIBRARY_PATH=${{MKLROOT}}/lib/intel64:$LD_LIBRARY_PATH

    module load igmm/mpi/gcc/mpich/3.1.4
    mpirun -np 8 /exports/igmm/eddie/haley-lab/Eilidh/tests/dissect/src/dissect  --reml \
      --bgen {params.bgen} \
      --pheno {input.pheno}  \
      --pheno-col 1 \
      --extract {input.region} \
      --fit-snp-grm \
      --fit-haplotype-grm \
      --min-snps-region 2 \
      --min-overlap-snps 0 \
      --reml-maxit 30 \
      --run-all-models \
      --out {params.output}
 
    """
