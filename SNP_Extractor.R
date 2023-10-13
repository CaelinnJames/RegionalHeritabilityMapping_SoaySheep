args <- commandArgs(trailingOnly = TRUE)
print(args)
Blocks <- read.table(paste0("Haps/Haps_",args[2],".txt"),header=F)
Hap <- Blocks[args[1],1:3]
SNPs <- read.table("sheep_geno_imputed_Plates1to97_20220627.bim")
SNPs <- SNPs[which(SNPs$V1 == Hap$V1),]
SNPs <- SNPs[-which(SNPs$V4 < Hap$V2 | SNPs$V4 > Hap$V3),]
write.table(SNPs[,2],paste0("SNPs/Chr",args[2],"/regional_snps_",args[2],"_",args[1],".txt"),col.names = FALSE,row.names = FALSE,quote = FALSE)

