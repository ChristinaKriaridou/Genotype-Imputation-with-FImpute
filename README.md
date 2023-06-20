# Genotype-Imputation-with-FImpute
Imputation for genomic selection and GWAS practical as part of the AquaIMPACT training course on Genomic Innovations for Aquaculture Breeding (27th June 2023 - Block 2).

## 1. DATA
For this tutorial, we are going to use the first chromosome of an Atlantic salmon dataset previously published here: https://doi.org/10.1186/s12864-015-2117-9.

## 2. PRE-IMPUTATION FILTERING OF STUDY GENOTYPES
Before you perform an imputation run with your study genotypes, you should filter the data to remove low-quality variants and individuals, as these can lower the accuracy of the results. Standard GWAS quality control filters are usually sufficient to prepare a dataset for imputation.
The dataset in this practical was filtered with PLINK v.1.9 (Purcell et al., 2007). Individuals with just one of their two parents genotyped or > 20% missing genotypes were excluded from the analysis. SNPs with > 10% missing genotypes; significant deviation from Hardy–Weinberg Equilibrium (P-value < 10−6); MAF < 0.05; or Mendelian error rates > 10% were also excluded from subsequent imputation analyses.
After filtering 606 individuals and 78,035 SNPs remained in the dataset. The first chromosome that we are going to use consists of 4,424 SNPs.

## 3. MASK A CERTAIN NUMBER OF SNPs PROPORTIONAL TO CHR LENGTH AND EQUALLY SPACED
To mask a certain number of SNPs on chromosome 1 (that we are going to later impute), we will use the code in this repository:
https://github.com/Roslin-Aquaculture/Select-SNPs-to-generate-low-density-panels.

The SNPs selected are equally distributed across the genome and proportionally to chromosome length. Additionally, this script always selects the first and the last SNP of each chromosome.
We are going to run the script to mask the SNPs of just chromosome 1 in this practical, but you can find out how many SNPs you need to mask to achieve a certain low density panel for all your chromosomes in the tutorial in the above link. 

```
#Check your working directory path
getwd()
#Set your working directory path (where you downloaded the "Masking" folder of this practical)
setwd("C:/Users/.../Imputation_tutorial/Masking")

#Load the R packages below
library(DescTools)
library(Siccuracy)
library(data.table)
library(tidyverse)
library(FNN)

#Number of SNPs to select from chromosome 1 for a given low-density panel
#Read the file with the chromosome start and end bp position. This file has the information just for
#chromosome 1. You can modify the script and add information for the rest of the chromosomes in the
#file according to your species.
chr_length<- read.table("scot_salmon_chr_length.txt", header = TRUE)
head(chr_length)
#Calculate total map length
total_length<- sum(chr_length$total_chr_distance)
#Divide the desired number of SNPs (e.g. 71) in the LD panel with the total map length
#The number of SNPs (71) was calculated from a file containing information for all chromosomes for a
#low density panel of 1000 SNPs equally spaced throughout the genome and according to chromosome length.
index<- 71/total_length
#Multiply the length of each chr with the index above and round 
n_snp_per_chr<- round(index*chr_length$total_chr_distance, 0)
sum(n_snp_per_chr)    #should sum up to the desired SNP number we want to keep for each chromosome  

# load SNP ID, chromosome and position file
snpmap <-read.table("scot_salmon_snp_info.txt",header=T,sep="\t")
head(snpmap)

#split snpmap to chromosomes (the output from split() is a list)
chr_list = split(snpmap, snpmap$Chr)

#Take each chromosome and divide it's length into equally distanced parts (equidistant positions),
#according to the number of snps we want to select (LD)
#These equidistant theoretical positions will be then used to find the nearest real position on the map  
b_chr<- list()
for (i in 1:length(chr_length$chr)) {
  b_chr[[i]]<- round(seq(chr_length$first_SNP[i], chr_length$last_SNP[i], length.out = n_snp_per_chr[i]), digits = 0)
}

#Check that the length of the list for each chromosome equals the number of SNPs we want to keep
#for this chromosome
for (i in 1:length(chr_length$chr)){
  print(length(b_chr[[i]]))
} 

######### CHROMOSOME 1 ##########

#Extract the BPPos and SNP ID column from the dataframes list 
chr_bppos = lapply(chr_list, "[", , "BPPos")
chr_snpid = lapply(chr_list, "[", , "SNPID")

#Select the SNPs we are going to keep (all the other SNPs will be set to missing)
#Find the nearest row of the theoritical positions to the real positions on the map file
s1<- c()
for (i in 1:length(b_chr[[1]])) {
  s1[i]<- knnx.index(chr_bppos[["1"]], b_chr[[1]][i], k=1)
  chr_bppos[["1"]][s1[i]] <- 0
  s1<- sort(s1)
  print(length(unique(s1)))
}


#The genotype files must be split per chromosome (just genotypes in the format of 0/1/2) and individuals
#should be reordered so that all the parents are at the beginning (rows 1-87)
#Read the genotype files per chromosome for all individuals
geno_file_names <- list.files(pattern= glob2rx("salmon_chr*"))
geno_file_names
salmon_pf_chr <- lapply (geno_file_names, read.table, sep=" ", header=F, row.names=NULL)
#number of individuals
nrow(salmon_pf_chr[[1]])
#number of SNPs
ncol(salmon_pf_chr[[1]])

#Mask all the SNPs of the offspring (rows 86-606 of the genotype file) that are not specified in s1
#and write them in a separate file for each chromosome. 
masking <- list()
masked_salmon_chr<- list()
m<- list()
for (i in 1:length(chr_length$chr)) {
  masking[[i]] = list(86:606, setdiff(1:length(chr_bppos[[i]]), s1))
  
  snpfile <- tempfile()
  salmon_geno<- data.matrix(salmon_pf_chr[[i]])
  class(salmon_geno)
  write.snps(salmon_geno, snpfile)
  fn <- tempfile()
  
  # This function runs through a single file and masks specified columns of specified rows. 
  masked_salmon_chr[[i]] <- mask_SNPs(snpfile, fn, masking=masking[[i]], na=5)
  m[[i]] <- read.snps(fn)
  #change the path according to where you want to save the masked file
  myfile <- file.path("C:/Users/.../Imputation_tutorial", paste0("masked_salmon_chr", i, ".txt"))
  write.table(m[[i]], file = myfile, sep = " ",
              row.names = FALSE, col.names = FALSE)
}
```
After we masked the SNPs in the first chromosome of the Atlantic salmon dataset (71 SNPs remained unmasked out of 4424), we remove the spaces between the genotypes. Then we add the individual IDs, the chip info for each individual and the header to create the inpute file for FImpute.
```
sed 's/ //g' masked_salmon_chr1.txt | paste ids_chip_606.txt - | sed '1i ID Chip Genotypes' > masked_salmon_chr1_fimpute.txt
```

## 4. FIMPUTE INPUT FILE FORMATS
### i. Map file
The map file contains information for each SNP present in the dataset. The columns of this file are: (1) SNP name, (2) chromosome number, (3) bp position, (4) chip position. For example, chip 2 may be a lower-density panel and SNP number 1 will start at the equivalent of SNP 20 in chip 1. This will change based on the common SNPs you have on the panels.

![image](https://github.com/ChristinaKriaridou/Genotype-Imputation-with-FImpute/assets/74717500/f6697f3d-05d5-4d84-9d6d-162952ab5510)

### ii. Genotype file
The genotype file contains the SNPs in columns and individuals in rows. The columns are: (1) ID, (2) Chip Number, (3) Genotypes (coded as allele dosages: 0,1,2,5=missing). There are no spaces present in the genotypes (saves space).
![image](https://github.com/ChristinaKriaridou/Genotype-Imputation-with-FImpute/assets/74717500/3dd41fc4-c7cd-407a-9d61-55831947eb21)

### iii. Pedigree file
The pedigree file provides family information for each individual in the genotype file. You need a file with (1) Animal ID, (2) Sire ID, (3) Dam ID, (4) Sex (M/F/U).

![image](https://github.com/ChristinaKriaridou/Genotype-Imputation-with-FImpute/assets/74717500/cf8006af-b79f-4779-997f-7ddfa249cfca)

### iv. Parameter file settings
This is a simple control file with options and file names. There are many more options in the documentation.
![image](https://github.com/ChristinaKriaridou/Genotype-Imputation-with-FImpute/assets/74717500/f005244a-007c-42fd-9a66-c997c27d688b)

## 5. RUNNING THE APPLICATION
FImpute [control filename] -o

If control file name is not specified, the program will prompt the user to enter it. Option –o forces the program to overwrite output folder if it already exists. 
```
FImpute control_file.txt -o
```
## 6. OUTPUT FILES
Some of the important output files you can check are:

**report.txt**: detailed report on the steps carried out by the software.

**genotypes_imp.txt**: this file contains the imputed genotypes for all individuals.

**parentage_test.txt**: this file is produced by the parentage test option in the control file and checks for parentage errors. It contains information about individual call rate, sire call rate, dam call rate, no. Mendelian inconsistencies, no. loci compared and a possible match for the parents.

**stat_snp.txt**: reports statistics on SNPs: SNP ID, chromosome number, positions, call frequencies, missing rate and minor allele frequency. 

**stat_snp_imp.txt**: reports statistics on SNPs after imputation.

**stat_anim.txt**: Reports statistics on individuals' genotypes: ID, chip number, call frequencies, homozygosity and missing rate. Missing calls are ignored for statistics on homozygosity and calls 0, 1 and 2. 

**stat_anim_imp.txt**: Reports statistics on individuals' genotypes after imputation

## 7. IMPUTATION ACCURACY
To calculate imputation accuracy (i.e. correlation between two genotype matrices) we are going to use a package called Siccuracy. Load the package in R with: 
```
library(Siccuracy)
```

If you haven’t installed the Siccuracy package use devtools to install directly from github:
```
library(devtools)
install_github(repo=’stefanedwards/Siccuracy’, ref = “master”)
```

After imputation we remove the header, the parents and add spaces between the imputed genotypes to prepare the file for imputation accuracy calculation. I already did that for this tutorial so you can skip this step.
```
cut -f 3 genotypes_imp.txt | sed '1d' | sed 's/./& /g' | awk '{print NR,$0}' | sed '1,86d' > genotypes_imp_for_accur_calc.txt
```
For the next steps in R, we will need the imputed file with the offspring (genotypes_imp_for_accur_calc.txt) and the true genotypes (genotypes_true_for_accur_calc.txt).
```
#Set your working directory with the path to the folder where you have downloaded the files of "Imputation_accuracy_calc"
setwd("C:/Users/.../Imputation_tutorial/Imputation_accuracy_calc")

#Calculate correlation between real and imputed data after imputing the first chromosome
#of the salmon genotypes with FImpute v.3
# standardized method
accur_fimpute <- imputation_accuracy('genotypes_true_for_accur_calc.txt', 'genotypes_imp_for_accur_calc.txt', na= 5)
round(accur_fimpute[["matcor"]], digits = 3)

# Pearson correlation between true and imputed genotypes per animal and SNP.
corr_animals<- as.data.frame(accur_fimpute[["animals"]][["cors"]])
write.table(corr_animals, file = "salmon_corr_per_animal.txt", sep = " ",
            row.names = TRUE, col.names =FALSE, quote = FALSE)
colnames(corr_animals)[1]  <- "correlation"

corr_snps<- as.data.frame(accur_fimpute[["snps"]][["cors"]])
write.table(corr_snps, file = "testsalmon_corr_per_snp.txt", sep = " ",
            row.names = TRUE, col.names =FALSE , quote=FALSE)
colnames(corr_snps)[1]  <- "correlation"
```

## 8. PLOT IMPUTATION ACCURACY PER SNP AND INDIVIDUAL
Plot correlation per SNP:
```
#Load the ggplot library
library(ggplot2)

#Read the "snp_info.txt" file from the  "FImpute_output_files" folder (change the path in the command below accordingly)
snp_info <- read.delim("C:/Users/S1899268/Desktop/Imputation_tutorial/FImpute_output_files/snp_info.txt")
#Add the correlation column to the snp_info data frame
snp_info$correlation <- corr_snps$correlation
#save separately the imputed and non imputed SNPs
non_imputed_snps<- snp_info[snp_info$chip_2 > 0,]
imputed_snps<- snp_info[snp_info$chip_2 == 0,]

#Plot the correlation for the imputed (black dots) and non imputed SNPs (blue dots) with ggplot
theme_set(theme_bw())
figure <- ggplot(imputed_snps, aes(BPPos, correlation)) + #geom_smooth(method = "loess", span=0.025, se=FALSE, size=1) + 
  geom_point(size=0.5) + ggtitle(paste("Imputation of chromosome 1 with FImpute v.3")) + xlab ("Position (bp)") + ylab ("Correlation")+
  theme(plot.title = element_text(hjust = 0.5))
figure2 <- figure + geom_point(data=non_imputed_snps, col="blue", size=1.5) + 
  ylim(c(0, 1.01))

#Save the graph in a pdf
pdf("plot_correlation_per_snp.pdf", width=10,height=7)
print(figure2)
dev.off()
```


## 9. POST-IMPUTATION FILTERING
It is standard practice to perform additional filtering once imputation run has completed. This is mainly to remove poorly imputed variants that might behave badly in association tests. Here, after imputation we filtered for MAF (< 0.05) with PLINK v.1.9 (Purcell et al., 2007).


