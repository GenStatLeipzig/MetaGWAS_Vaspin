# MetaGWAS Vaspin
Details of the Mendelian Randomization Analysis

**Last Updated: 08/06/2023**

Supporting code for the following draft:

* **Working title**: Genetic dissection of serum vaspin highlights its causal role in lipid metabolism
* **Short title**: MetaGWAS Vaspin

We conducted a meta-analysis of genome-wide association studies for serum-vaspin from six independent cohorts (N=7,446). Potential functional variants of vaspin were included in Mendelian Randomization (MR) analyses to assess possible causative chains between vaspin and HOMA or lipid traits. To further validate the MR analyses we analyzed data from GTEx, treated db/db mice with vaspin and measured serum lipids.   

For more information, please contact Jana Breitfeld (jana.breitfeld@medizin.uni-leipzig.de)

# Source File

If you want to reproduce our results, you will need to customize a source file, indicating
* path to R library (please use [R Version 4.x](https://cran.r-project.org/), all necessary packages are listed in the source file)
* path to data (summary statistics will be available after publication on zenodo)
* path to disease reference data sets (see below for links we used)

# Mendelian Randomization scripts

1) **Preparation of SNP data for MR**: extraction of relevant SNP data from VASPIN MetaGWAS results and available reference data for investigated diseases
2) **Adjusted MR with five SNPs as instruments**: covariance matrix adjusted MR of five SNPs (including two pairs of low LD partners)
3) **Unadjusted MR with three SNPs as instruments**: unadjusted MR as sensitivity analysis to support findings of adjusted MR
not included) **Calculation of covariance matrix for adjusted MR**: genetic data of LIFE-Adult used, can not be shared

# Tables
1) **Main Table 3** showing results of the adjusted MR with five SNPs (--> see script 02_MR_adjusted_fiveSNPs.R)

# Supplement Material
1) **MR plots** showing results of several MR methods for comparison (--> see script 02_MR_adjusted_fiveSNPs.R)
2) **Table MR single SNP statistics** (--> see script 02_MR_adjusted_fiveSNPs.R)
3) **Table MR results for adjusted and unadjusted MR** (--> see script 03_MR_unadjusted_threeSNPs.R)

# Links to disease reference data used 
* HOMA_IR: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005179/harmonised/ (PMID: 20081858)
* Triglycerides, Cholesterol and LDL-Cholesterol: https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/ (PMID: 34887591)
