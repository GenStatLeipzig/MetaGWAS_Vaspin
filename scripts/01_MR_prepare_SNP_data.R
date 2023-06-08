# In this script the data for the Mendelian Randomization is prepared.
# At first the summary statistics of the MetaGWAS of VASPIN are filtered for the five SNPs that will be the instruments of the MR analysis.
# Two SNPs needed to be added because they were missing in the HOMA reference data set.
# As a second step the SNP data for the above mentioned seven SNPs is added from four different disease reference data sets 
# (HOMA, Triglycerides, total Cholesterol and LDL-Cholesterol).


#####
# 0. Setting
#####
rm(list=ls())
source("../SourceFile.R")
setwd(basicpath_scripts)


#####
# 1. define list of my SNPs and get their data from VASPIN MetaGWAS
#####
# This list contains the five top-SNPs from chromosome 14 locus of the MetaGWAS.
# Two of this SNP are not available in the HOMA data, so two proxy SNPs of those were added to this list for the HOMA analysis.
mySNPs = c("rs7141073:94993744:C:T", "rs1956709:94967969:A:G", "rs4905216:95000960:G:C", "rs61978267:94980896:C:T", 
           "rs73338689:94973878:G:C", "rs12436152:94967748:C:T", "rs17094914:94975583:C:G") 
mySNPs.rs = c("rs7141073", "rs1956709", "rs4905216", "rs61978267", "rs73338689", "rs12436152", "rs17094914")

snpData = fread("../sumStats/MetaGWAs_VASPIN_ALL_summary_stats.txt.gz")
snpData = snpData[is.element(rsID, mySNPs), ]
cols2keep = c("rsID", "chromosome", "position", "N", "effect_allele", "other_allele", "EAF", "MAF", "beta", "SE", "logP")
colsOut = setdiff(colnames(snpData), cols2keep)
snpData[, get("colsOut") := NULL]
setnames(snpData, c("rsID", "chromosome", "position", "N", "effect_allele", "other_allele", "EAF", "MAF", "beta", "SE", "logP"), 
         c("SNP", "chr", "pos", "N", "effect_allele", "other_allele", "EAF", "MAF", "beta.vaspin", "SE.vaspin", "logP.vaspin"))
dummy = strsplit(snpData[, SNP], split = ":")
snpData[, rsID := sapply(dummy, function(x) return(x[1]))]
setcolorder(snpData, c("SNP", "rsID", "chr", "pos", "N", "EAF", "MAF", "effect_allele", "other_allele", "beta.vaspin", "SE.vaspin", "logP.vaspin"))
snpData


#####
# 2. add data of disease reference data sets
#####
# Homa
# Discovery sample description: 37,037 European ancestry individuals
dat.homa = fread("../reference_data/20081858-GCST005179-EFO_0004501.h.tsv.gz")
table(is.element(mySNPs.rs, dat.homa[, variant_id]))

#proxy-SNPs used for the missing two SNPs 
#rs61978267 -> rs12436152, distance , D' 0.9832, R^2 0.9354, Allele: C=C, T=T -> position 3 in Credible Set
#rs73338689 -> rs17094914, distance , D' 1.0, R^2 0.977, Allele: G=C C=G -> position 2 in Credible Set

matched = match(mySNPs.rs, dat.homa[, variant_id])
dat.homa = dat.homa[matched,]
setnames(dat.homa, c("hm_rsid", "hm_effect_allele","hm_other_allele", "hm_effect_allele_frequency", "beta", "standard_error", "p_value"), 
         c("rsID.homa", "EA.homa", "OA.homa", "EAF.homa", "beta.homa", "SE.homa", "P.homa"))
cols2keep = c("rsID.homa", "EAF.homa", "EA.homa", "OA.homa", "beta.homa", "SE.homa", "P.homa")
colsOut = setdiff(colnames(dat.homa), cols2keep)
dat.homa[, get("colsOut") := NULL]
setcolorder(dat.homa, cols2keep)
dat.homa[4, rsID.homa := "rs61978267"]
dat.homa[5, rsID.homa := "rs73338689"]
dat.homa[, rsID.homa2 := rsID.homa]

snpData = merge(x = snpData, y = dat.homa, by.x = "rsID", by.y = "rsID.homa", all = T)
snpData[, logP.homa := -log10(P.homa)]

# Triglycerides
dat.tri = fread("../reference_data/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")
dat.tri = dat.tri[is.element(rsID, mySNPs.rs), ]
cols2keep = c("rsID.tri", "EAF.tri", "EA.tri", "OA.tri", "beta.tri", "SE.tri", "logP.tri")
setnames(dat.tri, c("rsID", "POOLED_ALT_AF", "ALT", "REF", "EFFECT_SIZE", "SE", "pvalue_neg_log10"), cols2keep)
colsOut = setdiff(colnames(dat.tri), cols2keep)
dat.tri[, get("colsOut") := NULL]
setcolorder(dat.tri, cols2keep)
dat.tri[, rsID.tri2 := rsID.tri]
dat.tri

snpData = merge(x = snpData, y = dat.tri, by.x = "rsID", by.y = "rsID.tri", )
snpData

# total Cholesterol
dat.chol = fread("../reference_data/TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")
dat.chol = dat.chol[is.element(rsID, mySNPs.rs), ]
cols2keep = c("rsID.chol", "EAF.chol", "EA.chol", "OA.chol", "beta.chol", "SE.chol", "logP.chol")
setnames(dat.chol, c("rsID", "POOLED_ALT_AF", "ALT", "REF", "EFFECT_SIZE", "SE", "pvalue_neg_log10"), cols2keep)
colsOut = setdiff(colnames(dat.chol), cols2keep)
dat.chol[, get("colsOut") := NULL]
setcolorder(dat.chol, cols2keep)
dat.chol[, rsID.chol2 := rsID.chol]
dat.chol

snpData = merge(x = snpData, y = dat.chol, by.x = "rsID", by.y = "rsID.chol", )
snpData

# LDL
dat.ldl = fread("../reference_data/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")
dat.ldl = dat.ldl[is.element(rsID, mySNPs.rs), ]
cols2keep = c("rsID.ldl", "EAF.ldl", "EA.ldl", "OA.ldl", "beta.ldl", "SE.ldl", "logP.ldl")
setnames(dat.ldl, c("rsID", "POOLED_ALT_AF", "ALT", "REF", "EFFECT_SIZE", "SE", "pvalue_neg_log10"), cols2keep)
colsOut = setdiff(colnames(dat.ldl), cols2keep)
dat.ldl[, get("colsOut") := NULL]
setcolorder(dat.ldl, cols2keep)
dat.ldl[, rsID.ldl2 := rsID.ldl]
dat.ldl

snpData = merge(x = snpData, y = dat.ldl, by.x = "rsID", by.y = "rsID.ldl", )
snpData


#####
# 3. check data
#####
#check for allele switches
table(snpData[, effect_allele] == snpData[, EA.homa])
table(snpData[, other_allele] == snpData[, OA.homa])

table(snpData[, effect_allele] == snpData[, EA.tri])
table(snpData[, other_allele] == snpData[, OA.tri])
snpData[which(snpData[, effect_allele] != snpData[, EA.tri]), ]

table(snpData[, effect_allele] == snpData[, EA.chol])
table(snpData[, other_allele] == snpData[, OA.chol])
snpData[which(snpData[, effect_allele] != snpData[, EA.chol]), ]

table(snpData[, effect_allele] == snpData[, EA.ldl])
table(snpData[, other_allele] == snpData[, OA.ldl])
snpData[which(snpData[, effect_allele] != snpData[, EA.ldl]), ]

#switch alleles of rs4905216 for Triglyceride, LDL and total Cholesterol
snpData[rsID == "rs4905216", EAF.tri := 1 - EAF.tri]
snpData[rsID == "rs4905216", EA.tri := "C"]
snpData[rsID == "rs4905216", OA.tri := "G"]
snpData[rsID == "rs4905216", beta.tri := (-1) * beta.tri]

snpData[rsID == "rs4905216", EAF.chol := 1 - EAF.chol]
snpData[rsID == "rs4905216", EA.chol := "C"]
snpData[rsID == "rs4905216", OA.chol := "G"]
snpData[rsID == "rs4905216", beta.chol := (-1) * beta.chol]

snpData[rsID == "rs4905216", EAF.ldl := 1 - EAF.ldl]
snpData[rsID == "rs4905216", EA.ldl := "C"]
snpData[rsID == "rs4905216", OA.ldl := "G"]
snpData[rsID == "rs4905216", beta.ldl := (-1) * beta.ldl]

table(snpData[, rsID] == snpData[, rsID.homa2])
table(snpData[, rsID] == snpData[, rsID.tri2])
table(snpData[, rsID] == snpData[, rsID.chol2])
table(snpData[, rsID] == snpData[, rsID.ldl2])

snpData[, c("rsID.homa2", "rsID.tri2", "rsID.chol2", "rsID.ldl2") := NULL]
snpData

#####
# 4. save SNP data object
#####
save(snpData, file = "../data_MR/MR_snpData.RData")
