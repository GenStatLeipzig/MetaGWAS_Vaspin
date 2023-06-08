#'---
#'title: "Adjusted Mendel Rando of VASPIN and three disease data sets"
#'author: "Katrin Horn"
#'date: "`r Sys.Date()`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#'---

#' # presettings
rm(list=ls())
source("../SourceFile_aman.R")
setwd(basicpath_scripts)

#' # Load data and prepare Mendel Rando
# load prepared SNP data
load("../data_MR/MR_snpData.RData")
snpData

# load covariance matrix (due to low LD SNP partners)
load("../data_MR/MR_covMatrix.RData")
covMatrix

# define SNPs to use for different phenotypes
homa.SNPs = c("rs7141073:94993744:C:T", "rs1956709:94967969:A:G", "rs4905216:95000960:G:C", "rs12436152:94967748:C:T", "rs17094914:94975583:C:G")
other.SNPs = c("rs7141073:94993744:C:T", "rs1956709:94967969:A:G", "rs4905216:95000960:G:C", "rs61978267:94980896:C:T", "rs73338689:94973878:G:C")

# sort snpData by position
is.data.table(snpData)
setkey(snpData, "pos")

# check SNP order in data objects
table(colnames(covMatrix) == snpData[, SNP])

#' # Mendelian Randomization
#' ## Vaspin and HOMA adjusted with covariance matrix
m.SNP = match(homa.SNPs, snpData[, SNP])
m.cor = match(homa.SNPs, colnames(covMatrix))

MR.homa = mr_input(bx = snpData[m.SNP, beta.vaspin], bxse = snpData[m.SNP, SE.vaspin],
                   by = snpData[m.SNP, beta.homa], byse = snpData[m.SNP, SE.homa],
                   snps = snpData[m.SNP, rsID], exposure = "Vaspin", outcome = "HOMA", correlation = covMatrix[m.cor, m.cor])
MR.homa.stat = mr_allmethods(MR.homa, method = "all")
MR.homa.stat

result.homa = MR.homa.stat@Values[4,]
mr_plot(MR.homa, interactive = T, labels = T)

#' ## Vaspin and Triglycerides 
m.SNP = match(other.SNPs, snpData[, SNP])
m.cor = match(other.SNPs, colnames(covMatrix))

MR.tri = mr_input(bx = snpData[m.SNP, beta.vaspin], bxse = snpData[m.SNP, SE.vaspin],
                   by = snpData[m.SNP, beta.tri], byse = snpData[m.SNP, SE.tri],
                   snps = snpData[m.SNP, rsID], exposure = "Vaspin", outcome = "Triglycerides", correlation = covMatrix[m.cor, m.cor])
MR.tri.stat = mr_allmethods(MR.tri, method = "all")
MR.tri.stat

result.tri = MR.tri.stat@Values[4,]
mr_plot(MR.tri, interactive = T, labels = T)

#' ## Vaspin and Cholesterol adjusted
m.SNP = match(other.SNPs, snpData[, SNP])
m.cor = match(other.SNPs, colnames(covMatrix))

MR.chol = mr_input(bx = snpData[m.SNP, beta.vaspin], bxse = snpData[m.SNP, SE.vaspin],
                   by = snpData[m.SNP, beta.chol], byse = snpData[m.SNP, SE.chol],
                   snps = snpData[m.SNP, rsID], exposure = "Vaspin", outcome = "Cholesterol", correlation = covMatrix[m.cor, m.cor])
MR.chol.stat = mr_allmethods(MR.chol, method = "all")
MR.chol.stat

result.chol = MR.chol.stat@Values[4,]
mr_plot(MR.chol, interactive = T, labels = T)

#' ## Vaspin and LDL adjusted
m.SNP = match(other.SNPs, snpData[, SNP])
m.cor = match(other.SNPs, colnames(covMatrix))

MR.ldl = mr_input(bx = snpData[m.SNP, beta.vaspin], bxse = snpData[m.SNP, SE.vaspin],
                   by = snpData[m.SNP, beta.ldl], byse = snpData[m.SNP, SE.ldl],
                   snps = snpData[m.SNP, rsID], exposure = "Vaspin", outcome = "Cholesterol", correlation = covMatrix[m.cor, m.cor])
MR.ldl.stat = mr_allmethods(MR.ldl, method = "all")
MR.ldl.stat

result.ldl = MR.ldl.stat@Values[4,]
mr_plot(MR.ldl, interactive = T, labels = T)

#' # Prepare results table for paper draft
result = rbindlist(list(result.homa, result.tri, result.chol, result.ldl))
result[, Outcome := c("HOMA", "TG", "Chol", "LDL-Chol")]
setcolorder(result, c("Outcome", "Method", "Estimate", "Std Error", "95% CI ", " ", "P-value"))
setnames(result, c("Estimate", "Std Error"), c("Beta", "SE"))
WriteXLS(result, ExcelFileName = "../tables/Table_3_Results_of_Mendelian_Randomization.xlsx", SheetNames = "Result_Mendel_Rando") 


#' # Prepare supplementary table with single SNP statistics
snpData[, pValue.vaspin := 10^(logP.vaspin*(-1))]
snpData[, pValue.tri := 10^(logP.tri*(-1))]
snpData[, pValue.chol := 10^(logP.chol*(-1))]
snpData[, pValue.ldl := 10^(logP.ldl*(-1))]
colsOut = c("rsID", "EAF", "EAF.homa", "EA.homa", "OA.homa", "EAF.tri", "EA.tri", "OA.tri", "EAF.chol", "EA.chol", "OA.chol", 
            "EAF.ldl", "EA.ldl", "OA.ldl", "logP.vaspin", "logP.homa", "logP.tri", "logP.chol", "logP.ldl")
snpData[, get("colsOut") := NULL]

# calculate F statistics for vaspin
snpData[, FStat.vaspin := (beta.vaspin/SE.vaspin)^2]

# update column order
setcolorder(snpData, c("SNP", "chr", "pos", "N", "MAF", "effect_allele", "other_allele",
                       "beta.vaspin", "SE.vaspin", "FStat.vaspin", "pValue.vaspin", "beta.homa", "SE.homa", "P.homa", 
                       "beta.tri", "SE.tri", "pValue.tri", "beta.chol", "SE.chol", "pValue.chol", "beta.ldl", "SE.ldl", "pValue.ldl"))
snpData[, N := round(N)]

# remove unused SNP data entries 
setNA.other = c("beta.tri", "SE.tri", "pValue.tri", "beta.chol", "SE.chol", "pValue.chol", "beta.ldl", "SE.ldl", "pValue.ldl")
setNA.homa = c("beta.homa", "SE.homa", "P.homa")
snpData[SNP == "rs12436152:94967748:C:T", get("setNA.other") := NA]
snpData[SNP == "rs17094914:94975583:C:G", get("setNA.other") := NA]
snpData[SNP == "rs61978267:94980896:C:T", get("setNA.homa") := NA]
snpData[SNP == "rs73338689:94973878:G:C", get("setNA.homa") := NA]

# save SNP data in excel file
WriteXLS(snpData, ExcelFileName = "../tables/ST_MR_single_SNP_statistics.xlsx", SheetNames = "MR_SNPs")

#' # Make supplementary plot for different MR methods
plotA.l = mr_plot(MR.homa, interactive = F, labels = F)
plotA.r = mr_plot(mr_allmethods(MR.homa, method = "main", iterations = 50), orientate = TRUE)
plotB.l = mr_plot(MR.tri, interactive = F, labels = F)
plotB.r = mr_plot(mr_allmethods(MR.tri, method = "main", iterations = 50), orientate = TRUE)
plotC.l = mr_plot(MR.chol, interactive = F, labels = F)
plotC.r = mr_plot(mr_allmethods(MR.chol, method = "main", iterations = 50), orientate = TRUE)
plotD.l = mr_plot(MR.ldl, interactive = F, labels = F)
plotD.r = mr_plot(mr_allmethods(MR.ldl, method = "main", iterations = 50), orientate = TRUE)

# make the plots look nicer and with matching text sizes
huebschScatter = function(plot1, x, y, title){
  plot1 + theme_light(base_size = 10) + xlab(x) + ylab(y)+ ggtitle(title)
}
plotA.l2 = huebschScatter(plotA.l, "Effect sizes vaspin","Effect sizes HOMA", "A")
plotB.l2 = huebschScatter(plotB.l, "Effect sizes vaspin","Effect sizes TG", "B")
plotC.l2 = huebschScatter(plotC.l, "Effect sizes vaspin","Effect sizes Chol", "C")
plotD.l2 = huebschScatter(plotD.l, "Effect sizes vaspin","Effect sizes LDL-Chol", "D")

huebschOther = function(plot2, x, y){
  plot2 + theme_light(base_size = 10) + xlab(x) + ylab(y)
}
plotA.r2 = huebschOther(plotA.r, "Effect sizes vaspin","Effect sizes HOMA")
plotB.r2 = huebschOther(plotB.r, "Effect sizes vaspin","Effect sizes TG")
plotC.r2 = huebschOther(plotC.r, "Effect sizes vaspin","Effect sizes Chol")
plotD.r2 = huebschOther(plotD.r, "Effect sizes vaspin","Effect sizes LDL-Chol")

tiff(filename = "../figures/SF_MR_results_several_methods.tiff", width = 1600, height = 2000, res=250, compression = 'lzw')
plot_grid(plotA.l2, plotA.r2,
          plotB.l2, plotB.r2,
          plotC.l2, plotC.r2,
          plotD.l2, plotD.r2,
          ncol=2, nrow=4, align ="h", rel_widths = c(0.8,1))
dev.off()
