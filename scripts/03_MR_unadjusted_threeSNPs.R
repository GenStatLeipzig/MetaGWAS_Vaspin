#'---
#'title: "Unadjusted Mendelian Randomization of VASPIN and four disease data sets"
#'subtitle: "A little sensitivity analysis"
#'author: "Katrin Horn"
#'date: "`r Sys.Date()`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#'---

#' # Settings
rm(list=ls())
source("../SourceFile.R")
setwd(basicpath_scripts)


#' # Load data and prepare Mendel Rando
load("../data_MR/MR_snpData.RData")
snpData

# define SNPs to use for different phenotypes
homa.SNPs = c("rs7141073:94993744:C:T", "rs1956709:94967969:A:G", "rs12436152:94967748:C:T")
other.SNPs = c("rs7141073:94993744:C:T", "rs1956709:94967969:A:G", "rs61978267:94980896:C:T")

# sort snpData by position
is.data.table(snpData)
setkey(snpData, "pos")

#' # Mendelian Randomizations
#' ## Vaspin and HOMA
matched = match(homa.SNPs, snpData[,SNP])
MR.homa = mr_input(bx = snpData[matched, beta.vaspin], bxse = snpData[matched, SE.vaspin],
                   by = snpData[matched, beta.homa], byse = snpData[matched, SE.homa],
                   snps = snpData[matched, rsID], exposure = "Vaspin", outcome = "HOMA")
MR.homa.stat = mr_allmethods(MR.homa, method = "all")
MR.homa.stat

res.homa = MR.homa.stat@Values[4,]
mr_plot(MR.homa, interactive = T, labels = T)

#' ## Vaspin and Triglycerides
matched = match(other.SNPs, snpData[, SNP])
MR.tri = mr_input(bx = snpData[matched, beta.vaspin], bxse = snpData[matched, SE.vaspin],
                  by = snpData[matched, beta.tri], byse = snpData[matched, SE.tri],
                  snps = snpData[matched, rsID], exposure = "Vaspin", outcome = "Triglycerides")
MR.tri.stat = mr_allmethods(MR.tri, method = "all")
MR.tri.stat

res.tri = MR.tri.stat@Values[4,]
mr_plot(MR.tri, interactive = T, labels = T)

#' ## Vaspin and Cholesterol
MR.chol = mr_input(bx = snpData[matched, beta.vaspin], bxse = snpData[matched, SE.vaspin],
                   by = snpData[matched, beta.chol], byse = snpData[matched, SE.chol],
                   snps = snpData[matched, rsID], exposure = "Vaspin", outcome = "Cholesterol")
MR.chol.stat = mr_allmethods(MR.chol, method = "all")
MR.chol.stat

res.chol = MR.chol.stat@Values[4,]
mr_plot(MR.chol, interactive = T, labels = T)

#' ## Vaspin and LDL
MR.ldl = mr_input(bx = snpData[matched, beta.vaspin], bxse = snpData[matched, SE.vaspin],
                   by = snpData[matched, beta.ldl], byse = snpData[matched, SE.ldl],
                   snps = snpData[matched, rsID], exposure = "Vaspin", outcome = "Cholesterol")
MR.ldl.stat = mr_allmethods(MR.ldl, method = "all")
MR.ldl.stat

res.ldl = MR.ldl.stat@Values[4,]
mr_plot(MR.ldl, interactive = T, labels = T)

#' ## save results 
result.3SNPs = rbindlist(list(res.homa, res.tri, res.chol, res.ldl))
result.3SNPs[, Outcome := c("HOMA", "TG", "Chol", "LDL-Chol")]
setcolorder(result.3SNPs, c("Outcome", "Method", "Estimate", "Std Error", "95% CI ", " ", "P-value"))
setnames(result.3SNPs, c("Estimate", "Std Error"), c("Beta", "SE"))

# combine with results from adjustet MR analysis
result.5SNPs = read_xlsx(path = "../tables/Table_3_Results_of_Mendelian_Randomization.xlsx", sheet = 1)
result.5SNPs = as.data.table(result.5SNPs)
setnames(result.5SNPs, c("...6"), c(" "))

result.3SNPs[, Analysis := "correlation adjusted with 5 independent SNPs"]
result.5SNPs[, Analysis := "unadjusted with 3 independent SNPs"]

result = rbindlist(list(result.3SNPs, result.5SNPs), use.names = F)
setcolorder(result, c("Analysis", "Outcome", "Method", "Beta", "SE", "95% CI ", " ", "P-value"))
WriteXLS(result, ExcelFileName = "../tables/ST_MR_results.xlsx", SheetNames = "MR_results")
