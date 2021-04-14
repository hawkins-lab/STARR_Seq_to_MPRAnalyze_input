library(MPRAnalyze)
setwd("/net/hawkins/vol1/home/chrhsu/proj/STARR_Seq_to_MPRAnalyze/pipeline_output/13_final_mpranalyze_input")
colAnnot_rs <- read.csv("colAnnot_rs_dnagr5_rnagr1_barcodeallelic.csv", row.names = 1)

dnaCount_rs <- read.csv("DNACounts_rs_dnagrt5_rnagrt1.csv",row.names = 1)

rnaCount_rs <- read.csv("RNACounts_rs_dnagrt5_rnagrt1.csv", row.names = 1)

DNA_counts_rs <- as.matrix(dnaCount_rs)

RNA_counts_rs <- as.matrix(rnaCount_rs)

Th1_Rs_allele <- MpraObject(dnaCounts=DNA_counts_rs, rnaCounts=RNA_counts_rs, , colAnnot = colAnnot_rs, controls = NA_integer_, BPPARAM = NULL)

objrs_depth <- estimateDepthFactors(Th1_Rs_allele, lib.factor = "batch", which.lib = "dna", depth.estimator = "totsum")
objrs_depth <- estimateDepthFactors(Th1_Rs_allele, lib.factor = "batch", which.lib = "rna", depth.estimator = "totsum")

objrs_allele <- analyzeComparative(objrs_depth, dnaDesign = ~ batch + barcode_allelic + allele, rnaDesign = ~ batch + allele, reducedDesign = ~ allele)
objrs_barcodeallelic <- analyzeComparative(objrs_depth, dnaDesign = ~batch + barcode_allelic + allele, rnaDesign = ~ allele, reducedDesign = ~ allele)

#objrs_batch <- analyzeComparative(objrs_depth, dnaDesign = ~ batch + barcode + allele, rnaDesign = ~ batch + allele, reducedDesign = ~ batch)
results_rs_barcodeallele <- testLrt(objrs_allele)
results_rs_barcodeallelic <- testLrt(objrs_barcodeallelic)
#results_rs_batch <- testLrt(objrs_batch)
object_allele <- saveRDS(objrs_allele, "objrsdna5rna1_alleledepthcorrectedbarcodeallelic.rds")
object_barcodeallelic <- saveRDS(objrs_barcodeallelic, "objrsdna5rna1_alleledepthcorrectedbarcodeallelicrnamodelallele.rds")
#object_batch <- saveRDS(objrs_batch, "objrsdna5rna1_batchdepthcorrected.rds")
write.csv(results_rs_barcodeallele, "./Th1dna5rna1_rsmpranalyzereddesignallelewithdepthnormalizationbarcodeallelic.csv", row.names = TRUE)
write.csv(results_rs_barcodeallelic, "./Th1dna5rna1_rsmpranalyzereddesignallelewithdepthnormalizationbarcodeallelicrnamodelallele.csv", row.names = TRUE)
#write.csv(results_rs_batch, "./Th1dna5rna1_rsmpranalyzereddesignbatchwithdepthnormalization.csv", row.names = TRUE)
obj_scale <- analyzeComparative(objrs_depth, rnaDesign = ~ allele, reducedDesign = ~ 1, mode = "scale")
res.scale <- testLrt(obj_scale)
write.csv(res.scale , "./resscale.csv", row.names=TRUE)
