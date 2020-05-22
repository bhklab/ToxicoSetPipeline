 drugmat <- drugMatrix
 drugmat@drug
 
saveRDS(drugmat, "./R/drugMatrixDMSO.rds")
 drugmat <- readRDS("./R/drugMatrixDMSO.rds")
 
View(drugmat@molecularProfiles$rna@phenoData@data)

## DIFFENTIAL GENE EXPRESION ANALYSIS
library(limma)

#### ALL DRUGS AT THE SAME TIME

# Average replicate values for each drug gene pair
# eset <- summarize_eset_replicates(tSet)
eset <- drugmat@molecularProfiles$rna
eset <- eset[, eset$dose_level %in% c("Control","High") & eset$duration == "24"]

# Samples, drugs, dose
targets <- pData(eset)[, c("samplename", "drugid", "dose_level")]

# Create combine design condition
drug <- factor(targets$drugid)
relevel(drug, "DMSO")
#dose <- factor(targets$dose_level, levels = c("Control", "High"))

dose #low and middle doses not showing? 

# Create design matrix for experimental design
design <- model.matrix(~drug)
#relevel 

head(design)

# Fit a linear model based on design matrix
fit <- lmFit(eset, design)

# Predict coefficients using emperical Bayes moderation of SE
# Generate a t-stat, moderated F-stat and log-odds of differential expression
stats <- eBayes(fit)

View(stats)

# Generate DE tables
resultsList <- lapply(colnames(stats), function(coef) {
  topTable(stats, coef)
})
table <- topTable(stats, coef = "drugXolegel", number= nrow(exprs(eset)), adjust.method = "BH")

volcanoplot(stats, coef = "drugXolegel", style = "p-value", highlight = 0, names = stats$genes, hl.col = "blue",
            xlab = "Log2 Fold Change", ylab = NULL, pch = 16, cex = 0.35)

?volcanoplot
#checkvignettes for stats 
?volcanoplot 

View(stats[,"drugXolegel"])

?volcanoplot
#tryconnectivitymapping 
#casestudy try, compare to tggates (volcanoplot) 

rownames(stats[,"drugXolegel"])

x <- as.vector(table$logFC)
y <- as.vector (table$P.Value)


# Evaluate hypothesis test for each drug gene pair
results <- decideTests(stats, coef = "drugXolegel")

View(results@.Data)

results

#GSEA from drugmatrix 










#### TWO DRUGS AT A TIME

eset2 <- eset[, eset$drugid %in% c("Xolegel", "DMSO")]

# Samples, drugs, dose
targets2 <- pData(eset2)[, c("samplename", "drugid", "dose_level")]

# Create combine design condition
drug2 <- factor(targets2$drugid)
dose2 <- factor(targets2$dose_level, levels = c("Control", "High"))

# Create design matrix for experimental design
design2 <- model.matrix(~0+drug2+drug2:dose2)

# Fit a linear model based on design matrix
fit2 <- lmFit(eset2, design2)

# Predict coefficients using emperical Bayes moderation of SE
# Generate a t-stat, moderated F-stat and log-odds of differential expression
stats2 <- eBayes(fit2)

stats2@.Data

# Generate DE tables
topTable(stats2, coef = "drug2Xolegel:dose2High")
#Error in fit$coefficients[, coef] : subscript out of bounds
topTable

# Evaluate hypothesis test for each drug gene pair
results2 <- decideTests(stats2, adjust.method = "fdr")


