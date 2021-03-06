library(Biobase)

#Load tSet
tset <- readRDS("data/TGGATES.rds")

#Quality Check 3 : Effect of azathioprine on the NRF2-associated gene module(#144)
##########################################################################################
drug <- "azathioprine"
conc <- 72.8

# drug <- "diclofenac"
# conc <- 400
# 
# drug <- "flutamide"
# conc <- 50
# 
# drug <- "isoniazid"
# conc <- 10000

#Genes of interest
genes <-c("ENSG00000159231","ENSG00000125037", "ENSG00000102393","ENSG00000181019", "ENSG00000109854", "ENSG00000164220")

#apply function to extract exprs matrix
values <- lapply(genes, function(gene){
#subset pehnodata for desired drugs
drug_subset <- subset(pData(tset@molecularProfiles$rna),drugid == drug,select=c(samplename, dose_level, individual_id,concentration))
#subset for ony high conc
drug_subset_high <- subset(drug_subset, concentration == conc)
#extracting exprs
assay <- exprs(tset@molecularProfiles$rna)
#subsetting exprs matrix
drug_subset$expression <- assay[gene,as.character(drug_subset$samplename)]
drug_subset_high$expression <- assay[gene,as.character(drug_subset_high$samplename)]
#ctrl rep
ctrlA <- na.omit(drug_subset$expression[select=c(drug_subset$dose_level == "Control" & drug_subset$individual_id=="1")])
ctrlB <- na.omit(drug_subset$expression[select=c(drug_subset$dose_level == "Control" & drug_subset$individual_id=="2")])

highA <- na.omit(drug_subset_high$expression[select=c(drug_subset_high$dose_level == "High" & drug_subset_high$individual_id=="1")])
highB <- na.omit(drug_subset_high$expression[select=c(drug_subset_high$dose_level == "High" & drug_subset_high$individual_id=="2")])

ctrl <- rowMeans(cbind(ctrlA, ctrlB))
high <- rowMeans(cbind(highA, highB))


normalised_vehicle <- (high-ctrl)*100
return(normalised_vehicle)
})

values <- as.data.frame(do.call(rbind,values))
colnames(values) <- c(2,8,24)
rownames(values) <- genes


time <- c(2,8,24)
legendnames <- c('CBR3','EMC3','GLA','NQO1','HTAT1P2','F2RL2')
colours <- c("purple", "red","green","violet","orange", "turquoise")

png("results/qc3_plot.png", width = 800, height = 600)
matplot(x = time, y = t(values)+100, col=colours, 
        pch=rep(21,ncol(values)), type=c("b"), lty=rep(1,ncol(values)), lwd=rep(5,ncol(values)),
        bg=colours,
        xlim=range(0,4,8,12,16,20,24),ylim=range(0,100,200,300,400,500),main="Azathioprine (Mod 325)",
        xlab="Time",ylab="mRNA Level (% Vehicle)")
par(font=2) 
legend(23.1,520,legend=legendnames, 
       col = colours,
       pch=rep(21,ncol(values)), pt.bg = colours,
       text.col = 'black',
       lty=rep(1,ncol(values)),lwd = rep(5, ncol(values)), cex=0.75,xjust = 0.5,bty = "n", adj = 0.25)
dev.off()
