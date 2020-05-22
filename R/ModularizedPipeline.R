library(ToxicoGx)
library(Biobase)
library(affy)
library(affyio)
library(BiocManager)
library(rat2302rnensgcdf)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(data.table)
#devtools::install_github("bhklab/ToxicoGx", ref = "master")
#library(ToxicoGx)
#drugMatrix
#saveRDS(drugMatrix, "/Users/parwaiznijrabi/desktop/BHKlab/drugMatrix121219.rds")

create_phenoData <- function(species=c("Rat"), verbose = TRUE) { 
  if (verbose) {message("Creating phenoData object...")} 
  phenoData <- read.csv("./data/phenodataDMSOcontrols.csv", stringsAsFactors = FALSE)
  phenoData$Dose_Unit <- "μM"
  if(verbose) {message("phenoData object created!")}
  return(phenoData)
}
phenoData <- create_phenoData("Rat")


create_exprsData <- function(species=c("Rat"), phenoData, verbose = TRUE) {
  if (verbose) {message("Creating eset object")}
  if(species == "Rat"){
  }
  #esetNorm <- just.rma(filenames = celFiles, verbose = TRUE, cdfname = "rat2302rnensgcdf")
  esetNorm <- readRDS("./esetNorm.rds")
  eset <- esetNorm
  storageMode(eset)<-"environment"
  eset <-subset(eset, substr(rownames(eset@assayData$exprs), 0, 4) != "AFFX")
  storageMode(eset)<-"lockedEnvironment"
  annotation(eset)<-"rna"
  if (verbose) {message("eset object created!")}
  return(eset)
}
eset <- create_exprsData("Rat")


create_featureData <- function(species=c("Rat"), eset, verbose = TRUE){
  if (verbose) {message("Creating featureData object...")}
  if (species == "Rat"){
    ensembl <- useMart("ensembl")
    datasets <- listDatasets(ensembl)
    ensembl = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host="uswest.ensembl.org",ensemblRedirect = FALSE)
    storageMode(eset) <- "environment"
    affxrows <- rownames(eset@assayData$exprs)
    rownames(eset@assayData$exprs) <- substr(rownames(eset@assayData$exprs), 1, nchar(affxrows)-3)
    #saveRDS(affxrows, file = "./data/probesRatVitro1.rds") 
    CELgenes <- affxrows
    CELgenes1 <- gsub(".at", " ", CELgenes)
    results <-getBM(attributes=c("external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id","external_transcript_name","ensembl_transcript_id"), filters = "ensembl_gene_id",values=CELgenes1, mart=ensembl,checkFilters = TRUE)
    uniqueB <- results[!duplicated(results$ensembl_gene_id),]
    CELnotB <- unique(CELgenes1) [!unique(CELgenes1) %in% uniqueB$ensembl_gene_id]
    names(uniqueB) <- c("gene_name", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id")
    finalFeature <- uniqueB
    
    finalFeature$BEST <- NA
    names(finalFeature) <- c("Symbol", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id", "BEST")
    rownames(finalFeature) <- finalFeature$gene_id
    finalFeature$gene_id
    geneid1 <- finalFeature$gene_id
    
    for (i in 1:length(geneid1)) {
      geneid1[i] = paste(geneid1[i], "at", sep="_")
    }
    geneid1
    finalFeature$gene_id <- geneid1
    finalFeature$gene_id
    finalFeature[,1]
    rownames(finalFeature) = finalFeature$gene_id
    
    if(verbose) {message("featureData object created!")}
    return(finalFeature)
    
  }
}
create_featureData("Rat")


create_Expressionset <- function(species=c("Rat"), eset, verbose = TRUE){
  if (verbose) {message("Creating expressionset...")}
  if (species == "Rat"){
phenoData 
featureData <- readRDS("./featureData.rds")
pData(eset)<-phenoData
fData(eset)<-featureData
storageMode(eset)<-"lockedEnvironment"
saveRDS(eset1, file="/Users/parwaiznijrabi/Desktop/ExpressionSet.rds")
  }
}


######################toxicoset constructor function ######################

tSet <- ToxicoSet("drugmatrix_hepatocyte",
                  molecularProfiles=list("rna"=ExpressionSet),
                  cell=cell,
                  drug=drug,
                  sensitivityInfo=NULL,
                  sensitivityRaw=NULL,
                  sensitivityProfiles=NULL,
                  curationDrug=curationDrug,
                  curationCell=curationCell,
                  curationTissue=curationTissue,
                  datasetType=c("perturbation"),
                  verify = TRUE)

#all(rownames(fData(esetFINAL))%in% rownames(exprs(esetFINAL))) 
#length(intersect(pData(esetFINAL)[,"cellid"] , cell$cellid))
#length(intersect(pData(esetFINAL)[,"cellid"] , cell$cellid))
#length(intersect(pData(esetFINAL)[,"cellid"] , cell$cellid))
#View(esetFINAL@phenoData@data)
#View(esetFINAL@assayData$exprs)
#View(esetFINAL@featureData@data)

# TESTS #
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(unique(sensitivityInfo$cellid) , cell$cellid))