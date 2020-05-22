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

create_phenodata_DM <- function(species=c("Rat"), verbose = TRUE) { 
  if (verbose) {message("Creating phenodata_DM object...")} 
  phenodata_DM <- read.csv("data_DM/phenodataDMSOcontrols.csv", stringsAsFactors = FALSE)
  phenodata_DM$Dose_Unit <- "μM"
  
  if(verbose) {message("phenodata_DM object created!")}
  return(phenodata_DM)
}
phenodata_DM <- create_phenodata_DM("Rat")
rownames(phenodata_DM) <- phenodata_DM$samplename
phenodata_DM$duration <- floor(phenodata_DM$duration)
                                               
create_exprsdata_DM <- function(species=c("Rat"), phenodata_DM, verbose = TRUE) {
  if (verbose) {message("Creating eset object")}
  if(species == "Rat"){
  }
  #celFiles <- paste("/Users/sisira/Desktop/DrugMatrix_raw_files/drugmatrix.rathepatocyte.celfiles","/", phenodata_DM[,"celfilename"], sep="")
  
  #esetNorm <- just.rma(filenames = celFiles, verbose = TRUE, cdfname = "rat2302rnensgcdf")
  #saveRDS(esetNorm, "data/eset_DM.rds")
  eset <- readRDS("data/eset_DM.rds")
  
  storageMode(eset)<-"environment"
  eset <-subset(eset, substr(rownames(eset@assayData$exprs), 0, 4) != "AFFX")
  storageMode(eset)<-"lockedEnvironment"
  annotation(eset)<-"rna"
  if (verbose) {message("eset object created!")}
  return(eset)
}
eset <- create_exprsdata_DM("Rat", phenodata_DM)
colnames(eset) <- gsub(".CEL", "", colnames(eset))

create_featuredata_DM <- function(species=c("Rat"), eset, verbose = TRUE){
  if (verbose) {message("Creating featuredata_DM object...")}
  if (species == "Rat"){
    ensembl <- useMart("ensembl")
    datasets <- listDatasets(ensembl)
    ensembl = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host="uswest.ensembl.org",ensemblRedirect = FALSE)
    storageMode(eset) <- "environment"
    affxrows <- rownames(eset@assayData$exprs)
    rownames(eset@assayData$exprs) <- substr(rownames(eset@assayData$exprs), 1, nchar(affxrows)-3)
    #saveRDS(affxrows, file = "./data_DM/probesRatVitro1.rds") 
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
    
    if(verbose) {message("featuredata_DM object created!")}
    return(finalFeature)
    
  }
}
featuredata_DM <- create_featuredata_DM("Rat", eset)


create_Expressionset <- function(species=c("Rat"), eset, verbose = TRUE){
  if (verbose) {message("Creating expressionset...")}
  if (species == "Rat"){

pData(eset) <- phenodata_DM
fData(eset) <- featuredata_DM
storageMode(eset) <- "lockedEnvironment"
return(eset)
#saveRDS(eset, file="data_DM/ExpressionSet.rds")
  }
}

ExpressionSet <- create_Expressionset("Rat", eset)

all(rownames(phenodata_DM) == colnames(exprs(ExpressionSet)))

######################toxicoset constructor function ######################
create_curationDrug <- function(phenodata_DMData, verbose = TRUE){
  curationDrug <- unique(subset(phenodata_DM, select=c(drugid, dataset_drugid)))
  rownames(curationDrug) <- curationDrug$drugid
  names(curationDrug) <- c("unique.drugid", "dataset_drugid")
  
  return(curationDrug)
}

curationDrug <- create_curationDrug(phenodata_DMData)

create_curationCell <- function(phenodata_DM, verbose = TRUE){
  curationCell <- unique(subset(phenodata_DM, select=c(cellid)))
  curationCell$dataset_cellid <- curationCell$cellid
  names(curationCell) <- c("unique.cellid", "dataset_cellid")
  rownames(curationCell) <- curationCell$unique.cellid
  
  return(curationCell)
}
curationCell <- create_curationCell(phenodata_DM)


create_curationTissue <- function(phenodata_DM, verbose = TRUE){
  curationTissue <- unique(subset(phenodata_DM, select=c(organ_id)))
  curationTissue$dataset_tissueid <- "Liver"
  names(curationTissue)[1] <- "unique.tissueid"
  rownames(curationTissue) <- "Hepatocyte"
  
  return(curationTissue)
}
curationTissue <- create_curationTissue(phenodata_DM)


drug <- readRDS("data_DM/drug.rds")
rownames(drug)[1] <- "1-Naphthyl isothiocyanate" 
rownames(drug)[122] <- "Valproic acid" 
rownames(drug)[75] <- "Mercaptopurine" 
drug$drugid <- rownames(drug)


create_cell <- function(phenodata_DM, verbose = TRUE){
  cell <- unique(subset(phenodata_DM,select=c(cellid, organ_id, species, test_type)))
  names(cell)<-c("cellid","tissueid", "species","testType")
  cell$tissueid<-"Liver"
  rownames(cell) <- cell$cellid
  
  return(cell)
}
cell <- create_cell(phenodata_DM)

drugMatrix_new <- ToxicoSet("drugMatrix",
                  molecularProfiles=list("rna"=ExpressionSet),
                  cell=cell,
                  drug=drug,
                  sensitivityInfo=NULL,
                  sensitivityRaw=NULL,
                  sensitivityProfiles=NULL,
                  curationDrug=curationDrug,
                  curationCell=curationCell,
                  curationTissue=curationTissue,
                  datasetType = c("perturbation"),
                  verify = TRUE)

saveRDS(drugMatrix_new, "Updated tsets/drugMatrix.rds")
