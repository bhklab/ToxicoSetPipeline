drugcuration <- read.csv(file = "/Users/parwaiznijrabi/ToxicoGx/TGX/drugmappingnoCID.csv")
drugcuration

saveRDS(drugcuration, file ="/Users/parwaiznijrabi/Desktop/curationDrug3.rds")
readRDS("/Users/parwaiznijrabi/Desktop/curationDrug3.rds")


curationDrug <- readRDS("/Users/parwaiznijrabi/Desktop/curationDrug3.rds")
curationDrug

#125 compounds for hepatocyte 
#added loratadine, but may not be meant to be there. possiblt may be 124 compouds only 
curationdrugUpdate <- read.csv(file = "/Users/parwaiznijrabi/Desktop/curationDrug5.csv")
curationdrugUpdate
setdiff(rownames(drug), rownames(curationdrugUpdate))
curationdrugUpdate
drug #from drug1.R 
curationdrugUpdate
rownames(curationdrugUpdate) <- rownames(drug)
curationdrugUpdate

saveRDS(curationdrugUpdate, file ="/Users/parwaiznijrabi/Desktop/curationDrug8.rds")
