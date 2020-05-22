metadata <- read.csv("/Users/parwaiznijrabi/ToxicoGx/TGX/s_Hepatocyte.csv", header = TRUE, sep = "\t")
library(gsubfn)
library(textclean)
library(dplyr)
as.data.frame(metadata$Comment.chEMBL.ID.)
chemblID <- as.data.frame(unique(metadata$Comment.chEMBL.ID.))
chemblID

chemblID <- chemblID[-c(22), ]
chemblID1 <- as.data.frame(chemblID)

ncol(metadata)
metadata$Comment.chEMBL.ID.

which(colnames(metadata)== "Factor.Value.Compound.")
which(colnames(metadata)== "Comment.chEMBL.ID.")
which(colnames(metadata)== "Term.Source.REF.5")
which(colnames(metadata)== "Term.Accession.Number.5")
which(colnames(metadata)== "CHEBI")


metadata$CHEBI<- paste(metadata$Term.Source.REF.5, metadata$Term.Accession.Number.5)
metadata

drugcuration <- unique(metadata[,c(24,28,53)])
drugcuration

drugcuration[,order(colnames(drugcuration))]
drugcuration[c(2,1,3)]

drugcuration

drugcuration <- as.data.frame(drugcuration)

drugcuration$CHEBI <- mgsub(drugcuration$CHEBI, c(" "), c(""))

drugcuration

names(drugcuration) <- c("drug.id.dixa", "CHEMBL.ID", "CHEBI.ID")
drugcuration
nrow(drugcuration) #125
drugcuration <- drugcuration[-c(126), ]
nrow(curationDrug) #125- added loratadine 


drugcuration

as.data.frame(drugcuration$drug.id.dixa)
as.data.frame(curationDrug$drug.id.dixa)




#desloratadine repeats in drugcuration, on s_hepatocyte its shown as desloratadine but maps to loratadine as well?
#drugcuration has n=125, curationdrug (drug curation object) maps to 124. desloratadine is the dupliate
#after this, make cell object 
#then make molecular profiles object (check screenshot)
#compile all into toxicoset and check object
#start on quality control ()

nrow(drugcuration)

drugcuration
ncol(drugcuration)


drug <- read.csv(file = "/Users/parwaiznijrabi/Desktop/drug1csv.csv")


colnames(drug) <- c("drugmatrix.drugid", "CHEMBL.ID", "CHEBI.ID")
rownames(drug) <- c("1-Naphthyl Isothiocyanate", "Ethinyl estradiol", "Estradiol", 
                    "3-Methylcholanthrene", "4,4'-Methylenedianiline", "Acetaminophen", 
                    "Tretinoin", "Allyl alcohol", "Amantadine", "Amineptine", "Amiodarone", 
                    "Amitriptyline", "Antimycin A", "Atorvastatin", "Azathioprine", 
                    "Busulfan", "Buthionine sulfoxamine", "Cadmium Dichloride", "Carbon tetrachloride", 
                    "Carmustine", "Carvedilol", "Celecoxib", "N-Tetracosanoylphytosphingosine", 
                    "Cerivastatin", "Chlorambucil", "Chloroquine", "Chlorpromazine", 
                    "Cisplatin", "Citalopram", "Clofibrate", "Clomipramine", "Clotrimazole", 
                    "Clozapine", "Colchicine", "Coralgil", "Cyclophosphamide", "Cyclosporin A", 
                    "Cytochalasin B", "D-Galactosamine", "Danazol", "Daunorubicin", 
                    "Desloratadine", "Dexamethasone", "Didanosine", "Diethylstilbestrol", 
                    "Diphenhydramine", "Zonalon", "Doxofylline", "Doxorubicin", "Erlotinib", 
                    "Ethanol", "Etoposide", "Fenofibrate", "Fluoxetine", "Fluphenazine", 
                    "Gabapentin", "Gefitinib", "Gemfibrozil", "Gentamicin", "Griseofulvin", 
                    "Haloperidol", "Idarubicin", "Ifosfamide", "Imatinib", "Isotretinoin", 
                    "Itraconazole", "Kanamycin", "Xolegel", "Labetalol", "Lansoprazole", 
                    "Lithocholic Acid", "Lomustine", "Loratadine", "Marimastat", 
                    "Marcaptopurine", "Methapyrilene", "Microcystin-LR", "Tandutinib", 
                    "Monocrotaline", "N-Nitrosodiethylamine", "N,N-Dimethylformamide", 
                    "Nevirapine", "Nimesulide", "Norethindrone", "Nortriptyline", 
                    "Olanzapine", "Omeprazole", "Pantoprazole", "Paraquat dichloride", 
                    "Paroxetine", "Pemoline", "Perhexiline", "Phalloidin", "Prinomastat", 
                    "Progesterone", "Quetiapine", "Rabeprazole", "Ramipril", "Rifampicin", 
                    "Risperidone", "Rofecoxib", "Sertraline", "Sildenafil", "Simvastatin", 
                    "Sirolimus", "Sotalol", "Sparteine", "Sporidesmin", "Staurosporine", 
                    "Stavudine", "Streptozocin", "Sulconazole", "Sulindac", "Sulpiride", 
                    "Tacrolimus", "Tamoxifen", "Tenidap", "Tetracycline", "TGF-Beta-1, human recombinant", 
                    "Valdecoxib", "Valeric acid", "Valproic Acid", "Venlafaxine", 
                    "Vincaleukoblastine", "Zidovudine")


drug
UMCTRL <- data.frame("unmapped_control", "0", "0")
colnames(UMCTRL) <- c("drugmatrix.drugid", "CHEMBL.ID", "CHEBI.ID")
rownames(UMCTRL) <- "unmapped_control"
drug <- rbind(drug, UMCTRL)
drug

drug

drugslist <- curationdrugUpdate2$unique.drugid
dput(as.character(drugslist))

setdiff(rownames(drug), curationDrug$unique.drugid)

rownames(drug)
saveRDS(drug, file ="/Users/parwaiznijrabi/Desktop/drug4.rds")
drug1 <- readRDS("/Users/parwaiznijrabi/Desktop/tsetworking/drug4.rds")

setdiff(rownames(drug1), rownames(curationDrug))



curationDrug
