library(PharmacoGx)
library(Biobase)


drug.info <- readRDS("/pfs/processFIMM/drug.info.rds")
cell.info <- readRDS("/pfs/processFIMM/cell.info.rds")

curationCell <- readRDS("/pfs/processFIMM/curationCell.rds")
curationDrug <- readRDS("/pfs/processFIMM/curationDrug.rds")
curationTissue <- readRDS("/pfs/processFIMM/curationTissue.rds")

sens.info <- readRDS("/pfs/processFIMM/sens.info.rds")
sens.prof <- readRDS("/pfs/processFIMM/sens.prof.rds")

sens.raw <- readRDS("/pfs/processFIMM/sens.raw.rds")

profiles <- get(load(file.path("pfs", "SliceAssemble/profiles.RData")))
profiles <- profiles[rownames(sens.info),]
profiles <- apply(profiles, c(1,2), as.numeric)

sens.prof <- cbind(sens.prof, profiles)

sens.prof <- sens.prof[,-c(1, 2)]

message("aac correlations are")

message(cor(sens.prof[,"aac_published"], sens.prof[,"aac_recomputed"], use="pairwise.complete"))

dummyRNA <- Biobase::ExpressionSet()
pData(dummyRNA)$cellid <- character()
pData(dummyRNA)$batchid <- character()

fData(dummyRNA)$BEST <- vector()
fData(dummyRNA)$Symbol <- character()

annotation(dummyRNA) <- "rna"

cellsPresent <- sort(unique(sens.info$cellid))
cell.info <- cell.info[cellsPresent,]


drugsPresent <- sort(unique(sens.info$drugid))

drug.info <- drug.info[drugsPresent,]


drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
drug_all <- drug_all[which(!is.na(drug_all[ , "FIMM.drugid"])),]
drug_all <- drug_all[ , c("unique.drugid", "FIMM.drugid","smiles","inchikey","cid","FDA")]
rownames(drug_all) <- drug_all[ , "unique.drugid"]

drug_all <- drug_all[rownames(drug.info),]
drug.info[,c("smiles","inchikey","cid","FDA")] <- drug_all[,c("smiles","inchikey","cid","FDA")]

colnames(cell.info)[which(names(cell.info) == "cellid")] <- "unique.cellid"
colnames(drug.info)[which(names(drug.info) == "drugid")] <- "unique.drugid"

message("Making PSet")


FIMM <-   PharmacoSet(molecularProfiles=list("rna"=dummyRNA),
                      name="FIMM", 
                      cell=cell.info, 
                      drug=drug.info, 
                      sensitivityInfo=sens.info, 
                      sensitivityRaw=sens.raw, 
                      sensitivityProfiles=sens.prof, 
                      sensitivityN=NULL, 
                      curationCell=curationCell, 
                      curationDrug=curationDrug, 
                      curationTissue=curationTissue, 
                      datasetType="sensitivity")


message("Saving")


saveRDS(FIMM, file="/pfs/out/FIMM.rds", version=2)

