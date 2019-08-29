library(PharmacoGx)
library(Biobase)


drug.info <- readRDS("/pfs/FIMMdata/drug.info.rds")
cell.info <- readRDS("/pfs/FIMMdata/cell.info.rds")

curationCell <- readRDS("/pfs/FIMMdata/curationCell.rds")
curationDrug <- readRDS("/pfs/FIMMdata/curationDrug.rds")
curationTissue <- readRDS("/pfs/FIMMdata/curationTissue.rds")

sens.info <- readRDS("/pfs/FIMMdata/sens.info.rds")
sens.prof <- readRDS("/pfs/FIMMdata/sens.prof.rds")

sens.raw <- readRDS("/pfs/FIMMdata/sens.raw.rds")

profiles <- get(load(file.path("pfs", "FIMMProfiles/profiles.RData")))
profiles <- profiles[rownames(sensitivity_info),]


sens.prof <- cbind(sens.prof[,"aac_published"], profiles)

message("aac correlations are")

message(cor(sens.prof[,"aac_published"], sens.prof[,"AAC"]))

dummyRNA <- Biobase::ExpressionSet()
pData(dummyRNA)$cellid <- character()
pData(dummyRNA)$batchid <- character()

fData(dummyRNA)$BEST <- vector()
fData(dummyRNA)$Symbol <- character()

annotation(dummyRNA) <- "rna"


message("Making PSet")


FIMM <-   PharmacoSet(molecularProfiles=list("rna"=dummyRNA),
                      name="FIMM", 
                      cell=cell.info, 
                      drug=drug.info, 
                      sensitivityInfo=sens.info, 
                      sensitivityRaw=sens.raw.array, 
                      sensitivityProfiles=sens.prof, 
                      sensitivityN=NULL, 
                      curationCell=curationCell, 
                      curationDrug=curationDrug, 
                      curationTissue=curationTissue, 
                      datasetType="sensitivity")


message("Saving")


save(FIMM, file="FIMM.RData", version=2)

