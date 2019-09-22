library(PharmacoGx)
library(Biobase)


drug.info <- readRDS("/pfs/FIMMraw/drug.info.rds")
cell.info <- readRDS("/pfs/FIMMraw/cell.info.rds")

curationCell <- readRDS("/pfs/FIMMraw/curationCell.rds")
curationDrug <- readRDS("/pfs/FIMMraw/curationDrug.rds")
curationTissue <- readRDS("/pfs/FIMMraw/curationTissue.rds")

sens.info <- readRDS("/pfs/FIMMraw/sens.info.rds")
sens.prof <- readRDS("/pfs/FIMMraw/sens.prof.rds")

sens.raw <- readRDS("/pfs/FIMMraw/sens.raw.rds")

profiles <- get(load(file.path("pfs", "FIMMProfiles/profiles.RData")))
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

