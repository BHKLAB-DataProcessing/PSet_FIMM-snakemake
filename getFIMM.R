library(PharmacoGx)
library(Biobase)
library(SummarizedExperiment)

args = commandArgs(trailingOnly=TRUE)
ORCESTRA_ID <- args

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

emptyE <- ExpressionSet()
pData(emptyE)$cellid <- character()
pData(emptyE)$batchid <- character()
fData(emptyE)$BEST <- vector()
fData(emptyE)$Symbol <- character()
annotation(emptyE) <- "FIMM contains no molecular profiles of cell lines. This SE is empty placeholder."

emptySE <- SummarizedExperiment::SummarizedExperiment(
  ## TODO:: Do we want to pass an environment for better memory efficiency?
  assays=S4Vectors::SimpleList(as.list(Biobase::assayData(emptyE))
  ),
  # Switch rearrange columns so that IDs are first, probes second
  rowData=S4Vectors::DataFrame(Biobase::fData(emptyE),
                               rownames=rownames(Biobase::fData(emptyE)) 
  ),
  colData=S4Vectors::DataFrame(Biobase::pData(emptyE),
                               rownames=rownames(Biobase::pData(emptyE))
  ),
  metadata=list("experimentData" = emptyE@experimentData, 
                "annotation" = Biobase::annotation(emptyE), 
                "protocolData" = Biobase::protocolData(emptyE)
  )
)


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

curationCell <- curationCell[rownames(cell.info),]
curationDrug <- curationDrug[rownames(drug.info),]
curationTissue <- curationTissue[rownames(cell.info),]

message("Making PSet")


FIMM <-   PharmacoGx::PharmacoSet(molecularProfiles=list("rna"=emptySE),
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
dataset <- "FIMM"
#output ORCESTRA_ID and Pachyderm commit id
write.table(dataset, file="/pfs/out/dataset.txt", row.names = F ,quote = F, sep = "\t", col.names = F)
write.table(ORCESTRA_ID, file="/pfs/out/orcestra_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)				   
pach_commit_id <- Sys.getenv("PACH_OUTPUT_COMMIT_ID")
write.table(pach_commit_id, file="/pfs/out/commit_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)
