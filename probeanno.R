library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library("HGNChelper")
library("mitch")

# 450K
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
myann <- data.frame(anno[,c("UCSC_RefGene_Name","UCSC_RefGene_Group","Islands_Name","Relation_to_Island")])
gp <- myann[,"UCSC_RefGene_Name",drop=FALSE]
gp2 <- strsplit(gp$UCSC_RefGene_Name,";")
names(gp2) <- rownames(gp)
gp2 <- lapply(gp2,unique)
gt1 <- stack(gp2)
colnames(gt1) <- c("gene","probe")
gt1$probe <- as.character(gt1$probe)
dim(gt1)

if ( file.exists("new.hgnc.table.rds")) {
  new.hgnc.table <- readRDS("new.hgnc.table.rds")  
} else {
  new.hgnc.table <- getCurrentHumanMap()
  saveRDS(new.hgnc.table,"new.hgnc.table.rds")
}

fix <- checkGeneSymbols(gt1$gene,map=new.hgnc.table)
fix2 <- fix[which(fix$x != fix$Suggested.Symbol),]
length(unique(fix2$x))
gt1$gene <- fix$Suggested.Symbol
head(gt1)
saveRDS(gt1,"450K.rds")

# EPIC
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
myann <- data.frame(anno[,c("UCSC_RefGene_Name","UCSC_RefGene_Group","Islands_Name","Relation_to_Island")])
gp <- myann[,"UCSC_RefGene_Name",drop=FALSE]
gp2 <- strsplit(gp$UCSC_RefGene_Name,";")
names(gp2) <- rownames(gp)
gp2 <- lapply(gp2,unique)
gt2 <- stack(gp2)
colnames(gt2) <- c("gene","probe")
gt2$probe <- as.character(gt2$probe)
fix <- checkGeneSymbols(gt2$gene,map=new.hgnc.table)
fix2 <- fix[which(fix$x != fix$Suggested.Symbol),]
length(unique(fix2$x))
gt2$gene <- fix$Suggested.Symbol
saveRDS(gt2,"EPIC.rds")
