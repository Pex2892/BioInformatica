library("limma")
library("Biobase")
library("gplots")
library("data.table")

setwd('/Users/giuseppesgroi/Desktop/BioInformatica')

path.file.miRNA <- "THCA/miRNASeq/THCA__mir_HiSeq.hg19.mirbase20__tissueTypeAll__20180430100324.txt"
# THCA/miRNASeq/THCA__mir_HiSeq.hg18__tissueTypeAll__20180430100150.txt
# THCA/miRNASeq/THCA__mir_HiSeq.hg19.mirbase20__tissueTypeAll__20180430100324.txt

path.file.RNA <- "THCA/RNASeq/THCA__gene.normalized_RNAseq__tissueTypeAll__20180430141622.txt"
# THCA/RNASeq/THCA__gene_RNAseq__tissueTypeAll__20180430153443.txt
# THCA/RNASeq/THCA__gene.normalized_RNAseq__tissueTypeAll__20180430141622.txt

# =========== SETTINGS PARAMETERS ===========
p.value.max <- 0.05
logFoldChange.min <- 1.5

# =========== LOAD BIOSPECIMEN CLINICAL ===========
col.df.1 <- c("bcr_patient_uuid", "bcr_patient_barcode", "tumor_status", "ajcc_pathologic_tumor_stage", "histologic_diagnosis", "gender", "vital_status")
df.patient <- read.table("THCA/BiospecimenClinicalData/nationwidechildrens.org_clinical_patient_thca.txt", header = T, sep = "\t", quote="")[col.df.1]
df.patient <- df.patient[-c(1, 2),]

# =========== LOAD miRNA Seq ===========
df.miRNA <- read.table(path.file.miRNA, header = T, sep = "\t", check.names = FALSE, row.names = 1)
df.miRNA <- data.matrix(df.miRNA[-1, seq(from = 2, to = ncol(df.miRNA), by = 2)])

# =========== LOAD RNA Seq ===========
df.RNA <- read.table(path.file.RNA, header = T, sep = "\t", check.names = FALSE, row.names = 1)
df.RNA <- data.matrix(df.RNA)

# =========== INTEGRATION DATA ===========
rownames(df.patient) <- as.vector(df.patient$bcr_patient_barcode)

df.patient <- df.patient[-which(
  df.patient$ajcc_pathologic_tumor_stage == "[Not Available]" |
  df.patient$ajcc_pathologic_tumor_stage == ""
),] #Rimuovo tutti i levels '[Not Available]'

df.patient$ajcc_pathologic_tumor_stage <- droplevels(df.patient)$ajcc_pathologic_tumor_stage #rimuovo i livelli non utilizzati

col.patient.barcode <- as.vector(unique(df.patient$bcr_patient_barcode))

idx.col.miRNA <- array(unlist(lapply(col.patient.barcode, function(x) which(colnames(df.miRNA) %like% as.character(x)))))

idx.col.RNA <- array(unlist(lapply(col.patient.barcode, function(x) which(colnames(df.RNA) %like% as.character(x)))))

df.miRNA <- as.data.frame(df.miRNA[, idx.col.miRNA])
names(df.miRNA) <- substr(names(df.miRNA), -1, 12)

df.RNA <- as.data.frame(df.RNA[, idx.col.RNA])
names(df.RNA) <- substr(names(df.RNA), -1, 12)

col.union <- intersect(colnames(df.RNA), colnames(df.miRNA))
df.union <- data.matrix(rbind(df.miRNA[, col.union], df.RNA[, col.union]))

df.patient <- df.patient[match(colnames(df.union), df.patient$bcr_patient_barcode), ]

design.original <- model.matrix(~0+ajcc_pathologic_tumor_stage, data=df.patient)
colnames(design.original) <- c("S1", "S2", "S3", "S4", "S4a", "S4c")

# =========== ANALYSIS ===========
normalized.expressions <- voom(df.union, design.original)
initial.model <- lmFit(normalized.expressions, design.original)
limma.model <- eBayes(initial.model)

#ogni volta usate lo stadio precedente come controllo per vedere quali geni differiscono
contrasts <- makeContrasts(S2-S1, S3-S2, S4-S3, S4a-S4, S4c-S4a, levels=design.original)
contrasts.model <- eBayes(contrasts.fit(limma.model, contrasts))

results.all <- topTableF(contrasts.model, number=nrow(df.union), adjust.method="BH", p.value=p.value.max, lfc=logFoldChange.min)

#save.image(file="prova.RData")

results.S2.vs.S1 <- topTable(contrasts.model, coef=1, number=nrow(df.union), adjust.method="BH", p.value=p.value.max, lfc=logFoldChange.min)
results.S3.vs.S2 <- topTable(contrasts.model, coef=2, number=nrow(df.union), adjust.method="BH", p.value=p.value.max, lfc=logFoldChange.min)
results.S4.vs.S3 <- topTable(contrasts.model, coef=3, number=nrow(df.union), adjust.method="BH", p.value=p.value.max, lfc=logFoldChange.min)
results.S4a.vs.S4 <- topTable(contrasts.model, coef=4, number=nrow(df.union), adjust.method="BH", p.value=p.value.max, lfc=logFoldChange.min)
results.S4c.vs.S4a <- topTable(contrasts.model, coef=5, number=nrow(df.union), adjust.method="BH", p.value=p.value.max, lfc=logFoldChange.min)

# for (i in 1:5) {
#   name.file <- paste(paste("filename_0", i, sep="_"), ".txt", sep="")
#   write.table(results.all[i], name.file, sep="\t" ,row.names=TRUE, col.names = FALSE, quote = FALSE)
# }

#write.table(results.S4a.vs.S3[1], "filename.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)

# if(nrow(results.all) > 1) { #evito di stampare la heatmap, se il dataframe result è composto di una sola riga
#   de.genes <- rownames(results.all)
#   de.expressions <- normalized.expressions$E[de.genes,]
#   heatmap.2(de.expressions, col=redgreen(100), scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none")
# }

# gene <- "hsa-mir-4999"
# exp.values <- as.vector(normalized.expressions$E[gene, rownames(df.patient)])
# boxplot.data <- data.frame(Sample.Category=df.patient$ajcc_pathologic_tumor_stage, Sample.Value=exp.values, row.names=rownames(df.patient))
# boxplot(Sample.Value~Sample.Category, data=boxplot.data, main=paste0("Evaluation of ", gene), xlab="Sample Category", ylab="Expression Value")
