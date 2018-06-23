#source("Module_A.R")
#source("Module_B.R")
library("limma")
library("Biobase")
library("gplots")
#typeCancer <- "THYM"
#savePath <- "../TCGA-THYM/mirna_seq"

setwd('/Users/giuseppesgroi/Desktop/BioInformatica')

# =========== DOWNLOAD DATASET ===========
# DownloadRNASeqData(typeCancer, assayPlatform = NULL, tissueType = NULL, saveFolderName = savePath, inputPatientIDs = NULL)
# DownloadCPTACData(typeCancer, assayPlatform = NULL, tissueType = NULL, saveFolderName = savePath, inputPatientIDs = NULL)
# DownloadBiospecimenClinicalData(typeCancer, saveFolderName = savePath)
# mi_RNA<-DownloadmiRNASeqData(typeCancer, assayPlatform = NULL, tissueType = NULL, saveFolderName = "miRNA_seq", inputPatientIDs = NULL)

# =========== LOAD BIOSPECIMEN CLINICAL ===========
col.df.1 <- c("bcr_patient_uuid","bcr_patient_barcode","tumor_status","masaoka_stage","histologic_diagnosis","gender","vital_status")
#col.df.2 <- c("bcr_sample_barcode","bcr_sample_uuid","sample_type","vial_number","sample_type_id")
#c_df2=c("bcr_sample_barcode","bcr_analyte_barcode","bcr_analyte_uuid","analyte_type","subportion_sequence")
col.df.3 <- c("bcr_patient_uuid","bcr_sample_barcode","bcr_aliquot_barcode","bcr_aliquot_uuid","biospecimen_barcode_bottom")

df.1 <- read.table("THYM/BiospecimenClinical/nationwidechildrens.org_clinical_patient_thym.txt", header=T, sep="\t")[col.df.1]
df.1 <- df.1[-c(1,2),]
# df1 <- read.table("THYM/BiospecimenClinical/nationwidechildrens.org_biospecimen_sample_thym.txt", header=T, sep="\t")[c_df1]
# df1<-df1[-1,]
#df2 <- read.table("THYM/BiospecimenClinical/nationwidechildrens.org_biospecimen_analyte_thym.txt", header=T, sep="\t")[c_df2]
# df2<-df2[-1,]
df.3 <- read.table("THYM/BiospecimenClinical/nationwidechildrens.org_biospecimen_aliquot_thym.txt", header=T, sep="\t")[col.df.3]
df.3 <- df.3[-1,]

# =========== LOAD RNA Seq ===========
#df.RNA.1 <- read.table("THYM/RNASeq/THYM__gene_RNAseq__tissueTypeAll__20180421104414.txt", header = T, sep = "\t", check.names=FALSE)
#df.RNA.2 <- read.table("THYM/RNASeq/THYM__gene.normalized_RNAseq__tissueTypeAll__20180421104232.txt", header = T, sep = "\t", check.names=FALSE)

# df_RNA <- read.table("RNASeq/THYM__exon_RNAseq__tissueTypeAll__20180421110230.txt", header = T, sep = "\t")[1:3,]
# df_RNA1 <- read.table("RNASeq/THYM__exonJunction_RNAseq__tissueTypeAll__20180421111800.txt", header = T, sep = "\t")[1:3,]
# df_RNA4 <- read.table("RNASeq/THYM__isoform_RNAseq__tissueTypeAll__20180421105100.txt", header = T, sep = "\t")[1:3,]
# df_RNA5 <- read.table("RNASeq/THYM__isoform.normalized_RNAseq__tissueTypeAll__20180421104702.txt", header = T, sep = "\t")[1:3,]

# =========== LOAD miRNA Seq ===========
df.miRNA.1 <- read.table("THYM/miRNASeq/THYM__mir_HiSeq.hg18__tissueTypeAll__20180427152414.txt", header = T, sep = "\t", check.names=FALSE)
df.miRNA.1 <- df.miRNA.1[-1,]

df.miRNA.2 <- read.table("THYM/miRNASeq/THYM__mir_HiSeq.hg19.mirbase20__tissueTypeAll__20180427152542.txt", header = T, sep = "\t", check.names=FALSE)
df.miRNA.2 <- df.miRNA.2[-1,]


# =========== INTEGRATION DATA ===========
df.patient <- merge(x = df.1, y = df.3[,1:3], by = "bcr_patient_uuid", all.x = TRUE)
rownames(df.patient) <- as.vector(df.patient[,9]) # 9 si riferisce alla colonna aliquot_barcode

df.patient <- df.patient[-which(df.patient$masaoka_stage == "[Not Available]"),] #Rimuovo tutti i levels '[Not Available]'
df.patient$masaoka_stage <- droplevels(df.patient)$masaoka_stage #rimuovo i livelli non utilizzati

design.original <- model.matrix(~0+masaoka_stage, data=df.patient)
colnames(design.original) <- c("S1", "S2a", "S2b", "S3", "S4a","S4b") 

rownames(df.miRNA.1) <- as.vector(df.miRNA.1[,1])
df.miRNA.1 <- df.miRNA.1[,-1]
indx <- sapply(df.miRNA.1, is.factor)
df.miRNA.1[indx] <- lapply(df.miRNA.1[indx], function(x) as.numeric(as.character(x)))


#bisogna cercare le colonne di df_miRNA in design_original e salvarle
design.miRNA.1 <- design.original[row.names(design.original) %in% colnames(df.miRNA.1), ] #fa il matching tra le righe di design.original e df.miRNA.1
#faccio la trasposta di df.miRNA.1
df.miRNA.1 <- as.data.frame(t(df.miRNA.1))
#in df.miRNA.1 mantengo solo le righe che hanno almeno un matching
df.miRNA.1 <- as.data.frame(t(df.miRNA.1[row.names(df.miRNA.1) %in% row.names(design.miRNA.1), ]))
#prelevo solo i pazienti che hanno il matching
df.patient.1 <- df.patient[row.names(df.patient) %in% colnames(df.miRNA.1), ]

# =========== ANALYSIS ===========
normalized.expressions <- voom(df.miRNA.1, design.miRNA.1)
initial.model <- lmFit(normalized.expressions, design.miRNA.1)
limma.model <- eBayes(initial.model)
#ogni volta usate lo stadio precedente come controllo per vedere quali geni differiscono
contrasts <- makeContrasts(S2a-S1, S2b-S2a, S3-S2b, S4a-S3, S4b-S4a, levels=design.miRNA.1)
contrasts.model <- eBayes(contrasts.fit(limma.model, contrasts))
results.all <- topTableF(contrasts.model, number=nrow(df.miRNA.1), adjust.method="BH", p.value=0.2, lfc=1.5)
#results.S1.vs.N <- topTable(contrasts.model, coef=1, number=nrow(df_miRNA), adjust.method="BH", p.value=0.9, lfc=1.5)

de.genes <- rownames(results.all)
de.expressions <- normalized.expressions$E[de.genes,]
heatmap.2(de.expressions, col=redgreen(100), scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none")

# gene <- "hsa-mir-7-3"
# exp.values <- as.vector(normalized.expressions$E[gene, rownames(df.patient.1)])
# boxplot.data <- data.frame(Sample.Category=df.patient.1$masaoka_stage, Sample.Value=exp.values, row.names=rownames(df.patient.1))
# boxplot(Sample.Value~Sample.Category, data=boxplot.data, main=paste0("Evaluation of ", gene), xlab="Sample Category", ylab="Expression Value")
