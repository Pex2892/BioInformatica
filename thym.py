import os
import rpy2.robjects as robjects
import json

def analysisTHYM(idxAnalysis, file_miRNA, file_RNA, pvalue, foldchange, contrasts, geneSelected):
    root = "datasetTCGA" + os.sep + "THYM" + os.sep

    path_file_miRNA = root + "miRNASeq" + os.sep + file_miRNA
    path_file_RNA = root + "RNASeq" + os.sep + file_RNA

    name_file_rdata = "static/tempAnalysis" + os.sep + "THYM_workspace.RData"
    name_file_json = "static/tempAnalysis" + os.sep + "THYM_" + pvalue + "_" + foldchange + "_" + contrasts + ".json"
    name_file_heatmap = "static/tempAnalysis" + os.sep + "THYM_" + pvalue + "_" + foldchange + "_" + contrasts + "_heatmap.jpg"
    name_file_boxplot = "static/tempAnalysis" + os.sep + "THYM_" + pvalue + "_" + foldchange + "_" + contrasts + "_" + geneSelected + "_boxplot.jpg"
    name_file_mithril_1 = "static/tempAnalysis" + os.sep + "input_mithril_THYM_" + contrasts
    name_file_mithril_2 = "static/tempAnalysis" + os.sep + "input_mithril_THYM_" + contrasts + ".txt" #questo è per i contrasti che non sono ALL

    robjects.r('''
        library("limma")
        library("Biobase")
        library("gplots")
        library("data.table")
        library("jsonlite")
        
        if({0} == 4) {{ #Creazione File di input per mithril
            load(file="{3}") 

            if({4} == 0) {{
                for (i in 1:5) {{
                    name.file <- paste(paste("{11}", i, sep="_"), ".txt", sep="")
                    write.table(results[i], name.file, sep="\t" ,row.names=TRUE, col.names = FALSE, quote = FALSE)
                }} 
            }} else {{
                write.table(results[1], "{12}", sep="\t" ,row.names=TRUE, col.names = FALSE, quote = FALSE)                            
            }}
            
        }} else if({0} == 3) {{ #Creazione boxplot per singolo gene
            load(file="{3}")
             
            gene <- "{9}"
            exp.values <- as.vector(normalized.expressions$E[gene, rownames(df.patient)])
            boxplot.data <- data.frame(Sample.Category=df.patient$masaoka_stage, Sample.Value=exp.values, row.names=rownames(df.patient))
            
            jpeg(file="{10}")
            boxplot(Sample.Value~Sample.Category, data=boxplot.data, main=paste0("Evaluation of ", gene), xlab="Sample Category", ylab="Expression Value")
            dev.off()
            
        }} else {{
            if({0} == 1) {{ #Analisi per l'estrazione dei biomarcatori
                # =========== LOAD BIOSPECIMEN CLINICAL ===========
                col.df.1 <- c("bcr_patient_uuid","bcr_patient_barcode","tumor_status","masaoka_stage","histologic_diagnosis","gender","vital_status")
                df.patient <- read.table("datasetTCGA/THYM/BiospecimenClinicalData/nationwidechildrens.org_clinical_patient_thym.txt", header=T, sep="\t")[col.df.1]
                df.patient <- df.patient[-c(1,2),]
                
                # =========== LOAD miRNA Seq ===========
                df.miRNA <- read.table("{1}", header = T, sep = "\t", check.names=FALSE, row.names = 1)
                df.miRNA <- data.matrix(df.miRNA[-1, seq(from=2, to=ncol(df.miRNA), by=2)])
                
                # =========== LOAD RNA Seq ===========
                df.RNA <- read.table("{2}", header = T, sep = "\t", check.names=FALSE, row.names = 1)
                df.RNA <- data.matrix(df.RNA)
                
                # =========== INTEGRATION DATA ===========
                rownames(df.patient) <- as.vector(df.patient$bcr_patient_barcode)

                df.patient <- df.patient[-which(df.patient$masaoka_stage == "[Not Available]"),] #Rimuovo tutti i levels '[Not Available]'
                df.patient$masaoka_stage <- droplevels(df.patient)$masaoka_stage #rimuovo i livelli non utilizzati
                
                col.patient.barcode <- as.vector(unique(df.patient$bcr_patient_barcode))
                
                idx.col.miRNA <- array(unlist(lapply(col.patient.barcode, function(x) which(colnames(df.miRNA) %like% as.character(x)))))
                
                idx.col.RNA <- array(unlist(lapply(col.patient.barcode, function(x) which(colnames(df.RNA) %like% as.character(x)))))
                
                df.miRNA <- as.data.frame(df.miRNA[, idx.col.miRNA])
                names(df.miRNA) <- substr(names(df.miRNA), -1, 12)
                
                df.RNA <- as.data.frame(df.RNA[, idx.col.RNA])
                names(df.RNA) <- substr(names(df.RNA), -1, 12)
                
                col.union <- intersect(colnames(df.RNA), colnames(df.miRNA))
                df.union <- data.matrix(rbind(df.miRNA[,col.union], df.RNA[,col.union]))
                
                df.patient <- df.patient[match(colnames(df.union), df.patient$bcr_patient_barcode), ]
                
                design.original <- model.matrix(~0+masaoka_stage, data=df.patient)
                colnames(design.original) <- c("S1", "S2a", "S2b", "S3", "S4a","S4b") 
                
                # =========== ANALYSIS ===========
                normalized.expressions <- voom(df.union, design.original)
                initial.model <- lmFit(normalized.expressions, design.original)
                limma.model <- eBayes(initial.model)
                
                #ogni volta usate lo stadio precedente come controllo per vedere quali geni differiscono
                contrasts <- makeContrasts(S2a-S1, S2b-S2a, S3-S2b, S4a-S3, S4b-S4a, levels=design.original)
                contrasts.model <- eBayes(contrasts.fit(limma.model, contrasts))
                
                save.image(file="{3}")
            }} else if({0} == 2) {{ #Analisi per l'estrazione dei biomarcatori (Custom)
                load(file="{3}") 
            }}
            
            if({4} == 0) {{
                results <- topTableF(contrasts.model, number=nrow(df.union), adjust.method="BH", p.value={5}, lfc={6})
            }} else {{
                results <- topTable(contrasts.model, coef={4}, number=nrow(df.union), adjust.method="BH", p.value={5}, lfc={6})
            }}
             
            #print(nrow(results))
            
            if(is.data.frame(results) && nrow(results)!=0) {{   
                json <- toJSON(results, pretty = T)
                #cat(x) #serve a stamparlo
                write(json, "{7}")
                
                de.genes <- rownames(results)
                de.expressions <- normalized.expressions$E[de.genes,]
                jpeg(file="{8}")
                heatmap.2(de.expressions, col=redgreen(100), scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none")
                dev.off()
            }}
        }}
        '''.format(idxAnalysis, path_file_miRNA, path_file_RNA, name_file_rdata, contrasts, pvalue,
                   foldchange, name_file_json, name_file_heatmap, geneSelected, name_file_boxplot, name_file_mithril_1, name_file_mithril_2))

    nGenes = 0
    if os.path.isfile(name_file_json) and (idxAnalysis == "1" or idxAnalysis == "2"):

        #apro il json appena creato
        with open(name_file_json, 'r') as f:
            datastore = json.load(f)

        htmlGenes = '<option value="">Seleziona...</option>'
        htmlAllGenes = ''

        for idx,value in enumerate(datastore):
            htmlGenes += '<option value="'+value['_row']+'">'+value['_row']+'</option>'
            htmlAllGenes += value['_row']+"<br />"
            nGenes += 1
        # questo html, mi serve per inserirlo nella lista dei geni

        return htmlGenes, htmlAllGenes, nGenes, 0
    elif idxAnalysis == "3" or idxAnalysis == "4": #è importante questa condizione, perchè allora causerà un errore sulla creazione del boxplot o nell'estrazione dei file di input
        return '', '', nGenes, 0

    return '', '', nGenes, 1
