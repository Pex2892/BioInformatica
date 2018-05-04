import os
import rpy2.robjects as robjects
import json

def analysisTHYM(idxAnalysis, data_seq, path, pvalue, foldchange, contrasts, geneSelected):
    mess = ""

    path_with_file = "datasetTCGA" + os.sep + "THYM" + os.sep + data_seq + "Seq" + os.sep + path
    name_file_json = "static/tempAnalysis" + os.sep + "THYM_" + data_seq + pvalue + foldchange + contrasts + ".json"
    name_file_heatmap = "static/tempAnalysis" + os.sep + "THYM_" + data_seq + pvalue + foldchange + contrasts + "_heatmap.jpg"
    name_file_rdata = "static/tempAnalysis" + os.sep + "THYM_workspace.RData"
    name_file_boxplot = "static/tempAnalysis" + os.sep + "THYM_" + geneSelected + "_boxplot.jpg"

    name_file_mithril_1 = "static/tempAnalysis" + os.sep + "input_mithril_THYM_" + data_seq + "_" + contrasts
    name_file_mithril_2 = "static/tempAnalysis" + os.sep + "input_mithril_THYM_" + data_seq + "_" + contrasts +".txt" #questo è per i contrasti che non sono ALL

    robjects.r('''
        library("limma")
        library("Biobase")
        library("gplots")
        library("jsonlite")
        
        if({7} == 4) {{
            load(file="{6}") 
            
            if({5} == 0) {{
                for (i in 1:5) {{
                    name.file <- paste(paste("{10}", i, sep="_"), ".txt", sep="")
                    write.table(results[i], name.file, sep="\t" ,row.names=TRUE, col.names = FALSE)
                }} 
            }}
                
            if({5} > 0) 
                write.table(results[1], "{11}", sep="\t" ,row.names=TRUE, col.names = FALSE)                            
        }} else if({7} == 3) {{
            load(file="{6}") 
            gene <- "{8}"
            exp.values <- as.vector(normalized.expressions$E[gene, rownames(df.patient.1)])
            boxplot.data <- data.frame(Sample.Category=df.patient.1$masaoka_stage, Sample.Value=exp.values, row.names=rownames(df.patient.1))
            jpeg(file="{9}")
            boxplot(Sample.Value~Sample.Category, data=boxplot.data, main=paste0("Evaluation of ", gene), xlab="Sample Category", ylab="Expression Value")
            dev.off()
        }} else {{
            if({7} == 1) {{
                # =========== LOAD BIOSPECIMEN CLINICAL ===========
                col.df.1 <- c("bcr_patient_uuid","bcr_patient_barcode","tumor_status","masaoka_stage","histologic_diagnosis","gender","vital_status")
                df.biospecimen.1 <- read.table("datasetTCGA/THYM/BiospecimenClinicalData/nationwidechildrens.org_clinical_patient_thym.txt", header=T, sep="\t")[col.df.1]
                df.biospecimen.1 <- df.biospecimen.1[-c(1,2),]
        
                col.df.2 <- c("bcr_patient_uuid","bcr_sample_barcode","bcr_aliquot_barcode","bcr_aliquot_uuid","biospecimen_barcode_bottom")
                df.biospecimen.2 <- read.table("datasetTCGA/THYM/BiospecimenClinicalData/nationwidechildrens.org_biospecimen_aliquot_thym.txt", header=T, sep="\t")[col.df.2]
                df.biospecimen.2 <- df.biospecimen.2[-1,]
        
                # =========== LOAD Seq ===========
                df.load <- read.table("{0}", header = T, sep = "\t", check.names=FALSE)[-1,]
                
                # =========== INTEGRATION DATA ===========
                df.patient <- merge(x = df.biospecimen.1, y = df.biospecimen.2[,1:3], by = "bcr_patient_uuid", all.x = TRUE)
                rownames(df.patient) <- as.vector(df.patient[,9]) # 9 si riferisce alla colonna aliquot_barcode
        
                df.patient <- df.patient[-which(df.patient$masaoka_stage == "[Not Available]"),] #Rimuovo tutti i levels '[Not Available]'
                df.patient$masaoka_stage <- droplevels(df.patient)$masaoka_stage #rimuovo i livelli non utilizzati
        
                design.original <- model.matrix(~0+masaoka_stage, data=df.patient)
                colnames(design.original) <- c("S1", "S2a", "S2b", "S3", "S4a","S4b")
        
                rownames(df.load) <- as.vector(df.load[,1]) #rinomino l'indice della riga col l'id del gene
                df.load <- df.load[,-1] #rimuovo la colonna relativa l'id del gene
                indx <- sapply(df.load, is.factor)
                df.load[indx] <- lapply(df.load[indx], function(x) as.numeric(as.character(x))) #converto il df da factor a numeric
        
                #bisogna cercare le colonne di df_miRNA in design_original e salvarle
                design.load <- design.original[row.names(design.original) %in% colnames(df.load), ] #fa il matching tra le righe di design.original e df.miRNA.1
                #faccio la trasposta di df.miRNA.1
                df.load <- as.data.frame(t(df.load))
                #in df.miRNA.1 mantengo solo le righe che hanno almeno un matching
                df.load <- as.data.frame(t(df.load[row.names(df.load) %in% row.names(design.load), ]))
                #prelevo solo i pazienti che hanno il matching
                df.patient.1 <- df.patient[row.names(df.patient) %in% colnames(df.load), ]
                        
                # =========== ANALYSIS ===========
                normalized.expressions <- voom(df.load, design.load)
                initial.model <- lmFit(normalized.expressions, design.load)
                limma.model <- eBayes(initial.model)
        
                #ogni volta usate lo stadio precedente come controllo per vedere quali geni differiscono
                contrasts <- makeContrasts(S2a-S1, S2b-S2a, S3-S2b, S4a-S3, S4b-S4a, levels=design.load)
                contrasts.model <- eBayes(contrasts.fit(limma.model, contrasts))
                
                save.image(file="{6}")
            }} else if({7} == 2) {{
                load(file="{6}") 
            }}
            
            if({5} == 0)
                results <- topTableF(contrasts.model, number=nrow(df.load), adjust.method="BH", p.value={1}, lfc={2})
            
            if({5} > 0)
                results <- topTable(contrasts.model, coef={5}, number=nrow(df.load), adjust.method="BH", p.value={1}, lfc={2})
                
            print(nrow(results))
    
            json <- toJSON(results, pretty = T)
            #cat(x) #serve a stamparlo
            write(json, "{3}")
            
            de.genes <- rownames(results)
            de.expressions <- normalized.expressions$E[de.genes,]
            jpeg(file="{4}")
            heatmap.2(de.expressions, col=redgreen(100), scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none")
            dev.off()
        }}
        
        '''.format(path_with_file, float(pvalue), float(foldchange), name_file_json, name_file_heatmap, contrasts, name_file_rdata, idxAnalysis, geneSelected, name_file_boxplot, name_file_mithril_1, name_file_mithril_2))

    if idxAnalysis == "1" or idxAnalysis == "2":
        #apro il json appena creato
        with open(name_file_json, 'r') as f:
            datastore = json.load(f)

        htmlGenes = '<option value="">Seleziona...</option>'
        htmlAllGenes = ''
        for idx,value in enumerate(datastore):
            htmlGenes += '<option value="'+value['_row']+'">'+value['_row']+'</option>'
            htmlAllGenes += value['_row']+"<br />"
        # questo html, mi serve per inserirlo nella lista dei geni

        return htmlGenes, htmlAllGenes

    return '', ''