import os
import rpy2.robjects as robjects
import json

def analysisLIHC(idxAnalysis, file_miRNA, file_RNA, pvalue, foldchange, contrasts, geneSelected):
    root = "datasetTCGA" + os.sep + "LIHC" + os.sep

    path_file_miRNA = root + "miRNASeq" + os.sep + file_miRNA
    path_file_RNA = root + "RNASeq" + os.sep + file_RNA

    name_file_rdata = "static/tempAnalysis" + os.sep + "LIHC_workspace.RData"
    name_file_json = "static/tempAnalysis" + os.sep + "LIHC_" + pvalue + "_" + foldchange + "_" + contrasts + ".json"
    name_file_heatmap = "static/tempAnalysis" + os.sep + "LIHC_" + pvalue + "_" + foldchange + "_" + contrasts + "_heatmap.jpg"
    name_file_boxplot = "static/tempAnalysis" + os.sep + "LIHC_" + pvalue + "_" + foldchange + "_" + contrasts + "_" + geneSelected + "_boxplot.jpg"
    name_file_mithril_1 = "static/tempAnalysis" + os.sep + "input_mithril_LIHC_" + contrasts
    name_file_mithril_2 = "static/tempAnalysis" + os.sep + "input_mithril_LIHC_" + contrasts + ".txt" #questo è per i contrasti che non sono ALL