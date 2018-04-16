import os
import rpy2.robjects as robjects
import rpy2.interactive as r

def downloadData(typeTumor, biocli, rnaseq, mirnaseq):
    mess = ""
    rootDir = os.getcwd() + os.sep

    os.chdir("datasetTCGA")

    if not os.path.isdir(typeTumor):
        os.makedirs(typeTumor)  # creo cartella con typeTumor

    os.chdir(typeTumor) # accedo alla directory appena creata

    robjects.r.source(rootDir + "R_library/Module_A.R") #importo libreria
    robjects.r.source(rootDir + "R_library/Module_B.R") #importo libreria

    #installare tutte le dipendenze richieste --> httr, rjson, RCurl, HGNChelper

    if biocli == 'true':
        if not os.path.isdir("BiospecimenClinical"):
            os.makedirs("BiospecimenClinical")  # creo cartella per la prima opzione

        os.chdir("BiospecimenClinical")  # accedo alla directory appena creata

        # scarico i dati inserendo codice R
        robjects.r('DownloadBiospecimenClinicalData(cancerType = "{0}")'.format(typeTumor))

        os.chdir("..")


    if rnaseq == 'true':
        if not os.path.isdir("RNASeq"):
            os.makedirs("RNASeq")  # creo cartella per la seconda opzione

        os.chdir("RNASeq")  # accedo alla directory appena creata

        # scarico i dati inserendo codice R
        robjects.r('DownloadRNASeqData(cancerType = "{0}", assayPlatform = NULL, tissueType = NULL, inputPatientIDs = NULL)'.format(typeTumor))

        os.chdir("..")


    if mirnaseq == 'true':
        if not os.path.isdir("miRNASeq"):
            os.makedirs("miRNASeq")  # creo cartella per la seconda opzione

        os.chdir("miRNASeq")  # accedo alla directory appena creata

        # scarico i dati inserendo codice R
        robjects.r('DownloadmiRNASeqData(cancerType = "{0}", assayPlatform = NULL, tissueType = NULL, inputPatientIDs = NULL)'.format(typeTumor))

        os.chdir("..")

    os.chdir("..")
    os.chdir("..")