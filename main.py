from flask import Flask, render_template, request, jsonify
from dataset import *
from thym import *
from thca import *
from lihc import *
import json
import subprocess as sp
from subprocess import Popen, PIPE, STDOUT

# ***** SETTINGS *****
if not os.path.isdir('datasetTCGA'):
    os.makedirs('datasetTCGA')


def emptyCache():
    # svuoto la cartella tempAnalysis
    filelist = [f for f in os.listdir('static/tempAnalysis')]
    for f in filelist:
        os.remove(os.path.join('static/tempAnalysis', f))


if not os.path.isdir('static/tempAnalysis'):
    os.makedirs('static/tempAnalysis')
else:
    emptyCache()

app = Flask(__name__)


# ***** INDEX *****
@app.route('/')
def index():
    return render_template('index.html')


# ***** BIOMARKERS *****
@app.route('/biomarkers')
def biomarkers():
    #emptyCache()
    return render_template('biomarkers.html')


# ***** MITHrIL *****
@app.route('/MITHrIL')
def mithril():
    return render_template('MITHrIL.html')


# ***** DOWNLOADER *****
@app.route('/downloader', methods=['POST'])
def process():
    typeTumor = request.form['typeTumor']
    biocli = request.form['biocli']
    rnaseq = request.form['rnaseq']
    mirnaseq = request.form['mirnaseq']

    if typeTumor and biocli and rnaseq and mirnaseq:
        mess = downloadData(typeTumor, biocli, rnaseq, mirnaseq)
        return jsonify({'status': 1, 'type': 'success', 'mess': ''})

    return jsonify({'status': 0, 'type': 'error', 'mess': 'Missing data!'})


# ***** ANALYSIS BIOMARKERS *****
@app.route('/analysisBiomarkers', methods=['POST'])
def process2():
    nameTumor = request.form['tumor']
    idxAnalysis = request.form['idxAnalysis']
    file_miRNA = request.form['file_miRNA']
    file_RNA = request.form['file_RNA']
    pvalue = request.form['pvalue']
    foldchange = request.form['foldchange']
    contrasts = request.form['contrasts']
    geneSelected = request.form['geneSelected']

    if nameTumor and idxAnalysis and file_miRNA and file_RNA and pvalue and foldchange and contrasts:
        if nameTumor == 'THYM':
            htmlListGenes, htmlAllGenes, nGenes, emptyDF = analysisTHYM(idxAnalysis, file_miRNA, file_RNA, pvalue, foldchange, contrasts, geneSelected)
        elif nameTumor == 'THCA':
            htmlListGenes, htmlAllGenes, nGenes, emptyDF = analysisTHCA(idxAnalysis, file_miRNA, file_RNA, pvalue, foldchange, contrasts, geneSelected)
        elif nameTumor == 'LIHC':
            htmlListGenes, htmlAllGenes, nGenes, emptyDF = analysisLIHC(idxAnalysis, file_miRNA, file_RNA, pvalue, foldchange, contrasts, geneSelected)

        return jsonify({'status': 1, 'type': 'success', 'listGene': htmlListGenes, 'AllListGene': htmlAllGenes, 'nGenes': nGenes, 'emptyDF': emptyDF})

    return jsonify({'status': 0, 'type': 'error', 'mess': 'Missing data!'})


# ***** ANALYSIS MITHRIL *****
@app.route('/analysisMithril', methods=['POST'])
def process3():
    pathtempAnalysis = 'static/tempAnalysis/'

    inputFileName = request.form['fileName']
    outPertubationsName = pathtempAnalysis + request.form['outPertubationsName']
    outMainFile = pathtempAnalysis + request.form['outMainFile']

    rootDir = os.getcwd() + os.sep

    jsonHTML1 = []
    #jsonHTML2 = []

    if inputFileName:
        if os.path.exists(pathtempAnalysis + inputFileName):
            try:
                sp.call(['java', '-jar', rootDir + 'library/java/MITHrIL2.jar', 'merged-mithril', '-verbose', '-i',
                         pathtempAnalysis + inputFileName, '-o', outMainFile, '-p', outPertubationsName])

                with open(rootDir + outMainFile, "r") as f:
                    next(f)
                    for rows in f:
                        row = rows.split("\t")
                        jsonTMP = {}
                        jsonTMP['PathwayId'] = row[0]
                        jsonTMP['PathwayName'] = row[1]
                        jsonTMP['RawAccumulator'] = row[2]
                        jsonTMP['ImpactFactor'] = row[3]
                        jsonTMP['probability_pi'] = row[4]
                        jsonTMP['total_perturbation'] = row[5]
                        jsonTMP['corrected_accumulator'] = row[6]
                        jsonTMP['pValue'] = row[7]
                        jsonTMP['adjusted_pValue'] = row[8]
                        jsonHTML1.append(jsonTMP)

                    jsonHTML1 = json.dumps(jsonHTML1)

                    # with open(rootDir + outPertubationsName, "r") as f:
                    #     next(f)
                    #     for rows in f:
                    #         row = rows.split("\t")
                    #         jsonTMP = {}
                    #         jsonTMP['PathwayId'] = row[0]
                    #         jsonTMP['PathwayName'] = row[1]
                    #         jsonTMP['GeneId'] = row[2]
                    #         jsonTMP['GeneName'] = row[3]
                    #         jsonTMP['perturbation'] = row[4]
                    #         jsonTMP['accumulator'] = row[5]
                    #         jsonTMP['pValue'] = row[6]
                    #         jsonHTML2.append(jsonTMP)
                    #
                    #     jsonHTML2 = json.dumps(jsonHTML2)
                    #
                    #     print(jsonHTML2)

                return jsonify({'status': 1, 'type': 'success', 'tableHtml': jsonHTML1, 'mess': 'I file output di <b>MITHrIL</b>, sono stati generati con successo!'})
            except:
                return jsonify({'status': 0, 'type': 'error', 'mess': 'Si è verificato un problema durante l\'esecuzione del jar'})
        else:
            return jsonify({'status': 0, 'type': 'error', 'mess': 'Il file di input non è stato trovato, potrebbe non essere stato generato!'})
    else:
        return jsonify({'status': 0, 'type': 'error', 'mess': 'Si è verificato un problema durante l\'esecuzione!'})


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=2892, debug=True)
