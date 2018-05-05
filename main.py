from flask import Flask, render_template, request, jsonify
from dataset import *
from thym import *
import subprocess as sp

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
    data_seq = request.form['data_seq']
    path = request.form['path']
    pvalue = request.form['pvalue']
    foldchange = request.form['foldchange']
    contrasts = request.form['contrasts']
    idxAnalysis = request.form['idxAnalysis']
    geneSelected = request.form['geneSelected']

    if nameTumor and data_seq and path and pvalue and foldchange and contrasts and idxAnalysis:
        if nameTumor == 'THYM':
            htmlListGenes, htmlAllGenes = analysisTHYM(idxAnalysis, data_seq, path, pvalue, foldchange, contrasts,
                                                       geneSelected)

        return jsonify({'status': 1, 'type': 'success', 'listGene': htmlListGenes, 'AllListGene': htmlAllGenes})

    return jsonify({'status': 0, 'type': 'error', 'mess': 'Missing data!'})


# ***** ANALYSIS MITHRIL *****
@app.route('/analysisMithril', methods=['POST'])
def process3():
    fileName = request.form['fileName']

    rootDir = os.getcwd() + os.sep

    print(fileName)

    # nome del file errato, perchè non c'è la differenza fra mirna e rna
    # aggiungere controllo se il file cercato esiste nella dir tempAnalysis
    if fileName:
        try:
            sp.call(
                ['java', '-jar', rootDir + 'library/R/MITHrIL2.jar', 'merged-mithril', '-verbose', '-i', fileName, '-o', 'out.txt', '-p',
                 'outPertubations.txt'], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
            return jsonify({'status': 1, 'type': 'success', 'mess': 'Il file output di <b>MITHrIL</b>, sono stati con successo!'})
        except:
            return jsonify({'status': 0, 'type': 'error', 'mess': 'ERRORE JAVA!'})


if __name__ == '__main__':
    app.run(port=2892, debug=True)
