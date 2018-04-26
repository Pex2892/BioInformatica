from flask import Flask, render_template, request, jsonify
from dataset import *

# ***** SETTINGS *****
if not os.path.isdir('datasetTCGA'):
    os.makedirs('datasetTCGA')

app = Flask(__name__)


# ***** INDEX *****
@app.route('/')
def index():
	return render_template('index.html')







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


if __name__ == '__main__':
    app.run(port=2892, debug=True)