# Software engineering project
# Authors : DENET Lola, DESQUERRE Émilie, HUI Tongyuxuan, MEGUERDITCHIAN Caroline, NIU Wenli, PRATX Julie, VU Thao Uyen
# M2 Bioinformatics - Bordeaux university
# December 14, 2020

# flask_app.py
# Used to make the link between the front-end and back-end part of the program.
# Allows user-program interaction through a web interface linked to the Python program using Flask.


from flask import Flask, render_template, request
from datetime import timedelta
import main

app = Flask(__name__)

app.config['DEBUG'] = True
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = timedelta(seconds=1)

@app.route('/')
def home():
    return render_template('home.html', selectedMenu="Home")

@app.route('/infos')
def infos():
    return render_template('infos.html', selectedMenu="Infos")

@app.route('/sequences')
def sequences():
    res = main.id_list
    names = main.name_gene
    return render_template('sequences.html', idlist=res,namegene=names, selectedMenu="Sequences")


@app.route('/sequences',methods=['GET','POST'])
def getSequences():
    if request.method == "POST":
        if request.form["submit"] == 'submit':
            id_list = request.form.getlist('getSequences')
            main.dirName = main.get_fasta(id_list)
    res=main.dirName
    return render_template('alignment.html', res=res, selectedMenu="Alignment")

@app.route('/alignment', methods=['POST'])
def alignmentAndPhylogeny():
    if request.method == "POST":
        if request.form["fAli"] == "Clustal":
            main.clustal_alignment("multifasta.fasta", "alignment.fasta")
            type = 'clustal'
        elif request.form["fAli"] == "Muscle":
            main.muscle_alignment("multifasta.fasta", "alignment.fasta")
            type = 'fasta'

        if request.form["tree"] == "Neighbor Joining":
            main.NJ_tree("alignment.fasta", type)
        elif request.form["tree"] == "Maximum Likelihood":
            main.ML_tree("alignment.fasta", "multiple_sequences_alignment", type)
        res=main.dirName
        return render_template('phylogeny.html', res=res, selectedMenu="Phylogeny")


if __name__ == '__main__':
    app.debug = True
    app.run(host='0.0.0.0', port=8080)
