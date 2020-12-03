from flask import Flask, render_template, request
from datetime import timedelta
import main

app = Flask(__name__)

app.config['DEBUG'] = True
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = timedelta(seconds=1)

@app.route('/')
def mainAcceuil():
    return render_template('Acceuil.html', selectedMenu="Accueil")

@app.route('/infos')
def infos_caract():
    return render_template('infos.html', selectedMenu="Infos")

@app.route('/bdd')
def geneId():
    res = main.id_list
    return render_template('bdd.html', idlist=res, selectedMenu="Base de données")


@app.route('/bdd',methods=['GET','POST'])
def choixGene():
    if request.method == "POST":
        if request.form["submit"] == 'submit':
            id_list = request.form.getlist('choixGene')
            main.get_fasta(id_list)

    return render_template('Alignement.html')

# @app.route('/Alignement')
# def Alignement():
#     return render_template('Alignement.html', selectedMenu="Alignement")


@app.route('/Alignement', methods=['POST'])
def Alignement():
    if request.method == "POST":
        if request.form["fAli"] == "clustal_alignment":
            main.clustal_alignment("multifasta.fasta", "obtenu.fasta")
            type = 'clustal'
        elif request.form["fAli"] == "Muscle_alignement":
            main.muscle_alignment("multifasta.fasta", "obtenu.fasta")
            type = 'fasta'

        if request.form["tree"] == "NJ_tree":
            main.NJ_tree("obtenu.fasta", type)
            # return render_template('NJ_tree.html', selectedMenu="tree")
        elif request.form["tree"] == "ML_tree":
            main.ML_tree("obtenu.fasta", "msa_muscle", type)

        return render_template('tree.html')
        # if request.form["fAli"] == "clustal_alignment":
        #         main.clustal_alignment("multifasta.fasta","msa_clustal.fasta")
        #         main.NJ_tree("msa_clustal.fasta", 'clustal')
        #         return render_template('NJ_tree.html', selectedMenu="tree")
        #
        # elif request.form["fAli"] == "clustal_alignment":
        #         main.clustal_alignment("multifasta.fasta","msa_clustal.fasta")
        #         main.ML_tree("msa_clustal.fasta", "msa_muscle", 'clustal')
        #         return render_template('ML_tree.html', selectedMenu="tree")

        # elif request.form["fAli"] == "Muscle_alignement":
        #     if request.form["tree"] == "NJ_tree":
        #         main.muscle_alignment("multifasta.fasta","msa_muscle.fasta")
        #         main.NJ_tree("msa_muscle.fasta", 'fasta')
        #         return render_template('NJ_tree.html', selectedMenu="tree")
        #
        # elif request.form["fAli"] == "Muscle_alignement":
        #     if request.form["tree"] == "ML_tree":
        #         main.muscle_alignment("multifasta.fasta","msa_muscle.fasta")
        #         main.ML_tree("msa_muscle.fasta", "msa_muscle", 'fasta')
        #         return render_template('ML_tree.html', selectedMenu="tree")

if __name__ == '__main__':
    app.debug = True
    app.run(host='0.0.0.0', port=8080)
