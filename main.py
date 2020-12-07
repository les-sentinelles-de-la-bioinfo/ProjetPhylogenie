import sys
import time
import random
import string
from Bio import Entrez, AlignIO, SeqIO, Phylo
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import PhyloXML
from Bio.Phylo.Applications import PhymlCommandline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

def get_random_string(length):
    """Générer une chaîne aléatoire de longueur fixe"""
    str = string.ascii_lowercase
    return ''.join(random.choice(str) for i in range(length))

def get_fasta(id_list):
    # write a fasta file for each gene
    filename_list = []
    for gene_id in id_list:
        records = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
        filename = (gene_id + '.fasta').format(records)
        filename_list.append("static/data/sauvegardes/" + dirName + filename)
        print('Writing:{}'.format(filename))
        time.sleep(1)
        with open('static/data/sauvegardes/' + dirName + filename, 'w') as f:
            f.write(records.read())
            f.close()

    # combine files in multifasta file
    with open('static/data/sauvegardes/'+ dirName + 'multifasta.fasta', 'w') as outfile:
        for fname in filename_list:
            with open(fname) as infile:
                outfile.write(infile.read())


def clustal_alignment(infile, outfile):
    # create an alignment file with clustal omega
    if (user_OS == 'darwin'):
        clustal_exe = "static/tools/MacOS/clustal-omega-1.2.3-macosx"
    if (user_OS == 'linux'):
        clustal_exe = "static/tools/Linux/clustalo-1.2.4-Ubuntu-x86_64"
    if (user_OS == 'win32'):
        clustal_exe = current_path + "/static/tools/Windows/clustal-omega-1.2.2-win64/clustalo.exe"

    cline = ClustalwCommandline(clustal_exe, infile="static/data/sauvegardes/" + dirName + infile,
                                outfile="static/data/sauvegardes/" + dirName + outfile)
    stdout, stderr = cline()

def muscle_alignment(infile, outfile):
    #create an alignment file with muscke
    if(user_OS == 'darwin'):
        muscle_exe = "static/tools/MacOS/muscle3.8.31_i86darwin64"
    if(user_OS == 'linux'):
        muscle_exe = "static/tools/Linux/muscle3.8.31_i86linux64"
    if(user_OS == 'win32'):
        muscle_exe = current_path + "/static/tools/Windows/muscle3.8.31_i86win32.exe"

    in_file = "static/data/sauvegardes/" + dirName + infile
    out_file = "static/data/sauvegardes/" + dirName + outfile
    muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
    stdout, stderr = muscle_cline()

def NJ_tree(infile, file_type):
    #Tree creation with neighbor-joining
    filename = "static/data/sauvegardes/" + dirName + infile
    aln = AlignIO.read(filename, file_type) #clustal si alignement clustal, fasta si alignement fasta
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, 'nj') # nj ou UPGMA
    tree = constructor.build_tree(aln)
    # print(tree)
    #display a tree on terminal
    #Phylo.draw_ascii(tree)
    tree.ladderize()
    Phylo.draw(tree, do_show=False)
    Phylo.write(tree, 'static/data/sauvegardes/' + dirName + 'tree.txt', "newick")
    foo = current_path + "static/data/sauvegardes/" + dirName + 'tree.png'
    plt.savefig(foo)


def ML_tree(infile, outfile, file_type):
    # Tree creation with maximum-likelihood algorithm (phyML)
    # input : infile = .fasta alignment file that the user can import or paste, outfile = name of output file, file_type = clustal is the clustal too has been used, fasta if muscle tool has been used
    # output : .newick file and .png picture to display
    # phylogeny page should allow to choose maximum likelihood method

    # convert file to phylip
    records = SeqIO.parse("static/data/sauvegardes/" + dirName + infile, file_type)  # clustal <-> fasta
    count = SeqIO.write(records, "static/data/sauvegardes/" + dirName + outfile + ".phylip", "phylip")
    print("Converted %i records" % count)

    if (user_OS == 'darwin'):
        cmd = PhymlCommandline(cmd='static/tools/MacOS/PhyML-3.1/PhyML-3.1_macOS-MountainLion',
                               input='static/data/sauvegardes/' + dirName + outfile + '.phylip')
    if (user_OS == 'linux'):
        cmd = PhymlCommandline(cmd='static/tools/Linux/PhyML-3.1/PhyML-3.1_linux64',
                               input='static/data/sauvegardes/' + dirName + outfile + '.phylip')
    if (user_OS == 'win32'):
        cmd = PhymlCommandline(cmd= current_path + '/static/tools/Windows/PhyML-3.1/PhyML-3.1_win32.exe',
                               input='static/data/sauvegardes/' + dirName + outfile + '.phylip')

    out_log, err_log = cmd()
    tree = Phylo.read('static/data/sauvegardes/' + dirName + outfile + '.phylip_phyml_tree.txt', 'newick')
    Phylo.draw(tree, do_show=False)
    Phylo.write(tree, 'static/data/sauvegardes/' + dirName + 'tree.txt', "newick")
    foo = current_path + 'static/data/sauvegardes/' + dirName + 'tree.png'
    plt.savefig(foo)

##################################### MAIN ##############################################################
#Contact address
Entrez.email = "sentinelles.bioinfo@gmail.com"

# Define user's operating system
user_OS = sys.platform

current_path = os.path.dirname(__file__)
parent_path = os.path.dirname(current_path)

# list of genes
id_list = ["AY158636.1","AY158639.1","AY159811.1","AY159808.1","AY159809.1","AY158637.1","AY159810.1"]
name_gene = ["Vipera berus Pla2Vb", "Vipera berus AmtI2", "Vipera berus AmtI1", "Vipera aspis AmtI1", "Vipera aspis AmtI1", "Vipera aspis (AmtI2)", "Vipera aspis zinnikeri AmtI1"]

# Create a random name for each user session and create corresponding directories in order to allow multiple simultaneous uses without loss of data
dirName = get_random_string(10) + "/"
saveDir = "static/data/sauvegardes/"
os.makedirs(saveDir + dirName, exist_ok=True)
os.makedirs("static/figure/" + dirName, exist_ok=True)

## Function calls to use the program in the terminal without the user interface
#get_fasta(id_list)
#clustal_alignment("multifasta.fasta","msa_clustal.fasta")
#muscle_alignment("multifasta.fasta","msa_muscle.fasta")
#NJ_tree("msa_clustal.fasta", "clustal")
#ML_tree("msa_clustal.fasta", "msa_muscle", "clustal")
