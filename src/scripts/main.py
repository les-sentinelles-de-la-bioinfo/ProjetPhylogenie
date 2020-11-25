import sys
import time
from Bio import Entrez, AlignIO, SeqIO, Phylo
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import PhyloXML
from Bio.Phylo.Applications import PhymlCommandline
from matplotlib import pyplot as plt


#Contact address
Entrez.email = "sentinelles.bioinfo@gmail.com"

#Define user's operating system
user_OS = sys.platform

#list of genes
id_list_1 = ["MT298507.1", "HQ954792.1", "MK013995.1","KJ128666.1","JX218056.1","KP403684.1","HQ536294.1","MK013988.1"]
id_list_2 = ["AY158636.1","AY158639.1","AY159811.1","AY159808.1","AY159809.1","AY158637.1","AY159810.1"]

def get_fasta(id_list):
    #write a fasta file for each gene
    #input : list of gene id, the user can choose which list to use
    #output : .fasta multifasta file which is displayed on the website and can be downloaded
    
    filename_list = []
    for gene_id in id_list:
        records = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
        filename = (gene_id + '.fasta').format(records)
        filename_list.append("../../data/sauvegardes/"+filename)
        print('Writing:{}'.format(filename))
        time.sleep(1)
        with open('../../data/sauvegardes/'+filename, 'w') as f:
            f.write(records.read())
            f.close()

    #combine files in multifasta file
    with open('../../data/sauvegardes/multifasta.fasta', 'w') as outfile:
        for fname in filename_list:
            with open(fname) as infile:
                outfile.write(infile.read())


def clustal_alignment(infile, outfile):
    #create an alignment file with clustal omega
    #input : infile = multifasta file from the function get_fasta that the user can import or paste, outfile = name of the output file
    #output : .fasta alignment file, displayed on the website and can be downloaded
    #alignement page should allow to choose clustal tool 

    if(user_OS == 'darwin'):
        clustal_exe = "../../etc/tools/MacOS/clustal-omega-1.2.3-macosx"
    if(user_OS == 'linux'):
        clustal_exe = "../../etc/tools/Linux/clustalo-1.2.4-Ubuntu-x86_64"
    if(user_OS == 'win32'):
        clustal_exe = "../../etc/tools/Windows/clustal-omega-1.2.2-win64/clustalo.exe"

    cline = ClustalwCommandline(clustal_exe, infile="../../data/sauvegardes/"+infile, outfile="../../data/sauvegardes/" + outfile)
    stdout, stderr = cline()


def muscle_alignment(infile, outfile):
    #create an alignment file with muscle
    #input : infile = file from function get_fasta that the user can import or paste, outfile = name of the output file
    #output : .fasta alignment file, displayed on the website and can be downloaded
    #alignment page should allow to choose muscle tool 

    if(user_OS == 'darwin'):
        muscle_exe = "../../etc/tools/MacOS/muscle3.8.31_i86darwin64"
    if(user_OS == 'linux'):
        muscle_exe = "../../etc/tools/Linux/muscle3.8.31_i86linux64"
    if(user_OS == 'win32'):
        muscle_exe = "../../etc/tools/Windows/muscle3.8.31_i86win32.exe"

    in_file = "../../data/sauvegardes/" + infile
    out_file = "../../data/sauvegardes/" + outfile
    muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
    stdout, stderr = muscle_cline()


def NJ_tree(infile, file_type):
    #Tree creation with neighbor-joining
    #input : infile = .fasta alignment file that the user can import or paste, file_type = clustal if the clustal too has been used, fasta if muscle tool has been used
    #output : .png picture to display
    #phylogeny page should allow to choose neighbor-joining method

    filename = "../../data/sauvegardes/" + infile
    aln = AlignIO.read(filename, file_type) #clustal if clustal alignment, fasta if fasta alignment
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, 'nj') # nj or UPGMA
    tree = constructor.build_tree(aln)
    #print(tree)
    #display a tree on terminal
    #Phylo.draw_ascii(tree)
    tree.ladderize()
    Phylo.draw(tree, do_show=False)
    plt.savefig('../../data/images/NJ_tree.png')


def ML_tree(infile, outfile, file_type):
    #Tree creation with maximum-likelihood algorithm (phyML)
    #input : infile = .fasta alignment file that the user can import or paste, outfile = name of output file, file_type = clustal is the clustal too has been used, fasta if muscle tool has been used
    #output : .newick file and .png picture to display
    #phylogeny page should allow to choose maximum lieklihood method

    #convert file to phylip
    records = SeqIO.parse("../../data/sauvegardes/" + infile, file_type) # clustal <-> fasta
    count = SeqIO.write(records, "../../data/sauvegardes/" + outfile + ".phylip", "phylip")
    print("Converted %i records" % count)

    if(user_OS == 'darwin'):
        cmd = PhymlCommandline(cmd='../../etc/tools/MacOS/PhyML-3.1/PhyML-3.1_macOS-MountainLion', input='../../data/sauvegardes/' + outfile + '.phylip')
    if(user_OS == 'linux'):
        cmd = PhymlCommandline(cmd='../../etc/tools/Linux/PhyML-3.1/PhyML-3.1_linux64', input='../../data/sauvegardes/' + outfile + '.phylip')
    if(user_OS == 'win32'):
        cmd = PhymlCommandline(cmd='../../etc/tools/Windows/PhyML-3.1/PhyML-3.1_win32.exe', input='../../data/sauvegardes/' + outfile + '.phylip')
    
    out_log, err_log = cmd()
    tree = Phylo.read('../../data/sauvegardes/' + outfile + '.phylip_phyml_tree.txt', 'newick')
    Phylo.draw(tree, do_show=False)
    plt.savefig('../../data/images/ML_tree.png')




############### MAIN ###############


#get_fasta(id_list_1)
#clustal_alignment("multifasta.fasta","msa_clustal.fasta")
#muscle_alignment("multifasta.fasta","msa_muscle.fasta")
#NJ_tree("msa_muscle.fasta", 'fasta')
#ML_tree("msa_muscle.fasta", "msa_muscle", "fasta")
