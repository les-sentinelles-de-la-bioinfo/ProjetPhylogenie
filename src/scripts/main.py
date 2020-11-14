import os
from Bio import Entrez
from Bio import AlignIO
from Bio import SeqIO
from Bio import Phylo
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo import PhyloXML
from Bio.Phylo.Applications import PhymlCommandline
import time


#donne email en cas de problème
Entrez.email = "sentinelles.bioinfo@gmail.com"

#liste des gènes à récupérer
id_list = ["MT298507.1", "HQ954792.1", "MK013995.1","KJ128666.1","JX218056.1","KP403684.1","HQ536294.1","MK013988.1"]
filename_list = []

#écrit un fichier fasta pour chaque gène
n = 0
for gene_id in id_list:
    n = n + 1
    records = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
    filename = (gene_id + '.fasta').format(records)
    filename_list.append("../../data/sauvegardes/"+filename)
    print('Writing:{}'.format(filename))
    time.sleep(1)
    with open('../../data/sauvegardes/'+filename, 'w') as f:
        f.write(records.read())
        f.close()

#lis la séquence
#with open('test.gb', 'r') as f2:
#    data = f2.read()
#    #print(data)

#combine les fichiers
with open('../../data/sauvegardes/multifasta.fasta', 'w') as outfile:
    for fname in filename_list:
        with open(fname) as infile:
            outfile.write(infile.read())


#crée un fichier fasta d'alignement avec clustal
clustal_exe = "../../etc/tools/MacOS/clustal-omega-1.2.3-macosx"
#clustal_exe = "../../etc/tools/Linux/clustalo-1.2.4-Ubuntu-x86_64"
#clustal_exe = "../../etc/tools/Windows/clustal-omega-1.2.2-win64/clustalo.exe"

cline = ClustalwCommandline(clustal_exe, infile="../../data/sauvegardes/multifasta.fasta", outfile="../../data/sauvegardes/msa_clustal.clustal")
stdout, stderr = cline()


#crée un fichier fasta d'alignement avec muscle
muscle_exe = "../../etc/tools/MacOS/muscle3.8.31_i86darwin64"
#muscle_exe = "../../etc/tools/Linux/muscle3.8.31_i86linux64"
#muscle_exe = "../../etc/tools/Windows/muscle3.8.31_i86win32.exe"

in_file = "../../data/sauvegardes/multifasta.fasta"
out_file = "../../data/sauvegardes/msa_muscle.fasta"
muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
stdout, stderr = muscle_cline()


#création de l'arbre avec NJ ou UPGMA + clustal alignment
filename = "../../data/sauvegardes/msa_clustal.clustal"
aln = AlignIO.read(filename, 'clustal') #changer clustal par fasta si alignement muscle
#for alignment in AlignIO.parse(filename, "clustal"):
#    print("Alignment of length %i" % alignment.get_alignment_length())
print(aln)
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
print(dm)
constructor = DistanceTreeConstructor(calculator, 'nj') # nj ou UPGMA
tree = constructor.build_tree(aln)
print(tree)
#dessine un arbre sur le terminal
#Phylo.draw_ascii(tree)
tree.ladderize()
Phylo.draw(tree)


#création de l'arbre avec ML
#convert file to phylip
records = SeqIO.parse("../../data/sauvegardes/msa_clustal.clustal", "clustal") # clustal <-> fasta
count = SeqIO.write(records, "../../data/sauvegardes/msa_clustal.phylip", "phylip")
print("Converted %i records" % count)
cmd = PhymlCommandline(cmd='../../etc/tools/MacOS/PhyML-3.1/PhyML-3.1_macOS-MountainLion', input='../data/sauvegardes/msa_clustal.phylip')
#cmd = PhymlCommandline(cmd='../../etc/tools/Linux/PhyML-3.1/PhyML-3.1_linux64', input='../data/sauvegardes/msa_clustal.phylip')
#cmd = PhymlCommandline(cmd='../../etc/tools/Windows/PhyML-3.1/PhyML-3.1_win32.exe', input='../data/sauvegardes/msa_clustal.phylip')
out_log, err_log = cmd()
