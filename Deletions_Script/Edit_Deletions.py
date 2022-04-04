from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys

#Global variable
MolHandle = re.compile("AGCGGGAGACCGGGGTCTCTGAGC")
DelSize=3

#Definitions
def readInFile(inFileName):
    seqs=[]
    names=[]
    with open(inFileName) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqs.append(record.seq)
            names.append(record.id)
            names=[name.replace("insertion","deletion") for name in names]
    return(names,seqs)

def editSeqs(seqs):
    matches=[MolHandle.search(str(seq[i])) for i in range(len(seq))]
    boundaries=[(match.start(0),match.end(0)) for match in matches]
    newSeqs=[str(seq[i][0:boundaries[i][0]]+seq[i][(boundaries[i][1]+DelSize):len(seq[i])]) for i in range(len(seq))]
    return(newSeqs)

def outputSeqs(newSeqs,names,OFname):

    sequences=[SeqRecord(Seq(newSeqs[i]),id=names[i],description="") for i in range(len(newSeqs))]

    with open(OFname, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

#Main
#"./Oligopool_with_ATG_inclusion/EV-D8_Fermon_Insertion_library_300_oligo_len_MOD_TO_INCLUDE_ATG/Capsid_DI_Oligos.fasta"
if __name__ == "__main__":
    names, seq = readInFile(sys.argv[1])
    newSeqs=editSeqs(seq)
    outputFileName=sys.argv[2] 
    outputSeqs(newSeqs,names,outputFileName)
