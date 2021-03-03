# import click
from numpy import fabs
from Bio import SeqIO
import pandas as pd
from os import path, system
from sys import argv

seq_to_insert = {}
with open("Trimmed.fasta", "w") as fout:
    for rec in SeqIO.parse(argv[1], "fasta"):
        seq_id = rec.id.split("_")[0]
        seq = str(rec.seq).strip("N")
        seq_to_insert[seq_id] = seq
        fout.write(">%s\n%s\n"%(seq_id, seq))

def most_frequent(olst):
    lst = []
    for k in olst:
        if k in ["N", '-']:
            continue
        lst.append(k)
    if not lst:
        return 'N'
    
    return max(set(lst), key = lst.count) 

refs = {}
for rec in SeqIO.parse(argv[2], "fasta"):
    refs[rec.id] = list(str(rec.seq).upper())
refs = pd.DataFrame(refs)
refs = "".join(refs.apply(lambda x: most_frequent(list(x.values)), axis=1))

with open("ref.fasta", "w") as fout:
    fout.write(f">ref\n{refs}\n")
# print(refs)




system(f"blat -noHead ref.fasta Trimmed.fasta Trimmed.psl")
matches = pd.read_table("Trimmed.psl", header=None).sort_values(0, ascending=False).drop_duplicates(9)[[0,1,3,9,10,11,12,13,14,15,16]]

# matches
# matches = matches.loc[matches[9]==matches[13], [0,1,3,9,10,11,12,13,14,15,16]]

original_seq = {}
for rec in SeqIO.parse(argv[2], "fasta"):
    original_seq[rec.id] = str(rec.seq)


with open(argv[3], "w") as fout:
    reported = []
    for _, row in matches.iterrows():
        reported.append(row[9])
        try:
            if row[10]==row[12]:
                if row[11] == 0:
                    fout.write(">%s\n%s\n%s\n%s\n" % (row[9], original_seq[row[9]][:row[15]],seq_to_insert[row[9]],original_seq[row[9]][row[16]:]))
                else:
                    fout.write(">%s\n%s\n%s\n%s\n" % (row[9], original_seq[row[9]][:row[16]-row[10]],seq_to_insert[row[9]],original_seq[row[9]][row[16]:]))
            else:
                fout.write(">%s\n%s\n%s\n%s\n" % (row[9], original_seq[row[9]][:row[15]],seq_to_insert[row[9]],original_seq[row[9]][row[15]+row[10]:])) 
        except KeyError:
            continue
    for k in original_seq:
        if k in reported:
            continue
        fout.write(">%s\n%s\n"%(k, original_seq[k]))


system("rm Trimmed.* ref.fasta")













# def sanger_merge():
#     pass


# @click.command()
# @click.option("-a", help="Original Alignemt file, along with reference")
# # @click.option("-r", help="Reference Ids in the alignemnt")
# @click.option("-s", help="Sanger Fasta")
# @click.option("-o", help="Output fasta")
# @click.option("-b", help="blat path")

# def run(a,s,o,b):
#     '''
#     This code insert the sanger sequence in alignment or given sequences.
#     Sequences must overlapping region to sanger sequences. And sanger sequences should be unique in the genome.
#     blat must be in system path. If not installed, please download from USCS (https://www.youtube.com/watch?v=u5hS8vzmQXI). 
#     '''
#     if not path.exists(s):
#         exit("Sanger Sequence file path is not correct. Exiting . . . . ")
#     try:
#         pass
#     except:
#         exit("Given file is not in fasta")
        

#     pass

# if __name__ == '__main__':
#     run()






