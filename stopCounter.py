#!/usr/bin/env python

import click
from Bio import Seq, SeqIO
import pandas as pd
from os import path


@click.command()
@click.option("-a", help="Alignment file including reference", type=str,default="run69_aln_ref.fasta", show_default=True)
@click.option("-g", help="Genbank file for reference", type=str,default="../dnds/NC_045512.2-sib.gb", show_default=True)
@click.option("-o", help="Output table file", type=str,default="codonCount.csv", show_default=True)
@click.option("-b", help="Buggy sequence file", type=str, default="Buggy_sequences.fasta", show_default=True)
@click.option("-r", help="Reference ID", type=str, default="NC_045512.3", show_default=True)
def run(a, o, g, b, r):
    # Created by Anmol
    """
    Considering that in the alignment file, reference doesn't have any gap.
    """
    if not a:
        exit("Alignemt file not given. Exiting ....")
    elif not g:
        exit("Genbank file not given. Exiting . . . .")
    elif not o:
        exit("Output file path not given. Exiting . . . .")
    else:
        if not path.isfile(a):
            exit("Given alignment path is not a file or doesn't exists. Exiting . . . . ")
        if not path.isfile(g):
            exit("Given genbank file path is not a file or doesn't exists. Exiting . . . . ")
    
    sequences = {}
    for rec in SeqIO.parse(a, "fasta"):
        sequences[rec.id.split("_re")[0]] = list(str(rec.seq).upper())
    if r not in sequences:
        exit("Reference not in alignment file. . . . . ")
        
    try:
        sequences = pd.DataFrame(sequences)
        if '-' in sequences[r].values:
            seqt = sequences[sequences[r]=="-"].T
            ids2remove = []
            for col in seqt:
                ids2remove += list(seqt[seqt[col] != "-"].index)
            ids2remove = list(set(ids2remove))
            with open(b, "w") as fout:
                for id_ in ids2remove:
                    tseq = "".join(sequences[id_].values)
                    fout.write(">%s\n%s\n"%(id_, tseq))
                    print("%s written in %s for recheck" % (id_, b))
            sequences = sequences.drop(ids2remove, axis=1).drop(seqt.columns) # Recheck this part to avoid any error

        sequences[sequences=="-"]="N"
    except:
        exit("Given file is not alignment. Please correct file")
    
    codon_pos = {"pos":[], "codon":[]}
    #Convert for negative strand


    stop_codon_pos = []
    
    for seq_record in SeqIO.parse(g, "genbank"):
        for feature in seq_record.features:
            if feature.type == "CDS":
                for p in feature.location.parts:
                    ttab = sequences.loc[p.start.position: p.end.position-1]
                    for i in range(0, ttab.shape[0], 3):
                        tcodon = ttab.iloc[i:i+3]
                        amino = tcodon.apply(lambda x: ''.join(x.values)).apply(lambda x:str(Seq.Seq(x).translate()))
                        if "*" in amino.values:
                            db_pos = tcodon.index[0]
                            amino = amino.reset_index().rename(columns={0:db_pos}).T.drop("index")
                            stop_codon_pos.append(amino)

    stop_codon_pos = pd.concat(stop_codon_pos)
    stop_codon_pos.columns = sequences.columns
    stop_codon_pos.loc["StopCount"]=stop_codon_pos.apply(lambda x: list(x.values).count('*'))
    stop_codon_pos.to_csv(o)

if __name__ == '__main__':
    run()
