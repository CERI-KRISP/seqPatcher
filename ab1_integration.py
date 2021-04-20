#!/usr/bin/env python

# Author: Anmol Kiran
# Affiliation: University of Liverpool. UK and Malawi-Liverpool-Wellcome Trust, Malawi
# V0.1: 20/04/2021



import click
from Bio import SeqIO, AlignIO, Seq
from glob import glob
from collections import defaultdict
from click.decorators import command
import pandas as pd
import numpy as np
from os import makedirs, path, system
import tempfile as tmf
import subprocess as sb
from ruffus import * # 2.8.4
# from abx import *




def cmd(command):
    """Runs all the command using Popen communicate."""
    sb.Popen(command, stdout=sb.DEVNULL).communicate()



# def most_frequent(olst):
#     """Returns most frequent occuring value from given list. Ignoring ambiguous and -s. """
#     lst = []
#     for k in olst:
#         if k in "ACGT":
#             lst.append(k)
#     if not lst:
#         return 'N'
#     return max(set(lst), key = lst.count)



# def dummy_ref(alnfile, dummy_ref_file):
#     """Generate dummy reference sequence based most frequent occuring nucleotide in at the position in the alignment"""
#     refs = {}
#     for rec in SeqIO.parse(alnfile, "fasta"):
#         refs[rec.id] = list(str(rec.seq).upper())
#     refs = pd.DataFrame(refs)
#     refs = "".join(refs.apply(lambda x: most_frequent(list(x.values)), axis=1))
#     with open(dummy_ref_file, "w") as fout:
#         fout.write(f">ref\n{refs}\n")




# def integrate_seq(infiles, outfile):
#     """Inserts the sanger sequences in a given alignments"""
#     tmp_fold = infiles[2]
    
    
#     # Aligned sequences
#     original_seq = {}
#     for rec in SeqIO.parse(infiles[0], "fasta"):
#         original_seq[rec.id] = str(rec.seq)

#     if not original_seq: # Generate empty fasta
#         command=["touch", outfile]
#         cmd(command)
#         return
#     else:
        
#         flb = path.split(infiles[0])[1].split(".")[0]
#         f_dummy = f"{tmp_fold}/{flb}_dummy.fasta"
#         dummy_ref(infiles[0], f_dummy)

#     # Cleaned sanger sequences
#     seq_to_insert = {}
#     for rec in SeqIO.parse(infiles[1], "fasta"):
#         seq_to_insert[rec.id] = rec.seq
    
    
#     # mapping sanger sequences
#     p_dummy = f"{tmp_fold}/{flb}_dummy.psl"
#     command = ["blat", "-noHead", f_dummy, infiles[1], p_dummy]
#     cmd(command)
#     matches = pd.read_table(p_dummy, header=None).sort_values(0, ascending=False).drop_duplicates(9)[[0,1,3,9,10,11,12,13,14,15,16]] # Considering that all the sequences are on positive strand
#     matches["q_dummy"] = matches[9].apply(lambda x: x.split("_")[0])
#     matches = matches.sort_values('q_dummy')
    
    

#     # matches
#     # matches = matches.loc[matches[9]==matches[13], [0,1,3,9,10,11,12,13,14,15,16]]
#     # TODO: Create dummy insert ids for sanger data

    
    
#     # TODO: Reverse complement value for peak
#     # TODO: Add Strand and reverse strand position in the sequences, after trimming
    
    
    
    


#     with open(outfile, "w") as fout:
#         reported = []
#         for q_dummy in set(matches['q_dummy'].values):
#             if q_dummy in original_seq:
#                 reported.append(q_dummy)
#             else:
#                 continue
#             # TODO: Add temperory seq variable
#             tseq = original_seq[q_dummy]
#             for _, row in matches[matches['q_dummy']==q_dummy].iterrows():
#                 print(row)
#                 try:
#                     t_ins = str(seq_to_insert[row[9]]) if row[9]=="+" else str(seq_to_insert[row[9]].reverse_complement()) 
#                     if row[10]==row[12]:
#                         if row[11] == 0:
                            
#                             tseq = tseq[:row[15]]+ t_ins+ tseq[row[16]:]
#                             # fout.write(">%s\n%s\n%s\n%s\n" % (row["q_dummy"], original_seq[row["q_dummy"]][:row[15]],seq_to_insert[row["q_dummy"]],original_seq[row["q_dummy"]][row[16]:]))
#                         else:
#                             tseq = tseq[:row[16]-row[10]]+t_ins+tseq[row[16]:]
#                             # fout.write(">%s\n%s\n%s\n%s\n" % (row["q_dummy"], original_seq[row["q_dummy"]][:row[16]-row[10]],seq_to_insert[row["q_dummy"]],original_seq[row["q_dummy"]][row[16]:]))
#                     else:
#                         tseq = tseq[:row[15]]+t_ins+tseq[row[15]+row[10]:]
#                         # fout.write(">%s\n%s\n%s\n%s\n" % (row["q_dummy"], original_seq[row["q_dummy"]][:row[15]],seq_to_insert[row["q_dummy"]],original_seq[row["q_dummy"]][row[15]+row[10]:]))
#                 except KeyError:
#                     continue
#             fout.write(f">{q_dummy}\n{tseq}\n")
#         for k in original_seq:
#             if k in reported:
#                 continue
#             fout.write(">%s\n%s\n"%(k, original_seq[k]))
            
#     #TODO: Create Dummy files in case it fails with a warning output


#     # system("rm Trimmed.* ref.fasta")


# def reverse_check():
#     """To Check whether the sequences are opposite strand."""





def useful_short_range(lst, gap=10):
    def ranges(p):
        q = sorted(p)
        i = 0
        for j in range(1,len(q)):
            if q[j] > 1+q[j-1]:
                yield (q[i],q[j-1])
                i = j
        yield (q[i], q[-1])
    trange =  list(ranges(lst))
    # print(trange)
    
    myrange = []
    for rng in trange:
        d = rng[1]-rng[0]+1
        if (d<gap) and d%3==0:
            myrange.append(rng)
    return myrange




def useful_range(lst, gap=10):
    # print(gap)
    def ranges(p, gap):
        q = sorted(p)
        i = 0
        for j in range(1,len(q)):
            if q[j] > gap+q[j-1]:
                yield (q[i],q[j-1])
                i = j
        yield (q[i], q[-1])
    trange =  list(ranges(lst, gap))
    i = 0
    myrange = trange[0]
    for rng in trange:
        d = rng[1]-rng[0]+1
        if i < d:
            i = d
            myrange = rng
    return myrange


# def final_base(ref, fd, rv=None, bsc=None):
#     """"Returns final base for distinct base sites."""
#     if ref == "-":
#         return "" # TODO: Check for multiple of 3 inserts
#     elif fd=="-" or rv=="-":
#         return ref
#     elif ((fd==rv) or (rv==None)) and fd in "ACGT":
#         return fd
#     elif fd in amb_bases and rv in "ACGT":
#         return rv
#     else:
#         return fd








                
#     #             fout.write(f">{flb}\n{nuc_seq}\n")
#     #     for rec in SeqIO.parse("tmpfold/CoVid_S_Gene.fasta", "fasta"):
#     #             fout.write(f">ref\n{rec.seq}")
#     # command =["muscle", "-in", f"{tmp_fold}/{flb}.fasta", "-out", f"{tmp_fold}/{flb}.aln.fasta"]
#     # cmd(command)
#     # tseq = {}
#     # for rec in SeqIO.parse(f"{tmp_fold}/{flb}.fasta", "fasta"):
#     #     tseq[rec.id] = list(str(rec.seq).upper())
#     # tseq = pd.DataFrame(tseq)
#     # tseq["gap"] = tseq.apply(lambda x : '-' if '-' in x.values else '0', axis=1)
#     # tseq["alnpos"] = list(range(len(tseq)))
#     # core_aln = np.where(tseq.apply(lambda x : '-' if '-' in x.values else '0', axis=1).values!='-')[0]
#     # my_range = useful_range(core_aln)
#     # tseq = tseq.loc[my_range[0]:my_range[1]] # Reverse the values in case of reverse complement
#     # return tseq # TODO: Need to be retrun some other values as well for the final correction 


def df_reverse_complement(df):
    """Reverse Complement of the nucleotide dataframe."""
    df = df.loc[::-1]
    def rvc(x):
        return  str(Seq.Seq(x).complement())
    def rvc_dict(x):
        xc = {}
        for k in x:
            xc[rvc(k)] = x[k]
        return xc
    df["nuc"] = df["nuc"].apply(rvc)
    df["peak"] = df["peak"].apply(rvc_dict)
    return df





def ab1seq(fl, tmp_fold):
    """"ab1 to seq trimmed based on reference."""
    channels = ["DATA9", # G
                "DATA10", # A
                "DATA11", # T
                "DATA12" # C
                ]
    amb_bases = { # All ambiguous nucleotides
            "N":set(["A", "C", "G", "T"]),
            "R":set(["A", "G"]),
            "Y":set(["T", "C"]),
            "K":set(["G", "T"]),
            "M":set(["A", "C"]),
            "S":set(["G", "C"]),
            "W":set(["A", "T"]),
            "B":set(["C", "G", "T"]),
            "D":set(["A", "G", "T"]),
            "H":set(["A", "C", "T"]),
            "V":set(["A", "C", "G"]),
        }

    bases = {"DATA9":"G", "DATA10":"A", "DATA11":"T", "DATA12":"C"}
    
    record = SeqIO.read(fl, "abi")
    trace = defaultdict(list)

    for c in channels:
        trace[bases[c]] = record.annotations["abif_raw"][c]

    nuc_df = {"nuc":[], "peak":[]}
    
    for c in zip(record.annotations["abif_raw"]["PBAS1"],record.annotations["abif_raw"]["PLOC1"]):
        base = chr(c[0])
        nuc_df["nuc"].append(base)
        if base in amb_bases:
            td = {}
            for b in amb_bases[base]:
                td[b]=trace[b][c[1]]
            nuc_df["peak"].append(td)
        else:
            nuc_df["peak"].append({base:trace[base][c[1]]})
    nuc_df = pd.DataFrame(nuc_df)
    peak = nuc_df["peak"].apply(lambda x: np.mean(list(x.values()))).values

    correct_peak_pos = np.where(peak>100)[0] 
    min_pos, max_pos = np.min(correct_peak_pos), np.max(correct_peak_pos) # 1 is added as last position is not included
    nuc_df = nuc_df.loc[min_pos:max_pos]

    nuc_seq = "".join(nuc_df["nuc"].values)
    
    flb = path.split(fl)[1].split(".ab1")[0] # Consitering that the sample name is the first part
    # print(nuc_seq, flb)
    
    with open(f"{tmp_fold}/{flb}.fasta", "w")  as fout:
        fout.write(f">{flb}\n{nuc_seq}\n")
    
    command =["blat", "-noHead", f"{tmp_fold}/CoVid_S_Gene.fasta", f"{tmp_fold}/{flb}.fasta", f"{tmp_fold}/{flb}.psl"]# Mapping against reference
    cmd(command)
    df = pd.read_table(f"{tmp_fold}/{flb}.psl", header=None).sort_values(0, ascending=False).drop_duplicates(9)
    if len(df) != 1:
        exit(f"Looks like {fl} has multiple sequences. Please check. Exiting . . . . . ")
    
    for _, row in df.iterrows():
        if row[8]=="-":
            nuc_df = df_reverse_complement(nuc_df)
    return nuc_df





def aln_df_with_ref(seq_dict, flb, tmp_fold):
    """Alignment dataframe"""
    inf = f"{tmp_fold}/{flb}.in.fasta"
    otf = f"{tmp_fold}/{flb}.out.fasta"
    with open(inf, "w") as fout:
        for k in seq_dict:
            fout.write(f">{k}\n{seq_dict[k]}\n")
    command = ["muscle", "-in", inf, "-out", otf]
    cmd(command)
    sequences = {}
    for rec in SeqIO.parse(otf, "fasta"):
        sequences[rec.id] = list(str(rec.seq).upper())
    return pd.DataFrame(sequences)


def merge_base_peak(df, peak_dict):
    
    df["idx"] = list(df.index)
    

    # df = enumerate_columns(df)
    # print(df)
    for k in peak_dict:
        peak_dict[k][f"{k}_idx"] = list(df[df[f"{k}"]!="-"].index)#list(range(len(peak_dict[k]))
        df = df.merge(peak_dict[k], left_on="idx", right_on=f"{k}_idx", how="outer")
        del df[str(k)]
    # print(df)
    return df







def aln_clean(df):
    one_side = ["1"]
    if df.shape[1] == 4:
        idx = df[df["1_nuc"]!="-"].index
    else:
        one_side.append("2")
        idx = df[~((df["1_nuc"]=="-") | (df["2_nuc"]=="-"))].index
    u_range = useful_range(list(idx), 10) # Add the importnat gap. Default: 10
    for col in df.columns:
        if col=="ref":continue
        df.loc[~df.index.isin(list(range(u_range[0],u_range[1]+1))), col] = "-"
        
    def rep_paired_base(lst):
        # TODO: Can be improved later selective query here
        if lst["ref"] == "-":
            return "-"
        elif ((lst["1_nuc"] == lst["2_nuc"] == "-") or
              (lst["1_nuc"] != lst["2_nuc"] == lst["ref"]) or
              (lst["2_nuc"] != lst["1_nuc"] == lst["ref"]) or
              (lst["2_nuc"] == lst["1_nuc"] == lst["ref"]) or
              (lst["2_nuc"] != lst["1_nuc"] != lst["ref"]) or
              ((lst["2_nuc"] not in 'ACGT-') and  (lst["1_nuc"] not in 'ACGT-'))):
            return lst["ref"]
        elif lst["1_nuc"] == lst["2_nuc"] != lst["ref"]:
            return lst["1_nuc"]
        else:
            if (lst["2_nuc"] not in 'ACGT-') and  (lst["1_nuc"] in 'ACGT'):
                return lst["1_nuc"]
            elif (lst["1_nuc"] not in 'ACGT-') and  (lst["2_nuc"] in 'ACGT'):
                 return lst["2_nuc"]
        
    def rep_single_base(lst):
        if lst["ref"] == "-":
            return "-"
        elif lst["1_nuc"] == "-":
            return lst["ref"]
        elif lst["1_nuc"] not in 'ACGT-':
            pass
            # TODO: Set statistical parameters
#         else:
    df["rep"] = df.apply(rep_paired_base, axis=1)
    return df, u_range
            
        

        
        # pass
    # df["rep"] = ""
    # print(df)
    # df.apply(rep_paired_base, axis=1)
    # print(df)
    # pass


def ab2fasta(infiles, outfile, bc="neigh"): # Base criteria, max, neighbors, mixed # Inputfiles paired and none paired
    """Retains fasta and converts ab1 to fasta"""
    tmp_fold= path.split(outfile)[0]
    is_fasta = False
    
    output_file = open(outfile, "w")

    for fls in infiles:
        print(fls)
        if fls[0].endswith(".fasta"):
            for rec in SeqIO.parse(fl, "fasta"):
                    sequences[rec.id] = rec.seq
            is_fasta = True

        else:
            ab1seq_dfs = {}
            tsequences = {}
            flb = path.split(fls[0])[1].split("_")[0]
            # if flb != "K008703":
            #     continue
            for i, fl in enumerate(fls, start=1):
                ab1seq_dfs[i] = ab1seq(fl, tmp_fold)
                tsequences[i] = ''.join(ab1seq_dfs[i]['nuc'].values)
            for rec in SeqIO.parse(f"{tmp_fold}/CoVid_S_Gene.fasta", "fasta"):
                tsequences[rec.id] = str(rec.seq)
            
            for k in ab1seq_dfs:
                # print(ab1seq_dfs[k].columns)
                ab1seq_dfs[k].columns = [f"{k}_{col}" for col in ab1seq_dfs[k].columns]

            aln_with_peak = merge_base_peak(aln_df_with_ref(tsequences, flb, tmp_fold), ab1seq_dfs)
            xx = aln_with_peak.loc[1294:1299]
            # print(xx[pd.isnull(xx["1_nuc"])], "SSSSSSSS")
            for n in ["1", "2"]:
                cl = f"{n}_nuc"
                if cl in aln_with_peak.columns:
                    print(cl)
                    print(aln_with_peak.loc[pd.isnull(aln_with_peak[cl])].loc[1294:1299])
                    aln_with_peak.loc[pd.isnull(aln_with_peak[cl]),cl] = "-"
                    print(aln_with_peak.loc[1294:1299])
            # print("".join(aln_with_peak["ref"].values))
            # print(aln_with_peak[["ref","2_nuc","1_nuc"]].apply(set, axis=1))
            aln_with_peak, u_range = aln_clean(aln_with_peak)
            # print(aln_with_peak.loc[1294:1299], u_range)
            usr = useful_short_range(aln_with_peak[aln_with_peak[["ref", "1_nuc", "2_nuc"]].apply(lambda x: True if '-' in list(x.values) else False, axis=1)].index, 10)
            # print(usr)
            if usr:
                for us in usr:
                    # print(aln_with_peak.loc[us[0]:us[1], "1_nuc"].values)
                    if "".join(aln_with_peak.loc[us[0]:us[1], "1_nuc"].values) == "".join(aln_with_peak.loc[us[0]:us[1], "2_nuc"].values):
                        aln_with_peak.loc[us[0]:us[1], "rep"]= aln_with_peak.loc[us[0]:us[1], "1_nuc"].values
                        # print("Kiran", aln_with_peak.loc[us[0]:us[1]])
                # exit(f"Contact Anmol with {fls}")
            # print(aln_with_peak[aln_with_peak[["ref", "1_nuc", "2_nuc"]].apply(lambda x: True if '-' in list(x.values) else False, axis=1)].index)
            print(aln_with_peak, u_range)
            seq = "".join(aln_with_peak["rep"].values)
            print(seq)
            output_file.write(f">{flb} {u_range[0]} {u_range[1]}\n{seq}\n")
    output_file.close()
            # return , u_range
            # print("Breaking out", aln_with_peak)
            # exit()
        
            # print(ab1seq_dfs)
        # print(sequences)
        # exit()
#     # TODO: Clean Alignments, insert ref in the dict 
#     if is_fasta:
#         pass
#     else:
#         pass

#                 for k in files_dict_list:
#                     sanger_seq = open(f"{k}.fasta", "w")
#                     for fl in files_dict_list[k]:
#                         record = SeqIO.read(fl, "abi")

#                         trace = defaultdict(list)
#                         for c in channels:
#                             trace[bases[c]] = record.annotations["abif_raw"][c]

#                         peak = []
#                         nuc_seq = []
#                         for c in zip(record.annotations["abif_raw"]["PBAS1"],record.annotations["abif_raw"]["PLOC1"]):
#                             # print(i, chr(c[0]), c[1])
#                             base = chr(c[0])
#                             nuc_seq.append(base) #  Replace base with a value
#                             if base in amb_bases:
#                                 base_values = []
#                                 for b in amb_bases[base]:
#                                     base_values.append(trace[b][c[1]])
#                                 peak.append(np.mean(base_values))
#                             else:
#                                 peak.append(trace[base][c[1]])
#                         peak = np.array(peak)
#                         correct_peak_pos = np.where(peak>100)[0]
#                         min_pos, max_pos = np.min(correct_peak_pos), np.max(correct_peak_pos)+1 # 1 is added as last position is not included
#                         nuc_seq = "".join(nuc_seq)[min_pos:max_pos]
#                         flb = path.split(fl)[1].split(".")[0]
#                         sanger_seq.write(f">{flb}\n{nuc_seq}\n")
#                     sanger_seq.close()
#                     system(f"blat -noHead CoVid_S_Gene.fasta {k}.fasta {k}.psl")
#                     sequences = {}
#                     for rec in SeqIO.parse(f"{k}.fasta", "fasta"):
#                         sequences[rec.id] = rec.seq
#                     blat_map = pd.read_table(f"{k}.psl", header=None)
#                     blat_map = blat_map[blat_map[8]=='-']
#                     for _, row in blat_map.iterrows():
#                         sequences[row[9]] = sequences[row[9]].reverse_complement()
#                     sanger_seq = open(f"{k}.fasta", "w")
#                     for kn in sequences:
#                         sanger_seq.write(f">{kn}\n{sequences[kn]}\n")


#                     for rec in SeqIO.parse("CoVid_S_Gene.fasta", "fasta"):
#                         sanger_seq.write(f">ref\n{rec.seq}")
#                     sanger_seq.close()
#                     system(f"muscle -in {k}.fasta -out {k}.aln.fasta")
#                     tseq = {}
#                     for rec in SeqIO.parse(f"{k}.aln.fasta", "fasta"):
#                         tseq[rec.id] = list(str(rec.seq).upper())
#                     tseq = pd.DataFrame(tseq)
#                     tseq["gap"] = tseq.apply(lambda x : '-' if '-' in x.values else '0', axis=1)
#                     tseq["alnpos"] = list(range(len(tseq)))

#                     core_aln = where(tseq.apply(lambda x : '-' if '-' in x.values else '0', axis=1).values!='-')[0]

#                     core_aln_min, core_aln_max = np.min(core_aln), np.max(core_aln)

#                     xx = core_aln[1:]- core_aln[:-1] # Choose the longest block and correct that
#                     xx = np.insert(xx, 0, 1)
                    
#                     print(k, core_aln_min, core_aln[xx>10], core_aln_max)
#                     # TODO: For paired
#                     my_range = useful_range(core_aln)
#                     tseq = tseq.loc[my_range[0]:my_range[1]]
#                     with open(f"{k}.aln.fasta", "w") as fout:
#                         for c in tseq.columns:
#                             if c in ["gap","alnpos"]:
#                                 continue
#                             fout.write(f">{c}\n{''.join(tseq[c].values)}\n")




                # trace = defaultdict(list)
                # for c in channels:
                #     trace[c] = record.annotations["abif_raw"][c]
                # nucs = {"pos":[], "nuc":[], "pos_color":[]}
                # for i, c in enumerate(zip(record.annotations["abif_raw"]["PBAS1"],record.annotations["abif_raw"]["PLOC1"])):
                #     nucs["pos"].append(i)
                #     nucs["nuc"].append(chr(c[0]).upper())
                #     nucs["pos_color"].append(c[1])
                # nucs = pd.DataFrame(nucs)
                # nucs = nucs.iloc[100:-100, ]
                # # TODO: Find a way to report multiple of 3 nucleotide gaps
                # print(nucs, "Not working")
                # if not bc:
                #     pass
                # else:
                #     nucs["alt"] =  nucs["nuc"].values
                #     nucs_ambi = nucs[~nucs["nuc"].isin(["A", "C", "G", "T"])]
                #     for i, row in nucs_ambi.iterrows():
                #         # max
                #         peak_values = [trace["DATA9"][row["pos_color"]],
                #                     trace["DATA10"][row["pos_color"]],
                #                     trace["DATA11"][row["pos_color"]],
                #                     trace["DATA12"][row["pos_color"]]]
                #         if bc=="ref":
                #             sg_seq = ''.join(nucs['nuc'].values)
                #             with open("tmp/sanger_seq.fasta", "w") as fout:
                #                 fout.write(f">tseq\n{sg_seq}")
                #             # TODO: Replace with subprocess 
                #             command = ["blat", "-noDead", ref, {sanger_seq}, {sanger_seq.psl}]
                #             cmd(command)
                            
                #             for rec in SeqIO.parse("ref.fasta", "fasta"):
                #                 ref_seq = str(rec.seq).upper()
                #             mapping = pd.read_table("sanger_seq.psl", header=None).sort_values(0, ascending=False).head(1)
                #             if not mapping:
                #                 exit("Given sequence didn't match with the reference. Please check your input sequence. Exiting . . . . .")
                                
                #             for _, row in mapping.iterrows():
                #                 if row[0] < row[10]*0.9:
                #                     print("sequence match is less that 90%. Please check your sequence")
                #                 if row[17] == 1:
                #                     # TODO: Might need to check end and start of the mapping in case it doesn't match at the ends
                                    
                #                     pass
                #                 else:
                #                     pass
                #                 continue
                #             # mapping = 
                            
                #         elif bc=="max":
                #             max_index = peak_values.index(max(peak_values))
                #             base = "GATC"[max_index]
                #             nucs.loc[i,"alt"] = base
                #     # print(nucs)
                #             # print(base)
                #         elif bc=="neigh":
                #             if i == 0:
                #                 continue
                #             pre_nuc = nucs.loc[i-1, "nuc"]
                #             if pre_nuc in amb_bases[row["nuc"]]:
                #                 nucs.loc[i,"alt"] = pre_nuc
                #             else:
                #                 max_index = peak_values.index(max(peak_values))
                #                 base = "GATC"[max_index]
                #                 nucs.loc[i,"alt"] = base
    #             sequences[path.split(fl)[1].split(".")[0]] = ''.join(nucs['nuc'].values)
    #             # TODO: Return merged sequence file
    #         except:
    #             print("The file is neither ab1 or Fasta. Ignoring . . . ..")
    # with open(outfile, "w") as fout:
    #     for sec in sequences:
    #         fout.write(f">{sec}\n{sequences[sec]}\n")
        




# def sanger_seq(inff, outfile, is_fold=False):
#     """Generate readable nucleotide in fasta format from sanger ab1 files."""
#     print(inff, outfile, is_fold)
#     # Impute with alternative or max or ref
#     if is_fold:
#         files = glob(f"{inff}/*.ab1")
#         if not files:
#             exit(f"Given ab1 folder {inff} doen't contain any ab1 file. Exiting . . . . ")
#         else:
#             ab2fasta(files, outfile, "ref")
#     else:
#         command = ["cp", inff, outfile]
#     cmd(command)



# def cmd_check():
#     # TODO: "Check the listed commands which will be used by this program"
#     pass



@click.command()
@click.option("-sa_ab1", help="ab1 folder or sanger sequence file", type=str, default="ab1", show_default=True) # Convert this to folder
# @click.option("-fa", help="Fasta output file. If not given, only sequences will be printed in terminal", 
#               type=str, default=None, show_default=True)

@click.option("-asf", help="Assemblies folder containing fasta files", 
              type=str, default="assemblies", show_default=True)
# @click.option("-rf", help="Reference fasta", 
#               type=bool, default=False, show_default=True)
@click.option("-bs", help="Base replacement methos", 
              type=click.Choice([None, "max", "neigh", "ref"]), default=None, show_default=True)
@click.option("-outd", help="Output Folder", 
              type=str, default="Results", show_default=True)

def run(sa_ab1, asf, bs, outd):# , fa, asb, al, bs
    """
    Convert ab1 to Fasta based given parameters. must contain original sequence name., Must me N trimmed at the end
    work when either ab1 or sa parameter given. Must have ab1 files with ab1 suffix with sequence name as prefix\n
    .\n
    ├── ab1\n
    │   ├── x.ab1\n
    │   ├── y.ab1\n
    │   └── z.ab1\n
    ├── assembly.fasta\n
    ├── ref.fasta\n
    └── sanger.fasta\n
    Will generate `merged_output.fasta` in current folder.
    
    max: max peak
    neigh: earlier peak if ambiguous nucleotide contain earlier base
    ref: Use reference base in case of ambiguity
    """

    tmp_fold = "tmp"
    if not path.exists(tmp_fold):
        makedirs(tmp_fold)
    elif path.isfile(tmp_fold):
        exit("temperory folder name alaready exits in file format. exiting . . . . .")

        
    if not asf:
        exit("Assembly folder not given. Exiting . . . . . . .")
    elif not path.exists(asf) or not path.isdir(asf):
        exit("Given assembly folder doesn't exist or path is not a folder. Exiting . . . . .")
    else:
        assemblies = glob(f"{asf}/*.fasta")
    
    # TODO: Use below sequence to generate initial assembly file. Later it can be replace with other genes
    S_gene_seq = """ATGATATGATTTTATCTCTTCTTAGTAAAGGTAGACTTATAATTAGAGAAAACAACAGAGTTGTTATTTCTAGTG
ATGTTCTTGTTAACAACTAAACGAACAATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTT
AATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAA
GTTTTCAGATCCTCAGTTTTACATTCAACTCAGGACTTGTTCTTACCTTTCTTTTCCAATGTTACTTGGTTCCAT
GCTATACATGTCTCTGGGACCAATGGTACTAAGAGGTTTGATAACCCTGTCCTACCATTTAATGATGGTGTTTAT
TTTGCTTCCACTGAGAAGTCTAACATAATAAGAGGCTGGATTTTTGGTACTACTTTAGATTCGAAGACCCAGTCC
CTACTTATTGTTAATAACGCTACTAATGTTGTTATTAAAGTCTGTGAATTTCAATTTTGTAATGATCCATTTTTG
GGTGTTTATTACCACAAAAACAACAAAAGTTGGATGGAAAGTGAGTTCAGAGTTTATTCTAGTGCGAATAATTGC
ACTTTTGAATATGTCTCTCAGCCTTTTCTTATGGACCTTGAAGGAAAACAGGGTAATTTCAAAAATCTTAGGGAA
TTTGTGTTTAAGAATATTGATGGTTATTTTAAAATATATTCTAAGCACACGCCTATTAATTTAGTGCGTGATCTC
CCTCAGGGTTTTTCGGCTTTAGAACCATTGGTAGATTTGCCAATAGGTATTAACATCACTAGGTTTCAAACTTTA
CTTGCTTTACATAGAAGTTATTTGACTCCTGGTGATTCTTCTTCAGGTTGGACAGCTGGTGCTGCAGCTTATTAT
GTGGGTTATCTTCAACCTAGGACTTTTCTATTAAAATATAATGAAAATGGAACCATTACAGATGCTGTAGACTGT
GCACTTGACCCTCTCTCAGAAACAAAGTGTACGTTGAAATCCTTCACTGTAGAAAAAGGAATCTATCAAACTTCT
AACTTTAGAGTCCAACCAACAGAATCTATTGTTAGATTTCCTAATATTACAAACTTGTGCCCTTTTGGTGAAGTT
TTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGATTATTCT
GTCCTATATAATTCCGCATCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTCTGC
TTTACTAATGTCTATGCAGATTCATTTGTAATTAGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGA
AAGATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAATCTT
GATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGA
GATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCT
TTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTTGAA
CTTCTACATGCACCAGCAACTGTTTGTGGACCTAAAAAGTCTACTAATTTGGTTAAAAACAAATGTGTCAATTTC
AACTTCAATGGTTTAACAGGCACAGGTGTTCTTACTGAGTCTAACAAAAAGTTTCTGCCTTTCCAACAATTTGGC
AGAGACATTGCTGACACTACTGATGCTGTCCGTGATCCACAGACACTTGAGATTCTTGACATTACACCATGTTCT
TTTGGTGGTGTCAGTGTTATAACACCAGGAACAAATACTTCTAACCAGGTTGCTGTTCTTTATCAGGATGTTAAC
TGCACAGAAGTCCCTGTTGCTATTCATGCAGATCAACTTACTCCTACTTGGCGTGTTTATTCTACAGGTTCTAAT
GTTTTTCAAACACGTGCAGGCTGTTTAATAGGGGCTGAACATGTCAACAACTCATATGAGTGTGACATACCCATT
GGTGCAGGTATATGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAATCC
ATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTACTCTAATAACTCTATTGCCATACCCACA
AATTTTACTATTAGTGTTACCACAGAAATTCTACCAGTGTCTATGACCAAGACATCAGTAGATTGTACAATGTAC
ATTTGTGGTGATTCAACTGAATGCAGCAATCTTTTGTTGCAATATGGCAGTTTTTGTACACAATTAAACCGTGCT
TTAACTGGAATAGCTGTTGAACAAGACAAAAACACCCAAGAAGTTTTTGCACAAGTCAAACAAATTTACAAAACA
CCACCAATTAAAGATTTTGGTGGTTTTAATTTTTCACAAATATTACCAGATCCATCAAAACCAAGCAAGAGGTCA
TTTATTGAAGATCTACTTTTCAACAAAGTGACACTTGCAGATGCTGGCTTCATCAAACAATATGGTGATTGCCTT
GGTGATATTGCTGCTAGAGACCTCATTTGTGCACAAAAGTTTAACGGCCTTACTGTTTTGCCACCTTTGCTCACA
GATGAAATGATTGCTCAATACACTTCTGCACTGTTAGCGGGTACAATCACTTCTGGTTGGACCTTTGGTGCAGGT
GCTGCATTACAAATACCATTTGCTATGCAAATGGCTTATAGGTTTAATGGTATTGGAGTTACACAGAATGTTCTC
TATGAGAACCAAAAATTGATTGCCAACCAATTTAATAGTGCTATTGGCAAAATTCAAGACTCACTTTCTTCCACA
GCAAGTGCACTTGGAAAACTTCAAGATGTGGTCAACCAAAATGCACAAGCTTTAAACACGCTTGTTAAACAACTT
AGCTCCAATTTTGGTGCAATTTCAAGTGTTTTAAATGATATCCTTTCACGTCTTGACAAAGTTGAGGCTGAAGTG
CAAATTGATAGGTTGATCACAGGCAGACTTCAAAGTTTGCAGACATATGTGACTCAACAATTAATTAGAGCTGCA
GAAATCAGAGCTTCTGCTAATCTTGCTGCTACTAAAATGTCAGAGTGTGTACTTGGACAATCAAAAAGAGTTGAT
TTTTGTGGAAAGGGCTATCATCTTATGTCCTTCCCTCAGTCAGCACCTCATGGTGTAGTCTTCTTGCATGTGACT
TATGTCCCTGCACAAGAAAAGAACTTCACAACTGCTCCTGCCATTTGTCATGATGGAAAAGCACACTTTCCTCGT
GAAGGTGTCTTTGTTTCAAATGGCACACACTGGTTTGTAACACAAAGGAATTTTTATGAACCACAAATCATTACT
ACAGACAACACATTTGTGTCTGGTAACTGTGATGTTGTAATAGGAATTGTCAACAACACAGTTTATGATCCTTTG
CAACCTGAATTAGACTCATTCAAGGAGGAGTTAGATAAATATTTTAAGAATCATACATCACCAGATGTTGATTTA
GGTGACATCTCTGGCATTAATGCTTCAGTTGTAAACATTCAAAAAGAAATTGACCGCCTCAATGAGGTTGCCAAG
AATTTAAATGAATCTCTCATCGATCTCCAAGAACTTGGAAAGTATGAGCAGTATATAAAATGGCCATGGTACATT
TGGCTAGGTTTTATAGCTGGCTTGATTGCCATAGTAATGGTGACAATTATGCTTTGCTGTATGACCAGTTGCTGT
AGTTGTCTCAAGGGCTGTTGTTCTTGTGGATCCTGCTGCAAATTTGATGAAGACGACTCTGAGCCAGTGCTCAAA
GGAGTCAAATTACATTACACATAA"""
    
    with open(f"{tmp_fold}/CoVid_S_Gene.fasta", "w") as fout:
        fout.write(f">ref\n{S_gene_seq}\n")
    sanger_files = glob(f"{sa_ab1}/*") #  Either fasta or ab1 files
    # NOTE: Consider fasta file doen't ends with fasta amd ab1 file with ab1
    file_groups = {}
    for fl in sanger_files:
        flx = path.split(fl)[1]
        if flx.endswith(".fasta"):
            n_seq = 0
            for _ in SeqIO.parse(flx, "fasta"):
                n_seq += 1
            if n_seq != 1:
                exit(f"{flx} has {n_seq} sequences. Expecting 1. Exiting . . . . . .")
            flb = flx.split(".fasta")[0]
        elif flx.endswith(".ab1"):
            flb = flx.split("_")[0] # Fist part should be the name
        else:
            exit("File doesn't have fasta or ab1 extendion")
        if flb not in file_groups:
            file_groups[flb] = []
        file_groups[flb].append(fl)
    for k in file_groups:
        if len(file_groups[k]) > 2:
            exit(f"{len(file_groups[k])} files are associated with id: {k}. Expecting 1 or 2 only")
    # print(file_groups)
    for k in file_groups:
        ft = []
        for fl in file_groups[k]:
            ft.append(fl.split(".")[-1])
        if len(set(ft)) != 1:
            exit("Multiple fileformat for {k}. Exiting . . . . . . ")
    sanger_outputs = list(file_groups.values())

    @merge(sanger_outputs, f"{tmp_fold}/sanger.fasta") # Apply merge it. Ignore  if it fails to include
    def sanger_seq_r(inputfiles, outputfile):
        ab2fasta(inputfiles, outputfile)
    
    
    
    
    @transform(assemblies, formatter(".+/(?P<filebase>\w+)"),add_inputs(sanger_seq_r), "%s/{filebase[0]}.fasta" % outd)
    def alignments(infiles, outfile):
        # print(infiles,outfiles)
        sang_file = infiles[1]
        assembly = infiles[0]
        sanger_seq = {}
        sanger_seq_desc = {}
        for rec in SeqIO.parse(sang_file, "fasta"):
            sanger_seq_desc[rec.id] = rec.description
            sanger_seq[rec.id] = str(rec.seq)
            
        org_seq = {}
        for rec in SeqIO.parse(assembly, "fasta"):
            org_seq[rec.id]=str(rec.seq)

        flb = path.split(assembly)[1].split(".fasta")[0]
        command = ["blat", "-noHead", assembly, sang_file, f"{tmp_fold}/{flb}.psl"]
        # cmd(command)
        df = pd.read_table(f"{tmp_fold}/{flb}.psl", header=None)
        # print(df[df[9]==df[13]].sort_values(0, ascending=False))
        df = df[df[9]==df[13]].sort_values(0, ascending=False).drop_duplicates(9)
        print(df)
        for _, row in df.iterrows():
            # print(sanger_seq_desc[row[9]].split())
            start, end = list(map(int, sanger_seq_desc[row[9]].split()[1:]))
            if row[11] == 0:
                
                org_seq[row[9]] = org_seq[row[9]][:row[15]+start]+sanger_seq[row[9]][start:end+1]+org_seq[row[9]][row[15]+end+1:]
                
            elif row[12] == row[10]:
                org_seq[row[9]] = org_seq[row[9]][:row[16]-(row[10]-start)]+sanger_seq[row[9]][start:end+1]+org_seq[row[9]][row[16]-(row[10]-end-1):]
            else:
                exit("Contact Anmol")# TODO: Update
            pass
        with open(outfile, "w") as fout:
            for k in org_seq:
                print(k, len(org_seq[k]))
                fout.write(f">{k}\n{org_seq[k]}\n")
            
        #     df
        # print(infile, outfile)
        """Generates alignments if  sequences are not aligned"""
        # try:
        #     for _ in AlignIO.parse(infile, "fasta"):
        #         continue
        #     command = ["cp", infile, outfile]

        #     # TODO: Copy the file in desire folder
        # except:
        #     try:
        #         _ = SeqIO.parse(infile, "fasta")
        #         command = ["muscle", "-in", infile, "-out", outfile]
        #         # TODO: Align the sequence in desired folder
        #     except:
        #         print("Warning: Given file is not fasta")
        #         command = ["touch", outfile]
        # cmd(command)

    # @mkdir(outd)
    # @transform(alignments, formatter(".+/(?P<filebase>\w+)"), add_inputs(sanger_seq_r, tmp_fold), "%s/{filebase[0]}.fasta" % outd)       
    # def integration(infiles, outfile):
    #     print(infiles, outfile, "hello")
    #     integrate_seq(infiles, outfile)
        # pass

    pipeline_run(verbose=0)
    pipeline_printout_graph("flowchart.svg", "svg", no_key_legend=True )


if __name__ == '__main__':
    run()
