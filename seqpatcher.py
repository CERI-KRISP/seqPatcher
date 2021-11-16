#!/usr/bin/env python

# Author: Anmol Kiran
# Affiliation: University of Liverpool. UK
# and Malawi-Liverpool-Wellcome Trust, Malawi
# V0.0.1: 20/04/2021

from Bio import Seq, SeqIO
import pandas as pd
import numpy as np
import click

# from collections import ChainMap
from shutil import copyfile

# from functools import partial
import tempfile as tmf
from os import makedirs, path  # , access
from glob import glob

# from multiprocessing import cpu_count, Pool
from collections import defaultdict

# import sys
import subprocess as sb

# import os
import warnings

warnings.filterwarnings("ignore")
# warnings.simplefilter(action="ignore", category=FutureWarning)


__version__ = "0.0.1"

_s_gene_seq = """ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTT
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


def cmd(command):
    """Runs all the command using Popen communicate.
    :command: in form of list ie ['ls', '-l']
    :return: None
    """
    sb.Popen(command, stdout=sb.DEVNULL, stderr=sb.DEVNULL).communicate()


def min_max(values):
    """Returns minimum and maximum values in the given list/array.

    :values: list/array
    :returns: touple of min and max

    """
    return min(values), max(values)


def end_correction(aln_df):
    """In case of paired sequences Nucleotide will be replaced with apropriate nucleotide for simple processing

    :aln_df: Unprocessed alignment dataframe
    :returns: Processed data frame
    :TODO: Write what is being done here

    """

    pass


#           2       20


def trim(mmcount, length, aln_df):
    """Trim the sequence at the end based on mismatched in given length of

    :mmcount: Mismatch count in given range
    :length: Length at both end to explore
    :returns: trimmed alignment dataframe

    """
    # aln_len = length(aln)
    for col in aln_df.columns:
        if col == "ref":
            continue
        min_, max_ = min_max(aln_df[aln_df[col] != "-"].index)
        mismach_locations = aln_df[aln_df[col].isin("-NRYKMSWBDHV")].index
        mismach_locations = mismach_locations[
            (mismach_locations >= min_) & (mismach_locations <= max_)
        ]
        start_mismatch_location = mismach_locations[mismach_locations < length]
        if len(start_mismatch_location) >= mmcount:
            min_ = start_mismatch_location[-1] + 1
        end_mismatch_locations = mismach_locations[
            mismach_locations > (length(aln_df) - length)
        ]
        if len(end_mismatch_locations) >= mmcount:
            max_ = end_mismatch_locations[0] - 1
        aln_df.loc[:min_, col] = "-"
        aln_df.loc[max_:, col] = "-"

    return aln_df


def codon_aln(aln_df):
    """Correct alignment around codon

    :aln_df: normal alignment dataframe
    :returns: colodon alinment dataframe

    """
    df_shape = aln_df.shape[1]
    non_ref_seq = [seq_id for seq_id in df_shape.columns if seq_id != "ref"]
    if df_shape == 2:
        aln_min, aln_max = min_max(aln_df[aln_df[non_ref_seq[0]] != "-"].index)
        # TODO: Check whether is a 0 position if not change the locatio in no indel at next two positions
        #

        pass
    elif df_shape == 3:
        # TODO: Accept nucleotide at overhang and then try to change the codon alignmentt

        pass
    else:
        pass

    return


# TODO: How many codon should be allow to be remove without replacing with N
# TODO: Remove text from last mismatch - I don't think it will make any difference -  But allow to extend the ends for better mapping


def ranges(lst, given_gap=0):
    """
    A generator returns list of range based on given numbers
    [1,2,3,4,5, 10, 11, 12, 13, 14] => [[1,5], [10,14]]
    """
    lst_sorted = sorted(lst)
    init = 0
    for num in range(1, len(lst_sorted)):
        # reported_gap = lst_sorted[num] - lst_sorted[num - 1] - 1
        if lst_sorted[num] > given_gap + 1 + lst_sorted[num - 1]:
            # or (reported_gap % 3 != 0):
            # gap=0 means overlapping
            yield (lst_sorted[init], lst_sorted[num - 1])
            init = num
    yield (lst_sorted[init], lst_sorted[-1])


def useful_range(lst, gap=10):  # TODO: Add in command
    """Returns first occurance maximum length block range."""
    # TODO: Check with other that if they are happy with it
    # or they want other fragment to be selected

    trange = list(ranges(lst, gap))
    pre_frag_len = 0
    myrange = trange[0]
    for rng in trange:
        frag_len = rng[1] - rng[0] + 1
        if pre_frag_len < frag_len:
            pre_frag_len = frag_len
            myrange = rng
    return myrange


def df_reverse_complement(dataframe):
    """Reverse Complement of the nucleotide dataframe."""
    dataframe = dataframe.loc[::-1]

    def rvc(nuc):
        """Return complement base."""
        return str(Seq.Seq(nuc).complement())

    def rvc_dict(nuc_dict):
        """Returns complement dictionary."""
        temp_nuc_dict = {}
        for nuc in nuc_dict:
            temp_nuc_dict[rvc(nuc)] = nuc_dict[nuc]
        return temp_nuc_dict

    dataframe["nuc"] = dataframe["nuc"].apply(rvc)
    dataframe["peak"] = dataframe["peak"].apply(rvc_dict)
    return dataframe


def orient(seqfile, ref, tmp_fold):
    """Returns orientation of the sequence"""
    if seqfile.endswith(".ab1"):
        flb = path.split(seqfile)[1].rsplit(".", 1)[0]
        record = SeqIO.read(seqfile, "abi")
        seq = "".join(
            [chr(ascii_val)
             for ascii_val in record.annotations["abif_raw"]["PBAS1"]]
        )
        with open(f"{tmp_fold}/{flb}.fasta", "w") as fout:
            fout.write(f">tmp\n{seq}\n")
        seqfile = f"{tmp_fold}/{flb}.fasta"

    flb = path.split(seqfile)[1].rsplit(".", 1)[0]
    command = [
        "blat",
        "-noHead",
        ref,
        seqfile,
        f"{tmp_fold}/{flb}.psl",
    ]  # Mapping against reference
    cmd(command)
    blat_df = (
        pd.read_table(f"{tmp_fold}/{flb}.psl", header=None)
        .sort_values(0, ascending=False)
        .drop_duplicates(9)
    )
    if len(blat_df) != 1:
        print(f"No match found of {seqfile}. Ignoring")
        return None
    for _, row in blat_df.iterrows():
        if row[8] == "-":
            return "R"
        else:
            return "F"


def ab1seq(infile, tmp_fold, peak_selection=None):
    """ab1 to seq trimmed based on reference."""

    amb_bases = {  # All ambiguous nucleotides
        "N": set(["A", "C", "G", "T"]),
        "R": set(["A", "G"]),
        "Y": set(["T", "C"]),
        "K": set(["G", "T"]),
        "M": set(["A", "C"]),
        "S": set(["G", "C"]),
        "W": set(["A", "T"]),
        "B": set(["C", "G", "T"]),
        "D": set(["A", "G", "T"]),
        "H": set(["A", "C", "T"]),
        "V": set(["A", "C", "G"]),
    }

    bases = {"DATA9": "G", "DATA10": "A", "DATA11": "T", "DATA12": "C"}

    record = SeqIO.read(infile, "abi")
    trace = defaultdict(list)

    for channel in bases:
        trace[bases[channel]] = record.annotations["abif_raw"][channel]

    nuc_df = {"nuc": [], "peak": []}

    for channel in zip(
        record.annotations["abif_raw"]["PBAS1"],
        record.annotations["abif_raw"]["PLOC1"],
    ):
        ambi_base = chr(channel[0])
        if ambi_base in amb_bases:
            # if not peak_selection:
            nuc_df["nuc"].append(ambi_base)
            td = {}  # TODO: Please check what does td reprensts
            for base in amb_bases[ambi_base]:
                td[base] = trace[base][channel[1]]
            nuc_df["peak"].append(td)
            # elif peak_selection == "max":
            # tb = ""
            # td = 0
            # for base in amb_bases[ambi_base]:
            # if td < trace[base][channel[1]]:
            # td = trace[base][channel[1]]
            # tb = base
            # nuc_df["nuc"].append(tb)
            # nuc_df["peak"].append({tb: td})
        else:
            nuc_df["nuc"].append(ambi_base)
            nuc_df["peak"].append({ambi_base: trace[ambi_base][channel[1]]})
    nuc_df = pd.DataFrame(nuc_df)
    # In case of ambigious nucleotides
    peak = nuc_df["peak"].apply(lambda x: np.mean(list(x.values()))).values

    correct_peak_pos = np.where(peak > 100)[0]
    # TODO: Autothreshold
    min_pos, max_pos = np.min(correct_peak_pos), np.max(
        correct_peak_pos
    )  # 1 is added as last position is not included
    nuc_df = nuc_df.loc[min_pos:max_pos]

    if infile.endswith(".R.ab1"):
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


def merge_base_peak(nuc_df, peak_dict):
    """Merge the peak related information in nucleotide dataframe."""
    # NOTE: Useful only when there is ambigiuity in the base call, rest remove to save memory

    nuc_df["idx"] = list(nuc_df.index)

    # df = enumerate_columns(df)
    # print(df)
    to_drop = ["idx"]
    for nuc in peak_dict:
        # Expecting sequence will not change while aligning, except insertion of gaps
        peak_dict[nuc][f"{nuc}_idx"] = list(
            nuc_df[nuc_df[f"{nuc}"] != "-"].index
        )  # list(range(len(peak_dict[k]))
        nuc_df = nuc_df.merge(
            peak_dict[nuc], left_on="idx", right_on=f"{nuc}_idx", how="outer"
        )
        to_drop += [f"{nuc}_idx", f"{nuc}_nuc"]

    nuc_df = nuc_df.drop(to_drop, axis=1)
    return nuc_df


def rep_paired_base(lst):
    """Selecting representative base in presence of forward and
    reverse sanger sequencing data."""
    if lst["F"] == "-":
        if lst["R"] == "-":
            return "-"
        else:
            bmax = 0
            base = "A"
            peaks = lst["R_peak"]
            for bs in peaks:
                if peaks[bs] > bmax:
                    bmax = peaks[bs]
                    base = bs
            return base
    else:
        if lst["R"] == "-":
            bmax = 0
            base = "A"
            peaks = lst["F_peak"]
            for bs in peaks:
                if peaks[bs] > bmax:
                    bmax = peaks[bs]
                    base = bs
            return base
        else:
            # remove lst['ref'], in case that is present
            common = list((set(lst["R_peak"]) & set(
                lst["F_peak"])) - set([lst["ref"]]))
            if len(common) == 1:
                return common[0]
            else:
                return lst["ref"]


def aln_clean(aln_df, gap=10):
    """Clearning alihnment dataframe to remove unneccessary gaps and
    alignements."""
    sang_type = None
    if "F" in aln_df and "R" in aln_df:
        idx = aln_df[~((aln_df["F"] == "-") & (aln_df["R"] == "-"))].index
        sang_type = "P"  # Fore paired
    else:
        if "F" in aln_df:
            idx = aln_df[aln_df["F"] != "-"].index
            sang_type = "F"
        else:
            idx = aln_df[aln_df["R"] != "-"].index
            sang_type = "R"

    u_range = list(useful_range(list(idx), gap))
    for col in aln_df.columns:

        if col == "ref":
            continue
        aln_df.loc[: u_range[0] - 1, col] = "-"
        aln_df.loc[u_range[1] + 1:, col] = "-"

    if sang_type in "FR":
        aln_df["consensus"] = aln_df[sang_type].values
        ambi_indexes = aln_df[aln_df[sang_type].isin("NRYKMSWBDHV")].index
        for ambi_index in ambi_indexes:
            bmax = 0
            base = "A"
            peaks = aln_df.loc[ambi_index, f"{sang_type}_peak"]
            for bs in peaks:
                if peaks[bs] > bmax:
                    bmax = peaks[bs]
                    base = bs
            aln_df.loc[ambi_index, "consensus"] = base

        insert_ranges = aln_df.loc[u_range[0]: u_range[1]]
        insert_ranges = insert_ranges[insert_ranges["ref"] == "-"].index
        if insert_ranges.any():
            insert_ranges = ranges(insert_ranges)
            for insert_range in insert_ranges:
                if (insert_range[1] - insert_range[0] + 1) % 3 == 0:
                    continue
                else:
                    if (insert_range[1] - insert_range[0] + 1) < 3:
                        aln_df.loc[insert_range[0]
                            : insert_range[1], "consensus"] = "-"
                    else:
                        # TODO: Talk to san one more time
                        aln_df.loc[insert_range[0]
                            : insert_range[1], "consensus"] = "N"
        del_ranges = aln_df.loc[u_range[0]: u_range[1]]
        del_ranges = del_ranges[del_ranges[sang_type] == "-"].index
        if del_ranges.any():
            del_ranges = ranges(del_ranges)
            for del_range in del_ranges:
                if (del_range[1] - del_range[0] + 1) % 3 == 0:
                    continue
                else:
                    if (del_range[1] - del_range[0] + 1) < 3:
                        aln_df.loc[
                            del_range[0]: del_range[1], "consensus"
                        ] = aln_df.loc[del_range[0]: del_range[1], "ref"].values
                    else:
                        aln_df.loc[del_range[0]
                            : del_range[1], "consensus"] = "N"

    else:
        aln_df["consensus"] = aln_df.apply(rep_paired_base, axis=1)
        insert_ranges = aln_df.loc[u_range[0]: u_range[1]]
        insert_ranges = insert_ranges[
            (insert_ranges["ref"] == "-")
            & ((insert_ranges["F"] != "-") & (insert_ranges["R"] != "-"))
        ].index
        if insert_ranges.any():
            insert_ranges = ranges(insert_ranges)
            for insert_range in insert_ranges:
                if (insert_range[1] - insert_range[0] + 1) % 3 == 0:
                    continue
                else:
                    if (insert_range[1] - insert_range[0] + 1) < 3:
                        aln_df.loc[insert_range[0]
                            : insert_range[1], "consensus"] = "-"
                    else:
                        # TODO: Talk to san one more time
                        aln_df.loc[insert_range[0]
                            : insert_range[1], "consensus"] = "N"

        del_ranges = aln_df.loc[u_range[0]: u_range[1]]
        del_ranges = del_ranges[
            (del_ranges["ref"] != "-")
            & ((del_ranges["F"] == "-") & (del_ranges["R"] == "-"))
        ].index
        if del_ranges.any():
            del_ranges = ranges(del_ranges)
            for del_range in del_ranges:
                if (del_range[1] - del_range[0] + 1) % 3 == 0:
                    continue
                else:
                    if (del_range[1] - del_range[0] + 1) < 3:
                        aln_df.loc[
                            del_range[0]: del_range[1], "consensus"
                        ] = aln_df.loc[del_range[0]: del_range[1], "ref"].values
                    else:
                        aln_df.loc[del_range[0]
                            : del_range[1], "consensus"] = "N"

    return aln_df, u_range


def fasta_map2ref(infile, gap, tmp_fold, n3, idb):
    # TODO: Use some part for arranging the sequence
    """Integrates Sanger fasta to refgene

    :args: infile, outfile, tmp_fold, idb, gap


    """
    sequences = {}
    for rec in SeqIO.parse(infile, "fasta"):
        if infile.endswith(".R.fasta"):  # Generates revese complement
            sequences[rec.id] = rec.seq.reverse_complement()
        else:
            sequences[rec.id] = rec.seq

    for rec in SeqIO.parse(f"{tmp_fold}/ref.fasta", "fasta"):
        sequences["ref"] = rec.seq

    flb = path.split(infile)[1].split(".")[0]

    aln_df = aln_df_with_ref(sequences, flb, tmp_fold)
    # print(aln_df)
    mapped_index = aln_df[aln_df[flb] != "-"].index
    u_range = useful_range(mapped_index, gap)
    aln_df["concensus"] = "-"
    aln_df.loc[: u_range[0], "concensus"] = aln_df.loc[: u_range[0], "ref"]
    aln_df.loc[u_range[0]: u_range[1], "concensus"] = aln_df.loc[
        u_range[0]: u_range[1], flb
    ]
    aln_df.loc[u_range[1]:, "concensus"] = aln_df.loc[u_range[1]:, "ref"]

    if idb in ["del", "both"]:
        del_sites = aln_df[aln_df[flb] == "-"].index.values
        if len(del_sites):
            del_ranges = ranges(del_sites)
            del_ranges = [
                rng for rng in del_ranges if (rng[0] > u_range[0] & rng[1] < u_range[1])
            ]
            # Del codon is accepted when codon deletion is allowed else deletions are filled with gaps
            for rng in del_ranges:
                if n3:
                    if (rng[1] - rng[0] + 1) % 3 != 0:
                        aln_df.loc[rng[0]: rng[1], "concensus"] = aln_df.loc[
                            rng[0]: rng[1], "ref"
                        ]
                else:
                    aln_df.loc[rng[0]: rng[1], "concensus"] = aln_df.loc[
                        rng[0]: rng[1], "ref"
                    ]

    if idb in ["ins", "both"]:
        ins_sites = aln_df[aln_df["ref"] == "-"].index.values
        if len(ins_sites):
            ins_ranges = ranges(ins_sites)
            ins_ranges = [
                rng for rng in ins_ranges if (rng[0] > u_range[0] & rng[1] < u_range[1])
            ]
            for rng in ins_ranges:
                if n3:
                    if (rng[1] - rng[0] + 1) % 3 != 0:
                        aln_df.loc[rng[0]: rng[1], "concensus"] = aln_df.loc[
                            rng[0]: rng[1], "ref"
                        ]
                else:
                    aln_df.loc[rng[0]: rng[1], "concensus"] = aln_df.loc[
                        rng[0]: rng[1], "ref"
                    ]

    seq = "".join(aln_df["concensus"]).replace("-", "N")
    outfile = f"{tmp_fold}/sanger_converted_fasta/{flb}.fasta"

    with open(outfile, "w") as fout:
        fout.write(f">{flb} {u_range[0]} {u_range[1]}\n{seq}\n")


def ab1_2seq_map2ref(infiles, gap, tmp_fold):
    """TODO: Docstring for ab1_2seq.

    :infiles: TODO
    :returns: TODO

    """
    ab1seq_dfs = {}
    tsequences = {}
    flb = path.split(infiles[0])[1].split(".")[0]
    for fl in infiles:
        if fl.endswith(".F.ab1"):
            ab1seq_dfs["F"] = ab1seq(fl, tmp_fold)  # TODO: add the condition
            tsequences["F"] = "".join(ab1seq_dfs["F"]["nuc"].values)
        else:
            ab1seq_dfs["R"] = ab1seq(fl, tmp_fold)
            tsequences["R"] = "".join(ab1seq_dfs["R"]["nuc"].values)
    # TODO: Keep the ref name as ref
    for rec in SeqIO.parse(f"{tmp_fold}/ref.fasta", "fasta"):
        tsequences[rec.id] = str(rec.seq)

    for k in ab1seq_dfs:
        ab1seq_dfs[k].columns = [f"{k}_{col}" for col in ab1seq_dfs[k].columns]

    aln_with_peak = merge_base_peak(
        aln_df_with_ref(tsequences, flb, tmp_fold), ab1seq_dfs
    )
    # for n in ["F", "R"]:
    # cl = f"{n}_nuc"
    # if cl in aln_with_peak:
    # aln_with_peak.loc[pd.isnull(aln_with_peak[cl]), cl] = "-"
    # print(aln_with_peak)
    # exit(0)
    aln_with_peak, u_range = aln_clean(aln_with_peak)
    # aln_with_peak.to_csv("testxx.csv")
    # exit(0)
    seq = "".join(list(aln_with_peak["consensus"].values))
    # TODO: Generate sequence and exit
    outfile = path.split(infiles[0])[1].split(".")[0]
    outfile = f"{tmp_fold}/sanger_converted_fasta/{outfile}.fasta"

    output_file = open(outfile, "w")
    output_file.write(f">{flb} {u_range[0]} {u_range[1]}\n{seq}\n")
    output_file.close()


def ab2fasta(
    sang_dict, tmp_fold, gap, key, n3, idb  # , bc="neigh"
):  # Base criteria, max, neighbors, mixed # Inputfiles paired and none paired
    # sanger_outputs, tmp_fold, gap
    """Retains fasta and converts ab1 to fasta"""
    # print(key, sang_dict)
    infiles = sang_dict[key]

    if len(infiles) == 1 and infiles[0].endswith(".fasta"):
        fasta_map2ref(infiles[0], gap, tmp_fold, n3, idb)

    else:
        ab1_2seq_map2ref(infiles, gap, tmp_fold)


def files_and_groups(sanger_files):
    """List fasta and ab1 sequences as dictionary.

    :sanger_files: List of files in the folder
    :returns: dictionary of files as paired/single and fasta/ab1

    """
    file_groups = {}
    for file_ in sanger_files:
        flx = path.split(file_)[1]
        flb = flx.split(".")[0]
        if flb not in file_groups:
            file_groups[flb] = []
        file_groups[flb].append(file_)
    return file_groups


def non_overlapping_ids(asseblies, ab1s):
    """Check for ovelapping and non-ovelapping ids and generates csv table

    :asseblies: Fasta assembly containing folder
    :ab1s: Sanger generated ab1 or fasta files
    :returns: Pandas dataframe

    """
    # Assembly IDS
    assembly_ids = []
    for fl in glob(f"{asseblies}/*.fasta"):
        for rec in SeqIO.parse(fl, "fasta"):
            assembly_ids.append(rec.id)

    # Sanger sequences IDs
    sanger_fasta = []
    for fl in glob(f"{ab1s}/*.fasta"):
        for rec in SeqIO.parse(fl, "fasta"):
            sanger_fasta.append(rec.id)
    sanger_fasta_missing_assembly = ",".join(
        set(sanger_fasta) - set(assembly_ids))

    # Sanger forward ab1 IDs
    sanger_ab1_f = []
    for fl in glob(f"{ab1s}/*.F.ab1"):
        sanger_ab1_f.append(path.split(fl)[1].split(".F.ab1")[0])
    sanger_ab1_f_missing_assembly = ",".join(
        set(sanger_fasta) - set(sanger_ab1_f))

    # Sanger Reverse ab1 IDs
    sanger_ab1_r = []
    for fl in glob(f"{ab1s}/*.R.ab1"):
        sanger_ab1_r.append(path.split(fl)[1].split(".R.ab1")[0])
    sanger_ab1_r_missing_assembly = ",".join(
        set(sanger_fasta) - set(sanger_ab1_r))

    data_frame = {"assembly": [], "ab1_Forward": [],
                  "ab1_Reverse": [], "fasta": []}
    for assembly_id in assembly_ids:
        data_frame["assembly"].append(assembly_id)

        if assembly_id in sanger_ab1_f:
            data_frame["ab1_Forward"].append(1)
        else:
            data_frame["ab1_Forward"].append(0)

        if assembly_id in sanger_ab1_r:
            data_frame["ab1_Reverse"].append(1)
        else:
            data_frame["ab1_Reverse"].append(0)

        if assembly_id in sanger_fasta:
            data_frame["fasta"].append(1)
        else:
            data_frame["fasta"].append(0)

    deduct = False

    if (
        sanger_ab1_f_missing_assembly
        or sanger_ab1_r_missing_assembly
        or sanger_fasta_missing_assembly
    ):
        deduct = True
        data_frame["assembly"].append("No Assembly")
        data_frame["ab1_Forward"].append(sanger_ab1_f_missing_assembly)
        data_frame["ab1_Reverse"].append(sanger_ab1_r_missing_assembly)
        data_frame["fasta"].append(sanger_fasta_missing_assembly)

    data_frame = pd.DataFrame(data_frame)

    # Check for overlap
    if deduct:
        is_overlap = (
            data_frame.iloc[:-1][["ab1_Forward",
                                  "ab1_Reverse", "fasta"]].sum().sum()
        )
    else:
        is_overlap = data_frame[["ab1_Forward",
                                 "ab1_Reverse", "fasta"]].sum().sum()

    if not is_overlap:
        return pd.DataFrame()
    return data_frame


def integrate_in_assembly(outputfold, tmp_fold, sample_id):
    """Mergre sange sequences in NGS assemblies
    :outputfold: Final output folder
    :tempfold: Intermediate files generated by other part of the script
    :sample_id: Sample with NGS assembly and sange sequencing

    """
    # TODO: create psl folder in temp folder

    assembly = f"{tmp_fold}/assemblies/{sample_id}.fasta"
    sanger = f"{tmp_fold}/sanger_converted_fasta/{sample_id}.fasta"
    psl_file = f"{tmp_fold}/tmp/{sample_id}.psl"
    command = ["blat", "-noHead", assembly, sanger, psl_file]
    cmd(command)

    sanger_seq = {}
    sanger_seq_desc = {}
    for rec in SeqIO.parse(sanger, "fasta"):
        sanger_seq_desc[rec.id] = rec.description
        sanger_seq[rec.id] = str(rec.seq)
    org_seq = {}
    for rec in SeqIO.parse(assembly, "fasta"):
        org_seq[rec.id] = str(rec.seq)
    blat_df = pd.read_table(psl_file, header=None)

    blat_df = (
        blat_df[blat_df[9] == blat_df[13]]
        .sort_values(0, ascending=False)
        .drop_duplicates(9)
    )

    for _, row in blat_df.iterrows():
        start, end = list(map(int, sanger_seq_desc[row[9]].split()[1:]))
        if row[11] == 0:
            org_seq[row[9]] = (
                org_seq[row[9]][: row[15] + start]
                + sanger_seq[row[9]][start: end + 1]
                + org_seq[row[9]][row[15] + end + 1:]
            )

        elif row[12] == row[10]:
            org_seq[row[9]] = (
                org_seq[row[9]][: row[16] - (row[10] - start)]
                + sanger_seq[row[9]][start: end + 1]
                + org_seq[row[9]][row[16] - (row[10] - end - 1):]
            )
        else:
            if start >= row[11]:
                org_seq[row[9]] = (
                    org_seq[row[9]][: row[15] + start]
                    + sanger_seq[row[9]][start: end + 1]
                    + org_seq[row[9]][row[15] + end + 1:]
                )
            elif end < row[10]:
                org_seq[row[9]] = (
                    org_seq[row[9]][: row[16] - (row[10] - start)]
                    + sanger_seq[row[9]][start: end + 1]
                    + org_seq[row[9]][row[16] - (row[10] - end - 1):]
                )
            else:

                print("Please report at a bug at")
                print("https://github.com/krisp-kwazulu-natal/" "seqPatcher/issues")
    with open(f"{outputfold}/{sample_id}.fasta", "w") as fout:
        for k in org_seq:
            fout.write(f">{k}\n{org_seq[k]}\n")


@click.command()
@click.option(
    "-s",
    "--sanger-ab1",
    "sa_ab1",
    help="ab1 folder or sanger sequence file",
    type=str,
    default="sanger_ab1",
    show_default=True,
)  # Convert this to folder
# "/home/devil/Documents/San/Corona/Merging/Sanger/12April2021"
# @click.option("-fa", help="Fasta output file.
# If not given, only sequences will be printed in terminal",
#               type=str, default=None, show_default=True)
@click.option(
    "-a",
    "-assemblies-foder",
    "asf",
    help="Assemblies folder containing fasta files",
    type=str,
    default="unitest/assemblies",
    show_default=True,
)
@click.option(
    "-o",
    "--out-dir",
    "outd",
    help="Output Folder",
    type=str,
    default="Results",
    show_default=True,
)
@click.option(
    "-t",
    "--tab",
    "tab",
    help="CSV file for overlapping assemblies and sanger ids." " If not given, stdout.",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "-O",
    "--output-fasta",
    "ss",
    help="Sanger output fasta from ab1",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "-R",
    "--ref-gene-fasta-file",
    "rf",
    help="Refence gene file in fasta format",
    type=str,
    default=None,
    show_default=True,
)
# @click.option(
# "-n",
# "--cpu",
# "cpu",
# help="Number of CPU to use",
# type=int,
# default=1,
# show_default=True,
# )
@click.option(
    "-c",
    "--clean-intermediate",
    "ci",
    help="Remove intermediate file",
    type=bool,
    default=True,
    show_default=True,
)
@click.option(
    "-g",
    "--gap-allowed",
    "gap",
    help="Gap Allowed between aligned fragment to consider the region continuous",
    type=int,
    default=10,
    show_default=True,
)
@click.option(
    "-3",
    "--only-3-nuc",
    "n3",
    help="Allow  3 nucleotide InDels else replace with reference nucleotides",
    type=bool,
    default=True,
    show_default=True,
)
@click.option(
    "-x",
    "--indel-selection",
    "idb",
    help="Replace Insertion, Deletion or Both",
    type=click.Choice(["del", "ins", "both"]),
    default="del",
    show_default=True,
    multiple=False,
)
@click.version_option(__version__)
def run(sa_ab1, asf, outd, tab, ss, rf, ci, gap, n3, idb):  # , fa, asb, al, bscpu,
    # print(sa_ab1, asf, outd, tab, ss, rf, cpu, ci, gap, n3, idb)
    """
    Reports nucleotide sequence from Sanger chromatogram data based on user
    provided parameters and integrate that in assembly generated using NGS
    data"""
    # TODO: Integrate multi core system
    # if cpu < 1:
    # print("Number of CPU use given is < 1.")
    # print("Using default CPU value of 1")
    # cpu = 1
    # elif cpu > cpu_count() - 1:
    # print("Given cpu usage is more or equal to cpus avalable on system")
    # print(f"Setting CPU usage to {cpu_count() - 1 }")
    # cpu = cpu_count() - 1
    # else:
    # pass

    # pool = Pool(cpu)

    if not sa_ab1:
        exit("Sanger data folder is not given. Exiting . . . .")
    if not path.exists(sa_ab1) or not path.isdir(sa_ab1):
        exit(
            f"Given sanger data folder {sa_ab1} doesn't exist or path is not a folder."
            " Exiting . . . . ."
        )

    if not asf:
        exit("Assembly folder not given. Exiting . . . . . . .")
    if not path.exists(asf) or not path.isdir(asf):
        exit(
            f"Given assembly folder {asf} doesn't exist or path is not a folder."
            " Exiting . . . . ."
        )

    if not rf:
        print(
            "Reference sequence file is not given."
            " Considering sars-cov-2 spike protein sequence as reference"
        )
    elif not path.exists(rf) or not path.isfile(rf):
        print(
            f"Given reference file {rf} doesn't exist or path is not a file.")

        print("Considering sars-cov-2 spike protein sequence as reference")
        rf = None

    # tmp_fold = "tmp"
    tmp_fold = tmf.mkdtemp()

    # ----------Housekeeping----------------
    # Creating temporary  files and folders

    ref_path = f"{tmp_fold}/ref.fasta"
    makedirs(outd, exist_ok=True)
    for folder in [
        "assemblies",
        "sanger_raw",
        "sanger_converted_fasta",
        "sanger_final_fasta",
        "tmp",
    ]:
        makedirs(f"{tmp_fold}/{folder}", exist_ok=True)

    # Copying ref fasta
    with open(ref_path, "w") as fout:
        if not rf:
            fout.write(f">ref\n{_s_gene_seq}\n")
        else:
            seq_count = 0
            for rec in SeqIO.parse(rf, "fasta"):
                seq_count += 1
                if seq_count > 1:
                    exit(
                        f"{rf} contains more than 1 sequence. "
                        "Expect only one."
                        " Exiting."
                    )
                seq_desc = rec.description.split()
                seq = rec.seq
                if len(seq_desc) == 3 and seq_desc[1] == "CDS":
                    seq = seq[int(seq_desc[2]):]

                fout.write(f">{rec.id}\n{seq}\n")

            if not seq_count:
                exit(
                    f"{rf} contains 0 (zero) sequence. " "Expect only one." " Exiting."
                )

    sanger_files = glob(f"{sa_ab1}/*")
    if not sanger_files:
        exit(f"No file found in {sa_ab1} folder. Exiting . . . .")
    sanger_names = []
    for fl in sanger_files:
        if fl.endswith(".fasta"):
            for rec in SeqIO.parse(fl, "fasta"):
                if rec.id not in sanger_names:
                    sanger_names.append(rec.id)
        if fl.endswith(".ab1"):
            flb = path.split(fl)[1].split(".")[0]
            if flb not in sanger_names:
                sanger_names.append(flb)
    if ss:
        for fl in glob(f"{sa_ab1}/*"):
            if fl.endswith(".fasta"):
                for rec in SeqIO.parse(fl, "fasta"):
                    with open(f"{tmp_fold}/tmp/{rec.id}.fasta", "w") as fout:
                        fout.write(f">{rec.id}\n{rec.seq}\n")

                    l_r = orient(
                        f"{tmp_fold}/tmp/{rec.id}.fasta",
                        ref_path,
                        f"{tmp_fold}/tmp",
                    )

                    with open(
                        f"{tmp_fold}/sanger_raw/{rec.id}.{l_r}.fasta", "w"
                    ) as fout:
                        fout.write(f">{rec.id}\n{rec.seq}\n")
            if fl.endswith(".ab1"):
                fl_e = path.split(fl)[1]
                flb = fl_e.split(".")[0]
                l_r = orient(fl, ref_path, f"{tmp_fold}/tmp")
                copyfile(fl, f"{tmp_fold}/sanger_raw/{flb}.{l_r}.ab1")

        sanger_outputs = files_and_groups(glob(f"{tmp_fold}/sanger_raw/*"))
        for k in sanger_outputs:
            ab2fasta(sanger_outputs, tmp_fold, gap, k, n3, idb)

        with open(ss, "w") as fout:
            for fl in glob(f"{tmp_fold}/sanger_converted_fasta/*"):
                for rec in SeqIO.parse(fl, "fasta"):
                    coors = rec.description.split()[1:]
                    coors = int(coors[0]), int(coors[1])
                    fout.write(f">{rec.id}\n{rec.seq[coors[0]:coors[1]]}\n")

    assembly_files = glob(f"{asf}/*.fasta")
    if not assembly_files:
        exit(f"No file found in {asf} folder. Exiting . . . .")

    assembly_names = []
    for fl in assembly_files:
        for rec in SeqIO.parse(fl, "fasta"):
            assembly_names.append(rec.id)
    # else:
    # exit(f"No file is assembly folder {asf}. Exiting . . . .")
    if not assembly_names:
        exit("No assembly sequence found. Exiting . . . . .")

    common_ids = set(assembly_names) & set(sanger_names)
    if not common_ids:
        exit(
            "Genome assembly and sanger sequencing data doesn't have common"
            " id(s). Exiting..."
        )

    # Copying assembly to tmp folder
    for fl in glob(f"{asf}/*.fasta"):
        for rec in SeqIO.parse(fl, "fasta"):
            if rec.id in common_ids:
                with open(f"{tmp_fold}/assemblies/{rec.id}.fasta", "w") as fout:
                    fout.write(f">{rec.id}\n{rec.seq}\n")

    # Copying sanger files
    if not ss:  # Else already done earlier
        for fl in glob(f"{sa_ab1}/*"):
            if fl.endswith(".fasta"):
                for rec in SeqIO.parse(fl, "fasta"):
                    if rec.id in common_ids:  # TODO: Do something smart
                        with open(f"{tmp_fold}/tmp/{rec.id}.fasta", "w") as fout:
                            fout.write(f">{rec.id}\n{rec.seq}\n")

                        l_r = orient(
                            f"{tmp_fold}/tmp/{rec.id}.fasta",
                            ref_path,
                            f"{tmp_fold}/tmp",
                        )

                        with open(
                            f"{tmp_fold}/sanger_raw/{rec.id}.{l_r}.fasta", "w"
                        ) as fout:
                            fout.write(f">{rec.id}\n{rec.seq}\n")
            if fl.endswith(".ab1"):
                fl_e = path.split(fl)[1]
                flb = fl_e.split(".")[0]
                if flb in common_ids:

                    l_r = orient(fl, ref_path, f"{tmp_fold}/tmp")
                    copyfile(fl, f"{tmp_fold}/sanger_raw/{flb}.{l_r}.ab1")

    seq_id_df = non_overlapping_ids(
        f"{tmp_fold}/assemblies", f"{tmp_fold}/sanger_raw")

    seq_id_df = seq_id_df[
        ~(
            (seq_id_df["assembly"] == "No Assembly")
            | (
                (seq_id_df["ab1_Forward"] == 0)
                & (seq_id_df["ab1_Reverse"] == 0)
                & (seq_id_df["fasta"] == 0)
            )
        )
    ]

    if tab:
        seq_id_df.to_csv(tab, index=False)
    else:
        print(seq_id_df.to_csv(index=False))

    seq_id_df = seq_id_df[
        (
            ((seq_id_df["ab1_Forward"] == 1) | (seq_id_df["ab1_Reverse"] == 1))
            & (seq_id_df["fasta"] == 0)
        )
        | (
            ((seq_id_df["ab1_Forward"] == 0) & (seq_id_df["ab1_Reverse"] == 0))
            & (seq_id_df["fasta"] == 1)
        )
    ]
    print("The patcher executed for ..")
    print(",".join(seq_id_df["assembly"]))

    assemblies = glob(f"{tmp_fold}/assemblies/*.fasta")

    sanger_outputs = files_and_groups(glob(f"{tmp_fold}/sanger_raw/*"))
    for k in sanger_outputs:
        ab2fasta(sanger_outputs, tmp_fold, gap, k, n3, idb)

    for id_ in sanger_outputs:
        integrate_in_assembly(outd, tmp_fold, id_)

    if ci:
        command = ["rm", "-rf", tmp_fold]
        cmd(command)


if __name__ == "__main__":
    run()
