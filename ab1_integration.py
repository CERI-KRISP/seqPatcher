#!/usr/bin/env python

# Author: Anmol Kiran
# Affiliation: University of Liverpool. UK
# and Malawi-Liverpool-Wellcome Trust, Malawi
# V0.1: 20/04/2021


import subprocess as sb
import sys
from collections import defaultdict
from glob import glob
from os import makedirs, path

import click
import numpy as np
import pandas as pd
from Bio import Seq, SeqIO
from ruffus import mkdir  # 2.8.4
from ruffus import (
    add_inputs,
    formatter,
    merge,
    pipeline_printout_graph,
    pipeline_run,
    transform,
)


def cmd(command):
    """Runs all the command using Popen communicate."""
    sb.Popen(command, stdout=sb.DEVNULL).communicate()


# , stderr=sb.DEVNULL


def ranges(lst, given_gap=0):
    """
    A generator returns list of range based on given numbers
    [1,2,3,4,5, 10, 11, 12, 13, 14] => [[1,5], [10,14]]
    """
    lst_sorted = sorted(lst)
    init = 0
    for num in range(1, len(lst_sorted)):
        reported_gap = lst_sorted[num] - lst_sorted[num - 1] - 1
        if (lst_sorted[num] > given_gap + 1 + lst_sorted[num - 1]) or (
            reported_gap % 3 != 0
        ):
            # gap=0 means overlapping
            yield (lst_sorted[init], lst_sorted[num - 1])
            init = num
    yield (lst_sorted[init], lst_sorted[-1])


def useful_range(lst, gap=10):
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


def ab1seq(infile, tmp_fold):
    """ab1 to seq trimmed based on reference."""

    channels = ["DATA9", "DATA10", "DATA11", "DATA12"]  # G  # A  # T  # C
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

    for channel in channels:
        trace[bases[channel]] = record.annotations["abif_raw"][channel]

    nuc_df = {"nuc": [], "peak": []}

    for channel in zip(
        record.annotations["abif_raw"]["PBAS1"],
        record.annotations["abif_raw"]["PLOC1"],
    ):
        ambi_base = chr(channel[0])
        nuc_df["nuc"].append(ambi_base)
        if ambi_base in amb_bases:
            td = {}  # TODO: Please check what does td reprensts
            for base in amb_bases[ambi_base]:
                td[base] = trace[base][channel[1]]
            nuc_df["peak"].append(td)
        else:
            nuc_df["peak"].append({ambi_base: trace[ambi_base][channel[1]]})
    nuc_df = pd.DataFrame(nuc_df)
    peak = nuc_df["peak"].apply(lambda x: np.mean(list(x.values()))).values

    correct_peak_pos = np.where(peak > 100)[0]
    min_pos, max_pos = np.min(correct_peak_pos), np.max(
        correct_peak_pos
    )  # 1 is added as last position is not included
    nuc_df = nuc_df.loc[min_pos:max_pos]

    nuc_seq = "".join(nuc_df["nuc"].values)

    flb = path.split(infile)[1].split(".ab1")[
        0
    ]  # Consitering that the sample name is the first part

    with open(f"{tmp_fold}/{flb}.fasta", "w") as fout:
        fout.write(f">{flb}\n{nuc_seq}\n")

    command = [
        "blat",
        "-noHead",
        f"{tmp_fold}/CoVid_S_Gene.fasta",
        f"{tmp_fold}/{flb}.fasta",
        f"{tmp_fold}/{flb}.psl",
    ]  # Mapping against reference
    cmd(command)
    blat_df = (
        pd.read_table(f"{tmp_fold}/{flb}.psl", header=None)
        .sort_values(0, ascending=False)
        .drop_duplicates(9)
    )
    if len(blat_df) != 1:
        print(
            f"Looks like {infile} has multiple sequences. \
                    Please check. Exiting . . . . . "
        )
        sys.exit(0)

    for _, row in blat_df.iterrows():
        if row[8] == "-":
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

    nuc_df["idx"] = list(nuc_df.index)

    # df = enumerate_columns(df)
    # print(df)
    for nuc in peak_dict:
        peak_dict[nuc][f"{nuc}_idx"] = list(
            nuc_df[nuc_df[f"{nuc}"] != "-"].index
        )  # list(range(len(peak_dict[k]))
        nuc_df = nuc_df.merge(
            peak_dict[nuc], left_on="idx", right_on=f"{nuc}_idx", how="outer"
        )
        del nuc_df[str(nuc)]
    return nuc_df


def aln_clean(aln_df):
    """Clearning alihnment dataframe to remove unneccessary gaps and
    alignements."""
    one_side = ["1"]
    if aln_df.shape[1] == 5:
        idx = aln_df[aln_df["1_nuc"] != "-"].index
    else:
        one_side.append("2")
        idx = aln_df[
            ~((aln_df["1_nuc"] == "-") | (aln_df["2_nuc"] == "-"))
        ].index
    u_range = list(useful_range(list(idx), 10))
    # Add the importnat gap. Default: 10
    if len(one_side) == 1:
        # Avoid erroneous alignments at the end
        temp_aln_df = aln_df.loc[u_range[0]: u_range[0] + 8]
        if list(temp_aln_df["1_nuc"].values) != list(
            temp_aln_df["ref"].values
        ):
            aln_df.loc[u_range[0]: u_range[0] + 8, "1_nuc"] = "-"
            u_range[0] += 9
        temp_aln_df = aln_df.loc[u_range[1] - 8: u_range[1]]
        if list(temp_aln_df["1_nuc"].values) != list(
            temp_aln_df["ref"].values
        ):
            aln_df.loc[u_range[1] - 8: u_range[1], "1_nuc"] = "-"
            u_range[1] -= 9

    for col in aln_df.columns:
        if col == "ref":
            continue
        aln_df.loc[
            ~aln_df.index.isin(list(range(u_range[0], u_range[1] + 1))), col
        ] = "-"

    def rep_paired_base(lst):
        """Selecting reprentative base in presence of forward and
        reverse sanger sequencing data."""
        # TODO: Can be improved later selective query here
        if lst["ref"] == "-":
            return "-"
        if (
            (lst["1_nuc"] == lst["2_nuc"] == "-")
            or (lst["1_nuc"] != lst["2_nuc"] == lst["ref"])
            or (lst["2_nuc"] != lst["1_nuc"] == lst["ref"])
            or (lst["2_nuc"] == lst["1_nuc"] == lst["ref"])
            or (lst["2_nuc"] != lst["1_nuc"] != lst["ref"])
            or (
                (lst["2_nuc"] not in "ACGT-") and (lst["1_nuc"] not in "ACGT-")
            )
        ):
            return lst["ref"]
        if lst["1_nuc"] == lst["2_nuc"] != lst["ref"]:
            return lst["1_nuc"]
        if (lst["2_nuc"] not in "ACGT-") and (lst["1_nuc"] in "ACGT"):
            return lst["1_nuc"]
        if (lst["1_nuc"] not in "ACGT-") and (lst["2_nuc"] in "ACGT"):
            return lst["2_nuc"]
        # TODO: What should be returned when nothing works??

    def rep_single_base(lst):
        if lst["ref"] == "-":
            return "-"
        if (lst["1_nuc"] == "-") or (lst["1_nuc"] == lst["ref"]):
            return lst["ref"]
        if lst["1_nuc"] not in "ACGT-":
            return lst["ref"]
        return lst["1_nuc"]

    if len(one_side) > 1:
        aln_df["rep"] = aln_df.apply(rep_paired_base, axis=1)
    else:
        aln_df["rep"] = aln_df.apply(rep_single_base, axis=1)
    # print("".join(aln_df["rep"].values))
    return aln_df, u_range


def fasta_map2ref(infile, outfile, tmp_fold):
    """TODO: Docstring for fasta_map2ref.

    :arg1: TODO: Considering the the file has only on sequence
    :returns: TODO

    """
    fout = path.split(infile)[1] + ".psl"
    fout = f"{tmp_fold}/{fout}"
    command = [
        "blat",
        "-noHead",
        f"{tmp_fold}/CoVid_S_Gene.fasta",
        infile,
        fout,
    ]

    cmd(command)
    blat_df = pd.read_table(fout, header=None)
    blat_df = blat_df.sort_values(0, ascending=False)
    blat_df = blat_df.drop_duplicates(9).head(
        1
    )  # NOTE:Considering only best match

    sequences = {}
    for rec in SeqIO.parse(infile, "fasta"):
        sequences[rec.id] = rec.seq
    for rec in SeqIO.parse(f"{tmp_fold}/CoVid_S_Gene.fasta", "fasta"):
        sequences[rec.id] = rec.seq
    # TODO: Collect the details based

    output_file = open(outfile, "w")

    for _, row in blat_df.iterrows():
        if row[8] == "-":
            sequences[row[9]] = sequences[row[9]].reverse_complement()

        if row[17] == 1:
            seq = (
                sequences["ref"][: row[15]]
                + sequences[row[9]][row[11]: row[12]]
                + sequences["ref"][row[16]:]
            )
        else:

            tstarts = np.array(list(map(int, str(row[20])[:-1].split(","))))
            qstarts = np.array(list(map(int, str(row[19])[:-1].split(","))))
            blocks = np.array(list(map(int, str(row[18])[:-1].split(","))))

            qends = qstarts + blocks
            tends = tstarts + blocks
            begins = sequences["ref"][: row[15]]
            ends = sequences["ref"][row[16]:]
            mseq = ""
            for i in range(row[17] - 1):
                qgap = qstarts[i + 1] - qends[i]
                tgap = tstarts[i + 1] - tends[i]
                mseq += sequences[row[9]][qstarts[i]: qends[i]]

                if qgap:
                    if (qgap > 9) or (qgap % 3):
                        mseq += sequences["ref"][tends[i]: tstarts[i + 1]]
                    else:
                        mseq += "-" * qgap
                else:
                    if (tgap > 9) or (tgap % 3):
                        mseq += "-" * tgap
                    else:
                        mseq += sequences[row[9]][qends[i]: qstarts[i + 1]]
            seq = begins + mseq + ends
        output_file.write(f">{row[9]} {row[15]} {row[16]}\n{seq}\n")
    output_file.close()


def ab1_2seq_map2ref(infiles, outfile, tmp_fold):
    """TODO: Docstring for ab1_2seq.

    :infiles: TODO
    :returns: TODO

    """
    ab1seq_dfs = {}
    tsequences = {}
    flb = path.split(infiles[0])[1].split(".")[0]
    for i, fl in enumerate(infiles, start=1):
        ab1seq_dfs[i] = ab1seq(fl, tmp_fold)  # TODO: add the condition
        tsequences[i] = "".join(ab1seq_dfs[i]["nuc"].values)
    for rec in SeqIO.parse(f"{tmp_fold}/CoVid_S_Gene.fasta", "fasta"):
        tsequences[rec.id] = str(rec.seq)

    for k in ab1seq_dfs:
        ab1seq_dfs[k].columns = [f"{k}_{col}" for col in ab1seq_dfs[k].columns]

    aln_with_peak = merge_base_peak(
        aln_df_with_ref(tsequences, flb, tmp_fold), ab1seq_dfs
    )
    for n in ["1", "2"]:
        cl = f"{n}_nuc"
        if cl in aln_with_peak.columns:
            aln_with_peak.loc[pd.isnull(aln_with_peak[cl]), cl] = "-"
    aln_with_peak, u_range = aln_clean(aln_with_peak)
    if "2_nuc" in aln_with_peak.columns:
        usr = ranges(
            aln_with_peak[
                aln_with_peak[["ref", "1_nuc", "2_nuc"]].apply(
                    lambda row: True if "-" in list(row.values) else False,
                    axis=1,
                )
            ].index,
            10,
        )  # TODO Idiot fix it
    else:
        usr = ranges(
            aln_with_peak[
                aln_with_peak[["ref", "1_nuc"]].apply(
                    lambda row: True if "-" in list(row.values) else False,
                    axis=1,
                )
            ].index,
            10,
        )  # TODO: Idiot fix it
    if usr:
        for us in usr:
            if "2_nuc" in aln_with_peak.columns:
                if "".join(
                    aln_with_peak.loc[us[0]: us[1], "1_nuc"].values
                ) == "".join(aln_with_peak.loc[us[0]: us[1], "2_nuc"].values):
                    aln_with_peak.loc[
                        us[0]: us[1], "rep"
                    ] = aln_with_peak.loc[us[0]: us[1], "1_nuc"].values
            else:
                aln_with_peak.loc[us[0]: us[1], "rep"] = aln_with_peak.loc[
                    us[0]: us[1], "ref"
                ].values
    seq = "".join(list(aln_with_peak["rep"].values))

    output_file = open(outfile, "w")
    output_file.write(f">{flb} {u_range[0]} {u_range[1]}\n{seq}\n")
    output_file.close()


def ab2fasta(
    infiles, outfile  # , bc="neigh"
):  # Base criteria, max, neighbors, mixed # Inputfiles paired and none paired
    """Retains fasta and converts ab1 to fasta"""
    tmp_fold = path.split(outfile)[0]

    if len(infiles) == 1:
        if infiles[0].endswith(".fasta"):
            fasta_map2ref(infiles[0], outfile, tmp_fold)

    else:
        ab1_2seq_map2ref(infiles, outfile, tmp_fold)


def files_and_groups(sanger_files):
    """List fasta and ab1 sequences as dictionary.

    :sanger_files: List of files in the folder
    :returns: dictionary of files as paired/single and fasta/ab1

    """
    file_groups = {}
    for file_ in sanger_files:
        flx = path.split(file_)[1]
        if flx.endswith(".fasta"):
            flb = flx.split(".fasta")[0]
        elif flx.endswith(".ab1"):
            flb = flx.split(".")[0]  # Fist part should be the name
        else:
            print(f"{file_} doesn't have fasta or ab1 extention. Ignoring")
            continue
        if flb not in file_groups:
            file_groups[flb] = []
        file_groups[flb].append(file_)
    for group in file_groups:
        if len(file_groups[group]) > 2:
            print(
                f"{len(file_groups[group])} files are associated with id: {group}. \
                        Expecting 1 or 2 only"
            )
            sys.exit(0)
    for group in file_groups:
        filetype = []
        for file_ in file_groups[group]:
            filetype.append(file_.split(".")[-1])
        if len(set(filetype)) != 1:
            print("Multiple fileformat for {k}. Exiting . . . . . . ")
            sys.exit(0)
    sanger_outputs = list(file_groups.values())
    return sanger_outputs


@click.command()
@click.option(
    "-sa_ab1",
    help="ab1 folder or sanger sequence file",
    type=str,
    default="unitest/ab1_paired/ab1",
    show_default=True,
)  # Convert this to folder
# "/home/devil/Documents/San/Corona/Merging/Sanger/12April2021"
# @click.option("-fa", help="Fasta output file.
# If not given, only sequences will be printed in terminal",
#               type=str, default=None, show_default=True)
@click.option(
    "-asf",
    help="Assemblies folder containing fasta files",
    type=str,
    default="assemblies",
    show_default=True,
)  # /home/devil/Documents/San/Corona/Merging/Sanger/
# @click.option("-rf", help="Reference fasta",
#               type=bool, default=False, show_default=True)
# @click.option("-bs", help="Base replacement methos",
#               type=click.Choice([None, "max", "neigh", "ref"]),
# default=None, show_default=True)
@click.option(
    "-outd",
    help="Output Folder",
    type=str,
    default="Results",
    show_default=True,
)
def run(sa_ab1, asf, outd):  # , fa, asb, al, bs
    """
    Convert ab1 to Fasta based given parameters. must contain original
    sequence name., Must me N trimmed at the end
    work when either ab1 or sa parameter given. Must have ab1 files with ab1
    suffix with sequence name as prefix\n
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
    makedirs(tmp_fold, exist_ok=True)

    if not asf:
        print("Assembly folder not given. Exiting . . . . . . .")
        sys.exit(0)
    elif not path.exists(asf) or not path.isdir(asf):
        print(
            "Given assembly folder doesn't exist or path is not a folder. \
            Exiting . . . . ."
        )
        sys.exit(0)
    else:
        assemblies = glob(f"{asf}/*.fasta")

    # Later it can be replace with other genes
    s_gene_seq = """ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTT
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
        fout.write(f">ref\n{s_gene_seq}\n")

    sanger_outputs = files_and_groups(glob(f"{sa_ab1}/*"))

    @transform(
        sanger_outputs,
        formatter(r".+/(?P<filebase>\w+).(?P<tame>\w+)"),
        "%s/{filebase[0]}.fasta" % tmp_fold,
    )
    def ab1_2_fasta_r(inputfiles, outputfile):
        """TODO: Converts ab1 to fasta.

        :inp: TODO
        :returns: TODO

        """
        ab2fasta(inputfiles, outputfile)

    @merge(
        ab1_2_fasta_r, f"{tmp_fold}/sanger.fasta"
    )  # Apply merge it. Ignore  if it fails to include
    def sanger_seq_r(inputfiles, outputfile):
        with open(outputfile, "w") as fout:
            for fl in inputfiles:
                for rec in SeqIO.parse(fl, "fasta"):
                    fout.write(f">{rec.description}\n{rec.seq}\n")
        # print(inputfiles)

    @mkdir(outd)
    @transform(
        assemblies,
        formatter(r".+/(?P<filebase>\w+)"),  # NOTE: Anmol you have added an r
        add_inputs(sanger_seq_r),
        "%s/{filebase[0]}.fasta" % outd,
    )
    def alignments(infiles, outfile):
        # TODO: Add this function to run_pipeline function
        sang_file = infiles[1]
        assembly = infiles[0]
        sanger_seq = {}
        sanger_seq_desc = {}
        for rec in SeqIO.parse(sang_file, "fasta"):
            sanger_seq_desc[rec.id] = rec.description
            sanger_seq[rec.id] = str(rec.seq)

        org_seq = {}
        for rec in SeqIO.parse(assembly, "fasta"):
            org_seq[rec.id] = str(rec.seq)

        flb = path.split(assembly)[1].split(".fasta")[0]
        command = [
            "blat",
            "-noHead",
            assembly,
            sang_file,
            f"{tmp_fold}/{flb}.psl",
        ]
        cmd(command)
        blat_df = pd.read_table(f"{tmp_fold}/{flb}.psl", header=None)
        # return
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

                    print("You Idiot Contact Anmol")
        with open(outfile, "w") as fout:
            for k in org_seq:
                fout.write(f">{k}\n{org_seq[k]}\n")

    pipeline_run(verbose=9)
    pipeline_printout_graph("flowchart.svg", "svg", no_key_legend=True)


if __name__ == "__main__":
    run()
