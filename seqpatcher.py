#!/usr/bin/env python

# Author: Anmol Kiran
# Affiliation: University of Liverpool. UK
# and Malawi-Liverpool-Wellcome Trust, Malawi
# V0.1: 20/04/2021


import os
import subprocess as sb
import sys
from collections import defaultdict
from multiprocessing import cpu_count
from glob import glob
from os import makedirs, path, access
import tempfile as tmf

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
    """Runs all the command using Popen communicate.
    :command: in form of list ie ['ls', '-l']
    :return: None
    """
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


def ab1seq(infile, tmp_fold):
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
    # TODO: Autothreshold
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
        f"{tmp_fold}/ref.fasta",
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
        # NOTE: Considering both sequences should overlap
        # idx = aln_df[~((aln_df["1_nuc"] == "-") | (aln_df["2_nuc"] == "-"))].index
        # NOTE: Considering either sequence overlapping regions
        idx = aln_df[~((aln_df["1_nuc"] == "-") &
                       (aln_df["2_nuc"] == "-"))].index
        # TODO: Use below code if other team members allow
        # Checking with gap with common 3 nuc length gaps
        idx3 = aln_df[(aln_df["1_nuc"] == "-") &
                      (aln_df["2_nuc"] == "-")].index
        u3_range = list(useful_range(list(idx3), 1))
        # TODO: Select ranges with length of 3 and push them in the sequences

    u_range = list(useful_range(list(idx), 10))
    # Add the importnat gap. Default: 10
    if len(one_side) == 1:
        # Avoid erroneous alignments at the end
        temp_aln_df = aln_df.loc[u_range[0]: u_range[0] + 8]
        if list(temp_aln_df["1_nuc"].values) != list(temp_aln_df["ref"].values):
            aln_df.loc[u_range[0]: u_range[0] + 8, "1_nuc"] = "-"
            u_range[0] += 9
        temp_aln_df = aln_df.loc[u_range[1] - 8: u_range[1]]
        if list(temp_aln_df["1_nuc"].values) != list(temp_aln_df["ref"].values):
            aln_df.loc[u_range[1] - 8: u_range[1], "1_nuc"] = "-"
            u_range[1] -= 9

    for col in aln_df.columns:
        if col == "ref":
            continue
        aln_df.loc[
            ~aln_df.index.isin(list(range(u_range[0], u_range[1] + 1))), col
        ] = "-"

    def rep_paired_base(lst):
        """Selecting representative base in presence of forward and
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
            or ((lst["2_nuc"] not in "ACGT-") and (lst["1_nuc"] not in "ACGT-"))
        ):
            return lst["ref"]
        if lst["1_nuc"] == lst["2_nuc"] != lst["ref"]:
            return lst["1_nuc"]
        # if (lst["2_nuc"] not in "ACGT-") and (lst["1_nuc"] in "ACGT"):
        if (lst["2_nuc"] not in "ACGT") and (lst["1_nuc"] in "ACGT"):
            # NOTE: alowing gap in one base
            return lst["1_nuc"]
        # if (lst["1_nuc"] not in "ACGT-") and (lst["2_nuc"] in "ACGT"):
        if (lst["1_nuc"] not in "ACGT") and (lst["2_nuc"] in "ACGT"):
            # NOTE: alowing gap in one base
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


def fasta_map2ref(*args):
    # TODO: Use some part for arranging the sequence
    """Generate ouput file, sanger fasta integrated in refgne with range is
    description

    :args: infile, outfile, tmp_fold, idb, gap

    :returns: None

    """
    infile, outfile, tmp_fold, idb, gap = args
    fout = path.split(infile)[1] + ".psl"
    fout = f"{tmp_fold}/{fout}"
    command = [
        "blat",
        "-noHead",
        f"{tmp_fold}/ref.fasta",
        infile,
        fout,
    ]

    cmd(command)
    blat_df = pd.read_table(fout, header=None)
    blat_df = blat_df.sort_values(0, ascending=False)
    blat_df = blat_df.drop_duplicates(9).head(
        1)  # NOTE:Considering only best match

    sequences = {}
    flb = ""
    for rec in SeqIO.parse(infile, "fasta"):
        sequences[rec.id] = rec.seq
        flb = rec.id
    for rec in SeqIO.parse(f"{tmp_fold}/ref.fasta", "fasta"):
        sequences["ref"] = rec.seq

    for _, row in blat_df.iterrows():
        if row[8] == "-":
            sequences[row[9]] = sequences[row[9]].reverse_complement()
    aln_df = aln_df_with_ref(sequences, flb, tmp_fold)
    mapped_index = aln_df[aln_df[flb] != "-"].index
    u_range = useful_range(mapped_index, gap)
    ref = "".join(aln_df["ref"])
    seq = "".join(aln_df[flb])
    seq = ref[: u_range[0]] + seq[u_range[0]: u_range[1]] + ref[u_range[1]:]

    with open(outfile, "w") as fout:
        fout.write(f">{flb} {u_range[0]} {u_range[1]}\n{seq}\n")


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
    for rec in SeqIO.parse(f"{tmp_fold}/ref.fasta", "fasta"):
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
                if "".join(aln_with_peak.loc[us[0]: us[1], "1_nuc"].values) == "".join(
                    aln_with_peak.loc[us[0]: us[1], "2_nuc"].values
                ):
                    aln_with_peak.loc[us[0]: us[1], "rep"] = aln_with_peak.loc[
                        us[0]: us[1], "1_nuc"
                    ].values
            else:
                aln_with_peak.loc[us[0]: us[1], "rep"] = aln_with_peak.loc[
                    us[0]: us[1], "ref"
                ].values
    seq = "".join(list(aln_with_peak["rep"].values))
    # TODO: Generate sequence and exit

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
        # TODO: Auto detect forward and reverse
        elif flx.endswith(".ab1"):
            if flx.endswith("_F.ab1"):
                flb = flx.split("_F.ab1")[0]  # Fist part should be the name
            elif flx.endswith("_R.ab1"):
                flb = flx.split("_F.ab1")[0]
            else:
                print(f"{file_} file is not registered as forward or reverse.")
            # TODO: Perform pairwise corrections
        else:
            print(f"{file_} doesn't have fasta or ab1 extention. Ignoring")
            continue
        if flb not in file_groups:
            file_groups[flb] = []
        file_groups[flb].append(file_)
    for group in file_groups:
        if len(file_groups[group]) > 2:
            print(
                f"{len(file_groups[group])} files ({file_groups[group]}) are"
                f" associated with id: {group}. Expecting 1 fasta or 2 ab1"
                " files only"
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


def non_overlapping_ids(asseblies, ab1s):
    """Check for ovelapping and non-ovelapping ids and generates csv table

    :asseblies: Fasta assembly containing folder
    :ab1s: Sanger generated ab1 or fasta files
    :returns: Pandas dataframe

    """
    assembly_ids = []
    for fl in glob(f"{asseblies}/*.fasta"):
        for rec in SeqIO.parse(fl, "fasta"):
            assembly_ids.append(rec.id)

    sanger_fasta = []
    for fl in glob(f"{ab1s}/*.fasta"):
        for rec in SeqIO.parse(fl, "fasta"):
            sanger_fasta.append(rec.id)
    sanger_fasta_missing_assembly = ",".join(
        set(sanger_fasta) - set(assembly_ids))
    sanger_ab1_f = []

    for fl in glob(f"{ab1s}/*.ab1"):
        sanger_ab1_f.append(path.split(fl)[1].split("_F.ab1")[0])
    sanger_ab1_f_missing_assembly = ",".join(
        set(sanger_fasta) - set(sanger_ab1_f))
    sanger_ab1_r = []
    for fl in glob(f"{ab1s}/*_R.ab1"):
        sanger_ab1_r.append(path.split(fl)[1].split("_R.ab1")[0])
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

    data_frame["assembly"].append("No Assembly")
    data_frame["ab1_Forward"].append(sanger_ab1_f_missing_assembly)
    data_frame["ab1_Reverse"].append(sanger_ab1_r_missing_assembly)
    data_frame["fasta"].append(sanger_fasta_missing_assembly)
    return pd.DataFrame(data_frame)


@click.command()
@click.option(
    "-s",
    "--sanger-ab1",
    "sa_ab1",
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
    "-a",
    "-assemblies-foder",
    "asf",
    help="Assemblies folder containing fasta files",
    type=str,
    default="unitest/assemblies",
    show_default=True,
)  # /home/devil/Documents/San/Corona/Merging/Sanger/
# @click.option("-rf", help="Reference fasta",
#               type=bool, default=False, show_default=True)
# @click.option("-bs", help="Base replacement methos",
#               type=click.Choice([None, "max", "neigh", "ref"]),
# default=None, show_default=True)
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
    "-r",
    "--unpaired-ids",
    "mf",
    help="Report none ovelapping and IDs missing sanger seq",
    type=bool,
    default=True,
    show_default=True,
)
@click.option(
    "-t",
    "--missmatch-table",
    "mff",
    help="Mismatch id table csv file",
    type=str,
    default="mmf.csv",
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
@click.option(
    "-n",
    "--cpu",
    "cpu",
    help="Number of CPU to use",
    type=int,
    default=1,
    show_default=True,
)
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
    help="Allow gap of 3 nucleotide (only if supported by both forward and reverse)",
    type=bool,
    default=True,
    show_default=True,
)
@click.option(
    "-x",
    "--indel-selction",
    "idb",
    help="Replace Insertion, Deletion or Both",
    type=click.Choice(["del", "ins", "both"]),
    default="del",
    show_default=True,
    multiple=False,
)
def run(sa_ab1, asf, outd, mf, mff, ss, rf, cpu, ci, gap, n3, idb):  # , fa, asb, al, bs
    # TODO: Allowing insertion details only
    """
    Convert ab1 to Fasta based given parameters. must contain original
    sequence name., Must me N trimmed at the end
    work when either ab1 or sa parameter given. Must have ab1 files with ab1
    suffix with sequence name as prefix\n
    Will generate `merged_output.fasta` in current folder.

    max: max peak
    neigh: earlier peak if ambiguous nucleotide contain earlier base
    ref: Use reference base in case of ambiguity
    """
    if cpu < 1:
        print("Number of CPU use given is < 1.")
        print("Using default CPU value of 1")
        cpu = 1
    elif cpu > cpu_count() - 1:
        print("Given cpu usage is more or equal to cpus avalable on system")
        print(f"Setting CPU usage to {cpu_count() - 1 }")
        cpu = cpu_count() - 1
    else:
        pass

    seq_id_df = pd.DataFrame()  # Empty DF
    if mf:
        seq_id_df = non_overlapping_ids(asf, sa_ab1)
        if mff:
            seq_id_df.to_csv(mff, index=False)
        else:
            print(seq_id_df)

    tmp_fold = "/tmp"
    if access(tmp_fold, os.W_OK):
        tmp_fold = tmf.mkdtemp()
    else:
        print("/tmp is not writable. Trying current folder for temperory" " files")
        try:
            tmp_fold = tmf.mkdtemp(dir=".")
        except:
            exit("Current folder is not writable. Exiting . . . .")

    # TODO: Replace with a folder in /tmp or in current
    # makedirs(tmp_fold, exist_ok=True)
    # TODO: Separate all the assembled sequences here and list the in array
    if not asf:
        print("Assembly folder not given. Exiting . . . . . . .")
        sys.exit(0)
    elif not path.exists(asf) or not path.isdir(asf):
        print(
            f"Given assembly folder {asf} doesn't exist or path is not a folder."
            " Exiting . . . . ."
        )
        sys.exit(0)
    else:
        makedirs(f"{tmp_fold}/asseblies", exist_ok=True)
        print("Sepating assemblies . . . . .")
        # NOTE:Different file, same ids will result in overwriting the sequences
        for fl in glob(f"{asf}/*.fasta"):
            for rec in SeqIO.parse(fl, "fasta"):
                with open(f"{tmp_fold}/asseblies/{rec.id}.fasta", "w") as asbfile:
                    asbfile.write(f">{rec.id}\n{rec.seq}\n")
        assemblies = glob(f"{tmp_fold}/asseblies/*.fasta")

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

    with open(f"{tmp_fold}/ref.fasta", "w") as fout:
        if not rf:
            fout.write(f">ref\n{s_gene_seq}\n")
        else:
            if not path.exists(rf) or not path.isfile(rf):
                exit(f"{rf} path or file doesn't exist. Exiting . . . ")
            # TODO: Check if ref file exist
            seq_count = 0
            for rec in SeqIO.parse(rf, "fasta"):
                seq_count += 1
                if seq_count > 1:
                    exit(
                        f"{rf} contains more than 1 sequence. "
                        "Expect only one."
                        " Exiting."
                    )
                fout.write(f">{rec.id}\n{rec.seq}\n")
            if not seq_count:
                exit(
                    f"{rf} contains 0 (zero) sequence. " "Expect only one." " Exiting."
                )

    sanger_outputs = files_and_groups(glob(f"{sa_ab1}/*"))

    @transform(
        sanger_outputs,
        formatter(r".+/(?P<filebase>\w+).(?P<tame>\w+)"),
        "%s/{filebase[0]}.fasta" % tmp_fold,
        extras=[tmp_fold, gaps])
    def ab1_2_fasta_r(inputfiles, outputfile):
        """TODO: Converts ab1 to fasta.

        :input: ab1 file
        :outputfile: fasta file
        :returns: None

        """
        ab2fasta(inputfiles, outputfile)

    @merge(
        ab1_2_fasta_r, f"{tmp_fold}/sanger.fasta", extras=[ss]
    )  # Apply merge it. Ignore  if it fails to include
    def sanger_seq_r(inputfiles, outputfile, extras):
        # print(extras, "testing Anmol")
        extras_file = open(extras, "w") if extras else None
        with open(outputfile, "w") as fout:
            for fl in inputfiles:
                for rec in SeqIO.parse(fl, "fasta"):
                    fout.write(f">{rec.description}\n{rec.seq}\n")
                    if extras_file:
                        coors = rec.description.split()[1:]
                        coors = int(coors[0]), int(coors[1])
                        extras_file.write(
                            f">{rec.id}\n{rec.seq[coors[0]:coors[1]]}\n")
        if extras_file:
            extras_file.close()
        # TODO: if mergeing is allow, generate a local copy
        # TODO: Get input from extras

    @mkdir(outd)
    @transform(
        ab1_2_fasta_r,
        formatter(r".+/(?P<filebase>\w+)"),
        # TODO: Try to infer assemblies based on this
        "%s/{filebase[0]}.fasta" % outd,
    )
    def alignments(sanger_file, outfile):
        """Map sanger file to generated assemblies and insert sanger sequence."""

        assembly = f"{tmp_fold}/asseblies/{path.split(sanger_file)[1]}"
        sanger_seq = {}
        sanger_seq_desc = {}
        for rec in SeqIO.parse(sanger_file, "fasta"):
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
            sanger_file,
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

                    print("Please report at a bug at")
                    print(
                        "https://github.com/krisp-kwazulu-natal/"
                        "sars-cov-2-sequencing-merge-sanger/issues"
                    )
        with open(outfile, "w") as fout:
            for k in org_seq:
                fout.write(f">{k}\n{org_seq[k]}\n")

    pipeline_run(verbose=9, multithread=cpu)
    # pipeline_printout_graph("flowchart.svg", "svg", no_key_legend=True)


if __name__ == "__main__":
    run()