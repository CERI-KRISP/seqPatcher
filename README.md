# SARS-CoV-2 Sequencing: Merge Sanger

The script seqpatcher.py integrates the sanger sequenced SARS-CoV-2 S-gene segment into the corresponding HTS based SARS-CoV-2 genome assembly.

# Dependencies

## External (non-Python)

- [MUSCLE](https://www.drive5.com/muscle/downloads.htm)
- [BLAT](https://hgdownload.soe.ucsc.edu/admin/exe/)

Both the tools can be either be downloaded from provided links or installed using bioconda or ubuntu repositories. Tools must be in system path, accessible via terminal.

## Python modules

- pandas==1.2.4
- biopython==1.78
- click==7.1.2

### installation

- Install dependencies using `pip install -r requirements.txt`, download `seqpatcher.py` in your working directory and execute as instructed in [Execution](#execution) section.
- Install using `pip install git+https://github.com/krisp-kwazulu-natal/seqPatcher` or `git clone https://github.com/krisp-kwazulu-natal/seqPatcher && cd seqPatcher && python setup.py install`. `seqpatcher.py` should be installed in system path.

## File and sequence naming

1. File name or sample name must not have special characters.
2. Assembly files should be written as <sample_name>.fasta. Sequence's name should be the same as <sample_name>.
3. If assembly files contain multiple sequences, assembly file will be generated for each sample ID as explained in point 1.
4. If all the assemblies are in one file, name can be anything with fasta extension. However, the sequences' names be same as <sample_name>.
5. ab1 file should be named as <sample_name>.<xyz>.ab1. Where xyz can be any string.
6. If user already has as sanger fasta it should be named as <sample_name>.<xyz>.fasta. Where xyz can be any string.

<a name="execution" />
# Execution

- Execute `python seqpatcher.py --help` for help
- Execute `seqpatcher.py --help` for help, if installed in system path.
- To execute in current folder you must have the following folders
  - ab1: Containing ab1 files single/paired or pre-compiled fasta.
  - assemblies: Assembly files generated using HTS data. These can be merged in one file.
  - Sanger sequences merged fasta will be saved in `Result` Folder.
  - Some files will be generated user provided path.
    - Merged assemblies on one file.
    - Merge ab1 to fasta converted samples.

# Base selection

- **Note**: Below cases are valid when InDel are smaller than 10 and not multiple of three

## Paired ab1

InDels muliple of 3 nucleotides can be permitted else below rule will be applied.

| Ref    | Forward    | Reverse    | Final Outcome |
| ------ | ---------- | ---------- | ------------- |
| -      | Any base/- | Any base/- | -             |
| Base A | Base B     | Base B     | Base B        |
| Base A | Base A     | Base B     | Base A        |
| Base A | Base B     | Base A     | Base A        |
| Base A | Base B     | -          | Base A        |
| Base A | Base A     | -          | Base A        |
| Base A | -          | -          | Base A        |
| Base A | Base B     | Base Ambi  | Base B        |
| Base A | Base B     | Base C     | Base A        |

## Single ab1

| Reference | Forward/Reverse | Selected     |
| --------- | --------------- | ------------ |
| -         | Any base/-      | -            |
| Base A    | -               | Base A       |
| Base A    | Base B          | Base B       |
| Base A    | Ambi            | highest peak |

## Single Fasta

| Reference | FastaNucleotides | Selected |
| --------- | ---------------- | -------- |
| -         | Any base         | -        |
| Base A    | -                | Base A   |
| Base A    | Base B           | Base B   |
| Base A    | Ambi             | Base A   |

## Additional Features

1. Generate table of overlapping IDs of assemblies and sanger sequence data.
2. Converts sanger ab1 to fasta based on above criteria.

# License

GPLv3
