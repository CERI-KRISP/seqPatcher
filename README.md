# SARS-CoV-2 Sequencing: Merge Sanger

The script seqpatcher.py integrates the sanger sequenced SARS-CoV-2 S-gene segment into the corresponding HTS based SARS-CoV-2 genome assembly.

# Dependencies

## External (non-Python)

- [MUSCLE](https://www.drive5.com/muscle/downloads.htm)
- [BLAT](https://hgdownload.soe.ucsc.edu/admin/exe/)

Both the tools can be either be downloaded from provided links or installed using bioconda or ubuntu repositories. Tools must be in system the path, accessible via terminal.

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
3. If assembly files contain multiple sequences, assembly file will be generated for each sample ID as explained in point 2.
4. If all the assemblies are in one file, name can be anything with fasta extension. However, the sequences' names be same as <sample_name>.
5. ab1 file should be named as <sample_name>.\<xyz>.ab1. Where xyz can be any string.
6. If user already has as sanger fasta it should be named as <sample_name>.\<xyz>.fasta. Where xyz can be any string.

<a name='cmdoptions'>

# Commandline options

```
Usage: seqpatcher3.py [OPTIONS]

  Reports nucleotide sequence from Sanger chromatogram data based on user
  provided parameters and integrate that in assembly generated using NGS
  data

Options:
  -s, --sanger-ab1 TEXT           ab1 folder or sanger sequence file
                                  [default: sanger_ab1]

  -a, -assemblies-foder TEXT      Assemblies folder containing fasta files
                                  [default: unitest/assemblies]

  -o, --out-dir TEXT              Output Folder  [default: Results]
  -t, --tab TEXT                  CSV file for overlapping assemblies and
                                  sanger ids. If not given, stdout.

  -O, --output-fasta TEXT         Sanger output fasta from ab1
  -R, --ref-gene-fasta-file TEXT  Refence gene file in fasta format
  -c, --clean-intermediate BOOLEAN
                                  Remove intermediate file  [default: True]
  -g, --gap-allowed INTEGER       Gap Allowed between aligned fragment to
                                  consider the region continuous  [default:
                                  10]

  -3, --only-3-nuc BOOLEAN        Allow  3 nucleotide InDels else replace with
                                  reference nucleotides  [default: True]

  -x, --indel-selction [del|ins|both]
                                  Replace Insertion, Deletion or Both
                                  [default: del]

  --version                       Show the version and exit.
  --help                          Show this message and exit.
```

<a name="execution" />

# Execution

- Execute `python seqpatcher.py --help` for help
- Execute `seqpatcher.py --help` for help, if installed in system path.
- To execute in current folder you must have the following folders
  - ab1: Containing ab1 files single/paired or pre-compiled fasta.
  - assemblies: The folder contains assembly files generated using HTS data. These can be merged in one file.
  - Sanger sequences merged fasta will be saved in `Result` Folder.
  - Some files will be generated to user provided path.
    - Merged assemblies on one file.
    - Merged ab1 to fasta converted files.

# Base selection

- **Note**: Below cases are valid when InDel are smaller than 10 and not multiple of three

## Paired ab1

InDels muliple of 3 nucleotides can be permitted else below rule will be applied.

| Ref    | Forward    | Reverse    | Final Outcome                                                   |
| ------ | ---------- | ---------- | --------------------------------------------------------------- |
| -      | Any base/- | Any base/- | -                                                               |
| Base A | Base B     | Base B     | Base B                                                          |
| Base A | Base A     | Base B     | Base A                                                          |
| Base A | Base B     | Base A     | Base A                                                          |
| Base A | Base B     | -          | Base A                                                          |
| Base A | Base A     | -          | Base A                                                          |
| Base A | -          | -          | Base A                                                          |
| Base A | Base B     | Base Ambi  | Base B if (Base B is in Base Ambi) else Base A                  |
| Base A | Base Ambi  | Base Ambi  | Single common base in Ambi bases (excluding Base A) else Base A |
| Base A | Base B     | Base C     | Base A                                                          |

## Single ab1

| Reference | Forward/Reverse | Selected     |
| --------- | --------------- | ------------ |
| -         | Base A          | Base A       |
| Base A    | -               | Base A       |
| Base A    | Base B          | Base B       |
| Base A    | Ambi            | highest peak |

## Single Fasta

| Reference | FastaNucleotides | Selected |
| --------- | ---------------- | -------- |
| -         | Base A           | Base A   |
| Base A    | -                | Base A   |
| Base A    | Base B           | Base B   |
| Base A    | Ambi             | Base A   |

## Additional Features

1. Generate table of overlapping IDs of assemblies and sanger sequence data.
2. Converts sanger ab1 to fasta based on above criteria.

# License

GPLv3
