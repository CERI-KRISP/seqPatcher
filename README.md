# SARS-CoV-2 Sequencing: Merge Sanger

The script seqpatcher.py integrates the sanger sequenced SARS-CoV-2 S-gene segment into the corresponding HTS based SARS-CoV-2 genome assembly.


# Dependencies

## External (non-Python)

- [MUSCLE](https://www.drive5.com/muscle/downloads.htm)
- [BLAT](https://hgdownload.soe.ucsc.edu/admin/exe/)

Both the tools can be either be downloaded from provided links or installed using bioconda or ubuntu repositories.  Tools must be in system path, accessible via terminal.

## Python modules

- pandas==1.2.4
- biopython==1.78
- click==7.1.2
- ruffus==2.8.4

### installation
- Install dependencies using `pip install -r requirements.txt` , download `seqpatcher.py` in your working directory and execute as instructed in [Execution](#execution) section.
- Install using `pip install git+https://github.com/krisp-kwazulu-natal/sars-cov-2-sequencing-merge-sanger` or `git clone https://github.com/krisp-kwazulu-natal/sars-cov-2-sequencing-merge-sanger && cd sars-cov-2-sequencing-merge-sanger && python setup.py install`.  `seqpatcher.py` should be installed in system path.



## File and sequence naming

- Assembly files should be written as <sample_name>.fasta. Sequences' names be same as <sample_name>
- If all the assemblies are in one file, name can be anything with fasta
extension. However, the sequences' names be same as <sample_name>

- ab1 file should be named as <sample_name>.<F/R>.ab1

- If user already has as sanger fasta it should be named as
<sample_name>.U.fasta. Where U can be anything not containing '.'.
    - Sequence name should be as <sample_name>



<a name="execution" />
# Execution

- execute `python seqpatcher.py --help` for help
- execute `seqpatcher.py --help` for help, if installed in system path.
- To execute in current folder you must have following folders
  - ab1: Containing ab1 files paired
  - assemblies: Assembly files generated using HTS data. These can be merged in one file.
  - Result will in `Result` Folder


# Base selection

- **Note**: Below cases are valid when InDel are smaller than 10 and not multiple of three
## Paired ab1

|Ref|Forward|Reverse|Final Outcome|
|---|---|---|---|
|-|Any base/-|Any base/-|-|
|Base A|Base B|Base B|Base B|
|Base A|Base A|Base B|Base A|
|Base A|Base B|Base A|Base A|
|Base A|Base B|-|Base A|
|Base A|Base A|-|Base A|
|Base A|-|-|Base A|
|Base A|Base B|Base Ambi|Base B|
|Base A|Base B|Base C|Base A|




## Single ab1

|Reference|Forward/Reverse|Selected|
|---------|---------------|--------|
|-|Any base/-|-|
|Base A|-|Base A|
|Base A|Base B|Base B|
|Base A|Ambi|highest peak|


## Single Fasta
|Reference|FastaNucleotides|Selected|
|---------|---------------|--------|
|-|Any base|-|
|Base A|-|Base A|
|Base A|Base B|Base B|
|Base A|Ambi|Base A|




# TODOs

- [x] Testing
- [ ] Clean intermediate files
- [x] Start over option
- [x] fix bs option after discussion with San

# License

GPLv3
