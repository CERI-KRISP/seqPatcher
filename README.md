# SARS-CoV-2 Sequencing: Merge Sanger

The script integrates the sanger sequenced SARS-CoV-2 S-gene into HTS based generated SARS-CoV-2 genome assemblies.


# Dependencies

## External (non-Python)

- MUSCLE
- BLAT

## Python modules

- Pandas
- Biopython
- click
- Ruffus

## File and sequence naming

- Assembly files should be written as <sample_name>.fasta
- If all the assemblies are in one file, name can be anything with fasta
extension.

- ab1 file should be named as <sample_name>.<F/R>.ab1

- If user already has as sanger fasta it should be named as
<sample_name>.U.fasta. Where U can be anything not containing '.'.
    - Sequence name should be as <sample_name>




# Execution

- execute `python ab1_integration.py --help` for help
- To execute in current folder you must have following folders
  - ab1: Containing ab1 files paired
  - assemblies: generated assembly files using HTS data
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

- [ ] Code cleaning
- [ ] Testing
- [ ] Clean intermediate files
- [ ] Start over option
- [ ] Bug Fixes
- [ ] fix bs option after discussion with San

# License

GPLv3
