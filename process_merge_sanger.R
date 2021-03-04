library(readr)
library(dplyr)
library(writexl)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--nc_sanger_file", type="character", default = "./example_data/nc_sanger_seqs.csv",
                    help="Provide file containing nextclade results for sanger sequences")
parser$add_argument("--nc_merged_file", type="character", default = "./example_data/nc_runNN_seqs_aln_withSanger.csv",
                    help="Provide file containing nextclade results for merged file")
args <- parser$parse_args()
if (args$nc_sanger_file == "no_input_return_error" ) { stop("file not found and is a required argument.", call.=FALSE)}
if (args$nc_merged_file == "no_input_return_error" ) { stop("file not found and is a required argument.", call.=FALSE)}


## path to nextclade output of sanger sequences
#nc_sanger <- "sanger/compile/nextclade_sanger_all.csv"
path_sanger_nc_file <- args$nc_sanger_file
## path to merged sequences
#nc_merged_seqs <- "~/temp/Collaborations/covid/assembly_auto/merge_sanger/sanger/compile/nextclade_all_seqs_withSanger2.csv"
path_merged_seqs_nc_file <-  args$nc_merged_file

## specify work directory if different from current directory
wd <- dirname(path_merged_seqs_nc_file)

## Read in files
nc_sanger_seqs <- read_delim(path_sanger_nc_file, ";", escape_double = FALSE, trim_ws = TRUE)
nc_merged_seqs <- read_delim(path_merged_seqs_nc_file, ";", escape_double = FALSE, trim_ws = TRUE)

## Trim sequence name suffix
nc_sanger_seqs$seqName <- gsub("_cons_seq", "", nc_sanger_seqs$seqName)

## Get just the sequence names and aa subtitutions
nc_sanger_seqs_sub <- nc_sanger_seqs %>%  select(seqName, aaSubstitutions)
nc_merged_seqs_sub <- nc_merged_seqs %>% select(seqName, aaSubstitutions)

## Merge the files to see if all mutations on the left are present on the right
nc_sanger_seqs_compare <- inner_join(nc_sanger_seqs_sub, nc_merged_seqs_sub, by = "seqName")
nc_sanger_seqs_compare_full <- full_join(nc_sanger_seqs_sub, nc_merged_seqs_sub, by = "seqName")

## Save report to excel
write_xlsx(nc_sanger_seqs_compare, paste0(wd, "/nc_merge_report.xlsx"))
write_xlsx(nc_sanger_seqs_compare_full, paste0(wd, "/nc_merge_report_full.xlsx"))