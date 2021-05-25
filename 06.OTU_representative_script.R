#!/usr/bin/Rscript

# Script for the selection of representative sequence within OTU
# Input file = table with 6 columns:
#    1     id = Sequence ID
#    2   size = Number of reads
#    3    seq = DNA sequences
#    4    len = Sequences length
#    5   ACTG = Percentage of ACTG bases
#    6      N = Percentage of ambiguous (non-ACTG) bases
# Table should be pre-sorted by number of ambiguous nucleotides and number of reads

## Usage:
# ./06.OTU_representative_script.R --input "Seqs_tab.txt" --output "Rep_seq.fasta"


###############################
############################### Parse input parameters
###############################

## Script variables that we need to parse
#   INPUTFILE - input file name
#   OUTFILE   - name of the output file with results

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Input file name"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output file"),
  make_option(c("-v", "--verbose"), action="store", default=FALSE, type='logical', help="Print info for debugging")
)
opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$input)){
  cat("Input file is not specified.\n", file=stderr())
  stop()
}
if(is.na(opt$output)){
  cat("Output file is not specified.\n", file=stderr())
  stop()
}

## Assign variables
INPUTFILE <- opt$input
OUTFILE <- opt$output

NCORES_DT <- 1             # Number of CPU cores for data.table

## Log assigned variables
if(opt$verbose == TRUE){
  cat(paste("Input file: ", INPUTFILE, "\n", sep=""))
  cat(paste("Output file: ", OUTFILE, "\n", sep=""))
}


###############################
############################### Load packages
###############################

suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(Biostrings) )

## Set the number of threads that data.table should use
setDTthreads(threads = NCORES_DT)


###############################
############################### Main workflow
###############################


## Load sequence table
sb <- fread(INPUTFILE, header = TRUE, sep = "\t")

## Estimate median sequence length
med <- median(rep(x = sb$len, times = sb$size))

## If the number of elements is even - there would be no exact len match
## Select the nearest value
ok_len <- sb$len[ which.min(abs(sb$len - med)) ]

## Subset sequences to the desired length 
## and get the first one (with the least number of ambiguous nucleotides and largest number of reads)
ok <- sb[ len == ok_len, ][1]

## Export sequence to fasta
seq <- DNAStringSet(ok$seq)
names(seq) <- ok$id
writeXStringSet(x = seq, filepath = OUTFILE, format="fasta", width = 9999)

