## Demultiplex reads

export NCORES=6
export MINSCORE=93

export DATAFLD="./00_Input_FASTQ/"
export RESFLD="./01_Demux_FASTQ/"

mkdir -p "$RESFLD"

## Demultiplex with Lima
lima_demux () {
  # input = folder path with one fastq.gz file with sequences and one fasta file with barcodes

  RUNNAME=$(basename $1)

  RESPATH=$(echo "$RESFLD""$RUNNAME")
  mkdir -p "$RESPATH"
 
  INPF=$(find $1 -type f -name "*.fastq.gz")
  INPB=$(find $1 -type f -name "*_barcodes.fasta")

  lima \
    --same --ccs \
    --single-side -W 250 \
    --min-length 10 \
    --min-score "$MINSCORE" \
    --split-named \
    --num-threads "$NCORES" \
    --log-level INFO --log-file "$RESPATH"/_log.txt \
    "$INPF" \
    "$INPB" \
    "$RESPATH"/lima.fq.gz

}
export -f lima_demux

## Process all sequencing runs (each should be in a separate folder)
find "$DATAFLD" -maxdepth 1 -mindepth 1 -type d \
  | parallel -j 1 lima_demux {}

## Rename files
cd "$RESFLD"

rename --filename \
  's/^lima.//g; s/--.*$/.fq.gz/' \
  $(find . -name "*.fq.gz")


## Compress logs
find . -name "lima.lima.report" | parallel -j 8 "gzip -8 {}"
