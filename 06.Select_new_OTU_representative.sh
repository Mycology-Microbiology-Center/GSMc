## Script for the selection of alternative representative sequences
## In some OTUs, representative sequence could be a chimeric (too long) sequence
## Take the sequences with median length and highest number of reads

# 1. Extract sequences within OTU cluster (with ripgrep)
# 2. Remove old representative sequence from the cluster
# 3. Estimate sequence length
# 4. Select new representative sequence.
#    In case of several sequences with the same length (== median), pick the one with more reads
# 5. Re-annotate new representatives with BLAST


NCORES=8

## Function to extract new representative sequences
extract_seqs (){
  
  # $1 = input = OTU ID

  ## Databases
  UC="GL_OTUs.uc.gz"                       # VSEARCH-mapping file
  DB="GL_for_clustering.fasta.gz"          # Sequences used for OTU clustering
  RSCRIPT="06.OTU_representative_script.R" # R-script for representative sequence selection

  TMPDIR=$(mktemp -d -t "otu_XXXXXXXX")

  ## Extract OTU composition
  rg -z "$1" "$UC" > "$TMPDIR"/OTU.uc

  ## Get sequence IDs + remove old representative
  awk '{print $9}' "$TMPDIR"/OTU.uc \
    | sed -r 's/;size=[0-9]+//g' \
    | sort | uniq \
    | grep -v "$1" \
    > "$TMPDIR"/Seqs_in_OTU.txt

  ## Check if there are any sequences in this OTU (excluding old representative seq)
  if [[ $(head "$TMPDIR"/Seqs_in_OTU.txt | wc -l) -eq 0 ]]
  then
    echo "No sequences to choose from. Exit."
    rm -r "$TMPDIR"
    return 1
  fi

  ## Extract sequences
  rg -z -A 1 -f "$TMPDIR"/Seqs_in_OTU.txt --context-separator "" "$DB" | sed '/^$/d' > "$TMPDIR"/Seqs.fasta

  ## Estimate sequence length and number of ambiguous nucleotides, Extract number of reads
  ## Sort sequences by number of ambiguous nucleotides and number of reads
  seqkit fx2tab "$TMPDIR"/Seqs.fasta \
    --length --base-content ACTG --base-content NRYSWKMBDHV \
  | sed -r 's:\t+:\t:g' | sed 's/\t$//g' \
  | sed 's/;size=/\t/g' \
  | csvtk add-header -t -n id,size,seq,len,ACTG,N \
  | csvtk sort -t -k N:n -k size:nr \
  > "$TMPDIR"/Seqs_tab.txt

  ## Choose new representative
  "$RSCRIPT" --input "$TMPDIR"/Seqs_tab.txt --output "$TMPDIR"/"$1".fasta

  ## Compress results
  cat "$TMPDIR"/"$1".fasta | gzip > "$1".fasta.gz

  ## Copy cluster members (excluding old representative seq)
  cat "$TMPDIR"/Seqs.fasta | gzip > cluster_"$1".fasta.gz

  ## Clean up
  rm -r "$TMPDIR"
}

export -f extract_seqs

## Example:
# extract_seqs "0a1d6809390e863a55133f29d66462a98b533b40"


## Batch selection of new representative sequences.
# Input = text file with sequence IDs (`OTU_IDs.txt`), one ID per line
parallel -j "$NCORES" -a OTU_IDs.txt --progress 'extract_seqs {}'


