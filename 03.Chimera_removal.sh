## Script for sample-wise chimera removal
##   1. reference-based chimera removal
##   2. de novo chimera search. Sequences identified as putative chimeras will have lower priority in clustering

## Input = FASTQ, Output = FASTQ

cat > chimera_rm.sh <<'EOT'
# #!/bin/bash

  # $1 = input file name
  # $2 = output prefix
  # $3 = log file
  # $4 = reference database path

  echo -e "Start " "$1" > "$3"
  TMPDIR=$(mktemp -d -t "chimera_XXXXXXX")

  echo -e "\n### Dereplication\n" >> "$3"

  ## Dereplicate
  vsearch \
    --derep_fulllength "$1" \
    --output - \
    --fasta_width 0 \
    --threads 1 \
    --sizein --sizeout \
    2>> "$3" \
  | gzip -2 > "$TMPDIR"/derep.fa.gz


  echo -e "\n### Reference-based chimera removal\n" >> "$3"

  ## Reference-based
  vsearch \
    --uchime_ref "$TMPDIR"/derep.fa.gz \
    --db "$4" \
    --fasta_width 0 \
    --sizein --sizeout \
    --threads 1 \
    --chimeras - \
    --nonchimeras "$TMPDIR"/NoChimera_Ref.fasta \
    2>> "$3" \
  | gzip -2 > "$TMPDIR"/Ref_Chimeras.fasta.gz

  echo -e "\n### De novo chimera removal\n" >> "$3"

  ## De-novo
  numseqs=$(grep -c "^>" "$TMPDIR"/derep.fa.gz)
  if [[ $numseqs -gt 3 ]] ; then

    vsearch \
      --uchime_denovo "$TMPDIR"/derep.fa.gz \
      --chimeras "$TMPDIR"/Chimeras_Denovo.fasta \
      --nonchimeras - \
      --fasta_width 0 \
      --sizein --xsize \
      2>> "$3" \
    | gzip -2 > "$TMPDIR"/NoChimera_DeNovo.fasta.gz

  else
    echo "The number of sequences is too small. Do not perform de novo search." >> "$3"
  fi


  echo -e "\n### Sequence extraction\n" >> "$3"

  ## Get sequences
  if [ -f "$TMPDIR"/NoChimera_DeNovo.fasta.gz ]; then
    zcat "$TMPDIR"/NoChimera_DeNovo.fasta.gz | awk '$0 !~ /^>/ {print toupper($0)}' |  sort > "$TMPDIR"/tmp_ok.txt
  elif [ -f "$TMPDIR"/NoChimera_Ref.fasta ]; then
    awk '$0 !~ /^>/ {print toupper($0)}' "$TMPDIR"/NoChimera_Ref.fasta | sort > "$TMPDIR"/tmp_ok.txt
  else
    touch "$TMPDIR"/tmp_ok.txt
  fi

  if [ -f "$TMPDIR"/Chimeras_Denovo.fasta ]; then
    awk '$0 !~ /^>/ {print toupper($0)}' "$TMPDIR"/Chimeras_Denovo.fasta | sort > "$TMPDIR"/tmp_chim_denovo.txt
  else
    touch "$TMPDIR"/tmp_chim_denovo.txt
  fi

  if [ -f "$TMPDIR"/Ref_Chimeras.fasta.gz ]; then
    zcat "$TMPDIR"/Ref_Chimeras.fasta.gz | awk '$0 !~ /^>/ {print toupper($0)}' |  sort > "$TMPDIR"/chim_ref.txt
  else
    touch "$TMPDIR"/chim_ref.txt
  fi


  ## Remove de-novo chimeras identified with reference database
  comm -13 "$TMPDIR"/chim_ref.txt "$TMPDIR"/tmp_ok.txt > "$TMPDIR"/ok.txt

  ## Remove seqs identified as reference chimeras
  comm -13 "$TMPDIR"/chim_ref.txt "$TMPDIR"/tmp_chim_denovo.txt > "$TMPDIR"/chim_denovo.txt

  echo -e "Number of unique non-chimeric sequences: " $(wc -l < "$TMPDIR"/ok.txt) >> "$3"
  echo -e "Number of unique reference-based chimeric sequences: " $(wc -l < "$TMPDIR"/chim_ref.txt) >> "$3"
  echo -e "Number of unique de novo chimeric sequences: " $(wc -l < "$TMPDIR"/chim_denovo.txt) >> "$3"

  #### To improve speed, split patterns into chunks (with 100 lines each), then ripgrep

  ## Good sequences
  if [[ -s "$TMPDIR"/ok.txt ]]; then

      split -l 100 "$TMPDIR"/ok.txt "$TMPDIR"/patt.split.
      for CHUNK in "$TMPDIR"/patt.split.* ; do
        rg -z -B 1 -A 2 -x -f "$CHUNK" "$1" >> "$TMPDIR"/tmp_NoChimera.fq
      done
      rm "$TMPDIR"/patt.split.*

      sed '/^--$/d' "$TMPDIR"/tmp_NoChimera.fq | gzip -6 > "$2"_NoChimera.fq.gz
  fi

  ## De novo chimeras
  if [[ -s "$TMPDIR"/chim_denovo.txt ]]; then

    split -l 100 "$TMPDIR"/chim_denovo.txt "$TMPDIR"/patt.split.
    for CHUNK in "$TMPDIR"/patt.split.* ; do
      rg -z -B 1 -A 2 -x -f "$CHUNK" "$1" >> "$TMPDIR"/tmp_ChimDeNov.fq
    done
    rm "$TMPDIR"/patt.split.*

    sed '/^--$/d' "$TMPDIR"/tmp_ChimDeNov.fq | gzip -6 > "$2"_ChimDeNov.fq.gz
  fi

  ## Reference-based chimeras
  if [[ -s "$TMPDIR"/chim_ref.txt ]]; then
    split -l 100 "$TMPDIR"/chim_ref.txt "$TMPDIR"/patt.split.
    for CHUNK in "$TMPDIR"/patt.split.* ; do
      rg -z -B 1 -A 2 -x -f "$CHUNK" "$1" >> "$TMPDIR"/tmp_ChimRef.fq
    done
    rm "$TMPDIR"/patt.split.*

    sed '/^--$/d' "$TMPDIR"/tmp_ChimRef.fq  | gzip -6 > "$2"_ChimRef.fq.gz
  fi


  ## Clean up
  rm -r "$TMPDIR"

  echo -e "\nDone " "$1" >> "$3"

EOT

chmod +x chimera_rm.sh


## Example:
# ./chimera_rm.sh Sample.fq.gz Samp Samp.log Reference_chimera_DB.fasta.gz
