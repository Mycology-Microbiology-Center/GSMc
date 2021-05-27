## Prepare PacBio data for clustering:
## - trim primers and sort sequences to full-length and partial
##   full-length seqs will have higher priority in OTU clustering
## - extract ITS region with ITSxpress
##   full-ITS will have higher priority in OTU clustering
## - trim SSU and LSU from partial sequence -- lower priority

###################################################
################################################### Cutadapt - trim primers
###################################################

mkdir -p 02_PrimerCut

export NCORES=1
export MINLEN=80
export PRIMER_ERROR=3
export PRIMER_OVERLAP=16

## Function to convert IUPAC codes in primers. Author Sten Anslan, developed for PipeCraft
convert_IUPAC () {
  echo $1 | \
  if grep -q -E "R|Y|S|W|K|M|B|D|H|V|N|I"; then
      #define IUPAC codes
      R=$"[AG]"
      Y=$"[CT]"
      S=$"[GC]"
      W=$"[AT]"
      K=$"[GT]"
      M=$"[AC]"
      B=$"[CGT]"
      D=$"[AGT]"
      H=$"[ACT]"
      V=$"[ACG]"
      N=$"[ATGC]"
      I=$"[ATGC]"
      #replace IUPAC codes
      primer=$(echo $1 | \
      sed -e "s/R/$R/g; s/Y/$Y/g; \
      s/S/$S/g; s/W/$W/g; s/K/$K/g; \
      s/M/$M/g; s/B/$B/g; s/D/$D/g; \
      s/H/$H/g; s/V/$V/g; s/N/$N/g; \
      s/I/$I/g")
      #return convered primer
      echo $primer
  else
      #return original primer when no IUPAC codes were detected
      echo $1
  fi
}
export -f convert_IUPAC

## Function to construct reverse-comlement sequence
RC () {
  echo "$1" | tr "[ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv]" "[TACGAtacgaNnRrYySsWwMmKkVvHhDdBb]" | rev
}
export -f RC

export PRIMER_F="GTACACACCGCCCGTCG"
export PRIMER_R="CCTSCSCTTANTDATATGC"
export PRIMER_Fr=$( RC "$PRIMER_F")
export PRIMER_Rr=$( RC "$PRIMER_R")


## Trim primers and correct read orientation
primer_trim () {

  # $1 = input
  # $2 = output prefix
  # $3 = log

  echo -e "Start " "$1" > "$3"
  TMPDIR=$(mktemp -d -t "cutadapt_XXXXXXX")

  echo -e "Full-length sequence search (both primers required)\n" >> "$3"

  ## Both primers in correct orientation
  ## However, multiple primer occurrences possible (-> use `times`)
  cutadapt \
    -g "${PRIMER_F};required"..."${PRIMER_Rr};required" \
    -e "$PRIMER_ERROR" \
    --overlap "$PRIMER_OVERLAP" \
    --minimum-length "$MINLEN" \
    --times 4 \
    --revcomp \
    --cores "$NCORES" \
    --output "$TMPDIR"/trimmed_tmp.fq.gz \
    --untrimmed-output "$TMPDIR"/untrimmed_tmp.fq.gz \
    "$1" >> "$3"
  

  echo -e "\nRemove full sequences with additional sequences present \n" >> "$3"

  ## Remove sequences with incorrect primer orientation
  cutadapt \
    -g "$PRIMER_F" -g "$PRIMER_Fr" -g "$PRIMER_R" -g "$PRIMER_Rr" \
    -e "$PRIMER_ERROR" \
    --overlap "$PRIMER_OVERLAP" \
    --times 4 \
    --cores "$NCORES" \
    --discard-trimmed \
    --output "$2"_trimmed_Full.fq.gz \
    "$TMPDIR"/trimmed_tmp.fq.gz \
    >> "$3"

  echo -e "\nPrepare data for parial sequence trimming \n" >> "$3"

  
  ## Trim sequences to the first and last 50 bp
  seqkit subseq -r 1:50   "$TMPDIR"/untrimmed_tmp.fq.gz | seqkit replace -p "\s.+" | gzip -2 > "$TMPDIR"/utr_L.fq.gz
  seqkit subseq -r -50:-1 "$TMPDIR"/untrimmed_tmp.fq.gz | seqkit replace -p "\s.+" | gzip -2 > "$TMPDIR"/utr_R.fq.gz
  cat "$TMPDIR"/utr_L.fq.gz "$TMPDIR"/utr_R.fq.gz > "$TMPDIR"/utr_LR.fq.gz 
  
  ## Find primers
  rg -z $(convert_IUPAC $PRIMER_F) -B 1 --context-separator "" "$TMPDIR"/utr_LR.fq.gz | sed '/^$/d' | sed 's/@/>/' | seqkit seq --name > "$TMPDIR"/utr_F_names_tmp.txt
  rg -z $(convert_IUPAC $PRIMER_Fr) -B 1 --context-separator "" "$TMPDIR"/utr_LR.fq.gz | sed '/^$/d' | sed 's/@/>/' | seqkit seq --name >> "$TMPDIR"/utr_F_names_tmp.txt
  
  rg -z $(convert_IUPAC $PRIMER_R) -B 1 --context-separator "" "$TMPDIR"/utr_LR.fq.gz | sed '/^$/d' | sed 's/@/>/' | seqkit seq --name > "$TMPDIR"/utr_R_names_tmp.txt
  rg -z $(convert_IUPAC $PRIMER_Rr) -B 1 --context-separator "" "$TMPDIR"/utr_LR.fq.gz | sed '/^$/d' | sed 's/@/>/' | seqkit seq --name >> "$TMPDIR"/utr_R_names_tmp.txt
  
  ## Exclude duplicated reads
  cat "$TMPDIR"/utr_F_names_tmp.txt | sort | uniq > "$TMPDIR"/utr_F_names.txt
  comm -13 "$TMPDIR"/utr_F_names.txt <(sort "$TMPDIR"/utr_R_names_tmp.txt) > "$TMPDIR"/utr_R_names.txt
  
  ## Extract reads
  seqkit grep -f "$TMPDIR"/utr_F_names.txt "$TMPDIR"/untrimmed_tmp.fq.gz -o "$TMPDIR"/utr_F.fq.gz
  seqkit grep -f "$TMPDIR"/utr_R_names.txt "$TMPDIR"/untrimmed_tmp.fq.gz -o "$TMPDIR"/utr_R.fq.gz
  
  echo -e "\nParial sequence trimming - F at 5'-end\n" >> "$3"

  ## If F primer is found within first 50 bases
  # fix read orientation, optionally remove (multiple) R primers, 
  # remove short sequences, and remove reads with duplicated F primers
  cutadapt \
      -g "$PRIMER_F" \
      -e "$PRIMER_ERROR" \
      --overlap "$PRIMER_OVERLAP" \
      --revcomp \
      --discard-untrimmed \
      "$TMPDIR"/utr_F.fq.gz \
      2> "$TMPDIR"/F_log_01.log \
    | cutadapt \
      -a "$PRIMER_R" -a "$PRIMER_Rr" \
      -e "$PRIMER_ERROR" \
      --overlap "$PRIMER_OVERLAP" \
      --times 4 \
      - \
      2> "$TMPDIR"/F_log_02.log \
    | cutadapt \
      -g "${PRIMER_F}" -a "${PRIMER_Fr}" \
      -e 2 \
      --overlap 17 \
      --times 4 \
      --minimum-length "$MINLEN" \
      --output "$TMPDIR"/F_trimmed.fq.gz \
      - > "$TMPDIR"/F_log_03.log

  ## Combine logs
  echo -e "\n...Trim F (required)\n" >> "$3"
  cat "$TMPDIR"/F_log_01.log >> "$3"
  echo -e "\n...Trim R and Rr (optional)\n" >> "$3"
  cat "$TMPDIR"/F_log_02.log >> "$3"
  echo -e "\n...Trim duplicated 5'-F or 3'-Fr (optional)\n" >> "$3"
  cat "$TMPDIR"/F_log_03.log >> "$3"


  echo -e "\nParial sequence trimming - R at 5'-end\n" >> "$3"

  ## If R primer is found within first 50 bases (excluding sequences with F in the first 50)
  # fix read orientation, optionally remove (multiple) F primers, 
  # remove short sequences, and remove reads with duplicated R primers
  cutadapt \
      -g "$PRIMER_R" \
      -e "$PRIMER_ERROR" \
      --overlap "$PRIMER_OVERLAP" \
      --revcomp \
      --discard-untrimmed \
      "$TMPDIR"/utr_R.fq.gz \
      2> "$TMPDIR"/R_log_01.log \
    | cutadapt \
      -a "$PRIMER_F" -a "$PRIMER_Fr" \
      -e "$PRIMER_ERROR" \
      --overlap "$PRIMER_OVERLAP" \
      --times 6 \
      - \
      2> "$TMPDIR"/R_log_02.log \
    | cutadapt \
      -g "${PRIMER_R}" -a "${PRIMER_Rr}" \
      -e 2 \
      --overlap 17 \
      --times 6 \
      --minimum-length "$MINLEN" \
      --output "$TMPDIR"/R_trimmed.fq.gz \
      - > "$TMPDIR"/R_log_03.log
  
  ## Combine logs
  echo -e "\n...Trim R (required)\n" >> "$3"
  cat "$TMPDIR"/R_log_01.log >> "$3"
  echo -e "\nTrim F and Fr (optional)\n" >> "$3"
  cat "$TMPDIR"/R_log_02.log >> "$3"
  echo -e "\n...Trim duplicated 5'-R or 3'-Rr (optional)\n" >> "$3"
  cat "$TMPDIR"/R_log_03.log >> "$3"

  ## Combine partial sequences
  cat "$TMPDIR"/F_trimmed.fq.gz "$TMPDIR"/R_trimmed.fq.gz > "$2"_trimmed_Part.fq.gz

  ## Clean temp files
  rm -r "$TMPDIR"

  echo -e "Done " "$1" >> "$3"
}
export -f primer_trim


find "01_Demux" -name "*.fq.gz" | rush -j 1 \
  " mkdir -p ./02_PrimerCut/{/%}; \
    mkdir -p ./02_PrimerCut/logs; \
    primer_trim {} ./02_PrimerCut/{/%}/{%:} ./02_PrimerCut/logs/{%:}.log
  "


###################################################
################################################### Prepare data for ITS extraction
###################################################

mkdir -p Primer_trim

NCORES=1

## Combine full-length sequences, add sample ID to header
find ./02_PrimerCut/ -name "*_trimmed_Full.fq.gz" | rush -j "$NCORES" \
  "seqkit replace -p '\s.+' -r ';sample={%:}' {} -o ./Primer_trim/{%:}_renamed.fq.gz"

cat ./Primer_trim/*.fq.gz > GL_trimmed_Full_Merg.fq.gz


## Combine partial sequences, add sample ID to header
find ./02_PrimerCut/ -name "*_trimmed_Part.fq.gz" | rush -j "$NCORES" \
  "seqkit replace -p '\s.+' -r ';sample={%:}' {} -o ./Primer_trim/{%:}_renamed.fq.gz"

cat ./Primer_trim/*.fq.gz > GL_trimmed_Part_Merg.fq.gz




###################################################
################################################### Run ITSxpress on HPC
###################################################

## Optionally, split data into parts to run ITSxpress in parallel


## ITSxpress - extract full-ITS region
itsxpress \
  --fastq trimmed_reads.fq.gz \
  --cluster_id 1.0 \
  --single_end \
  --region "ALL" \
  --taxa "All" \
  --outfile Extracted_ITS_full.fq.gz \
  --log ITSxpress_Full.log \
  --threads 10

## ITSxpress - extract ITS1 region
itsxpress \
  --fastq trimmed_reads.fq.gz \
  --cluster_id 1.0 \
  --single_end \
  --region "ITS1" \
  --taxa "All" \
  --outfile Extracted_ITS1.fq.gz \
  --log ITSxpress_ITS1.log \
  --threads 10

## ITSxpress - extract ITS2 region
itsxpress \
  --fastq trimmed_reads.fq.gz \
  --cluster_id 1.0 \
  --single_end \
  --region "ITS2" \
  --taxa "All" \
  --outfile Extracted_ITS2.fq.gz \
  --log ITSxpress_ITS2.log \
  --threads 10

