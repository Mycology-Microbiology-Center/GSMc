## OTU clustering

NCORES=8

mkdir -p /tmp/DEREP/

## Remove sequences (after ITSx) shorter than:
export MINLEN_FULL=250    # full-length ITS

# Function for quality filtering and dereplication, short sequence removal
filter_and_derep () {

  # $1 = input file (fastq) or stdin
  # $2 = min length

  # samp_name="$(basename $1 .fq.gz)"

  vsearch \
    --fastq_filter "$1" \
    --fastq_minlen "$2" \
    --fastq_maxee 2 \
    --fastq_maxns 1 \
    --fastq_qmax 93 \
    --fastq_ascii 33 \
    --relabel_sha1 \
    --threads 1 \
    --fastaout - \
  | vsearch \
    --derep_fulllength - \
    --output - \
    --fasta_width 0 \
    --threads 1 \
    --sizein --sizeout
}
export -f filter_and_derep


## Full ITS, non-chimeric
zcat *_NoChimera.fq.gz \
  | filter_and_derep - "$MINLEN_FULL" \
  | gzip -2 > /tmp/DEREP/GSMc_Full_NoChim.fa.gz

## Full ITS putative de novo chimera
zcat *_ChimDeNov.fq.gz \
  | filter_and_derep - "$MINLEN_FULL" \
  | gzip -2 > /tmp/DEREP/GSMc_ChimDeNov.fa.gz

## Function to sort by abundance and sequence length
sort_seqs (){
  seqkit fx2tab --length "$1" \
    | sed -r 's:\t+:\t:g' | sed 's/\t$//g' \
    | awk '{x=$1; gsub("^[[:alnum:]]+;size=", "", x); print $0 "\t" x}' \
    | csvtk sort -H -t -k 4:nr -k 3:nr \
    | awk '{print $1 "\t" $2 }' \
    | seqkit tab2fx -w 0 \
    | gzip -6 --rsyncable
}

sort_seqs /tmp/DEREP/GSMc_Full_NoChim.fa.gz > GSMc_Full_NoChim.fa.gz
sort_seqs /tmp/DEREP/GSMc_ChimDeNov.fa.gz   > GSMc_ChimDeNov.fa.gz

rm -r /tmp/DEREP/


## Concatenate sequences for OTU clustering

# 1) trimmed, high-quality UNITE-INSDc sequences (UNITE_Full.fa.gz)
# 2) trimmed, high-quality GSMc sequences (GSMc_Full_NoChim.fa.gz)
# 2) untrimmed, potentially partial sequences (UNITE_Part.fa.gz and GSMc_Part_NoChim.fa.gz)
# 3) GSMc sequences marked as putatively de novo chimeras (GSMc_ChimDeNov.fa.gz)

cat \
 UNITE_Full.fa.gz \
 GSMc_Full_NoChim.fa.gz \
 UNITE_Part.fa.gz \
 GSMc_Part_NoChim.fa.gz \
 GSMc_ChimDeNov.fa.gz \
 > GL_for_clustering.fasta.gz


## OTU clustering on HPC
vsearch \
  --cluster_smallmem GL_for_clustering.fasta.gz \
  --id 0.98 \
  --iddef 2 \
  --sizein --sizeout \
  --strand both \
  --usersort \
  --relabel_sha1 \
  --threads 40 \
  --centroids GL_OTUs_centr.fasta \
  --consout GL_OTUs_cons.fasta \
  --uc GL_OTUs.uc


## Compress results
gzip -6 --rsyncable GL_OTUs_centr.fasta &
gzip -6 --rsyncable GL_OTUs_cons.fasta &
gzip -6 --rsyncable GL_OTUs.uc



############################################
############################################ Mapping with VSEARCH
############################################

mkdir -p /tmp/DEREP_Samp/

export NCORES=12
export MINLEN=40

## Dereplicate at sample level, add sample ID to the header
derep_rename () {
  samp_name="$(basename $1 .fq.gz)"

  seqkit seq --min-len "$MINLEN" -w 0 "$1" \
  | vsearch \
    --derep_fulllength - \
    --output - \
    --fasta_width 0 \
    --threads 1 \
    --relabel_sha1 \
    --sizein --sizeout \
  | sed 's/>.*/&;sample='"$samp_name"';/' \
  | gzip -4 > /tmp/DEREP_Samp/"$samp_name"_derep.fa.gz

}
export -f derep_rename

find . -type f -name "*.fq.gz" | parallel -j "$NCORES" --progress "derep_rename {}"


## Concatenate samples for mapping (single file)
cat /tmp/DEREP_Samp/*.fa.gz > PacBio_Derep_ForMapping.fa.gz

rm -r /tmp/DEREP_Samp


## Map reads to OTUs on HPC
vsearch \
  --usearch_global PacBio_Derep_ForMapping.fa.gz \
  --db "GL_OTUs_centr.fasta.gz" \
  --id 0.98 \
  --strand both \
  --qmask none \
  --dbmask none \
  --sizein --sizeout \
  --fasta_width 0 \
  --otutabout "OTU_table.txt" \
  --threads 40

## To improve mapping speed, it could be performed in parallel for groups of samples.
## The obtained tables should be merged then into a single OTU table.
