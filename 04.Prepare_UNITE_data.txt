## Prepare UNITE+INSDc data for clustering
## - ITS extraction
## - sorting sequences into full-length and partial sequences
## - prefix dereplication


#########################################################
######################################################### Extract ITS
#########################################################

NCORES=8

ITSx \
  -i UNITE.fasta \
  --complement T \
  --save_regions all \
  --graphical F \
  --positions T \
  -E 1e-1 \
  -t all \
  --cpu "$NCORES" \
  --preserve T \
  -o ITSX



### Combine partial sequences

# Workflow:
#    1. full ITS
#  + 2. _ITS1 + 5.8S + ITS2_   [partial, not landed in 'full']
#  + 3.  ITS1 + 5.8S
#  + 4.         5.8S + ITS2
#  + 5.  ITS1
#  + 6.                ITS2
#    7. non_detections


## Prepare tables for ID matching
seqkit fx2tab *.full.fasta | sed 's/\t$//g' | csvtk add-header -t -n id,Full > tmp_0_FullITS.txt
seqkit fx2tab *.ITS1.fasta | sed 's/\t$//g' | csvtk add-header -t -n id,ITS1 > tmp_1_ITS1.txt
seqkit fx2tab *.5_8S.fasta | sed 's/\t$//g' | csvtk add-header -t -n id,58S  > tmp_1_s58.txt
seqkit fx2tab *.ITS2.fasta | sed 's/\t$//g' | csvtk add-header -t -n id,ITS2 > tmp_1_ITS2.txt

## Join ITS fragments
csvtk join -t -f "id" tmp_1_ITS1.txt tmp_1_s58.txt tmp_1_ITS2.txt   > tmp_2_ITS1_58S_ITS2.txt
csvtk join -t -f "id" tmp_1_ITS1.txt tmp_1_s58.txt                  > tmp_3_ITS1_58S.txt
csvtk join -t -f "id"                  tmp_1_s58.txt tmp_1_ITS2.txt > tmp_4_58S_ITS2.txt


## Remove duplicated sequences
R

## Load data
s_0_FullITS       <- read.delim("tmp_0_FullITS.txt")
s_1_ITS1          <- read.delim("tmp_1_ITS1.txt")
s_1_s58           <- read.delim("tmp_1_s58.txt")
s_1_ITS2          <- read.delim("tmp_1_ITS2.txt")
s_2_ITS1_58S_ITS2 <- read.delim("tmp_2_ITS1_58S_ITS2.txt")
s_3_ITS1_58S      <- read.delim("tmp_3_ITS1_58S.txt")
s_4_58S_ITS2      <- read.delim("tmp_4_58S_ITS2.txt")

## Concatenate ITS parts
concat_its <- function(z){
  clz <- c("ITS1", "X58S", "ITS2") %in% colnames(s_2_ITS1_58S_ITS2)[-1]
  if(sum(clz) == 3){                                           # _ITS1 + 5.8S + ITS2_
    z$Full <- paste(z$ITS1, z$X58S, z$ITS2, sep = "")
  } else {
    if(clz[1] == TRUE & clz[2] == TRUE & clz[3] == FALSE){     # ITS1 + 5.8S
      z$Full <- paste(z$ITS1, z$X58S, sep = "")
    }
    if(clz[1] == FALSE & clz[2] == TRUE & clz[3] == TRUE){     # 5.8S + ITS2
      z$Full <- paste(z$X58S, z$ITS2, sep = "")
    }
  }

  z <- z[, c("id", "Full")]
  return(z)
}

## Progressively add sequences
RES <- data.frame(s_0_FullITS, Type = "0_Full_length")

tmp <- concat_its(s_2_ITS1_58S_ITS2)   # _ITS1 + 5.8S + ITS2_
tmp <- tmp[!tmp$id %in% RES$id, ]
tmp <- data.frame(tmp, Type = "1_ITS1_5.8S_ITS2")
RES <- rbind(RES, tmp)
rm(tmp)

tmp <- concat_its(s_3_ITS1_58S)        # ITS1 + 5.8S
tmp <- tmp[!tmp$id %in% RES$id, ]
tmp <- data.frame(tmp, Type = "2_ITS1_5.8S")
RES <- rbind(RES, tmp)
rm(tmp)

tmp <- concat_its(s_3_ITS1_58S)        # 5.8S + ITS2
tmp <- tmp[!tmp$id %in% RES$id, ]
tmp <- data.frame(tmp, Type = "3_5.8S_ITS2")
RES <- rbind(RES, tmp)
rm(tmp)

tmp <- s_1_ITS1                        # ITS1
colnames(tmp)[2] <- "Full"
tmp <- tmp[!tmp$id %in% RES$id, ]
tmp <- data.frame(tmp, Type = "4_ITS1")
RES <- rbind(RES, tmp)
rm(tmp)

tmp <- s_1_ITS2                        # ITS2
colnames(tmp)[2] <- "Full"
tmp <- tmp[!tmp$id %in% RES$id, ]
tmp <- data.frame(tmp, Type = "5_ITS2")
RES <- rbind(RES, tmp)
rm(tmp)

write.table(x = RES, file = "tmp_Res.txt", quote = F, sep = "\t", row.names = F, col.names = F)

q("no")


## Convert table back to fasta,
## Remove leading and trailing Ns
awk '{ print $1 "\t" $2 }' tmp_Res.txt \
  | seqkit tab2fx -w 0  \
  | seqkit replace -p "^n+|n+$" -r "" -is -w 0 \
  | gzip -6 > tmp_Res.fasta.gz

## Clean up
rm Res_UNITE_uniq.part_*.5_8S.fasta
rm Res_UNITE_uniq.part_*.SSU.fasta
rm Res_UNITE_uniq.part_*.LSU.fasta
rm Res_UNITE_uniq.part_*.ITS2.fasta
rm Res_UNITE_uniq.part_*.ITS1.fasta
rm Res_UNITE_uniq.part_*.full.fasta
rm Res_UNITE_uniq.part_*.chimeric.fasta
rm Res_UNITE_uniq.part_*_no_detections.fasta

rm tmp_4_58S_ITS2.txt
rm tmp_3_ITS1_58S.txt
rm tmp_2_ITS1_58S_ITS2.txt
rm tmp_1_s58.txt
rm tmp_1_ITS2.txt
rm tmp_1_ITS1.txt
rm tmp_0_FullITS.txt



#########################################################
######################################################### Prefix dereplication and sorting
#########################################################

## Prefix dereplication
vsearch \
  --derep_prefix tmp_Res.fasta.gz \
  --sizein --sizeout \
  --output - \
  --uc UNITE_PrefDerep.uc \
  | gzip -6 > UNITE_PrefDerep.fasta.gz


## Sort sequences by number of ambiguous nucleotides and length
seqkit fx2tab UNITE_PrefDerep.fasta.gz \
  --length --base-content ACTG --base-content NRYSWKMBDHV \
  | sed -r 's:\t+:\t:g' | sed 's/\t$//g' \
  | sed 's/;size=/\t/g' \
  | csvtk add-header -t -n id,size,seq,len,ACTG,N \
  | csvtk sort -t -k N:n -k len:nr \
  > UNITE_PrefDerep_table.txt


## Then:
## - Remove short sequences and sequences with large number of ambiguous nucleotides
## - Prepare data for sequence mapping (`UNITE_for_mapping.fasta`)
##   each sequences should have `size` and `sample` annotations in USEARCH format

