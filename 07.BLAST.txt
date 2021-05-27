## Taxonomy annotation of sequences

## MegaBLAST
blastn \
  -task megablast \
  -outfmt=6 \
  -strand both \
  -query input.fasta \
  -db "/path/to/BLAST/database" \
  -max_target_seqs 10 \
  -max_hsps 1 \
  -out megablast_output.txt \
  -num_threads 8


## BLASTn
blastn \
  -task blastn \
  -outfmt=6 \
  -strand both \
  -query input.fasta \
  -db "/path/to/BLAST/database" \
  -max_target_seqs 10 \
  -max_hsps 1 \
  -out blastn_output.txt \
  -num_threads 8

