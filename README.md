## The Global Soil Mycobiome consortium dataset

This repository contains the code associated with the paper:

Tedersoo L, Mikryukov V, Anslan S, Bahram M, Khalid AN, Corrales A, Agan A, Vasco-Palacios AM, Saitta A, Antonelli A, Rinaldi AC, Verbeken A, Sulistyo BP, Tamgnoue B, Furneaux B, Duarte Ritter C, Nyamukondiwa C, Sharp C, Marín C, Dai DQ, Gohar D, Sharmah D, Biersma EM, Cameron EK, De Crop E, Otsing E, Davydov EA, Albornoz FA, Brearley FQ, Buegger F, Gates G, Zahn G, Bonito G, Hiiesalu I, Hiiesalu I, Zettur I, Barrio IC, Pärn J, Heilmann-Clausen J, Ankuda J, Kupagme JY, Sarapuu J, Maciá-Vicente JG, Fovo JD, Geml J, Alatalo JM, Alvarez-Manjarrez J, Monkai J, Põldmaa K, Runnel K, Adamson K, Bråthen KA, Pritsch K, Tchan KI, Armolaitis K, Hyde KD, Newsham KK, Panksep K, Adebola LA, Lamit LJ, Saba M, da Silva Cáceres ME, Tuomi M, Gryzenhout M, Bauters M, Bálint M, Wijayawardene N, Hagh-Doust N, Yorou NS, Kurina O, Mortimer PE, Meidl P, Nilsson RH, Puusepp R, Casique-Valdés R, Drenkhan R, Garibay-Orijel R, Godoy R, Alfarraj S, Rahimlou S, Põlme S, Dudov SV, Mundra S, Ahmed T, Netherway T, Henkel TW, Roslin T, Fedosov VE, Onipchenko VG, Yasanthika WAE, Lim YW, Piepenbring M, Klavina D, Kõljalg U, and Abarenkov K. **The Global Soil Mycobiome consortium dataset for boosting fungal diversity research** // *Fungal Diversity* 111, 573-588 (2021) [DOI:10.1007/s13225-021-00493-7](https://link.springer.com/article/10.1007/s13225-021-00493-7).


## Contents


<table>
  <tr>
   <td>
<ul>

<li>01.Demultiplex.sh
</li>
</ul>
   </td>
   <td>Functions used to demultiplex PacBio sequencing runs
   </td>
  </tr>
  <tr>
   <td>
<ul>

<li>02.Extract_ITS.sh
</li>
</ul>
   </td>
   <td>Functions used for trimming the primers and extracting ITS region
   </td>
  </tr>
  <tr>
   <td>
<ul>

<li>03.Chimera_removal.sh
</li>
</ul>
   </td>
   <td>Sample-wise reference-based and <em>de novo</em> chimera removal
   </td>
  </tr>
  <tr>
   <td>
<ul>

<li>04.Prepare_UNITE_data.sh
</li>
</ul>
   </td>
   <td>Preparation of UNITE+INSDc data for clustering
   </td>
  </tr>
  <tr>
   <td>
<ul>

<li>05.Clustering.sh
</li>
</ul>
   </td>
   <td>OTU clustering and sequence mapping
   </td>
  </tr>
  <tr>
   <td>
<ul>

<li>06.OTU_representative_script.R

<li>06.Select_new_OTU_representative.sh
</li>
</ul>
   </td>
   <td>Scripts for selecting alternative representative sequences
   </td>
  </tr>
  <tr>
   <td>
<ul>

<li>07.BLAST.sh
</li>
</ul>
   </td>
   <td>Functions used for taxonomic annotation
   </td>
  </tr>
</table>



## Data availability

The results of the analysis are available from the PlutoF data repository [DOI 10.15156/BIO/2263453](https://doi.org/10.15156/BIO/2263453) and include an OTU table with corresponding sample and OTU (taxonomic and functional) metadata in spreadsheet and Biological Observation Matrix (BIOM) formats.

UNITE 9.01 beta dataset (used for reference-based chimera identification) available at https://doi.org/10.15156/BIO/1444285

Sequence database used for BLAST-based identification at the kingdom level: https://doi.org/10.15156/BIO/1444347

## Dependencies

These scripts require a shell/Linux computing environment and R version 4.0.5.

The following software was used in the analysis:

[VSEARCH](https://github.com/torognes/vsearch/) v.2.17.0 (Rognes et al., 2016)

[cutadapt](https://github.com/marcelm/cutadapt/) v.3.4 (Martin 2011)

[seqkit](https://github.com/shenwei356/seqkit/) v.0.16.0 (Shen et al., 2016)

[ITSxpress](https://github.com/USDA-ARS-GBRU/itsxpress/) v.1.8.0 (Rivers et al. 2018)

[ripgrep](https://github.com/BurntSushi/ripgrep) v.12.1.1

[rush](https://github.com/shenwei356/rush) v.0.4.2

[LIMA](https://github.com/pacificbiosciences/barcoding/) v.2.0.0 (PacBio)

[csvtk](https://github.com/shenwei356/csvtk) v.0.23.0

[BLAST+](https://github.com/ncbi/blast_plus_docs) v.2.11.0 (Camacho et al. 2009)

[GNU bash](https://www.gnu.org/software/bash/) v.5.0.17

[GNU parallel](https://www.gnu.org/software/parallel/) v.20210422 (Tange, 2021)

[GNU awk](https://www.gnu.org/software/gawk/) v.5.1.0

[GNU find](https://www.gnu.org/software/findutils/) v.4.7.0

[GNU sed](https://www.gnu.org/software/sed/) v.4.8

[R](https://cran.r-project.org/) v.4.0.5 (R Core Team, 2021)

[data.table](https://github.com/Rdatatable/data.table) v.1.14.0 (Dowle and Srinivasan, 2021)

[Biostrings](https://github.com/Bioconductor/Biostrings) v.2.6.0 (Pagès et al., 2021)

[plyr](https://github.com/hadley/plyr) v.1.8.6 (Wickham, 2011)


## Acknowledgements

The bulk of this work was supported by the Estonian Science Foundation (grants PRG632, PSG136, MOBTP198, PUT1170), Norway-Baltic financial mechanism (grant EMP442) and Novo Nordisk Fonden (Silva Nova).


## Contact

For all code-related questions, please file a [GitHub issue](https://github.com/Mycology-Microbiology-Center/GSMc/issues).

Please email leho.tedersoo@ut.ee or vladimir.mikryukov@ut.ee for any additional questions about the analytical methods used in this paper. All other relevant data are available from the authors upon request.
