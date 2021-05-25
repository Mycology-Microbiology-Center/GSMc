## The Global Soil Mycobiome consortium dataset

This repository contains the code associated with the paper:

Leho Tedersoo, Vladimir Mikryukov, Sten Anslan, Mohammad Bahram, Abdul Nasir Khalid, Adriana Corrales, Ahto Agan, Aída-Marcela Vasco-Palacios, Alessandro Saitta, Alexandre Antonelli, Andrea C. Rinaldi, Annemieke Verbeken, Bobby P. Sulistyo, Boris Tamgnoue, Brendan Furneaux, Camila Duarte Ritter, Casper Nyamukondiwa, Cathy Sharp, César Marín, D Q. Dai, Daniyal Gohar, Dipon Sharmah, Elisabeth Machteld Biersma, Erin K. Cameron, Eske De Crop, Eveli Otsing, Evgeny A. Davydov, Felipe E. Albornoz, Francis Q. Brearley, Franz Buegger, Genevieve Gates, Geoffrey Zahn, Gregory Bonito, Indrek Hiiesalu, Inga Hiiesalu, Irma Zettur, Isabel C. Barrio, Jaan Pärn, Jacob Heilmann-Clausen, Jelena Ankuda, John Y. Kupagme, Joosep Sarapuu, Jose G. Maciá-Vicente, Joseph Djeugap Fovo, József Geml, Juha M. Alatalo, Julieta Alvarez-Manjarrez, Jutamart Monkai, Kadri Põldmaa, Kadri Runnel, Kalev Adamson, Kari A. Bråthen, Karin Pritsch, Kassim I. Tchan, Kęstutis Armolaitis, Kevin D. Hyde, Kevin K. Newsham, Kristel Panksep, Lateef A. Adebola, Louis J. Lamit, Malka Saba, Marcela E. da Silva Cáceres, Maria Tuomi, Marieka Gryzenhout, Marijn Bauters, Miklós Bálint, Nalin Wijayawardene, Niloufar Hagh-Doust, Nourou S. Yorou, Olavi Kurina, Peter E. Mortimer, Peter Meidl, R. Henrik Nilsson, Rasmus Puusepp, Rebeca Casique-Valdés, Rein Drenkhan, Roberto Garibay-Orijel, Roberto Godoy, Saleh Alfarraj, Saleh Rahimlou, Sergei Põlme, Sergey V. Dudov, Sunil Mundra, Talaat Ahmed, Tarquin Netherway, Terry W. Henkel, Tomas Roslin, Vladimir E. Fedosov, Vladimir G. Onipchenko, W. A. Erandi Yasanthika, Young Woon Lim, Urmas Kõljalg, Kessy Abarenkov. The Global Soil Mycobiome consortium dataset for boosting fungal diversity research.

The manuscript is currently under review.


## Contents


<table>
  <tr>
   <td>
<ul>

<li>01.Demultiplex.txt
</li>
</ul>
   </td>
   <td>Functions used to demultiplex PacBio sequencing runs
   </td>
  </tr>
  <tr>
   <td>
<ul>

<li>02.Extract_ITS.txt
</li>
</ul>
   </td>
   <td>Functions used for trimming the primers and extracting ITS region
   </td>
  </tr>
  <tr>
   <td>
<ul>

<li>03.Chimera_removal.txt
</li>
</ul>
   </td>
   <td>Sample-wise reference-based and <em>de novo</em> chimera removal
   </td>
  </tr>
  <tr>
   <td>
<ul>

<li>04.Prepare_UNITE_data.txt
</li>
</ul>
   </td>
   <td>Preparation of UNITE+INSDc data for clustering
   </td>
  </tr>
  <tr>
   <td>
<ul>

<li>05.Clustering.txt
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

<li>06.Select_new_OTU_representative.txt
</li>
</ul>
   </td>
   <td>Scripts for selecting alternative representative sequences
   </td>
  </tr>
  <tr>
   <td>
<ul>

<li>07.BLAST.txt
</li>
</ul>
   </td>
   <td>Functions used for taxonomic annotation
   </td>
  </tr>
</table>



## Data availability

The results of the analysis are available from the PlutoF data repository DOI and include an OTU table with corresponding sample and OTU (taxonomic and functional) metadata in spreadsheet and Biological Observation Matrix (BIOM) formats.


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

This bulk of this work was supported by the Estonian Science Foundation (grants PRG632, MOBTP198) and Norway-Baltic financial mechanism (grant EMP442).


## Contact

For all code-related questions, please file a [GitHub issue](https://github.com/Mycology-Microbiology-Center/GSMc/issues).

Please email leho.tedersoo@ut.ee or vladimir.mikryukov@ut.ee for any additional questions about the analytical methods used in this paper. All other relevant data are available from the authors upon request.