# OfftargetsPop
Potential off-targets for sgRNAs can be found when an in silico

In order to verify how many other potential off-targets for sgRNAs can be found when an in silico 
analysis is conducted not only against the reference human genome, but considering variants from 
specific populations. For this, we used Crispritz, a tool that can be used for this pourpose and 
also can also be adjusted to predict off-target effects considering a certain number of mismatch 
or RNA/DNA bulges in these potential loci compared to the sgRNA.

We considered variants from:

* (1) the whole population sampled in gnomad v3 (71,702 genomes); 
* (2) the individual populations represented in gnomad v3; 
* (3) the Brazilian population sampled in AbraOM (1200 genomes)

Note: We have run all of our analysis in hg38 annotation 

## Requeriments

* CRISPRitz: https://github.com/pinellolab/CRISPRitz
* R Bioconduct Packages: stringr, dplyr, ggplot2, VennDiagram and RColorBrewer
* Storage local: 100 GB

## Instalation

CRISPRitz Installation and Usage from https://github.com/InfOmics/CRISPRitz

```bash
# crispritz only works with python 3.8
conda install python=3.8 
conda uninstall -c bioconda crispritz
conda install crispritz
```

Install vcflib from https://github.com/vcflib/vcflib

## Run the Scripts

In order to execute the process in backround, write the script in .sh file and then use the command:

```bash
sh Gene_therapy_script.sh > loginfo.out &
```

Analysis results off-targets GT
Note: Before running, change the line below in the .R file.

````
setwd("C:/Users/Karina Griesi/Hospital Albert Einstein/PROADI-AF - General/Anemia_Falciforme/Anemia_Falciforme_compartilhado_20200917/TRI NIO 18-20/off-targets/GT_analysis")
````

```bash
R --file=GT_offtargets_final_analysis2.r
```

## Author's contributions 

D.C.T and KGO designed the study, conducted the bioinformatics analysis, produced
the figures and tables and contributed to writing the manuscript. RW, KTM and ALS
designed the study and contributed to writing the manuscript. All authors read and 
approved the final manuscript.

## Acknowledgements

We acknowledge Murilo de Castro Cervato and Renato David Puga
all the support for bioinformatic analysis and troubleshooting and Rodrigo
Ferreira de Carvalho for the technical support related to IT infrastrucutre.
