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
