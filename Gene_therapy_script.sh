######Gene Therapy 2 - All pops########

#In order to execute the process in backround, write the script in .sh file and then use the command:
#bash execute.sh > loginfo.out &

#####Description######
#In order to verify how many other potential off-targets for sgRNAs can be found when an in silico 
#analysis is conducted not only against the reference human genome, but considering variants from 
#specific populations. For this, we used Crispritz, a tool that can be used for this pourpose and 
#also can also be adjusted to predict off-target effects considering a certain number of mismatch 
#or RNA/DNA bulges in these potential loci compared to the sgRNA.
#We considered variants from: (1) the whole population sampled in gnomad v3 (71,702 genomes); 
#(2) the individual populations represented in gnomad v3; 
#(3) the Brazilian population sampled in AbraOM (1200 genomes)
#We have run all of our analysis in hg38 annotation 


#####Instalation#####
#CRISPRitz Installation and Usage from https://github.com/InfOmics/CRISPRitz
#conda install python=3.8 #crispritz only works with python 3.8
#conda uninstall -c bioconda crispritz
#conda install crispritz
#Install vcflib from https://github.com/vcflib/vcflib


######Files to be downloaded and necessary processing######
#create the following directories:
mkdir hg38_ref
mkdir abraom
mkdir vcf_abraom_filtered
mkdir gnomadv3
mkdir vcf_all_filtered
mkdir vcf_afr_filtered
mkdir vcf_amr_filtered
mkdir vcf_asj_filtered
mkdir vcf_eas_filtered
mkdir vcf_fin_filtered
mkdir vcf_nfe_filtered
mkdir vcf_sas_filtered
mkdir pam
mkdir annotation

######Reference chromosomes#####
#Download reference chromosomes from: 
#http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/
#download the fasta file of each chromosome (individually) and store them in a same folder (hg_38/). The file must be labeled as “chrN.fa” in order to be correctly identified by Crispritz
#Decompressing fasta files
cd hg38_ref
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
gunzip chr$value.fa.gz 
done
cd ..

######Vcf files with population variants######
#download vcf files from gnomad v.3 from: https://gnomad.broadinstitute.org/downloads
#store the files in gnomadv3/
cd gnomadv3
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
wget
https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.chr$value.vcf.bgz
mv gnomad.genomes.r3.0.sites.chr$value.vcf.bgz gnomad.genomes.r3.0.sites.chr$value.vcf.gz
done

#Decompressing gnomad vcf files
for value in 15 16 17 18 19 20 21 22 X Y
do
gunzip -c gnomad.genomes.r3.0.sites.chr$value.vcf.gz > gnomad.genomes.r3.0.sites.chr$value.vcf
done
cd ..

#download vcf files from AbraOM (this was obtained directly from the authors; store at abraom/)
cd abraom/
mv abraom_genomes.txt abraom_genomes.vcf

#the vcf files must be separated by chromosome and stored all in a same folder (separate gnomad and abraom vcfs in distinct folders). Also the vcf files must contain the label “.chrN.” in order to be correctly identified by Crispritz


######Filtering vcf file based on variant frequency >=0.01#######
#filtering gnomad vcf for entire population
cd gnomadv3/
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
vcffilter -o -f "AF > 0.01" -f "AF = 0.01" gnomadgenomes.r3.0.sites.chr$value.vcf > /data/vcf_all_filtered/gnomadgenomes.chr$value.r3.0.filtered01.vcf
done


#filtering gnomad vcf for african population
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
vcffilter -o -f "AF_afr > 0.01" -f "AF_afr = 0.01" gnomad.genomes.r3.0.sites.chr$value.vcf > /data/vcf_afr_filtered/gnomadgenomes.chr$value.r3.0.filtered01_afr.vcf
done


#filtering gnomad vcf for latin population (amr)
for value in 1 2 3 4 5 6 
do
vcffilter -o -f "AF_amr > 0.01" -f "AF_amr = 0.01" gnomad.genomes.r3.0.sites.chr$value.vcf > /data/vcf_amr_filtered/gnomadgenomes.chr$value.r3.0.filtered01_amr.vcf
done


#filtering gnomad vcf for Ashkenazi Jewish population
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
vcffilter -o -f "AF_asj > 0.01" -f "AF_asj = 0.01" gnomad.genomes.r3.0.sites.chr$value.vcf > /data/vcf_asj_filtered/gnomadgenomes.chr$value.r3.0.filtered01_asj.vcf
done


#filtering gnomad vcf for East Asian population
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
vcffilter -o -f "AF_eas > 0.01" -f "AF_eas = 0.01" gnomad.genomes.r3.0.sites.chr$value.vcf > /data/vcf_eas_filtered/ gnomadgenomes.chr$value.r3.0.filtered01_eas.vcf
done


#filtering gnomad vcf for Finnish population
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
vcffilter -o -f "AF_fin > 0.01" -f "AF_fin = 0.01" gnomad.genomes.r3.0.sites.chr$value.vcf > /data/vcf_fin_filtered/ gnomadgenomes.chr$value.r3.0.filtered01_fin.vcf
done

#filtering gnomad vcf for non-Finnish European population
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
vcffilter -o -f "AF_nfe > 0.01" -f "AF_nfe = 0.01" gnomad.genomes.r3.0.sites.chr$value.vcf > /data/vcf_nfe_filtered/ gnomadgenomes.chr$value.r3.0.filtered01_nfe.vcf
done


#filtering gnomad vcf for South Asian population
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
vcffilter -o -f "AF_sas > 0.01" -f "AF_sas = 0.01" gnomad.genomes.r3.0.sites.chr$value.vcf > /data/vcf_sas_filtered/ gnomadgenomes.chr$value.r3.0.filtered01_sas.vcf
done


#Remove decompressed vcfs
rm *.vcf
cd ..

#filtering ABraOM vcf (this vcf is not originally separated by chromosomes; so lets filter the 
#whole genome first and then separate the file into individual chromosomes)
cd abraom/
awk -F "\t" 'BEGIN {OFS="\t"}; { if($8 >= 0.01) { print }}' abraom_genomes.vcf | sed '1d'   > abraom_genomes_int1.vcf
awk -F "\t" 'BEGIN {OFS="\t"}; {print $1,$2,$3,$4,$5,$6,".", "."}' abraom_genomes_int1.vcf > abraom_genomes_filtrado_0.01.vcf
#separate chromosomes
awk 'BEGIN {OFS="\t"}; { if($1 == "chr1") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}'> ../vcf_abraom_filtered/abraom.chr1.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr2") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr2>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}'> ../vcf_abraom_filtered/abraom.chr2.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr3") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr3>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr3.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr4") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr4>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr4.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr5") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr5>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr5.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr6") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr6>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr6.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr7") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr7>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr7.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr8") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr8>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr8.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr9") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr9>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr9.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr10") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr10>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr10.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr11") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr11>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr11.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr12") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr12>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr12.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr13") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr13>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr13.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr14") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr14>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr14.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr15") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr15>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr15.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr16") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr16>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr16.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr17") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr17>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr17.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr18") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr18>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr18.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr19") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr19>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr19.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr20") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr20>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr20.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr21") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr21>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr21.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chr22") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chr22>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chr22.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chrX") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chrX>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chrX.filtrado_0.01.vcf
awk 'BEGIN {OFS="\t"}; { if($1 == "chrY") { print }}' abraom_genomes_filtrado_0.01.vcf | awk 'BEGIN {print "##fileformat=VCFv4.3\n##contig=<ID=chrY>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print}' > ../vcf_abraom_filtered/abraom.chrY.filtrado_0.01.vcf

#bgzip abraom.chr1.filtrado_0.01.vcf
#tabix -p vcf abraom.chr1.filtrado_0.01.vcf.gz
#vcf-validator abraom.chr1.filtrado_0.01.vcf.gz

cd ..

######Compress all filtered vcf files##########
#an update in crispritz requirews now that the vcf files are compressed

#compressing gnomad vcf for entire population
cd vcf_all_filtered
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
gzip gnomadgenomes.chr$value.r3.0.filtered01.vcf
done
cd ..

#compressing gnomad vcf for african population
cd vcf_afr_filtered
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
gzip gnomadgenomes.chr$value.r3.0.filtered01_afr.vcf
done
cd ..

#compressing gnomad vcf for latin population (amr)
cd vcf_amr_filtered
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
gzip gnomadgenomes.chr$value.r3.0.filtered01_amr.vcf
done
cd ..

#compressing gnomad vcf for Ashkenazi Jewish population
cd vcf_asj_filtered
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
gzip gnomadgenomes.chr$value.r3.0.filtered01_asj.vcf
done
cd ..

#compressing gnomad vcf for East Asian population
cd vcf_eas_filtered
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
gzip gnomadgenomes.chr$value.r3.0.filtered01_eas.vcf
done
cd ..

#compressing gnomad vcf for Finnish population
cd vcf_fin_filtered
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
gzip gnomadgenomes.chr$value.r3.0.filtered01_fin.vcf
done
cd ..

#compressing gnomad vcf for non-Finnish European population
cd vcf_nfe_filtered
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
gzip gnomadgenomes.chr$value.r3.0.filtered01_nfe.vcf
done
cd ..

#compressing gnomad vcf for South Asian population
cd vcf_sas_filtered
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
gzip gnomadgenomes.chr$value.r3.0.filtered01_sas.vcf
done
cd ..

#compressing abraom vcf
cd vcf_abraom_filtered
for value in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
gzip abraom.chr$.filtrado_0.01.vcf
done
cd ..


######Install annovar for annotation#####
!wget https://github.com/Varstation/T1-2020/raw/master/annovar/annovar.zip

!unzip annovar.zip
!rm annovar.zip

#download hg38 build version
perl annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene annovar/humandb/

#download gnomad_v3 
perl annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad30_genome annovar/humandb/
#this need to be unziped
!unzip hg38_gnomad30_genome.txt.zip

######Annotation file (previous tests)#######
#Download annotation files from https://genome.ucsc.edu/cgi-bin/hgTables; use GENCODE v32
#check if the assembly is correctly assigned to GRCh38
#for exons
#group: Genes and gene predictions; track: Gencode v32; table: known gene; region: genome; output format: BED; file type returned: plain text
#click in get output
#select exons plus 0 bases at each end
#click on get bed

#for introns
#group: Genes and gene predictions; track: Gencode v32; table: known gene; region: genome; output format: BED; file type returned: plain text
#click in get output
#select exons plus 100 bases at each end
#click on get bed

#for promoters
#group: Genes and gene predictions; track: Gencode v32; table: known gene; region: genome; output format: BED; file type returned: plain text
#click in get output
#select upstream by 1000 bases
#click on get bed

#for DNase clusters
#group: Regulation; track: DNase Clusters; table: wgEncodeRegDnaseClustered; region: genome; output format: BED; file type returned: plain text
#click in get output
#select whole gene
#click on get bed

#for TF clusters
#group: Regulation; track: TF Clusters; table: encRegTfbsClustered; region: genome; output format: BED; file type returned: plain text
#click in get output
#select whole gene
#click on get bed

#Note: The annotation of the sequences of the first 100 base pairs of the introns alone cannot be 
#generated in UCSC. Thus, we generate an annotation file that contains the exon sequence +/- 100bp. 
#In this file, we annotate these sequences as introns, while in the exons bed file, the sequences 
#will be annotated as exons. In this way, off-targets in the exons will appear repeated twice in the 
#table, once annotated as an exon and once as an intron. This second annotation as an intron will 
#drop out when we remove the repeated lines (check if the exon annotation always comes first than 
#the intron annotation and if it is the one that is maintained, in fact).

#To generate a unique file containing all the bed annotation files:
#remove the extra columns
cd annotation
#Compiling annotation terms
awk -F"\t" -v OFS="\t" '{ print $1,$2,$3,"exon" }' exons_hg38.bed > temp_exons_hg38.bed
awk -F"\t" -v OFS="\t" '{ print $1,$2,$3,"TF" }' TF_clusters_hg38.bed > temp_TF_clusters_hg38.bed
awk -F"\t" -v OFS="\t" '{ print $1,$2,$3,"DNase" }' DNase_clusters_hg38.bed > temp_DNase_clusters_hg38.bed
awk -F"\t" -v OFS="\t" '{ print $1,$2,$3,"promoter" }' promoters_hg38.bed > temp_promoters_hg38.bed
awk -F"\t" -v OFS="\t" '{ print $1,$2,$3,"exon100" }' exons100bp_hg38.bed > temp_exon100.bed
awk -F"\t" -v OFS="\t" '{ print $1,$2,$3,"intron_full" }' introns_hg38.bed > temp_introns.bed
awk -F"\t" -v OFS="\t" '{ print $1,$2,$3,"5UTR" }' 5UTR_hg38.bed > temp_5UTR.bed
awk -F"\t" -v OFS="\t" '{ print $1,$2,$3,"3UTR" }' 3UTR_hg38.bed > temp_3UTR.bed
cat temp_* | sort -k1,1 -k2,2n > hg38_annotation_full.bed
rm temp_*
cd ..

######Running Crispritz#######
######Step 1 – Add-variants#######
#Add variants to the reference genome fasta file, creating an enriched fasta genome 
#input: directory with reference chromosome fasta files; directory with vcf files
#output: directory containing the enriched fasta files, one folder for SNP enriched chromosomes and one with indel enriched chromosomes
#model: crispritz.py add-variants vcf_directory/ hg38_ref_directory/


#gnomad all
crispritz.py add-variants vcf_all_filtered/ hg38_ref/
mv variants_genome variants_genome_gnomad_all

#gnomad afr
crispritz.py add-variants vcf_afr_filtered/ hg38_ref/ 
mv variants_genome variants_genome_gnomad_afr

#gnomad amr
crispritz.py add-variants vcf_amr_filtered/ hg38_ref/ 
mv variants_genome variants_genome_gnomad_amr

#gnomad asj
crispritz.py add-variants vcf_asj_filtered/ hg38_ref/ 
mv variants_genome variants_genome_gnomad_asj

#gnomad eas
crispritz.py add-variants vcf_eas_filtered/ hg38_ref/ 
mv variants_genome variants_genome_gnomad_eas

#gnomad fin
crispritz.py add-variants vcf_fin_filtered/ hg38_ref/ 
mv variants_genome variants_genome_gnomad_fin

#gnomad nfe
crispritz.py add-variants vcf_nfe_filtered/ hg38_ref/ 
mv variants_genome variants_genome_gnomad_nfe

#gnomad sas
crispritz.py add-variants vcf_sas_filtered/ hg38_ref/ 
mv variants_genome variants_genome_gnomad_sas

#abraom
crispritz.py add-variants vcf_abraom_filtered/ hg38_ref/
mv variants_genome variants_genome_abraom


#######Step 2 - Index########
#This tool is created to generate an index genome (similar to the bwa-index step). This step is time 
#consuming (from 30 to 60 minutes) but helps to save a lot of execution time while searching with lot 
#of guides and with the support of bulges (DNA and RNA).
#input: name of the genome to be created; directory containing the enriched fasta files; a text file containing the pam and a space separated number indicating the length of the PAM sequence (e.g. Cas9 PAM is NNNNNNNNNNNNNNNNNNNNNGG 3). The sequence is composed by 20 Ns and NGG, followed by number 3, representing the length of the PAM sequence; Number of bulges to include in the database to perform the following search (i.e. the max number bulges allowed for DNA and RNA when searching on the database)
#output: Directory containing an index genome in .bin format, separated into single chromosome files, containing all the candidate targets for a selected PAM, adding also characters to perform bulge search
#model: crispritz.py index-genome name_index_genome_toBeCreated/ enriched_fastas_directory/ pam_directory/pam_file.txt -bMax 2


#gnomad all
crispritz.py index-genome index_SNPs_gnomad_all variants_genome_gnomad_all/SNPs_genome/hg38_ref_enriched pam/pam.txt -bMax 2

#gnomad afr
crispritz.py index-genome index_SNPs_gnomad_afr variants_genome_gnomad_afr/SNPs_genome/hg38_ref_enriched pam/pam.txt -bMax 2

#gnomad amr
crispritz.py index-genome index_SNPs_gnomad_amr variants_genome_gnomad_amr/SNPs_genome/hg38_ref_enriched pam/pam.txt -bMax 2

#gnomad asj
crispritz.py index-genome index_SNPs_gnomad_asj variants_genome_gnomad_asj/SNPs_genome/hg38_ref_enriched pam/pam.txt -bMax 2

#gnomad eas
crispritz.py index-genome index_SNPs_gnomad_eas variants_genome_gnomad_eas/SNPs_genome/hg38_ref_enriched pam/pam.txt -bMax 2

#gnomad fin
crispritz.py index-genome index_SNPs_gnomad_fin variants_genome_gnomad_fin/SNPs_genome/hg38_ref_enriched pam/pam.txt -bMax 2

#gnomad nfe
crispritz.py index-genome index_SNPs_gnomad_nfe variants_genome_gnomad_nfe/SNPs_genome/hg38_ref_enriched pam/pam.txt -bMax 2

#gnomad sas
crispritz.py index-genome index_SNPs_gnomad_sas variants_genome_gnomad_sas/SNPs_genome/hg38_ref_enriched pam/pam.txt -bMax 2

#abraom
crispritz.py index-genome index_SNPs_abraom variants_genome_abraom/SNPs_genome/hg38_ref_enriched pam/pam.txt -bMax 2
 

#ref
crispritz.py index-genome index_SNPs_ref hg38_ref/ pam/pam.txt -bMax 2


#####Step 3 - Search##### 
#Performs the search for off-targets considering mismatches and bulges, in a certain number as 
#established by the user
#input: directory containing the index genome; text file containing the pam sequence; text file containing the guides; name of output file; tag to activate index search; number of mismatches; size of DNA and/or RNA bulges; output type (-r off-targets list only, -p profile only, -t everything); command to activate score generation (-scores followed by the directory of the fasta genome) – gives a score probability of those mismatches occur in that off-target locus
#output: targets file, containing all genomic targets for the guides set; profile file, containing a matrix-like representation of guides behavior (bp/mm, total on-/off- target, targets per mismatch threshold); extended profile file, containing the motif matrix for each guide and each mismatch threshold, useful to create visual analysis of the guides behaviour
#model: crispritz.py search index_genome_directory/ pam_directory/ pam_file.txt guides_directory/ guides_file.txt name_output_file -index -mm 4 -bDNA 2 -bRNA 2 -th 4 -var -t 


#gnomad all:
crispritz.py search genome_library/NGG_2_index_SNPs_gnomad_all/ pam/pam.txt guides/GeneTherapy.txt gnomad_all_SNPs_search_GT -index -mm 4 -bDNA 1 -bRNA 2 -th 4 -var -t 

#gnomad afr:
crispritz.py search genome_library/NGG_2_index_SNPs_gnomad_afr/ pam/pam.txt guides/GeneTherapy.txt gnomad_afr_SNPs_search_GT -index -mm 4 -bDNA 1 -bRNA 2 -th 4 -var -t 

#gnomad amr
crispritz.py search genome_library/NGG_2_index_SNPs_gnomad_amr/ pam/pam.txt guides/GeneTherapy.txt gnomad_amr_SNPs_search_GT -index -mm 4 -bDNA 1 -bRNA 2 -th 4 -var -t 

#gnomad asj
crispritz.py search genome_library/NGG_2_index_SNPs_gnomad_asj/ pam/pam.txt guides/GeneTherapy.txt gnomad_asj_SNPs_search_GT -index -mm 4 -bDNA 1 -bRNA 2 -th 4 -var -t 

#gnomad eas
crispritz.py search genome_library/NGG_2_index_SNPs_gnomad_eas/ pam/pam.txt guides/GeneTherapy.txt gnomad_eas_SNPs_search_GT -index -mm 4 -bDNA 1 -bRNA 2 -th 4 -var -t 

#gnomad fin
crispritz.py search genome_library/NGG_2_index_SNPs_gnomad_fin/ pam/pam.txt guides/GeneTherapy.txt gnomad_fin_SNPs_search_GT -index -mm 4 -bDNA 1 -bRNA 2 -th 4 -var -t 

#gnomad nfe
crispritz.py search genome_library/NGG_2_index_SNPs_gnomad_nfe/ pam/pam.txt guides/GeneTherapy.txt gnomad_nfe_SNPs_search_GT -index -mm 4 -bDNA 1 -bRNA 2 -th 4 -var -t 

#gnomad sas
crispritz.py search genome_library/NGG_2_index_SNPs_gnomad_sas/ pam/pam.txt guides/GeneTherapy.txt gnomad_sas_SNPs_search_GT -index -mm 4 -bDNA 1 -bRNA 2 -th 4 -var -t 

#abraom
crispritz.py search genome_library/NGG_2_index_SNPs_abraom/ pam/pam.txt guides/GeneTherapy.txt abraom_SNPs_search_GT -index -mm 4 -bDNA 1 -bRNA 2 -th 4 -var -t 

#reference
crispritz.py search genome_library/NGG_2_index_SNPs_ref/ pam/pam.txt guides/GeneTherapy.txt ref_search_GT -mm 4 -bDNA 1 -bRNA 2 -th 4 -var -t -scores hg38_ref/


######Step 4 – Filter tables for targets with 3 or less variants and generate scores#######
#If you have too many variants in the same target, as the score function will calculate the CFD 
#for each of combination, if you have a string of 20nt with many variants, the function will try 
#to generate all the possible combinations. To solve this, let's filter targets with 3 variants at most


#gnomad_all SNPs
awk -F"\t" -v fld=3 '{print $0"\t"gsub(/Y/,"",$fld)}' gnomad_all_SNPs_search_GT.targets.txt  | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/R/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/S/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/W/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/K/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/M/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/B/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/D/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/H/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/V/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/y/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/r/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/s/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/w/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/k/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/m/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/b/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/d/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/h/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/v/,"",$fld)}' > count_variants.txt


awk '$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30<4' count_variants.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > gnomad_all_SNPs_GT_filteredVar.targets.txt

crispritz.py scores gnomad_all_SNPs_GT_filteredVar.targets.txt pam/pam.txt guides/GeneTherapy.txt variants_genome_gnomad_all/SNPs_genome/hg38_ref_enriched


#gnomad_afr
awk -F"\t" -v fld=3 '{print $0"\t"gsub(/Y/,"",$fld)}' gnomad_afr_SNPs_search_GT.targets.txt  | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/R/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/S/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/W/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/K/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/M/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/B/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/D/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/H/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/V/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/y/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/r/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/s/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/w/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/k/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/m/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/b/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/d/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/h/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/v/,"",$fld)}' > count_variants_afr.txt


awk '$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30<4' count_variants_afr.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > gnomad_afr_SNPs_GT_filteredVar.targets.txt

crispritz.py scores gnomad_afr_SNPs_GT_filteredVar.targets.txt pam/pam.txt guides/GeneTherapy.txt variants_genome_gnomad_afr/SNPs_genome/hg38_ref_enriched

#gnomad amr
awk -F"\t" -v fld=3 '{print $0"\t"gsub(/Y/,"",$fld)}' gnomad_amr_SNPs_search_GT.targets.txt  | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/R/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/S/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/W/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/K/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/M/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/B/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/D/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/H/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/V/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/y/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/r/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/s/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/w/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/k/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/m/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/b/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/d/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/h/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/v/,"",$fld)}' > count_variants_amr.txt


awk '$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30<4' count_variants_amr.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > gnomad_amr_SNPs_GT_filteredVar.targets.txt

crispritz.py scores gnomad_amr_SNPs_GT_filteredVar.targets.txt pam/pam.txt guides/GeneTherapy.txt variants_genome_gnomad_amr/SNPs_genome/hg38_ref_enriched

#gnomad asj
awk -F"\t" -v fld=3 '{print $0"\t"gsub(/Y/,"",$fld)}' gnomad_asj_SNPs_search_GT.targets.txt  | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/R/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/S/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/W/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/K/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/M/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/B/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/D/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/H/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/V/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/y/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/r/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/s/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/w/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/k/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/m/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/b/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/d/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/h/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/v/,"",$fld)}' > count_variants_asj.txt


awk '$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30<4' count_variants_asj.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > gnomad_asj_SNPs_GT_filteredVar.targets.txt

crispritz.py scores gnomad_asj_SNPs_GT_filteredVar.targets.txt pam/pam.txt guides/GeneTherapy.txt variants_genome_gnomad_asj/SNPs_genome/hg38_ref_enriched

#gnomad eas
awk -F"\t" -v fld=3 '{print $0"\t"gsub(/Y/,"",$fld)}' gnomad_eas_SNPs_search_GT.targets.txt  | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/R/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/S/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/W/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/K/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/M/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/B/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/D/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/H/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/V/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/y/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/r/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/s/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/w/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/k/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/m/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/b/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/d/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/h/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/v/,"",$fld)}' > count_variants_eas.txt


awk '$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30<4' count_variants_eas.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > gnomad_eas_SNPs_GT_filteredVar.targets.txt

crispritz.py scores gnomad_eas_SNPs_GT_filteredVar.targets.txt pam/pam.txt guides/GeneTherapy.txt variants_genome_gnomad_eas/SNPs_genome/hg38_ref_enriched

#gnomad fin 
awk -F"\t" -v fld=3 '{print $0"\t"gsub(/Y/,"",$fld)}' gnomad_fin_SNPs_search_GT.targets.txt  | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/R/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/S/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/W/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/K/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/M/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/B/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/D/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/H/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/V/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/y/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/r/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/s/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/w/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/k/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/m/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/b/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/d/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/h/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/v/,"",$fld)}' > count_variants_fin.txt


awk '$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30<4' count_variants_fin.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > gnomad_fin_SNPs_GT_filteredVar.targets.txt

crispritz.py scores gnomad_fin_SNPs_GT_filteredVar.targets.txt pam/pam.txt guides/GeneTherapy.txt variants_genome_gnomad_fin/SNPs_genome/hg38_ref_enriched

#gnomad nfe
awk -F"\t" -v fld=3 '{print $0"\t"gsub(/Y/,"",$fld)}' gnomad_nfe_SNPs_search_GT.targets.txt  | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/R/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/S/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/W/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/K/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/M/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/B/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/D/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/H/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/V/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/y/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/r/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/s/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/w/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/k/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/m/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/b/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/d/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/h/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/v/,"",$fld)}' > count_variants_nfe.txt


awk '$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30<4' count_variants_nfe.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > gnomad_nfe_SNPs_GT_filteredVar.targets.txt

crispritz.py scores gnomad_nfe_SNPs_GT_filteredVar.targets.txt pam/pam.txt guides/GeneTherapy.txt variants_genome_gnomad_nfe/SNPs_genome/hg38_ref_enriched

#gnomad sas
awk -F"\t" -v fld=3 '{print $0"\t"gsub(/Y/,"",$fld)}' gnomad_sas_SNPs_search_GT.targets.txt  | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/R/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/S/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/W/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/K/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/M/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/B/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/D/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/H/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/V/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/y/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/r/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/s/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/w/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/k/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/m/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/b/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/d/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/h/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/v/,"",$fld)}' > count_variants_sas.txt


awk '$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30<4' count_variants_sas.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > gnomad_sas_SNPs_GT_filteredVar.targets.txt

crispritz.py scores gnomad_sas_SNPs_GT_filteredVar.targets.txt pam/pam.txt guides/GeneTherapy.txt variants_genome_gnomad_sas/SNPs_genome/hg38_ref_enriched


#abraom
awk -F"\t" -v fld=3 '{print $0"\t"gsub(/Y/,"",$fld)}' abraom_SNPs_search_GT.targets.txt  | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/R/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/S/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/W/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/K/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/M/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/B/,"",$fld)}'| awk -F"\t" -v fld=3 '{print $0"\t"gsub(/D/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/H/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/V/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/y/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/r/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/s/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/w/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/k/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/m/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/b/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/d/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/h/,"",$fld)}' | awk -F"\t" -v fld=3 '{print $0"\t"gsub(/v/,"",$fld)}' > count_variants_abraom.txt


awk '$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30<4' count_variants_abraom.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > abraom_SNPs_GT_filteredVar.targets.txt

crispritz.py scores abraom_SNPs_GT_filteredVar.targets.txt pam/pam.txt guides/GeneTherapy.txt variants_genome_abraom/SNPs_genome/hg38_ref_enriched


######Step 5 - Filter the targets of interest#######
#Take the targets only with bulges
#gnomad_all SNPs
awk -F "\t" 'BEGIN {OFS="\t"}; {if($8 == 0) { print }}' gnomad_all_GT_SNPs_filteredVar.targets.CFD.txt > gnomad_all_GT_SNPs_bulges.txt 

awk -F "\t" 'BEGIN {OFS="\t"}; {if($10 < 5) { print }}' gnomad_all_GT_SNPs_filteredVar.targets.CFD.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {if($11 >= 0.2 ) { print }}' > gnomad_all_GT_SNPs_mismatches.txt 
cat gnomad_all_GT_SNPs_bulges.txt gnomad_all_GT_SNPs_mismatches.txt > gnomad_all_GT_SNPs_FinalTargets.txt


#gnomad_afr SNPs
awk -F "\t" 'BEGIN {OFS="\t"}; {if($8 == 0) { print }}' gnomad_afr_GT_SNPs_filteredVar.targets.CFD.txt > gnomad_afr_GT_SNPs_bulges.txt

awk -F "\t" 'BEGIN {OFS="\t"}; {if($10 < 5) { print }}' gnomad_afr_GT_SNPs_filteredVar.targets.CFD.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {if($11 >= 0.2 ) { print }}' > gnomad_afr_GT_SNPs_mismatches.txt 
cat gnomad_afr_GT_SNPs_bulges.txt gnomad_afr_GT_SNPs_mismatches.txt > gnomad_afr_GT_SNPs_FinalTargets.txt

#gnomad amr
awk -F "\t" 'BEGIN {OFS="\t"}; {if($8 == 0) { print }}' gnomad_amr_GT_SNPs_filteredVar.targets.CFD.txt > gnomad_amr_GT_SNPs_bulges.txt

awk -F "\t" 'BEGIN {OFS="\t"}; {if($10 < 5) { print }}' gnomad_amr_GT_SNPs_filteredVar.targets.CFD.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {if($11 >= 0.2 ) { print }}' > gnomad_amr_GT_SNPs_mismatches.txt 
cat gnomad_amr_GT_SNPs_bulges.txt gnomad_amr_GT_SNPs_mismatches.txt > gnomad_amr_GT_SNPs_FinalTargets.txt

#gnomad asj
awk -F "\t" 'BEGIN {OFS="\t"}; {if($8 == 0) { print }}' gnomad_asj_GT_SNPs_filteredVar.targets.CFD.txt > gnomad_asj_GT_SNPs_bulges.txt

awk -F "\t" 'BEGIN {OFS="\t"}; {if($10 < 5) { print }}' gnomad_asj_GT_SNPs_filteredVar.targets.CFD.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {if($11 >= 0.2 ) { print }}' > gnomad_asj_GT_SNPs_mismatches.txt 
cat gnomad_asj_GT_SNPs_bulges.txt gnomad_asj_GT_SNPs_mismatches.txt > gnomad_asj_GT_SNPs_FinalTargets.txt

#gnomad eas
awk -F "\t" 'BEGIN {OFS="\t"}; {if($8 == 0) { print }}' gnomad_eas_GT_SNPs_filteredVar.targets.CFD.txt > gnomad_eas_GT_SNPs_bulges.txt

awk -F "\t" 'BEGIN {OFS="\t"}; {if($10 < 5) { print }}' gnomad_eas_GT_SNPs_filteredVar.targets.CFD.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {if($11 >= 0.2 ) { print }}' > gnomad_eas_GT_SNPs_mismatches.txt 
cat gnomad_eas_GT_SNPs_bulges.txt gnomad_eas_GT_SNPs_mismatches.txt > gnomad_eas_GT_SNPs_FinalTargets.txt

#gnomad fin
awk -F "\t" 'BEGIN {OFS="\t"}; {if($8 == 0) { print }}' gnomad_fin_GT_SNPs_filteredVar.targets.CFD.txt > gnomad_fin_GT_SNPs_bulges.txt

awk -F "\t" 'BEGIN {OFS="\t"}; {if($10 < 5) { print }}' gnomad_fin_GT_SNPs_filteredVar.targets.CFD.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {if($11 >= 0.2 ) { print }}' > gnomad_fin_GT_SNPs_mismatches.txt 
cat gnomad_fin_GT_SNPs_bulges.txt gnomad_fin_GT_SNPs_mismatches.txt > gnomad_fin_GT_SNPs_FinalTargets.txt

#gnomad nfe
awk -F "\t" 'BEGIN {OFS="\t"}; {if($8 == 0) { print }}' gnomad_nfe_GT_SNPs_filteredVar.targets.CFD.txt > gnomad_nfe_GT_SNPs_bulges.txt

awk -F "\t" 'BEGIN {OFS="\t"}; {if($10 < 5) { print }}' gnomad_nfe_GT_SNPs_filteredVar.targets.CFD.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {if($11 >= 0.2 ) { print }}' > gnomad_nfe_GT_SNPs_mismatches.txt 
cat gnomad_nfe_GT_SNPs_bulges.txt gnomad_nfe_GT_SNPs_mismatches.txt > gnomad_nfe_GT_SNPs_FinalTargets.txt

#gnomad sas
awk -F "\t" 'BEGIN {OFS="\t"}; {if($8 == 0) { print }}' gnomad_sas_GT_SNPs_filteredVar.targets.CFD.txt > gnomad_sas_GT_SNPs_bulges.txt

awk -F "\t" 'BEGIN {OFS="\t"}; {if($10 < 5) { print }}' gnomad_sas_GT_SNPs_filteredVar.targets.CFD.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {if($11 >= 0.2 ) { print }}' > gnomad_sas_GT_SNPs_mismatches.txt 
cat gnomad_sas_GT_SNPs_bulges.txt gnomad_sas_GT_SNPs_mismatches.txt > gnomad_sas_GT_SNPs_FinalTargets.txt


#abraom SNPs
awk -F "\t" 'BEGIN {OFS="\t"}; {if($8 == 0) { print }}' abraom_GT_SNPs_filteredVar.targets.CFD.txt > abraom_GT_SNPs_bulges.txt

awk -F "\t" 'BEGIN {OFS="\t"}; {if($10 < 5) { print }}' abraom_GT_SNPs_filteredVar.targets.CFD.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {if($11 >= 0.2 ) { print }}' > abraom_GT_SNPs_mismatches.txt 
cat abraom_GT_SNPs_bulges.txt abraom_GT_SNPs_mismatches.txt > abraom_GT_SNPs_FinalTargets.txt


#reference
awk -F "\t" 'BEGIN {OFS="\t"}; {if($8 == 0) { print }}' ref_search_GT.targets.CFD.txt > reference_GT_bulges.txt

awk -F "\t" 'BEGIN {OFS="\t"}; {if($10 < 5) { print }}' ref_search_GT.targets.CFD.txt | awk -F "\t" 'BEGIN {OFS="\t"}; {if($11 >= 0.2 ) { print }}' > reference_GT_mismatches.txt 
cat reference_GT_bulges.txt reference_GT_mismatches.txt > reference_GT_FinalTargets.txt

#Step 6 – annotate off-targets - previous tests
#input: Targets file, containing all genomic targets for the guides set; Bed file containing the annotations; Name of output file; Samples ID file, containing the list of samples with their associated Population and Superpopulation (Optional)
#output: Targets file with annotation (identical file as the targets file in input) with an added column containing the annotations); One summary file, counting all the annotations per mismatch number.
#model: crispritz.py annotate-results targetFile_from_previous_step.txt annotations.bed name_output_file.annotated


#gnomad all
crispritz.py annotate-results gnomad_all_GT_SNPs_FinalTargets.txt hg38_annotation_full.bed gnomad_all_GT_SNPs_FinalTargets.annotated

crispritz.py annotate-results gnomad_all_GT_INDELs_FinalTargets.txt hg38_annotation_full.bed gnomad_all_GT_INDELs_FinalTargets.annotated

#gnomad afr
crispritz.py annotate-results gnomad_afr_GT_SNPs_FinalTargets.txt hg38_annotation_full.bed gnomad_afr_GT_SNPs_FinalTargets.annotated.txt
crispritz.py annotate-results gnomad_afr_GT_INDELs_FinalTargets.txt hg38_annotation_full.bed gnomad_afr_GT_INDELs_FinalTargets.annotated

#abraom
crispritz.py annotate-results abraom_GT_SNPs_FinalTargets.txt hg38_annotation_full.bed abraom_GT_SNPs_FinalTargets.annotated
crispritz.py annotate-results abraom_GT_INDELs_FinalTargets.txt hg38_annotation_full.bed abraom_GT_INDELs_FinalTargets.annotated

#reference
crispritz.py annotate-results ref_GT_SNPs_FinalTargets.txt hg38_annotation_full.bed ref_GT_SNPs_FinalTargets.annotated


#from here, the script was written in R (file: GT_offtarget_final_analysis.R)

