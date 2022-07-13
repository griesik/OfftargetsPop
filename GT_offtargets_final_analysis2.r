#Analysis results off-targets GT
setwd("C:/Users/Karina Griesi/Hospital Albert Einstein/PROADI-AF - General/Anemia_Falciforme/Anemia_Falciforme_compartilhado_20200917/TRIÊNIO 18-20/off-targets/GT_analysis")
load("GT_offtargets_final_analysis_env2.RData")
library(stringr)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)

#####load sgRNAs table####
sgRNAs=read.delim("sgRNAs.txt", header=FALSE)
sgRNAs$V1=as.factor(sgRNAs$V1)
colnames(sgRNAs)=c("crRNANoBulges","gene")

#####managing the tables to get unique variants and to prepare vcf files for annotation#####
#####abraom####
#read file
abraom = read.delim ("abraom_SNPs_GT_FinalTargets.txt", header= FALSE)

#copy crRNA column in a new column to don´t loose bulge info and 
#remove "-"
abraom [,12]= gsub("-","",abraom[,2])
#generate a new column, concatenating crRNA, chromosome and cluster position
abraom [,13]= paste (abraom[,12],"_",abraom[,4],":",abraom[,6], sep="")
colnames(abraom)=c("Bulge_type","crRNA", "DNA", "Chromosome", "Position", "Cluster_Position", 
                   "Direction", "Mismatches", "Bulge_Size", "Total", "CFD", "crRNANoBulges", "crRNAPos")
#sort dataframe by this last created column and by CFD in descending order
abraom_unique = abraom[order(abraom$crRNAPos,-abraom$CFD),]
#take only unique values
abraom_unique = distinct(abraom_unique, crRNAPos, .keep_all = TRUE)
#create a column only with chromosome numbers
abraom_unique[,14]=gsub("chr","",abraom_unique[,4])
#create a column with the first base of the crRNA
abraom_unique[,15]=paste(str_sub(abraom_unique$crRNA, start=1, end=1))
#create a column with fake variants
abraom_unique[,16]="NA"
  
for (j in c(1:nrow(abraom_unique)))
{
  if (abraom_unique$V15[j]=="A") abraom_unique$V16[j]="G"
  if (abraom_unique$V15[j]=="C") abraom_unique$V16[j]="T"
  if (abraom_unique$V15[j]=="G") abraom_unique$V16[j]="A"
  if (abraom_unique$V15[j]=="T") abraom_unique$V16[j]="C"
  
}
abraom_unique[,17]="abraom"

vcf_abraom= abraom_unique[,c(14,6)]
vcf_abraom[,2]=as.numeric(vcf_abraom[,2])
vcf_abraom[,3]="."
vcf_abraom[,c(4,5)]=abraom_unique[,c(15,16)]
vcf_abraom[,c(6:9)]="."
vcf_abraom = vcf_abraom[order(vcf_abraom$V14,vcf_abraom$Cluster_Position),]



write.table (vcf_abraom, file="vcf_abraom_toVEP.txt", row.names = FALSE, sep="\t",
             quote = FALSE, col.names = FALSE)

#need to complete this part

#vep_abraom = ensemblVEP(vcf_abraom, param)

#####Gnomad_afr####
#read file
gnomad_afr = read.delim ("gnomad_afr_SNPs_GT_FinalTargets.txt")
#copy crRNA column in a new column to don´t loose bulge info
#remove "-"
gnomad_afr [,12]= gsub("-","",gnomad_afr[,2])
#generate a new column, concatenating crRNA, chromosome and cluster position
gnomad_afr [,13]= paste (gnomad_afr[,12],"_",gnomad_afr[,4],":",gnomad_afr[,6], sep="")
colnames(gnomad_afr)=c("Bulge_type","crRNA", "DNA", "Chromosome", "Position", "Cluster_Position", 
                   "Direction", "Mismatches", "Bulge_Size", "Total", "CFD", "crRNANoBulges", "crRNAPos")

#sort dataframe by this last created column and by CFD in descending order
gnomad_afr_unique = gnomad_afr[order(gnomad_afr$crRNAPos,-gnomad_afr$CFD),]
#take only unique values
gnomad_afr_unique = distinct(gnomad_afr_unique, crRNAPos, .keep_all = TRUE)
#create a column only with chromosome numbers
gnomad_afr_unique[,14]=gsub("chr","",gnomad_afr_unique[,4])
#create a column with the first base of the crRNA
gnomad_afr_unique[,15]=paste(str_sub(gnomad_afr_unique$crRNA, start=1, end=1))
#create a column with fake variants
gnomad_afr_unique[,16]="NA"

for (j in c(1:nrow(gnomad_afr_unique)))
{
  if (gnomad_afr_unique$V15[j]=="A") gnomad_afr_unique$V16[j]="G"
  if (gnomad_afr_unique$V15[j]=="C") gnomad_afr_unique$V16[j]="T"
  if (gnomad_afr_unique$V15[j]=="G") gnomad_afr_unique$V16[j]="A"
  if (gnomad_afr_unique$V15[j]=="T") gnomad_afr_unique$V16[j]="C"
  
}
gnomad_afr_unique[,17]="gnomad_afr"

vcf_gnomad_afr= gnomad_afr_unique[,c(14,6)]
vcf_gnomad_afr[,2]=as.numeric(vcf_gnomad_afr[,2])
vcf_gnomad_afr[,3]="."
vcf_gnomad_afr[,c(4,5)]=gnomad_afr_unique[,c(15,16)]
vcf_gnomad_afr[,c(6:9)]="."
vcf_gnomad_afr = vcf_gnomad_afr[order(vcf_gnomad_afr$V14,vcf_gnomad_afr$Cluster),]

write.table (vcf_gnomad_afr, file="vcf_gnomad_afr_toVEP.txt", row.names = FALSE, sep="\t",
             quote = FALSE, col.names = FALSE)


#####Gnomad_all####
#read file
gnomad_all = read.delim ("gnomad_all_SNPs_GT_FinalTargets.txt")
#copy crRNA column in a new column to don´t loose bulge info
#remove "-"
gnomad_all [,12]= gsub("-","",gnomad_all[,2])
#generate a new column, concatenating crRNA, chromosome and cluster position
gnomad_all [,13]= paste (gnomad_all[,12],"_",gnomad_all[,4],":",gnomad_all[,6], sep="")
colnames(gnomad_all)=c("Bulge_type","crRNA", "DNA", "Chromosome", "Position", "Cluster_Position", 
                       "Direction", "Mismatches", "Bulge_Size", "Total", "CFD", "crRNANoBulges", "crRNAPos")

#sort dataframe by this last created column and by CFD in descending order
gnomad_all_unique = gnomad_all[order(gnomad_all$crRNAPos,-gnomad_all$CFD),]
#take only unique values
gnomad_all_unique = distinct(gnomad_all_unique, crRNAPos, .keep_all = TRUE)
#create a column only with chromosome numbers
gnomad_all_unique[,14]=gsub("chr","",gnomad_all_unique[,4])
#create a column with the first base of the crRNA
gnomad_all_unique[,15]=paste(str_sub(gnomad_all_unique$crRNA, start=1, end=1))
#create a column with fake variants
gnomad_all_unique[,16]="NA"

for (j in c(1:nrow(gnomad_all_unique)))
{
  if (gnomad_all_unique$V15[j]=="A") gnomad_all_unique$V16[j]="G"
  if (gnomad_all_unique$V15[j]=="C") gnomad_all_unique$V16[j]="T"
  if (gnomad_all_unique$V15[j]=="G") gnomad_all_unique$V16[j]="A"
  if (gnomad_all_unique$V15[j]=="T") gnomad_all_unique$V16[j]="C"
  
}
gnomad_all_unique[,17]="gnomad_all"

vcf_gnomad_all= gnomad_all_unique[,c(14,6)]
vcf_gnomad_all[,2]=as.numeric(vcf_gnomad_all[,2])
vcf_gnomad_all[,3]="."
vcf_gnomad_all[,c(4,5)]=gnomad_all_unique[,c(15,16)]
vcf_gnomad_all[,c(6:9)]="."
vcf_gnomad_all = vcf_gnomad_all[order(vcf_gnomad_all$V14,vcf_gnomad_all$Cluster),]

write.table (vcf_gnomad_all, file="vcf_gnomad_all_toVEP.txt", row.names = FALSE, sep="\t",
             quote = FALSE, col.names = FALSE)


#####Gnomad_amr####
#read file
gnomad_amr = read.delim ("gnomad_amr_SNPs_GT_FinalTargets.txt")
#copy crRNA column in a new column to don´t loose bulge info
#remove "-"
gnomad_amr [,12]= gsub("-","",gnomad_amr[,2])
#generate a new column, concatenating crRNA, chromosome and cluster position
gnomad_amr [,13]= paste (gnomad_amr[,12],"_",gnomad_amr[,4],":",gnomad_amr[,6], sep="")
colnames(gnomad_amr)=c("Bulge_type","crRNA", "DNA", "Chromosome", "Position", "Cluster_Position", 
                       "Direction", "Mismatches", "Bulge_Size", "Total", "CFD", "crRNANoBulges", "crRNAPos")

#sort datamrame by this last created column and by CFD in descending order
gnomad_amr_unique = gnomad_amr[order(gnomad_amr$crRNAPos,-gnomad_amr$CFD),]
#take only unique values
gnomad_amr_unique = distinct(gnomad_amr_unique, crRNAPos, .keep_all = TRUE)
#create a column only with chromosome numbers
gnomad_amr_unique[,14]=gsub("chr","",gnomad_amr_unique[,4])
#create a column with the first base of the crRNA
gnomad_amr_unique[,15]=paste(str_sub(gnomad_amr_unique$crRNA, start=1, end=1))
#create a column with fake variants
gnomad_amr_unique[,16]="NA"

for (j in c(1:nrow(gnomad_amr_unique)))
{
  if (gnomad_amr_unique$V15[j]=="A") gnomad_amr_unique$V16[j]="G"
  if (gnomad_amr_unique$V15[j]=="C") gnomad_amr_unique$V16[j]="T"
  if (gnomad_amr_unique$V15[j]=="G") gnomad_amr_unique$V16[j]="A"
  if (gnomad_amr_unique$V15[j]=="T") gnomad_amr_unique$V16[j]="C"
  
}

gnomad_amr_unique[,17]="gnomad_amr"

vcf_gnomad_amr= gnomad_amr_unique[,c(14,6)]
vcf_gnomad_amr[,2]=as.numeric(vcf_gnomad_amr[,2])
vcf_gnomad_amr[,3]="."
vcf_gnomad_amr[,c(4,5)]=gnomad_amr_unique[,c(15,16)]
vcf_gnomad_amr[,c(6:9)]="."
vcf_gnomad_amr = vcf_gnomad_amr[order(vcf_gnomad_amr$V14,vcf_gnomad_amr$Cluster),]

write.table (vcf_gnomad_amr, file="vcf_gnomad_amr_toVEP.txt", row.names = FALSE, sep="\t",
             quote = FALSE, col.names = FALSE)


#####Gnomad_asj####
#read file
gnomad_asj = read.delim ("gnomad_asj_SNPs_GT_FinalTargets.txt")
#copy crRNA column in a new column to don´t loose bulge info
#remove "-"
gnomad_asj [,12]= gsub("-","",gnomad_asj[,2])
#generate a new column, concatenating crRNA, chromosome and cluster position
gnomad_asj [,13]= paste (gnomad_asj[,12],"_",gnomad_asj[,4],":",gnomad_asj[,6], sep="")
colnames(gnomad_asj)=c("Bulge_type","crRNA", "DNA", "Chromosome", "Position", "Cluster_Position", 
                       "Direction", "Mismatches", "Bulge_Size", "Total", "CFD", "crRNANoBulges", "crRNAPos")

#sort dataframe by this last created column and by CFD in descending order
gnomad_asj_unique = gnomad_asj[order(gnomad_asj$crRNAPos,-gnomad_asj$CFD),]
#take only unique values
gnomad_asj_unique = distinct(gnomad_asj_unique, crRNAPos, .keep_all = TRUE)
#create a column only with chromosome numbers
gnomad_asj_unique[,14]=gsub("chr","",gnomad_asj_unique[,4])
#create a column with the first base of the crRNA
gnomad_asj_unique[,15]=paste(str_sub(gnomad_asj_unique$crRNA, start=1, end=1))
#create a column with fake variants
gnomad_asj_unique[,16]="NA"

for (j in c(1:nrow(gnomad_asj_unique)))
{
  if (gnomad_asj_unique$V15[j]=="A") gnomad_asj_unique$V16[j]="G"
  if (gnomad_asj_unique$V15[j]=="C") gnomad_asj_unique$V16[j]="T"
  if (gnomad_asj_unique$V15[j]=="G") gnomad_asj_unique$V16[j]="A"
  if (gnomad_asj_unique$V15[j]=="T") gnomad_asj_unique$V16[j]="C"
  
}

gnomad_asj_unique[,17]="gnomad_asj"

vcf_gnomad_asj= gnomad_asj_unique[,c(14,6)]
vcf_gnomad_asj[,2]=as.numeric(vcf_gnomad_asj[,2])
vcf_gnomad_asj[,3]="."
vcf_gnomad_asj[,c(4,5)]=gnomad_asj_unique[,c(15,16)]
vcf_gnomad_asj[,c(6:9)]="."
vcf_gnomad_asj = vcf_gnomad_asj[order(vcf_gnomad_asj$V14,vcf_gnomad_asj$Cluster),]

write.table (vcf_gnomad_asj, file="vcf_gnomad_asj_toVEP.txt", row.names = FALSE, sep="\t",
             quote = FALSE, col.names = FALSE)



#####Gnomad_eas####
#read file
gnomad_eas = read.delim ("gnomad_eas_SNPs_GT_FinalTargets.txt")
#copy crRNA column in a new column to don´t loose bulge info
#remove "-"
gnomad_eas [,12]= gsub("-","",gnomad_eas[,2])
#generate a new column, concatenating crRNA, chromosome and cluster position
gnomad_eas [,13]= paste (gnomad_eas[,12],"_",gnomad_eas[,4],":",gnomad_eas[,6], sep="")
colnames(gnomad_eas)=c("Bulge_type","crRNA", "DNA", "Chromosome", "Position", "Cluster_Position", 
                       "Direction", "Mismatches", "Bulge_Size", "Total", "CFD", "crRNANoBulges", "crRNAPos")

#sort dateasame by this last created column and by CFD in descending order
gnomad_eas_unique = gnomad_eas[order(gnomad_eas$crRNAPos,-gnomad_eas$CFD),]
#take only unique values
gnomad_eas_unique = distinct(gnomad_eas_unique, crRNAPos, .keep_all = TRUE)
#create a column only with chromosome numbers
gnomad_eas_unique[,14]=gsub("chr","",gnomad_eas_unique[,4])
#create a column with the first base of the crRNA
gnomad_eas_unique[,15]=paste(str_sub(gnomad_eas_unique$crRNA, start=1, end=1))
#create a column with fake variants
gnomad_eas_unique[,16]="NA"

for (j in c(1:nrow(gnomad_eas_unique)))
{
  if (gnomad_eas_unique$V15[j]=="A") gnomad_eas_unique$V16[j]="G"
  if (gnomad_eas_unique$V15[j]=="C") gnomad_eas_unique$V16[j]="T"
  if (gnomad_eas_unique$V15[j]=="G") gnomad_eas_unique$V16[j]="A"
  if (gnomad_eas_unique$V15[j]=="T") gnomad_eas_unique$V16[j]="C"
  
}

gnomad_eas_unique[,17]="gnomad_eas"

vcf_gnomad_eas= gnomad_eas_unique[,c(14,6)]
vcf_gnomad_eas[,2]=as.numeric(vcf_gnomad_eas[,2])
vcf_gnomad_eas[,3]="."
vcf_gnomad_eas[,c(4,5)]=gnomad_eas_unique[,c(15,16)]
vcf_gnomad_eas[,c(6:9)]="."
vcf_gnomad_eas = vcf_gnomad_eas[order(vcf_gnomad_eas$V14,vcf_gnomad_eas$Cluster),]

write.table (vcf_gnomad_eas, file="vcf_gnomad_eas_toVEP.txt", row.names = FALSE, sep="\t",
             quote = FALSE, col.names = FALSE)


#####Gnomad_fin####
#read file
gnomad_fin = read.delim ("gnomad_fin_SNPs_GT_FinalTargets.txt")
#copy crRNA column in a new column to don´t loose bulge info
#remove "-"
gnomad_fin [,12]= gsub("-","",gnomad_fin[,2])
#generate a new column, concatenating crRNA, chromosome and cluster position
gnomad_fin [,13]= paste (gnomad_fin[,12],"_",gnomad_fin[,4],":",gnomad_fin[,6], sep="")
colnames(gnomad_fin)=c("Bulge_type","crRNA", "DNA", "Chromosome", "Position", "Cluster_Position", 
                       "Direction", "Mismatches", "Bulge_Size", "Total", "CFD", "crRNANoBulges", "crRNAPos")

#sort datfiname by this last created column and by CFD in descending order
gnomad_fin_unique = gnomad_fin[order(gnomad_fin$crRNAPos,-gnomad_fin$CFD),]
#take only unique values
gnomad_fin_unique = distinct(gnomad_fin_unique, crRNAPos, .keep_all = TRUE)
#create a column only with chromosome numbers
gnomad_fin_unique[,14]=gsub("chr","",gnomad_fin_unique[,4])
#create a column with the first base of the crRNA
gnomad_fin_unique[,15]=paste(str_sub(gnomad_fin_unique$crRNA, start=1, end=1))
#create a column with fake variants
gnomad_fin_unique[,16]="NA"

for (j in c(1:nrow(gnomad_fin_unique)))
{
  if (gnomad_fin_unique$V15[j]=="A") gnomad_fin_unique$V16[j]="G"
  if (gnomad_fin_unique$V15[j]=="C") gnomad_fin_unique$V16[j]="T"
  if (gnomad_fin_unique$V15[j]=="G") gnomad_fin_unique$V16[j]="A"
  if (gnomad_fin_unique$V15[j]=="T") gnomad_fin_unique$V16[j]="C"
  
}

gnomad_fin_unique[,17]="gnomad_fin"

vcf_gnomad_fin= gnomad_fin_unique[,c(14,6)]
vcf_gnomad_fin[,2]=as.numeric(vcf_gnomad_fin[,2])
vcf_gnomad_fin[,3]="."
vcf_gnomad_fin[,c(4,5)]=gnomad_fin_unique[,c(15,16)]
vcf_gnomad_fin[,c(6:9)]="."
vcf_gnomad_fin = vcf_gnomad_fin[order(vcf_gnomad_fin$V14,vcf_gnomad_fin$Cluster),]

write.table (vcf_gnomad_fin, file="vcf_gnomad_fin_toVEP.txt", row.names = FALSE, sep="\t",
             quote = FALSE, col.names = FALSE)


#####Gnomad_nfe####
#read file
gnomad_nfe = read.delim ("gnomad_nfe_SNPs_GT_FinalTargets.txt")
#copy crRNA column in a new column to don´t loose bulge info
#remove "-"
gnomad_nfe [,12]= gsub("-","",gnomad_nfe[,2])
#generate a new column, concatenating crRNA, chromosome and cluster position
gnomad_nfe [,13]= paste (gnomad_nfe[,12],"_",gnomad_nfe[,4],":",gnomad_nfe[,6], sep="")
colnames(gnomad_nfe)=c("Bulge_type","crRNA", "DNA", "Chromosome", "Position", "Cluster_Position", 
                       "Direction", "Mismatches", "Bulge_Size", "Total", "CFD", "crRNANoBulges", "crRNAPos")

#sort dataframe by this last created column and by CFD in descending order
gnomad_nfe_unique = gnomad_nfe[order(gnomad_nfe$crRNAPos,-gnomad_nfe$CFD),]
#take only unique values
gnomad_nfe_unique = distinct(gnomad_nfe_unique, crRNAPos, .keep_all = TRUE)
#create a column only with chromosome numbers
gnomad_nfe_unique[,14]=gsub("chr","",gnomad_nfe_unique[,4])
#create a column with the first base of the crRNA
gnomad_nfe_unique[,15]=paste(str_sub(gnomad_nfe_unique$crRNA, start=1, end=1))
#create a column with fake variants
gnomad_nfe_unique[,16]="NA"

for (j in c(1:nrow(gnomad_nfe_unique)))
{
  if (gnomad_nfe_unique$V15[j]=="A") gnomad_nfe_unique$V16[j]="G"
  if (gnomad_nfe_unique$V15[j]=="C") gnomad_nfe_unique$V16[j]="T"
  if (gnomad_nfe_unique$V15[j]=="G") gnomad_nfe_unique$V16[j]="A"
  if (gnomad_nfe_unique$V15[j]=="T") gnomad_nfe_unique$V16[j]="C"
  
}

gnomad_nfe_unique[,17]="gnomad_nfe"

vcf_gnomad_nfe= gnomad_nfe_unique[,c(14,6)]
vcf_gnomad_nfe[,2]=as.numeric(vcf_gnomad_nfe[,2])
vcf_gnomad_nfe[,3]="."
vcf_gnomad_nfe[,c(4,5)]=gnomad_nfe_unique[,c(15,16)]
vcf_gnomad_nfe[,c(6:9)]="."
vcf_gnomad_nfe = vcf_gnomad_nfe[order(vcf_gnomad_nfe$V14,vcf_gnomad_nfe$Cluster),]

write.table (vcf_gnomad_nfe, file="vcf_gnomad_nfe_toVEP.txt", row.names = FALSE, sep="\t",
             quote = FALSE, col.names = FALSE)


#####Gnomad_sas####
#read file
gnomad_sas = read.delim ("gnomad_sas_SNPs_GT_FinalTargets.txt")
#copy crRNA column in a new column to don´t loose bulge info
#remove "-"
gnomad_sas [,12]= gsub("-","",gnomad_sas[,2])
#generate a new column, concatenating crRNA, chromosome and cluster position
gnomad_sas [,13]= paste (gnomad_sas[,12],"_",gnomad_sas[,4],":",gnomad_sas[,6], sep="")
colnames(gnomad_sas)=c("Bulge_type","crRNA", "DNA", "Chromosome", "Position", "Cluster_Position", 
                       "Direction", "Mismatches", "Bulge_Size", "Total", "CFD", "crRNANoBulges", "crRNAPos")

#sort dataframe by this last created column and by CFD in descending order
gnomad_sas_unique = gnomad_sas[order(gnomad_sas$crRNAPos,-gnomad_sas$CFD),]
#take only unique values
gnomad_sas_unique = distinct(gnomad_sas_unique, crRNAPos, .keep_all = TRUE)
#create a column only with chromosome numbers
gnomad_sas_unique[,14]=gsub("chr","",gnomad_sas_unique[,4])
#create a column with the first base of the crRNA
gnomad_sas_unique[,15]=paste(str_sub(gnomad_sas_unique$crRNA, start=1, end=1))
#create a column with fake variants
gnomad_sas_unique[,16]="NA"

for (j in c(1:nrow(gnomad_sas_unique)))
{
  if (gnomad_sas_unique$V15[j]=="A") gnomad_sas_unique$V16[j]="G"
  if (gnomad_sas_unique$V15[j]=="C") gnomad_sas_unique$V16[j]="T"
  if (gnomad_sas_unique$V15[j]=="G") gnomad_sas_unique$V16[j]="A"
  if (gnomad_sas_unique$V15[j]=="T") gnomad_sas_unique$V16[j]="C"
  
}

gnomad_sas_unique[,17]="gnomad_sas"

vcf_gnomad_sas= gnomad_sas_unique[,c(14,6)]
vcf_gnomad_sas[,2]=as.numeric(vcf_gnomad_sas[,2])
vcf_gnomad_sas[,3]="."
vcf_gnomad_sas[,c(4,5)]=gnomad_sas_unique[,c(15,16)]
vcf_gnomad_sas[,c(6:9)]="."
vcf_gnomad_sas = vcf_gnomad_sas[order(vcf_gnomad_sas$V14,vcf_gnomad_sas$Cluster),]

write.table (vcf_gnomad_sas, file="vcf_gnomad_sas_toVEP.txt", row.names = FALSE, sep="\t",
             quote = FALSE, col.names = FALSE)


#####Reference####
#read file
ref = read.delim ("reference_GT_FinalTargets.txt")
#copy crRNA column in a new column to don´t loose bulge info
#remove "-"
ref [,12]= gsub("-","",ref[,2])
#generate a new column, concatenating crRNA, chromosome and cluster position
ref [,13]= paste (ref[,12],"_",ref[,4],":",ref[,6], sep="")
colnames(ref)=c("Bulge_type","crRNA", "DNA", "Chromosome", "Position", "Cluster_Position", 
                       "Direction", "Mismatches", "Bulge_Size", "Total", "CFD", "crRNANoBulges", "crRNAPos")

#sort dataframe by this last created column and by CFD in descending order
ref_unique = ref[order(ref$crRNAPos,-ref$CFD),]
#take only unique values
ref_unique = distinct(ref_unique, crRNAPos, .keep_all = TRUE)
#create a column only with chromosome numbers
ref_unique[,14]=gsub("chr","",ref_unique[,4])
#create a column with the first base of the crRNA
ref_unique[,15]=paste(str_sub(ref_unique$crRNA, start=1, end=1))
#create a column with fake variants
ref_unique[,16]="NA"

for (j in c(1:nrow(ref_unique)))
{
  if (ref_unique$V15[j]=="A") ref_unique$V16[j]="G"
  if (ref_unique$V15[j]=="C") ref_unique$V16[j]="T"
  if (ref_unique$V15[j]=="G") ref_unique$V16[j]="A"
  if (ref_unique$V15[j]=="T") ref_unique$V16[j]="C"
  
}

ref_unique[,17]="ref"

vcf_ref= ref_unique[,c(14,6)]
vcf_ref[,2]=as.numeric(vcf_ref[,2])
vcf_ref[,3]="."
vcf_ref[,c(4,5)]=ref_unique[,c(15,16)]
vcf_ref[,c(6:9)]="."
vcf_ref = vcf_ref[order(vcf_ref$V14,vcf_ref$Cluster),]

write.table (vcf_ref, file="vcf_ref_toVEP.txt", row.names = FALSE, sep="\t",
             quote = FALSE, col.names = FALSE)


#####Combine all sets (before annotation)####
totalTableAll=rbind(gnomad_all_unique, gnomad_afr_unique,gnomad_amr_unique,
                 gnomad_asj_unique, gnomad_eas_unique, gnomad_fin_unique,
                 gnomad_nfe_unique, gnomad_sas_unique, ref_unique)

totalTableAll= totalTableAll[order(totalTableAll$crRNAPos,-totalTableAll$CFD),]
totalTableAll[,18]=seq(1:27035)

colnames (totalTableAll)[c(14:18)]=c("Chr","REF","ALTfake","higherCFD","lineTotalTableAll")


totalTable=distinct(totalTableAll,crRNAPos, .keep_all=TRUE)



presence = as.data.frame(gnomad_all_unique$crRNAPos)
presence[,2]="gnomad_all"
colnames(presence)[1]="crRNAPos"

totalTable = left_join(totalTable, presence, by="crRNAPos")

presence = as.data.frame(gnomad_afr_unique$crRNAPos)
presence[,2]="gnomad_afr"
colnames(presence)[1]="crRNAPos"

totalTable = left_join(totalTable, presence, by="crRNAPos")

presence = as.data.frame(gnomad_amr_unique$crRNAPos)
presence[,2]="gnomad_amr"
colnames(presence)[1]="crRNAPos"

totalTable = left_join(totalTable, presence, by="crRNAPos")


presence = as.data.frame(gnomad_asj_unique$crRNAPos)
presence[,2]="gnomad_asj"
colnames(presence)[1]="crRNAPos"

totalTable = left_join(totalTable, presence, by="crRNAPos")

presence = as.data.frame(gnomad_eas_unique$crRNAPos)
presence[,2]="gnomad_eas"
colnames(presence)[1]="crRNAPos"

totalTable = left_join(totalTable, presence, by="crRNAPos")


presence = as.data.frame(gnomad_fin_unique$crRNAPos)
presence[,2]="gnomad_fin"
colnames(presence)[1]="crRNAPos"

totalTable = left_join(totalTable, presence, by="crRNAPos")

presence = as.data.frame(gnomad_nfe_unique$crRNAPos)
presence[,2]="gnomad_nfe"
colnames(presence)[1]="crRNAPos"

totalTable = left_join(totalTable, presence, by="crRNAPos")


presence = as.data.frame(gnomad_sas_unique$crRNAPos)
presence[,2]="gnomad_sas"
colnames(presence)[1]="crRNAPos"

totalTable = left_join(totalTable, presence, by="crRNAPos")

presence = as.data.frame(ref_unique$crRNAPos)
presence[,2]="reference"
colnames(presence)[1]="crRNAPos"

totalTable = left_join(totalTable, presence, by="crRNAPos")


colnames(totalTable)[c(19:27)]=c("gnomad_all","gnomad_afr","gnomad_amr",
                                 "gnomad_asj", "gnomad_eas", "gnomad_fin",
                                 "gnomad_nfe", "gnomad_sas", "reference")


totalTable$gnomad_all[which(totalTable$gnomad_all%in%NA)]="0" 
totalTable$gnomad_afr[which(totalTable$gnomad_afr%in%NA)]="0"
totalTable$gnomad_amr[which(totalTable$gnomad_amr%in%NA)]="0"
totalTable$gnomad_asj[which(totalTable$gnomad_asj%in%NA)]="0"
totalTable$gnomad_eas[which(totalTable$gnomad_eas%in%NA)]="0"
totalTable$gnomad_fin[which(totalTable$gnomad_fin%in%NA)]="0"
totalTable$gnomad_nfe[which(totalTable$gnomad_nfe%in%NA)]="0"
totalTable$gnomad_sas[which(totalTable$gnomad_sas%in%NA)]="0"
totalTable$reference[which(totalTable$reference%in%NA)]="0" 

#two CCR5 sgRNAs and one for EB have less than 20 bases; 
#2 PDCD1 sgRNAs have 21 bases; the results are not precise for these

totalTable=subset(totalTable, totalTable$crRNANoBulges!="GACTATGCTGCCGCCCAGTNNNN" &
                    totalTable$crRNANoBulges!= "GCAGAAGGGGACAGTAAGANNNN" &
                    totalTable$crRNANoBulges!= "GTCCGCAGCTTTCTCGANNNNNN" &
                    totalTable$crRNANoBulges!="GGCAGTTGTGTGACACGGAANNN" &
                    totalTable$crRNANoBulges!="GGCGTGACTTCCACATGAGCNNN")

totalTable=left_join(totalTable, sgRNAs, by="crRNANoBulges")
totalTable=totalTable[,c(1,2,3,28,4:27)]


####Annotate variants#####
#create a file to run in annovar
#this file will be suitable for annotation of the regions

annovar=totalTable[,c(14,5,5,15,16)]
write.table (annovar, file="annovar_input.avinput", row.names = FALSE, 
           sep="\t", col.names= FALSE, quote = FALSE)

#this is the command to be run in perl:
#perl annovar/table_annovar.pl annovar_input.avinput annovar/humandb/ --dot2underline -buildver hg38 -out annot_targets -remove -protocol refGene -operation g -nastring "."

#load annotated table
totalTableAnn=read.delim("annot_targets.hg38_multianno.txt")
totalTableAnn=cbind(totalTable,totalTableAnn[,c(6,7,8)])

CFDhigherNotRef=totalTableAnn[(totalTableAnn$higherCFD !="gnomad_all"),]
CFDhigherNotRef=CFDhigherNotRef[(CFDhigherNotRef$reference=="reference"),]
write.csv(CFDhigherNotRef,file="CFDhigherNotRef.csv", row.names = FALSE)

#totalTableAnn_noRef_exc=subset(totalTableAnn, totalTableAnn$gnomad_all != 0 |
 #                                totalTableAnn$gnomad_afr != 0 |
  #                               totalTableAnn$gnomad_asj != 0 |
   #                              totalTableAnn$gnomad_amr != 0 |
    #                             totalTableAnn$gnomad_eas != 0 |
     #                            totalTableAnn$gnomad_fin != 0 |
      #                           totalTableAnn$gnomad_nfe != 0 |
       #                          totalTableAnn$gnomad_sas != 0)
#


write.csv (totalTableAnn, file = "FinalTargets_concatenated.csv", row.names = FALSE)

######Exclusives####
exc_gnomad_all= subset (totalTableAnn, totalTableAnn$gnomad_all =="gnomad_all" &
                          totalTableAnn$gnomad_afr == "0" &
                          totalTableAnn$gnomad_amr =="0" &
                          totalTableAnn$gnomad_asj =="0" &
                          totalTableAnn$gnomad_eas =="0" &
                          totalTableAnn$gnomad_fin =="0" &
                          totalTableAnn$gnomad_nfe =="0" &
                          totalTableAnn$gnomad_sas =="0" &
                          totalTableAnn$reference == "0")

exc_gnomad_afr= subset (totalTableAnn, totalTableAnn$gnomad_afr =="gnomad_afr" &
                          totalTableAnn$gnomad_all == "0" &
                          totalTableAnn$gnomad_amr =="0" &
                          totalTableAnn$gnomad_asj =="0" &
                          totalTableAnn$gnomad_eas =="0" &
                          totalTableAnn$gnomad_fin =="0" &
                          totalTableAnn$gnomad_nfe =="0" &
                          totalTableAnn$gnomad_sas =="0" &
                          totalTableAnn$reference == "0")

exc_gnomad_amr= subset (totalTableAnn, totalTableAnn$gnomad_amr =="gnomad_amr" &
                          totalTableAnn$gnomad_afr == "0" &
                          totalTableAnn$gnomad_all =="0" &
                          totalTableAnn$gnomad_asj =="0" &
                          totalTableAnn$gnomad_eas =="0" &
                          totalTableAnn$gnomad_fin =="0" &
                          totalTableAnn$gnomad_nfe =="0" &
                          totalTableAnn$gnomad_sas =="0" &
                          totalTableAnn$reference == "0")

exc_gnomad_asj= subset (totalTableAnn, totalTableAnn$gnomad_asj =="gnomad_asj" &
                          totalTableAnn$gnomad_afr == "0" &
                          totalTableAnn$gnomad_amr =="0" &
                          totalTableAnn$gnomad_all =="0" &
                          totalTableAnn$gnomad_eas =="0" &
                          totalTableAnn$gnomad_fin =="0" &
                          totalTableAnn$gnomad_nfe =="0" &
                          totalTableAnn$gnomad_sas =="0" &
                          totalTableAnn$reference == "0")

exc_gnomad_eas= subset (totalTableAnn, totalTableAnn$gnomad_eas =="gnomad_eas" &
                          totalTableAnn$gnomad_afr == "0" &
                          totalTableAnn$gnomad_amr =="0" &
                          totalTableAnn$gnomad_asj =="0" &
                          totalTableAnn$gnomad_all =="0" &
                          totalTableAnn$gnomad_fin =="0" &
                          totalTableAnn$gnomad_nfe =="0" &
                          totalTableAnn$gnomad_sas =="0" &
                          totalTableAnn$reference == "0")

exc_gnomad_fin= subset (totalTableAnn, totalTableAnn$gnomad_fin =="gnomad_fin" &
                          totalTableAnn$gnomad_afr == "0" &
                          totalTableAnn$gnomad_amr =="0" &
                          totalTableAnn$gnomad_asj =="0" &
                          totalTableAnn$gnomad_eas =="0" &
                          totalTableAnn$gnomad_all =="0" &
                          totalTableAnn$gnomad_nfe =="0" &
                          totalTableAnn$gnomad_sas =="0" &
                          totalTableAnn$reference == "0")

exc_gnomad_nfe= subset (totalTableAnn, totalTableAnn$gnomad_nfe =="gnomad_nfe" &
                          totalTableAnn$gnomad_afr == "0" &
                          totalTableAnn$gnomad_amr =="0" &
                          totalTableAnn$gnomad_asj =="0" &
                          totalTableAnn$gnomad_eas =="0" &
                          totalTableAnn$gnomad_fin =="0" &
                          totalTableAnn$gnomad_all =="0" &
                          totalTableAnn$gnomad_sas =="0" &
                          totalTableAnn$reference == "0")

exc_gnomad_sas= subset (totalTableAnn, totalTableAnn$gnomad_sas =="gnomad_sas" &
                          totalTableAnn$gnomad_afr == "0" &
                          totalTableAnn$gnomad_amr =="0" &
                          totalTableAnn$gnomad_asj =="0" &
                          totalTableAnn$gnomad_eas =="0" &
                          totalTableAnn$gnomad_fin =="0" &
                          totalTableAnn$gnomad_nfe =="0" &
                          totalTableAnn$gnomad_all =="0" &
                          totalTableAnn$reference == "0")

exc_ref= subset (totalTableAnn, totalTableAnn$reference =="reference" &
                          totalTableAnn$gnomad_afr == "0" &
                          totalTableAnn$gnomad_amr =="0" &
                          totalTableAnn$gnomad_asj =="0" &
                          totalTableAnn$gnomad_eas =="0" &
                          totalTableAnn$gnomad_fin =="0" &
                          totalTableAnn$gnomad_nfe =="0" &
                          totalTableAnn$gnomad_sas =="0" &
                          totalTableAnn$gnomad_all == "0")

exc_notInReference= subset (totalTableAnn, totalTableAnn$reference =="0")

vcf_exc_notInReference= exc_notInReference[,c(15,7)]
vcf_exc_notInReference[,2]=as.numeric(vcf_exc_notInReference[,2])
vcf_exc_notInReference[,3]="."
vcf_exc_notInReference[,c(4,5)]=exc_notInReference[,c(16,17)]
vcf_exc_notInReference[,c(6:9)]="."
vcf_exc_notInReference= vcf_exc_notInReference[order(vcf_exc_notInReference$Chr,vcf_exc_notInReference$Cluster_Position),]
write.table (vcf_exc_notInReference, file="vcf_excNotInRef_toVEP.txt", row.names = FALSE, sep="\t",
             quote = FALSE, col.names = FALSE)


write.table (exc_gnomad_afr, file = "exclusive_gnomad_afr.txt", row.names = FALSE, 
             sep="\t", quote = FALSE)
write.table (exc_gnomad_asj, file = "exclusive_gnomad_asj.txt", row.names = FALSE, 
             sep="\t", quote = FALSE)
write.table (exc_gnomad_eas, file = "exclusive_gnomad_eas.txt", row.names = FALSE, 
             sep="\t", quote = FALSE)
write.table (exc_gnomad_fin, file = "exclusive_gnomad_fin.txt", row.names = FALSE, 
             sep="\t", quote = FALSE)
write.table (exc_gnomad_nfe, file = "exclusive_gnomad_nfe.txt", row.names = FALSE, 
             sep="\t", quote = FALSE)
write.table (exc_gnomad_sas, file = "exclusive_gnomad_sas.txt", row.names = FALSE, 
             sep="\t", quote = FALSE)
write.table (exc_ref, file = "exclusive_reference.txt", row.names = FALSE, 
             sep="\t", quote = FALSE)








#####Incorporating the frequency information to the final Table####
#create a file only with targets not found in reference to annotate
#the frequency of variants in them
#this is the command to be run in perl for that:
#perl annovar/table_annovar.pl annovar_input_freq.avinput annovar/humandb/ --dot2underline -buildver hg38 -out annot_targets -remove -protocol refGene,gnomad30_genome -operation g,f -nastring "."


#annovarNotInRef = exc_notInReference [,c(15,5,5,15,16)]
annovarNotInRef = exc_notInReference [,c(15,6,6,16,17)]
annovarNotInRef[,6]=seq(1:234)
annovarNotInRef[,7]=str_count(exc_notInReference$crRNA, "-")#how many DNA bulges
annovarNotInRef[,8]=str_count(exc_notInReference$DNA, "-")#how many RNA bulges
annovarNotInRef[,9]=23+(annovarNotInRef[,7]-annovarNotInRef[,8])
annovarNotInRefF=as.data.frame(matrix(nrow=1,ncol=6))
colnames(annovarNotInRefF)=c("V1","V2","V3","V4","V5","V6")

for (j in (1:nrow(annovarNotInRef)))
{
  x=annovarNotInRef[j,2]
  y=annovarNotInRef$V9[j]
  seq=seq(x+1,x+y, by=1)
  newmatrix=as.data.frame(matrix(nrow=y,ncol=5))
  newmatrix$V1=annovarNotInRef[j,1]
  newmatrix$V2=seq
  newmatrix$V3=newmatrix$V2
  newmatrix$V4="0"
  newmatrix$V5="G"
  newmatrix$V6=annovarNotInRef[j,6]
  annovarNotInRefF=rbind(annovarNotInRefF,newmatrix)
}



annovarNotInRefF=annovarNotInRefF[-1,]
annovarNotInRefF=annovarNotInRefF[order(annovarNotInRefF$V1,
                                        annovarNotInRefF$V2),]

annovarNotInRefC=annovarNotInRefF
annovarNotInRefC$V5="C"

annovarNotInRefA=annovarNotInRefF
annovarNotInRefA$V5="A"

annovarNotInRefT=annovarNotInRefF
annovarNotInRefT$V5="T"

annovarNotInRefF2=rbind(annovarNotInRefF,annovarNotInRefA,
                        annovarNotInRefC,annovarNotInRefT)
annovarNotInRefF2=annovarNotInRefF2[order(annovarNotInRefF2$V1,
                                          annovarNotInRefF2$V2,
                                          annovarNotInRefF2$V6),]



write.table (annovarNotInRefF2, file = "annovar_input_freqNew.avinput", row.names = FALSE, 
             sep="\t", quote = FALSE, col.names = FALSE)




#load the table generated by annovar
annovarFreq=read.delim("annot_targets_freq.hg38_multiannoNew.txt")
annovarFreq=cbind(annovarFreq,annovarNotInRefF2$V6)
annovarFreq=subset(annovarFreq, annovarFreq$AF != ".")
annovarFreqCheck=annovarFreq

#remove frequencies in log10 base, which is being interpreted as character, not number
for (j in (1:nrow(annovarFreq)))
  {if (str_detect(annovarFreq$AF[j], "e")==TRUE) annovarFreq$AF[j]=0
  if (str_detect(annovarFreq$AF_afr[j], "e")==TRUE) annovarFreq$AF_afr[j]=0
  if (str_detect(annovarFreq$AF_amr[j], "e")==TRUE) annovarFreq$AF_amr[j]=0
  if (str_detect(annovarFreq$AF_asj[j], "e")==TRUE) annovarFreq$AF_asj[j]=0
  if (str_detect(annovarFreq$AF_eas[j], "e")==TRUE) annovarFreq$AF_eas[j]=0
  if (str_detect(annovarFreq$AF_fin[j], "e")==TRUE) annovarFreq$AF_fin[j]=0
  if (str_detect(annovarFreq$AF_nfe[j], "e")==TRUE) annovarFreq$AF_nfe[j]=0
  if (str_detect(annovarFreq$AF_sas[j], "e")==TRUE) annovarFreq$AF_sas[j]=0}


annovarFreq$AF=as.numeric(annovarFreq$AF)
annovarFreq$AF_afr=as.numeric(annovarFreq$AF_afr)
annovarFreq$AF_amr=as.numeric(annovarFreq$AF_amr)
annovarFreq$AF_asj=as.numeric(annovarFreq$AF_asj)
annovarFreq$AF_eas=as.numeric(annovarFreq$AF_eas)
annovarFreq$AF_fin=as.numeric(annovarFreq$AF_fin)
annovarFreq$AF_nfe=as.numeric(annovarFreq$AF_nfe)
annovarFreq$AF_sas=as.numeric(annovarFreq$AF_sas)


annovarFreq=subset(annovarFreq, annovarFreq$AF >=0.01 |
                     annovarFreq$AF_afr >=0.01 |
                     annovarFreq$AF_amr >=0.01 |
                     annovarFreq$AF_asj >=0.01 |
                     annovarFreq$AF_eas >=0.01 |
                     annovarFreq$AF_fin >=0.01 |
                     annovarFreq$AF_nfe >=0.01 |
                     annovarFreq$AF_sas >=0.01)
annovarFreq=annovarFreq[,c(1:5,11:24,6:10)]
exc_notInReference$line=seq(1:234)
colnames(annovarFreq)[19]="line"
annovarFreqF=left_join(exc_notInReference,annovarFreq, by="line")
annovarFreqF=annovarFreqF[,-c(15:17,33,34,37:39,41,47,49:53)]
colnames(annovarFreqF)[29]="Var.Pos"




hist(annovarFreqF$AF, col = "blue", 
     main = paste("Variants frequency distribution\nGnomad All"),
     ylim = range(0:250),
     xlab = "Population requency", ylab= "# of variants")
hist(annovarFreqF$AF_afr, col = "blue", 
     main = paste("Variants frequency distribution\nGnomad afr"),
     ylim = range(0:250),
     xlab = "Population requency", ylab= "# of variants")
hist(annovarFreqF$AF_amr, col = "blue", 
     main = paste("Variants frequency distribution\nGnomad amr"),
     ylim = range(0:250),
     xlab = "Population requency", ylab= "# of variants")
hist(annovarFreqF$AF_asj, col = "blue", 
     main = paste("Variants frequency distribution\nGnomad asj"),
     ylim = range(0:250),
     xlab = "Population requency", ylab= "# of variants")
hist(annovarFreqF$AF_eas, col = "blue", 
     main = paste("Variants frequency distribution\nGnomad eas"),
     ylim = range(0:250),
     xlab = "Population requency", ylab= "# of variants")
hist(annovarFreqF$AF_fin, col = "blue", 
     main = paste("Variants frequency distribution\nGnomad fin"),
     ylim = range(0:250),
     xlab = "Population requency", ylab= "# of variants")
hist(annovarFreqF$AF_nfe, col = "blue", 
     main = paste("Variants frequency distribution\nGnomad nfe"),
     ylim = range(0:250),
     xlab = "Population requency", ylab= "# of variants")
hist(annovarFreqF$AF_sas, col = "blue", 
     main = paste("Variants frequency distribution\nGnomad sas"),
     ylim = range(0:250),
     xlab = "Population requency", ylab= "# of variants")




write.csv(annovarFreqF, file="exc_notInRef_annovar_freqNew.csv",
            row.names = FALSE)


#variants that matters for off-target detection were manually curated
#the file below contain only those variants
variantsMatters=read.csv("allele_freq_only_important_variants.csv")
hist(variantsMatters$AF, col = "blue", 
     main = paste("Variants frequency distribution\nGnomad all"),
     ylim = range(0:250),
     xlab = "Population frequency", ylab= "# of variants")


#####Venn diagram Total targets x Ref####

offtargetsPop=setdiff(totalTableAnn,exc_ref)
offtargetsPop=offtargetsPop$crRNAPos

offtargetsRef=subset(totalTableAnn,totalTableAnn$reference=="reference")
offtargetsRef=offtargetsRef$crRNAPos

myCol = brewer.pal(2, "Pastel2")
venn.diagram(x = list(offtargetsRef,offtargetsPop),
             category.names = c("Reference" , "Total"),
             filename = 'venn_diagramm.png',
             output=TRUE,
             lwd=2,  
             lty = 'blank',
               fill = myCol[c(1,2)],
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans"
             )




write.table(offtargetsPop, file="offtargetsPop_venndiag.txt", 
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(offtargetsRef, file="offtargetsRef_venndiag.txt", 
            row.names = FALSE, quote = FALSE, col.names = FALSE)



#####Number of off-targets per sgRNA####

offTargetsCounts=totalTableAll[(totalTableAll$CFD !=1),]
#ref
offTargetsCounts_ref =table((offTargetsCounts[(offTargetsCounts$higherCFD=="ref"),])$crRNANoBulges)
offTargetsCounts_ref=as.data.frame(offTargetsCounts_ref)


#gnomad_all
offTargetsCounts_all = table((offTargetsCounts[(offTargetsCounts$higherCFD=="gnomad_all"),])$crRNANoBulges)
offTargetsCounts_all=as.data.frame(offTargetsCounts_all)

#per pop
offTargetsCounts_afr = table((offTargetsCounts[(offTargetsCounts$higherCFD=="gnomad_afr"),])$crRNANoBulges)
offTargetsCounts_afr=as.data.frame(offTargetsCounts_afr)

offTargetsCounts_amr = table((offTargetsCounts[(offTargetsCounts$higherCFD=="gnomad_amr"),])$crRNANoBulges)
offTargetsCounts_amr=as.data.frame(offTargetsCounts_amr)

offTargetsCounts_asj = table((offTargetsCounts[(offTargetsCounts$higherCFD=="gnomad_asj"),])$crRNANoBulges)
offTargetsCounts_asj=as.data.frame(offTargetsCounts_asj)

offTargetsCounts_eas = table((offTargetsCounts[(offTargetsCounts$higherCFD=="gnomad_eas"),])$crRNANoBulges)
offTargetsCounts_eas=as.data.frame(offTargetsCounts_eas)

offTargetsCounts_fin = table((offTargetsCounts[(offTargetsCounts$higherCFD=="gnomad_fin"),])$crRNANoBulges)
offTargetsCounts_fin=as.data.frame(offTargetsCounts_fin)

offTargetsCounts_nfe = table((offTargetsCounts[(offTargetsCounts$higherCFD=="gnomad_nfe"),])$crRNANoBulges)
offTargetsCounts_nfe=as.data.frame(offTargetsCounts_nfe)

offTargetsCounts_sas = table((offTargetsCounts[(offTargetsCounts$higherCFD=="gnomad_sas"),])$crRNANoBulges)
offTargetsCounts_sas=as.data.frame(offTargetsCounts_sas)

offTargetsCounts=cbind(offTargetsCounts_ref,offTargetsCounts_all[,2],
                       offTargetsCounts_afr[,2],offTargetsCounts_amr[,2],
                       offTargetsCounts_asj[,2],offTargetsCounts_eas[,2],
                       offTargetsCounts_fin[,2],offTargetsCounts_nfe[,2],
                       offTargetsCounts_sas[,2])

offTargetsCounts=subset(offTargetsCounts, offTargetsCounts$Var1!="GACTATGCTGCCGCCCAGTNNNN" &
                          offTargetsCounts$Var1!= "GCAGAAGGGGACAGTAAGANNNN" &
                          offTargetsCounts$Var1!= "GTCCGCAGCTTTCTCGANNNNNN" &
                          offTargetsCounts$Var1!="GGCAGTTGTGTGACACGGAANNN" &
                          offTargetsCounts$Var1!="GGCGTGACTTCCACATGAGCNNN")

#whole data set
totalTableNoOn=totalTable[(totalTable$CFD !=1),]
offTargetsCounts_whole = table(totalTableNoOn$crRNANoBulges)
offTargetsCounts_whole=as.data.frame(offTargetsCounts_whole)

offTargetsCounts=cbind(offTargetsCounts,offTargetsCounts_whole[,2])
colnames(offTargetsCounts)[1]="crRNANoBulges"


offTargetsCounts=left_join(offTargetsCounts,sgRNAs, by="crRNANoBulges")
offTargetsCounts=offTargetsCounts[,c(1,12,2,11,3:10)]

colnames(offTargetsCounts)=c("sgRNA","gene","ref", "count_total",
                             "gnomad_all",
                             "gnomad_afr",
                             "gnomad_amr","gnomad_asj", 
                             "gnomad_eas", "gnomad_fin",
                             "gnomad_nfe", "gnomad_sas")

offTargetsCounts$percentage=(offTargetsCounts$count_total-offTargetsCounts$ref)/offTargetsCounts$ref

write.table(offTargetsCounts, file="offtargets_sgRNA_counts.txt", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)


#####Number of off-targets per pop#####
#Count the total number of variants found in each population and not
#found in reference and not found in gnomad_all
#the question behind this here is: if you only consider the variants 
#with freq >=1% in gnomad_all, you loose specific variants in a certain
#population
#combine the number of exclusive off targets per pop

#create a matrix with this data

offtargets_notIngnomadAll = as.data.frame(matrix(nrow=3,ncol=8))
colnames(offtargets_notIngnomadAll)=c("gnomad_all", "gnomad_afr",
                                      "gnomad_amr","gnomad_asj", 
                                      "gnomad_eas", "gnomad_fin",
                                      "gnomad_nfe", "gnomad_sas")
rownames(offtargets_notIngnomadAll)=c("Not in Reference", 
                                      "Not in gnomad_all",
                                      "Exclusive")


offtargets_notIngnomadAll$gnomad_all[1]=nrow(subset(totalTable, 
                                                    totalTable$gnomad_all=="gnomad_all" &
                                                      totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_afr[1]=nrow(subset(totalTable, 
                                                     totalTable$gnomad_afr=="gnomad_afr" &
                                                       totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_amr[1]=nrow(subset(totalTable, 
                                                     totalTable$gnomad_amr=="gnomad_amr" &
                                                       totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_asj[1]=nrow(subset(totalTable, 
                                                     totalTable$gnomad_asj=="gnomad_asj" &
                                                       totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_eas[1]=nrow(subset(totalTable, 
                                                     totalTable$gnomad_eas=="gnomad_eas" &
                                                       totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_fin[1]=nrow(subset(totalTable, 
                                                     totalTable$gnomad_fin=="gnomad_fin" &
                                                       totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_nfe[1]=nrow(subset(totalTable, 
                                                     totalTable$gnomad_nfe=="gnomad_nfe" &
                                                       totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_sas[1]=nrow(subset(totalTable, 
                                                     totalTable$gnomad_sas=="gnomad_sas" &
                                                       totalTable$reference == "0"))


offtargets_notIngnomadAll$gnomad_afr[2]=nrow(subset(totalTable, 
                                                    totalTable$gnomad_afr=="gnomad_afr" &
                                                      totalTable$gnomad_all =="0" &
                                                      totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_amr[2]=nrow(subset(totalTable, 
                                                    totalTable$gnomad_amr=="gnomad_amr" &
                                                      totalTable$gnomad_all =="0" &
                                                      totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_asj[2]=nrow(subset(totalTable, 
                                                    totalTable$gnomad_asj=="gnomad_asj" &
                                                      totalTable$gnomad_all =="0" &
                                                      totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_eas[2]=nrow(subset(totalTable, 
                                                    totalTable$gnomad_eas=="gnomad_eas" &
                                                      totalTable$gnomad_all =="0" &
                                                      totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_fin[2]=nrow(subset(totalTable, 
                                                    totalTable$gnomad_fin=="gnomad_fin" &
                                                      totalTable$gnomad_all =="0" &
                                                      totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_nfe[2]=nrow(subset(totalTable, 
                                                    totalTable$gnomad_nfe=="gnomad_nfe" &
                                                      totalTable$gnomad_all =="0" &
                                                      totalTable$reference == "0"))
offtargets_notIngnomadAll$gnomad_sas[2]=nrow(subset(totalTable, 
                                                    totalTable$gnomad_sas=="gnomad_sas" &
                                                      totalTable$gnomad_all =="0" &
                                                      totalTable$reference == "0"))

offtargets_notIngnomadAll$gnomad_afr[3]=nrow(exc_gnomad_afr)
offtargets_notIngnomadAll$gnomad_amr[3]=nrow(exc_gnomad_amr)
offtargets_notIngnomadAll$gnomad_asj[3]=nrow(exc_gnomad_asj)
offtargets_notIngnomadAll$gnomad_eas[3]=nrow(exc_gnomad_eas)
offtargets_notIngnomadAll$gnomad_fin[3]=nrow(exc_gnomad_fin)
offtargets_notIngnomadAll$gnomad_nfe[3]=nrow(exc_gnomad_nfe)
offtargets_notIngnomadAll$gnomad_sas[3]=nrow(exc_gnomad_sas)



offtargets_notIngnomadAll=t(offtargets_notIngnomadAll)
offtargets_notIngnomadAll=as.data.frame(offtargets_notIngnomadAll)
write.table(offtargets_notIngnomadAll, file="offtargets_Pops_counts.txt", 
            row.names = TRUE, quote = FALSE, col.names = TRUE)



  
######CFD distribution#####
cfd=totalTable[,12]
hist(cfd, col = "blue", 
     main = paste("CFD Distribution\nAll off-targets"),
     ylim = range(0:1200),
     xlab = "CFD", ylab= "# of off-targets")

cfd2=exc_notInReference[,12]
hist(cfd2, col = "blue", 
     main = paste("CFD Distribution\nOff-targets not in Reference"),
     ylim = range(1:150),
     xlab = "CFD", ylab= "# of off-targets")

cfd3=exc_gnomad_afr[,12]
hist(cfd3, col = "blue", 
     main = paste("CFD Distribution\nGnomad Afr exclusive off-targets"),
     ylim = range(1:120),
     xlab = "CFD", ylab= "# of off-targets")


#####Function distribution#####
genomefunc=as.factor(totalTableAnn[,27])
genomefunc=(table(genomefunc))
names=names(genomefunc)
pal=brewer.pal(11,"Set3")

barplot(genomefunc, cex.names = 0.4, ylim = range(0:1200),
        ylab = "# of off-targets", col = pal, xlim= range(0:12))

genomefunc=as.data.frame(genomefunc)
numbers=genomefunc[,2]
Functions=paste(genomefunc[,1]," ",
                "(",numbers,")")
ggplot(genomefunc, aes(x="", y=Freq, fill= Functions)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+ theme_void()+
  scale_fill_brewer(palette="Paired")
        
#only exc_notInRef
genomefunc=as.factor(exc_notInReference[,29])
genomefunc=(table(genomefunc))
names=names(genomefunc)
pal=brewer.pal(11,"Set3")

barplot(genomefunc, cex.names = 0.4, ylim = range(0:1200),
        ylab = "# of off-targets", col = pal, xlim= range(0:12))

genomefunc=as.data.frame(genomefunc)
numbers=genomefunc[,2]
Functions=paste(genomefunc[,1]," ",
                "(",numbers,")")
ggplot(genomefunc, aes(x="", y=Freq, fill= Functions)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+ theme_void()+
  scale_fill_brewer(palette="Paired")



#####Annotate variants with vep####
#Load vcf files generated above directly into VEP website
#abraom
vep_conseq_abraom=read.delim("VEP_abraom_conseq.txt")
vep_conseq_abraom=vep_conseq_abraom[,c(2:4)]
vep_abraom=read.delim("VEP_abraom.txt")
vep_abraom=vep_abraom[,c(2,6,11,12)]
vep_abraom = distinct(vep_abraom, Location, .keep_all = TRUE)#this will make some gene annotations to disapear
abraom_unique[,18]=paste(abraom_unique$V15,":",abraom_unique$Cluster,
                         "-",abraom_unique$Cluster, sep="")
abraom_unique = abraom_unique[order(abraom_unique$V15,abraom_unique$Cluster),]
colnames(abraom_unique)[18]="Location"
abraom_unique_ann=left_join(abraom_unique, vep_conseq_abraom,by="Location" )
abraom_unique_ann=left_join(abraom_unique_ann,vep_abraom, by="Location")
abraom_unique_ann=distinct(abraom_unique_ann,crRNAPos, .keep_all=TRUE)


#gnomad_afr
vep_conseq_gnomad_afr=read.delim("VEP_gnomad_afr_conseq.txt")
vep_conseq_gnomad_afr=vep_conseq_gnomad_afr[,c(2:4)]
vep_gnomad_afr=read.delim("VEP_gnomad_afr.txt")
vep_gnomad_afr=vep_gnomad_afr[,c(2,6,11,12)]
vep_gnomad_afr = distinct(vep_gnomad_afr, Location, .keep_all = TRUE)#this will make some gene annotations to disapear
gnomad_afr_unique[,18]=paste(gnomad_afr_unique$V15,":",gnomad_afr_unique$Cluster,
                             "-",gnomad_afr_unique$Cluster, sep="")
gnomad_afr_unique = gnomad_afr_unique[order(gnomad_afr_unique$V15,gnomad_afr_unique$Cluster),]
colnames(gnomad_afr_unique)[18]="Location"
gnomad_afr_unique_ann=left_join(gnomad_afr_unique, vep_conseq_gnomad_afr,by="Location" )
gnomad_afr_unique_ann=left_join(gnomad_afr_unique_ann,vep_gnomad_afr, by="Location")
gnomad_afr_unique_ann=distinct(gnomad_afr_unique_ann,crRNAPos, .keep_all=TRUE)


#gnomad_all
vep_conseq_gnomad_all=read.delim("VEP_gnomad_all_conseq.txt")
vep_conseq_gnomad_all=vep_conseq_gnomad_all[,c(2:4)]
vep_gnomad_all=read.delim("VEP_gnomad_all.txt")
vep_gnomad_all=vep_gnomad_all[,c(2,6,11,12)]
vep_gnomad_all = distinct(vep_gnomad_all, Location, .keep_all = TRUE)#this will make some gene annotations to disapear
gnomad_all_unique[,18]=paste(gnomad_all_unique$V15,":",gnomad_all_unique$Cluster,
                             "-",gnomad_all_unique$Cluster, sep="")
gnomad_all_unique = gnomad_all_unique[order(gnomad_all_unique$V15,gnomad_all_unique$Cluster),]
colnames(gnomad_all_unique)[18]="Location"
gnomad_all_unique_ann=left_join(gnomad_all_unique, vep_conseq_gnomad_all,by="Location" )
gnomad_all_unique_ann=left_join(gnomad_all_unique_ann,vep_gnomad_all, by="Location")
gnomad_all_unique_ann=distinct(gnomad_all_unique_ann,crRNAPos, .keep_all=TRUE)

#reference
vep_conseq_ref=read.delim("VEP_reference_conseq.txt")
vep_conseq_ref=vep_conseq_ref[,c(2:4)]
vep_ref=read.delim("VEP_reference.txt")
vep_ref=vep_ref[,c(2,6,11,12)]
vep_ref = distinct(vep_ref, Location, .keep_all = TRUE)#this will make some gene annotations to disapear
ref_unique[,18]=paste(ref_unique$V15,":",ref_unique$Cluster,
                      "-",ref_unique$Cluster, sep="")
ref_unique = ref_unique[order(ref_unique$V15,ref_unique$Cluster),]
colnames(ref_unique)[18]="Location"
ref_unique_ann=left_join(ref_unique, vep_conseq_ref,by="Location" )
ref_unique_ann=left_join(ref_unique_ann,vep_ref, by="Location")
ref_unique_ann=distinct(ref_unique_ann,crRNAPos, .keep_all=TRUE)


#####Analyzing only SCD directed therapies####
SCDTableAnn=subset(totalTableAnn, totalTableAnn$gene=="HBG_197"|
                     totalTableAnn$gene=="HBG_196"|
                     totalTableAnn$gene=="HBG_115"|
                     totalTableAnn$gene=="HBG_195"|
                     totalTableAnn$gene=="HBB_DEVER"|
                     totalTableAnn$gene=="HBB_PARK")
write.table(SCDTableAnn, file="SCD_FinalTargets_concatenated.txt", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)


SCDannovarFreqF=subset(annovarFreqF, annovarFreqF$gene=="HBG_197"|
                         annovarFreqF$gene=="HBG_196"|
                         annovarFreqF$gene=="HBG_115"|
                         annovarFreqF$gene=="HBG_195"|
                         annovarFreqF$gene=="HBB_DEVER"|
                         annovarFreqF$gene=="HBB_PARK")
write.table(SCDannovarFreqF, file="SCD_exc_notInRef_annovar_freq.txt", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

hist(SCDannovarFreqF$AF, col = "blue", 
     main = paste("SCD genes variant frequencies distribution\nGnomad All"),
     ylim = range(0:250),
     xlab = "Population frequency", ylab= "# of variants")
hist(annovarFreqF$AF_afr, col = "blue", 
     main = paste("SCD genes variant frequencies distribution\nGnomad afr"),
     ylim = range(0:250),
     xlab = "Population frequency", ylab= "# of variants")



######SCD CFD distribution#####
cfd=SCDTableAnn[,12]
hist(cfd, col = "blue", 
     main = paste("CFD Distribution\nSCD All off-targets"),
     ylim = range(0:800),
     xlab = "CFD", ylab= "# of off-targets")

cfd2=exc_notInReference[,12]
hist(cfd2, col = "blue", 
     main = paste("CFD Distribution\nSCD Off-targets not in Reference"),
     ylim = range(1:150),
     xlab = "CFD", ylab= "# of off-targets")

cfd3=exc_gnomad_afr[,11]
hist(cfd3, col = "blue", 
     main = paste("CFD Distribution\nGnomad Afr exclusive off-targets"),
     ylim = range(1:120),
     xlab = "CFD", ylab= "# of off-targets")


#####Function distribution#####
genomefunc=as.factor(SCDTableAnn[,27])
genomefunc=(table(genomefunc))
names=names(genomefunc)
pal=brewer.pal(11,"Set3")

barplot(genomefunc, cex.names = 0.4, ylim = range(0:1200),
        ylab = "# of off-targets", col = pal, xlim= range(0:12))

genomefunc=as.data.frame(genomefunc)
numbers=genomefunc[,2]
Functions=paste(genomefunc[,1]," ",
                "(",numbers,")")
ggplot(genomefunc, aes(x="", y=Freq, fill= Functions)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+ theme_void()+
  scale_fill_brewer(palette="Paired")

#####SCD Number of off-targets per pop#####
#Count the total number of variants found in each population and not
#found in reference and not found in gnomad_all
#the question behind this here is: if you only consider the variants 
#with freq >=1% in gnomad_all, you loose specific variants in a certain
#population
#combine the number of exclusive off targets per pop

#create a matrix with this data

SCD_offtargets_notIngnomadAll = as.data.frame(matrix(nrow=3,ncol=8))
colnames(SCD_offtargets_notIngnomadAll)=c("gnomad_all", "gnomad_afr",
                                      "gnomad_amr","gnomad_asj", 
                                      "gnomad_eas", "gnomad_fin",
                                      "gnomad_nfe", "gnomad_sas")
rownames(SCD_offtargets_notIngnomadAll)=c("Not in Reference", 
                                      "Not in gnomad_all",
                                      "Exclusive")


SCD_offtargets_notIngnomadAll$gnomad_all[1]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_all=="gnomad_all" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_afr[1]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_afr=="gnomad_afr" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_amr[1]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_amr=="gnomad_amr" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_asj[1]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_asj=="gnomad_asj" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_eas[1]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_eas=="gnomad_eas" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_fin[1]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_fin=="gnomad_fin" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_nfe[1]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_nfe=="gnomad_nfe" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_sas[1]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_sas=="gnomad_sas" &
                                                      SCDTableAnn$reference == "0"))


SCD_offtargets_notIngnomadAll$gnomad_afr[2]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_afr=="gnomad_afr" &
                                                      SCDTableAnn$gnomad_all =="0" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_amr[2]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_amr=="gnomad_amr" &
                                                      SCDTableAnn$gnomad_all =="0" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_asj[2]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_asj=="gnomad_asj" &
                                                      SCDTableAnn$gnomad_all =="0" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_eas[2]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_eas=="gnomad_eas" &
                                                      SCDTableAnn$gnomad_all =="0" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_fin[2]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_fin=="gnomad_fin" &
                                                      SCDTableAnn$gnomad_all =="0" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_nfe[2]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_nfe=="gnomad_nfe" &
                                                      SCDTableAnn$gnomad_all =="0" &
                                                      SCDTableAnn$reference == "0"))
SCD_offtargets_notIngnomadAll$gnomad_sas[2]=nrow(subset(SCDTableAnn, 
                                                    SCDTableAnn$gnomad_sas=="gnomad_sas" &
                                                      SCDTableAnn$gnomad_all =="0" &
                                                      SCDTableAnn$reference == "0"))

SCD_exc_gnomad_afr=subset(exc_gnomad_afr,exc_gnomad_afr$gene=="HBG_197"|
                            exc_gnomad_afr$gene=="HBG_196"|
                            exc_gnomad_afr$gene=="HBG_195"|
                            exc_gnomad_afr$gene=="HBG_115"|
                            exc_gnomad_afr$gene=="HBB_DEVER"|
                            exc_gnomad_afr$gene=="HBG_PARK")


SCD_offtargets_notIngnomadAll$gnomad_afr[3]=nrow(SCD_exc_gnomad_afr)
SCD_offtargets_notIngnomadAll$gnomad_amr[3]=0
SCD_offtargets_notIngnomadAll$gnomad_asj[3]=5
SCD_offtargets_notIngnomadAll$gnomad_eas[3]=7
SCD_offtargets_notIngnomadAll$gnomad_fin[3]=7
SCD_offtargets_notIngnomadAll$gnomad_nfe[3]=2
SCD_offtargets_notIngnomadAll$gnomad_sas[3]=7



SCD_offtargets_notIngnomadAll=t(SCD_offtargets_notIngnomadAll)
SCD_offtargets_notIngnomadAll=as.data.frame(SCD_offtargets_notIngnomadAll)
write.table(SCD_offtargets_notIngnomadAll, file="SCD_offtargets_Pops_counts.txt", 
            row.names = TRUE, quote = FALSE, col.names = TRUE)

