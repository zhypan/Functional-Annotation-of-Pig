
#This script from Colin, and zhangyuan make some modify
#1 calculate the H3K27ac depth for each of the enhancers using the deepTools multiBamSummary command.

cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AAGs
sbatch -p high -t 1-0 -c 24 /group/zhougrp/FAANG/Final_Paper_Analysis/signal_count.sh E6_Gs.bed H3K27ac /group/zhougrp/zhangyuan/Pig_all/Aligned_Reads/H3K27ac_*.bam

#signal_count.sh:
#output is H3K27ac_counts.tsv
#!/bin/bash
#module load conda3
#source activate /group/zhougrp/FAANG/bosTau9/.snakemake/conda/ebbe0aed
module load deepTools
multiBamSummary BED-file --BED $1 --bamfiles ${@:3} -o ${2}_counts.npz -p 24 --outRawCounts ${2}_counts.tsv -e 200



#2 Next we have to normalize the counts in the same way that the RNA-seq counts are normalized. 
#First you have to use this script to convert the table into a different format for edgeR to do the normalizing:
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AAGs
mv H3K27ac_counts.tsv Pig_H3K27ac_counts.tsv
python prepare_table.py Pig_H3K27ac_counts.tsv > Pig_H3K27ac_counts_prepared.tsv
#python prepare_table.py: 
import sys
with open(sys.argv[1]) as f:
    samples = [x.strip("'") for x in f.readline().split()[3:]]
    print(','.join(sorted(samples)))
    for line in f:
        cols = line.split()
        name = (cols[0] + ":" + cols[1] + "-" + cols[2])
        line = {}
        for i, num in enumerate(cols[3:]):
            line[samples[i]] = num
        print(','.join([name] + [line[s] for s in sorted(samples)]))




#3 Then you can run the R script to do the normalizing:
Rscript /group/zhougrp/FAANG/Final_Paper_Analysis/edger_norm.R Pig_H3K27ac_counts_prepared.tsv Pig_H3K27ac_counts_cpms.csv
#edger_norm.R: 
library(edgeR)
args = commandArgs(trailingOnly=TRUE)
data <- read.table(args[1], sep=",", row.names=1, header=T)
RG <- DGEList(counts=data, group=rep(1,ncol(data)))
RG <- calcNormFactors(RG)
RG <- cpm(RG)
write.table(RG, file=args[2], sep=",", quote=F)





#4 Now you can finally run the target gene prediction. use all the CPU cores of a bigmem node, logging into a bigmem node. 
srun -p bigmemh -t1-0 -c62 --pty bash

#5 Then you can run the script. 
#You have to give it the Annotation_Expression_EdgeR files including TMM from EdgeR, then the H3K27ac_counts_cpms.csv you just created, 
#then the TSS for pig genome, 
#the TAD file creat from CTCF by FIMO, and finally the name of the output file to create:
python /group/zhougrp/FAANG/Final_Paper_Analysis/egcorr.py /group/zhougrp/zhangyuan/pig_RNA_1/Tables/Annotation_Expression_EdgeR_{P348,P350}.tsv Pig_H3K27ac_counts_cpms.csv /group/zhougrp/zhangyuan/genome/pig/part_change_esemb100/TSS_esemble100_colin_1.bed /group/zhougrp/FAANG/Final_Paper_Analysis/Pig_CTCF_TADs.bed Target_gene/output.tsv


#The output will have these columns:
1 - Enhancer region
2 - Spearman correlation of H3K27ac and RNA-seq across samples
3 - p-value for correlation
4 - Gene ID
5 - Distance of enhancer from TSS
6 - q-value for correlation

#6 used a q-value <= 0.05 to filter the confident predictions. You can do that with this command:
awk '{ if ($6 <= 0.05) print }' output.tsv > output_confident.tsv