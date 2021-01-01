Clustering of epigenomes with state

#1.devide the genome to 200bp window
#!/bin/bash
cd /group/zhougrp/zhangyuan/genome/pig
bedtools makewindows -g susScr11.chrom.sizes_nochrUn.txt -w 200 > 1.txt
sort -k1,1 -k2,2n 1.txt > susScr11.chrom.sizes_nochrUn_windows_200bp.bed

#2 bam to bed
#!/bin/bash
cd /group/zhougrp/zhangyuan/Pig_all/Aligned_Reads
ls *.bam | while read id;
do
echo $id
bedtools bamtobed -i $id | sort -k1,1 -k2,2n > ../Bed_file/${id%%.*}.bed
done

#3. get the raw count in each 200bp window of all assay
#!/bin/bash
cd /group/zhougrp/zhangyuan/Pig_all/Bed_file
ls H3K27ac*.bed H3K27me3*.bed H3K4me1*.bed H3K4me3*.bed Input*.bed  | while read id;
do
echo ${id}
bedtools intersect -a /group/zhougrp/zhangyuan/genome/pig/susScr11.chrom.sizes_nochrUn_windows_200bp.bed -b ${id} -c -sorted >  "../Bamcount/"${id%%.*}"_per_200bp.bedg"
done


#4. calculate –log10(Poisson P value) by R
#!/usr/bin/Rscript
setwd('/group/zhougrp/zhangyuan/Pig_all/Bamcount') 
## Load counts per window in chip sample
filename=read.table('Sample_list.txt')
filename=as.matrix(filename)
makers =c("H3K27ac","H3K27me3","H3K4me3","H3K4me1")
sink("hello.txt")
for(sample in filename)
{
    for (maker in makers){
    bedg = paste(maker, "_", sample, "_per_200bp.bedg", sep="")
    input_1 = paste( "Input_", sample, "_per_200bp.bedg", sep="")
    out_put=paste("../Poisson_P/", maker, "_", sample, "_per_200bp_2.txt", sep="")
    print(bedg)
    print(input_1)
    print(out_put)
    chip.bedg <- read.table(bedg)
    names(chip.bedg) <- c("chrom", "start", "end","counts")
    print(head(chip.bedg)) 
    ## Load counts per window in the input sample
    input.bedg <- read.table(input_1)
    names(input.bedg) <- c("chrom", "start", "end","counts")
    head(input.bedg)

    ## The fourth column of the bedgraph file contains the counts
    chip.counts <- chip.bedg[,4] 
    input.counts <- input.bedg[,4] 
    head(chip.counts)

    ## Count total counts
    chip.total <- sum(chip.counts)
    input.total <- sum(input.counts)


    ## Option 1: mormalize the data based on the library sum
    norm.input.libsum <- (input.counts+1) * chip.total / input.total
    head(norm.input.libsum)
    norm.input <- norm.input.libsum


    ## Choose a window size (200 or 50)
    window.size <- 200

    # a.  Number of reads in the test (FNR ChIP-seq)
    # b.  Number of reads in the input
    # c.  Normalized number of reads in the input (norm.input)
    read.stats <- cbind(chip.bedg, input.counts, norm.input)

    # d.  Reads per base in the normalized input
    read.stats$input.norm.rpb <- read.stats$norm.input/window.size
      
    # e.  Fold enrichment (ratio between test and normalized input)
    read.stats$fold <- read.stats$counts/read.stats$norm.input

    # f.  Log-fold enrichment (log(10) of the fold enrichment)
    read.stats$log.fold <- log(base=10, read.stats$fold)

    # sum(read.stats$norm.input != read.stats$input.counts*m/n) == 0
    ## Compute the Poisson p-value
    read.stats$ppois <- ppois(q=read.stats$counts-1, 
                              lambda=read.stats$norm.input, 
                              lower.tail=FALSE, log=FALSE)


    read.stats$log.ppois <- -log(base=10, read.stats$ppois)

    write.table(as.data.frame(read.stats), file=out_put, sep="\t",quote=F, row.names = FALSE)

}
}
sink()

less Sample_list.txt
Adipose_P348
Adipose_P350
Cecum_P348
Cecum_P350
Cerebellum_P348
Cerebellum_P350
Colon_P348
Colon_P350
Cortex_P348
Cortex_P350
Duodenum_P348
Duodenum_P350
Hypothalamus_P348
Hypothalamus_P350
Ileum_P348
Ileum_P350
Jejunum_P348
Jejunum_P350
Liver_P348
Liver_P350
Lung_P348
Lung_P350
Muscle_P348
Muscle_P350
Spleen_P348
Spleen_P350
Stomach_P348
Stomach_P350




#5 extract mark –log10(Poisson P value) for state. For H3K4me1, H3K27ac  we used state EnhA, for H3K4me3 state TssA; for H3K27me3 state ReprPC; 

#!/bin/bash
cd /group/zhougrp/zhangyuan/Pig_all/Poisson_P
ls *H3K4me1*_2.txt| while read id;
do
echo ${id}
bedtools intersect -wa -a <(sort -k1,1 -k2,2n ${id}) -b /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AAGs/E6_Gs.bed  -sorted >  "/group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig/"$(basename $id "_per_200bp_2.txt")"_me1_EnhA.txt"
done

#!/bin/bash
cd /group/zhougrp/zhangyuan/Pig_all/Poisson_P
ls *H3K27ac*_2.txt| while read id;
do
echo ${id}
bedtools intersect -wa -a <(sort -k1,1 -k2,2n ${id}) -b /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AAGs/E6_Gs.bed  -sorted >  "/group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig/"$(basename $id "_per_200bp_2.txt")"_ac_EnhA.txt"
done


H3K4me3

#!/bin/bash
cd /group/zhougrp/zhangyuan/Pig_all/Poisson_P
ls *H3K4me3*_2.txt| while read id;
do
echo ${id}
bedtools intersect -wa -a <(sort -k1,1 -k2,2n ${id}) -b /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AAGs/E1_Gs.bed  -sorted >  "/group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig/"$(basename $id "_per_200bp_2.txt")"_K4_TssA.txt"
done


 H3K27me3 

#!/bin/bash
cd /group/zhougrp/zhangyuan/Pig_all/Poisson_P
ls *H3K27me3*_2.txt| while read id;
do
echo ${id}
bedtools intersect -wa -a <(sort -k1,1 -k2,2n ${id}) -b /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AAGs/E13_Gs.bed  -sorted >  "/group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig/"$(basename $id "_per_200bp_2.txt")"_K27_Repr.txt"
done


#H3K4me1 in TssA
#!/bin/bash
cd /group/zhougrp/zhangyuan/Pig_all/Poisson_P
ls *H3K4me1*_2.txt| while read id;
do
echo ${id}
bedtools intersect -wa -a <(sort -k1,1 -k2,2n ${id}) -b /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AAGs/E1_Gs.bed  -sorted >  "/group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig/"$(basename $id "_per_200bp_2.txt")"_me1_TssA.txt"
done




#6 combine the –log10(Poisson P value) of smaples for state

cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig/
ls H3K4me1*me1_EnhA.txt | while read id;
do
  echo $id
  cat ${id} | cut -f 11 | paste -s >> 2.txt
done
awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' 2.txt > 3.txt
rm 2.txt
cut  -f 1-3 H3K4me1_Liver_P348_me1_EnhA.txt > 5.txt
paste 5.txt 3.txt |sed 's/ /\t/g' > all_tissue_me1_EnhA.csv

#获取组织序列
ls H3K4me1*me1_EnhA.txt | sed 's/H3K4me1_//g' | sed 's/_me1_EnhA.txt//g' | perl -p -e 's/\n/ /g' >  header.txt 
vim header
cat header.txt all_tissue_me1_EnhA.csv |sed 's/ /\t/g' > all_tissue_me1_EnhA_last.csv


#ac_EnhA
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig/
ls H3K27ac*ac_EnhA.txt | while read id;
do
  echo $id
  cat ${id} | cut -f 11 | paste -s >> 2.txt
done
awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' 2.txt > 3.txt
rm 2.txt
cut  -f 1-3 H3K27ac_Adipose_P348_ac_EnhA.txt > 5.txt
paste 5.txt 3.txt |sed 's/ /\t/g' > all_tissue_ac_EnhA.csv

cat header.txt all_tissue_ac_EnhA.csv |sed 's/ /\t/g' > all_tissue_ac_EnhA_last.csv



#K4_TssA.txt
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig/
ls H3K4me3*_K4_TssA.txt | while read id;
do
  echo $id
  cat ${id} | cut -f 11 | paste -s >> 2.txt
done
awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' 2.txt > 3.txt
rm 2.txt
cut  -f 1-3 H3K4me3_Adipose_P348_K4_TssA.txt > 5.txt
paste 5.txt 3.txt |sed 's/ /\t/g' > all_tissue_K4_TssA.csv

cat header.txt all_tissue_K4_TssA.csv |sed 's/ /\t/g' > all_tissue_K4_TssA_last.csv



#K27_Repr.txt
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig/
ls H3K27me3*_K27_Repr.txt | while read id;
do
  echo $id
  cat ${id} | cut -f 11 | paste -s >> 2.txt
done
awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' 2.txt > 3.txt
rm 2.txt
cut  -f 1-3 H3K27me3_Adipose_P348_K27_Repr.txt > 5.txt
paste 5.txt 3.txt |sed 's/ /\t/g' > all_tissue_K27_Repr.csv

cat header.txt all_tissue_K27_Repr.csv |sed 's/ /\t/g' > all_tissue_K27_Repr_last.csv


#H3K4me1 in TssA
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig/
ls H3K4me1*me1_TssA.txt| while read id;
do
  echo $id
  cat ${id} | cut -f 11 | paste -s >> 2.txt
done
awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' 2.txt > 3.txt
rm 2.txt
cut  -f 1-3 H3K4me1_Adipose_P348_me1_TssA.txt > 5.txt
paste 5.txt 3.txt |sed 's/ /\t/g' > all_tissue_me1_TssA.csv

cat header.txt all_tissue_me1_TssA.csv |sed 's/ /\t/g' > all_tissue_me1_TssA_last.csv




#7. hierarchical clustering in R
ac_EnhA

#!/usr/bin/Rscript
setwd('/group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig') 
data <- read.csv('all_tissue_ac_EnhA_last.csv', sep = "\t", header = TRUE) 
head(data)
data1 = data[,c(4:31)]
head(data1)
data1 <- t(data1)
df2 <- as.matrix((scale(data1)))

library(ggplot2)
library(dendextend)
library("factoextra")

result <- dist(df2, method = "euclidean")
result_hc <- hclust(d = result, method = "ward.D2")


pdf("all_tissue_ac_EnhA_color_1.pdf")
fviz_dend(result_hc, k = 4, 
          cex = 0.5, 
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, 
          rect = TRUE,
          horiz = TRUE,
          rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          rect_fill = TRUE         
)
dev.off()

require("igraph")
pdf("all_tissue_ac_EnhA_color_3.pdf")
fviz_dend(result_hc, k = 4, k_colors = "jco",
type = "phylogenic", repel = TRUE)
dev.off()




K27_Repr

#!/usr/bin/Rscript
setwd('/group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig') 
data <- read.csv('all_tissue_K27_Repr_last.csv', sep = "\t", header = TRUE) 
head(data)
data1 = data[,c(4:31)]
head(data1)
data1 <- t(data1)
df2 <- as.matrix((scale(data1)))

library(ggplot2)
library(dendextend)
library("factoextra")

result <- dist(df2, method = "euclidean")
result_hc <- hclust(d = result, method = "ward.D2")


pdf("all_tissue_K27_Repr_color_1.pdf")
fviz_dend(result_hc, k = 4, 
          cex = 0.5, 
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, 
          rect = TRUE,
          horiz = TRUE,
          rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          rect_fill = TRUE         
)
dev.off()

require("igraph")
pdf("all_tissue_K27_Repr_color_3.pdf")
fviz_dend(result_hc, k = 4, k_colors = "jco",
type = "phylogenic", repel = TRUE)
dev.off()




K4_TssA.txt

#!/usr/bin/Rscript
setwd('/group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig') 
data <- read.csv('all_tissue_K4_TssA_last.csv', sep = "\t", header = TRUE) 
head(data)
data1 = data[,c(4:31)]
head(data1)
data1 <- t(data1)
df2 <- as.matrix((scale(data1)))

library(ggplot2)
library(dendextend)
library("factoextra")

result <- dist(df2, method = "euclidean")
result_hc <- hclust(d = result, method = "ward.D2")


pdf("all_tissue_K4_TssA_color_1.pdf")
fviz_dend(result_hc, k = 4, 
          cex = 0.5, 
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, 
          rect = TRUE,
          horiz = TRUE,
          rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          rect_fill = TRUE         
)
dev.off()


me1_TssA.txt

#!/usr/bin/Rscript
setwd('/group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig') 
data <- read.csv('all_tissue_me1_TssA_last.csv', sep = "\t", header = TRUE) 
head(data)
data1 = data[,c(4:31)]
head(data1)
data1 <- t(data1)
df2 <- as.matrix((scale(data1)))

library(ggplot2)
library(dendextend)
library("factoextra")

result <- dist(df2, method = "euclidean")
result_hc <- hclust(d = result, method = "ward.D2")


pdf("all_tissue_me1_TssA_color_1.pdf")
fviz_dend(result_hc, k = 4, 
          cex = 0.5, 
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, 
          rect = TRUE,
          horiz = TRUE,
          rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          rect_fill = TRUE         
)
dev.off()


require("igraph")
pdf("all_tissue_me1_TssA_color_3.pdf")
fviz_dend(result_hc, k = 4, k_colors = "jco",
type = "phylogenic", repel = TRUE)
dev.off()


#try complete
e1_EnhA   complete

#!/usr/bin/Rscript
setwd('/group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/AA_Cluster_epig') 
data <- read.csv('all_tissue_me1_EnhA_last.csv', sep = "\t", header = TRUE) 
head(data)
data1 = data[,c(4:31)]
head(data1)
data1 <- t(data1)
df2 <- as.matrix((scale(data1)))

library(ggplot2)
library(dendextend)
library("factoextra")

result <- dist(df2, method = "euclidean")
result_hc <- hclust(d = result, method = "complete")


pdf("all_tissue_me1_EnhA_complete.pdf")
fviz_dend(result_hc, k = 4, 
          cex = 0.5, 
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, 
          rect = TRUE,
          horiz = TRUE,
          rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          rect_fill = TRUE         
)
dev.off()
