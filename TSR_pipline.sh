
#For enhancer E6


#!/bin/bash
####get each state for each tissues
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify
ls *15_segments.bed | while read id;
do
  echo $id
  b=$(basename $id "_15_segments.bed")
  echo $b
  for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
  do
     grep -w $state $id | sort -k1,1 -k2,2n > "state_variability/"$b"_"$state".bed"

  done
done  


#Gs
####get Gs: total regions of each state (Gs) across 14 tissues 
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability
for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
cat *$i".bed" |sort -k1,1 -k2,2n > 1.bed
bedtools merge -i 1.bed > "AAGs/"$i"_Gs.bed"
done


#E6_Gs
## E6 is EnhA enhancer, count the E6 bins of each tissue in E6 Gs regions
#!/bin/bash
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability
ls *E6*bed | while read id;
do
echo $id
bedtools intersect -a <(sort -k1,1 -k2,2n /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AAGs/E6_Gs.bed) -b <(sort -k1,1 -k2,2n ${id}) -c -sorted >  "AARegulatory_module/"${id%%.*}"_Gs_enhancer.bed"
done




#paste all tissues together
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AARegulatory_module
ls *E6_TssA_EnhA_Gs.bed  | while read id;
do
  echo $id
  cat ${id} | cut -f 4 | paste -s >> 2.txt
done
awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' 2.txt > 3.txt
rm 2.txt
cut  -f 1-3 Cecum_E6_TssA_EnhA_Gs.bed > 5.txt
paste 5.txt 3.txt |sed 's/ /\t/g' >6.txt
cat header.txt 6.txt |sed 's/ /\t/g' > all_E6_TssA_EnhA_Gs.csv


#normalized the one count and count the number for each region
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AARegulatory_module
cat all_E6_Gs_enhancer.csv | cut -f 4-17 | sed 's/10/1/g'| sed 's/11/1/g'| sed 's/12/1/g'| sed 's/13/1/g'| sed 's/14/1/g'| sed 's/15/1/g' | sed 's/16/1/g'| sed 's/17/1/g'| sed 's/2/1/g'| sed 's/3/1/g'| sed 's/4/1/g'| sed 's/5/1/g' | sed 's/6/1/g'| sed 's/7/1/g'|  sed 's/8/1/g'| sed 's/9/1/g' > 3.txt
cat 3.txt | awk '{for(i=1;i<=NF;i++){a[NR]+=$i}print $0,a[NR]}' > 4.txt
cut  -f 1-3 all_E6_Gs_enhancer.csv > 5.txt
paste  5.txt 4.txt |sed 's/ /\t/g' > all_E6_Gs_enhancer_one_count.csv



#extract 17 module
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AARegulatory_module
cat all_E6_Gs_enhancer_one_count.csv | awk '{if (($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==14) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_All_common_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if (($5+$7+$9+$11+$12)==5 && ($4+$6+$8+$10+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_gut_common_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if (($6+$8+$10)==3 && ($4+$5+$7+$9+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Brain_common_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($4==1 && ($5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Adipose_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($5==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Cecum_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($6==1 && ($4+$5+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Cerebellum_E6_.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($7==1 && ($4+$6+$5+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Colon_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($8==1 && ($4+$6+$7+$5+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Cortex_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($9==1 && ($4+$6+$7+$8+$5+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Duodenum_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($10==1 && ($4+$6+$7+$8+$9+$5+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Hypothalamus_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($11==1 && ($4+$6+$7+$8+$9+$10+$5+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Ileum_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($12==1 && ($4+$6+$7+$8+$9+$10+$11+$5+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Jejunum_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($13==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$5+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Liver_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($14==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$5+$15+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Lung_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($15==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$5+$16+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Muscle_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($16==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$5+$17)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Spleen_E6.txt
cat all_E6_Gs_enhancer_one_count.csv | awk '{if ($17==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$5)==0) print $0}' > AA_TSR_E6_Gs_enhancer/TSR_Stomach_E6.txt







#For promoter E1

#!/bin/bash
####get each state for each tissues
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify
ls *15_segments.bed | while read id;
do
  echo $id
  b=$(basename $id "_15_segments.bed")
  echo $b
  for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
  do
     grep -w $state $id | sort -k1,1 -k2,2n > "state_variability/"$b"_"$state".bed"

  done
done  


#Gs
####get Gs: total regions of each state (Gs) across 14 tissues 
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability
for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
cat *$i".bed" |sort -k1,1 -k2,2n > 1.bed
bedtools merge -i 1.bed > "AAGs/"$i"_Gs.bed"
done


#E1_Gs
## E1 is TssA promoter, count the E1 bins of each tissue in E1 Gs regions
#!/bin/bash
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability
ls *E1.bed | while read id;
do
echo $id
bedtools intersect -a <(sort -k1,1 -k2,2n /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AAGs/E1_Gs.bed) -b <(sort -k1,1 -k2,2n ${id}) -c -sorted >  "AARegulatory_module/"${id%%.*}"_Gs_promoter.bed"
done



#paste all tissues together
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AARegulatory_module
ls *E1_Gs_promoter.bed  | while read id;
do
  echo $id
  cat ${id} | cut -f 4 | paste -s >> 2.txt
done
awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' 2.txt > 3.txt
rm 2.txt
cut  -f 1-3 Cecum_E1_Gs_promoter.bed > 5.txt
paste 5.txt 3.txt |sed 's/ /\t/g' >6.txt
cat header.txt 6.txt |sed 's/ /\t/g' > all_E1_Gs_promoter.csv


#normalized the one count and count the number for each region
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AARegulatory_module
cat all_E1_Gs_promoter.csv | cut -f 4-17 | sed 's/10/1/g'| sed 's/11/1/g'| sed 's/12/1/g'| sed 's/13/1/g'| sed 's/14/1/g'| sed 's/15/1/g' | sed 's/16/1/g'| sed 's/17/1/g'| sed 's/2/1/g'| sed 's/3/1/g'| sed 's/4/1/g'| sed 's/5/1/g' | sed 's/6/1/g'| sed 's/7/1/g'|  sed 's/8/1/g'| sed 's/9/1/g' > 3.txt
cat 3.txt | awk '{for(i=1;i<=NF;i++){a[NR]+=$i}print $0,a[NR]}' > 4.txt
cut  -f 1-3 all_E1_Gs_promoter.csv> 5.txt
paste  5.txt 4.txt |sed 's/ /\t/g' > all_E1_Gs_promoter_one_count.csv


#extract 17 module
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AARegulatory_module
cat all_E1_Gs_promoter_one_count.csv | awk '{if (($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==14) print $0}' > AA_TSR_E1/TSR_All_common_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if (($5+$7+$9+$11+$12)==5 && ($4+$6+$8+$10+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_gut_common_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if (($6+$8+$10)==3 && ($4+$5+$7+$9+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Brain_common_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($4==1 && ($5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Adipose_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($5==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Cecum_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($6==1 && ($4+$5+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Cerebellum_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($7==1 && ($4+$6+$5+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Colon_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($8==1 && ($4+$6+$7+$5+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Cortex_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($9==1 && ($4+$6+$7+$8+$5+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Duodenum_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($10==1 && ($4+$6+$7+$8+$9+$5+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Hypothalamus_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($11==1 && ($4+$6+$7+$8+$9+$10+$5+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Ileum_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($12==1 && ($4+$6+$7+$8+$9+$10+$11+$5+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Jejunum_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($13==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$5+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Liver_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($14==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$5+$15+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Lung_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($15==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$5+$16+$17)==0) print $0}' > AA_TSR_E1/TSR_Muscle_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($16==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$5+$17)==0) print $0}' > AA_TSR_E1/TSR_Spleen_E1_promoter.txt
cat all_E1_Gs_promoter_one_count.csv  | awk '{if ($17==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$5)==0) print $0}' > AA_TSR_E1/TSR_Stomach_E1_promoter.txt











#For repressor E13

#!/bin/bash
####get each state for each tissues
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify
ls *15_segments.bed | while read id;
do
  echo $id
  b=$(basename $id "_15_segments.bed")
  echo $b
  for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
  do
     grep -w $state $id | sort -k1,1 -k2,2n > "state_variability/"$b"_"$state".bed"

  done
done  


#Gs
####get Gs: total regions of each state (Gs) across 14 tissues 
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability
for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
cat *$i".bed" |sort -k1,1 -k2,2n > 1.bed
bedtools merge -i 1.bed > "AAGs/"$i"_Gs.bed"
done




#E13_Gs
#!/bin/bash
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability
ls *E13.bed | while read id;
do
echo $id
bedtools intersect -a <(sort -k1,1 -k2,2n /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AAGs/E13_Gs.bed) -b <(sort -k1,1 -k2,2n ${id}) -c -sorted >  "AARegulatory_module/"${id%%.*}"_Gs_repress.bed"
done



#paste all tissues together
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AARegulatory_module
ls *E13_Gs_repress.bed  | while read id;
do
  echo $id
  cat ${id} | cut -f 4 | paste -s >> 2.txt
done
awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' 2.txt > 3.txt
rm 2.txt
cut  -f 1-3 Cecum_E13_Gs_repress.bed > 5.txt
paste 5.txt 3.txt |sed 's/ /\t/g' >6.txt
cat header.txt 6.txt |sed 's/ /\t/g' > all_E13_Gs_repress.csv


#normalized the one count and count the number for each region
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AARegulatory_module
cat all_E13_Gs_repress.csv | cut -f 4-17 | sed 's/10/1/g'| sed 's/11/1/g'| sed 's/12/1/g'| sed 's/13/1/g'| sed 's/14/1/g'| sed 's/15/1/g' | sed 's/16/1/g'| sed 's/17/1/g'| sed 's/2/1/g'| sed 's/3/1/g'| sed 's/4/1/g'| sed 's/5/1/g' | sed 's/6/1/g'| sed 's/7/1/g'|  sed 's/8/1/g'| sed 's/9/1/g' > 3.txt
cat 3.txt | awk '{for(i=1;i<=NF;i++){a[NR]+=$i}print $0,a[NR]}' > 4.txt
cut  -f 1-3 all_E13_Gs_repress.csv > 5.txt
paste  5.txt 4.txt |sed 's/ /\t/g' > all_E13_Gs_repress_one_count.csv



#extract 17 module
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AARegulatory_module
cat all_E13_Gs_repress_one_count.csv | awk '{if (($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==14) print $0}' > AA_TSR_E13/TSR_All_common_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if (($5+$7+$9+$11+$12)==5 && ($4+$6+$8+$10+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_gut_common_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if (($6+$8+$10)==3 && ($4+$5+$7+$9+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Brain_common_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($4==1 && ($5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Adipose_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($5==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Cecum_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($6==1 && ($4+$5+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Cerebellum_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($7==1 && ($4+$6+$5+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Colon_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($8==1 && ($4+$6+$7+$5+$9+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Cortex_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($9==1 && ($4+$6+$7+$8+$5+$10+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Duodenum_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($10==1 && ($4+$6+$7+$8+$9+$5+$11+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Hypothalamus_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($11==1 && ($4+$6+$7+$8+$9+$10+$5+$12+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Ileum_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($12==1 && ($4+$6+$7+$8+$9+$10+$11+$5+$13+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Jejunum_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($13==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$5+$14+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Liver_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($14==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$5+$15+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Lung_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($15==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$5+$16+$17)==0) print $0}' > AA_TSR_E13/TSR_Muscle_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($16==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$5+$17)==0) print $0}' > AA_TSR_E13/TSR_Spleen_E13_repr.txt
cat all_E13_Gs_repress_one_count.csv  | awk '{if ($17==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$5)==0) print $0}' > AA_TSR_E13/TSR_Stomach_E13_repr.txt










