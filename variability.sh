

#This script try to do Chromatin state variability
#1.1 acquire each state for each tissues
#!/bin/bash
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


#1.2 acquire Gs, All tissues of each state, merge the region 
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability
for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
cat *$i".bed" |sort -k1,1 -k2,2n > 1.bed
bedtools merge -i 1.bed > "AAGs/"$i"_Gs.bed"
done


#gs_number.txt
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AAGs
ls *_Gs.bed | while read id;
do
echo $id >> 1.txt 
cat $id| awk '{print ($3-$2)}' |awk '{sum+=$1} END {print sum}' >> 2.txt

done
paste 1.txt 2.txt > gs_number.txt
rm 1.txt 2.txt




#origin Tissue_state_lenth/Total_state_lenth
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability
rm AAAA.txt
for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
  ls *_${i}.bed | while read id;
  do
   C=$(basename $id "_${i}.bed")
   echo $C
   A=$(cat $id| awk '{print ($3-$2)}' |awk '{sum+=$1} END {print sum}')
   B=$(cat AAGs/${i}_Gs.bed | awk '{print ($3-$2)}' |awk '{sum+=$1} END {print sum}')
   echo ${i} $C $A $B $(bc <<< "scale=4; $A/$B ") >> AAAA.txt
  done
done
sed -i 's/ /\t/g' AAAA.txt



#prepare E1.txt-----E15.txt
for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
grep -w ${i} AAAA.txt | cut -f 5 | paste -s >> 2.txt
done
awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' 2.txt > 3.txt
rm 2.txt
paste tissues1.txt 3.txt |sed 's/ /\t/g' > Relative_genome_coverage1.csv
rm 3.txt



#chr_state_var.csv
#!/bin/bash
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability
for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
rm AAstatechange/${i}_state.txt
cat ${i}.txt |while read id;
  do echo ${id}
    cat $id| sort -k1,1 -k2,2n > 1.bed
    bedtools merge -i 1.bed > 2.bed
    A=$(cat 2.bed| awk '{print ($3-$2)}' |awk '{sum+=$1} END {print sum}')
    B=$(cat AAGs/${i}_Gs.bed | awk '{print ($3-$2)}' |awk '{sum+=$1} END {print sum}')
    echo $A $B $(bc <<< "scale=4; $A/$B ") >> AAstatechange/${i}_state.txt
  done
done

cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/state_variability/AAstatechange
for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
cat ${i}_state.txt | cut -d " " -f 3 | paste -s >> 2.txt
done
awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' 2.txt | sed 's/ /\t/g' > chr_state_var.csv
rm 2.txt











