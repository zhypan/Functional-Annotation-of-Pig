#!/bin/bash
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify
echo OUTPUTSAMPLE_pig_five_maker_15 > fold_enrichment_region_Ajust_TSS.txt
ls *segments.bed | while read id
do
  ls /group/zhougrp/zhangyuan/pig_RNA_1/expressiondir/Average/Ajust_TSS/*.bed | while read bedfile
  do
     echo $bedfile
     ls -lh $bedfile
     echo state $id ${bedfile##*/} $(basename $bedfile .bed) >> fold_enrichment_region_Ajust_TSS.txt
     for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
     do
     A=$(grep -w $i $id | awk '{print ($3-$2) }' |awk '{sum+=$1} END {print sum}')
     B=$(cat $bedfile | awk '{print ($3-$2) }' |awk '{sum+=$1} END {print sum}')
     C=$(bedtools intersect  -a <(grep -w $i $id) -b $bedfile |awk '{print ($3-$2) }' |awk '{sum+=$1} END {print sum}')
     D=2478444698
     if [ "$C" = "" ]; 
     then 
     let C=0
     echo $C
     fi
     echo $i $id ${bedfile##*/} $(bc <<< "scale=10;($C/$A)/($B/$D)") >> fold_enrichment_region_Ajust_TSS.txt
     done
  done
done
sed 's/ /\t/g' fold_enrichment_region_Ajust_TSS.txt > fold_enrichment_region_Ajust_TSS.csv



touch fold_enrichment.txt
ls *segments.bed | while read id
do
  echo $id
  grep -w $id fold_enrichment_region_Ajust_TSS.csv > 1.txt
  ls /group/zhougrp/zhangyuan/pig_RNA_1/expressiondir/Average/Ajust_TSS/*.bed | while read bedfile
  do
    grep -w ${bedfile##*/} 1.txt | cut -f 4 | paste -s >> 2.txt
  done
  awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' 2.txt > $(basename $id ".bed")_fold_enrich1.txt
  rm 2.txt
  paste state.txt $(basename $id ".bed")_fold_enrich1.txt |sed 's/ /\t/g' > "enrichment/"$(basename $id ".bed")_fold_enrich_Ajust_TSS.txt
  rm $(basename $id ".bed")_fold_enrich1.txt
done

