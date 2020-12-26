Large-scale_chromatin_structure

#get 2M fregments

grep chrUn_NW -v susScr11.chrom.sizes > susScr11.chrom.sizes_nochrUn.txt
#!/bin/bash 
rm 1.txt 2.txt 4.txt
cat susScr11.chrom.sizes_nochrUn.txt | while read id;
do
echo $id
arr=($id)
size=${arr[1]}
chr=${arr[0]}
sum=0
i=0
square=0
  if [ "$size" -gt 2000000 ];
  then
    i=2000000

      while(( i <= $size ))
      do
      echo $chr >>4.txt
      echo $square >>1.txt
      echo $i>>2.txt
      let "square=i+1"   
      let "sum+=i"
      let "i += 2000000"
      done
      echo $chr >>4.txt
      echo $square >>1.txt
      echo  $size >> 2.txt

  else
    echo $chr >>4.txt
    echo $i >>1.txt
    echo  $size >> 2.txt 
  fi

done
paste 4.txt 1.txt 2.txt | sort -k1,1 -k2,2n > 2M_size_genome.txt


#state frecancy state bin/total bin
#!/bin/bash 
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify
ls *15_segments.bed | while read tissue
do
    echo $tissue
    rm "Large_scale_chromatin_stru/"$(basename $tissue "_15_segments.bed")"_15_Large_scale.txt"
    cat /group/zhougrp/zhangyuan/genome/pig/2M_size_genome.txt | while read id;
    do
    echo $id
    arr=($id)
    chr=${arr[0]}  
    start=${arr[1]}
    end=${arr[2]}
    echo ${chr} 
    echo ${start} 
    echo ${end}

    echo "cat "$tissue"| awk '{if (\$1==\""$chr"\" && \$2>="$start" && \$3<= "$end" ) print}' > segment.txt" > ${dd}.sh; bash ${dd}.sh
         for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
         do
         A=$(grep -w $i segment.txt | wc -l  | awk '{print $1}')
         B=$(cat segment.txt| wc -l| awk '{print $1}')
         if [ $B = 0 ]||[ $A = 0 ]; 
         then
         let A=0 
         let B=1
         echo $B
         fi
         echo ${tissue} ${chr} ${start} ${end} ${i}  $(bc <<< "scale=10;$A/$B") >> "Large_scale_chromatin_stru/"$(basename $tissue "_15_segments.bed")"_15_Large_scale.txt"
         done
    done
done


#average state frecancy of 14 tissues
#!/bin/bash 
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_pig_five_all_maker_15_last_modify/Large_scale_chromatin_stru
ls *_15_Large_scale.txt | while read id;
do
  cat $id | cut -d " " -f 6 | paste -s >> 2.txt
done
awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' 2.txt > 3.txt
rm 2.txt
cat 3.txt | awk '{print ($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14)/14 }' > 4.txt
cut -d " " -f 1-5 Jejunum_15_Large_scale.txt > 5.txt
paste 5.txt 4.txt |sed 's/ /\t/g' > All_15_Large_scale.csv

for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
  grep -w $i All_15_Large_scale.csv | cut -f 6  | paste -s >> 2.txt
done
awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' 2.txt > 3.txt
rm 2.txt
grep -w E1 Jejunum_15_Large_scale.txt | cut -d " " -f 1-4 >5.txt
paste 5.txt 3.txt |sed 's/ /\t/g' > 6.txt
cat -n 6.txt >All_15_Large_scale_last.csv



