#!/bin/bash
bash /group/zhougrp/FAANG/Hi-C/juicer2.sh -d /group/zhougrp/zhangyuan/pig_HIC/$1 -q bigmemh -l bigmemh -g susScr11 -a '$1' -p /group/zhougrp/Genomes/susScr11/susScr11.chrom.sizes -y /group/zhougrp/FAANG/Hi-C/susScr11_HindIII.txt -z /group/zhougrp/Genomes/susScr11/susScr11.fa -D /home/mmhalste/Hi-C/
