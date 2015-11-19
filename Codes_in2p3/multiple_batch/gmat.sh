#!/bin/bash


lieu="$PWD"

for i in LHTMORPH_exp2 LHTMORPH_SI2 morph_exp morph_SI LHT_exp LHT_SI 
  do 
    
  mkdir "$PWD/${i}_results"
  lieu1="$PWD/${i}"
  lieu2="$PWD/"
  for j in `seq 20`
    do
    echo "iteration_${j}"
    qsub  -l ct=46:00:00 -l vmem=2G -l fsize=50M -j y -o "${lieu1}_results"/qsout.txt $PWD/R2.sh ${i} ${j} $lieu1 $lieu2
    done
  done 