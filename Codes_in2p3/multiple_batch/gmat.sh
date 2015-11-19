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
  
#   exemple if you want to run  a different number of jobs for somme .R 
#   
#   for i in  morph_SI LHT_exp LHT_SI # LHTMORPH_exp2
#   do 
#     
#   mkdir "$PWD/${i}_results"
#   lieu1="$PWD/${i}"
#   lieu2="$PWD/"
#   for j in `seq 5`
#     do
#     echo "iteration_${j}"
#     qsub  -l ct=46:00:00 -l vmem=2G -l fsize=50M -j y -o "${lieu1}_results"/qsout.txt $PWD/R2.sh ${i} ${j} $lieu1 $lieu2
#     done
#   done 