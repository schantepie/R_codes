#!/bin/bash

nbjob=300
time=46:00:00
memory=3G
taille=100M

echo qsub -l ct=${time} -l vmem=${memory} -l fsize=${taille} -j y -o $PWD/qsout.txt $PWD/Rsim.sh $PWD/ $"##_index_+1##" > submitfile.txt

/afs/in2p3.fr/home/c/calvat/public/jobsubmit/create_job_list.py -i ${nbjob} -t $PWD/submitfile.txt -o $PWD/task.txt

qsub -l demon=1 -q demon /afs/in2p3.fr/home/c/calvat/public/jobsubmit/jobs_submit.py -p $PWD/task.txt
