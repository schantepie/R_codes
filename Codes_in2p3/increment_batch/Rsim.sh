#!/bin/bash

cd $TMPDIR	#emplacement de la machine qui d'acceuil temporaire
lieu1=$1	#recupere l'argument 2 = lieu1  
numjob=$2

rsync -av --exclude '*.sh' $lieu1/ ./ #importe les datas sur la machine d'acceuil

ln -s /afs/in2p3.fr/home/throng/mnhn/bin/R ./R	#cre un lanceur R dans le temporaire
cat *.R | R --slave --args $numjob			#lance ton script dans R et fais passer  l argument 1 = le numero du job

mkdir "${lieu1}results"				#cre un dossier results
mv *.R* "${lieu1}results"				#copie Rdata et script dans results