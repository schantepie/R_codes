#!/bin/bash

cd $TMPDIR	#fichier temporaire de la machine qui t'acceuil
lieu1=$3	#recupere l argument 2 = lieu1  
lieu2=$4


#cp  $lieu1/* .				#copie ton dossier dans le temporaire
#cp -r $lieu2/ !"(*.R|*.sh)" .
rsync -av --exclude '*.R' --exclude '*.sh' $lieu2/ .
cp "${lieu1}.R" .


ln -s /afs/in2p3.fr/home/throng/mnhn/bin/R ./R	#cre un lanceur R dans le temporaire
cat *.R | R --slave --args $2			#lance ton script dans R et fais passer  l argument 1 = le numero de ton job

# mkdir "${lieu1}_results"				#cre un dossier results
mv *.R* "${lieu1}_results"				#copie Rdata et script dans results



