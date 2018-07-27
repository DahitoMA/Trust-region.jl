#!/bin/bash
# Pour lancer le fichier tu fais bash script_julia_SGE.sh
# Tu peux aussi faire chmod +x script_julia_SGE.sh une seule fois puis faire ./script_julia_SGE

for (( i=1; i <= 105; i++ ))
# -b pour dire si tu lances un binaire ou un script (réponse y ou n)
# -cwd pour dire que tu travailles dans le répertoire courant (Current Work Directory)
# -V pour faire passer les variables d'environnement
# -e pour designer le fichier qui contient le STDERR
# -o pour désigner le fichier qui contient le STDOUT
do
  hostname
  qsub -V -cwd -e STDERR/$i.txt -o STDOUT/$i.txt -l h_rt=1:30:00 -b y julia launch_CUTEst_pbs.jl $i
  sleep 50
done
# Attention les paramètres que tu fais passer à Julia sont dans la variable ARGS et de type String!
