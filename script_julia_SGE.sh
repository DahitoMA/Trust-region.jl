#!/bin/bash
# Pour lancer le fichier tu fais bash script_julia_SGE.sh
# Tu peux aussi faire chmod +x script_julia_SGE.sh une seule fois puis faire ./script_julia_SGE

for (( i=1; i <= 24; i++ ))
# -b pour dire si tu lances un binaire ou un script (réponse y ou n)
# -cwd pour dire que tu travailles dans le répertoire courant (Current Work Directory)
# -V pour faire passer les variables d'environnement
# -e pour designer le fichier qui contient le STDERR
# -o pour désigner le fichier qui contient le STDOUT
do
  hostname
	qsub -V -cwd -b y julia launch_CUTEstLS_pbs.jl $i
	# qsub -V -cwd -b y julia somme.jl $i 1000
done
# Attention les paramètres que tu fais passer à Julia sont dans la variable ARGS et de type String!
# Bon courage
