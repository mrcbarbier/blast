import os, subprocess,shlex

#Envoie tous les jobs au serveur simultanement

for rec in open('dirlist.dat'):
	directory=rec[:-1]
	subprocess.Popen(shlex.split('qsub ./'+directory+directory[:-1]+'.pbs'))
