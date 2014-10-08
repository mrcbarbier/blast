import os , subprocess,shlex


#Genere pour chaque valeur du parametre a echantilloner:
# - un dossier
# - spec.h approprie
# - une copie de la configuration initiale
# - xxx.out     l'executable
# - xxx.pbs     le job a envoyer sur le cluster


codepath1="../freec3clean.cpp" #Chemin vers le code source
codepath2="../fonctions.cpp" 
initconfig="../initconfig.dat" #Chemin vers la configuration initiale


#Valeurs fixes (recouvre valeurs par defaut ci-dessous)

val={}
val["nu"]=50000
val["ncoll"]=100000
val["dim"]=2
val["rdens"]=0.2
val["alpha"]=.7
val["printpas"]= 50000
val["tstep"]=-1


truebasedir=os.getcwd()+'/'
basedir=""
dirs=[]
param={}

#Parametre a echantillonner

param["masslaw"]=[]
indexlist=[]
for i in range (5) :
	indexlist.append( str(1*i-2))
	if len(indexlist[-1])==1 :
		indexlist[-1]=''+indexlist[-1]
for index in indexlist:
	dirs.append(""+str(index)+"/")
	param["masslaw"].append(int(index))


#Valeurs par defaut

define = "#define"

defaultspecs = '''
//Physically relevant
#define ncoll 500			// Max # of collisions
#define tmax 10.0			// Max rescaled time (whichever smallest)
#define dim 2			// Dimension of space
#define rdens 0.2			// Reduced density
#define size 1			// volume of the box
#define denslaw 1			// TODO: initial density distribution (1=homogeneous, other=power law exponent)
#define masslaw 0			// keeping number density uniform and changing particle masses with distance (0=homogeneous)
#define walltype 0			// boundary condition (0 periodic, 1 specular, 2 sticky)

#define alpha 0.8			// Restitution coefficient for inelastic collisions


//Precision and averaging
#define nu 2000			// Total particle number
#define nsys 1			// Number of systems for averaging
#define nslim 2			// Max number of systems for which all positions and velocities are printed to file

#define cutoff 0.0000000000000005			// Inelastic collapse regularization

//Initial configuration
#define partinit .3		// Initial energetic particles (if >1: nb of particles, if <1: radius/width)
#define blastMach -200		//Mach number of particles in blast (if >0 finite temperature outside, if<0 infinite Mach)
#define planar 0		// 0= radial shock, 1= planar shock
#define radinit 0			// 1= Initial velocity of blast particles purely radial
#define reload 0			// TODO: Reload config (name in " " or 0 = no)
#define nbstart 0			// Starting number for output files (useful if reloading)

//Outputs
#define dobs 0			// Dimension of observation for moments of the velocity distribution
#define evolpas 10		// # of collisions b/w frequent output (a few quantities whose time evolution we want to follow)
#define colwait -1			// # collisions before starting measures

#define stopr 0			// Blocks all printf except errors
#define PI 3.14159265358979323846'''


terms=[]
dflt=[]
comment=[]
headers=[]
for line in defaultspecs.split('\n'):
    if not line or line[:2]=='//':
        continue
    lsp1=line.split('//')
    if len(lsp1)>1:
        comment.append(lsp1[1])
    else:
        comment.append("")
    lsp2=lsp1[0].split()
    define,term,df=lsp2
    terms.append(term)
    dflt.append(df)

default ={}
com={}
for i in xrange(len(terms)):
	default[terms[i]]=dflt[i]
	if comment[i] != "" :
		com[terms[i]]="			// "+comment[i]

fn=open('dirlist.dat','w')
for direc in dirs :
	fn.write(str(direc)+'\n')
fn.close()


for z in xrange(len(dirs)) :
	directory=dirs[z]
	filename = basedir + directory + "specs.h"
	
	d = os.path.dirname(filename)
	if not os.path.exists(d):
		os.makedirs(d)
	fn = open(filename,'w')
	
	monjob = directory[:-1]
	job = open(basedir+directory+monjob+".pbs",'w')
	subname=monjob+".pbs"
	monjob+='.out'
#	#PBS -l cput=00:01:00,mem=10mb,nodes=1:ppn=1
	job.write('''
#!/bin/bash
#PBS -V
#PBS -q q-1sem

# On cree un repertoire temporaire sur le disque local du noeud de calcul
TMP=/tmp/$USER.$PBS_JOBID
mkdir $TMP

#on copie le programe + les fichiers d'entree
cd /'''+truebasedir[1:]+directory+'''
cp ./'''+monjob+''' $TMP
cp ./initconfig.dat $TMP

cd $TMP

#on lance le calcul
time  ./'''+monjob+''' > resultats

#on rapatrie les resultats dans le repertoire de donnee
mv -f $TMP/* '''+truebasedir+directory+'''
mv -f $TMP '''+truebasedir+'trash/')
	
	iline=0
	for line in defaultspecs.split('\n'):
		if not line or line.split()[0]!=define:
			fn.write(line+'\n')
			continue
		i=terms[iline]
		iline+=1
		fn.write(define +" "+ i+ " ")
		if i in param.keys() :
			fn.write(str(param[i][z]))
		elif i in val.keys():
			fn.write(str(val[i]))
		else :
			fn.write(default[i])
		if i in com.keys() :
			fn.write(com[i])
		fn.write("\n")

	print monjob
	#cmd =["g++", ' ~/Hydro/freec3b.cpp', "-O2","-lm",'-I', './040', "-o", './040/'+monjob]
	os.popen('cp {}  ./{}'.format(initconfig,directory))
	subprocess.Popen(shlex.split("g++ -o ./{}{} {} {}  -O2 -lm -I ./".format(directory,monjob,codepath1,codepath2) ),cwd='./')
	#cmd =["g++", ' ~/Hydro/hello.cpp']
	#subprocess.Popen(cmd)
	#subprocess.call('cd '+directory+';pwd; g++ -o ./'+monjob+' ~/Hydro/freec3b.cpp -O2 -lm -I ./' ,shell=True)
	#os.execv('/usr/bin/g++',['foo',"-o ./"+directory+monjob," ~/Hydro/freec3b.cpp"," -O2"," -lm ","-I ./"])
	#subprocess.call('ssh cluster ; qsub '+truebasedir+directory+subname+';exit',shell=True)
	
	fn.close()


