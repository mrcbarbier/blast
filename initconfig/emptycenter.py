
#Takes the configuration "initconfig.dat" and removes all particles within a radius rad, output into "initconfig2.dat"

fin = open("initconfig.dat","r")

fout=open("initconfig2.dat","w")

rad=0.03

for i in fin :
    tmp=i.split()
    if float(tmp[0])**2+float(tmp[1])**2> rad**2 :
        fout.write(i)

fin.close()
fout.close()
