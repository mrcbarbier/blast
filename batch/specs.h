#define ncoll 5000000			// Max # of collisions
#define tmax 100.0			// Max time
#define len 100			// Resolution of the velocity distribution function (number of bins in the histogram)
#define dim 2			// Dimension of space
#define rdens 0.2			// Reduced density
#define size 1			// volume of the box
#define denslaw 1			// initial density distribution (1=homogeneous, other=power law)
#define masslaw 0			// same thing but keeping particle number constant and changing masses (0=homogeneous)
#define coldgas 0			// artificial cold gas (realign velocities)
#define walltype 0			// boundary condition (0 periodic, 1 specular, 2 sticky)

#define nu 50000			// Total particle number
#define nsys 1			// Number of systems for averaging
#define nslim 2			// Max number of systems for which all positions and velocities are printed

#define alpha 0.7			// Restitution coefficient for inelastic collisions
#define rescal 1			// Final rescaling of the velocities if alpha < 1.0 ?
#define dobs 0			// Dimension of observation for moments
#define cutoff 0.0000000000000005			// Inelastic collapse regularization

#define partinit 1			// Initial energetic particles (nb >1, radius < 1, 0 = all,-1=blast non zero temperature,-2 = planar shock)
#define radinit 0			// Initial velocity chosen radial
#define boxtyp 1			// Box type 0 square,1 polar
#define NBOX {50,1000}			// Box #
#define Gnbox 128			// bin # for correlation function
#define printpas 50000			// # collisions b/w full print
#define tstep -1			// time step b/w full print (>0 -> overrides prev line)
#define evolpas 10			// Nb of collisions b/w small print
#define colwait 0			// # collisions before starting measures

#define reload 0			// Reload config (name in " " or 0 = no)
#define nbstart 0			// Starting number for prints

#define stopr 0			// Blocks all printf except errors
#define PI 3.14159265358979323846
