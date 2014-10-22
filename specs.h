//Physically relevant
#define ncoll 5000			// Max # of collisions
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
#define ideal_ambient 1		//For temperature outside, is the ambient gas ideal? (i.e. no collisions between non-blast)
#define planar 0		// 0= radial shock, 1= planar shock
#define radinit 0			// 1= Initial velocity of blast particles purely radial
#define reload 0			// TODO: Reload config (name in " " or 0 = no)
#define nbstart 0			// Starting number for output files (useful if reloading)

//Outputs
#define dobs 0			// Dimension of observation for moments of the velocity distribution
#define printpas 100			// # of collisions b/w full output  (position and velocity for all particles if nsys<nslim)
#define tstep -1			// time step b/w full output (>0 -> overrides prev line)
#define evolpas 10		// # of collisions b/w frequent output (a few quantities whose time evolution we want to follow)
#define colwait -1			// # collisions before starting measures

#define stopr 0			// Blocks all printf except errors

#define PI 3.14159265358979323846
