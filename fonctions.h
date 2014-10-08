
extern long int overlap ;
extern double epsilon;

int fbox(double * r, double *r0, double Rmax); //position to box id

double xbox(int b1, int k, double *r0, double Rmax); //box id to box center position

////////////////////////////////// CHAMPS BOITES ET RADIAUX //////////////////////////////////////////////////////////:

void r0calc(double ** r, double ** v, double * r0); // calcul de la position du centre de masse des particules en mouvement

double er(int j, double * r, double * r0);

void r0calc(double ** r, double ** v, double * r0);

void convrad(double * r, double *r0);

double rand01();

double timeval(double * r1, double * r2, double * v1, double * v2, double sigma);
