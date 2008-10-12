/* LIFTED AS A BACKSTOP FROM snpMatrix 1.5.4, as 1.5.5 has SLOW snp.rhs.tests */
int wcenter(const double *y, int n, const double *weight, const int *stratum, 
	    int nstrata, int resid, double *ynew);

int wresid(const double *y, int n, const double *weight, const double *x, 
	   double *ynew);

double wssq(const double *y, int n, const double *weight);

double wsum(const double *y, int n, const double *weight);

double wspr(const double *y, const double *x, int n, const double *weight);
