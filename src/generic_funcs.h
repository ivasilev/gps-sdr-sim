#ifndef _GENERIC_FUNCS_H
#define _GENERIC_FUNCS_H


void subVect(double *y, const double *x1, const double *x2);
double normVect(const double *x);
double dotProd(const double *x1, const double *x2);

void xyz2llh(const double *xyz, double *llh);
void llh2xyz(const double *llh, double *xyz);
void ltcmat(const double *llh, double t[3][3]);
void ecef2neu(const double *xyz, double t[3][3], double *neu);
void neu2azel(double *azel, const double *neu);

int replaceExpDesignator(char *str, int len);

#endif // _GENERIC_FUNCS_H
