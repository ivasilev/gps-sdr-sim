//============================================================================
//
// ionoutc.h - Ionospheric delay handling
//
//============================================================================

#ifndef _IONOUTC_H_
#define _IONOUTC_H_

class Ionoutc
{
public:
	int enable;
	int vflg;
	double alpha0,alpha1,alpha2,alpha3;
	double beta0,beta1,beta2,beta3;
	double A0,A1;
	int dtls,tot,wnt;
	int dtlsf,dn,wnlsf;
};

#endif // _IONOUTC_H_
