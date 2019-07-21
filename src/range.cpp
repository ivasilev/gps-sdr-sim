//============================================================================
//
// range.cpp - range handling
//
//============================================================================

#include <math.h>

#include "range.h"
#include "ephemeris.h"
#include "ionoutc.h"
#include "gpstime.h"
#include "generic_funcs.h"

/*! \brief Compute range between a satellite and the receiver
 *  \param[out] rho The computed range
 *  \param[in] eph Ephemeris data of the satellite
 *  \param[in] g GPS time at time of receiving the signal
 *  \param[in] xyz position of the receiver
 */
void Range::Compute(const Ephemeris& eph, const Ionoutc& ionoutc, const GpsTime& paramTime, const double xyz[])
{
	double pos[3],vel[3],clk[2];
	double los[3];
	double tau;
	double localRange,localRate;
	double xrot,yrot;

	double llh[3],neu[3];
	double tmat[3][3];
	
	// SV position at time of the pseudorange observation.
	eph.Satpos(paramTime, pos, vel, clk);

	// Receiver to satellite vector and light-time.
	subVect(los, pos, xyz);
	tau = normVect(los)/SPEED_OF_LIGHT;

	// Extrapolate the satellite position backwards to the transmission time.
	pos[0] -= vel[0]*tau;
	pos[1] -= vel[1]*tau;
	pos[2] -= vel[2]*tau;

	// Earth rotation correction. The change in velocity can be neglected.
	xrot = pos[0] + pos[1]*OMEGA_EARTH*tau;
	yrot = pos[1] - pos[0]*OMEGA_EARTH*tau;
	pos[0] = xrot;
	pos[1] = yrot;

	// New observer to satellite vector and satellite range.
	subVect(los, pos, xyz);
	localRange = normVect(los);
	d = localRange;

	// Pseudorange.
	range = localRange - SPEED_OF_LIGHT*clk[0];

	// Relative velocity of SV and receiver.
	localRate = dotProd(vel, los)/localRange;

	// Pseudorange rate.
	rate = localRate; // - SPEED_OF_LIGHT*clk[1];

	// Time of application.
	g = paramTime;

	// Azimuth and elevation angles.
	xyz2llh(xyz, llh);
	ltcmat(llh, tmat);
	ecef2neu(los, tmat, neu);
	neu2azel(azel, neu);

	// Add ionospheric delay
	iono_delay = ionosphericDelay(ionoutc, g, llh, azel);
	range += iono_delay;

	return;
}

double Range::ionosphericDelay(const Ionoutc& ionoutc, const GpsTime& g, const double * const llh, const double * const azel)
{
	double iono_delay = 0.0;
	double E,phi_u,lam_u,F;

	if (ionoutc.enable==FALSE)
		return (0.0); // No ionospheric delay

	E = azel[1]/PI;
	phi_u = llh[0]/PI;
	lam_u = llh[1]/PI;

	// Obliquity factor
	F = 1.0 + 16.0*pow((0.53 - E),3.0);

	if (ionoutc.vflg==FALSE)
		iono_delay = F*5.0e-9*SPEED_OF_LIGHT;
	else
	{
		double t,psi,phi_i,lam_i,phi_m,phi_m2,phi_m3;
		double AMP,PER,X,X2,X4;

		// Earth's central angle between the user position and the earth projection of
		// ionospheric intersection point (semi-circles)
		psi = 0.0137/(E + 0.11) - 0.022;
		
		// Geodetic latitude of the earth projection of the ionospheric intersection point
		// (semi-circles)
		phi_i = phi_u + psi*cos(azel[0]);
		if(phi_i>0.416)
			phi_i = 0.416;
		else if(phi_i<-0.416)
			phi_i = -0.416;

		// Geodetic longitude of the earth projection of the ionospheric intersection point
		// (semi-circles)
		lam_i = lam_u + psi*sin(azel[0])/cos(phi_i*PI);

		// Geomagnetic latitude of the earth projection of the ionospheric intersection
		// point (mean ionospheric height assumed 350 km) (semi-circles)
		phi_m = phi_i + 0.064*cos((lam_i - 1.617)*PI);
		phi_m2 = phi_m*phi_m;
		phi_m3 = phi_m2*phi_m;

		AMP = ionoutc.alpha0 + ionoutc.alpha1*phi_m
			+ ionoutc.alpha2*phi_m2 + ionoutc.alpha3*phi_m3;
		if (AMP<0.0)
			AMP = 0.0;

		PER = ionoutc.beta0 + ionoutc.beta1*phi_m
			+ ionoutc.beta2*phi_m2 + ionoutc.beta3*phi_m3;
		if (PER<72000.0)
			PER = 72000.0;

		// Local time (sec)
		t = SECONDS_IN_DAY/2.0*lam_i + g.sec;
		while(t>=SECONDS_IN_DAY)
			t -= SECONDS_IN_DAY;
		while(t<0)
			t += SECONDS_IN_DAY;

		// Phase (radians)
		X = 2.0*PI*(t - 50400.0)/PER;

		if(fabs(X)<1.57)
		{
			X2 = X*X;
			X4 = X2*X2;
			iono_delay = F*(5.0e-9 + AMP*(1.0 - X2/2.0 + X4/24.0))*SPEED_OF_LIGHT;
		}
		else
			iono_delay = F*5.0e-9*SPEED_OF_LIGHT;
	}

	return (iono_delay);
}

