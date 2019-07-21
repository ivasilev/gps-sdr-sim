//============================================================================
//
// range.h - representation of a satelite range
//
//============================================================================

#ifndef _RANGE_H_
#define _RANGE_H_

#include "gpstime.h"

class Ephemeris;
class Ionoutc;
class GpsTime;

class Range
{
public:
	GpsTime g;
	double range; // pseudorange
	double rate;
	double d; // geometric distance
	double azel[2];
	double iono_delay;

public:
    void Compute(const Ephemeris& eph, const Ionoutc& ionoutc, const GpsTime& g, const double xyz[]);

private:
    static double ionosphericDelay(const Ionoutc& ionoutc, const GpsTime& g, const double * const llh, const double * const azel);
};

#endif  // _RANGE_H_
