//============================================================================
//
// range.h - representation of a satelite range
//
//============================================================================

#ifndef _RANGE_H_
#define _RANGE_H_

#include "gpstime.h"

class Range
{
public:
	GpsTime g;
	double range; // pseudorange
	double rate;
	double d; // geometric distance
	double azel[2];
	double iono_delay;
};

#endif  // _RANGE_H_
