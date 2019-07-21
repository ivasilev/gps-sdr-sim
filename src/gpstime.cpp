//============================================================================
//
// time.cpp - timing-related functionality, implements DateTime and GpsTime
//
//============================================================================


#include "gpstime.h"
#include <math.h>

GpsTime::GpsTime() : week(0), sec(0) {} 

/*! \brief Convert a UTC date into a GPS date
 *  \param[in] t input date in UTC form
 *  \param[out] g output date in GPS form
 */
GpsTime::GpsTime(const DateTime& t)
{
	int doy[12] = {0,31,59,90,120,151,181,212,243,273,304,334};
	int ye;
	int de;
	int lpdays;

	ye = t.y - 1980;

	// Compute the number of leap days since Jan 5/Jan 6, 1980.
	lpdays = ye/4 + 1;
	if ((ye%4)==0 && t.m<=2)
		lpdays--;

	// Compute the number of days elapsed since Jan 5/Jan 6, 1980.
	de = ye*365 + doy[t.m-1] + t.d + lpdays - 6;

	// Convert time to GPS weeks and seconds.
	week = de / 7;
	sec = (double)(de%7)*SECONDS_IN_DAY + t.hh*SECONDS_IN_HOUR 
		+ t.mm*SECONDS_IN_MINUTE + t.sec;

	return;
}

GpsTime::GpsTime(const GpsTime& g, const double dt)
{
	week = g.week;
	sec = g.sec + dt;

	sec = round(sec*1000.0)/1000.0; // Avoid rounding error

	while (sec>=SECONDS_IN_WEEK)
	{
		sec -= SECONDS_IN_WEEK;
		week++;
	}

	while (sec<0.0)
	{
		sec += SECONDS_IN_WEEK;
		week--;
	}
}

double GpsTime::Sub(const GpsTime& g)
{
	double dt;

	dt = sec - g.sec;
	dt += (double)(week - g.week) * SECONDS_IN_WEEK;

	return(dt);
}

DateTime::DateTime() : y(0), m(0), d(0), hh(0), mm(0), sec(0) {}

DateTime::DateTime(const GpsTime& g)
{
	// Convert Julian day number to calendar date
	int c = (int)(7*g.week + floor(g.sec/86400.0)+2444245.0) + 1537;
	int d = (int)((c-122.1)/365.25);
	int e = 365*d + d/4;
	int f = (int)((c-e)/30.6001);

	d = c - e - (int)(30.6001*f);
	m = f - 1 - 12*(f/14);
	y = d - 4715 - ((7 + m)/10);

	hh = ((int)(g.sec/3600.0))%24;
	mm = ((int)(g.sec/60.0))%60;
	sec = g.sec - 60.0*floor(g.sec/60.0);

	return;
}

