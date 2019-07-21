//============================================================================
//
// time.h - timing-related functionality
//
//============================================================================

#ifndef _GPSTIME_H_
#define _GPSTIME_H_

#include "constants.h"

// Forward declarations
class DateTime;
class GpsTime;

/*! \brief Structure representing GPS time */
class GpsTime
{
public:
	int week;	/*!< GPS week number (since January 1980) */
	double sec; 	/*!< second inside the GPS \a week */

public:
    GpsTime();
    GpsTime(const DateTime& t);
    GpsTime(const GpsTime& g, const double dt);
    double Sub(const GpsTime& g);
};

/*! \brief Structure repreenting UTC time */
class DateTime
{
public:
	int y; 		/*!< Calendar year */
	int m;		/*!< Calendar month */
	int d;		/*!< Calendar day */
	int hh;		/*!< Calendar hour */
	int mm;		/*!< Calendar minutes */
	double sec;	/*!< Calendar seconds */
public:
    DateTime();
    DateTime(const GpsTime& g);
};

#endif // _GPSTIME_H_
