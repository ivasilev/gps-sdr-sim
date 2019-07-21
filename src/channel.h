//============================================================================
//
// channel.h - single channel handling
//
//============================================================================

#ifndef _CHANNEL_H_
#define _CHANNEL_H_

#include "constants.h"
#include "gpstime.h"
#include "range.h"

class Ephemeris;
class Ionoutc;

/*! \brief Structure representing a Channel */
class Channel
{
public:
	int prn;	/*< PRN Number */
	int ca[CA_SEQ_LEN]; /*< C/A Sequence */
	double f_carr;	/*< Carrier frequency */
	double f_code;	/*< Code frequency */
	double carr_phase;
	double code_phase; /*< Code phase */
	GpsTime g0;	/*!< GPS time at start */
	unsigned long sbf[5][N_DWRD_SBF]; /*!< current subframe */
	unsigned long dwrd[N_DWRD]; /*!< Data words of sub-frame */
	int iword;	/*!< initial word */
	int ibit;	/*!< initial bit */
	int icode;	/*!< initial code */
	int dataBit;	/*!< current data bit */
	int codeCA;	/*!< current C/A code */
	double azel[2];
	Range rho0;

public:
    void ComputeCodePhase(const Range& rho1, const double dt);
    int GenerateNavMsg(const GpsTime&  g, const int init);
    void Eph2sbf(const Ephemeris& eph, const Ionoutc& ionoutc);

private:
    static unsigned long computeChecksum(unsigned long source, int nib);
    static unsigned long countBits(unsigned long v);

};

#endif // _CHANNEL_H_

