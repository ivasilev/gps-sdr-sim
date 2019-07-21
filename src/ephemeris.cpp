#include "ephemeris.h"
#include <math.h>

#include "generic_funcs.h"

/*! \brief Compute Satellite position, velocity and clock at given time
 *  \param[in] eph Ephemeris data of the satellite
 *  \param[in] g GPS time at which position is to be computed
 *  \param[out] pos Computed position (vector)
 *  \param[out] vel Computed velociy (vector)
 *  \param[clk] clk Computed clock
 */
void Ephemeris::Satpos(const GpsTime& g, double *pos, double *vel, double *clk) const
{
    // Computing Satellite Velocity using the Broadcast Ephemeris
    // http://www.ngs.noaa.gov/gps-toolbox/bc_velo.htm

    double tk;
    double mk;
    double ek;
    double ekold;
    double ekdot;
    double cek,sek;
    double pk;
    double pkdot;
    double c2pk,s2pk;
    double uk;
    double ukdot;
    double cuk,suk;
    double ok;
    double sok,cok;
    double ik;
    double ikdot;
    double sik,cik;
    double rk;
    double rkdot;
    double xpk,ypk;
    double xpkdot,ypkdot;

    double relativistic, OneMinusecosE, tmp;

    tk = g.sec - toe.sec;

    if(tk>SECONDS_IN_HALF_WEEK)
        tk -= SECONDS_IN_WEEK;
    else if(tk<-SECONDS_IN_HALF_WEEK)
        tk += SECONDS_IN_WEEK;

    mk = m0 + n*tk;
    ek = mk;
    ekold = ek + 1.0;

    OneMinusecosE = 0; // Suppress the uninitialized warning.
    while(fabs(ek-ekold)>1.0E-14)
    {
        ekold = ek;
        OneMinusecosE = 1.0-ecc*cos(ekold);
        ek = ek + (mk-ekold+ecc*sin(ekold))/OneMinusecosE;
    }

    sek = sin(ek);
    cek = cos(ek);

    ekdot = n/OneMinusecosE;

    relativistic = -4.442807633E-10*ecc*sqrta*sek;

    pk = atan2(sq1e2*sek,cek-ecc) + aop;
    pkdot = sq1e2*ekdot/OneMinusecosE;

    s2pk = sin(2.0*pk);
    c2pk = cos(2.0*pk);

    uk = pk + cus*s2pk + cuc*c2pk;
    suk = sin(uk);
    cuk = cos(uk);
    ukdot = pkdot*(1.0 + 2.0*(cus*c2pk - cuc*s2pk));

    rk = A*OneMinusecosE + crc*c2pk + crs*s2pk;
    rkdot = A*ecc*sek*ekdot + 2.0*pkdot*(crs*c2pk - crc*s2pk);

    ik = inc0 + idot*tk + cic*c2pk + cis*s2pk;
    sik = sin(ik);
    cik = cos(ik);
    ikdot = idot + 2.0*pkdot*(cis*c2pk - cic*s2pk);

    xpk = rk*cuk;
    ypk = rk*suk;
    xpkdot = rkdot*cuk - ypk*ukdot;
    ypkdot = rkdot*suk + xpk*ukdot;

    ok = omg0 + tk*omgkdot - OMEGA_EARTH*toe.sec;
    sok = sin(ok);
    cok = cos(ok);

    pos[0] = xpk*cok - ypk*cik*sok;
    pos[1] = xpk*sok + ypk*cik*cok;
    pos[2] = ypk*sik;

    tmp = ypkdot*cik - ypk*sik*ikdot;

    vel[0] = -omgkdot*pos[1] + xpkdot*cok - tmp*sok;
    vel[1] = omgkdot*pos[0] + xpkdot*sok + tmp*cok;
    vel[2] = ypk*cik*ikdot + ypkdot*sik;

    // Satellite clock correction
    tk = g.sec - toc.sec;

    if(tk>SECONDS_IN_HALF_WEEK)
        tk -= SECONDS_IN_WEEK;
    else if(tk<-SECONDS_IN_HALF_WEEK)
        tk += SECONDS_IN_WEEK;

    clk[0] = af0 + tk*(af1 + tk*af2) + relativistic - tgd;
    clk[1] = af1 + 2.0*tk*af2;

    return;
}

int Ephemeris::CheckSatVisibility(const GpsTime& g, double *xyz, double elvMask, double *azel) const
{
    double llh[3],neu[3];
    double pos[3],vel[3],clk[3],los[3];
    double tmat[3][3];

    if (vflg != 1)
        return (-1); // Invalid

    xyz2llh(xyz,llh);
    ltcmat(llh, tmat);

    Satpos(g, pos, vel, clk);
    subVect(los, pos, xyz);
    ecef2neu(los, tmat, neu);
    neu2azel(azel, neu);

    if (azel[1]*R2D > elvMask)
        return (1); // Visible
    // else
    return (0); // Invisible
}

