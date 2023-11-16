#ifndef _NZMG_H
#define _NZMG_H

/* Define the parameters for the International Ellipsoid
   used for the NZGD49 datum (and hence for NZMG) */

#define NZMG_A  6378388.0
#define NZMG_RF 297.0

/* Routines to convert NZMG to latitude and longitude
   and vice versa.  Northing (n) and Easting (e) are in
   metres, Latitude (lt) and Longitude (ln) are in
   radians */

void nzmg_geod( double n, double e, double *lt, double *ln );
void geod_nzmg( double lt, double ln, double *n, double *e );


#endif