#ifndef _NZTM_H
#define _NZTM_H

/* Define the parameters for the International Ellipsoid
   used for the NZGD2000 datum (and hence for NZTM) */

#define NZTM_A  6378137 
#define NZTM_RF 298.257222101

#define NZTM_CM    173.0 
#define NZTM_OLAT    0.0 
#define NZTM_SF      0.9996 
#define NZTM_FE    1600000.0 
#define NZTM_FN    10000000.0

/* Routines to convert NZTM to latitude and longitude
   and vice versa.  Northing (n) and Easting (e) are in
   metres, Latitude (lt) and Longitude (ln) are in
   radians */

void nztm_geod( double n, double e, double *lt, double *ln );
void geod_nztm( double lt, double ln, double *n, double *e );


#endif
