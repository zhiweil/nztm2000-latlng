#include "nztm.h"

#include <math.h>

#define TEST_NZTM
#ifdef TEST_NZTM
#include <stdio.h>
#endif

/* Defines PI (from Abramowitz and Stegun Table 1.1) */
#define PI 3.1415926535898
#define TWOPI (2.0 * PI)
#define rad2deg (180 / PI)

/* Structure used to define a TM projection */

typedef struct
{
   double meridian; /* Central meridian */
   double scalef;   /* Scale factor */
   double orglat;   /* Origin latitude */
   double falsee;   /* False easting */
   double falsen;   /* False northing */
   double utom;     /* Unit to metre conversion */

   double a, rf, f, e2, ep2; /* Ellipsoid parameters */
   double om;                /* Intermediate calculation */
} tmprojection;

static double meridian_arc(tmprojection *tm, double lt);

/* Initiallize the TM structure  */

static void define_tmprojection(tmprojection *tm, double a, double rf,
                                double cm, double sf, double lto, double fe, double fn, double utom)
{

   double f;

   tm->meridian = cm;
   tm->scalef = sf;
   tm->orglat = lto;
   tm->falsee = fe;
   tm->falsen = fn;
   tm->utom = utom;
   if (rf != 0.0)
      f = 1.0 / rf;
   else
      f = 0.0;
   tm->a = a;
   tm->rf = rf;
   tm->f = f;
   tm->e2 = 2.0 * f - f * f;
   tm->ep2 = tm->e2 / (1.0 - tm->e2);

   tm->om = meridian_arc(tm, tm->orglat);
}

/***************************************************************************/
/*                                                                         */
/*  meridian_arc                                                           */
/*                                                                         */
/*  Returns the length of meridional arc (Helmert formula)                 */
/*  Method based on Redfearn's formulation as expressed in GDA technical   */
/*  manual at http://www.anzlic.org.au/icsm/gdatm/index.html               */
/*                                                                         */
/*  Parameters are                                                         */
/*    projection                                                           */
/*    latitude (radians)                                                   */
/*                                                                         */
/*  Return value is the arc length in metres                               */
/*                                                                         */
/***************************************************************************/

static double meridian_arc(tmprojection *tm, double lt)
{
   double e2 = tm->e2;
   double a = tm->a;
   double e4;
   double e6;
   double A0;
   double A2;
   double A4;
   double A6;

   e4 = e2 * e2;
   e6 = e4 * e2;

   A0 = 1 - (e2 / 4.0) - (3.0 * e4 / 64.0) - (5.0 * e6 / 256.0);
   A2 = (3.0 / 8.0) * (e2 + e4 / 4.0 + 15.0 * e6 / 128.0);
   A4 = (15.0 / 256.0) * (e4 + 3.0 * e6 / 4.0);
   A6 = 35.0 * e6 / 3072.0;
   return a * (A0 * lt - A2 * sin(2 * lt) + A4 * sin(4 * lt) - A6 * sin(6 * lt));
}

/*************************************************************************/
/*                                                                       */
/*   foot_point_lat                                                      */
/*                                                                       */
/*   Calculates the foot point latitude from the meridional arc          */
/*   Method based on Redfearn's formulation as expressed in GDA technical*/
/*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html            */
/*                                                                       */
/*   Takes parameters                                                    */
/*      tm definition (for scale factor)                                 */
/*      meridional arc (metres)                                          */
/*                                                                       */
/*   Returns the foot point latitude (radians)                           */ /*                                                                       */
/*************************************************************************/

static double foot_point_lat(tmprojection *tm, double m)
{
   double f = tm->f;
   double a = tm->a;
   double n;
   double n2;
   double n3;
   double n4;
   double g;
   double sig;
   double phio;

   n = f / (2.0 - f);
   n2 = n * n;
   n3 = n2 * n;
   n4 = n2 * n2;

   g = a * (1.0 - n) * (1.0 - n2) * (1 + 9.0 * n2 / 4.0 + 225.0 * n4 / 64.0);
   sig = m / g;

   phio = sig + (3.0 * n / 2.0 - 27.0 * n3 / 32.0) * sin(2.0 * sig) + (21.0 * n2 / 16.0 - 55.0 * n4 / 32.0) * sin(4.0 * sig) + (151.0 * n3 / 96.0) * sin(6.0 * sig) + (1097.0 * n4 / 512.0) * sin(8.0 * sig);

   return phio;
}

/***************************************************************************/
/*                                                                         */
/*   tmgeod                                                                */
/*                                                                         */
/*   Routine to convert from Tranverse Mercator to latitude and longitude. */
/*   Method based on Redfearn's formulation as expressed in GDA technical  */
/*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html              */
/*                                                                         */
/*   Takes parameters                                                      */
/*      input easting (metres)                                             */
/*      input northing (metres)                                            */
/*      output latitude (radians)                                          */
/*      output longitude (radians)                                         */
/*                                                                         */
/***************************************************************************/

static void tm_geod(tmprojection *tm,
                    double ce, double cn, double *ln, double *lt)
{
   double fn = tm->falsen;
   double fe = tm->falsee;
   double sf = tm->scalef;
   double e2 = tm->e2;
   double a = tm->a;
   double cm = tm->meridian;
   double om = tm->om;
   double utom = tm->utom;
   double cn1;
   double fphi;
   double slt;
   double clt;
   double eslt;
   double eta;
   double rho;
   double psi;
   double E;
   double x;
   double x2;
   double t;
   double t2;
   double t4;
   double trm1;
   double trm2;
   double trm3;
   double trm4;

   cn1 = (cn - fn) * utom / sf + om;
   fphi = foot_point_lat(tm, cn1);
   slt = sin(fphi);
   clt = cos(fphi);

   eslt = (1.0 - e2 * slt * slt);
   eta = a / sqrt(eslt);
   rho = eta * (1.0 - e2) / eslt;
   psi = eta / rho;

   E = (ce - fe) * utom;
   x = E / (eta * sf);
   x2 = x * x;

   t = slt / clt;
   t2 = t * t;
   t4 = t2 * t2;

   trm1 = 1.0 / 2.0;

   trm2 = ((-4.0 * psi + 9.0 * (1 - t2)) * psi + 12.0 * t2) / 24.0;

   trm3 = ((((8.0 * (11.0 - 24.0 * t2) * psi - 12.0 * (21.0 - 71.0 * t2)) * psi + 15.0 * ((15.0 * t2 - 98.0) * t2 + 15)) * psi + 180.0 * ((-3.0 * t2 + 5.0) * t2)) * psi + 360.0 * t4) / 720.0;

   trm4 = (((1575.0 * t2 + 4095.0) * t2 + 3633.0) * t2 + 1385.0) / 40320.0;

   *lt = fphi + (t * x * E / (sf * rho)) * (((trm4 * x2 - trm3) * x2 + trm2) * x2 - trm1);

   trm1 = 1.0;

   trm2 = (psi + 2.0 * t2) / 6.0;

   trm3 = (((-4.0 * (1.0 - 6.0 * t2) * psi + (9.0 - 68.0 * t2)) * psi + 72.0 * t2) * psi + 24.0 * t4) / 120.0;

   trm4 = (((720.0 * t2 + 1320.0) * t2 + 662.0) * t2 + 61.0) / 5040.0;

   *ln = cm - (x / clt) * (((trm4 * x2 - trm3) * x2 + trm2) * x2 - trm1);
}

/***************************************************************************/
/*                                                                         */
/*   geodtm                                                                */
/*                                                                         */
/*   Routine to convert from latitude and longitude to Transverse Mercator.*/
/*   Method based on Redfearn's formulation as expressed in GDA technical  */
/*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html              */
/*   Loosely based on FORTRAN source code by J.Hannah and A.Broadhurst.    */
/*                                                                         */
/*   Takes parameters                                                      */
/*      input latitude (radians)                                           */
/*      input longitude (radians)                                          */
/*      output easting  (metres)                                           */
/*      output northing (metres)                                           */
/*                                                                         */
/***************************************************************************/

static void geod_tm(tmprojection *tm,
                    double ln, double lt, double *ce, double *cn)
{
   double fn = tm->falsen;
   double fe = tm->falsee;
   double sf = tm->scalef;
   double e2 = tm->e2;
   double a = tm->a;
   double cm = tm->meridian;
   double om = tm->om;
   double utom = tm->utom;
   double dlon;
   double m;
   double slt;
   double eslt;
   double eta;
   double rho;
   double psi;
   double clt;
   double w;
   double wc;
   double wc2;
   double t;
   double t2;
   double t4;
   double t6;
   double trm1;
   double trm2;
   double trm3;
   double gce;
   double trm4;
   double gcn;

   dlon = ln - cm;
   while (dlon > PI)
      dlon -= TWOPI;
   while (dlon < -PI)
      dlon += TWOPI;

   m = meridian_arc(tm, lt);

   slt = sin(lt);

   eslt = (1.0 - e2 * slt * slt);
   eta = a / sqrt(eslt);
   rho = eta * (1.0 - e2) / eslt;
   psi = eta / rho;

   clt = cos(lt);
   w = dlon;

   wc = clt * w;
   wc2 = wc * wc;

   t = slt / clt;
   t2 = t * t;
   t4 = t2 * t2;
   t6 = t2 * t4;

   trm1 = (psi - t2) / 6.0;

   trm2 = (((4.0 * (1.0 - 6.0 * t2) * psi + (1.0 + 8.0 * t2)) * psi - 2.0 * t2) * psi + t4) / 120.0;

   trm3 = (61 - 479.0 * t2 + 179.0 * t4 - t6) / 5040.0;

   gce = (sf * eta * dlon * clt) * (((trm3 * wc2 + trm2) * wc2 + trm1) * wc2 + 1.0);
   *ce = gce / utom + fe;

   trm1 = 1.0 / 2.0;

   trm2 = ((4.0 * psi + 1) * psi - t2) / 24.0;

   trm3 = ((((8.0 * (11.0 - 24.0 * t2) * psi - 28.0 * (1.0 - 6.0 * t2)) * psi + (1.0 - 32.0 * t2)) * psi - 2.0 * t2) * psi + t4) / 720.0;

   trm4 = (1385.0 - 3111.0 * t2 + 543.0 * t4 - t6) / 40320.0;

   gcn = (eta * t) * ((((trm4 * wc2 + trm3) * wc2 + trm2) * wc2 + trm1) * wc2);
   *cn = (gcn + m - om) * sf / utom + fn;

   return;
}

/* Define a static implementation of tmprojection */
/* Note: for some implementations it may be better to create this
   dynamically and develop modified versions of the transformation
   functions to take this as a parameter */

static tmprojection nztm_projection;
static int initiallized = 0;

static tmprojection *get_nztm_projection(void)
{
   if (!initiallized)
   {
      define_tmprojection(&nztm_projection, NZTM_A, NZTM_RF,
                          NZTM_CM / rad2deg, NZTM_SF, NZTM_OLAT / rad2deg, NZTM_FE, NZTM_FN,
                          1.0);
      initiallized = 1;
   }
   return &nztm_projection;
}

/* Functions implementation the TM projection specifically for the
   NZTM coordinate system
*/

void nztm_geod(double n, double e, double *lt, double *ln)
{
   tmprojection *nztm = get_nztm_projection();
   tm_geod(nztm, e, n, ln, lt);
}

void geod_nztm(double lt, double ln, double *n, double *e)
{
   tmprojection *nztm = get_nztm_projection();
   geod_tm(nztm, ln, lt, e, n);
}

#ifdef TEST_NZTM
int main(int argc, char *argv[])
{
   double e, n, lt, ln, e1, n1;
   while (1)
   {
      printf("Enter NZTM easting, northing: ");
      if (scanf("%lf%lf", &e, &n) != 2)
         break;
      nztm_geod(n, e, &lt, &ln);
      geod_nztm(lt, ln, &n1, &e1);

      printf("\nInput NZTM e,n:  %12.3lf %12.3lf\n", e, n);
      printf("Output Lat/Long: %12.6lf %12.6lf\n", lt * rad2deg, ln * rad2deg);
      printf("Output NZTM e,n: %12.3lf %12.3lf\n", e1, n1);
      printf("Difference:      %12.3lf %12.3lf\n\n", e1 - e, n1 - n);
   }
   return 0;
}
#endif
