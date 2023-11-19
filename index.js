const Decimal = require("decimal.js");

/**
 * Define the parameters for the International Ellipsoid
 * used for the NZGD2000 datum (and hence for NZTM).
 * */
const NZTM_A = Decimal("6378137");
const NZTM_RF = Decimal("298.257222101");
const NZTM_CM = Decimal("173.0");
const NZTM_OLAT = Decimal("0.0");
const NZTM_SF = Decimal("0.9996");
const NZTM_FE = Decimal("1600000.0");
const NZTM_FN = Decimal("10000000.0");

const PI = Decimal("3.1415926535898");
const TWOPI = PI.times(Decimal("2.0"));
const rad2deg = Decimal("180").div(PI);

// six decimal points for all values
const FIXED_DP = 6;
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
const meridian_arc = (tm, lt) => {
  const e2 = tm.e2;
  const a = tm.a;
  const e4 = e2.times(e2);
  const e6 = e4.times(e2);

  const A0 = Decimal("1")
    .minus(e2.div(Decimal("4.0")))
    .minus(Decimal("3.0").times(e4).div(Decimal("64.0")))
    .minus(Decimal("5.0").times(e6).div(Decimal("256.0")));
  const A2 = Decimal("3.0")
    .div(Decimal("8.0"))
    .times(
      e2
        .plus(e4.div(Decimal("4.0")))
        .plus(Decimal("15.0").times(e6).div(Decimal("128.0")))
    );
  const A4 = Decimal("15.0")
    .div(Decimal("256.0"))
    .times(e4.plus(Decimal("3.0").times(e6).div(Decimal("4.0"))));
  const A6 = Decimal("35.0").times(e6).div(Decimal("3072.0"));

  let tmp = A0.times(lt);
  tmp = tmp.minus(A2.times(Decimal("2.0").times(lt).sin()));
  tmp = tmp.plus(A4.times(Decimal("4.0").times(lt).sin()));
  tmp = tmp.minus(A6.times(Decimal("6.0").times(lt).sin()));
  return a.times(tmp);
};

/* Initiallize the TM structure  */
const define_tmprojection = (a, rf, cm, sf, lto, fe, fn, utom) => {
  const tm = {};
  let f;
  tm.meridian = cm;
  tm.scalef = sf;
  tm.orglat = lto;
  tm.falsee = fe;
  tm.falsen = fn;
  tm.utom = utom;
  if (rf.toNumber() != 0.0) {
    f = Decimal("1.0").div(rf);
  } else {
    f = Decimal("0.0");
  }
  tm.a = a;
  tm.rf = rf;
  tm.f = f;
  tm.e2 = f.times(Decimal("2.0")).minus(f.times(f));
  tm.ep2 = tm.e2.div(Decimal("1.0").minus(tm.e2));
  tm.om = meridian_arc(tm, tm.orglat);

  return tm;
};

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
const foot_point_lat = (tm, m) => {
  const f = tm.f;
  const a = tm.a;

  const n = f.div(Decimal("2.0").minus(f));
  const n2 = n.times(n);
  const n3 = n2.times(n);
  const n4 = n2.times(n2);

  let tmp = Decimal("1.0")
    .plus(Decimal("9.0").times(n2).div(Decimal("4.0")))
    .plus(Decimal("225.0").times(n4).div(Decimal("64.0")));
  const g = a
    .times(Decimal("1.0").minus(n))
    .times(Decimal("1.0").minus(n2))
    .times(tmp);
  const sig = m.div(g);

  let phio = sig;
  phio = phio.plus(
    Decimal("3.0")
      .times(n)
      .div(Decimal("2.0"))
      .minus(Decimal("27.0").times(n3).div(Decimal("32.0")))
      .times(Decimal("2.0").times(sig).sin())
  );
  phio = phio.plus(
    Decimal("21.0")
      .times(n2)
      .div(Decimal("16.0"))
      .minus(Decimal("55.0").times(n4).div(Decimal("32.0")))
      .times(Decimal("4.0").times(sig).sin())
  );
  phio = phio.plus(
    Decimal("151.0")
      .times(n3)
      .div(Decimal("96.0"))
      .times(Decimal("6.0").times(sig).sin())
  );
  phio = phio.plus(
    Decimal("1097.0")
      .times(n4)
      .div(Decimal("512.0"))
      .times(Decimal("8.0").times(sig).sin())
  );

  return phio;
};

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
const tm_geod = (tm, ce, cn) => {
  const fn = tm.falsen;
  const fe = tm.falsee;
  const sf = tm.scalef;
  const e2 = tm.e2;
  const a = tm.a;
  const cm = tm.meridian;
  const om = tm.om;
  const utom = tm.utom;

  const cn1 = cn.minus(fn).times(utom).div(sf).plus(om);
  const fphi = foot_point_lat(tm, cn1);
  const slt = fphi.sin();
  const clt = fphi.cos();

  const eslt = Decimal("1.0").minus(e2.times(slt).times(slt));
  const eta = a.div(eslt.sqrt());
  const rho = eta.times(Decimal("1.0").minus(e2)).div(eslt);
  const psi = eta.div(rho);

  const E = ce.minus(fe).times(utom);
  const x = E.div(eta.times(sf));
  const x2 = x.times(x);

  const t = slt.div(clt);
  const t2 = t.times(t);
  const t4 = t2.times(t2);

  let trm1 = Decimal("1.0").div(Decimal("2.0"));
  let trm2 = Decimal("-4.0")
    .times(psi)
    .plus(Decimal("9.0").times(Decimal("1").minus(t2)));
  trm2 = trm2.times(psi);
  trm2 = trm2.plus(Decimal("12.0").times(t2)).div(Decimal("24.0"));
  let trm3 = Decimal("8.0")
    .times(Decimal("11.0").minus(Decimal("24.0").times(t2)))
    .times(psi);
  trm3 = trm3
    .minus(
      Decimal("12.0").times(Decimal("21.0").minus(Decimal("71.0").times(t2)))
    )
    .times(psi);
  trm3 = trm3
    .plus(
      Decimal("15.0").times(
        Decimal("15.0").times(t2).minus(Decimal("98.0")).times(t2)
      )
    )
    .times(psi);
  trm3 = trm3
    .plus(
      Decimal("180.0").times(
        Decimal("-3.0").times(t2).plus(Decimal("5.0").times(t2))
      )
    )
    .times(psi);
  trm3 = trm3.plus(Decimal("360.0").times(t4)).div(Decimal("720.0"));

  let trm4 = Decimal("1575.0")
    .times(t2)
    .plus(Decimal("4095.0"))
    .times(t2)
    .plus(Decimal("3633.0"))
    .times(t2)
    .plus(Decimal("1385.0"))
    .div(Decimal("720.0"));
  const latitude1 = t.times(x).times(E).div(sf.times(rho));
  const latitude2 = trm4
    .times(x2)
    .minus(trm3)
    .times(x2)
    .plus(trm2)
    .times(x2)
    .minus(trm1);
  const latitude = fphi.plus(latitude1.times(latitude2));

  trm1 = Decimal("1.0");
  trm2 = psi.plus(Decimal("2.0").times(t2)).div(Decimal("6.0"));
  trm3 = Decimal("-4.0")
    .times(Decimal("1.0").minus(Decimal("6.0").times(t2)))
    .times(psi);
  trm3 = trm3.plus(Decimal("9.0").minus(Decimal("68.0").times(t2))).times(psi);
  trm3 = trm3.plus(Decimal("72.0").times(t2)).times(psi);
  trm3 = trm3.plus(Decimal("24.0").times(t4)).div(Decimal("120.0"));
  trm4 = Decimal("720.0")
    .times(t2)
    .plus(Decimal("1320.0"))
    .times(t2)
    .plus(Decimal("662.0"))
    .times(t2)
    .plus(Decimal("61.0"))
    .div(Decimal("5040.0"));

  const longitude1 = x.div(clt);
  const longitude2 = trm4
    .times(x2)
    .minus(trm3)
    .times(x2)
    .plus(trm2)
    .times(x2)
    .minus(trm1);
  const longitude = cm.minus(longitude1.times(longitude2));

  return { latitude, longitude };
};

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
const geod_tm = (tm, ln, lt) => {
  const fn = tm.falsen;
  const fe = tm.falsee;
  const sf = tm.scalef;
  const e2 = tm.e2;
  const a = tm.a;
  const cm = tm.meridian;
  const om = tm.om;
  const utom = tm.utom;

  let dlon = ln.minus(cm);
  while (dlon.toNumber() > PI.toNumber()) {
    dlon = dlon.minus(TWOPI);
  }
  while (dlon.toNumber() < PI.neg().toNumber()) {
    dlon = dlon.plus(TWOPI);
  }
  const m = meridian_arc(tm, lt);
  const slt = lt.sin();
  const eslt = Decimal("1.0").minus(e2.times(slt.times(slt)));
  const eta = a.div(eslt.sqrt());
  const rho = eta.times(Decimal("1.0").minus(e2)).div(eslt);
  const psi = eta.div(rho);

  const clt = lt.cos();
  const w = dlon.times(Decimal("1.0"));

  const wc = clt.times(w);
  const wc2 = wc.times(wc);

  const t = slt.div(clt);
  const t2 = t.times(t);
  const t4 = t2.times(t2);
  const t6 = t2.times(t4);

  let trm1, trm2, trm3, trm4;
  trm1 = psi.minus(t2).div(Decimal("6.0"));
  trm2 = Decimal("4.0")
    .times(Decimal("1.0").minus(Decimal("6.0").times(t2)))
    .times(psi);
  trm2 = trm2.plus(Decimal("1.0")).plus(Decimal("8.0").times(t2)).times(psi);
  trm2 = trm2.minus(Decimal("2.0").times(t2)).times(psi);
  trm2 = trm2.plus(t4).div(Decimal("120.0"));
  trm3 = Decimal("61")
    .minus(Decimal("479.0").times(t2))
    .plus(Decimal("179.0").times(t4))
    .minus(t6)
    .div(Decimal("5040.0"));

  const gce1 = sf.times(eta).times(dlon).times(clt);
  const gce2 = trm3
    .times(wc2)
    .plus(trm2)
    .times(wc2)
    .plus(trm1)
    .times(wc2)
    .plus(Decimal("1.0"));
  const easting = gce1.times(gce2).div(utom).plus(fe);

  trm1 = Decimal("1.0").div(Decimal("2.0"));
  trm2 = Decimal("4.0")
    .times(psi)
    .plus(Decimal("1"))
    .times(psi)
    .minus(t2)
    .div(Decimal("24.0"));
  trm3 = Decimal("8.0")
    .times(Decimal("11.0").minus(Decimal("24.0").times(t2)))
    .times(psi);
  trm3 = trm3
    .minus(
      Decimal("28.0").times(Decimal("1.0").minus(Decimal("6.0").times(t2)))
    )
    .times(psi);
  trm3 = trm3.plus(Decimal("1.0").minus(Decimal("32.0").times(t2))).times(psi);
  trm3 = trm3.minus(Decimal("2.0").times(t2)).times(psi);
  trm3 = trm3.plus(t4).div(Decimal("720.0"));
  trm4 = Decimal("1385.0")
    .minus(Decimal("3111.0").times(t2))
    .plus(Decimal("543.0").times(t4))
    .minus(t6)
    .div(Decimal("40320.0"));
  const gcn1 = eta.times(t);
  const gcn2 = trm4
    .times(wc2)
    .plus(trm3)
    .times(wc2)
    .plus(trm2)
    .times(wc2)
    .plus(trm1)
    .times(wc2);
  const gcn = gcn1.times(gcn2);
  const northing = gcn.plus(m).minus(om).times(sf).div(utom).plus(fn);

  return { northing, easting };
};

/* Global that will be initialsed just one for efficiency */
let nztm_projection = null;
const get_nztm_projection = () => {
  if (!nztm_projection) {
    nztm_projection = define_tmprojection(
      NZTM_A,
      NZTM_RF,
      NZTM_CM.div(rad2deg),
      NZTM_SF,
      NZTM_OLAT.div(rad2deg),
      NZTM_FE,
      NZTM_FN,
      Decimal("1.0")
    );
  }
  return nztm_projection;
};

/* Functions implementation the TM projection specifically for the
   NZTM coordinate system
*/

/**
 *
 * @param {northing in meters} - n
 * @param {easting  in meters} - e
 * @returns {longitude, latitude} - latitude and Longitude are in degree
 */
exports.nztm_geod = (e, n) => {
  let decN = Decimal(`${n}`);
  let decE = Decimal(`${e}`);

  const tm = get_nztm_projection();
  const result = tm_geod(tm, decE, decN);

  return {
    northing: n,
    easting: e,
    longitude: Number(result.longitude.times(rad2deg).toFixed(FIXED_DP)),
    latitude: Number(result.latitude.times(rad2deg).toFixed(FIXED_DP)),
  };
};

/**
 *
 * @param {latitude in degree} lt
 * @param {longitude in degree} ln
 * @returns {northing, easting} - northing and easting in meters
 */
exports.geod_nztm = (lt, ln) => {
  let decLT = Decimal(`${lt}`).div(rad2deg);
  let decLN = Decimal(`${ln}`).div(rad2deg);
  const tm = get_nztm_projection();
  const result = geod_tm(tm, decLN, decLT);
  return {
    latitude: lt,
    longitude: ln,
    northing: result.northing.toFixed(0),
    easting: result.easting.toFixed(0),
  };
};
