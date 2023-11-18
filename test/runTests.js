const nzmg = require("../index");

/**
 * The following output is what returns by running the original C code.
 * Tests in this file will use this in tests validation.
 *
 * echo "1783295.000  5868193.000" | ./nztm
 *
 * Enter NZTM easting, northing:
 * Input NZTM e,n:   1783295.000  5868193.000
 * Output Lat/Long:   -37.314852   175.068489
 * Output NZTM e,n:  1783295.000  5868193.000
 * Difference:             0.000       -0.000
 *
 *
 * echo "1375175 5086098" | ./nztm
 * Enter NZTM easting, northing:
 * Input NZTM e,n:   1375175.000  5086098.000
 * Output Lat/Long:   -44.343561   170.179492
 * Output NZTM e,n:  1375175.000  5086098.000
 * Difference:            -0.000       -0.000
 */
console.log("Test started.");
const northingEasting = [
  {
    easting: 1783295,
    northing: 5868193,
  },
  {
    easting: 1375175,
    northing: 5086098,
  },
];

for (let ne of northingEasting) {
  console.log(
    "-------------------------------------------------------------------------"
  );
  console.log(`Original: (easting, northing) =(${ne.easting}, ${ne.northing})`);
  const resp = nzmg.nztm_geod(ne.easting, ne.northing);
  console.log(
    `Converted: (latitude, longitude) =(${resp.latitude}, ${resp.longitude})`
  );
  const resp1 = nzmg.geod_nztm(resp.latitude, resp.longitude);
  console.log(
    `Converted: (easting, northing) = (${resp1.easting}, ${resp1.northing})`
  );

  if (resp.easting == resp1.easting && resp.northing == resp1.northing) {
    console.log(
      `NZTM coordinate { easting: ${resp.easting}, northing: ${resp.northing}} has been coverted successfully.`
    );
  } else {
    console.log(
      `NZTM coordinate { easting: ${resp.easting}, northing: ${resp.northing}} has failed.`
    );
  }
}
console.log(
  "-------------------------------------------------------------------------"
);
console.log("Test Completed");
