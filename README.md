# nztm2000-latlng

Convert values between NZTM2000 and latitude/longitude. The implementation is a direct translation of the C source code found on [this LINZ website](https://www.linz.govt.nz/products-services/geodetic/geodetic-software-downloads#nztm2000). The precision of this module is also comparable to that of the origional C implementation.

The original C files can be found in the "original" folder of [this GIT repository](https://www.npmjs.com/package/nztm2000-latlng).

## 1. Methods

### 1.1 tm_geod

Routine to convert from New Zeland Tranverse Mercator (NZTM2000) to latitude and longitude. Six decimal point places are kept for the converted latititude and langitude.

**Parameters**

```javascript
tm_geod(easting, northing); // easting and northig are in meters;
```

**Return**

```json
{
  "easting": 1783295, // in meters, whole number
  "northing": 5868193, // in meters, whole number
  "latitude": -37.314852, // in degree, six decimal places
  "longitude": 175.068489 // in degree, six decimal places
}
```

### 1.2 geod_tm

Routine to convert from latitude and longitude to New Zeland Tranverse Mercator (NZTM2000). The converted easting and northing are whole numbers.

**Parameters**

```javascript
geod_tm(latitude, longitude); // latitude and longitue are in degrees
```

**Return**

```json
{
  "easting": 1783295, // in meters, whole number
  "northing": 5868193, // in meters, whole number
  "latitude": -37.314852, // in degree, six decimal places
  "longitude": 175.068489 // in degree, six decimal places
}
```

## 2. Test

The module replicated the C testing code as well by the following npm script.

```bash
npm run test
```

The conversion can be checked by [the LINZ online converter](https://www.geodesy.linz.govt.nz/concord/index.cgi?IS=NZTM&IH=-&PN=N&IC=H&IO=EN&ID=S&OS=NZGD2000&OH=-&OC=D&OO=EN&OD=S&YEAR=now&OP=2&do_entry=1).

## 3. Run C implementation

Please install gcc or any C compiler on the development host machine.

Complie the C implementation:

```bash
cd original
gcc -o nztm nztm.c
```

Run the code:

```bash
echo "1783295.000  5868193.000" | ./nztm
```
