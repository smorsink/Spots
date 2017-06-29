/***************************************************************************************/
/*                                      Atmo.h

    This is the header file for Atmo.cpp, which holds ComputeCurve and atmosphere routines
    called in Spot.cpp and probably elsewhere.

    Was split from Chi.h

*/
/***************************************************************************************/

#define NDIM 5  //
#define MPTS 6  //

#define MAX_NUMBINS 256 // REMEMBER TO CHANGE THIS IN STRUCT.H AS WELL!!
#define NCURVES 301      // REMEMBER TO CHANGE THIS IN STRUCT.H AS WELL!! number of different light curves that it will calculate





// Calculates the light curve, when given all the angles
// MOVE TO atmo.h
class LightCurve ComputeCurve( class LightCurve* angles );


// Find value in array
int Find(double val, std::vector<double> array);

// Linear Interpolation
double Linear(double x,double Xa, double Ya, double Xb, double Yb);

// Linear Interpolation in log-log space
double LogLinear(double x,double Xa, double Ya, double Xb, double Yb);

// Linear Interpolation in an array
double Interpolate(double X_INT, std::vector<double> X, std::vector<double> Y);

// Linear Interpolation in am array and in log-log space
double LogInterpolate(double X_INT, std::vector<double> X, std::vector<double> Y);

// Round value to the nearest value in an array
int Round(int n, double z, std::vector<double> v);

// Reading NSATMOS hydrogen tables
void Read_NSATMOS(double T, double M, double R);

// Hydrogen
double Hydrogen(double E, double cos_theta);

// Hydrogen
double Hydrogen2(int E_dex, double cos_theta);

// Reading old NSX helium tables
void Read_NSX(double T, double M, double R);

// Helium
double Helium(double E, double cos_theta);

// Helium
double Helium2(int E_dex, double cos_theta);

// Reading NSX hydrogen table
void Read_NSXH(double T, double M, double R);

// NSX Hydrogen
double NSXH(double E, double cos_theta);

// NSX Hydrogen
double NSXH2(int E_dex, double cos_theta);

// Reading McPHAC Hydrogen table
void Read_McPHAC(double T, double M, double R);

// NSX Hydrogen
double McPHAC(double E, double cos_theta);

// NSX Hydrogen
double McPHAC2(int E_dex, double cos_theta);

// Reading new NSX helium table
void Read_NSXHe(double T, double M, double R);

// NSX Helium
double NSXHe(double E, double cos_theta);

// NSX Helium
double NSXHe2(int E_dex, double cos_theta);


// Cole's McPHAC
double McPHACC(double E, double cos_theta);

// Cole's McPHAC
double McPHACC2(int E_dex, double cos_theta);


double NSXHnew(double E, double cos_theta, double T, double lgrave, class LightCurve mexmcc);


// Cole's McPHAC in spotMex
double McPHACC4(int E_dex, double cos_theta, double T, double M, double R, class LightCurve mexmcc);


double Line( double T, double E , double E1, double E2);

double LineBandFlux( double T, double E1, double E2, double L1, double L2 );



// flux from a specific energy band, for atmosphere models
double AtmosEBandFlux( unsigned int model, double cos_theta, double E1, double E2 );

// flux from a specific energy band, for helium model at log-spaced energy points
double AtmosEBandFlux2( unsigned int model, double cos_theta, double E1, double E2 );

double AtmosEBandFlux4new( unsigned int model, double cos_theta, double T, double lgrav, double E1, double E2, class LightCurve mexmcc);


