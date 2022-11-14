/***************************************************************************************/
/*                                      Chi.h

    This is the header file for Chi.cpp, which holds the computational methods/functions 
    called in Spot.cpp and probably elsewhere.
    
    And now it computes chi^2 too!
*/
/***************************************************************************************/

#define NDIM 5  //
#define MPTS 6  //

#define MAX_NUMBINS 128 // REMEMBER TO CHANGE THIS IN STRUCT.H AS WELL!!
//#define NCURVES 1      // REMEMBER TO CHANGE THIS IN STRUCT.H AS WELL!! number of different light curves that it will calculate
#define NCURVES 700


// Calculates chi^2
double ChiSquare( class DataStruct* obsdata, class LightCurve* curve);

double BackChi ( class DataStruct* obsdata, class LightCurve* curve, class LightCurve* backcurve);
double BandChi ( class DataStruct* obsdata, class LightCurve* curve, double back, unsigned int band);

class LightCurve SpotShape( int pieces, int p, int numtheta, double theta_1, double rho, class LightCurve* incurve,  class OblModelBase* modptr);

double SpotIntegrand( double rho, double zeta, class LightCurve* curve, class OblModelBase* model);

class LightCurve ReBinCurve(class DataStruct* obsdata, class LightCurve* curve);


// Calculates angles
class LightCurve ComputeAngles( class LightCurve* incurve,
				                class OblDeflectionTOA* defltoa );

class LightCurve Bend ( class LightCurve* incurve,
			class OblDeflectionTOA* defltoa);

class LightCurve ReadBend ( class LightCurve* incurve,
			    char *bend_file);


class LightCurve ShiftCurve( class LightCurve* angles, double phishift);



//STAY IN CHI.CPP
// Normalizes the light curve flux to 1
class LightCurve Normalize1( double Flux[NCURVES][MAX_NUMBINS], unsigned int numbins );

class LightCurve Normalize2( double Flux[NCURVES][MAX_NUMBINS], unsigned int numbins );


//Calculates Legendre polynomial P2 for equation 8, MLCB
double LegP2( double costheta );



// Calculates Legendre polynomial P4 for equation 8, MLCB
double LegP4( double costheta );

