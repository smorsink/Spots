/***************************************************************************************/
/*                                      Ism.h

    This is the header file for Ism.cpp 

*/
/***************************************************************************************/

#define NDIM 5  //
#define MPTS 6  //

//#define MAX_NUMBINS 256 // REMEMBER TO CHANGE THIS IN STRUCT.H AS WELL!!
//#define NCURVES 1      // REMEMBER TO CHANGE THIS IN STRUCT.H AS WELL!! number of different light curves that it will calculate
//#define NCURVES 300
//

class LightCurve Wabs (class LightCurve* incurve, unsigned int attenuation, double nh);

class LightCurve Attenuate (class LightCurve* incurve, class ISM* ism);

void ReadTBNEW(double nh, class ISM* ism);
