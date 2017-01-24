/***************************************************************************************/
/*                                      Instru.h

    This is the header file for Instru.cpp, which contains attenuation tables
    pre-calculated for NICER and instrument response curve. 

*/
/***************************************************************************************/

#define NDIM 5  //
#define MPTS 6  //

#define MAX_NUMBINS 128 // REMEMBER TO CHANGE THIS IN STRUCT.H AS WELL!!
#define NCURVES 50      // REMEMBER TO CHANGE THIS IN STRUCT.H AS WELL!! number of different light curves that it will calculate

//

class LightCurve Attenuate (class LightCurve* incurve, unsigned int attenuation);

class LightCurve Inst_Res (class LightCurve* incurve, unsigned int inst_curve);

class LightCurve Background_list (class LightCurve* incurve, char *background_file);

