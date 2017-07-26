/***************************************************************************************/
/*                                      Instru.h

    This is the header file for Instru.cpp, which contains attenuation tables
    pre-calculated for NICER and instrument response curve. 

*/
/***************************************************************************************/

#define NDIM 5  //
#define MPTS 6  //

#define MAX_NUMBINS 256 // REMEMBER TO CHANGE THIS IN STRUCT.H AS WELL!!
#define NCURVES 351      // REMEMBER TO CHANGE THIS IN STRUCT.H AS WELL!! number of different light curves that it will calculate

//

class LightCurve Attenuate (class LightCurve* incurve, unsigned int attenuation, double nh);

class LightCurve Inst_Res (class LightCurve* incurve, unsigned int inst_curve);

class LightCurve Inst_Res2 (class LightCurve* incurve, unsigned int inst_curve);


class LightCurve Background_list (class LightCurve* incurve, char *background_file);

class LightCurve AGN_Background (class LightCurve* incurve, double agnbackground, double nh);

class LightCurve Sky_Background (class LightCurve* incurve, double skybackground);

class LightCurve PowerLaw_Background (class LightCurve* incurve, double agnbackground, double nh);
