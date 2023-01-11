/***************************************************************************************/
/*                                      Instru.h

    This is the header file for Instru.cpp, which accesses
    NICER response and area matrices.

*/
/***************************************************************************************/

#define NDIM 5  //
#define MPTS 6  //

#define MAX_NUMBINS 256 // REMEMBER TO CHANGE THIS IN STRUCT.H AS WELL!!
//#define NCURVES 1      // REMEMBER TO CHANGE THIS IN STRUCT.H AS WELL!! number of different light curves that it will calculate
#define NCURVES 300
//


class LightCurve Inst_Res (class LightCurve* incurve, unsigned int inst_curve);

class LightCurve Inst_Res2 (class LightCurve* incurve, unsigned int inst_curve, unsigned long int *start, double** response);

class LightCurve Inst_Res3 (class LightCurve* incurve, unsigned int inst_curve, unsigned long int *start, double** response);

