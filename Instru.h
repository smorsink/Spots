/***************************************************************************************/
/*                                      Instru.h

    This is the header file for Instru.cpp, which accesses
    NICER response and area matrices.

*/
/***************************************************************************************/

#define NDIM 5  //
#define MPTS 6  //



class LightCurve Inst_Res (class LightCurve* incurve, unsigned int inst_curve);

class LightCurve Inst_Res2 (class LightCurve* incurve, unsigned int inst_curve, unsigned long int *start, double** response);

class LightCurve Inst_Res3 (class LightCurve* incurve, unsigned int inst_curve, unsigned long int *start, double** response);

void ReadResponse(class Instrument* nicer);

class LightCurve ConvertEnergyChannels(class LightCurve* incurve, class Instrument* nicer);

class LightCurve ApplyResponse(class LightCurve* incurve, class Instrument* nicer);
