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

double Attenuate (unsigned int p, double flux_before, unsigned int attenuation);

double Inst_Res (unsigned int p, double flux_before, unsigned int inst_curve);

double Background_list (unsigned int p, double flux_before, char *background_file);

