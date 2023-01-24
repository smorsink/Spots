/***************************************************************************************/
/*                                    Instru.cpp

    This holds the instrument-specifc routines used in Spot.cpp. Contains 
    NICER response and area matrix routines.  

*/
/***************************************************************************************/

#include "matpack.h"
#include <exception>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unistd.h>
#include "Exception.h"
#include "Units.h"
#include "Struct.h"
#include "time.h"
#include "Instru.h"
#include "nrutil.h"
#include <stdio.h>
using namespace std;


// Read in the combined ARF and RMF response matrix
// Up to date NICER response for 2023

void ReadResponse(class Instrument* nicer){

  nicer->elo = dvector(0,NUM_NICER_CHANNELS);
  nicer->ehi = dvector(0,NUM_NICER_CHANNELS);
  nicer->start = lvector(0,NUM_NICER_CHANNELS);
  nicer->response = dmatrix(0,NUM_NICER_CHANNELS,0,NUM_NICER_CHANNELS);
  
  std::ifstream file;
  file.open("instrument/NICERarray50_rsp.txt");
	
  if(file.is_open()) {
    for (unsigned int p(0);p<NUM_NICER_CHANNELS;p++){ // was NCURVES
      file >> nicer->elo[p];
      file >> nicer->ehi[p];
      file >> nicer->start[p];

      /* if (p < 2)
	std::cout << std::endl << " p = " << p 
		  << " Photon Energy Range: " << nicer->elo[p]
		  << " to " << nicer->ehi[p] << " keV" << std::endl;*/
	    
	    
      for (unsigned int j(0); j<=300; j++){
	file >> nicer->response[p][j];
	/*if (p < 2)
	  std::cout << " r[" << p << "," << j << "]= " << nicer->response[p][j] ;*/
	      
      }
    }
  }else{
    throw( Exception( "instrument response curve file is not found" ));
  }
  file.close();
}

class LightCurve ApplyResponse(class LightCurve* incurve, class Instrument* nicer){

  class LightCurve curve;

  std::cout << "Entered ApplyResponse!" << std::endl;
  
  curve = (*incurve);

  // Interpolate computed flux in each energy band to find the flux in the nicer bands

  // Apply response matrix to the interpolated flux

  std::cout << "NCURVES = " << NCURVES << std::endl;
  
  for (unsigned int p(0); p<NCURVES; p++){

    //std::cout << "ApplyResponse: p = " << p << std::endl;

    if (p==2)
      std::cout
      << "curve.elo[" << p <<"] = " << curve.elo[p]
      << " nicer.elo[p] = " << nicer->elo[p]
	      << std::endl;


  }


  return curve;
  
  
}




class LightCurve Inst_Res2 (class LightCurve* incurve, unsigned int inst_curve, unsigned long int* start, double** response){

  class LightCurve curve, newcurve;
	unsigned int  numbands, numbins;


	//Pour numbands and numbins into newcurve
	curve = (*incurve);
	
	numbands = curve.numbands;
	numbins = curve.numbins;
	newcurve.numbands = numbands;
	newcurve.numbins = numbins;
	

	int newindex;

    for (unsigned int p =0; p<NCURVES; p++)
	for (unsigned int i=0; i<numbins;i++)
	 newcurve.f[p][i] = 0.0;


    for (unsigned int p = 0; p < NCURVES; p++){
      //if (p==299) std::cout << "p=" << p << "start[p] = " << start[p] << " curve[p][0]=" << curve.f[p][0] << std::endl;
      for (unsigned int j=0; j<=76; j++){   
	newindex = j + start[p] - 1;   
	//std::cout << "j="<< j << " newindex=" << newindex << std::endl;
	if ( newindex < NCURVES)
	  for (unsigned int i = 0; i < numbins; i++){
	    newcurve.f[newindex][i] += curve.f[p][i] * response[p][j];
	  }
      }

    }

    for(unsigned int p=0; p<NCURVES; p++){
      for (unsigned int i=0; i<numbins; i++){
	curve.f[p][i] = newcurve.f[p][i];
      }
    }

    return curve;
}

//class LightCurve Inst_Res3 (class LightCurve* incurve, unsigned int inst_curve){
class LightCurve Inst_Res3 (class LightCurve* incurve, unsigned int inst_curve, unsigned long int* start, double** response){

  class LightCurve curve, newcurve;
	unsigned int  numbands, numbins;

	//Pour numbands and numbins into newcurve
	curve = (*incurve);
	
	numbands = curve.numbands;
	numbins = curve.numbins;
	newcurve.numbands = numbands;
	newcurve.numbins = numbins;
	
	int newindex;

	
    for (unsigned int p =0; p<NCURVES; p++)
	for (unsigned int i=0; i<numbins;i++)
	 newcurve.f[p][i] = 0.0;


    for (unsigned int p = 0; p < NCURVES; p++){
      for (unsigned int j=0; j<=start[p]; j++){   
	  for (unsigned int i = 0; i < numbins; i++){
	    newcurve.f[j][i] += curve.f[p][i] * response[p][j];
	  }
	  /*if (p==0)
	    std::cout << "j=" << j
		      << " response[0][j] = " << response[p][j]
		      << " newflux[j][0] = " << newcurve.f[j][0]
		      << std::endl;*/
      }

    }

    for(unsigned int p=0; p<NCURVES; p++){
      for (unsigned int i=0; i<numbins; i++){
	curve.f[p][i] = newcurve.f[p][i];
      }
    }

    return curve;
}



