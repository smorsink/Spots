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

  // NUM_NICER_CHANNELS is the number of NICER energy channels that will be read in.
  // There are more channels than this in the response matrix

  
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

      if ( nicer->start[p] != 0 )
	std::cout << "WARNING: start[" << p << "] = " << nicer->start[p] << std::endl;
      

	    
      // The number 300 comes from the structure of the response matrix file
      
      for (unsigned int j(0); j<=300; j++){
	file >> nicer->response[p][j];

	      
      }
    }
  }else{
    throw( Exception( "instrument response curve file is not found" ));
  }
  file.close();
}



// Read in a light curve computed for numbands different energy bands.
// Interpolate to find the same light curve computed for the NICER energy channels
class NICERCurve ConvertEnergyChannels(class LightCurve* incurve, class Instrument* nicer){


  class NICERCurve newcurve;
  double energy, energyfactor;

  
      for (unsigned int p(0);p<NUM_NICER_CHANNELS;p++){ // loop through the NICER energy channels

	energy = nicer->elo[p];

	double factor = (energy - incurve->elo[0])/(incurve->ehi[0] - incurve->elo[0]);
	int index = factor;


	
	// Interpolate to find the flux in the NICER energy channel

	newcurve.elo[p] = nicer->elo[p];
	newcurve.ehi[p] = nicer->ehi[p];

	energyfactor = (energy - incurve->elo[index])/(incurve->elo[index+1]-incurve->elo[index]);

	for (unsigned int i(0); i<=incurve->numbins; i++){
	  // linear interpolation for each timebin	 	  
	  newcurve.f[p][i] = incurve->f[index][i] + (incurve->f[index+1][i] - incurve->f[index][i]) * energyfactor;

	  newcurve.f[p][i] *= (1.0)/(incurve->elo[index+1]-incurve->elo[index]);

	  	  
	}	
      }

      // newcurve now has NUM_NICER_CHANNELS energy bands!
      newcurve.numbands = NUM_NICER_CHANNELS;
      
      return newcurve;
  
}


class NICERCurve ApplyResponse(class NICERCurve* incurve, class Instrument* nicer){

  class NICERCurve curve, newcurve;
  double factor;
  
  // std::cout << "Entered ApplyResponse!" << std::endl;
  
  curve = (*incurve);

  // Interpolate computed flux in each energy band to find the flux in the nicer bands

  // Apply response matrix to the interpolated flux
  // It is important that the number of energy channels in "incurve" be the same as the number of channels used in the response matrix.



  for (unsigned int q(0);q<300;q++)
    for (unsigned int i(0);i<curve.numbins;i++)
      newcurve.f[q][i] = 0;
  

  
  
  for (unsigned int p(0); p<NUM_NICER_CHANNELS; p++){

    // p refers to the index of the original photon.
    // This photon gets mapped to channels labelled q
     
    for (unsigned int q(0); q < 300; q++){

      factor = nicer->response[p][q] * (nicer->ehi[p]-nicer->elo[p]);
      if (factor != 0){
	for( unsigned int i(0);i<curve.numbins;i++){	
	  newcurve.f[q][i] += curve.f[p][i] * factor;
	}
      }
    }
  }

  
  newcurve.numbands = 300;
  newcurve.numbins = curve.numbins;
  
  return newcurve;
  
  
}

