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
      
      /*if (p > 301)
	std::cout << std::endl << " p = " << p 
		  << " Photon Energy Range: " << nicer->elo[p]
		  << " to " << nicer->ehi[p] << " keV" << std::endl;*/
	    
      // The number 300 comes from the structure of the response matrix file
      
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



// Read in a light curve computed for numbands different energy bands.
// Interpolate to find the same light curve computed for the NICER energy channels
class LightCurve ConvertEnergyChannels(class LightCurve* incurve, class Instrument* nicer){

  std::cout << "ConvertEnergyChannels: We computed " << incurve->numbands << " number of bands" << std::endl;
  std::cout << "                       starting with " << incurve->elo[0] << " keV at intervals Delta(E) = "
	    << incurve->ehi[0] - incurve->elo[0] << " keV" << std::endl;

  class LightCurve newcurve;
  double energy, energyfactor;

  newcurve = (*incurve);
   
      for (unsigned int p(0);p<NUM_NICER_CHANNELS;p++){ // loop through the NICER energy channels

	/*if (p<3){
	  std::cout
	    << "Computed curve.elo[" << p <<"] = " << incurve->elo[p]
	    << " nicer.elo[p] = " << nicer->elo[p]
	    << std::endl;
	  std::cout
	    << "Computed curve.elhi[" << p <<"] = " << incurve->ehi[p]
	    << " nicer.ehi[p] = " << nicer->ehi[p]
	    << std::endl;	  
	    }*/

	energy = nicer->elo[p];

	double factor = (energy - incurve->elo[0])/(incurve->ehi[0] - incurve->elo[0]);
	int index = factor;

	/*	if (p>300){
	  
	  std::cout
	    <<" channel = " << p
	    << " NICER channel energy = " << energy 
	    << " factor = " << factor
	    << " index = " << index
	    << " closest energy[index] = " << incurve->elo[index]
	    << std::endl;
	    }*/
	  
	
	// Interpolate to find the flux in the NICER energy channel

	newcurve.elo[p] = nicer->elo[p];
	newcurve.ehi[p] = nicer->ehi[p];

	energyfactor = (energy - incurve->elo[index])/(incurve->elo[index+1]-incurve->elo[index]);

	for (unsigned int i(0); i<=newcurve.numbins; i++){
	  // linear interpolation for each timebin	 	  
	  newcurve.f[p][i] = incurve->f[index][i] + (incurve->f[index+1][i] - incurve->f[index][i]) * energyfactor;

	  newcurve.f[p][i] *= (1.0)/(incurve->elo[index+1]-incurve->elo[index]);

	  
	  /* if(i==0 && p<3)
	    std::cout
	    << "f1 = " << incurve->f[index][i] << " f2 = " << incurve->f[index+][i]*/
	  
	}	
      }

      // newcurve now has NUM_NICER_CHANNELS energy bands!
      newcurve.numbands = NUM_NICER_CHANNELS;
      
      return newcurve;
  
}


class LightCurve ApplyResponse(class LightCurve* incurve, class Instrument* nicer){

  class LightCurve curve, newcurve;
  double factor;
  
  std::cout << "Entered ApplyResponse!" << std::endl;
  
  curve = (*incurve);

  // Interpolate computed flux in each energy band to find the flux in the nicer bands

  // Apply response matrix to the interpolated flux
  // It is important that the number of energy channels in "incurve" be the same as the number of channels used in the response matrix.
  // std::cout << "NCURVES = " << NCURVES << std::endl;
  //std::cout << "curve.numbands = " << curve.numbands << std::endl;
  //std::cout << " elo = " << nicer->elo[0] << " ehi = " << nicer->ehi[0] << std::endl;

  
  for (unsigned int p(0); p<NUM_NICER_CHANNELS; p++){

    // p refers to the index of the original photon.
    // This photon gets mapped to channels labelled q

    //  if (p<10)
    //std::cout << " p = " << p << " flux[p][0] = " << curve.f[p][0] << std::endl;
    
    for (unsigned int q(0); q < 300; q++){

      factor = nicer->response[p][q] * (nicer->ehi[p]-nicer->elo[p]);
      //if ( q == 30 && p<10)
	//	std::cout << " p = " << p << "factor = " << factor << std::endl;
      if (factor != 0){
	for(unsigned int i(0);i<curve.numbins;i++){	
	  newcurve.f[q][i] += curve.f[p][i] * factor;
	}
      }
    }

  }
  newcurve.numbands = 300;
  newcurve.numbins = curve.numbins;

  //  std::cout << "Response[0][30] = " << nicer->response[0][30] << std::endl;
  // std::cout << "After Response: f[30][0] = " << newcurve.f[30][0] << std::endl;
  
  return newcurve;
  
  
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



