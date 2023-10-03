/***************************************************************************************/
/*                                    Ism.cpp

    Routines that take into account the attenuation by the ISM

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
#include "Atmo.h"
#include "OblDeflectionTOA.h"
#include "OblModelBase.h"
#include "PolyOblModelNHQS.h"
#include "Exception.h"
#include "Units.h"
#include "Struct.h"
#include "time.h"
#include "Ism.h"
#include "nrutil.h"
#include <stdio.h>
using namespace std;




class LightCurve Attenuate (class LightCurve* incurve, class ISM* tbnew){

	class LightCurve curve;
	unsigned int numbins;
	unsigned int numbands;

	double factor;
	int index;

	//	double attenuation;
	
	curve = (*incurve);
	numbins = curve.numbins;
	numbands = curve.numbands;

	// The light curve is computed at NCURVES observed energy values defined by
	// curve.elo[p]

	// The attenuation is defined for energies tbnew->energy

	//std::cout << "numbands = " << numbands << std::endl;

	for (unsigned int p = 0; p < numbands; p++){

	  //std::cout << "Attenuate: p = " << p << " computed energy = " << curve.elo[p] << std::endl;

	  factor = curve.elo[p]*1e2 - 10.0 + 0.5; //The 0.5 takes care of the rounding from double to int
	  index = factor;

	  //std::cout << "factor = " << factor
	  //<< " index = " << index << " NH energy = " << tbnew->energy[index]
	  //	    << " atten = " << tbnew->attenuation[index]
	  //	    << std::endl;
	  
	  
	  
	  for (unsigned int i = 0; i < numbins; i++){
        
	    curve.f[p][i] *= tbnew->attenuation[index]; 
        	
	  }
	}
	return curve;
}



void ReadTBNEW(double nh, class ISM* tbnew ){

  double **atten;
  atten = dmatrix(0,400,0,1191);

  int inh;

  //double *ee, *aa;

  tbnew->energy = dvector(0,1191);
  tbnew->attenuation = dvector(0,1191);

  // TBNEW file has entries for 1190 photon energies!

  
  // Read in the attenuation file

      std::cout << "Read in TBNEW!" << std::endl;
      std::ifstream ism;       // input stream
      ism.open("ISM/tbnew_full.txt");
      char line[265]; // line of the data file being read in

      std::cout << "nh = " << nh << "e18 cm^2" << std::endl;
      inh = nh;
      //std::cout << "inh = " << inh  << std::endl;
      
      double get_nh, get_e, get_att;
      for (unsigned int k(1); k< 400; k++){ // k is a label for the value of nH
	for (unsigned int i(0); i<1191; i++){ // i is a label for the photon energy
	  ism.getline(line,265);
	  //std::cout << "line = " << line << std::endl;
	  sscanf( line, "%lf %lf %lf", &get_nh, &get_e, &get_att );
	  if (k==1){
	    tbnew->energy[i] = get_e;
	    if (inh < 10 && inh > 0)
	      tbnew->attenuation[i] = get_att;  
	  }
	  else{

	    if (inh >= 10 && k == 2+(nh-10)/5)
	      tbnew->attenuation[i] = get_att;
	    
	  }
	  
	  atten[k][i] = get_att;
	
	  
	}
      }
     
     
   
      //   std::cout << "energy 0 = " << tbnew->energy[0] << " attenuation = " << tbnew->attenuation[0] << std::endl;
      //std::cout << "energy 1 = " << tbnew->energy[1] << " attenuation = " << tbnew->attenuation[1] << std::endl;

     

	// 20230117 FIX UP Energies

      free_dmatrix(atten,0,400,0,1191);
     
}
   
   

