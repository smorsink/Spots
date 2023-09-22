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


class LightCurve Wabs (class LightCurve* incurve, unsigned int attenuation, double nh){
  // OLD - Not used
	class LightCurve curve, newcurve;
	double temp;
	char cwd[1024];
	std::vector<double> atten_ener, atten_width, atten_factors;
	unsigned int atten_size(0), numbands, numbins;

	//Pour numbands and numbins into newcurve

	curve = (*incurve);
	numbands = curve.numbands;
	numbins = curve.numbins;
	newcurve.numbands = numbands;
	newcurve.numbins = numbins;

    //Read in ism attenuation factors
    //getcwd(cwd, sizeof(cwd));
    //sprintf(ismdir,"%s/ISM",cwd);
    //chdir(ismdir);

	//cout << "Reading Attenuation File" << std::endl;
	ifstream file;
	switch (attenuation) {
		case 1: // j0030 WABS
			file.open("ISM/j0030_wabs_1.8e20.txt");
			break;
		case 2: // j0030 TBABS
			file.open("ISM/j0030_tbabs_1.8e20.txt");
			break;
		case 3: // j0437 WABS
			file.open("ISM/j0437_wabs_4e19.txt");
			break;
		case 4: // j0437 TBABS
			file.open("ISM/j0437.tbabs_4e19.txt");
			break;
		case 5: //
			file.open("ISM/tbnew_4e19.txt");
			break;
	}

	if(file.is_open()){
        while (file >> temp) {
            atten_ener.push_back(temp);
            file >> temp;
            atten_width.push_back(temp);
            file >> temp;
            atten_factors.push_back(temp);
            atten_size += 1;
            //cout << temp << endl;
        }
    }
    else{
        cout << "attenuation file is not found" << endl;
    }
    file.close();
	chdir(cwd);

	for (unsigned int p = 0; p < numbands; p++){
        for (unsigned int i = 0; i < numbins; i++){
        
	  newcurve.f[p][i] = curve.f[p][i] * pow(atten_factors[p],nh);
        	

        	/*
        	if (newcurve.f[p][i] > incurve.f[p][i]){
				cout << "attenuation factor greater than 1!" << endl;
        	}
        	*/
        }
	}
	return newcurve;
}


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
	  // std::cout << "Attenuate: p = " << p << " computed energy = " << curve.elo[p] << std::endl;

	  factor = curve.elo[p]*1e2 - 10.0 + 0.5; //The 0.5 takes care of the rounding from double to int
	  index = factor;

	  /* std::cout << "factor = " << factor
	    << " index = " << index << " NH energy = " << tbnew->energy[index]
		    << " atten = " << tbnew->attenuation[index]
		    << std::endl;*/
	  
	  
	  
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

  tbnew->energy = dvector(0,1191);
  tbnew->attenuation = dvector(0,1191);

  // TBNEW file has entries for 1190 photon energies!

  
  // Read in the attenuation file

      std::cout << "Read in TBNEW!" << std::endl;
      std::ifstream ism;       // input stream
      ism.open("ISM/tbnew_full.txt");
      char line[265]; // line of the data file being read in
      double get_nh, get_e, get_att;
      for (unsigned int k(1); k< 400; k++){ // k is a label for the value of nH
	for (unsigned int i(0); i<1191; i++){ // i is a label for the photon energy
	  ism.getline(line,265);
	  //std::cout << "line = " << line << std::endl;
	  sscanf( line, "%lf %lf %lf", &get_nh, &get_e, &get_att );
	  tbnew->energy[i] = get_e;
	  atten[k][i] = get_att;
	  /*if (i==0 && k < 300)
	    std::cout
	      << " k = " << k 
	      << " nh = " << get_nh
	      << " energy = " << get_e << " atten = " << get_att <<std::endl;*/
	  
	}
      }
      std::cout << "nh = " << nh << "e18 cm^2" << std::endl;
      inh = nh;
      std::cout << "inh = " << inh  << std::endl;
     
      // Read in appropriate column of matrix into vector (depending on value of NH)
      if (inh < 10 && inh > 0)
	tbnew->attenuation = atten[1];      
      if (inh >= 10){
	inh = 2+(nh-10)/5;
	std::cout << "inh = " << inh  << std::endl;
	tbnew->attenuation = atten[inh];
      }
      std::cout << "energy 0 = " << tbnew->energy[0] << " attenuation = " << tbnew->attenuation[0] << std::endl;
      std::cout << "energy 1 = " << tbnew->energy[1] << " attenuation = " << tbnew->attenuation[1] << std::endl;

     

	// 20230117 FIX UP Energies

      free_dmatrix(atten,0,400,0,1191);
     
}
   
   

