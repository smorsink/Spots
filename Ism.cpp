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


class LightCurve Attenuate (class LightCurve* incurve, double* tbnew){

	class LightCurve curve;
	unsigned int numbins;

	curve = (*incurve);
	numbins = curve.numbins;

	for (unsigned int p = 0; p < NCURVES; p++){
	  //std::cout << "p = " << p << "atten = " << tbnew[p] << std::endl;
	  for (unsigned int i = 0; i < numbins; i++){
        
	    curve.f[p][i] *= tbnew[p]; 
        	
	  }
	}
	return curve;
}



void ReadTBNEW(double nh, double *tbnew ){

  double **atten;
  atten = dmatrix(0,400,0,1191);
  double *tbnew_nh, *tbnew_energy;
  tbnew_nh = dvector(0,1191);  // attenuation for each photon energy
  tbnew_energy = dvector(0,1191); // value of photon energy 
  int inh;

  // TBNEW file has entries for 1190 photon energies!

  
  // Read in the attenuation file

      std::cout << "Read in TBNEW!" << std::endl;
      std::ifstream ism;       // input stream
      ism.open("ISM/tbnew_full.txt");
      char line[265]; // line of the data file being read in
      double get_nh, get_e, get_att;
      for (unsigned int k(1); k<= 400; k++){
	for (unsigned int i(0); i<1191; i++){
	  ism.getline(line,265);
	  //std::cout << "line = " << line << std::endl;
	  sscanf( line, "%lf %lf %lf", &get_nh, &get_e, &get_att );
	  tbnew_energy[i] = get_e;
	  atten[k][i] = get_att;
	}
      }
      std::cout << "nh = " << nh << "e18 cm^2" << std::endl;
      inh = nh;
      //std::cout << "inh = " << inh  << std::endl;
     
     
      if (inh < 10 && inh > 0)
	tbnew_nh = atten[1];      
      if (inh >= 10){
	inh = 2+(nh-10)/5;
	tbnew_nh = atten[inh];
      }
      

      for (unsigned int p(0); p<NCURVES; p++){

	if (p%2==0){ //p is even
	  tbnew[p] = 0.25*(3.0*tbnew_nh[p/2] + tbnew_nh[p/2+1]);
	}
	else{ //p is odd
	  tbnew[p] = 0.25 * (tbnew_nh[p/2] + 3.0*tbnew_nh[p/2+1]);
	}	
      }
      
      //std::cout << "tbnew[0] = " << tbnew[0] << std::endl;

	// 20230117 FIX UP Energies

      free_dmatrix(atten,0,400,0,1191);
      //free_dvector(tbnew_nh,0,1191);
}
   
   

