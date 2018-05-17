/***************************************************************************************/
/*                                    Instru.cpp

    This holds the instrument-specifc routines used in Spot.cpp. Contains attenuation tables
    pre-calculated for NICER and instrument response curve. 

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
#include "instru.h"
#include <stdio.h>
using namespace std;


class LightCurve Attenuate (class LightCurve* incurve, unsigned int attenuation, double nh, double* tbnew){

	class LightCurve curve, newcurve;
	unsigned int numbands, numbins;

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



class LightCurve Inst_Res2 (class LightCurve* incurve, unsigned int inst_curve){

  class LightCurve curve, newcurve;
	unsigned int  numbands, numbins;
	double temp;

	//Pour numbands and numbins into newcurve
	curve = (*incurve);
	
	numbands = curve.numbands;
	numbins = curve.numbins;
	newcurve.numbands = numbands;
	newcurve.numbins = numbins;
	
	double elow,ehigh;
	int newindex;

    for (unsigned int p =0; p<NCURVES; p++)
	for (unsigned int i=0; i<numbins;i++)
	 newcurve.f[p][i] = 0.0;


    for (unsigned int p = 0; p < NCURVES; p++){
      //if (p==299) std::cout << "p=" << p << "start[p] = " << start[p] << " curve[p][0]=" << curve.f[p][0] << std::endl;
      for (unsigned int j=0; j<=76; j++){   
	newindex = j + curve.start[p] - 1;   
	//std::cout << "j="<< j << " newindex=" << newindex << std::endl;
	if ( newindex < NCURVES)
	  for (unsigned int i = 0; i < numbins; i++){
	    newcurve.f[newindex][i] += curve.f[p][i] * curve.response[p][j];
	    /*if (newindex==1 && i==0) std::cout << " j=" << j << " newindex=" << newindex 
	    				       << " f["<<p<<"]=" << curve.f[p][i] 
	    				       << " response["<<p<<"]=" << curve.response[p][j]  
	    				       << " flux = " << newcurve.f[newindex][i] << std::endl;*/
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

class LightCurve Read_Background_Guess (class LightCurve* incurve, char *background_file){

	class LightCurve newcurve, curve;
	unsigned int background_size(0), numbands, numbins;
	double temp;
	char line[1024];
	//char cwd[1024], bgddir[1024];
	std::vector<double> background_list;

	FILE *file;

	//Pour numbands and numbins into newcurve
	curve = (*incurve);
	numbands = curve.numbands;
	numbins = curve.numbins;
	newcurve.numbands = numbands;
	newcurve.numbins = numbins;


	file = fopen(background_file,"r");	
	for (unsigned int p = 0; p < numbands-1; p++){

	  fgets(line,20,file);
	  //std::cout << "line = " << line << std::endl;
	  
	  sscanf(line,"%lf;",&temp);
		  
	  //std::cout << "p = " << p << " back = " << temp << std::endl;

	  newcurve.f[p][0] = temp;

	}
	
	fclose(file);
       
	return newcurve;
}



class LightCurve Background_list (class LightCurve* incurve, char *background_file){

	class LightCurve newcurve, curve;
	unsigned int background_size(0), numbands, numbins;
	double temp;
	char cwd[1024], bgddir[1024];
	std::vector<double> background_list;

	//Pour numbands and numbins into newcurve
	curve = (*incurve);
	numbands = curve.numbands;
	numbins = curve.numbins;
	newcurve.numbands = numbands;
	newcurve.numbins = numbins;

    //Read in background numbers
    getcwd(cwd, sizeof(cwd));
    sprintf(bgddir,"%s/Background",cwd);
    chdir(bgddir);

	ifstream file;
	//std::cout << "background_file = " << background_file << std::endl;
	file.open(background_file);

	if(file.is_open()){
        while (file >> temp) {
            background_list.push_back(temp);
            background_size += 1;
        }
	}
	else{
	  throw( Exception( "background list file is not found" ));
	}
	file.close();
	chdir(cwd);


	for (unsigned int p = 0; p < numbands; p++){
        for (unsigned int i = 0; i < numbins; i++){
        	if (p >= background_size){
        		newcurve.f[p][i] = curve.f[p][i];
        	} 
        	else {
        		newcurve.f[p][i] = curve.f[p][i] + background_list[p];
        	}
        }
	}
	return newcurve;
}


class LightCurve AGN_Background (class LightCurve* incurve, double agnbackground, double nh){

	class LightCurve newcurve, curve;
	unsigned int numbands, numbins;
	double temp, gamma, n(0.0), A;
	char cwd[1024], bgddir[1024];
	std::vector<double> nicer_area, ism;

	//Pour numbands, numbins, and nh into newcurve
	curve = (*incurve);
	numbands = curve.numbands;
	numbins = curve.numbins;
	//nh = curve.nh;
	newcurve.numbands = numbands;
	newcurve.numbins = numbins;
	//newcurve.nh = nh;

	//Read in NICER effective area
    getcwd(cwd, sizeof(cwd));
    sprintf(bgddir,"%s/Area",cwd);
    chdir(bgddir);

    ifstream file;
	file.open("NICER_fine_area.txt");

	if(file.is_open()){
        while (file >> temp) {
        	file >> temp;
            nicer_area.push_back(temp);
        }
	}
	else{
	  throw( Exception( "nicer_area_file is not found" ));
	}
	file.close();
	chdir(cwd);

    //Read in ISM factors
    getcwd(cwd, sizeof(cwd));
    sprintf(bgddir,"%s/ISM",cwd);
    chdir(bgddir);

    ifstream file1;
	file1.open("tbnew_4e19.txt");

	if(file1.is_open()){
        while (file1 >> temp) {
        	file1 >> temp;
        	file1 >> temp;
            ism.push_back(temp);
        }
	}
	else{
	  throw( Exception( "ISM file is not found" ));
	}
	file1.close();
	chdir(cwd);

	//Calculate Normalization

	gamma = -2.59;

	for (unsigned int i = 0; i < 301; i++){
		n += nicer_area[i] * pow((i/100+0.1),gamma) * 0.01 * pow(ism[i],nh);
	}
	A = agnbackground/n;


	// Apply Background
	for (unsigned int p = 0; p < numbands; p++){
        for (unsigned int i = 0; i < numbins; i++){
        	newcurve.f[p][i] = curve.f[p][i] + nicer_area[p] * pow((p/100+0.1),gamma) * A * 0.01 * pow(ism[p],nh) / numbins;        	
        }
	}

    return newcurve;
}

/* Phase-independent power law background */

class LightCurve PowerLaw_Background (class LightCurve* incurve, double normalization, double power){

  class LightCurve curve;
  curve = (*incurve);
	 
  double E_diff = (curve.para.E_band_upper_1 - curve.para.E_band_lower_1)/curve.numbands;

  // Apply Background
  for (unsigned int p = 0; p < NCURVES; p++){
    double  E0 = (curve.para.E_band_lower_1+(p+0.5)*E_diff);
    for (unsigned int i = 0; i < curve.numbins; i++){
      curve.f[p][i] = normalization * pow(E0,power);        	
    }
  }

  std::cout << "Background flux[0][0] = " << curve.f[0][0] << std::endl;

  return curve;
}




class LightCurve Sky_Background (class LightCurve* incurve, double skybackground){

	class LightCurve newcurve, curve;
	unsigned int numbands, numbins;
	double temp;
	char cwd[1024], bgddir[1024];
	std::vector<double> skyback;

	//Pour numbands and numbins into newcurve
	curve = (*incurve);
	numbands = curve.numbands;
	numbins = curve.numbins;
	//nh = curve.nh;
	newcurve.numbands = numbands;
	newcurve.numbins = numbins;
	//newcurve.nh = nh;

    //Read in background numbers
    getcwd(cwd, sizeof(cwd));
    sprintf(bgddir,"%s/Background",cwd);
    chdir(bgddir);

	ifstream file;
	//std::cout << "background_file = " << background_file << std::endl;
	file.open("skyback_1191.txt");

	if(file.is_open()){
        while (file >> temp) {
        	file >> temp;
            skyback.push_back(temp);
        }
	}
	else{
	  throw( Exception( "Diffuse sky background file is not found" ));
	}
	file.close();
	chdir(cwd);



	for (unsigned int p = 0; p < numbands; p++){
        for (unsigned int i = 0; i < numbins; i++){
        		newcurve.f[p][i] = curve.f[p][i] + skyback[p]*skybackground/0.3*16/numbins;
				//This list is normalized to 0.3 counts per second with 16 phse bins, calculated using May 2014 NICER effective area	
        }
	}
	return newcurve;
}

