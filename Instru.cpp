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
#include <stdio.h>
using namespace std;

double Attenuate (unsigned int p, double flux_before, unsigned int attenuation){

	double flux_after(0), temp;
	char cwd[1024], ismdir[1024];
	std::vector<double> atten_ener, atten_width, atten_factors;
	unsigned int atten_size(0);

    //Read in ism attenuation factors
    getcwd(cwd, sizeof(cwd));
    sprintf(ismdir,"%s/ISM",cwd);
    chdir(ismdir);

	//cout << "Reading Attenuation File" << std::endl;
	ifstream file;
	switch (attenuation) {
		case 1: // j0030 WABS
			file.open("j0030_wabs_1.8e20.txt");
			break;
		case 2: // j0030 TBABS
			file.open("j0030_tbabs_1.8e20.txt");
			break;
		case 3: // j0437 WABS
			file.open("j0437.wabs_4e19.txt");
			break;
		case 4: // j0437 TBABS
			file.open("j0437.tbabs_4e19.txt");
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

	if (p >= atten_size){
		flux_after = flux_before;
	}else{
		flux_after = flux_before * atten_factors[p];
	}

	if (flux_before < flux_after){
		cout << "attenuation factor greater than 1!" << endl;
	}

	return flux_after;
}

double Inst_Res (unsigned int p, double flux_before, unsigned int inst_curve){

	unsigned int response_size(0);
	double flux_after(0),temp;
	char cwd[1024], resdir[1024];
	std::vector<double> response_list;

	getcwd(cwd, sizeof(cwd));
	sprintf(resdir,"%s/Area",cwd);
	chdir(resdir);

	ifstream file;
	file.open("niceravgarea.txt");

	if(file.is_open()) {
		while (file >> temp) {
			file >> temp;
			response_list.push_back(temp);
			response_size += 1;
		}
	}else{
        throw( Exception( "instrument response curve file is not found" ));
    }
    file.close();
	chdir(cwd);

	if (p >= response_size){

		 throw( Exception( "instrument response curve doesn't have this band" ));
		 return -1;

	}else{

		flux_after = flux_before * response_list[p];
		return flux_after;
	}
	
}


double Background_list (unsigned int p, double flux_before, char *background_file){
	
	unsigned int background_size(0);
	double flux_after(0),temp;
	char cwd[1024], bgddir[1024];
	std::vector<double> background_list;


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

	if (p >= background_size){

		flux_after = flux_before;

	}else{

	  //std::cout << "p = " << p << " back[p] = " << background_list[p] << " flux = " << flux_before << "+"
	  //	    << background_list[p] << std::endl;

		flux_after = flux_before + background_list[p];
	}

	return flux_after;
}
