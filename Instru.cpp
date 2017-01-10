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

double Attenuate (double flux_before, unsigned int attenuation){

	double flux_after(0), temp;
	char cwd[1024], ismdir[1024];
	std::vector<double> atten_ener, atten_width, atten_factors;

    //Read in hydrogen atmosphere parameters
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
        }
    }
    else{
        cout << "attenuation file is not found" << endl;
    }
    file.close();
	chdir(cwd);

	/*

	some search and interpolate routine to find the right attenuation factor (fact)
	flux_after = flux_before * fact;

	*/

	return flux_after;
}

double Inst_Res (double flux_before, unsigned int inst_curve){

	double flux_after(0);
	return flux_after;

}
