/***************************************************************************************/
/*                                   Response.cpp   */
/* Reads in a light curve and a response matrix and applies response to model */

   
/***************************************************************************************/



// INCLUDE ALL THE THINGS! 
// If you do not get this reference, see 
//       http://hyperboleandahalf.blogspot.ca/2010/06/this-is-why-ill-never-be-adult.html
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <exception>
#include <vector>
#include <string>
#include "Struct.h"
#include "nrutil.h"
#include <unistd.h>
#include <string.h>

// MAIN
int main ( int argc, char** argv ) try {  // argc, number of cmd line args; 
                                          // argv, actual character strings of cmd line args

  /*********************************************/
  /* VARIABLE DECLARATIONS AND INITIALIZATIONS */
  /*********************************************/
    
  std::ofstream out;      // output stream; printing information to the output file

  std::ifstream in;

  
  unsigned long int *start;                    // Starting channels for Instrument Response 
  double *offaxis, *arf;
  double **response;          // Instrument Response Curve

  double *energy;
  
  start = lvector(0,NCURVES);
  offaxis = dvector(0,NCURVES);
  arf = dvector(0,NCURVES);

  response = dmatrix(0,NCURVES,0,400);
  energy = dvector(0,NCURVES);
    
  unsigned int NS_model(1),       // Specifies oblateness (option 3 is spherical)
    spectral_model(0),    // Spectral model choice (initialized to blackbody)
    beaming_model(0),     // Beaming model choice (initialized to isotropic)
    numbins(MAX_NUMBINS), // Number of time or phase bins for one spin period; Also the number of flux data points
    databins(MAX_NUMBINS),   // Number of phase bins in the data
    numphi(1),            // Number of azimuthal (projected) angular bins per spot
    numtheta(1),          // Number of latitudinal angular bins per spot
    numtheta_in(1),
    spotshape(0), 		  // Spot shape; 0=standard
    numbands(NCURVES),    // Number of energy bands;
    attenuation(0),       // Attenuation flag, specific to NICER targets with implemented factors
    inst_curve(0);		  // Instrument response flag, 1 = NICER response curve


  int NlogTeff, Nlogg, NlogE, Nmu, Npts;

  char out_file[256] = "flux.txt",    // Name of file we send the output to; unused here, done in the shell script
       in_file[256] = "No File Name Specified!", 
       T_mesh_file[100],              // Input file name for a temperature mesh, to make a spot of any shape
       data_file[256],                // Name of input file for reading in data
    filenameheader[256]="Run",
       background_file[256];

         
 
    		
  // Create LightCurve data structure
  class LightCurve curve, newcurve;  // variables curve and normalized curve, of type LightCurve
  //class LightCurve curve;
  //class LightCurve *flxcurve;
  //class DataStruct obsdata;           // observational data as read in from a file


  /*********************************************************/
  /* READING IN PARAMETERS FROM THE COMMAND LINE ARGUMENTS */
  /*********************************************************/
    
    for ( int i(1); i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {
	            
	    case 'i': // Name of input file
	            	sscanf(argv[i+1], "%s", in_file);
	            	break;
	                	          
	    case 'n':  // Number of phase or time bins
	                sscanf(argv[i+1], "%u", &databins);
					if ( databins < MIN_NUMBINS) {
			  			numbins = MIN_NUMBINS;
					}
					else
			  			numbins = databins;		 
	                break;
	               
	    case 'o':  // Name of output file
	                sscanf(argv[i+1], "%s", out_file);
	                break;
	          

	    case 'R':  // Instrument Response Curve
	    			sscanf(argv[i+1], "%u", &inst_curve);
	    			break;
	      	          	               
            } // end switch	
        } // end if
    } // end for


 

    /*************************/
    /* OPENING THE INPUT FILE */
    /*************************/
   
      std::cout << "Opening input file" << std::endl;	
      in.open( in_file );  // opening the file with observational data


      /****************************************/
      /* READING IN FLUXES FROM THE DATA FILE */
      /****************************************/
    

    if (in.is_open()){
      std::cout << "reading data" << std::endl;
      std::cout << "number of data bins " << databins << " number of bands " << numbands << std::endl;
    	double temp;
	for (unsigned int k = 0; k < numbands; k++){
	  for (unsigned int j = 0; j < databins; j++){
	    in >> temp;
	    if (j==0) energy[k] = temp;	    
	    in >> temp;
	    curve.t[j] = temp;

	    in >> temp;
	    curve.f[k][j] = temp;
	  }
    	}
    	in.close();
    }
    // Finished reading in the data file
	

    /******************************************/
    /*      OPEN INSTRUMENT RESPONSE CURVE    */
    /******************************************/


     //Version 1.02 (March 2018) of response matrix

      	std::ifstream file;
	std::cout << "Opening v1.02 NICER Response Matrix" << std::endl;
	file.open("Area/rmf_v1.02_newformat.txt");
	double elow,ehigh,junk;
	int channel;
	// start[p] holds the number of non-zero elements
	
	if(file.is_open()) {
	  for (unsigned int p(0);p<NCURVES;p++){
	    file >> channel;
	    file >> elow;
	    file >> ehigh;
	    file >> start[p];
	    if (p==0){
	      std::cout << "Channel #" << channel
			<< " Elo = " << elow
			<< " Ehi = " << ehigh
			<< " Number of Nonzeros = " << start[p]
			<< std::endl;	      
	    }
	    
	    for (unsigned int j(0); j<start[p]; j++){
	      file >> response[p][j];
	      if (p==0)
		std::cout << response[p][j] << " ";
	    }
	    if (p==0)
	      std::cout << std::endl;
	  }	  
	}
        file.close();

	std::cout << "Off-Axis Correction Area v1.02" << std::endl;
	file.open("Area/offset_correction_chans25to355_noheader.txt");
	if(file.is_open()) {
	  for (unsigned int p(0);p<24;p++){
	    offaxis[p] = 0.0;
	  }	  
	  for (unsigned int p(24);p<NCURVES;p++){
	    file >> channel;
	    file >> junk;
	    file >> junk;
	    file >> offaxis[p];	    
	  }	  
	}	
        file.close();

	std::cout << "On-Axis Area v1.02" << std::endl;
	file.open("Area/ni_xrcall_onaxis_v1.02_arf_noheader.txt");
	if(file.is_open()) {	  
	  for (unsigned int p(0);p<NCURVES;p++){
	    file >> channel;
	    file >> junk;
	    file >> junk;
	    file >> arf[p];
	    file >> junk;
	    file >> junk;
	    file >> junk;
	    file >> junk;
	    file >> junk; 	    
	  }	  
	}	
        file.close();

	
	
	// Apply the Area Response (ARF)

	for (unsigned int p =0; p<NCURVES; p++)
	  for (unsigned int i=0; i<numbins;i++)
	    curve.f[p][i] *= arf[p];

       
	// Apply the Response Matrix
	
	for (unsigned int p =0; p<NCURVES; p++)
	  for (unsigned int i=0; i<numbins;i++)
	    newcurve.f[p][i] = 0.0;


	for (unsigned int p = 0; p < NCURVES; p++){
	  for (unsigned int j=0; j<=start[p]; j++){   
	    for (unsigned int i = 0; i < numbins; i++){
	      newcurve.f[j][i] += curve.f[p][i] * response[p][j];
	    }
	 
	  }

	}




 
 
    /***************************************************/
    /* WRITING COLUMN HEADINGS AND DATA TO OUTPUT FILE */
    /***************************************************/
  
 
    //std::ofstream out;
      out.open(out_file, std::ios_base::trunc);
      out.precision(10);
      if ( out.bad() || out.fail() ) {
	std::cerr << "Couldn't open output file: " << out_file << " Exiting." << std::endl;
	return -1;
      }
      else
	std::cout << "Opening "<< out_file << " for printing " << std::endl;
 

      for ( unsigned int p(25); p < 300; p++ ) {
	for ( unsigned int i(0); i < numbins; i++ ) {
	  if (numbands !=1)
	    
	  
	  out << p << "\t";
	  out << curve.t[i]<< "\t";		
	  //	  out << curve.f[p][i] << "\t";
	  //out << newcurve.f[p][i] * 1e6 << "\t";
	  out << newcurve.f[p][i] * offaxis[p] * 1e6/32.0 << "\t";
	  out    << std::endl;
	}
      }


      

    
      out.close();
    




    return 0;
} 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}

