/***************************************************************************************/
/*                                   MeasuredSpot.cpp

    This code reads in a pulse profile in time/energy space, computed with no ISM absorption.
    ISM absorption is added, and the instrument response file is applied.
    
    Based on code written by Coire Cadeau and modified by Sharon Morsink and 
    Abigail Stevens.
    
    PGxx refers to equation xx in Poutanen and Gierlinski 2003, arxiv: 0303084v1
    MLCBxx refers to equation xx in Morsink, Leahy, Cadeau & Braga 2007, arxiv: 0703123v2
    
    (C) Coire Cadeau, 2007; Source (C) Coire Cadeau 2007, all rights reserved.
*/
/***************************************************************************************/

// Changes

// 2016-08-05 - SMM: Added a function to Chi.cpp called Bend. Bend computes the look-up table for the bending angles. Removed the section of code in Spot.cpp where this is computed. 
// 2016-08-05 - SMM: Bend is now called one time each latitude
// 2016-09-19 - SMM: This version should do large spots correctly! However 2nd spot hasn't been changed.
// 2016-11-13 - SMM: When spot goes over the pole, changed the order that pieces are computed.
// 2016-12-01 - SMM: Adding dynamic memory allocation
// 2023-01-24 - SMM: Attenuation due to the ISM (TBNEW) is applied to the signal before the instrumental response is added.

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
#include "OblDeflectionTOA.h"
//#include "Chi.h"
//#include "Atmo.h"
//#include "Hydrogen.h"
//#include "TimeDelays.h"
#include "Instru.h"
#include "Ism.h"
//#include "PolyOblModelNHQS.h"
//#include "PolyOblModelCFLQS.h"
//#include "SphericalOblModel.h"
//#include "OblModelBase.h"
#include "Units.h"
#include "Exception.h"
#include "Struct.h"
#include "time.h"
#include "interp.h"
#include "nrutil.h"
#include <unistd.h>
#include <string.h>

// MAIN
int main ( int argc, char** argv ) try {  // argc, number of cmd line args; 
                                          // argv, actual character strings of cmd line args


  std::cout << "Hello!" << std::endl;
  
  /*********************************************/
  /* VARIABLE DECLARATIONS AND INITIALIZATIONS */
  /*********************************************/
    
  std::ofstream out;      // output stream; printing information to the output file

  double incl_1(90.0),          // Inclination angle of the observer, in degrees
    incl_2(90.0),               // PI - incl_1; needed for computing flux from second hot spot, since cannot have a theta greater than 
    theta_1(90.0),              // Emission angle (latitude) of the first upper spot, in degrees, down from spin pole
    theta_2(90.0),              // Emission angle (latitude) of the second lower spot, in degrees, up from spin pole (180 across from first spot)
    d_theta_2(0.0),
    d_incl_2(0.0),
    mass,                       // Mass of the star, in M_sun
    rspot(0.0),                 // Radius of the star at the spot, in km
    mass_over_req,              // Dimensionless mass divided by radius ratio
    omega,                      // Frequency of the spin of the star, in Hz
    req,                        // Radius of the star at the equator, in km
    phaseshift(0.0),                    // Phase shift of spot (in radians)
    spot_temperature(0.0),      // Temperature of first spot, in the star's frame, in Kelvin
    spot2_temperature(0.0),		// Temperature of second spot
    rho(0.0),                   // Angular radius of the first spot, in radians
    rho2(0.0),
    dphi(1.0),                  // Each chunk of azimuthal angle projected onto equator, when broken up into the bins (see numphi)
    phishift,
    mu_1(1.0),                  // = cos(theta_1), unitless
    mu_2(1.0),                  // = cos(theta_2), unitless
    cosgamma,                   // Cos of the angle between the radial vector and the vector normal to the surface; defined in equation 13, MLCB
    //Flux[NCURVES][MAX_NUMBINS], // Array of fluxes; Each curve gets its own vector of fluxes based on photons in bins.
    //Temp[NCURVES][MAX_NUMBINS],
    E_band_lower_1(2.0),        // Lower bound of first energy band to calculate flux over, in keV.
    E_band_upper_1(3.0),        // Upper bound of first energy band to calculate flux over, in keV.
    chisquared(1.0),            // The chi^2 of the data; only used if a data file of fluxes is inputed
    distance(3.0857e20),        // Distance from earth to the NS, in meters; default is 10kpc
    obstime(1.0),               // Length of observation (in seconds)
    phase_2(0.5),				// Phase of second spot, 0 < phase_2 < 1
    nh(0.0);					// real nh = nh*e18

 
  
  double SurfaceArea(0.0);
  double Eobs;

  double E0, DeltaE(0.0), DeltaLogE(0.0), rot_par;
    
  unsigned int NS_model(1),       // Specifies oblateness (option 3 is spherical)
    spectral_model(0),    // Spectral model choice (initialized to blackbody)
    beaming_model(0),     // Beaming model choice (initialized to isotropic)
    numbins(MAX_NUMBINS), // Number of time or phase bins for one spin period; Also the number of flux data points
    databins(MAX_NUMBINS),   // Number of phase bins in the data
    numphi(1),            // Number of azimuthal (projected) angular bins per spot
    numtheta(1),          // Number of latitudinal angular bins per spot
    numtheta_in(1),
    spotshape(0), 		  // Spot shape; 0=standard
    numbands(NCURVES),    // Number of energy bands that will be computed;
    inst_curve(0);		  // Instrument response flag, 0 = No instrument response; 1 = NICER response curve

  // Note: space will be allocated for a total of NCURVES different energy bands
  // We will compute only numbands different energy bands

  
  int NlogTeff, Nlogg, NlogE, Nmu, Npts;

  char out_file[256] = "flux.txt",    // Name of file we send the output to; unused here, done in the shell script
    test_file[256] = "test.txt",
       bend_file[256] = "No File Name Specified!", 
       data_file[256],                // Name of input file for reading in data
    filenameheader[256]="Run";
     
         
  // flags!
  bool incl_is_set(false),         // True if inclination is set at the command line (inclination is a necessary variable)
    	 theta_is_set(false),        // True if theta is set at the command line (theta is a necessary variable)
    	 mass_is_set(false),         // True if mass is set at the command line (mass is a necessary variable)
    	 rspot_is_set(false),        // True if rspot is set at the command line (rspot is a necessary variable)
    	 omega_is_set(false),        // True if omega is set at the command line (omega is a necessary variable)
    	 model_is_set(false),        // True if NS model is set at the command line (NS model is a necessary variable)
    	 datafile_is_set(false),     // True if a data file for inputting is set at the command line
    kelvin(false),                    // True if Temperature is in Kelvin; Otherwise in keV
    logEflag(false),
    	 ignore_time_delays(false),  // True if we are ignoring time delays
         bend_file_is_set(false),
    two_spots(false);           // True if we are modelling a NS with two antipodal hot spots
    	
   
    		
  // Create LightCurve data structure
  class LightCurve curve;  // variables curve and normalized curve, of type LightCurve
  // class LightCurve *flxcurve, *flxcurve2;
  //class DataStruct obsdata;           // observational data as read in from a file


  
  /*********************************************************/
  /* READING IN PARAMETERS FROM THE COMMAND LINE ARGUMENTS */
  /*********************************************************/
    
    for ( int i(1); i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {

        case 'A': // ISM column density, in multiples of 1e18cm^2
        			sscanf(argv[i+1], "%lf", &nh);
        			break;
	            
	    
	                
	 
	            
	    case 'I': // Name of input file
	            	sscanf(argv[i+1], "%s", data_file);
	            	datafile_is_set = true;
	            	break;
	                
	          	    

	 
	          
	    case 'n':  // Number of phase or time bins
	                sscanf(argv[i+1], "%u", &databins);
			//	if ( databins < MIN_NUMBINS) {
			//			numbins = MIN_NUMBINS;
			//		}
			//		else
			numbins = databins;		 
	                break;

	   
			
	    case 'o':  // Name of output file
	                sscanf(argv[i+1], "%s", out_file);
	                break;
	          
	   

	    case 'R':  // Instrument Response Curve
	    			sscanf(argv[i+1], "%u", &inst_curve);
	    			break;
	      	          
	   
			
	    case 'S': // Number of Energy bands
	      			sscanf(argv[i+1], "%u", &numbands);
	      			break;

	   
	            
	    case 'u': // Lower limit of first energy band, in keV
	            	sscanf(argv[i+1], "%lf", &E_band_lower_1);
	            	//E_band_lower_1_set = true;
	            	break;
	            
	    case 'U': // Upper limit of first energy band, in keV
	            	sscanf(argv[i+1], "%lf", &E_band_upper_1);
	            	//E_band_upper_1_set = true;
	            	break;
	                
	            
	    case 'x': // Part of funny NICER line
	            	sscanf(argv[i+1], "%lf", &E0);
	            	break;
	            
	    case 'X': // Part of funny NICER line
	            	sscanf(argv[i+1], "%lf", &DeltaE);
	            	break;
	               	
	    case 'Z': // Observation time (in seconds)
	      sscanf(argv[i+1], "%lf", &obstime);
	      break;
	            
	            
                case 'h': default: // Prints help
		  std::cout << "\n\nMeasuredSpot help:  -flag description [default value]\n" << std::endl      	            		  
			    << "-A ISM column density, in multiples of base value [0]" << std::endl
			    << "      base value is 1e18 cm^-2 " << std::endl
			    << "-b Bending Angle File" << std::endl
			    << "-B Phase of second spot, 0 < phase_2 < 1 [0.5]" << std::endl
			    << "-C Temperature of second spot, in log(K) [0.0]" << std::endl
			    << "-d Size of second spot. [0.0]" << std::endl
			    << "-D Distance from earth to star, in kpc. [10]" << std::endl
			    << "-e * Latitudinal location of emission region, in degrees, between 0 and 90." << std::endl
			    << "-E Offset latitudinal angle of second spot (from antipodal), in degrees [0.0]" << std::endl
			    << "-f * Spin frequency of star, in Hz." << std::endl
			    << "-g Atmosphere beaming model [0]:" << std::endl
			    << "      0 for BB, no beaming" << std::endl
			    << "      11 for Wynn's full NSXH table" << std::endl
			    << "-i * Inclination of observer, in degrees, between 0 and 90." << std::endl
			    << "-I Input filename." << std::endl
		                     
			    << "-K Kelvin or keV" << std::endl
			    << "-l Time shift (or phase shift), in seconds." << std::endl
			    << "-L Offset inclination of observer, in degrees, between 0 and 90." << std::endl
			    << "-m * Mass of star in Msun." << std::endl          
			    << "-n Number of phase or time bins. [128]" << std::endl
      	  	               
			    << "-o Output filename." << std::endl
		                      
			    << "-p Angular radius of spot, rho, in radians. [0.0]" << std::endl
			    << "-P Hot spot shape model [0]:" << std::endl
			    << "      0 for standard circular" << std::endl
			    << "      1 for circular in rotating frame" << std::endl
			    << "-q * Oblateness Model of star: [1]" << std::endl
			    << "      1 for Neutron/Hybrid quark star poly model" << std::endl
			    << "      2 for CFL quark star poly model" << std::endl
			    << "      3 for spherical model" << std::endl
			    << "-R Instrument response curve [0]:" << std::endl
			    << "      0 No instrument response" << std::endl
			    << "      1 NICER combined ARF & RMF response matrix" << std::endl
			    << "-r * Radius of star (at the equator), in km." << std::endl
			    << "-s Spectral model of radiation: [0]" << std::endl
			    << "      0 for monochromatic." << std::endl
			    << "      1 for NICER fake lines" << std::endl
			    << "      2 for BB integrated" << std::endl
			    << "      3 for Atmosphere integrated" << std::endl
			    << "-S Number of energy bands: [NCURVES]" << std::endl
			    << "-t Number of theta bins for large spots. Must be < 30. [1]" << std::endl
			    << "-T Temperature of the spot, in keV or Kelvin. [0]" << std::endl
			    << "-u Energy bands' lower limit, in keV. [2]" << std::endl
			    << "-U Energy bands' upper limit, in keV. [3]" << std::endl
			    << "-x NICER funny line E0" << std::endl
			    << "-X NICER funny line DeltaE" << std::endl
			    << "-Z Observation time in seconds [1.0]" << std::endl
			    << "-2 Flag for calculating two hot spots. Using this sets it to true. [false]" << std::endl
			    << " Note: '*' next to description means required input parameter." << std::endl
			    << std::endl;
		  // return 0;
            } // end switch	
        } // end if
    } // end for


    
    
 
 

    
    // Widths of energy intervals
    if (logEflag) { // logarithmic intervals
      DeltaLogE = (log10(E_band_upper_1) - log10(E_band_lower_1))/(numbands);
      //std::cout << "DeltaLogE = " << DeltaLogE << std::endl;
    }
    else{ // linear intervals
      DeltaE = (E_band_upper_1 - E_band_lower_1)/(numbands);
      //std::cout << "DeltaE = " << DeltaE << std::endl;
    }



    
 
    /**********************************/
    /* PASS VALUES INTO THE STRUCTURE */
    /**********************************/    	
    
    curve.para.mass = mass;
    curve.para.mass_over_r = mass_over_req;
    curve.para.omega = omega;
    curve.para.omega_bar_sq = rot_par;
    curve.para.radius = req;
    curve.para.req = req;
    curve.para.theta = theta_1;
    curve.para.theta_c = theta_1;
    curve.para.incl = incl_1;
    curve.para.temperature = spot_temperature;
    //curve.para.ts = ts;
    curve.para.E_band_lower_1 = E_band_lower_1;
    curve.para.E_band_upper_1 = E_band_upper_1;

    curve.para.DeltaE = DeltaE;
    curve.para.DeltaLogE = DeltaLogE;
    curve.para.E0 = E0;
    curve.para.distance = distance;
    curve.numbins = numbins;
    curve.para.phaseshift=phaseshift;

    curve.flags.ignore_time_delays = ignore_time_delays;
    curve.flags.spectral_model = spectral_model;
    curve.flags.beaming_model = beaming_model;
    curve.flags.ns_model = NS_model;
    curve.flags.bend_file = bend_file_is_set;
    curve.flags.inst_curve = inst_curve;
    curve.numbands = numbands; // The number of energy bands computed
    curve.tbands = NCURVES; // The number of energy bands allocated into memory
    curve.cbands = numbands; // The number of energy bands computed
    curve.flags.kelvin = kelvin;
    curve.flags.logEflag = logEflag;
    curve.flags.spotshape = spotshape;


    /***************************/
    /* Read in TBNew ISM File  */
    /***************************/ 

    class ISM ism;   
    // Read in TBNew using the correct value of NH.
    if (nh != 0.0)
      ReadTBNEW(nh,&ism);
    

   

  
    /******************************************/
    /*      OPEN INSTRUMENT RESPONSE CURVE    */
    /******************************************/
    class Instrument nicer;
    std::cout << "Read in the Instrumental Response: ints_curve = " << inst_curve << std::endl;

    if (inst_curve == 1 ){ // READ in the response matrix
 
      ReadResponse( &nicer);
      
    }
    

    /**************************************/
    /* Initialize time, flux and energies */
    /**************************************/

    
    //flxcurve = &normcurve;


    for ( unsigned int i(0); i < numbins; i++ ) {
      curve.t[i] = i / (1.0 * numbins); // + phaseshift/(2.0*Units::PI);  // defining the time used in the lightcurves
        for ( unsigned int p(0); p < curve.numbands; p++ ) {
            curve.f[p][i] = 0.0;
	}
    } 

    // Define photon energies for computation

    for (unsigned int p(0); p <= numbands; p++){

      curve.elo[p] = E_band_lower_1 + p*DeltaE;
      curve.ehi[p] = E_band_lower_1 + (p+1)*DeltaE;

      if (p == 0)
	std::cout << "elo[0] = " << curve.elo[p] << " ehi[0] = " << curve.ehi[p] << std::endl;

    }


    // Read in the time/energy waveform

      std::ifstream in;       // input stream
      in.open(data_file);
      char line[265]; // line of the data file being read in
      double get_e, get_t, get_f, junk1, junk2;

      for (unsigned p(0);p<numbands;p++){ // loop through the energy bands
	for (unsigned j(0);j<numbins;j++){ // loop through the time bins

	  in.getline(line,265);  
	  sscanf( line, "%lf %lf %lf %lf %lf", &get_e, &get_t, &get_f, &junk1, &junk2 );
	  curve.f[p][j] = get_f;
	  //std::cout << line << std::endl;

	}
      }

      std::cout << "Input to MeasuredSpot: F[0][0] = " << curve.f[0][0] << std::endl;
      

   


    /**********************************/
    /*       APPLYING ATTENUATION     */
    /**********************************/

    if (nh != 0){
      std::cout << "Apply ISM to the Signal!" << std::endl;
      std::cout << " Before ISM flux[0][0] = " << curve.f[0][0] << std::endl;
      
      curve = Attenuate(&curve,&ism);

      std::cout << " After ISM flux[0][0] = " << curve.f[0][0] << std::endl;

    }

 

    
    /***************************************/
    /* START OF INSTRUMENT EFFECT ROUTINES */
    /***************************************/

    // Interpolate to create all the other energy bands
    std::cout << "Number of energy bands computed = " << numbands << std::endl;

    curve = ConvertEnergyChannels(&curve, &nicer);
    numbands = curve.numbands;

        std::cout << "Number of NICER Energy Bands = " << numbands << std::endl;


    /***************************************************/
    /* WRITING COLUMN HEADINGS AND DATA TO OUTPUT FILE */
    /***************************************************/
    
    sprintf(test_file,"Test/detect3.txt");
    out.open(test_file, std::ios_base::trunc);
    out.precision(10);
    if ( out.bad() || out.fail() ) {
      std::cerr << "Couldn't open output file: " << out_file << " Exiting." << std::endl;
      return -1;
    }
    else
      std::cout << "Opening "<< out_file << " for printing " << std::endl;

    //out << "#Spot before ISM absorption" << std::endl;
    
    for ( unsigned int p(0); p < numbands; p++ ) {
      for ( unsigned int i(0); i < numbins; i++ ) {

	Eobs = curve.elo[p];
	  
	out << Eobs << "\t";
	out << curve.t[i]<< "\t";		
	out << curve.f[p][i] << "\t";
	out << p << std::endl;
      }
    }
    out.close();
    


   

    
    /******************************************/
    /*     MULTIPLYING OBSERVATION TIME       */
    /******************************************/
    for ( unsigned int i(0); i< databins; i++){
      for (unsigned int p(0);p<numbands-1;p++){
	curve.f[p][i] *= obstime;
      }
    }



    


    
    /******************************************/
    /*  APPLYING INSTRUMENT RESPONSE CURVE    */
    /******************************************/

    //    std::cout << "Apply Instrument Response to Spot: ints_curve = " << curve.flags.inst_curve << std::endl;
    if (curve.flags.inst_curve > 0){
      std::cout << "Applying Instrument Response" << std::endl;
      curve = ApplyResponse(&curve,&nicer);
    } // Finished Applying the Response Matrix

    // double spotcounts = 0.0;
    //double bkgcounts = 0.0;

   


 
    
    /***************************************************/
    /* WRITING COLUMN HEADINGS AND DATA TO OUTPUT FILE */
    /***************************************************/

    //sprintf(out_file,out_file);
      out.open(out_file, std::ios_base::trunc);
      out.precision(10);
      if ( out.bad() || out.fail() ) {
	std::cerr << "Couldn't open output file: " << out_file << " Exiting." << std::endl;
	return -1;
      }
      else
	std::cout << "Opening "<< out_file << " for printing " << std::endl;

      //out << "#Response matrix applied to light curve" << std::endl;
      
      for ( unsigned int p(30); p < curve.numbands; p++ ) {
	for ( unsigned int i(0); i < curve.numbins; i++ ) {

	  Eobs = curve.elo[p];
	  
	  // out << Eobs << "\t";
	  out << p << "\t";
	  out << i << "\t";
	  // out << curve.t[i]<< "\t";		
	  out << curve.f[p][i] << std::endl;
	  //out << p << std::endl;
	}
      }
      out.close();
    
      //delete model;
    return 0;
} 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}

