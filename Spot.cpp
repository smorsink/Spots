/***************************************************************************************/
/*                                   Spot.cpp

    This code produces a pulse profile once a set of parameter describing the star, 
    spectrum, and hot spot have been inputed.
    
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
#include "Chi.h"
#include "Atmo.h"
#include "TimeDelays.h"
#include "Instru.h"    
#include "PolyOblModelNHQS.h"
#include "PolyOblModelCFLQS.h"
#include "SphericalOblModel.h"
#include "OblModelBase.h"
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
  std::ofstream allflux;      // Prints information about the star's surface area
  std::ofstream testout;  // testing output stream;
  std::ofstream param_out;// piping out the parameters and chisquared

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
    ts(0.0),                    // Phase shift or time off-set from data; Used in chi^2 calculation
    spot_temperature(0.0),      // Temperature of first spot, in the star's frame, in keV
    spot2_temperature(0.0),		// Temperature of second spot
    rho(0.0),                   // Angular radius of the first spot, in radians
    rho2(0.0),
    dphi(1.0),                  // Each chunk of azimuthal angle projected onto equator, when broken up into the bins (see numphi)
    phishift,
    aniso(0.586),               // Anisotropy parameter
    Gamma1(2.0),                // Spectral index
    Gamma2(2.0),                // Spectral index
    Gamma3(2.0),                // Spectral index
    mu_1(1.0),                  // = cos(theta_1), unitless
    mu_2(1.0),                  // = cos(theta_2), unitless
    cosgamma,                   // Cos of the angle between the radial vector and the vector normal to the surface; defined in equation 13, MLCB
    //Flux[NCURVES][MAX_NUMBINS], // Array of fluxes; Each curve gets its own vector of fluxes based on photons in bins.
    Temp[NCURVES][MAX_NUMBINS],
    E_band_lower_1(2.0),        // Lower bound of first energy band to calculate flux over, in keV.
    E_band_upper_1(3.0),        // Upper bound of first energy band to calculate flux over, in keV.
    E_band_lower_2(5.0),        // Lower bound of second energy band to calculate flux over, in keV.
    E_band_upper_2(6.0),        // Upper bound of second energy band to calculate flux over, in keV.
    background(0.0),	        // One background value for all bands.
    dsbackground(0.0),			// Diffuse sky background normalization
    agnbackground(0.0),			// AGN background normalization
    chisquared(1.0),            // The chi^2 of the data; only used if a data file of fluxes is inputed
    distance(3.0857e20),        // Distance from earth to the NS, in meters; default is 10kpc
    obstime(1.0),               // Length of observation (in seconds)
    phase_2(0.5),				// Phase of second spot, 0 < phase_2 < 1
    powerlaw(0.0),
    nh(1.0);					// real nh = nh*4e19

  unsigned long int *start;                    // Starting channels for Instrument Response 
  double *arf, *offaxis;
  double **response;          // Instrument Response Curve

  start = lvector(0,NCURVES);
  arf = dvector(0,NCURVES);
  offaxis = dvector(0,NCURVES);
  response = dmatrix(0,NCURVES,0,400);
  
  double SurfaceArea(0.0);

  double E0, L1, L2, DeltaE, rot_par;
    
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
       bend_file[256] = "No File Name Specified!", 
       T_mesh_file[100],              // Input file name for a temperature mesh, to make a spot of any shape
       data_file[256],                // Name of input file for reading in data
    filenameheader[256]="Run",
       background_file[256];

         
  // flags!
  bool incl_is_set(false),         // True if inclination is set at the command line (inclination is a necessary variable)
    	 theta_is_set(false),        // True if theta is set at the command line (theta is a necessary variable)
    	 mass_is_set(false),         // True if mass is set at the command line (mass is a necessary variable)
    	 rspot_is_set(false),        // True if rspot is set at the command line (rspot is a necessary variable)
    	 omega_is_set(false),        // True if omega is set at the command line (omega is a necessary variable)
    	 model_is_set(false),        // True if NS model is set at the command line (NS model is a necessary variable)
    	 datafile_is_set(false),     // True if a data file for inputting is set at the command line
    	 ignore_time_delays(false),  // True if we are ignoring time delays
         bend_file_is_set(false),
    	 T_mesh_in(false),           // True if we are varying spot shape by inputting a temperature mesh for the hot spot's temperature
   	 two_spots(false),           // True if we are modelling a NS with two antipodal hot spots
    	 //only_second_spot(false),    // True if we only want to see the flux from the second hot spot (does best with normalize_flux = false)
    	 pd_neg_soln(false),
    background_file_is_set(false);
    		
  // Create LightCurve data structure
  class LightCurve curve, normcurve;  // variables curve and normalized curve, of type LightCurve
  //class LightCurve curve;
  class LightCurve *flxcurve;
  class DataStruct obsdata;           // observational data as read in from a file

 

  
  /*********************************************************/
  /* READING IN PARAMETERS FROM THE COMMAND LINE ARGUMENTS */
  /*********************************************************/
    
    for ( int i(1); i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {

        case 'a': // Attenuation flag, corresponds to four attenuation files for NICER targets
        			sscanf(argv[i+1], "%u", &attenuation);
        			// 0 = nothing happens, light curve produced as normal
        			// 1 = NICER J0030 WABS
        			// 2 = NICER J0030 TBABS
        			// 3 = NICER J0437 WABS
        			// 4 = NICER J0437 TBABS
        			// 5 = NICER J0437 TBnew (with fine bands)
        			break;

        case 'A': // ISM column density, in multiples of 1e18cm^2
        			sscanf(argv[i+1], "%lf", &nh);
        			break;
	            
	    case 'b': // Bending Angle File
	            	sscanf(argv[i+1], "%s", bend_file);	
					bend_file_is_set = true;
	            	break;

	    case 'B': // Phase of second spot, 0 < phase_2 < 1 (antipodal = 0.5)
	    			sscanf(argv[i+1], "%lf", &phase_2);
				phase_2 *= -1.0;
	    			break;

	    case 'C': // Temperature of second spot
	    			sscanf(argv[i+1], "%lf", &spot2_temperature);
	    			break;

	    case 'd': // Size of second spot
	    			sscanf(argv[i+1], "%lf", &rho2);
	    			break;	
	            
	    case 'D':  // Distance to NS in kpc
	            	sscanf(argv[i+1], "%lf", &distance);
			distance *= 3.08567758149e19; // Convert to metres -- Correct factor
			//distance *= 3.1e19;
	            	break;
	                
	    case 'e':  // Emission angle of the spot (degrees), latitude
	                sscanf(argv[i+1], "%lf", &theta_1);
			theta_1 *= Units::PI/180.0;
	                theta_is_set = true;
	                break;

		case 'E':  // Offset emission angle of the second spot (degrees), latitude
	                sscanf(argv[i+1], "%lf", &d_theta_2);
	                break;      

	    case 'f':  // Spin frequency (Hz)
	                sscanf(argv[i+1], "%lf", &omega);
	                omega_is_set = true;
	                break;

	    case 'g':  // Spectral Model, beaming (graybody factor)
	                sscanf(argv[i+1],"%u", &beaming_model);
	                // 0 = BB, no beaming
	                // 2 = Fake Spectral Line
	                // 3 = NSATMOS Hydrogen
	                // 4 = NSX Helium
	                // 5 = NSX Hydrogen
					// 7 = Hopf Function
			// 10 = McPhac 
			// 11 = NSX Hydrogen (New version by Wynn Ho)
			// 12 = Blackbody Test
			// 13 = Hopf function Test
			// 14 = BB + Hopf Test
	                break;
	                
	    case 'i':  // Inclination angle of the observer (radians)
	                sscanf(argv[i+1], "%lf", &incl_1);
			incl_1 *= Units::PI/180.0;
	                incl_is_set = true;
	                break;
	            
	    case 'I': // Name of input file
	            	sscanf(argv[i+1], "%s", data_file);
	            	datafile_is_set = true;
	            	break;
	                
	    case 'j': // Diffuse Sky Background
	                sscanf(argv[i+1], "%lf", &dsbackground);
	                break;

	    case 'J': // AGN Background
	    	        sscanf(argv[i+1], "%lf", &agnbackground);
	    	        break;
	          	          
	    case 'l':  // Time shift, phase shift.
	                sscanf(argv[i+1], "%lf", &ts);
	                //ts *= -1;
	                break;

	    case 'L': // Powerlaw for Powerlaw Background
	      sscanf(argv[i+1], "%lf", &powerlaw);
	      break;

	    case 'k': // Background in low energy band (between 0 and 1)
	      			sscanf(argv[i+1], "%lf", &background);
	      			break;

	    case 'K': // Background file specified for each band (between 0 and 1)
	      			sscanf(argv[i+1], "%s", background_file);
	      			background_file_is_set = true;
				std::cout << "Background file is set to : " << background_file << std::endl;
	      			break;
	          	          
	    case 'm':  // Mass of the star (solar mass units)
	                sscanf(argv[i+1], "%lf", &mass);
	                mass_is_set = true;
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
	          
	    case 'p':  // Angular Radius of spot (radians)
	                sscanf(argv[i+1], "%lf", &rho);
	                break;

	    case 'P':  // Spot shape model 0=default=standard; 1=circular in rotating frame; 2=other
	                sscanf(argv[i+1], "%u", &spotshape);
	                break;	          

	    case 'q':  // Oblateness model (default is 1)
	                sscanf(argv[i+1], "%u", &NS_model);
	                model_is_set = true;
	                break;

	    case 'R':  // Instrument Response Curve
	    			sscanf(argv[i+1], "%u", &inst_curve);
	    			break;
	      	          
	    case 'r':  // Radius of the star at the equator(km)
	                sscanf(argv[i+1], "%lf", &req);
	                rspot_is_set = true;
	                break;

	    case 's':  // Spectral Model
	                sscanf(argv[i+1], "%u", &spectral_model);
					// 0 = blackbody
					// 1 = thin line for nicer
	        		// 2 = *slow and old* integrated flux for energy bands
	        		// 3 = *fast and correct* integrated flux for energy bands
					break;
			
	    case 'S': // Number of Energy bands
	      			sscanf(argv[i+1], "%u", &numbands);
	      			break;

	    case 't':  // Number of theta bins 
	                sscanf(argv[i+1], "%u", &numtheta);
			numtheta_in = numtheta;
	                break;
	          
	    case 'T':  // Temperature of the spot, in the star's frame, in keV
	                sscanf(argv[i+1], "%lf", &spot_temperature);
	                break;
	            
	    case 'u': // Lower limit of first energy band, in keV
	            	sscanf(argv[i+1], "%lf", &E_band_lower_1);
	            	//E_band_lower_1_set = true;
	            	break;
	            
	    case 'U': // Upper limit of first energy band, in keV
	            	sscanf(argv[i+1], "%lf", &E_band_upper_1);
	            	//E_band_upper_1_set = true;
	            	break;
	                
	    case 'v': // NICER funny line E1
	            	sscanf(argv[i+1], "%lf", &L1);
	            	break;
	            	
	    case 'V': // NICER funny line E2
	            	sscanf(argv[i+1], "%lf", &L2);
	            	break;
	            
	    case 'x': // Part of funny NICER line
	            	sscanf(argv[i+1], "%lf", &E0);
	            	break;
	            
	    case 'X': // Part of funny NICER line
	            	sscanf(argv[i+1], "%lf", &DeltaE);
	            	break;
	            	
	    case 'z': // Input file for temperature mesh
	            	sscanf(argv[i+1], "%s", T_mesh_file);
	            	T_mesh_in = true;
	            	break;
			
	    case 'Z': // Observation time (in seconds)
	      			sscanf(argv[i+1], "%lf", &obstime);
	      			break;

	            	
	    case '2': // If the user want two spots
	            	two_spots = true;
	            	break;
	            	
	            case '3': // Header for file name
	            	sscanf(argv[i+1],"%s", filenameheader);
	            	break;
	            
	            case '8': // Param_degen gave a negative solution
	            	pd_neg_soln = true;
	            	break;
	            
                case 'h': default: // Prints help
      	            std::cout << "\n\nSpot help:  -flag description [default value]\n" << std::endl
       	            		  
                              << "-a Attenuation Flag [0]:" << std::endl
                              << "      1 for NICER J0030 WABS model" << std::endl
                              << "      2 for NICER J0030 TBABS model" << std::endl
                              << "      3 for NICER J0437 WABS model" << std::endl
                              << "      4 for NICER J0437 TBABS model" << std::endl
                              << "      5 for NICER J0437 'fine channels' TBnew model" << std::endl
                              << "-A ISM column density, in multiples of base value [0]" << std::endl
                              << "      base value is 1e18 " << std::endl
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
                              << "      1 for BB + graybody" << std::endl
                              << "      2 for NICER fake spectral lines" << std::endl
                              << "      3 for NSATMOS (Khaled's)" << std::endl
                              << "      4 for *Old* NSX Helium" << std::endl
                              << "      5 for *limited T/g* NSXH" << std::endl
                              << "      6 for BB * (1 - cosbeta*coseta^2)" << std::endl
                              << "      7 for BB Hopf Function" << std::endl
                              << "      8 for *limited T/g* Slavko's McPHAC" << std::endl
                              << "      9 for *limited T/g* NSX Helium" << std::endl
                              << "      10 for Cole's McPHAC correct E/T version" << std::endl
                              << "      11 for Wynn's full NSXH table" << std::endl
		                      << "-i * Inclination of observer, in degrees, between 0 and 90." << std::endl
                              << "-I Input filename." << std::endl
		                      << "-j Diffuse Sky Background, normalized to expected counts/second in NICER [0.0]" << std::endl
		                      << "-J AGN Background, normalized to expected counts/second in NICER [0.0]" << std::endl
		                      << "-k Constant background for all bands" << std::endl
		                      << "-K Listed backgrounds for each band" << std::endl
		                      << "-l Time shift (or phase shift), in seconds." << std::endl
		                      << "-L Offset inclination of observer, in degrees, between 0 and 90." << std::endl
		                      << "-m * Mass of star in Msun." << std::endl          
      	  	                  << "-n Number of phase or time bins. [128]" << std::endl
      	  	                  << "-N Flag for normalizing the flux. Using this sets it to true. [false]" << std::endl
		                      << "-o Output filename." << std::endl
		                      << "-O Name for second output file, file has alternative format" << std::endl
		                      << "-p Angular radius of spot, rho, in radians. [0.0]" << std::endl
		                      << "-P Hot spot shape model [0]:" << std::endl
		                      << "      0 for standard circular" << std::endl
		                      << "      1 for circular in rotating frame" << std::endl
		                      << "-q * Oblateness Model of star: [1]" << std::endl
		                      << "      1 for Neutron/Hybrid quark star poly model" << std::endl
		                      << "      2 for CFL quark star poly model" << std::endl
		                      << "      3 for spherical model" << std::endl
		                      << "-R Instrument response curve [0]:" << std::endl
		                      << "      0 nothing's done as if all channels are 1 cm^2" << std::endl
		                      << "      1 NICER 2014 'fine channels' effective areas" << std::endl
		                      << "-r * Radius of star (at the equator), in km." << std::endl
		                      << "-s Spectral model of radiation: [0]" << std::endl
		                      << "      0 for monochromatic." << std::endl
		                      << "      1 for NICER fake lines" << std::endl
		                      << "      2 for BB integrated" << std::endl
		                      << "      3 for Atmosphere integrated" << std::endl
		                      << "-S Number of energy bands: [NCURVES]" << std::endl
		                      << "-t Number of theta bins for large spots. Must be < 30. [1]" << std::endl
		                      << "-T Temperature of the spot, in keV. [0]" << std::endl
		                      << "-u Energy bands' lower limit, in keV. [2]" << std::endl
		                      << "-U Energy bands' upper limit, in keV. [3]" << std::endl
		                      << "-v NICER funny line L1" << std::endl
		                      << "-V NICER funny line L2" << std::endl
		                      << "-x NICER funny line E0" << std::endl
		                      << "-X NICER funny line DeltaE" << std::endl
		                      << "-z Input file name for temperature mesh." << std::endl
		                      << "-Z Observation time in seconds [1.0]" << std::endl
		                      << "-2 Flag for calculating two hot spots. Using this sets it to true. [false]" << std::endl
		                      << " Note: '*' next to description means required input parameter." << std::endl
		                      << std::endl;
	                return 0;
            } // end switch	
        } // end if
    } // end for


    
    
    if (bend_file_is_set){ // Read in table of bending angles for all M/R

      // ReadBend allocates memory for the look up tables
      // Reads in tables for b, psi, d(cosalpha)/d(cospsi), toa
      // All depend on M/R
      // Memory must be freed at the end of the program!!!!

      curve = ReadBend(&curve,bend_file); 

    }

    /***********************************************/
    /* CHECKING THAT THE NECESSARY VALUES WERE SET */
    /***********************************************/
    
    if( !( incl_is_set && theta_is_set
	    && mass_is_set && rspot_is_set
	    && omega_is_set && model_is_set ) ) {
        throw( Exception(" Not all required parameters were specified. Exiting.\n") );
        return -1;
    }
    
   
    
    /******************************************/
    /* SENSIBILITY CHECKS ON INPUT PARAMETERS */
    /******************************************/
    
    if ( numtheta < 1 ) {
        throw( Exception(" Illegal number of theta bins. Must be positive. Exiting.\n") );
        return -1;
    }
    if ( spot_temperature < 0.0 || rho < 0.0 || mass < 0.0 || rspot < 0.0 || omega < 0.0 ) {
        throw( Exception(" Cannot have a negative temperature, spot size, NS mass, NS radius, spin frequency. Exiting.\n") );
        return -1;
    }
    if ( rho == 0.0 && numtheta != 1 ) {
        std::cout << "\nWarning: Setting theta bin to 1 for a trivially sized spot.\n" << std::endl;
        numtheta = 1;
    }
   
    if ( numbins > MAX_NUMBINS || numbins <= 0 ) {
    	throw( Exception(" Illegal number of phase bins. Must be between 1 and MAX_NUMBINS, inclusive. Exiting.\n") );
    	return -1;
    }
    
 

    /*****************************************************/
    /* UNIT CONVERSIONS -- MAKE EVERYTHING DIMENSIONLESS */
    /*****************************************************/

    std::cout << "Spin Frequency = " << omega << " Hz" << std::endl;
    std::cout << "Mass = " << mass << " MSUN" << std::endl;
    std::cout << "R_eq = " << req << " km" << std::endl;
    std::cout << "Distance = " << distance << "m" << std::endl;

    //mass_over_req = 1.0/6.1138;

    mass_over_req = mass/(req) * Units::GMC2 * 1;
    std::cout << "GM/(R_eqc^2) = " << mass_over_req << std::endl;
    std::cout << "R/M = " << 1.0/mass_over_req << std::endl;
 
    //incl_1 *= (Units::PI / 180.0);  // radians
    //d_incl_2 *= (Units::PI / 180.0);  // radians
    //theta_1 *= (Units::PI / 180.0); // radians
    //d_theta_2 *= (Units::PI / 180.0);  // radians
    theta_2 = theta_1+d_theta_2; // radians
    mu_1 = cos( theta_1 );
    mu_2 = mu_1; 
    mass = Units::cgs_to_nounits( mass*Units::MSUN, Units::MASS );
    req = Units::cgs_to_nounits( req*1.0e5, Units::LENGTH );
   

    std::cout << "R/M = " << req/mass << std::endl;
    //mass_over_req = mass/req;

    omega = Units::cgs_to_nounits( 2.0*Units::PI*omega, Units::INVTIME );
    distance = Units::cgs_to_nounits( distance*100, Units::LENGTH );
    rot_par = pow(omega*req,2)/mass_over_req;

 
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
    curve.para.aniso = aniso;
    curve.para.Gamma1 = Gamma1;
    curve.para.Gamma2 = Gamma2;
    curve.para.Gamma3 = Gamma3;
    curve.para.temperature = spot_temperature;
    curve.para.ts = ts;
    curve.para.E_band_lower_1 = E_band_lower_1;
    curve.para.E_band_upper_1 = E_band_upper_1;
    curve.para.E_band_lower_2 = E_band_lower_2;
    curve.para.E_band_upper_2 = E_band_upper_2;
    curve.para.E0 = E0;
    curve.para.distance = distance;
    curve.numbins = numbins;

    curve.flags.ignore_time_delays = ignore_time_delays;
    curve.flags.spectral_model = spectral_model;
    curve.flags.beaming_model = beaming_model;
    curve.flags.ns_model = NS_model;
    curve.flags.bend_file = bend_file_is_set;
    curve.flags.attenuation = attenuation;
    curve.flags.inst_curve = inst_curve;
    curve.numbands = numbands;
    curve.fbands = FBANDS;
    curve.tbands = NCURVES;
    curve.cbands = CBANDS;

    curve.flags.spotshape = spotshape;

    double **atten;
    atten = dmatrix(0,400,0,1191);
    double *tbnew, *tbnew_nh;
    tbnew_nh = dvector(0,1191);
    tbnew = dvector(0,NCURVES);
    int inh;

    // Read in the attenuation file
    if (curve.flags.attenuation == 5){
      std::cout << "Read in ISM" << std::endl;
      std::ifstream ism;       // input stream
      ism.open("ISM/tbnew_full.txt");
      char line[265]; // line of the data file being read in
      double get_nh, get_e, get_att;
      for (unsigned int k(1); k<= 400; k++){
	for (unsigned int i(0); i<1191; i++){
	  ism.getline(line,265);
	  //std::cout << "line = " << line << std::endl;
	  sscanf( line, "%lf %lf %lf", &get_nh, &get_e, &get_att );
	  atten[k][i] = get_att;
	}
      }
      std::cout << "nh = " << nh << "e18 cm^2" << std::endl;
      inh = nh;
      std::cout << "inh = " << inh  << std::endl;
      if (inh == 0)
	curve.flags.attenuation = 0;
      else{
      
	//	std::cout << "index = " << 2+(inh-10)/5 << std::endl;
      
	if (inh < 10 && inh > 0)
	  tbnew_nh = atten[1];      
	if (inh >= 10){
	  inh = 2+(nh-10)/5;
	  tbnew_nh = atten[inh];
	}
      }

      for (unsigned int p(0); p<NCURVES; p++){

	if (p%2==0){ //p is even
	  tbnew[p] = 0.25*(3.0*tbnew_nh[p/2] + tbnew_nh[p/2+1]);
	}
	else{ //p is odd
	  tbnew[p] = 0.25 * (tbnew_nh[p/2] + 3.0*tbnew_nh[p/2+1]);
	}	
      }
      
	std::cout << "tbnew[0] = " << tbnew[0] << std::endl;
 

      free_dmatrix(atten,0,400,0,1191);
    }
    //    free_dvector(tbnew_nh,0,1191);      
   
    // Read in the attenuation file
    if (curve.flags.attenuation == 6){
      std::cout << "Read in ISM (approximate version) " << std::endl;
      std::ifstream ism;       // input stream
      ism.open("ISM/tbnew/tbnew0.004.txt");
      char line[265]; // line of the data file being read in
      double get_nh, get_e, get_att;
      int m=0;
	for (unsigned int i(0); i<1191; i++){
	  ism.getline(line,265);
	  //std::cout << "line = " << line << std::endl;
	  sscanf( line, "%lf %lf %lf", &get_nh, &get_e, &get_att );
	  atten[m][i] = get_att;
	    //pow(get_att,nh);
	}

      std::cout << "nh = " << nh << " x 0.4 x e20 cm^2" << std::endl;
      inh = nh; 
      std::cout << "inh = " << inh  << std::endl;
      if (inh == 0)
	curve.flags.attenuation = 0;
      else{
      
	for (unsigned int p(0); p<NCURVES; p++){

	  if (p%2==0){ //p is even
	    tbnew[p] = 0.25*(3.0*atten[0][p/2] + atten[0][p/2+1]);
	  }
	  else{ //p is odd
	    tbnew[p] = 0.25 * (atten[0][p/2] + 3.0*atten[0][p/2+1]);
	  }
	  
	  tbnew[p] *= pow(tbnew[p],nh-1);

	 

	}
	std::cout << "tbnew[0] = " << tbnew[0] << std::endl;
      
	free_dmatrix(atten,0,400,0,1191);
      }
    //    free_dvector(tbnew_nh,0,1191);      
    }

   // Define the Atmosphere Model
    if (curve.flags.beaming_model == 10){ // *cole* McPHAC Hydrogen Atmosphere

      char atmodir[1024], cwd[1024];
      getcwd(cwd, sizeof(cwd));
		
      std::cout << "Reading in McPhac Binary File " << std::endl;

      //sprintf(atmodir,"%s/atmosphere",cwd); Commented this out to see if this is a problem
      chdir(atmodir);
      FILE *Hspecttable;
      Hspecttable=fopen("Hatm8000dT0.05.bin","rb");
      curve.mccangl = dvector(0,51);
      curve.mcloget = dvector(0,101);
      curve.mccinte = dvector(0,1595001);
      int e_index(0);
      double junk(0.0), junk2(0.0), junk3(0.0), junk4(0.0);
      for (int i = 0; i < 50; i++){
	fread(&junk,sizeof(double),1,Hspecttable);
	fread(&junk2,sizeof(double),1,Hspecttable);
	fread(&junk3,sizeof(double),1,Hspecttable);
	fread(&curve.mccangl[i],sizeof(double),1,Hspecttable);    		
	fread(&curve.mccinte[i],sizeof(double),1,Hspecttable);

	//setprecision(10);
		
	/*	std::cout << "i = " << i 
		<< " junk = " << junk
		<< " junk2 = " << junk2
		<< " junk3 = " << junk3
		<< " costheta=" << curve.mccangl[i]
		<< " thing = " << curve.mccinte[i]
		<< std::endl;*/
		
	if (i%50 == 0){
	  curve.mcloget[e_index] = junk3;
	  //	  std::cout << "e_index = " << e_index << " log(e/t)=" << curve.mcloget[e_index] << std::endl;
	  e_index ++;
	}
      }
      
      for (int i = 50; i < 1595001; i++){
	fread(&junk,sizeof(double),1,Hspecttable);
	fread(&junk2,sizeof(double),1,Hspecttable);
	fread(&junk3,sizeof(double),1,Hspecttable);
	fread(&junk4,sizeof(double),1,Hspecttable);    		
	fread(&curve.mccinte[i],sizeof(double),1,Hspecttable);
	
	if (i%50==0 && e_index < 100){
	  curve.mcloget[e_index] = junk3;
	  //std::cout << "e_index = " << e_index << " log(e/t)=" << curve.mcloget[e_index] << std::endl;
	  e_index ++;
	}

      }
      fclose(Hspecttable);
      chdir(cwd);
      std::cout << "finished reading Cole's McPHAC" << std::endl;
    } // End Spectrum Option 10


    // NSX 
    
      if (curve.flags.beaming_model == 11){ // Wynn Ho's NSX-H atmosphere

	double loget,mu,logt,logg,logi;

	
	std::cout << "Reading in NSX version 200804" << std::endl;
	
	//	std::ifstream Hspecttable; 
	//	Hspecttable.open( "atmosphere/nsx_H_v171019.out" );  // current version of NSX

      FILE *Hspecttable;
      Hspecttable=fopen("../atmosphere/nsx_H_v200804.out","r");
      std::cout << "Done opening file" << std::endl;
	  
	NlogTeff = 35;
	curve.mclogTeff = dvector(0,NlogTeff);
	for (int i=0;i<NlogTeff;i++){
	  curve.mclogTeff[i] = 5.10 + i*0.05;	  
	}

	Nlogg = 14; // updated 20221116
	curve.mclogg = dvector(0,Nlogg);
	for (int i=0;i<Nlogg;i++){
	  curve.mclogg[i] = 13.7 + 0.1*i;
	}

	NlogE = 166;
	curve.mcloget = dvector(0,NlogE);
	for (int i=0;i<NlogE;i++){
	  curve.mcloget[i] = -1.30 + i*0.02;
	}

	Nmu = 67;
	Npts =  (NlogTeff*Nlogg*NlogE);

	curve.NlogTeff = NlogTeff;
	curve.Nlogg = Nlogg;
	curve.NlogE = NlogE;
	curve.Nmu = Nmu;
	curve.Npts = Npts;

	std::cout << "Npts = " << Npts << std::endl;

	curve.mccangl = dvector(0,Nmu);
	curve.mccinte = dvector(0,Npts*Nmu);

	std::cout << "Finished allocating memory " << std::endl;
	
    	for (int i = 0; i < Npts*Nmu; i++){
	  
	  fscanf(Hspecttable,"%lf %lf %lf %lf %lf\n",
		 &loget, &mu,
		 &logi, &logt, &logg);	
 
	  curve.mccinte[i] = logi;
	  //std::cout << " logi = " << curve.mccinte[i] << std::endl;
		    
	  if ( i < Nmu ) {
	    curve.mccangl[i] = mu;
	    //std::cout << " mu = " << curve.mccangl[i] << std::endl;
	  }
	  
	  if (logt == 5.1 && mu == 0.5 && loget == 1.0 && logg == 15)
	    std::cout
	    << " i = " << i
	    << " i%(Nmu) = " << i%(Nmu)
	    << " logT = " << logt << " = " << curve.mclogTeff[i/(Nlogg*NlogE*Nmu)]
	    << " logg = " << logg << " = " << curve.mclogg[i/(NlogE*Nmu*NlogTeff)]
	    << " logE = " << loget << " = " << curve.mcloget[i/(Nmu*Nlogg*NlogTeff)]
	    << " cos(theta) = " << curve.mccangl[i%(Nmu)]
		    << " I = " << curve.mccinte[i]
		    << std::endl;
	  
	}
       
	fclose(Hspecttable);
    	std::cout << "finished reading Wynn Ho's NSX-H" << std::endl;

      } // End Spectrum Option 11


    /*************************/
    /* OPENING THE DATA FILE */
    /*************************/
   
    if ( datafile_is_set ) {
      std::cout << "Opening data file" << std::endl;	
      std::ifstream data; //(data_file);      // the data input stream
      data.open( data_file );  // opening the file with observational data

      if ( data.fail() || data.bad() || !data ) {
      std::cout << "fail in loading data" << std::endl;
	throw( Exception("Couldn't open data file."));
	return -1;
      }
      // need to do the following to use arrays of pointers (as defined in Struct)
      for (unsigned int y(0); y < 300; y++) {
	//for (unsigned int y(0); y < numbands; y++) {
      	//std::cout << "setting pointers" << std::endl;
	obsdata.t = new double[databins];
	obsdata.f[y] = new double[databins];
	obsdata.err[y] = new double[databins];
      }
 
      /****************************************/
      /* READING IN FLUXES FROM THE DATA FILE */
      /****************************************/
    
    if (data.is_open()){
      std::cout << "reading data" << std::endl;
      std::cout << "number of data bins " << databins << " number of bands " << numbands << std::endl;
    	double temp;
	for (unsigned int k = 25; k < 300; k++){
	  //	  	for (unsigned int k = 0; k < numbands; k++){
	  for (unsigned int j = 0; j < databins; j++){
	    data >> temp;
    		
	    //std::cout << "band = " << k; 
	    data >> temp;
	    obsdata.t[j] = temp;
	    //	std::cout << "data t["<<j <<"] = " << obsdata.t[j] ;

		data >> temp;
		obsdata.f[k][j] = temp;
		//	std::cout << " f[k][j] = " << temp << std::endl;
		//data >> temp;
		obsdata.err[k][j] = sqrt(temp);
		//std::cout << " err[k][j] = " << temp << std::endl;
	  }
	  //data >> temp;
    	}
    	data.close();
    }

  
       
      // Read in data file to structure "obsdata" (observed data)
      // f[1][i] flux in low energy band
      // f[2][i] flux in high energy band
      obsdata.numbins = databins;
      obsdata.numbands = 300;

      data.close();
      //ts = obsdata.t[0]; // Don't do this if you want manually setting ts to do anything!!
      //obsdata.shift = obsdata.t[0];
      obsdata.shift = ts;

      //std::cout << "Finished reading data from " << data_file << ". " << std::endl;
    } // Finished reading in the data file
	

    /******************************************/
    /*      OPEN INSTRUMENT RESPONSE CURVE    */
    /******************************************/

    std::cout << "Read in the Instrumental Response: ints_curve = " << curve.flags.inst_curve << std::endl;


    if (curve.flags.inst_curve == 2){ // May 2014 version of response matrix

      	std::ifstream file;
	file.open("Area/NICER_May2014_rsp.txt");
	double elow,ehigh;
	
	if(file.is_open()) {
	  for (unsigned int p(0);p<NCURVES;p++){
	    file >> elow;
	    file >> ehigh;
	    file >> start[p];
	    for (unsigned int j(0); j<=76; j++){
	      file >> response[p][j]; 
	    }
	  }
	}else{
        throw( Exception( "instrument response curve file is not found" ));
        }
        file.close();


	
    }

        if (curve.flags.inst_curve == 3){ //Version 1.02 (March 2018) of response matrix

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
	}else{
        throw( Exception( "instrument response curve file is not found" ));
        }
        file.close();

	//std::ifstream file;
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
	else{
        throw( Exception( "off-axis response file is not found" ));
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
    }



    /***************************/
    /* START SETTING THINGS UP */
    /***************************/ 

  

    /*********************************************************************************/
    /* Set up model describing the shape of the NS; oblate, funky quark, & spherical */
    /*********************************************************************************/
	
    OblModelBase* model;
    if ( NS_model == 1 ) { // Oblate Neutron Hybrid Quark Star model
        // Default model for oblate neutron star

      std::cout << " Oblate Neutron Star" << std::endl;
      model = new PolyOblModelNHQS( req,
		   		    mass_over_req,
				    rot_par );

      std::cout << " r_pole = " <<  model->R_at_costheta(1.0) << std::endl;

       
    }
    else if ( NS_model == 2 ) { // Oblate Colour-Flavour Locked Quark Star model
        // Alternative model for quark stars (not very different)
        model = new PolyOblModelCFLQS( req,
				     mass_over_req,
				     rot_par );
        //printf("Oblate Colour-Flavour Locked Quark Star. ");
    }
    else if ( NS_model == 3 ) { // Standard spherical model
        // Spherical neutron star
      //req = rspot;
      rspot = req;
        model = new SphericalOblModel( rspot );
        std::cout << "Spherical Model. " << std::endl;
    }
    else {
        throw(Exception("\nInvalid NS_model parameter. Exiting.\n"));
        return -1;
    }

   
    /****************************/
    /* Initialize time and flux */
    /****************************/

    //std::cout << "Number of Time bins = " << numbins << " " << "Number of Energy Bands = " << numbands << std::endl;
    //std::cout << "NCURVES=" << NCURVES << std::endl;
    flxcurve = &normcurve;

    for ( unsigned int i(0); i < numbins; i++ ) {
        curve.t[i] = i / (1.0 * numbins);// + ts;  // defining the time used in the lightcurves
        for ( unsigned int p(0); p < curve.tbands; p++ ) {
            curve.f[p][i] = 0.0;
	}
    } 

  

    /*********************************/
    /* FIRST HOT SPOT - STANDARD CASE*/
    /*********************************/

		
    if ( T_mesh_in ) {
    	std::cout << "WARNING: code can't handle a spot asymmetric over the pole with a temperature mesh." << std::endl;
    	spot_temperature = 2;
    }
    curve.para.temperature = spot_temperature;

    int pieces;
    //Does the spot go over the pole?
    if ( rho > theta_1){ // yes
    	pieces=2;
    }
    else{ //no
    	pieces=1;
    }
    // If spot is in 2 pieces, p=0 is the crescent; p=1 is the symmetric part over the pole

    if (numbands != 1)
      numbands = NCURVES;

    if ( rho <= 1e-2 )
      numtheta = 1;

    //std::cout << "A: flxcurvef = " << flxcurve->f[0][0] << std::endl;
    
    for (unsigned int p(0);p<pieces;p++){
      //    for (unsigned int p(0);p<1;p++){
      	curve = SpotShape(pieces,p,numtheta,theta_1,rho, &curve, model);
      	double deltatheta(0.0);

	// Looping through the mesh of the spot
      	for (unsigned int k(0); k < numtheta; k++) { // Loop through the circles

	  //  for (unsigned int k(1); k < 2; k++) { // Loop through the circles

	  std::cout << "Theta_spot = " << theta_1 *180/Units::PI
		    << " k = " << k
		    << " theta_k = " << curve.para.theta_k[k]
		    << std::endl;

	  
	  deltatheta = curve.para.dtheta[k];
	  double thetak = curve.para.theta_k[k];
	  double phi_edge;
	  phi_edge = curve.para.phi_k[k];

	  dphi = 2.0*Units::PI/(numbins*1.0);
	  
	  mu_1 = cos(thetak);
	  if (fabs(mu_1) < DBL_EPSILON) mu_1 = 0.0;

	  if (NS_model != 3)
	    rspot = model->R_at_costheta(mu_1);
	  else
	    rspot = req;

	  // Values we need in some of the formulas.
	  cosgamma = model->cos_gamma(mu_1);   // model is pointing to the function cos_gamma
	  curve.para.cosgamma = cosgamma;

	  curve.para.radius = rspot; // load rspot into structure
	  curve.para.mass_over_r = mass_over_req * req/rspot;

	  OblDeflectionTOA* defltoa = new OblDeflectionTOA(model, mass, curve.para.mass_over_r , rspot); 
	  curve = Bend(&curve,defltoa);

	    numphi = 2.0*phi_edge/dphi;
	    phishift = 2.0*phi_edge - numphi*dphi;
	  
	  curve.para.dS = pow(rspot,2) * sin(thetak) * deltatheta * dphi;

	  if (numphi==0){
	    numphi = 1;
	    phishift = 0.0;
	    curve.para.dS = pow(rspot,2) * sin(thetak) * deltatheta * (2.0*phi_edge);
	  }

	  if (spotshape!=2)
	    curve.para.dS *= curve.para.gamma_k[k];

	  curve.para.theta = thetak;

	  SurfaceArea += pow(rspot,2) * sin(thetak) * deltatheta * 2.0*Units::PI / curve.para.cosgamma;

	  if (numtheta==1){  //For a spot with only one theta bin (used for small spot)
	    numphi=1;
	    phi_edge=0.0;
	    dphi=0.0;
	    phishift = 0.0;
	    curve.para.dS = 2.0*Units::PI * pow(rspot,2) * (1.0 - cos(rho));
	    if ( spotshape == 1 ) curve.para.dS /= curve.para.gamma_k[k];
	    if ( spotshape == 0 ) curve.para.dS *= curve.para.gamma_k[k];
	  }
       
    
	  if ( NS_model != 3 ) curve.para.dS /= curve.para.cosgamma;

	  //std::cout << "numphi = " << numphi << std::endl;
	  
	  for ( unsigned int j(0); j < 0.5*numphi ; j++ ) {// looping through the phi divisions
	  //for ( int j(0); j < 10 ; j++ ) {// looping through the phi divisions

	    //curve.para.phi_0 =  -phi_edge + (j+0.5)*dphi;			

	    curve.para.phi_0 = 0.0;
	   
	    
	    //Heart of spot, calculate curve for the first phi bin - otherwise just shift
	    if ( j==0){
	      //	      std::cout << "starting ComputeAngles" << std::endl;

	      curve = ComputeAngles(&curve, defltoa); 	
	      //std::cout << "Spot1: starting ComputeCurve " << std::endl;
	      curve = ComputeCurve(&curve);
	      //std::cout << "Spot1: starting TimeDelays k = " << k << std::endl;
	      curve = TimeDelays(&curve);
	      //std::cout << "Spot1: finished TimeDelays" << std::endl;
	    }


	  
	
	  

	    // Add curves, load into Flux array
	    for ( int i(0); i < numbins; i++ ) {
	      int q(i+j);
	      int qq(i-j);
	      if (q>=numbins){
		//std::cout << "??? q= " << q << std::endl;
		q+=-numbins;
	      }
	      if (qq<0) qq += numbins;
	      for ( unsigned int pp(0); pp < curve.cbands; pp++ ) {
		flxcurve->f[pp][i] += curve.f[pp][q];
		if (j!=0) flxcurve->f[pp][i] += curve.f[pp][qq];
	      }
	    } // ending Add curve
	  } // end for-j-loop

	  
	  // Add in the missing bit.
	  
	  if (phishift != 0.0 ){ // Add light from last bin, which requires shifting
	    for ( unsigned int i(0); i < numbins; i++ ) {
	      int q(i+numphi-1);
	      if (q>=numbins) q+=-numbins;
	      for ( unsigned int pp(0); pp < curve.cbands; pp++ ) {
		Temp[pp][i] = curve.f[pp][q];
	      }
	    }
	    for
	      (unsigned int pp(0); pp < curve.cbands; pp++ )
	      for ( unsigned int i(0); i < numbins; i++ ) 
	    	curve.f[pp][i] = Temp[pp][i];	  	
    
	    curve = ShiftCurve(&curve,phishift);
	    
	    for( unsigned int pp(0); pp < curve.cbands; pp++ )
	      for ( unsigned int i(0); i < numbins; i++ ) {
		//flxcurve->f[pp][i] +=  (Temp[pp][i]+curve.f[pp][i])*phishift/dphi * 0.5 ;
		//flxcurve->f[pp][i] +=  2.0*(curve.f[pp][i])*phishift/dphi ;
			flxcurve->f[pp][i] +=  1.0*(curve.f[pp][i])*phishift/dphi ;
	      }
	  } //end of last bin
	  
	  delete defltoa;
	  
	  //free memory here
	  free_dvector(curve.defl.psi_b,0,301);
	  free_dvector(curve.defl.b_psi,0,301);
	  free_dvector(curve.defl.dcosa_dcosp_b,0,301);
	  free_dvector(curve.defl.toa_b,0,301);
      	} // closing for loop through theta divisions
    } // End Standard Case of first spot

    std::cout << "Finished the first spot!\n" << std::endl;
    
    /**********************************************************/
    /* SECOND HOT SPOT -- Mirroring of first hot spot         */
    /**********************************************************/

    //Setting new parameters for second spot
    if ( two_spots ) {
      numtheta = numtheta_in;
    	std::cout << "Starting Spot 2" << std::endl;
	incl_2 = incl_1;
    	curve.para.incl = incl_2;
 	theta_2 = theta_1 + d_theta_2; //d_theta_2 = 7.78 degrees, for 56+131.78 - 180
	mu_2 = cos(theta_2);
    	curve.para.theta = theta_2;  // keeping theta the same, but changing inclination
    	cosgamma = model->cos_gamma(mu_2);
    	curve.para.cosgamma = cosgamma;

    	if ( T_mesh_in ) {
      		std::cout << "WARNING: code can't handle a spot asymmetric over the pole with a temperature mesh." << std::endl;
      		spot_temperature = 2;
    	}
    	curve.para.temperature = spot2_temperature;
    	int pieces;
    	//Does the spot go over the pole?
    	if ( rho2 > theta_2){ // yes
	  pieces=2;
	}
    	else{ //no
	  pieces=1;
	}
		
	for (unsigned int p(0);p<pieces;p++){
	  curve = SpotShape(pieces,p,numtheta,theta_2,rho2, &curve, model);

	  double deltatheta(0.0);

	  // Looping through the mesh of the spot
	  for (unsigned int k(0); k < numtheta; k++) { // Loop through the circles
	 
	    deltatheta = curve.para.dtheta[k];

	    double thetak = curve.para.theta_k[k];
	    double phi_edge = (Units::PI * 2 * phase_2)-curve.para.phi_k[k]; // for phase difference of 0.5626 cycles

	    //std::cout << "Spot 2: k=" << k << " thetak = " << thetak << std::endl;
	    
	    dphi = 2.0*Units::PI/(numbins*1.0);

	    mu_2 = cos(thetak);
	    if (fabs(mu_2) < DBL_EPSILON) mu_2 = 0.0;

	    if (NS_model != 3)
	      rspot = model->R_at_costheta(mu_2);
	    else
	      rspot = req;

	    // Values we need in some of the formulas.
	    cosgamma = model->cos_gamma(mu_2);
	    curve.para.cosgamma = cosgamma;

	    curve.para.radius = rspot; // load rspot into structure
	    curve.para.mass_over_r = mass_over_req * req/rspot;

		

	    //std::cout << " Entering defltoa M/R = " << curve.para.mass_over_r << std::endl; 

	    OblDeflectionTOA* defltoa = new OblDeflectionTOA(model, mass, curve.para.mass_over_r , rspot); 

	    //std::cout << " Entering Bend M/R = " << curve.para.mass_over_r << std::endl; 


	    curve = Bend(&curve,defltoa);

	    //numphi = 2.0*(Units::PI-phi_edge)/dphi;
	    numphi = 2.0*((Units::PI * 2 * phase_2)-phi_edge)/dphi;
	    //std::cout << numphi << " " << phi_edge << " " << dphi << std::endl;
	    //phishift = 2.0*(Units::PI-phi_edge)- numphi*dphi;
	    phishift = 2.0*((Units::PI * 2 * phase_2)-phi_edge)- numphi*dphi;
	    //phishift = 2.0*(Units::PI-phi_edge)+ numphi*dphi;

	    curve.para.dS = pow(rspot,2) * sin(thetak) * deltatheta * dphi * curve.para.gamma_k[k];
	    curve.para.theta = thetak;

	    if (numtheta==1){  //For a spot with only one theta bin (used for small spot)
	      numphi=1;
	      //phi_edge=-1*Units::PI;
	      phi_edge=-1*(Units::PI * 2 * phase_2);
	      dphi=0.0;
	      phishift = 0.0;
	      curve.para.dS = 2.0*Units::PI * pow(rspot,2) * (1.0 - cos(rho2)) ;
	      if ( spotshape == 1 ) curve.para.dS /= curve.para.gamma_k[k];
	      if ( spotshape == 0 ) curve.para.dS *= curve.para.gamma_k[k];
	    }
	    //std::cout << numphi << " " << phi_edge << " " << dphi << " " << phishift << " " << curve.para.dS << std::endl;
       
	    if ( NS_model != 3 ) curve.para.dS /= curve.para.cosgamma;

	    for ( unsigned int j(0); j < numphi; j++ ) {// looping through the phi divisions
	      curve.para.phi_0 =  phi_edge + (j+0.5)*dphi;			
	      //std::cout << phi_edge << " " << dphi << " " << curve.para.phi_0 << " " << j << " " << numphi << std::endl;
	      //Heart of spot, calculate curve for the first phi bin - otherwise just shift
	      if ( j==0){
			 
		curve = ComputeAngles(&curve, defltoa); 
		curve = ComputeCurve(&curve);
		curve = TimeDelays(&curve);
	      }
	
	      if ( curve.para.temperature == 0.0 ) {// Flux is zero for parts with zero temperature
		for ( unsigned int i(0); i < numbins; i++ ) {
		  for ( unsigned int pp(0); pp < curve.tbands; pp++ ) curve.f[pp][i] = 0.0;
		}
	      }

	      // Add curves, load into Flux array
	      for ( unsigned int i(0); i < numbins; i++ ) {
		int q(i+j);
		if (q>=numbins) q+=-numbins;
		for ( unsigned int pp(0); pp < curve.cbands; pp++ ) {
				  //Flux[p][i] += curve.f[p][q];
				  flxcurve->f[pp][i] += curve.f[pp][q];
		}
	      } // ending Add curves
	    } // end for-j-loop

	    // Add in the missing bit.   
	       
	    if (phishift != 0.0){ // Add light from last bin, which requires shifting
	      for ( unsigned int i(0); i < numbins; i++ ) {
		int q(i+numphi-1);
		if (q>=numbins) q+=-numbins;
		for ( unsigned int pp(0); pp < curve.cbands; pp++ ) {
		  Temp[pp][i] = curve.f[pp][q];
		}
	      }
	      for (unsigned int pp(0); pp < curve.cbands; pp++ )
		for ( unsigned int i(0); i < numbins; i++ ) curve.f[pp][i] = Temp[pp][i];	  
	    
	      curve = ShiftCurve(&curve,phishift);
	      
	      for( unsigned int pp(0); pp < curve.cbands; pp++ )
		for ( unsigned int i(0); i < numbins; i++ ) {
		  //Flux[pp][i] += curve.f[pp][i]*phishift/dphi      ;
		  //flxcurve->f[p][i] += curve.f[pp][i]*phishift/dphi      ;
		  flxcurve->f[pp][i] +=  (Temp[pp][i]+curve.f[pp][i])*phishift/dphi * 0.5 ;
		}
	    } //end of last bin  
	    
	    delete defltoa;

	    //insert free memory here
	    free_dvector(curve.defl.psi_b,0,301);
	    free_dvector(curve.defl.b_psi,0,301);
	    free_dvector(curve.defl.dcosa_dcosp_b,0,301);
	    free_dvector(curve.defl.toa_b,0,301);
	  } // closing for loop through theta divisions
    	} // End Standard Case of second spot
    } // closing second spot

    // Flxcurve holds the current version of the waveform

    // Integrate over the phase bins

    /*for (unsigned int i(0);i<numbins;i++){
	std::cout << "C1. Flux in Bin 0: = " << flxcurve->f[0][i] << std::endl;
	}*/
      
    int binfactor (numbins/databins);
    std::cout << "binfactor = " << binfactor << std::endl;

    for (unsigned int p(0);p<numbands-1;p++){
      for (unsigned int i(0);i<databins;i++){
	curve.f[p][i] = 0.0;
	int index, offset(0);
	
	index = 0 + i*binfactor + offset;
	if (index>=numbins) index -= numbins;
	if (index<0) index += numbins;
	curve.f[p][i] += flxcurve->f[p][index]; // Average over the phase bins


	index++;
	if (index>=numbins) index -= numbins;
	if (index<0) index += numbins;
	curve.f[p][i] += 4.0*flxcurve->f[p][index]; // Average over the phase bins

	index++;
	if (index>=numbins) index -= numbins;
	if (index<0) index += numbins;
	  curve.f[p][i] += 2.0*flxcurve->f[p][index]; // Average over the phase bins

	index++;
	if (index>=numbins) index -= numbins;
	if (index<0) index += numbins;	  
	curve.f[p][i] += 4.0*flxcurve->f[p][index]; // Average over the phase bins

	index++;
	if (index>=numbins) index -= numbins;
	if (index<0) index += numbins;	  
	curve.f[p][i] += flxcurve->f[p][index];
	 
	curve.f[p][i] *= 1.0/(12.0);

	if (p==0)
	  curve.t[i] = curve.t[i*binfactor];

	/*	if (p==0)
	  std::cout << " i=" << i
		    << " t[0]=" << flxcurve->t[i]
		    << " newt[0]=" << curve.t[i]
		    << " f[0]=" << flxcurve->f[p][binfactor*i+0] 
		    << " f[1]=" << flxcurve->f[p][binfactor*i+1]
		    << " new f[0] = " << curve.f[p][i]
		    << std::endl;
	*/
      }
    }
    numbins = databins;
        
	for ( unsigned int i(0); i< databins; i++){
	  for (unsigned int p(0);p<numbands-1;p++){
	    flxcurve->f[p][i] = curve.f[p][i];
	  }
	}

 

	
    /***************************************/
    /* START OF INSTRUMENT EFFECT ROUTINES */
    /***************************************/

    // Interpolate to create all the other energy bands
    
    unsigned int q = 0;
    unsigned int index = 0;
    unsigned int factor = curve.tbands/curve.cbands;

    std::cout << "tbands = " << curve.tbands
	      << " numbins = " << numbins
	      << std::endl;
   
    for (unsigned int p = 0; p < curve.tbands; p++){
      q = (p/factor);
      index = p - (q*factor);
      //std::cout << " p = " << p 
      //	<< "index = " << index << std::endl;
      for (unsigned int i = 0; i < numbins; i++){

	if (index == 0)
	  curve.f[p][i] = flxcurve->f[q][i];
	else{
	  // linear interp
	  curve.f[p][i] = flxcurve->f[q][i] * (factor - index)/(factor*1.0) + flxcurve->f[q+1][i] * index/(factor*1.0); 
	} 
      }
    }


    /**********************************/
    /*       APPLYING ATTENUATION     */
    /**********************************/

    if (curve.flags.attenuation != 0){
      std::cout << "Apply ISM to the Signal!" << std::endl;
      //if (curve.flags.attenuation >= 5)
      //curve = Attenuate(&curve,curve.flags.attenuation,nh,tbnew);
	//if (curve.flags.attenuation == 3)
	//curve = Wabs(&curve, curve.flags.attenuation, nh);
    }
    
      
   
        


    /******************************************/
    /*         ADDING BACKGROUND              */
    /******************************************/

    /* Create Phase-independent Powerlaw Background */
    
    if (powerlaw!=0){
      std::cout << "Create the powerlaw background!" << std::endl;
      (*flxcurve) = PowerLaw_Background(&curve,1.0,powerlaw);
      // Add attenuation to the background
      //if (curve.flags.attenuation == 5)
      //(*flxcurve) = Attenuate(flxcurve,curve.flags.attenuation,nh,tbnew);
    }
      
    if (background_file_is_set){
      //curve = Background_list(&curve, background_file);      
    } 

    else {
      if (background != 0)
	// Add a constant background in all bands
	for (unsigned int p = 0; p < numbands; p++){	   
	  for (unsigned int i = 0; i < numbins; i++){	     
	    curve.f[p][i] += background;
	    if (p==0)
	      std::cout << "Added background: i="<<i << " flux = " << curve.f[p][i] << std::endl;
	  }
	}
    }
	
    if (agnbackground > 0){
      std::cout << "AGN Background" << std::endl;		
      curve = AGN_Background(&curve, agnbackground, nh);
    }

    if (dsbackground > 0){
      std::cout << "Sky Background" << std::endl;
      curve = Sky_Background(&curve, dsbackground);
    }


    // If databins < numbins then rebin the theoretical curve down 
    std::cout << "databins = " << databins << ", numbins = " << numbins << std::endl;
    std::cout << " Rebin the data? " << std::endl;
    if (databins < numbins) {
        obsdata.numbins = databins;
	std::cout << " Rebin the data! " << std::endl;
        curve = ReBinCurve(&obsdata,&curve);
	(*flxcurve) = ReBinCurve(&obsdata,flxcurve);
        numbins = databins;
    }

    //    if (numbands != 1)
    //numbands = curve.fbands;
    if (numbands !=1)
      numbands = curve.tbands;

    std::cout << "numbands = " << numbands << std::endl;

    /******************************************/
    /*  APPLYING INSTRUMENT RESPONSE CURVE    */
    /******************************************/

    //    std::cout << "Apply Instrument Response to Spot: ints_curve = " << curve.flags.inst_curve << std::endl;
    if (curve.flags.inst_curve > 0){
      std::cout << "Applying Instrument Response" << std::endl;


           
      //if (curve.flags.inst_curve == 2)
      //curve = Inst_Res2(&curve, curve.flags.inst_curve,start,response);

      if (curve.flags.inst_curve == 3){ // Apply response to both signal and background
	std::cout << "Apply the Response Matrix " << std::endl;

	// Apply the ARF
	for (unsigned int p = 0; p < numbands; p++){
	  for (unsigned int i = 0; i < numbins; i++){          	   
	    curve.f[p][i] *= arf[p];
	    flxcurve->f[p][i] *= arf[p];
	  }
	}	

	// Redistribution Matrix
	//curve = Inst_Res3(&curve, curve.flags.inst_curve,start,response); // Signal
	//(*flxcurve) = Inst_Res3(flxcurve, curve.flags.inst_curve,start,response); // Background

	// Offaxis Response
	for (unsigned int p = 0; p < 355; p++){
	  for (unsigned int i = 0; i < numbins; i++){          	   
	    curve.f[p][i] *= offaxis[p];
	    flxcurve->f[p][i] *= offaxis[p];
	  }
	}	
      }     
    } // Finished Applying the Response Matrix

    double spotcounts = 0.0;
    double bkgcounts = 0.0;

    if (powerlaw != 0.0){

    /******************************************/
    /*     MULTIPLYING OBSERVATION TIME       */
    /******************************************/
    for (unsigned int p = 0; p < 300; p++){
      //std::cout << "band " << p << std::endl; 
      //std::cout << "spotcounts = " << spotcounts << " f=" << curve.f[p][0] << std::endl;
        for (unsigned int i = 0; i < numbins; i++){          
	  curve.f[p][i] *= obstime/(databins);  
	  //curve.f[p][i] *= obstime;  
	  spotcounts += curve.f[p][i];
	  bkgcounts += flxcurve->f[p][i];
        }
    }
    }


    
    
    //numbands = 300;

    // bkgcounts *= 1.0011;

    if (powerlaw != 0.0){
      double totalcounts = 0.0;
      for (unsigned int p = 0; p < numbands; p++){  
        for (unsigned int i = 0; i < numbins; i++){           
	  //curve.f[p][i] *= 1e6/spotcounts;
	  //curve.f[p][i] += flxcurve->f[p][i] * obstime / 6.114e5;;
	  curve.f[p][i] += flxcurve->f[p][i] * 1e6 / bkgcounts;
	  totalcounts += curve.f[p][i];
        }
      }
      std::cout << "Total Counts = " << totalcounts << std::endl;
    }
    

    /************************************************************/
    /* If data file is set, calculate chi^2 fit with simulation */
    /************************************************************/
    if ( datafile_is_set ) {
      curve.numbands = 300; // For synthetic 2 spot test
    	std::cout << "calculating chi squared" << std::endl;
    	chisquared = ChiSquare ( &obsdata, &curve );
	
	if (background_file_is_set){
	  //normcurve = Read_Background_Guess(&curve, background_file);   
	   chisquared = BackChi ( &obsdata, &curve, &normcurve );	   
	}	
    }
    
    std::cout << "Spot: m = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN 
	      << " Msun, r = " << Units::nounits_to_cgs(req, Units::LENGTH )*1.0e-5 
	      << " km, f = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) 
	      << " Hz, i = " << incl_1 * 180.0 / Units::PI 
	      << ", e = " << theta_1 * 180.0 / Units::PI 
	      << ", X^2 = " << chisquared 
	      << std::endl;    

    std::cout.precision(10);
    std::cout << "Log(Likelihood) = " <<  chisquared << std::endl;



    std::cout << "Spot Counts = " << spotcounts << std::endl;
    std::cout << "Background Counts = " << bkgcounts << std::endl;
    std::cout << "Spot: numbands = " << numbands<< std::endl;

 
    
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
 


      double E_diff, Eobs;
      E_diff = (E_band_upper_1 - E_band_lower_1)/numbands;
      Eobs = E0;

      
      for ( unsigned int p(0); p < numbands; p++ ) {
	for ( unsigned int i(0); i < numbins; i++ ) {
	  if (numbands !=1)
	    Eobs = curve.para.E_band_lower_1 + p*E_diff;
	  out << Eobs << "\t";
	  out << curve.t[i]<< "\t";		
	  out << curve.f[p][i] << "\t";
	  out << p << std::endl;
	}
      }
      	
      /*
      for ( unsigned int p(25); p < 300; p++ ) {
	for ( unsigned int i(0); i < numbins; i++ ) {
	  if (numbands !=1)
	    Eobs = curve.para.E_band_lower_1 + p*E_diff;
	  
	  out << p << "\t"; // NICER Channels
	  out << curve.t[i]<< "\t";	
	  out << curve.f[p][i] - flxcurve->f[p][i] * 1e6 / bkgcounts << "\t" ;	
	  out << curve.f[p][i] << "\t"; // Signal from Spots, with ISM and off-axis response + Background
	  //	  out << flxcurve->f[p][i] *1e6/(bkgcounts) << "\t"; //Background, with off-axis response
	  //      out << flxcurve->f[p][i] << "\t";
	  out << std::endl;
	}
      }

      */
      

    
      out.close();
    



 
    // Free previously allocated memory

     if (curve.flags.beaming_model == 10 || curve.flags.beaming_model >= 12  ){ // *cole* McPHAC Hydrogen Atmosphere        
       free_dvector(curve.mccangl,0,51);
       free_dvector(curve.mcloget,0,101);
       free_dvector(curve.mccinte,0,1595001);
     }
     if (curve.flags.beaming_model == 11){ // New NSX-H model
	free_dvector(curve.mclogTeff,0,NlogTeff);
	free_dvector(curve.mclogg,0,Nlogg);
	free_dvector(curve.mcloget,0,NlogE);
	free_dvector(curve.mccangl,0,Nmu);
	free_dvector(curve.mccinte,0,Npts*Nmu);
     }



    delete model;
    return 0;
} 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}

