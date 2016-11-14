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
#include "PolyOblModelNHQS.h"
#include "PolyOblModelCFLQS.h"
#include "SphericalOblModel.h"
#include "OblModelBase.h"
#include "Units.h"
#include "Exception.h"
#include "Struct.h"
#include "time.h"
#include <string.h>

// MAIN
int main ( int argc, char** argv ) try {  // argc, number of cmd line args; 
                                          // argv, actual character strings of cmd line args

  /*********************************************/
  /* VARIABLE DECLARATIONS AND INITIALIZATIONS */
  /*********************************************/
    
  std::ofstream out;      // output stream; printing information to the output file
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
    bbrat(1.0),                 // Ratio of blackbody to Compton scattering effects, unitless
    ts(0.0),                    // Phase shift or time off-set from data; Used in chi^2 calculation
    spot_temperature(0.0),      // Inner temperature of the spot, in the star's frame, in keV
    rho(0.0),                   // Angular radius of the inner bullseye part of the spot, in degrees (converted to radians)
    dphi(1.0),                  // Each chunk of azimuthal angle projected onto equator, when broken up into the bins (see numphi)
    phishift,
    aniso(0.586),               // Anisotropy parameter
    Gamma1(2.0),                // Spectral index
    Gamma2(2.0),                // Spectral index
    Gamma3(2.0),                // Spectral index
    mu_1(1.0),                  // = cos(theta_1), unitless
    mu_2(1.0),                  // = cos(theta_2), unitless
    cosgamma,                   // Cos of the angle between the radial vector and the vector normal to the surface; defined in equation 13, MLCB
    Flux[NCURVES][MAX_NUMBINS], // Array of fluxes; Each curve gets its own vector of fluxes based on photons in bins.
    Temp[NCURVES][MAX_NUMBINS],
    E_band_lower_1(2.0),        // Lower bound of first energy band to calculate flux over, in keV.
    E_band_upper_1(3.0),        // Upper bound of first energy band to calculate flux over, in keV.
    E_band_lower_2(5.0),        // Lower bound of second energy band to calculate flux over, in keV.
    E_band_upper_2(6.0),        // Upper bound of second energy band to calculate flux over, in keV.
    background[NCURVES],
    T_mesh[30][30],             // Temperature mesh over the spot; same mesh as theta and phi bins, assuming square mesh
    chisquared(1.0),            // The chi^2 of the data; only used if a data file of fluxes is inputed
    distance(3.0857e22),        // Distance from earth to the NS, in meters; default is 10kpc
    B;                          // from param_degen/equations.pdf 2
   

  double E0, E1, E2, DeltaE, rot_par;
    
  unsigned int NS_model(1),       // Specifies oblateness (option 3 is spherical)
    spectral_model(0),    // Spectral model choice (initialized to blackbody)
    beaming_model(0),     // Beaming model choice (initialized to isotropic)
    numbins(MAX_NUMBINS), // Number of time or phase bins for one spin period; Also the number of flux data points
    numphi(1),            // Number of azimuthal (projected) angular bins per spot
    numtheta(1),          // Number of latitudinal angular bins per spot
    spotshape(0), // Spot shape; 0=standard
    numbands(NCURVES); // Number of energy bands;

  



  char out_file[256] = "flux.txt",    // Name of file we send the output to; unused here, done in the shell script
         out_dir[80],                   // Directory we could send to; unused here, done in the shell script
         T_mesh_file[100],              // Input file name for a temperature mesh, to make a spot of any shape
         data_file[256],                // Name of input file for reading in data
	  //testout_file[256] = "test_output.txt", // Name of test output file; currently does time and two energy bands with error bars
         filenameheader[256]="Run";

         
  // flags!
  bool incl_is_set(false),         // True if inclination is set at the command line (inclination is a necessary variable)
    	 theta_is_set(false),        // True if theta is set at the command line (theta is a necessary variable)
    	 mass_is_set(false),         // True if mass is set at the command line (mass is a necessary variable)
    	 rspot_is_set(false),        // True if rspot is set at the command line (rspot is a necessary variable)
    	 omega_is_set(false),        // True if omega is set at the command line (omega is a necessary variable)
    	 model_is_set(false),        // True if NS model is set at the command line (NS model is a necessary variable)
    	 datafile_is_set(false),     // True if a data file for inputting is set at the command line
    	 ignore_time_delays(false),  // True if we are ignoring time delays
    	 T_mesh_in(false),           // True if we are varying spot shape by inputting a temperature mesh for the hot spot's temperature
    	 normalize_flux(false),      // True if we are normalizing the output flux to 1
    	 //E_band_lower_2_set(false),  // True if the lower bound of the second energy band is set
    	 //E_band_upper_2_set(false),  // True if the upper bound of the second energy band is set
    	 two_spots(false),           // True if we are modelling a NS with two antipodal hot spots
    	 only_second_spot(false),    // True if we only want to see the flux from the second hot spot (does best with normalize_flux = false)
    	 pd_neg_soln(false);
		
  // Create LightCurve data structure
  class LightCurve curve, normcurve;  // variables curve and normalized curve, of type LightCurve
  class DataStruct obsdata;           // observational data as read in from a file


  /*********************************************************/
  /* READING IN PARAMETERS FROM THE COMMAND LINE ARGUMENTS */
  /*********************************************************/
    
    for ( int i(1); i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {
	    case 'a':  // Anisotropy parameter
	                sscanf(argv[i+1], "%lf", &aniso);
	                break;
	            
	    case 'b': // Blackbody ratio
	            	sscanf(argv[i+1], "%lf", &bbrat);
	            	break;
	            
	    case 'D':  // Distance to NS
	            	sscanf(argv[i+1], "%lf", &distance);
	            	break;
	                
	    case 'e':  // Emission angle of the spot (degrees), latitude
	                sscanf(argv[i+1], "%lf", &theta_1);
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
	                // 1 = BB + graybody
	                // 2 = Fake Spectral Line
	                // 3 = Hydrogen
	                // 4 = Helium
	                break;
	                
	    case 'i':  // Inclination angle of the observer (degrees)
	                sscanf(argv[i+1], "%lf", &incl_1);
	                incl_is_set = true;
	                break;

	    // Offset inclination angle set as case L
	            
	    case 'I': // Name of input file
	            	sscanf(argv[i+1], "%s", data_file);
	            	datafile_is_set = true;
	            	break;
	                
	    case 'j':  // Flag for calculating only the second (antipodal) hot spot
	            	only_second_spot = true;
	            	break;
	          	          
	    case 'l':  // Time shift, phase shift.
	                sscanf(argv[i+1], "%lf", &ts);
	                break;

	    case 'L':  // Offset inclination angle of the observer (degrees)
	                sscanf(argv[i+1], "%lf", &d_incl_2);
	                break;

	    case 'k': // Background in low energy band (between 0 and 1)
	      sscanf(argv[i+1], "%lf", &background[2]);
	      break;

	    case 'K': // Background in high energy band (between 0 and 1)
	      sscanf(argv[i+1], "%lf", &background[3]);
	      break;
	          	          
	    case 'm':  // Mass of the star (solar mass units)
	                sscanf(argv[i+1], "%lf", &mass);
	                mass_is_set = true;
	                break;
	          
	    case 'n':  // Number of phase or time bins
	                sscanf(argv[i+1], "%u", &numbins);
	                break;
	                
	    case 'N': // Flag for not normalizing the flux output
	            	normalize_flux = true;
	            	break;

	    case 'o':  // Name of output file
	                sscanf(argv[i+1], "%s", out_file);
	                break;

	    case 'O':  case 'A': // This option controls the name of the output directory
	                sscanf(argv[i+1], "%s", out_dir);
	                break;
	          
	    case 'p':  // Angular Radius of spot (radians)
	                sscanf(argv[i+1], "%lf", &rho);
	                break;

	    case 'P':  // Spot shape model 0=default=standard 
	                sscanf(argv[i+1], "%u", &spotshape);
	                break;	          

	    case 'q':  // Oblateness model (default is 1)
	                sscanf(argv[i+1], "%u", &NS_model);
	                model_is_set = true;
	                break;
	      	          
	    case 'r':  // Radius of the star at the equator(km)
	                sscanf(argv[i+1], "%lf", &req);
	                rspot_is_set = true;
	                break;

	    case 's':  // Spectral Model
	                sscanf(argv[i+1], "%u", &spectral_model);
			// 0 = blackbody
			// 1 = thin line for nicer
	        // 2 = integrated flux for energy bands
			break;
			
	    case 'S': // Number of Energy bands
	      sscanf(argv[i+1], "%u", &numbands);
	      break;

	    case 't':  // Number of theta bins 
	                sscanf(argv[i+1], "%u", &numtheta);
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
	            	sscanf(argv[i+1], "%lf", &E1);
	            	break;
	            	
	    case 'V': // NICER funny line E2
	            	sscanf(argv[i+1], "%lf", &E2);
	            	break;
	            
	    case 'x': // Scattering radius, in kpc
	            	sscanf(argv[i+1], "%lf", &E0);
	            	break;
	            
	    case 'X': // Scattering intensity, units not specified
	            	sscanf(argv[i+1], "%lf", &DeltaE);
	            	break;
	            	
	    case 'z': // Input file for temperature mesh
	            	sscanf(argv[i+1], "%s", T_mesh_file);
	            	T_mesh_in = true;
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
                              << "-a Anisotropy parameter. [0.586]" << std::endl
                              << "-b Ratio of blackbody flux to comptonized flux. [1.0]" << std::endl
                              << "-d Ignores time delays in output (see source). [0]" << std::endl
                              << "-D Distance from earth to star, in meters. [~10kpc]" << std::endl
                              << "-e * Latitudinal location of emission region, in degrees, between 0 and 90." << std::endl
                              << "-f * Spin frequency of star, in Hz." << std::endl
                              << "-g Graybody factor of beaming model: 0 = isotropic, 1 = Gray Atmosphere. [0]" << std::endl
		                      << "-i * Inclination of observer, in degrees, between 0 and 90." << std::endl
                              << "-I Input filename." << std::endl
		                      << "-j Flag for computing only the second (antipodal) hot spot. [false]" << std::endl
		                      << "-l Time shift (or phase shift), in seconds." << std::endl
		                      << "-m * Mass of star in Msun." << std::endl          
      	  	                  << "-n Number of phase or time bins. [128]" << std::endl
      	  	                  << "-N Flag for normalizing the flux. Using this sets it to true. [false]" << std::endl
		                      << "-o Output filename." << std::endl
		                      << "-O Name of the output directory." << std::endl
		                      << "-p Angular radius of spot, rho, in radians. [0.0]" << std::endl
		                      << "-q * Model of star: [3]" << std::endl
		                      << "      1 for Neutron/Hybrid quark star poly model" << std::endl
		                      << "      2 for CFL quark star poly model" << std::endl
		                      << "      3 for spherical model" << std::endl
		                      << "-r * Radius of star (at the equator), in km." << std::endl
		                      << "-s Spectral model of radiation: [0]" << std::endl
		                      << "      0 for bolometric light curve." << std::endl
		                      << "      1 for blackbody in monochromatic energy bands (must include T option)." << std::endl
		                      << "-t Number of theta bins for large spots. Must be < 30. [1]" << std::endl
		                      << "-T Temperature of the spot, in keV. [2]" << std::endl
		                      << "-u Low energy band, lower limit, in keV. [2]" << std::endl
		                      << "-U Low energy band, upper limit, in keV. [3]" << std::endl
		                      << "-v High energy band, lower limit, in keV. [5]" << std::endl
		                      << "-V High energy band, upper limit, in keV. [6]" << std::endl
		                      << "-x Scattering radius, in kpc." << std::endl
		                      << "-X Scattering intensity, units unspecified." << std::endl
		                      << "-z Input file name for temperature mesh." << std::endl
		                      << "-2 Flag for calculating two hot spots, on both magnetic poles. Using this sets it to true. [false]" << std::endl
		                      << "-3 File name header, for use with Ferret and param_degen." << std::endl
		                      << "-8 Sets a flag to indicate that param_degen gave a negative solution." << std::endl
		                      << " Note: '*' next to description means required input parameter." << std::endl
		                      << std::endl;
	                return 0;
            } // end switch	
        } // end if
    } // end for

    /***********************************************/
    /* CHECKING THAT THE NECESSARY VALUES WERE SET */
    /***********************************************/
    
    if( !( incl_is_set && theta_is_set
	    && mass_is_set && rspot_is_set
	    && omega_is_set && model_is_set ) ) {
        throw( Exception(" Not all required parameters were specified. Exiting.\n") );
        return -1;
    }
    
   
    /***************************************************************/
    /* READING TEMPERATURE VALUES FROM AN INPUT FILE, IF SPECIFIED */
    /***************************************************************/
	
    if ( T_mesh_in ) {
    	std::ifstream inStream(T_mesh_file);
    	std::cout << "** Tmeshfile name = " << T_mesh_file << std::endl;
    	unsigned int n(0);

    	if ( !inStream ) {
        	throw( Exception( "Temperature mesh input file could not be opened. \nCheck that the actual input file name is the same as the -z flag parameter value. \nExiting.\n") );
        	return -1;
        }
        std::string inLine;
        while ( std::getline(inStream, inLine) )  // getting numtheta from the file
        	n++;
        if (n != numtheta) {
        	throw( Exception( "Numtheta from mesh file does not match numtheta from command line input. \nExiting.\n") );
        	return -1;
        }
        numtheta = numphi = n;
        // go back to the beginning of the file
        inStream.clear();
        inStream.seekg (0, std::ios::beg);
        for ( unsigned int k(0); k < numtheta; k++ ) {
    		for ( unsigned int j(0); j < numphi; j++ ) {
    			inStream >> T_mesh[k][j];
    		}
    	}
    	inStream.close();
    	
    	// Printing the temperature mesh, for checking that it was read in correctly
    	std::cout << "Temperature Mesh of Spot:" << std::endl;
    	for ( unsigned int i(0); i < numtheta; i++ ) {
    		for ( unsigned int j(0); j < numphi; j++ ) {
    	    	std::cout << T_mesh[i][j] << "  ";
    		}
    		std::cout << std::endl;
    	}
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

    mass_over_req = mass/(req) * Units::GMC2;
    incl_1 *= (Units::PI / 180.0);  // radians
    if ( only_second_spot ) incl_1 = Units::PI - incl_1; // for doing just the 2nd hot spot
    theta_1 *= (Units::PI / 180.0); // radians
    theta_2 = theta_1+d_theta_2; // radians
    //rho *= (Units::PI / 180.0);  // rho is input in radians
    mu_1 = cos( theta_1 );
    mu_2 = mu_1; 
    mass = Units::cgs_to_nounits( mass*Units::MSUN, Units::MASS );
    req = Units::cgs_to_nounits( req*1.0e5, Units::LENGTH );
   
    omega = Units::cgs_to_nounits( 2.0*Units::PI*omega, Units::INVTIME );
    distance = Units::cgs_to_nounits( distance*100, Units::LENGTH );
    rot_par = pow(omega*req,2)/mass_over_req;
	
    //std::cout << "Dimensionless: Mass/Radius = " << mass_over_req  << " M/R = " << mass/req << std::endl; 

    /**********************************/
    /* PASS VALUES INTO THE STRUCTURE */
    /**********************************/    	
    
    curve.para.mass = mass;
    curve.para.mass_over_r = mass_over_req;
    curve.para.omega = omega;
    curve.para.radius = req;
    curve.para.req = req;
    curve.para.theta = theta_1;
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
    curve.para.distance = distance;
    curve.numbins = numbins;
    //curve.para.rsc = r_sc;
    //curve.para.Isc = I_sc;

    //numphi = numtheta; // code currently only handles a square mesh over the hotspot
  
    curve.flags.ignore_time_delays = ignore_time_delays;
    curve.flags.spectral_model = spectral_model;
    curve.flags.beaming_model = beaming_model;
    curve.flags.ns_model = NS_model;
    curve.numbands = numbands;
    curve.flags.spotshape = spotshape;


   // Define the Observer's Spectral Model

    if (curve.flags.spectral_model == 0){ // NICER: Monochromatic Obs at E0=1keV
      curve.para.E0 = 1.0;
      curve.numbands = 1;
    }
    if (curve.flags.spectral_model == 1){ // NICER Line
      curve.para.E0 = E0; // Observed Energy in keV
      curve.para.E1 = E1; // Lowest Energy in keV in Star's frame
      curve.para.E2 = E2; // Highest Energy in keV
      curve.para.DeltaE = DeltaE; // Delta(E) in keV
    }


    // initialize background
    for ( unsigned int p(0); p < numbands; p++ ) background[p] = 0.0;

    /*************************/
    /* OPENING THE DATA FILE */
    /*************************/
   
    if ( datafile_is_set ) {
      std::cout << "setting data file" << std::endl;	
      std::ifstream data; //(data_file);      // the data input stream
      std::cout << "opening data file" << std::endl;
      data.open( data_file );  // opening the file with observational data
      char line[265]; // line of the data file being read in
      unsigned int numLines(0), i(0);
      if ( data.fail() || data.bad() || !data ) {
      std::cout << "fail in loading data" << std::endl;
	throw( Exception("Couldn't open data file."));
	return -1;
      }
      // need to do the following to use arrays of pointers (as defined in Struct)
      for (unsigned int y(0); y < numbands; y++) {
      	std::cout << "setting pointers" << std::endl;
	obsdata.t = new double[numbins];
	obsdata.f[y] = new double[numbins];
	obsdata.err[y] = new double[numbins];
      }
 
      /****************************************/
      /* READING IN FLUXES FROM THE DATA FILE */
      /****************************************/
    	
      while ( data.getline(line,265) ) {
	i = numLines;
	double get_t;
	double get_f1;
	double get_err1;
	double get_f2;
	double get_err2;
			
	sscanf( line, "%lf %lf %lf %lf %lf", &get_t, &get_f1, &get_err1, &get_f2, &get_err2 );
	//std::cout << "allocate to matrix" << std::endl;	
	obsdata.t[i] = get_t;	
	obsdata.f[0][i] = get_f1;
	obsdata.err[0][i] = get_err1;
	obsdata.f[1][i] = get_f2;
	obsdata.err[1][i] = get_err2;
	// PG1808 data has 2 sigma error bars -- divide by 2!
	//obsdata.err[1][i] *= 2.0; // to account for systematic errors
	//obsdata.err[2][i] *= 2.0; // to account for systematic errors
	numLines++;
      } 

      //if ( numLines != numbins ) {
	//throw (Exception( "Numbins indicated in command-line not equal to numbins in data file."));
	std::cout << "Warning! Numbins from command-line not equal to numbins in data file." << std::endl;
	std::cout << "Command-line numbins = " << numbins <<", data file numbins = " << numLines << std::endl;
	std::cout << "\t! Setting numbins = numbins from data file." << std::endl;
	numbins = numLines;
	curve.numbins = numLines;
	//return -1;
      //}
       
      // Read in data file to structure "obsdata" (observed data)
      // f[1][i] flux in low energy band
      // f[2][i] flux in high energy band
      obsdata.numbins = numbins;

      data.close();
      //ts = obsdata.t[0]; // Don't do this if you want manually setting ts to do anything!!
      //obsdata.shift = obsdata.t[0];
      obsdata.shift = ts;

      //std::cout << "Finished reading data from " << data_file << ". " << std::endl;
    } // Finished reading in the data file
		
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

      //std::cout << " r_pole = " <<  model->R_at_costheta(1.0) << std::endl;

       
    }
    else if ( NS_model == 2 ) { // Oblate Colour-Flavour Locked Quark Star model
        // Alternative model for quark stars (not very different)
        model = new PolyOblModelCFLQS( req,
				     PolyOblModelBase::zetaparam(mass,rspot),
				     PolyOblModelBase::epsparam(omega, mass, rspot) );
        //printf("Oblate Colour-Flavour Locked Quark Star. ");
    }
    else if ( NS_model == 3 ) { // Standard spherical model
        // Spherical neutron star
      //req = rspot;
      rspot = req;
        model = new SphericalOblModel( rspot );
        //printf("Spherical Model. ");
    }
    else {
        throw(Exception("\nInvalid NS_model parameter. Exiting.\n"));
        return -1;
    }

   
    /****************************/
    /* Initialize time and flux */
    /****************************/
	
    for ( unsigned int i(0); i < numbins; i++ ) {
        curve.t[i] = i / (1.0 * numbins);// + ts;  // defining the time used in the lightcurves
        for ( unsigned int p(0); p < numbands; p++ ) {
            Flux[p][i] = 0.0;                   // initializing flux to 0
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
    else //no
    	pieces=1;

    // If spot is in 2 pieces, p=0 is the crescent; p=1 is the symmetric part over the pole



    for (unsigned int p(0);p<pieces;p++){

      curve = SpotShape(pieces,p,numtheta,theta_1,rho, &curve);

      double deltatheta(0.0);

    // Looping through the mesh of the spot
      	for (unsigned int k(0); k < numtheta; k++) { // Loop through the circles
	 
	  deltatheta = curve.para.dtheta[k];

	  double thetak = curve.para.theta_k[k];

	  double phi_edge = curve.para.phi_k[k];

	  dphi = 2.0*Units::PI/(numbins*1.0);

	  mu_1 = cos(thetak);
	  if (fabs(mu_1) < DBL_EPSILON) mu_1 = 0.0;

	  if ( mu_1 < 0.0){
	    std::cout << "Southern Hemisphere! mu=" << mu_1 << std::endl;
			  //  mu_1 = fabs(mu_1);
			  //thetak = Units::PI - thetak;
			  //curve.para.incl = Units::PI - incl_1;
	  }

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
	  curve.para.theta = thetak;

	  if (numtheta==1){  //For a spot with only one theta bin (used for small spot)
	    numphi=1;
	    phi_edge=0.0;
	    dphi=0.0;
	    phishift = 0.0;
	    curve.para.dS = 2.0*Units::PI * pow(rspot,2) * (1.0 - cos(rho));
	  }
       
	  if ( NS_model != 3 ) curve.para.dS /= curve.para.cosgamma;

	  for ( unsigned int j(0); j < numphi; j++ ) {// looping through the phi divisions
	    
	    curve.para.phi_0 =  -phi_edge + (j+0.5)*dphi;			

	    //Heart of spot, calculate curve for the first phi bin - otherwise just shift
	    if ( j==0){ 
	      curve = ComputeAngles(&curve, defltoa); 
	      curve = ComputeCurve(&curve);
	    }
	
	    if ( curve.para.temperature == 0.0 ) {// Flux is zero for parts with zero temperature
	      for ( unsigned int i(0); i < numbins; i++ ) {
		for ( unsigned int p(0); p < numbands; p++ ) curve.f[p][i] = 0.0;
	      }
	    }

	    // Add curves, load into Flux array
	    for ( unsigned int i(0); i < numbins; i++ ) {
	      int q(i+j);
	      if (q>=numbins) q+=-numbins;
	      for ( unsigned int p(0); p < numbands; p++ ) {
		Flux[p][i] += curve.f[p][q];
	      }
	    } // ending Add curves
	  } // end for-j-loop

	  // Add in the missing bit.      
	  if (phishift != 0.0){ // Add light from last bin, which requires shifting
	    for ( unsigned int i(0); i < numbins; i++ ) {
	      int q(i+numphi-1);
	      if (q>=numbins) q+=-numbins;
	      for ( unsigned int pp(0); pp < numbands; pp++ ) {
		Temp[pp][i] = curve.f[pp][q];
	      }
	    }
	    for (unsigned int pp(0); pp < numbands; pp++ )
	      for ( unsigned int i(0); i < numbins; i++ ) curve.f[pp][i] = Temp[pp][i];	  
	    
	    curve = ShiftCurve(&curve,phishift);
       
	    for( unsigned int pp(0); pp < numbands; pp++ )
	      for ( unsigned int i(0); i < numbins; i++ ) Flux[pp][i] += curve.f[pp][i]*phishift/dphi      ;
	  } //end of last bin  
	  delete defltoa;
      	} // closing for loop through theta divisions
    } // End Standard Case of first spot

    /**********************************************************/
    /* SECOND HOT SPOT -- Mirroring of first hot spot         */
    /**********************************************************/

    //Setting new parameters for second spot
    if ( two_spots ) {
    	incl_2 = Units::PI - incl_1 + d_incl_2; // keeping theta the same, but changing inclination
    	curve.para.incl = incl_2;
    	theta_2 = theta_1 + d_theta_2;
    	curve.para.theta = theta_2;  // keeping theta the same, but changing inclination
    	cosgamma = model->cos_gamma(mu_2);
    	curve.para.cosgamma = cosgamma;
    		    		
    	if ( T_mesh_in ) {
      		std::cout << "WARNING: code can't handle a spot asymmetric over the pole with a temperature mesh." << std::endl;
      		spot_temperature = 2;
    	}
    	curve.para.temperature = spot_temperature;

    	int pieces;
    	//Does the spot go over the pole?
    	if ( rho > theta_2){ // yes
      		pieces=2;
   		}
    	else //no
 			pieces=1;

    	for (unsigned int p(0);p<pieces;p++){

    		double deltatheta = 2.0*rho/numtheta;
    		if (pieces==2){
				if (p==0) deltatheta = (rho-theta_2)/numtheta;
				else deltatheta = (2.0*theta_2)/numtheta;
    		}      
    		//std::cout << "numtheta=" << numtheta << " theta = " << theta_2 << " delta(theta) = " << deltatheta << std::endl;
     
      		// Looping through the mesh of the spot
    		for (unsigned int k(0); k < numtheta; k++) { // Loop through the circles
				//std::cout << "k=" << k << std::endl;

				double thetak = theta_2 - rho + (k+0.5)*deltatheta; 
				double phi_edge, phij;

				if (pieces==2){
	  				if (p==0){
	    				thetak = (k+0.5)*deltatheta;
	    				phi_edge = Units::PI;
	  				}
	  				else thetak = rho - theta_2 + (k+0.5)*deltatheta;  	
				}

				//std::cout << "k=" << k << " thetak = " << thetak << std::endl;
				dphi = 2.0*Units::PI/(numbins*1.0);

				// What is the value of radius at this angle?
				// For spherical star, rspot = req;
				rspot = req; // Spherical star
				curve.para.radius = rspot; // load rspot into structure
				curve.para.mass_over_r = mass_over_req * req/rspot;

				OblDeflectionTOA* defltoa = new OblDeflectionTOA(model, mass, curve.para.mass_over_r , rspot); 
				//std::cout << "Now Compute the bending angles by entering Bend" << std::endl;
				curve = Bend(&curve,defltoa);

				if ( (pieces==2 && p==1) || (pieces==1)){ //crescent-shaped 2nd piece, or the one circular piece if spot doesn't go over pole
	  				double cos_phi_edge = (cos(rho) - cos(theta_2)*cos(thetak))/(sin(theta_2)*sin(thetak)) * -1;
	  				//std::cout << cos_phi_edge << std::endl;
	  				if (  cos_phi_edge > 1.0 || cos_phi_edge < -1.0 ) cos_phi_edge = 1.0;
					if ( fabs( sin(theta_2) * sin(thetak) ) > 0.0) { // checking for a divide by 0
	    				phi_edge = acos( cos_phi_edge );   // value of phi (a.k.a. azimuth projected onto equatorial plane) at the edge of the circular spot at some latitude thetak
	    				//std::cout << phi_edge << std::endl;
	  				}
	  				else {  // trying to divide by zero
	    				throw( Exception(" Tried to divide by zero in calculation of phi_edge for spot 2. Likely, thetak = 0. Exiting.") );
	    				return -1;
	  				}
				} 
				numphi = 2.0*(Units::PI-phi_edge)/dphi;
				phishift = 2.0*(Units::PI-phi_edge)- numphi*dphi;

				curve.para.dS = pow(rspot,2) * sin(thetak) * deltatheta * dphi;
				curve.para.theta = thetak;

				if (numtheta==1){
	  				numphi=1;
	  				phi_edge=-1*Units::PI;
	  				dphi=0.0;
	  				phishift = 0.0;
	  				curve.para.dS = 2.0*Units::PI * pow(rspot,2) * (1.0 - cos(rho));
				}
       

    			for ( unsigned int j(0); j < numphi; j++ ) {   // looping through the phi divisions
					phij = phi_edge + (j+0.5)*dphi;
					curve.para.phi_0 = phij;

					if ( NS_model == 1 || NS_model == 2 ) curve.para.dS /= curve.para.cosgamma;
					if ( j==0){ // Only do computation if its the first phi bin - otherwise just shift
	  					curve = ComputeAngles(&curve, defltoa); 
	  					curve = ComputeCurve(&curve);
					}
	
					if ( curve.para.temperature == 0.0 ) {
	  					for ( unsigned int i(0); i < numbins; i++ ) {
	    					for ( unsigned int p(0); p < numbands; p++ )  curve.f[p][i] = 0.0;
	  					}
					}

					// Add curves, load into Flux array
					for ( unsigned int i(0); i < numbins; i++ ) {
	 					int q(i+j);
	  					if (q>=numbins) q+=-numbins;
	  					for ( unsigned int p(0); p < numbands; p++ ) Flux[p][i] += curve.f[p][q];  
					} // ending Add curves
    			} // end for-j-loop

      			// Add in the missing bit.      
      			if (phishift != 0.0){ // Add light from last bin, which requires shifting
      				for ( unsigned int i(0); i < numbins; i++ ) {
	  					int q(i+numphi-1);
	  					if (q>=numbins) q+=-numbins;
	  					for ( unsigned int p(0); p < numbands; p++ ) Temp[p][i] = curve.f[p][q];	  
					}

					for (unsigned int p(0); p < numbands; p++ ){
	  					for ( unsigned int i(0); i < numbins; i++ ) curve.f[p][i] = Temp[p][i];	  
   	     			}
					curve = ShiftCurve(&curve,phishift);
       
					for( unsigned int p(0); p < numbands; p++ ){
	  					for ( unsigned int i(0); i < numbins; i++ ) Flux[p][i] += curve.f[p][i]*phishift/dphi      ;
					}
      			} // closing the missing bit.      
    		} // closing for loop through theta divisions
    	} // closing for loop for calculating each piece, either the symmetric or crescent
	} // closing second spot
            
    // You need to be so super sure that ignore_time_delays is set equal to false.
    // It took almost a month to figure out that that was the reason it was messing up.
     
    
    /*******************************/
    /* NORMALIZING THE FLUXES TO 1 */
    /*******************************/
	    
    // Normalizing the flux to 1 in low energy band. 
    if ( normalize_flux ) {
      normcurve = Normalize1( Flux, numbins );

      // Add background to normalized flux
      for ( unsigned int i(0); i < numbins; i++ ) {
	for ( unsigned int p(0); p < numbands; p++ ) {
	  Flux[p][i] = normcurve.f[p][i] + background[p];
	}
      } 
      
      // Renormalize to 1.0
      normcurve = Normalize1( Flux, numbins );

      // fun little way around declaring the Chi.cpp method Normalize as returning a matrix! so that we can still print out Flux
      for ( unsigned int i(0); i < numbins; i++ ) {
	for ( unsigned int p(0); p < numbands; p++ ) {
	  curve.f[p][i] = normcurve.f[p][i];
	}
      }  
    } // Finished Normalizing
    else{ // Curves are not normalized
      for ( unsigned int i(0); i < numbins; i++ ) {
	for ( unsigned int p(0); p < numbands; p++ ) {
	  curve.f[p][i] = Flux[p][i];
	}
      }  
    }
     
    /************************************************************/
    /* If data file is set, calculate chi^2 fit with simulation */
    /************************************************************/
	
    if ( datafile_is_set ) {
    	std::cout << "calculating chi squared" << std::endl;
    	chisquared = ChiSquare ( &obsdata, &curve );
    }
    
    std::cout << "Spot: m = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN 
	      << " Msun, r = " << Units::nounits_to_cgs(req, Units::LENGTH )*1.0e-5 
	      << " km, f = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) 
	      << " Hz, i = " << incl_1 * 180.0 / Units::PI 
	      << ", e = " << theta_1 * 180.0 / Units::PI 
	      << ", X^2 = " << chisquared 
	      << std::endl;    


    /*******************************/
    /* CALCULATING PULSE FRACTIONS */
    /*******************************/
    
    double avgPulseFraction(0.0); // average of pulse fractions across all energy bands
    double sumMaxFlux(0.0), sumMinFlux(0.0);
    double overallPulseFraction(0.0);
    for ( unsigned int j(0); j < numbands; j++ ) {
      sumMaxFlux += curve.maxFlux[j];
      sumMinFlux += curve.minFlux[j];
      avgPulseFraction += curve.pulseFraction[j];
    }
    overallPulseFraction = (sumMaxFlux - sumMinFlux) / (sumMaxFlux + sumMinFlux);
    avgPulseFraction /= (NCURVES);
	
    double U(0.0), Q(0.0), A(0.0);   // as defined in PG19 and PG20
    U = (1 - 2 * mass / rspot) * sin(incl_1) * sin(theta_1); // PG20
    Q = 2 * mass / rspot + (1 - 2 * mass / rspot) * cos(incl_1) * cos(theta_1); // PG20
    A = U / Q; // PG19
    B = rspot * sin(incl_1) * sin(theta_1) * omega / ( sqrt( 1 - ((2*mass) / rspot) ) ); // param_degen/equations.pdf 2    

    /********************************************/
    /* WRITING THE SIMULATION TO AN OUTPUT FILE */
    /********************************************/ 
    	
    out.open(out_file, std::ios_base::trunc);
    if ( out.bad() || out.fail() ) {
        std::cerr << "Couldn't open output file. Exiting." << std::endl;
        return -1;
    }

    if ( rho == 0.0 ) rho = Units::PI/180.0;

    out << "% Photon Flux for hotspot on a NS. \n"
        << "% R_sp = " << Units::nounits_to_cgs(rspot, Units::LENGTH )*1.0e-5 << " km; "
        << "% R_eq = " << Units::nounits_to_cgs(req, Units::LENGTH )*1.0e-5 << " km; "
        << "% M = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN << " Msun; "
        << "% Spin = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) << " Hz \n"
        << "% Gravitational redshift 1+z = " << 1.0/sqrt(1.0 - 2*mass/rspot) << "\n"
        << "% Inclination angle = " << incl_1 * 180.0/Units::PI << " degrees \n"
        << "% Emission angle = " << theta_1 * 180.0/Units::PI << " degrees \n"
        << "% Angular radius of spot = " << rho << " radians \n"
        << "% Number of spot bins = " << numtheta << " \n"
        << "% Distance to NS = " << Units::nounits_to_cgs(distance, Units::LENGTH)*.01 << " m \n" 
        << "% Phase shift or time lag = " << ts << " \n"
      //        << "# Pulse fractions: bolo = " << curve.pulseFraction[0] <<", mono = " << curve.pulseFraction[1] << ", \n"
      //<< "#                  low E band = " << curve.pulseFraction[2] << ", high E band = " << curve.pulseFraction[3] << " \n"
      //<< "# Rise/Fall Asymm: bolo = " << curve.asym[0] <<", mono = " << curve.asym[1] << ", \n"
      //<< "#                  low E band = " << curve.asym[2] << ", high E band = " << curve.asym[3] << " \n"
      //<< "# B = " << B << " \n"
      //<< "# B/(1+z) = " << B * sqrt(1.0 - 2*mass/rspot) << "\n"
      //<< "# U = " << U << " \n"
      //<< "# Q = " << Q << " \n"
      //<< "# A = U/Q = " << A << " \n"
      //<< "# %(Amp-Bolo) = " << (A - curve.pulseFraction[0])/A * 100.0 << " \n"
      //<< "# %(Bolo-low) = " << (curve.pulseFraction[0]-curve.pulseFraction[2])/curve.pulseFraction[0] * 100.0 << " \n"
        << std::endl;
        
    if (datafile_is_set)
    	out << "% Data file " << data_file << ", chisquared = " << chisquared << std::endl;

    if ( T_mesh_in )
    	out << "% Temperature mesh input: "<< T_mesh_file << std::endl;
    else
    	out << "% Spot temperature, (star's frame) kT = " << spot_temperature << " keV " << std::endl;
	if ( NS_model == 1)
    	out << "% Oblate NS model " << std::endl;
    else if (NS_model == 3)
    	out << "% Spherical NS model " << std::endl;
    if ( beaming_model == 0 )
        out << "% Isotropic Emission " << std::endl;
    else
        out << "% Limb Darkening for Gray Atmosphere (Hopf Function) " << std::endl;
    if ( normalize_flux )
    	out << "% Flux normalized to 1 " << std::endl;
    else
    	out << "% Flux not normalized " << std::endl;
    
    /***************************************************/
    /* WRITING COLUMN HEADINGS AND DATA TO OUTPUT FILE */
    /***************************************************/
  
    out << "%\n"
    	<< "% Column 1: phase bins (0 to 1)\n";

    if (curve.flags.spectral_model==0){
      out << "% Column 2: Monochromatic Number flux (photons/(cm^2 s keV) measured at energy (at infinity) of " << curve.para.E0 << " keV\n";
      out << "%"
	  << std::endl;
      for ( unsigned int i(0); i < numbins; i++ ) {
        out << curve.t[i]<< "\t" ;
	for ( unsigned int p(0); p < numbands; p++ ) { 
            out << curve.f[p][i] << "\t" ;
        }
	out << i;
        out << std::endl;
      }
    }

    if (curve.flags.spectral_model==1){
      out << "% Column 2: Photon Energy (keV) in Observer's frame" << std::endl;
      out << "% Column 3: Number flux (photons/(cm^2 s) " << std::endl;
      out << "% " << std::endl;

      for ( unsigned int p(0); p < numbands; p++ ) { 
	for ( unsigned int i(0); i < numbins; i++ ) {
	  out << curve.t[i]<< "\t" ;
	  out << curve.para.E0 + p*curve.para.DeltaE << "\t";
	  out << curve.f[p][i] << "\t" << std::endl;
	  //out << i;
	  //out << std::endl;
	}

      }
    }


    if (curve.flags.spectral_model==2){
    	double E_diff;
		E_diff = (E_band_upper_1 - E_band_lower_1)/numbands;
    	for (unsigned int k(0); k < numbands; k++){
      		out << "% Column" << k+2 << ": Integrated Number flux (photons/(cm^2 s) measured between energy (at infinity) of " << curve.para.E_band_lower_1+k*E_diff << " keV and " << curve.para.E_band_lower_1+(k+1)*E_diff << " keV\n";    		
    	}
      	out << "%" << std::endl;
      	for ( unsigned int i(0); i < numbins; i++ ) {
        	out << curve.t[i]<< "\t" ;
			for ( unsigned int p(0); p < numbands; p++ ) { 
            	out << curve.f[p][i] << "\t" ;
        	}
			out << i;
        	out << std::endl;
      	}
    }


    if (curve.flags.spectral_model==3){
    	double E_diff;
		E_diff = (E_band_upper_1 - E_band_lower_1)/numbands;
    	for (unsigned int k(0); k < numbands; k++){
      		out << "% Column " << k+2 << ": Integrated Number flux (photons/(cm^2 s) measured between energy (at infinity) of " << curve.para.E_band_lower_1+k*E_diff << " keV and " << curve.para.E_band_lower_1+(k+1)*E_diff << " keV\n";    		
    	}
      	out << "%" << std::endl;
      	for ( unsigned int i(0); i < numbins; i++ ) {
        	out << curve.t[i]<< "\t" ;
			for ( unsigned int p(0); p < numbands; p++ ) { 
            	out << curve.f[p][i] << "\t" ;
            	out << 1e-6 << "\t" ; // Fake errors.  
        	}
			//out << i;
        	out << std::endl;
      	}
    }



      //<< "# Column 4: Number flux (photons/(cm^2 s)) in the energy band " << E_band_lower_1 << " keV to " << E_band_upper_1 << " keV \n"
      //<< "# Column 5: Number flux (photons/(cm^2 s)) in the energy band " << E_band_lower_2 << " keV to " << E_band_upper_2 << " keV \n"

    
  


    out.close();
    
    //delete defltoa;
    delete model;
    return 0;
} 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
