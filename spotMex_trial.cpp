/***************************************************************************************/
/*                                   SpotMex.cpp

    This code produces a pulse profile once a set of parameter describing the star, 
    spectrum, and hot spot have been inputed.
    
    This is the matlab executable version of Spot.cpp, by Abigail Stevens (2012/2013)
    
*/
/***************************************************************************************/

// INCLUDE ALL THE THINGS! 
// If you do not get this reference, see 
//       http://hyperboleandahalf.blogspot.ca/2010/06/this-is-why-ill-never-be-adult.html
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <exception>
#include <vector>
#include <string>
#include <string.h>
#include <stdio.h>
#include "OblDeflectionTOA.h"
#include "Chi.h"
#include "PolyOblModelNHQS.h"
#include "PolyOblModelCFLQS.h"
#include "SphericalOblModel.h"
#include "OblModelBase.h"
#include "Units.h"
#include "Exception.h"
#include "Struct.h"
#include "time.h"
#include "io64.h"
#include "mex.h"
#include <sstream>
#include "Atmo.h"    
#include <string.h>

// MAIN

void mexFunction ( int numOutputs, mxArray *theOutput[], int numInputs, const mxArray *theInput[]) {
//theOutput[] is an array containing Null pointers.  We need to allocate the memory when returning the results
//numOutputs is the expected number of output pointers, i.e., the number of variables to be returned (unlike C, MatLab can have more than one)
//theInput[] is an array containing pointers to data of mxArray type
//numInputs is the number of input pointers, i.e., the number of variables passed to the function

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
    B,                          // from param_degen/equations.pdf 2
    *curveOut,
    *chiOut;

  double E0, E1, E2, DeltaE, rot_par;
    
  unsigned int NS_model(1),       // Specifies oblateness (option 3 is spherical)
    spectral_model(0),    // Spectral model choice (initialized to blackbody)
    beaming_model(0),     // Beaming model choice (initialized to isotropic)
    numbins(MAX_NUMBINS), // Number of time or phase bins for one spin period; Also the number of flux data points
    numphi(1),            // Number of azimuthal (projected) angular bins per spot
    numtheta(1),          // Number of latitudinal angular bins per spot
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


    // Setting up the output parameters
    int dimSize[2];
    dimSize[0] = 1;
    dimSize[1] = 1;

    // output needs to be an array so set it up as a [1,1] array holding one value
    theOutput[0] = mxCreateNumericArray(2, dimSize, mxDOUBLE_CLASS, mxREAL);

    chiOut = mxGetPr(theOutput[0]);


	/********************************************************/
    /* READ INFORMATION PASSED BY MEX FUNCTION              */
    /* MAKE SURE CORRECT NUMBER OF PARAMETERS ARE PASSED IN */
    /********************************************************/
   
	if ( numInputs != 23 ) { // need to go through and figure out how many inputs we should have
		//mexPrintf("numInputs = %d, expecting 23", numInputs);
		mexErrMsgTxt("Wrong number of inputs.");
	}
    printf("right number of inputs %\n");
	mass = mxGetScalar(theInput[0]); // double
	req = mxGetScalar(theInput[1]); // radius at the equator; double
	omega = mxGetScalar(theInput[2]); // spin frequency; double	
	incl_1 = mxGetScalar(theInput[3]); // inclination angle; double
	theta_1 = mxGetScalar(theInput[4]); // emission angle; double
	ts = mxGetScalar(theInput[5]); // phase shift/time shift; double
	numbins = mxGetScalar(theInput[6]); // int
	NS_model = mxGetScalar(theInput[7]); // int
	rho = mxGetScalar(theInput[8]); // spot size in radian; double 
    spot_temperature = mxGetScalar(theInput[9]); // keV; double
	distance = mxGetScalar(theInput[10]); // in meters; double
	numtheta = mxGetScalar(theInput[11]); // int
	spectral_model = mxGetScalar(theInput[12]); // 3 for new integrated atmosphere scheme; int
	numbands = mxGetScalar(theInput[13]); // int
	E_band_lower_1 = mxGetScalar(theInput[14]); // keV; double 
	E_band_upper_1 = mxGetScalar(theInput[15]); // keV; double
	beaming_model = mxGetScalar(theInput[16]); // 3 for hydrogen, 4 for helium; int			
	int spots_2 = mxGetScalar(theInput[17]); // int
	if (spots_2 == 1)
		two_spots = false;
	if (spots_2 == 2) 
		two_spots = true;
    std::cout << two_spots << std::endl;
	//E0 = mxGetScalar(theInput[18]); // monochromatic band energy; double
	obsdata.t = mxGetPr(theInput[18]); // array of double
	// To be changed to a for loop that depends on number of bands
	obsdata.f[0] = mxGetPr(theInput[19]); // array of double
    obsdata.f[1] = mxGetPr(theInput[20]); // array of double
    obsdata.err[0] = mxGetPr(theInput[21]); // array of double
    obsdata.err[1] = mxGetPr(theInput[22]); // array of double


    

    /***********************************************/
    /* CHECKING THAT THE NECESSARY VALUES WERE SET */
    /***********************************************/
    /*
    if( !( incl_is_set && theta_is_set
	    && mass_is_set && rspot_is_set
	    && omega_is_set && model_is_set ) ) {
        throw( Exception(" Not all required parameters were specified. Exiting.\n") );
        return -1;
    }
    */
   
    /*******************************************************************/
    /* READING TEMPERATURE VALUES FROM AN INPUT FILE, NOT USED FOR NOW */
    /*******************************************************************/
	/*
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
    */
    /******************************************/
    /* SENSIBILITY CHECKS ON INPUT PARAMETERS */
    /******************************************/
    /*
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
    
 
 	*/
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

    //ts = obsdata.t[0]; // Don't do this if you want manually setting ts to do anything!!
    //obsdata.shift = obsdata.t[0];
    obsdata.shift = ts;
    obsdata.numbins = numbins;

    /*************************/
    /* OPENING THE DATA FILE */
    /*************************/
    /*************************/
    /* ALREADY READ IN!!!    */
    /*************************/

    /*
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
    */
      /****************************************/
      /* READING IN FLUXES FROM THE DATA FILE */
      /****************************************/
    /*   
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
    //std::cout << "Warning! Numbins from command-line not equal to numbins in data file." << std::endl;
    std::cout << "Command-line numbins = " << numbins <<", data file numbins = " << numLines << std::endl;
    //std::cout << "\t! Setting numbins = numbins from data file." << std::endl;
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
    */    
		
    /***************************/
    /* START SETTING THINGS UP */
    /***************************/ 

    // Change to computation of rspot!!!!

    // Calculate the Equatorial Radius of the star.
    //req = calcreq( omega, mass, theta_1, rspot );  // implementation of MLCB11

    //rspot = req;

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
    	/*
        throw(Exception("\nInvalid NS_model parameter. Exiting.\n"));
        return -1;
        */
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

    for (unsigned int p(0);p<pieces;p++){

    	double deltatheta = 2.0*rho/numtheta;

      	if (pieces==2){
			if (p==0) deltatheta = (rho-theta_1)/numtheta; //symmetric over pole
			else deltatheta = (2.0*theta_1)/numtheta; //crescent shaped
      	}     
      	//std::cout << "numtheta=" << numtheta << " theta = " << theta_1 << " delta(theta) = " << deltatheta << std::endl;
     
    // Looping through the mesh of the spot
      	for (unsigned int k(0); k < numtheta; k++) { // Loop through the circles
	  //std::cout << "k=" << k << std::endl;

			double thetak = theta_1 - rho + (k+0.5)*deltatheta; 
			double phi_edge, phij;

			if (pieces==2){
	  			if (p==0){
	    			thetak = (k+0.5)*deltatheta;
	    			phi_edge = Units::PI;
	  			}
	  			else {
	    			thetak = rho - theta_1 + (k+0.5)*deltatheta;
	  			}
			}


			dphi = 2.0*Units::PI/(numbins*1.0);

			// What is the value of radius at this angle?
			// For spherical star, rspot = req;

			mu_1 = cos(thetak);
			if (fabs(mu_1) < DBL_EPSILON) mu_1 = 0.0;

			//std::cout << "k=" << k << " thetak = " << thetak << " mu=" << mu_1 << std::endl;

			if ( mu_1 < 0.0){
			  std::cout << "Southern Hemisphere!" << std::endl;
			  mu_1 = fabs(mu_1);
			  thetak = Units::PI - thetak;
			  curve.para.incl = Units::PI - incl_1;
			}

			if (NS_model != 3)
			  rspot = model->R_at_costheta(mu_1);
			else
			  rspot = req;

			//std::cout << "Spot: req = " << req 
			//	  <<" rspot = " << rspot << std::endl;

			// Values we need in some of the formulas.
			cosgamma = model->cos_gamma(mu_1);   // model is pointing to the function cos_gamma
			curve.para.cosgamma = cosgamma;

			curve.para.radius = rspot; // load rspot into structure
			curve.para.mass_over_r = mass_over_req * req/rspot;

			OblDeflectionTOA* defltoa = new OblDeflectionTOA(model, mass, curve.para.mass_over_r , rspot); 
			curve = Bend(&curve,defltoa);

			if ( (pieces==2 && p==1) || (pieces==1)){ //crescent-shaped 2nd piece, or the one circular piece if spot doesn't go over pole
	  			double cos_phi_edge = (cos(rho) - cos(theta_1)*cos(thetak))/(sin(theta_1)*sin(thetak));
	  			if (  cos_phi_edge > 1.0 || cos_phi_edge < -1.0 ) cos_phi_edge = 1.0;
				if ( fabs( sin(theta_1) * sin(thetak) ) > 0.0) { // checking for a divide by 0
	    			phi_edge = acos( cos_phi_edge );   // value of phi (a.k.a. azimuth projected onto equatorial plane) at the edge of the circular spot at some latitude thetak
	    			//std::cout << phi_edge << std::endl;
	  			}
	  			else {  // trying to divide by zero
	  				/*
	    			throw( Exception(" Tried to divide by zero in calculation of phi_edge for spot 1. Likely, thetak = 0. Exiting.") );
	  	  			return -1;
	  	  			*/
	  			}
			}
      
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
       
			if ( NS_model == 1 || NS_model == 2 ) curve.para.dS /= curve.para.cosgamma;

			for ( unsigned int j(0); j < numphi; j++ ) {// looping through the phi divisions
			  	phij = -phi_edge + (j+0.5)*dphi;
				curve.para.phi_0 = phij;

			

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
	  				for ( unsigned int p(0); p < numbands; p++ ) {
	    				Temp[p][i] = curve.f[p][q];
	  				}
				}
				for (unsigned int p(0); p < numbands; p++ )
	  				for ( unsigned int i(0); i < numbins; i++ ) curve.f[p][i] = Temp[p][i];	  
   	     
				curve = ShiftCurve(&curve,phishift);
       
				for( unsigned int p(0); p < numbands; p++ )
	  				for ( unsigned int i(0); i < numbins; i++ ) Flux[p][i] += curve.f[p][i]*phishift/dphi      ;
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
	  					/*
	    				throw( Exception(" Tried to divide by zero in calculation of phi_edge for spot 2. Likely, thetak = 0. Exiting.") );
	    				return -1;
	    				*/
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
     
    /***************************************************/
    /* Data file should already be available           */
    /***************************************************/
	
    chisquared = ChiSquare ( &obsdata, &curve );
    
    
    std::cout << "Spot: m = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN 
	      << " Msun, r = " << Units::nounits_to_cgs(req, Units::LENGTH )*1.0e-5 
	      << " km, f = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) 
	      << " Hz, i = " << incl_1 * 180.0 / Units::PI 
	      << ", e = " << theta_1 * 180.0 / Units::PI 
	      << ", X^2 = " << chisquared 
	      << std::endl;    
 

	if (std::isnan(chisquared)) {
		mexPrintf("Chisquared is nan! Setting it to 1,000,000.\n");
		chisquared = 10000000.0;
	}
	else {
		mexPrintf("X^2 = %f\n\n", chisquared);
	}
	
	if (ts > 1.0) ts -= 1.0; // to make sure ts is in [0:1], for printing out
    //std::cout << "Pre printing to everything." << std::endl;

	/*****************************************************************/
    /* DUMPING DATA INTO MATLAB OUTPUT                               */
    /* For filewriting, need to use fopen, fprintf, and fclose (etc) */
    /*****************************************************************/
    //std::cout << "Pre output." << std::endl;
	
    chiOut[0] = chisquared; // saved this to theOutput[0] at top just after declarations
    //std::cout << "Output 1." << std::endl;
    dimSize[1] = 3; // number of columns (1: time, 2: flux in 1st energy band, 3: flux in 2nd energy band)
    //std::cout << "Output 2." << std::endl;

    dimSize[0] = (int)numbins; // number of rows ( = numbins)
    //std::cout << "Output 3." << std::endl;

    theOutput[1] = mxCreateNumericArray(2, dimSize, mxDOUBLE_CLASS, mxREAL); // formatting/setting up the output
    //std::cout << "Output 4." << std::endl;

    curveOut = mxGetPr(theOutput[1]); // matlab and mex love pointers
	

 
    delete model;
}
