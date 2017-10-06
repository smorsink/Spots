/***************************************************************************************/
/*                                   SpotMex_new.cpp

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
#include "instru.h"    
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
#include "interp.h"
#include "nrutil.h"
#include <unistd.h>
#include <string.h>

// MAIN
void mexFunction ( int numOutputs, mxArray *theOutput[], int numInputs, const mxArray *theInput[]) {

std::cout << "Hello World!" << std::endl;

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
    background[NCURVES],	        // One background value for all bands.
    dsbackground(0.0),			// Diffuse sky background normalization
    agnbackground(0.0),			// AGN background normalization
    plbackground(0.0),			// Phase-dependent power law background normalization
    T_mesh[30][30],             // Temperature mesh over the spot; same mesh as theta and phi bins, assuming square mesh
    chisquared(1.0),            // The chi^2 of the data; only used if a data file of fluxes is inputed
    distance(3.0857e20),        // Distance from earth to the NS, in meters; default is 10kpc
    obstime(1.0),               // Length of observation (in seconds)
    phase_2(0.5),				// Phase of second spot, 0 < phase_2 < 1
    nh(1.0);					// real nh = nh*4e19
   
  double SurfaceArea(0.0);

  double E0, L1, L2, DeltaE, rot_par;
    
  unsigned int NS_model(1),       // Specifies oblateness (option 3 is spherical)
    spectral_model(0),    // Spectral model choice (initialized to blackbody)
    beaming_model(0),     // Beaming model choice (initialized to isotropic)
    numbins(MAX_NUMBINS), // Number of time or phase bins for one spin period; Also the number of flux data points
    databins(MAX_NUMBINS),   // Number of phase bins in the data
    numphi(1),            // Number of azimuthal (projected) angular bins per spot
    numtheta(1),          // Number of latitudinal angular bins per spot
    spotshape(0), 		  // Spot shape; 0=standard
    numbands(NCURVES),    // Number of energy bands;
    attenuation(0),       // Attenuation flag, specific to NICER targets with implemented factors
    inst_curve(0);		  // Instrument response flag, 1 = NICER response curve

  unsigned int energybandfactor(1); //We compute NCURVES energy bands, and interpolate to get NCURVESxenergybandfactor final energy bands.

  int NlogTeff, Nlogg, NlogE, Nmu, Npts;

  char out_file[256] = "flux.txt",    // Name of file we send the output to; unused here, done in the shell script
       out_file1[256] = "flux2.txt",    // Name of file we send the output to; unused here, done in the shell script
       bend_file[256] = "No File Name Specified!", 
       //out_dir[80],                   // Directory we could send to; unused here, done in the shell script
       T_mesh_file[100],              // Input file name for a temperature mesh, to make a spot of any shape
       data_file[256],                // Name of input file for reading in data
       //testout_file[256] = "test_output.txt", // Name of test output file; currently does time and two energy bands with error bars
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
    	 normalize_flux(false),      // True if we are normalizing the output flux to 1
    	 //E_band_lower_2_set(false),  // True if the lower bound of the second energy band is set
    	 //E_band_upper_2_set(false),  // True if the upper bound of the second energy band is set
    	 two_spots(false),           // True if we are modelling a NS with two antipodal hot spots
    	 //only_second_spot(false),    // True if we only want to see the flux from the second hot spot (does best with normalize_flux = false)
    	 pd_neg_soln(false),
    	 background_file_is_set(false),
    	 out_file1_is_set(false);	 // True to produce second output file with alternate format
		
  // Create LightCurve data structure
  class LightCurve curve, normcurve;  // variables curve and normalized curve, of type LightCurve
  //class LightCurve curve;
  class LightCurve *flxcurve;
  class DataStruct obsdata;           // observational data as read in from a file


	double *curveOut, *chiOut;

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
   
    //std::cout << "number of inputs is " << numInputs << std::endl;
	mass = mxGetScalar(theInput[0]); // double
	req = mxGetScalar(theInput[1]); // radius at the equator; double
	omega = mxGetScalar(theInput[2]); // spin frequency; double	
	incl_1 = mxGetScalar(theInput[3]); // inclination angle; double
	theta_1 = mxGetScalar(theInput[4]); // emission angle; double
	ts = mxGetScalar(theInput[5]); // phase shift/time shift; double
	databins = mxGetScalar(theInput[6]); // int
    if (databins < MIN_NUMBINS){
        numbins = MIN_NUMBINS;
    }
    else{
        numbins = databins;
    }
	NS_model = mxGetScalar(theInput[7]); // int
	rho = mxGetScalar(theInput[8]); // spot size in radian; double 
    spot_temperature = mxGetScalar(theInput[9]); // keV; double
	distance = mxGetScalar(theInput[10]); // in kpc; double
    distance *= 3.0857e19; // convert distance to meters
	numtheta = mxGetScalar(theInput[11]); // int
	spectral_model = mxGetScalar(theInput[12]); // 3 for new integrated atmosphere scheme; int
	numbands = mxGetScalar(theInput[13]); // int
	E_band_lower_1 = mxGetScalar(theInput[14]); // keV; double 
	E_band_upper_1 = mxGetScalar(theInput[15]); // keV; double
	beaming_model = mxGetScalar(theInput[16]); // 3 for hydrogen, 4 for helium; int			
	int spots_2 = mxGetScalar(theInput[17]); // int
	if (spots_2 == 1){
		two_spots = false;
        //std::cout << "one spot" << std::endl;
	}
    if (spots_2 == 2){ 
		two_spots = true;
        std::cout << "two_spots" << std::endl;
	}
    
    
    std::cout << "Spot: m = " << mass
	      << " Msun, r = " << req
	      << " km, f = " << omega 
	      << " Hz, i = " << incl_1 
	      << ", e = " << theta_1  
	      << ", ts = " << ts
                << ", rho = " << rho
                << ", T = " << spot_temperature
               << ", ObsTime = " << obstime << "s" 
            << ", Distance = " << distance
	      << std::endl;   
    

	//double *datatime = mxGetPr(theInput[18]);
    //std::cout << "Read in the time vector! " << std::endl;
    
    obsdata.t = mxGetPr(theInput[18]); // array of double
	//for (int i=0;i<databins;i++){
	//obsdata.t[i] = datatime[i];
	//std::cout << "i = " << i << " time = " << obsdata.t[i] << std::endl;
	//}



	int bend_file_is = mxGetScalar(theInput[19]);

    if (bend_file_is == 1){
        bend_file_is_set = true;
        curve.flags.bend_file = true;
        curve.defl.mr = mxGetPr(theInput[20]);
        curve.defl.num_mr = 1000;
        double *bend_data_b = mxGetPr(theInput[21]);
        double *bend_data_psi = mxGetPr(theInput[22]);
        double *bend_data_dcosa = mxGetPr(theInput[23]);
        double *bend_data_toa = mxGetPr(theInput[24]);

        curve.defl.psi = dmatrix(0,1001,0,301);
        curve.defl.b = dmatrix(0,1001,0,301);
        curve.defl.dcosa = dmatrix(0,1001,0,301);
        curve.defl.toa = dmatrix(0,1001,0,301);
              
        for (int j = 0; j < 1001; j++){
            for (int i = 0; i < 301; i++){
                curve.defl.b[j][i] = bend_data_b[j*301+i];
                curve.defl.psi[j][i] = bend_data_psi[j*301+i];
                curve.defl.dcosa[j][i] = bend_data_dcosa[j*301+i];
                curve.defl.toa[j][i] = bend_data_toa[j*301+i];
            }
        }
       
    }

    //std::cout << "Successfully read in the Bend file " << std::endl;


    spotshape = mxGetScalar(theInput[25]);
    
    obstime = mxGetScalar(theInput[26]); // in Mega-seconds
    //obstime *= 1e6; // convert to seconds
   
       
       
    inst_curve = mxGetScalar(theInput[27]);
    attenuation = mxGetScalar(theInput[28]);
    double *atmodata_mccinte = mxGetPr(theInput[29]);
    curve.mccinte = dvector(0,1595001);
    for (int i = 0; i < 1595001; i++){
        curve.mccinte[i] = atmodata_mccinte[i];
    }

   double *atmodata_mccangl = mxGetPr(theInput[30]);
    curve.mccangl = dvector(0,51);
    for (int i = 0; i < 51; i++){
        curve.mccangl[i] = atmodata_mccangl[i];
    }
    curve.mcloget = dvector(0,101);
    curve.mcloget = mxGetPr(theInput[31]);
 
    
  
    for (int i = 0; i < numbands-1; i++){
    obsdata.f[i] = mxGetPr(theInput[32+2*i]); // array of double
    //std::cout << "i=" << i << " flux = " << obsdata.f[i][15] << std::endl;
    //obsdata.f[i] = mxGetPr(theInput[31+3*i]); // array of double
    //obsdata.err[i] = mxGetPr(theInput[32+3*i]); // array of double
    background[i] = mxGetScalar(theInput[33+2*i]);
       //std::cout << " bg["<<i<<"]=" << background[i] ;
    }
      // std::cout <<"We Read in the flux!!!" << std::endl;
    /*for( int i = 0; i < databins ; i++){   
    std::cout << "i = " << i  << " Flux[0][i]  = " << obsdata.f[0][i] << " background = " << background[0] << std::endl;
    }*/
    
    /*****************************************************/
    /* UNIT CONVERSIONS -- MAKE EVERYTHING DIMENSIONLESS */
    /*****************************************************/

    mass_over_req = mass/(req) * Units::GMC2;
    //std::cout << "GM/(R_eqc^2) = " << mass_over_req << std::endl;
 
    incl_1 *= (Units::PI / 180.0);  // radians
    //d_incl_2 *= (Units::PI / 180.0);  // radians
    //if ( only_second_spot ) incl_1 = Units::PI - incl_1; // for doing just the 2nd hot spot
    theta_1 *= (Units::PI / 180.0); // radians
    //d_theta_2 *= (Units::PI / 180.0);  // radians
    //theta_2 = theta_1+d_theta_2; // radians
    //rho *= (Units::PI / 180.0);  // rho is input in radians
    mu_1 = cos( theta_1 );
    //mu_2 = mu_1; 
    mass = Units::cgs_to_nounits( mass*Units::MSUN, Units::MASS );
    req = Units::cgs_to_nounits( req*1.0e5, Units::LENGTH );
   
    omega = Units::cgs_to_nounits( 2.0*Units::PI*omega, Units::INVTIME );
    distance = Units::cgs_to_nounits( distance*100, Units::LENGTH );
    rot_par = pow(omega*req,2)/mass_over_req;




    /**********************************/
    /* PASS VALUES INTO THE STRUCTURE */
    /**********************************/  
    
   // std::cout << "Now load it all into the structure!" << std::endl;
    
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

    //std::cout << "ints_curve = " << curve.flags.inst_curve << std::endl;


   // Define the Observer's Spectral Model

   
    /***************************/
    /* START SETTING THINGS UP */
    /***************************/ 

  

    /*********************************************************************************/
    /* Set up model describing the shape of the NS; oblate, funky quark, & spherical */
    /*********************************************************************************/
	
    OblModelBase* model;
    if ( NS_model == 1 ) { // Oblate Neutron Hybrid Quark Star model
        // Default model for oblate neutron star

      //std::cout << " Oblate Neutron Star" << std::endl;
      model = new PolyOblModelNHQS( req,
		   		    mass_over_req,
				    rot_par );

      //std::cout << " r_pole = " <<  model->R_at_costheta(1.0) << std::endl;

       
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
            flxcurve->f[p][i] = 0.0;
	}
    } 

  

    /*********************************/
    /* FIRST HOT SPOT - STANDARD CASE*/
    /*********************************/
    
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

    numbands = NCURVES;


    for (unsigned int p(0);p<pieces;p++){

      	curve = SpotShape(pieces,p,numtheta,theta_1,rho, &curve, model);
      	double deltatheta(0.0);

	// Looping through the mesh of the spot
      	for (unsigned int k(0); k < numtheta; k++) { // Loop through the circles
	  //   		for (unsigned int k(9); k < 10; k++) { // Loop through the circles
	 
	  deltatheta = curve.para.dtheta[k];

	  double thetak = curve.para.theta_k[k];
	  	
	  double phi_edge = curve.para.phi_k[k];

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

	  for ( unsigned int j(0); j < numphi; j++ ) {// looping through the phi divisions
	    
	    curve.para.phi_0 =  -phi_edge + (j+0.5)*dphi;			

	    //Heart of spot, calculate curve for the first phi bin - otherwise just shift
	    if ( j==0){
	      	      //std::cout << "starting ComputeAngles" << std::endl;

	      curve = ComputeAngles(&curve, defltoa); 	
	      //std::cout << "starting ComputeCurve " << std::endl;
	      curve = ComputeCurve(&curve);
	      //std::cout << "starting TimeDelays" << std::endl;
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
		//if (p==0 && i==0) std::cout << "before adding: flux[0][0] =" << curve.f[0][0] << std::endl;

		flxcurve->f[pp][i] += curve.f[pp][q];
		//if (p==0 && i==0) std::cout << "after adding j=" << j << ": flux[0][0] =" << flxcurve->f[0][0] << std::endl;
	      }
	    } // ending Add curves
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
	    for (unsigned int pp(0); pp < curve.cbands; pp++ )
	      for ( unsigned int i(0); i < numbins; i++ ) 
		curve.f[pp][i] = Temp[pp][i];	  	
    
	    curve = ShiftCurve(&curve,phishift);
	    
	    for( unsigned int pp(0); pp < curve.cbands; pp++ )
	      for ( unsigned int i(0); i < numbins; i++ ) {
		//Flux[pp][i] += curve.f[pp][i]*phishift/dphi      ;
		//	if (pp==0 && i==0) std::cout << "before fractional shift: flux[0][0] =" << flxcurve->f[0][0] << std::endl;
		flxcurve->f[pp][i] +=  curve.f[pp][i]*phishift/dphi ;
		//if (pp==0 && i==0) std::cout << "fractional shift: flux[0][0] =" << flxcurve->f[0][0] << std::endl;
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

    /**********************************************************/
    /* SECOND HOT SPOT -- Mirroring of first hot spot         */
    /**********************************************************/

    //Setting new parameters for second spot
    if ( two_spots ) {
    	std::cout << "Starting Spot 2" << std::endl;
    	incl_2 = Units::PI - incl_1 + d_incl_2; // keeping theta the same, but changing inclination
    	curve.para.incl = incl_2;
    	theta_2 = theta_1 - d_theta_2; //d_theta_2 = 7.78 degrees, for 56+131.78 - 180
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

	    //double phi_edge = Units::PI-curve.para.phi_k[k]; // from our convention for second spot, true phi_edge should be pi - SpotShape:phi_k 
	    double phi_edge = (Units::PI * 2 * phase_2)-curve.para.phi_k[k]; // for phase difference of 0.5626 cycles

	    dphi = 2.0*Units::PI/(numbins*1.0);

	    mu_2 = cos(thetak);
	    if (fabs(mu_2) < DBL_EPSILON) mu_2 = 0.0;

	    if ( mu_2 < 0.0){
	      //std::cout << "Southern Hemisphere! mu=" << mu_2 << std::endl;
	      //mu_2 = fabs(mu_2);
	      //thetak = Units::PI - thetak;
	      //curve.para.incl = Units::PI - incl_2;
	    }

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
		  flxcurve->f[p][i] += curve.f[pp][i]*phishift/dphi      ;
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

    
    /***************************************/
    /* START OF INSTRUMENT EFFECT ROUTINES */
    /***************************************/

    // Interpolate to create all the other energy bands
    
    unsigned int q = 0;
    unsigned int index = 0;
    unsigned int factor = curve.tbands/curve.cbands;

    double tvec[5], fvec[5], err;
    int npt=2;
    int zip=0;

    for (unsigned int p = 0; p < curve.tbands; p++){
      q = (p/factor);
      index = p - (q*factor);
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
      std::cout << "ISM" << std::endl;
      curve = Attenuate(&curve,curve.flags.attenuation,nh);
    }
    
    /******************************************/
    /*         ADDING BACKGROUND              */
    /******************************************/
/*
    if (background_file_is_set){
      curve = Background_list(&curve, background_file);      
    } 

    else {
      if (background != 0)
	// Add a constant background in all bands
	for (unsigned int p = 0; p < numbands; p++){	   
	  for (unsigned int i = 0; i < numbins; i++){	     
	    curve.f[p][i] += background;
	    //	    if (p==0)
	    //std::cout << "Added background: i="<<i << " flux = " << curve.f[p][i] << std::endl;
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
	
*/
    /******************************************/
    /*  APPLYING INSTRUMENT RESPONSE CURVE    */
    /******************************************/
    if (curve.flags.inst_curve > 0){
      std::cout << "Applying Instrument Response" << std::endl;
      curve = Inst_Res2(&curve, curve.flags.inst_curve);
    }



    /* Create Phase-independent Powerlaw Background */

    (*flxcurve) = PowerLaw_Background(&curve,1.0,-2.0);    
    if (curve.flags.inst_curve > 0){
      //std::cout << "Applying Instrument Response to the Powerlaw Background" << std::endl;
      (*flxcurve) = Inst_Res2(flxcurve, curve.flags.inst_curve);
    }

   


    // If databins < numbins then rebin the theoretical curve down 
    //std::cout << "databins = " << databins << ", numbins = " << numbins << std::endl;
    //std::cout << " Rebin the data? " << std::endl;
    if (databins < numbins) {
        obsdata.numbins = databins;
	//  std::cout << " Rebin the data! " << std::endl;
        curve = ReBinCurve(&obsdata,&curve);
	    (*flxcurve) = ReBinCurve(&obsdata,flxcurve);
        numbins = databins;
    }


    numbands = curve.fbands;
    // Count the photons!
    double spotcounts = 0.0;
    double bkgcounts = 0.0;
    /******************************************/
    /*     MULTIPLYING OBSERVATION TIME       */
    /******************************************/
    for (unsigned int p = 0; p < numbands; p++){
       for (unsigned int i = 0; i < numbins; i++){          
	  //curve.f[p][i] *= obstime/(databins);  
	  curve.f[p][i] *= obstime;  
	  spotcounts += curve.f[p][i];
	  bkgcounts += flxcurve->f[p][i];
        }
    }
    
    numbands = curve.tbands;

    //double totalcounts = 0.0;
    for (unsigned int p = 0; p < numbands; p++){  
        for (unsigned int i = 0; i < numbins; i++){           
	  //curve.f[p][i] *= 1e4/spotcounts;
	  //curve.f[p][i] += background[p];
	  curve.f[p][i] += flxcurve->f[p][i] * 1e4 / bkgcounts;
	  //curve.f[p][i] = flxcurve->f[p][i] * 1e6 / bkgcounts;
	  //totalcounts += curve.f[p][i];
        }
    }
    //std::cout << "Total Counts = " << totalcounts << std::endl;
    
           


    /************************************************************/
    /* If data file is set, calculate chi^2 fit with simulation */
    /************************************************************/
	//std::cout << obsdata.f[0][7] << " " << curve.f[0][7] << " " << std::endl;
    //if ( datafile_is_set ) {
    	std::cout << "calculating chi squared" << std::endl;
    	chisquared = ChiSquare ( &obsdata, &curve );
    //}
    
    std::cout << "Spot: m = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN 
	      << " Msun, r = " << Units::nounits_to_cgs(req, Units::LENGTH )*1.0e-5 
	      << " km, f = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) 
	      << " Hz, i = " << incl_1 * 180.0 / Units::PI 
	      << ", e = " << theta_1 * 180.0 / Units::PI 
	      << ", X^2 = " << chisquared 
	      << std::endl;    


    //std::cout << "Chi^2[0] = " << obsdata.chi[0] << std::endl;

	if (std::isnan(chisquared)||std::isinf(chisquared)) {
		mexPrintf("Chisquared is nan! Setting it to 10^11.\n");
		chisquared = 100000000000.0;
	}
	else {
		mexPrintf("X^2 = %f\n\n", chisquared);
	}
	
	if (ts > 1.0) ts -= 1.0; // to make sure ts is in [0:1], for printing out
    //std::cout << "Pre printing to everything." << std::endl;
    
    
    std::cout << "Final value of flux[0][0] = " << curve.f[0][0] << std::endl;
 
  
   



 
    // Free previously allocated memory

   /*  if (curve.flags.beaming_model == 10 || curve.flags.beaming_model >= 12  ){ // *cole* McPHAC Hydrogen Atmosphere        
       std::cout << "Free mccangle!" << std::endl;
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
     }*/

	/*****************************************************************/
    /* DUMPING DATA INTO MATLAB OUTPUT                               */
    /* For filewriting, need to use fopen, fprintf, and fclose (etc) */
    /*****************************************************************/
    //std::cout << "Pre output." << std::endl;
	
    chiOut[0] = chisquared; // saved this to theOutput[0] at top just after declarations
   // std::cout << "Output 1." << std::endl;
    dimSize[1] = 3; // number of columns (1: time, 2: flux in 1st energy band, 3: flux in 2nd energy band)
    //std::cout << "Output 2." << std::endl;

    dimSize[0] = (int)numbins; // number of rows ( = numbins)
    //std::cout << "Output 3." << std::endl;

    theOutput[1] = mxCreateNumericArray(2, dimSize, mxDOUBLE_CLASS, mxREAL); // formatting/setting up the output
    //std::cout << "Output 4." << std::endl;

    //curveOut = mxGetPr(theOutput[1]); // matlab and mex love pointers
	
        free_dmatrix(curve.defl.psi,0,1001,0,301);
        free_dmatrix(curve.defl.b,0,1001,0,301);
        free_dmatrix(curve.defl.dcosa,0,1001,0,301);
        free_dmatrix(curve.defl.toa,0,1001,0,301);
        

	free_dvector(curve.mccinte,0,1595001);
    //free_dvector(curve.mccangl,0,51);
    //free_dvector(curve.mcloget,0,101);
	std::cout << "Freed the memory" << std::endl;
    
    //delete model;

  
}



