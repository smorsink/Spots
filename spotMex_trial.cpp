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
#include "nrutil.h"
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
  std::ofstream allflux;  // Prints information about the star's surface area
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

  double SurfaceArea(0.0);
  double E0, L1, L2, DeltaE, rot_par;
    
  unsigned int NS_model(1),       // Specifies oblateness (option 3 is spherical)
    spectral_model(0),    // Spectral model choice (initialized to blackbody)
    beaming_model(0),     // Beaming model choice (initialized to isotropic)
    numbins(MAX_NUMBINS), // Number of time or phase bins for one spin period; Also the number of flux data points
    databins(MAX_NUMBINS),
    numphi(1),            // Number of azimuthal (projected) angular bins per spot
    numtheta(1),          // Number of latitudinal angular bins per spot
    spotshape(0), // Spot shape; 0=standard
    numbands(NCURVES); // Number of energy bands;


  char out_file[256] = "flux.txt",    // Name of file we send the output to; unused here, done in the shell script
       bend_file[256] = "angles100.txt", 
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
    	 bend_file_is_set(false),
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
        std::cout << "one spot" << std::endl;
	}
    if (spots_2 == 2){ 
		two_spots = true;
        std::cout << "two_spots" << std::endl;
	}
    //E0 = mxGetScalar(theInput[18]); // monochromatic band energy; double
	obsdata.t = mxGetPr(theInput[18]); // array of double
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
        //std::cout << "reorganizing bending data" << std::endl;

        curve.defl.psi = dmatrix(0,1001,0,301);
        curve.defl.b = dmatrix(0,1001,0,301);
        curve.defl.dcosa = dmatrix(0,1001,0,301);
        curve.defl.toa = dmatrix(0,1001,0,301);
        curve.defl.psi_b = dvector(0,301);
        curve.defl.b_psi = dvector(0,301);
        curve.defl.dcosa_dcosp_b = dvector(0,301);
        curve.defl.toa_b = dvector(0,301);

       
        
        for (int j = 0; j < 1001; j++){
            for (int i = 0; i < 301; i++){
                curve.defl.b[j][i] = bend_data_b[j*301+i];
                curve.defl.psi[j][i] = bend_data_psi[j*301+i];
                curve.defl.dcosa[j][i] = bend_data_dcosa[j*301+i];
                curve.defl.toa[j][i] = bend_data_toa[j*301+i];
            }
        }
        //std::cout << "M/R[5] = " << curve.defl.mr[5] << std::endl;
        //std::cout << "b[5][10] = " << curve.defl.b[5][10] << std::endl;
    }

    spotshape = mxGetScalar(theInput[25]);
        //std::cout << "spotshape = " << spotshape << std::endl;
    double ObsTime;
    ObsTime = mxGetScalar(theInput[26]); // in Mega-seconds
    ObsTime *= 1e6; // convert to seconds
        //std::cout << "ObsTime = " << ObsTime << std::endl;
    
       std::cout << "Spot: m = " << mass
	      << " Msun, r = " << req
	      << " km, f = " << omega 
	      << " Hz, i = " << incl_1 
	      << ", e = " << theta_1  
	      << ", ts = " << ts
               << ", ObsTime = " << ObsTime << "s" 
	      << std::endl;   
    
    
       
    for (int i = 0; i < numbands; i++){
    obsdata.f[i] = mxGetPr(theInput[27+3*i]); // array of double
    obsdata.err[i] = mxGetPr(theInput[28+3*i]); // array of double
    background[i] = mxGetScalar(theInput[29+3*i]);
        std::cout << " bg["<<i<<"]=" << background[i] ;
    }
       std::cout << std::endl;
    

   
   
 
 	
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
    curve.numbands = numbands;
    //curve.para.rsc = r_sc;
    //curve.para.Isc = I_sc;

    //numphi = numtheta; // code currently only handles a square mesh over the hotspot
  
    curve.flags.ignore_time_delays = ignore_time_delays;
    curve.flags.spectral_model = spectral_model;
    curve.flags.beaming_model = beaming_model;
    curve.flags.ns_model = NS_model;
    curve.flags.spotshape = spotshape;



   // Define the Observer's Spectral Model

    if (curve.flags.spectral_model == 0){ // NICER: Monochromatic Obs at E0=1keV
      curve.para.E0 = 1.0;
      curve.numbands = 1;
    }
    if (curve.flags.spectral_model == 1){ // NICER Line
      curve.para.E0 = E0; // Observed Energy in keV
      curve.para.L1 = L1; // Lowest Energy in keV in Star's frame
      curve.para.L2 = L2; // Highest Energy in keV
      curve.para.DeltaE = DeltaE; // Delta(E) in keV
    }

 
    obsdata.shift = ts;
    obsdata.numbins = databins;


   
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
        
      //vstd::cout << " r_pole = " <<  model->R_at_costheta(1.0) << std::endl;

       
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
        printf("Spherical Model. ");
    }
    else {
    	/*
        throw(Exception("\nInvalid NS_model parameter. Exiting.\n"));
        return -1;
        */
    }

   
    if (NS_model != 3)
        rspot = model->R_at_costheta(cos(theta_1));
      else
        rspot = req;
    /*
    std::cout << "Spot Temperature = " << spot_temperature << "keV" << std::endl;
    std::cout << "Obs Temperature = " << spot_temperature*sqrt(1-2.0*mass_over_req*req/rspot) << std::endl;
    std::cout << "Correct Temperature = " << 2.0*sqrt(3.0/5.0) << std::endl;
    
    if ( fabs(spot_temperature*sqrt(1-2.0*mass_over_req*req/rspot) - 2.0*sqrt(3.0/5.0)) > 0.05 ){
            chisquared = 1e10;
            std::cout << " Wrong temperature!!!!!!!!!" << std::endl;
    }
    else{
            // Correct Temperature!
    */
    
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

      curve = SpotShape(pieces,p,numtheta,theta_1,rho, &curve, model);

      double deltatheta(0.0);

    // Looping through the mesh of the spot
        for (unsigned int k(0); k < numtheta; k++) { // Loop through the circles
     
      deltatheta = curve.para.dtheta[k];

      double thetak = curve.para.theta_k[k];
      // std::cout << "k = " << k 
      //        << "theta = " << thetak
      //        << std::endl;

      double phi_edge = curve.para.phi_k[k];

      dphi = 2.0*Units::PI/(numbins*1.0);

      mu_1 = cos(thetak);
      if (fabs(mu_1) < DBL_EPSILON) mu_1 = 0.0;

      //if ( mu_1 < 0.0){
        //std::cout << "Southern Hemisphere! mu=" << mu_1 << std::endl;
              //mu_1 = fabs(mu_1);
              //thetak = Units::PI - thetak;
              //curve.para.incl = Units::PI - incl_1;
      //}

      if (NS_model != 3)
        rspot = model->R_at_costheta(mu_1);
      else
        rspot = req;

      // Values we need in some of the formulas.
      cosgamma = model->cos_gamma(mu_1);   // model is pointing to the function cos_gamma
      curve.para.cosgamma = cosgamma;

      curve.para.radius = rspot; // load rspot into structure
      curve.para.mass_over_r = mass_over_req * req/rspot;

      //std::cout << " Entering defltoa M/R = " << curve.para.mass_over_r << std::endl; 
      OblDeflectionTOA* defltoa = new OblDeflectionTOA(model, mass, curve.para.mass_over_r , rspot); 
      //std::cout << " Entering bend M/R = " << curve.para.mass_over_r << std::endl; 
      curve = Bend(&curve,defltoa);
      //std::cout << " Max b/R = " << curve.defl.b_R_max << curve.defl.b_psi[3*NN] << std::endl; 


      numphi = 2.0*phi_edge/dphi;
      phishift = 2.0*phi_edge - numphi*dphi;

      curve.para.dS = pow(rspot,2) * sin(thetak) * deltatheta * dphi;
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
       
      //std::cout << numphi << " " << phi_edge << " " << dphi << " " << phishift << " " << curve.para.dS << std::endl;
       
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
        else{ //no
            pieces=1;
        }

      

        
    for (unsigned int p(0);p<pieces;p++){

      curve = SpotShape(pieces,p,numtheta,theta_2,rho, &curve, model);

      double deltatheta(0.0);

    // Looping through the mesh of the spot
        for (unsigned int k(0); k < numtheta; k++) { // Loop through the circles
     
      deltatheta = curve.para.dtheta[k];

      double thetak = curve.para.theta_k[k];

      double phi_edge = Units::PI-curve.para.phi_k[k]; // from our convention for second spot, true phi_edge should be pi - SpotShape:phi_k 

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

      numphi = 2.0*(Units::PI-phi_edge)/dphi;
      phishift = 2.0*(Units::PI-phi_edge)- numphi*dphi;
      //phishift = 2.0*(Units::PI-phi_edge)+ numphi*dphi;

      curve.para.dS = pow(rspot,2) * sin(thetak) * deltatheta * dphi * curve.para.gamma_k[k];
      curve.para.theta = thetak;

      if (numtheta==1){  //For a spot with only one theta bin (used for small spot)
        numphi=1;
        phi_edge=-1*Units::PI;
        dphi=0.0;
        phishift = 0.0;
        curve.para.dS = 2.0*Units::PI * pow(rspot,2) * (1.0 - cos(rho)) ;
        if ( spotshape == 1 ) curve.para.dS /= curve.para.gamma_k[k];
        if ( spotshape == 0 ) curve.para.dS *= curve.para.gamma_k[k];
      }
      //std::cout << numphi << " " << phi_edge << " " << dphi << " " << phishift << " " << curve.para.dS << std::endl;
       
      if ( NS_model != 3 ) curve.para.dS /= curve.para.cosgamma;

      for ( unsigned int j(0); j < numphi; j++ ) {// looping through the phi divisions
        
        curve.para.phi_0 =  phi_edge + (j+0.5)*dphi;            
        //std::cout << phi_edge << " " << dphi << " " << curve.para.phi_0 << std::endl;
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
    } // End Standard Case of second spot
    

    } // closing second spot
            
    // You need to be so super sure that ignore_time_delays is set equal to false.
    // It took almost a month to figure out that that was the reason it was messing up.
     
    /*******************************/
    /* ADDING BACKGROUND TO CURVE  */
    /*******************************/

    for (unsigned int p = 0; p < numbands; p++){
        //std::cout << "band " << p << std::endl; 
        for (unsigned int i = 0; i < numbins; i++){
            //std::cout << Flux[p][i] << std::endl;
            Flux[p][i] += background[p];
                Flux[p][i] *= ObsTime;  
        }
    }

    /*******************************/
    /* NORMALIZING THE FLUXES TO 1 */
    /*******************************/
	    
    // Normalizing the flux to 1 in low energy band. 
    if ( normalize_flux ) {
      normcurve = Normalize1( Flux, numbins );
      std::cout << "Normalize!!" << std::endl;

      // Add background to normalized flux
      // background implemented in previous section.
    /*  
      for ( unsigned int i(0); i < numbins; i++ ) {
	for ( unsigned int p(0); p < numbands; p++ ) {
	  Flux[p][i] = normcurve.f[p][i] + background[p];
	}
      }
    */

      
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

    if (databins < numbins)
      curve = ReBinCurve(&obsdata,&curve);
     
    /***************************************************/
    /* Data file should already be available           */
    /***************************************************/
	std::cout << "numbins =" << numbins << " numbands =" << numbands << std::endl;

    chisquared = ChiSquare ( &obsdata, &curve );
    
    chisquared = 0.0;
    
    for (unsigned int i(0);i<30;i++){
     chisquared += obsdata.chi[i];   
    }
    
    
    //std::cout << "Warning chi^2 is only for band 0+1+2 " << std::endl;
    //chisquared = obsdata.chi[0] + obsdata.chi[1] + obsdata.chi[2];
    
    //}
    
    std::cout << "Spot: m = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN 
	      << " Msun, r = " << Units::nounits_to_cgs(req, Units::LENGTH )*1.0e-5 
	      << " km, f = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) 
	      << " Hz, i = " << incl_1 * 180.0 / Units::PI 
	      << ", e = " << theta_1 * 180.0 / Units::PI 
	      << ", X^2 = " << chisquared 
	      << std::endl;   
    
   
 

	if (std::isnan(chisquared)||std::isinf(chisquared)) {
		mexPrintf("Chisquared is nan! Setting it to 1,000,000.\n");
		chisquared = 100000000000.0;
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
   // std::cout << "Output 1." << std::endl;
    dimSize[1] = 3; // number of columns (1: time, 2: flux in 1st energy band, 3: flux in 2nd energy band)
    //std::cout << "Output 2." << std::endl;

    dimSize[0] = (int)numbins; // number of rows ( = numbins)
    //std::cout << "Output 3." << std::endl;

    theOutput[1] = mxCreateNumericArray(2, dimSize, mxDOUBLE_CLASS, mxREAL); // formatting/setting up the output
    //std::cout << "Output 4." << std::endl;

    //curveOut = mxGetPr(theOutput[1]); // matlab and mex love pointers
	

    /*********************************/
    /* Writing lightcurve output     */
    /*********************************/

/*        
    out.open("output_spotMex.txt");

    out << "Spot: m = " << Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN 
	      << " Msun, r = " << Units::nounits_to_cgs(req, Units::LENGTH )*1.0e-5 
	      << " km, f = " << Units::nounits_to_cgs(omega, Units::INVTIME)/(2.0*Units::PI) 
	      << " Hz, i = " << incl_1 * 180.0 / Units::PI 
	      << ", e = " << theta_1 * 180.0 / Units::PI 
                << ", ts = " << ts
                << ", ObsTime = " << ObsTime 
	      << ", X^2 = " << chisquared 
	      << std::endl;    
 
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
                //out << 1e-3 << "\t" ; // Fake errors. 
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
                //out << 1e-6 << "\t" ; // Fake errors.  
            }
            //out << i;
            out << std::endl;
        }
    }



      //<< "# Column 4: Number flux (photons/(cm^2 s)) in the energy band " << E_band_lower_1 << " keV to " << E_band_upper_1 << " keV \n"
      //<< "# Column 5: Number flux (photons/(cm^2 s)) in the energy band " << E_band_lower_2 << " keV to " << E_band_upper_2 << " keV \n"

    
  


    out.close();

*/
     free_dmatrix(curve.defl.psi,0,1001,0,301);
        free_dmatrix(curve.defl.b,0,1001,0,301);
        free_dmatrix(curve.defl.dcosa,0,1001,0,301);
        free_dmatrix(curve.defl.toa,0,1001,0,301);
        
        free_dvector(curve.defl.psi_b,0,301);
        free_dvector(curve.defl.b_psi,0,301);
        free_dvector(curve.defl.dcosa_dcosp_b,0,301);
        free_dvector(curve.defl.toa_b,0,301);
    
    
    delete model;
}
