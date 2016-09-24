/***************************************************************************************/
/*                                     Atmo.cpp

    This holds the ComputeCurve routine used in Spot.cpp and atmosphere routines that are 
    used by ComputeCurve.

	Was split from Chi.cpp, the large files with every computational routines.

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



/**************************************************************************************/
/* ComputeCurve:                                                                      */
/*              computes the flux of each light curve                                 */
/*																					  */
/* pass: angles = all the angles necessary to compute the flux correctly;             */
/*                computed in the routine/method/function above [radians or unitless] */
/**************************************************************************************/
class LightCurve ComputeCurve( class LightCurve* angles ) {
	
	/******************************************/
	/* VARIABLE DECLARATIONS FOR ComputeCurve */
	/******************************************/
	
    class LightCurve curve;

    std::ifstream input;
    std::ofstream output;
    std::ofstream ttt;


    bool infile_is_set;   // If the user has specified an input file
    double 
           mass,               // Mass of the star, in M_sun
           radius,             // Radius of the star at the spot, in km
      mass_over_r,
           temperature,        // Temperature of the spot, in keV
      //           E_mono,             // Energy observed, in keV, for each monochromatic energy
           E_band_lower_1,     // Lower bound of energy band for flux integration, in keV
           E_band_upper_1,     // Upper bound of energy band for flux integration, in keV
           E_band_lower_2,     // Lower bound of energy band for flux integration, in keV
           E_band_upper_2,     // Upper bound of energy band for flux integration, in keV
           redshift,           // Gravitational redshift = 1 + z = (1-2M/R)^{-1/2}
           bolo,               // Bolometric flux; bolo = sigma T^4/pi
           gray(1.0);          // Graybody factor (when = 1, not effective)
        
    double E0, E1, E2, DeltaE, E_obs;

    unsigned int numbins(MAX_NUMBINS);  // Time bins of light curve (usually 128)
    unsigned int numbands(NCURVES);  // Number of Energy Bands
        
    std::vector< double > totflux(MAX_NUMBINS, 0.0); // integrated flux

    std::vector< bool > nullcurve(NCURVES,true); // true means that the curve is zero everywhere

    //std::vector< double > softbb(MAX_NUMBINS, 0.0);  // blackbody soft flux
    //std::vector< double > softcm(MAX_NUMBINS, 0.0);  // compton soft flux

    // double softbbave(0.0), softcmave(0.0), hardave(0.0);  // average value of softbb and softcm and high energy band

    /*********************/
    /* SETTING THINGS UP */
    /*********************/
    
    // One monochromatic energy, hardwired value, in keV
    //    E_mono = 1.0;

    curve = (*angles);
    mass = curve.para.mass;                     // unitless
    radius = curve.para.radius;                 // unitless
    mass_over_r = curve.para.mass_over_r;
    temperature = curve.para.temperature;       // T in keV 
    infile_is_set = curve.flags.infile_is_set;
    numbins = curve.numbins;
    numbands = curve.numbands;
    E_band_lower_1 = curve.para.E_band_lower_1;     // in keV
    E_band_upper_1 = curve.para.E_band_upper_1;     // in keV
    E_band_lower_2 = curve.para.E_band_lower_2;     // in keV
    E_band_upper_2 = curve.para.E_band_upper_2;     // in keV
    E0 = curve.para.E0;
    E1 = curve.para.E1;
    E2 = curve.para.E2;
    DeltaE = curve.para.DeltaE;
   
    redshift = 1.0 / sqrt( 1 - 2.0 * mass_over_r);

    //bolo = 2.0e12/15.0 * pow(Units::H_PLANCK,-3) * pow(Units::C,-2) * pow( Units::PI * Units::EV * temperature, 4);
    //bolo *= 1.0e-3/Units::EV;

    bolo = 2.0e9 * 2.404 * Units::C * pow(temperature * Units::EV/(Units::H_PLANCK * Units::C) , 3); // use this one! probably!
    // the e9 in the beginning is for changing T^3 from keV to eV
    // 2.404 comes from evaluating Bradt equation 6.17 (modified, for photon number count units), using the Riemann zeta function for z=3

    if (curve.flags.beaming_model == 3){ // Hydrogen Atmosphere
        Read_NSATMOS(curve.para.temperature, curve.para.mass, curve.para.radius); // Reading NSATMOS FILES Files
        cout << "Using hydrogen atmosphere" << endl;
    }

    if (curve.flags.beaming_model == 4){ // Hydrogen Atmosphere
        Read_NSX(curve.para.temperature, curve.para.mass, curve.para.radius); // Reading NSATMOS FILES Files
        cout << "Using helium atmosphere" << endl;
    }

    for ( unsigned int i(0); i < numbins; i++ ) { // Compute flux for each phase bin

    if ( curve.dOmega_s[i] != 0.0 ) {

	if (std::isnan(gray) || gray == 0) std::cout << "gray = " << gray << std::endl;
	if (std::isnan(curve.eta[i]) || curve.eta[i] == 0) std::cout << "eta[i="<<i<<"] = " << curve.eta[i] << std::endl;
	if (std::isnan(redshift) || redshift == 0) std::cout << "redshift = " << redshift << std::endl;

    /*******************************************************************/
    /* COMPUTING LIGHT CURVE FOR MONOCHROMATIC ENERGY, p = 0           */
    /*      First computes in [erg/(s cm^2 Hz), converts to            */
    /*      photons/(s cm^2 keV)                                       */
    /*******************************************************************/
	if (curve.flags.spectral_model == 0){ // Monochromatic Observation
    
        if (curve.flags.beaming_model == 0){ // Blackbody, no beaming
	       // Moonochromatic light curve in energy flux erg/(s cm^2 Hz)
            curve.f[0][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * BlackBody(temperature,E0*redshift/curve.eta[i]); 
	       // Units: erg/(s cm^2 Hz)
	       //Convert to photons/(s cm^2 keV)
	        curve.f[0][i] *= (1.0 / ( E0 * Units::H_PLANCK )); // Units: photons/(s cm^2 keV)
            if (curve.f[0][i] != 0.0) nullcurve[0] = false;
	    }

	if (curve.flags.beaming_model == 1 || curve.flags.beaming_model == 5 || curve.flags.beaming_model == 6 ){ // Blackbody, graybody factor
	  //            gray = Gray(curve.cosbeta[i]*curve.eta[i]); // gets the graybody factor

	  if (curve.flags.beaming_model == 1)
	    gray = Gray(curve.cosbeta[i]*curve.eta[i]); // gets the graybody factor
	  if (curve.flags.beaming_model == 5)
	    gray = pow( curve.cosbeta[i]*curve.eta[i], 2.0);
	  if (curve.flags.beaming_model == 6)
	    gray = 1.0 - pow( curve.cosbeta[i]*curve.eta[i], 2.0);


	  curve.f[0][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * BlackBody(temperature,E0*redshift/curve.eta[i]); 
	  curve.f[0][i] *= (1.0 / ( E0 * Units::H_PLANCK )); // Units: photons/(s cm^2 keV)
	  if (curve.f[0][i] != 0.0) nullcurve[0] = false;
        }

        if (curve.flags.beaming_model == 2){ // Funny Line Emission, not calculated
        }

        if (curve.flags.beaming_model == 3){ // Hydrogen Atmosphere
            curve.f[0][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * Hydrogen(E0 * redshift/curve.eta[i], curve.cosbeta[i]*curve.eta[i]);
            curve.f[0][i] *= (1.0 / ( E0 * Units::H_PLANCK ));
        }

        if (curve.flags.beaming_model == 4){ // Helium Atmosphere
            curve.f[0][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * Helium(E0 * redshift/curve.eta[i], curve.cosbeta[i]*curve.eta[i]);
            curve.f[0][i] *= (1.0 / ( E0 * Units::H_PLANCK ));
        }
    }

    /***********************/
    /* FUNNY LINE EMISSION */
    /***********************/
	if (curve.flags.spectral_model == 1){ // Funny Line Emission for NICER

        for (unsigned int p=0; p<numbands; p++){
	    
	    E_obs = E0 + p*DeltaE;
	    //curve.f[p][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * Line(temperature,E_obs*redshift/curve.eta[i],E1,E2); 
	    // Units: erg/(s cm^2 Hz)
	    //Convert to photons/(s cm^2 keV)
	    //curve.f[p][i] *= (1.0 / ( E_obs * Units::H_PLANCK )); // Units: photons/(s cm^2 keV)
	    // Multiply by the width of the energy band;
	    //curve.f[p][i] *= DeltaE; // Units: photons/(s cm^2)

	    //curve.eta[i] = 1.0;

	    curve.f[p][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * LineBandFlux(temperature, (E_obs-0.5*DeltaE)*redshift/curve.eta[i], (E_obs+0.5*DeltaE)*redshift/curve.eta[i], E1, E2); // Units: photon/(s cm^2)

	    if (curve.f[p][i] != 0.0) nullcurve[0] = false;
        }
	}

    /*******************************************************************/
    /* COMPUTING BLACKBODY LIGHT CURVE FOR INTEGRATED FLUX             */
    /* Units: photons/(cm^2 s)                                         */
    /*******************************************************************/

    if (curve.flags.spectral_model == 2){ // Integrated Flux of Energy Bands
        double E_diff = (E_band_upper_1 - E_band_lower_1)/numbands;

        for (unsigned int p = 0; p<numbands; p++){
            if (curve.flags.beaming_model == 0){ //blackbody
                curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * EnergyBandFlux(temperature, (E_band_lower_1+p*E_diff)*redshift/curve.eta[i], (E_band_lower_1+(p+1)*E_diff)*redshift/curve.eta[i]); // Units: photon/(s cm^2)
            }

            if (curve.flags.beaming_model == 1){ //blackbody with graybody
                gray = Gray(curve.cosbeta[i]*curve.eta[i]); // gets the graybody factor
                curve.f[p][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * EnergyBandFlux(temperature, (E_band_lower_1+p*E_diff)*redshift/curve.eta[i], (E_band_lower_1+(p+1)*E_diff)*redshift/curve.eta[i]); // Units: photon/(s cm^2)
            }
            if (curve.flags.beaming_model == 3){ //hydrogen
                curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * AtmosEBandFlux(curve.flags.beaming_model, curve.cosbeta[i]*curve.eta[i], (E_band_lower_1+p*E_diff)*redshift/curve.eta[i], (E_band_lower_1+(p+1)*E_diff)*redshift/curve.eta[i]); // Units: photon/(s cm^2)        
            }
            if (curve.flags.beaming_model == 4){ //helium
                curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * AtmosEBandFlux(curve.flags.beaming_model, curve.cosbeta[i]*curve.eta[i], (E_band_lower_1+p*E_diff)*redshift/curve.eta[i], (E_band_lower_1+(p+1)*E_diff)*redshift/curve.eta[i]); // Units: photon/(s cm^2)        
            }
        }
    }

    if (curve.flags.spectral_model == 3){ // *exactly* Integrated Flux of Energy Bands
        double E_diff = (E_band_upper_1 - E_band_lower_1)/numbands;
        for (unsigned int p = 0; p<numbands; p++){
            if (curve.flags.beaming_model == 4){ //helium
                curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * AtmosEBandFlux2(curve.flags.beaming_model, curve.cosbeta[i]*curve.eta[i], (E_band_lower_1+p*E_diff)*redshift/curve.eta[i], (E_band_lower_1+(p+1)*E_diff)*redshift/curve.eta[i]); // Units: photon/(s cm^2)        
            }
        }
    }

	if (curve.flags.spectral_model == 7) { // Not Used Right Now.
			
	  /*******************************************/
	  /* COMPUTING BOLOMETRIC LIGHT CURVE, p = 0 */
	  /* Units: photons/(cm^2 s)                 */
	  /*******************************************/
    		
	  // Bolometric Light Curve for energy flux erg/(cm^2 s)
	  curve.f[0][i] = bolo * gray * curve.dOmega_s[i] * pow(curve.eta[i],5) * pow(redshift,-4); // Units: erg/(cm^2 s)

	  // Bolometric Light Curve for photon number flux photons/(cm^2 s)
	  curve.f[0][i] = bolo * gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3); // Units: photons/(cm^2 s)
                                     
	  /***************************************************/
	  /* COMPUTING PHOTON NUMBER FLUX FOR AN ENERGY BAND */
	  /*   	1 < p < NCURVES-1                          */
	  /*      Units: photons/(s cm^2)                    */
	  /***************************************************/
    		
	  // First energy band
	  curve.f[NCURVES-2][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * EnergyBandFlux(temperature, E_band_lower_1*redshift/curve.eta[i], E_band_upper_1*redshift/curve.eta[i]); // Units: photon/(s cm^2)
	  // Second energy band
	  curve.f[NCURVES-1][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * EnergyBandFlux(temperature, E_band_lower_2*redshift/curve.eta[i], E_band_upper_2*redshift/curve.eta[i]); // Units: photon/(s cm^2)
	}		
      }
      else { // if curve.dOmega_s[i] == 0.0
	for ( unsigned int p(0); p < numbands; p++) {
	  curve.f[p][i] = 0.0;
	}
      }

    } // ending the for(i) loop
    

	
    /***********************************************************/
    /* DEALING WITH THE TIME DELAYS, REBINNING THE LIGHT CURVE */
    /* This is where the jumpy problems tend to be.            */
    /***********************************************************/

    curve.flags.ignore_time_delays = false;
    		
    if ( !curve.flags.ignore_time_delays ) { // if we are not ignoring the time delays        
      int k(0), j(0); // index placeholders; approximately, k is i+1 and j is i-1
      // but time delays mean that j isn't always i-1
      // used in the linear interpolation
        
      /********************************/
      /* LOOP THROUGH THE LIGHTCURVES */
      /********************************/
    	
      for (unsigned int p(0); p < numbands; p++) {
      // for (unsigned int p(1); p < 2; p++) {

	if (!nullcurve[p]) {
	

	/*********************/
	/* MORE DECLARATIONS */
	/*********************/
    		
	  std::vector< double > newflux(MAX_NUMBINS, 0.0);                  // rebinned flux (to account for photon travel time)
	  unsigned int imax(0), imin(0);                                    // index of element with maximum flux value, minimum flux value
	  double max_discrete_flux(0.0), min_discrete_flux(curve.f[p][0]);  // maximum flux value and minimum flux value assigned to a grid point (discrete)
	  double tx(0), ta(0), tb(0), tc(0);                                // a,b,c: three-point interpolation on a parabola for min and max areas (time a, time b, time c)
	  double fa(0), fb(0), fc(0);                                       // corresponding fluxes for three-point interpolation
	  int ia(0), ic(0);                                                 // indices in array for where flux is fa, fc
	  double temporary1(0.0), temporary2(0.0);                          // makes math easier
	  double maximum(0.0), minimum(100000.0);                           // true (continuous) maximum and minimum flux values
	  double tmin(0.0);                                                 // value of t at the true minimum

	  int ec1(0), ec2(0),j1,j2;
	  double te1, te2, tt1, tt2, slope1, slope2;

	  /**********************************************************/
	  /* START WITH FINDING THE MAXIMUM FLUX OF THE LIGHT CURVE */
	  /**********************************************************/
			
	  /********************************************************/
	  /* FINDING THE DISCRETE MAXIMUM FLUX OF THE LIGHT CURVE */
	  /********************************************************/
    		
	  for ( unsigned int i(0); i < numbins; i++ ) {
	    if ( curve.f[p][i] > max_discrete_flux ) { // tells you where the maximum is
	      imax = i;  
	      max_discrete_flux = curve.f[p][i];
	    }
	  }
	        
	  /******************************************/
	  /* FINDING THE TRUE MAXIMUM VALUE OF FLUX */
	  /******************************************/
    		
	  if ( imax == 0 ) {
	    ia = numbins - 1;
	    ta = curve.t_o[ia] - 1.0;
	  }
	  else { // else imax == -1
	    ia = imax - 1;
	    ta = curve.t_o[ia];
	  }
	  tb = curve.t_o[imax];
	  if ( imax == numbins - 1 ) {
	    ic = 0;
	    tc = curve.t_o[ic] + 1.0;
	  }
	  else { //else imax == 0
	    ic = imax + 1;
	    tc = curve.t_o[ic];
	  }
	  fa = curve.f[p][ia];
	  fb = curve.f[p][imax];
	  fc = curve.f[p][ic];
            
	  /**********************************************/
	  /* NUMERICAL RECIPES, PARABOLIC INTERPOLATION */
	  /* Equation 10.3.1 (our t is their x)         */
	  /**********************************************/
    		
	  if ( ( (tb-ta)*(fb-fc) - (tb-tc)*(fb-fa) ) != 0.0 ) { // otherwise you get a big fat NAN for a flux
	    tx = tb - 0.5 * (pow(tb-ta,2)*(fb-fc) - pow(tb-tc,2)*(fb-fa)) / ((tb-ta)*(fb-fc) - (tb-tc)*(fb-fa));
	    temporary1 =  (fa-fc)/( pow(tc-tx,2) - pow(ta-tx,2)) ;
	  }
	  if ( ( (tb-ta)*(fb-fc) - (tb-tc)*(fb-fa) ) == 0.0 ) { // to avoid dividing by zero
	    /********/      tx = 0.0; // is this what it should be?
	  }
	  if ( ( pow(tc-tx,2) - pow(ta-tx,2) ) != 0.0 ) {
	    temporary1 =  (fa-fc)/( pow(tc-tx,2) - pow(ta-tx,2)) ;
	  }
	  if ( ( pow(tc-tx,2) - pow(ta-tx,2) ) == 0.0 ) {
	    temporary1 = 0.0;
	  }
	  maximum = fb - temporary1*pow(tb-tx,2);
	  // maximum is the maximum flux for the specific "p" light curve
	  // tx is the time of maximum 

	  /************************************************/
	  /* NOW FIND THE MINIMUM FLUX OF THE LIGHT CURVE */
	  /* Note: real light curves won't have eclipses  */
	  /* If eclipsed, then min = 0                    */
	  /************************************************/
    		
	  if ( !curve.eclipse ){// && !curve.ingoing ) {   
            
	    /********************************************************/
	    /* FINDING THE DISCRETE MINIMUM FLUX OF THE LIGHT CURVE */
	    /* If not eclipsed                                      */
	    /********************************************************/
    			
	    for ( unsigned int i(0); i < numbins; i++) {
	      if (curve.f[p][i] < min_discrete_flux) {
		imin = i;
		min_discrete_flux = curve.f[p][i];
	      }
	    }
	        
	    //std::cout << "imin = " << imin << "minflux discrete = " << min_discrete_flux << std::endl;
    
	    /******************************************/
	    /* FINDING THE TRUE MINIMUM VALUE OF FLUX */
	    /* If not eclipsed                        */
	    /******************************************/
    			
	    if ( imin == 0 ) {
	      ia = numbins - 1;
	      ta = curve.t_o[ia] - 1.0;
	    }
	    else {  //else imin == 1
	      ia = imin - 1;
	      ta = curve.t_o[ia];
	    }
	    tb = curve.t_o[imin];
	    if ( imin == numbins-1 ) {
	      ic = 0;
	      tc = curve.t_o[ic]+1.0;
	    }
	    else {  //else imin == 2
	      ic = imin+1;
	      tc = curve.t_o[ic];
	    }
	    fa = curve.f[p][ia];
	    fb = curve.f[p][imin];
	    fc = curve.f[p][ic];
                	
	    /**********************************************/
	    /* NUMERICAL RECIPES, PARABOLIC INTERPOLATION */
	    /* Equation 10.3.1 (our t is their x)         */
	    /**********************************************/

	    if ( ( (tb-ta)*(fb-fc) - (tb-tc)*(fb-fa) ) != 0.0 ) { // otherwise you get a big fat NAN for a flux
	      tmin = tb - 0.5*(pow(tb-ta,2)*(fb-fc) - pow(tb-tc,2)*(fb-fa)) / ((tb-ta)*(fb-fc) - (tb-tc)*(fb-fa));  
	    }
	    if ( ( pow(tc-tmin,2) - pow(ta-tmin,2) ) != 0.0 ) {
	      temporary2 =  (fa-fc)/( pow(tc-tmin,2) - pow(ta-tmin,2)) ;
	    }
	    if ( ( (tb-ta)*(fb-fc) - (tb-tc)*(fb-fa) ) == 0.0) {  // what should happen for dividing by zero?
	      tmin = 0.0;
	    }
	    if ( ( pow(tc-tmin,2) - pow(ta-tmin,2) ) == 0.0 ) { // what should happen for dividing by zero?
	      temporary2 = 0.0;
	    }
           	 	
	    minimum = fb - temporary2 * pow(tb-tmin,2);

	    /*	    std::cout << "imin = " << imin 
		      << "tmin = " << tmin
		      << "True minflux = " << minimum << std::endl;*/
           	 
            } // ending "not eclipsed, not ingoing" section

	  // minimum is the minimum value of flux
	  // tmin is the time of the minimum
	  // only true of not eclipsed

            
	  /*******************************************/
	  /* FOR ECLIPSING LIGHT CURVES, MINIMUM = 0 */
	  /*******************************************/
    		
	  else if ( curve.eclipse ) {
	    minimum = 0.0;

	    // Find the first eclipsed bin; set ec1=i;
	    for (unsigned int i(0);i<numbins;i++){
	      if (curve.f[p][i]==0.0){
		ec1=i;
		break;
	      }
	    }
	    
	    //	    std::cout << "eclipse begins at bin" << ec1 << std::endl;
	    // Extrapolate to find beginning of eclipse (linear extrapolation)
	    j1 = ec1-2;
	    j2 = ec1-1;

	    if (ec1==0){
	      j1+=numbins;
	      j2+=numbins;
	    }
	    if (ec1==1){
	      j1+=numbins;
	    }
	    tt1 = curve.t_o[j1];
	    if (ec1<2)
	      tt1 +=-1.0;
	    tt2 = curve.t_o[j2];
	    if (ec1==0)
	      tt2+=-1.0;
	    


	    slope1 = (curve.f[p][j1]-curve.f[p][j2])/(tt1-tt2);
	
	    if (slope1 != 0.0)
	      te1  = - (curve.f[p][j2])/slope1 + tt2 ;
	    else
	      te1 = tt1;


	    // After the eclipse, find the first bin with non-zero flux
	    for (unsigned int i(ec1);i<numbins;i++){
	      if (curve.f[p][i]>0.0){
		ec2=i;
		break;
	      }
	    }
	    //std::cout << "eclipse ends at bin" << ec2 << std::endl;
	    j1=ec2;
	    j2=ec2+1;
	    if (ec2==numbins-1){
	      j2+=-numbins;
	    }
	    tt1=curve.t_o[j1];
	    tt2=curve.t_o[j2];
	    if (ec2==numbins-1)
	      tt2+=1.0;

	    slope2 = (curve.f[p][j1]-curve.f[p][j2])/(tt1-tt2);
	    te2 = tt1 - curve.f[p][j1]/slope2;
	    //std::cout << "Eclipse ends at time = " << te2 << std::endl;


	  } // ending "yes eclipsed
            
            /****************************************************************************/
	  /* FOR NOT ECLIPSED AND YES INGOING PHOTONS, POLITELY SET TO ZERO AND CRASH */
	  /****************************************************************************/
    		
	  else {    //if there is no eclipse and the photon is intially ingoing
	    throw( Exception(" Photons are ingoing and not eclipsing. This is an issue that the code cannot handle.") );
	    minimum = 0.0;
	    maximum = 0.0;
	    break;
	  } // ending "not eclipsed, yes ingoing" section
            
	  /****************************************************/
	  /* COMPUTING THE PULSE FRACTION FOR THE LIGHT CURVE */
	  /****************************************************/
   			
	  curve.minFlux[p] = minimum;
	  curve.maxFlux[p] = maximum;
	  curve.pulseFraction[p] = (curve.maxFlux[p] - curve.minFlux[p]) / (curve.maxFlux[p] + curve.minFlux[p]);

	  curve.asym[p] = (tmin - tx) - 0.5;
	  if (curve.asym[p] < 0.0) curve.asym[p]+=1.0;


	  // Initializing totflux
	  for ( unsigned int i(0); i < numbins; i++ )
	    totflux.at(i) = 0.0;
			
		           		
	  /**************************************************************/
	  /* ADDING FLUXES FROM ALL THE PHASE BINS AND OTHER FUN THINGS */
	  /**************************************************************/

	  for ( unsigned int i(0); i < numbins; i++ ) {  // for-i-loop, looping through the phase bins
            
	    k = i + 1;
	    j = i;		 
	    if ( k == static_cast<int>(numbins) ) 
	      k = 0;
	 		
	    /* if (curve.eclipse) {  // begin eclipse section
	      if (curve.f[p][k] == 0.0) {
		if (curve.f[p][j] != 0.0) {
		  k = i;
		  j = i-1;
		}
		else if (curve.f[p][j-1] != 0.0) {
		  k = j - 1;
		  j = j - 2;
		}
	      }
	      else {
		if ( curve.f[p][j] == 0.0 ){
		  j = k;
		  k = k + 1;
		}
	      } 
	      } // finished eclipse section*/

	      // Check to see if we're near the maximum
	    if ( (i <= imax + 1 && i >= imax -1) || 
		 (imax == 0 && (i <= 1 || i == numbins-1)) || 
		 (imax == numbins-1 && (i==0 || i >= numbins-2))) { // parabolic interpolation near the maximum
	      if ( imax == numbins-1 && i == 0 ) {
		newflux.at(i) = maximum - temporary1 * pow(curve.t[i] + 1.0 - tx,2);  // intermediate value of flux
	      }
	      else {
		if ( imax == 0 && i == numbins-1 )
		  newflux.at(i) = maximum - temporary1 * pow(curve.t[i] - 1.0 - tx,2);
		else
		  newflux.at(i) = maximum - temporary1 * pow(curve.t[i] - tx,2);
	      }
	      //std::cout << "MAX i = " << i << " time = " << curve.t[i] << " flux  = " << newflux.at(i) << std::endl;
	    }
	            
	    // Not near the maximum
	    else {	            
	      	     
		double t1, t2;
		// Find point to the left of the "ith" point 
		// Time delays mean that j isn't always i-1

		// Check to see if t_o[0] > t_o[i]

		if ( curve.t_o[0] > curve.t[i]){ // if thing
		  if ( i==0){
		    if (curve.t_o[numbins-1] > 1.0){
		      j = numbins-2;
		      t1 = curve.t_o[j]-1.0;

		      k = numbins -1;
		      t2 = curve.t_o[k]-1.0;
		    }
		    else{
		      j = numbins-1;
		      t1 = curve.t_o[j]-1.0;
		      k = 0;
		      t2 = curve.t_o[k];
		    }
		  } // end if i==0
		  else{ //i>0
		    j = i-2;
		    if (j<0) j += numbins;
		    t1 = curve.t_o[j];
		    if (t1 > 1.0) t1 += -1.0;
		    k = j+1;
		    if (k>=numbins) k += -numbins;
		    t2 = curve.t_o[k];
		  } 
		}
		else { //Default
		  j = 0;
		  while ( (curve.t_o[j] <= curve.t[i]) && (j < static_cast<int>(numbins)) ) 
		    j++;
		  j--;
		  if ( j < 0 ) // because otherwise the computer gives us garbage and the flux looks ridiculous
		    j = numbins - abs(j);
		  t1 = curve.t_o[j]; // time to the left of the time we're interested in

		  if ( j == static_cast<int>(numbins) - 1 ) {
		    k = 0;
		  }
		  else { // else5; effectively, if i != 0 because that would make j != numbins-1
		    k = j + 1;
		    while (curve.t_o[k] <= curve.t[i] && k < static_cast<int>(numbins)) 
		      k++;
		  }
		  t2 = curve.t_o[k]; // time to the right of the point we're interested in
		  if (k==0) t2 += 1.0;

		}

		/*if (i==0){
		  std::cout << "bin 0 " << std::endl;
		  std::cout << "j=" << j << " t_o[j}=" << curve.t_o[j] << " t1 = " << t1 << std::endl;
		  std::cout << "k=" << k << " t_o[k}=" << curve.t_o[k] << " t2 = " << t2 << std::endl;
		  }*/

		newflux.at(i) = curve.f[p][j] + (curve.f[p][k]-curve.f[p][j])/(t2-t1) * (curve.t[i]-t1); // linear interpolation!

		if (curve.eclipse){
		  if (curve.t[i] < te1 && curve.t[i] > te1 - 1.0/numbins){
		    //std::cout << "time just before eclipse at te1=" << te1 <<"! i=" << i << " t = " << curve.t[i] << std::endl;
		    newflux.at(i) = slope1 * (curve.t[i] - te1);
		  }
		  if (curve.t[i] < te1 + 1.0/numbins && curve.t[i] > te1)
		    newflux.at(i) = 0.0;

		  if (curve.t[i] > te2 && curve.t[i] < te2 + 1.0/numbins){
		    // std::cout << "eclipse ends at te2=" << te2 <<"! i=" << i << " t = " << curve.t[i] << std::endl;
		    newflux.at(i) = slope2 * (curve.t[i] - te2);
		  }

		  if (curve.t[i] > te1 && curve.t[i] < te2)
		    newflux.at(i) = 0.0;
		}

		if (minimum != 0.0)
		  if ( fabs(newflux.at(i) - minimum)/minimum <= 0.000001) {

		    //std::cout << "near minimum!!!! " << std::endl;
		    if ( i == numbins-1 ) 
		      newflux.at(i) = minimum - temporary2 * pow(curve.t[i] - 1 - tmin,2); 
		    else {
		      newflux.at(i) = minimum - temporary2 * pow(curve.t[i] - tmin,2);
		    }
		    //std::cout << "MIN i = " << i << " time = " << curve.t[i] << " flux  = " << newflux.at(i) << std::endl;



		  }

		//std::cout << "i = " << i << " j = " << j << " k = " << k 
		//	  << "t=" << curve.t[i] << " f=" << newflux.at(i) 
		//	  <<std::endl;
		//else
		//newflux.at(i) = curve.f[p][i];

		// }// end else-not-near the minimum
	    }
	
	    /* newflux vs t_e corresponds to the new re-binned light curve.
	       It corresponds to the same light curve as bolflux vs t_o */
	    
	    if ( newflux.at(i) < 0.0 )
	      newflux.at(i) = 0.0;
	  
	  } // closes the for(i) loop, going through the curves.
    	
	  /****************************************************************/
	  /* SUMMING THE FLUXES AT EACH PHASE BIN ACROSS ALL LIGHT CURVES */
	  /****************************************************************/
    	
	  for ( unsigned int i(0); i < numbins; i++ ) {
            totflux.at(i) += newflux.at(i);   // setting totflux = newflux
	  }
	  //if (p==0)
	  //ttt.open("time.txt", std::ios_base::trunc);

	  for ( unsigned int i(0); i < numbins; i++ ) {
	    /*if ( p==0 ) 
	      ttt << curve.t_o[i] << " " << curve.f[p][i] << " " << curve.t[i] << " " << totflux.at(i) 
	    	  << " " << i
	    	  << std::endl;*/
	    
            curve.f[p][i] = totflux.at(i);
	  }
	  // only plotting versus evenly spaced time -- not using t_o
	}
      }// end for-p-loop
    } // end time delay section
    
    return curve;

} // end ComputeCurve






/**************************************************************************************/
/* Hydrogen:                                                                          */
/*           computes the monochromatic blackbody flux in units of erg/cm^2           */
/*                                                                                    */
/*                                                                                    */
/**************************************************************************************/

// Bi-section search for a value in an array
int Find(double val, std::vector<double> array){
    
    bool found = false;
    int low = 0;
    int high = array.size();
    int mid = 0;
    
    while(low <= high){
        mid  = (low + high) / 2;   // midpoint
        if(array[mid] <= val && array[mid+1] >= val){
            found = true;
            break;
        }else if(val < array[mid]){
            high = mid - 1; // decrease the higher edge 1 step
        }else{
            low = mid + 1;  // increase the lower edge 1 step
        }
    }
    return mid;
}

// Linear interpolation
double Linear(double x,double Xa, double Ya, double Xb, double Yb){
    
    return  (Ya + ((Yb - Ya) * (x - Xa) / (Xb - Xa)));
}

// Linear interpolation in log-log space
double LogLinear(double x,double Xa, double Ya, double Xb, double Yb){
    
    return  (exp(log(Ya) + ((log(Yb) - log(Ya)) * (log(x) - log(Xa)) / (log(Xb) - log(Xa)))));
}

// Find and calculate Linear Interpolation
double Interpolate(double X_INT, std::vector<double> X, std::vector<double> Y){
    
    int RowNum = Find(X_INT,X);
    double d_X1 = X_INT - X[RowNum];
    double d_X2 = X[RowNum+1] - X[RowNum];
    double d_Y = Y[RowNum+1] - Y[RowNum];
    
    return  (Y[RowNum] + (d_Y * d_X1 / d_X2));
}

// Find and calculate Linear Interpolation in log-log space
double LogInterpolate(double X_INT, std::vector<double> X, std::vector<double> Y){
    
    int RowNum = Find(X_INT,X);
    double Ya = Y[RowNum];
    double Yb = Y[RowNum+1];
    double Xa = X[RowNum];
    double Xb = X[RowNum+1];
    return  (exp(log(Ya) + ((log(Yb) - log(Ya)) * (log(X_INT) - log(Xa)) / (log(Xb) - log(Xa)))));
}


// Round value to the nearest value in an array
int Round(int n, double z, std::vector<double> v){
    double mid = (v[n] + v[n+1])/2.0;
    if (z >= mid) {
        return n + 1;
    }else{
        return n;
    }
}

// Set up hydrogen/helium arrays and dummy variables
std::vector<double> Es,F,FF,FFF,FFFF,I,II,III,IIII;
double X, Y, X1, Y1, X2, Y2;

// Read the four hydrogen intensity files to peform the four point interpolation
void Read_NSATMOS(double T, double M, double R){
    
    double delta, lgrav, lt, temp, dump;
    int i_lgrav, i_lt(0), n_lgrav, n_lt, size_logt(10), size_lsgrav(11), size_mu(12);
    char s1[40],s2[40],s3[40],s4[40], atmodir[1024], cwd[1024];
    std::vector<double> logt, lsgrav; 

    //Read in hydrogen atmosphere parameters
    getcwd(cwd, sizeof(cwd));
    sprintf(atmodir,"%s/atmosphere",cwd);
    chdir(atmodir);
    ifstream H_atmo_para;
    H_atmo_para.open("nsatmos-info.txt");

    if(H_atmo_para.is_open()){
        //logt
        for (int i = 1; i <= size_logt; i++){
            H_atmo_para >> temp;
            logt.push_back(temp);  
        }

        //lgrav
        for (int i = 1; i <= size_lsgrav; i++){
            H_atmo_para >> temp;
            lsgrav.push_back(temp);
        }

        //discarding mu choices
        for (int i = 1; i <= size_mu; i++){
            H_atmo_para >> dump;
        }
    H_atmo_para.close();
    }

    //Find correct logt and lgrav paramter choice
    M = Units::nounits_to_cgs(M, Units::MASS);
    R = Units::nounits_to_cgs(R, Units::LENGTH);
    delta = 1 / sqrt(1 - (2 * Units::G * M / (R * Units::C * Units::C)));
    lgrav = log10(delta * Units::G * M / (R * R));
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));
    
    i_lt = Find(lt,logt);
    i_lgrav = Find(lgrav,lsgrav);

    n_lgrav = i_lgrav + 1;
    n_lt = i_lt + 1;

    //Load hydrogen atmosphere files
    sprintf(s1,"nsatmos_edd4.emerg%02d%02d01",n_lt,n_lgrav); 
    sprintf(s2,"nsatmos_edd4.emerg%02d%02d01",n_lt+1,n_lgrav);
    sprintf(s3,"nsatmos_edd4.emerg%02d%02d01",n_lt,n_lgrav+1);
    sprintf(s4,"nsatmos_edd4.emerg%02d%02d01",n_lt+1,n_lgrav+1); 
    ifstream H_table1;
    H_table1.open(s1);
    ifstream H_table2;
    H_table2.open(s2);
    ifstream H_table3;
    H_table3.open(s3);
    ifstream H_table4;
    H_table4.open(s4);
    
    if(H_table1.is_open()){
        while (H_table1 >> temp) {
            F.push_back(temp);
            H_table1 >> temp;
            I.push_back(temp);
        }
    }else{
        cout << "NSATMOS files is not found  " << s1 << endl;
    }
    H_table1.close();
    
    if(H_table2.is_open()){
        while (H_table2 >> temp) {
            FF.push_back(temp);
            H_table2 >> temp;
            II.push_back(temp);
        }
    }else{
        cout << "NSATMOS files is not found  " << s2 << endl;
    }
    H_table2.close();
    
    if(H_table3.is_open()){
        while (H_table3 >> temp) {
            FFF.push_back(temp);
            H_table3 >> temp;
            III.push_back(temp);
        }
    }else{
        cout << "NSATMOS files is not found  " << s3 << endl;
    }
    H_table3.close();
    
    if(H_table4.is_open()){
        while (H_table4 >> temp) {
            FFFF.push_back(temp);
            H_table4 >> temp;
            IIII.push_back(temp);
        }
    }else{
        cout << "NSATMOS files is not found  " << s4 << endl;
    }
    H_table4.close();
    chdir(cwd);
    
    X = lt;
    Y = lgrav;
    
    X1 = logt[i_lt];
    Y1 = lsgrav[i_lgrav];
    X2 = logt[i_lt + 1];
    Y2 = lsgrav[i_lgrav + 1];
}

// Calculate the final interpolated intensity
double Hydrogen(double E, double cos_theta){
    double freq, P, temp, dump;
    double I_int[8],Q[4],R[2];
    int i_mu(0), n_mu, down, up,  size_logt(10), size_lsgrav(11), size_mu(12);
    char atmodir[1024], cwd[1024];
    std::vector<double> mu, F_temp,FF_temp,FFF_temp,FFFF_temp,I_temp,II_temp,III_temp,IIII_temp;

    //Convert energy point to frequency
    freq = 1E3 * E * Units::EV / Units::H_PLANCK;

    //Read in atmosphere parameters
    getcwd(cwd, sizeof(cwd));
    sprintf(atmodir,"%s/atmosphere",cwd);
    chdir(atmodir);
    ifstream H_atmo_para;
    H_atmo_para.open("nsatmos-info.txt");

    if(H_atmo_para.is_open()){
        //discard logt and lgrav values
        for (int i = 1; i <= size_logt+size_lsgrav; i++){
            H_atmo_para >> dump;
        }

        //mu
        for (int i = 1; i <= size_mu; i++){
            H_atmo_para >> temp;
            mu.push_back(temp);
        }
    H_atmo_para.close();
    }
    
    //Find proper mu choice
    for (int m = 0; m < size_mu; ++m) {
        if (cos_theta <= mu[m]) {
            i_mu = m;
        }
    }
    n_mu = i_mu + 1;
    
    //Read and interpolate to proper frequency 
   for (int i = 0; i < 2; ++i) {
       down = (i_mu + i) * 128;
       up = down + 128;
       //cout << down << endl;
       
        for (int j = down; j < up; ++j) {
            F_temp.push_back(F[j]);
            I_temp.push_back(I[j]);
        }
        I_int[i] = LogInterpolate(freq,F_temp,I_temp);
    }
                    
    for (int i = 0; i < 2; ++i) {
        down = (i_mu + i) * 128;
        up = down + 128;
        
        for (int j = down; j < up; ++j) {
            FF_temp.push_back(FF[j]);
            II_temp.push_back(II[j]);
        }
        I_int[i+2] = LogInterpolate(freq,FF_temp,II_temp);
    }
    
    for (int i = 0; i < 2; ++i) {
        down = (i_mu + i) * 128;
        up = down + 128;

        for (int j = down; j < up; ++j) {
            FFF_temp.push_back(FFF[j]);
            III_temp.push_back(III[j]);
        }
        I_int[i+4] = LogInterpolate(freq,FFF_temp,III_temp);
    }
    
    for (int i = 0; i < 2; ++i) {
        down = (i_mu + i) * 128;
        up = down + 128;

        for (int j = down; j < up; ++j) {
            FFFF_temp.push_back(FFFF[j]);
            IIII_temp.push_back(IIII[j]);
        }
        I_int[i+6] = LogInterpolate(freq,FFFF_temp,IIII_temp);
    }
    chdir(cwd);


    // Perform interpolation to correct mu (cos_theta)
    Q[0] = LogLinear(cos_theta,mu[i_mu],I_int[0],mu[i_mu+1],I_int[1]);
    Q[1] = LogLinear(cos_theta,mu[i_mu],I_int[2],mu[i_mu+1],I_int[3]);
    Q[2] = LogLinear(cos_theta,mu[i_mu],I_int[4],mu[i_mu+1],I_int[5]);
    Q[3] = LogLinear(cos_theta,mu[i_mu],I_int[6],mu[i_mu+1],I_int[7]); 

    // Interpolation to local gravity
    R[0] = LogLinear(Y,Y1,Q[0],Y2,Q[2]);
    R[1] = LogLinear(Y,Y1,Q[1],Y2,Q[3]);

    // Interpolation to temperature
    P = LogLinear(X,X1,R[0],X2,R[1]);

    // Set to zero at small angle
    if (cos_theta < 0.005) P = 0;

    return P;
}

/**************************************************************************************/
/* Helium:                                                                            */
/*           mirroring of hydrogen routines                                           */
/*                                                                                    */
/**************************************************************************************/

//double X, Y, X1, Y1, X2, Y2;
//std::vector<double> F, FF, FFF, FFFF, I, II, III, IIII;

// Read the helium intensity file, set to 4 sets of intensity based on combinations of parameters of nearest 2 logt and 2 lgrav,
// which can be used to peform the four point interpolation.
//const std::vector<double> mu_he {1,0.99995,0.9998,0.99875,0.996195,0.984808,0.965926,0.939693,0.906308,0.866025,0.819152,0.766044,0.707107,0.642787,0.573577,0.499998,0.422622,0.342021,0.258816,0.173652,0.0871556,0.00999616,9.63268e-05};
//const std::vector<double> logt_he {5.50515,5.60206,5.69897,5.79934,5.89763,6,6.11394,6.20412,6.30103,6.39794,6.50515,6.60206,6.69897};
//const std::vector<double> lsgrav_he {13.6021,13.7993,14,14.2041,14.3802,14.6021,14.7993,14.8976};
void Read_NSX(double T, double M, double R){
 
    double delta, lt, lgrav, temp, dump;
    int size_logt(13), size_lsgrav(8), size_mu(23), i_lt, i_lgrav, n_lt, n_lgrav, size_ener(125), size_set, skip_to, skip_two;
    char atmodir[1024], cwd[1024];
    std::vector<double> freq, logt, lsgrav;

    //setting values of lt and lgrav based on input T, M, and R. Also sets to load using first mu value.   
    M = Units::nounits_to_cgs(M, Units::MASS);
    R = Units::nounits_to_cgs(R, Units::LENGTH);
    delta = 1 / sqrt(1 - (2 * Units::G * M / (R * Units::C * Units::C)));
    lgrav = log10(delta * Units::G * M / (R * R));
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));

    //obtain helium atmosphere parameters
    getcwd(cwd, sizeof(cwd));
    sprintf(atmodir,"%s/atmosphere",cwd);
    chdir(atmodir);
    ifstream He_atmo_para;
    He_atmo_para.open("nsx-info.txt");

    if(He_atmo_para.is_open()){
        //logt
        for (int i = 1; i <= size_logt; i++){
            He_atmo_para >> temp;
            logt.push_back(temp); 
        }

        //lgrav
        for (int i = 1; i <= size_lsgrav; i++){
            He_atmo_para >> temp;
            lsgrav.push_back(temp);
        }

        //discarding mu choices
        for (int i = 1; i <= size_mu; i++){
            He_atmo_para >> dump;
        }
    He_atmo_para.close();
    }

    i_lt = Find(lt,logt);
    i_lgrav = Find(lgrav,lsgrav);

    n_lt = i_lt+1;
    n_lgrav = i_lgrav+1;
 

    //Open helium spectra file    
    ifstream He_table1;
    He_table1.open("nsx_He_2.out");


    if(He_table1.is_open()){
        //discarding logt choices
        He_table1 >> size_logt;
        for (int i = 1; i <= size_logt; i++){
            He_table1 >> dump;
         }

        //discarding local gravity choices
        He_table1 >> size_lsgrav;
        for (int i = 1; i <= size_lsgrav; i++){
            He_table1 >> dump;
        }

        //discarding mu choices
        He_table1 >> size_mu;
        for (int i = 1; i <= size_mu; i++){
            He_table1 >> dump;
         }

        //loading energy points
        He_table1 >> size_ener;
        for (int i = 1; i <= size_ener; i++){
            He_table1 >> temp;
            Es.push_back(temp);
            temp = temp * 1E3 * Units::EV / Units::H_PLANCK;
            freq.push_back(temp);
        }

        F = freq;
        FF = freq;
        FFF = freq;
        FFFF = freq;

        //calculate size of each parameter set
        //should be equal to size_ener*size_mu = 125*23 = 2875
        size_set = size_ener*size_mu;

  
        //Calculate how many values to skip until desired parameter set
        skip_to = i_lt*size_lsgrav*size_mu*size_ener+i_lgrav*size_mu*size_ener;
        skip_two = (size_lsgrav-2)*size_mu*size_ener;
  
        //Skipping and loading helium spectrum with chosen parameters/condition.
        int k = 0;
        while (k < skip_to) {
            He_table1 >> dump;
            //cout << "Value " << dump << " has been dumped." << endl;
            k++;
        }
        for (int i = 1; i <= size_set; i++) {
            He_table1 >> temp;
            I.push_back(temp);
        }

        //loading second set
        for (int i = 1; i <= size_set; i++) {
            He_table1 >> temp;
            III.push_back(temp);
        }

        //loading third set
        k = 0;
        while (k < skip_two) {
            He_table1 >> dump;
            k++;
        }
        for (int i = 1; i <= size_set; i++) {
            He_table1 >> temp;
            II.push_back(temp);
        }

        //loading fourth set
        for (int i = 1; i <= size_set; i++) {
            He_table1 >> temp;
            IIII.push_back(temp);
        }
    }else{
        cout << "NSX_file is not found" << endl;
    }
    He_table1.close();
    chdir(cwd);
    
    X = lt;
    Y = lgrav;
    
    X1 = logt[i_lt];
    Y1 = lsgrav[i_lgrav];
    X2 = logt[i_lt + 1];
    Y2 = lsgrav[i_lgrav + 1];  

}

// Calculate the final interpolated intensity
double Helium(double E, double cos_theta){
    double freq, P, size_logt(13), size_lsgrav(8), size_mu(23), temp, dump;
    double I_int[8],Q[4],R[2],mu[23];
    int i_mu(0), down, mid, up;
    char atmodir[1024], cwd[1024];
    std::vector<double> F_temp,FF_temp,FFF_temp,FFFF_temp,I_temp,II_temp,III_temp,IIII_temp,Iv_temp,IIv_temp,IIIv_temp,IIIIv_temp;

    // Convert energy to frequency
    freq = 1E3 * E * Units::EV / Units::H_PLANCK;

    // Read in helium atmosphere parameter choices
    getcwd(cwd, sizeof(cwd));
    sprintf(atmodir,"%s/atmosphere",cwd);
    chdir(atmodir);
    ifstream He_atmo_para;
    He_atmo_para.open("nsx-info.txt");

    if(He_atmo_para.is_open()){
        //discard logt and lsgrav choices
        for (int i = 1; i <= size_logt+size_lsgrav; i++){
            He_atmo_para >> dump;
        }

        //mu
        for (int i = 1; i <= size_mu; i++){
            He_atmo_para >> temp;
            mu[i-1] = temp;
        }
    He_atmo_para.close();
    }

    //finding mu value 
    for (int m = 0; m < size_mu; ++m) {
        if (cos_theta <= mu[m]) {
            i_mu = m;
        }
    }
 
    down = i_mu * 125;
    mid = down + 125;
    up = mid + 125;       

    //Read and interpolate to correct frequency
    for (int j = down; j < mid; ++j) {
        I_temp.push_back(I[j]);
    }
    for (int j = mid; j < up; ++j){
        Iv_temp.push_back(I[j]);
    }
    F_temp = F;
    I_int[0] = LogInterpolate(freq,F_temp,I_temp);
    I_int[1] = LogInterpolate(freq,F_temp,Iv_temp);

    for (int j = down; j < mid; ++j) {
        II_temp.push_back(II[j]);
    }
    for (int j = mid; j < up; ++j){
        IIv_temp.push_back(II[j]);
    }
    FF_temp = F;
    I_int[2] = LogInterpolate(freq,FF_temp,II_temp);
    I_int[3] = LogInterpolate(freq,FF_temp,IIv_temp);   


    for (int j = down; j < mid; ++j) {
        III_temp.push_back(III[j]);
    }
    for (int j = mid; j < up; ++j) {
        IIIv_temp.push_back(III[j]);
    }
    FFF_temp = F;
    I_int[4] = LogInterpolate(freq,FFF_temp,III_temp);
    I_int[5] = LogInterpolate(freq,FFF_temp,IIIv_temp);
    

    for (int j = down; j < mid; ++j) {
        IIII_temp.push_back(IIII[j]);
    }
    for (int j = mid; j < up; ++j) {
        IIIIv_temp.push_back(IIII[j]);
    }
    FFFF_temp = F;
    I_int[6] = LogInterpolate(freq,FFFF_temp,IIII_temp);
    I_int[7] = LogInterpolate(freq,FFFF_temp,IIIIv_temp);
    chdir(cwd);

    // Interpolate to chosen mu
    Q[0] = LogLinear(cos_theta,mu[i_mu],I_int[0],mu[i_mu+1],I_int[1]);
    Q[1] = LogLinear(cos_theta,mu[i_mu],I_int[2],mu[i_mu+1],I_int[3]);
    Q[2] = LogLinear(cos_theta,mu[i_mu],I_int[4],mu[i_mu+1],I_int[5]);
    Q[3] = LogLinear(cos_theta,mu[i_mu],I_int[6],mu[i_mu+1],I_int[7]); 

    // Interpolate to chosen local gravity
    R[0] = LogLinear(Y,Y1,Q[0],Y2,Q[2]);
    R[1] = LogLinear(Y,Y1,Q[1],Y2,Q[3]);

    // Interpolate to chosen temperature
    P = LogLinear(X,X1,R[0],X2,R[1]);

    // Set to zero at small angle
    if (cos_theta < 9.63268e-05) P = 0;

    return P;
}

// Calculate the final interpolated intensity
// Given which energy point to call (0-124) at arbitrary angle, interpolates intensity
// Interpolates to angle, local gravity, and temperature
double Helium2(int E_dex, double cos_theta){
    double P, size_logt(13), size_lsgrav(8), size_mu(23), temp, dump;
    double Q[4],R[2],mu[23];
    int i_mu(0), down, mid, up;
    char atmodir[1024], cwd[1024];
    std::vector<double> I_temp,Iv_temp,II_temp,IIv_temp,III_temp,IIIv_temp,IIII_temp,IIIIv_temp;

    // Read in helium atmosphere parameter choices
    getcwd(cwd, sizeof(cwd));
    sprintf(atmodir,"%s/atmosphere",cwd);
    chdir(atmodir);
    ifstream He_atmo_para;
    He_atmo_para.open("nsx-info.txt");

    if(He_atmo_para.is_open()){
        //discard logt and lsgrav choices
        for (int i = 1; i <= size_logt+size_lsgrav; i++){
            He_atmo_para >> dump;
        }

        //mu
        for (int i = 1; i <= size_mu; i++){
            He_atmo_para >> temp;
            mu[i-1] = temp;
        }
    He_atmo_para.close();
    }

    //finding mu value 
    for (int m = 0; m < size_mu; ++m) {
        if (cos_theta <= mu[m]) {
            i_mu = m;
        }
    }

    down = i_mu * 125;
    mid = down + 125;
    up = mid + 125;       

    //Read and interpolate to correct frequency
    for (int j = down; j < mid; ++j) {
        I_temp.push_back(I[j]);
    }
    for (int j = mid; j < up; ++j){
        Iv_temp.push_back(I[j]);
    } 
    for (int j = down; j < mid; ++j) {
        II_temp.push_back(II[j]);
    }
    for (int j = mid; j < up; ++j){
        IIv_temp.push_back(II[j]);
    }  
    for (int j = down; j < mid; ++j) {
        III_temp.push_back(III[j]);
    }
    for (int j = mid; j < up; ++j) {
        IIIv_temp.push_back(III[j]);
    }
    for (int j = down; j < mid; ++j) {
        IIII_temp.push_back(IIII[j]);
    }
    for (int j = mid; j < up; ++j) {
        IIIIv_temp.push_back(IIII[j]);
    }
    chdir(cwd);

    // Interpolate to chosen mu
    Q[0] = LogLinear(cos_theta,mu[i_mu],I_temp[E_dex],mu[i_mu+1],Iv_temp[E_dex]);
    Q[1] = LogLinear(cos_theta,mu[i_mu],II_temp[E_dex],mu[i_mu+1],IIv_temp[E_dex]);
    Q[2] = LogLinear(cos_theta,mu[i_mu],III_temp[E_dex],mu[i_mu+1],IIIv_temp[E_dex]);
    Q[3] = LogLinear(cos_theta,mu[i_mu],IIII_temp[E_dex],mu[i_mu+1],IIIIv_temp[E_dex]); 

    // Interpolate to chosen local gravity
    R[0] = LogLinear(Y,Y1,Q[0],Y2,Q[2]);
    R[1] = LogLinear(Y,Y1,Q[1],Y2,Q[3]);

    // Interpolate to chosen temperature
    P = LogLinear(X,X1,R[0],X2,R[1]);

    // Set to zero at small angle
    if (cos_theta < 9.63268e-05) P = 0;

    return P;
}


/**************************************************************************************/
/* Blackbody:                                                                         */
/*           computes the monochromatic blackbody flux in units of erg/cm^2			  */
/*																					  */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E = monochromatic energy in keV * redshift / eta                             */
/**************************************************************************************/
double BlackBody( double T, double E ) {   // Blackbody flux in units of erg/cm^2
    return ( 2.0e9 / pow(Units::C * Units::H_PLANCK, 2) * pow(E * Units::EV, 3) / (exp(E/T) - 1) ); // shouldn't it have a pi?
    // the e9 is to switch E from keV to eV; Units::EV gets it from eV to erg, since it's first computed in erg units.
    // the switch from erg units to photon count units happens above just after this is called.
} // end Blackbody


/**************************************************************************************/
/* Line:                                                                         */
/*           computes the monochromatic blackbody flux in units of erg/cm^2*/
/*           for a small range of emission frame energies. */
/*																					  */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E = monochromatic energy in keV * redshift / eta                             */
/**************************************************************************************/
double Line( double T, double E , double E1, double E2) {   // Blackbody flux in units of erg/cm^2
  // If  E1 <= E <= E2 then blackbody flux. Zero otherwise  
  // E, E1, E2 are all as locally measured in the rest frame of the star
  if ( E <= E2 && E >= E1)
    return ( 2.0e9 / pow(Units::C * Units::H_PLANCK, 2) * pow(E * Units::EV, 3) / (exp(E/T) - 1) ); 

  else
    return 0.0;
   
} // end line

/**************************************************************************************/
/* EnergyBandFlux:                                                                    */
/*                computes the blackbody flux in units of erg/cm^2 using trapezoidal  */
/*                rule for approximating an integral                                  */
/*                variant of Bradt equation 6.6                                       */
/*                T, E1, E2 put into eV in this routine                               */
/*																					  */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E1 = lower bound of energy band in keV * redshift / eta                      */
/*       E2 = upper bound of energy band in keV * redshift / eta                      */
/*       L1 = lower limit of emitted energy in star's frame                           */
/*       L2 = lower limit of emitted energy in star's frame                           */
/**************************************************************************************/
double LineBandFlux( double T, double E1, double E2, double L1, double L2 ) {
  T *= 1e3; // from keV to eV
  // x = E / T
  E1 *= 1e3; // from keV to eV
  E2 *= 1e3; // from keV to eV
  L1 *= 1e3;
  L2 *= 1e3;

  /********************************************/
  /* VARIABLE DECLARATIONS FOR LINEBandFlux */
  /********************************************/
	
	// a, b, x, n, h as defined by Mathematical Handbook eqn 15.16 (Trapezoidal rule to approximate definite integrals)
	double a = E1 / T;          // lower bound of integration
	double b = E2 / T;          // upper bound of integration
	double current_x(0.0);      // current value of x, at which we are evaluating the integrand; x = E / T; unitless
	unsigned int current_n(0);  // current step
	unsigned int n_steps(400); // total number of steps
	double h = (b - a) / n_steps;     // step amount for numerical integration; the size of each step
	double integral_constants = 2.0 * pow(T*Units::EV,3) / pow(Units::C,2) / pow(Units::H_PLANCK,3); // what comes before the integral when calculating flux using Bradt eqn 6.6 (in units of photons/cm^2/s)
	double flux(0.0);           // the resultant energy flux density; Bradt eqn 6.17
	
	// begin trapezoidal rule
	current_x = a + h * current_n;
	if ( L1/T <= current_x && current_x <= L2/T)
	  flux = Bradt_flux_integrand(current_x);

	for ( current_n = 1; current_n < n_steps-1; current_n++ ) {
		current_x = a + h * current_n;
		if ( L1/T <= current_x && current_x <= L2/T)
		  flux += 2.0 * Bradt_flux_integrand(current_x);
	}
	
	current_x = a + h * current_n;
	if ( L1/T <= current_x && current_x <= L2/T)
	  flux += Bradt_flux_integrand(current_x);
	
	flux *= h/2.0;	
	// end trapezoidal rule; numerical integration complete!
	
	flux *= integral_constants;

	return flux;
} // end LINEBandFlux

/**************************************************************************************/
/* EnergyBandFlux:                                                                    */
/*                computes the blackbody flux in units of erg/cm^2 using trapezoidal  */
/*                rule for approximating an integral                                  */
/*                variant of Bradt equation 6.6                                       */
/*                T, E1, E2 put into eV in this routine                               */
/*																					  */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E1 = lower bound of energy band in keV * redshift / eta                      */
/*       E2 = upper bound of energy band in keV * redshift / eta                      */
/**************************************************************************************/
double EnergyBandFlux( double T, double E1, double E2 ) {
	T *= 1e3; // from keV to eV
	// x = E / T
	E1 *= 1e3; // from keV to eV
	E2 *= 1e3; // from keV to eV
	
	/********************************************/
   	/* VARIABLE DECLARATIONS FOR EnergyBandFlux */
    /********************************************/
	
	// a, b, x, n, h as defined by Mathematical Handbook eqn 15.16 (Trapezoidal rule to approximate definite integrals)
	double a = E1 / T;          // lower bound of integration
	double b = E2 / T;          // upper bound of integration
	double current_x(0.0);      // current value of x, at which we are evaluating the integrand; x = E / T; unitless
	unsigned int current_n(0);  // current step
	unsigned int n_steps(3000); // total number of steps
	double h = (b - a) / n_steps;     // step amount for numerical integration; the size of each step
	double integral_constants = 2.0 * pow(T*Units::EV,3) / pow(Units::C,2) / pow(Units::H_PLANCK,3); // what comes before the integral when calculating flux using Bradt eqn 6.6 (in units of photons/cm^2/s)
	double flux(0.0);           // the resultant energy flux density; Bradt eqn 6.17
	
	// begin trapezoidal rule
	current_x = a + h * current_n;
	flux = Bradt_flux_integrand(current_x);

	for ( current_n = 1; current_n < n_steps-1; current_n++ ) {
		current_x = a + h * current_n;
		flux += 2.0 * Bradt_flux_integrand(current_x);
	}
	
	current_x = a + h * current_n;
	flux += Bradt_flux_integrand(current_x);
	
	flux *= h/2.0;	
	// end trapezoidal rule; numerical integration complete!
	
	flux *= integral_constants;

	return flux;
} // end EnergyBandFlux

/**************************************************************************************/
/* AtmosEBandFlux:                                                                    */
/*                computes the flux of atmosphere models in units of counts/s/cm^2    */
/*                using Simpson's rule for approximating an integral                  */
/*                See Numerical Recipes eq 4.1.4                                      */
/*                requires usage of Hydrogen/Helium routines                          */
/*                lT, lgrav are already set by which intensity files to load          */
/*                cos_theta is used in Hydrogen/Helium routines                       */
/*                E1, E2 are in keV, will be converted                                */
/*                                                                                    */
/* pass: model = atmospheric model (0 = Hydrogen, 1 = Helium)                         */
/*       cos_theta = angle at one particular time bin                                 */
/*       E1 = lower bound of energy band in keV * redshift / eta                      */
/*       E2 = upper bound of energy band in keV * redshift / eta                      */
/**************************************************************************************/
double AtmosEBandFlux( unsigned int model, double cos_theta, double E1, double E2 ) {

    /********************************************/
    /* VARIABLE DECLARATIONS FOR AtmosEBandFlux */
    /********************************************/
    
    int n_steps(3000);       // total number of steps
    double step_size;           // step size in energy scale
    double current_e;           // central energy of current step
    double e_l, e_m, e_u;       // lower, middle, and upper energy values of Simpson's rule
    double current_n;           // integrated flux in current step
    double flux(0.0);           // total integrated flux

    
    step_size = (E2-E1)/n_steps;

    if (model == 3){ //Hydrogen
        for (int j = 1; j <= n_steps; j++){
            
            current_e = step_size*(j-1+1/2) + E1;
            e_l = (current_e-step_size/2);
            e_m = (current_e);
            e_u = (current_e+step_size/2);
            current_n = step_size * (1/3*Hydrogen(e_l,cos_theta)/e_l + 4/3*Hydrogen(e_m,cos_theta)/e_m + 1/3*Hydrogen(e_u,cos_theta)/e_u);
            flux += current_n;
            
        }
    }
    if (model == 4){ //Helium
        for (int j = 1; j <= n_steps; j++){
            current_e = step_size*(j-1+1/2) + E1;
            e_l = (current_e-step_size/2);
            e_m = (current_e);
            e_u = (current_e+step_size/2);
            current_n = step_size * (1/3*Helium(e_l,cos_theta)/e_l + 4/3*Helium(e_m,cos_theta)/e_m + 1/3*Helium(e_u,cos_theta)/e_u);
            flux += current_n;
        }
    }

    flux = flux/Units::H_PLANCK;
    return flux;
}



/**************************************************************************************/
/* AtmosEBandFlux2: (ONLY FOR HELIUM FOR NOW)                                         */
/*                computes the flux of atmosphere models in units of counts/s/cm^2    */
/*                using Simpson's rule for approximating an integral                  */
/*                See Numerical Recipes eq 4.1.3-5, 4.1.13, 4.1.14                    */
/*                requires usage of Helium2 routine                                   */
/*                lt, lgrav are already set by which intensity files to load          */
/*                cos_theta is used in Hydrogen/Helium routines                       */
/*                E1, E2 are in keV, will be converted for internal processing        */
/*                E_dex needs to be found in this routine and passed into Helium2     */
/*                                                                                    */
/* pass: model = beaming model (helium = 4)                                           */
/*       cos_theta = angle at one particular time bin                                 */
/*       E1 = lower bound of energy band in keV * redshift / eta                      */
/*       E2 = upper bound of energy band in keV * redshift / eta                      */
/**************************************************************************************/
double AtmosEBandFlux2( unsigned int model, double cos_theta, double E1, double E2 ) {

    /********************************************/
    /* VARIABLE DECLARATIONS FOR AtmosEBandFlux */
    /********************************************/
    
    int e1_dex(0);              // largest energy index that's smaller than E1
    int e2_dex(0);              // largest energy index that's smaller than E2
    int n_steps(0);             // number of energy points within band
    int ener_size(125);         // number of energy choices
    //int e_dex;                  // energy index of current step
    //double current_e;           // central energy of current step
    //double current_n;           // integrated flux in current step
    double flux(0.0);           // total integrated flux

    for (int m = 0; m < ener_size; m++) {
        if (E1 >= Es[m]) {
            e1_dex = m;
        }
        if (E2 >= Es[m]) {
            e2_dex = m;
        }
    }


    //calculate number of energy points within band  
    n_steps = e2_dex - e1_dex;

    if (n_steps == 0){ // zero energy points within bandwidth: (4.1.3) one trapzoid
        cout << "0 steps" << endl;
        flux = (E2 - E1) / 2 * (Helium(E1,cos_theta) / E1 + Helium(E2,cos_theta) / E2);
    }
    if (n_steps == 1){ // one energy points within bandwidth: (4.1.3) two trapzoids
        cout << "1 steps" << endl;
        int e_dex = e1_dex+1; // index of the energy point
        double e_m = F[e_dex] * Units::H_PLANCK / Units::EV / 1E3; // energy point in keV
        flux = (e_m - E1) / 2 * (Helium(E1,cos_theta) / E1 + Helium2(e_dex,cos_theta) / e_m);  // first trapezoid
        flux += (E2 - e_m) / 2 * (Helium2(e_dex,cos_theta) / e_m + Helium(E2,cos_theta) / E2); // second trapezoid
    }
    if (n_steps == 2){ // two energy points within bandwidth: (4.1.3) three trapzoids
        cout << "2 steps" << endl;
        int el_dex = e1_dex+1; // index of the first energy point within bandwidth
        int eu_dex = e2_dex;   // index of the last energy point within bandwidth
        double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
        double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
        double h = log(e_u)-log(e_l); // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (Helium(E1,cos_theta) / E1 + Helium2(el_dex,cos_theta) / e_l);  // first trapezoid
        flux += h * (Helium2(el_dex,cos_theta) + Helium2(eu_dex,cos_theta)) / 2;                // middle trapezoid, exact
        flux += (E2 - e_u) / 2 * (Helium2(eu_dex,cos_theta) / e_u + Helium(E2,cos_theta) / E2); // last trapezoid
    }
    if (n_steps == 3){ // three energy points within bandwidth: (4.1.3, 4.1.4) Simpson's + two trapezoids
        cout << "3 steps" << endl;
        int el_dex = e1_dex+1; // index of the first energy point within bandwidth
        int em_dex = e1_dex+2; // index of the middle energy point within bandwidth
        int eu_dex = e2_dex;   // index of the last energy point within bandwidth
        double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
        double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
        double h = (log(e_u)-log(e_l)) / 2; // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (Helium(E1,cos_theta) / E1 + Helium2(el_dex,cos_theta) / e_l);  // first trapezoid
        flux += h * (Helium2(el_dex,cos_theta) + Helium2(em_dex,cos_theta) * 4 + Helium2(eu_dex,cos_theta)) / 3; // Simpson's, exact
        flux += (E2 - e_u) / 2 * (Helium2(eu_dex,cos_theta) / e_u + Helium(E2,cos_theta) / E2); // second trapezoid
    }
    if (n_steps == 4){ // four energy points within bandwidth: (4.1.3, 4.1.5) 8/3 Simpson's + two trapezoids
        cout << "4 steps" << endl;
        int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
        int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
        int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
        int eu_dex = e2_dex;    // index of the last energy point within bandwidth
        double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
        double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
        double h = (log(e_u)-log(e_l)) / 3; // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (Helium(E1,cos_theta) / E1 + Helium2(el_dex,cos_theta) / e_l);  // first trapezoid
        flux += h * (Helium2(el_dex,cos_theta) + Helium2(em1_dex,cos_theta) * 3 + Helium2(em2_dex,cos_theta) * 3 + Helium2(eu_dex,cos_theta)) * 3 / 8; // Simpson's 3/8, exact             
        flux += (E2 - e_u) / 2 * (Helium2(eu_dex,cos_theta) / e_u + Helium(E2,cos_theta) / E2); // second trapezoid
    }
    if (n_steps == 5){ // five energy points within bandwidth: (4.1.3, 4.1.13) Simpson's "4,2" + two trapezoids  
        cout << "5 steps" << endl;
        int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
        int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
        int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
        int em3_dex = e1_dex+4; // index of the fourth energy point within bandwidth        
        int eu_dex = e2_dex;    // index of the last energy point within bandwidth
        double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
        double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
        double h = (log(e_u)-log(e_l)) / 4; // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (Helium(E1,cos_theta) / E1 + Helium2(el_dex,cos_theta) / e_l);  // first trapezoid
        flux += h * (Helium2(el_dex,cos_theta) + Helium2(em1_dex,cos_theta) * 4 + Helium2(em2_dex,cos_theta) * 2 + Helium2(em3_dex,cos_theta) * 4 + Helium2(eu_dex,cos_theta)) / 3; // Simpson's "4,2", exact             
        flux += (E2 - e_u) / 2 * (Helium2(eu_dex,cos_theta) / e_u + Helium(E2,cos_theta) / E2); // second trapezoid
    }
    if (n_steps >= 6){ // six or more energy points within bandwidth: (4.1.3, 4.1.14) Simpson's cubic + two trapezoids
        cout << "6 steps" << endl;
        int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
        int el1_dex = e1_dex+2; // index of the second energy point within bandwidth
        int el2_dex = e1_dex+3; // index of the third energy point within bandwidth        
        int eu_dex = e2_dex;    // index of the last energy point within bandwidth
        int eu1_dex = e2_dex-1; // index of the second last energy point within bandwidth
        int eu2_dex = e2_dex-2; // index of the third last energy point within bandwidth
        double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
        double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
        double h = (log(e_u)-log(e_l)) / (n_steps-1); // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (Helium(E1,cos_theta) / E1 + Helium2(el_dex,cos_theta) / e_l);  // first trapezoid
        // Simpson's cubic, exact
        flux += h * (Helium2(el_dex,cos_theta) * 9 + Helium2(el1_dex,cos_theta) * 28 + Helium2(el2_dex,cos_theta) * 23) / 24; // first three coefficients
        for (int m = 1; m <= n_steps-6; m++) flux += h * (Helium2(el_dex+m+3,cos_theta)); // middle coefficients
        flux += h * (Helium2(eu2_dex,cos_theta) * 23 + Helium2(eu1_dex,cos_theta) * 28 + Helium2(eu_dex,cos_theta) * 9) / 24; // last three coefficients
        // end Simpson's cubic
        flux += (E2 - e_u) / 2 * (Helium2(eu_dex,cos_theta) / e_u + Helium(E2,cos_theta) / E2); // second trapezoid                
    }


    
    flux = flux/Units::H_PLANCK;
    return flux;
}


/**************************************************************************************/
/* Bradt_flux_integrand:                                                              */
/*                      integrand of Bradt eqn 6.6 when integrating over nu, modified */
/*                      so the exponent is 2 not 3, so that it comes out as photon    */
/*                      number flux, instead of erg flux                              */
/*																					  */
/* pass: x = current_x from above routine                                             */
/**************************************************************************************/
double Bradt_flux_integrand( double x ) {
	return ( pow(x,2) / (exp(x) - 1) );  // 2 (not 3) for photon number flux
} // end Bradt_flux_integrand

/**************************************************************************************/
/* Gray:																			  */
/*		computes and returns the limb darkening factors for a Gray electron -         */
/*		scattering atmosphere (Hopf function) using values from Chandrasekhar's       */
/*		"Radiative Transfer", table 13												  */
/*																					  */
/* pass: cosine = curve.cosbeta[i] * curve.eta[i]                                     */
/**************************************************************************************/
double Gray( double cosine ) {
	
	/**********************************/
   	/* VARIABLE DECLARATIONS FOR Gray */
    /**********************************/
	
    double F[11],   // limb darkening values for electron-scattering atmosphere
           mu[11];  // a table of cos(theta) values
    int i;          // loop variable

    for ( i = 0; i <= 10; i++ ) {
        mu[i] = 0.0 + i*0.1;
    }
    /*
    // Value from Mihalas's "Stellar Atmospheres"
    F[0] = 0.4330;
    F[1] = 0.5401;
    F[2] = 0.6280;
    F[3] = 0.7112;
    F[4] = 0.7921;
    F[5] = 0.8716;
    F[6] = 0.9501;
    F[7] = 1.0280;
    F[8] = 1.1053;
    F[9] = 1.1824;
    F[10] = 1.2591;
    */
  
    // Values from Chandrasekhar's "Radiative Transfer", Table 13
    F[0] = 1.0;
    F[1] = 1.2647;
    F[2] = 1.48009;
    F[3] = 1.68355;
    F[4] = 1.88105;
    F[5] = 2.07496;
    F[6] = 2.2665;
    F[7] = 2.45639;
    F[8] = 2.64503;
    F[9] = 2.83274;
    F[10] = 3.01973;
  
    for ( i = 0; i <= 10; i++ ) {
        F[i] *= (0.5/1.194);
    }

    int index(0);
    //i=1;

    if (cosine == 1.0)
    	index = 10;
    else {
        while ( cosine >= mu[index]) 
        	index++;
    }
    // cosine lies in the range   mu[index-1] <= cosine <= mu[index]

    if (cosine == mu[index-1]) 
    	return F[index-1];
    else {
        if (cosine == mu[index]) 
        	return F[index];
        else {   
            return ( F[index-1] + (cosine - mu[index-1]) * (F[index] - F[index-1])/(mu[index]-mu[index-1]) );
        }
    }
} // end Gray
