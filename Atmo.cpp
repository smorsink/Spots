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
#include "interp.h"
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


    std::ofstream ttt;

    double 
      mass_over_r,
           temperature,        // Temperature of the spot, in keV
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
        
    //  std::vector< double > totflux(MAX_NUMBINS, 0.0); // integrated flux

    /*********************/
    /* SETTING THINGS UP */
    /*********************/
    
    // One monochromatic energy, hardwired value, in keV
    //    E_mono = 1.0;
    //std::cout << "starting computeangles" << std::endl;
    curve = (*angles);

    mass_over_r = curve.para.mass_over_r;
    temperature = curve.para.temperature;       // T in keV 

    numbins = curve.numbins;
    numbands = curve.numbands;
    E_band_lower_1 = curve.para.E_band_lower_1;     // in keV
    E_band_upper_1 = curve.para.E_band_upper_1;     // in keV
    E_band_lower_2 = curve.para.E_band_lower_2;     // in keV
    E_band_upper_2 = curve.para.E_band_upper_2;     // in keV
    E0 = curve.para.E0;
    E1 = curve.para.L1;
    E2 = curve.para.L2;
    DeltaE = curve.para.DeltaE;
    //cout << curve.mccinte[0] << endl;
    //std::cout << "beaming model = " << curve.flags.beaming_model << std::endl;
   
    redshift = 1.0 / sqrt( 1 - 2.0 * mass_over_r);

    double M = Units::nounits_to_cgs(curve.para.mass, Units::MASS);
    double R = Units::nounits_to_cgs(curve.para.req, Units::LENGTH);
    double delta = 1 / sqrt(1 - 2.0*curve.para.mass/curve.para.req);
    double cos_theta = cos(curve.para.theta);
    double obl_approx = 1 + (-0.791 + 0.776 * curve.para.mass/curve.para.req) * curve.para.omega_bar_sq * (1.0-pow(cos_theta,2)) 
		  + (1.138 - 1.431 * curve.para.mass/curve.para.req) * curve.para.omega_bar_sq *  pow(cos_theta,2);
    std::cout << "obl_approx = " << obl_approx << std::endl;
    double lgrav = log10(delta * Units::G * M / R / R * obl_approx);
    std::cout << "log(g) = " << lgrav << std::endl;
    std::cout << curve.flags.beaming_model << std::endl;



    bolo = 2.0e9 * 2.404 * Units::C * pow(temperature * Units::EV/(Units::H_PLANCK * Units::C) , 3); // use this one! probably!
    // the e9 in the beginning is for changing T^3 from keV to eV
    // 2.404 comes from evaluating Bradt equation 6.17 (modified, for photon number count units), using the Riemann zeta function for z=3

    //std::cout << " numbins = " << numbins << std::endl;

    // numbins=0;

    for ( unsigned int i(0); i < numbins; i++ ) { // Compute flux for each phase bin

      //std::cout << "i = " << i << std::endl;
      //std::cout << "beaming model = " << curve.flags.beaming_model << std::endl;


      if ( curve.dOmega_s[i] != 0.0 ) {

	// Default value is gray = 1
	gray = 1.0; // beaming_model == 0

	// Calculate Beaming for Modified Blackbodies

	if (curve.flags.beaming_model == 1)
	    gray = Gray(curve.cosbeta[i]*curve.eta[i]); // Chandrasekhar limb-darkening	  
	if (curve.flags.beaming_model == 5)
	    gray = pow( curve.cosbeta[i]*curve.eta[i], 2.0);
	if (curve.flags.beaming_model == 6)
	    gray = 1.0 - pow( curve.cosbeta[i]*curve.eta[i], 2.0);
	if (curve.flags.beaming_model == 7) // Hopf Function
	    gray = 0.42822+0.92236*curve.cosbeta[i]*curve.eta[i]-0.085751*pow(curve.cosbeta[i]*curve.eta[i],2);
	     
	/*******************************************************************/
	/* COMPUTING LIGHT CURVE FOR MONOCHROMATIC ENERGY, p = 0           */
	/*      First computes in [erg/(s cm^2 Hz), converts to            */
	/*      photons/(s cm^2 keV)                                       */
	/*******************************************************************/


	if (curve.flags.spectral_model == 0){ // Monochromatic Observation of a modified blackbody
	  	double E_diff = (E_band_upper_1 - E_band_lower_1)/numbands;

		if ( curve.flags.beaming_model == 0 || curve.flags.beaming_model == 1 || curve.flags.beaming_model == 6 || curve.flags.beaming_model == 7){

		  curve.f[0][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * BlackBody(temperature,E0*redshift/curve.eta[i]); 
		  curve.f[0][i] *= (1.0 / ( E0 * Units::H_PLANCK )); // Units: photons/(s cm^2 keV)
           
	  	}

	  	if (curve.flags.beaming_model == 2){ // Funny Line Emission, not calculated
	  	}

	  	if (curve.flags.beaming_model == 3){ // NSATMOS Hydrogen Atmosphere
		  for (unsigned int p = 0; p<numbands; p++){
		    curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * Hydrogen((E_band_lower_1+p*E_diff)*redshift/curve.eta[i], curve.cosbeta[i]*curve.eta[i]);
		    curve.f[p][i] *= (1.0 / ( (E_band_lower_1+(p+0.5)*E_diff) * Units::H_PLANCK ));
		  }
	  	}

	  	if (curve.flags.beaming_model == 4){ // NSX Helium Atmosphere
		  for (unsigned int p = 0; p<numbands; p++){
		    curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * Helium((E_band_lower_1+p*E_diff)*redshift/curve.eta[i], curve.cosbeta[i]*curve.eta[i]);
		    curve.f[p][i] *= (1.0 / ( (E_band_lower_1+(p+0.5)*E_diff) * Units::H_PLANCK ));
		  }
	  	}

	  	if (curve.flags.beaming_model == 5){ // NSXH Atmosphere
		  for (unsigned int p = 0; p<numbands; p++){
		    curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * NSXH((E_band_lower_1+(p+0.5)*E_diff)*redshift/curve.eta[i], curve.cosbeta[i]*curve.eta[i]);
		    curve.f[p][i] *= (1.0 / ( (E_band_lower_1+(p+0.5)*E_diff) * Units::H_PLANCK ));
	    	
		  }
	  	}

	  	if (curve.flags.beaming_model == 8){ // *slavko* McPHAC
		  for (unsigned int p = 0; p<numbands; p++){
		    curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * McPHAC((E_band_lower_1+(p+0.5)*E_diff)*redshift/curve.eta[i], curve.cosbeta[i]*curve.eta[i]);
		    curve.f[p][i] *= (1.0 / ( (E_band_lower_1+(p+0.5)*E_diff) * Units::H_PLANCK ));
		    curve.f[p][i] *= E_diff;
		  }
	  	}

	  	if (curve.flags.beaming_model == 9){ // *new* NSX Helium Atmosphere
		  for (unsigned int p = 0; p<numbands; p++){
		    curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * NSXHe((E_band_lower_1+(p+0.5)*E_diff)*redshift/curve.eta[i], curve.cosbeta[i]*curve.eta[i]);
		    curve.f[p][i] *= (1.0 / ( (E_band_lower_1+(p+0.5)*E_diff) * Units::H_PLANCK ));
		  }
	  	}
	 	
	 
	        	
	  	if (curve.flags.beaming_model == 10){ // *cole* McPHACC3
		  for (unsigned int p = 0; p<numbands; p++){
		    //cout << p << endl;
		    curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * McPHACC3new((E_band_lower_1+(p+0.5)*E_diff)*redshift/curve.eta[i], curve.cosbeta[i]*curve.eta[i], curve.para.temperature, lgrav, curve);
		    curve.f[p][i] *= (1.0 / ( (E_band_lower_1+(p+0.5)*E_diff) * Units::H_PLANCK ));
		    
curve.f[p][i] *= E_diff; // Fake Integration
		   
		  }
	  	}


		if (curve.flags.beaming_model == 11){ // New NSX-H
		  //cout << "Calling NSXHnew!" << endl;
		  for (unsigned int p = 0; p<numbands; p++){
		    //cout << p << endl;
		    curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * NSXHnew((E_band_lower_1+(p+0.5)*E_diff)*redshift/curve.eta[i], curve.cosbeta[i]*curve.eta[i], curve.para.temperature, curve.para.mass, curve.para.radius, curve);
		    curve.f[p][i] *= (1.0 / ( (E_band_lower_1+(p+0.5)*E_diff) * Units::H_PLANCK ));
		    //if (isnan(curve.f[p][i])) cout << "curve at " << p << " " << " is nan!" << endl;
		  }
	  	}
	
	  	
	  	
	  	


	} // Spectral_model == 0 

	/***********************/
	/* FUNNY LINE EMISSION */
	/***********************/
	if (curve.flags.spectral_model == 1){ // Funny Line Emission for NICER

	  for (unsigned int p=0; p<numbands; p++){
	    
	    E_obs = E0 + p*DeltaE;
	    curve.f[p][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * LineBandFlux(temperature, (E_obs-0.5*DeltaE)*redshift/curve.eta[i], (E_obs+0.5*DeltaE)*redshift/curve.eta[i], E1, E2); // Units: photon/(s cm^2)


	  }
	}

	/*******************************************************************/
	/* COMPUTING BLACKBODY LIGHT CURVE FOR INTEGRATED FLUX             */
	/* Units: photons/(cm^2 s)                                         */
	/*******************************************************************/

	if (curve.flags.spectral_model == 2){ // Integrated Flux for Modified Blackbodies
	  double E_diff = (E_band_upper_1 - E_band_lower_1)/numbands;

	  for (unsigned int p = 0; p<numbands; p++){
           
	    curve.f[p][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) 
	    	 * EnergyBandFlux(temperature, (E_band_lower_1+p*E_diff)*redshift/curve.eta[i], (E_band_lower_1+(p+1)*E_diff)*redshift/curve.eta[i]); // Units: photon/(s cm^2)

            }          
        }
      

	if (curve.flags.spectral_model == 3){ // *exactly* Integrated Flux of Energy Bands
	  	double E_diff = (E_band_upper_1 - E_band_lower_1)/numbands;
	  	for (unsigned int p = 0; p<numbands; p++){
            if (curve.flags.beaming_model == 3 || curve.flags.beaming_model == 4 || curve.flags.beaming_model == 5 || curve.flags.beaming_model == 8 || curve.flags.beaming_model == 9){
	      		curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * AtmosEBandFlux2(curve.flags.beaming_model, curve.cosbeta[i]*curve.eta[i], (E_band_lower_1+p*E_diff)*redshift/curve.eta[i], (E_band_lower_1+(p+1)*E_diff)*redshift/curve.eta[i]); // Units: photon/(s cm^2)        	      		
	      	}
	    if (curve.flags.beaming_model == 10){ // Cole's McPhac File
	      	       	      		
	      curve.f[p][i] = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) * AtmosEBandFlux3new(curve.flags.beaming_model, curve.cosbeta[i]*curve.eta[i], curve.para.temperature,lgrav, (E_band_lower_1+p*E_diff)*redshift/curve.eta[i], (E_band_lower_1+(p+1)*E_diff)*redshift/curve.eta[i], curve); // Units: photon/(s cm^2)        
	      	}
	   
	  	}
	}
	
      	}
      	else { // if curve.dOmega_s[i] == 0.0
			for ( unsigned int p(0); p < numbands; p++) {
	  			curve.f[p][i] = 0.0;
			}
      	}
    //cout << "phase " << i << " done." << endl;
    } // ending the for(i) loop
    

    //std::cout << "McPhac (Old) Flux = " << McPHAC(1.0,0.5) << std::endl;
    

    //    std::cout << "McPhac Flux = " << 
    //McPHACC3new(1.0,0.5, curve.para.temperature, curve.para.mass, curve.para.radius, curve)
    //	      << std::endl;


    //std::cout << "NSXH Flux = " <<
    //NSXH(1.0, 0.5)
    //	      << std::endl;
    
    



	
    /***********************************************************/
    /* DEALING WITH THE TIME DELAYS, REBINNING THE LIGHT CURVE */
    /* This is where the jumpy problems tend to be.            */
    /***********************************************************/

    curve.flags.ignore_time_delays = false;
    		
    if ( !curve.flags.ignore_time_delays ) { // if we are not ignoring the time delays        
      int k(0), j(0); // index placeholders; approximately, k is i+1 and j is i-1
      // but time delays mean that j isn't always i-1
      // used in the linear interpolation
        
      //std::cout << "Taking care of time delays!" << std::endl;


      int ecl1(0), ecl2(0), j1, j2, k1, k2;
      double tt1,tt2, ttt1, ttt2;

      if (curve.eclipse){

	// If curve has an eclipse, find out where 
	// Find the first eclipsed bin; set ec1=i;
	for (unsigned int i(0);i<numbins;i++){
	  if (curve.dOmega_s[i]==0.0){
	    ecl1=i;
	    break;
	  }
	}
	j1 = ecl1-2;
	j2 = ecl1-1;
	
	if (ecl1==0){
	  j1+=numbins;
	  j2+=numbins;
	}
	if (ecl1==1){
	  j1+=numbins;
	}
	tt1 = curve.t_o[j1];
	if (ecl1<2)
	  tt1 +=-1.0;
	tt2 = curve.t_o[j2];
	if (ecl1==0)
	  tt2+=-1.0;

	// After the eclipse, find the first bin with non-zero flux
	for (unsigned int i(ecl1);i<numbins;i++){
	  if (curve.dOmega_s[i]>0.0){
	    ecl2=i;
	    break;
	  }
	}
	//std::cout << "eclipse ends at bin" << ecl2 << std::endl;
	k1=ecl2;
	k2=ecl2+1;
	if (ecl2==numbins-1){
	  k2+=-numbins;
	}
	ttt1=curve.t_o[k1];
	ttt2=curve.t_o[k2];
	if (ecl2==numbins-1)
	  ttt2+=1.0;


      } // End of Eclipse finder

      /********************************/
      /* LOOP THROUGH THE LIGHTCURVES */
      /********************************/
       
      for (unsigned int p(0); p < numbands; p++) {
      // for (unsigned int p(1); p < 2; p++) {

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
	  
	  double te1, te2, slope1, slope2;

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

	  //std::cout << "maximum takes place at time=" << tx << std::endl;

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

	    slope1 = (curve.f[p][j1]-curve.f[p][j2])/(tt1-tt2);
	
	    if (slope1 != 0.0)
	      te1  = - (curve.f[p][j2])/slope1 + tt2 ;
	    else
	      te1 = tt1;

	    slope2 = (curve.f[p][k1]-curve.f[p][k2])/(ttt1-ttt2);
	    te2 = ttt1 - curve.f[p][k1]/slope2;

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
   			
	  //curve.minFlux[p] = minimum;
	  //curve.maxFlux[p] = maximum;
	  //curve.pulseFraction[p] = (curve.maxFlux[p] - curve.minFlux[p]) / (curve.maxFlux[p] + curve.minFlux[p]);

	  //curve.asym[p] = (tmin - tx) - 0.5;
	  //if (curve.asym[p] < 0.0) curve.asym[p]+=1.0;


	  // Initializing totflux
	  //for ( unsigned int i(0); i < numbins; i++ )
	  //totflux.at(i) = 0.0;
			
		           		
	  /**************************************************************/
	  /* ADDING FLUXES FROM ALL THE PHASE BINS AND OTHER FUN THINGS */
	  /**************************************************************/

	  for ( unsigned int i(0); i < numbins; i++ ) {  // for-i-loop, looping through the phase bins
            
	    k = i + 1;
	    j = i;		 
	    if ( k == static_cast<int>(numbins) ) 
	      k = 0;
	 		
	   

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
		    
		  }
	    }
	
	    /* newflux vs t_e corresponds to the new re-binned light curve.
	       It corresponds to the same light curve as bolflux vs t_o */
	    
	    if ( newflux.at(i) < 0.0 )
	      newflux.at(i) = 0.0;
	  
	  } // closes the for(i) loop, going through the curves.
    	
	  /****************************************************************/
	  /* SUMMING THE FLUXES AT EACH PHASE BIN ACROSS ALL LIGHT CURVES */
	  /****************************************************************/
    	
	  // for ( unsigned int i(0); i < numbins; i++ ) {
	  //totflux.at(i) += newflux.at(i);   // setting totflux = newflux
	  //}
	  //if (p==0)
	  //ttt.open("newtime.txt", std::ios_base::trunc);

	  for ( unsigned int i(0); i < numbins; i++ ) {
	    //if ( p==0 ) 
	    //ttt << curve.t_o[i] << " " 
	    //	  << curve.f[p][i] << " " << curve.t[i] << " " << totflux.at(i) 
	    //	  << " " << i
	    //	  << std::endl;
	    
	    //            curve.f[p][i] = totflux.at(i);
            curve.f[p][i] = newflux.at(i);
	  }
	  // only plotting versus evenly spaced time -- not using t_o
	
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
    /*
	if (abs(x - Xa) < 0.0001*abs(x)){
		//cout << x << " " << Xa << " " << abs(x-Xa) << " " << 0.001*abs(x) << endl;
		return Ya;

	}else{
    	return  (Ya + ((Yb - Ya) * (x - Xa) / (Xb - Xa)));
	}
	*/
    return  (Ya + ((Yb - Ya) * (x - Xa) / (Xb - Xa)));
}

// Linear interpolation in log-log space
double LogLinear(double x,double Xa, double Ya, double Xb, double Yb){
    /*
    if (abs(log(x) - log(Xa)) < 0.0001*log(x)){
    	//cout << log(x) << " " << log(Xa) << " " << abs(log(x) - log(Xa)) << " " << 0.01*log(x) << endl;
    	return Ya;
    }else{
	    return  (exp(log(Ya) + ((log(Yb) - log(Ya)) * (log(x) - log(Xa)) / (log(Xb) - log(Xa)))));
    }
    */
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
std::vector<double> mu_2, Es,F,FF,FFF,FFFF,I,II,III,IIII;
double X, Y, X1, Y1, X2, Y2;

// Read the four hydrogen intensity files to peform the four point interpolation
void Read_NSATMOS(double T, double M, double R){
    
    double delta, lgrav, lt, temp;
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
            H_atmo_para >> temp;
            mu_2.push_back(temp);
        }
    H_atmo_para.close();
    }

    //Find correct logt and lgrav paramter choice
    M = Units::nounits_to_cgs(M, Units::MASS);
    R = Units::nounits_to_cgs(R, Units::LENGTH);
    delta = 1 / sqrt(1 - (2 * Units::G * M / (R * Units::C * Units::C)));
    lgrav = log10(delta * Units::G * M / (R * R));
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));
    cout << "temperature in log(K) is " << lt << endl;
    cout << "specific gravity in log(g) is " << lgrav << endl;
    
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
            Es.push_back(temp);
            temp = temp * 1E3 * Units::EV / Units::H_PLANCK;
            //cout << temp << endl;
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
            temp = temp * 1E3 * Units::EV / Units::H_PLANCK;
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
            temp = temp * 1E3 * Units::EV / Units::H_PLANCK;
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
            temp = temp * 1E3 * Units::EV / Units::H_PLANCK;
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
    double freq(0), P, ener_spacing, first_ener, freq_index; 
    double I_int[8],Q[4],R[2];
    int i_mu(0), n_mu, i_f(0), n_f, down, up, size_mu(12);
    //char atmodir[1024], cwd[1024];
    std::vector<double> mu,I_temp,II_temp,III_temp,IIII_temp,Iv_temp,IIv_temp,IIIv_temp,IIIIv_temp;

    //Convert energy point to frequency
    freq = 1E3 * E * Units::EV / Units::H_PLANCK;

    //Read in atmosphere parameters
    mu = mu_2;
    
    //Find proper mu choice
    for (int m = 0; m < size_mu; ++m) {
        if (cos_theta <= mu[m]) {
            i_mu = m;
        }
    }
    n_mu = i_mu + 1;

    //Find proper freqency choice
    ener_spacing = 1.0471;
    first_ener = 0.050119;
    freq_index = log(E / first_ener) / log(ener_spacing);
    i_f = (int) freq_index;
    n_f = i_f + 1;

    down = i_mu * 125;
    up = down + 125;
    for (int j = down; j < up; ++j) {
        I_temp.push_back(I[j]);
    }
    I_int[0] = LogLinear(E, Es[i_f], I_temp[i_f], Es[n_f], I_temp[n_f]);

    down = n_mu * 125;
    up = down + 125;
    for (int j = down; j < up; ++j) {
        Iv_temp.push_back(I[j]);
    }
    I_int[1] = LogLinear(E, Es[i_f], Iv_temp[i_f], Es[n_f], Iv_temp[n_f]);

    down = i_mu * 125;
    up = down + 125;
    for (int j = down; j < up; ++j) {
        II_temp.push_back(II[j]);
    }
    I_int[2] = LogLinear(E, Es[i_f], II_temp[i_f], Es[n_f], II_temp[n_f]);

    down = n_mu * 125;
    up = down + 125;
    for (int j = down; j < up; ++j) {
        IIv_temp.push_back(II[j]);
    }
    I_int[3] = LogLinear(E, Es[i_f], IIv_temp[i_f], Es[n_f], IIv_temp[n_f]);
                    
    down = i_mu * 125;
    up = down + 125;
    for (int j = down; j < up; ++j) {
        III_temp.push_back(III[j]);
    }
    I_int[4] = LogLinear(E, Es[i_f], III_temp[i_f], Es[n_f], III_temp[n_f]);

    down = n_mu * 125;
    up = down + 125;
    for (int j = down; j < up; ++j) {
        IIIv_temp.push_back(III[j]);
    }
    I_int[5] = LogLinear(E, Es[i_f], III_temp[i_f], Es[n_f], IIIv_temp[n_f]);
    
    down = i_mu * 125;
    up = down + 125;
    for (int j = down; j < up; ++j) {
        IIII_temp.push_back(IIII[j]);
    }
    I_int[6] = LogLinear(E, Es[i_f], IIII_temp[i_f], Es[n_f], IIII_temp[n_f]);

    down = n_mu * 125;
    up = down + 125;
    for (int j = down; j < up; ++j) {
        IIIIv_temp.push_back(IIII[j]);
    }
    I_int[7] = LogLinear(E, Es[i_f], IIIIv_temp[i_f], Es[n_f], IIIIv_temp[n_f]);

    //chdir(cwd);


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

double Hydrogen2(int E_dex, double cos_theta){
    double P, size_mu(12);
    double Q[4],R[2];
    int i_mu(0), down, mid, up;
    //char atmodir[1024], cwd[1024];
    std::vector<double> mu, I_temp,Iv_temp,II_temp,IIv_temp,III_temp,IIIv_temp,IIII_temp,IIIIv_temp;

    // Read in helium atmosphere parameter choices
    mu = mu_2;

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
    //chdir(cwd);

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
            He_atmo_para >> temp;
            mu_2.push_back(temp);
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
    double freq, P, size_mu(23), ener_spacing, first_ener, freq_index;
    double I_int[8],Q[4],R[2];
    int i_mu(0), down, mid, up, i_f, n_f;
    char atmodir[1024], cwd[1024];
    std::vector<double> mu, F_temp,FF_temp,FFF_temp,FFFF_temp,I_temp,II_temp,III_temp,IIII_temp,Iv_temp,IIv_temp,IIIv_temp,IIIIv_temp;

    // Convert energy to frequency
    freq = 1E3 * E * Units::EV / Units::H_PLANCK;

    // Read in helium atmosphere parameter choices
    getcwd(cwd, sizeof(cwd));
    sprintf(atmodir,"%s/atmosphere",cwd);
    chdir(atmodir);
    ifstream He_atmo_para;
    He_atmo_para.open("nsx-info.txt");

    mu = mu_2;
    //finding mu value 
    for (int m = 0; m < size_mu; ++m) {
        if (cos_theta <= mu[m]) {
            i_mu = m;
        }
    }

    //Find proper freqency choice
    ener_spacing = 1.0471;
    first_ener = 0.050119;
    freq_index = log(E / first_ener) / log(ener_spacing);
    i_f = (int) freq_index;
    n_f = i_f + 1;
 
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
    I_int[0] = LogLinear(E, Es[i_f], I_temp[i_f], Es[n_f], I_temp[n_f]);
    I_int[1] = LogLinear(E, Es[i_f], Iv_temp[i_f], Es[n_f], Iv_temp[n_f]);

    for (int j = down; j < mid; ++j) {
        II_temp.push_back(II[j]);
    }
    for (int j = mid; j < up; ++j){
        IIv_temp.push_back(II[j]);
    }
    I_int[2] = LogLinear(E, Es[i_f], II_temp[i_f], Es[n_f], II_temp[n_f]);
    I_int[3] = LogLinear(E, Es[i_f], IIv_temp[i_f], Es[n_f], IIv_temp[n_f]);


    for (int j = down; j < mid; ++j) {
        III_temp.push_back(III[j]);
    }
    for (int j = mid; j < up; ++j) {
        IIIv_temp.push_back(III[j]);
    }
    I_int[4] = LogLinear(E, Es[i_f], III_temp[i_f], Es[n_f], III_temp[n_f]);
    I_int[5] = LogLinear(E, Es[i_f], IIIv_temp[i_f], Es[n_f], IIIv_temp[n_f]);
    

    for (int j = down; j < mid; ++j) {
        IIII_temp.push_back(IIII[j]);
    }
    for (int j = mid; j < up; ++j) {
        IIIIv_temp.push_back(IIII[j]);
    }
    I_int[6] = LogLinear(E, Es[i_f], IIII_temp[i_f], Es[n_f], IIII_temp[n_f]);
    I_int[7] = LogLinear(E, Es[i_f], IIIIv_temp[i_f], Es[n_f], IIIIv_temp[n_f]);
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
    double P, size_mu(23);
    double Q[4],R[2];
    int i_mu(0), down, mid, up;
    char atmodir[1024], cwd[1024];
    std::vector<double> mu, I_temp,Iv_temp,II_temp,IIv_temp,III_temp,IIIv_temp,IIII_temp,IIIIv_temp;

    // Read in helium atmosphere parameter choices
    getcwd(cwd, sizeof(cwd));
    sprintf(atmodir,"%s/atmosphere",cwd);
    chdir(atmodir);
    ifstream He_atmo_para;
    He_atmo_para.open("nsx-info.txt");

    mu = mu_2;

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
/* NSX Hydrogen:                                                                      */
/*    This version works for the first NSXH files from Wynn.                          */
/*                                                                                    */
/**************************************************************************************/
void Read_NSXH(double T, double M, double R){
    
    double delta, lgrav, lt, temp, dump;
    char atmodir[1024], cwd[1024];
    std::vector<double> logt, lsgrav; 

    //Read in hydrogen atmosphere parameters
    getcwd(cwd, sizeof(cwd));
    sprintf(atmodir,"%s/atmosphere",cwd);
    chdir(atmodir);

    //Find correct logt and lgrav paramter choice
    M = Units::nounits_to_cgs(M, Units::MASS);
    R = Units::nounits_to_cgs(R, Units::LENGTH);
    delta = 1 / sqrt(1 - (2 * Units::G * M / (R * Units::C * Units::C)));
    lgrav = log10(delta * Units::G * M / (R * R));
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));
    cout << "temperature in log(K) is " << lt << endl;
    cout << "gravity in log(cgs units) is " << lgrav << endl;

    //Load hydrogen atmosphere files
    ifstream H_table1;
    if(lt <= 6.1 && lt >= 6){
    	H_table1.open("nsx_spint0_6.05g1425_nrp11.out");
    	cout << "loading log(T) = 6.05 file" << endl;
	}

	if(lt <= 5.7 && lt >= 5.6){
    	H_table1.open("nsx_spint0_5.69g1425_nrp1.out");
    	cout << "loading log(T) = 5.69 file" << endl;
	}

	if(lt <= 6.35 && lt >= 6.25){
    	H_table1.open("nsx_spint0_6.3g1425_nrp1.out");
    	cout << "loading log(T) = 6.3 file" << endl;
	}
    
    if(H_table1.is_open()){
    	for (int j = 1; j <= 100; j++) {
    		H_table1 >> temp;
            F.push_back(temp);
            temp = temp / 1E3 / Units::EV * Units::H_PLANCK;
            Es.push_back(temp);
            //cout << temp << endl;
            H_table1 >> dump;
            H_table1 >> temp;
            I.push_back(temp);

            for (int i = 1; i <= 255; i++) {
            	H_table1 >> dump;
            	H_table1 >> dump;
            	H_table1 >> temp;
            	I.push_back(temp);
        	}
    	}
    }else{
        cout << "NSXH file is not found" << endl;
    }
    H_table1.close();

    ifstream H_table2;
    if(lt <= 6.1 && lt >= 6){
    	H_table2.open("nsx_spint0_6.05g1425_nrp11.out");
    	cout << "loading log(T) = 6.05 file" << endl;
	}

	if(lt <= 5.7 && lt >= 5.6){
    	H_table2.open("nsx_spint0_5.69g1425_nrp1.out");
    	cout << "loading log(T) = 5.69 file" << endl;
	}

	if(lt <= 6.35 && lt >= 6.25){
    	H_table2.open("nsx_spint0_6.3g1425_nrp1.out");
    	cout << "loading log(T) = 6.3 file" << endl;
	}
    
    if(H_table2.is_open()){
        for (int i = 1; i <= 256; i++) {
            H_table2 >> dump;
            H_table2 >> temp;
            mu_2.push_back(temp);
            H_table2 >> dump;
        }
    }else{
        cout << "NSXH file is not found (while in second reading stage) " << endl;
    }
    //cout << I[0] << " " << F[0] << " " << Es[0] << " " << mu_2[0] << endl;
    H_table2.close();
    chdir(cwd);
}

// Calculate the final interpolated intensity
double NSXH(double E, double cos_theta){
    double freq, P, mu_spacing, theta, mu_index, freq_spacing, first_freq, freq_index;
    double I_int[2];
    int i_mu(0), n_mu, i_f(0), n_f, e_size(100);
    std::vector<double> mu, I_temp,Iv_temp;

    //cout << "starting NSXH" << endl;
    //Convert energy point to frequency
    freq = 1E3 * E * Units::EV / Units::H_PLANCK;

    mu = mu_2;
    
    //Find proper mu choice
    mu_spacing = ((Units::PI/2) - 0.0047) / 255;
    theta = acos (cos_theta);
    mu_index = ((Units::PI/2) - theta) / mu_spacing;
    i_mu = (int) mu_index;
    n_mu = i_mu + 1;

    
    //Read and interpolate to proper frequency
    for (int i = 0; i <= 99; i++){
    	I_temp.push_back(I[i*256+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	Iv_temp.push_back(I[i*256+n_mu]);   	
    }
    
    //Find proper freqency choice
    freq_spacing = exp(log(F[e_size-1]/F[0])/(e_size-1));
    first_freq = F[0];
    freq_index = log(freq / first_freq) / log(freq_spacing);
    i_f = (int) freq_index;
    n_f = i_f + 1;
    //cout << F[e_size-1] << " " << F[0] << endl;
    //cout << freq_spacing << " " << first_freq << " " << freq_index << " " << i_f << " " << n_f << endl;

    I_int[0] = LogLinear(freq, F[i_f], I_temp[i_f], F[n_f], I_temp[n_f]);
    I_int[1] = LogLinear(freq, F[i_f], Iv_temp[i_f], F[n_f], Iv_temp[n_f]);
    

    // Perform interpolation to correct mu (cos_theta)
    P = Linear(cos_theta,mu[i_mu],I_int[0],mu[n_mu],I_int[1]);
    cout << P << endl;

    return P;
}


double NSXH2(int E_dex, double cos_theta){
    double P, mu_spacing, theta, mu_index;
    int i_mu(0), n_mu(0);
    std::vector<double> mu, I_temp,Iv_temp,II_temp,IIv_temp,III_temp,IIIv_temp,IIII_temp,IIIIv_temp;

    //Read in atmosphere parameters

    mu = mu_2;
    //Find proper mu choice
    mu_spacing = ((Units::PI/2) - 0.0047) / 255;
    theta = acos (cos_theta);
    mu_index = ((Units::PI/2) - theta) / mu_spacing;

    i_mu = (int) mu_index;
    n_mu = i_mu + 1;
    //cout << mu_index << " " << i_mu << " " << n_mu << endl;

    //Read in the two sets intensities for two angles
    for (int i = 0; i <= 99; i++){
    	I_temp.push_back(I[i*256+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	Iv_temp.push_back(I[i*256+n_mu]);   	
    }
                       

    // Interpolate to chosen mu
    P = Linear(cos_theta,mu[i_mu],I_temp[E_dex],mu[i_mu+1],Iv_temp[E_dex]);

    return P;
}


/**************************************************************************************/
/* McPHAC Hydrogen:                                                                   */
/*    This version works for the McPHAC files from Slavko                             */
/*                                                                                    */
/**************************************************************************************/
void Read_McPHAC(double T, double M, double R){
    
    double delta, lgrav, lt, temp, dump;
    char atmodir[1024], cwd[1024];
    std::vector<double> logt, lsgrav; 

    //Read in hydrogen atmosphere parameters
    getcwd(cwd, sizeof(cwd));
    sprintf(atmodir,"%s/atmosphere",cwd);
    chdir(atmodir);

    //Find correct logt and lgrav paramter choice
    M = Units::nounits_to_cgs(M, Units::MASS);
    R = Units::nounits_to_cgs(R, Units::LENGTH);
    delta = 1 / sqrt(1 - (2 * Units::G * M / (R * Units::C * Units::C)));
    lgrav = log10(delta * Units::G * M / (R * R));
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));
    cout << "temperature in log(K) is " << lt << endl;
    cout << "gravity in log(cgs units) is " << lgrav << endl;

    //Load hydrogen atmosphere files
    ifstream H_table1;
    if(lt <= 6.1 && lt >= 6){
    	H_table1.open("McPHAC_logT6.05_EmergentSpectrum.4401.7.dat");
    	cout << "loading log(T) = 6.05 file" << endl;
	}

	if(lt <= 5.7 && lt >= 5.6){
    	H_table1.open("McPHAC_logT5.69_EmergentSpectrum.4401.7.dat");
    	cout << "loading log(T) = 5.69 file" << endl;
	}

	if(lt <= 6.35 && lt >= 6.25){
    	H_table1.open("McPHAC_logT6.3_EmergentSpectrum.4401.7.dat");
    	cout << "loading log(T) = 6.3 file" << endl;
	}
    
    if(H_table1.is_open()){
    	for (int j = 1; j <= 100; j++) {
    		H_table1 >> temp;
            F.push_back(temp);
            temp = temp / 1E3 / Units::EV * Units::H_PLANCK;
            Es.push_back(temp);
            H_table1 >> dump;
            H_table1 >> temp;
            I.push_back(temp);

            for (int i = 1; i <= 255; i++) {
            	H_table1 >> dump;
            	H_table1 >> dump;
            	H_table1 >> temp;
            	I.push_back(temp);
        	}
    	}
    }else{
        cout << "McPHAC file is not found" << endl;
    }
    H_table1.close();

    ifstream H_table2;
    if(lt <= 6.1 && lt >= 6){
    	H_table2.open("McPHAC_logT6.05_EmergentSpectrum.4401.7.dat");
    	cout << "loading log(T) = 6.05 file" << endl;
	}

	if(lt <= 5.7 && lt >= 5.6){
    	H_table2.open("McPHAC_logT5.69_EmergentSpectrum.4401.7.dat");
    	cout << "loading log(T) = 5.69 file" << endl;
	}

	if(lt <= 6.35 && lt >= 6.25){
    	H_table2.open("McPHAC_logT6.3_EmergentSpectrum.4401.7.dat");
    	cout << "loading log(T) = 6.3 file" << endl;
	}
    
    if(H_table2.is_open()){
        for (int i = 1; i <= 256; i++) {
            H_table2 >> dump;
            H_table2 >> temp;
            mu_2.push_back(temp);
            H_table2 >> dump;
        }
    }else{
        cout << "McPHAC file is not found (while in second reading stage) " << endl;
    }
    //cout << I[0] << " " << F[0] << " " << Es[0] << " " << mu_2[0] << endl;
    H_table2.close();
    chdir(cwd);
}

// Calculate the final interpolated intensity
double McPHAC(double E, double cos_theta){
    double freq, P, P1, mu_spacing, theta, mu_index, freq_spacing, first_freq, freq_index;
    double I_int[2];
    int i_mu(0), n_mu, i_f(0), n_f, e_size(100);
    std::vector<double> mu, I_temp,Iv_temp;

    //cout << "starting NSXH" << endl;
    //Convert energy point to frequency
    freq = 1E3 * E * Units::EV / Units::H_PLANCK;

    mu = mu_2;
    
    //Find proper mu choice
    mu_spacing = ((Units::PI/2) - 0.0047) / 255;
    theta = acos (cos_theta);
    mu_index = ((Units::PI/2) - theta) / mu_spacing;
    i_mu = (int) mu_index;
    n_mu = i_mu + 1;
    //cout << i_mu << " " << n_mu << endl;

    
    //Read and interpolate to proper frequency
    for (int i = 0; i <= 99; i++){
    	I_temp.push_back(I[i*256+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	Iv_temp.push_back(I[i*256+n_mu]);   	
    }
    
    //Find proper freqency choice
    freq_spacing = exp(log(F[e_size-1]/F[0])/(e_size-1));
    first_freq = F[0];
    freq_index = log(freq / first_freq) / log(freq_spacing);
    i_f = (int) freq_index;
    n_f = i_f + 1;
    //cout << F[e_size-1] << " " << F[0] << endl;
    //cout << freq_spacing << " " << first_freq << " " << freq_index << " " << i_f << " " << n_f << endl;

    I_int[0] = LogLinear(freq, F[i_f], I_temp[i_f], F[n_f], I_temp[n_f]);
    I_int[1] = LogLinear(freq, F[i_f], Iv_temp[i_f], F[n_f], Iv_temp[n_f]);
    
    //cout << E << endl;

    // Perform interpolation to correct mu (cos_theta)
    //P = LogLinear(cos_theta,mu[i_mu],I_int[0],mu[n_mu],I_int[1]);
    P1 = LogLinear(cos_theta,cos((255-i_mu)*mu_spacing),I_int[0],cos((255-n_mu)*mu_spacing),I_int[1]);
    
    cout << "Intensity = " << P1 << endl;

    return P1;
}


double McPHAC2(int E_dex, double cos_theta){
    double P, P1, mu_spacing, theta, mu_index;
    int i_mu(0), n_mu(0);
    std::vector<double> mu, I_temp,Iv_temp,II_temp,IIv_temp,III_temp,IIIv_temp,IIII_temp,IIIIv_temp;

    //Read in atmosphere parameters

    mu = mu_2;
    //Find proper mu choice
    mu_spacing = ((Units::PI/2) - 0.0047) / 255;
    theta = acos (cos_theta);
    mu_index = ((Units::PI/2) - theta) / mu_spacing;

    i_mu = (int) mu_index;
    n_mu = i_mu + 1;
    //cout << mu_index << " " << i_mu << " " << n_mu << endl;

    //Read in the two sets intensities for two angles
    for (int i = 0; i <= 99; i++){
    	I_temp.push_back(I[i*256+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	Iv_temp.push_back(I[i*256+n_mu]);   	
    }
                       

    // Interpolate to chosen mu
    //P = LogLinear(cos_theta,mu[i_mu],I_temp[E_dex],mu[i_mu+1],Iv_temp[E_dex]);
    P1 = LogLinear(cos_theta,cos((255-i_mu)*mu_spacing),I_temp[E_dex],cos((255-n_mu)*mu_spacing),Iv_temp[E_dex]);
    return P1;
}





/**************************************************************************************/
/* NSX Helium (new):                                                                  */
/*    This version works for the first NSXHe files from Wynn.                         */
/*                                                                                    */
/**************************************************************************************/
void Read_NSXHe(double T, double M, double R){
    
    double delta, lgrav, lt, temp, dump;
    char atmodir[1024], cwd[1024];
    std::vector<double> logt, lsgrav; 

    //Read in hydrogen atmosphere parameters
    getcwd(cwd, sizeof(cwd));
    sprintf(atmodir,"%s/atmosphere",cwd);
    chdir(atmodir);

    //Find correct logt and lgrav paramter choice
    M = Units::nounits_to_cgs(M, Units::MASS);
    R = Units::nounits_to_cgs(R, Units::LENGTH);
    delta = 1 / sqrt(1 - (2 * Units::G * M / (R * Units::C * Units::C)));
    lgrav = log10(delta * Units::G * M / (R * R));
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));
    cout << "temperature in log(K) is " << lt << endl;
    cout << "gravity in log(cgs units) is " << lgrav << endl;

    //Load hydrogen atmosphere files
    ifstream H_table1;
    if(lt <= 6.1 && lt >= 6){
    	H_table1.open("nsx_spint0_6.05g1425He_nrp1.out");
    	cout << "loading log(T) = 6.05 file" << endl;
	}

	if(lt <= 5.7 && lt >= 5.6){
    	H_table1.open("nsx_spint0_5.69g1425He_nrp1.out");
    	cout << "loading log(T) = 5.69 file" << endl;
	}

	if(lt <= 6.35 && lt >= 6.25){
    	H_table1.open("nsx_spint0_6.3g1425He_nrp1.out");
    	cout << "loading log(T) = 6.3 file" << endl;
	}
    
    if(H_table1.is_open()){
    	for (int j = 1; j <= 100; j++) {
    		H_table1 >> temp;
            F.push_back(temp);
            temp = temp / 1E3 / Units::EV * Units::H_PLANCK;
            Es.push_back(temp);
            H_table1 >> dump;
            H_table1 >> temp;
            I.push_back(temp);

            for (int i = 1; i <= 255; i++) {
            	H_table1 >> dump;
            	H_table1 >> dump;
            	H_table1 >> temp;
            	I.push_back(temp);
        	}
    	}
    }else{
        cout << "NSXH file is not found" << endl;
    }
    H_table1.close();

    ifstream H_table2;
    if(lt <= 6.1 && lt >= 6){
    	H_table2.open("nsx_spint0_6.05g1425He_nrp1.out");
    	cout << "loading log(T) = 6.05 file" << endl;
	}

	if(lt <= 5.7 && lt >= 5.6){
    	H_table2.open("nsx_spint0_5.69g1425He_nrp1.out");
    	cout << "loading log(T) = 5.69 file" << endl;
	}

	if(lt <= 6.35 && lt >= 6.25){
    	H_table2.open("nsx_spint0_6.3g1425He_nrp1.out");
    	cout << "loading log(T) = 6.3 file" << endl;
	}
    
    if(H_table2.is_open()){
        for (int i = 1; i <= 256; i++) {
            H_table2 >> dump;
            H_table2 >> temp;
            mu_2.push_back(temp);
            H_table2 >> dump;
        }
    }else{
        cout << "NSXH file is not found (while in second reading stage) " << endl;
    }
    //cout << I[0] << " " << F[0] << " " << Es[0] << " " << mu_2[0] << endl;
    H_table2.close();
    chdir(cwd);
}

// Calculate the final interpolated intensity
double NSXHe(double E, double cos_theta){
    double freq, P, mu_spacing, theta, mu_index, freq_spacing, first_freq, freq_index;
    double I_int[2];
    int i_mu(0), n_mu, i_f(0), n_f, e_size(100);
    std::vector<double> mu, I_temp,Iv_temp;

    //cout << "starting NSXH" << endl;
    //Convert energy point to frequency
    freq = 1E3 * E * Units::EV / Units::H_PLANCK;

    mu = mu_2;
    
    //Find proper mu choice
    mu_spacing = ((Units::PI/2) - 0.0047) / 255;
    theta = acos (cos_theta);
    mu_index = ((Units::PI/2) - theta) / mu_spacing;
    i_mu = (int) mu_index;
    n_mu = i_mu + 1;

    
    //Read and interpolate to proper frequency
    for (int i = 0; i <= 99; i++){
    	I_temp.push_back(I[i*256+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	Iv_temp.push_back(I[i*256+n_mu]);   	
    }
    
    //Find proper freqency choice
    freq_spacing = exp(log(F[e_size-1]/F[0])/(e_size-1));
    first_freq = F[0];
    freq_index = log(freq / first_freq) / log(freq_spacing);
    i_f = (int) freq_index;
    n_f = i_f + 1;
    //cout << F[e_size-1] << " " << F[0] << endl;
    //cout << freq_spacing << " " << first_freq << " " << freq_index << " " << i_f << " " << n_f << endl;

    I_int[0] = LogLinear(freq, F[i_f], I_temp[i_f], F[n_f], I_temp[n_f]);
    I_int[1] = LogLinear(freq, F[i_f], Iv_temp[i_f], F[n_f], Iv_temp[n_f]);
    

    // Perform interpolation to correct mu (cos_theta)
    P = LogLinear(cos_theta,mu[i_mu],I_int[0],mu[n_mu],I_int[1]);
    //cout << P << endl;

    return P;
}


double NSXHe2(int E_dex, double cos_theta){
    double P, mu_spacing, theta, mu_index;
    int i_mu(0), n_mu(0);
    std::vector<double> mu, I_temp,Iv_temp,II_temp,IIv_temp,III_temp,IIIv_temp,IIII_temp,IIIIv_temp;

    //Read in atmosphere parameters

    mu = mu_2;
    //Find proper mu choice
    mu_spacing = ((Units::PI/2) - 0.0047) / 255;
    theta = acos (cos_theta);
    mu_index = ((Units::PI/2) - theta) / mu_spacing;

    i_mu = (int) mu_index;
    n_mu = i_mu + 1;
    //cout << mu_index << " " << i_mu << " " << n_mu << endl;

    //Read in the two sets intensities for two angles
    for (int i = 0; i <= 99; i++){
    	I_temp.push_back(I[i*256+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	Iv_temp.push_back(I[i*256+n_mu]);   	
    }
                       

    // Interpolate to chosen mu
    P = LogLinear(cos_theta,mu[i_mu],I_temp[E_dex],mu[i_mu+1],Iv_temp[E_dex]);

    return P;
}



/**************************************************************************************/
/* McPHAC Hydrogen:                                                                   */
/*    This version works for the original McPHAC files from Cole                      */
/*                                                                                    */
/**************************************************************************************/
void Read_McPHACC(double T, double M, double R){
    
    double delta, lt, lgrav, temp, dump, real_T;
    int size_logt(13), size_lsgrav(8), size_mu(23), i_lt, i_lgrav, n_lt, n_lgrav, size_ener(125), size_set, skip_to, skip_two;
    char s1[40], s2[40], s3[40], s4[40], atmodir[1024], cwd[1024];
    std::vector<double> freq, logt, lsgrav;

    //setting values of lt and lgrav based on input T, M, and R. Also sets to load using first mu value.   
    M = Units::nounits_to_cgs(M, Units::MASS);
    R = Units::nounits_to_cgs(R, Units::LENGTH);
    delta = 1 / sqrt(1 - (2 * Units::G * M / (R * Units::C * Units::C)));
    lgrav = log10(delta * Units::G * M / (R * R));
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));
    cout << "temperature in log(K) is " << lt << endl;
    cout << "gravity in log(cgs units) is " << lgrav << endl;

    //obtain helium atmosphere parameters
    getcwd(cwd, sizeof(cwd));
    sprintf(atmodir,"%s/atmosphere/mcphacc",cwd);
    chdir(atmodir);

    i_lt = (lt-5.1)/0.05; //if we need to load 1st temperature, i_lt = 0. this is discrete math
    i_lgrav = (lgrav-13.7)/0.1;

    n_lt = i_lt+1;
    n_lgrav = i_lgrav+1;

    sprintf(s1,"mcphacc_T%d_g%d.txt",n_lt,n_lgrav);
    sprintf(s2,"mcphacc_T%d_g%d.txt",n_lt+1,n_lgrav);  
    sprintf(s3,"mcphacc_T%d_g%d.txt",n_lt,n_lgrav+1); 
    sprintf(s4,"mcphacc_T%d_g%d.txt",n_lt+1,n_lgrav+1);

    cout << s1 << endl;
    cout << s2 << endl;
    cout << s3 << endl;
    cout << s4 << endl; 

    real_T = pow(10.0,lt);

    ifstream H_table1;
    H_table1.open(s1);
    ifstream H_table2;
    H_table2.open(s2);
    ifstream H_table3;
    H_table3.open(s3);
    ifstream H_table4;
    H_table4.open(s4);

    if(H_table1.is_open()){
    	for (int i = 1; i <= 5000; i++) {
    		H_table1 >> dump;
    		H_table1 >> dump;
    		H_table1 >> temp;
    		temp = temp*pow(real_T,3);
        	I.push_back(temp);
        }
        H_table1.close();    	
    }else{
    	cout << s1 << " is not found." << endl;
    }

    if(H_table2.is_open()){
    	for (int i = 1; i <= 5000; i++) {
    		H_table2 >> dump;
    		H_table2 >> dump;
    		H_table2 >> temp;
    		temp = temp*pow(real_T,3);
        	II.push_back(temp);
        }
        H_table2.close();    	
    }else{
    	cout << s2 << " is not found." << endl;
    }     

    if(H_table3.is_open()){
    	for (int i = 1; i <= 5000; i++) {
    		H_table3 >> dump;
    		H_table3 >> dump;
    		H_table3 >> temp;
    		temp = temp*pow(real_T,3);
        	III.push_back(temp);
        }
        H_table3.close();    	
    }else{
    	cout << s3 << " is not found." << endl;
    }

    if(H_table4.is_open()){
    	for (int i = 1; i <= 5000; i++) {
    		H_table4 >> dump;
    		H_table4 >> dump;
    		H_table4 >> temp;
    		temp = temp*pow(real_T,3);
        	IIII.push_back(temp);
        }
        H_table4.close();    	
    }else{
    	cout << s4 << " is not found." << endl;
    }

    //cout << i_lt << " " << n_lt << " " << i_lgrav << " " << n_lgrav << endl;
	/* 
    ifstream H_table1;
    H_table1.open("cole_mcphac_table2.txt");

    if(H_table1.is_open()){

    	skip_to = (i_lt * 11 + i_lgrav)*5000; //skipping unwanted lines.

    	for (int j = 1; j <= skip_to; j++) {
     		H_table1 >> dump;	//skipping temperature
     		H_table1 >> dump;   //skipping specific gravity
     		H_table1 >> dump;	//skipping energy
     		H_table1 >> dump;	//skipping angle
     		H_table1 >> dump;	//skipping intensity
    	}

    	for (int j = 1; j <= 5000; j++) {
    		H_table1 >> temp;	//skipping temperature
    		X1 = temp;
     		H_table1 >> temp;   //skipping specific gravity
     		Y1 = temp;
     		H_table1 >> dump;	//skipping energy
     		H_table1 >> dump;	//skipping angle
    		H_table1 >> temp;
    		temp = temp*pow(real_T,3);
            I.push_back(temp);
    	}

    	for (int j = 1; j <= 5000; j++) {
    		H_table1 >> dump;	//skipping temperature
     		H_table1 >> temp;   //skipping specific gravity
     		Y2 = temp;
     		H_table1 >> dump;	//skipping energy
     		H_table1 >> dump;	//skipping angle
    		H_table1 >> temp;
    		temp = temp*pow(real_T,3);
            III.push_back(temp);
    	}

    	skip_two = (11-2)*5000;

    	for (int j = 1; j <= skip_two; j++) {
     		H_table1 >> dump;	//skipping temperature
     		H_table1 >> dump;   //skipping specific gravity
     		H_table1 >> dump;	//skipping energy
     		H_table1 >> dump;	//skipping angle
     		H_table1 >> dump;	//skipping intensity
    	}

    	for (int j = 1; j <= 5000; j++) {
    		H_table1 >> temp;	//skipping temperature
    		X2 = temp;
     		H_table1 >> dump;   //skipping specific gravity
     		H_table1 >> dump;	//skipping energy
     		H_table1 >> dump;	//skipping angle
    		H_table1 >> temp;
    		temp = temp*pow(real_T,3);
            II.push_back(temp);
    	}

    	for (int j = 1; j <= 5000; j++) {
    		H_table1 >> dump;	//skipping temperature
     		H_table1 >> dump;   //skipping specific gravity
     		H_table1 >> dump;	//skipping energy
     		H_table1 >> dump;	//skipping angle
    		H_table1 >> temp;
    		temp = temp*pow(real_T,3);
            IIII.push_back(temp);
    	}
    	

    }else{
        cout << "McPHAC (Cole's) file is not found" << endl;
    }
    H_table1.close();
    */

    cout << "finished loading Cole's McPHAC intensities" << endl; 

    ifstream H_table5;
    H_table5.open("mcphacc_T1_g1.txt");
    
    if(H_table5.is_open()){
        H_table5 >> temp;
        temp = pow(10.0,temp);
        Es.push_back(temp);	//recording first energy
        //cout << temp << endl;
        temp = temp * 1E3 * Units::EV / Units::H_PLANCK;
        F.push_back(temp);
        for (int i = 1; i <= 49; i++){
        	H_table5 >> temp;
        	mu_2.push_back(temp);	//recording first 49 angles
    		H_table5 >> dump;	//skipping intensity
     		H_table5 >> dump;   //skipping energy        	
        }
        H_table5 >> temp;
        mu_2.push_back(temp);	//recording 50th angle;
        H_table5 >> dump;	//skipping intensity

        //entering second energy group
        for (int j = 2; j <= 100; j++){
        	H_table5 >> temp;
        	temp = pow(10.0,temp);
        	Es.push_back(temp);	//recording second and up to last energy
        	//cout << temp << endl;
        	temp = temp * 1E3 * Units::EV / Units::H_PLANCK;
        	F.push_back(temp);
        	for (int i = 1; i <= 49; i++){
    			H_table5 >> dump;	//skipping angle
     			H_table5 >> dump;   //skipping intensity
     			H_table5 >> dump;	//skipping energy
        	}
        	H_table5 >> dump;	//skipping angle
     		H_table5 >> dump;   //skipping intensity      	
        }
    }else{
        cout << "McPHAC (Cole's) file is not found (while in second reading stage) " << endl;
    }
    //cout << I[0] << " " << F[0] << " " << Es[0] << " " << mu_2[0] << endl;
    H_table5.close();   

    X = pow(10.0,lt);
    Y = pow(10.0,lgrav);
    X1 = pow(10.0,5.1+0.05*i_lt);
    X2 = pow(10.0,5.1+0.05*n_lt);
    Y1 = pow(10.0,13.7+0.1*i_lgrav);
    Y2 = pow(10.0,13.7+0.1*n_lgrav);
    //cout << T << " " << lt << endl;
    chdir(cwd);
	/*
    for (int i = 0; i <= 99; i++){
    	cout << Es[i] << endl;
    }
    */
    
    
}

// Calculate the final interpolated intensity
double McPHACC(double E, double cos_theta){
    double freq, P, mu_spacing, theta, mu_index, ener_spacing, first_ener, freq_index;
    double I_int[8],Q[4],R[2];
    int i_mu(0), n_mu, i_f(0), n_f, e_size(100);
    std::vector<double> mu, I_temp, Iv_temp, II_temp, IIv_temp, III_temp, IIIv_temp, IIII_temp, IIIIv_temp;

    //cout << "starting NSXH" << endl;
    //Convert energy point to frequency
    //freq = 1E3 * E * Units::EV / Units::H_PLANCK;

    mu = mu_2;
    
    //Find proper mu choice
    mu_spacing = ((Units::PI/2) - 0.024084) / 49;
    theta = acos (cos_theta);
    mu_index = ((Units::PI/2) - theta) / mu_spacing;
    i_mu = (int) mu_index;
    n_mu = i_mu + 1;
    //cout << i_mu << " " << n_mu << endl;

    //Find proper freqency choice
    ener_spacing = pow(10.0,0.0338);
    first_ener = Es[0];
    freq_index = log(E / first_ener) / log(ener_spacing);
    //cout << E << " " << first_ener << " " << freq_index << endl;
    i_f = (int) freq_index;
    n_f = i_f + 1;
    
    //Read and interpolate to proper frequency
    for (int i = 0; i <= 99; i++){
    	I_temp.push_back(I[i*50+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	Iv_temp.push_back(I[i*50+n_mu]);   	
    }
    //cout << I_temp[i_f] << endl;
    I_int[0] = LogLinear(E, Es[i_f], I_temp[i_f], Es[n_f], I_temp[n_f]);
    I_int[1] = LogLinear(E, Es[i_f], Iv_temp[i_f], Es[n_f], Iv_temp[n_f]);

    for (int i = 0; i <= 99; i++){
    	II_temp.push_back(II[i*50+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	IIv_temp.push_back(II[i*50+n_mu]);   	
    }

    I_int[2] = LogLinear(E, Es[i_f], II_temp[i_f], Es[n_f], II_temp[n_f]);
    I_int[3] = LogLinear(E, Es[i_f], IIv_temp[i_f], Es[n_f], IIv_temp[n_f]);
    
    for (int i = 0; i <= 99; i++){
    	III_temp.push_back(III[i*50+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	IIIv_temp.push_back(III[i*50+n_mu]);   	
    }

    I_int[4] = LogLinear(E, Es[i_f], III_temp[i_f], Es[n_f], III_temp[n_f]);
    I_int[5] = LogLinear(E, Es[i_f], IIIv_temp[i_f], Es[n_f], IIIv_temp[n_f]);

    for (int i = 0; i <= 99; i++){
    	IIII_temp.push_back(IIII[i*50+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	IIIIv_temp.push_back(IIII[i*50+n_mu]);   	
    }

    I_int[6] = LogLinear(E, Es[i_f], IIII_temp[i_f], Es[n_f], IIII_temp[n_f]);
    I_int[7] = LogLinear(E, Es[i_f], IIIIv_temp[i_f], Es[n_f], IIIIv_temp[n_f]);

    // Interpolate to chosen mu
    Q[0] = Linear(theta,acos(mu[i_mu]),I_int[0],acos(mu[i_mu+1]),I_int[1]);
    Q[1] = Linear(theta,acos(mu[i_mu]),I_int[2],acos(mu[i_mu+1]),I_int[3]);
    Q[2] = Linear(theta,acos(mu[i_mu]),I_int[4],acos(mu[i_mu+1]),I_int[5]);
    Q[3] = Linear(theta,acos(mu[i_mu]),I_int[6],acos(mu[i_mu+1]),I_int[7]);


    // Set to highest mu values for emission normal to surface
    if (cos_theta > 0.999710) {
    	Q[0] = I_int[0];
    	Q[1] = I_int[2];
    	Q[2] = I_int[4];
    	Q[3] = I_int[6]; 
    	//cout << Q[0] << " " << Q[1] << " " << Q[2] << " " << Q[3] << endl;   	
    }

    // Interpolate to chosen local gravity
    R[0] = LogLinear(Y,Y1,Q[0],Y2,Q[2]);
    R[1] = LogLinear(Y,Y1,Q[1],Y2,Q[3]);

    // Interpolate to chosen temperature
    P = LogLinear(X,X1,R[0],X2,R[1]);

    // Set to zero at small angle
    if (cos_theta < 0.015629) P = 0;

    //cout << P << endl;

    //if (isnan(P)) cout << Q[0] << " " << R[0] << endl;
    //if (isnan(P)) cout << i_mu << " " << mu[i_mu] << " " << i_mu+1 << " " << mu[i_mu] << endl;
    //if (isnan(P)) cout << "P is nan!" << endl;
    //cout << P << endl;

    return P;
}


double McPHACC2(int E_dex, double cos_theta){
    double P, P1, mu_spacing, theta, mu_index;
    double Q[4],R[2];
    int i_mu(0), n_mu(0);
    std::vector<double> mu, I_temp,Iv_temp,II_temp,IIv_temp,III_temp,IIIv_temp,IIII_temp,IIIIv_temp;

    mu = mu_2;
    //Find proper mu choice
    mu_spacing = ((Units::PI/2) - 0.024084) / 49;
    theta = acos (cos_theta);
    mu_index = ((Units::PI/2) - theta) / mu_spacing;
    i_mu = (int) mu_index;
    n_mu = i_mu + 1;


    for (int i = 0; i <= 99; i++){
    	I_temp.push_back(I[i*50+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	Iv_temp.push_back(I[i*50+n_mu]);   	
    }

    for (int i = 0; i <= 99; i++){
    	II_temp.push_back(II[i*50+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	IIv_temp.push_back(II[i*50+n_mu]);   	
    }
    
    for (int i = 0; i <= 99; i++){
    	III_temp.push_back(III[i*50+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	IIIv_temp.push_back(III[i*50+n_mu]);   	
    }

    for (int i = 0; i <= 99; i++){
    	IIII_temp.push_back(IIII[i*50+i_mu]);
    	//cout << I[i*100+i_mu] << endl;
     	IIIIv_temp.push_back(IIII[i*50+n_mu]);   	
    }
 
    // Interpolate to chosen mu
    Q[0] = Linear(theta,acos(mu[i_mu]),I_temp[E_dex],acos(mu[n_mu]),Iv_temp[E_dex]);
    Q[1] = Linear(theta,acos(mu[i_mu]),II_temp[E_dex],acos(mu[n_mu]),IIv_temp[E_dex]);
    Q[2] = Linear(theta,acos(mu[i_mu]),III_temp[E_dex],acos(mu[n_mu]),IIIv_temp[E_dex]);
    Q[3] = Linear(theta,acos(mu[i_mu]),IIII_temp[E_dex],acos(mu[n_mu]),IIIIv_temp[E_dex]); 

    // Set to highest mu values for emission normal to surface
    if (cos_theta > 0.999710) {
    	Q[0] = I_temp[E_dex];
    	Q[1] = II_temp[E_dex];
    	Q[2] = III_temp[E_dex];
    	Q[3] = IIII_temp[E_dex]; 
    	//cout << cos_theta << " " << theta << " " << acos(mu[i_mu]) << " " << I_temp[E_dex] << " " << acos(mu[n_mu]) << endl;   	
    }

    // Interpolate to chosen local gravity
    R[0] = LogLinear(Y,Y1,Q[0],Y2,Q[2]);
    R[1] = LogLinear(Y,Y1,Q[1],Y2,Q[3]);

    // Interpolate to chosen temperature
    P = LogLinear(X,X1,R[0],X2,R[1]);

    // Set to zero at small angle
    if (cos_theta < 0.015629) P = 0;

    //if (isnan(P)) cout << Q[0] << " " << R[0] << endl;

    return P;
}


// Calculate the final interpolated intensity
double McPHACC3(double E, double cos_theta, double T, double M, double R, class LightCurve mexmcc){
	double delta, obl_approx, lgrav, lt, th_spacing, theta, th_index, ener_spacing, first_ener, ener_index;
	double e0, e1, th0, th1, grav0, grav1, t0, t1;
	double I_temp[16], I_int[8], J[4], K[2], L(0.0);
	int i_f, n_f, i_lt, i_lgrav, n_lt, n_lgrav, i_mu, n_mu, first_inte;


    //setting values of lt and lgrav based on input T, M, and R. Also sets to load using first mu value.   
    M = Units::nounits_to_cgs(M, Units::MASS);
    R = Units::nounits_to_cgs(R, Units::LENGTH);
    delta = 1 / sqrt(1 - (2 * Units::G * M / (R * Units::C * Units::C)));
    obl_approx = (1 + (-0.791 + 0.776 * mexmcc.para.mass_over_r) * pow(sin(mexmcc.para.omega_bar_sq),2) + (1.138 - 1.431 * mexmcc.para.mass_over_r) * pow(cos(mexmcc.para.omega_bar_sq),2));
    lgrav = log10(delta * Units::G * M / (R * R));
    //lgrav = log10(delta * Units::G * M / R / R * obl_approx);
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));
    //cout << "temperature in log(K) is " << lt << endl;
    //cout << "gravity in log(cgs units) is " << lgrav << endl;
    //cout << "cos_theta is " << cos_theta << endl;
    //cout << "energy is " << E << endl;

    i_lt = (lt-5.1)/0.05; //if we need to load 1st temperature, i_lt = 0. this is discrete math
    n_lt = i_lt+1;
    t0 = 5.1+0.05*i_lt;
    t1 = t0+0.05;
    //cout << i_lt << " " << n_lt << " " << t0 << " " << t1 << endl;
    //cout << t0 << " " << t1 << endl;

    i_lgrav = (lgrav-13.7)/0.1;
    n_lgrav = i_lgrav+1;
    grav0 = 13.7+0.1*i_lgrav;
    grav1 = grav0+0.1;
    //cout << grav0 << " " << grav1 << endl;
    
    //Find proper mu choice
    n_mu = 1;
    while (cos_theta > mexmcc.mccangl[n_mu] && n_mu < 49){
    	n_mu += 1;
    }

    i_mu = n_mu - 1;
    th0 = acos(mexmcc.mccangl[i_mu+1]);
    th1 = acos(mexmcc.mccangl[n_mu+1]);
    theta = acos(cos_theta);
    //cout << mexmcc.mccangl[i_mu] << " " << mexmcc.mccangl[n_mu] << " " << cos_theta << endl;
    //cout << th0 << " " << th1 << " " << acos(cos_theta) << endl;

    //Find proper freqency choice
    ener_spacing = pow(10.0,0.0338);
	first_ener = 0.004969469;
    ener_index = log(E / first_ener) / log(ener_spacing);
    i_f = (int) ener_index;
    n_f = i_f + 1;
    e0 = first_ener*pow(ener_spacing,i_f);
    e1 = e0*ener_spacing;
    //cout << e0 << " " << e1 << " " << E << endl;


    first_inte = (i_lt*11 + i_lgrav) * 5000 + i_f * 50 + i_mu +1;
    //cout << i_lt << " " << i_lgrav << " " << i_f << " " << i_mu << " " << first_inte << endl;
    I_temp[0] = mexmcc.mccinte[first_inte]*pow(10.0,t0*3.0);
    double freq = 1E3 * e0 * Units::EV / Units::H_PLANCK;
    //cout << mexmcc.mccinte[first_inte] << " " << pow(10.0,t0*3.0) << " " << t0 << endl;
    //cout << I_temp[0] << " " << first_inte << " " << t0 << " " << grav0 << " " << freq << " " << th0 << " " << mexmcc.mccangl[i_mu+1] << endl;
    //cout << mexmcc.mccangl[0] << " " << mexmcc.mccangl[49] << endl;
    I_temp[1] = mexmcc.mccinte[first_inte+1]*pow(10.0,t0*3.0);
    I_temp[2] = mexmcc.mccinte[first_inte+50]*pow(10.0,t0*3.0);
    I_temp[3] = mexmcc.mccinte[first_inte+51]*pow(10.0,t0*3.0);

    I_temp[4] = mexmcc.mccinte[first_inte+5000]*pow(10.0,t0*3.0);
    I_temp[5] = mexmcc.mccinte[first_inte+5001]*pow(10.0,t0*3.0);
    I_temp[6] = mexmcc.mccinte[first_inte+5050]*pow(10.0,t0*3.0);
    I_temp[7] = mexmcc.mccinte[first_inte+5051]*pow(10.0,t0*3.0);

    I_temp[8] = mexmcc.mccinte[first_inte+55000]*pow(10.0,t1*3.0);
    I_temp[9] = mexmcc.mccinte[first_inte+55001]*pow(10.0,t1*3.0);
    I_temp[10] = mexmcc.mccinte[first_inte+55050]*pow(10.0,t1*3.0);
    I_temp[11] = mexmcc.mccinte[first_inte+55051]*pow(10.0,t1*3.0);

    I_temp[12] = mexmcc.mccinte[first_inte+60000]*pow(10.0,t1*3.0);
    I_temp[14] = mexmcc.mccinte[first_inte+60001]*pow(10.0,t1*3.0);
    I_temp[13] = mexmcc.mccinte[first_inte+60050]*pow(10.0,t1*3.0);
    I_temp[15] = mexmcc.mccinte[first_inte+60051]*pow(10.0,t1*3.0);

    I_int[0] = LogLinear(E, e0, I_temp[0], e1, I_temp[2]); //t0, grav0, th0
    I_int[1] = LogLinear(E, e0, I_temp[1], e1, I_temp[3]); //t0, grav0, th1
    I_int[2] = LogLinear(E, e0, I_temp[4], e1, I_temp[6]); //t0, grav1, th0
    I_int[3] = LogLinear(E, e0, I_temp[5], e1, I_temp[7]); //t0, grav1, th1
    I_int[4] = LogLinear(E, e0, I_temp[8], e1, I_temp[10]);//t1, grav0, th0
    I_int[5] = LogLinear(E, e0, I_temp[9], e1, I_temp[11]);//t1, grav0, th1
    I_int[6] = LogLinear(E, e0, I_temp[12], e1, I_temp[14]);//t1, grav1, th0
    I_int[7] = LogLinear(E, e0, I_temp[13], e1, I_temp[15]);//t1, grav1, th1
    
    //cout << I_int[0] << " " << I_int[1] << " " << I_int[2] << " " << I_int[3] << endl;

    // Interpolate to chosen mu
    J[0] = Linear(theta,th0,I_int[0],th1,I_int[1]); //t0, grav0
    J[1] = Linear(theta,th0,I_int[2],th1,I_int[3]); //t0, grav1
    J[2] = Linear(theta,th0,I_int[4],th1,I_int[5]); //t1, grav0
    J[3] = Linear(theta,th0,I_int[6],th1,I_int[7]); //t1, grav1

    // Interpolate to chosen local gravity
    K[0] = pow(10.0,Linear(lgrav,grav0,log10(J[0]),grav1,log10(J[1]))); //t0
    K[1] = pow(10.0,Linear(lgrav,grav0,log10(J[2]),grav1,log10(J[3]))); //t1

    // Interpolate to chosen temperature
    L = pow(10.0,Linear(lt,t0,log10(K[0]),t1,log10(K[1])));

    // Set to zero at small angle
    if (cos_theta < 0.015629) L = 0;
    //cout << P << endl;

    //cout << I_temp[0] << " " << I_int[0] << " " << J[0] << " " << K[0] << " " << L << endl;
    /*
    if (isnan(L)) {
    	cout << theta << " " << th0 << " " << th1 << " " << i_mu << " " << n_mu << endl;
    	cout << cos_theta << " " << mexmcc.mccangl[i_mu] << " " << mexmcc.mccangl[n_mu] << endl;
    	cout << theta << " " << I_int[0] << " " << J[0] << " " << J[1] << " " << J[2] << " " << J[3] << endl;
    }
    */
    //cout << L << endl;

    return L;
}

// Calculate the final interpolated intensity
// This "new" version takes into account that the energy is really the ratio: E/kT
double McPHACC3new(double E, double cos_theta, double T, double lgrav, class LightCurve mexmcc){
	double delta, obl_approx, lt, th_spacing, theta, th_index, ener_spacing, first_ener, ener_index;
	double e0, e1, th0, th1, grav0, grav1, t0, t1;
	double I_temp[16], I_int[8], J[4], K[2], L(0.0);
	int i_f, n_f, i_lt, i_lgrav, n_lt, n_lgrav, i_mu, n_mu, first_inte;
	double mu0, mu1;

       
	


    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));

    // This part is purely for testing against the old mcphac file.
    /*
    lt = 6.3;
    T = pow(10,lt) * 1e-3/(Units::EV / Units::K_BOLTZ);
    lgrav = 14.25;
    */

    // Find the correct temperature range
    i_lt = (lt-5.1)/0.05; //if we need to load 1st temperature, i_lt = 0. this is discrete math
    if (i_lt < 1) i_lt = 1;


    i_lgrav = (lgrav-13.7)/0.1;
    if (i_lgrav < 1) i_lgrav = 1;

    //Find proper mu choice
    n_mu = 1;
    while (cos_theta > mexmcc.mccangl[n_mu] && n_mu < 50){
    	n_mu += 1;
    }
    i_mu = n_mu - 1;

    //    mu0 = mexmcc.mccangl[i_mu+1];
    //mu1 = mexmcc.mccangl[n_mu+1];

    //Find proper freqency choice

    ener_spacing = pow(10.0,0.0338);
    ener_index = (log10(E/T) + 1.30369)/0.0338;
    i_f = (int) ener_index;  

    // Now do a 4pt interpolation
    double evec[4];
    double ivec[4][4][4][4];
    double err;
    double muvec[4];
    double gvec[4];
    double tvec[4];
    
    int ii_mu(i_mu-1);
    if (ii_mu < 0)
      ii_mu = 0;
    if (ii_mu > 46)
      ii_mu = 46;

    for (int r(0); r<4; r++){
      tvec[r] =  5.1+0.05*(i_lt-1+r);
      //std::cout << "tvec[r]=" << tvec[r] << std::endl;
      for (int q(0); q<4; q++){
	gvec[q] = 13.7+0.1*(i_lgrav-1+q);
	//std::cout << "logg = " << gvec[q] << std::endl;
	for (int k(0); k<4; k++){
	  muvec[k] =  mexmcc.mccangl[ii_mu+k];
	  //std::cout << "muvec[k] = " << muvec[k] << std::endl;
	  // Interpolate over Energy for fixed Teff, gravity, and mu
	  for( int j(0); j<4; j++){
	    //evec[j] = pow(10,mexmcc.mcloget[i_f-1+j]);
	    evec[j] = mexmcc.mcloget[i_f-1+j];
	    first_inte = ((i_lt-1+r)*11 + i_lgrav-1+q) * 5000 + (i_f-1+j) * 50 + (ii_mu + k) +1 -1;
	    t0 = tvec[r];
	    ivec[r][q][k][j] = mexmcc.mccinte[first_inte]*pow(10.0,t0*3.0);
      
	    //std::cout << "evec[j] = " << evec[j] 
	    //	      << " ivec[j] = " << ivec[r][q][k][j]
	    //	      << std::endl;
	  }
	  I_int[k] = polint(evec,ivec[r][q][k],4,log10(E/T),&err);
	  //cout << "0 mu[k] = " << muvec[k]
	  //   <<" E=" << E/T << " New 4pt Interpolated: I = " << I_int[k] 
	  //   << " err = " << err << endl;
	}
	// Intepolate over mu for fixed Teff, gravity
	J[q] = polint(muvec,I_int,4,cos_theta,&err);
	//cout << " costheta = " << cos_theta << " Interpolated I = " << J[q]
	//   << " err = " << err << std::endl;
      }
      // Interpolate over logg for fixed Teff
      K[r] = polint(gvec,J,4,lgrav,&err);
      //cout << " logg = " << lgrav 
      //   << " Interpolated I = " 
      //   <<  K[r] 
      //   << " err = " << err << endl;      
    }

    L = polint(tvec,K,4,lt,&err);

    //std::cout << " log(T_eff) = " << lt 
    //	      << " Interpolated I = " << L
    //	      << " err = " << err << std::endl;


  
    //  cout << L << endl;
    
    if (isnan(L)) {
    	cout << theta << " " << th0 << " " << th1 << " " << i_mu << " " << n_mu << endl;
    	cout << cos_theta << " " << mexmcc.mccangl[i_mu] << " " << mexmcc.mccangl[n_mu] << endl;
    	cout << theta << " " << I_int[0] << " " << J[0] << " " << J[1] << " " << J[2] << " " << J[3] << endl;
    }

        return L;
}

// Calculate the final interpolated intensity
// This "new" version takes into account that the energy is really the ratio: E/kT
double NSXHnew(double E, double cos_theta, double T, double M, double R, class LightCurve mexmcc){
	double delta, obl_approx, lgrav, lt, th_spacing, theta, th_index, ener_spacing, first_ener, ener_index;
	double e0, e1, th0, th1, grav0, grav1, t0, t1;
	double I_temp[16], I_int[8], J[4], K[2], L(0.0);
	int i_f, n_f, i_lt, i_lgrav, n_lt, n_lgrav, i_mu, n_mu, first_inte;
	double mu0, mu1;

    //setting values of lt and lgrav based on input T, M, and R. Also sets to load using first mu value.   
    M = Units::nounits_to_cgs(M, Units::MASS);
    R = Units::nounits_to_cgs(R, Units::LENGTH);
    delta = 1 / sqrt(1 - (2 * Units::G * M / (R * Units::C * Units::C)));
    obl_approx = (1 + (-0.791 + 0.776 * mexmcc.para.mass_over_r) * pow(sin(mexmcc.para.omega_bar_sq),2) + (1.138 - 1.431 * mexmcc.para.mass_over_r) * pow(cos(mexmcc.para.omega_bar_sq),2));
    lgrav = log10(delta * Units::G * M / (R * R));
    //lgrav = log10(delta * Units::G * M / R / R * obl_approx);
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));
/*
    cout << endl; 
    cout << "temperature in log(K) is " << lt << endl;
    cout << "gravity in log(cgs units) is " << lgrav << endl;
    cout << "cos_theta is " << cos_theta << endl;
    cout << "energy is " << E << " keV" << endl;
    cout << "temperature is " << T  << " keV" << endl;
*/

    i_lt = (lt-5.1)/0.05; //if we need to load 1st temperature, i_lt = 0. this is discrete math
    n_lt = i_lt+1;
    t0 = 5.1+0.05*i_lt;
    t1 = t0+0.05;
    //    cout << i_lt << " " << n_lt << " " << t0 << " " << t1 << endl;
    //cout << t0 << " " << t1 << endl;

    i_lgrav = (lgrav-13.7)/0.1;
    n_lgrav = i_lgrav+1;
    grav0 = 13.7+0.1*i_lgrav;
    grav1 = grav0+0.1;
    //cout << grav0 << " " << grav1 << endl;
    
    //Find proper mu choice
    n_mu = 1;
    while (cos_theta > mexmcc.mccangl[n_mu] && n_mu < 67){
    	n_mu += 1;
    }

    i_mu = n_mu - 1;
    //    th0 = acos(mexmcc.mccangl[i_mu+1]);
    //th1 = acos(mexmcc.mccangl[n_mu+1]);
    //theta = acos(cos_theta);

    mu0 = mexmcc.mccangl[i_mu];
    mu1 = mexmcc.mccangl[n_mu];


    //cout << mexmcc.mccangl[i_mu] << " " << mexmcc.mccangl[n_mu] << " " << cos_theta << endl;
    //cout << th0 << " " << th1 << " " << acos(cos_theta) << endl;

    //Find proper freqency choice

    ener_spacing = pow(10.0,0.02);
    ener_index = (log10(E) + 1.32 )/0.02;

    i_f = (int) ener_index;
    n_f = i_f + 1;
    //    e0 = first_ener*pow(ener_spacing,i_f);
    //e1 = e0*ener_spacing;
    
    e0 = pow(10,mexmcc.mcloget[i_f]);
    e1 = pow(10,mexmcc.mcloget[n_f]);


    //cout << "log(e0) = " << mexmcc.mcloget[i_f] << std::endl;
    //  cout << "log(e0)=" << log10(e0) << " " << log10(e1) << " " << log10(E/T) << endl;
	  //cout << "e0=" << e0 << " " << e1 << " " << E/T << endl;
    
    first_inte = (i_lt*11 + i_lgrav) * 5000 + i_f * 50 + i_mu +1 -1;
    

    I_temp[0] = mexmcc.mccinte[first_inte];
    //cout << endl;
    //cout << "0: costheta = " <<  mexmcc.mccangl[i_mu] << " log(E0/kT) = " << log10(e0) << " I = " << I_temp[0] << std::endl;
    I_temp[2] = mexmcc.mccinte[first_inte+50];
    //cout << "2: costheta = " <<  mexmcc.mccangl[i_mu] << " log(E1/kT) = " << log10(e1) << " I = " << I_temp[2] << std::endl;




    I_temp[1] = mexmcc.mccinte[first_inte+1];
    //cout << "1: costheta = " <<  mexmcc.mccangl[n_mu] << " log(E0/kT) = " << log10(e0) << " I = " << I_temp[1] << std::endl;
    I_temp[3] = mexmcc.mccinte[first_inte+51];
    //cout << "3: costheta = " <<  mexmcc.mccangl[n_mu] << " log(E1/kT) = " << log10(e1) << " I = " << I_temp[3] << std::endl;

   


    I_temp[4] = mexmcc.mccinte[first_inte+5000];
    I_temp[5] = mexmcc.mccinte[first_inte+5001];
    I_temp[6] = mexmcc.mccinte[first_inte+5050];
    I_temp[7] = mexmcc.mccinte[first_inte+5051];

    I_temp[8] = mexmcc.mccinte[first_inte+55000];
    I_temp[9] = mexmcc.mccinte[first_inte+55001];
    I_temp[10] = mexmcc.mccinte[first_inte+55050];
    I_temp[11] = mexmcc.mccinte[first_inte+55051];

    I_temp[12] = mexmcc.mccinte[first_inte+60000];
    I_temp[14] = mexmcc.mccinte[first_inte+60001];
    I_temp[13] = mexmcc.mccinte[first_inte+60050];
    I_temp[15] = mexmcc.mccinte[first_inte+60051];

    I_int[0] = LogLinear(E, e0, I_temp[0], e1, I_temp[2]); //t0, grav0, th0
    //cout << "0 E Interpolated: I = " << I_int[0] << endl;

    I_int[1] = LogLinear(E, e0, I_temp[1], e1, I_temp[3]); //t0, grav0, th1
    //cout << "1 E Interpolated: I = " << I_int[1] << endl;
    

    I_int[2] = LogLinear(E, e0, I_temp[4], e1, I_temp[6]); //t0, grav1, th0
    I_int[3] = LogLinear(E, e0, I_temp[5], e1, I_temp[7]); //t0, grav1, th1
    I_int[4] = LogLinear(E, e0, I_temp[8], e1, I_temp[10]);//t1, grav0, th0
    I_int[5] = LogLinear(E, e0, I_temp[9], e1, I_temp[11]);//t1, grav0, th1
    I_int[6] = LogLinear(E, e0, I_temp[12], e1, I_temp[14]);//t1, grav1, th0
    I_int[7] = LogLinear(E, e0, I_temp[13], e1, I_temp[15]);//t1, grav1, th1
    
    //cout << I_int[0] << " " << I_int[1] << " " << I_int[2] << " " << I_int[3] << endl;

    // Interpolate to chosen mu
    J[0] = Linear(cos_theta,mu0,I_int[0],mu1,I_int[1]); //t0, grav0
    //cout << "J mu interpolated = " << J[0] << std::endl;
    J[1] = Linear(cos_theta,mu0,I_int[2],mu1,I_int[3]); //t0, grav1
    J[2] = Linear(cos_theta,mu0,I_int[4],mu1,I_int[5]); //t1, grav0
    J[3] = Linear(cos_theta,mu0,I_int[6],mu1,I_int[7]); //t1, grav1

    // Interpolate to chosen local gravity
    //cout << "grav0 = " << grav0 << " J0 = " << J[0] << endl;
    // cout << "grav1 = " << grav1 << " J1 = " << J[1] << endl;

    //lgrav = 14.25;

     // cout << " Interp to grav = " << lgrav << endl;

    K[0] = pow(10.0,Linear(lgrav,grav0,log10(J[0]),grav1,log10(J[1]))); //t0
    K[1] = pow(10.0,Linear(lgrav,grav0,log10(J[2]),grav1,log10(J[3]))); //t1

    //cout << "t0 = " << t0 << " K0 = " << K[0] << endl;
    //cout << "t1 = " << t1 << " K1= " << K[1] << endl;

    //lt = 6.3;

    // cout << "Interp to t=" << lt << endl;


    // Interpolate to chosen temperature
    L = pow(10.0,Linear(lt,t0,log10(K[0]),t1,log10(K[1])));

    // Set to zero at small angle
    if (cos_theta < 0.015629) L = 0;
    //cout << L << endl;

    //cout << I_temp[0] << " " << I_int[0] << " " << J[0] << " " << K[0] << " " << L << endl;
    /*
    if (isnan(L)) {
    	cout << theta << " " << th0 << " " << th1 << " " << i_mu << " " << n_mu << endl;
    	cout << cos_theta << " " << mexmcc.mccangl[i_mu] << " " << mexmcc.mccangl[n_mu] << endl;
    	cout << theta << " " << I_int[0] << " " << J[0] << " " << J[1] << " " << J[2] << " " << J[3] << endl;
    }
    */
        return L;
} // End of New NSX-H from Wynn Ho

// Calculate the final interpolated intensity
double McPHACC4(int E_dex, double cos_theta, double T, double M, double R, class LightCurve mexmcc){
	double delta, obl_approx, lgrav, lt, th_spacing, theta, th_index, ener_spacing, first_ener, ener_index;
	double th0, th1, grav0, grav1, t0, t1;
	double I_temp[16], I_int[8], J[4], K[2], L(0.0);
	int i_f, n_f, i_lt, i_lgrav, n_lt, n_lgrav, i_mu, n_mu, first_inte;

	//cout << "McPHACC4" << endl;
    //setting values of lt and lgrav based on input T, M, and R. Also sets to load using first mu value.   
    M = Units::nounits_to_cgs(M, Units::MASS);
    R = Units::nounits_to_cgs(R, Units::LENGTH);
    delta = 1 / sqrt(1 - (2 * Units::G * M / (R * Units::C * Units::C)));
    obl_approx = (1 + (-0.791 + 0.776 * mexmcc.para.mass_over_r) * pow(sin(mexmcc.para.omega_bar_sq),2) + (1.138 - 1.431 * mexmcc.para.mass_over_r) * pow(cos(mexmcc.para.omega_bar_sq),2));
    //lgrav = log10(delta * Units::G * M / (R * R));
    lgrav = log10(delta * Units::G * M / R / R * obl_approx);
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));
    //cout << "temperature in log(K) is " << lt << endl;
    //cout << "gravity in log(cgs units) is " << lgrav << endl;
    //cout << "cos_theta is " << cos_theta << endl;
    //cout << "energy is " << E << endl;
    //cout << mexmcc.para.omega << " " << mexmcc.para.theta << endl;

    i_lt = (lt-5.1)/0.05; //if we need to load 1st temperature, i_lt = 0. this is discrete math
    n_lt = i_lt+1;
    t0 = 5.1+0.05*i_lt;
    t1 = t0+0.05;
    //cout << t0 << " " << t1 << endl;

    i_lgrav = (lgrav-13.7)/0.1;
    n_lgrav = i_lgrav+1;
    grav0 = 13.7+0.1*i_lgrav;
    grav1 = grav0+0.1;
    //cout << grav0 << " " << grav1 << endl;
    
    //Find proper mu choice
    //cout << "finding mu" << endl;
    n_mu = 1;
    while (cos_theta > mexmcc.mccangl[n_mu] && n_mu < 49){
    	n_mu += 1;
    }

    i_mu = n_mu - 1;
    th0 = acos(mexmcc.mccangl[i_mu]);
    th1 = acos(mexmcc.mccangl[n_mu]);
    theta = acos(cos_theta);
    //cout << mexmcc.mccangl[i_mu] << " " << mexmcc.mccangl[n_mu] << " " << cos_theta << endl;
    //cout << th0 << " " << th1 << " " << acos(cos_theta) << endl;
    //cout << mexmcc.mccangl[i_mu] << " " << mexmcc.mccangl[n_mu] << " " << cos_theta << endl;


    first_inte = (i_lt*11 + i_lgrav) * 5000 + E_dex * 50 + i_mu;
    //cout << i_lt << " " << i_lgrav << " " << E_dex << " " << i_mu << " " << first_inte << endl;

    I_int[0] = mexmcc.mccinte[first_inte]*pow(10.0,t0*3.0); //t0, grav0, th0
    I_int[1] = mexmcc.mccinte[first_inte+1]*pow(10.0,t0*3.0); //t0, grav0, th1
    I_int[2] = mexmcc.mccinte[first_inte+5000]*pow(10.0,t0*3.0); //t0, grav1, th0
    I_int[3] = mexmcc.mccinte[first_inte+5001]*pow(10.0,t0*3.0); //t0, grav1, th1
    I_int[4] = mexmcc.mccinte[first_inte+55000]*pow(10.0,t1*3.0);//t1, grav0, th0
    I_int[5] = mexmcc.mccinte[first_inte+55001]*pow(10.0,t1*3.0);//t1, grav0, th1
    I_int[6] = mexmcc.mccinte[first_inte+60000]*pow(10.0,t1*3.0);//t1, grav1, th0
    I_int[7] = mexmcc.mccinte[first_inte+60001]*pow(10.0,t1*3.0);//t1, grav1, th1
    
    // Interpolate to chosen mu
    J[0] = Linear(theta,th0,I_int[0],th1,I_int[1]); //t0, grav0
    J[1] = Linear(theta,th0,I_int[2],th1,I_int[3]); //t0, grav1
    J[2] = Linear(theta,th0,I_int[4],th1,I_int[5]); //t1, grav0
    J[3] = Linear(theta,th0,I_int[6],th1,I_int[7]); //t1, grav1

    // Interpolate to chosen local gravity
    K[0] = pow(10.0,Linear(lgrav,grav0,log10(J[0]),grav1,log10(J[1]))); //t0
    K[1] = pow(10.0,Linear(lgrav,grav0,log10(J[2]),grav1,log10(J[3]))); //t1

    // Interpolate to chosen temperature
    L = pow(10.0,Linear(lt,t0,log10(K[0]),t1,log10(K[1])));

    // Set to zero at small angle
    if (cos_theta < 0.015629) L = 0;

    //cout << I_temp[0] << " " << I_int[0] << " " << J[0] << " " << K[0] << " " << L << endl;
    /*
    if (isnan(L)) {
    	cout << theta << " " << th0 << " " << th1 << " " << i_mu << " " << n_mu << endl;
    	cout << cos_theta << " " << mexmcc.mccangl[i_mu] << " " << mexmcc.mccangl[n_mu] << endl;
    	cout << theta << " " << I_int[0] << " " << J[0] << " " << J[1] << " " << J[2] << " " << J[3] << endl;
    }
    */
    //cout << L << endl;

    return L;
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
	//unsigned int n_steps(900); // total number of steps
	unsigned int n_steps(800);
	// This number of steps (100) is optimized for Delta(E) = 0.3 keV
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

    if (model == 5){ //Hydrogen
        for (int j = 1; j <= n_steps; j++){            
            current_e = step_size*(j-1+1/2) + E1;
            e_l = (current_e-step_size/2);
            e_m = (current_e);
            e_u = (current_e+step_size/2);
            current_n = step_size * (1/3*NSXH(e_l,cos_theta)/e_l + 4/3*NSXH(e_m,cos_theta)/e_m + 1/3*NSXH(e_u,cos_theta)/e_u);
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
    int ener_size(0);         // number of energy choices
    //int e_dex;                  // energy index of current step
    //double current_e;           // central energy of current step
    //double current_n;           // integrated flux in current step
    double flux(0.0);           // total integrated flux

    if (model == 5 || model == 8 || model == 9 || model == 10){
    	ener_size = 100;
    	//cout << "integrating ener_size = " << ener_size << endl;
    }
    if (model == 3 || model == 4){
    	ener_size = 125;
    }

    for (int m = 0; m < ener_size; m++) {
        if (E1 >= Es[m]) {
            e1_dex = m;
        }
        if (E2 >= Es[m]) {
            e2_dex = m;
        }
		//cout << Es[m] << " ";
    }

    //cout << endl;

    //calculate number of energy points within band  
    n_steps = e2_dex - e1_dex;
    //cout << E1 << " " << E2 << " " << e1_dex << " " << e2_dex << " " << n_steps << endl;
    
    /*
    if (n_steps >= 4){
    	cout << E1 << " " << E2 << " " << e1_dex << " " << e2_dex << " " << n_steps << endl;
    }
    */


    if (model == 3){ //Hydrogen
        if (n_steps == 0){ // zero energy points within bandwidth: (4.1.3) one trapzoid
            //cout << "0 steps" << endl;
            flux = (E2 - E1) / 2 * (Hydrogen(E1,cos_theta) / E1 + Hydrogen(E2,cos_theta) / E2);
        }
        if (n_steps == 1){ // one energy points within bandwidth: (4.1.3) two trapzoids
            //cout << "1 steps" << endl;
            int e_dex = e1_dex+1; // index of the energy point
            double e_m = F[e_dex] * Units::H_PLANCK / Units::EV / 1E3; // energy point in keV
            flux = (e_m - E1) / 2 * (Hydrogen(E1,cos_theta) / E1 + Hydrogen2(e_dex,cos_theta) / e_m);  // first trapezoid
            flux += (E2 - e_m) / 2 * (Hydrogen2(e_dex,cos_theta) / e_m + Hydrogen(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 2){ // two energy points within bandwidth: (4.1.3) three trapzoids
            //cout << "2 steps" << endl;
            int el_dex = e1_dex+1; // index of the first energy point within bandwidth
            int eu_dex = e2_dex;   // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = log(e_u)-log(e_l); // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (Hydrogen(E1,cos_theta) / E1 + Hydrogen2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (Hydrogen2(el_dex,cos_theta) + Hydrogen2(eu_dex,cos_theta)) / 2;                // middle trapezoid, exact
            flux += (E2 - e_u) / 2 * (Hydrogen2(eu_dex,cos_theta) / e_u + Hydrogen(E2,cos_theta) / E2); // last trapezoid
        }
        if (n_steps == 3){ // three energy points within bandwidth: (4.1.3, 4.1.4) Simpson's + two trapezoids
            //cout << "3 steps" << endl;
            int el_dex = e1_dex+1; // index of the first energy point within bandwidth
            int em_dex = e1_dex+2; // index of the middle energy point within bandwidth
            int eu_dex = e2_dex;   // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 2; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (Hydrogen(E1,cos_theta) / E1 + Hydrogen2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (Hydrogen2(el_dex,cos_theta) + Hydrogen2(em_dex,cos_theta) * 4 + Hydrogen2(eu_dex,cos_theta)) / 3; // Simpson's, exact
            flux += (E2 - e_u) / 2 * (Hydrogen2(eu_dex,cos_theta) / e_u + Hydrogen(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 4){ // four energy points within bandwidth: (4.1.3, 4.1.5) 8/3 Simpson's + two trapezoids
            //cout << "4 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 3; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (Hydrogen(E1,cos_theta) / E1 + Hydrogen2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (Hydrogen2(el_dex,cos_theta) + Hydrogen2(em1_dex,cos_theta) * 3 + Hydrogen2(em2_dex,cos_theta) * 3 + Hydrogen2(eu_dex,cos_theta)) * 3 / 8; // Simpson's 3/8, exact             
            flux += (E2 - e_u) / 2 * (Hydrogen2(eu_dex,cos_theta) / e_u + Hydrogen(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 5){ // five energy points within bandwidth: (4.1.3, 4.1.13) Simpson's "4,2" + two trapezoids  
            //cout << "5 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int em3_dex = e1_dex+4; // index of the fourth energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 4; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (Hydrogen(E1,cos_theta) / E1 + Hydrogen2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (Hydrogen2(el_dex,cos_theta) + Hydrogen2(em1_dex,cos_theta) * 4 + Hydrogen2(em2_dex,cos_theta) * 2 + Hydrogen2(em3_dex,cos_theta) * 4 + Hydrogen2(eu_dex,cos_theta)) / 3; // Simpson's "4,2", exact             
            flux += (E2 - e_u) / 2 * (Hydrogen2(eu_dex,cos_theta) / e_u + Hydrogen(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps >= 6){ // six or more energy points within bandwidth: (4.1.3, 4.1.14) Simpson's cubic + two trapezoids
            //cout << "6 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int el1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int el2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            int eu1_dex = e2_dex-1; // index of the second last energy point within bandwidth
            int eu2_dex = e2_dex-2; // index of the third last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / (n_steps-1); // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (Hydrogen(E1,cos_theta) / E1 + Hydrogen2(el_dex,cos_theta) / e_l);  // first trapezoid
            // Simpson's cubic, exact
            flux += h * (Hydrogen2(el_dex,cos_theta) * 9 + Hydrogen2(el1_dex,cos_theta) * 28 + Hydrogen2(el2_dex,cos_theta) * 23) / 24; // first three coefficients
            for (int m = 1; m <= n_steps-6; m++) flux += h * (Hydrogen2(el_dex+m+3,cos_theta)); // middle coefficients
            flux += h * (Hydrogen2(eu2_dex,cos_theta) * 23 + Hydrogen2(eu1_dex,cos_theta) * 28 + Hydrogen2(eu_dex,cos_theta) * 9) / 24; // last three coefficients
            // end Simpson's cubic
            flux += (E2 - e_u) / 2 * (Hydrogen2(eu_dex,cos_theta) / e_u + Hydrogen(E2,cos_theta) / E2); // second trapezoid                
        }
    }



    if (model == 4){ //Helium
        if (n_steps == 0){ // zero energy points within bandwidth: (4.1.3) one trapzoid
            //cout << "0 steps" << endl;
            flux = (E2 - E1) / 2 * (Helium(E1,cos_theta) / E1 + Helium(E2,cos_theta) / E2);
        }
        if (n_steps == 1){ // one energy points within bandwidth: (4.1.3) two trapzoids
            //cout << "1 steps" << endl;
            int e_dex = e1_dex+1; // index of the energy point
            double e_m = F[e_dex] * Units::H_PLANCK / Units::EV / 1E3; // energy point in keV
            flux = (e_m - E1) / 2 * (Helium(E1,cos_theta) / E1 + Helium2(e_dex,cos_theta) / e_m);  // first trapezoid
            flux += (E2 - e_m) / 2 * (Helium2(e_dex,cos_theta) / e_m + Helium(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 2){ // two energy points within bandwidth: (4.1.3) three trapzoids
            //cout << "2 steps" << endl;
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
            //cout << "3 steps" << endl;
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
            //cout << "4 steps" << endl;
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
            //cout << "5 steps" << endl;
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
            //cout << "6 steps" << endl;
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
    }
    

    if (model == 5){ //NSX Hydrogen
        if (n_steps == 0){ // zero energy points within bandwidth: (4.1.3) one trapzoid
            //cout << "0 steps" << endl;
            flux = (E2 - E1) / 2 * (NSXH(E1,cos_theta) / E1 + NSXH(E2,cos_theta) / E2);
        }
        if (n_steps == 1){ // one energy points within bandwidth: (4.1.3) two trapzoids
            //cout << "1 steps" << endl;
            int e_dex = e1_dex+1; // index of the energy point
            double e_m = F[e_dex] * Units::H_PLANCK / Units::EV / 1E3; // energy point in keV
            flux = (e_m - E1) / 2 * (NSXH(E1,cos_theta) / E1 + NSXH2(e_dex,cos_theta) / e_m);  // first trapezoid
            flux += (E2 - e_m) / 2 * (NSXH2(e_dex,cos_theta) / e_m + NSXH(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 2){ // two energy points within bandwidth: (4.1.3) three trapzoids
            //cout << "2 steps" << endl;
            int el_dex = e1_dex+1; // index of the first energy point within bandwidth
            int eu_dex = e2_dex;   // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = log(e_u)-log(e_l); // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (NSXH(E1,cos_theta) / E1 + NSXH2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (NSXH2(el_dex,cos_theta) + NSXH2(eu_dex,cos_theta)) / 2;                // middle trapezoid, exact
            flux += (E2 - e_u) / 2 * (NSXH2(eu_dex,cos_theta) / e_u + NSXH(E2,cos_theta) / E2); // last trapezoid
        }
        if (n_steps == 3){ // three energy points within bandwidth: (4.1.3, 4.1.4) Simpson's + two trapezoids
            //cout << "3 steps" << endl;
            int el_dex = e1_dex+1; // index of the first energy point within bandwidth
            int em_dex = e1_dex+2; // index of the middle energy point within bandwidth
            int eu_dex = e2_dex;   // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 2; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (NSXH(E1,cos_theta) / E1 + NSXH2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (NSXH2(el_dex,cos_theta) + NSXH2(em_dex,cos_theta) * 4 + NSXH2(eu_dex,cos_theta)) / 3; // Simpson's, exact
            flux += (E2 - e_u) / 2 * (NSXH2(eu_dex,cos_theta) / e_u + NSXH(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 4){ // four energy points within bandwidth: (4.1.3, 4.1.5) 8/3 Simpson's + two trapezoids
            //cout << "4 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 3; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (NSXH(E1,cos_theta) / E1 + NSXH2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (NSXH2(el_dex,cos_theta) + NSXH2(em1_dex,cos_theta) * 3 + NSXH2(em2_dex,cos_theta) * 3 + NSXH2(eu_dex,cos_theta)) * 3 / 8; // Simpson's 3/8, exact             
            flux += (E2 - e_u) / 2 * (NSXH2(eu_dex,cos_theta) / e_u + NSXH(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 5){ // five energy points within bandwidth: (4.1.3, 4.1.13) Simpson's "4,2" + two trapezoids  
            //cout << "5 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int em3_dex = e1_dex+4; // index of the fourth energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 4; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (NSXH(E1,cos_theta) / E1 + NSXH2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (NSXH2(el_dex,cos_theta) + NSXH2(em1_dex,cos_theta) * 4 + NSXH2(em2_dex,cos_theta) * 2 + NSXH2(em3_dex,cos_theta) * 4 + NSXH2(eu_dex,cos_theta)) / 3; // Simpson's "4,2", exact             
            flux += (E2 - e_u) / 2 * (NSXH2(eu_dex,cos_theta) / e_u + NSXH(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps >= 6){ // six or more energy points within bandwidth: (4.1.3, 4.1.14) Simpson's cubic + two trapezoids
            //cout << "6 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int el1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int el2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            int eu1_dex = e2_dex-1; // index of the second last energy point within bandwidth
            int eu2_dex = e2_dex-2; // index of the third last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / (n_steps-1); // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (NSXH(E1,cos_theta) / E1 + NSXH2(el_dex,cos_theta) / e_l);  // first trapezoid
            // Simpson's cubic, exact
            flux += h * (NSXH2(el_dex,cos_theta) * 9 + NSXH2(el1_dex,cos_theta) * 28 + NSXH2(el2_dex,cos_theta) * 23) / 24; // first three coefficients
            for (int m = 1; m <= n_steps-6; m++) flux += h * (NSXH2(el_dex+m+3,cos_theta)); // middle coefficients
            flux += h * (NSXH2(eu2_dex,cos_theta) * 23 + NSXH2(eu1_dex,cos_theta) * 28 + NSXH2(eu_dex,cos_theta) * 9) / 24; // last three coefficients
            // end Simpson's cubic
            flux += (E2 - e_u) / 2 * (NSXH2(eu_dex,cos_theta) / e_u + NSXH(E2,cos_theta) / E2); // second trapezoid                
        }
    }


    if (model == 8){ //NICER's McPHAC
        if (n_steps == 0){ // zero energy points within bandwidth: (4.1.3) one trapzoid
            //cout << "0 steps" << endl;
            flux = (E2 - E1) / 2 * (McPHAC(E1,cos_theta) / E1 + McPHAC(E2,cos_theta) / E2);
        }
        if (n_steps == 1){ // one energy points within bandwidth: (4.1.3) two trapzoids
            //cout << "1 steps" << endl;
            int e_dex = e1_dex+1; // index of the energy point
            double e_m = F[e_dex] * Units::H_PLANCK / Units::EV / 1E3; // energy point in keV
            flux = (e_m - E1) / 2 * (McPHAC(E1,cos_theta) / E1 + McPHAC2(e_dex,cos_theta) / e_m);  // first trapezoid
            flux += (E2 - e_m) / 2 * (McPHAC2(e_dex,cos_theta) / e_m + McPHAC(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 2){ // two energy points within bandwidth: (4.1.3) three trapzoids
            //cout << "2 steps" << endl;
            int el_dex = e1_dex+1; // index of the first energy point within bandwidth
            int eu_dex = e2_dex;   // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = log(e_u)-log(e_l); // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (McPHAC(E1,cos_theta) / E1 + McPHAC2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (McPHAC2(el_dex,cos_theta) + McPHAC2(eu_dex,cos_theta)) / 2;                // middle trapezoid, exact
            flux += (E2 - e_u) / 2 * (McPHAC2(eu_dex,cos_theta) / e_u + McPHAC(E2,cos_theta) / E2); // last trapezoid
        }
        if (n_steps == 3){ // three energy points within bandwidth: (4.1.3, 4.1.4) Simpson's + two trapezoids
            //cout << "3 steps" << endl;
            int el_dex = e1_dex+1; // index of the first energy point within bandwidth
            int em_dex = e1_dex+2; // index of the middle energy point within bandwidth
            int eu_dex = e2_dex;   // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 2; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (McPHAC(E1,cos_theta) / E1 + McPHAC2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (McPHAC2(el_dex,cos_theta) + McPHAC2(em_dex,cos_theta) * 4 + McPHAC2(eu_dex,cos_theta)) / 3; // Simpson's, exact
            flux += (E2 - e_u) / 2 * (McPHAC2(eu_dex,cos_theta) / e_u + McPHAC(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 4){ // four energy points within bandwidth: (4.1.3, 4.1.5) 8/3 Simpson's + two trapezoids
            //cout << "4 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 3; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (McPHAC(E1,cos_theta) / E1 + McPHAC2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (McPHAC2(el_dex,cos_theta) + McPHAC2(em1_dex,cos_theta) * 3 + McPHAC2(em2_dex,cos_theta) * 3 + McPHAC2(eu_dex,cos_theta)) * 3 / 8; // Simpson's 3/8, exact             
            flux += (E2 - e_u) / 2 * (McPHAC2(eu_dex,cos_theta) / e_u + McPHAC(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 5){ // five energy points within bandwidth: (4.1.3, 4.1.13) Simpson's "4,2" + two trapezoids  
            //cout << "5 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int em3_dex = e1_dex+4; // index of the fourth energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 4; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (McPHAC(E1,cos_theta) / E1 + McPHAC2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (McPHAC2(el_dex,cos_theta) + McPHAC2(em1_dex,cos_theta) * 4 + McPHAC2(em2_dex,cos_theta) * 2 + McPHAC2(em3_dex,cos_theta) * 4 + McPHAC2(eu_dex,cos_theta)) / 3; // Simpson's "4,2", exact             
            flux += (E2 - e_u) / 2 * (McPHAC2(eu_dex,cos_theta) / e_u + McPHAC(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps >= 6){ // six or more energy points within bandwidth: (4.1.3, 4.1.14) Simpson's cubic + two trapezoids
            //cout << "6 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int el1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int el2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            int eu1_dex = e2_dex-1; // index of the second last energy point within bandwidth
            int eu2_dex = e2_dex-2; // index of the third last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / (n_steps-1); // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (McPHAC(E1,cos_theta) / E1 + McPHAC2(el_dex,cos_theta) / e_l);  // first trapezoid
            // Simpson's cubic, exact
            flux += h * (McPHAC2(el_dex,cos_theta) * 9 + McPHAC2(el1_dex,cos_theta) * 28 + McPHAC2(el2_dex,cos_theta) * 23) / 24; // first three coefficients
            for (int m = 1; m <= n_steps-6; m++) flux += h * (McPHAC2(el_dex+m+3,cos_theta)); // middle coefficients
            flux += h * (McPHAC2(eu2_dex,cos_theta) * 23 + McPHAC2(eu1_dex,cos_theta) * 28 + McPHAC2(eu_dex,cos_theta) * 9) / 24; // last three coefficients
            // end Simpson's cubic
            flux += (E2 - e_u) / 2 * (McPHAC2(eu_dex,cos_theta) / e_u + McPHAC(E2,cos_theta) / E2); // second trapezoid                
        }
    }


    if (model == 9){ //new NSX Helium
        if (n_steps == 0){ // zero energy points within bandwidth: (4.1.3) one trapzoid
            //cout << "0 steps" << endl;
            flux = (E2 - E1) / 2 * (NSXHe(E1,cos_theta) / E1 + NSXHe(E2,cos_theta) / E2);
        }
        if (n_steps == 1){ // one energy points within bandwidth: (4.1.3) two trapzoids
            //cout << "1 steps" << endl;
            int e_dex = e1_dex+1; // index of the energy point
            double e_m = F[e_dex] * Units::H_PLANCK / Units::EV / 1E3; // energy point in keV
            flux = (e_m - E1) / 2 * (NSXHe(E1,cos_theta) / E1 + NSXHe2(e_dex,cos_theta) / e_m);  // first trapezoid
            flux += (E2 - e_m) / 2 * (NSXHe2(e_dex,cos_theta) / e_m + NSXHe(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 2){ // two energy points within bandwidth: (4.1.3) three trapzoids
            //cout << "2 steps" << endl;
            int el_dex = e1_dex+1; // index of the first energy point within bandwidth
            int eu_dex = e2_dex;   // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = log(e_u)-log(e_l); // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (NSXHe(E1,cos_theta) / E1 + NSXHe2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (NSXHe2(el_dex,cos_theta) + NSXHe2(eu_dex,cos_theta)) / 2;                // middle trapezoid, exact
            flux += (E2 - e_u) / 2 * (NSXHe2(eu_dex,cos_theta) / e_u + NSXHe(E2,cos_theta) / E2); // last trapezoid
        }
        if (n_steps == 3){ // three energy points within bandwidth: (4.1.3, 4.1.4) Simpson's + two trapezoids
            //cout << "3 steps" << endl;
            int el_dex = e1_dex+1; // index of the first energy point within bandwidth
            int em_dex = e1_dex+2; // index of the middle energy point within bandwidth
            int eu_dex = e2_dex;   // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 2; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (NSXHe(E1,cos_theta) / E1 + NSXHe2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (NSXHe2(el_dex,cos_theta) + NSXHe2(em_dex,cos_theta) * 4 + NSXHe2(eu_dex,cos_theta)) / 3; // Simpson's, exact
            flux += (E2 - e_u) / 2 * (NSXHe2(eu_dex,cos_theta) / e_u + NSXHe(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 4){ // four energy points within bandwidth: (4.1.3, 4.1.5) 8/3 Simpson's + two trapezoids
            //cout << "4 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 3; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (NSXHe(E1,cos_theta) / E1 + NSXHe2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (NSXHe2(el_dex,cos_theta) + NSXHe2(em1_dex,cos_theta) * 3 + NSXHe2(em2_dex,cos_theta) * 3 + NSXHe2(eu_dex,cos_theta)) * 3 / 8; // Simpson's 3/8, exact             
            flux += (E2 - e_u) / 2 * (NSXHe2(eu_dex,cos_theta) / e_u + NSXHe(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 5){ // five energy points within bandwidth: (4.1.3, 4.1.13) Simpson's "4,2" + two trapezoids  
            //cout << "5 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int em3_dex = e1_dex+4; // index of the fourth energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 4; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (NSXHe(E1,cos_theta) / E1 + NSXHe2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (NSXHe2(el_dex,cos_theta) + NSXHe2(em1_dex,cos_theta) * 4 + NSXHe2(em2_dex,cos_theta) * 2 + NSXHe2(em3_dex,cos_theta) * 4 + NSXHe2(eu_dex,cos_theta)) / 3; // Simpson's "4,2", exact             
            flux += (E2 - e_u) / 2 * (NSXHe2(eu_dex,cos_theta) / e_u + NSXHe(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps >= 6){ // six or more energy points within bandwidth: (4.1.3, 4.1.14) Simpson's cubic + two trapezoids
            //cout << "6 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int el1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int el2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            int eu1_dex = e2_dex-1; // index of the second last energy point within bandwidth
            int eu2_dex = e2_dex-2; // index of the third last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / (n_steps-1); // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (NSXHe(E1,cos_theta) / E1 + NSXHe2(el_dex,cos_theta) / e_l);  // first trapezoid
            // Simpson's cubic, exact
            flux += h * (NSXHe2(el_dex,cos_theta) * 9 + NSXHe2(el1_dex,cos_theta) * 28 + NSXHe2(el2_dex,cos_theta) * 23) / 24; // first three coefficients
            for (int m = 1; m <= n_steps-6; m++) flux += h * (NSXHe2(el_dex+m+3,cos_theta)); // middle coefficients
            flux += h * (NSXHe2(eu2_dex,cos_theta) * 23 + NSXHe2(eu1_dex,cos_theta) * 28 + NSXHe2(eu_dex,cos_theta) * 9) / 24; // last three coefficients
            // end Simpson's cubic
            flux += (E2 - e_u) / 2 * (NSXHe2(eu_dex,cos_theta) / e_u + NSXHe(E2,cos_theta) / E2); // second trapezoid                
        }
    }

    
    if (model == 10){ //Cole's McPHAC Spot routines
        if (n_steps == 0){ // zero energy points within bandwidth: (4.1.3) one trapzoid
            //cout << "0 steps" << endl;
            flux = (E2 - E1) / 2 * (McPHACC(E1,cos_theta) / E1 + McPHACC(E2,cos_theta) / E2);
        }
        if (n_steps == 1){ // one energy points within bandwidth: (4.1.3) two trapzoids
            //cout << "1 steps" << endl;
            int e_dex = e1_dex+1; // index of the energy point
            double e_m = F[e_dex] * Units::H_PLANCK / Units::EV / 1E3; // energy point in keV
            flux = (e_m - E1) / 2 * (McPHACC(E1,cos_theta) / E1 + McPHACC2(e_dex,cos_theta) / e_m);  // first trapezoid
            flux += (E2 - e_m) / 2 * (McPHACC2(e_dex,cos_theta) / e_m + McPHACC(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 2){ // two energy points within bandwidth: (4.1.3) three trapzoids
            //cout << "2 steps" << endl;
            int el_dex = e1_dex+1; // index of the first energy point within bandwidth
            int eu_dex = e2_dex;   // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = log(e_u)-log(e_l); // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (McPHACC(E1,cos_theta) / E1 + McPHACC2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (McPHACC2(el_dex,cos_theta) + McPHACC2(eu_dex,cos_theta)) / 2;                // middle trapezoid, exact
            flux += (E2 - e_u) / 2 * (McPHACC2(eu_dex,cos_theta) / e_u + McPHACC(E2,cos_theta) / E2); // last trapezoid
        }
        if (n_steps == 3){ // three energy points within bandwidth: (4.1.3, 4.1.4) Simpson's + two trapezoids
            //cout << "3 steps" << endl;
            int el_dex = e1_dex+1; // index of the first energy point within bandwidth
            int em_dex = e1_dex+2; // index of the middle energy point within bandwidth
            int eu_dex = e2_dex;   // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 2; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (McPHACC(E1,cos_theta) / E1 + McPHACC2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (McPHACC2(el_dex,cos_theta) + McPHACC2(em_dex,cos_theta) * 4 + McPHACC2(eu_dex,cos_theta)) / 3; // Simpson's, exact
            flux += (E2 - e_u) / 2 * (McPHACC2(eu_dex,cos_theta) / e_u + McPHACC(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 4){ // four energy points within bandwidth: (4.1.3, 4.1.5) 8/3 Simpson's + two trapezoids
            //cout << "4 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 3; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (McPHACC(E1,cos_theta) / E1 + McPHACC2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (McPHACC2(el_dex,cos_theta) + McPHACC2(em1_dex,cos_theta) * 3 + McPHACC2(em2_dex,cos_theta) * 3 + McPHACC2(eu_dex,cos_theta)) * 3 / 8; // Simpson's 3/8, exact             
            flux += (E2 - e_u) / 2 * (McPHACC2(eu_dex,cos_theta) / e_u + McPHACC(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps == 5){ // five energy points within bandwidth: (4.1.3, 4.1.13) Simpson's "4,2" + two trapezoids  
            //cout << "5 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int em3_dex = e1_dex+4; // index of the fourth energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / 4; // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (McPHACC(E1,cos_theta) / E1 + McPHACC2(el_dex,cos_theta) / e_l);  // first trapezoid
            flux += h * (McPHACC2(el_dex,cos_theta) + McPHACC2(em1_dex,cos_theta) * 4 + McPHACC2(em2_dex,cos_theta) * 2 + McPHACC2(em3_dex,cos_theta) * 4 + McPHACC2(eu_dex,cos_theta)) / 3; // Simpson's "4,2", exact             
            flux += (E2 - e_u) / 2 * (McPHACC2(eu_dex,cos_theta) / e_u + McPHACC(E2,cos_theta) / E2); // second trapezoid
        }
        if (n_steps >= 6){ // six or more energy points within bandwidth: (4.1.3, 4.1.14) Simpson's cubic + two trapezoids
            //cout << "6 steps" << endl;
            int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
            int el1_dex = e1_dex+2; // index of the second energy point within bandwidth
            int el2_dex = e1_dex+3; // index of the third energy point within bandwidth        
            int eu_dex = e2_dex;    // index of the last energy point within bandwidth
            int eu1_dex = e2_dex-1; // index of the second last energy point within bandwidth
            int eu2_dex = e2_dex-2; // index of the third last energy point within bandwidth
            double e_l = F[el_dex] * Units::H_PLANCK / Units::EV / 1E3; // first energy point in keV
            double e_u = F[eu_dex] * Units::H_PLANCK / Units::EV / 1E3; // last energy point in keV
            double h = (log(e_u)-log(e_l)) / (n_steps-1); // difference between log-spaced energy points
            flux = (e_l - E1) / 2 * (McPHACC(E1,cos_theta) / E1 + McPHACC2(el_dex,cos_theta) / e_l);  // first trapezoid
            // Simpson's cubic, exact
            flux += h * (McPHACC2(el_dex,cos_theta) * 9 + McPHACC2(el1_dex,cos_theta) * 28 + McPHACC2(el2_dex,cos_theta) * 23) / 24; // first three coefficients
            for (int m = 1; m <= n_steps-6; m++) flux += h * (McPHACC2(el_dex+m+3,cos_theta)); // middle coefficients
            flux += h * (McPHACC2(eu2_dex,cos_theta) * 23 + McPHACC2(eu1_dex,cos_theta) * 28 + McPHACC2(eu_dex,cos_theta) * 9) / 24; // last three coefficients
            // end Simpson's cubic
            flux += (E2 - e_u) / 2 * (McPHACC2(eu_dex,cos_theta) / e_u + McPHACC(E2,cos_theta) / E2); // second trapezoid                
        }
    }
    

    flux = flux/Units::H_PLANCK;
    return flux;
}



// This is version makes use of Cole's version of McPhac
double AtmosEBandFlux3( unsigned int model, double cos_theta, double T, double M, double R, double E1, double E2, class LightCurve mexmcc){

    int e1_dex(0);              // largest energy index that's smaller than E1
    int e2_dex(0);              // largest energy index that's smaller than E2
    int n_steps(0);             // number of energy points within band
    int ener_size(0);         // number of energy choices
    //int e_dex;                  // energy index of current step
    //double current_e;           // central energy of current step
    //double current_n;           // integrated flux in current step
    double flux(0.0);           // total integrated flux

    ener_size = 125;

	double ener_spacing = pow(10.0,0.0338);
    double first_ener = 0.004969469;
    double ener_index = log(E1 / first_ener) / log(ener_spacing);
    e1_dex = (int) ener_index;
    ener_index = log(E2 / first_ener) / log(ener_spacing);
    e2_dex = (int) ener_index;

    n_steps = e2_dex - e1_dex;

    if (n_steps == 0){ // zero energy points within bandwidth: (4.1.3) one trapzoid
        //cout << "0 steps" << endl;
        flux = (E2 - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2);
    }
    if (n_steps == 1){ // one energy points within bandwidth: (4.1.3) two trapzoids
        //cout << "1 steps" << endl;
        int e_dex = e1_dex+1; // index of the energy point
        double e_m = first_ener*pow(ener_spacing,e_dex); // energy point in keV
        flux = (e_m - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC4(e1_dex,cos_theta, T, M, R, mexmcc) / e_m);  // first trapezoid
        flux += (E2 - e_m) / 2 * (McPHACC4(e_dex,cos_theta, T, M, R, mexmcc) / e_m + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2); // second trapezoid
    }
    if (n_steps == 2){ // two energy points within bandwidth: (4.1.3) three trapzoids
        //cout << "2 steps" << endl;
        int el_dex = e1_dex+1; // index of the first energy point within bandwidth
        int eu_dex = e2_dex;   // index of the last energy point within bandwidth
        double e_l = first_ener*pow(ener_spacing,el_dex); // first energy point in keV
        double e_u = first_ener*pow(ener_spacing,eu_dex); // last energy point in keV
        double h = log(e_u)-log(e_l); // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) / e_l);  // first trapezoid
        flux += h * (McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) + McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc)) / 2;                // middle trapezoid, exact
        flux += (E2 - e_u) / 2 * (McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc) / e_u + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2); // last trapezoid
    }
    if (n_steps == 3){ // three energy points within bandwidth: (4.1.3, 4.1.4) Simpson's + two trapezoids
        //cout << "3 steps" << endl;
        int el_dex = e1_dex+1; // index of the first energy point within bandwidth
        int em_dex = e1_dex+2; // index of the middle energy point within bandwidth
        int eu_dex = e2_dex;   // index of the last energy point within bandwidth
        double e_l = first_ener*pow(ener_spacing,el_dex); // first energy point in keV
        double e_u = first_ener*pow(ener_spacing,eu_dex); // last energy point in keV
        double h = (log(e_u)-log(e_l)) / 2; // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) / e_l);  // first trapezoid
        flux += h * (McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) + McPHACC4(em_dex,cos_theta, T, M, R, mexmcc) * 4 + McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc)) / 3; // Simpson's, exact
        flux += (E2 - e_u) / 2 * (McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc) / e_u + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2); // second trapezoid
    }
    if (n_steps == 4){ // four energy points within bandwidth: (4.1.3, 4.1.5) 8/3 Simpson's + two trapezoids
        //cout << "4 steps" << endl;
        int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
        int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
        int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
        int eu_dex = e2_dex;    // index of the last energy point within bandwidth
        double e_l = first_ener*pow(ener_spacing,el_dex); // first energy point in keV
        double e_u = first_ener*pow(ener_spacing,eu_dex); // last energy point in keV
        double h = (log(e_u)-log(e_l)) / 3; // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) / e_l);  // first trapezoid
        flux += h * (McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) + McPHACC4(em1_dex,cos_theta, T, M, R, mexmcc) * 3 + McPHACC4(em2_dex,cos_theta, T, M, R, mexmcc) * 3 + McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc)) * 3 / 8; // Simpson's 3/8, exact             
        flux += (E2 - e_u) / 2 * (McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc) / e_u + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2); // second trapezoid
    }
    if (n_steps == 5){ // five energy points within bandwidth: (4.1.3, 4.1.13) Simpson's "4,2" + two trapezoids  
        //cout << "5 steps" << endl;
        int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
        int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
        int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
        int em3_dex = e1_dex+4; // index of the fourth energy point within bandwidth        
        int eu_dex = e2_dex;    // index of the last energy point within bandwidth
        double e_l = first_ener*pow(ener_spacing,el_dex); // first energy point in keV
        double e_u = first_ener*pow(ener_spacing,eu_dex); // last energy point in keV
        double h = (log(e_u)-log(e_l)) / 4; // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) / e_l);  // first trapezoid
        flux += h * (McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) + McPHACC4(em1_dex,cos_theta, T, M, R, mexmcc) * 4 + McPHACC4(em2_dex,cos_theta, T, M, R, mexmcc) * 2 + McPHACC4(em3_dex,cos_theta, T, M, R, mexmcc) * 4 + McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc)) / 3; // Simpson's "4,2", exact             
        flux += (E2 - e_u) / 2 * (McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc) / e_u + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2); // second trapezoid
    }
    if (n_steps >= 6){ // six or more energy points within bandwidth: (4.1.3, 4.1.14) Simpson's cubic + two trapezoids
        //cout << "6 steps" << endl;
        int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
        int el1_dex = e1_dex+2; // index of the second energy point within bandwidth
        int el2_dex = e1_dex+3; // index of the third energy point within bandwidth        
        int eu_dex = e2_dex;    // index of the last energy point within bandwidth
        int eu1_dex = e2_dex-1; // index of the second last energy point within bandwidth
        int eu2_dex = e2_dex-2; // index of the third last energy point within bandwidth
        double e_l = first_ener*pow(ener_spacing,el_dex); // first energy point in keV
        double e_u = first_ener*pow(ener_spacing,eu_dex); // last energy point in keV
        double h = (log(e_u)-log(e_l)) / (n_steps-1); // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) / e_l);  // first trapezoid
        // Simpson's cubic, exact
        flux += h * (McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) * 9 + McPHACC4(el1_dex,cos_theta, T, M, R, mexmcc) * 28 + McPHACC4(el2_dex,cos_theta, T, M, R, mexmcc) * 23) / 24; // first three coefficients
        for (int m = 1; m <= n_steps-6; m++) flux += h * (McPHACC4(el_dex+m+3,cos_theta, T, M, R, mexmcc)); // middle coefficients
        flux += h * (McPHACC4(eu2_dex,cos_theta, T, M, R, mexmcc) * 23 + McPHACC4(eu1_dex,cos_theta, T, M, R, mexmcc) * 28 + McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc) * 9) / 24; // last three coefficients
        // end Simpson's cubic
        flux += (E2 - e_u) / 2 * (McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc) / e_u + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2); // second trapezoid                
    }

    flux = flux/Units::H_PLANCK;
    return flux;
}



// This is version makes use of Cole's version of McPhac
double AtmosEBandFlux3new( unsigned int model, double cos_theta, double T, double lgrav, double E1, double E2, class LightCurve mexmcc){

    int e1_dex(0);              // largest energy index that's smaller than E1
    int e2_dex(0);              // largest energy index that's smaller than E2
    int n_steps(0);             // number of energy points within band
    double flux(0.0);           // total integrated flux

 


    double ener_spacing = pow(10.0,0.0338);
    double first_ener = 0.04969469;
    double ener_index = log10(E1/T / first_ener) / log10(ener_spacing);
    e1_dex = (int) ener_index;
    ener_index = log10(E2/T / first_ener) / log10(ener_spacing);
    e2_dex = (int) ener_index;

    n_steps = 1* (e2_dex - e1_dex);

    // n_steps = 4;


    if (n_steps == 0){ // zero energy points within bandwidth: (4.1.3) one trapzoid
      //cout << "0 steps" << endl;
      flux = (E2 - E1) / 2.0 * (McPHACC3new(E1,cos_theta, T, lgrav, mexmcc) / E1 + McPHACC3new(E2,cos_theta, T, lgrav, mexmcc) / E2);
    }
    if (n_steps == 1){ // one energy points within bandwidth: (4.1.3) two trapzoids
      //cout << "1 step" << endl;
        int e_dex = e1_dex+1; // index of the energy point
        double e_m = first_ener*pow(ener_spacing,e_dex)*T; // energy point in keV
	
	e_m = 0.5*(E2+E1);

	double counts_m = McPHACC3new(e_m,cos_theta, T,lgrav, mexmcc) / e_m;

        flux = (e_m - E1) / 2 * (McPHACC3new(E1,cos_theta, T, lgrav, mexmcc) / E1 
				 + counts_m);  // first trapezoid
        flux += (E2 - e_m) / 2 * (counts_m + McPHACC3new(E2,cos_theta, T,lgrav, mexmcc) / E2); // second trapezoid
    }
    if (n_steps >= 2){ // two energy points within bandwidth: (4.1.3) three trapzoids
      //cout << "n steps= " << n_steps << endl;
     
       double e_step = (E2-E1)/(n_steps+1.0);
       double e_l, e_u;
       double counts_l, counts_u;

       flux = e_step * 0.5 * McPHACC3new(E1,cos_theta, T, lgrav, mexmcc) / E1;
       //cout << "flux = " << flux << endl;

       for (int i(1); i<=n_steps; i++){

	 e_l = E1 + i*e_step;
	 counts_l = McPHACC3new(e_l,cos_theta, T, lgrav, mexmcc) / e_l;
	 flux += e_step * (counts_l);
	 //cout << "flux = " << flux << endl;

       }

       flux += e_step * 0.5 *  McPHACC3new(E2,cos_theta, T,lgrav, mexmcc) / E2;
       // cout << "flux = " << flux << endl;

    }
   
    flux = flux/Units::H_PLANCK;
    return flux;
} // new version of energy integration





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
