/***************************************************************************************/
/*                                     Chi.cpp

    This holds the computational methods/functions called in Spot.cpp and probably
    elsewhere.
    
    And now it computes chi^2 too!
    
    PGxx means equation xx in Poutanen & Gierlinski 2003, MNRAS 343
    MLCBxx means equation xx in Morsink, Leahy, Cadeau & Braga 2007, ApJ 663 
*/
/***************************************************************************************/

// Changes
// 2016-08-05 - SMM: Function Bend creates the look-up table of bending angles
// re-upload of 2016-09-12 commitment


#include "matpack.h"
#include <exception>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unistd.h>
#include "Chi.h"
#include "OblDeflectionTOA.h"
#include "OblModelBase.h"
#include "PolyOblModelNHQS.h"
#include "Exception.h"
#include "Units.h"
#include "Struct.h"
#include "time.h"
#include "nrutil.h"
#include "interp.h"
#include <stdio.h>
using namespace std;

/**************************************************************************************/
/* ChiSquare:                                                                         */
/*           computes the chi^2 fit of data vs a simulation                           */
/*																					  */
/* pass: obsdata = light curve fluxes from observational data                         */
/*       ts = time shift or phase shift (normalized)                                  */
/**************************************************************************************/
double ChiSquare ( class DataStruct* obsdata, class LightCurve* curve) {
        
    /***************************************/
    /* VARIABLE DECLARATIONS FOR ChiSquare */
    /***************************************/
    
    unsigned int numbins,  // Number of phase (or time) bins the light curve is cut into
                 numbands;   

    int k,      // Array index variable
      n;      // Array index variable
    
    double ts,                              // time shift, so the phase of the simulation matches the phase of the data
         chisquare(0.0),                  // Computed chi^2
         tempflux[NCURVES][MAX_NUMBINS],  // Temporary array to store the flux
         min_location;                    // 
    
    cout << "starting chi squared calculation" << endl;
    numbins = obsdata->numbins;
    numbands = curve->numbands;
    ts = curve->para.ts;

    for ( unsigned int z(1); z<=1 ; z++ ) { // for different epochs

        while ( ts < 0.0 ) {
            ts += 1.0;
        }
        while ( ts > 1.0 ) {
            ts -= 1.0;
        }
        
        min_location = ts * numbins; // real version of the bin with the num
        int new_b = min_location;
        double new_shift = (min_location-new_b)/(numbins*1.0);
       
        // Rebinning the data and store shifted data back in Flux        
        for ( unsigned int i(0); i < numbins; i++ ) {
            k = i - new_b; //May changed
            if (k > static_cast<int>(numbins)-1) k -= numbins;
            if (k < 0) k += numbins;
            
	    for (unsigned int p(0);p<numbands;p++){
	      tempflux[p][i] = curve->f[p][k]; // putting things from curve->f's k bin into tempflux's i bin           
	    }
        }
        //cout << "integer shift complete" << endl;

        for ( unsigned int i(0); i < numbins; i++ ) {
            n = i - 1;
            if ( n < 0 ) n += numbins;
            // n = i+1;
            if ( n > static_cast<int>(numbins) - 1 ) n -= numbins;
            
	    for (unsigned int p(0); p < numbands; p++){
	      curve->f[p][i] = tempflux[p][i] + (tempflux[p][n]-tempflux[p][i]) * new_shift * numbins;         
	    }
        }
        //cout << "fractional shift complete complete" << endl;
        
        // Compute chisquare for shifted data
        
	chisquare = 0.0;
	for ( unsigned int j(0); j<numbands; j++){
	  obsdata->chi[j] = 0.0;
	  for ( unsigned int i(0); i < numbins; i++ ) {
	    obsdata->chi[j] += pow( (obsdata->f[j][i] - curve->f[j][i])/obsdata->err[j][i], 2);
            //chisquare += pow( (obsdata->f[j][i] - curve->f[j][i])/obsdata->err[j][i], 2);
	  }
	  std::cout << " chi^2[" << j << "] = " << obsdata->chi[j] ;
	  chisquare += obsdata->chi[j]; 
	}   

	std::cout << std::endl;
	std::cout << "Total Chi^2 = " << chisquare << std::endl;

    } // End epoch loop
    

    return chisquare;
    
}

class LightCurve SpotShape( int pieces, int p, int numtheta, double theta_1, double rho, class LightCurve* incurve, class OblModelBase* model){
  // If the spot is in one piece, then p=0
  // If the spot is in two pieces, then:
  // p=0 --> crescent shape
  // p=1 --> polar cap

  class LightCurve curve;
  curve = *incurve;

  double r_eq = model->R_at_costheta(0.0); // Equatorial radius

  if (p==1){ // polar cap

    double deltatheta = (rho-theta_1)/numtheta; //symmetric over pole
    for(int k(0); k < numtheta; k++){
      curve.para.dtheta[k] = deltatheta;
      curve.para.theta_k[k] = (k+0.5)*curve.para.dtheta[k];
      curve.para.phi_k[k] = Units::PI;
    }
  }
  if (p==0){

    if (curve.flags.spotshape!=1){
      double deltatheta(2.0*rho/numtheta);
      double deltamu(0.0);
      if (pieces==2){
	  deltatheta = 2.0*theta_1/numtheta; // crescent 
	  deltamu = 2.0*sin(rho)*sin(theta_1)/numtheta; 
      }
      for (int k(0); k < numtheta; k++){
	curve.para.dtheta[k] = deltatheta;
	double thetak = theta_1 - rho + (k+0.5)*deltatheta;  
	double muk(0.0);
	if (pieces==2){	  
	  muk = cos(rho-theta_1) - (k+0.5)*deltamu;
	  thetak = acos(muk);
	  curve.para.dtheta[k] = deltamu / sin(thetak);
	  //thetak = rho - theta_1 + (k+0.5)*curve.para.dtheta[k]; // midpoint of interval
	  //thetak = rho - theta_1 + (k+1.0)*deltatheta;	  // lower edge of interval
	  //thetak = rho - theta_1 + (k)*deltatheta;	  // top edge of interval
	}   
	curve.para.theta_k[k] = thetak;
		    
	  double cos_phi_edge = (cos(rho) - cos(theta_1)*cos(thetak))/(sin(theta_1)*sin(thetak));	
	  if (  cos_phi_edge > 1.0 || cos_phi_edge < -1.0 ) cos_phi_edge = 1.0;
	  if ( fabs( sin(theta_1) * sin(thetak) ) > 0.0) { // checking for a divide by 0
	    curve.para.phi_k[k] = acos( cos_phi_edge );   
	    // value of phi (a.k.a. azimuth projected onto equatorial plane) at the edge of the circular spot at some latitude thetak
	    // std::cout << "k="<< k << " theta_k="<< thetak << " phi_edge=" << curve.para.phi_k[k] << std::endl;

	}
      }
    } // end spotshape=0 case

  if (curve.flags.spotshape==1){

    if (p==0){ // case of just one spot that doesn't go over the pole; OR crescent shaped part of spot

   
	double dzeta(Units::PI/(numtheta));
	double zeta_0(0.0);
	double zeta_k(0.0);
     

      double north(1.0*rho); // Length of the segment joining the spot centre to the North-most part of the spot
      // Need to integrate to find north on an oblate star
      int k;
      double dtheta = fabs(rho)/(numtheta);

      north = 0.0;
      for (k=1; k<=numtheta; k++){ // Integrate from centre of spot to northernmost edge

	double thetak = theta_1 + dtheta*k; 
	double costhetak = cos(thetak);
	north += dtheta * model->R_at_costheta(costhetak)/r_eq;	
      }
      //std::cout << "North length = " << north << std::endl;


      //for (k=31; k < 32; k++){
      for (k=0; k < numtheta; k++){
	zeta_k = (0.5+k)*dzeta + zeta_0;

	// Trapezoidal Rule
	double a = 0.0;          // lower bound of integration
	double b = north;          // upper bound of integration
	double current_rho(0.0);      // current value of x, at which we are evaluating the integrand
	unsigned int current_n(0);  // current step
	unsigned int n_steps(10); // total number of steps
	double h = (b - a) / n_steps;     // step amount for numerical integration; the size of each step
	
	double length(0.0);           // the integrated length
	double oldlength(0.0);
	double rho_star(0.0);

	// begin trapezoidal rule
	current_rho = a + h * current_n;
	length = SpotIntegrand(current_rho,zeta_k, &curve, model);
	//length = 0.0;

	/*	std::cout << "rho=" << current_rho
		<< " length=" << length*h/2.0 << std::endl;*/

	for ( current_n = 1; current_n <= n_steps-1; current_n++ ) {
		current_rho = a + h * current_n;
		oldlength = length;
		length += 2.0 * SpotIntegrand(current_rho,zeta_k,&curve,model);
		/*	std::cout << "rho=" << current_rho
			<< " length=" << length*h/2.0 << std::endl;*/
		if (length >= north*2.0/h){
		  /* std::cout << " north = " << north
			    << " curr_rho = " << current_rho
			    << " length = " << length * h/2.0
			    << " old_rho = " << current_rho - h
			    << " oldlength = " << oldlength * h/2.0
			    << std::endl;*/
		  rho_star = current_rho - h + (north*2.0/h - oldlength)*h/(length-oldlength);
		  //std::cout << " interpol rho = " << rho_star << endl;
		}
	}
	
	current_rho = a + h * current_n;
	length += SpotIntegrand(current_rho,zeta_k,&curve,model);
	/*std::cout << "rho=" << current_rho
	  << " length=" << length*h/2.0 << std::endl;*/
	rho_star = current_rho;

		if (length >= north*2.0/h){
		  /*std::cout << " north = " << north
			    << " curr_rho = " << current_rho
			    << " length = " << length * h/2.0
			    << " old_rho = " << current_rho - h
			    << " oldlength = " << oldlength * h/2.0
			    << std::endl;*/
		  rho_star = current_rho - h + (north*2.0/h - oldlength)*h/(length-oldlength);
		  //std::cout << " interpol rho = " << rho_star << endl;
		}

	length *= h/2.0;	
	// end trapezoidal rule; numerical integration complete!

	/*	std::cout << "k = " << k 
		  << " zeta = " << zeta_k 
		  << " rho = " << rho_star
		  << " length = " << north
		  << std::endl;*/

	double theta_0(curve.para.theta_c);

	double costheta(cos(theta_0)*cos(rho_star) + sin(theta_0)*sin(rho_star)*cos(zeta_k));

	/*	std::cout << "k=" << k
		  << " cos(theta) = " << costheta
		  << " R(costheta)/R_eq = " << model->R_at_costheta(costheta)/model->R_at_costheta(0.0)
		  << std::endl;*/

	curve.para.theta_k[k] = acos(costheta);

	double denom = (sin(theta_0)*cos(rho_star) - cos(theta_0)*sin(rho_star)*cos(zeta_k));
	double tanphi = sin(rho_star)*sin(zeta_k)/(denom);

	if (denom>0.0)
	  curve.para.phi_k[k] = atan(tanphi);
	else
	  curve.para.phi_k[k] = Units::PI + atan(tanphi);

	/*	std::cout << "k = " << k
		  << " theta_k = " << curve.para.theta_k[k]
		  << " phi_k = " << curve.para.phi_k[k]
		  << std::endl;*/


      } // end for-k-loop

      if (pieces == 1){ // spot does not cover the pole
	k=0;
	curve.para.dtheta[k] = curve.para.theta_k[0] - (curve.para.theta_c - rho);
	for (k=1; k<numtheta-1; k++){
	  curve.para.dtheta[k] = 0.5 * (curve.para.theta_k[k+1] - curve.para.theta_k[k-1]);
	}
	k=numtheta-1;
	curve.para.dtheta[k] =  (curve.para.theta_c + rho) - curve.para.theta_k[k];
      }
      if (pieces == 2 && p==0){ //crescent -shaped bit

	k=0;
	curve.para.dtheta[k] = curve.para.theta_k[0] - (-curve.para.theta_c + rho);
	for (k=1; k<numtheta-1; k++){
	  curve.para.dtheta[k] = 0.5 * (curve.para.theta_k[k+1] - curve.para.theta_k[k-1]);
	}
	k=numtheta-1;
	curve.para.dtheta[k] =  (curve.para.theta_c + rho) - curve.para.theta_k[k];
      }


    }
  }// end spotshape=1 case
  } //p=0 case

  /*  std::cout << "pieces = " << pieces 
	    << " p = " << p
	    << " curve.para.radius = " << curve.para.req
	    << std::endl;*/

  for ( int k(0); k<numtheta; k++){
    double rtheta(curve.para.req); // value of r at present value of theta
    //double theta_test = curve.para.theta_k[k] + 0.5*curve.para.dtheta[k];
    double theta_test = curve.para.theta_k[k];
    rtheta = model->R_at_costheta( cos(theta_test));

    // What is radius versus rtheta

    double mass_over_r = curve.para.mass_over_r * curve.para.req/rtheta;
    double red = 1.0/sqrt(1.0-2.0*mass_over_r);
    double omega = curve.para.omega;
    //omega=0.0;
    double speed = omega*rtheta*sin(theta_test)*red; // MLCB34
    curve.para.gamma_k[k] = sqrt(1.0/(1.0-pow(speed,2))); 

    /*   std::cout << "k=" << k
	      << " curve.para.theta = " << theta_test
	      << " rtheta = " << rtheta 
	      << " red = " << red
	      << " speed = " << speed
	      << " gamma = " << curve.para.gamma_k[k]
	      << std::endl;
    */
  }


  return curve;

}

double SpotIntegrand( double rho, double zeta, class LightCurve* incurve, class OblModelBase* model){
  class LightCurve curve;
  curve = *incurve;

  double theta_0(curve.para.theta_c);
  double costheta(cos(theta_0)*cos(rho) + sin(theta_0)*sin(rho)*cos(zeta));
  double sintheta( sqrt(1.0 - pow(costheta,2)) );
  double r_eq(model->R_at_costheta(0.0));
  double rtheta(model->R_at_costheta(costheta));// value of r at present value of theta
 

  double mass_over_r = curve.para.mass_over_r * r_eq/rtheta; // change for oblate star
  double red = 1.0/sqrt(1.0-2.0*mass_over_r);
  double omega = curve.para.omega;
  //omega=0.0;
  double speed = omega*rtheta*sintheta*red; // MLCB34
  double gammasq = 1.0/(1.0-pow(speed,2)); 
 
  double integrand = rtheta/r_eq * sqrt( 1.0 + gammasq * pow( omega*rtheta*red * sin(theta_0)*sin(zeta),2));

  return integrand;
}




/**************************************************************************************/
/* ComputeAngles:                                                                     */
/*              computes all angles necessary to create the x-ray light curve         */
/*																					  */
/* pass: incurve = contains needed values like mass, radius, inclination, theta_0     */
/*       defltoa =                                                                    */
/**************************************************************************************/
class LightCurve ComputeAngles ( class LightCurve* incurve,
				                 class OblDeflectionTOA* defltoa) {

	/*******************************************/
	/* VARIABLE DECLARATIONS FOR ComputeAngles */
	/*******************************************/
	
    class LightCurve curve;
    curve = *incurve;

    double theta_0,                   // Emission angle (latitude) of the spot, in radians          
           incl,                      // inclination angle of the observer, in radians
           mass,                      // Mass of the star, in M_sun
           radius,                    // Radius of the star (at whatever little bit we're evaluating at)
      req,                            // Radius at the equator
      mass_over_r,
           omega,                     // Spin frequency of the neutron star, in Hz
           cosgamma,                  // Cos of the angle between the radial vector and the surface normal vector
           shift_t,                   // Time off-set from data
           b_guess(0.0),              // Impact parameter; starting off with reasonable guess then refining it
           mu(1.0),                   // = cos(theta_0), unitless
           speed(0.0),                // Velocity of the spot, defined in MLCB34
           alpha(0.0),                // Zenith angle, in radians
           sinalpha(0.0),             // Sin of zenith angle, defined in MLCB19
           cosalpha(1.0),             // Cos of zenith angle, used in MLCB30
           b(0.0),                    // Photon's impact parameter
      b_R, // b/radius
           toa_val(0.0),              // Time of arrival, MLCB38
           phi_0,                     // Azimuthal location of the spot, in radians
           dS,                        // Surface area of the spot, defined in MLCB2; computed in Spot.cpp
           distance;                  // Distance from Earth to NS, inputted in meters

    double red;
           
    unsigned int numbins(MAX_NUMBINS);// Number of phase bins the light curve is split into; same as number of flux data points
    numbins = curve.numbins;


    bool ingoing(false);
    bool infile_is_set(false);

    std::vector< double > phi_em(numbins, 0.0);   // Azimuth as measured from the star's emission frame; one for each phase bin
  

    // These are calculated in the second loop.
    //std::vector< double > dcosalpha_dcospsi(numbins, 0.0);    // Used in MLCB30
    std::vector< double > cosdelta(numbins, 0.0);             // 
    std::vector< double > cosxi(numbins, 0.0);                // Used in Doppler boost factor, defined in MLCB35

    // vectors for 4-point interpolation
    std::vector< double > psi_k(4, 0.0);       // Bending angle, sub k?
    std::vector< double > b_k(4, 0.0);         // Impact parameter, sub k?
    std::vector< double > dp_k(4, 0.0);      
    std::vector< double > t_k(4, 0.0);          

    /************************************************************************************/
    /* SETTING THINGS UP - keep in mind that these variables are in dimensionless units */
    /************************************************************************************/
    
    dS = curve.para.dS;
    theta_0 = curve.para.theta;
    incl = curve.para.incl;
    phi_0 = curve.para.phi_0;
    mass = curve.para.mass;
    radius = curve.para.radius; //present location on the star
    req =     curve.para.req;  // equatorial radius of the star
    mass_over_r = curve.para.mass_over_r;
    red = 1.0/sqrt(1.0-2.0*mass_over_r);
    omega = curve.para.omega;
    cosgamma = curve.para.cosgamma;  // for spherical, cosgamma = 0
    distance = curve.para.distance;
    shift_t = curve.para.ts;
    infile_is_set = curve.flags.infile_is_set;
    //speed = omega*radius*sin(theta_0) / sqrt( 1.0 - 2.0*mass_over_r ); // MLCB34
    speed = omega*radius*sin(theta_0)*red; // MLCB34
    //std::cout << "Speed = " << speed << std::endl;
    mu = cos(theta_0);

    double singamma(0.0);
    singamma = sqrt( 1.0 - pow( cosgamma, 2.0 ));

    //std::cout << "ComputeAngles:" << std::endl;
    //std::cout << "ComputeAngles: b_R_max = " << curve.defl.b_psi[3*NN] << curve.defl.b_R_max << std::endl;

    if (mu < 0.0){
      /* std::cout << "ComputeAngles: Southern Hemisphere!"
		<< " mu = " << mu 
		<< std::endl;*/
      singamma *= -1.0;
    }

    //initial assumptions
    curve.eclipse = false;
    curve.ingoing = false;
    curve.problem = false;
	
    /********************************************************/
    /* COMPUTE EMISSION TIME, PHASE BINS, AND LIGHT BENDING */
    /********************************************************/

    for ( unsigned int i(0); i < numbins; i++ ) { // opening For-Loop-1
        // SMM: Added an offset of phi_0
        // SMM: Time is normalized to the spin period so 0 < t_e < 1 
        // curve.t[i] = i/(1.0*numbins) + shift_t; (This is computed in Spot.cpp)
        phi_em.at(i) = phi_0 + (2.0*Units::PI) * curve.t[i]; // phi_em is the same thing as "phi" in PG; changes with t
	curve.cospsi[i] = cos(incl)*cos(theta_0) + sin(incl)*sin(theta_0)*cos(phi_em.at(i));
	curve.psi[i] = acos(curve.cospsi[i]);
       
    } // closing For-Loop-1

    int j(0);
    double b_R_min(0), psimin(0);

    // Check to see if this is an oblate star
    // If it is oblate compute the values of b_min, psi_max_in
    if (curve.flags.ns_model != 3){
      b_R_min = defltoa->b_R_min_ingoing(radius, cos(theta_0));
      psimin = defltoa->psi_ingoing(b_R_min,curve.defl.b_R_max, curve.defl.psi_max, radius,&curve.problem);
    }
    

    for ( unsigned int i(0); i < numbins; i++ ) { // opening For-Loop-2	

        int sign(0);
        double b_R_val(0.0);
        bool result(false);
        double b1(0.0), b2(0.0), psi1(0.0), psi2(0.0);
	double dc1(0.0), dc2(0.0), t1(0.0), t2(0.0);
        double xb(0.0);
        int k(0);
        //j=0;
        /**************************************************************************/
	/* TEST FOR VISIBILITY FOR EACH VALUE OF b, THE PHOTON'S IMPACT PARAMETER */
	/**************************************************************************/

	//	std::cout << "i = " << i 
	//	  << " psi = " << curve.psi[i] << std::endl; 

        if ( curve.psi[i] < curve.defl.psi_max ) {
            if (curve.psi[i] > curve.defl.psi_b[j] )
	            while ( curve.psi[i] > curve.defl.psi_b[j] ) {
	                j++;      
		    }
            else {
	            while ( curve.psi[i] < curve.defl.psi_b[j] ) {
	              j--;
	            }
	           j++;
            }
      
            b1 = curve.defl.b_psi[j-1];
            b2 = curve.defl.b_psi[j];
            psi1 = curve.defl.psi_b[j-1];
            psi2 = curve.defl.psi_b[j];
	    dc1 = curve.defl.dcosa_dcosp_b[j-1];
	    dc2 = curve.defl.dcosa_dcosp_b[j];
	    t1 = curve.defl.toa_b[j-1];
	    t2 = curve.defl.toa_b[j];
	   
            k = j - 2;      
            if ( j == 1 ) k = 0;
            if ( j == 3 * NN ) k = 3 * NN - 3;

	   
            for ( j = 0; j < 4; j++ ) {
	            b_k.at(j) = curve.defl.b_psi[k+j];
	            psi_k.at(j) = curve.defl.psi_b[k+j];
	            dp_k.at(j) = curve.defl.dcosa_dcosp_b[k+j];
	            t_k.at(j) = curve.defl.toa_b[k+j];
            }
  
            // 4-pt interpolation to find the correct value of b given psi.
	    xb = curve.psi[i];
            b_guess = (xb-psi_k.at(1))*(xb-psi_k.at(2))*(xb-psi_k.at(3))*b_k.at(0)/
	                  ((psi_k.at(0)-psi_k.at(1))*(psi_k.at(0)-psi_k.at(2))*(psi_k.at(0)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(2))*(xb-psi_k.at(3))*b_k.at(1)/
	                  ((psi_k.at(1)-psi_k.at(0))*(psi_k.at(1)-psi_k.at(2))*(psi_k.at(1)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(1))*(xb-psi_k.at(3))*b_k.at(2)/
	                  ((psi_k.at(2)-psi_k.at(0))*(psi_k.at(2)-psi_k.at(1))*(psi_k.at(2)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(1))*(xb-psi_k.at(2))*b_k.at(3)/
                	  ((psi_k.at(3)-psi_k.at(0))*(psi_k.at(3)-psi_k.at(1))*(psi_k.at(3)-psi_k.at(2)));
			 

	    curve.dcosalpha_dcospsi[i] = (xb-psi_k.at(1))*(xb-psi_k.at(2))*(xb-psi_k.at(3))*dp_k.at(0)/
	                  ((psi_k.at(0)-psi_k.at(1))*(psi_k.at(0)-psi_k.at(2))*(psi_k.at(0)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(2))*(xb-psi_k.at(3))*dp_k.at(1)/
	                  ((psi_k.at(1)-psi_k.at(0))*(psi_k.at(1)-psi_k.at(2))*(psi_k.at(1)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(1))*(xb-psi_k.at(3))*dp_k.at(2)/
	                  ((psi_k.at(2)-psi_k.at(0))*(psi_k.at(2)-psi_k.at(1))*(psi_k.at(2)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(1))*(xb-psi_k.at(2))*dp_k.at(3)/
                	  ((psi_k.at(3)-psi_k.at(0))*(psi_k.at(3)-psi_k.at(1))*(psi_k.at(3)-psi_k.at(2)));

            toa_val = (xb-psi_k.at(1))*(xb-psi_k.at(2))*(xb-psi_k.at(3))*t_k.at(0)/
	                  ((psi_k.at(0)-psi_k.at(1))*(psi_k.at(0)-psi_k.at(2))*(psi_k.at(0)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(2))*(xb-psi_k.at(3))*t_k.at(1)/
	                  ((psi_k.at(1)-psi_k.at(0))*(psi_k.at(1)-psi_k.at(2))*(psi_k.at(1)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(1))*(xb-psi_k.at(3))*t_k.at(2)/
	                  ((psi_k.at(2)-psi_k.at(0))*(psi_k.at(2)-psi_k.at(1))*(psi_k.at(2)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(1))*(xb-psi_k.at(2))*t_k.at(3)/
                	  ((psi_k.at(3)-psi_k.at(0))*(psi_k.at(3)-psi_k.at(1))*(psi_k.at(3)-psi_k.at(2)));

	    // b_guess = (b2-b1)/(psi2-psi1) * (curve.psi[i] - psi2) + b2;
	    // curve.dcosalpha_dcospsi[i] = (dc2-dc1)/(psi2-psi1) * (curve.psi[i] - psi2) + dc2;
	    // toa_val = (t2-t1)/(psi2-psi1) * (curve.psi[i] - psi2) + t2;

        } // ending psi.at(i) < curve.defl.psi_max

	b_guess *= radius;
	b2 *= radius;

        /***********************************************/
	/* FINDING IF A SOLUTION EXISTS, SETTING FLAGS */
	/***********************************************/

	//Change so that bval is really b/radius;

        result = defltoa->b_from_psi( curve.psi[i], radius, mu, b_R_val, sign, curve.defl.b_max, 
				      curve.defl.psi_max, b_R_min*radius,psimin, b_guess,
				      &curve.problem );

        if ( result == false && i == 0) { 
            curve.visible[i] = false;
            curve.t_o[i] = 0.0;
            curve.dOmega_s[i] = 0.0;
            curve.eclipse = true;
	     
        }
        else if ( result == false ) { 
            curve.visible[i] = false;
            curve.t_o[i] = curve.t[i] + curve.t_o[i-1] - curve.t[i-1];
            curve.dOmega_s[i] = 0.0;
            curve.eclipse = true;
        }
        else { // there is a solution
	  //b = bval;
	    b_R = b_R_val;
            if ( sign < 0 ) { // if the photon is initially ingoing (only a problem in oblate models)
	            ingoing = true;
	            curve.ingoing = true;
	            //std::cout << "ingoing!"<< std::endl;
            }
            else if ( sign > 0 ) {
	            ingoing = false;
            }
            else {
	            throw( Exception("Chi.cpp: sign not returned as + or - with success.") ); // used to say "ObFluxApp.cpp"
            }
            
	    double b_maximum = radius/sqrt(1.0 - 2.0*mass_over_r);
	    if ( (fabs(b-b_maximum) < 1e-7) && (b > 0.0) && (b > b_maximum) ) { 
	      // this corrects for b being ever so slightly over bmax, which yields all kinds of errors in OblDeflectionTOA
	      std::cout << "Setting b = b_max." << std::endl;
	      b = b_maximum - DBL_EPSILON;
	    }
            curve.b[i] = b/radius;
            
            /*******************************************************/
	    /* IF A SOLUTION EXISTS, CALCULATING ANGLES:  sinalpha */
	    /*                    						  cosalpha */
	    /*                    						  alpha    */
	    /*                   						  cosdelta */
	    /*                							  cosbeta  */
	    /*******************************************************/


            sinalpha =  b_R * sqrt( 1.0 - 2.0 * mass_over_r );  // PG4, reconfigured
            cosalpha = sqrt(fabs( 1.0 - sinalpha * sinalpha )); 
	    //   alpha    = asin( sinalpha );
            

	    if (sinalpha > 1)
	      alpha = 0.5*Units::PI;
	    else
	      alpha    = asin( sinalpha );


	    if ( sign < 0 ) { // alpha is greater than pi/2
	      //std::cout << "i=" << i << " sign < 0 " << std::endl;
	      alpha = Units::PI - alpha;
	      cosalpha *= -1;
	    }



            cosdelta.at(i) =  (cos(incl) - cos(theta_0)*curve.cospsi[i]) / (sin(theta_0)*sin(curve.psi[i])) ;

            if ( theta_0 == 0.0 )  // cosdelta needs to be redone if theta_0 = 0
            	cosdelta.at(i) = sqrt( 1 - pow( sin(incl) * sin(phi_em.at(i)) / sin(curve.psi[i]) ,2) );
	    // law of sines from MLCB Fig 1 and in MLCB17 and just above
     
            if ( (cos(theta_0) < 0) && (cos(incl) < 0 ) ) {
	            cosdelta.at(i) *= -1.0;
            }

            if ( sin(curve.psi[i]) != 0.0 ) {
	            curve.cosbeta[i] = cosalpha * cosgamma + sinalpha * singamma * cosdelta.at(i);
            }
            else {
	            curve.cosbeta[i] = cosalpha * cosgamma;
            }

	    if ( fabs(cosalpha) < 0.01 ) {
	            curve.cosbeta[i] = (cosalpha*cosgamma +  singamma * cosdelta.at(i));
		    //std::cout << "small cosbeta=" << curve.cosbeta[i] << std::endl;
		    }


            if ( curve.cosbeta[i] < 0.0 || curve.cosbeta[i] > 1.0 ) { // cosbeta > 1.0 added by Abbie, Feb 22 2013
	      /* std::cerr << "i = " << i
	        	          << ", Not visible at phi_em = " << 180.0 * phi_em.at(i) / Units::PI 
	                   	  << ", cos(beta) = " << curve.cosbeta[i] 
		                  << ", cos(alpha) = " << cosalpha
		                  << ", cos(gamma) = " << cosgamma
		                  << ", cos(delta) = " << cosdelta.at(i)
		                  << " (visibility condition)." << std::endl << std::endl;*/
	            curve.visible[i] = false;
            }
            else {
	            curve.visible[i] = true;
            }
			
			/********************************************************/
			/* IF THE PHOTON IS VISIBLE, COMPUTE: dpsi_db           */     
			/*                                    toa               */
			/*									  cosxi             */
			/*									  eta               */
			/*									  t_o               */
			/*									  psi               */
			/*									  R_dpsi_db         */
			/*									  dcosalpha_dcospsi */
			/*									  dOmega_s          */
			/********************************************************/
			
            if ( curve.visible[i] ) { // visible 
            	if (alpha == 0.0 && curve.psi[i] == 0.0 & phi_em.at(i) == 0.0) 
            		cosxi.at(i) = 0.0; // to avoid NAN errors from dividing by 0; appears when incl = theta at i=0
	            else 
	            	cosxi.at(i) = - sinalpha * sin(incl) * sin(phi_em.at(i)) / sin(curve.psi[i]);  // PG11
	            curve.eta[i] = sqrt( 1.0 - speed*speed ) / (1.0 - speed*cosxi.at(i) ); // Doppler boost factor, MLCB33
	            
		    

	            if ((1.0 - speed*cosxi.at(i)) == 0.0) 
	            	std::cout << "dividing by zero" << std::endl;
	            if((std::isnan(speed) || speed == 0.0) && theta_0 != 0.0 )
	            	std::cout << "speed = " << speed << " at i = " << i << std::endl;

	            if ( ingoing ) {
		      b = b_R * radius;
		      toa_val = defltoa->toa_ingoing( b, radius, mu, &curve.problem );
	            }
                
		    toa_val += (req/radius - 1.0);
		    toa_val += -  2.0 * mass_over_r * log( (radius/req) * (1.0 - 2.0 * mass_over_r)/(1.0 - 2.0 * mass_over_r*radius/req)) ;

		    //Change toa_val so it is dimensionless. Then multiply by radius to give correct dimensions!!!!		    
		    // Correct for emission from location different from equator.

	            curve.t_o[i] = curve.t[i] + (omega * toa_val * radius) / (2.0 * Units::PI);	           
         

	            curve.dOmega_s[i] = (dS / (distance * distance)) 
	                               * (1.0 / (1.0 - 2.0 * mass_over_r)) 
	                               * curve.cosbeta[i] 
	                               * curve.dcosalpha_dcospsi[i];  // PG8

	        

            } // end visible
            
            /**************************************************************/
			/* IF THE PHOTON IS NOT VISIBLE, POLITELY CRASH OR GO TO ZERO */
			/**************************************************************/
      
            else { // not visible; we think that it shouldn't matter if it's not visible at i=0
	      curve.t_o[i] = curve.t[i] + (curve.t_o[i-1] - curve.t[i-1]) ; // t_o is not defined properly, so we'll set it to emission time
	      //curve.t_o[i] = curve.t[i] ;
	      curve.dOmega_s[i] = 0.0;    // don't see the spot, so dOmega = 0
	      curve.cosbeta[i] = 0.0;     // doesn't matter, doesn't enter into calculation
	      curve.eta[i] = 1.0;	        // doesn't matter, doesn't enter into calculation
	      // std::cout << "toa = " << (curve.t_o[i-1] - curve.t[i-1]) << std::endl;
	            //}
            } // end not visible
        } // end "there is a solution"
    }  // closing For-Loop-2

    //  std::cout << "End of ComputeAngles" << std::endl; 

    return curve;

} // End ComputeAngles

/**************************************************************************************/
/* ReadBend:                                                                          */
/*         Reads in the table of alpha, b, psi, dpsi_db for different M/R values      */
/* pass: incurve = contains needed values like mass, radius, inclination, theta_0     */
/*       defltoa =                                                                    */
/**************************************************************************************/
class LightCurve ReadBend ( class LightCurve* incurve,
				                 char *bend_file) {
  std::ifstream in;     // Bending Angle File

  class LightCurve curve;
  curve = *incurve;

  char line[265]; // line of the data file being read in
  double get_alpha, get_bR, get_psi, get_dcos, get_toa, get_mr;
  unsigned int get_i, numangles, num_mr;
  
  in.open(bend_file);
  // Read in value of NN
  in.getline(line,265);      
  sscanf( line, "%u", &numangles);
  std::cout << "numangles = " << numangles << std::endl;
  if ( numangles != NN)
    std::cerr << "Something is really wrong here: numangles != NN "<< std::endl;

  // Read in value of num_mr
  in.getline(line,265);      
  sscanf( line, "%u", &num_mr);
  curve.defl.num_mr = num_mr;
  std::cout << "num_mr = " << num_mr << std::endl;
  if ( num_mr > MR)
    std::cerr << "Problem: num_mr > MR... check struct.h " << std::endl;
 
  curve.defl.mr = dvector(0,num_mr);
  curve.defl.psi = dmatrix(0,num_mr,0,3*numangles+1);
  curve.defl.b = dmatrix(0,num_mr,0,3*numangles+1);
  curve.defl.dcosa = dmatrix(0,num_mr,0,3*numangles+1);
  curve.defl.toa = dmatrix(0,num_mr,0,3*numangles+1);
  
     // Loop through the M/R values
      for (unsigned int j(0);j<=num_mr; j++){
	//in.getline(line,265);      
	//sscanf( line, "% M/R= %lf", &get_mr );
	//std::cout << "line = " << line << std::endl;
	//curve.defl.mr[j] = get_mr;
	//in.getline(line,265);  

	for (unsigned int i(0);i<=3*numangles; i++){

	  in.getline(line,265);  
	  sscanf( line, "%lf %lf %lf %lf %lf %lf %u", &get_alpha, &get_bR, &get_psi, &get_dcos, &get_toa, &get_mr, &get_i );

	  if (i==0){
	    curve.defl.mr[j] = get_mr;
	    //std::cout << "M/R = " << curve.defl.mr[j] << std::endl;
	  }

	  curve.defl.b[j][i] = get_bR;
	  curve.defl.psi[j][i] = get_psi;
	  curve.defl.dcosa[j][i] = get_dcos;
	  curve.defl.toa[j][i] = get_toa;
	  
	} // end-for-i loop

      } // end-for-j loop     

      in.close();
    
      //std::cout << "ReadBend: test psi = " << curve.defl.psi[10][10] << std::endl;

      return curve;

}




/**************************************************************************************/
/* Bend:                                                                          */
/*        Interpolates to find bending angles for specific M/R value              */
/* pass: incurve = contains needed values like mass, radius, inclination, theta_0     */
/*       defltoa =                                                                    */
/**************************************************************************************/
class LightCurve Bend ( class LightCurve* incurve,
				                 class OblDeflectionTOA* defltoa) {

	/*******************************************/
	/* VARIABLE DECLARATIONS FOR Bend          */
	/*******************************************/
	
    class LightCurve curve;
    curve = *incurve;

    std::ofstream bend;      // output stream; printing information to the output file
    std::ifstream angles;   // input stream; read in the bending angles from a file

    int printflag(0); // Set to 1 if you want to print out stuff to a file.
    int readflag(0);
    int computeflag(0); // Set to 1 if you want to compute the bending angles.
    int interpflag(1);

    double 
      eps,
      mass,                      // Mass of the star, in M_sun
      radius,                    // Radius of the star (at whatever little bit we're evaluating at)
      mass_over_r,
      alpha(0.0),                // Zenith angle, in radians
      sinalpha(0.0),             // Sin of zenith angle, defined in MLCB19
      cosalpha(1.0),             // Cos of zenith angle, used in MLCB30
      b_R,
      psi(0.0),
      toa_val(0.0),              // Time of arrival, MLCB38
      dpsi_db_val(0.0),          // Derivative of MLCB20 with respect to b
      dcosa_dcosp;
          
  
    /* std::cout << "Entering Bend: M/R = " 
	      << curve.para.mass_over_r
	      << std::endl;

	      std::cout << "Bend: test psi = " << curve.defl.psi[10][10] << std::endl;*/
   
    /************************************************************************************/
    /* SETTING THINGS UP - keep in mind that these variables are in dimensionless units */
    /************************************************************************************/
    
    mass = curve.para.mass;
    radius = curve.para.radius;
    mass_over_r = curve.para.mass_over_r;

    if (computeflag){ // Compute Lookup Table


   /**********************************************************/
    /* Compute maximum deflection for purely outgoing photons */
    /**********************************************************/
	
    double  b_mid;  // the value of b, the impact parameter, at 90% of b_max
    curve.defl.b_max =  defltoa->bmax_outgoing(radius); // telling us the largest value of b
    curve.defl.b_R_max = curve.defl.b_max/radius;
    curve.defl.psi_max = defltoa->psi_max_outgoing_u(curve.defl.b_max/radius,&curve.problem); // telling us the largest value of psi

    /********************************************************************/
    /* COMPUTE b VS psi LOOKUP TABLE, GOOD FOR THE SPECIFIED M/R AND mu */
    /********************************************************************/
	
    b_mid = curve.defl.b_R_max * 0.9; 
    // since we want to split up the table between coarse and fine spacing. 
    // 0 - 90% is large spacing, 90% - 100% is small spacing. b_mid is this value at 90%.
    curve.defl.b_psi[0] = 0.0; // definitions of the actual look-up table for b when psi = 0
    curve.defl.psi_b[0] = 0.0; // definitions of the actual look-up table for psi when b = 0
    dpsi_db_val = defltoa->dpsi_db_outgoing_u( 0.0, &curve.problem );
    curve.defl.dcosa_dcosp_b[0] = fabs( (1.0 - 2.0 * mass_over_r) / (dpsi_db_val));
    curve.defl.toa_b[0] = toa_val = defltoa->toa_outgoing_u( 0.0, &curve.problem );

    // computation for b < b_mid
    for ( unsigned int i(1); i < NN+1; i++ ) { /* compute table of b vs psi points */
        curve.defl.b_psi[i] = b_mid * i / (NN * 1.0);
        curve.defl.psi_b[i] = defltoa->psi_outgoing_u(curve.defl.b_psi[i], curve.defl.b_R_max, curve.defl.psi_max, &curve.problem); 
	dpsi_db_val = defltoa->dpsi_db_outgoing_u( curve.defl.b_psi[i], &curve.problem );
	curve.defl.toa_b[i] = defltoa->toa_outgoing_u(curve.defl.b_psi[i], &curve.problem );

	sinalpha = curve.defl.b_psi[i] * sqrt(1.0 - 2.0*mass_over_r);
	cosalpha = sqrt(1.0 - pow(sinalpha,2));

	curve.defl.dcosa_dcosp_b[i] = fabs( sinalpha/cosalpha * sqrt(1.0 - 2.0*mass_over_r) / (sin(curve.defl.psi_b[i]) *dpsi_db_val));
    }

    // Compute slope of d(cosalpha)/d(psi) near cos(alpha)=0
      eps = 1e-2;
      double x1, x2, y1, y2, yslope, yy;
  
      x1 = Units::PI/2.0 - eps; 
      x2 = x1+eps*0.1;

      b_R = sin(x1) / sqrt( 1.0 - 2.0 * mass_over_r );
      psi = defltoa->psi_outgoing_u( b_R, curve.defl.b_R_max, curve.defl.psi_max, &curve.problem);
      dpsi_db_val = defltoa->dpsi_db_outgoing_u( b_R, &curve.problem );
      y1 = fabs( sin(x1)/cos(x1) * sqrt(1.0 - 2.0*mass_over_r) / (sin(fabs(psi)) * dpsi_db_val));

      b_R = sin(x2) / sqrt( 1.0 - 2.0 * mass_over_r );
      psi = defltoa->psi_outgoing_u( b_R, curve.defl.b_R_max, curve.defl.psi_max, &curve.problem);
      dpsi_db_val = defltoa->dpsi_db_outgoing_u( b_R,&curve.problem );
      y2 = fabs( sin(x2)/cos(x2) * sqrt(1.0 - 2.0*mass_over_r) / (sin(fabs(psi)) * dpsi_db_val));

      yslope = (y2-y1)/(x2-x1);


    // computation for b > b_mid
    for ( unsigned int i(NN+1); i < 3*NN; i++ ) { /* compute table of b vs psi points */
        curve.defl.b_psi[i] = b_mid + (curve.defl.b_R_max - b_mid) / 2.0 * (i - NN) / (NN * 1.0); 
        curve.defl.psi_b[i] = defltoa->psi_outgoing_u(curve.defl.b_psi[i], curve.defl.b_R_max, curve.defl.psi_max, &curve.problem); 
	sinalpha = curve.defl.b_psi[i] * sqrt(1.0 - 2.0*mass_over_r);
	curve.defl.toa_b[i] = defltoa->toa_outgoing_u(curve.defl.b_psi[i], &curve.problem );


	if (cosalpha >= cos(x1)){
	  dpsi_db_val = defltoa->dpsi_db_outgoing_u( curve.defl.b_psi[i], &curve.problem );
	  cosalpha = sqrt(1.0 - pow(sinalpha,2));
	  curve.defl.dcosa_dcosp_b[i] = fabs( sinalpha/cosalpha * sqrt(1.0 - 2.0*mass_over_r) / (sin(curve.defl.psi_b[i]) * dpsi_db_val));
	}
	else{
	  alpha = asin(sinalpha);
	  yy = y2 + yslope*(alpha - x2);
	  curve.defl.dcosa_dcosp_b[i] = yy;
	}
    }
    
    curve.defl.b_psi[3*NN] = curve.defl.b_R_max;   // maximums
    curve.defl.psi_b[3*NN] = curve.defl.psi_max;
    curve.defl.toa_b[3*NN] = defltoa->toa_outgoing_u(curve.defl.b_max/radius, &curve.problem );

    alpha = Units::PI/2.0;
    yy = y2 + yslope*(alpha - x2);
    curve.defl.dcosa_dcosp_b[3*NN]  = yy;
    // Finished computing lookup table
    }    

    if (readflag){ // Read in values instead of computing them 

      // Allocate Memory -- Look up table for specific M/R value
      curve.defl.psi_b = dvector(0,3*NN+1);
      curve.defl.b_psi = dvector(0,3*NN+1);
      curve.defl.dcosa_dcosp_b = dvector(0,3*NN+1);
      curve.defl.toa_b = dvector(0,3*NN+1);

      angles.open("bend.txt");
      char line[265]; // line of the data file being read in
      //unsigned int numLines(0);
      double get_alpha, get_bR, get_psi, get_dcos, get_toa;

      // std::cout << "Reading in a file is not yet implemented" << std::endl;
 
      for (unsigned int i(1);i<=5;i++){
	angles.getline(line,265);
	std::cout << "line=" << line << std::endl;
	}

      for (unsigned int i(0);i<=3*NN;i++){
	angles.getline(line,265);
	sscanf( line, "%lf %lf %lf %lf %lf", &get_alpha, &get_bR, &get_psi, &get_dcos, &get_toa );
	//std::cout << "line=" << line << std::endl;
	//std::cout << "radius = " << radius << " b/R=" << get_bR << std::endl;
	curve.defl.b_psi[i] = get_bR ;
	curve.defl.psi_b[i] = get_psi;
	curve.defl.toa_b[i] = get_toa; // toa is now dimensionless
	curve.defl.dcosa_dcosp_b[i] = get_dcos;
      }

      curve.defl.b_R_max =  curve.defl.b_psi[3*NN];
      curve.defl.b_max = curve.defl.b_R_max * radius;
      curve.defl.psi_max =  curve.defl.psi_b[3*NN];
    } // End readflag==1

    //printflag = 0;
    if (printflag){

      //initial assumptions
      curve.problem = false;
	
      // Open output file and print out header

      bend.open("bend.txt", std::ios_base::trunc);
      bend.precision(10);
      bend << "# Mass = " <<  Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN << " Msun "    << std::endl;
      bend << "# R_sp = " << Units::nounits_to_cgs(radius, Units::LENGTH )*1.0e-5 << " km " << std::endl;
      bend << "# M/R = " << curve.para.mass_over_r << std::endl;
      bend << "# R/M = " << 1.0/curve.para.mass_over_r << std::endl;
      bend << "#alpha #b/R #psi #dcosalpha/dcospsi #toa*C/R" << std::endl;

      for (unsigned int i(0); i <= 3*NN; i++ ) { 

	//b = curve.defl.b_psi[i] * radius;
	sinalpha = curve.defl.b_psi[i] * sqrt( 1.0 - 2.0 * mass_over_r );
	alpha = asin(sinalpha);

	psi = curve.defl.psi_b[i];
	dcosa_dcosp = curve.defl.dcosa_dcosp_b[i];
	toa_val = curve.defl.toa_b[i];
				
	    bend 
	      << alpha << " " 
	      << curve.defl.b_psi[i] << " " 
	      << psi << " " 
	      << dcosa_dcosp << " " 
	      << toa_val << " "
	      << std::endl;	       		
       	   
      }
    } // End printflag==1

    if (interpflag){

    // Interpolate to find bending angles for a particular m/r value

      double mr_lo, mr_hi, mr_step;
      unsigned int num_mr(curve.defl.num_mr);
      unsigned int numangles(NN);

      mr_lo = curve.defl.mr[0];
      mr_hi = curve.defl.mr[num_mr];
      mr_step = (mr_hi-mr_lo)/num_mr;

      int jlo, jhi;
    
      if (mass_over_r < mr_hi){
	jlo = (mass_over_r - mr_lo)/mr_step;
	jhi = jlo+1;
      }
      else{
	jlo = num_mr-1;
	jhi = num_mr;
      }

      /*     std::cout << "mass_over_r = " << mass_over_r <<
	" jlo=" << jlo << " m/r[jlo]=" << curve.defl.mr[jlo] << std::endl;
      std::cout << "mass_over_r = " << mass_over_r <<
	" jhi=" << jhi << " m/r[jhi]=" << curve.defl.mr[jhi] << std::endl;
      */

  

    curve.defl.psi_b = interplin( curve.defl.mr, curve.defl.psi, numangles, mass_over_r, &jlo);
    curve.defl.b_psi =  interplin( curve.defl.mr, curve.defl.b, numangles, mass_over_r, &jlo);
    curve.defl.dcosa_dcosp_b =  interplin( curve.defl.mr, curve.defl.dcosa, numangles, mass_over_r, &jlo);
    curve.defl.toa_b =  interplin( curve.defl.mr, curve.defl.toa, numangles, mass_over_r, &jlo);
    
    curve.defl.psi_max = curve.defl.psi_b[3*numangles];
    curve.defl.b_R_max = curve.defl.b_psi[3*numangles];
    curve.defl.b_max =  curve.defl.b_R_max * curve.para.radius;

    // std::cout << "b/R_max = " << curve.defl.b_R_max << std::endl;
    /*
      bend.open("test.txt", std::ios_base::trunc);
      bend.precision(10);
      bend << "# Mass = " <<  Units::nounits_to_cgs(mass, Units::MASS)/Units::MSUN << " Msun "    << std::endl;
      bend << "# R_sp = " << Units::nounits_to_cgs(radius, Units::LENGTH )*1.0e-5 << " km " << std::endl;
      bend << "# M/R = " << curve.para.mass_over_r << std::endl;
      bend << "# R/M = " << 1.0/curve.para.mass_over_r << std::endl;
      bend << "#alpha #b/R #psi #dcosalpha/dcospsi #toa*C/R" << std::endl;

      for (unsigned int i(0); i <= 3*NN; i++ ) { 

	//b = curve.defl.b_psi[i] * radius;
	sinalpha = curve.defl.b_psi[i] * sqrt( 1.0 - 2.0 * mass_over_r );
	alpha = asin(sinalpha);
	if (i==3*NN)
	  alpha = Units::PI/2.0;

	psi = curve.defl.psi_b[i];
	dcosa_dcosp = curve.defl.dcosa_dcosp_b[i];
	toa_val = curve.defl.toa_b[i];
				
	    bend 
	      << alpha << " " 
	      << curve.defl.b_psi[i] << " " 
	      << psi << " " 
	      << dcosa_dcosp << " " 
	      << toa_val << " "
	      << std::endl;	       		
       	   
      }

    */

    } // End interpflag==1


    //std::cout << "Bend: TEST psi[10] = " << curve.defl.psi_b[10] << std::endl;

    return curve;

} // End Bend



/**************************************************************************************/
/* ShiftCurve:                                                                      */
/*              Shifts curve by azimuthal angle phi                                 */
/*																					  */
/* pass: angles = all the angles necessary to compute the flux correctly;             */
/*                computed in the routine/method/function above [radians or unitless] */
/**************************************************************************************/
class LightCurve ShiftCurve( class LightCurve* angles, double phishift) {
		
    class LightCurve curve;

    std::ifstream input;
    std::ofstream output;
    std::ofstream ttt;

    unsigned int numbins(MAX_NUMBINS);  // Time bins of light curve (usually 128)
    unsigned int numbands(NCURVES);  // Number of Energy Bands
        
    std::vector< double > totflux(MAX_NUMBINS, 0.0); // integrated flux
    std::vector< bool > nullcurve(NCURVES,false); // true means that the curve is zero everywhere

    double timeshift;

    /*********************/
    /* SETTING THINGS UP */
    /*********************/

    curve = (*angles);
    numbins = curve.numbins;
    numbands = curve.numbands;

    timeshift = phishift/(2.0*Units::PI);

    //std::cout << "TimeShift = " << timeshift <<std::endl;

    // Shift light curve by angle phishift<(2pi)/numbins
  
    int k(0), j(0); // index placeholders; approximately, k is i+1 and j is i-1

        
      /********************************/
      /* LOOP THROUGH THE LIGHTCURVES */
      /********************************/
    	
    for (unsigned int p(0); p < numbands; p++) {
  
	if (!nullcurve[p]) {
	

	/*********************/
	/* MORE DECLARATIONS */
	/*********************/
    		
	  std::vector< double > newflux(MAX_NUMBINS, 0.0);                  // rebinned flux (to account for photon travel time)
	  double minimum(100000.0);                           // true (continuous) maximum and minimum flux values
	 
	  int ec1(0), ec2(0),j1,j2;
	  double te1, te2, tt1, tt2, slope1, slope2;

  	                   
	  /*******************************************/
	  /* FOR ECLIPSING LIGHT CURVES, MINIMUM = 0 */
	  /*******************************************/
    		
	  if ( curve.eclipse ) {
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
	    tt1 = curve.t[j1];
	    if (ec1<2)
	      tt1 +=-1.0;
	    tt2 = curve.t[j2];
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
	    tt1=curve.t[j1];
	    tt2=curve.t[j2];
	    if (ec2==numbins-1)
	      tt2+=1.0;

	    slope2 = (curve.f[p][j1]-curve.f[p][j2])/(tt1-tt2);
	    te2 = tt1 - curve.f[p][j1]/slope2;
	    //std::cout << "Eclipse ends at time = " << te2 << std::endl;
	  } // ending "yes eclipsed
            
		           	   
	  for ( unsigned int i(0); i < numbins; i++ ) {  // for-i-loop, looping through the phase bins
            
		double t1, t2;
		// Find point to the left of the "ith" point; call it j	       
		if ( i==0){		     
		  j = numbins-1;
		  t1 = curve.t[j]-1.0;
		  k = 0;
		  t2 = curve.t[k];		    
		} // end if i==0		 
		else { //Default
		  j = i-1;
		  t1 = curve.t[j]; // time to the left of the time we're interested in
		  k = i;
		  t2 = curve.t[k]; // time to the right of the point we're interested in
		}
		newflux.at(i) = curve.f[p][k] - (curve.f[p][k]-curve.f[p][j])/(t2-t1) * (timeshift); // linear interpolation!

		if (curve.eclipse){
		  if (curve.t[i] < te1 && curve.t[i] > te1 - 1.0/numbins){
		    //std::cout << "time just before eclipse at te1=" << te1 <<"! i=" << i << " t = " << curve.t[i] << std::endl;
		    newflux.at(i) = -slope1 * (timeshift) + curve.f[p][k];
		  }
		  if (curve.t[i] < te1 + 1.0/numbins && curve.t[i] > te1)
		    newflux.at(i) = 0.0;

		  if (curve.t[i] > te2 && curve.t[i] < te2 + 1.0/numbins){
		    // std::cout << "eclipse ends at te2=" << te2 <<"! i=" << i << " t = " << curve.t[i] << std::endl;
		    newflux.at(i) = -slope2 * (timeshift) + curve.f[p][k];
		  }

		  if (curve.t[i] > te1 && curve.t[i] < te2)
		    newflux.at(i) = 0.0;
		}

		
	  
	
		/* newflux vs t corresponds to the new shifted light curve.
		   It corresponds to the old curve is curve.f */
	    
		if ( newflux.at(i) < 0.0 )
		  newflux.at(i) = 0.0;
	  
	  } // closes the for(i) loop, going through the curves.
    	
	 
	  // if (p==0)
	  //ttt.open("shift.txt", std::ios_base::trunc);

	  for ( unsigned int i(0); i < numbins; i++ ) {
	    /*if ( p==0 ) 
	      ttt << curve.t[i] << " " << curve.f[p][i] << " " << curve.t[i] << " " << newflux.at(i) 
	    	  << " " << i
	    	  << std::endl;
	    */
            curve.f[p][i] = newflux.at(i);
	  }
	 

	} // end NOT NULL
    }// end for-p-loop
   
    
    return curve;

} // end ShiftCurve


/**************************************************************************************/
/* Normalize1:                                                                         */
/*           normalizes the fluxes in each energy band to 1 by dividing by 
	     the average flux of each light curve													  */
/*																					  */
/* pass: Flux = the flux for each part of the light curve [units don't matter]        */
/*       numbins = number of phase bins the light curve is divided into [Unitless]    */
/**************************************************************************************/
class LightCurve Normalize1( double Flux[NCURVES][MAX_NUMBINS], unsigned int numbins ) {

  /***************************************/
  /* VARIABLE DECLARATIONS FOR Normalize */
  /***************************************/
	
    class LightCurve newcurve;  // new curve structure to hold the normalized lightcurve flux values
    unsigned int i, p;          // loop variables
    double norm[NCURVES];       // normalization factor, one for each curve

    for ( p = 0; p < NCURVES; p++ ) {
        norm[p] = 0.0;
        for ( i = 0; i < numbins; i++ ) {
            norm[p] += Flux[p][i]; 
            newcurve.f[p][i] = Flux[p][i];
        }
        if ( norm[p] != 0.0) norm[p] /= (numbins*1.0); // makes norm the average value for each curve
    }


    for ( p = 0; p < NCURVES; p++ ) {
        for ( i = 0; i < numbins; i++ ) {
            if ( norm[p] != 0.0)  newcurve.f[p][i] /= norm[p]; 
            else newcurve.f[p][i] = 1.0;
        }
    }
    
    return newcurve;
} // end Normalize1

/**************************************************************************************/
/* Normalize2:                                                                         */
/*           normalizes the fluxes in low energy band to 1 by dividing by 
	     the average flux of the low energy light curve		       	  */
/*									       	  */
/* pass: Flux = the flux for each part of the light curve [units don't matter]        */
/*       numbins = number of phase bins the light curve is divided into [Unitless]    */
/**************************************************************************************/
class LightCurve Normalize2( double Flux[NCURVES][MAX_NUMBINS], unsigned int numbins ) {

  /***************************************/
  /* VARIABLE DECLARATIONS FOR Normalize */
  /***************************************/
	
    class LightCurve newcurve;  // new curve structure to hold the normalized lightcurve flux values
    unsigned int i, p;          // loop variables
    double norm[NCURVES];       // normalization factor, one for each curve

    for ( p = 0; p < NCURVES; p++ ) {
        norm[p] = 0.0;
        for ( i = 0; i < numbins; i++ ) {
            norm[p] += Flux[p][i]; 
            newcurve.f[p][i] = Flux[p][i];
        }
        if ( norm[p] != 0.0) norm[p] /= (numbins*1.0); // makes norm the average value for each curve
    }


    for ( p = 0; p < NCURVES; p++ ) {
        for ( i = 0; i < numbins; i++ ) {
            if ( norm[2] != 0.0)  newcurve.f[p][i] /= norm[2]; 
            else newcurve.f[p][i] = 1.0;
        }
    }
    
    return newcurve;
} // end Normalize2



#define TINY 1.0e-10
#define NMAX 300
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define NDIM 5
#define MPTS 6








/**************************************************************************************/
/* LegP2:																			  */
/*       calculating a Legendre polynomial P2										  */
/*       called in the routine/method/function recalc								  */
/*																					  */
/* pass: costheta = cos of angle of the hotspot down from geom. north pole [unitless] */
/**************************************************************************************/
double LegP2( double costheta ) {
    double P2;
    P2 = 0.5 * ( 3.0 * pow( costheta, 2 ) - 1.0 );
    return P2;
} // end LegP2

/**************************************************************************************/
/* LegP4:																			  */
/*       calculating a Legendre polynomial P4										  */
/*       called in the routine/method/function recalc								  */
/*																					  */
/* pass: costheta = cos of angle of the hotspot down from geom. north pole [unitless] */
/**************************************************************************************/
double LegP4( double costheta ) {
    double P4;
    P4 = 0.125 * ( 35 * pow( costheta, 4 ) - 30 * pow( costheta, 2 ) + 3);
    return P4;
} // end LegP4


// end Chi.cpp
