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
    
    unsigned int numbins;  // Number of phase (or time) bins the light curve is cut into
    
    int k,      // Array index variable
    	n;      // Array index variable
    
    double ts,                              // time shift, so the phase of the simulation matches the phase of the data
    	   chisquare(0.0),                  // Computed chi^2
    	   tempflux[NCURVES][MAX_NUMBINS],  // Temporary array to store the flux
    	   min_location;                    // 
    
    cout << "starting chi squared calculation" << endl;
    numbins = obsdata->numbins;
    ts = curve->para.ts;
    cout << "pointers set" << endl;
    cout << "timeshift is " << ts << endl;
    for ( unsigned int z(1); z<=1 ; z++ ) { // for different epochs
        cout << obsdata->f[0][0] << " " << curve->f[0][0] << " " << obsdata->err[0][0] << endl;
        
        while ( ts < 0.0 ) {
            ts += 1.0;
        }
        while ( ts > 1.0 ) {
            ts -= 1.0;
        }
        
        min_location = ts * numbins; // real version of the bin with the num
        int new_b = min_location;
        double new_shift = (min_location-new_b)/(numbins*1.0);
        //cout << numbins << endl;
        cout << ts << " " << min_location << " " << new_shift << " " << new_b << endl;

        // Rebinning the data and store shifted data back in Flux
        
        for ( unsigned int i(0); i < numbins; i++ ) {
            k = i - new_b; //May changed
            if (k > static_cast<int>(numbins)-1) k -= numbins;
            if (k < 0) k += numbins;
            
            unsigned int p = 0;
        	unsigned int q = 0;
            tempflux[p][i] = curve->f[q][k]; // putting things from curve->f's k bin into tempflux's i bin
            p = 1;
        	q = 1;
        	tempflux[p][i] = curve->f[q][k]; // putting things from curve->f's k bin into tempflux's i bin
        }
        cout << "integer shift complete" << endl;
        for ( unsigned int i(0); i < numbins; i++ ) {
            n = i - 1;
            if ( n < 0 ) n += numbins;
            // n = i+1;
            if ( n > static_cast<int>(numbins) - 1 ) n -= numbins;
            
            unsigned int p = 0;
        	unsigned int q = 0;
            curve->f[q][i] = tempflux[p][i] + (tempflux[p][n]-tempflux[p][i]) * new_shift * numbins;
            p = 1;
        	q = 1;
        	curve->f[q][i] = tempflux[p][i] + (tempflux[p][n]-tempflux[p][i]) * new_shift * numbins;
        }
        cout << "fractional shift complete complete" << endl;
        
        // Compute chisquare for shifted data
        
        // energy band 1
        unsigned int p = 0;
        unsigned int q = 0;
        //cout << obsdata->f[0][i] << " " << curve->f[0][i] << " " << obsdata->err[0][i] << endl;
        for ( unsigned int i(0); i < numbins; i++ ) {
            chisquare += pow( (obsdata->f[p][i] - curve->f[q][i])/obsdata->err[p][i], 2);
            // using different 'p' and 'q' because we're not comparing the exact same columns of obsdata.f and Flux
        }
        cout << "chisquare for band 1 is " << chisquare << endl;
		// energy band 2
        p = 1;
        q = 1;
        for ( unsigned int i(0); i < numbins; i++ ) {
            chisquare += pow( (obsdata->f[p][i] - curve->f[q][i])/obsdata->err[p][i], 2);
            // using different 'p' and 'q' because we're not comparing the exact same columns of obsdata.f and Flux
        }

        cout << "chisquare for band 1+2 is " << chisquare << endl;

    }
    
    //chisquare += pow( (bbrat - 0.32)/(0.064) , 2);
    
    //std::cout << "M = " << Units::nounits_to_cgs(curve->para.mass, Units::MASS)/Units::MSUN << ", R = " << Units::nounits_to_cgs(curve->para.radius, Units::LENGTH )*1.0e-5 << ", i = " 
    //<< curve->para.incl * 180.0 / Units::PI <<", e = " << curve->para.theta * 180.0 / Units::PI <<", ts = " << ts << std::endl;

    //std::cout << "\nChisquare = " << chisquare << std::endl;

    return chisquare;
    
}

class LightCurve SpotShape( int pieces, int p, int numtheta, double theta_1, double rho, class LightCurve* incurve){

  class LightCurve curve;
  curve = *incurve;

  if (curve.flags.spotshape==0){

  double deltatheta(2.0*rho/numtheta);

  if (pieces==2)
    if (p==0)
      deltatheta = 2.0*theta_1/numtheta; // crescent or single piece
    else
      deltatheta = (rho-theta_1)/numtheta; //symmetric over pole

  for (int k(0); k < numtheta; k++){

    curve.para.dtheta[k] = deltatheta;

    double thetak = theta_1 - rho + (k+0.5)*deltatheta; 

    if (pieces==2){
      if (p==1){
	thetak = (k+0.5)*curve.para.dtheta[k];
	curve.para.phi_k[k] = Units::PI;
      }
      else {
	thetak = rho - theta_1 + (k+0.5)*deltatheta;
      }
    }
   
    curve.para.theta_k[k] = thetak;

    if (p==0){ //crescent-shaped piece, or the one circular piece if spot doesn't go over pole
	    
      double cos_phi_edge = (cos(rho) - cos(theta_1)*cos(thetak))/(sin(theta_1)*sin(thetak));
	   
      if (  cos_phi_edge > 1.0 || cos_phi_edge < -1.0 ) cos_phi_edge = 1.0;
      if ( fabs( sin(theta_1) * sin(thetak) ) > 0.0) { // checking for a divide by 0
	curve.para.phi_k[k] = acos( cos_phi_edge );   
	// value of phi (a.k.a. azimuth projected onto equatorial plane) at the edge of the circular spot at some latitude thetak
      }	   
    }



  }
  }
  return curve;

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
           toa_val(0.0),              // Time of arrival, MLCB38
           dpsi_db_val(0.0),          // Derivative of MLCB20 with respect to b
           phi_0,                     // Azimuthal location of the spot, in radians
           dS,                        // Surface area of the spot, defined in MLCB2; computed in Spot.cpp
           distance;                  // Distance from Earth to NS, inputted in meters

    double eps, epspsi,dcosa_dcosp, red;
           
    unsigned int numbins(MAX_NUMBINS);// Number of phase bins the light curve is split into; same as number of flux data points
    numbins = curve.numbins;


    bool ingoing(false);
    bool infile_is_set(false);

    std::vector< double > phi_em(numbins, 0.0);   // Azimuth as measured from the star's emission frame; one for each phase bin
    std::vector< double > psi(numbins, 0.0);      // Bending angle, as defined in MLCB20

    // These are calculated in the second loop.
    std::vector< double > dcosalpha_dcospsi(numbins, 0.0);    // Used in MLCB30
    std::vector< double > cosdelta(numbins, 0.0);             // 
    std::vector< double > cosxi(numbins, 0.0);                // Used in Doppler boost factor, defined in MLCB35

    // vectors for 4-point interpolation
    std::vector< double > psi_k(4, 0.0);       // Bending angle, sub k?
    std::vector< double > b_k(4, 0.0);         // Impact parameter, sub k?
    
    /************************************************************************************/
    /* SETTING THINGS UP - keep in mind that these variables are in dimensionless units */
    /************************************************************************************/
    
    dS = curve.para.dS;
    theta_0 = curve.para.theta;
    incl = curve.para.incl;
    phi_0 = curve.para.phi_0;
    mass = curve.para.mass;
    radius = curve.para.radius;
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

    //std::cout << "ComputeAngles: radius = " << radius << std::endl;


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
        psi.at(i) = acos(cos(incl)*cos(theta_0) + sin(incl)*sin(theta_0)*cos(phi_em.at(i))); // PG1; this theta_0 is the theta_0 from spot.cpp
       
    } // closing For-Loop-1

    int j(0);
    double bmin(0), psimin(0);

    // Check to see if this is an oblate star
    // If it is oblate compute the values of b_min, psi_max_in
    if (curve.flags.ns_model != 3){
      bmin = defltoa->bmin_ingoing(radius, cos(theta_0));
      psimin = defltoa->psi_ingoing(bmin,curve.defl.b_max, curve.defl.psi_max, cos(theta_0),&curve.problem);
       std::cout << "ComputeAngles: Oblate! b_min_ingoing = " << bmin 
		<< " psi_min_ingoing = " << psimin
		<< " b_max_outgoing = " << curve.defl.b_max
		<< " psi_max = " << curve.defl.psi_max
		<< std::endl;
    }
    

	
    for ( unsigned int i(0); i < numbins; i++ ) { // opening For-Loop-2
        int sign(0);
        double bval(0.0);
        bool result(false);
        double b1(0.0), b2(0.0), psi1(0.0), psi2(0.0);
        double xb(0.0);
        int k(0);
        //j=0;
        /**************************************************************************/
	/* TEST FOR VISIBILITY FOR EACH VALUE OF b, THE PHOTON'S IMPACT PARAMETER */
	/**************************************************************************/

        if ( psi.at(i) < curve.defl.psi_max ) {
            if ( psi.at(i) > curve.defl.psi_b[j] )
	            while ( psi.at(i) > curve.defl.psi_b[j] ) {
	                j++;      
                }
            else {
	            while ( psi.at(i) < curve.defl.psi_b[j] ) {
	              j--;
	            }
	           j++;
            }
      
            b1 = curve.defl.b_psi[j-1];
            b2 = curve.defl.b_psi[j];
            psi1 = curve.defl.psi_b[j-1];
            psi2 = curve.defl.psi_b[j];
            k = j - 2;
      
            if ( j == 1 ) k = 0;
            if ( j == 3 * NN ) k = 3 * NN - 3;

            for ( j = 0; j < 4; j++ ) {
	            b_k.at(j) = curve.defl.b_psi[k+j];
	            psi_k.at(j) = curve.defl.psi_b[k+j];
            }
  
            // 4-pt interpolation to find the correct value of b given psi.
            xb = psi.at(i);
            b_guess = (xb-psi_k.at(1))*(xb-psi_k.at(2))*(xb-psi_k.at(3))*b_k.at(0)/
	                  ((psi_k.at(0)-psi_k.at(1))*(psi_k.at(0)-psi_k.at(2))*(psi_k.at(0)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(2))*(xb-psi_k.at(3))*b_k.at(1)/
	                  ((psi_k.at(1)-psi_k.at(0))*(psi_k.at(1)-psi_k.at(2))*(psi_k.at(1)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(1))*(xb-psi_k.at(3))*b_k.at(2)/
	                  ((psi_k.at(2)-psi_k.at(0))*(psi_k.at(2)-psi_k.at(1))*(psi_k.at(2)-psi_k.at(3)))
	                  +(xb-psi_k.at(0))*(xb-psi_k.at(1))*(xb-psi_k.at(2))*b_k.at(3)/
                	  ((psi_k.at(3)-psi_k.at(0))*(psi_k.at(3)-psi_k.at(1))*(psi_k.at(3)-psi_k.at(2)));
        } // ending psi.at(i) < curve.defl.psi_max

        /***********************************************/
	/* FINDING IF A SOLUTION EXISTS, SETTING FLAGS */
	/***********************************************/

	//std::cout << "Compute Angles: i = " << i << " psi = " << psi.at(i) << std::endl;

        result = defltoa->b_from_psi( fabs(psi.at(i)), radius, mu, bval, sign, curve.defl.b_max, 
				      curve.defl.psi_max, bmin,psimin, b_guess, fabs(psi.at(i)), b2, fabs(psi.at(i))-psi2, 
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
            b = bval;
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


            sinalpha =  b * sqrt( 1.0 - 2.0 * mass_over_r ) / radius;  // PG4, reconfigured
            cosalpha = sqrt( 1.0 - sinalpha * sinalpha ); 
            alpha    = asin( sinalpha );
            

            if ( sign < 0 ) { // alpha is greater than pi/2
	            alpha = Units::PI - alpha;
            }

            cosdelta.at(i) =  (cos(incl) - cos(theta_0)*cos(psi.at(i))) / (sin(theta_0)*sin(psi.at(i))) ;
            if ( theta_0 == 0.0 )  // cosdelta needs to be redone if theta_0 = 0
            	cosdelta.at(i) = sqrt( 1 - pow( sin(incl) * sin(phi_em.at(i)) / sin(psi.at(i)) ,2) ); // law of sines from MLCB Fig 1 and in MLCB17 and just above
     
            if ( (cos(theta_0) < 0) && (cos(incl) < 0 ) ) {
	            cosdelta.at(i) *= -1.0;
            }
            if ( sin(psi.at(i)) != 0.0 ) {
	            curve.cosbeta[i] = cosalpha * cosgamma + sinalpha * sqrt( 1.0 - pow( cosgamma, 2.0 )) * cosdelta.at(i);
	            //if( std::isnan(curve.cosbeta[i]) ) std::cout << "cosdelta.at(i="<<i<<") = " << cosdelta.at(i) << std::endl;
            }
            else {
	            curve.cosbeta[i] = cosalpha * cosgamma;
            }

            if ( cosalpha < 0.01 ) {
	            curve.cosbeta[i] = (Units::PI/2.0 - alpha + sqrt(2) * sqrt(1.0-cosgamma) * cosdelta.at(i));
  			}
            if ( curve.cosbeta[i] < 0.0 || curve.cosbeta[i] > 1.0 ) { // cosbeta > 1.0 added by Abbie, Feb 22 2013
	            std::cerr << "i = " << i
	        	          << ", Not visible at phi_em = " << 180.0 * phi_em.at(i) / Units::PI 
	                   	  << ", cos(beta) = " << curve.cosbeta[i] 
		                  << ", cos(alpha) = " << cosalpha
		                  << ", cos(gamma) = " << cosgamma
		                  << ", cos(delta) = " << cosdelta.at(i)
		                  << " (visibility condition)." << std::endl << std::endl;
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
            	if (alpha == 0.0 && psi.at(i) == 0.0 & phi_em.at(i) == 0.0) 
            		cosxi.at(i) = 0.0; // to avoid NAN errors from dividing by 0; appears when incl = theta at i=0
	            else 
	            	cosxi.at(i) = - sinalpha * sin(incl) * sin(phi_em.at(i)) / sin(fabs(psi.at(i)));  // PG11
	            curve.eta[i] = sqrt( 1.0 - speed*speed ) / (1.0 - speed*cosxi.at(i) ); // Doppler boost factor, MLCB33
	            
	            if ((1.0 - speed*cosxi.at(i)) == 0.0) 
	            	std::cout << "dividing by zero" << std::endl;
	            if((std::isnan(speed) || speed == 0.0) && theta_0 != 0.0 )
	            	std::cout << "speed = " << speed << " at i = " << i << std::endl;
	            //if(std::isnan(cosxi.at(i)) || cosxi.at(i) == 0.0) 
	            //	std::cout << "cosxi(i="<<i<<") = " << cosxi.at(i) << std::endl;
	            

	            if ( ingoing ) {
	              //std::cout << "ComputeAngles:Ingoing b = " << b << std::endl;
		      dpsi_db_val = defltoa->dpsi_db_ingoing_u( b, radius, mu, &curve.problem );
		      toa_val = defltoa->toa_ingoing( b, radius, mu, &curve.problem );
	            }
                else {

		  if (b != curve.defl.b_max ){
	                dpsi_db_val = defltoa->dpsi_db_outgoing_u( b, radius, &curve.problem );
	                //if (i == 0) std::cout << "b = " << b <<", dpsi_db = " << dpsi_db_val << std::endl;

	                toa_val = defltoa->toa_outgoing_u( b, radius, &curve.problem );
	                //std::cout << "dpsi_db_val = " << dpsi_db_val << ", toa_val = " << toa_val << std::endl;
		  }
		  else{
		    toa_val = defltoa->toa_outgoing_u( b, radius, &curve.problem );
		    	eps = 1e-6;
			b = curve.defl.b_max * sqrt(1.0 - eps);
			epspsi = defltoa->psi_outgoing( b, radius, curve.defl.b_max, curve.defl.psi_max, &curve.problem);
			dpsi_db_val = defltoa->dpsi_db_outgoing_u( b, radius, &curve.problem );
			dcosa_dcosp = sqrt(1.0-2*mass_over_r) / (sqrt(eps) * sin(fabs(epspsi)) * radius * dpsi_db_val) * sqrt(1.0-eps);

		  }


		}

	            //std::cout << "Done computing TOA " << std::endl;
	            curve.t_o[i] = curve.t[i] + (omega * toa_val) / (2.0 * Units::PI);
	            curve.psi[i] = psi.at(i);
	            curve.R_dpsi_db[i] = dpsi_db_val * radius;

		    if (b != curve.defl.b_max){
		    if ( psi.at(i) == 0 && alpha == 0 ) 
		      curve.dcosalpha_dcospsi[i] = fabs( (1.0 - 2.0 * mass_over_r) / curve.R_dpsi_db[i]);
		    //if (psi.at(i) == 0 && alpha == 0 ) curve.dcosalpha_dcospsi[i] = 0.0;
	            else 
		      curve.dcosalpha_dcospsi[i] = fabs( sinalpha/cosalpha * sqrt(1.0 - 2.0*mass_over_r) / (sin(fabs(psi.at(i))) * curve.R_dpsi_db[i]) );
		    }
		    else
		      curve.dcosalpha_dcospsi[i] = dcosa_dcosp; 
		   

	            curve.dOmega_s[i] = (dS / (distance * distance)) 
	                               * (1.0 / (1.0 - 2.0 * mass_over_r)) 
	                               * curve.cosbeta[i] 
	                               * curve.dcosalpha_dcospsi[i];  // PG8
		    //std::cout << "t = " << curve.t[i] << " dOmega = " << curve.dOmega_s[i] << std::endl;

	            /***********************************/
				/* FLAGS IF A VALUE IS NAN OR ZERO */
				/***********************************/
	            if (std::isnan(dpsi_db_val) || dpsi_db_val == 0) std::cout << "dpsi_db_val = " << dpsi_db_val << "at i = " << i << std::endl;
				if (std::isnan(psi.at(i))) std::cout << "psi.at(i="<<i<<") = " << psi.at(i) << std::endl;
				if (std::isnan(curve.dOmega_s[i])) std::cout << "dOmega is NAN at i = " << i << std::endl;
				if (std::isnan(curve.cosbeta[i])) std::cout << "cosbeta is NAN at i = " << i << std::endl;
				if (std::isnan(curve.dcosalpha_dcospsi[i])) std::cout << "dcosalpha_dcospsi is NAN at i = " << i<< std::endl;
				if (std::isnan(sinalpha)) std::cout << "sinalpha is NAN at i = " << i << std::endl;
				if (std::isnan(cosalpha)) std::cout << "cosalpha is NAN at i = " << i << std::endl;
				if (std::isnan(psi.at(i))) std::cout << "psi is NAN at i = " << i << std::endl;
				if (std::isnan(mass)) std::cout << "mass is NAN at i = " << i << std::endl;
				if (std::isnan(radius)) std::cout << "radius is NAN at i = " << i << std::endl;
				

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
	      std::cout << "toa = " << (curve.t_o[i-1] - curve.t[i-1]) << std::endl;
	            //}
            } // end not visible
        } // end "there is a solution"
    }  // closing For-Loop-2

    return curve;

} // End ComputeAngles

/**************************************************************************************/
/* Bend:                                                                              */
/*              Creates a Table of alpha, b, psi, dpsi_db                             */
/*																					  */
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

    int printflag(0); // Set to 1 if you want to print out stuff to a file.

    double 
      eps,
      mass,                      // Mass of the star, in M_sun
      radius,                    // Radius of the star (at whatever little bit we're evaluating at)
      mass_over_r,
      alpha(0.0),                // Zenith angle, in radians
      sinalpha(0.0),             // Sin of zenith angle, defined in MLCB19
      cosalpha(1.0),             // Cos of zenith angle, used in MLCB30
      b(0.0),                    // Photon's impact parameter
      psi(0.0),
      toa_val(0.0),              // Time of arrival, MLCB38
      dpsi_db_val(0.0),          // Derivative of MLCB20 with respect to b
      dcosa_dcosp;
          
    unsigned int numbins(MAX_NUMBINS);// Number of phase bins the light curve is split into; same as number of flux data points
    numbins = 100;


    //std::cout << "Entering Bend" << std::endl;

    /************************************************************************************/
    /* SETTING THINGS UP - keep in mind that these variables are in dimensionless units */
    /************************************************************************************/
    
    mass = curve.para.mass;
    radius = curve.para.radius;
    mass_over_r = curve.para.mass_over_r;

   /**********************************************************/
    /* Compute maximum deflection for purely outgoing photons */
    /**********************************************************/
	
    double  b_mid;  // the value of b, the impact parameter, at 90% of b_max
    curve.defl.b_max =  defltoa->bmax_outgoing(radius); // telling us the largest value of b
    curve.defl.psi_max = defltoa->psi_max_outgoing_u(curve.defl.b_max,radius,&curve.problem); // telling us the largest value of psi

    /********************************************************************/
    /* COMPUTE b VS psi LOOKUP TABLE, GOOD FOR THE SPECIFIED M/R AND mu */
    /********************************************************************/
	
    b_mid = curve.defl.b_max * 0.9; 
    // since we want to split up the table between coarse and fine spacing. 
    // 0 - 90% is large spacing, 90% - 100% is small spacing. b_mid is this value at 90%.
    curve.defl.b_psi[0] = 0.0; // definitions of the actual look-up table for b when psi = 0
    curve.defl.psi_b[0] = 0.0; // definitions of the actual look-up table for psi when b = 0
    
    for ( unsigned int i(1); i < NN+1; i++ ) { /* compute table of b vs psi points */
        curve.defl.b_psi[i] = b_mid * i / (NN * 1.0);
        curve.defl.psi_b[i] = defltoa->psi_outgoing_u(curve.defl.b_psi[i], radius, curve.defl.b_max, curve.defl.psi_max, &curve.problem); 
	// calculates the integral
    }

    // For arcane reasons, the table is not evenly spaced.
    for ( unsigned int i(NN+1); i < 3*NN; i++ ) { /* compute table of b vs psi points */
        curve.defl.b_psi[i] = b_mid + (curve.defl.b_max - b_mid) / 2.0 * (i - NN) / (NN * 1.0); 
	// spacing for the part where the points are closer together
        curve.defl.psi_b[i] = defltoa->psi_outgoing_u(curve.defl.b_psi[i], radius, curve.defl.b_max, curve.defl.psi_max, &curve.problem); // referenced same as above
    }
    
    curve.defl.b_psi[3*NN] = curve.defl.b_max;   // maximums
    curve.defl.psi_b[3*NN] = curve.defl.psi_max;
    // Finished computing lookup table
    
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

      for (unsigned int i(0); i < numbins; i++ ) { 

	alpha = i/(numbins-1.0) * Units::PI/2.0; 
	sinalpha = sin(alpha);
	cosalpha = sqrt( 1.0 - sinalpha * sinalpha ); 

	b = sinalpha * radius / sqrt( 1.0 - 2.0 * mass_over_r );
	psi = defltoa->psi_outgoing_u( b, radius, curve.defl.b_max, curve.defl.psi_max, &curve.problem);
	dpsi_db_val = defltoa->dpsi_db_outgoing_u( b, radius, &curve.problem );
	toa_val = defltoa->toa_outgoing_u( b, radius, &curve.problem );

	if (psi==0 && alpha == 0)
	  dcosa_dcosp =  fabs( (1.0 - 2.0 * mass_over_r) / (radius * dpsi_db_val));
      	//dcosa_dcosp =  fabs( (1.0 - 2.0 * mass_over_r));
	else
	  dcosa_dcosp = fabs( sinalpha/cosalpha * sqrt(1.0 - 2.0*mass_over_r) / (sin(fabs(psi)) * radius*dpsi_db_val));

	if (cosalpha == 0){
	  eps = 5e-7;
	  b = curve.defl.b_max * sqrt(1.0 - eps);
	  psi = defltoa->psi_outgoing_u( b, radius, curve.defl.b_max, curve.defl.psi_max, &curve.problem);
	  dpsi_db_val = defltoa->dpsi_db_outgoing_u( b, radius, &curve.problem );
	  dcosa_dcosp = sqrt(1.0-2*mass_over_r) / (sqrt(eps) * sin(fabs(psi)) * radius * dpsi_db_val) * sqrt(1.0-eps);

	  b = curve.defl.b_max;
	  psi = defltoa->psi_outgoing_u( b, radius, curve.defl.b_max, curve.defl.psi_max, &curve.problem);

	}
		
	bend 
	  //<< i << " " 
	  << alpha << " " 
	  << b/radius << " " 
	  //	<< cosalpha << " "
	  << psi << " " 
	  //<< cos(psi) << " "
	  //<< 1.0/dcosa_dcosp << " " 
	  << dcosa_dcosp << " " 
	  //<< Units::nounits_to_cgs(toa_val,Units::TIME)  << " " 
	  << toa_val / radius 
	  << std::endl;
	       			   
      }
    }

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
