/***************************************************************************************/
/*                                     Atmo.cpp

    This holds the ComputeCurve routine used in Spot.cpp and atmosphere routines that are 
    used by ComputeCurve.

	Was split from Chi.cpp, the large files with every computational routines.

*/
/********************************************************************************NS*******/

#include "matpack.h"
#include <exception>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unistd.h>
#include "Atmo.h"
#include "McPhac.h"
#include "BlackBody.h"
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
    mass_over_r, mass_over_req,
    temperature,        // Temperature of the spot, in keV
    E_band_lower_1,     // Lower bound of energy band for flux integration, in keV
    E_band_upper_1,     // Upper bound of energy band for flux integration, in keV
    E_band_lower_2,     // Lower bound of energy band for flux integration, in keV
    E_band_upper_2,     // Upper bound of energy band for flux integration, in keV
    redshift,           // Gravitational redshift = 1 + z = (1-2M/R)^{-1/2}
    bolo,               // Bolometric flux; bolo = sigma T^4/pi
    gray(1.0);          // Graybody factor (when = 1, not effective)
  
  double E0, E1, E2, DeltaE;

  unsigned int numbins(MAX_NUMBINS);  // Time bins of light curve (usually 128)
  unsigned int numbands(NCURVES);  // Number of Energy Bands
        
  //  std::vector< double > totflux(MAX_NUMBINS, 0.0); // integrated flux




  /*********************/
  /* SETTING THINGS UP */
  /*********************/
    
  // One monochromatic energy, hardwired value, in keV
  //    E_mono = 1.0;
  std::cout << "ComputeCurve: Starting ComputeCurve" << std::endl;
  curve = (*angles);

  mass_over_r = curve.para.mass_over_r;
  mass_over_req = mass_over_r * curve.para.radius/curve.para.req;

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
  std::cout << "beaming model = " << curve.flags.beaming_model << std::endl;
   
  redshift = 1.0 / sqrt( 1 - 2.0 * mass_over_r);

  //std::cout << "redshift = " << redshift << std::endl;

  // double M = Units::nounits_to_cgs(curve.para.mass, Units::MASS);
  double R = Units::nounits_to_cgs(curve.para.req, Units::LENGTH);
  //   double delta = 1 / sqrt(1 - 2.0*curve.para.mass/curve.para.req);
  double delta = 1.0 / sqrt(1.0 - 2.0*mass_over_req);

  double cos_theta = cos(curve.para.theta);
  double obl_approx = 1.0
    + ((-0.791 + 0.776 * mass_over_req) * curve.para.omega_bar_sq  
       + (-1.315 * mass_over_req + 2.431 * pow(mass_over_req,2)) * pow(curve.para.omega_bar_sq,2.0)
       + (-1.172 * mass_over_req) * pow(curve.para.omega_bar_sq,3.0)) * (1.0-pow(cos_theta,2)) 
    + ((1.138 - 1.431 * mass_over_req) * curve.para.omega_bar_sq 
       + (0.653 * mass_over_req - 2.864 * pow(mass_over_req,2)) * pow(curve.para.omega_bar_sq,2.0)
       + (0.975 * mass_over_req) * pow(curve.para.omega_bar_sq,3.0))  *  pow(cos_theta,2)
    + (13.47 * mass_over_req - 27.13 * pow(mass_over_req,2)) *  pow(curve.para.omega_bar_sq,2.0) * cos_theta * (1.0 - cos_theta);
  double lgrav = log10(delta * mass_over_req * pow(Units::C,2)/R * obl_approx);
  //std::cout << "theta = " << curve.para.theta * 180.0/Units::PI << " cos(theta) = " << cos_theta << " g = " << delta * mass_over_req * pow(Units::C,2)/R * obl_approx
  //	      << " log(g) = " << lgrav << std::endl;


  std::cout << "ComputeCurve: Gravity log(g) = " << lgrav << std::endl;

  
  int i_lgrav = (lgrav-13.7)/0.1;
  if (i_lgrav < 1) i_lgrav = 1;

  // New on Nov 2

  //i_lgrav--;

  //
  
  double gvec[5];
  for (int q(0); q<4; q++){
    gvec[q+1] = 13.7+0.1*(i_lgrav-1+q);

    std::cout << "Gravity lgrav=" << lgrav << " gvec["<<q <<"]="<<gvec[q+1] <<std::endl;
    
  }

  bolo = 2.0e9 * 2.404 * Units::C * pow(temperature * Units::EV/(Units::H_PLANCK * Units::C) , 3); // use this one! probably!
  // the e9 in the beginning is for changing T^3 from keV to eV
  // 2.404 comes from evaluating Bradt equation 6.17 (modified, for photon number count units), using the Riemann zeta function for z=3

  //std::cout << "ATMO: Number of Energy bands = " << numbands << std::endl;
  //std::cout << "Eband_upper = " << E_band_upper_1 << std::endl;
  //std::cout << "Eband_lower = " << E_band_lower_1 << std::endl;

  double E_diff = 1.0;
  if ( numbands == 1){
    E0 = curve.para.E0;
  }
  else{
    E_diff = (E_band_upper_1 - E_band_lower_1)/numbands;
    //std::cout << "Final Energy band width = " << E_diff << std::endl;

    // Compute new energy band width

    numbands = curve.cbands;
    //std::cout << "Number of bands to be computed = " << numbands << std::endl;
    E_diff *= (curve.tbands*1.0)/(numbands*1.0);
    //std::cout << "New Energy band width = " << E_diff << std::endl;
  }

  //std::cout << "curve.flags.spectral_model = " << curve.flags.spectral_model << std::endl;
  //std::cout << "curve.flags.beaming_model = " << curve.flags.beaming_model << std::endl;
  
  
  for ( unsigned int i(0); i < numbins; i++ ) { // Compute flux for each phase bin

    if ( curve.dOmega_s[i] != 0.0 ) {

      // Default value is gray = 1
      gray = 1.0; // beaming_model == 0

      // Calculate Beaming for Modified Blackbodies

     if (curve.flags.beaming_model == 6)
	gray = pow( curve.cosbeta[i]*curve.eta[i], 2.0);
      
      if (curve.flags.beaming_model == 7)
	gray = 1.0 - pow( curve.cosbeta[i]*curve.eta[i], 2.0);
      if (curve.flags.beaming_model == 8) // Hopf Function
	gray = 0.42822+0.92236*curve.cosbeta[i]*curve.eta[i]-0.085751*pow(curve.cosbeta[i]*curve.eta[i],2);
	     
      /*******************************************************************/
      /* COMPUTING LIGHT CURVE FOR MONOCHROMATIC ENERGY, p = 0           */
      /*      First computes in [erg/(s cm^2 Hz), converts to            */
      /*      photons/(s cm^2 keV)                                       */
      /*******************************************************************/


      if (curve.flags.spectral_model == 0){ // Monochromatic Observation 


	
	
	if ( curve.flags.beaming_model == 0 || curve.flags.beaming_model == 6 || curve.flags.beaming_model == 7){

	  for (unsigned int p = 0; p<numbands; p++){
	    if (numbands != 1)
	      E0 = (E_band_lower_1+p*E_diff);
	    //if (i==0) std::cout << "p = " << p << " E0bs = " << E0 << " Ediff = " << E_diff << std::endl;
	    curve.f[p][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3)
	      * BlackBody(temperature,E0*redshift/curve.eta[i])*2.0e9 / pow(Units::C * Units::H_PLANCK, 2) * pow(temperature * Units::EV, 3); 
	    curve.f[p][i] *= (1.0 / ( E0 * Units::H_PLANCK )); // Units: photons/(s cm^2 keV)
	    curve.f[p][i] *= E_diff; // Units: photons/(s cm^2)
	    // In case of numbands == 1 E_diff = 1.0, and the units are photon/(s cm^2 keV)
	  }
	}
	      	        	
	if (curve.flags.beaming_model == 10){ // *cole* McPHACC3

	  double cos_theta = curve.cosbeta[i]*curve.eta[i];
	  int theta_index = th_index( cos_theta, &curve);
	  double solidangle = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3);

	  for (unsigned int p = 0; p<numbands; p++){
	    //	    E0 = (E_band_lower_1+(p*curve.tbands/curve.cbands+0.5)*E_diff);
	    E0 = E_band_lower_1 + (p+0.5)*E_diff;
	    //std::cout << "p = " << p << " E = " << E0 << " E_diff=" << E_diff << std::endl;
	    curve.f[p][i] = solidangle
	      * McPHACC3new(E0*redshift/curve.eta[i], 
			    cos_theta, theta_index, curve.para.temperature, lgrav, i_lgrav, gvec, &curve);
	    curve.f[p][i] *= (1.0 / ( E0 * Units::H_PLANCK )); 
	    curve.f[p][i] *= E_diff; // Fake Integration
	    //std::cout << "p = " << p << " i = " << i << " flux = " << curve.f[p][i] << std::endl;
	  }
	}


	/*************************************************/
	if (curve.flags.beaming_model == 11){ // New NSX-H
	  //cout << "Calling NSXHnew!" << endl;	  
	  double cos_theta = curve.cosbeta[i]*curve.eta[i];
	  int theta_index = th_index_nsx( cos_theta, &curve);
	  
	  double solidangle = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3);
	  double Elo, Ehi, Emid;
	  	  
	  for (unsigned int p = 0; p<numbands; p++){
	    
	      /* E0 = E_band_lower_1 + (p + 0.5)*E_diff;
	      curve.f[p][i] = solidangle
	      * NSXHnew(E0*redshift/curve.eta[i], 
			    cos_theta, theta_index, curve.para.temperature, lgrav, i_lgrav, gvec, &curve);
	    curve.f[p][i] *= (1.0 / ( (E0) * Units::H_PLANCK ));	    
	    curve.f[p][i] *= E_diff; // Fake Integration

	    std::cout << "1. Atmo.cpp:  E0 = " << E0 <<  " flux[0][0] = " << curve.f[p][i] << std::endl;
	      */

	    E0 = E_band_lower_1 + (p + 0.5)*E_diff;

	    //Testing
	    if (p==0 && i==0)
	      curve.f[p][i] = solidangle
	      * NSXHnew(E0*redshift/curve.eta[i], 
			cos_theta, theta_index, curve.para.temperature, lgrav, &curve);
	    curve.f[p][i] *= (1.0 / ( (E0) * Units::H_PLANCK ));	    
	    curve.f[p][i] *= E_diff; // Fake Integration

	    /* if (  p==0){
	      std::cout << "2. Atmo.cpp:  E0 = " << E0 <<  " flux[0][0] = " << curve.f[p][i] << std::endl;
	      }*/
	    /*
	     Elo = E_band_lower_1 + (p + 0)*E_diff;
	    Emid = E_band_lower_1 + (p + 0.5)*E_diff;
	    Ehi = E_band_lower_1 + (p + 1)*E_diff;

	    curve.f[p][i] =  NSXHnew(Elo*redshift/curve.eta[i], 
				     cos_theta, theta_index, curve.para.temperature, lgrav, i_lgrav, gvec, &curve);

	    curve.f[p][i] +=  NSXHnew(Emid*redshift/curve.eta[i], 
				     cos_theta, theta_index, curve.para.temperature, lgrav, i_lgrav, gvec, &curve);

	    curve.f[p][i] += NSXHnew(Ehi*redshift/curve.eta[i], 
				     cos_theta, theta_index, curve.para.temperature, lgrav, i_lgrav, gvec, &curve);
	    
				     curve.f[p][i] *= 1.0/3.0 * solidangle * E_diff / ( Emid * Units::H_PLANCK);
	    */

	  }
	}
      } // Spectral_model == 0 

      //std::cout << "Finished with Hydrogen!" << std::endl;

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

	double cos_theta = curve.cosbeta[i]*curve.eta[i];
	int theta_index = th_index( cos_theta, &curve);
	double solidangle = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3);

	for (unsigned int p = 0; p<numbands; p++){
	  
	  if (curve.flags.beaming_model == 10){ // Cole's McPhac File
	    curve.f[p][i] = solidangle
	      * AtmosEBandFlux3new(curve.flags.beaming_model,cos_theta, theta_index,
				   curve.para.temperature,lgrav, i_lgrav, gvec,
				   (E_band_lower_1+p*E_diff)*redshift/curve.eta[i],
				   (E_band_lower_1+(p+1)*E_diff)*redshift/curve.eta[i], curve); // Units: photon/(s cm^2)    
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
  





    
    return curve;

} // end ComputeCurve






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






// NSX - Hydrogen Atmosphere Computed by Wynn Ho
// Calculate the final interpolated intensity
// This "new" version takes into account that the energy is really the ratio: E/kT
double NSXHnew(double E, double cos_theta, int theta_index, double T, double lgrav,  class LightCurve* curve){

  //std::cout << "Welcome to NSX! E=" << E << std::endl;

  class LightCurve mexmcc;
  mexmcc = (*curve);

  double lt, ener_index;
	
  double I_int[5], J[5], K[5], L(0.0);
  int i_f, i_lt,  first_inte;  
  int ii_mu = theta_index;

  // CALCULATE LOG(T) AND CORRECT INDEX
  lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));

  // CALCULATE LOG(E/T) 
   double logET = log10(E/T);

   
  //TEST
   /*
  lt = 5.1;
  lgrav = 15;
  cos_theta = 0.5;
  ii_mu = 43;
  logET = 1.0;
   */
  // End TEST

  
  // Find the correct temperature range
    i_lt = (lt-5.1)/0.05 ; //if we need to load 1st temperature, i_lt = 0. this is discrete math
    if (i_lt < 1) i_lt = 0;
    if (i_lt > 35) i_lt = 35;


    // Find the correct gravity index
    int i_lgrav = (lgrav-13.7)/0.1;
    if (i_lgrav > 2) i_lgrav--;
    if (i_lgrav > 10) i_lgrav=10;
    
    /*
    std::cout << "log(T) = " << lt << " i_lt = " << i_lt << " log(T) = " << mexmcc.mclogTeff[i_lt] << std::endl;
    std::cout << "log(g) = " << lgrav << " i_lg=" << i_lgrav << " log(g) = " << mexmcc.mclogg[i_lgrav] << std::endl;
    std::cout << "mu = " << cos_theta << " imu=" << ii_mu << " mu = " << mexmcc.mccangl[ii_mu] << std::endl;
    */
    
    double ivec[5][5][5][5];
    double err, err1;
    double tvec[5], gvec[5];
    
    int npt(4), tpt(4), gpt(4), mpt(4);
    
    if ( npt > 4)
      L = 0.0;
    else{



      
    //Find proper freqency choice
    ener_index = (logET + 1.30)/0.02;
    i_f = (int) ener_index;
    std::cout << "ener_index = " << ener_index << " i_f = " << i_f << std::endl;
    // SMM 2022 if (npt==2 || npt==4) i_f +=1;
    if (i_f < 0) i_f = 1;
    // Original version
    if (i_f > 162) i_f=162;
    // SMM: 20221117 Change to:
    if (i_f > 164) i_f=164;
    /*std::cout << "log(E/T) = " << log10(E/T)
	      << " i_f = " << i_f
	      << " log(E/T) = " << mexmcc.mcloget[i_f]
	      << std::endl;*/
    // int ii_mu = th_index( cos_theta, &mexmcc);

  
    
    
    for (int r(0); r<tpt; r++){
      //std::cout << "r = " << r << std::endl << std::endl;
      tvec[r+1] =  5.1+0.05*(i_lt-1+r+1);
      for (int q(0); q<gpt; q++){
	gvec[q+1] =  mexmcc.mclogg[i_lgrav+q];
	//std::cout << "q = " << q << std::endl;
	for (int k(0); k<mpt; k++){
	  //std::cout << "k = " << k << std::endl;
	  for( int j(0); j<npt; j++){
	   	    
	    first_inte = (i_lt + r) * mexmcc.NlogE * mexmcc.Nmu * mexmcc.Nlogg
	      + (i_lgrav+q) * mexmcc.NlogE * mexmcc.Nmu
	      + (i_f + j) * mexmcc.Nmu
	      + (ii_mu + k) ;
	    
	    // linear interpolation
	    ivec[r][q][k][j+1] = mexmcc.mccinte[first_inte];

	    // logarithmic interpolation
	    //ivec[r][q][k][j+1] = log10(mexmcc.mccinte[first_inte]);
	    /*
	      std::cout << "j = " << j
			<< "first_inte = " << first_inte
			<< " tvec = " << tvec[r+1]
			<< " gvec = " << gvec[q+1]
			<< " costheta = " << mexmcc.mccangl[ii_mu+k]
		      << " evec = " << mexmcc.mcloget[i_f+j] 
		      << " ivec = " << ivec[r][q][k][j+1] << std::endl;*/
	  }
	  
	  // linear interpolation

	  // E=T*100.0;// TEST
	  
	  I_int[k+1] = polint(&mexmcc.mcloget[i_f-1],ivec[r][q][k],npt,logET,&err1);
	  /*if (r == 0 && q == 0 ){
	    I_int[k+1] = polint(&mexmcc.mcloget[i_f-1],ivec[r][q][k],npt,log10(E/T),&err1);		  
	    std::cout << "mu = " << mexmcc.mccangl[ii_mu+k]
	      << " Interpolated value = " << I_int[k+1] << std::endl;
	      }*/
	  //std::cout << "mu = " << mexmcc.mccangl[ii_mu+k]
	  // << " Interpolated value = " << I_int[k+1] << std::endl;
	}
	// Intepolate over mu for fixed Teff, gravity
	J[q+1] = polint(&mexmcc.mccangl[ii_mu-1],I_int,mpt,cos_theta,&err);

	/*if (r == 0){
	  J[q+1] = polint(&mexmcc.mccangl[ii_mu-1],I_int,mpt,cos_theta,&err);
	  std::cout << "logg = " << mexmcc.mclogg[i_lgrav+q] << " Interp = " << J[q+1] << std::endl;
	  }*/
	//std::cout << "logg = " << mexmcc.mclogg[i_lgrav+q] << " Interp = " << J[q+1] << std::endl;
      }
      // Interpolate over logg for fixed Teff
      //K[r+1] = polint(gvec,J,gpt,lgrav,&err);
      
      K[r+1] = polint(gvec,J,gpt,lgrav,&err);
	
      //std::cout << " Interp logT=" << tvec[r+1] << " IntK= " << K[r+1] << std::endl;
      //}
    }

    // Interpolate over log(Teff)
    L = polint(tvec,K,tpt,lt,&err);

    //std::cout << "Interpolated I = " << L << std::endl;
    
    
    }
   
    /* std::cout << "Log(I/T^3) = " << L << " I/T^3 = " << pow(10.0,L)
	      << " I = " << pow(10.0,L+lt*3.0)
	      << std::endl;*/
   
    
    return pow(10.0,L+lt*3.0);
	
}














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



