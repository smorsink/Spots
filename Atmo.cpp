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
#include "Hydrogen.h"
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
  
  double E0, E1, E2, DeltaE, DeltaLogE;

  unsigned int numbins(MAX_NUMBINS);  // Time bins of light curve (usually 128)
  unsigned int numbands(0);  // Number of Energy Bands

  bool logEflag(false);

  /*********************/
  /* SETTING THINGS UP */
  /*********************/
    
  // One monochromatic energy, hardwired value, in keV
  //    E_mono = 1.0;
  //std::cout << "ComputeCurve: Starting ComputeCurve" << std::endl;
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
  DeltaLogE = curve.para.DeltaLogE;
  logEflag = curve.flags.logEflag;


  // DeltaE = (E_band_upper_1 - E_band_lower_1)/(numbands-1.0);
  /*     std::cout << "COMPUTE CURVE: Lowest energy = " << E_band_lower_1
		<< " = " << curve.elo[0] 
		<< " Highest energy = " << E_band_upper_1
		<< " numbands = curve.numbands = " << numbands
		<< " Delta(E) = " << DeltaE
		<< std::endl;*/

  redshift = 1.0 / sqrt( 1 - 2.0 * mass_over_r);
 
  double R = Units::nounits_to_cgs(curve.para.req, Units::LENGTH);
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


  //std::cout << "ComputeCurve: Gravity log(g) = " << lgrav << std::endl;

  
  int i_lgrav = (lgrav-13.7)/0.1;
  if (i_lgrav < 1) i_lgrav = 1;


  
  double gvec[5];
  for (int q(0); q<4; q++){
    gvec[q+1] = 13.7+0.1*(i_lgrav-1+q);    
  }

  bolo = 2.0e9 * 2.404 * Units::C * pow(temperature * Units::EV/(Units::H_PLANCK * Units::C) , 3); // use this one! probably!
  // the e9 in the beginning is for changing T^3 from keV to eV
  // 2.404 comes from evaluating Bradt equation 6.17 (modified, for photon number count units), using the Riemann zeta function for z=3

  
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

      double E_diff = DeltaE;

      if (curve.flags.spectral_model == 0){ // Monochromatic Observation 


	// Check units for Temperature!!!!
	
	if ( curve.flags.beaming_model == 0 || curve.flags.beaming_model == 6 || curve.flags.beaming_model == 7){

	  if (curve.flags.kelvin){
	    temperature *= Units::K_BOLTZ/Units::EV*1e-3 ; // Convert T to keV
	  }
	  
	  for (unsigned int p = 0; p<numbands; p++){
	    if (numbands != 1)
	      E0 = (E_band_lower_1+p*E_diff);

	    if (logEflag){ // logarithmic
	      E0 = curve.para.E_band_lower_1 * pow(10,p*DeltaLogE);
	    }
	    else{
	      E0 = curve.para.E_band_lower_1 + (p+0.5)*E_diff;
	    }
	    
	   
	    curve.f[p][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3)
	      * BlackBody(temperature,E0*redshift/curve.eta[i])*2.0e9 / pow(Units::C * Units::H_PLANCK, 2) * pow(temperature * Units::EV, 3);
	    
	    //curve.f[p][i] *= (1.0 / ( E0 * Units::H_PLANCK )); // Units: photons/(s cm^2 keV)
	    //curve.f[p][i] *= E_diff; // Units: photons/(s cm^2)
	    // In case of numbands == 1 E_diff = 1.0, and the units are photon/(s cm^2 keV)

	    if (logEflag){ // Logarithmic
	      if (p==0 && i==0)
		std::cout << " Logarithmic! " << std::endl;
	      curve.f[p][i] *= (1.0 / ( E0 * Units::H_PLANCK ));	    
	      curve.f[p][i] *=  0.01 ; // Fake Integration	  
	    }
	    else{ // Linear
	      curve.f[p][i] *= (1.0 / ( (E0) * Units::H_PLANCK ));	    
	      curve.f[p][i] *= DeltaE; // Fake Integration	   
	    }

	    
	  }
	}
	      	        	

   
	/*************************************************/
	if (curve.flags.beaming_model == 11){ // New NSX-H
	  //std::cout << "Calling NSXHnew!" << endl;	  
	  double cos_theta = curve.cosbeta[i]*curve.eta[i];
	  int theta_index = th_index_nsx( cos_theta, &curve);
	  
	  double solidangle = curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3);

	  //std::cout << "Temperature = " << curve.para.temperature << std::endl;
	  
	  for (unsigned int p = 0; p<=numbands; p++){

	    // E0 is the observed photon energy in keV
	    E0 = curve.elo[p]; // Use the lower limit of the energy band;

	    // NSXHnew returns specific intensity in units of erg (s Hz)^{-1} cm^-2 str^-1
	    
	    
	    curve.f[p][i] = solidangle
	      * NSXHnew(E0*redshift/curve.eta[i], 
			cos_theta, theta_index, curve.para.temperature, lgrav, &curve);

	    if (logEflag){ // Logarithmic
	      if (p==0 && i==0)
		std::cout << " Logarithmic! " << std::endl;
	      curve.f[p][i] *= (1.0 / ( E0 * Units::H_PLANCK ));	    
	      curve.f[p][i] *=  0.01 ; // Fake Integration	  
	    }
	    else{ // Linear
	      curve.f[p][i] *= (1.0 / ( (E0) * Units::H_PLANCK ));	    
	      curve.f[p][i] *= DeltaE/numbins; // Fake Integration
	      //curve.f[p][i] *= (1.0 / (E0 * 1e3 * Units::EV )); // divide by photon energy in erg
	      // f now has units of number of photons/(cm^2)
	      // multiply by delta t x delta f = time bin x photon frequency bin
	      //curve.f[p][i] *= (DeltaE * 1e3 * Units::EV / Units::H_PLANCK)/(numbins);
	    }

	    /*if (p==0 && i==0)
	      std::cout << "ATMO (H): Omega = " << curve.dOmega_s[i]
			<< " solid angle = " << solidangle
			<< " E0 = " << E0
			<< " DeltaLogE = " << DeltaLogE
			<< " DeltaE = " << DeltaE
			<< " flux = " << curve.f[p][i]
			<< " log(10) = " << log(10.0)
			<< std::endl;*/
	    
	  }
	}
      } // Spectral_model == 0 

      //std::cout << "Finished with Hydrogen!" << std::endl;

      /*******************************************************************/
      /* COMPUTING BLACKBODY LIGHT CURVE FOR INTEGRATED FLUX             */
      /* Units: photons/(cm^2 s)                                         */
      /*******************************************************************/

      if (curve.flags.spectral_model == 2){ // Integrated Flux for Modified Blackbodies
	
	double Elo, Ehi;
	
	  if (curve.flags.kelvin){
	    temperature *= Units::K_BOLTZ/Units::EV*1e-3 ; // Convert T to keV
	  }

	for (unsigned int p = 0; p<numbands; p++){

	    if (logEflag){ // logarithmic
	      Elo = curve.para.E_band_lower_1 * pow(10,(p-0.5)*DeltaLogE);
	      Ehi = curve.para.E_band_lower_1 * pow(10,(p+0.5)*DeltaLogE);
	    }
	    else{
	      Elo = curve.para.E_band_lower_1 + (p+0.5)*E_diff;
	      Elo = curve.para.E_band_lower_1 + (p+1.5)*E_diff;
	    }

	   
	  curve.f[p][i] = gray * curve.dOmega_s[i] * pow(curve.eta[i],4) * pow(redshift,-3) 
	    * EnergyBandFlux(temperature, (Elo)*redshift/curve.eta[i], (Ehi)*redshift/curve.eta[i]); // Units: photon/(s cm^2)
	  
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



