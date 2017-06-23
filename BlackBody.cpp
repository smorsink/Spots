// Routines for Computing BlackBody 

#include "matpack.h"
#include <exception>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unistd.h>
#include "BlackBody.h"
#include "Exception.h"
#include "Units.h"
#include "Struct.h"
#include "time.h"
#include "interp.h"
#include <stdio.h>
using namespace std;

/**************************************************************************************/
/* Blackbody:                                                                         */
/*           computes the monochromatic blackbody flux in units of erg/cm^2			  */
/*																					  */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E = monochromatic energy in keV * redshift / eta                             */
/**************************************************************************************/
double BlackBody( double T, double E ) {   // Blackbody flux in units of erg/cm^2
  //return ( 2.0e9 / pow(Units::C * Units::H_PLANCK, 2) * pow(E * Units::EV, 3) / (exp(E/T) - 1) ); 

     return (pow(E/T,3.0)/(exp(E/T) - 1));

    // the e9 is to switch E from keV to eV; Units::EV gets it from eV to erg, since it's first computed in erg units.
    // the switch from erg units to photon count units happens above just after this is called.
} // end Blackbody


/**************************************************************************************/
/* BlackBodyTab:                                                                      */
/* Reads in the tabulated blackbody file                                              */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E = monochromatic energy in keV * redshift / eta                             */
/**************************************************************************************/
double BlackBodyTabLogLog( double T, double E, class LightCurve mexmcc ) {   // Blackbody flux in units of erg/cm^2

  int i_f, npt;
  
  npt =  4;

  double loget( log10(E/T) );
  double intensity, err;

  double  ener_index = ( loget + 1.30369)/0.0338;
  i_f = (int) ener_index;  
    if (npt==2 || npt==4) i_f +=1;
    if (i_f < 1) i_f = 1;
    if (i_f > 95) i_f=95;

    double evec[7];
    double ivec[7];
    // int first_inte;

    for( int j(0); j<npt; j++){
      evec[j+1] = mexmcc.mcloget[i_f-1+j];
      ivec[j+1] = log10(mexmcc.mccinte[i_f-1+j]);       
    }

    intensity =  pow(10.0,polint( evec, ivec, npt, loget, &err )); 

     return (intensity);

     // return ( 2.0e9 / pow(Units::C * Units::H_PLANCK, 2) * pow(T * Units::EV, 3) *  intensity); 
    // the e9 is to switch E from keV to eV; Units::EV gets it from eV to erg, since it's first computed in erg units.
    // the switch from erg units to photon count units happens above just after this is called.
} // end Blackbody

/**************************************************************************************/
/* BlackBodyTab:                                                                      */
/* Reads in the tabulated blackbody file                                              */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E = monochromatic energy in keV * redshift / eta                             */
/**************************************************************************************/
double BlackBodyTabLogLinear( double T, double E, class LightCurve mexmcc ) {   // Blackbody flux in units of erg/cm^2

  int i_f, npt;
  
  npt =  4;

  double loget( log10(E/T) );
  double intensity, err;

  double  ener_index = ( loget + 1.30369)/0.0338;
  i_f = (int) ener_index;  
    if (npt==2 || npt==4) i_f +=1;
    if (i_f < 1) i_f = 1;
    if (i_f > 95) i_f=95;

    double evec[7];
    double ivec[7];
    // int first_inte;

    for( int j(0); j<npt; j++){
      evec[j+1] = mexmcc.mcloget[i_f-1+j];
      ivec[j+1] = mexmcc.mccinte[i_f-1+j];       
    }

    intensity =  polint( evec, ivec, npt, loget, &err ); 

     return (intensity);

     // return ( 2.0e9 / pow(Units::C * Units::H_PLANCK, 2) * pow(T * Units::EV, 3) *  intensity); 
    // the e9 is to switch E from keV to eV; Units::EV gets it from eV to erg, since it's first computed in erg units.
    // the switch from erg units to photon count units happens above just after this is called.
} // end Blackbody




/**************************************************************************************/
/* HopfTab:                                                                      */
/* Reads in the tabulated blackbody file                                              */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E = monochromatic energy in keV * redshift / eta                             */
/**************************************************************************************/
double HopfTab( double cosalpha, class LightCurve mexmcc ) {   // Blackbody flux in units of erg/cm^2

  int i_mu,n_mu, npt;
  
  npt =  4;

  //Find proper mu choice
  n_mu = 1;
  while (cosalpha > mexmcc.mccangl[n_mu] && n_mu < 50){
    n_mu += 1;
  }
  i_mu = n_mu - 1;

  if (i_mu==0)
    i_mu = 1;

  double intensity, err;

  double muvec[7];
  double ivec[7];
  //int first_inte;

  for( int j(0); j<npt; j++){
    muvec[j+1] = mexmcc.mccangl[i_mu-1+j];
    ivec[j+1] = mexmcc.mccinte[i_mu-1+j];  
    //if (cosalpha < 0.01)
    //std::cout << "i_mu=" << i_mu << " j = " << j << " mu = " << muvec[j+1] << " ivec = " << ivec[j+1] << std::endl;
  }

  // if (cosalpha < 0.01)
  //intensity =  printpolint( muvec, ivec, npt, cosalpha, &err ); 
  //else
    intensity =  polint( muvec, ivec, npt, cosalpha, &err ); 

  return (intensity);

  //     return ( 2.0e9 / pow(Units::C * Units::H_PLANCK, 2) * pow(T * Units::EV, 3) *  intensity); 
    // the e9 is to switch E from keV to eV; Units::EV gets it from eV to erg, since it's first computed in erg units.
    // the switch from erg units to photon count units happens above just after this is called.
} // end Blackbody


/**************************************************************************************/
/* BlackBodyTab:                                                                      */
/* Reads in the tabulated blackbody file                                              */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E = monochromatic energy in keV * redshift / eta                             */
/**************************************************************************************/
double BlackBodyHopfTab( double T, double E, double cosalpha, class LightCurve mexmcc ) {   // Blackbody flux in units of erg/cm^2

  int i_f, npt;
  int i_mu,n_mu;
  
  npt =  4;


  
  double loget( log10(E/T) );
  double intensity, err;

  double  ener_index = ( loget + 1.30369)/0.0338;
  i_f = (int) ener_index;  
    if (npt==2 || npt==4) i_f +=1;
    if (i_f < 1) i_f = 1;
    if (i_f > 95) i_f=95;

  //Find proper mu choice
  n_mu = 1;
  while (cosalpha > mexmcc.mccangl[n_mu] && n_mu < 50){
    n_mu += 1;
  }
  i_mu = n_mu - 1;

  if (i_mu==0)
    i_mu = 1;


    double evec[7];
    double muvec[7];
    double ivec[7][7];

    double I[7];
    
    int index;

    for( int k(0); k<npt; k++){
      muvec[k+1] = mexmcc.mccangl[i_mu-1+k];
      //std::cout << "k=" << k << " mu = " <<  mexmcc.mccangl[i_mu-1+k] << std::endl;
      for( int j(0); j<npt; j++){
	evec[j+1] = mexmcc.mcloget[i_f-1+j];
	index = (i_f+j-1)*50 + (i_mu-1+k);
	ivec[k+1][j+1] = log10(mexmcc.mccinte[index]);   
	/*std::cout << "j="<< j
		  << " evec[j+1] = " << evec[j+1]
		  << " index= " << index
		  << " intensity = " <<  mexmcc.mccinte[index]
		  << std::endl;*/
      }
      I[k+1] =  pow(10.0,polint( evec, ivec[k+1], npt, loget, &err )); 
    }
    intensity = polint(muvec,I,npt,cosalpha, &err);

    return (intensity);

    //return ( 2.0e9 / pow(Units::C * Units::H_PLANCK, 2) * pow(T * Units::EV, 3) *  intensity); 
    // the e9 is to switch E from keV to eV; Units::EV gets it from eV to erg, since it's first computed in erg units.
    // the switch from erg units to photon count units happens above just after this is called.
} // end Blackbody




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

double Hopf( double cosalpha) {

  return (0.42822+0.92236*cosalpha-0.085751*pow(cosalpha,2)) ;

}
