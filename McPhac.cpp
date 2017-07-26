/***************************************************************************************/
/*                                     McPhac.cpp

    This holds the routines used in Spot.cpp and atmosphere routines that are 
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
#include "McPhac.h"
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

int th_index(double cos_theta, class LightCurve* mexmcc){

  int n_mu, i_mu;

    //Find proper mu choice
    n_mu = 1;
    while (cos_theta > mexmcc->mccangl[n_mu] && n_mu < 50){
    	n_mu += 1;
    }
    i_mu = n_mu - 1;
    // n_mu += -1;

    int ii_mu(i_mu-1);
    if (ii_mu < 0)
      ii_mu = 0;
    if (ii_mu > 46)
      ii_mu = 46;


    return ii_mu;

}



// Calculate the final interpolated intensity
// This "new" version takes into account that the energy is really the ratio: E/kT
double McPHACC3new(double E, double cos_theta, int theta_index, double T, double lgrav, double gvec[4], class LightCurve* curve){


  class LightCurve mexmcc;
  mexmcc = (*curve);

	double lt, ener_index;
	//double e0, e1, th0, th1, t0, t1;
	
	double I_int[4], J[4], K[4], L(0.0);
	int i_f, i_lt, i_lgrav, first_inte;
	std::ofstream out1, out2;      // output stream; printing information to the output file
  
      	int ii_mu = theta_index;


    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));

    // Find the correct temperature range
    i_lt = (lt-5.1)/0.05 ; //if we need to load 1st temperature, i_lt = 0. this is discrete math
    if (i_lt < 1) i_lt = 1;
    if (i_lt > 26) i_lt = 26;


    i_lgrav = (lgrav-13.7)/0.1;
    if (i_lgrav < 1) i_lgrav = 1;

 

    // Now do a npt interpolation
    //double evec[5];
    double ivec[4][4][4][4];
    double err, err1;
    double tvec[4];

    
    int npt(3);
    
    //Find proper freqency choice
    ener_index = (log10(E/T) + 1.30369)/0.0338;
    i_f = (int) ener_index;  
    if (npt==2 || npt==4) i_f +=1;
    if (i_f < 1) i_f = 1;
    if (i_f > 95) i_f=95;

    // int ii_mu = th_index( cos_theta, &mexmcc);

    for (int r(0); r<3; r++){
      tvec[r+1] =  5.1+0.05*(i_lt-1+r);
      //t0 = tvec[r+1];
      //std::cout << "tvec[r]=" << tvec[r+1] << std::endl;
      for (int q(0); q<3; q++){
	//gvec[q+1] = 13.7+0.1*(i_lgrav-1+q);
	//std::cout << "logg = " << gvec[q+1] << std::endl;
	for (int k(0); k<3; k++){
	  //muvec[k+1] =  mexmcc.mccangl[ii_mu+k];
	  for( int j(0); j<npt; j++){
	    //evec[j] = pow(10,mexmcc.mcloget[i_f-1+j]);
	    //evec[j+1] = mexmcc.mcloget[i_f-1+j];
	    first_inte = ((i_lt-1+r)*11 + i_lgrav-1+q) * 5000 + (i_f-1+j) * 50 + (ii_mu + k);

	    // linear interpolation
	    ivec[r][q][k][j+1] = mexmcc.mccinte[first_inte];

	    // logarithmic interpolation
	    //ivec[r][q][k][j+1] = log10(mexmcc.mccinte[first_inte]);	    
	  }
	  
	  // linear interpolation
	  I_int[k+1] = polint(&mexmcc.mcloget[i_f-2],ivec[r][q][k],npt,log10(E/T),&err1);
	 
	  // logarithmic interpolation
	  //I_int[k+1] = pow(10.0,polint(&mexmcc.mcloget[i_f-2],ivec[r][q][k],npt,log10(E/T),&err1));	
	}
	// Intepolate over mu for fixed Teff, gravity
	J[q+1] = polint(&mexmcc.mccangl[ii_mu-1],I_int,3,cos_theta,&err);
      }
      // Interpolate over logg for fixed Teff
      K[r+1] = polint(gvec,J,3,lgrav,&err);   
    }

    // Interpolate over log(Teff)
    L = polint(tvec,K,3,lt,&err)*pow(10.0,lt*3.0);

    //if (cos_theta < 0.015629) L = 0;
    
    return L;
}








// This is version makes use of Cole's version of McPhac
double AtmosEBandFlux3new( unsigned int model, double cos_theta, int theta_index, double T, double lgrav, double gvec[4], double E1, double E2, class LightCurve mexmcc){

    int e1_dex(0);              // largest energy index that's smaller than E1
    int e2_dex(0);              // largest energy index that's smaller than E2
    int n_steps(0);             // number of energy points within band
    double flux(0.0);           // total integrated flux

    //double gvec[4];


    double ener_spacing = pow(10.0,0.0338);
    double first_ener = 0.04969469;
    double ener_index = log10(E1/T / first_ener) / log10(ener_spacing);
    e1_dex = (int) ener_index;
    ener_index = log10(E2/T / first_ener) / log10(ener_spacing);
    e2_dex = (int) ener_index;

    n_steps = 1* (e2_dex - e1_dex);

    // int theta_index = th_index(cos_theta,&mexmcc);
    n_steps = 4;


    if (n_steps == 0){ // zero energy points within bandwidth: (4.1.3) one trapzoid
      //cout << "0 steps" << endl;
      flux = (E2 - E1) / 2.0 * (McPHACC3new(E1,cos_theta,theta_index, T, lgrav, gvec, &mexmcc) / E1 + McPHACC3new(E2,cos_theta,theta_index, T, lgrav, gvec, &mexmcc) / E2);
    }
    if (n_steps == 1){ // one energy points within bandwidth: (4.1.3) two trapzoids
      //cout << "1 step" << endl;
        int e_dex = e1_dex+1; // index of the energy point
        double e_m = first_ener*pow(ener_spacing,e_dex)*T; // energy point in keV
	
	e_m = 0.5*(E2+E1);

	
	double counts_m = McPHACC3new(e_m,cos_theta,theta_index, T,lgrav, gvec, &mexmcc) / e_m;

        flux = (e_m - E1) / 2 * (McPHACC3new(E1,cos_theta,theta_index, T, lgrav, gvec, &mexmcc) / E1 
				 + counts_m);  // first trapezoid
        flux += (E2 - e_m) / 2 * (counts_m + McPHACC3new(E2,cos_theta, theta_index, T,lgrav,gvec, &mexmcc) / E2); // second trapezoid
    }
    if (n_steps >= 2){ // two energy points within bandwidth: (4.1.3) three trapzoids
      //cout << "n steps= " << n_steps << endl;
     
       double e_step = (E2-E1)/(n_steps+1.0);
       double e_l;
       double counts_l;

       flux = e_step * 0.5 * McPHACC3new(E1,cos_theta, theta_index, T, lgrav,gvec, &mexmcc) / E1;
       //cout << "flux = " << flux << endl;

       for (int i(1); i<=n_steps; i++){

	 e_l = E1 + i*e_step;
	 counts_l = McPHACC3new(e_l,cos_theta,theta_index, T, lgrav,gvec, &mexmcc) / e_l;
	 flux += e_step * (counts_l);
	 //cout << "flux = " << flux << endl;

       }

       flux += e_step * 0.5 *  McPHACC3new(E2,cos_theta, theta_index, T,lgrav,gvec, &mexmcc) / E2;
       // cout << "flux = " << flux << endl;

    }
   
    flux = flux/Units::H_PLANCK;
    return flux;
} // new version of energy integration



