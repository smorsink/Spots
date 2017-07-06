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



















// Calculate the final interpolated intensity
// This "new" version takes into account that the energy is really the ratio: E/kT
double McPHACC3new(double E, double cos_theta, double T, double lgrav, class LightCurve mexmcc){
	double lt, ener_index;
	//double e0, e1, th0, th1, t0, t1;
	double t0;
	double I_int[9], J[5], K[3], L(0.0);
	int i_f, i_lt, i_lgrav, i_mu, n_mu, first_inte;
	//double mu0, mu1;
	std::ofstream out1, out2;      // output stream; printing information to the output file
  
       
	
	out2.open("InterpolationTests/interp.txt", std::ios_base::app);
	out2 << "#log(E/T) E (keV) I/T^3       I       T (keV)  logT   logg   cosa   Interp  Error(I/T^3)" << std::endl;



    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));

    // This part is purely for testing against the old mcphac file.
    /*    
    lt = 5.8;
    T = pow(10,lt) * 1e-3/(Units::EV / Units::K_BOLTZ);
    lgrav = 14.1;
    */

    //    std::cout << "McPhac logt = " << lt << std::endl;

    // Find the correct temperature range
    i_lt = (lt-5.1)/0.05 ; //if we need to load 1st temperature, i_lt = 0. this is discrete math
    if (i_lt < 1) i_lt = 1;
    if (i_lt > 26) i_lt = 26;

    // std::cout << "log(T) = " << lt << std::endl;
    //std::cout << "i_lt = " << i_lt << std::endl;

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



    // Now do a npt interpolation
    double evec[7];
    double ivec[7][7][7][7];
    double err, err1;
    double muvec[7];
    double gvec[7];
    double tvec[7];

    // for (int npt(2);npt<7;npt++){
    
    int npt(4);
    
    //Find proper freqency choice
    ener_index = (log10(E/T) + 1.30369)/0.0338;
    i_f = (int) ener_index;  
    if (npt==2 || npt==4) i_f +=1;
    if (i_f < 1) i_f = 1;
    if (i_f > 95) i_f=95;

    int ii_mu(i_mu-1);
    if (ii_mu < 0)
      ii_mu = 0;
    if (ii_mu > 46)
      ii_mu = 46;

    for (int r(0); r<4; r++){
      tvec[r+1] =  5.1+0.05*(i_lt-1+r);
      t0 = tvec[r+1];
      //std::cout << "tvec[r]=" << tvec[r+1] << std::endl;
      for (int q(0); q<4; q++){
	gvec[q+1] = 13.7+0.1*(i_lgrav-1+q);
	//std::cout << "logg = " << gvec[q+1] << std::endl;
	for (int k(0); k<4; k++){
	  muvec[k+1] =  mexmcc.mccangl[ii_mu+k];
	  //muvec[k+1] =  acos(mexmcc.mccangl[ii_mu+k]);
	  //std::cout << "muvec[k] = " << muvec[k+1] << std::endl;
	  // Interpolate over Energy for fixed Teff, gravity, and mu
	  //std::cout << "E = " << E << " T = " << T << " E/T =" << E/T << std::endl;  
	  //std::cout << "Log10(E/T) = " << log10(E/T) << std::endl;

	  /* if ( fabs(tvec[r+1] - 5.8) < 1e-10 && (gvec[q+1]==14.1) && (muvec[k+1]==.2323)){

	    out1.open("InterpolationTests/LowT.txt", std::ios_base::trunc);

	    for (int i(0); i<100; i++){

	      
	      double log_e_t = mexmcc.mcloget[1-1+i];
	      double energy = pow(10,log_e_t) * pow(10,tvec[r+1]) * 1e-3/(Units::EV / Units::K_BOLTZ);
	      first_inte = ((i_lt-1+r)*11 + i_lgrav-1+q) * 5000 + (1-1+i) * 50 + (ii_mu + k) +1 -1;
	      double iovert3 = mexmcc.mccinte[first_inte];

	      out1 << "#log(E/T) E (keV)    I/T^3        I    T (keV)  logT   logg   cosa index" << std::endl;
	      out1 << log_e_t << "  " 
			<< energy << " "
			<< iovert3 << " "
			<< iovert3 * pow(10.0,t0*3.0) << " "
			<<  pow(10,tvec[r+1]) * 1e-3/(Units::EV / Units::K_BOLTZ) << " "
			<< tvec[r+1] << "    "
			<< gvec[q+1] << "   "
			<< muvec[k+1] << " "
			<< i
			<< std::endl;


	    }
	    out1.close();*/
	    
	


	  for( int j(0); j<npt; j++){
	    //evec[j] = pow(10,mexmcc.mcloget[i_f-1+j]);
	    evec[j+1] = mexmcc.mcloget[i_f-1+j];
	    first_inte = ((i_lt-1+r)*11 + i_lgrav-1+q) * 5000 + (i_f-1+j) * 50 + (ii_mu + k) +1 -1;

	    //  ivec[r][q][k][j+1] = mexmcc.mccinte[first_inte]*pow(10.0,t0*3.0);

	    ivec[r][q][k][j+1] = mexmcc.mccinte[first_inte];

	    //ivec[r][q][k][j+1] = log10(mexmcc.mccinte[first_inte]);

	    /* std::cout << "j = " << j<< std::endl;
	    std::cout << "E = " << E << " T = " << T << " E/T =" << E/T << std::endl;  
	    std::cout << "Log10(E/T) = " << log10(E/T) << std::endl;
	    std::cout << "evec[j] = " << evec[j+1] 
	    	      << " ivec[j] = " << ivec[r][q][k][j+1]
		      << " mcint = " << mexmcc.mccinte[first_inte]
	    	      << std::endl;*/
	    
	  }
	  
	  //	  I_int[k+1] = polint(evec,ivec[r][q][k],npt,log10(E/T),&err1)*pow(10.0,t0*3.0);
	  I_int[k+1] = polint(evec,ivec[r][q][k],npt,log10(E/T),&err1);
	
	  /* if ( fabs(tvec[r+1] - 6.4) < 1e-10 && (gvec[q+1]==14.1) && (muvec[k+1]==.2323)){

	     I_int[k+1] = polint(evec,ivec[r][q][k],npt,log10(E/T),&err1)*pow(10.0,t0*3.0);
	    //I_int[k+1] = printpolint(evec,ivec[r][q][k],npt,log10(E/T),&err1);

	    // out2 << "#log(E/T) E (keV) I/T^3       I       T (keV)  logT   logg   cosa   Interp  Error(I/T^3)" << std::endl;
	    
	    out2 
		    << log10(E/T) << " "
		    << E << "   "
		    << polint(evec,ivec[r][q][k],npt,log10(E/T),&err1) << "  "
		    << I_int[k+1] << " "
		    << T << " "
		    << tvec[r+1] << "    "
		    << gvec[q+1] << "   "
		    << muvec[k+1] << " "
		    << npt << "       " 
		    << err1
		    << endl;
	  }
	  else
	    I_int[k+1] = polint(evec,ivec[r][q][k],npt,log10(E/T),&err1)*pow(10.0,t0*3.0);*/

	  //if (isnan(I_int[k])){
	  /* 
	    cout << evec[1] << " " << evec[2] << " " << evec[3] << " " << evec[4] << " " << log10(E/T) << " " << ivec[r][q][k][1] << " " << ivec[r][q][k][2] << " " << ivec[r][q][k][3] << " " << ivec[r][q][k][4] << endl;
	  	I_int[k] = printpolint(evec,ivec[r][q][k],4,log10(E/T),&err);
	  	cout << "Interpolated Intensity = " << I_int[k]
	  	  << "err = " << err << endl;
	  */
		//}
	  /*  cout << "0 mu[k] = " << muvec[k+1]
	     <<" E=" << E/T << " New 4pt Interpolated: I = " << I_int[k+1] 
	     << " err = " << err << endl;*/
	}
	// Intepolate over mu for fixed Teff, gravity
	J[q+1] = polint(muvec,I_int,4,cos_theta,&err);
	/*	cout 
		  << " logg = " << gvec[q+1]
		  << " costheta = " << cos_theta << " Interpolated I = " << J[q+1]
		  << " err = " << err << std::endl;*/
      }
      // Interpolate over logg for fixed Teff
      K[r+1] = polint(gvec,J,4,lgrav,&err);
      /*cout << " logg = " << lgrav 
         << " Interpolated I = " 
         <<  K[r+1] 
         << " err = " << err << endl;   */   
    }

    L = polint(tvec,K,4,lt,&err)*pow(10.0,lt*3.0);
    // L = polint(tvec,K,4,lt,&err);


    
	    out2 
		    << log10(E/T) << " "
		    << E << "   "
		    << L * pow(10.0,-lt*3.0) << "  "
		    << L << " "
		    << T << " "
		    << lt << "    "
		    << lgrav << "   "
		    << cos_theta << " "
		    << npt << "       " 
		    << err1
		    << endl;



    
    
    /*    std::cout << " log(T_eff) = " << lt 
    	      << " Interpolated I = " << L
    	      << " err = " << err << std::endl;
    */

  
    //  cout << L << endl;
    if (cos_theta < 0.015629) L = 0;
    
    if (isnan(L)) {
    	cout << cos_theta << " " << mexmcc.mccangl[i_mu] << " " << mexmcc.mccangl[n_mu] << endl;
    	cout << I_int[0] << " " << J[0] << " " << J[1] << " " << J[2] << " " << J[3] << endl;
    }

        return L;
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
       double e_l;
       double counts_l;

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



