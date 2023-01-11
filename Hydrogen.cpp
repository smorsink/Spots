/***************************************************************************************/
/*                                     Hydrogen.cpp

    This holds the routines used to load and access the Hydrogen Atmosphere Model
    computed by Wynn Ho. 

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
#include "Hydrogen.h"
#include "Exception.h"
#include "Units.h"
#include "Struct.h"
#include "time.h"
#include "interp.h"
#include "nrutil.h"
#include <stdio.h>
using namespace std;



// Computes the correct theta index for NSX
int th_index_nsx(double cos_theta, class LightCurve* mexmcc){

  int n_mu;

  //std::cout << "th_index_nsx: cos_theta = " << cos_theta << std::endl;  
  //Find proper mu choice
    n_mu = 1;
    while (cos_theta < mexmcc->mccangl[n_mu] && n_mu < mexmcc->Nmu){
    	n_mu += 1;
    }
    /*std::cout << "th_index_nsx: cos_theta = " << cos_theta << " nmu = " << n_mu
	      << " mu[i] = " << mexmcc->mccangl[n_mu]
	      << " mu[i-1] = " << mexmcc->mccangl[n_mu-1]
	      << " mu[i-2] = " << mexmcc->mccangl[n_mu-2]
	      << std::endl;*/
    return (n_mu-1);
}



// Read in Wynn Ho's Hydrogen Table
void ReadNSXHnew(class LightCurve* curve){

	double loget,mu,logt,logg,logi;
	
	int NlogTeff, Nlogg, NlogE, Nmu, Npts;

	
	std::cout << "Reading in NSX version 200804" << std::endl;
	
	//	std::ifstream Hspecttable; 
	//	Hspecttable.open( "atmosphere/nsx_H_v171019.out" );  // current version of NSX

      FILE *Hspecttable;
      Hspecttable=fopen("../atmosphere/nsx_H_v200804.out","r");
      std::cout << "Done opening file" << std::endl;
	  
	NlogTeff = 35;
	curve->mclogTeff = dvector(0,NlogTeff);
	for (int i=0;i<NlogTeff;i++){
	  curve->mclogTeff[i] = 5.10 + i*0.05;	  
	}

	Nlogg = 14; // updated 20221116
	curve->mclogg = dvector(0,Nlogg);
	for (int i=0;i<Nlogg;i++){
	  curve->mclogg[i] = 13.7 + 0.1*i;
	}

	NlogE = 166;
	curve->mcloget = dvector(0,NlogE);
	for (int i=0;i<NlogE;i++){
	  curve->mcloget[i] = -1.30 + i*0.02;
	}

	Nmu = 67;
	Npts =  (NlogTeff*Nlogg*NlogE);

	curve->NlogTeff = NlogTeff;
	curve->Nlogg = Nlogg;
	curve->NlogE = NlogE;
	curve->Nmu = Nmu;
	curve->Npts = Npts;

	std::cout << "Npts = " << Npts << std::endl;

	curve->mccangl = dvector(0,Nmu);
	curve->mccinte = dvector(0,Npts*Nmu);

	std::cout << "Finished allocating memory " << std::endl;
	
    	for (int i = 0; i < Npts*Nmu; i++){
	  
	  fscanf(Hspecttable,"%lf %lf %lf %lf %lf\n",
		 &loget, &mu,
		 &logi, &logt, &logg);	
 
	  curve->mccinte[i] = logi;
	  //std::cout << " logi = " << curve.mccinte[i] << std::endl;
		    
	  if ( i < Nmu ) {
	    curve->mccangl[i] = mu;
	    //std::cout << " mu = " << curve.mccangl[i] << std::endl;
	  }
	  
	  /* if (logt == 5.1 && mu == 0.5 && loget == 1.0 && logg == 15)
	    std::cout
	    << " i = " << i
	    << " i%(Nmu) = " << i%(Nmu)
	    << " logT = " << logt << " = " << curve.mclogTeff[i/(Nlogg*NlogE*Nmu)]
	    << " logg = " << logg << " = " << curve.mclogg[i/(NlogE*Nmu*NlogTeff)]
	    << " logE = " << loget << " = " << curve.mcloget[i/(Nmu*Nlogg*NlogTeff)]
	    << " cos(theta) = " << curve.mccangl[i%(Nmu)]
		    << " I = " << curve.mccinte[i]
		    << std::endl;*/
	  
	}
       
	fclose(Hspecttable);
    	std::cout << "finished reading Wynn Ho's NSX-H" << std::endl;



  

}







// NSX - Hydrogen Atmosphere Computed by Wynn Ho
// Calculate the final interpolated intensity
// This "new" version takes into account that the energy is really the ratio: E/kT
double NSXHnew(double E, double cos_theta, int theta_index, double T, double lgrav,  class LightCurve* curve){

  //std::cout << "Welcome to NSX! E=" << E << std::endl;

  class LightCurve mexmcc;
  mexmcc = (*curve);

  double lt, ener_index;
  double logET;
	
  double I_int[5], J[5], K[5], L(0.0);
  int i_f, i_lt,  first_inte;  
  int ii_mu = theta_index;

  if ( ii_mu > 62)
    ii_mu = 62;

  // CALCULATE LOG(T) AND CORRECT INDEX
  if (curve->flags.kelvin){
    lt = log10(T);
    logET = log10(E/(T*Units::K_BOLTZ/Units::EV*1e-3)); // Convert T to keV
    //std::cout << "Temperature in Kelvin! log(T) = " << lt << std::endl;
  }
  else{
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));
      // CALCULATE LOG(E/T) 
    logET = log10(E/T);

  }

   
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
    
    //    int npt(4), tpt(4), gpt(4), mpt(4);
    int npt(2), tpt(2), gpt(2), mpt(2);
    
    if ( npt > 4)
      L = 0.0;
    else{



      
    //Find proper freqency choice
    ener_index = (logET + 1.30)/0.02;
    i_f = (int) ener_index;
    //std::cout << "ener_index = " << ener_index << " i_f = " << i_f << std::endl;
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

	    /* if ( cos_theta < 0.2 && r==0 && q==0 && j==0)
	      std::cout << "j = " << j
			<< "first_inte = " << first_inte
			<< " tvec = " << tvec[r+1]
			<< " gvec = " << gvec[q+1]
			<< " costheta = " << mexmcc.mccangl[ii_mu+k]
		      << " evec = " << mexmcc.mcloget[i_f+j] 
		      << " ivec = " << ivec[r][q][k][j+1] << std::endl;
	    */
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










