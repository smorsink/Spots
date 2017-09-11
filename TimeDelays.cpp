#include "matpack.h"
#include <exception>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unistd.h>
#include "TimeDelays.h"
#include "Exception.h"
#include "Units.h"
#include "Struct.h"
#include "time.h"
#include "interp.h"
#include <stdio.h>
using namespace std;


/**************************************************************************************/
/* TimeDelays:                                                                      */
/*              computes the flux of each light curve                                 */
/*																					  */
/* pass: angles = all the angles necessary to compute the flux correctly;             */
/*                computed in the routine/method/function above [radians or unitless] */
/**************************************************************************************/
class LightCurve TimeDelays( class LightCurve* angles ) {
	
  class LightCurve curve;
  curve = (*angles);

  int numbins = curve.numbins;
  int numbands = curve.numbands;

  double tvec[4], fvec[4], err;

  std::ofstream ttt;

    /***********************************************************/
    /* DEALING WITH THE TIME DELAYS, REBINNING THE LIGHT CURVE */
    /* This is where the jumpy problems tend to be.            */
    /***********************************************************/

    curve.flags.ignore_time_delays = false;
    		
    if ( !curve.flags.ignore_time_delays ) { // if we are not ignoring the time delays        
      int k(0), j(0); // index placeholders; approximately, k is i+1 and j is i-1
      // but time delays mean that j isn't always i-1
      // used in the linear interpolation
        
      std::cout << "Taking care of time delays!" << std::endl;


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
       
      for (unsigned int p(0); p < NCURVES; p++) {
      //for (unsigned int p(0); p < 1; p++) {

	//std::cout << "TimeDelays: p=" << p << std::endl;

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

	    // std::cout << "TimeDelays: i=" << i <<std::endl;
            
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

		int npt=3;
		int start=0;
		int index=0;
		j -= 1;

		for (int m=1; m<=npt; m++){
		  start=0;
		  
		  if (j+m+start >= numbins)
		    start -= numbins;
		  if (j+m+start < 0)
		    start += numbins;

		  tvec[m] = curve.t_o[j+m+start];
		  if (start < 0) tvec[m] += 1.0;
		  
		  if (curve.t[i] == 0.0){
		    tvec[m] -= 1.0;
		  }

		  fvec[m] = curve.f[p][j+m+start];
		  //std::cout << " START = " << start << " " 
		  //	    << "m=" << m << " t=" << tvec[m] << " f=" << fvec[m] << std::endl;

		}

		j+=1;

		//newflux.at(i) = curve.f[p][j] + (curve.f[p][k]-curve.f[p][j])/(t2-t1) * (curve.t[i]-t1); // linear interpolation!
		
		//double oldflux = newflux.at(i);

		//std::cout << "Old Interpolation time = " << curve.t[i] << "flux = " << newflux.at(i) << std::endl; 

		//std::cout << "do new interpolation!" << std::endl;
		newflux.at(i) = polint(tvec,fvec,npt,curve.t[i],&err);
		//std::cout << "New Interpolation time = " << curve.t[i] << "flux = " << newflux.at(i) << std::endl; 

		/*	if (p==0)
		  std::cout << "i = " << i  
			    << " time = " << curve.t[i]
		    << " ********* Flux diff = " << (oldflux - newflux.at(i)) 
			  << " oldflux = " << oldflux 
			    << " newflux = " << newflux.at(i)
			    << " percent diff = " << (oldflux-newflux.at(i))/oldflux * 100 
			  << std::endl;
		*/
	

		if (curve.eclipse){
		  if (curve.t[i] < te1 && curve.t[i] > te1 - 1.0/numbins){
		    /* if (p==0){
		      std::cout << "time just before eclipse at te1=" << te1 <<"! i=" << i 
				<< " t = " << curve.t[i] << std::endl;
		      std::cout << "before correction, flux= " << newflux.at(i) << std::endl;
		      }*/
		    newflux.at(i) = slope1 * (curve.t[i] - te1);
		    /*if (p==0) 
		      std::cout << "after correction, flux= " << newflux.at(i) << std::endl;*/
		  }
		  if (curve.t[i] < te1 + 1.0/numbins && curve.t[i] > te1)
		    newflux.at(i) = 0.0;

		  if (curve.t[i] > te2 && curve.t[i] < te2 + 1.0/numbins){
		    /* if (p==0)
		       std::cout << "eclipse ends at te2=" << te2 <<"! i=" << i << " t = " << curve.t[i] << std::endl;*/
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
	  /*if (p==0) {

	    std::cout << "Opening newtime.txt for export " << std::endl;

	    ttt.open("newtime512-3.txt", std::ios_base::trunc);
	    }*/
	  for ( unsigned int i(0); i < numbins; i++ ) {
	    /*if ( p==0 ) 
	      ttt << curve.t_o[i] << " " 
	    	  << curve.f[p][i] << " " << curve.t[i] << " " << newflux.at(i) 
	    	  << " " << i
	    	  << std::endl;*/
	    
	    //            curve.f[p][i] = totflux.at(i);
            curve.f[p][i] = newflux.at(i);
	  }
	  // only plotting versus evenly spaced time -- not using t_o
	
      }// end for-p-loop
    } // end time delay section


    
    return curve;

} // end ComputeCurve

