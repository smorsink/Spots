/***************************************************************************************/
/*                               OblDeflectionTOA.cpp

    This code calculates the relationships between impact parameter (b), deflection angle,
    and ray time of arrival (TOA), given for an oblate NS model.
    
    Based on code written by Coire Cadeau and modified by Sharon Morsink and 
    Abigail Stevens.
    
    (C) Coire Cadeau, 2007; Source (C) Coire Cadeau 2007, all rights reserved.
    Permission is granted for private use only, and not distribution, either verbatim or
    of derivative works, in whole or in part.
    This code is not thoroughly tested or guaranteed for any particular use.
*/
/***************************************************************************************/

#include "matpack.h"

#include <exception>
#include <cmath>
#include <iostream>
#include "OblDeflectionTOA.h"
#include "OblModelBase.h"
#include "Exception.h"
#include "Units.h"
#include "Struct.h" //Jan 21 (year? prior to 2012)
// Globals are bad, but there's no easy way to get around it
// for this particular case (need to pass a member function
// to a Matpack routine which does not have a signature to accomodate
// passing in the object pointer).
#define EPSILON 1E-15

class LightCurve curve; //Jan 21 (year? prior to 2012)

const OblDeflectionTOA* OblDeflectionTOA_object;
double OblDeflectionTOA_b_value;
double OblDeflectionTOA_b_over_r;
double OblDeflectionTOA_costheta_value;
double OblDeflectionTOA_psi_value;
double OblDeflectionTOA_b_max_value;
double OblDeflectionTOA_psi_max_value;
//double OblDeflectionTOA_b_guess;
//double OblDeflectionTOA_psi_guess;
double OblDeflectionTOA_rspot;

/**********************************************************/
/* OblDeflectionTOA_psi_integrand_wrapper:                */
/*                                                        */
/* Passed r (NS radius, unitless), *prob (problem toggle) */
/* Returns integrand if it is not nan.                    */
/**********************************************************/


double OblDeflectionTOA_psi_integrand_wrapper_u ( double u,  bool *prob ) {
  double integrand( OblDeflectionTOA_object->psi_integrand_u( OblDeflectionTOA_b_over_r, u)); 
  //OblDeflectionTOA_object is looking up psi_integrand by feeding it OblDeflectionTOA_b_over_r, which is defined in psi_outgoing
   if( std::isnan(integrand) ) {
    std::cout << "integrand is nan!" << std::endl;
    std::cout << "b/r=" << OblDeflectionTOA_b_over_r << std::endl;
    std::cout << "u = " << u  << std::endl;
    std::cout << "(u*b/r)^2)=" << pow(u*OblDeflectionTOA_b_over_r,2);
    std::cerr << "ERROR in OblDeflectionTOA_psi_integrand_wrapper_u(): returned NaN." << std::endl;
    *prob = true;
    integrand = -7888.0;
    return integrand;
    }
  return integrand;
}

double OblDeflectionTOA_psi_ingoing_integrand_wrapper_u ( double u,  bool *prob ) {
  double integrand( OblDeflectionTOA_object->psi_ingoing_integrand_u( OblDeflectionTOA_b_over_r, u)); 
 
  //OblDeflectionTOA_object is looking up psi_integrand by feeding it OblDeflectionTOA_b_over_r, which is defined in psi_outgoing
  if( std::isnan(integrand) ) {
    std::cout << "integrand is nan!" << std::endl;
    std::cout << "b/r=" << OblDeflectionTOA_b_over_r << std::endl;
    std::cout << "u = " << u  << std::endl;
    std::cout << "(u*b/r)^2)=" << pow(u*OblDeflectionTOA_b_over_r,2);
    std::cerr << "ERROR in OblDeflectionTOA_psi_ingoing_integrand_wrapper_u(): returned NaN." << std::endl;
    *prob = true;
    integrand = -7888.0;
    return integrand;
  }
  return integrand;
}



/************************************************************/
/* OblDeflectionTOA_dpsi_db_integrand_wrapper:              */
/*                                                          */
/* Passed r (NS radius, unitless), *prob (problem toggle)   */
/* Returns integrand if it is not nan.                      */
/************************************************************/
double OblDeflectionTOA_dpsi_db_integrand_wrapper_u ( double u, bool *prob ) {
	double integrand( OblDeflectionTOA_object->dpsi_db_integrand_u( OblDeflectionTOA_b_over_r, u ) );
  
  	if ( std::isnan(integrand) ) {
    	std::cout << "dpsi_db_u integrand is nan!" << std::endl;
    	std::cout << "u (km) = " << Units::nounits_to_cgs(u, Units::LENGTH)/1.0e5 << std::endl;
    	std::cerr << "ERROR in OblDeflectionTOA_dpsi_db_integrand_wrapper(): returned NaN." << std::endl;
    	*prob = true;
    	integrand = -7888.0;
    	return integrand;
  	}
  
  	return integrand;
}

/************************************************************/
/* OblDeflectionTOA_dpsi_db_integrand_wrapper:              */
/*                                                          */
/* Passed r (NS radius, unitless), *prob (problem toggle)   */
/* Returns integrand if it is not nan.                      */
/************************************************************/
double OblDeflectionTOA_dpsi_db_ingoing_integrand_wrapper_u ( double u, bool *prob ) {
	double integrand( OblDeflectionTOA_object->dpsi_db_ingoing_integrand_u( OblDeflectionTOA_b_over_r, u ) );
  
  	if ( std::isnan(integrand) ) {
    	std::cout << "dpsi_db_u integrand is nan!" << std::endl;
    	std::cout << "u (km) = " << Units::nounits_to_cgs(u, Units::LENGTH)/1.0e5 << std::endl;
    	std::cerr << "ERROR in OblDeflectionTOA_dpsi_db_integrand_wrapper(): returned NaN." << std::endl;
    	*prob = true;
    	integrand = -7888.0;
    	return integrand;
  	}
  
  	return integrand;
}


/***********************************************************/
/* OblDeflectionTOA_toa_integrand_minus_b0_wrapper:        */
/*                                                         */
/* Passed r (NS radius, unitless), *prob (problem toggle)  */
/* Returns integrand if it is not nan.                     */
/***********************************************************/
double OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u ( double u, bool *prob ) {
  	double integrand( OblDeflectionTOA_object->toa_integrand_minus_b0_u( OblDeflectionTOA_b_over_r, u ) );
  	if ( std::isnan(integrand) ) {
    	std::cout << "toa_minus_b0_u integrand is nan!" << std::endl;
    	std::cout << "r (km) = " << Units::nounits_to_cgs(u, Units::LENGTH)/1.0e5 << std::endl;
    	std::cerr << "ERROR in OblDeflectionTOA_toa_integrand_minus_b0_wrapper(): returned NaN." << std::endl;
    	*prob = true;
    	integrand = -7888.0;
    	return integrand;
  	}
  	return integrand;
}

double OblDeflectionTOA_toa_ingoing_integrand_minus_b0_wrapper_u ( double u, bool *prob ) {
  	double integrand( OblDeflectionTOA_object->toa_ingoing_integrand_minus_b0_u( OblDeflectionTOA_b_over_r, u ) );
  	if ( std::isnan(integrand) ) {
	  std::cout << "toa_minus_b0_u ingoing integrand is nan!" << std::endl;
	  std::cout << "r (km) = " << Units::nounits_to_cgs(u, Units::LENGTH)/1.0e5 << std::endl;
	  std::cerr << "ERROR in OblDeflectionTOA_toa_integrand_minus_b0_wrapper(): returned NaN." << std::endl;
    	*prob = true;
    	integrand = -7888.0;
    	return integrand;
  	}
  	return integrand;
}

/*****************************************************/
/*****************************************************/
double OblDeflectionTOA_rcrit_zero_func_wrapper ( double rc ) {
  	return double ( OblDeflectionTOA_object->rcrit_zero_func( rc, 
  	                	OblDeflectionTOA_b_value ) );
}

/***********************************************************/
/* OblDeflectionTOA_b_from_psi_ingoing_zero_func_wrapper:  */
/*                                                         */
/* Passed b (photon impact parameter)                      */
/* Returns                      */
/***********************************************************/
double OblDeflectionTOA_b_from_psi_ingoing_zero_func_wrapper ( double b_R ) {
  	return double ( OblDeflectionTOA_object->b_from_psi_ingoing_zero_func( b_R, 
									       OblDeflectionTOA_b_max_value/OblDeflectionTOA_rspot,
									       OblDeflectionTOA_psi_max_value,
						OblDeflectionTOA_rspot,
						OblDeflectionTOA_psi_value ) );
}


// End of global pollution.


// Defining constants
//const double OblDeflectionTOA::INTEGRAL_EPS = 1.0e-7;
const double OblDeflectionTOA::FINDZERO_EPS = 1.0e-6;

const double OblDeflectionTOA::RFINAL_MASS_MULTIPLE = 1.0e7;
//const double OblDeflectionTOA::DIVERGENCE_GUARD = 2.0e-2; // set to 0 to turn off
const double OblDeflectionTOA::DIVERGENCE_GUARD = 0.0;
const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_N_MAX = 100000;
const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_N_MAX_1 = 10000;
const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_N = 2000;
const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_N_1 = 2000;
const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_INGOING_N = 4000;
//const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_INGOING_N = 100000;
//const double OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_POWER = 4.0; 
const double OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_POWER = 4.0;

OblDeflectionTOA::OblDeflectionTOA ( OblModelBase* modptr, const double& mass_nounits, const double& mass_over_r_nounits, const double& radius_nounits ) 
  : r_final( RFINAL_MASS_MULTIPLE * mass_nounits ) {
  this->modptr = modptr;
  mass = mass_nounits;
  mass_over_r = mass_over_r_nounits;
  rspot = radius_nounits;
}

/************************************************************/
/* OblDeflectionTOA::bmax_outgoing                          */
/*														    */
/* Computes the maximal value of b for outgoing light rays  */
/* Passed rspot (radius of NS at spot, unitless)            */
/* Returns maximal value of b (in funny unitless units!)    */
/************************************************************/
double OblDeflectionTOA::bmax_outgoing ( const double& rspot ) const {

  double bmax_outgoing ( rspot / sqrt( 1.0 - 2.0 * get_mass_over_r() ) );
  	return bmax_outgoing;
}

/************************************************************/
/* OblDeflectionTOA::b_R_max_outgoing                          */
/* Uses value of M/R where R=rspot						       			      	    */
/* Computes the maximal value of b/R for outgoing light rays  */
/************************************************************/
double OblDeflectionTOA::b_R_max_outgoing ( const double& mass_over_r ) const {

  double b_R_max_outgoing ( 1.0 / sqrt( 1.0 - 2.0 * mass_over_r ) );
  	return b_R_max_outgoing;
}



/*******************************************************************************/
/* OblDeflectionTOA::bmin_ingoing                                              */
/*														                       */
/* Computes value of b corresponding to the most ingoing light ray allowed,    */
/*   i.e. the ray that just grazes the surface and moves toward axis with b=0  */
/* Passed cos_theta (cosine of theta (what is theta here?))                    */
/* Returns minimum value of b allowed                                          */
/*******************************************************************************/
double OblDeflectionTOA::bmin_ingoing ( const double& rspot, const double& cos_theta ) const { 

  	// Need to return the value of b corresponding to dR/dtheta of the surface.

  	double drdth ( modptr->Dtheta_R( cos_theta ) );
	double b ( rspot / sqrt( (1.0 - 2.0 * get_mass_over_r()) + pow( drdth / rspot , 2.0 ) ) );
 
  	// NaN check:
  	if ( std::isnan(b) ) {
    	std::cerr << "ERROR in OblDeflectionTOA::bmin_ingoing(): returned NaN." << std::endl;
    	b = -7888.0;
    	return b;
  	}

  	return b;    // returning b!!!
}

/*******************************************************************************/
/* OblDeflectionTOA::b_R_min_ingoing                                              */
/*														                       */
/* Computes value of b/r corresponding to the most ingoing light ray allowed,    */
/*   i.e. the ray that just grazes the surface and moves toward axis with b=0  */
/* Passed cos_theta (cosine of theta (what is theta here?))                    */
/* Returns minimum value of b allowed                                          */
/*******************************************************************************/
double OblDeflectionTOA::b_R_min_ingoing ( const double& rspot, const double& cos_theta ) const { 

  	// Need to return the value of b/R corresponding to dR/dtheta of the surface.

  	double drdth ( modptr->Dtheta_R( cos_theta ) );
	double b_R ( 1.0 / sqrt( (1.0 - 2.0 * get_mass_over_r()) + pow( drdth / rspot , 2.0 ) ) );
 
  	// NaN check:
  	if ( std::isnan(b_R) ) {
    	std::cerr << "ERROR in OblDeflectionTOA::bmin_ingoing(): returned NaN." << std::endl;
    	b_R = -7888.0;
    	return b_R;
  	}

  	return b_R;    // returning b!!!
}




/*************************************************************/
/* OblDeflectionTOA::ingoing_allowed                         */
/*														     */
/* Determines if an ingoing light ray is allowed             */
/*   for locations where dR/dtheta \neq 0                    */
/* Passed cos_theta (cosine of theta (what is theta here?))  */
/* Returns boolean (true if allowed)                         */
/*************************************************************/
bool OblDeflectionTOA::ingoing_allowed ( const double& cos_theta ) {
  	// ingoing rays can be a consideration for locations where dR/dtheta \neq 0
  	//costheta_check( cos_theta );
  	return bool ( modptr->Dtheta_R( cos_theta ) != 0.0 );
}

/************************************************************************/
/* OblDeflectionTOA::rcrit                                              */
/*														                */
/* Determines if an ingoing light ray is allowed                        */
/*   for locations where dR/dtheta \neq 0                               */
/* Passed b (photon impact parameter), cos_theta (cos of spot angle?),  */
/*        prob (problem toggle)                                         */
/* Returns                                    */
/************************************************************************/
double OblDeflectionTOA::rcrit ( const double& b, const double& rspot, bool *prob ) const {
  	double candidate;
       
  	//double r ( modptr->R_at_costheta( cos_theta ) );
 
  	if ( b == bmax_outgoing(rspot) ) {
    	return double(rspot);
  	}

  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_value = b;

	/*	candidate = MATPACK::FindZero(modptr->R_at_costheta(1.0),
				modptr->R_at_costheta(0.0),
				OblDeflectionTOA_rcrit_zero_func_wrapper);*/ // rcrit_guess

	candidate = MATPACK::FindZero(0.6*rspot,
				      rspot,
				      OblDeflectionTOA_rcrit_zero_func_wrapper); // rcrit_guess
 
 
  	return candidate;
}



/*****************************************************/
double OblDeflectionTOA::psi_outgoing_u ( const double& b_R, 
					  const double& b_R_max, const double& psi_max, 
					  bool *prob ) const{

  double dummy;

  if ( b_R > b_R_max || b_R < 0.0 ) {
    std::cerr << "ERROR in OblDeflectionTOA::psi_outgoing_u(): b out-of-range." << std::endl;
    std::cerr << "M/R = " << get_mass_over_r() 
	      << " b/R = " << b_R 
	      << " b_R_max = " << b_R_max
	      << std::endl;



    *prob = true;
    	dummy = -7888.0;
    	return dummy;
  }

  if ( fabs( b_R - b_R_max ) < 1e-8 ) { // essentially the same as b=b_max
    return psi_max;
  }
  else {  // set globals
    OblDeflectionTOA_object = this;
    OblDeflectionTOA_b_over_r = b_R;     // this will be fed to another function

    double psi(0.0);
    double split(0.90);  

    if ( b_R != 0.0 ) {

      if ( b_R > 0.9995*b_R_max){
	psi = Integration( 0.0, split, OblDeflectionTOA_psi_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
	psi += Integration( split, 1.0, OblDeflectionTOA_psi_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1*10 );
      }
      else{

	if ( b_R > 0.995*b_R_max){	
	  psi = Integration( 0.0, split, OblDeflectionTOA_psi_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
	  psi += Integration( split, 1.0, OblDeflectionTOA_psi_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1 );
	}
	else
	  psi = Integration( 0.0, 1.0, OblDeflectionTOA_psi_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
	// default case 
      }
    }				
    return psi;
  }
}





double OblDeflectionTOA::psi_max_outgoing_u ( const double& b_R, bool *prob ) const{
  	double dummy;
	double mass_over_r( get_mass_over_r());

  	if ( (b_R - b_R_max_outgoing(mass_over_r) > 1e-3) || b_R < 0.0 ) {
    	std::cerr << "ERROR in OblDeflectionTOA::psi_max_outgoing(): b out-of-range." << std::endl;

	std::cerr << "M/R = " << get_mass_over_r() 
	      << " b/R = " << b_R 
	      << " b_R_max = " <<  b_R_max_outgoing(mass_over_r)	      << std::endl;

    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}

  	// set globals
  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_over_r = b_R;

	double psi, split(0.9);  	

	psi = Integration( 0.0, split, OblDeflectionTOA_psi_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
	psi += Integration( split, 1.0, OblDeflectionTOA_psi_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1*10 );
					
  	return psi;
}

/*****************************************************/
/*****************************************************/
double OblDeflectionTOA::psi_ingoing ( const double& b_R, const double& b_R_max, const double& psimax, 
				       const double& rspot, bool *prob ) const {

  double psi, psi_in, rcrit;
	
  	// set globals
  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_value = b_R * rspot;
	OblDeflectionTOA_b_max_value = b_R_max * rspot;
	OblDeflectionTOA_psi_max_value = psimax;

	if ( fabs(b_R-b_R_max) < 1e-6){
	  psi_in = 0.0;
	}
	else{

	  rcrit = this->rcrit( b_R*rspot, rspot, &curve.problem );	 
	  OblDeflectionTOA_b_over_r = b_R*rspot/rcrit;
	  psi_in = 2.0*Integration( rcrit/rspot, 1.0, OblDeflectionTOA_psi_ingoing_integrand_wrapper_u,10000 );
	}
	
  	
   	OblDeflectionTOA_b_over_r = b_R;
	
	if ( b_R == b_R_max)
	  psi = psi_in +  OblDeflectionTOA_psi_max_value;
	  else
	  psi = psi_in + this->psi_outgoing_u( b_R, b_R_max,
						      OblDeflectionTOA_psi_max_value , &curve.problem );

	//	std::cout << "Total psi = " << psi << std::endl;

  	return psi;
}

/*****************************************************/
// Major Change to b_from_psi
// a guess for b is given to this routine
// if the photon is purely outgoing nothing is changed
// if the photon is initially ingoing a calculation is done
/*****************************************************/
bool OblDeflectionTOA::b_from_psi ( const double& psi, const double& rspot, const double& cos_theta, double& b_R,
				    int& rdot, const double& bmax_out, 
				    const double& psi_out_max,
				    const double& bmin, const double& psimin,
				    const double& b_guess, 
				    bool *prob ) {

  // return bool indicating whether a solution was found.
  // store result in b_R
  // store -1 in rdot if its an ingoing solution, +1 otherwise.

  	double bcand;    //
  	double dummyb;   //
  	double dummypsi; //

	OblDeflectionTOA_b_max_value = bmax_out;
	OblDeflectionTOA_psi_max_value = psi_out_max;

  	if ( std::isinf(psi_out_max) ) {
	  std::cerr << "ERROR in OblDeflectionTOA::b_from_psi(): psi_out_max = infinity" << std::endl;
	  *prob = true;
	  dummyb = -7888.0;
	  return dummyb;
  	}

  	if ( psi < 0 ) {
	  std::cerr << "ERROR in OblDeflectionTOA::b_from_psi: need psi >= 0" << std::endl;
	  *prob = true;
	  dummypsi = -7888.0;
	  return dummypsi;
  	}
  	else if ( psi == 0.0 ) {
	  b_R = 0.0;
	  rdot = 1;
	  return true;
  	}
  	else if ( psi == psi_out_max ) {
	  b_R = bmax_out/rspot;
	  rdot = 1;
	  return true;
  	}
  	else if ( psi < psi_out_max ) { // begin normal outgoing case, psi < psi_out_max

	  bcand = b_guess;

	  if ( fabs(bcand) <= EPSILON || fabs(bmax_out - bcand) <= EPSILON ) { // this indicates no soln
	    std::cerr << "ERROR in OblDeflectionTOA::b_from_psi(): outgoing returned no solution?" << std::endl;
	    *prob = true;
	    b_R = -7888.0;
	    return b_R;
	  }
	  else {
	    b_R = bcand/rspot;
	    rdot = 1;
	    return true;
	  }
  	} // end normal outgoing case, psi < psi_out_max
  	
  	else { // psi > psi_out_max 
	  // Test to see if ingoing photons are allowed
	  bool ingoing_allowed( this->ingoing_allowed(cos_theta) );
	 

	  if ( ingoing_allowed ) {

	    //std::cout << "b_from_psi: Ingoing!!!" << std::endl;


	    if ( psi > psimin ) {
	      //std::cout << "b_from_psi: psi > psimin" << std::endl;
	      return false;
	    }
	    else if ( psi == psimin ) {
	      b_R = bmin/rspot;
	      rdot = -1;
	      return true;
	    }
	    else { // psi_out_max < psi < psi_in_max
	      OblDeflectionTOA_object = this;
	      OblDeflectionTOA_costheta_value = cos_theta;
	      OblDeflectionTOA_psi_value = psi;
	      OblDeflectionTOA_rspot = rspot;

	      /*
	      bcand = MATPACK::FindZero(bmin, bmax_out,
					OblDeflectionTOA_b_from_psi_ingoing_zero_func_wrapper,
					OblDeflectionTOA::FINDZERO_EPS);
	      */

	      double bcand_R = MATPACK::FindZero(bmin/rspot, bmax_out/rspot,
					OblDeflectionTOA_b_from_psi_ingoing_zero_func_wrapper,
					OblDeflectionTOA::FINDZERO_EPS);


	      /*std::cout << "b_from_psi: bmin_R = " << bmin/rspot
			<< " bmax_R = " << bmax_out/rspot
			<< " bcand_R = " << bcand_R 
			<< std::endl;*/

 
	      if( fabs(bmin/rspot - bcand_R) <= EPSILON || fabs(bmax_out/rspot - bcand_R) <= EPSILON ) { // this indicates no soln
		std::cerr << "ERROR in OblDeflectionTOA::b_from_psi(): ingoing returned no solution?" << std::endl;
		*prob = true;
		b_R = -7888.0;
		return b_R;
	      }
	      else {
		b_R = bcand_R;
		rdot = -1;
		return true;
	      }
	    }
	  } // end ingoing photons are allowed
	  else { // ingoing photons not allowed.
	    //std::cerr << "Eclipse!!" << std::endl;
	    return false;
	  }
  	}
  	std::cerr << "ERROR in OblDeflectionTOA::b_from_psi(): reached end of function?" << std::endl;
  	*prob = true;
  	return false;
}


double OblDeflectionTOA::dpsi_db_outgoing_u( const double& b_R, bool *prob ) const {
  
  // dpsi_db returns a dimensionless quantity. You'll need to divide by radius to make it have dimensions of 1/length

  double mass_over_r( get_mass_over_r());
  double dummy; //

  double b_R_max = b_R_max_outgoing(mass_over_r);

  	if ( (b_R > b_R_max || b_R < 0.0) ) { 
    	std::cerr << "ERROR inOblDeflectionTOA::dpsi_db_outgoing(): b out-of-range." << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}
  
  	// set globals
  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_over_r = b_R;

	double dpsidb(0.0);
	double split(0.5);

	if (b_R < 0.95*b_R_max)
	  dpsidb = Integration( 0.0, 1.0, OblDeflectionTOA_dpsi_db_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N*10 );
	else
	  if ( b_R < 0.99995*b_R_max){
	      	
	    //dpsidb = Integration( 0.0, split, OblDeflectionTOA_dpsi_db_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N*10 );
	    dpsidb = Integration( 0.0, 1.0, OblDeflectionTOA_dpsi_db_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1*100 );
	  }
	  else{
	    if ( b_R < 0.999995*b_R_max){
	      	
	    dpsidb = Integration( 0.0, split, OblDeflectionTOA_dpsi_db_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
	    dpsidb += Integration( split, 1.0, OblDeflectionTOA_dpsi_db_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1*1000 );
	    }
	    else{
	    dpsidb = Integration( 0.0, split, OblDeflectionTOA_dpsi_db_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N );
	    dpsidb += Integration( split, 1.0, OblDeflectionTOA_dpsi_db_integrand_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1*1000 );
	    }

	  }
	  
	    
	//dpsidb *= 1.0/rspot;
  	return dpsidb;
}


double OblDeflectionTOA::dpsi_db_ingoing_u( const double& b, const double& rspot, const double& cos_theta, bool *prob ) {

  double dpsidb_in;
  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_value = b;

  
	double rcrit = this->rcrit( b, rspot, &curve.problem );
 

	//	double bmax(rspot/sqrt(1.0-2.0*get_mass_over_r()));

	OblDeflectionTOA_b_over_r = b/rcrit;

	dpsidb_in = 2.0/rcrit * Integration( rcrit/rspot, 1.0, OblDeflectionTOA_dpsi_db_ingoing_integrand_wrapper_u,100000 );
	//std::cout << "Approx 10000: toa_in = " << toa_in << std::endl;

	OblDeflectionTOA_b_over_r = b/rspot;
  	return double( dpsidb_in + this->dpsi_db_outgoing_u( b/rspot, &curve.problem ) );

}


double OblDeflectionTOA::toa_outgoing_u ( const double& b_R, bool *prob ) {
  // Returns a dimensionless toa!!!!
  // Must multiply by radius later on to proper dimensions!!!

  double mass_over_r(get_mass_over_r());

  double b_R_max=b_R_max_outgoing(mass_over_r);
  //double b_R=b/rspot;

  OblDeflectionTOA_object = this;
  OblDeflectionTOA_b_over_r = b_R;

  	// Note that for a b=0 ray, Delta T = int(1/(1-2*M/r), r),
  	// which is r + 2 M ln (r - 2M)

	double toa(0.0);
	double split(0.90);

	if (b_R < 0.95*b_R_max)
	  toa = Integration( 0.0, 1.0, OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u,TRAPEZOIDAL_INTEGRAL_N/10 );
	else
	  if ( b_R < 0.9999*b_R_max){
	      	
	    toa = Integration( 0.0, split, OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u,TRAPEZOIDAL_INTEGRAL_N/10 );
	    toa += Integration( split, 1.0, OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1/10);

	  }
	  else{
	    toa = Integration( 0.0, split, OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u,TRAPEZOIDAL_INTEGRAL_N/10 );
	    toa += Integration( split, 1.0, OblDeflectionTOA_toa_integrand_minus_b0_wrapper_u,TRAPEZOIDAL_INTEGRAL_N_1*100 );
	  }
	 
	//toa *= rspot;
	   
	//	double toa_b0_eqsurf = (rspot - req) + 2.0 * get_mass() * ( log( rspot - 2.0 * get_mass() )
	//	       				 - log( req - 2.0 * get_mass() ) );

	//	return double( toa - toa_b0_eqsurf);

	return (toa);

}


/*****************************************************/
/*****************************************************/
double OblDeflectionTOA::toa_ingoing ( const double& b, const double& rspot, const double& cos_theta, bool *prob ) {

  double toa_in;
  	OblDeflectionTOA_object = this;
  	OblDeflectionTOA_b_value = b;

  
	double rcrit = this->rcrit( b, rspot, &curve.problem );
	/*	std::cout << "b = " << b
		  << " rspot = " << rspot
		  << " rcrit = " << rcrit << std::endl;*/

	//	double bmax(rspot/sqrt(1.0-2.0*get_mass_over_r()));

	OblDeflectionTOA_b_over_r = b/rcrit;

	// test=1000;
	int test(1000);

	toa_in = 2.0*Integration( rcrit/rspot, 1.0, OblDeflectionTOA_toa_ingoing_integrand_minus_b0_wrapper_u,test );
	//std::cout << "Approx 10000: toa_in = " << toa_in << std::endl;

	OblDeflectionTOA_b_over_r = b/rspot;
  	return double( toa_in + this->toa_outgoing_u( b/rspot, &curve.problem ) );

}


/*****************************************************/
// Integrand with respect to u=R/r. b_R = b/R
/*****************************************************/
double OblDeflectionTOA::psi_integrand_u ( const double& b_R, const double& u) const { 
  // double integrand( b / (r * r * sqrt( 1.0 - pow( b / r, 2.0 ) * (1.0 - 2.0 * get_mass_over_r() * get_rspot()/r ) )) );
  double integrand ( b_R / sqrt( 1.0 - pow(u*b_R,2)*(1.0-2.0*get_mass_over_r()*u) ));
  	return integrand;
}

/*****************************************************/
// Integrand with respect to u=rc/r. b_R = bin/rc
/*****************************************************/
double OblDeflectionTOA::psi_ingoing_integrand_u ( const double& b_R, const double& u) const { 
  // double integrand( b / (r * r * sqrt( 1.0 - pow( b / r, 2.0 ) * (1.0 - 2.0 * get_mass_over_r() * get_rspot()/r ) )) );

  double x( 1.0 - pow( b_R, -2.0)); // x = 2M/rc

  double integrand ( b_R / sqrt( 1.0 - pow(u*b_R,2)*(1.0-x*u) ));
  	return integrand;
}

/*****************************************************/
// gives the integrand from equation 38, MLCB
/*****************************************************/
double OblDeflectionTOA::toa_ingoing_integrand_minus_b0_u ( const double& b_R, const double& u ) const {   
  double x( 1.0 - pow( b_R, -2.0)); // x = 2M/rc
  double integrand(  get_rspot()/(u*u) * (1.0 / (1.0 - x*u)) 
		    * ( (1.0 / sqrt( 1.0 - pow( b_R * u , 2.0 ) * (1.0 - x*u) )) - 1.0) );
  	return integrand;
}

/*****************************************************/
// integrand for the derivative of equation 20, MLCB
/*****************************************************/
double OblDeflectionTOA::dpsi_db_integrand_u ( const double& b_R, const double& u ) const { 
 
  double integrand (   pow( 1.0 - pow( b_R*u , 2.0 ) * (1.0 - 2.0*get_mass_over_r() * u) , -1.5));
  // divide by rspot after!!!!
  	return integrand;
}

/*****************************************************/
// integrand for the derivative of equation 20, MLCB
/*****************************************************/
double OblDeflectionTOA::dpsi_db_ingoing_integrand_u ( const double& b_R, const double& u ) const { 
  //Divide this by r_crit afterwards!!!!!  

  double x( 1.0 - pow( b_R, -2.0)); // x = 2M/rc

  double integrand (   pow( 1.0 - pow( b_R*u , 2.0 ) * (1.0 - x*u) , -1.5));

  	return integrand;
}






/*****************************************************/
// gives the integrand from equation 38, MLCB
/*****************************************************/
double OblDeflectionTOA::toa_integrand_minus_b0_u ( const double& b_R, const double& u ) const {   
  double integrand(  1.0/(u*u) * (1.0 / (1.0 - 2.0 * get_mass_over_r()*u)) 
		    * ( (1.0 / sqrt( 1.0 - pow( b_R * u , 2.0 ) * (1.0 - 2.0 * get_mass_over_r()*u) )) - 1.0) );

  // Multiply final value by rspot!!!!!
  	return integrand;
}

/*****************************************************/
// used in oblate shape of star
/*****************************************************/
double OblDeflectionTOA::rcrit_zero_func ( const double& rc, const double& b ) const { 
  //	return double( rc - b * sqrt( 1.0 - 2.0 * get_mass() / rc ) );
  return double( rc - b * sqrt( 1.0 - 2.0 * get_mass_over_r() *get_rspot() / rc ) );
}

/*****************************************************/
// used in oblate shape of star
/*****************************************************/
double OblDeflectionTOA::b_from_psi_ingoing_zero_func ( const double& b_R, const double& b_R_max, const double& psimax,
							const double& rspot, 
							const double& psi ) const { 
  return double( psi - this->psi_ingoing(b_R, b_R_max, psimax, rspot, &curve.problem) );

}

/*****************************************************/
/*****************************************************/
/*double OblDeflectionTOA::b_from_psi_outgoing_zero_func ( const double& b_R, 
							 const double& psi,
							 const double& b_max, 
							 const double& psi_max 
							  ) const {
  return double( psi - this->psi_outgoing_u( b_R, b_max, psi_max, &curve.problem ) );
}
*/

/*****************************************************/
/* OblDeflectionTOA::TrapezoidalInteg                */
/*                                                   */
/* Numerically integrates using the trapezoid rule.  */
/*****************************************************/
double OblDeflectionTOA::TrapezoidalInteg ( const double& a, const double& b, 
                                            double (*func)(double x, bool *prob),
					   						const long int& N ) {
  	if ( a == b ) return double(0.0);

  	double integral(0.0);
  	//const long int N( OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_N );
  	const double power( OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_POWER );

  	for ( int i(1); i <= N; i++ ) {
    	double left, right, mid, width, fmid;
    	left = TrapezoidalInteg_pt( a, b, power, N, i-1, &curve.problem );
    	right = TrapezoidalInteg_pt( a, b, power, N, i, &curve.problem );
    	mid = (left + right) / 2.0;
    	width = right - left;
    	fmid = (*func)(mid, &curve.problem); //Double Check
    
    	integral += (width * fmid);
  	}

  	return integral;
}

/*****************************************************/
/* OblDeflectionTOA::TrapezoidalInteg                */
/*                                                   */
/* Numerically integrates using the trapezoid rule.  */


/*******************************************************/
/* OblDeflectionTOA::TrapezoidalInteg_pt               */
/*                                                     */
/* Called in TrapezoidalInteg                          */
/* Power spacing of subdivisions for TrapezoidalInteg  */
/* i ranges from 0 to N                                */
/*******************************************************/
double OblDeflectionTOA::TrapezoidalInteg_pt ( const double& a, const double& b, 
					      					   const double& power, const long int& N,
					      					   const long int& i, bool *prob ) {
  	double dummy;
  	if ( N <= 0 ) {
    	std::cerr << "ERROR in OblDeflectionTOA::TrapezoidalInteg_pt: N <= 0." << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}
  	if ( i < 0 || i > N ) {
    	std::cerr << "ERROR in OblDeflectionTOA::TrapezoidalInteg_pt: i out of range." << std::endl;
    	*prob = true;
    	dummy = -7888.0;
    	return dummy;
  	}
  	return double( a + (b - a) * pow( 1.0 * i / N , power ) );
}
/*******************************************************/
/* OblDeflectionTOA::Integration                       */
/*                                                     */
/* Integrates! By calling TrapezoidalInteg             */
/* Here, N = TRAPEZOIDAL_INTEGRAL_N                    */
/*******************************************************/
inline double OblDeflectionTOA::Integration ( const double& a, const double& b, 
                                              double (*func)(double x, bool *prob),
					     					  const long int& N ) {
  	return TrapezoidalInteg( a, b, func, N );
  	//return MATPACK::AdaptiveSimpson( a, b, func, OblDeflectionTOA::INTEGRAL_EPS );
}

