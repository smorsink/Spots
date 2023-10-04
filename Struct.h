/***************************************************************************************/
/*                                    Struct.h

    This code creates the many data storage structures used in other pieces of code.
    
    Based on code written by Coire Cadeau and modified by Sharon Morsink and 
    Abigail Stevens.
    
    Does Coire still get to claim copyright? I don't think we personally own our research 
    at the University of Alberta. Anyways, here it is:
    (C) Coire Cadeau, 2007; Source (C) Coire Cadeau 2007, all rights reserved.
*/
/***************************************************************************************/

#ifndef STRUCT_H
#define STRUCT_H

#include <exception>
#include <vector>
#include <float.h>

#define NN 100            // lookup table for bending angle (deflection angle) calculation
#define MAX_NUMBINS 256  //  how many time bins the light curve is cut up into
#define MIN_NUMBINS 256   // We need a minimum number of bins since the curves won't be accurate if we use too few bins.
#define NCURVES 1500
//#define NCURVES 500
// NCURVES is the number of energy bands that the code will compute before the instrument response is applied.

//#define FBANDS 300        // Final number of energy channels
#define CBANDS 500         // Number of energy bands computed
#define MR 1000             // Maximum number of m/r values
#define NUM_NICER_CHANNELS 1500   // Number of NICER energy channels
#define NUM_NICER_BINS 32



struct Parameters {      // local bit of spot information
  double theta;          // Angle between the NS spin axis and the spot bin; in radians
  double theta_c;        // Angle between the NS spin axis and the centre of the spot; in radians
  double phi_0;          // Azimuthal angular location of the centre of the spot; in radians
  double dS;             // Area of grid bit
  double incl;           // Inclination angle (between NS spin axis and observer's line of sight); in radians
  double temperature;    // Temperature of the spot, in the spot's frame; in keV or Kelvin (depends on Temperature Flag)
  double mass;           // Mass of the star; unitless in here
  double radius;         // Radius of the spot; unitless in here
  double req;            // Radius at the equator
  double rspot;          // Radius at the spot
  double mass_over_r;    // Dimensionless mass divided by radius
  double omega;          // Spin frequency of the star; unitless in here
  double omega_bar_sq;	 // Rotation parameter, omega_bar_squared in AlGendy & Morsink 2014
  double cosgamma;       // Gamma is angle between true normal to surface and radial vector; probably in radians
  double bbrat;          // Ratio of blackbody-to-comptonization of spectrum
  double ts;             // Time shift between 0 and 1 (beginning of period, end of period)
  double E_band_lower_1; // Lower bound of energy band for flux calculation; in keV
  double E_band_upper_1; // Upper bound of energy band for flux calculation; in keV
  double E_band_lower_2; // Lower bound of energy band for flux calculation; in keV
  double E_band_upper_2; // Upper bound of energy band for flux calculation; in keV
  double distance;       // Distance from earth to NS; in meters
  double E0; // NICER
  double L1; // NICER
  double L2; // NICER
  double DeltaE; // width of energy band
  double DeltaLogE; //logarithmic width of energy band
  double theta_k[256]; // location of spot bins
  double phaseshift;   // Overall phase shift of spot
  double phi_k[256];   // phi location of edge of spot
  double dtheta[256]; // width of spot bins
  double gamma_k[256]; // value of lorentz gamma at the spot bin
};

struct Flags {
	double shift_t;               // for shifting of time, to match light curve phases
	bool infile_is_set;           // if an input file has been set
	bool ignore_time_delays;      // if we should ignore time delays
	bool bend_file;				  // yes means bend file has been read, passing to bend.
	unsigned int spectral_model;  // stating which model we're using -- definitions of models given elsewhere
	unsigned int beaming_model;   // stating which model we're using -- definitions of models given elsewhere
    unsigned int ns_model; // Shape model: 3=spherical; 1=oblate
    unsigned int attenuation;	  // attenuation flag for the four NICER target files
    unsigned int inst_curve;      // instrument response curve flag for NICER
    unsigned int spotshape;
  bool kelvin; // true means temperature in Kelvin; false means in keV
  bool logEflag; // true means energy bands evenly spaced in log space; false means linear space
};

class Defl {
	public:                       // allocates the memory for the lookup table -- not evenly spaced
	double *psi_b; // psi values for specific M/R
	double *b_psi; // b values for specific M/R
	double *dcosa_dcosp_b; // dcos(alpha)/dcos(psi) values for specific M/R
	double *toa_b; // toa values for specific M/R
	double psi_max;               // largest possible value of psi
	double b_max;                 // largest possible value of b
	double b_R_max;               // largest possible value of b/R
	int num_mr;              // number of M/R values
	double *mr;              // values of M/R
	double **b;        // values of b/R for different M/R
	double **psi;      // values of psi for different M/R
	double **dcosa;    // values of dcosa/dcospsi for different M/R
	double **toa;      // values of toa for different M/R

	//class OblDeflectionTOA defltoa;
};

class LightCurve {                     // Stores all the data about the light curve!
	public:
	double t[MAX_NUMBINS];                 // one-dimensional array that hold the value of time of emission; normalized between 0 and 1
	double f[CBANDS][MAX_NUMBINS];        // two-dimensional array of fluxes (one one-dimensional array for each energy curve)
	double elo[CBANDS];                   // Lowest energy included in energy band "p"
	double ehi[CBANDS];                   // Highest energy included in energy band "p"
	//bool visible[MAX_NUMBINS];             // is the spot visible at that point
	double t_o[MAX_NUMBINS];               // the time in the observer's frame; takes into account the light travel time
	double cosbeta[MAX_NUMBINS];           // as seen in MLCB17 (cos of zenith angle, between the norm vector and initial photon direction)
	double eta[MAX_NUMBINS];               // doppler shift factor; MLCB33
	double psi[MAX_NUMBINS];               // bending angle; MLCB15
	double cospsi[MAX_NUMBINS];
	//double R_dpsi_db[MAX_NUMBINS];         // derivative with respect to b of MLCB20 times the radius
	double b[MAX_NUMBINS];                 // impact parameter; defined in dimensionless units
	//double dcosalpha_dcospsi[MAX_NUMBINS]; // appears in MLCB30
	double dOmega_s[MAX_NUMBINS];          // solid angle weight factor; MLCB30
	struct Parameters para;                // parameters from above; para is like i, Parameters is like Integer
	struct Flags flags;                    // flags from above
	class Defl defl;                       // deflection from above
	double *mccinte;  						   // intensity values
	double *mccangl;
	double *mcloget;                        // log(photonenergy/temperature)
	double *mclogTeff;
	double *mclogg;
	int Nmu;
	int NlogTeff;
	int Nlogg;
	int NlogE;
	int Npts;
	unsigned int numbins;                  // Number of time or phase bins for one spin period; Also the number of flux data points
	unsigned int numbands;
	unsigned int fbands; // The Final Number of energy bands. Value = 300 right now
	unsigned int tbands; // The Total Number of energy bands that are required by the response matrix = 350 right now
	unsigned int cbands; // The Computed Number of energy bands in approximation. Either using 350 for no approx or 35. 
	bool eclipse;                          // True if an eclipse occurs
	bool ingoing;                          // True if one or more photons are ingoing
	bool problem;                          // True if a problem occurs
	unsigned int count;                    // for outputting command line args in Chisquare, chi.cpp

};

class NICERCurve {                     // Stores all the data about the light curve!
	public:
	double t[NUM_NICER_BINS];                 // one-dimensional array that hold the value of time of emission; normalized between 0 and 1
	double f[NCURVES][NUM_NICER_BINS];        // two-dimensional array of fluxes (one one-dimensional array for each energy curve)
	double elo[NCURVES];                   // Lowest energy included in energy band "p"
	double ehi[NCURVES];                   // Highest energy included in energy band "p"

	unsigned int numbins;                  // Number of time or phase bins for one spin period; Also the number of flux data points
	unsigned int numbands;

};


class DataStruct {             // if reading in data, this would be the experimental data
	public:
	double *t;                       // time
	double *f[NCURVES];              // flux; an array of double pointers
	double *err[NCURVES];            // error bars
	double chisquare;                // total chi squared 
	double chi[NCURVES];            // chi square for each energy band
	double shift;                    // if we need it; unused
	unsigned int numbins;            // Number of time or phase bins for one spin period; 
	//Also the number of flux data points
	unsigned int numbands;
	};


class Instrument { // Information about the instrument's response matrix, etc
 public:
  unsigned long *start;
  double **response; // Response matrix (combined ARF and RFM)
  double *elo;       // Lowest energy (in keV) photon in channel p
  double *ehi;       // Highest energy photon in channel p
};

class ISM { // Attenuation by the ISM
 public:
  double *attenuation; // Vector of attenuations for each photon energy
  double *energy;     // Vector of photon energies
};


#endif // STRUCT_H
