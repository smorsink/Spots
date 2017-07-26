

double McPHACC3new(double E, double cos_theta, int theta_index, double T, double lgrave, double ggvec[4], class LightCurve* mexmcc);



//double AtmosEBandFlux3new( unsigned int model, double cos_theta, double T, double lgrav, double E1, double E2, class LightCurve mexmcc);
double AtmosEBandFlux3new( unsigned int model, double cos_theta, int theta_index, double T, double lgrav, double gvec[4], double E1, double E2, class LightCurve mexmcc);


int th_index(double cos_theta, class LightCurve* mexmcc);
