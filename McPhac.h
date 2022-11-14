

double McPHACC3new(double E, double cos_theta, int theta_index, double T, double lgrave, int i_lgrav, double ggvec[4], class LightCurve* mexmcc);



//double AtmosEBandFlux3new( unsigned int model, double cos_theta, double T, double lgrav, double E1, double E2, class LightCurve mexmcc);
double AtmosEBandFlux3new( unsigned int model, double cos_theta, int theta_index, double T, double lgrav, int i_lgrav, double gvec[4], double E1, double E2, class LightCurve mexmcc);


int th_index(double cos_theta, class LightCurve* mexmcc);

int th_index_nsx(double cos_theta, class LightCurve* mexmcc);
