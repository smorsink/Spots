
double BlackBody( double T, double E );


double EnergyBandFlux( double T, double E1, double E2 );


double Bradt_flux_integrand( double x );

double BlackBodyTabLogLog( double T, double E, class LightCurve mexmcc );
double BlackBodyTabLogLinear( double T, double E, class LightCurve mexmcc );

double HopfTab(double cosalpha, class LightCurve mexmcc);

double BlackBodyHopfTab( double T, double E, double cosalpha, class LightCurve mexmcc );

double Hopf( double cosalpha);
