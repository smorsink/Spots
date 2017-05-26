2017 May 24 - WCG Ho

Format similar to that used for Xspec models NSX, NSMAXG, NSATMOS, and NSA.

First 8 rows:
Row 1: = 15, number of effective temperatures NlogTeff
Row 2: logTeff (K) = 5.1,5.2,...,6.5
Row 3: = 11, number of surface gravities Nlogg
Row 4: logg (cm s^-2) = 13.7,13.8,...,14.7
Row 5: = 137, number of photon energies NlogE
Row 6: logE (keV) = -1.32,-1.30,...,1.4
Row 7: = 67, number of angles Nmu [=cos(theta)]
Row 8: mu, ~ every 1.5 degrees in range [1,1.0e-6]

Remaining NlogTeff*Nlogg*NlogE = 15*11*137 = 22605 rows:
Row     9: (logTeff,logg)=(5.1,13.7) and logE=-1.32:
 67 values of specific intensity Inu (erg s^-1 cm^-2 Hz^-1 ster^-1)
 for mu=1, mu=0.999999995,..., mu=1.0e-6
Row    10: (logTeff,logg)=(5.1,13.7) and logE=-1.30: 67 values of Inu
...
Row   145: (logTeff,logg)=(5.1,13.7) and logE=+1.40: 67 values of Inu
Row   146: (logTeff,logg)=(5.1,13.8) and logE=-1.32: 67 values of Inu
...
Row  1515: (logTeff,logg)=(5.1,14.7) and logE=+1.40: 67 values of Inu
Row  1516: (logTeff,logg)=(5.2,13.7) and logE=-1.32: 67 values of Inu
...
Row 22613: (logTeff,logg)=(6.5,14.7) and logE=+1.40: 67 values of Inu

Note: For low logTeff models, calculations do not extend up to logE=1.4 since
very low Inu at these energies.  Therefore Inu is set to 1.0e-60 for energies
outside model ranges (see eg rows 145).
