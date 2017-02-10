This directory contains synthetic data, and auxiliary files, for NICER-like observations of a PSR J0437-like source.  

The files are:

1. J0437dataNICERplskyAGNNH2e20.dat
   The first column is the phase (in cycles, after folding on the rotational frequency).
   The second column is the centroid of the energy bin (in keV; there are 30 bins of width 0.1 keV each)
   The third column is the number of counts in each phase-energy bin.
   The counts are the sum of the counts from the two spots and from all other sources of X-rays, including the diffuse sky background, the NICER instrumental background, an unmodulated power law, and a nearby AGN.


2. niceravg.area and nicer_xti_CDR-BOM.area
   NICER effective area curves; Zaven kindly provided the source file nicer_xti_CDR-BOM.area.  niceravg.area is an average over 0.1 keV bins; note that the first bin averages over 0.1 through 0.15 keV, whereas other bins average over 0.16 through 0.25 keV, 0.26 through 0.35 keV, and so on.


3. deflection.bin
   This is a table, in binary format (this speeds up the I/O), which has five columns:
   1. R/M
   2. Initial angle of local ray from local radial direction, in radians
   3. Total angular distance traveled by ray to infinity, in radians
   4. Propagation time to infinity, minus propagation time to infinity for a radial ray, for a 1 solar mass star of the indicated R/M.
   5. Ratio of the solid angle subtended by an infinitesimal bundle of rays at the starting R/M, to the solid angle subtended by the same bundle at infinity.

To read this binary table we use the following C code:

   defltable=fopen("deflection.bin","rb");

   for (i=0; i<NRdefl; i++)
   {
      for (j=0; j<Nangdefl; j++)
      {
               fread(&Rdefl[i],sizeof(double),1,defltable);
               fread(&angrad[j],sizeof(double),1,defltable);
               fread(&deflangle[i][j],sizeof(double),1,defltable);
               fread(&defltime[i][j],sizeof(double),1,defltable);
               fread(&deflconv[i][j],sizeof(double),1,defltable);
      }
   }
   fclose(defltable);

   Here NRdefl=490 is the number of values of R/M (from 3.1 to 8 in steps of 0.01) and Nangdefl=1200 is the number of angles from the local radial direction (from 0 to 2.2 radians, uniformly spaced; we consider angles greater than pi/2 because we are treating oblate stars).


4. Hatm8000dT0.05.bin
   This is the H atmosphere table, in binary, which we generated using McPHAC with resolutions such that the bolometric flux is within 0.1% of \sigma_SB T_eff^4 at all temperatures.  It gives the specific intensity as a function of the effective temperature, surface gravity, photon energy, and angle from the normal.  The columns are:
   1. log_10 of the effective temperature in Kelvin
   2. log_10 of the surface gravity in cm s^{-2}
   3. log_10 of the photon energy in keV
   4. cosine of the angle to the surface normal (unitless)
   5. Specific intensity in cgs units divided by Teff^3 (this improves the interpolation properties)

To read this binary table we use the following C code:

   Hspecttable=fopen("Hatm8000dT0.05.bin","rb");
   for (i=0; i<NlogTeff; i++)
   {
      for (j=0; j<Nlogg; j++)
      {
         for (k=0; k<Nnu; k++)
         {
            for (k1=0; k1<Nmu; k1++)
            {
               fread(&logTefftab[i],sizeof(double),1,Hspecttable);
               fread(&loggtab[j],sizeof(double),1,Hspecttable);
               fread(&lognutab[k],sizeof(double),1,Hspecttable);
               fread(&mutab[k1],sizeof(double),1,Hspecttable);
               fread(&Inutab[k][k1][j][i],sizeof(double),1,Hspecttable);
            }
         }
      }
   }
   fclose(Hspecttable);

Here NlogTeff=29 (log_10 Teff from 5.1 through 6.5 in steps of 0.05), Nlogg=11 (log_10 g from 13.7 through 14.7 in steps of 0.1), Nnu=100 (log_10 E/keV = -1.30369 through 2.04273 in steps of 0.338) and Nmu=50 (mu=cos(theta)=0.015629 through 0.999710 in steps that are roughly uniform in theta (dtheta=0.031 radians); note that the angle gridding is determined within McPHAC rather than being set by the user.


5. skyback.txt
   Diffuse sky + Soyuz background file, provided by Slavko.  
   First column is photon energy in keV, second is count rate in counts per second


6. We also provide a background from a nearby AGN; Zaven recommended as a reference Wang et al. 1998, MNRAS, 293, 397

7. For ISM absorption we used the Morrison and McCammon 1983, ApJ 270, 119 cross sections (this is used in the wabs module of XSPEC).

8. We put in an unmodulated power law with a number index of 1.6, based on the Guillot et al. 2016 (arXiv:1512.03957) fit including NuSTAR data.

9. The total observation time we simulated was 10^6 (one million) seconds.

10. The parameters of the run were based on the Bogdanov 2013 fits:
For the two spots:
1. R=13 km
2. M=1.44 Msun (thus R/M=6.1138)
3. thetaobs=0.733 (all angles in radians; our prior here was uniform from 0.71 to 0.76)
4. thetac1=0.6283
5. dtheta1=0.01
6. kT1=0.231139 keV (in surface comoving frame)
7. thetac2=2.077
8. dtheta2=0.33
9. kT2=0.0577846 keV (in surface comoving frame)
10. dphi between spots: 0.5625 cycles
11. distance fixed at 156.3 pc
12. NH=2x10^{20} cm^{-2}
13. Rotational frequency as seen by a distant observer: 173.6 Hz

Total counts for spot 1: ~640,000; spot 2: ~370,000; power law: ~100,000; AGN: ~400,000; diffuse sky + Soyuz: ~300,000, all from 0.1 to 3 keV.