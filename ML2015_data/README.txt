This directory contains data for Figures 3a, 3b, 3c, and 3d from Miller and Lamb 2015, ApJ, 808, 31.

Note that:

1. We used the Morsink et al. 2007 (ApJ, 663, 1244) oblate Schwarzschild shape model, rather than the AlGendy and Morsink 2014 (ApJ, 791, 78) shape model.

2. The synthetic waveforms were computed and analyzed *without* the Lorentz factor that is now agreed to be correct.

3. For all except case 3d, the synthetic data are for a single, uniform, circular spot.  For case 3d, the synthetic data involve a gradient in temperature in the N-S direction for a circular spot, but in the analysis we make our standard assumption that the spot has uniform emission in addition to being circular.

4. We assumed that the spectrum was a blackbody as seen in the local comoving frame of the surface, but that the beaming pattern was the Hopf pattern: the intensity at an angle $\alpha^\prime$ from the normal, as seen in the local comoving frame, is given by 
   $0.42822+0.92236*\cos(\alpha^\prime)-0.085751*\cos^2(\alpha^\prime)$
times the blackbody intensity.

5. For all except case 3d we assumed that the we know the temperature of the emission as seen by the *distant* observer (rather than that we know the temperature as seen in the comoving frame).  In cases 3a, 3b, and 3c we generated the synthetic data assuming that kT=2.0 keV as seen in the local comoving frame.  Thus, for example, for R/M=5, we fix the temperature seen by a distant observer to $kT_{\rm distant}=2.0~{\rm keV}(1-2/5)^{1/2}$.

For case 3d, the temperature was a free parameter but the distance was assumed to be exactly 10 kpc in our analysis.

6. We assume that we know the rotational frequency and the duration of the observation, but that we do not know any other information about the systems.  For example, for all except case 3d we searched distances from 5-20 kpc (we used 10 kpc to generate the synthetic data, and assumed d=10 kpc in the special case 3d), and we searched spot angular radii from 0 degrees to 90 degrees, as well as searching over the full range of observer and spot center inclinations.

7. We aimed for roughly a million expected counts from the spot, and for normalization purposes we assumed a detector of area 1 cm^2 that is perfect; it has 100% response at all energies.

8. We do not have any interstellar absorption in the data; N_H=0.

9. We have several million counts of unmodulated background in the data, but we do not assume that we know anything about that background when we perform the analysis.

10. The data are in three columns: (1) phase (16 bins; when analyzing the data we do not assume that we know the phase of peak counts, but instead fit that as a free parameter), (2) energy in keV (the listed energy is the midpoint of the energy range; the bin widths are 0.3 keV), and (3) counts in that phase-energy bin.

These factors should be taken into account in the analysis of the data.  See Miller and Lamb 2015 for additional information.

Details about the individual data sets are:

1. Data for Figure 3a (dataML2015Fig3a.dat).  The rotational frequency is 600 Hz as seen by a distant observer, and R/M=5.0.  To reach roughly a million expected counts, we used a duration of 320209.39 seconds to go along with the 1.0 cm^2 effective area of our perfect detector.

2. Data for Figure 3b (dataML2015Fig3b.dat).  The rotational frequency is 300 Hz as seen by a distant observer, and R/M=5.0.  To reach roughly a million expected counts, we used a duration of 324223.7 seconds to go along with the 1.0 cm^2 effective area of our perfect detector.

3. Data for Figure 3c (dataML2015Fig3c.dat).  The rotational frequency is 600 Hz as seen by a distant observer, and R/M=6.35.  To reach roughly a million expected counts, we used a duration of 140382.8 seconds to go along with the 1.0 cm^2 effective area of our perfect detector.

4. Data for Figure 3d (dataML2015Fig3d.dat).  The rotational frequency is 600 Hz as seen by a distant observer, and R/M=5.  To reach roughly a million expected counts, we used a duration of 254967.4 seconds to go along with the 1.0 cm^2 effective area of our perfect detector.  As indicated above, for this special case we assumed d=10 kpc in our analysis, but allowed the spot temperature to be a free parameter.
