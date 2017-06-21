// oldstuff.c -- routines that aren't being used anymore.

/**************************************************************************************/
/* Gray:																			  */
/*		computes and returns the limb darkening factors for a Gray electron -         */
/*		scattering atmosphere (Hopf function) using values from Chandrasekhar's       */
/*		"Radiative Transfer", table 13												  */
/*																					  */
/* pass: cosine = curve.cosbeta[i] * curve.eta[i]                                     */
/**************************************************************************************/
double Gray( double cosine ) {
	
	/**********************************/
   	/* VARIABLE DECLARATIONS FOR Gray */
    /**********************************/
	
    double F[11],   // limb darkening values for electron-scattering atmosphere
           mu[11];  // a table of cos(theta) values
    int i;          // loop variable

    for ( i = 0; i <= 10; i++ ) {
        mu[i] = 0.0 + i*0.1;
    }
    /*
    // Value from Mihalas's "Stellar Atmospheres"
    F[0] = 0.4330;
    F[1] = 0.5401;
    F[2] = 0.6280;
    F[3] = 0.7112;
    F[4] = 0.7921;
    F[5] = 0.8716;
    F[6] = 0.9501;
    F[7] = 1.0280;
    F[8] = 1.1053;
    F[9] = 1.1824;
    F[10] = 1.2591;
    */
  
    // Values from Chandrasekhar's "Radiative Transfer", Table 13
    F[0] = 1.0;
    F[1] = 1.2647;
    F[2] = 1.48009;
    F[3] = 1.68355;
    F[4] = 1.88105;
    F[5] = 2.07496;
    F[6] = 2.2665;
    F[7] = 2.45639;
    F[8] = 2.64503;
    F[9] = 2.83274;
    F[10] = 3.01973;
  
    for ( i = 0; i <= 10; i++ ) {
        F[i] *= (0.5/1.194);
    }

    int index(0);
    //i=1;

    if (cosine == 1.0)
    	index = 10;
    else {
        while ( cosine >= mu[index]) 
        	index++;
    }
    // cosine lies in the range   mu[index-1] <= cosine <= mu[index]

    if (cosine == mu[index-1]) 
    	return F[index-1];
    else {
        if (cosine == mu[index]) 
        	return F[index];
        else {   
            return ( F[index-1] + (cosine - mu[index-1]) * (F[index] - F[index-1])/(mu[index]-mu[index-1]) );
        }
    }
} // end Gray



/**************************************************************************************/
/* McPHAC Hydrogen:                                                                   */
/*    This version works for the original McPHAC files from Cole                      */
/*                                                                                    */
/**************************************************************************************/
void Read_McPHACC(double T, double M, double R){
    
    double delta, lt, lgrav, temp, dump, real_T;
    //int size_logt(13), size_lsgrav(8), size_mu(23), i_lt, i_lgrav, n_lt, n_lgrav, size_ener(125), size_set, skip_to, skip_two;
    int i_lt, i_lgrav, n_lt, n_lgrav;
    char s1[40], s2[40], s3[40], s4[40], atmodir[1024], cwd[1024];
    std::vector<double> freq, logt, lsgrav;

    //setting values of lt and lgrav based on input T, M, and R. Also sets to load using first mu value.   
    M = Units::nounits_to_cgs(M, Units::MASS);
    R = Units::nounits_to_cgs(R, Units::LENGTH);
    delta = 1 / sqrt(1 - (2 * Units::G * M / (R * Units::C * Units::C)));
    lgrav = log10(delta * Units::G * M / (R * R));
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));
    cout << "temperature in log(K) is " << lt << endl;
    cout << "gravity in log(cgs units) is " << lgrav << endl;

    //obtain helium atmosphere parameters
    getcwd(cwd, sizeof(cwd));
    sprintf(atmodir,"%s/atmosphere/mcphacc",cwd);
    chdir(atmodir);

    i_lt = (lt-5.1)/0.05; //if we need to load 1st temperature, i_lt = 0. this is discrete math
    i_lgrav = (lgrav-13.7)/0.1;

    n_lt = i_lt+1;
    n_lgrav = i_lgrav+1;

    sprintf(s1,"mcphacc_T%d_g%d.txt",n_lt,n_lgrav);
    sprintf(s2,"mcphacc_T%d_g%d.txt",n_lt+1,n_lgrav);  
    sprintf(s3,"mcphacc_T%d_g%d.txt",n_lt,n_lgrav+1); 
    sprintf(s4,"mcphacc_T%d_g%d.txt",n_lt+1,n_lgrav+1);

    cout << s1 << endl;
    cout << s2 << endl;
    cout << s3 << endl;
    cout << s4 << endl; 

    real_T = pow(10.0,lt);

    ifstream H_table1;
    H_table1.open(s1);
    ifstream H_table2;
    H_table2.open(s2);
    ifstream H_table3;
    H_table3.open(s3);
    ifstream H_table4;
    H_table4.open(s4);

    if(H_table1.is_open()){
    	for (int i = 1; i <= 5000; i++) {
    		H_table1 >> dump;
    		H_table1 >> dump;
    		H_table1 >> temp;
    		temp = temp*pow(real_T,3);
        	I.push_back(temp);
        }
        H_table1.close();    	
    }else{
    	cout << s1 << " is not found." << endl;
    }

    if(H_table2.is_open()){
    	for (int i = 1; i <= 5000; i++) {
    		H_table2 >> dump;
    		H_table2 >> dump;
    		H_table2 >> temp;
    		temp = temp*pow(real_T,3);
        	II.push_back(temp);
        }
        H_table2.close();    	
    }else{
    	cout << s2 << " is not found." << endl;
    }     

    if(H_table3.is_open()){
    	for (int i = 1; i <= 5000; i++) {
    		H_table3 >> dump;
    		H_table3 >> dump;
    		H_table3 >> temp;
    		temp = temp*pow(real_T,3);
        	III.push_back(temp);
        }
        H_table3.close();    	
    }else{
    	cout << s3 << " is not found." << endl;
    }

    if(H_table4.is_open()){
    	for (int i = 1; i <= 5000; i++) {
    		H_table4 >> dump;
    		H_table4 >> dump;
    		H_table4 >> temp;
    		temp = temp*pow(real_T,3);
        	IIII.push_back(temp);
        }
        H_table4.close();    	
    }else{
    	cout << s4 << " is not found." << endl;
    }

    //cout << i_lt << " " << n_lt << " " << i_lgrav << " " << n_lgrav << endl;
	/* 
    ifstream H_table1;
    H_table1.open("cole_mcphac_table2.txt");

    if(H_table1.is_open()){

    	skip_to = (i_lt * 11 + i_lgrav)*5000; //skipping unwanted lines.

    	for (int j = 1; j <= skip_to; j++) {
     		H_table1 >> dump;	//skipping temperature
     		H_table1 >> dump;   //skipping specific gravity
     		H_table1 >> dump;	//skipping energy
     		H_table1 >> dump;	//skipping angle
     		H_table1 >> dump;	//skipping intensity
    	}

    	for (int j = 1; j <= 5000; j++) {
    		H_table1 >> temp;	//skipping temperature
    		X1 = temp;
     		H_table1 >> temp;   //skipping specific gravity
     		Y1 = temp;
     		H_table1 >> dump;	//skipping energy
     		H_table1 >> dump;	//skipping angle
    		H_table1 >> temp;
    		temp = temp*pow(real_T,3);
            I.push_back(temp);
    	}

    	for (int j = 1; j <= 5000; j++) {
    		H_table1 >> dump;	//skipping temperature
     		H_table1 >> temp;   //skipping specific gravity
     		Y2 = temp;
     		H_table1 >> dump;	//skipping energy
     		H_table1 >> dump;	//skipping angle
    		H_table1 >> temp;
    		temp = temp*pow(real_T,3);
            III.push_back(temp);
    	}

    	skip_two = (11-2)*5000;

    	for (int j = 1; j <= skip_two; j++) {
     		H_table1 >> dump;	//skipping temperature
     		H_table1 >> dump;   //skipping specific gravity
     		H_table1 >> dump;	//skipping energy
     		H_table1 >> dump;	//skipping angle
     		H_table1 >> dump;	//skipping intensity
    	}

    	for (int j = 1; j <= 5000; j++) {
    		H_table1 >> temp;	//skipping temperature
    		X2 = temp;
     		H_table1 >> dump;   //skipping specific gravity
     		H_table1 >> dump;	//skipping energy
     		H_table1 >> dump;	//skipping angle
    		H_table1 >> temp;
    		temp = temp*pow(real_T,3);
            II.push_back(temp);
    	}

    	for (int j = 1; j <= 5000; j++) {
    		H_table1 >> dump;	//skipping temperature
     		H_table1 >> dump;   //skipping specific gravity
     		H_table1 >> dump;	//skipping energy
     		H_table1 >> dump;	//skipping angle
    		H_table1 >> temp;
    		temp = temp*pow(real_T,3);
            IIII.push_back(temp);
    	}
    	

    }else{
        cout << "McPHAC (Cole's) file is not found" << endl;
    }
    H_table1.close();
    */

    cout << "finished loading Cole's McPHAC intensities" << endl; 

    ifstream H_table5;
    H_table5.open("mcphacc_T1_g1.txt");
    
    if(H_table5.is_open()){
        H_table5 >> temp;
        temp = pow(10.0,temp);
        Es.push_back(temp);	//recording first energy
        //cout << temp << endl;
        temp = temp * 1E3 * Units::EV / Units::H_PLANCK;
        F.push_back(temp);
        for (int i = 1; i <= 49; i++){
        	H_table5 >> temp;
        	mu_2.push_back(temp);	//recording first 49 angles
    		H_table5 >> dump;	//skipping intensity
     		H_table5 >> dump;   //skipping energy        	
        }
        H_table5 >> temp;
        mu_2.push_back(temp);	//recording 50th angle;
        H_table5 >> dump;	//skipping intensity

        //entering second energy group
        for (int j = 2; j <= 100; j++){
        	H_table5 >> temp;
        	temp = pow(10.0,temp);
        	Es.push_back(temp);	//recording second and up to last energy
        	//cout << temp << endl;
        	temp = temp * 1E3 * Units::EV / Units::H_PLANCK;
        	F.push_back(temp);
        	for (int i = 1; i <= 49; i++){
    			H_table5 >> dump;	//skipping angle
     			H_table5 >> dump;   //skipping intensity
     			H_table5 >> dump;	//skipping energy
        	}
        	H_table5 >> dump;	//skipping angle
     		H_table5 >> dump;   //skipping intensity      	
        }
    }else{
        cout << "McPHAC (Cole's) file is not found (while in second reading stage) " << endl;
    }
    //cout << I[0] << " " << F[0] << " " << Es[0] << " " << mu_2[0] << endl;
    H_table5.close();   

    X = pow(10.0,lt);
    Y = pow(10.0,lgrav);
    X1 = pow(10.0,5.1+0.05*i_lt);
    X2 = pow(10.0,5.1+0.05*n_lt);
    Y1 = pow(10.0,13.7+0.1*i_lgrav);
    Y2 = pow(10.0,13.7+0.1*n_lgrav);
    //cout << T << " " << lt << endl;
    chdir(cwd);
	/*
    for (int i = 0; i <= 99; i++){
    	cout << Es[i] << endl;
    }
    */
    
    
}






// This is version makes use of Cole's version of McPhac
double AtmosEBandFlux3( unsigned int model, double cos_theta, double T, double M, double R, double E1, double E2, class LightCurve mexmcc){

    int e1_dex(0);              // largest energy index that's smaller than E1
    int e2_dex(0);              // largest energy index that's smaller than E2
    int n_steps(0);             // number of energy points within band
    int ener_size(0);         // number of energy choices
    //int e_dex;                  // energy index of current step
    //double current_e;           // central energy of current step
    //double current_n;           // integrated flux in current step
    double flux(0.0);           // total integrated flux

    ener_size = 125;

	double ener_spacing = pow(10.0,0.0338);
    double first_ener = 0.004969469;
    double ener_index = log(E1 / first_ener) / log(ener_spacing);
    e1_dex = (int) ener_index;
    ener_index = log(E2 / first_ener) / log(ener_spacing);
    e2_dex = (int) ener_index;

    n_steps = e2_dex - e1_dex;

    if (n_steps == 0){ // zero energy points within bandwidth: (4.1.3) one trapzoid
        //cout << "0 steps" << endl;
        flux = (E2 - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2);
    }
    if (n_steps == 1){ // one energy points within bandwidth: (4.1.3) two trapzoids
        //cout << "1 steps" << endl;
        int e_dex = e1_dex+1; // index of the energy point
        double e_m = first_ener*pow(ener_spacing,e_dex); // energy point in keV
        flux = (e_m - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC4(e1_dex,cos_theta, T, M, R, mexmcc) / e_m);  // first trapezoid
        flux += (E2 - e_m) / 2 * (McPHACC4(e_dex,cos_theta, T, M, R, mexmcc) / e_m + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2); // second trapezoid
    }
    if (n_steps == 2){ // two energy points within bandwidth: (4.1.3) three trapzoids
        //cout << "2 steps" << endl;
        int el_dex = e1_dex+1; // index of the first energy point within bandwidth
        int eu_dex = e2_dex;   // index of the last energy point within bandwidth
        double e_l = first_ener*pow(ener_spacing,el_dex); // first energy point in keV
        double e_u = first_ener*pow(ener_spacing,eu_dex); // last energy point in keV
        double h = log(e_u)-log(e_l); // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) / e_l);  // first trapezoid
        flux += h * (McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) + McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc)) / 2;                // middle trapezoid, exact
        flux += (E2 - e_u) / 2 * (McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc) / e_u + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2); // last trapezoid
    }
    if (n_steps == 3){ // three energy points within bandwidth: (4.1.3, 4.1.4) Simpson's + two trapezoids
        //cout << "3 steps" << endl;
        int el_dex = e1_dex+1; // index of the first energy point within bandwidth
        int em_dex = e1_dex+2; // index of the middle energy point within bandwidth
        int eu_dex = e2_dex;   // index of the last energy point within bandwidth
        double e_l = first_ener*pow(ener_spacing,el_dex); // first energy point in keV
        double e_u = first_ener*pow(ener_spacing,eu_dex); // last energy point in keV
        double h = (log(e_u)-log(e_l)) / 2; // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) / e_l);  // first trapezoid
        flux += h * (McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) + McPHACC4(em_dex,cos_theta, T, M, R, mexmcc) * 4 + McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc)) / 3; // Simpson's, exact
        flux += (E2 - e_u) / 2 * (McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc) / e_u + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2); // second trapezoid
    }
    if (n_steps == 4){ // four energy points within bandwidth: (4.1.3, 4.1.5) 8/3 Simpson's + two trapezoids
        //cout << "4 steps" << endl;
        int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
        int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
        int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
        int eu_dex = e2_dex;    // index of the last energy point within bandwidth
        double e_l = first_ener*pow(ener_spacing,el_dex); // first energy point in keV
        double e_u = first_ener*pow(ener_spacing,eu_dex); // last energy point in keV
        double h = (log(e_u)-log(e_l)) / 3; // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) / e_l);  // first trapezoid
        flux += h * (McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) + McPHACC4(em1_dex,cos_theta, T, M, R, mexmcc) * 3 + McPHACC4(em2_dex,cos_theta, T, M, R, mexmcc) * 3 + McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc)) * 3 / 8; // Simpson's 3/8, exact             
        flux += (E2 - e_u) / 2 * (McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc) / e_u + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2); // second trapezoid
    }
    if (n_steps == 5){ // five energy points within bandwidth: (4.1.3, 4.1.13) Simpson's "4,2" + two trapezoids  
        //cout << "5 steps" << endl;
        int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
        int em1_dex = e1_dex+2; // index of the second energy point within bandwidth
        int em2_dex = e1_dex+3; // index of the third energy point within bandwidth        
        int em3_dex = e1_dex+4; // index of the fourth energy point within bandwidth        
        int eu_dex = e2_dex;    // index of the last energy point within bandwidth
        double e_l = first_ener*pow(ener_spacing,el_dex); // first energy point in keV
        double e_u = first_ener*pow(ener_spacing,eu_dex); // last energy point in keV
        double h = (log(e_u)-log(e_l)) / 4; // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) / e_l);  // first trapezoid
        flux += h * (McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) + McPHACC4(em1_dex,cos_theta, T, M, R, mexmcc) * 4 + McPHACC4(em2_dex,cos_theta, T, M, R, mexmcc) * 2 + McPHACC4(em3_dex,cos_theta, T, M, R, mexmcc) * 4 + McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc)) / 3; // Simpson's "4,2", exact             
        flux += (E2 - e_u) / 2 * (McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc) / e_u + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2); // second trapezoid
    }
    if (n_steps >= 6){ // six or more energy points within bandwidth: (4.1.3, 4.1.14) Simpson's cubic + two trapezoids
        //cout << "6 steps" << endl;
        int el_dex = e1_dex+1;  // index of the first energy point within bandwidth
        int el1_dex = e1_dex+2; // index of the second energy point within bandwidth
        int el2_dex = e1_dex+3; // index of the third energy point within bandwidth        
        int eu_dex = e2_dex;    // index of the last energy point within bandwidth
        int eu1_dex = e2_dex-1; // index of the second last energy point within bandwidth
        int eu2_dex = e2_dex-2; // index of the third last energy point within bandwidth
        double e_l = first_ener*pow(ener_spacing,el_dex); // first energy point in keV
        double e_u = first_ener*pow(ener_spacing,eu_dex); // last energy point in keV
        double h = (log(e_u)-log(e_l)) / (n_steps-1); // difference between log-spaced energy points
        flux = (e_l - E1) / 2 * (McPHACC3(E1,cos_theta, T, M, R, mexmcc) / E1 + McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) / e_l);  // first trapezoid
        // Simpson's cubic, exact
        flux += h * (McPHACC4(el_dex,cos_theta, T, M, R, mexmcc) * 9 + McPHACC4(el1_dex,cos_theta, T, M, R, mexmcc) * 28 + McPHACC4(el2_dex,cos_theta, T, M, R, mexmcc) * 23) / 24; // first three coefficients
        for (int m = 1; m <= n_steps-6; m++) flux += h * (McPHACC4(el_dex+m+3,cos_theta, T, M, R, mexmcc)); // middle coefficients
        flux += h * (McPHACC4(eu2_dex,cos_theta, T, M, R, mexmcc) * 23 + McPHACC4(eu1_dex,cos_theta, T, M, R, mexmcc) * 28 + McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc) * 9) / 24; // last three coefficients
        // end Simpson's cubic
        flux += (E2 - e_u) / 2 * (McPHACC4(eu_dex,cos_theta, T, M, R, mexmcc) / e_u + McPHACC3(E2,cos_theta, T, M, R, mexmcc) / E2); // second trapezoid                
    }

    flux = flux/Units::H_PLANCK;
    return flux;
}

// Calculate the final interpolated intensity
double McPHACC3(double E, double cos_theta, double T, double M, double R, class LightCurve mexmcc){
	double delta, obl_approx, lgrav, lt, theta, ener_spacing, first_ener, ener_index;
	double e0, e1, th0, th1, grav0, grav1, t0, t1;
	double I_temp[16], I_int[8], J[4], K[2], L(0.0);
	int i_f, n_f, i_lt, i_lgrav, n_lt, n_lgrav, i_mu, n_mu, first_inte;


    //setting values of lt and lgrav based on input T, M, and R. Also sets to load using first mu value.   
    M = Units::nounits_to_cgs(M, Units::MASS);
    R = Units::nounits_to_cgs(R, Units::LENGTH);
    delta = 1 / sqrt(1 - (2 * Units::G * M / (R * Units::C * Units::C)));
    obl_approx = (1 + (-0.791 + 0.776 * mexmcc.para.mass_over_r) * pow(sin(mexmcc.para.omega_bar_sq),2) + (1.138 - 1.431 * mexmcc.para.mass_over_r) * pow(cos(mexmcc.para.omega_bar_sq),2));
    lgrav = log10(delta * Units::G * M / (R * R));
    //lgrav = log10(delta * Units::G * M / R / R * obl_approx);
    lt = log10(1E3 * (T * Units::EV / Units::K_BOLTZ));
    //cout << "temperature in log(K) is " << lt << endl;
    //cout << "gravity in log(cgs units) is " << lgrav << endl;
    //cout << "cos_theta is " << cos_theta << endl;
    //cout << "energy is " << E << endl;

    i_lt = (lt-5.1)/0.05; //if we need to load 1st temperature, i_lt = 0. this is discrete math
    n_lt = i_lt+1;
    t0 = 5.1+0.05*i_lt;
    t1 = t0+0.05;
    //cout << i_lt << " " << n_lt << " " << t0 << " " << t1 << endl;
    //cout << t0 << " " << t1 << endl;

    i_lgrav = (lgrav-13.7)/0.1;
    n_lgrav = i_lgrav+1;
    grav0 = 13.7+0.1*i_lgrav;
    grav1 = grav0+0.1;
    //cout << grav0 << " " << grav1 << endl;
    
    //Find proper mu choice
    n_mu = 1;
    while (cos_theta > mexmcc.mccangl[n_mu] && n_mu < 49){
    	n_mu += 1;
    }

    i_mu = n_mu - 1;
    th0 = acos(mexmcc.mccangl[i_mu+1]);
    th1 = acos(mexmcc.mccangl[n_mu+1]);
    theta = acos(cos_theta);
    //cout << mexmcc.mccangl[i_mu] << " " << mexmcc.mccangl[n_mu] << " " << cos_theta << endl;
    //cout << th0 << " " << th1 << " " << acos(cos_theta) << endl;

    //Find proper freqency choice
    ener_spacing = pow(10.0,0.0338);
	first_ener = 0.004969469;
    ener_index = log(E / first_ener) / log(ener_spacing);
    i_f = (int) ener_index;
    n_f = i_f + 1;
    e0 = first_ener*pow(ener_spacing,i_f);
    e1 = e0*ener_spacing;
    //cout << e0 << " " << e1 << " " << E << endl;


    first_inte = (i_lt*11 + i_lgrav) * 5000 + i_f * 50 + i_mu +1;
    //cout << i_lt << " " << i_lgrav << " " << i_f << " " << i_mu << " " << first_inte << endl;
    I_temp[0] = mexmcc.mccinte[first_inte]*pow(10.0,t0*3.0);
    //double freq = 1E3 * e0 * Units::EV / Units::H_PLANCK;
    //cout << mexmcc.mccinte[first_inte] << " " << pow(10.0,t0*3.0) << " " << t0 << endl;
    //cout << I_temp[0] << " " << first_inte << " " << t0 << " " << grav0 << " " << freq << " " << th0 << " " << mexmcc.mccangl[i_mu+1] << endl;
    //cout << mexmcc.mccangl[0] << " " << mexmcc.mccangl[49] << endl;
    I_temp[1] = mexmcc.mccinte[first_inte+1]*pow(10.0,t0*3.0);
    I_temp[2] = mexmcc.mccinte[first_inte+50]*pow(10.0,t0*3.0);
    I_temp[3] = mexmcc.mccinte[first_inte+51]*pow(10.0,t0*3.0);

    I_temp[4] = mexmcc.mccinte[first_inte+5000]*pow(10.0,t0*3.0);
    I_temp[5] = mexmcc.mccinte[first_inte+5001]*pow(10.0,t0*3.0);
    I_temp[6] = mexmcc.mccinte[first_inte+5050]*pow(10.0,t0*3.0);
    I_temp[7] = mexmcc.mccinte[first_inte+5051]*pow(10.0,t0*3.0);

    I_temp[8] = mexmcc.mccinte[first_inte+55000]*pow(10.0,t1*3.0);
    I_temp[9] = mexmcc.mccinte[first_inte+55001]*pow(10.0,t1*3.0);
    I_temp[10] = mexmcc.mccinte[first_inte+55050]*pow(10.0,t1*3.0);
    I_temp[11] = mexmcc.mccinte[first_inte+55051]*pow(10.0,t1*3.0);

    I_temp[12] = mexmcc.mccinte[first_inte+60000]*pow(10.0,t1*3.0);
    I_temp[14] = mexmcc.mccinte[first_inte+60001]*pow(10.0,t1*3.0);
    I_temp[13] = mexmcc.mccinte[first_inte+60050]*pow(10.0,t1*3.0);
    I_temp[15] = mexmcc.mccinte[first_inte+60051]*pow(10.0,t1*3.0);

    I_int[0] = LogLinear(E, e0, I_temp[0], e1, I_temp[2]); //t0, grav0, th0
    I_int[1] = LogLinear(E, e0, I_temp[1], e1, I_temp[3]); //t0, grav0, th1
    I_int[2] = LogLinear(E, e0, I_temp[4], e1, I_temp[6]); //t0, grav1, th0
    I_int[3] = LogLinear(E, e0, I_temp[5], e1, I_temp[7]); //t0, grav1, th1
    I_int[4] = LogLinear(E, e0, I_temp[8], e1, I_temp[10]);//t1, grav0, th0
    I_int[5] = LogLinear(E, e0, I_temp[9], e1, I_temp[11]);//t1, grav0, th1
    I_int[6] = LogLinear(E, e0, I_temp[12], e1, I_temp[14]);//t1, grav1, th0
    I_int[7] = LogLinear(E, e0, I_temp[13], e1, I_temp[15]);//t1, grav1, th1
    
    //cout << I_int[0] << " " << I_int[1] << " " << I_int[2] << " " << I_int[3] << endl;

    // Interpolate to chosen mu
    J[0] = Linear(theta,th0,I_int[0],th1,I_int[1]); //t0, grav0
    J[1] = Linear(theta,th0,I_int[2],th1,I_int[3]); //t0, grav1
    J[2] = Linear(theta,th0,I_int[4],th1,I_int[5]); //t1, grav0
    J[3] = Linear(theta,th0,I_int[6],th1,I_int[7]); //t1, grav1

    // Interpolate to chosen local gravity
    K[0] = pow(10.0,Linear(lgrav,grav0,log10(J[0]),grav1,log10(J[1]))); //t0
    K[1] = pow(10.0,Linear(lgrav,grav0,log10(J[2]),grav1,log10(J[3]))); //t1

    // Interpolate to chosen temperature
    L = pow(10.0,Linear(lt,t0,log10(K[0]),t1,log10(K[1])));

    // Set to zero at small angle
    if (cos_theta < 0.015629) L = 0;
    //cout << P << endl;

    //cout << I_temp[0] << " " << I_int[0] << " " << J[0] << " " << K[0] << " " << L << endl;
    /*
    if (isnan(L)) {
    	cout << theta << " " << th0 << " " << th1 << " " << i_mu << " " << n_mu << endl;
    	cout << cos_theta << " " << mexmcc.mccangl[i_mu] << " " << mexmcc.mccangl[n_mu] << endl;
    	cout << theta << " " << I_int[0] << " " << J[0] << " " << J[1] << " " << J[2] << " " << J[3] << endl;
    }
    */
    //cout << L << endl;

    return L;
}

