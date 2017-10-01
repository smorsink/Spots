/***************************************************************************************/
/*                                   Format.cpp

This code reads in a data file and changes it into our format.

The program needs to know:
(1) number of time bins
(2) number of energy bins

*/
/***************************************************************************************/

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <exception>
#include <vector>
#include <string>
#include <string.h>

// MAIN
int main ( int argc, char** argv ) try {  // argc, number of cmd line args; 
                                          // argv, actual character strings of cmd line args

  /*********************************************/
  /* VARIABLE DECLARATIONS AND INITIALIZATIONS */
  /*********************************************/
    
  std::ofstream out;      // output stream; printing information to the output file
  std::ifstream in1;       // input stream
  std::ifstream in2, in3;

  char out_file[256] = "output.txt";
  char in_file1[256] = "No file specified";
  char in_file2[256] = "No file specified";
  char in_file3[256] = "No file Specified";
		
  int numbins(16), numbands(300);

  /*********************************************************/
  /* READING IN PARAMETERS FROM THE COMMAND LINE ARGUMENTS */
  /*********************************************************/
    
    for ( int i(1); i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {

	    case 't': // Number of time bins
	      sscanf(argv[i+1], "%u", &numbins);
	      break;

	    case 'i':  // Name of input file
	                sscanf(argv[i+1], "%s", in_file1);
	                break;

	    case 'j':  // Name of input file
	                sscanf(argv[i+1], "%s", in_file2);
	                break;

	    case 'k':  // Name of input file
	                sscanf(argv[i+1], "%s", in_file3);
	                break;

	    case 'e': // Number of energy bins
	      sscanf(argv[i+1], "%u", &numbands);
	      break;

	    case 'o':  // Name of output file
	                sscanf(argv[i+1], "%s", out_file);
	                break;
	                
            } // end switch	
        } // end if
    } // end for

    double time[numbins];
    double energy[numbands];
    double flux[numbands][numbins];
    double err[numbands][numbins];

    double diff;

      in1.open(in_file1);
      in2.open(in_file2);
      in3.open(in_file3);
      out.open(out_file);
      out.precision(10);

      char line[265]; // line of the data file being read in
      double get_t, get_e, get_f1, get_f2, get_f3, get_eem, get_et;
      double diffchisquare, diffloglikely, loglike3, loglike2;

      for (unsigned i(0);i<numbands;i++){ // loop through the energy bands

	for (unsigned j(0);j<numbins;j++){ // loop through the time bins

	  in1.getline(line,265);  
	  sscanf( line, "%lf %lf %lf", &get_t, &get_e, &get_f1 );

	  in2.getline(line,265);  
	  sscanf( line, "%lf %lf %lf", &get_t, &get_e, &get_f2);

	  in3.getline(line,265);  
	  sscanf( line, "%lf %lf %lf", &get_t, &get_e, &get_f3);
	  
	  if (i==0 && j==0)
	    std::cout << "f1 = " << get_f1
		      << " f2 = " << get_f2
		      << " f3 = " << get_f3
		      << std::endl;


	  diff = (get_f1-get_f2)/get_f1 * 100;

	  if(std::isnan(diff) )
	    diff= 0.0;

	  diffchisquare += get_f1 * ( pow( 1.0 - get_f2/get_f1,2) - pow( 1.0 - get_f3/get_f1,2) ) ;

	  loglike2 += get_f1 *  (log(get_f2)) - get_f2;

	  loglike3 += get_f1 *  (log(get_f3)) - get_f3;

	  diffloglikely += get_f1 *  (log(get_f3/get_f2)) - get_f3 + get_f2 ;

	
	  
	    out << get_t << " " 
		<< get_e << " "
		<< diff << " "
		<< get_f1 << " "
		<< get_f2 << " "
	      //<< get_eem << " "
	      //<< get_et
		<< std::endl;


	}
	energy[i] = get_e;
      }


      std::cout.precision(10);

      //std::cout << "Chi^2 = " << chisquare << std::endl;
      std::cout << "Chi^2/2 = " << 0.5*diffchisquare << std::endl;
      //std::cout << "dof = " << numbands * numbins << std::endl;
      std::cout << "logL2 = " << loglike2 << std::endl;
      std::cout << "logL3 = " << loglike3 << std::endl;

      std::cout << "Delta Log(L) = " << diffloglikely << std::endl;


      in1.close();
      in2.close();
      in3.close();
      out.close();

    return 0;
 } 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
