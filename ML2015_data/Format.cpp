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
  std::ifstream in;       // input stream

  char out_file[256] = "output.txt";
  char in_file[256] = "No file specified";
		
  int numbins(16), numbands(30);

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
	                sscanf(argv[i+1], "%s", in_file);
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

      in.open(in_file);

      char line[265]; // line of the data file being read in
      double get_t, get_e, get_f;

      for (unsigned i(0);i<numbands;i++){ // loop through the energy bands

	for (unsigned j(0);j<numbins;j++){ // loop through the time bins

	  in.getline(line,265);  
	  sscanf( line, "%lf %lf %lf", &get_t, &get_e, &get_f);

	  if (i==0)
	    time[j] = get_t;
	  
	  flux[i][j] = get_f;
	  err[i][j] = sqrt(get_f);
	}
	energy[i] = get_e;
      }

      in.close();
    
      out.open(out_file);

      for (unsigned j(0);j<numbins;j++){ // loop through the time bins

	out << time[j] << " " ;

	for (unsigned i(0);i<numbands;i++){ // loop through the energy bands

	  out << flux[i][j] << " "
	      << err[i][j] << " ";
	   
	}
	out << std::endl;

      }

      out.close();

    return 0;
 } 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
