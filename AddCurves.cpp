/***************************************************************************************/
/*                                   AddCurves.cpp

This code reads in 2 light curves and adds (or subtracts) them together

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
  std::ifstream in2;

  char out_file[256] = "output.txt";
  char in_file1[256] = "No file specified";
  char in_file2[256] = "No file specified";

		
  int numbins(128), numbands(1);

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
      out.open(out_file);
      out.precision(10);

      char line[265]; // line of the data file being read in
      double get_f1, get_f2, get_eem, get_et;
      double chisquare, loglikely;
      double get_t, get_e;
      int channel;

      for (unsigned i(0);i<numbands;i++){ // loop through the energy bands

	for (unsigned j(0);j<numbins;j++){ // loop through the time bins

	  in1.getline(line,265);  
	  sscanf( line, "%lf %lf %lf %d", &get_e, &get_t, &get_f1, &channel );
	  //sscanf( line, "%lf %lf", &get_t, &get_f1 );

	  if (i == 0 && j == 0)
	    std::cout << "line1 = " << line << std::endl;

	  in2.getline(line,265);
	  sscanf( line, "%lf %lf %lf %d", &get_e, &get_t, &get_f2, &channel );
		  
	  //sscanf( line, "%lf %lf", &get_t, &get_f2);

	  if (i == 0 && j == 0)
	    std::cout << "line2 = " << line << std::endl;


	  //if ( i==0)

	  diff = (get_f1-get_f2);

	  if(std::isnan(diff) || diff < 0.0 )
	    diff= 0.0;

	 
	  if (i == 0 && j == 0)
	    std::cout << "Diff = " << diff << std::endl;


	  
	    out << get_e << " " 
		<< get_t << " "
		<< diff << " "
		<< get_f1 << " "
		<< get_f2 << " "
	      //<< get_eem << " "
	      //<< get_et
		<< std::endl;


	}
	energy[i] = get_e;
      }


 

      in1.close();
      in2.close();
      out.close();

    return 0;
 } 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
