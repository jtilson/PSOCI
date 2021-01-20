/**************************************************************************************
 * Copyright (c) 2010 RENCI.
 * All rights reserved. This program and the accompanying materials
 * MAY BE available under the terms of the RENCI Open Source License
 * UNC at Chapel Hill which accompanies this distribution, and is available at
 * http://www.renci.org/resources/open-source-software-license

 * New implementation of PSOCI:

 Classes: 

 Description: 

 This code tests access to the Integral layer provided by the Columbus SIFS code
 
 History:

**************************************************************************************/
/**
 *   @file driverTestVector.C
 *
 */
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include <cmath>

#include <getopt.h>

// Now start testing for GA
#include <ga++.h>

#include <dra.h>
#define GA_DATA_TYPE C_DBL

//End ga++

#include "PsociTimer.hpp"
#include "PsociIntegrals.hpp"
#include "driver.integrals.hpp"

using namespace std;

void parseUserValues(int argc, char **argv)
{ 
  int c;
  while (1)
    {
      static struct option long_options[] =
        {
          {"defaults",            no_argument,       0, 'd'},
          {"filename",            required_argument, 0, 'f'},
          {"help",                no_argument,       0, 'h'},
          {0, 0, 0, 0}
        };
      int option_index = 0;
      c = getopt_long (argc, argv, "dhf:",
                       long_options, &option_index);

      if (c == -1)
        break;

      switch (c)
        {
        case 0:
          if (long_options[option_index].flag != 0)
            break;
          if (optarg)
          cout << endl;
          break;

        case 'd':
          cout << "option -d: Current list of argument values " << endl;
          cout << "root_filename \t\t" << root_filename << endl;
          exit(0);
          break;

        case 'f':
          root_filename = optarg;
          break;

        case 'h':
          cout << "Usage: driver [options] " << endl;
          cout << "option -d: Current list of argument values " << endl;
          cout << " --filename \t\t <string> filename of MOINTS data  " << endl;
          cout << " --defaults \t\t Report current settings of argument variables " << endl;
          cout << " --version \t\t Report version of the Slammr library " << endl;
          cout << " --help \t\t This message" << endl;
          exit(0);
          break;

        default:
          abort ();
        }
    }
  return;
}


// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//
//  Start driver program 
//
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#define GCOUT if (g_rank==0) cout

int main(int argc, char **argv)
{
  parseUserValues(argc, argv);

  // Set up global MPI fabric
  
  int g_size;
  int g_rank;
  
  // Set up Global Array fabric
  // heap and stack are per-core quantities - for collective operations
  // ALlocate in terms of DOUBLES for now

  const long wordSize = sizeof(double);
  const long OneGigaByte = 128*1024*1024 * wordSize; // A total of 1 Gig but Generally thinking interms of doubles (words)... 
  long maxMemPerCore = OneGigaByte;

  unsigned int heap=9000000, stack=9000000;
  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, 0);
  GA::setMemoryLimit( maxMemPerCore );
  g_size = GA::nodes();
  g_rank = GA::nodeid();
  if ( GA::usesMA() ) cout << "GA memory is coming from MA " << endl;
  if ( GA::usesFAPI() ) cout << "GA is using Fortran indexing " << endl;
  
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  GCOUT << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  
  GCOUT << "Parallel GA Run with " << g_size << " Total processors " << endl;
  
  string filename = root_filename;

// Open up MOINTS
   
  std::cout << "MAIN: Try to read header data out of a COlumbus file " << std::endl;

   Fint l_unit = 2;
   PsociIntegrals mos( 0, l_unit, filename );
   mos.printFilename();
   Fint unit = mos.fetchUnit();
   mos.OpenFile();

   if ( mos.sifrh1() != 0 ) {
     cerr << "sifrh1 not zero value aborting: " << endl;
     cerr << "probably incompatible 32bit/64bit builds " << endl;
     //GA::Terminate(); 
     //exit(1);
   }

   mos.printSifrh1();
   cout << endl << endl;

   int status = mos.sifrh2();
   mos.printSifrh2();

   status = mos.fetchOneElectronInts();
   int stat = mos.brdcstOneElectronInts();

   mos.set1ePrint();
   mos.printOneElectronInts();

   mos.fetchTwoElectronInts();
   int statb = mos.brdcstTwoElectronInts();

//   mos.set2ePrint();
   mos.printTwoElectronInts();

   GA::sync();


// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


  GA::Terminate();
  
}



