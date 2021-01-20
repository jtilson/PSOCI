/**************************************************************************************
 * Copyright (c) 2010 RENCI.
 * All rights reserved. This program and the accompanying materials
 * MAY BE available under the terms of the RENCI Open Source License
 * UNC at Chapel Hill which accompanies this distribution, and is available at
 * http://www.renci.org/resources/open-source-software-license

 * New implementation of PSOCI:

 Classes: 

 Description: 

 This code only tests the determinant generation steps. Vector restarts, subspace solutions
 disk-resident arrays, hamiltonian construcxtionand utils are in other tests
 
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
#include "PsociVector.hpp"
#include "PsociGArestart.hpp"
#include "PsociDRAservices.hpp"
#include "PsociDRArestart.hpp"
#include "PsociGAbasis.hpp"
#include "PsociConfigs.hpp"
#include "PsociDeterminants.hpp"

#include "driver.configs.hpp"


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
          cout << " --filename \t\t <string> filename of spatial configuration data  " << endl;
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
  PsociConfigs configs( filename ); //File is opened/closed at the read

  if ( g_rank == 0 ) {
    configs.printFilename();
    if ( configs.readConfigs() != 0 ) {
       cerr << "Failed to read configs at numfgs =" << endl;
       GA::Terminate();
    }
    configs.printParams();
    configs.printConfigurations();
  }
  configs.brdcstConfigs( 0 );

  if ( g_rank == 0 ) {
  vector<pair<int,string> > word;
  int maxspatials = configs.numTotalSpatials();
  cout << "Num total spatials is " << maxspatials << endl;
  int newsym = configs.fetchConfig( 1, maxspatials, word );
  cout << "new sym and string are " << newsym << endl;
  }

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


  PsociDeterminants deters( &configs );
  deters.printParams();
  int maxspatials = deters.maxSpatials(); 

// Process the Configs in chunks per processor
  pair<int,double> time;
  
  if ( g_rank == 0 ) {
  deters.fetchConfigs( 1, maxspatials, time );
  }

  cout << "I am " << time.first << " TIme to fetch cnfigs is " << time.second << endl;
 
 deters.printSpatialsData( 1 );
 deters.printSpatialsData( 2 );
 //deters.printSpatialsData( 1000 );

 deters.assembleGlobalMaxDetSef();

 cout << "local max ndeti per spatial is " << deters.fetchLocalMaxDetPerSpatial() << endl;
 cout << "global max ndeti per spatial is " << deters.fetchGlobalMaxDetPerSpatial() << endl;
 cout << "local max nseti per spatial is " << deters.fetchLocalMaxSefPerSpatial() << endl;
 cout << "global max nseti per spatial is " << deters.fetchGlobalMaxSefPerSpatial() << endl;


 cout << "Print spin parity is " << deters.fetchSpinParity() << endl;
 deters.assembleGlobalNumSpatials();
 cout << "num spatials readin by PsociConfigs method is " << deters.maxSpatials() << endl;
 cout << "Max num actually processed globally are " << deters.fetchGlobalMaxSpatials() << endl;
 cout << "Max ndet processed globally are " << deters.fetchGlobalMaxDet() << endl;
 cout << "Max nsef processed globally are " << deters.fetchGlobalMaxSef() << endl;




 deters.tearDownDeterminants();




//Generally do not want eveyone print out their arrays.....
//but for now go ahead

  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX finished XXXXXXXXXXXXXXXXXXX " << endl;
  
  GA::Terminate();
  
}



