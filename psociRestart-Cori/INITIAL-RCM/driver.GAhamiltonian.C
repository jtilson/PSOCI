/**************************************************************************************
 * Copyright (c) 2010,2011 RENCI.
 * All rights reserved. This program and the accompanying materials
 * MAY BE available under the terms of the RENCI Open Source License
 * UNC at Chapel Hill which accompanies this distribution, and is available at
 * http://www.renci.org/resources/open-source-software-license

 * New implementation of PSOCI:

 Classes: 

 Description: 
 
 History:

**************************************************************************************/
/**
 *   @file driver.GAhamiltonian.C
 *
 */
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include <cmath>
// Now start testing for GA
#include <ga++.h>

#include <mpi.h>

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
#include "PsociGADeterminants.hpp"
#include "PsociHamiltonian.hpp"
#include "PsociGAhamiltonian.hpp"


// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//
//  Start driver program 
//
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#define GCOUT if (g_rank==0) cout

// New declarations
long driverFetchConfigurations( PsociConfigs & configs );

int main(int argc, char **argv)
{
  // Set up global MPI fabric
  
  int g_size;
  int g_rank;
  
// Set up Global Array fabric
// heap and stack are per-core quantities - for collective operations
// ALlocate in terms of DOUBLES for now

// Set them up 4/2 MegaWords each
  const unsigned long heap = 4 * 1024 * 1024, stack = 2 * 1024 * 1024 ;

// Generally at most 25% memory/core. We want to leave 50% for O/S activities and 25% for C++ dynamic container allocations

// Franklin and HOPPER at NERSC
//   const unsigned long maxGAMemPerCore = 512 * 1024 * 1024; 

// Blueridge
 const unsigned long maxGAMemPerCore = 1 * 1024 * 1024 * 1024;

/* Regardless of the value of availableMemory() maxMemPerCore is the maximum
   size (per core) available for GA space ( excluding heap and stack )

   GA_DATA_TYPE  is used by ma_init and opnly applies to heap and stack
*/

  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, maxGAMemPerCore );

  g_size = GA::nodes();
  g_rank = GA::nodeid();
  if ( GA::usesMA() ) cout << "GA memory is coming from MA " << endl;
  if ( GA::usesFAPI() ) cout << "GA is using Fortran indexing " << endl;

  GCOUT << "heap and stack are" << heap << " " << stack << endl;
  GCOUT << "Specified size of " << maxGAMemPerCore << endl;
  GCOUT << "GA usesMA ? " << GA::usesMA() << endl;
  GCOUT << "GA memory limited ? " << GA::memoryLimited() << endl;
  GCOUT << "GA inquireMemory is " << GA::inquireMemory() << endl;
  GCOUT << "May not be GA array limit: GA memory avail " << GA::memoryAvailable() << endl;
  
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  GCOUT << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  GCOUT << "Parallel GA Run with " << g_size << " Total processors " << endl;
  
// Perform a restart like configs read
  string filename = "configs.test";
  PsociConfigs configs( filename ); //File is opened/closed at the read
  long configBytes = driverFetchConfigurations( configs );
  cout << "Current per-core Configs memory (MB) " << configBytes/1000000 << endl;
  cout << endl;

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//Construct determinants
// Need to fetch a subset of configs, create determinants and push those dets to GA space

  cout << "Build distributed GA_DET space " << endl;

  GA::GlobalArray * g_det; // declare global det data 
  GA::GlobalArray * g_det_num;
  GA::GlobalArray * g_nsef;

  PsociGADeterminants deters_det( g_det, g_det_num, g_nsef );
  deters_det.generateDistributedDeterminants( configs );

  cout << "Before nsef " << deters_det.fetchGlobalMaxSefs() << endl;
  vector<int> l_nsef;
  deters_det.assembleAndBrdcstNsef( l_nsef ); //don't need l_nsef anymore
  cout << "after nsef " << deters_det.fetchGlobalMaxSefs() << endl;

//Collective call
/* At this point alll determinants should have been generated AND pushed into GAdeterminant space.
   No subsequent locality is required
*/
   
   string filername="moints";
//   string filername="/home/jtilson/RuO+2-NewBasis-SDCI-PSOCIQA/moints";
//   string filername="./INTEGRALS/RuOSCI-B1-2.75au/moints";

   const int ROOT_READER=0;
   Fint l_unit = 2;
   PsociIntegrals mos( ROOT_READER, l_unit, filername );
   mos.fetchAndProcessIntegrals();

// Hamiltonian GA space is now wholely allocated and maintained within the GAhamiltonian class



// Need to differentiate between construction and restarts

  int maxsef= deters_det.fetchGlobalMaxSefs();

#ifdef SUPPLEMENTAL
  const int maxsparse=128;
  const int maxsparseBig=5000;
  const int maxwidth = 1 * maxsef; // Around 20% can be allocated for 'big' data sets 
#else
  const int maxsparse=5000;
#endif


  cout << "Setting GAH sparsity to " << maxsparse << endl;
  cout << "maxsef is " << maxsef << endl;
  GA::printStats();

  
#ifdef SUPPLEMENTAL
  cout << "supplemental maxwidth is " << maxwidth << endl;
  cout << "supplemental maxsparseBig is " << maxsparseBig << endl;
//  PsociGAhamiltonian hamilton(maxsparse, maxsef, g_cimat, g_icicol, g_number, g_diag_sef, maxsparseBig, maxwidth, g_cimat_supp, g_icicol_supp, g_number_supp, &deters_det, &mos );
  PsociGAhamiltonian hamilton(maxsparse, maxsef, maxsparseBig, maxwidth, &deters_det, &mos );
#else
  //PsociGAhamiltonian hamilton(maxsparse, maxsef, g_cimat, g_icicol, g_number, g_diag_sef, &deters_det, &mos );
  PsociGAhamiltonian hamilton(maxsparse, maxsef, &deters_det, &mos );
#endif


// Need this for the hamiltonian  but it comes from det processing and must be included on restarts

  int nsize = hamilton.computeNsefEntryPoints( l_nsef );

  //hamilton.printNsefEntryPoints();
/*
  cout << "nsize is " << nsize << endl;
  if ( GA::nodeid() == 0 ) {
  cout << "GLobal l_nsef ENTRY POINTS " << l_nsef.size() << endl;
  for(int i=0; i<deters.fetchGlobalMaxSpatials(); ++i ) {
     cout << GA::nodeid() << " " << l_nsef[i];
  }
  cout << endl;
  }
*/

  hamilton.printKsym();
  hamilton.printPerSpatialValues();
  hamilton.hamiltonianDistribution();

// Compute H
/* Not used anymore  applies only to the fetch2D code
  hamilton.SetGranularity( 1 );
  hamilton.PrintGranularity();
*/

  pair<int,double> Hinfo;
  long numH = hamilton.constructGlobalHamiltonian(Hinfo);
  //long numH = hamilton.constructGlobalHamiltonianSparseStats();

  cout << "Total number of elements above threshold  " << numH << endl;
  cout << "Timing: Me is " << Hinfo.first << " time is " << Hinfo.second << endl;

  GA::printStats();


  if ( g_rank == 0 ) hamilton.printDiags();

// Open up some 1D pretend arrays


/*
  cout  << "Print the CIMAT " << endl;
  //hamilton.dumpCimat();

*/

  cout << " I am " << GA::nodeid() << " at the final SYNC " << endl;
  GA::sync();
  GA::Terminate();
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX finished XXXXXXXXXXXXXXXXXXX " << endl;
  
}

// Fetch vectors from Disk if a restart is called.
// A collective call
// Returns current estimate of per-core memory use. Can be torn down befor H construction
long driverFetchConfigurations( PsociConfigs & configs )
{
  int g_rank = GA::nodeid();

  if ( g_rank == 0 ) {
    configs.printFilename();
    if ( configs.readConfigs() != 0 ) {
       cerr << "Failed to read configs at numfgs =" << endl;
       GA::Terminate();
    }
    configs.printParams();
    configs.printConfigurations(); 
  }
  configs.brdcstConfigs( 0 ); //now replicated to all nodes

  int maxspatials;
  int nbasis;
  vector<pair<int,string> > word;
  maxspatials = configs.numTotalSpatials();
  nbasis = configs.numTotalSpatials();
  cout << "Num total spatials is " << maxspatials << endl;
  return( maxspatials*nbasis*sizeof(int) );
}



