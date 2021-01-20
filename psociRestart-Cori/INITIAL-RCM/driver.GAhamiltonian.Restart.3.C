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
 *   @file driver.GAhamiltonian.Restart.C
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
#include "PsociGAvectors.hpp"
#include "PsociDRAservices.hpp"
#include "PsociDRArestart.hpp"
#include "PsociGAbasis.hpp"
#include "PsociConfigs.hpp"
#include "PsociDeterminants.hpp"
#include "PsociGADeterminants.hpp"
#include "PsociHamiltonian.hpp"
#include "PsociGAhamiltonian.hpp"
#include "PsociGAsubspace.hpp"
#include "PsociConstructHamiltonian.hpp"

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
  
  // Franklin or HOPPER
  //   const unsigned long maxGAMemPerCore = 512 * 1024 * 1024; 
  
  // Blueridge
  const unsigned long maxGAMemPerCore = 1 * 1024 * 1024 * 1024;
  
  /* Regardless of the value of availableMemory() maxMemPerCore is the maximum
     size (per core) available for GA space ( excluding heap and stack
     
     GA_DATA_TYPE  is used by ma_init and only applies to heap and stack
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

// Construct determinants
// Need to fetch a subset of configs, create determinants and push those dets to GA space

  cout << "Build distributed GA_DET space " << endl;
  
  // Pretend the need to do a restart set up as an enumeration
  //Set this up so we can choose to NOT store constructed data to disk
  
   //HSOURCE construct_hamiltonian = RESTART;
   HSOURCE construct_hamiltonian = CONSTRUCT;
  
  
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
/* At this point all determinants should have been generated AND pushed into GAdeterminant space.
   No subsequent locality is required
*/
  
   string filername="moints";
   //   string filername="/home/jtilson/RuO+2-NewBasis-SDCI-PSOCIQA/moints";
   //   string filername="./INTEGRALS/RuOSCI-B1-2.75au/moints";
   
   const int ROOT_READER=0;
   Fint l_unit = 2;
   PsociIntegrals mos( ROOT_READER, l_unit, filername );
   mos.fetchAndProcessIntegrals();

   double fcore =  mos.fetchCoreEnergy();


  int maxsef= deters_det.fetchGlobalMaxSefs();

/* 
  Build basic code to emulate a H construction/restart scenario. If restart simply
  attempt to read the data else try to create H and then store to disk. Add some checking
  so that if one attempts to Create a hamiltonian that already exists we do not by default
  destroy the current H.
*/


/* All Hamiltonian GA space is now handled within the GAhamiltonian method
 */

 
#ifdef SUPPLEMENTAL
// Do not use this

  const int maxsparse=128;
  const int maxsparseBig=5000;
  const int maxwidth = 1 * maxsef; // Around 20% can be allocated for 'big' data sets 
#else
  //const int maxsparse=3000;
    const int maxsparse=3000;
#endif

 cout << "Setting GAH sparsity to " << maxsparse << endl;
 cout << "maxsef is " << maxsef << endl;
 GA::printStats();
 
 /* Allocate DRAs: This is always needed: However we want to eliminate the chance that existing H data
    might get destroyed by accident.
 */  

// First simply start up the service: needed by all
 
  
 PsociDRAservices DRAservice;
 InitializeDRAservice( g_rank, DRAservice );
 
 //Now allocate the DRA arrays: currently as RW: BUt we may want to be more picky.
 //Declare the basis objects used by both restart modes
 
 string CIarrayname = "DRA_cimat";
 string CIfilename = "DRA_cimat";
 PsociDRArestart d_cimat( g_rank, CIarrayname, CIfilename );
 
 string ICIarrayname = "DRA_icicol";
 string ICIfilename = "DRA_icicol";
 PsociDRArestart d_icicol( g_rank, ICIarrayname, ICIfilename );
 
 string NUMarrayname = "DRA_number";
 string NUMfilename = "DRA_number";
 PsociDRArestart d_number( g_rank, NUMarrayname, NUMfilename );
 
 string DIAGarrayname = "DRA_diag";
 string DIAGfilename = "DRA_diag";
 PsociDRArestart d_diag_sef( g_rank, DIAGarrayname, DIAGfilename );
 
 // Now choose mode CONSTRUCT or RESTART ( or NONE )

 if ( construct_hamiltonian == CONSTRUCT ) {
   int ndims = 2;
   cout << "Preparing to Construct Hamiltonian " << endl;
   AllocateDRA( construct_hamiltonian, ndims, maxsef, maxsparse, d_cimat, d_icicol, d_number, d_diag_sef);

 } else if( construct_hamiltonian == RESTART ) { // Plan on fetching hamiltonian   
   cout << "Preparing to READ hamiltonian data from disk " << endl;
   SpecifyDRA(construct_hamiltonian, d_cimat,d_icicol,d_number,d_diag_sef);  
   
 } // if NONE then neither read nor dump 
 

 /* Still need this because the GA space needs to be initialized and the constructor actually
    creates the array space
 */
 
// SUPPLEMENTAL doesn;t yet work with DRA data layers
#ifdef SUPPLEMENTAL
 cout << "Current implementation does not support SUPPLEMENTAL " << endl;  
 GA::Terminate();
 
 cout << "supplemental maxwidth is " << maxwidth << endl;
 cout << "supplemental maxsparseBig is " << maxsparseBig << endl;
 //PsociGAhamiltonian hamilton(maxsparse, maxsef, g_cimat, g_icicol, g_number, g_diag_sef, maxsparseBig, maxwidth, g_cimat_supp, g_icicol_supp, g_number_supp, &deters_det, &mos );
 PsociGAhamiltonian hamilton(maxsparse, maxsef, maxsparseBig, maxwidth, &deters_det, &mos );
#else
 //PsociGAhamiltonian hamilton(maxsparse, maxsef, g_cimat, g_icicol, g_number, g_diag_sef, &deters_det, &mos );
 PsociGAhamiltonian hamilton(maxsparse, maxsef, &deters_det, &mos );
#endif

 // Extra bit neede by all but data comes from Det processing
 
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

 if ( construct_hamiltonian == CONSTRUCT ) {
   cout << "Construct distributed hamiltonian " << endl;
   
   pair<int,double> Hinfo;
   long numH = hamilton.constructGlobalHamiltonian(Hinfo);
   
   cout << "Total number of elements above threshold  " << numH << endl;
   cout << "Timing: Me is " << Hinfo.first << " time is " << Hinfo.second << endl;
   
   if ( construct_hamiltonian != NONE ) {
     cout << " Go to the dump " << endl;     
     pair<int,double> dumpTime;
     SendHamiltonianToDRA(construct_hamiltonian,d_cimat,d_icicol,d_number,d_diag_sef,hamilton,dumpTime);
     cout << dumpTime.first << "reports a dump time of " << dumpTime.second << endl;
   }
   
 } else if( construct_hamiltonian == RESTART ) { // Attempt to read hamiltonian from disk
   
   double timeread = psociTime(); 
   pair<int,double> readTime;
   FetchHamiltonianFromDRA(construct_hamiltonian,d_cimat,d_icicol,d_number,d_diag_sef,hamilton,readTime);
   cout << readTime.first << "reports a read time of " << readTime.second << endl;
 }
 
 // At this point we have the GA data in place We actually do not need DRAS anymore either.
 
 
 cout << "Closing up DRA as we have our data now " << endl;
 ExitDRASERVICE( DRAservice );
 
 
 // Move onto the iterative diagonalization
 // Remembere hamiltonian computes the matrix-vector products

 GA::printStats();
 
 if ( g_rank == 0 ) hamilton.printDiags();
 
 
 // Move onto the iterative diagonalization step.
 // Open up some 1D pretend arrays
 
 cout << "TEST out VECTOR processing " << endl;
 
 ///////////////////////////////////////////////////////////////
 /*
   PsociGArestart() -> create GA space and read vectors - or - guestimate if desired.
   PsociGAbasis ( handle from above )
   PsociGArestart() -> dump current data back to disk ( could have changed size )
   
   1) PsociGArestart()
*/

 int numSoughtRoots = 20; // more like 50 to 100 int maxsubspace = 4 * numSoughtRoots; // 4 expansions per eigenvector (

/* Choose the value for maxsubspace based on preliminary testing with RuO2+, VFCI wavefunctions
   at 40 roots. It was set to be 2* greater than the minimum post-compressed subspace
*/

 int maxsubspace = 4 * numSoughtRoots;
 int minvecs = 2; // in compression a multipland 

 int vdim=2;
 int dims[2],chunk[2];
 dims[0] = maxsef ;// maxsef upwars of 5,000,000;
 dims[1] = numSoughtRoots+maxsubspace;
//  chunk[0] = dims[0]; // Need to do performance tuning on this want to keep vectors together
//  chunk[1] = -1;
 
// Want distribution to be close to that of cimat and friends
 chunk[0] = -1;
 chunk[1] = numSoughtRoots+maxsubspace; // horizontal strips
 
 char * vector_title = "vector";
 char * Hvector_title = "Hvector";
 char * rnorm_title = "rnorm";
 GA::GlobalArray * g_v;
 GA::GlobalArray * g_Hv;
 GA::GlobalArray * g_rnorms;

 g_v  = GA::createGA( C_DBL, vdim, dims, (char *) vector_title, chunk);
 g_Hv = GA::createGA( C_DBL, vdim, dims, (char *) Hvector_title, chunk);
 
 cout << "driver sought is " << numSoughtRoots << endl;

// Only need enough rnorms for the sought vectors - of course one cannot reallocate after this

 dims[0] = maxsef ;// maxsef upwars of 5,000,000;
 dims[1] = numSoughtRoots;

// Want distribution to be close to that of cimat and friends
 g_rnorms = GA::createGA( C_DBL, vdim, dims, (char *) rnorm_title, chunk);
 

 PsociGAvectors gavectors(numSoughtRoots, g_v );
 gavectors.specifyDiagonalGuessVectors( g_v );
 //gavectors.printCurrentSoughtRoots( g_v );

 string vecname="vector_test.txt";
 gavectors.openVectorWrite( GA::nodeid(), vecname );
 gavectors.dumpVector(GA::nodeid(), g_v );
 gavectors.closeVectorFile();
 
 // Try to read thenm back in and see if the diagonal with augmentation stuff works.
 
 gavectors.openVectorRead( GA::nodeid(), vecname );
 cout << gavectors.getNumberVectors();
 cout << gavectors.getNumberVectorsSought();
 gavectors.readVector( GA::nodeid(), g_v );
 
 /* Note in the following one couldf also use PsociGAbasis to augment the system
    followed by manually resetting number of  roots requested usin ghte resizeVectors method
 */
 /*
   cout << "INCREASE VECS" << endl;
   numSoughtRoots += 4;
   gavectors.setNumberVectors(numSoughtRoots);
   gavectors.specifyDiagonalGuessVectors( g_v );
   cout << gavectors.getNumberVectors();
   cout << gavectors.getNumberVectorsSought();
   gavectors.closeVectorFile();
   
   gavectors.openVectorWrite( GA::nodeid(), vecname );
   gavectors.dumpVector(GA::nodeid(), g_v );
   gavectors.closeVectorFile();
 */
 
 //gavectors.printCurrentSoughtRoots( g_v );
 
 // Note I do not need gavectors anymore until the end when I want to dumpo solution vectors
 // works fine
 // Take set of vectors and incorporate intoPsociGAbasis
 
 /*
   2) PsociGAbasis
   
 */
 
  PsociGAhamiltonian * CIhandle = &hamilton;
  PsociGAbasis vector_set( numSoughtRoots, CIhandle , g_v, g_Hv );
  vector<double> vecnorm;
  //vector_set.fetchHandle()->print();
  
  vector_set.getAllNorms( vecnorm );
  GA::sync();
  for (int i=0; i<numSoughtRoots; ++i) {
    cout << "rank is " << GA::nodeid() << " VECNORM data are " << vecnorm[i] << endl;
  }
  cout << "Now normalize vectors and print norms again " << endl;
  pair<int,double> time;
  vector<double> vecnorm2;
  vector_set.getAllNorms( vecnorm2, time);
  cout << "Times rank is " << time.first << " Norm time is " << time.second << endl;
  
  vector_set.normAllVectors();
  vector_set.getAllNorms( vecnorm );
  
  int or_status = vector_set.orthonormalizeAllVectors();
  if ( or_status != 0 ) {
    cout << "A linear dependancy was found at i " << or_status << endl;
    cout << "This currently is a hard failure: Aborting " << endl;
    GA::Terminate();
  }
  GCOUT << "VECTORS are orthonormal "<< endl;
  
  GCOUT << GA::nodeid() << " Sought roots " << vector_set.CurrentNumSoughtRoots() << endl;
  GCOUT << GA::nodeid() << " Existing roots " << vector_set.CurrentExistingRoots() << endl;
  GCOUT << GA::nodeid() << " Max roots " << vector_set.CurrentMaxRoots() << endl;
  GCOUT << GA::nodeid() << " Req products " << vector_set.CurrentReqHVProducts();
  GCOUT << GA::nodeid() << " Existing products " << vector_set.CurrentExistingHVProducts();
  
  
  //Need to compute some H*v products to test subspace code
  
  int num = vector_set.computeHvProducts( 0, numSoughtRoots );
  GCOUT << "Computed " << num << " H*v products " << endl;
  
  //////////////////////////////////////////////////////////////////////////////
  /* 
     Take this current set of vectors and apply them to PsociGAsubspace
     It is there where the real subspace work occures.
  */
  
  //PsociGAbasis does the preconditioning and H*v work
  
  GCOUT << "INITIALIZE the SUBSPACE " << endl;
  
/* Create an rnorm matrix 
*/

// Try to converge the lowest root
   
// This "setpresent" call is here because we PRECOMPUTED some terms as partof toolkit testing: NOT DOING THIS ADDS A FW h*V TO THE CALC

// IMPORTANT
  vector_set.setpresent_matrixProducts( numSoughtRoots ); // tell the code some were already computed
  
  vector<double> evals( numSoughtRoots+maxsubspace ); // Need these for preconditioning - level shifting
  PsociGAbasis * sub_vector_set = &vector_set;
  
  PsociGAsubspace subspace( sub_vector_set, hamilton.fetchDIAG_SEF(), g_rnorms  );
  
  
  GCOUT << "maxsubspace is " << subspace.fetchMaxSubspace() << endl;
  GCOUT << "existing subspace " << subspace.fetchExistingSubspace() << endl;
  GCOUT << "subspace CI rows is " << subspace.numCIBasis() << endl;
  GCOUT << "subspace Vec rows is " << subspace.numVectorBasis() << endl;
  GCOUT << "subspace Vec cols is " << subspace.numRoots() << endl;
  
  bool converged=false;
  vector<double> rnorms; // Replicated list of total rnorms
  int count=0;
  int total_count=0;
  int numcont=0;
  
  int curroot=1;
   
// What is the right value for here?
  double rtol = MIN_PRECISION; // From PsociGAbasis: Tighter than this causes linear dependancies
   
   /* The current algorithm consists of ensuring ALL sought vectors are converged
      before declaring success.
   */
   
  vector<bool> vecs_converged( vector_set.CurrentNumSoughtRoots(), false ); 
  
// Converge the current set of vectors

  pair<int,double> diagTime;
  double diag = psociTime();
  
  double coarse_rtol = 0.001; //Initial coarse first-pass
  
/* For this current approach we first attempt to converge to a coarse levle of tolerance
   then repeat at the higher level of tolerance. Clearly much better ways are possible
   but that will need to come later

   At the start of the tighter tolerance cyclei we use as starting vectors the current set
   of converged vectors by calling the resize method. This reorthogonalizes the desired set
   and recomputes all H*v terms
*/

/* For now we do not have a decent way to reduce the space so only do one iteration for now
*/
  //int max_iterations_per_rtol = numSoughtRoots+maxsubspace ;
  int max_iterations_per_rtol = 1000; // picked out of a hat ! 
  double current_rtol = coarse_rtol;
  //double current_rtol = rtol;

  int subcount;
  count=0;
  converged=false;

  subcount = 0;
  double sumDiagTime = 0.0;
  double tempTime;

  int compressSubspace;
  int matrixVectorProducts=0;

  while( count < max_iterations_per_rtol && converged != true) {
    
    ++total_count;
    ++count;
    ++subcount;
    
    tempTime = psociTime();

/* We don't need this anymore
    subspace.generateSubspace();
    matrixVectorProducts += vector_set.CurrentExistingRoots();
*/
    if ( subcount == 1 ) { //trigger on subcount because we update every contraction
      subspace.generateSubspace();
      matrixVectorProducts += vector_set.CurrentExistingRoots();
      //subspace.generateIncrementalSubspace( vector_set.CurrentExistingRoots() ); // need this to build basis subspace
    } else {
      subspace.generateIncrementalSubspace( 1 );
      matrixVectorProducts += 1;
    }

    subspace.diagSubspace();
    subspace.computeRnorms( rnorms );
    //subspace.computeRnormsGA( rnorms );
    subspace.brdcstEvals( evals );
   
    sumDiagTime += psociTime() - tempTime;
    
    //GCOUT << count << " current SUBSPACE size is " << vector_set.CurrentExistingRoots() << endl;
    //GCOUT << "Maxcount and count " << maxsubspace<< " " << subcount << endl;
    if ( GA::nodeid() == 0 ) {
      cout << "Current roots (au) : rnorm" << endl;
      for(int i=0; i< numSoughtRoots; ++i ) {
        cout << evals[i] << " : " << rnorms[i] << endl;
      }
      subspace.printVectorContents( rnorms );
    }
    
    converged = true; //hope for the best
    
    for ( int i=0; i < numSoughtRoots; ++i ) { //Check all rnorms
      if ( rnorms[i] > current_rtol ) {
	curroot = i+1; // Not the best way: Start at the bottom and do one at a time
	converged = false;
	break;
      }
    }
    
    if( converged == false ) {
      double shift = evals[ curroot-1 ]; // Try to update vector == curroot 
      //GCOUT << "curroot is " << curroot << endl;
      int nump = vector_set.preconditionAndNewVectorIncore(curroot,shift, g_rnorms, hamilton.fetchDIAG_SEF() );
      
      // at this point g_scratch is populated and ready for appending.
      
      vector_set.appendScratchVector();
    }
    
    if ( converged == true && current_rtol == coarse_rtol ) {
      //GCOUT << "switch and compress to tight tolerance " << endl;
      converged = false;
      subcount = 0;
      current_rtol = rtol;
      curroot = 1; // Start over
      vector_set.expandAndReplaceVectors( subspace.fetchEvecsHandle(), minvecs *numSoughtRoots );
      ++compressSubspace;
      GCOUT << count << " NEW SUBSPACE size is " << vector_set.CurrentExistingRoots() << endl;
    }


// Expand current Evectors into g_vector ( keep subspace at least 3 * vectors in size )
    if ( vector_set.CurrentExistingRoots() > maxsubspace ) {
      cout << "COMPRESSING" << vector_set.CurrentExistingRoots()<<" "<< maxsubspace << endl;
      vector_set.expandAndReplaceVectors( subspace.fetchEvecsHandle(), minvecs *numSoughtRoots ); 
      ++compressSubspace;
      subcount = 0; // forces a complete rebuild of the H*v matrix
    }

/* Keeping 2 or 3 expansion vectors per sought vector was based on some initial testing
   using the VFCI RuO2+ configuration but is generally consistent with

   M. Crouzeix, B. Philippe, M. Sadkane, "The Davidson Method" 
   SIAM J. Sci. Comput. Vol 15, No. 1 pp 62-76, January 1994.

and also
 
   J.H. van Lenthe, P. Pulay, "A Space Saving Modification of Davidson's Eigenvector ALgorithm",
   J. Comp. Chem. Vol 11, No. 10, pp 1164-1168, 1990.
*/


  } // converged and/or out of cycles
  
  diagTime.second =  psociTime() - diag;
  diagTime.first = GA::nodeid();
  
  cout << GA::nodeid() <<  " sumDiagTime is " << sumDiagTime << endl;
  if ( GA::nodeid()== 0 ) {
    GCOUT << "Number of SUBSPACE contractions is " << compressSubspace << endl;
    GCOUT << "Number of total iterations is " << count << endl;
    GCOUT << "TOTAL MATRIX VECTOR PRODUCTS is " << matrixVectorProducts << endl;
    GCOUT << "Current results: fcore is " << fcore << " " << endl;
    for(int i=0; i< vector_set.CurrentNumSoughtRoots(); ++i ) {
      GCOUT << "Eigenvalue " << i << " Total energy (au) " << evals[i]+fcore <<  " Electronic (au) " << evals[i] << " rnorm " << rnorms[i] << endl;
    }
  }
  
  cout << "GA rank " << diagTime.first << " Diagonalization time (sec) is " << diagTime.second << endl;
  cout << "At end of iteration " << total_count << " with all converged = " << converged << endl;
  
  GCOUT << "Total cycles is " << total_count << endl;
  GCOUT << "Number of subspace contractions is " << numcont << endl;
  
  subspace.destroySubspace();
  
  
  // End of the test
  
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



