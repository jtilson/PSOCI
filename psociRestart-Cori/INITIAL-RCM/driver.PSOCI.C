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

#include <getopt.h>

#include <cmath>
#include <ga++.h>

#include <mpi.h>

#include <dra.h>
#define GA_DATA_TYPE C_DBL

#define GCOUT if (g_rank==0) cout

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
#include "PsociCIAnalysis.hpp"

#include "driver.PSOCI.hpp"

#ifdef RCRMONITOR
extern "C" {
#include "blackboard.h"
}
#endif


void parseUserValues(int argc, char **argv)
{
  int c;
  while (1)
    {
      static struct option long_options[] =
        {
          {"defaults (d)",            no_argument,       0, 'd'},
          {"numRoots (n)",            required_argument, 0, 'n'},
          {"moints_filename (m)",     required_argument, 0, 'm'},
          {"configs_filename (c)",    required_argument, 0, 'c'},
          {"vector_filename (v)",     required_argument, 0, 'v'},
          {"vector_restart (r)",      no_argument,       0, 'r'},
          {"hamiltonian_restart (H)", no_argument,       0, 'H'},
          {"hamiltonian_store (W)",   no_argument,       0, 'W'},
          {"subspace_min_mult (i)",   required_argument, 0, 'i'},
          {"subspace_max_mult (x)",   required_argument, 0, 'x'},
          {"subspace_iterations (I)", required_argument, 0, 'I'},
          {"coarse_rtol (C)",         required_argument, 0, 'C'},
          {"final_rtol (F)",          required_argument, 0, 'F'},
          {"maxsparse (s)",           required_argument, 0, 's'},
          {"nxtval_chunk (N)",         required_argument, 0, 'N'},
          {"showevals (S)",           no_argument,       0, 'S'},
          {"version (V)",             no_argument,       0, 'V'},
          {"help (h)",                no_argument,       0, 'h'},
          {0, 0, 0, 0}
        };
      int option_index = 0;
      c = getopt_long (argc, argv, "n:s:dm:c:v:ri:a:I:C:F:N:hHVSW:",
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
          cout << "moints_filename \t " << moints_filename << endl;
          cout << "vector_filename \t " << vector_filename << endl;
          cout << "configs_filename \t " << configs_filename << endl;
          cout << "subspace_min_mult \t " << subspace_min_mult << endl;
          cout << "subspace_max_mult \t " << subspace_max_mult << endl;
          cout << "subspace_iterations \t " << subspace_iterations << endl;
	  cout << "vector_restart \t\t " <<vector_restart << endl;
	  cout << "hamiltonian_restart \t " <<hamiltonian_restart << endl;
          cout << "hamiltonian_store \t " <<hamiltonian_store << endl;
	  cout << "coarse_rtol \t\t " << coarse_rtol << endl;
	  cout << "final_rtol \t\t " << final_rtol << endl;
	  cout << "maxsparse \t\t " << maxsparse << endl;
          cout << "nxtval_chunk \t\t " << nxtval_chunk << endl;
          cout << "showevals \t\t "<< showevals << endl;
          exit(0);
          break;
        case 'm':
          moints_filename = optarg;
          break;
        case 's':
          maxsparse = atoi(optarg);
          break;
        case 'S':
          showevals = true;
          break;
        case 'c':
          configs_filename = optarg;
          break;
        case 'v':
          vector_filename = optarg;
          break;
        case 'r':
          vector_restart = true;
          break;
        case 'H':
          hamiltonian_restart = true;
          break;
        case 'i':
          subspace_min_mult = atoi(optarg);
          break;
        case 'x':
          subspace_max_mult = atoi(optarg);
          break;
        case 'I':
          subspace_iterations = atoi(optarg);
          break;
        case 'C':
          coarse_rtol  = atof(optarg);
          break;
        case 'W':
          hamiltonian_store = true;
          break;
        case 'F':
          final_rtol  = atof(optarg);
          break;
        case 'N':
          nxtval_chunk = atoi(optarg); 
          if ( nxtval_chunk > 10 ) {
              //cout <<" Warning NXTVAL_CHUNK is large " << nxtval_chunk << endl;
          }
          break;
        case 'V':
          cout << "option -v: PSOCI version is currently " << PSOCI_VERSION << endl;
          exit(0);
          break;
        case 'n':
          numRoots  = atoi(optarg);
          break;
        case 'h':
          cout << "Usage: driver [options] " << endl;
          cout << "option -d: Current list of argument values " << endl;
          cout << "--moints_filename \t\t <string> filename of input moints file " << endl;
          cout << "--vector_filename \t\t <string> filename of input/output CI coefficient file " << endl;
          cout << "--configs_filename \t\t <string> filename of input configuration file " << endl;
          cout << "--subspace_min_mult \t\t (int) subspace_min_mult*numRoots specifies minimum subspace size " << endl;
          cout << "--subspace_max_mult \t\t (int) subspace_max_mult*numRoots specifies maximum subspace size " << endl;
          cout << "--subspace_iterations \t\t (int) maximum number of subspace iterations per tolerance setting " << endl;
          cout << "--vector_restart \t\t (bool) requests a restart - requires a valid CI coef file " << endl;
          cout << "--hamiltonian_restart \t\t (bool) requests a CI matrix restart - requires valid DRA files " << endl;
          cout << "--hamiltonian_store \t\t (bool) requests a CI dump to disk " << endl;
          cout << "--coarse_rtol \t\t\t (float) rnorm convergence for coarse step " << endl;
          cout << "--final_rtol \t\t\t (float) rnorm convergence for final vectors " << endl;
	  cout << "--maxsparse \t\t\t (int) Maximum anticipated number of CI elements per sef" << endl;
          cout << "--nxtval_chunk \t\t\t (int) Chunk size for nxtval task blocking " << endl;
          cout << "--showevals \t\t\t (bool) Print evals every subspace iteration" << endl;
          cout << "--numRoots \t\t\t (int) Number of eigenvectors to converge " << endl;
          cout << "--defaults \t\t Report current settings of argument variables " << endl;
          cout << "--version \t\t Report version of PSOCI " << endl;
          cout << "--help \t\t This message" << endl;
          exit(0);
          break;

        default:
          abort ();
        }
    }
  return;
}



// New declarations
long driverFetchConfigurations( PsociConfigs & configs );


// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//
//  Start driver program 
//
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#ifdef GATRACE
extern "C" {
    void trace_init_( long * num );
    void trace_end_( long * me );
    void trace_stime_( void );
    void trace_etime_( void );
}
#endif

int main(int argc, char **argv)
{
  parseUserValues(argc, argv);
  
  // Set up global MPI fabric
  
  int g_size;
  int g_rank;
  
  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, maxGAMemPerCore );
  
  g_size = GA::nodes();
  g_rank = GA::nodeid();
  if ( GA::usesMA() ) cout << "GA memory is coming from MA " << endl;
  if ( GA::usesFAPI() ) cout << "GA is using Fortran indexing " << endl;

  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  GCOUT << "RENCI Parallel Spin-Orbit Configuration Interaction program (PSOCI)" << endl;
  GCOUT << "Parallel (MPI) Version " << PSOCI_VERSION << endl;
  GCOUT << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
  GCOUT << endl;
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  
  GCOUT << endl;
  GCOUT << " -------------- setup ---------------------" << endl << endl;
  GCOUT << endl;
  GCOUT << "Current settings for this calculation are " << endl;
  GCOUT << "----------------------------------------- " << endl;


  GCOUT << "moints_filename \t" << moints_filename << endl;
  GCOUT << "vector_filename \t" << vector_filename << endl;
  GCOUT << "configs_filename \t" << configs_filename << endl;
  GCOUT << "subspace_min_mult \t" << subspace_min_mult << endl;
  GCOUT << "subspace_max_mult \t" << subspace_max_mult << endl;
  GCOUT << "subspace_iterations \t" << subspace_iterations << endl;
  GCOUT << "vector_restart \t\t" <<vector_restart << endl;
  GCOUT << "coarse_rtol \t\t" << coarse_rtol << endl;
  GCOUT << "final_rtol \t\t" << final_rtol << endl;
  GCOUT << "maxsparse \t\t" << maxsparse << endl;
  GCOUT << "nxtval_chunk \t\t"<< nxtval_chunk << endl;
  GCOUT << "numRoots \t\t" << numRoots << endl;
  GCOUT << " ------------------------------------------------------" << endl << endl;

  GCOUT << "Parallel Run with " << g_size << " Total processors " << endl;
  GCOUT << "Per core HEAP/STACK set to " << heap << " " << stack << endl;
  GCOUT << "Maximum GA space per core constrained to " << maxGAMemPerCore << endl;
  GCOUT << " ------------------------------------------------------" << endl << endl;

// Start the calculation
  
  GCOUT << "GA usesMA ? " << GA::usesMA() << endl;
  GCOUT << "GA memory limited ? " << GA::memoryLimited() << endl;
  GCOUT << "GA inquireMemory is " << GA::inquireMemory() << endl;
  GCOUT << "May not be GA array limit: GA memory avail " << GA::memoryAvailable() << endl;
  
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  GCOUT << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  GCOUT << "Parallel GA Run with " << g_size << " Total processors " << endl;
  
  GCOUT << "residual norm coarse convergence is  set to 0.001 " << endl;
  GCOUT << "residual norm final convergence is  set to " << MIN_PRECISION  << endl;
     
// Perform a restart like configs read
  PsociConfigs configs( configs_filename ); //File is opened/closed at the read
  long configBytes = driverFetchConfigurations( configs );
  GCOUT << "Current per-core Configs memory (MB) " << configBytes/1000000 << endl;
  GCOUT << endl;

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

// Construct determinants
// Need to fetch a subset of configs, create determinants and push those dets to GA space

  GCOUT << "Build distributed GA_DET space " << endl;
  
  //Set this up so we can choose to NOT store constructed data to disk
  
  //HSOURCE construct_hamiltonian = CONSTRUCT;
  HSOURCE construct_hamiltonian = NONE;

  if (hamiltonian_store) {
    GCOUT << "PSOCI will attempt a store of the DRA files " << endl;
    construct_hamiltonian = CONSTRUCT;
  }

  if (hamiltonian_restart) {
    GCOUT << "PSOCI will attempt a restart from the DRA files " << endl;
    construct_hamiltonian = RESTART;
  }

/////////////////////////////////////////////////
/* process configurations and subsequent generation of the determinants 
*/
  GA::GlobalArray * g_det; // declare global det data 
  GA::GlobalArray * g_det_num;
  GA::GlobalArray * g_nsef;
  GA::GlobalArray * g_orbmap;
  
/* NON-direct version 
  PsociGADeterminants deters_det( g_det, g_det_num, g_nsef );
  deters_det.generateDistributedDeterminants( configs );
*/
// DIRECT code follows
  PsociGADeterminants deters_det( configs, g_det, g_det_num, g_nsef, g_orbmap );
  deters_det.generateDistributedDeterminants();
  
  vector<int> l_nsef;
  deters_det.assembleAndBrdcstNsef( l_nsef ); //don't need l_nsef anymore
/* The next only for direct*/
//  deters_det.destroyJoutfgGA(); // Need these in analysis
  
/* At this point all determinants should have been generated AND pushed into GAdeterminant space.
*/
  
   //   string filername="moints";
   //   string filername="/home/jtilson/RuO+2-NewBasis-SDCI-PSOCIQA/moints";
   //   string filername="./INTEGRALS/RuOSCI-B1-2.75au/moints";
   
///////////////////////////////////////////////
/* Open up moints file, read it and brdcst all ints 
   The file open/close must be handled by the Fortran layer because
   in the Columbus sifrd2 layer (2e ints) it atttempts to open the 
   2e data file 
*/ 
   const int ROOT_READER=0;
   Fint l_unit = 2;
   PsociIntegrals mos( ROOT_READER, l_unit, moints_filename );
   mos.fetchAndProcessIntegrals();
   
   double fcore =  mos.fetchCoreEnergy(); // core energy
   int maxsef= deters_det.fetchGlobalMaxSefs(); // current problem size
   
   GCOUT << "Current global maxsef is " << maxsef << endl;

   if ( GA::nodeid() == 0 ) GA::printStats();
 
/////////////////////////////////////////////////
/* Begin the generation of the hamiltonian matrix
*/
 
   /* Allocate DRAs: This is always needed: However we want to eliminate the chance that existing H data
      might get destroyed by accident.
   */  
   // First simply start up the service: needed by all
   
   PsociDRAservices DRAservice;

   int multIOservers = 1;
   int numIOservers = multIOservers * GA::clusterNnodes();

   //InitializeDRAservice( g_rank, DRAservice, g_size );
   
   GCOUT << "Number of total IO servers is " << numIOservers << endl;

   InitializeDRAservice( g_rank, DRAservice, numIOservers );
   
   //Now allocate the DRA arrays: currently as RW: BUt we may want to be more picky.
   //Declare the basis objects used by both restart modes
   
   string CIarrayname = DRA_CIfilename;
   string CIfilename = DRA_CIfilename;
   PsociDRArestart d_cimat( g_rank, CIarrayname, CIfilename );
 
   string ICIarrayname = DRA_ICIfilename;
   string ICIfilename = DRA_ICIfilename;
   PsociDRArestart d_icicol( g_rank, ICIarrayname, ICIfilename );
   
   string NUMarrayname = DRA_NUMfilename;
   string NUMfilename = DRA_NUMfilename;
   PsociDRArestart d_number( g_rank, NUMarrayname, NUMfilename );
   
   string DIAGarrayname = DRA_DIAGfilename;
   string DIAGfilename = DRA_DIAGfilename;
   PsociDRArestart d_diag_sef( g_rank, DIAGarrayname, DIAGfilename );
 
   // Now choose mode CONSTRUCT or RESTART ( or NONE )
   
   int readSparse = 0;
   int retvalue=0;

   if ( construct_hamiltonian == CONSTRUCT ) {
     int ndims = 2;
     GCOUT << "Preparing to Construct Hamiltonian " << endl;
     AllocateDRA( construct_hamiltonian, ndims, maxsef, maxsparse, d_cimat, d_icicol, d_number, d_diag_sef);

   } else if( construct_hamiltonian == RESTART ) { // Plan on fetching hamiltonian   

     GCOUT << "Preparing to READ hamiltonian data from disk " << endl;
     retvalue = SpecifyDRA(construct_hamiltonian, d_cimat,d_icicol,d_number,d_diag_sef, readSparse);  
     GCOUT << "READ SPARSE was " << readSparse << endl;
     
   } // if NONE then neither read nor dump - NOT REALLY an option to the user yet
   
   /* Still need this because the GA space needs to be initialized and the constructor actually
      creates the array space
   */

   if( construct_hamiltonian == RESTART ) {
     GCOUT << "On Hamiltonian READ a values for sparsity were identified: Changing to match " << readSparse << endl;
     maxsparse = readSparse;
   }

#ifdef SUPPLEMENTAL
 cout << "NO SUPPLEMENTAL Current implementation of driver.PSOCI.x Use driver.PSOCI.splittest.x instead " << endl;  
 GA::Terminate();
#endif

 PsociGAhamiltonian hamilton(maxsparse, maxsef, &deters_det, &mos );
 GCOUT << "Resetting NXTVAL CHUNK to " << nxtval_chunk << endl;
 hamilton.setNxtvalChunk( nxtval_chunk );
 
 // Extra bit neede by all but data comes from Det processing
 
 int nsize = hamilton.computeNsefEntryPoints( l_nsef );
 //hamilton.printNsefEntryPoints();
 
 hamilton.printKsym();
 hamilton.printPerSpatialValues();
 hamilton.hamiltonianDistribution();

// Compute H

/* These now not used anymore  applies only to the fetch2D code
   hamilton.SetGranularity( 1 );
   hamilton.PrintGranularity();
*/

#ifdef GATRACE
   long maxnum = 10000;
   long * numevents = &maxnum;
   trace_init_( numevents );
#endif

 pair<int,double> Hinfo;
 if ( construct_hamiltonian == CONSTRUCT || construct_hamiltonian == NONE ) {

   
#ifdef DIRECT
   GCOUT << "Construct DIRECT distributed hamiltonian " << endl;
   long numH = hamilton.constructGlobalHamiltonianDirect(Hinfo);
#else
   GCOUT << "Construct distributed hamiltonian " << endl;
   long numH = hamilton.constructGlobalHamiltonian(Hinfo);
#endif

   GCOUT  << "Total number of elements above threshold  " << numH << endl;
   GCOUT  << "Timing: Me is " << Hinfo.first << " time is " << Hinfo.second << endl;
   
   if ( g_rank == 0 ){
      GA::printStats();
   }

   if ( construct_hamiltonian != NONE ) {
     pair<int,double> dumpTime;

     SendHamiltonianToDRA(construct_hamiltonian,d_cimat,d_icicol,d_number,d_diag_sef,hamilton,dumpTime);
     GCOUT << dumpTime.first << "reports a dump to DRA time of " << dumpTime.second << endl;
   }
   
 } else if( construct_hamiltonian == RESTART ) { // Attempt to read hamiltonian from disk
   
   double timeread = psociTime(); 
   pair<int,double> readTime;
   FetchHamiltonianFromDRA(construct_hamiltonian,d_cimat,d_icicol,d_number,d_diag_sef,hamilton,readTime);
   GCOUT  << readTime.first << "reports a read DRA time of " << readTime.second << endl;
 }
 
 // At this point we have the GA data in place We actually do not need DRAS anymore either.
 
 ExitDRASERVICE( DRAservice );

  GCOUT << "Finished with DRA " << endl;
#ifdef GATRACE
   long * gme = (long *)&g_rank;
   trace_end_( gme );
#endif

 
 // Move onto the iterative diagonalization
 // Remember hamiltonian computes the matrix-vector products
 
 // GA::printStats();
 
 if ( g_rank == 0 ) hamilton.printDiags();
 
 // Move onto the iterative diagonalization step.
 // Open up some 1D pretend arrays
 
 //////////////////////////////////////////////////////////////
 /*
   PsociGArestart() -> create GA space and read vectors - or - guestimate if desired.
   PsociGAbasis ( handle from above )
   PsociGArestart() -> dump current data back to disk ( could have changed size )
   
   1) PsociGArestart()
 */
 
 /* Choose the value for maxsubspace based on preliminary testing with RuO2+, VFCI wavefunctions
    at 40 roots. It was set to be 2* greater than the minimum post-compressed subspace
 */
 
 int numSoughtRoots = numRoots; // more like 50 to 100 int maxsubspace = 4 * numSoughtRoots; // 4 expansions per eigenvector (
 int maxsubspace =  subspace_max_mult * numSoughtRoots;
 int minvecs = subspace_min_mult; // in compression a multipland 

/* Keeping 2 or 3 expansion vectors per sought vector was based on some initial testing
   using the VFCI RuO2+ configuration but is generally consistent with

   M. Crouzeix, B. Philippe, M. Sadkane, "The Davidson Method" 
   SIAM J. Sci. Comput. Vol 15, No. 1 pp 62-76, January 1994.

and also
 
   J.H. van Lenthe, P. Pulay, "A Space Saving Modification of Davidson's Eigenvector ALgorithm",
   J. Comp. Chem. Vol 11, No. 10, pp 1164-1168, 1990.
*/

// Maximum number of iterations in the solution 
 int max_iterations_per_rtol = subspace_iterations;
 
  

 int vdim=2;
 int dims[2],chunk[2];
 dims[0] = maxsef ;
 dims[1] = numSoughtRoots+maxsubspace;
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
 g_v->zero();
 g_Hv = GA::createGA( C_DBL, vdim, dims, (char *) Hvector_title, chunk);
 g_Hv->zero();
 
 GCOUT << "NUMBER of converged sought vectors is " << numSoughtRoots << endl;

// Only need enough rnorms for the sought vectors - of course one cannot reallocate after this

 dims[0] = maxsef ;// maxsef upwards of 5,000,000;
 dims[1] = numSoughtRoots;

 // Want distribution to be close to that of cimat and friends
 g_rnorms = GA::createGA( C_DBL, vdim, dims, (char *) rnorm_title, chunk);
 g_rnorms->zero();
 
/* Get initial guess vectors and process initial H*v products
   Attempt to read vectors on restart. We could defer to diag in the event of a falure
   but let's trap it and exit
*/

 PsociGAvectors gavectors(numSoughtRoots, g_v );

 string resname=vector_filename;
 if (vector_restart ) {
       GCOUT << "Attempt a restart on preexisting CI coefs in the file " << resname << endl;
       gavectors.openVectorRead( GA::nodeid(), resname );
       gavectors.readVector( GA::nodeid(), g_v );
       int restartNumber = gavectors.getNumberVectors();
       GCOUT << "Vector Restart - current number of READ vectors is " << restartNumber << endl;
} else {
      GCOUT << "Build NEW CI coeficients from a diagonal guess " << resname << endl;
}
  gavectors.specifyDiagonalGuessVectors( g_v ); // handles both cases and appends/contracts as necessary
 //gavectors.printCurrentSoughtRoots( g_v );
  
 
 /* Note in the following one could also use PsociGAbasis to augment the system
    followed by manually resetting number of roots requested using the resizeVectors method
 */
 // Note I do not need gavectors anymore until the end when I want to dumpo solution vectors
 // and perform some analysis
 
 /*
   2) PsociGAbasis  
 */
 
  PsociGAhamiltonian * CIhandle = &hamilton;
  PsociGAbasis vector_set( numSoughtRoots, CIhandle , g_v, g_Hv );
  vector<double> vecnorm;
  
  vector_set.normAllVectors(); 
  int or_status = vector_set.orthonormalizeAllVectors();
  if ( or_status != 0 ) {
    cout << "A linear dependancy was found at i " << or_status << endl;
    GA::error(" linear dependancy: This currently is a hard failure: Aborting at ", or_status );
  }
  
  GCOUT << GA::nodeid() << " Sought roots " << vector_set.CurrentNumSoughtRoots() << endl;
  GCOUT << GA::nodeid() << " Existing roots " << vector_set.CurrentExistingRoots() << endl;
  GCOUT << GA::nodeid() << " Max roots " << vector_set.CurrentMaxRoots() << endl;
  GCOUT << GA::nodeid() << " Req products " << vector_set.CurrentReqHVProducts() << endl;
  GCOUT << GA::nodeid() << " Existing products " << vector_set.CurrentExistingHVProducts();
  
  //Need to compute some H*v products  for initial set of vectors
  
  int num = vector_set.computeHvProducts( 0, numSoughtRoots );
  vector_set.setpresent_matrixProducts( numSoughtRoots ); // tell the code some were already computed
// IMPORTANT: This "setpresent" call is here because we PRECOMPUTED some terms
// NOT DOING THIS ADDS A FW H*v TO THE diagonalization startup
  GCOUT << "Computed " << num << " H*v products " << endl;
  
  //////////////////////////////////////////////////////////////////////////////
  /* 
     Take this current set of vectors and apply them to PsociGAsubspace
     It is there where the real subspace work occures.
  */
  
  //PsociGAbasis does the preconditioning and H*v work
  
  GCOUT << "INITIALIZE the SUBSPACE " << endl;

  vector<double> evals( maxsubspace ); // Need these for preconditioning - level shifting


  PsociGAbasis * sub_vector_set = &vector_set;
  
  PsociGAsubspace subspace( sub_vector_set, hamilton.fetchDIAG_SEF(), g_rnorms  );
  
  GCOUT << "maxsubspace is " << subspace.fetchMaxSubspace() << endl;
  GCOUT << "existing subspace " << subspace.fetchExistingSubspace() << endl;
  GCOUT << "subspace CI rows is " << subspace.numCIBasis() << endl;
  GCOUT << "subspace Vec rows is " << subspace.numVectorBasis() << endl;
  GCOUT << "subspace Vec cols is " << subspace.numRoots() << endl;
  
  bool converged=false;
  vector<double> rnorms( numSoughtRoots, -1.0 ); // Replicated list of total rnorms
  int count=0;
  int total_count=0;
  int numcont=0;
  
  int curroot=1; //Start optimize from the first root
   
  // What is the right value for here?
  double rtol = final_rtol; 
  
  /* The current algorithm consists of ensuring ALL sought vectors are converged
     before declaring success.
  */
   
  vector<bool> vecs_converged( vector_set.CurrentNumSoughtRoots(), false ); 
  
  pair<int,double> diagTime;
  double diag = psociTime();
  
  double coarse_rtol = 0.001; //Initial coarse first-pass
  
  /* For this current approach we first attempt to converge to a coarse levle of tolerance
     then repeat at the higher level of tolerance. Clearly much better ways are possible
     but that will need to come later
     
     At the start of the tighter tolerance cycle we use as starting vectors the current set
     of converged vectors by calling the resize method. This reorthogonalizes the desired set
     and recomputes all H*v terms
  */
  
  /* For now we do not have a decent way to reduce the space so only do one iteration for now
   */

#ifdef RCRMONITOR
  if (GA::nodeid() == 0 ) cout << "Enabled RCR logging at the Subspace solution" << endl;
  char * logfilename = "loggerSubspace.txt";
  int stats =  RCRstartLogging( logfilename );

  char * location = "CnsrtSbspce";
  stats = RCRpushLocation(location);
#endif

  double current_rtol = coarse_rtol; //start low then increase later
  
  count=0;
  converged=false;


  int compressSubspace=0;
  int matrixVectorProducts=0;
  
  matrixVectorProducts += vector_set.CurrentExistingRoots(); // They were prebuilt in this example

  double tempTime2, tempTime, genTimes =0.0, rnormTimes=0.0, diagTimes=0.0, precondTimes=0.0, expandTimes=0.0, appendTimes=0.0; // debug performance issues in the solution steps
  double rnormvalue;

#ifdef NEWORTHO
  GCOUT << "Using NEW Orthonormalization procedure for appended vectors " << endl;
#endif
#ifdef NEWRNOMRS
  GCOUT << "USING NEW Residual norms CODE " << endl;
#endif
#ifdef DIAGGOP
  GCOUT << "USING diagGOP version " << endl;
#endif
#ifdef INCREMENTAL
  GCOUT << "USING INCREMENTAL GENERATE " << endl;
#endif
#ifdef NEWORTHOALL
  GCOUT << "USING alternative data-parallel orthoALL for MGS " << endl;
#endif
#ifdef NEWGETAPPEND
  GCOUT << "USING New getappend methods " << endl;
#endif
#ifdef  MATRIXPRODUCTGOP
 GCOUT << "USING alternative matrixVector product methods" << endl;
#endif

  int numUpdates = vector_set.CurrentExistingRoots();

  while( count < max_iterations_per_rtol && converged != true) {
    ++total_count;
    ++count;
    
    tempTime = psociTime();
    tempTime2 = tempTime;


#ifdef INCREMENTAL
    subspace.generateIncrementalSubspace( numUpdates );
#else
    subspace.generateSubspace();
#endif

    genTimes += psociTime() - tempTime;
    tempTime = psociTime();

#ifdef DIAGGOP
    subspace.diagSubspaceGOP();
#else
    subspace.diagSubspace(); 
#endif
    subspace.brdcstEvals( evals );

    diagTimes += psociTime() - tempTime;
    tempTime = psociTime();

/* This approach is a little sloppy. It doesn't go back and check if the 'Old' vectors 
   are getting poor again. However, for low rtol this doesn;t matter and for final_rtol 
   it doesn't seem to occur. Doing it this way pays big performance dividends
*/
#ifdef NEWRNOMRS
    converged = true;
    int itest = curroot; 
    while ( itest <= numSoughtRoots ) {
    int status = subspace.computeRnormSingle(itest, rnormvalue );
    rnorms[itest-1] = rnormvalue; // sometimes a root may be converged by chance ( and skipped )
    if ( rnormvalue > current_rtol ) {
    converged = false;
    curroot = itest; //So we could make this an array of rnorms to keep
    break;
    }
    ++itest;
    }
#else
    subspace.computeRnorms( rnorms );
    converged = true; //hope for the best
    for ( int i=0; i < numSoughtRoots; ++i ) { //Check all rnorms
      if ( rnorms[i] > current_rtol ) {
	curroot = i+1; // Not an increment simply choose current vector 
	converged = false;
	break;
      }
    }
#endif

    rnormTimes += psociTime() - tempTime;

    if ( showevals && GA::nodeid() == 0 ) {
     cout << "Current vector under optimization " << curroot<<" "<<rnorms[curroot-1]<<endl;
      cout << "Current roots (au) : rnorm" << endl;
      for(int i=0; i< numSoughtRoots; ++i ) {
        cout << evals[i] << " : " << rnorms[i] << endl;
      }
      subspace.printVectorContents( rnorms );
      cout << "Current solution timers report " << diagTimes<<" "<<" "<<genTimes<<" " <<rnormTimes<<" "<<expandTimes<<" "<<appendTimes << " " << precondTimes << endl;
    }
    
    if( converged == false ) {
      tempTime = psociTime();
      double shift = evals[ curroot-1 ]; // Try to update vector == curroot 
      int nump = vector_set.preconditionAndNewVectorIncore(curroot, shift, g_rnorms, hamilton.fetchDIAG_SEF() );
      precondTimes += psociTime() - tempTime;
      tempTime = psociTime();
      vector_set.appendScratchVector(); // at this point g_scratch is populated and ready for appending.
      matrixVectorProducts += 1;
      appendTimes += psociTime() - tempTime;
      numUpdates = 1;
    }
    
/* Note on vector_restart it is likely all vectors will be converged to the lower rtol 
   But the subspace has not been built up yet so we do not want to compress the space yet
*/
    if ( converged == true && current_rtol == coarse_rtol ) {
      tempTime = psociTime();
      converged = false;
      count=0; // it's maxcounts PER rtol
      current_rtol = rtol;
      curroot = 1; // Start over
      if ( vector_set.CurrentExistingRoots() >= minvecs * numSoughtRoots ) {
        vector_set.expandAndReplaceVectors( subspace.fetchEvecsHandle(), minvecs * numSoughtRoots );
        ++compressSubspace;
        matrixVectorProducts += vector_set.CurrentExistingRoots(); // All were recalculated 
      }
      expandTimes += psociTime() - tempTime;
      numUpdates = vector_set.CurrentExistingRoots();
    } 

// Expand current subspace vectors into g_vector ( keep subspace at least minvecs * vectors in size )
    if ( vector_set.CurrentExistingRoots() > maxsubspace ) {
      tempTime = psociTime();
      vector_set.expandAndReplaceVectors( subspace.fetchEvecsHandle(), minvecs * numSoughtRoots ); 
      ++compressSubspace;
      matrixVectorProducts += vector_set.CurrentExistingRoots(); // All are recalculated 
      expandTimes += psociTime() - tempTime;
      numUpdates = vector_set.CurrentExistingRoots();
    }

  diagTime.second +=  psociTime() - tempTime2;

  } // converged and/or out of cycles

/* COMPUTE the actual all-roots set of rnorms for printing out since the new
   high performance approach doesn't update old vectors all the time
*/
#ifdef NEWRNOMRS
// Can fail if not all 40 roots have been processed

  GCOUT << "Compute a last all-vector residual matrix " << endl;
  subspace.computeRnorms( rnorms );
#endif

  diagTime.first = GA::nodeid();

  for(int i=0; i< numSoughtRoots; ++i ) {
     if ( rnorms[i] <= 0 ) {
     GA::error("An rnorm was missed at i = ", i+1 );
     }
  }
// The final set of solutions must be re-expanded into the full space and compressed to
// numSoughtRoots.

  vector_set.expandAndReplaceVectors( subspace.fetchEvecsHandle(), vector_set.CurrentExistingRoots() );
  vector_set.resizeVectors( numSoughtRoots ); // This zeros out all higher order roots

#ifdef RCRMONITOR
  stats =  RCRpopLocation();
  stats =  RCRstopLogging();
#endif

// Dump current full set of vectors to disk for possible restart or post processing
  string vecname=vector_filename;
  GCOUT << "Writing current set of vector coefs to disk to the file " << vecname << endl;
  gavectors.openVectorWrite( GA::nodeid(), vecname );
  gavectors.dumpVector(GA::nodeid(), numSoughtRoots, g_v );
  gavectors.closeVectorFile();
  GCOUT << " ----------------------------------------------------" << endl << endl;
  
  GCOUT << " DiagTime is " << diagTime.second  << endl;
  if ( GA::nodeid()== 0 ) {
    GCOUT << "Hamiltonian construction time " <<  Hinfo.second << endl;
    GCOUT << "Number of SUBSPACE contractions is " << compressSubspace << endl;
    GCOUT << "Number of total iterations is " << total_count << endl;
    GCOUT << "TOTAL MATRIX VECTOR PRODUCTS is " << matrixVectorProducts << endl;
    GCOUT << "Current results: fcore is " << fcore << " " << endl;
    for(int i=0; i< vector_set.CurrentNumSoughtRoots(); ++i ) {
      GCOUT << "Eigenvalue " << i << " Total energy (au) " << setprecision(6) << evals[i]+fcore <<  " Electronic (au) " << evals[i] << " rnorm " << setprecision(10) << rnorms[i] << endl;
    }
  }
  
  GCOUT << "GA rank " << diagTime.first << " Diagonalization time (sec) is " << diagTime.second << endl;
  GCOUT << "At end of iteration " << total_count << " with all converged = " << converged << endl;
  GCOUT << "Total cycles is " << total_count << endl;
  GCOUT << "Number of subspace contractions is " << numcont << endl;
  
  subspace.destroySubspace();

// Need diags for the analysis part
  vector<double> diags;
  hamilton.fetchDiags( diags );

  hamilton.freeHamiltonian(); // Free up space 
  
  GCOUT << "Begin CI analysis for numRootsSought = " << vector_set.CurrentNumSoughtRoots() << endl;

  pair <double,double>  thresh;
  thresh = 0.0001;
  int num3 = analyzeEnergyContributions(thresh, fcore,l_nsef, evals, diags, g_v, vector_set, hamilton, deters_det);
  GA::sync();

  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX finished XXXXXXXXXXXXXXXXXXX " << endl;
  GA::Terminate();
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
      GA::error("Failed to read configs at numfgs =",-1);
    }
    configs.printParams();
    configs.printConfigurations(); 
  }
  configs.brdcstConfigs( 0 ); //now replicated to all nodes
  
  long maxspatials;
  long nbasis;
  vector<pair<int,string> > word;
  maxspatials = configs.numTotalSpatials();
  nbasis = configs.numTotalSpatials();
  GCOUT << "Num total spatials is " << maxspatials <<" Num total basis is " << nbasis<< endl;
  long memory = maxspatials*nbasis*sizeof(int);
  return( memory );
}



