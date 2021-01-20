/**************************************************************************************
 * Copyright (c) 2010,2011 RENCI.
 * All rights reserved. This program and the accompanying materials
 * MAY BE available under the terms of the RENCI Open Source License
 * UNC at Chapel Hill which accompanies this distribution, and is available at
 * http://www.renci.org/resources/open-source-software-license

 * New implementation of PSOCI:

 * Preliminary out-of-core RCM implementation
 
 Classes: 

 Description: 
 
 History:

**************************************************************************************/
/**
 *   @file driver.PSOCI.splittest.RCM.hpp
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

extern "C"
{
//#include <mpi.h>
#include <ga.h>
}

#include <dra.h>
#define GA_DATA_TYPE C_DBL

// We want to create subdirectories for storing H data
#include <sys/types.h>
#include <sys/stat.h>

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
//#include "PsociProcessMOcoefs.hpp"
#include "PsociNOpopulations.hpp"
#include "PerformRCM.hpp"


#include "driver.PSOCI.hpp"

void parseUserValues(int argc, char **argv)
{
  int c;
  while (1)
    {
      static struct option long_options[] =
        {
          {"defaults (d)",            no_argument,       0, 'd'},
          {"numRoots (n)",            required_argument, 0, 'n'},
          {"numNOroots (R)",          required_argument, 0, 'R'},
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
          {"vectorchunk (P)",         required_argument, 0, 'P'},
          {"supp_maxsparse (u)",      required_argument, 0, 'u'},
          {"percent_maxwidth (p)",    required_argument, 0, 'p'},
          {"nxtval_chunk (N)",        required_argument, 0, 'N'},
          {"granularity (g)",         required_argument, 0, 'g'},
          {"numJchunks (j)",          required_argument, 0, 'j'},
          {"scatter_products (T)",    no_argument,       0, 'T'},
          {"showevals (S)",           no_argument,       0, 'S'},
          {"version (V)",             no_argument,       0, 'V'},
          {"help (h)",                no_argument,       0, 'h'},
          {0, 0, 0, 0}
        };
      int option_index = 0;
      c = getopt_long (argc, argv, "TP:R:Wn:s:dm:c:g:v:ri:a:I:C:F:N:hHVSu:p:",
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
          cout << "vectorchunk \t\t " << l_vectorchunk << endl;
          cout << "supp_maxsparse \t\t " << supp_maxsparse << endl;
          cout << "percent_maxwidth \t " << percent_maxwidth << endl;
          cout << "granularity \t\t " << l_granularity << endl;
          cout << "numJchunks \t\t " << l_numJchunks << endl;
          cout << "nxtval_chunk \t\t " << nxtval_chunk << endl;
          cout << "showevals \t\t "<< showevals << endl;
          cout << "Number of CI roots \t " << numRoots << endl;
          cout << "Number of NO sets \t " << numNOroots << endl;
          cout << "scatter_products \t " << scatter_products << endl;
          exit(0);
          break;
        case 'm':
          moints_filename = optarg;
          break;
        case 'R':
          numNOroots = atoi(optarg);
          break;
        case 's':
          maxsparse = atoi(optarg);
          break;
        case 'S':
          showevals = true;
          break;
        case 'j':
          l_numJchunks = atoi(optarg);
          break;
        case 'c':
          configs_filename = optarg;
          break;
        case 'g':
          l_granularity = atoi(optarg);
          break;
        case 'P':
          l_vectorchunk = atoi(optarg);
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
        case 'W':
          hamiltonian_store = true;
          break;
        case 'T':
          scatter_products = true;
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
        case 'F':
          final_rtol  = atof(optarg);
          break;
        case 'u':
          supp_maxsparse  = atoi(optarg);
          break;
        case 'p':
          percent_maxwidth  = atoi(optarg);
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
          cout << "--vectorchunk \t\t\t (int) Number of chunks to break up vector in H*v (1-6 is opt) " << endl;
          cout << "--supp_maxsparse \t\t (int) Maximum anticipated number of supplemental CI elements per sef" << endl;
          cout << "--percent_maxwidth \t\t (int) Maximum (%) number of CI sefs in supplemental CI" << endl;
          cout << "--nxtval_chunk \t\t\t (int) Chunk size for nxtval task blocking " << endl;
          cout << "--granularity \t\t\t (int) Max width for blocked configuration fetching " << endl;
          cout << "--numJchunks \t\t\t (int) Number of J spatial chunks for H/CI build" << endl;
          cout << "--showevals \t\t\t (bool) Print evals every subspace iteration" << endl;
          cout << "--numRoots \t\t\t (int) Number of eigenvectors to converge " << endl;
          cout << "--numNOroots \t\t\t (int) Number of eigenvectors for populations" << endl;
          cout << "--scatter_products \t\t (bool) Choose ga_scatter method for H*v products " << endl;
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
  
/* Newer version of ARMCI prefer you to Init MPI yourself*/


  int g_size;
  int g_rank;

#ifndef ANLARMCI

//cout << "inside get configs" << endl;
#ifdef MPITHREADMULTIPLE
  int  mpiprovided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpiprovided );
//  cout << "provided is " << mpiprovided << endl;
#else
  MPI_Init(&argc, &argv);   
#endif

//  GA_Initialize_ltd(maxGAMemPerCore);
  GA_Initialize();

#ifdef MPITHREADMULTIPLE
  if ( mpiprovided != MPI_THREAD_MULTIPLE ) GA::error("MPI_THREAD_MULTIPLE not provided",mpiprovided);
#endif

  int memoryStackWords = stack/8;
  int memoryHeapWords  = heap/8;
  int maxGAMemPerCoreWords = maxGAMemPerCore/8;

  if ( GA_Uses_ma() ) {
    if ( GA::nodeid() == 0 ) cout << "INIT from within " << endl;
    MA_init(C_DBL, memoryStackWords, memoryHeapWords+maxGAMemPerCoreWords );
  } else {
    if ( GA::nodeid() == 0 ) cout << "INIT without GA space" << endl;
    MA_init(C_DBL, memoryStackWords, memoryHeapWords);
  }

#ifdef MPITHREADMULTIPLE
   if ( GA::nodeid() == 0 ) cout << "MPITHREADMULTIPLE Enabled" << endl;
#endif

#endif

#ifdef ANLARMCI

#ifdef GAEXP 
// Needed because the full init call causes assertion problems
  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, 0 );
#else
  if ( GA::nodeid() == 0 ) cout << "maxGAMemPerCore is " << maxGAMemPerCore << endl;
  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, maxGAMemPerCore );
#endif

#endif
  
  g_size = GA::nodes();
  g_rank = GA::nodeid();

// if(! MA_init(C_DBL, stack, heap) ) GA_Error("MA_init failed",stack+heap);

  if ( GA::usesMA() ) cout << "GA memory is coming from MA " << endl;
  if ( GA::usesFAPI() ) cout << "GA is using Fortran indexing " << endl;

  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  GCOUT << "RENCI Parallel Spin-Orbit Configuration Interaction program (PSOCI)" << endl;
  GCOUT << "Parallel (MPI) Version " << PSOCI_VERSION << endl;
  GCOUT << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
  GCOUT << endl;
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  GCOUT << "Enforces DRA write/read to perform RCM XXXXXXXXXXXXXXXXX" << endl << endl;
  
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
  GCOUT << "Scatter H*v products is " << scatter_products << endl;
  GCOUT << "maxsparse \t\t" << maxsparse << endl;
  if ( l_vectorchunk <= 0 ) {
     l_vectorchunk = GA::nodes();
     GCOUT << "resetting vectorchunk to g_size" << endl;
  }
  GCOUT << "vectorchunk \t\t" << l_vectorchunk << endl;
#ifdef SUPPLEMENTAL
  GCOUT << "supp_maxsparse \t\t" << supp_maxsparse << endl;
  GCOUT << "percent_maxwidth \t" << percent_maxwidth << endl;
#else
  GCOUT << "Supplemental not being used "<< endl;
  GCOUT << "supp_maxsparse and percent_maxwidth ignored " << endl;
#endif
  if ( nxtval_chunk > 10 ) {
    GCOUT << "The new ORBMAP Hamiltonian doesn't need large values of nxtval_chunk anymore. In fact it may impact performance" << endl;
    GCOUT << "We recommend no more than 10: For now we wil reset it for you  " << endl;
    nxtval_chunk = 10;
  }
  GCOUT << "nxtval_chunk \t\t"<< nxtval_chunk << endl;
  GCOUT << "granularity \t\t"<< l_granularity << endl;

  if ( l_numJchunks > 100 ) {
    GCOUT << "We recommend numJchunks no larger than 100  We will lreset it for you " << endl;
    l_numJchunks = 100;
  }
  GCOUT << "numJchunks \t\t" << l_numJchunks << endl;
  GCOUT << "numRoots \t\t" << numRoots << endl;
  if ( numNOroots < 0 ) numNOroots = numRoots;
  if ( numNOroots > numRoots ) numNOroots = numRoots;
  GCOUT << "numNOroots \t\t" << numNOroots << endl;
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
  
  GCOUT << "residual norm coarse convergence is  set to " << coarse_rtol << endl;
  GCOUT << "residual norm final convergence is  set to " << final_rtol   << endl;

#ifdef USEMPI
  GCOUT << "WIll use MPI collective in the matrix vector products" << endl;
#endif
     
#ifdef GETANDBRDCST
  GCOUT << "Using rank=0 + brdcst methods in matrixVector prods " << endl;
#else
  GCOUT << "Using parallel get+dgop methods in matrixVector prods " << endl;
#endif

// Perform a restart like configs read
  PsociConfigs configs( configs_filename ); //File is opened/closed at the read
  long configBytes = driverFetchConfigurations( configs );

  GCOUT << "Current per-core Configs memory (MB) " << configBytes/(1024*1024) << endl;
  GCOUT << endl;

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

// Construct determinants
// Need to fetch a subset of configs, create determinants and push those dets to GA space

  GCOUT << "Build distributed GA_DET space " << endl;

#ifdef CONFIGSCRAMBLE
  GCOUT << "USING SCRAMBLE methods to build GA_DET space " << endl;
#endif
  

#ifdef NEWACCUMULATEDSCATTER

#ifdef TACC
  GCOUT << " USING TACC LowLevel method " << endl;
#endif
#ifdef LOOPNPUT 
  GCOUT << " USING LOOPNPUT LowLevel method " << endl;
#endif
#ifdef SCATFLAT 
  GCOUT << " USING SCATFLAT LowLevel method " << endl;
#endif
#ifdef SCATNOFLAT 
  GCOUT << " USING SCATNOFLAT LowLevel method " << endl;
#endif

#endif

// This works fine for Hopper

GCOUT << "Number of cluster nodes is " <<  GA::clusterNnodes() << endl;;
GCOUT << "I am cluster node " << GA::clusterNodeid() << endl;

  // Pretend the need to do a restart set up as an enumeration
  //Set this up so we can choose to NOT store constructed data to disk
  
  //HSOURCE construct_hamiltonian = CONSTRUCT;

// OVERRIDE for RCM testing

   hamiltonian_store=true;
   GCOUT<< "OVERRIDE value for hamiltonian_store: set to TRUE"<<endl;
  

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

  pair<int, double> gendettime;

#ifdef DIRECT
  GCOUT << " DIRECT METHOD " << endl;
#else
  GCOUT << "NONDIRECT WAVEFUNCTION METHOD" << endl;
#endif

// For now must use SINGLENODECONFIGS
#ifdef DIRECT
#ifndef SINGLENODECONFIGS
  GA::error("DIRECT compiled code must also specify SINGLENODECONFIGS ",1);
#endif
   PsociGADeterminants deters_det( &configs,  g_det, g_det_num, g_nsef, g_orbmap);
//  deters_det.generateDistributedDeterminants(); // One needs to delete configs 
   deters_det.generateDistributedDeterminantsSmallMemoryGA(gendettime );
#else

   PsociGADeterminants deters_det( &configs, g_det, g_det_num, g_nsef, g_orbmap );
#ifdef SINGLENODECONFIGS
  GCOUT << "Build determinates with single core configs method: small memory footprint" << endl;
  deters_det.generateDistributedDeterminantsSmallMemoryGA(gendettime );
#else
  deters_det.generateDistributedDeterminants();
#endif


#endif

#ifdef SINGLENODECONFIGS
  GCOUT << "Time for distributed determinant computations: me is " << gendettime.first<<" time " << gendettime.second << endl;
#endif

  vector<int> l_nsef;
  deters_det.assembleAndBrdcstNsef( l_nsef ); //don't need l_nsef anymore
  
  int numBasisFtns = deters_det.fetchNumBasisFunctions(); // Need this for CIDEN work

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

#ifdef SUPPLEMENTAL
   const int numfiles = 7;
#else
   const int numfiles = 4;
#endif

 
   /* Allocate DRAs: This is always needed: However we want to eliminate the chance that existing H data
      might get destroyed by accident.
   */  
   // First simply start up the service: needed by all
   
   PsociDRAservices DRAservice;

   int multIOservers = 1;
   int numIOservers = max( 1, multIOservers * GA::nodes() ); // two per node

   //numIOservers = 1;
   // I have no idea if GA doesn anyting smart here about distribution
   numIOservers = GA::clusterNnodes();

//Whoops make these logicals
// But not yet.
// for big hopper run use this: numIOservers == GA::nodes() / 12;

   if ( numIOservers == GA::nodes() && GA::nodes() >= 16 ) {
       numIOservers == GA::nodes() / 16; // Im TACC and cluster nodes don;t work here
        if ( GA::nodeid() == 0 ) cout << "cluster nodes don't work: divide by 16 " << endl;
   }

   numIOservers /= numfiles; 
   
   GCOUT << "Number of total IO servers/files per GA is " << numIOservers << endl;

   InitializeDRAservice( g_rank, DRAservice, numIOservers );

GCOUT << "done with INIT DRA" << endl;
   
   //Now allocate the DRA arrays: currently as RW: BUt we may want to be more picky.
   //Declare the basis objects used by both restart modes
   
   // Create subdirectories to mitigate performance problem on largew core count runs
   // Also must account for num IO servers < g_size

#ifdef NOSUBDIR
// all files dumped into working dir
   string newsubdirname="./dra-0/";
//   string newsubdirname="./";
#else

  // GCOUT << "Done with initialization" << endl;
// how many DRA files per core

// choose some unique string to facilitate deleting these directories

   string base="dra-";
  // const int manual_max = 336; // keep this a multiple of 7 ( or 3) 
   const int manual_max = 840; // keep this a multiple of 7 ( or 3) 

/* Hopper 24 cores per node. 26 OSSs 156 OSTs
   72,000 / 156 -> 460 OSTs (pseudoDisk) ( ~460 cores) / dir -> 3220 files
*/
   //const int manual_max = 7200; // keep this a multiple of 7 ( or 3) 


   //const int manual_max = 2016 //288 cores

// dump files into a subdir. we want numsubdirs << numcores and deterministic

   const int maxfilesperdir=max(manual_max,numfiles); // needs some degree of optimization

   int numsubdirs = max(1,numfiles * g_size / maxfilesperdir);  // numdirs based on files
   GCOUT << "DRA numsubdirs set to " << numsubdirs << endl;

   int width = max(1, g_size/numsubdirs); // total g_ranks per dir

// Now figure out which of these I use
// this fails on the read. Seems to be a screwup in DRA regaridng the trailing number on the files

   int mydir=0;
   int lo=0, hi=0;

   for( int i=0; i< numsubdirs; ++i ) {
      lo = i * width; 
      hi = min(lo + width, g_size);

      if ( g_rank >= lo && g_rank < hi ) {
         mydir = i;
//cout << "in loop " << i<<" "<<width<<" "<<lo<<" "<<hi<< " "<<g_rank<<" "<<mydir<<endl;
      }
//cout << "data " << " "<<numsubdirs << " "<<width<<" "<<lo<<" "<<hi<< " "<<g_rank<<" "<<mydir<<endl;
   }

   std::string s;
   std::stringstream out;
   out << mydir;
   s = out.str();

   string newsubdirname="./"+base+s+"/";

   //cout << "my new dir name is " << newsubdirname << endl;
//   int dirstatus = mkdir(newsubdirname.c_str(),S_IRWXU );
#endif 
// end NOSUBDIR

   if ( GA::nodeid() == 0 ) {
      int dirstatus = mkdir(newsubdirname.c_str(),S_IRWXU );
   }



   string CIarrayname = DRA_CIfilename;
   string CIfilename = newsubdirname+DRA_CIfilename;
   PsociDRArestart d_cimat( g_rank, CIarrayname, CIfilename );
 
   string ICIarrayname = DRA_ICIfilename;
   string ICIfilename = newsubdirname+DRA_ICIfilename;
   PsociDRArestart d_icicol( g_rank, ICIarrayname, ICIfilename );
   
   string NUMarrayname = DRA_NUMfilename;
   string NUMfilename = newsubdirname+DRA_NUMfilename;
   PsociDRArestart d_number( g_rank, NUMarrayname, NUMfilename );
   
   GCOUT << GA::nodeid() <<" cifile " << CIfilename << endl;
#ifdef SUPPLEMENTAL

GCOUT << "SUPP DRA-1" << endl;

   string CIarrayname_supp = DRA_CIfilename_supp;
   string CIfilename_supp = newsubdirname+DRA_CIfilename_supp;
   PsociDRArestart d_cimat_supp( g_rank, CIarrayname_supp, CIfilename_supp );

GCOUT << "SUPP DRA-2" << endl;

   string ICIarrayname_supp = DRA_ICIfilename_supp;
   string ICIfilename_supp = newsubdirname+DRA_ICIfilename_supp;
   PsociDRArestart d_icicol_supp( g_rank, ICIarrayname_supp, ICIfilename_supp );

GCOUT << "SUPP DRA-3" << endl;

   string NUMarrayname_supp = DRA_NUMfilename_supp;
   string NUMfilename_supp = newsubdirname+DRA_NUMfilename_supp;
   PsociDRArestart d_number_supp( g_rank, NUMarrayname_supp, NUMfilename_supp );
#endif

   string DIAGarrayname = DRA_DIAGfilename;
   string DIAGfilename = newsubdirname+DRA_DIAGfilename;
   PsociDRArestart d_diag_sef( g_rank, DIAGarrayname, DIAGfilename );
 
   // Now choose mode CONSTRUCT or RESTART ( or NONE )
   
   int readSparse = 0;
   int readSparse_supp=0;
   int readWidth = 0;
   int ndims = 2;

   int maxwidth=0;

#ifdef SUPPLEMENTAL
   maxwidth = percent_maxwidth * (1+ (maxsef  / 100L)) ; // round upP - hack over the fact that 30*maxsef > 2GB
   GCOUT << " WIDTH is " << percent_maxwidth<<" % of maxsef = "<< maxsef<<" Computed real width is "<< maxwidth<<endl;
#endif

   if ( construct_hamiltonian == CONSTRUCT ) {
     GCOUT << "Preparing to Construct Hamiltonian " << endl;

#ifdef SUPPLEMENTAL
     AllocateDRA( construct_hamiltonian, ndims, maxsef, maxsparse, supp_maxsparse, maxwidth, d_cimat, d_icicol, d_number, d_diag_sef, d_cimat_supp, d_icicol_supp, d_number_supp);
#else
     AllocateDRA( construct_hamiltonian, ndims, maxsef, maxsparse, d_cimat, d_icicol, d_number, d_diag_sef);
#endif

   } else if( construct_hamiltonian == RESTART ) { // Plan on fetching hamiltonian
     GCOUT << "Preparing to READ hamiltonian data from disk " << endl;

#ifdef SUPPLEMENTAL
     int retvalue = SpecifyDRA(construct_hamiltonian, d_cimat,d_icicol,d_number,d_diag_sef, d_cimat_supp, d_icicol_supp, d_number_supp, readSparse, readSparse_supp, readWidth );
#else
     int retvalue = SpecifyDRA(construct_hamiltonian, d_cimat,d_icicol,d_number,d_diag_sef, readSparse);
#endif
     GCOUT << "READ SPARSE was " << readSparse << endl;
   } 
   
   if( construct_hamiltonian == RESTART ) {
     GCOUT << "On Hamiltonian READ values for sparsity and supp values were identified: Changing to match " << endl;
     maxsparse = readSparse;
     supp_maxsparse = min( maxsef, readSparse_supp);
     maxwidth = readWidth;
     GCOUT <<"maxsparse is "<<maxsparse<<" maxsparse_supp is "<< supp_maxsparse<<" maxwidth is "<<maxwidth<<endl;
   }

GCOUT << "Memory System before H " << endl;
if ( GA::nodeid() == 0 ) system ("cat /proc/meminfo | grep active");


#ifdef SUPPLEMENTAL
 PsociGAhamiltonian hamilton(maxsparse, maxsef, supp_maxsparse, maxwidth, &deters_det, &mos );
#else
 PsociGAhamiltonian hamilton(maxsparse, maxsef, &deters_det, &mos );
#endif

 GCOUT << "Setting NXTVAL CHUNK to " << nxtval_chunk << endl;
 hamilton.setNxtvalChunk( nxtval_chunk );
 
 // Extra bit neede by all but data comes from Det processing
 
 int nsize = hamilton.computeNsefEntryPoints( l_nsef );
 //hamilton.printNsefEntryPoints();
 
 hamilton.printKsym();
 hamilton.printPerSpatialValues();
 hamilton.hamiltonianDistribution();

// Compute H

/* These now not used anymore  applies only to the fetch2D code
   hamilton.SetGranularity( l_granularity );
   hamilton.PrintGranularity();
*/

#ifdef GATRACE
   long maxnum = 10000;
   long * numevents = &maxnum;
   trace_init_( numevents );
#endif

 pair<int,double> Hinfo;
 if ( construct_hamiltonian == CONSTRUCT || construct_hamiltonian == NONE ) {


// this doesn;t do anything anymiore

   GCOUT << "Construct Experimental FETCH2D distributed hamiltonian " << endl;
   hamilton.SetGranularity( l_granularity );
   hamilton.PrintGranularity();

   //GCOUT << "Construct Experimental 2D TWOWAY DIRECT distributed hamiltonian " << endl;
   //long numH = hamilton.constructGlobalHamiltonian2WAYDirect(Hinfo);

#ifdef DIRECT

#ifdef NEWORBMAP
GCOUT << "Using replicated orbmap method " << endl;
long numH = hamilton.constructGlobalHamiltonianAltOrbMapDirect( Hinfo );
#else
GCOUT << "Using NEW distributed ChunkJHmap method " << endl;
hamilton.SetJChunk( l_numJchunks );
hamilton.PrintJChunk();
long numH =  hamilton.constructGlobalHamiltonianAltOrbMapDirectChunkJmap( Hinfo );
//cout << "newer INNER TEST NEW H method " << endl;
//long numH =  hamilton.constructGlobalHamiltonianAltOrbMapDirectChunkJmapInner( Hinfo );
#endif

#else

//      long numH = hamilton.constructGlobalHamiltonian2WAY(Hinfo);
#ifdef NEWORBMAP
  GA::error("NEWORBMAP MUST NOT be specified for NONDIRECT" ,1);
#endif
  GCOUT << "Using NEW distributed NONDIRECT ChunkJHmap method " << endl;
hamilton.SetJChunk( l_numJchunks );
hamilton.PrintJChunk();
long numH =  hamilton.constructGlobalHamiltonianAltOrbMapChunkJmap( Hinfo );
GCOUT << "returned frim build H" << endl;

#endif

/*
#ifdef DIRECT
//     GCOUT << "Construct Experimental 2D DIRECT distributed hamiltonian " << endl;
//     long numH = hamilton.constructGlobalHamiltonianExpDirect(Hinfo);
     GCOUT << "Construct Experimental ORBMAP 2D TWOWAY DIRECT distributed hamiltonian " << endl;
     //deleted long numH = hamilton.constructGlobalHamiltonian2WAYDirect(Hinfo);
     long numH = constructGlobalHamiltonianAltOrbMapDirect( Hinfo );
#elif TWOWAY
      GCOUT << "Construct Experimental 2D TWOWAY distributed hamiltonian " << endl;
      long numH = hamilton.constructGlobalHamiltonian2WAY(Hinfo);
#else
     GCOUT << "Construct Experimental 2D distributed EXPALT hamiltonian " << endl;
//     long numH = hamilton.constructGlobalHamiltonianExp(Hinfo);
//    long numH = hamilton.constructGlobalHamiltonianExpAlternative(Hinfo);
#endif
#elif FETCHALL 
   GCOUT << "Construct Experimental get-ALL hamiltonian " << endl;
   hamilton.SetGranularity( -1 ); // -1 triggers a FULL block get one time
   hamilton.PrintGranularity();
//   long numH = hamilton.constructGlobalHamiltonianExpAll(Hinfo);

#ifdef DIRECT
   GCOUT << "Construct DIRECT distributed hamiltonian " << endl;
//   long numH = hamilton.constructGlobalHamiltonianDirect(Hinfo); // THis one works fine
    long numH = hamilton.constructGlobalHamiltonianAltOrbMapDirect(Hinfo);
#else
   GCOUT << "Construct distributed hamiltonian " << endl;
   long numH = hamilton.constructGlobalHamiltonian(Hinfo);
#endif
#endif
*/
   
   GCOUT  << "Total number of elements above threshold  " << numH << endl;
   GCOUT  << "Timing: Me is " << Hinfo.first << " time is " << Hinfo.second << endl;
   
   if ( g_rank == 0 ){
      GA::printStats();
   }

   if ( construct_hamiltonian != NONE ) {

   pair<int,double> dumpTime;
#ifdef SUPPLEMENTAL
   SendHamiltonianToDRA(construct_hamiltonian,d_cimat,d_icicol,d_number,d_diag_sef,d_cimat_supp,d_icicol_supp,d_number_supp, hamilton, dumpTime);
#else
   SendHamiltonianToDRA(construct_hamiltonian,d_cimat,d_icicol,d_number,d_diag_sef, hamilton, dumpTime);
#endif
   GCOUT << dumpTime.first << "DRA reports a dump time of " << dumpTime.second << endl;
   }
  }

 GCOUT << "RE-READ H data form DISK ALWAYS " << endl;

   pair<int,double> readTime;
   construct_hamiltonian = RESTART;
#ifdef SUPPLEMENTAL
   FetchHamiltonianFromDRA(construct_hamiltonian,d_cimat,d_icicol,d_number,d_diag_sef,d_cimat_supp,d_icicol_supp,d_number_supp, hamilton,readTime);
#else
   FetchHamiltonianFromDRA(construct_hamiltonian,d_cimat,d_icicol,d_number,d_diag_sef,hamilton,readTime);
#endif

   GCOUT  << readTime.first << "reports a read DRA time of " << readTime.second << endl;

// Perform RCM (direct): we can use hamilton with changes in-situ
//
// PsociGAhamiltonian hamilton(maxsparse, maxsef, supp_maxsparse, maxwidth, &deters_det, &mos );
   GCOUT << "Perform RCM transformation, send results back to DISK AND into GA space for solution" << endl;


   string testfilename="OrderedNodeRanks.txt";

   // PerformRCM OutOfCOreRCM( g_rank, ICIarrayname, testfilename );
 
   pair<int,double> readICICOLTime;
   construct_hamiltonian = RESTART;

//Warning: in-situ reuse of GA space: Be sure to read FINAL H before eigenanalysis

  //FAILSFetchICICOLFromDRA(construct_hamiltonian,d_icicol,d_icicol_supp, hamilton,readICICOLTime);
  FetchHamiltonianFromDRA(construct_hamiltonian,d_cimat,d_icicol,d_number,d_diag_sef,d_cimat_supp,d_icicol_supp,d_number_supp, hamilton,readTime);
  GA::sync();

   //hamilton.dumpCimat();

   //int junk = hamilton.generateNodeRanking();
   GCOUT << "RCM InMemory"<< endl;
   int junk = hamilton.generateNodeRankingInMemory();





// uncomment me when ready to do a real RCM
//  FetchHamiltonianFromDRA(construct_hamiltonian,d_cimat,d_icicol,d_number,d_diag_sef,d_cimat_supp,d_icicol_supp,d_number_supp, hamilton,readTime);

  GCOUT << "Finished with DRA " << endl;


#ifdef GATRACE
   long * gme = (long *)&g_rank;
   trace_end_( gme );
#endif

/*
   GCOUT << "Manually stop after H" << endl;
   GA::Terminate();
*/

 // Move onto the iterative diagonalization
 // Remember hamiltonian computes the matrix-vector products

 // GA::printStats();

 //if ( g_rank == 0 ) hamilton.printDiags();

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

 if ( numRoots <= 0 ) GA::error(" numRoots <=0 ",numRoots );
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

 GCOUT << " VECTOR distribution " << endl;
 g_v->printDistribution();


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
//       gavectors.readVector( GA::nodeid(), g_v );
       gavectors.readVectorScalable( GA::nodeid(), g_v );

       int restartNumber = gavectors.getNumberVectors();
       GCOUT << "Vector Restart - current number of READ vectors is " << restartNumber << endl;
} else {
      GCOUT << "Build NEW CI coefficients from a diagonal guess " << resname << endl;
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
  vector_set.setProducts( scatter_products );
  GCOUT << "H*v method scatter selection is " << vector_set.getProducts() << endl;

  if ( scatter_products ) {
#ifndef NEWACCUMULATEDSCATTER
     GA::error("Must build code with -DNEWACCUMULATEDSCATTER in order to use scatter_products",1);
#endif
  }

// For performance enhancement of matrixvector products

  if ( maxsef < l_vectorchunk ) {
     GCOUT << "Modification of l_vectorchunk : It is too big for maxsef" << endl;
     l_vectorchunk = maxsef; // rare but we want to be able to do tiny problems
  }

// Reset ( if necc) the value of vector_chunk

  GCOUT << GA::nodeid() << "NEW Original l_vectorchunk is " << l_vectorchunk << endl;
  hamilton.resetNchunk( l_vectorchunk ); // YES IT SEEMS REDUNDANT given we are about to potentially change it

  const int maxmempercore = 1024; // use 1024 for Edison
  l_vectorchunk = vector_set.resetVectorChunk( maxmempercore ) ;
  if ( GA::nodeid()==0 ) cout <<" NEW CHECK:setMatrixChunk " << l_vectorchunk<< endl;

  vector_set.setMatrixVectorChunk( l_vectorchunk );

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
  GCOUT << GA::nodeid() << " Existing products " << vector_set.CurrentExistingHVProducts() << endl;

  // Preread the local cimat, icicol, and number data for each core.

#ifdef PREFETCHHAMILTONIAN
  GCOUT << " PREFETCHING HAMILTONIAN " << endl;
#ifdef MATRIXBALANCE
#ifdef SQUAREHAMILTONIAN
  GA::error(" Mutually exclusive build configuration: -DSQUAREHAMILTONIAN and -DMATRIXBALANCE",1);
#endif
  GCOUT << " PREFETCH LIST CIMAT data " << endl;
  vector_set.generateListCIMAT();

//  vector_set.prefetchListCIMAT();
//  vector_set.destroyGACIMAT(); 


// need to destpoy memory allocate din generateList
#else
//  GCOUT << " PREFETCH CIMAT data " << endl;
//  vector_set.prefetchCIMAT();
#endif
#endif

// we can delete H at this point

  GCOUT << " pre-generate stagger chunk list current vector chuk is " << l_vectorchunk << endl;
  vector_set.generateStaggeredChunkList(l_vectorchunk );
  vector_set.allocateScratchCIMATSpace();
#ifdef NEWACCUMULATEDSCATTER
  vector_set.allocateScratchVectorSpaceScatterAccumulate();
#else
  vector_set.allocateScratchVectorSpace();
#endif

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

#ifdef GATHERHV
 GCOUT << "USING gather-scatter form of the matrixVector product methods" << endl;
#endif

#ifdef MODIFIED_EXPAND
 GCOUT << "USING experimental expandAndReplaceNew method: " << endl;
#endif

#ifdef STAGGERCHUNKLIST
 GCOUT << "MatrixVector routines will not use BRDCST on vector chunks " << endl;
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

// is it possible cores could get out of phase and be in different brdcsts?

    subspace.brdcstEvals( evals ); // rnormSIngle uses local space

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
#ifdef DETAILEDCHECK
    if ( evals.size() < itest ) GA::error("evals is too small ",evals.size() ) ;
#endif
    int status = subspace.computeRnormSingle(itest, evals, rnormvalue );
    rnorms[itest-1] = rnormvalue; // sometimes a root may be converged by chance ( and skipped )
    if ( rnormvalue > current_rtol ) {
    converged = false;
    curroot = itest; //So we could make this an array of rnorms to keep
    break;
    }
    ++itest;
    }
#else
    subspace.computeRnorms( evals, rnorms );
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

#ifdef MODIFIED_EXPAND
        vector_set.expandAndReplaceVectorsNew( subspace.fetchEvecsHandle(), minvecs * numSoughtRoots );
#else
        vector_set.expandAndReplaceVectors( subspace.fetchEvecsHandle(), minvecs * numSoughtRoots );
        matrixVectorProducts += vector_set.CurrentExistingRoots(); // All were recalculated 
#endif
        ++compressSubspace;
      }
      expandTimes += psociTime() - tempTime;
      numUpdates = vector_set.CurrentExistingRoots();
    }

// Expand current subspace vectors into g_vector ( keep subspace at least minvecs * vectors in size )
    if ( vector_set.CurrentExistingRoots() > maxsubspace ) {
      tempTime = psociTime();
#ifdef MODIFIED_EXPAND
      vector_set.expandAndReplaceVectorsNew( subspace.fetchEvecsHandle(), minvecs * numSoughtRoots );
#else
      vector_set.expandAndReplaceVectors( subspace.fetchEvecsHandle(), minvecs * numSoughtRoots );
      matrixVectorProducts += vector_set.CurrentExistingRoots();
#endif
      ++compressSubspace;
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
  subspace.computeRnorms( evals, rnorms );
#endif

  diagTime.first = GA::nodeid();

  for(int i=0; i< numSoughtRoots; ++i ) {
     if ( rnorms[i] <= 0 ) {
     GA::error("An rnorm was missed at i = ", i+1 );
     }
  }
// The final set of solutions must be re-expanded into the full space and compressed to
// numSoughtRoots.

#ifdef MODIFIED_EXPAND
  vector_set.expandAndReplaceVectorsNew( subspace.fetchEvecsHandle(), vector_set.CurrentExistingRoots() );
#else
  vector_set.expandAndReplaceVectors( subspace.fetchEvecsHandle(), vector_set.CurrentExistingRoots() );
  matrixVectorProducts += vector_set.CurrentExistingRoots();
  vector_set.resizeVectors( numSoughtRoots );
#endif

#ifdef RCRMONITOR
  stats =  RCRpopLocation();
  stats =  RCRstopLogging();
#endif

// Destroy local cimat, icicol, and number data

  GCOUT << " DESTROY PREFETCHED LOCAL H data " << endl;
  vector_set.destroyLocalCIMAT();

  GCOUT << " DESTROY STAGGERED data " << endl;
  vector_set.destroyStaggeredChunkList();

// Dump current full set of vectors to disk for possible restart or post processing
  string vecname=vector_filename;
  GCOUT << "Writing " << numSoughtRoots << " vector coefs to disk to the file " << vecname << endl;
  gavectors.openVectorWrite( GA::nodeid(), vecname );
  gavectors.dumpVectorScalable(GA::nodeid(), numSoughtRoots, g_v );
  //gavectors.dumpVector(GA::nodeid(), numSoughtRoots, g_v );
/*
  gavectors.closeVectorFile();
*/
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
  GCOUT << "Number of subspace contractions is " << compressSubspace << endl;

  //GCOUT << "Skip subspace destroy " << endl;
  subspace.destroySubspace();
  
// Need diags for the analysis part - don't bother if numNOroots = 0
 
  if ( numNOroots != 0 ) {

    vector<double> diags;
    hamilton.fetchDiags( diags );

    GCOUT << "Begin CI analysis for numRootsSought = " << vector_set.CurrentNumSoughtRoots() << endl;
    
    pair<double,double> thresh;
    thresh.first = 0.0001;
    thresh.second = 0.0;
#ifdef DUMPREFS
    thresh.second = 0.001; // a non-zero value is sufficient to trigger a ref dump
#endif
    int num3 = analyzeEnergyContributions(thresh, fcore,l_nsef, evals, diags, g_v, vector_set, hamilton, deters_det);
    
    GA::sync();
    GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX finished XXXXXXXXXXXXXXXXXXX " << endl;
    GCOUT << "Building 1e CIDEN for number of roots = "<< numSoughtRoots << endl;

///  vector_set.generateStaggeredChunkList(l_vectorchunk );
  vector_set.destroyScratchVectorSpace();

/*
 Since we are finished we can reuse the existing GA H allocations if we go ahead and zero them out.
 Do not delete the CI vectors or you'll need to reload them
 We should reread the previously dumed vectors to bootstrap the density building.
*/

    hamilton.flushH(); // This zeros out internal arrays. It maintains basic global parameters
    vector<vector<double> > densityCI;
    
// hamilton carries the current Jchunk value in
    hamilton.PrintJChunk();
    PsociNOpopulations population(numNOroots, g_v, &hamilton );
    string title4="test";
    string title5="NOCs";
    int nunit=5;
    int munit=7;
    string fname="mocoef";
    string fnname="nocoef";
 
    population.selectMOs(title4, munit, fname); 
    population.selectNOs(title5, nunit, numNOroots, fnname);

    int numr=population.readMOs();

// Maybe we should combine these: esp fetchOrbitalFileHeader with writeNOheaders
// NOTE only rank==0 does anything with these so arrays are NOT defined for rank <> 0

   population.fetchOrbitalFileHeader();
   string no_title="NaturalOrbitals";

// Write out with the actual coefs instead
// cout << "WRITE OUT THE HEADERS " << endl;
// population.writeNOheaders( no_title );

   pair<int,double> popT;
   int numH = population.populationDensityMatrix( popT );
   GCOUT  << "Total Population Analysis (CI+NOs+Mullikan): Timing: Me is " << popT.first << " time is " << popT.second << endl;

/* At this point we should have dumped ALL NO sets to individual files using the COlumbus library mowrit. Now, 
   the read in MOs and the to-write NOs share much information (nbf, symmetries, etc.) SO we can grab most of that 
   from the existing MOs. The only difference will be the need to add the occupations.
*/

    GA::sync();

    GCOUT << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    GCOUT << " Finished with NO and population analysis " << endl;

  } else { 
    GCOUT << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    GCOUT << "population analysis not performed,. numNOroots == 0 " << endl;
  }
  GA::Terminate(); // gets called by the destructor
}

// A collective call
// Returns current estimate of per-core memory use. Can be torn down befor H construction

long driverFetchConfigurations( PsociConfigs & configs )
{
  int g_rank = GA::nodeid();
  
#ifdef SINGLENODECONFIGS
    configs.readConfigsToGA();
    configs.brdcstParams( 0 ); // only brdcst parameters
    int howmany = configs.getLocalConfigsFromGA();
    configs.printDistributedConfigurations(); 
#else
  if ( g_rank == 0 ) {
    configs.printFilename();
    if ( configs.readConfigs() != 0 ) {
      GA::error("Failed to read configs at numfgs =",-1);
    }
    configs.printParams();
    configs.printConfigurations(); 
  }
  configs.brdcstConfigs( 0 ); //now replicated to all nodes
#endif
  
// GA::sync();
// cout << "done with printout " << endl;
// configs.printParams();

  long maxspatials;
  long nbasis;

  vector<pair<int,string> > word;
  maxspatials = configs.numTotalSpatials();
  nbasis = configs.numBasisFunctions();
  long memory = maxspatials*nbasis*sizeof(int);

#ifdef SINGLENODECONFIGS
  GCOUT << " rank=0 processing of configs " << endl; 
#endif
  GCOUT << " Memory is " << memory<< " Num total spatials is " << maxspatials <<" Num total basis is " << nbasis<< endl;
  return( memory );
}





