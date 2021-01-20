/**************************************************************************************
 * Copyright (c) 2010,2011 RENCI.
 * All rights reserved. This program and the accompanying materials
 * MAY BE available under the terms of the RENCI Open Source License
 * UNC at Chapel Hill which accompanies this distribution, and is available at
 * http://www.renci.org/resources/open-source-software-license

 * New implementation of PSOCI:

 Classes: 

 Description: 
 Estimate the sparsity parameters by initializing G to zero and simply running through all ther H construction
 without saving any of it.
 
 History:

**************************************************************************************/
/**
 *   @file driver.PSOCI.Hestimates.C
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


int main(int argc, char **argv)
{
  parseUserValues(argc, argv);
  
  int g_size;
  int g_rank;

  int memoryStackWords = stack/8;
  int memoryHeapWords  = heap/8;
  int maxGAMemPerCoreWords = maxGAMemPerCore/8;

  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, maxGAMemPerCore );

  if ( GA::nodeid() == 0 ) cout << "maxGAMemPerCore is " << maxGAMemPerCore << endl;
  
  // eed this for MPI-PR on Stampede2
  MPI_Comm comm = GA_MPI_Comm();

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
GCOUT << "get configs" << endl;
  PsociConfigs configs( configs_filename ); //File is opened/closed at the read
GCOUT << " start configs " << endl;
  long configBytes = driverFetchConfigurations( configs );

// Loose all nodes above here
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

   if ( numIOservers == GA::nodes() && GA::nodes() >= 32 ) {
       numIOservers == GA::nodes() / 32; // Im TACC and cluster nodes don;t work here
        if ( GA::nodeid() == 0 ) cout << "cluster nodes don't work: divide by 32 " << endl;
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
      }
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
   maxwidth = 1 + ( percent_maxwidth * (maxsef  / 100L) ) ; // round upP - hack over the fact that 30*maxsef > 2GB
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
     supp_maxsparse = readSparse_supp;
     maxwidth = readWidth;
     GCOUT <<"maxsparse is "<<maxsparse<<" maxsparse_supp is "<< supp_maxsparse<<" maxwidth is "<<maxwidth<<endl;
   }

GCOUT << "Memory System before H " << endl;
if ( GA::nodeid() == 0 ) system ("cat /proc/meminfo | grep active");

GCOUT << "Construct H" << endl;

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
//####long numH =  hamilton.constructGlobalHamiltonianAltOrbMapChunkJmap( Hinfo );
//####GCOUT << "returned frim build H" << endl;
GCOUT << "Get extimates" << endl;

hamilton.constructGlobalHamiltonianAltOrbMapChunkJmapEstimatesOnly();

#endif
 cout << "Exit" << endl;
  GA::Terminate(); // gets called by the destructor
}
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
  
  GA::SERVICES.sync();
  //GA::sync();
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





