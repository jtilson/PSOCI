//TODO delete the PsociHamiltonian objects in the destructor

/**************************************************************************************

* Copyright (c) 2010,2011,2012 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: This method controls access to distributed Hamiltonian matrices
              including their creation ( in parallel ) or read from DRA checkpoints.
              Generally, these are high level distributed operations that call 
              PsociHamiltonian:: to construct/process H blocks

              These are a mix of Local and Collective calls. 

 History:

              Removed many of the checks and balances becasue the GA calls were just
              too expensive on lots of cores

              Using the new preprocessing flag DETAILEDHAMCHECK to enable many
              checks on input data. This should only be needed when modifying the program

              Addey a new preprocessing flag  VISUALIZEHAM to dump non-zero values to local files
              generally this should NOT be used. Furthermore, it only is simple to use in conjunction
              with the STATIC preprocessing flag

              Added the hamiltonian method 2WAY (-DTWOWAY).  This process iconf in chunks of width
              l_nxtval_chunk. Do not go too wide as this will impact load balancing. Keep the value
              between 10 and 100. 

              Big changes, the H and CIDEN methods needed to screen out I-J terms prior to fetching J
              JOUTFG objects.  This was done by building a new GA that only contains the orbital maps (orbmap)
              (size=nbf). Then when an I is selected, the full set of J's that cannot interact are excluded. 
              We still need to do liots of GA gets ( to fetch orbnmaps) but they are small and can be 
              gotten as EAGER messages.As a result, chunking and granularoty are no longer useful.
              
**************************************************************************************/
/**
 *   @file PsociGAhamiltonian.C
 *
 */

// TODO refactor .....

#include <cmath> 
#include <list>


// Not much need for this but perhaps benchmarking against ga:readinc is warrented 

#ifdef TCGMSG 
#include <tcgmsg.h>
#endif

#include "PsociTaskManager.hpp"
#include "PsociTimer.hpp"
#include "PsociDeterminants.hpp"
#include "PsociHamiltonian.hpp"
#include "PsociGAhamiltonian.hpp"
#include "PsociBlasLapack.hpp"

#ifdef RCRMONITOR
extern "C" {
#include "blackboard.h"
}
#endif


extern "C" {
#include "mpi.h"
#ifdef USEMPI
#include "ga-mpi.h"
#endif
}
  

#ifdef GATRACE
extern "C" {
  void trace_init_( long * num );
  void trace_end_( long * me );
  void trace_stime_( void );
  void trace_etime_( void );
  void trace_genrec_( int * ga, int * ilo, int * ihi, int * jlo, int * jhi, int * op );
}
#endif

#ifdef TACC
extern "C" {
#include "typesf2c.h" // This defines Integer which is explicitely used on Hopper
#include "ga-papi.h"
}
#endif

#define NNDAF(i) (i * (i + 1)) / 2 

using namespace std;

/* New constructure that doesnot actually create any GA space: it is intended for the direct matrixvector products 
*/


PsociGAhamiltonian::PsociGAhamiltonian(string & testdir, int maxsparse, int maxsef,int maxsparseBig, int maxwidth, PsociGADeterminants * deters, PsociIntegrals * ints)
{
/*
  l_g_cimat = g_cimat;
  l_g_icicol = g_icicol;
  l_g_number = g_number;
  
  l_g_cimat_supp = g_cimat_supp;
  l_g_icicol_supp = g_icicol_supp;
  l_g_number_supp = g_number_supp;

  l_g_diag_sef = g_diag_sef;
*/

  l_nchunk = default_vector_chunk;

  l_maxsparse = maxsparse;
  l_maxsef = maxsef;
  
  l_maxsparse_supp = maxsparseBig;
  l_maxwidth = maxwidth;
  
  l_deters = deters; // Got GA dets
  
  l_nbf = l_deters->fetchNumBasisFunctions();
  l_nelec = l_deters->fetchNumElectrons();
  l_maxspatials = l_deters->fetchMaxSpatials();
  
  l_gaHamiltonian  = new PsociHamiltonian( l_nbf, l_nelec );
  l_gaHamiltonian->specifyIntegrals( ints ) ; // Got replicated ints
  
  l_gaHamiltonian->setKsym(l_deters->fetchGlobalSymmetry() );
  l_gaHamiltonian->setPerSpatialValues( l_deters->fetchGlobalMaxDetPerSpatial(), l_deters->fetchGlobalMaxSefPerSpatial() );

/* Set the default value for agglomerated fetches */
  l_granularity = 10; // Not demonstrated useful yet.
  l_numJchunksHbuild = 1;

// Need to find where to POPULATE this now. We should probably do it first time through
// We donot want to reengineer the maxtir-vector linear algebra
 

  createGAhamiltonianDiagOnly();

  char handleMessage4[] = "PsociGAhamiltonian:diag_sef";
  l_g_diag_sef->checkHandle( handleMessage4 );
  char handleMessage5[] = "PsociGAhamiltonian:matrix";

// need these as well

   local_dims[0] = maxsef;
   local_dims[1] = maxsparse;

    if (GA::nodeid() == 0 ) cout << "Maximum row prefetching (MAX_AGGREGATE) is "<< MAX_AGGREGATE << endl;
    if (GA::nodeid() == 0 ) cout << "MINIMAP_MAX paramaters is (MINIMAP_MAX) " <<  MINIMAP_MAX << endl;
    if (GA::nodeid() == 0 ) cout << "Max elements before deflation is (MAX_DEFLATE_ELEMENTS) " << MAX_DEFLATE_ELEMENTS << endl;

   if (GA::nodeid() == 0 ) cout << "EXPERIMENTAL: SIZES int versus Integer "<< sizeof(int)<<" "<<sizeof(Integer) << endl;
#ifdef TACC
   if ( sizeof(int) != sizeof(Integer) ) {
      cout << "LowLevel products: int <> Integer: try recompiling with the -DGAEXP flag " << endl;
#ifndef GAEXP 
      GA::error(" int <> Integer and -DGAEXP not specified ",-1);
#endif
   }
#endif
}
// deters MUST be preallocated and filled before getting here...it basically is the entry to GAdets

// Collective but not all need actually do any work.
PsociGAhamiltonian::PsociGAhamiltonian(int maxsparse, int maxsef, PsociGADeterminants * deters, PsociIntegrals * ints )
{
  /*
    l_g_cimat = g_cimat;
    l_g_icicol = g_icicol;
    l_g_number = g_number;
    l_g_diag_sef = g_diag_sef;
  */

  l_nchunk = default_vector_chunk;
  
  l_maxsparse = maxsparse;
  l_maxsef = maxsef;
  
  l_deters = deters; // Got GA dets
  
  l_nbf = l_deters->fetchNumBasisFunctions();
  l_nelec = l_deters->fetchNumElectrons();
  l_maxspatials = l_deters->fetchMaxSpatials();
  
  l_nxtval_chunk = NXTVAL_CHUNK_SIZE; // Only used id doing NXTVAL local balancing

  l_gaHamiltonian  = new PsociHamiltonian( l_nbf, l_nelec ); 
  l_gaHamiltonian->specifyIntegrals( ints ) ; // Got replicated ints
  
  l_gaHamiltonian->setKsym(l_deters->fetchGlobalSymmetry() );
  l_gaHamiltonian->setPerSpatialValues( l_deters->fetchGlobalMaxDetPerSpatial(), l_deters->fetchGlobalMaxSefPerSpatial() );
  
  l_granularity = 10;
  l_numJchunksHbuild = 1;
  
  if (GA::nodeid() == 0 ) cout << "createGAhamiltonian" << endl;
  createGAhamiltonian();
  
  char handleMessage[] = "PsociGAhamiltonian:matrix";
  l_g_cimat->checkHandle( handleMessage );
  
  char handleMessage2[] = "PsociGAhamiltonian:mapping";
  l_g_icicol->checkHandle( handleMessage2 );
  
  char handleMessage3[] = "PsociGAhamiltonian:number";
  l_g_number->checkHandle( handleMessage3 );
  
  char handleMessage4[] = "PsociGAhamiltonian:diag_sef";
  l_g_diag_sef->checkHandle( handleMessage4 );
  
  int type, idim;
  l_g_cimat->inquire( &type, &idim , local_dims );
  l_g_cimat->distribution( GA::nodeid(), local_cilo, local_cihi );
  if ( idim != 2 ) {
    GA::error("erroneous idim for cimat: aborting ",-1);
  } 
  if (GA::nodeid() == 0 ) cout << "finished cimat distribution" << endl;


#ifdef VISUALIZEHAM
/* Extra bit for visualizing H matrix
*/
    string tempString = "hamiltonian";
    int num = GA::nodeid();
    stringstream sout;
    sout << num;
    string scout = sout.str();
    scout.insert(0,"-");
    scout.insert(0,tempString);
    
    string Hname = scout;
    ohamiltonianfile.open( Hname.c_str(), std::ios::out) ;
#endif
}

void PsociGAhamiltonian::destroyGAspace()
{
  if ( GA::nodeid() == 0 ) cout << "destroyGAspace::DESTROY GA-BASED H/ICICOL/SUPP" << endl; 
  l_g_cimat->destroy();
  l_g_icicol->destroy();
  l_g_number->destroy();
  //l_g_diag_sef->destroy();//todo need this for rnorm processing keep me
#ifdef SUPPLEMENTAL
  l_g_cimat_supp->destroy();
  l_g_icicol_supp->destroy();
  l_g_number_supp->destroy();
#endif
}

/*
 PsociGAhamiltonian::PsociGAhamiltonian(int maxsparse, int maxsef, GA::GlobalArray * g_cimat, GA::GlobalArray * g_icicol, GA::GlobalArray * g_number, 
					GA::GlobalArray * g_diag_sef, int maxsparseBig, int maxwidth, GA::GlobalArray * g_cimat_supp, 
					GA::GlobalArray * g_icicol_supp, GA::GlobalArray * g_number_supp,
					PsociGADeterminants * deters, PsociIntegrals * ints )
*/
PsociGAhamiltonian::PsociGAhamiltonian(int maxsparse, int maxsef,int maxsparseBig, int maxwidth, PsociGADeterminants * deters, PsociIntegrals * ints)
{
/*
  l_g_cimat = g_cimat;
  l_g_icicol = g_icicol;
  l_g_number = g_number;
  
  l_g_cimat_supp = g_cimat_supp;
  l_g_icicol_supp = g_icicol_supp;
  l_g_number_supp = g_number_supp;

  l_g_diag_sef = g_diag_sef;
*/

  l_nchunk = default_vector_chunk;

  l_maxsparse = maxsparse;
  l_maxsef = maxsef;
  
  l_maxsparse_supp = maxsparseBig;
  l_maxwidth = maxwidth;
  
  l_deters = deters; // Got GA dets
  
  l_nbf = l_deters->fetchNumBasisFunctions();
  l_nelec = l_deters->fetchNumElectrons();
  l_maxspatials = l_deters->fetchMaxSpatials();
  
  l_gaHamiltonian  = new PsociHamiltonian( l_nbf, l_nelec );
  l_gaHamiltonian->specifyIntegrals( ints ) ; // Got replicated ints
  
  l_gaHamiltonian->setKsym(l_deters->fetchGlobalSymmetry() );
  l_gaHamiltonian->setPerSpatialValues( l_deters->fetchGlobalMaxDetPerSpatial(), l_deters->fetchGlobalMaxSefPerSpatial() );

/* Set the default value for agglomerated fetches */
  l_granularity = 10; // Not demonstrated useful yet.
  l_numJchunksHbuild = 1;
  
  //cout << "call create 2Array max spatials are  " << l_maxspatials <<" "<<  endl;

  if (GA::nodeid() == 0 ) cout << "2D createGAhamiltonian" << endl;
  createGAhamiltonian2Array();

  
  char handleMessage[] = "PsociGAhamiltonian:matrix";
  l_g_cimat->checkHandle( handleMessage );
  char handleMessage2[] = "PsociGAhamiltonian:mapping";
  l_g_icicol->checkHandle( handleMessage2 );
  char handleMessage3[] = "PsociGAhamiltonian:number";
  l_g_number->checkHandle( handleMessage3 );
  char handleMessage4[] = "PsociGAhamiltonian:diag_sef";
  l_g_diag_sef->checkHandle( handleMessage4 );
  char handleMessage5[] = "PsociGAhamiltonian:matrix";
  l_g_cimat_supp->checkHandle( handleMessage5 );
  char handleMessage6[] = "PsociGAhamiltonian:mapping";
  l_g_icicol_supp->checkHandle( handleMessage6 );
  char handleMessage7[] = "PsociGAhamiltonian:number";
  l_g_number_supp->checkHandle( handleMessage7 );


/* May be incomplete but SHOULD be okay. This just means that
   calls that need tio go to the split array may be off-chip
*/

  int type, idim;
  l_g_cimat->inquire( &type, &idim , local_dims );

  if (GA::nodeid() == 0 ) cout << "cimat distribution" << endl;
  l_g_cimat->distribution( GA::nodeid(), local_cilo, local_cihi );
  if ( idim != 2 ) {
    GA::error("erroneous idim for cimat: aborting ",-1);
  }
  
#ifdef VISUALIZEHAM
/* Extra bit for visualizing H matrix
*/
    string tempString = "hamiltonian";
    int num = GA::nodeid();
    stringstream sout;
    sout << num;
    string scout = sout.str();
    scout.insert(0,"-");
    scout.insert(0,tempString);

    string Hname = scout;
    ohamiltonianfile.open( Hname.c_str(), std::ios::out) ;
#endif

    if (GA::nodeid() == 0 ) cout << "Maximum row prefetching (MAX_AGGREGATE) is "<< MAX_AGGREGATE << endl;
    if (GA::nodeid() == 0 ) cout << "MINIMAP_MAX paramaters is (MINIMAP_MAX) " <<  MINIMAP_MAX << endl;
    if (GA::nodeid() == 0 ) cout << "Max elements before deflation is (MAX_DEFLATE_ELEMENTS) " << MAX_DEFLATE_ELEMENTS << endl;

   if (GA::nodeid() == 0 ) cout << "EXPERIMENTAL: SIZES int versus Integer "<< sizeof(int)<<" "<<sizeof(Integer) << endl;
#ifdef TACC
   if ( sizeof(int) != sizeof(Integer) ) {
      cout << "LowLevel products: int <> Integer: try recompiling with the -DGAEXP flag " << endl;
#ifndef GAEXP 
      GA::error(" int <> Integer and -DGAEXP not specified ",-1);
#endif
   }
#endif
}

/*
  PsociGAhamiltonian::PsociGAhamiltonian(int maxsparse, int maxsef, GA::GlobalArray * g_cimat, GA::GlobalArray * g_icicol, GA::GlobalArray * g_number, 
					GA::GlobalArray * g_diag_sef, int maxsparseBig, int maxwidth, GA::GlobalArray * g_cimat_supp, 
					GA::GlobalArray * g_icicol_supp, GA::GlobalArray * g_number_supp,
					PsociGADeterminants * deters, PsociIntegrals * ints )

  EXPERIMENTAL transition moment code
  NOTE: the LHS of the H construction is driven by the configs1 data set: THat way we do not need to change much of the
        index management code
*/
PsociGAhamiltonian::PsociGAhamiltonian(int maxsparse, int maxsef,int maxsparseBig, int maxwidth, PsociGADeterminants * deters, PsociGADeterminants * deters2, PsociIntegrals * ints)
{
/*
  l_g_cimat = g_cimat;
  l_g_icicol = g_icicol;
  l_g_number = g_number;
  
  l_g_cimat_supp = g_cimat_supp;
  l_g_icicol_supp = g_icicol_supp;
  l_g_number_supp = g_number_supp;

  l_g_diag_sef = g_diag_sef;
*/

  l_nchunk = default_vector_chunk;
  
  l_maxsparse = maxsparse;
  l_maxsef = maxsef;
  
  l_maxsparse_supp = maxsparseBig;
  l_maxwidth = maxwidth;
  
  l_deters = deters; // Got GA dets
  l_deters2 = deters2; // Got GA dets (J loop) 
  
  l_nbf = l_deters->fetchNumBasisFunctions();
  l_nelec = l_deters->fetchNumElectrons();
  l_maxspatials = l_deters->fetchMaxSpatials();
  
// Do we need to set ksym for both configs ?

  l_gaHamiltonian  = new PsociHamiltonian( l_nbf, l_nelec );
  l_gaHamiltonian->specifyIntegrals( ints ) ; // Got replicated ints
  
  l_gaHamiltonian->setKsym(l_deters->fetchGlobalSymmetry() );
  l_gaHamiltonian->setPerSpatialValues( l_deters->fetchGlobalMaxDetPerSpatial(), l_deters->fetchGlobalMaxSefPerSpatial() );

/* Set the default value for agglomerated fetches */
  l_granularity = 10; // Not demonstrated useful yet.

//  cout << "Performing a SPLIT-GA hamiltonian creation " << endl;
  
  if (GA::nodeid() == 0 ) cout << "Next createGAhamiltonian" << endl;
  createGAhamiltonian2Array(); // treats 2nd dimension has sparse ( and thus generic )
  
  char handleMessage[] = "PsociGAhamiltonian:matrix";
  l_g_cimat->checkHandle( handleMessage );
  char handleMessage2[] = "PsociGAhamiltonian:mapping";
  l_g_icicol->checkHandle( handleMessage2 );
  char handleMessage3[] = "PsociGAhamiltonian:number";
  l_g_number->checkHandle( handleMessage3 );
  char handleMessage4[] = "PsociGAhamiltonian:diag_sef";
  l_g_diag_sef->checkHandle( handleMessage4 );
  char handleMessage5[] = "PsociGAhamiltonian:matrix";
  l_g_cimat_supp->checkHandle( handleMessage5 );
  char handleMessage6[] = "PsociGAhamiltonian:mapping";
  l_g_icicol_supp->checkHandle( handleMessage6 );
  char handleMessage7[] = "PsociGAhamiltonian:number";
  l_g_number_supp->checkHandle( handleMessage7 );


/* May be incomplete but SHOULD be okay. This just means that
   calls that need tio go to the split array may be off-chip
*/

  int type, idim;
  l_g_cimat->inquire( &type, &idim , local_dims );
  if (GA::nodeid() == 0 ) cout << "next cimat distribution" << endl;
  l_g_cimat->distribution( GA::nodeid(), local_cilo, local_cihi );
  if ( idim != 2 ) {
    GA::error("erroneous idim for cimat: aborting ",-1);
  }

#ifdef VISUALIZEHAM
/* Extra bit for visualizing H matrix
*/
    string tempString = "hamiltonian";
    int num = GA::nodeid();
    stringstream sout;
    sout << num;
    string scout = sout.str();
    scout.insert(0,"-");
    scout.insert(0,tempString);

    string Hname = scout;
    ohamiltonianfile.open( Hname.c_str(), std::ios::out) ;
#endif
}

PsociGAhamiltonian::~PsociGAhamiltonian()
{
   delete l_gaHamiltonian ; 
}

void PsociGAhamiltonian::resetNchunk( int nchunk )
{
     l_nchunk = nchunk;
}

int PsociGAhamiltonian::fetchNchunk()
{
    return( l_nchunk );
}


/* Not always needed but empties (zeros) out existing arrays
   so that they may be reused. A typical use is computation of the
   1e CI DENSITY - in the current version for population analysis, however, 
   we do not actually need to store to GA so this is rather moot.
*/
void PsociGAhamiltonian::flushH()
{
  char handleMessage[] = "PsociGAhamiltonian:matrix";
  l_g_cimat->checkHandle( handleMessage );
  l_g_cimat->zero();
  
  char handleMessage2[] = "PsociGAhamiltonian:mapping";
  l_g_icicol->checkHandle( handleMessage2 );
  l_g_icicol->zero();
  
  char handleMessage3[] = "PsociGAhamiltonian:number";
  l_g_number->checkHandle( handleMessage3 );
  l_g_number->zero();
  
  char handleMessage4[] = "PsociGAhamiltonian:diag_sef";
  l_g_diag_sef->checkHandle( handleMessage4 );
  l_g_diag_sef->zero();
  
#ifdef SUPPLEMENTAL
  char handleMessage5[] = "PsociGAhamiltonian:matrix";
  l_g_cimat_supp->checkHandle( handleMessage5 );
  l_g_cimat_supp->zero();
  
  char handleMessage6[] = "PsociGAhamiltonian:mapping";
  l_g_icicol_supp->checkHandle( handleMessage6 );
  l_g_icicol_supp->zero();
  
  char handleMessage7[] = "PsociGAhamiltonian:number";
  l_g_number_supp->checkHandle( handleMessage7 );
  l_g_number_supp->zero();
#endif
}


/* DISABLED: A value > 0 that specifies how many configurations should be fetched at a time. We currently 
   set it to 1 since screening had a much bigger effect
*/

GA::GlobalArray * PsociGAhamiltonian::fetchCIMAT()
{
  return( l_g_cimat );
}

GA::GlobalArray * PsociGAhamiltonian::fetchICICOL()
{
  return( l_g_icicol );
}

GA::GlobalArray * PsociGAhamiltonian::fetchNUMBER()
{
  return( l_g_number );
}

int PsociGAhamiltonian::fetchNELEC() 
{
  return( l_nelec );
}

GA::GlobalArray * PsociGAhamiltonian::fetchCIMAT_supp()
{
#ifndef SUPPLEMENTAL
  GA::error(" MUST compilte code with -DSUPPLEMENTAL",-1);
#endif
  return( l_g_cimat_supp );
}

GA::GlobalArray * PsociGAhamiltonian::fetchICICOL_supp()
{
#ifndef SUPPLEMENTAL
  GA::error(" MUST compilte code with -DSUPPLEMENTAL",-1);
#endif
  return( l_g_icicol_supp );
}

GA::GlobalArray * PsociGAhamiltonian::fetchNUMBER_supp()
{
#ifndef SUPPLEMENTAL
  GA::error(" MUST compilte code with -DSUPPLEMENTAL",-1);
#endif
  return( l_g_number_supp );
}


GA::GlobalArray * PsociGAhamiltonian::fetchDIAG_SEF()
{
  return( l_g_diag_sef );
}


void PsociGAhamiltonian::SetGranularity( int width ) 
{
  l_granularity = width;
}

void PsociGAhamiltonian::SetJChunk( int numchunks )
{
  l_numJchunksHbuild = numchunks;
}

void PsociGAhamiltonian::PrintJChunk()
{
  if ( GA::nodeid() != 0 ) return;
   cout << "l_numJchunksHbuild = " << l_numJchunksHbuild<< endl;
   return;
}

void PsociGAhamiltonian::PrintGranularity()
{
  if ( GA::nodeid() == 0 ) {
    cout << "Current determinant (maximum) fetch granularity is " << l_granularity << endl;
  }
}

//Local 
void PsociGAhamiltonian::printPerSpatialValues()
{
  //cout << "Node id is " << GA::nodeid() << endl;
  l_gaHamiltonian->printPerSpatialValues();
}

//Local 
void PsociGAhamiltonian::printKsym()
{
  //cout << "Node id is " << GA::nodeid() << endl;
  l_gaHamiltonian->printKsym();
}

//Timer wrapped version
/*
long PsociGAhamiltonian::constructGlobalHamiltonianDirect( pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalHamiltonianDirect();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( status );
}
*/

// Collective: construct Global Hamiltonian using DIRECT joutfg approach.
// Need to carve up configuration space in some reasonable way.

// This model uses a compute-as-needed joutfg method to eliminate calls to GA
// A node computes all contributions for the given i then stores them to GA.


// Collective - use -DNXTVAL or -DSTATIC for load balancing choices
// Set NXTVAL as the default as it is MUCH better
// use -DSUPPLEMENTAL to exploit split-GA storage for memory savings 

/* Build preprocessing flags: This is NOT a fast method

   -DDIRECT -DSUPPLEMENTAL (best)
   -DDIRECT  disables SUPPLEMENTAL storage 

    No need to use this anymore. It is possible if we introduce the orbmap screening
    in the futre this could be a good way to go.May, 2012.
*/

/*
long PsociGAhamiltonian::constructGlobalHamiltonianDirect()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  
  //cout << "I am " << GA::nodeid() << " at the hamiltonian " << endl;
  //GA::SERVICES.sync();
  
  int counter=0;
  long allcount=0;
  
  //  cout << "Using fetch2D to minimize gets " << endl;
  
#ifdef SUPPLEMENTAL
  // Initialize suppval;
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif
  
#ifdef STATIC
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposiiton is contiguous. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL:Direct based load balancing scheme being used " << endl;
  jobs.initNxtask();// Initialize jobs.nxtask
  long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
  long inexttask=-1;
  int ilo = 0;
  int ihi = l_maxspatials;
#endif
  
  pair<int,double> info1;
  
  int numElems=0;
  long totalElems=0;
  int numdone=0;

  int block_status;
  
  // Build Hamiltonian from scratch
  
  double temp, temp2;
  double jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;
  
  JOUTFG fetchi, fetchj;
  PsociDeterminants * tdeter = l_deters->fetchInternalDeters(); 
  tdeter->copyInternalMapToVector();

// Copy internal multimap object to a simpler vector object
// Try destroying the multimap here
  
  for(int iconf = ilo; iconf < ihi; ++iconf ) {  
#ifndef STATIC
    ++inexttask; 
    if ( inexttask == nexttask ) {
#endif
      
      tdeter->computeConfigs( iconf+1, fetchi );
 
      // One configuration at a time for now
      vector<vector<double> > cimat(fetchi.nsefi );
      vector<vector<int> > icicol(fetchi.nsefi );
      vector<int> number(fetchi.nsefi );
      
      temp = psociTime();
      numElems = 0;
      for(int jconf=0; jconf <= iconf; ++jconf)
	{ 
          temp2 = psociTime();
          tdeter->computeConfigs( jconf+1, fetchj );
          fetchTime += psociTime() - temp2;
	  
	  ++allcount;
          temp2 = psociTime() ;
	  block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj, sefsef );
          sociTime += psociTime() - temp2;
          
          temp2 = psociTime();
	  if ( block_status ) {
                   numElems += packSefsef( fetchi, fetchj , sefsef, cimat, icicol, number );
          }
          packTime += psociTime() - temp2;
	}
      if ( numElems > 0 ) { // For this to not be true would be somewhat rare
	totalElems += numElems;
#ifdef SUPPLEMENTAL
	temp2 = psociTime();
	uploadCIMATtoGASupplemental( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
	uploadTime += psociTime() - temp2;
#else
	uploadCIMATtoGA( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
#endif
      }
#ifndef STATIC
      nexttask = jobs.nxtask( g_size );
      jTime += psociTime() - temp; 
    }
#endif

  }
  // This stuff still has meaning to me in the stdout.
  cout << GA::nodeid() << "jTime is " << jTime<<endl;
  cout << GA::nodeid() << "packTime is " << packTime<<endl;
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "uploadTime is " << uploadTime<<endl;
  cout << GA::nodeid() << "fetchTime is " << fetchTime<<endl;
  
  // Print all of them for now since they indicate load balance
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  GA::SERVICES.sync();
  
  tdeter->tearDownDeterminants();// no longer need the replicated CONFIGS/VEC_CONFIGS data
  
#ifdef SUPPLEMENTAL
  //Destroy suppval;
  destroySuppVal();
#endif
  
#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif
  
#ifdef RCRMONITOR
  stats =  RCRpopLocation();
  stats =  RCRstopLogging();
#endif
  
  //return( totalElems ); // switched for performance measurment
  return( numberNonzeroElements( totalElems ) );
}
*/

//Timer wrapped version
/*
long PsociGAhamiltonian::constructGlobalHamiltonian( pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalHamiltonian();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( status );
} 
*/

//Timer wrapped version
/*
long PsociGAhamiltonian::constructGlobalHamiltonianExpDirect( pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalHamiltonianExpDirect();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( status );
} 
*/

// Collective: construct Global Hamiltonian using DIRECT joutfg approach.
// Need to carve up configuration space in some reasonable way.

// This model uses a compute-as-needed joutfg method to eliminate calls to GA
// A node computes all contributions for the given i then stores them to GA.


// Collective - use -DNXTVAL or -DSTATIC for load balancing choices
// Set NXTVAL as the default as it is MUCH better
// use -DSUPPLEMENTAL to exploit split-GA storage for memory savings 
// Build with FETCH2D to get blocks working
/* Build preprocessing flags

   -DDIRECT -DFETCH2D -DSUPPLEMENTAL (best)
   -DDIRECT -DFETCH2D disables SUPPLEMENTAL storage 

   NOTE THis has been disabled as an option from the test "driver" codes.
   We default to using the TWOWAY DIRECT instead. Preprocessing flags
   will not change this
*/

/*
long PsociGAhamiltonian::constructGlobalHamiltonianExpDirect()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  
  int counter=0;
  long allcount=0;
  
  if ( GA::nodeid() == 0 ) cout << "Using computeConfig2D: ExpDirect Hamiltonian " << endl;
  
#ifdef SUPPLEMENTAL
  // Initialize suppval;
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif
  
#ifdef STATIC
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposiiton is contiguous. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL:Direct 2D based load balancing scheme being used " << endl;
  jobs.initNxtask();// Initialize jobs.nxtask
  long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
  long inexttask=-1;
  int ilo = 0;
  int ihi = l_maxspatials;
#endif
  
  pair<int,double> info1;
  
  int numElems=0;
  long totalElems=0;
  int block_status;
  
  // Build Hamiltonian from scratch
  
  double temp, temp2;
  double jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;
  
  JOUTFG fetchi;
  vector<JOUTFG> fetchj;
  
  PsociDeterminants * tdeter = l_deters->fetchInternalDeters(); 
  tdeter->copyInternalMapToVector(); // use tearDown to delete these
  
  for(int iconf = ilo; iconf < ihi; ++iconf ) {  
    
#ifndef STATIC
    ++inexttask; 
    if ( inexttask == nexttask ) {
    // cout << GA::nodeid() << " Doing a nexttask I J are " << iconf <<  " " << inexttask << " " << nexttask << endl;
#endif
      
      tdeter->computeConfigs( iconf+1, fetchi );
 
      vector<vector<double> > cimat( fetchi.nsefi );
      vector<vector<int> > icicol( fetchi.nsefi );
      vector<int> number( fetchi.nsefi );
      
      temp = psociTime();
      numElems = 0;
      for(int jconf=0; jconf <= iconf; jconf += l_granularity)
	{ 
          int kfinal = min( l_granularity, iconf - jconf + 1 );
          temp2 = psociTime();
          fetchj.clear();
          tdeter->computeConfigs( jconf+1, jconf+kfinal, fetchj );
          fetchTime += psociTime() - temp2;

         for(int jconf_block=0; jconf_block < fetchj.size(); ++jconf_block )
          {
	  ++allcount;
          temp2 = psociTime() ;
	  block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj[jconf_block], sefsef );
          sociTime += psociTime() - temp2;
          
          temp2 = psociTime();
	  if ( block_status ) {
                   numElems += packSefsef(fetchi, fetchj[jconf_block] , sefsef, cimat, icicol, number );
          }
          packTime += psociTime() - temp2;
	}
       }
      
      if ( numElems > 0 ) { // For this to not be true would be somewhat rare
        totalElems += numElems;
#ifdef SUPPLEMENTAL
        temp2 = psociTime();
	uploadCIMATtoGASupplemental( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
        uploadTime += psociTime() - temp2;
#else
	uploadCIMATtoGA( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
#endif
      }
#ifndef STATIC
      nexttask = jobs.nxtask( g_size );
      jTime += psociTime() - temp;
    }
#endif
    
  }
  cout << GA::nodeid() << "jTime is " << jTime<<endl;
  cout << GA::nodeid() << "packTime is " << packTime<<endl;
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "uploadTime is " << uploadTime<<endl;
  cout << GA::nodeid() << "fetchTime is " << fetchTime<<endl;

  // Print all of them for now since they indicate load balance
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  GA::SERVICES.sync();
  
#ifdef SUPPLEMENTAL
  //Destroy suppval;
  destroySuppVal();
#endif
  
#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif
  
#ifdef RCRMONITOR
  stats =  RCRpopLocation();
  stats =  RCRstopLogging();
#endif
  
  //return( totalElems ); // switched for performance measurment
  return( numberNonzeroElements( totalElems ) );
}
*/

// Collective: construct Global Hamiltonian.
// Need to carve up configuration space in some reasonable way.

/* It is possible that a node has no config data. This should not cause an error BUT
   it should be flagged as a potential scalability problem.
   
   A fetch of a config is of a variable size depending on the alpha/beta list lengths etc. Those details
   are handled in the underlying layer.
   
   A node computes all contributions for the given i then stores them to GA. This is not a 
   bad compile time choice. THe only constraint is that one cannot is the SUPPLEMENTAL storage
   model.

   Collective - use -DNXTVAL or -DSTATIC for load balancing choices
   Set NXTVAL as the default as it is MUCH better
   use -DSUPPLEMENTAL to exploit split-GA storage for memory savings 

   This version has been superceded by thatg using orbmaps for prescreening. For now, 
   we simply comments out this code ( below) and inserted the replacement code subsequently.
   May 2012.
*/

/*
/////////////////////////////////////////////////////////////// REMOVE THE BELOW CODE EVENTUALLY

long PsociGAhamiltonian::constructGlobalHamiltonian()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  
#ifdef PREALLOCATE
  int l_maxlength = l_deters->fetchMaxLength();
  vector<COMPRESS_SIZE> buffer( l_maxlength );
#endif

  //cout << "I am " << GA::nodeid() << " at the hamiltonian " << endl;
  //GA::SERVICES.sync();
  
#ifdef RCRMONITOR
  if (GA::nodeid() == 0 ) cout << "Enabled RCR logging at the Hamiltionian " << endl;
  char * logfilename = "logger.txt";
  int stats =  RCRstartLogging( logfilename );
  
  char * location = "CnsrtH";
  stats = RCRpushLocation(location);
#endif
  
  int counter=0;
  long allcount=0;
  
  //  cout << "Using fetch2D to minimize gets " << endl;
  
#ifdef SUPPLEMENTAL
  // Initialize suppval;
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif
  
#ifdef STATIC
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposiiton is contiguous. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used " << endl;
  jobs.initNxtask();// Initialize jobs.nxtask
  long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
  long inexttask=-1;
  int ilo = 0;
  int ihi = l_maxspatials;
#endif
  
  if ( g_rank == 0 ) cout << " EXPERIMENTAL ORBMAP screening code " << endl;
  vector<COMPRESS_SIZE> imap( l_nbf+2);
  vector<COMPRESS_SIZE> jmap( l_nbf+2);
  int countskip=0, countrun=0; // currently only used for statistics

  pair<int,double> info1;
  
  //   Put the fetch configs process here.......
  //   for(int iconf=0; iconf< l_maxspatials; ++iconf ) {
  //   for(int iconf=g_rank; iconf< l_maxspatials; iconf += g_size ) { 
  
  int numElems=0;
  long totalElems=0;

  int block_status;
  
  // Build Hamiltonian from scratch
  // for(int iconf=ihi; iconf>= ilo; --iconf ) { --iconf is better for static but since we rarely use it ......
  // No need to block on I terms for now
  
  double temp, temp2;
  double mapTime=0.0, jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;

  JOUTFG fetchi, fetchj;
#ifdef PREALLOCATE
  l_deters->preallocateTempJoutfgObject( fetchi );
  l_deters->preallocateTempJoutfgObject( fetchj );
#endif

  int numorbmap;
  for(int iconf = ilo; iconf < ihi; ++iconf ) {  
    
#ifndef STATIC
    ++inexttask; 
    if ( inexttask == nexttask ) {
#endif
      
      l_deters->fetchOrbMap( iconf+1, imap); // for pair-wise comparison

#ifdef PREALLOCATE
      l_deters->fetchAndUnpackDeterminantDataPreallocate( iconf+1, info1, buffer, fetchi );
#else
      l_deters->fetchAndUnpackDeterminantData( iconf+1, info1, fetchi );
#endif
      // One configuration at a time for now
      vector<vector<double> > cimat(fetchi.nsefi );
      vector<vector<int> > icicol(fetchi.nsefi );
      vector<int> number(fetchi.nsefi );
      
      vector<JOUTFG>::iterator jit;
      temp = psociTime();
      numElems = 0;


      for(int jconf=0; jconf <= iconf; ++jconf)
	{ 
          temp2 = psociTime();
          l_deters->fetchOrbMap( jconf+1, jmap);
          numorbmap = l_deters->compareOrbMap( imap, jmap );
          mapTime += psociTime() - temp2;

          if ( numorbmap != 0 ) {
          ++countrun;

          temp2 = psociTime();
#ifdef PREALLOCATE
	  l_deters->fetchAndUnpackDeterminantDataPreallocate( jconf+1, info1, buffer, fetchj );
#else
          l_deters->fetchAndUnpackDeterminantData( jconf+1, info1, fetchj );
#endif

          fetchTime += psociTime() - temp2;

	  ++allcount;
          temp2 = psociTime() ;
	  block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj, sefsef );
          sociTime += psociTime() - temp2;
          
          temp2 = psociTime();
	  if ( block_status ) {
                   numElems += packSefsef( fetchi, fetchj , sefsef, cimat, icicol, number );
#ifdef VISUALIZEHAM
               ohamiltonianfile << fetchi.index<< " " <<fetchj.index<< " " <<numElems<< endl;
#endif
          }
          packTime += psociTime() - temp2;
          } else {
            ++countskip;
          }
	}
#ifdef VISUALIZEHAM
       ohamiltonianfile << endl;
#endif
      
      if ( numElems > 0 ) { // For this to not be true would be somewhat rare
        totalElems += numElems; 
#ifdef SUPPLEMENTAL
        temp2 = psociTime();
	uploadCIMATtoGASupplemental( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
        uploadTime += psociTime() - temp2;
#else
	uploadCIMATtoGA( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
#endif
      }
#ifndef STATIC
      nexttask = jobs.nxtask( g_size );
      jTime += psociTime() - temp;
    }
#endif
    
  }
  cout << GA::nodeid() << "countrun is " << countrun << endl;
  cout << GA::nodeid() << "countskip is " << countskip << endl;
  cout << GA::nodeid() << "jTime is " << jTime<<endl;
  cout << GA::nodeid() << "packTime is " << packTime<<endl;
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "uploadTime is " << uploadTime<<endl;
  cout << GA::nodeid() << "fetchTime is " << fetchTime<<endl;
  cout << GA::nodeid() << "maptime is " << mapTime<<endl;


  // Print all of them for now since they indicate load balance
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  GA::SERVICES.sync();
  
#ifdef SUPPLEMENTAL
  //Destroy suppval;
  destroySuppVal();
#endif
  
#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif
  
#ifdef RCRMONITOR
  stats =  RCRpopLocation();
  stats =  RCRstopLogging();
#endif
  
  //return( totalElems ); // switched for performance measurment
  return( numberNonzeroElements( totalElems ) );
}

////////////////////////////////////////////////////////////////// END REMOVAL OF OLDER H code
*/

								 
//Timer wrapped version
/*
long PsociGAhamiltonian::constructGlobalHamiltonianExp( pair<int,double> & info )
{ 
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalHamiltonianExp();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid(); 
  return( status );
} 
*/


// Collective: construct Global Hamiltonian.- REPLACEMENT version

/* It is possible that a node has no config data. This should not cause an error BUT
   it should be flagged as a potential scalability problem.
   
   A fetch of a config is of a variable size depending on the alpha/beta list lengths etc. Those details
   are handled in the underlying layer.
   
   A node computes all contributions for the given i then stores them to GA. This is not a 
   bad compile time choice. 

   Collective - use -DNXTVAL or -DSTATIC for load balancing choices
   Set NXTVAL as the default as it is MUCH better

   This version superceds by including orbmaps for prescreening. For now, 
   May 2012.

   Recommended compiler settings:
   -DNXTVAL -DPREALLOCATE
   Related sugestions include -DALTERNATIVEDENSITY

   -DRCRMONITOR is highly experimental at this time
*/

// Collective - use -DNXTVAL or -DSTATIC for load balancing choices
// Set NXTVAL as the default as it is MUCH better

/*
long PsociGAhamiltonian::constructGlobalHamiltonian()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  
#ifdef PREALLOCATE
  int l_maxlength = l_deters->fetchMaxLength();
  vector<COMPRESS_SIZE> buffer( l_maxlength );
#endif
  
#ifdef RCRMONITOR
  if (GA::nodeid() == 0 ) cout << "Enabled RCR logging at the Hamiltionian " << endl;
  char * logfilename = "logger.txt";
  int stats =  RCRstartLogging( logfilename );
  
  char * location = "CnsrtH";
  stats = RCRpushLocation(location);
#endif
  
  int counter=0;
  long allcount=0;
  
#ifdef SUPPLEMENTAL
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif
  
#ifdef STATIC
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposiiton is contiguous. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used " << endl;
  jobs.initNxtask();// Initialize jobs.nxtask
  long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
  long inexttask=-1;
  int ilo = 0;
  int ihi = l_maxspatials;
#endif
  
#ifdef PREALLOCATE
  if ( g_rank == 0 ) cout << " Using PREALLOCATED fetch methods " << endl;
#endif

  if ( g_rank == 0 ) cout << " EXPERIMENTAL ORBMAP screening code " << endl;

#ifndef NEWORBMAP
  vector<COMPRESS_SIZE> imap( l_nbf+2);
#endif

  int countskip=0, countrun=0; // currently only used for statistics

  pair<int,double> info1;
  
  //   Put the fetch configs process here.......
  //   for(int iconf=0; iconf< l_maxspatials; ++iconf ) {
  //   for(int iconf=g_rank; iconf< l_maxspatials; iconf += g_size ) { 
  
  int numElems=0;
  long totalElems=0;

  int block_status;
  
  // Build Hamiltonian from scratch
  // for(int iconf=ihi; iconf>= ilo; --iconf ) { --iconf is better for static but since we rarely use it ......
  // No need to block on I terms for now
  
  double temp, temp2;
  double mapTime=0.0, jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;

  JOUTFG fetchi, fetchj;
#ifdef PREALLOCATE
  l_deters->preallocateTempJoutfgObject( fetchi );
  l_deters->preallocateTempJoutfgObject( fetchj );
#endif

  int numorbmap;
  int countfetch=0;
  
#ifdef NEWORBMAP
  if ( GA::nodeid() ==0 ) cout << "USING NEW GOP based ORBMAP approach " << endl;
  vector<vector<pair<short,short> > > ijpairlist( l_maxspatials );
  int listsize = l_deters->fetchAllCompressedOrbMap( ijpairlist ); // Collective: performs a BRDCST underneath
   if ( GA::nodeid() ==0 ) cout << "replicated orbmap size is " << listsize << endl;
#endif
  
  for(int iconf = ilo; iconf < ihi; ++iconf ) {  
    
#ifndef STATIC
    ++inexttask; 
    if ( inexttask == nexttask ) {
#endif  

#ifndef NEWORBMAP
      l_deters->fetchOrbMap( iconf+1, imap); // for pair-wise comparison
#endif
      
#ifdef PREALLOCATE
      l_deters->fetchAndUnpackDeterminantDataPreallocate( iconf+1, info1, buffer, fetchi );
#else
      l_deters->fetchAndUnpackDeterminantData( iconf+1, info1, fetchi );
#endif
      ++countfetch;

      // One configuration at a time for now
      vector<vector<double> > cimat(fetchi.nsefi );
      vector<vector<int> > icicol(fetchi.nsefi );
      vector<int> number(fetchi.nsefi );
      
      vector<JOUTFG>::iterator jit;
      temp = psociTime();
      numElems = 0;
      
      temp2 = psociTime();
      vector<int> list;
#ifdef NEWORBMAP
      countrun +=  l_deters->compareAllOrbMapRemoveZeros(iconf+1, ijpairlist, list ); // new methods
#else
      countrun +=  l_deters->compareAllOrbMap(iconf+1, imap, list ); // new methods
#endif
      mapTime = psociTime() - temp2;
      allcount += iconf+1;
      
      vector<COMPRESS_SIZE>::iterator jit2;
      for(jit2=list.begin(); jit2!=list.end(); ++jit2) {
	int jconf = (*jit2);

	temp2 = psociTime();
#ifdef PREALLOCATE
	  l_deters->fetchAndUnpackDeterminantDataPreallocate( jconf+1, info1, buffer, fetchj );
#else
          l_deters->fetchAndUnpackDeterminantData( jconf+1, info1, fetchj );
#endif
          ++countfetch;

          fetchTime += psociTime() - temp2;
	  
          temp2 = psociTime() ;
	  block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj, sefsef );
          sociTime += psociTime() - temp2;
          
          temp2 = psociTime();
	  if ( block_status ) {
	    numElems += packSefsef( fetchi, fetchj , sefsef, cimat, icicol, number );
          }
          packTime += psociTime() - temp2;
      } 
      if ( numElems > 0 ) { // For this to not be true would be somewhat rare
        totalElems += numElems; 

#ifdef SUPPLEMENTAL
        temp2 = psociTime();
	uploadCIMATtoGASupplemental( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
        uploadTime += psociTime() - temp2;
#else
	uploadCIMATtoGA( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
#endif

	//#endif
      }
#ifndef STATIC
      nexttask = jobs.nxtask( g_size );
      jTime += psociTime() - temp;
    }
#endif
    
  }
  cout << GA::nodeid() << " countrun is " << countrun << endl;
  cout << GA::nodeid() << " countskip is " << countskip << endl;
  cout << GA::nodeid() << " allcount is " << allcount << endl;
  cout << GA::nodeid() << " countfetches is " << countfetch << endl;
  cout << GA::nodeid() << " jTime is " << jTime<<endl;
  cout << GA::nodeid() << " packTime is " << packTime<<endl;
  cout << GA::nodeid() << " sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << " uploadTime is " << uploadTime<<endl;
  cout << GA::nodeid() << " fetchTime is " << fetchTime<<endl;
  cout << GA::nodeid() << " maptime is " << mapTime<<endl;
  
  // Print all of them for now since they indicate load balance

  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  GA::SERVICES.sync();
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif
  
#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif
  
#ifdef RCRMONITOR
  stats =  RCRpopLocation();
  stats =  RCRstopLogging();
#endif
  
  return( numberNonzeroElements( totalElems ) );
}
*/


/* Another experimental version of H construction that was less than useful
   Keep the code for future reference.
*/
/*
long PsociGAhamiltonian::constructGlobalHamiltonianExp()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  
  cout << "I am " << GA::nodeid() << " at the exp SPLIT hamiltonian: This code is antiquated and should be used " << endl;
  GA::Terminate();

  GA::SERVICES.sync();
  
  int counter=0;
  long allcount=0;
  
  double temp, temp2;
  double jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;
  
#ifdef PREALLOCATE
  int l_maxlength = l_deters->fetchMaxLength();
  vector<COMPRESS_SIZE> buffer( l_granularity * l_maxlength );
#endif

#ifdef SUPPLEMENTAL
  // Initialize suppval;
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif
  
#ifdef STATIC
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposiiton is contiguous. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used " << endl;
  jobs.initNxtask();// Initialize jobs.nxtask
  long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
  long inexttask=-1;
  int ilo = 0;
  int ihi = l_maxspatials;
#endif

#ifdef PREALLOCATE
  if ( g_rank == 0 ) cout << " USING PREALLOCATED FETCH METHODS " << endl;
#endif
  
  pair<int,double> info1;
  
  // Preallocate space for fetching determinant information from GA
  
  JOUTFG fetchi;
  
#ifdef PREALLOCATE
  l_deters->preallocateTempJoutfgObject( fetchi );
  vector<JOUTFG> fetchj( l_granularity, fetchi );
#else
  vector<JOUTFG> fetchj;
#endif

  //   Put the fetch configs process here.......
  //   for(int iconf=0; iconf< l_maxspatials; ++iconf ) {
  //   for(int iconf=g_rank; iconf< l_maxspatials; iconf += g_size ) { 
  
  int numElems=0;
  long totalElems=0;
  int block_status;
  
  for(int iconf = ilo; iconf < ihi; ++iconf ) {  
    
//I will do the set of nxtval_chunk at a TIME if Dynamic

#ifndef STATIC
    ++inexttask;
    if ( inexttask == nexttask ) {
//    cout << GA::nodeid() << " Doing a nexttask I J are " << iconf <<  " " << inexttask << " " << nexttask << endl;
#endif
      
    //  fetchi.resize(1); 
      
#ifdef PREALLOCATE
      l_deters->fetchAndUnpackDeterminantDataPreallocate( iconf+1, info1, buffer, fetchi );
#else
      l_deters->fetchAndUnpackDeterminantData( iconf+1, info1, fetchi );
#endif
      
      // One configuration at a time for now
      vector<vector<double> > cimat(fetchi.nsefi );
      vector<vector<int> > icicol(fetchi.nsefi );
      vector<int> number(fetchi.nsefi );
      
      //Need to block these gets up into larger groups l_granularity in width
      
      numElems = 0;
      for(int jconf=0; jconf <= iconf; jconf += l_granularity )  
	{ 
	  temp2 = psociTime();
	  int kfinal = min( l_granularity, iconf - jconf + 1 );
#ifdef PREALLOCATE
	  l_deters->fetchAndUnpackDeterminantData2DPreallocate( jconf+1, jconf + kfinal, buffer, fetchj );
#else
	  fetchj.clear();
	  l_deters->fetchAndUnpackDeterminantData2D( jconf+1, jconf + kfinal, fetchj );
#endif
	  fetchTime += psociTime() - temp2;
	  
	  // We use kfinal because the Preallocated data sets are always the same SIZE
	  
#ifdef PREALLOCATE
	  for(int jconf_block=0; jconf_block < kfinal; ++jconf_block )
	    {
	      ++allcount;
	      temp2 = psociTime();
	      block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj[jconf_block], sefsef );
	      sociTime += psociTime() - temp2;
	      
	      if ( block_status ) {
		numElems += packSefsef( fetchi, fetchj[jconf_block] , sefsef, cimat, icicol, number );
	      }
	    }
#else
	  for(int jconf_block=0; jconf_block<fetchj.size(); ++jconf_block )
	    {
	      ++allcount;
	      temp2 = psociTime();
	      block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj[jconf_block], sefsef );
	      sociTime += psociTime() - temp2;
	      if ( block_status ) {
		numElems += packSefsef( fetchi, fetchj[jconf_block] , sefsef, cimat, icicol, number );
	      }
	    }
#endif
	}
      
      if ( numElems > 0 ) { // For this to not be true would be somewhat rare
	totalElems += numElems;
#ifdef SUPPLEMENTAL
	uploadCIMATtoGASupplemental( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
#else
	uploadCIMATtoGA( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
#endif
      }
#ifndef STATIC
      nexttask = jobs.nxtask( g_size );
    }
#endif
  }
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "fetchTime is " << fetchTime<<endl;
  
  // Print all of them for now since they indicate load balance
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  GA::SERVICES.sync();
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif
  
#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif
  
#ifdef RCRMONITOR
  stats =  RCRpopLocation();
  stats =  RCRstopLogging();
#endif
  
  return( numberNonzeroElements( totalElems ) );
}
*/

//Timer wrapped version
long PsociGAhamiltonian::constructGlobalHamiltonian2WAY( pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalHamiltonian2WAY();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid(); 
  return( status );
} 


/* A 2D H update scheme. Several more layers of indexing are required BUT at the benefit
   of many fewer calls to fetch
   
   BEST compiler options are : -DTWOWAY -DPREALLOCATE -DFETCH2D -DSUPPLEMENTAL - more local memory 
                2nd BEST are : -DTWOWAY -DFETCH2D -DSUPPLEMENTAL          
                     min are : -DTWOWAY -DFETCH2D - looses SPLIT-GA storage


   Not full integrated with the new code.It does work but doesn't perform the orbmazp screening. This needs to 
   be modified moving-forward, May 2012.

*/
/* DISABLED for now

long PsociGAhamiltonian::constructGlobalTransitionHamiltonian2WAY()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  if ( g_rank == 0 ) cout << "I am " << GA::nodeid() << " at the fetch2d 2WAY-BLOCKED SPLIT TRANSITION hamiltonian " << endl;
   
  int counter=0;
  long allcount=0;

  double temp1, temp2;
  double jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;
  
#ifdef PREALLOCATE
  int l_maxlength = l_deters->fetchMaxLength();
  vector<COMPRESS_SIZE> buffer( max(l_granularity,l_nxtval_chunk) * l_maxlength );
#endif
  
#ifdef SUPPLEMENTAL
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif

#ifdef STATIC
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposition is contiguous and inclusive. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used " << endl;
  jobs.initNxtask();
  long nexttask = jobs.nxtaskAlt( g_size ); // Keep J's agglomerated and dynamically balance I for now
  long inexttask=-1;
  int ilo = 0;
  int ihi = l_maxspatials;

#endif
#ifdef PREALLOCATE
  if ( g_rank == 0 ) cout << " USING PREALLOCATED FETCH METHODS " << endl;
#endif

#ifdef PREALLOCATE
  JOUTFG temp;
  l_deters->preallocateTempJoutfgObject( temp );
  l_deters2->preallocateTempJoutfgObject( temp );
  vector<JOUTFG> fetchi( l_nxtval_chunk, temp );
  vector<JOUTFG> fetchj( l_granularity, temp );
#else
  vector<JOUTFG> fetchi;
  vector<JOUTFG> fetchj;
#endif

  pair<int,double> info1;
  
  int numElems=0;
  long totalElems=0;
  int block_status;
  
  int totalFetches=0; // Used for Tuning studies
  
  // Build Hamiltonian from scratch
  
  //   Increment the loop by l_nxtval_chunk within the loop.  This way
  //   we allow nxttask to dynamically determine a staritng iconf for a 
  //   contiguous compute set
  //
  
// Just in case afterall we do an accumulate 
  flushH(); // this simply zeros out cimat,icicol, num and the supplimental spaces

// cout << "ilo and ihi values " << GA::nodeid() <<" " <<ilo<<" "<<ihi<<endl;

  int countme=0;
  
  for(int iconf = ilo; iconf < ihi; ++iconf ) {
    
#ifndef STATIC
    if ( nexttask == iconf ) { // in the loop stride is implemented 
      
#endif
      int ifinal = min( l_nxtval_chunk, ihi - iconf ); // ihi already adds a ONE
      vector<vector<vector<double> > > cimat(ifinal);
      vector<vector<vector<int> > > icicol( ifinal);
      vector<vector<int> > number(ifinal);
      
      temp2 = psociTime();
#ifdef PREALLOCATE
      l_deters->fetchAndUnpackDeterminantData2DPreallocate( iconf+1, iconf+ifinal, buffer, fetchi );
#else
      fetchi.clear();
      l_deters->fetchAndUnpackDeterminantData2D( iconf+1, iconf+ifinal, fetchi );
#endif

      fetchTime += psociTime() - temp2;
      ++totalFetches;

      vector<int> iindex( ifinal); // populate absolute I index array - added code but of fixed work;
      for(int ii=0; ii< ifinal; ++ii ) {
        iindex[ii]=iconf+ii;
      }

      vector<int> numElems(ifinal, 0);

      for(int iconf_block=0; iconf_block < ifinal; ++iconf_block )
        {
          int nsefi = fetchi[iconf_block].nsefi;
              if ( cimat[iconf_block].size() != nsefi ) {
                cimat[iconf_block].resize(nsefi );
                icicol[iconf_block].resize(nsefi );
                number[iconf_block].resize(nsefi );
              }
        }

      for(int jconf=0; jconf < iconf+ifinal ; jconf += l_granularity )
        {
          int kfinal = min( l_granularity, (iconf+ifinal-1) - jconf + 1  );
          temp2 = psociTime();
#ifdef PREALLOCATE
          l_deters2->fetchAndUnpackDeterminantData2DPreallocate( jconf+1, jconf + kfinal, buffer, fetchj );
#else
          fetchj.clear();
          l_deters2->fetchAndUnpackDeterminantData2D( jconf+1, jconf + kfinal, fetchj );
#endif
          fetchTime += psociTime() - temp2;
          ++totalFetches; // comment me out eventually
	  
          vector<int> jindex( kfinal); // populate absolute J index array
          for(int jj=0; jj< kfinal; ++jj ) {
            jindex[jj]=jconf+jj;
          }
	  
          //SOCI always clears sefsef
	  
          for(int iconf_block=0; iconf_block < ifinal; ++iconf_block ) //always ONE for now
            {
              for(int jconf_block=0; jconf_block < kfinal; ++jconf_block )
                {
                  if ( jindex[jconf_block] <= iindex[iconf_block] ) {
                    ++allcount;
                    temp2 = psociTime();
                    block_status = l_gaHamiltonian->generateSOCIblocks( fetchi[iconf_block], fetchj[jconf_block], sefsef );
                    sociTime += psociTime() - temp2;
		    
                    if ( block_status ) {
                      int numTest = packSefsef( fetchi[iconf_block], fetchj[jconf_block] , sefsef, cimat[iconf_block], icicol[iconf_block], number[iconf_block] );
                      numElems[iconf_block] += numTest;
                    }
                  }
                }
            }
        } // jconf
      
      // Update CIMAT. For a given I we need ALL contributing J values 
      //      int iconf_test=iconf;
      for(int iconf_block=0; iconf_block < ifinal; ++iconf_block )
        {
          if ( numElems[iconf_block] > 0 ) { // For this to not be true would be somewhat rare
	    //            ++iconf_test;
	    //            cout << "iconf " << iconf_test<<" "<<numElems[iconf_block];
	    
            totalElems += numElems[iconf_block];
#ifdef SUPPLEMENTAL
            uploadCIMATtoGASupplemental( l_nsef[ fetchi[iconf_block].index - 1 ], cimat[iconf_block], icicol[iconf_block], number[iconf_block]);
#else
            uploadCIMATtoGA( l_nsef[ fetchi[iconf_block].index - 1 ], cimat[iconf_block], icicol[iconf_block], number[iconf_block]);
#endif
          }
        }
#ifndef STATIC
      nexttask = jobs.nxtaskAlt( g_size );
      iconf = iconf + (l_nxtval_chunk - 1); // skip to next possible value from jobs.nxtask
    }
#endif
    
  } // iconf
  
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "fetchTime is " << fetchTime<<endl;
  cout << "allcount is " << allcount<<endl; //this is correct
  cout << "total Fetches is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif
  
#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif
  
  //GA::SERVICES.sync();
  //cout << "LEAVING H construction " << endl;
  return( numberNonzeroElements( totalElems ) );
}
*/





/* A 2D Split H update scheme. Several more layers of indexing are required. NOTE This is a modified version.
   Relative to the original, J granularity is no longer used ( it is ignored). nxtval_chunk is not used much
   but initial performance suggests it to be unneccessary.

   This version uses orbmap to prescreen J configurations for a given I. PREALLOCATE seems to work on blueridge
   but needs to be checked n other platforms

   For now we'll keep the -DFETCH2D requirement but moving forward it can be removed. May 2012.
   
   BEST compiler options are :  -DPREALLOCATE -DSUPPLEMENTAL - more local memory 
                2nd BEST are :  -DSUPPLEMENTAL          


   THIS is the CORRECT production H method to use.

*/

/* This is the only method suitable for getting details Hamiltonian images. That is because
   one MUST compile PSOCI as nondirect AND not using  SINGLENODECONFIGS . One could BUT
   having the spatial configurations blocked by symmetry is helpful

   SINCE this is not a recommend H builder, go ahead and simply add the clunky stuff needed
   to collect and store the data
*/
long PsociGAhamiltonian::constructGlobalHamiltonian2WAY()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  if ( g_rank == 0 ) cout << "I am " << GA::nodeid() << " at the ORBMAP prescreening fetch2d 2WAY-BLOCKED SPLIT hamiltonian " << endl;
  
#ifdef SQUAREHAMILTONIAN
  if ( g_rank == 0 ) cout << "I am " << GA::nodeid() << " SQUARE HAMILTONIAN " << endl; 
#endif

  int counter=0;
  long allcount=0;
  
  double temp1, temp2, temp3;
  double mapTime=0.0, jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;
  
#ifdef PREALLOCATE
  int l_maxlength = l_deters->fetchMaxLength();
  vector<COMPRESS_SIZE> buffer( l_nxtval_chunk * l_maxlength );
#endif
  
#ifdef SUPPLEMENTAL
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif

#ifdef STATIC
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposition is contiguous and inclusive. We add ONE so that subsequent loop argument can be '<' instead of '<='
#endif
#ifndef STATIC
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used " << endl;
  PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
  jobs.initNxtask();
//  long nexttask = jobs.nxtaskAlt( g_size ); // Keep J's agglomerated and dynamically balance I for now
long nexttask = jobs.nxtask( g_size );
  int ilo = 0;
  int ihi = l_maxspatials;
#endif
  
#ifdef PREALLOCATE
  if ( g_rank == 0 ) cout << " USING PREALLOCATED FETCH METHODS " << endl;
#endif
  
  JOUTFG fetchi;
  JOUTFG fetchj;
#ifdef PREALLOCATE
  l_deters->preallocateTempJoutfgObject( fetchi );
  l_deters->preallocateTempJoutfgObject( fetchj );
#endif


#ifdef DUMPHAMILTONIAN
  int num=GA::nodeid();
  stringstream sout;
  string hamfile="hamiltonian";
  string scount;
  sout << num;
  scount = sout.str();
  scount.insert(0,"-");
  hamfile.append( scount );

  ofstream ohamfile;
  ohamfile.open( hamfile.c_str() );
#endif

  pair<int,double> info1;
  
  int numElems=0;
  long totalElems=0;
  int block_status;
  
  int totalFetches=0; // Used for Tuning studies
  int orbmapfetch=0;
  
  // Build Hamiltonian from scratch
  // Just in case afterall we do an accumulate 

  flushH(); // this simply zeros out cimat,icicol, num and the supp spaces
  
  int countme=0;
  int countrun=0;
  
  vector<int>::iterator jit;
  
#ifdef NEWORBMAP
  if ( GA::nodeid() ==0 ) cout << "USING (2WAY-based) NEW GOP based ORBMAP approach " << endl;
  vector<vector<pair<short,short> > > ijpairlist( l_maxspatials );
  int listsize = l_deters->fetchAllCompressedOrbMap( ijpairlist ); // Collective: performs a BRDCST underneath
   if ( GA::nodeid() ==0 ) cout << "replicated orbmap size is " << listsize << endl;
#endif

  for(int iconf = ilo; iconf < ihi; ++iconf ) { 
    
#ifndef STATIC
    if ( nexttask == iconf ) { // in the loop stride is implemented 
#endif

      vector<vector<double> > cimat;
      vector<vector<int> > icicol;
      vector<int> number;
      
      temp2 = psociTime();
#ifdef PREALLOCATE
      l_deters->fetchAndUnpackDeterminantDataPreallocate( iconf+1, buffer, fetchi );
#else
      l_deters->fetchAndUnpackDeterminantData( iconf+1, fetchi );
#endif
      fetchTime += psociTime() - temp2;
      ++totalFetches;
      
      int numElems= 0;

      vector<int> list;
	  int nsefi = fetchi.nsefi;
	  if ( cimat.size() != nsefi ) {
	    cimat.resize(nsefi );
	    icicol.resize(nsefi );
	    number.resize(nsefi );
          }

	  // For each I compute orbmap list of J

          temp3 = psociTime();
#ifndef NEWORBMAP
	  vector<COMPRESS_SIZE> imap( l_nbf+2  );
	  l_deters->fetchOrbMap( iconf+1, imap); 
          ++orbmapfetch;
	  countrun +=  l_deters->compareAllOrbMap(iconf+1, imap, list); // fetch jlist for each i usually 100-300 entries
#else
          countrun +=  l_deters->compareAllOrbMapRemoveZeros(iconf+1, ijpairlist, list );
#endif
          mapTime += psociTime() - temp3;
	  ++allcount;
      
// Convert this to loops over the J lists returned from the orbmaps

	  for( jit=list.begin(); jit!=list.end();++jit) { // now grab actual list

            int jconf = (*jit);

            temp2 = psociTime();
#ifdef PREALLOCATE
	    l_deters->fetchAndUnpackDeterminantDataPreallocate( jconf+1, buffer, fetchj );
#else
	    l_deters->fetchAndUnpackDeterminantData( jconf+1, fetchj );
#endif
            fetchTime += psociTime() - temp2;
	    ++totalFetches;

	    temp2 = psociTime();
	    block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj, sefsef );
	    sociTime += psociTime() - temp2;
	    
            int numTest=0;
	    if ( block_status ) { // Still could be an empty interaction
	      numTest = packSefsef( fetchi, fetchj, sefsef, cimat, icicol, number );
	      numElems += numTest;
	    }
#ifdef DUMPHAMILTONIAN
  if ( numTest > 0 ) {
    ohamfile << iconf<< " " << jconf<< " " << numTest << endl;
  }
#endif

	  }
	  totalElems += numElems;
          if ( numElems > 0 ) {
#ifdef SUPPLEMENTAL
	  uploadCIMATtoGASupplemental( l_nsef[ fetchi.index - 1 ], cimat, icicol, number);
#else
	  uploadCIMATtoGA( l_nsef[ fetchi.index - 1 ], cimat, icicol, number);
#endif
        }

#ifndef STATIC
      nexttask = jobs.nxtask( g_size );
    }
#endif
  } // iconf
#ifdef SQUAREHAMILTONIAN
  cout << GA::nodeid() << " square hamiltonian times " << endl;
#endif
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "fetchTime is " << fetchTime<<endl;
  cout << GA::nodeid() << "mapTime is " << mapTime<<endl;
  cout << "allcount is " << allcount<<endl; //this is correct
  cout << "countrun is " << countrun << endl;
  cout << "total Fetches is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
  cout << "total orbmap fetches is " << orbmapfetch << endl;
  //cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif
  
#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif
  
#ifdef DUMPHAMILTONIAN
    ohamfile.close();
#endif

  GA::SERVICES.sync(); // important 
  return( numberNonzeroElements( totalElems ) ); // DGOP underneath
}


//Timer wrapped version
long PsociGAhamiltonian::constructGlobalTransitionHamiltonian2WAY( pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalTransitionHamiltonian2WAY();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( status );
}

/* A 2D H update scheme. Several more layers of indexing are required BUT at the benefit
   of many fewer calls to fetch
   
   BEST compiler options are : -DTWOWAY -DPREALLOCATE -DFETCH2D -DSUPPLEMENTAL - more local memory 
                2nd BEST are : -DTWOWAY -DFETCH2D -DSUPPLEMENTAL          
                     min are : -DTWOWAY -DFETCH2D - looses SPLIT-GA storage
*/

long PsociGAhamiltonian::constructGlobalTransitionHamiltonian2WAY()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  if ( g_rank == 0 ) cout << "I am " << GA::nodeid() << " at the fetch2d 2WAY-BLOCKED SPLIT TRANSITION hamiltonian " << endl;
  
  int counter=0;
  long allcount=0;
  
  double temp1, temp2;
  double jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;
  
#ifdef PREALLOCATE
  int l_maxlength = l_deters->fetchMaxLength();
  vector<COMPRESS_SIZE> buffer( max(l_granularity,l_nxtval_chunk) * l_maxlength );
#endif
  
#ifdef SUPPLEMENTAL
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif

#ifdef STATIC
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposition is contiguous and inclusive. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used " << endl;
  PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
  jobs.initNxtask();
  long nexttask = jobs.nxtaskAlt( g_size ); // Keep J's agglomerated and dynamically balance I for now
//  long inexttask=-1;
  int ilo = 0;
  int ihi = l_maxspatials;
#endif
  
#ifdef PREALLOCATE
  if ( g_rank == 0 ) cout << " USING PREALLOCATED FETCH METHODS " << endl;
#endif
  
#ifdef PREALLOCATE
  JOUTFG temp;
  l_deters->preallocateTempJoutfgObject( temp );
  l_deters2->preallocateTempJoutfgObject( temp );
  vector<JOUTFG> fetchi( l_nxtval_chunk, temp );
  vector<JOUTFG> fetchj( l_granularity, temp );
#else
  vector<JOUTFG> fetchi;
  vector<JOUTFG> fetchj;
#endif

  pair<int,double> info1;
  
  int numElems=0;
  long totalElems=0;
  int block_status;
  
  int totalFetches=0; // Used for Tuning studies
  
  // Build Hamiltonian from scratch
  
/* Increment the loop by l_nxtval_chunk within the loop.  This way
   we allow nxttask to dynamically determine a staritng iconf for a 
   contiguous compute set
*/
  
// Just in case afterall we do an accumulate 

  flushH(); // this simply zeros out cimat,icicol, num and the supplimental spaces

// cout << "ilo and ihi values " << GA::nodeid() <<" " <<ilo<<" "<<ihi<<endl;

  int countme=0;

  for(int iconf = ilo; iconf < ihi; ++iconf ) { 
    
#ifndef STATIC
    if ( nexttask == iconf ) { // in the loop stride is implemented 
      
#endif
      int ifinal = min( l_nxtval_chunk, ihi - iconf ); // ihi already adds a ONE
      vector<vector<vector<double> > > cimat(ifinal);
      vector<vector<vector<int> > > icicol( ifinal);
      vector<vector<int> > number(ifinal);
      
      temp2 = psociTime();
#ifdef PREALLOCATE
      l_deters->fetchAndUnpackDeterminantData2DPreallocate( iconf+1, iconf+ifinal, buffer, fetchi );
#else
      l_deters->fetchAndUnpackDeterminantData2D( iconf+1, iconf+ifinal, fetchi );
#endif
      
      fetchTime += psociTime() - temp2;
      ++totalFetches;
      
      vector<int> iindex( ifinal); // populate absolute I index array - added code but of fixed work;
      for(int ii=0; ii< ifinal; ++ii ) {
	iindex[ii]=iconf+ii;
      }
      
      vector<int> numElems(ifinal, 0); 
      
      for(int iconf_block=0; iconf_block < ifinal; ++iconf_block )
	{
	  int nsefi = fetchi[iconf_block].nsefi;
	  if ( cimat[iconf_block].size() != nsefi ) {
	    cimat[iconf_block].resize(nsefi );
	    icicol[iconf_block].resize(nsefi );
	    number[iconf_block].resize(nsefi );
	  }
	}
      
      for(int jconf=0; jconf < iconf+ifinal ; jconf += l_granularity )
	{
	  int kfinal = min( l_granularity, (iconf+ifinal-1) - jconf + 1  );
	  temp2 = psociTime();
#ifdef PREALLOCATE
	  l_deters2->fetchAndUnpackDeterminantData2DPreallocate( jconf+1, jconf + kfinal, buffer, fetchj );
#else
	  fetchj.clear();
	  l_deters2->fetchAndUnpackDeterminantData2D( jconf+1, jconf + kfinal, fetchj );
#endif
	  
	  fetchTime += psociTime() - temp2;
	  ++totalFetches; // comment me out eventually
	  
	  vector<int> jindex( kfinal); // populate absolute J index array
	  for(int jj=0; jj< kfinal; ++jj ) {
	    jindex[jj]=jconf+jj;
	  }
	  
	  //SOCI always clears sefsef
	  
	  for(int iconf_block=0; iconf_block < ifinal; ++iconf_block ) //always ONE for now
	    {
	      for(int jconf_block=0; jconf_block < kfinal; ++jconf_block )
		{
		  if ( jindex[jconf_block] <= iindex[iconf_block] ) {
		    ++allcount;
		    temp2 = psociTime();
		    block_status = l_gaHamiltonian->generateSOCIblocks( fetchi[iconf_block], fetchj[jconf_block], sefsef ); 
		    sociTime += psociTime() - temp2;
		    
		    if ( block_status ) {
		      int numTest = packSefsef( fetchi[iconf_block], fetchj[jconf_block] , sefsef, cimat[iconf_block], icicol[iconf_block], number[iconf_block] );
		      numElems[iconf_block] += numTest;
		    }
		  } 
		}
	    }
	} // jconf
      
      /* Update CIMAT. For a given I we need ALL contributing J values */
      //      int iconf_test=iconf;
      for(int iconf_block=0; iconf_block < ifinal; ++iconf_block )
	{
	  if ( numElems[iconf_block] > 0 ) { // For this to not be true would be somewhat rare
	    //            ++iconf_test;
	    //            cout << "iconf " << iconf_test<<" "<<numElems[iconf_block];
	    
	    totalElems += numElems[iconf_block];
#ifdef SUPPLEMENTAL
	    uploadCIMATtoGASupplemental( l_nsef[ fetchi[iconf_block].index - 1 ], cimat[iconf_block], icicol[iconf_block], number[iconf_block]);
#else
	    uploadCIMATtoGA( l_nsef[ fetchi[iconf_block].index - 1 ], cimat[iconf_block], icicol[iconf_block], number[iconf_block]);
#endif
	  }
	} 
#ifndef STATIC
      nexttask = jobs.nxtaskAlt( g_size );
      iconf = iconf + (l_nxtval_chunk - 1); // skip to next possible value from jobs.nxtask
    }
#endif
    
  } // iconf
  
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "fetchTime is " << fetchTime<<endl;
  cout << "allcount is " << allcount<<endl; //this is correct
  cout << "total Fetches is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif
  
#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif

  GA::SERVICES.sync(); 

  return( numberNonzeroElements( totalElems ) );
}

//Timer wrapped version
/*
long PsociGAhamiltonian::constructGlobalHamiltonianExpAlternative( pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalHamiltonianExpAlternative();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( status );
}
*/

//Direct are not fully tested
long PsociGAhamiltonian::constructGlobalHamiltonianAltOrbMapDirect(pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalHamiltonianAltOrbMapDirect();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( status );
}

/* A 2D H update scheme. Several more layers of indexing are required BUT at the benefit
   of many fewer calls to fetch
   
   In this modification,  DIRECT refers to on-the-fly computation of the determinants from 
   spatial configurations.  On a good networking machine, the non-fetch2D direct versus no-direct
   were about the same; they were bound by new operations. On a slow network (Ranger?) this may
   have a big effect
 
   Lastly, we retain the original GAdeterminants call  1) to compute GLobals and 2) for subsequent
   analyses which use GA.

   BEST compiler options are :-DDIRECT -DSUPPLEMENTAL -DSINGLENODECONFIGS - uses less GA local memory 
                2nd BEST are :-DDIRECT -DSINGLENODECONFIGS

   Significantly modified  for using ORBMAPS
*/

long PsociGAhamiltonian::constructGlobalHamiltonianAltOrbMapDirect()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  if ( g_rank == 0 ) cout << "I am " << GA::nodeid() << " at the fetch2d DIRECT SPLIT ORBMAP hamiltonian " << endl;

  int counter=0;
  long allcount=0;

  double temp1, temp2;
  double jTime=0.0, orbmapTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;

#ifdef SUPPLEMENTAL
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif

#ifdef SQUAREHAMILTONIAN
  if ( g_rank == 0 ) cout << "I am " << GA::nodeid() << " SQUARE DIRECT HAMILTONIAN " << endl;
#endif

#ifdef STATIC
  GA::error(" Cannot do DIRECT and STATIC since g_det is not created",1);
#endif
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used " << endl;
  PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
  jobs.initNxtask();
  long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
//  int ilo = 0;
//  int ihi = l_maxspatials;
  int ilo = 0;
  int ihi = l_maxspatials;


/*
#ifdef PREALLOCATE
  if ( g_rank == 0 ) cout << " USING PREALLOCATED FETCH METHODS " << endl;
#endif
*/

  JOUTFG fetchi; 
  JOUTFG fetchj;

  pair<int,double> info1;

  long totalElems=0;
  int block_status;


  double compareT=0.0, fetchT=0.0;

  int totalFetches=0; // Used for Tuning studies

  // Build Hamiltonian from scratch
  // Just in case afterall we do an accumulate 

  flushH(); // this simply zeros out cimat,icicol, num and the supp spaces

  int countme=0;
  int countrun=0;

#ifdef NEWORBMAP
  if ( GA::nodeid() ==0 ) cout << "USING (2WAY-based) NEW GOP based ORBMAP approach: " << endl;
  vector<vector<pair<short,short> > > ijpairlist( l_maxspatials );
  int listsize = l_deters->fetchAllCompressedOrbMap( ijpairlist ); // Collective: performs a BRDCST underneath
   if ( GA::nodeid() ==0 ) cout << "replicated orbmap size is " << listsize << endl;
#endif

  vector<int> list; // list of identically non-zero js per i

  vector<COMPRESS_SIZE> imap( l_nbf+2  ); // a spatial configuration
  vector<COMPRESS_SIZE> jmap( l_nbf+2  );

  for(int iconf = ilo; iconf < ihi; ++iconf ) {

    int numElems = 0;

    if ( nexttask == iconf ) { // in the loop stride is implemented 
      vector<vector<double> > cimat(1);
      vector<vector<int> > icicol(1);
      vector<int> number(1);

      temp2 = psociTime();
      list.clear();
#ifndef NEWORBMAP
          double temp3 = psociTime();
          l_deters->fetchOrbMap( iconf+1, imap);
          fetchT += psociTime() - temp3;
          temp3 = psociTime();
          countrun +=  l_deters->compareAllOrbMap(iconf+1, imap, list); // fetch jlist- usually 100-300 entries
          compareT += psociTime() - temp3;
#else
          countrun +=  l_deters->compareAllOrbMapRemoveZeros(iconf+1, ijpairlist, list );
#endif
     orbmapTime += psociTime() - temp2;

     temp2 = psociTime();
     l_deters->computeConfigsGA( iconf+1, imap, fetchi ); // also does a fetchOrb underneath
      fetchTime += psociTime() - temp2;
      ++totalFetches;

      int nsefi = fetchi.nsefi;
          if ( cimat.size() != nsefi ) {
            cimat.resize( nsefi );
            icicol.resize( nsefi );
            number.resize( nsefi );
          }

      vector<int>::iterator jit;
      for( jit=list.begin(); jit!=list.end();++jit) { 
          int jconf = (*jit);

          temp2 = psociTime();
          l_deters->computeConfigsGA( jconf+1, jmap, fetchj );
          fetchTime += psociTime() - temp2;
          ++totalFetches; // comment me out eventually

//SOCI always clears sefsef

             ++allcount;
             temp2 = psociTime();

              block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj, sefsef );
              sociTime += psociTime() - temp2;

              if ( block_status ) {
              int numTest = packSefsef( fetchi, fetchj , sefsef, cimat, icicol, number );
              numElems += numTest;
               }
        } // jit

      /* Update CIMAT. For a given I we need ALL contributing J values */
          if ( numElems > 0 ) { // For this to not be true would be somewhat rare

            totalElems += numElems;
#ifdef SUPPLEMENTAL
            uploadCIMATtoGASupplemental( l_nsef[ fetchi.index - 1 ], cimat, icicol, number);
#else
            uploadCIMATtoGA( l_nsef[ fetchi.index - 1 ], cimat, icicol, number);
#endif
          }
      nexttask = jobs.nxtask( g_size );
    }

cout << endl; // EXTRA SPACE FOR GNUPLIT
  } // iconf

#ifdef SQUAREHAMILTONIAN
  cout << GA::nodeid() << " square hamiltonian times " << endl;
#else
  cout << GA::nodeid() << " trianglular hamiltonian times " << endl;
#endif
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "computeConfigsTime is " << fetchTime<<endl;
  cout << "allcount is " << allcount<<endl; //this is correct
  cout << "total computes is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
  cout << "orbmap time is " << orbmapTime << endl;
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  cout << "fetcT is " << fetchT << " compreT is " << compareT << endl;

//  tdeter->tearDownDeterminants();// no longer need the replicated CONFIGS/VEC_CONFIGS data

#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif

  jobs.destroyNxtask(); //Close'er up

  return( numberNonzeroElements( totalElems ) );
}

/* drive non orbmap code over an outer loop: specialize the method to this one case
*/



long PsociGAhamiltonian::constructGlobalHamiltonianAltOrbMapDirectChunkJmap(pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalHamiltonianAltOrbMapDirectChunkJmap();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( status );
}


long PsociGAhamiltonian::constructGlobalHamiltonianAltOrbMapChunkJmap(pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalHamiltonianAltOrbMapChunkJmap();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( status );
}


/* Highly experimental method: Requires alternative (accumulation based) upload methods */
/* Current default approach 
   Now includes jphase information */
// 

long PsociGAhamiltonian::constructGlobalHamiltonianAltOrbMapDirectChunkJmap()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  if ( g_rank == 0 ) cout << "I am " << GA::nodeid() << " EXPERIMENTAL at the JMAP CHUNK DIRECT SPLIT ORBMAP hamiltonian " << endl;
  
  int counter=0;
  long allcount=0;
  
  double temp1, temp2;
  double jTime=0.0, orbmapTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;

#ifdef SUPPLEMENTAL
 // GA::error(" AltOrbMapDirectChunkJmap must not use -DSUPPLEMENTAL",-1);
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif
#ifdef SQUAREHAMILTONIAN
  if ( g_rank == 0 ) cout << "I am " << GA::nodeid() << " SQUARE DIRECT HAMILTONIAN " << endl;
#endif
#ifdef NEWORBMAP
  GA::error(" cannot use constructGlobalHamiltonianAltOrbMapDirectChunkJmap method with NEWORBMAP",-1);
#endif

#ifdef STATIC
  if ( g_rank == 0 ) cout << "STATIC based load balancing scheme" << endl;  
#endif

#ifdef REVERSE
  if ( g_rank == 0 ) cout << " Perform looping in reverse " << endl;
#endif


  //if ( GA::nodeid() ==0 ) cout << "FETCH jmaps as a subset " << endl;
  
  const int numjchunks=l_numJchunksHbuild; // carve up jmap space into 10 chunks
  
  vector<int> list; // list of identically non-zero js per i
  
  int chunk_width;
  int l_maxspatials = l_deters->fetchMaxSpatials();
  
  (l_maxspatials % numjchunks==0)? chunk_width = l_maxspatials / numjchunks: chunk_width = 1 + (l_maxspatials / numjchunks);
  if ( GA::nodeid() == 0 ) cout << "Jchunk stats num: "<<numjchunks<<" width " << chunk_width << endl;
  int jchunk_index=0; // these track the actual real index numbers (0,maxefs-1)
  
  double jfetchTime = 0.0; 
  
  pair<int,double> info1;
  
  long totalElems=0;
  int block_status;
  
  double compareT=0.0, fetchT=0.0;
  double assemTime=0.0, compareTime=0.0;
 
  double temp3, totalCompute=0.0;

  int didblock=0;
  long totalblock=0;
  int excitblock=0;
  int phaseblock=0;
  int possibleblock=0;

  int totalFetches=0; // Used for Tuning studies
  
#ifndef STATIC
  if ( g_rank == 0 ) cout << " JSUBSET CHUNK: NXTVAL based load balancing scheme being used " << endl;
#endif

#ifdef REPLICATEDETVALUE
  if ( g_rank == 0 ) cout << "USING NEW JPHASE code for H build"<<endl; ;
// doesn't do anything
   vector<char> * rep_phases = l_deters->fetchReplicatedPhasesHandle();
 //  cout << "H build rep phase is size " << rep_phases->size()<<endl; // returns iopen
#endif
  
  // Build Hamiltonian from scratch
  // Just in case afterall we do an accumulate 

  flushH(); // this simply zeros out cimat,icicol, num and the supp spaces
  
//#ifndef REPLICATEDETVALUE
    //if ( g_rank == 0 ) cout << "Outer loop over prefetched J chunks " << endl;

    for(int jchunk=0; jchunk < numjchunks; ++jchunk ) { // data parallel

    if ( g_rank == 0 ) cout << " J loop = " << jchunk <<" "<< numjchunks << endl;

    double tempj = psociTime();
    vector<COMPRESS_SIZE> jsubsetlist; // jmap subset buffer gets sized in fetchJsubsetlist 
    vector<int> jindexlist;

//collectives could be used here

    int jnumlist = fetchJsubsetlist( jchunk, jindexlist, jsubsetlist, chunk_width, l_maxspatials ); //GA on the inside. resize width inside to be perfect
    jfetchTime += psociTime() - tempj;

    //GA::SERVICES.sync();
    //if ( g_rank == 0 )  cout << "jnumlist is " << GA::nodeid() <<" " << jnumlist << " " << chunk_width << " subset " << jsubsetlist.size() <<  endl;

//#endif
    
    // start old code and NXTVAL driver
#ifndef STATIC
    PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
#ifdef REVERSE
    jobs.initNxtask( l_maxspatials );
#else
    jobs.initNxtask();
#endif


#ifdef REVERSE
    long nexttask = jobs.nxtaskRev( g_size, l_maxspatials );
#else
    long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
#endif

#endif

    int ilo = 0;
    int ihi = l_maxspatials;
    
    JOUTFG fetchi; 
    JOUTFG fetchj;
    
    // Build Hamiltonian from scratch
    // Just in case afterall we do an accumulate 
    // flushH(); // this simply zeros out cimat,icicol, num and the supp spaces

    int countme=0;
    int countrun=0;
    
    vector<COMPRESS_SIZE> imap( l_nbf+2  ); // a spatial configuration
    vector<COMPRESS_SIZE> jmap( l_nbf+2  );

#ifndef STATIC
#ifdef REVERSE
    for( int iconf = ihi-1; iconf >= ilo; --iconf ) {
#else
    for( int iconf = ilo; iconf < ihi; ++iconf ) {
#endif
#endif

#ifdef STATIC
#ifdef REVERSE
//untested
    for( int iconf = ihi-1-g_rank; iconf >=ilo; iconf -= g_size) {
#else
    for( int iconf = g_rank; iconf < ihi; iconf += g_size) {
#endif
#endif
    
    int numElems = 0;
    totalblock += (iconf+1);

#ifndef STATIC
    if ( nexttask == iconf ) { // in the loop stride is implemented 
#endif

      // do this above
      vector<vector<double> > cimat(1);
      vector<vector<int> > icicol(1);
      vector<int> number(1);
      
      temp2 = psociTime();
      list.clear();

      int newJs=1;

// a compute step already gets imap
      temp2 = psociTime();
      l_deters->computeConfigsGA( iconf+1, imap, fetchi ); // also does a fetchOrb underneath
      fetchTime += psociTime() - temp2;
      ++totalFetches;
      
// THis does lower the number of blocks
#ifdef REPLICATEDETVALUE
//       temp2 = psociTime();
//       int iphase = rep_phases->at(fetchi.index-1);
//       vector<int> jphases;
//       newJs = assemblePhasePossibles(iconf, iphase, jphases, rep_phases );
//       phaseblock += jphases.size();
//       assemTime += psociTime() - temp2;
#endif

//       if ( newJs > 0 ) {

/*
#ifdef REPLICATEDETVALUE
       int currentrun =  l_deters->compareAllOrbMapMatrixVectorProducts(iconf+1, imap, jphases, list );
       countrun += currentrun;
       jphases.clear();
       excitblock += list.size();
       //if ( GA::nodeid() == 0 ) cout << "curr excit lsit size is " << currentrun << " " << list.size() << " emp J " << jphases.size() << endl;
       orbmapTime += psociTime() - temp2;
#endif
*/

//#ifndef REPLICATEDETVALUE
      temp2 = psociTime();
      //int currentrun =  l_deters->compareLocalJSubsetOrbMap(jnumlist, iconf+1, imap, jindexlist, jsubsetlist, list); // 

      int currentrun =  l_deters->compareLocalJSubsetOrbMap(jnumlist, iconf+1, rep_phases, imap, jindexlist, jsubsetlist, list); // 
      countrun += currentrun;
      phaseblock=-1;
      excitblock += list.size();
      orbmapTime += psociTime() - temp2;
//#endif

      // new countrun could be zero so we need to account for it.
      
      if ( currentrun != 0 ) {

	int nsefi = fetchi.nsefi;
	if ( cimat.size() != nsefi ) {
	  cimat.resize( nsefi );
	  icicol.resize( nsefi );
	  number.resize( nsefi );
	}

        temp3=psociTime();

	vector<int>::iterator jit;
	for( jit=list.begin(); jit!=list.end();++jit) { 

//cout << "start list " << list.size() << " " << (*jit) << endl;

          ++didblock;
          int jconf = (*jit);
	  
          temp2 = psociTime();
          l_deters->computeConfigsGA( jconf+1, jmap, fetchj );
          fetchTime += psociTime() - temp2;
          ++totalFetches; // comment me out eventually
	  
	  //SOCI always clears sefsef
	  
	  ++allcount;
	  temp2 = psociTime();
	  
//cout << "SOCI block " << endl;
	  block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj, sefsef );
	  sociTime += psociTime() - temp2;
	  

	  if ( block_status ) {
//cout << "block_status " << endl;
	    int numTest = packSefsef( fetchi, fetchj , sefsef, cimat, icicol, number );
	    numElems += numTest;
	  }
        } // jit
        totalCompute += psociTime() - temp3;
	
	/* Update CIMAT. For a given I we need ALL contributing J values */
	if ( numElems > 0 ) { // For this to not be true would be somewhat rare
	  totalElems += numElems;
	  
	  /* since we repeat loop over iconf we must be able to accumulate into H and not just put
	   */
#ifdef SUPPLEMENTAL
	  temp2 = psociTime();
	  uploadCIMATtoGASupplementalAcc( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
	  uploadTime += psociTime() - temp2;
#else
	  uploadCIMATtoGAAcc( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
#endif
	  //cant use this            uploadCIMATtoGAAcc( l_nsef[ fetchi.index - 1 ], cimat, icicol, number);
	}
      } // currentrun != 0
     //} // newjs > 0

#ifndef STATIC
#ifdef REVERSE
      nexttask = jobs.nxtaskRev( g_size, l_maxspatials );
#else
      nexttask = jobs.nxtask( g_size );
#endif
    }
#endif
    
  } // iconf

#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif

  } // jchunk
    
// accumulate blocks

  if ( GA::nodeid() == 0 ) {
#ifdef SQUAREHAMILTONIAN
  cout << GA::nodeid() << " square hamiltonian times " << endl;
#else 
  cout << GA::nodeid() << " trianglular hamiltonian times " << endl;
#endif
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "computeConfigsTime is " << fetchTime<<endl;
  cout << GA::nodeid() << "jfetchTime is " << jfetchTime<<endl;
  cout << "allcount is " << allcount<<endl; //this is correct
  cout << "total computes is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
  cout << "orbmap time is " << orbmapTime << endl;
  cout << "totalCompute exclude upload " << totalCompute << endl;
 cout << "upload time is " << uploadTime << endl;
  cout << "summed didblock phaseblock excitblock totalblock " <<didblock<<" "<< phaseblock<<" "<<excitblock<<" "<<totalblock<<endl;
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  //cout << "fetchT is " << fetchT << " compreT is " << compareT << endl;
  }
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif   

  char nop[] = "+";
  long gone=1;
//  GA::igop( &phaseblock, gone, nop);
//  if ( GA::nodeid() == 0 ) cout << "SUMMED phaseblock " << phaseblock<< endl;
/*
  GA::igop( &ldidblock, gone, nop);
  if ( GA::nodeid() == 0 ) cout << "SUMMED didblock " <<didblock<< endl;
  GA::igop( &ltotalFetches, gone, nop);
  if ( GA::nodeid() == 0 ) cout << "SUMMED fetches " << totalFetches<< endl;
*/

  return( numberNonzeroElements( totalElems ) );
}

// New nondirect H version
long PsociGAhamiltonian::constructGlobalHamiltonianAltOrbMapChunkJmap()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  if ( g_rank == 0 ) cout << "I am " << GA::nodeid() << "EXPERIMENTAL at the JMAP CHUNK NONDIRECT SPLIT ORBMAP hamiltonian " << endl;
  
  int counter=0;
  long allcount=0;
  
  double temp1, temp2;
  double jTime=0.0, orbmapTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;

#ifdef SUPPLEMENTAL
 // GA::error(" AltOrbMapDirectChunkJmap must not use -DSUPPLEMENTAL",-1);
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif
#ifdef NEWORBMAP
  GA::error(" cannot use constructGlobalHamiltonianAltOrbMapChunkJmap method with NEWORBMAP",-1);
#endif
#ifdef STATIC
  if ( g_rank == 0 ) cout << "STATIC non direct"<< endl;
#endif
#ifdef REVERSE
  if ( g_rank == 0 ) cout << "NONDIRECT: Using reverse iconf indexing " << endl;
#endif

  //if ( GA::nodeid() ==0 ) cout << "FETCH jmaps as a subset " << endl;
  
  const int numjchunks=l_numJchunksHbuild; // carve up jmap space into 10 chunks
  
  vector<int> list; // list of identically non-zero js per i
  
  int chunk_width;
  int l_maxspatials = l_deters->fetchMaxSpatials();
  
  (l_maxspatials % numjchunks==0)? chunk_width = l_maxspatials / numjchunks: chunk_width = 1 + (l_maxspatials / numjchunks);
  if ( GA::nodeid() ==0 ) cout << "Jchunk stats num: "<<numjchunks<<" width " << chunk_width << endl;
  int jchunk_index=0; // these track the actual real index numbers (0,maxefs-1)
  
  double jfetchTime = 0.0; 
  
  pair<int,double> info1;
  
  long totalElems=0;
  int block_status;
  
  double compareT=0.0, fetchT=0.0;
  double assemTime=0.0, compareTime=0.0;
 
  double temp3, totalCompute=0.0;

  int didblock=0;
  long totalblock=0;
  int excitblock=0;
  int phaseblock=0;
  int possibleblock=0;

  int totalFetches=0; // Used for Tuning studies
  
  if ( g_rank == 0 ) cout << " JSUBSET CHUNK: NXTVAL based load balancing scheme being used " << endl;

#ifdef REPLICATEDETVALUE
  if ( g_rank == 0 ) cout << "USING NEW JPHASE code for H build"<<endl; ;
 // doesn't do anything
   vector<char> * rep_phases = l_deters->fetchReplicatedPhasesHandle();
 //  cout << "H build rep phase is size " << rep_phases->size()<<endl; // returns iopen
#endif
  
  // Build Hamiltonian from scratch
  // Just in case afterall we do an accumulate 

  flushH(); // this simply zeros out cimat,icicol, num and the supp spaces
  
//#ifndef REPLICATEDETVALUE
    if ( g_rank == 0 ) cout << "Outer loop over prefetched J chunks " << endl;

    for(int jchunk=0; jchunk < numjchunks; ++jchunk ) { // data parallel
    double tempj = psociTime();
    vector<COMPRESS_SIZE> jsubsetlist; // jmap subset buffer gets sized in fetchJsubsetlist 
    vector<int> jindexlist;
    int jnumlist = fetchJsubsetlist( jchunk, jindexlist, jsubsetlist, chunk_width, l_maxspatials ); //GA on the inside. resize width inside to be perfect
    jfetchTime += psociTime() - tempj;
//#endif
    
// Can we set this up to loop in reverse ? 

    // start old code and NXTVAL driver
    // maybe need ot restart here

#ifndef STATIC
    PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );

#ifdef REVERSE
    jobs.initNxtask( l_maxspatials );
#else
    jobs.initNxtask();
#endif

#ifdef REVERSE
    long nexttask = jobs.nxtaskRev( g_size, l_maxspatials );
#else
    long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
#endif
#endif

    int ilo = 0;
    int ihi = l_maxspatials;
//#endif
    
    JOUTFG fetchi; 
    JOUTFG fetchj;
    
    // Build Hamiltonian from scratch
    // Just in case afterall we do an accumulate 
    // flushH(); // this simply zeros out cimat,icicol, num and the supp spaces

    int countme=0;
    int countrun=0;
    
    vector<COMPRESS_SIZE> imap( l_nbf+2  ); // a spatial configuration
    vector<COMPRESS_SIZE> jmap( l_nbf+2  );
  
// Can we loop in reverse ?

#ifndef STATIC
#ifdef REVERSE
  for( int iconf = ihi-1; iconf >= ilo; --iconf ) {
#else
  for( int iconf = ilo; iconf < ihi; ++iconf ) {
#endif
#endif

#ifdef STATIC
#ifdef REVERSE
  for( int iconf = ihi-1-g_rank; iconf >= ilo; iconf-=g_size ) {
#else
  for( int iconf = ilo+g_rank; iconf < ihi; iconf+=g_size ) {
#endif
#endif

    int numElems = 0;
    totalblock += (iconf+1);

    //cout << GA::nodeid() <<" iconf " << iconf << " my nexttask is " << nexttask << endl;
#ifndef STATIC
    if ( nexttask == iconf ) { // in the loop stride is implemented 
#endif

     //cout << GA::nodeid() <<" " << ilo<<" "<<ihi<< " iconf : doing some work " << iconf<<" "<<nexttask << endl;

      // do this above
      vector<vector<double> > cimat(1);
      vector<vector<int> > icicol(1);
      vector<int> number(1);
      
      temp2 = psociTime();
      list.clear();

      int newJs=1;

// a compute step already gets imap

      temp2 = psociTime();
      l_deters->fetchAndUnpackDeterminantData( iconf+1, info1, fetchi );
      fetchTime += psociTime() - temp2;
      ++totalFetches;
      
      temp2 = psociTime();
      //int currentrun =  l_deters->compareLocalJSubsetOrbMap(jnumlist, iconf+1, imap, jindexlist, jsubsetlist, list); // 

      l_deters->fetchOrbMap( iconf+1, imap);
      int currentrun =  l_deters->compareLocalJSubsetOrbMap(jnumlist, iconf+1, rep_phases, imap, jindexlist, jsubsetlist, list); // 
      countrun += currentrun;
      phaseblock=-1;
      excitblock += list.size();
      orbmapTime += psociTime() - temp2;

      // new countrun could be zero so we need to account for it.
      
      if ( currentrun != 0 ) {

	int nsefi = fetchi.nsefi;
	if ( cimat.size() != nsefi ) {
	  cimat.resize( nsefi );
	  icicol.resize( nsefi );
	  number.resize( nsefi );
	}

        temp3=psociTime();

	vector<int>::iterator jit;
       // cout << "jit list size is " << list.size() << endl;

	for( jit=list.begin(); jit!=list.end();++jit) { 
          ++didblock;
          int jconf = (*jit);
	  
          temp2 = psociTime();
          l_deters->fetchAndUnpackDeterminantData( jconf+1, info1, fetchj );
          fetchTime += psociTime() - temp2;
          ++totalFetches; // comment me out eventually
	  
	  //SOCI always clears sefsef
	  
	  ++allcount;
	  temp2 = psociTime();
	  
	  block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj, sefsef );
	  sociTime += psociTime() - temp2;

	  if ( block_status ) {
	    int numTest = packSefsef( fetchi, fetchj , sefsef, cimat, icicol, number );
	    numElems += numTest;
	  }
        } // jit

        totalCompute += psociTime() - temp3;
	
	/* Update CIMAT. For a given I we need ALL contributing J values */
	if ( numElems > 0 ) { // For this to not be true would be unlikely (but possible with disjoint Jblocks) 
	  totalElems += numElems;
	  
	  /* since we repeat loop over iconf we must be able to accumulate into H and not just put
	   */
#ifdef SUPPLEMENTAL
	  temp2 = psociTime();
	  uploadCIMATtoGASupplementalAcc( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
	  uploadTime += psociTime() - temp2;
#else
	  uploadCIMATtoGAAcc( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
#endif
	  //cant use this            uploadCIMATtoGAAcc( l_nsef[ fetchi.index - 1 ], cimat, icicol, number);
	}
      } // currentrun != 0
     //} // newjs > 0

#ifndef STATIC
#ifdef REVERSE
      nexttask = jobs.nxtaskRev( g_size, l_maxspatials );
#else
      nexttask = jobs.nxtask( g_size );
#endif
    }
#endif
  } // iconf


#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif

//ifndef REPLICATEDETVALUE
  } // jchunk
//#endif
    
// accumulate blocks

  if ( GA::nodeid() == 0 ) {
  cout << GA::nodeid() << " trianglular hamiltonian times " << endl;
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "fetchConfigsTime is " << fetchTime<<endl;
  cout << GA::nodeid() << "jfetchTime is " << jfetchTime<<endl;
  cout << "allcount is " << allcount<<endl; //this is correct
  cout << "total computes is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
  cout << "orbmap time is " << orbmapTime << endl;
  cout << "totalCompute exclude upload " << totalCompute << endl;
 cout << "upload time is " << uploadTime << endl;
  cout << "summed didblock phaseblock excitblock totalblock " <<didblock<<" "<< phaseblock<<" "<<excitblock<<" "<<totalblock<<endl;
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  //cout << "fetchT is " << fetchT << " compreT is " << compareT << endl;
  }
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif   

  char nop[] = "+";
  long gone=1;
//  GA::igop( &phaseblock, gone, nop);
//  if ( GA::nodeid() == 0 ) cout << "SUMMED phaseblock " << phaseblock<< endl;
/*
  GA::igop( &ldidblock, gone, nop);
  if ( GA::nodeid() == 0 ) cout << "SUMMED didblock " <<didblock<< endl;
  GA::igop( &ltotalFetches, gone, nop);
  if ( GA::nodeid() == 0 ) cout << "SUMMED fetches " << totalFetches<< endl;
*/

  return( numberNonzeroElements( totalElems ) );
}
// Crude version to siomply report sparsitry values
void PsociGAhamiltonian::constructGlobalHamiltonianAltOrbMapChunkJmapEstimatesOnly()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  if ( g_rank == 0 ) cout << "I am " << GA::nodeid() << "EXPERIMENTAL at the JMAP CHUNK NONDIRECT SPLIT ORBMAP Estimatyes only: hamiltonian " << endl;
  
  int counter=0;
  long allcount=0;
  
  double temp1, temp2;
  double jTime=0.0, orbmapTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;

#ifdef SUPPLEMENTAL
 // GA::error(" AltOrbMapDirectChunkJmap must not use -DSUPPLEMENTAL",-1);
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif
#ifdef NEWORBMAP
  GA::error(" cannot use constructGlobalHamiltonianAltOrbMapChunkJmap method with NEWORBMAP",-1);
#endif
#ifdef STATIC
  if ( g_rank == 0 ) cout << "STATIC non direct"<< endl;
#endif
#ifdef REVERSE
  if ( g_rank == 0 ) cout << "NONDIRECT: Using reverse iconf indexing " << endl;
#endif

  //if ( GA::nodeid() ==0 ) cout << "FETCH jmaps as a subset " << endl;
  
  const int numjchunks=l_numJchunksHbuild; // carve up jmap space into 10 chunks
  
  vector<int> list; // list of identically non-zero js per i
  
  int chunk_width;
  int l_maxspatials = l_deters->fetchMaxSpatials();
  
  (l_maxspatials % numjchunks==0)? chunk_width = l_maxspatials / numjchunks: chunk_width = 1 + (l_maxspatials / numjchunks);
  if ( GA::nodeid() ==0 ) cout << "Jchunk stats num: "<<numjchunks<<" width " << chunk_width << endl;
  int jchunk_index=0; // these track the actual real index numbers (0,maxefs-1)
  
  double jfetchTime = 0.0; 
  
  pair<int,double> info1;
  
  long totalElems=0;
  int block_status;
  
  double compareT=0.0, fetchT=0.0;
  double assemTime=0.0, compareTime=0.0;
 
  double temp3, totalCompute=0.0;

  int didblock=0;
  long totalblock=0;
  int excitblock=0;
  int phaseblock=0;
  int possibleblock=0;

  int totalFetches=0; // Used for Tuning studies
  
  if ( g_rank == 0 ) cout << " JSUBSET CHUNK: NXTVAL based load balancing scheme being used " << endl;

#ifdef REPLICATEDETVALUE
  if ( g_rank == 0 ) cout << "USING NEW JPHASE code for H build"<<endl; ;
 // doesn't do anything
   vector<char> * rep_phases = l_deters->fetchReplicatedPhasesHandle();
 //  cout << "H build rep phase is size " << rep_phases->size()<<endl; // returns iopen
#endif
  
  // Build Hamiltonian from scratch
  // Just in case afterall we do an accumulate 

  flushH(); // this simply zeros out cimat,icicol, num and the supp spaces
  
//#ifndef REPLICATEDETVALUE
    if ( g_rank == 0 ) cout << "Outer loop over prefetched J chunks " << endl;

    for(int jchunk=0; jchunk < numjchunks; ++jchunk ) { // data parallel
    double tempj = psociTime();
    vector<COMPRESS_SIZE> jsubsetlist; // jmap subset buffer gets sized in fetchJsubsetlist 
    vector<int> jindexlist;
    int jnumlist = fetchJsubsetlist( jchunk, jindexlist, jsubsetlist, chunk_width, l_maxspatials ); //GA on the inside. resize width inside to be perfect
    jfetchTime += psociTime() - tempj;
//#endif
    
// Can we set this up to loop in reverse ? 

    // start old code and NXTVAL driver
    // maybe need ot restart here

#ifndef STATIC
    PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );

#ifdef REVERSE
    jobs.initNxtask( l_maxspatials );
#else
    jobs.initNxtask();
#endif

#ifdef REVERSE
    long nexttask = jobs.nxtaskRev( g_size, l_maxspatials );
#else
    long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
#endif
#endif

    int ilo = 0;
    int ihi = l_maxspatials;
//#endif
    
    JOUTFG fetchi; 
    JOUTFG fetchj;
    
    // Build Hamiltonian from scratch
    // Just in case afterall we do an accumulate 
    // flushH(); // this simply zeros out cimat,icicol, num and the supp spaces

    int countme=0;
    int countrun=0;
    
    vector<COMPRESS_SIZE> imap( l_nbf+2  ); // a spatial configuration
    vector<COMPRESS_SIZE> jmap( l_nbf+2  );
  
// Can we loop in reverse ?

#ifndef STATIC
#ifdef REVERSE
  for( int iconf = ihi-1; iconf >= ilo; --iconf ) {
#else
  for( int iconf = ilo; iconf < ihi; ++iconf ) {
#endif
#endif

#ifdef STATIC
#ifdef REVERSE
  for( int iconf = ihi-1-g_rank; iconf >= ilo; iconf-=g_size ) {
#else
  for( int iconf = ilo+g_rank; iconf < ihi; iconf+=g_size ) {
#endif
#endif

    int numElems = 0;
    totalblock += (iconf+1);

    //cout << GA::nodeid() <<" iconf " << iconf << " my nexttask is " << nexttask << endl;
#ifndef STATIC
    if ( nexttask == iconf ) { // in the loop stride is implemented 
#endif

     //cout << GA::nodeid() <<" " << ilo<<" "<<ihi<< " iconf : doing some work " << iconf<<" "<<nexttask << endl;

      // do this above
      vector<vector<double> > cimat(1);
      vector<vector<int> > icicol(1);
      vector<int> number(1);
      
      temp2 = psociTime();
      list.clear();

      int newJs=1;

// a compute step already gets imap

      temp2 = psociTime();
      l_deters->fetchAndUnpackDeterminantData( iconf+1, info1, fetchi );
      fetchTime += psociTime() - temp2;
      ++totalFetches;
      
      temp2 = psociTime();
      //int currentrun =  l_deters->compareLocalJSubsetOrbMap(jnumlist, iconf+1, imap, jindexlist, jsubsetlist, list); // 

      l_deters->fetchOrbMap( iconf+1, imap);
      int currentrun =  l_deters->compareLocalJSubsetOrbMap(jnumlist, iconf+1, rep_phases, imap, jindexlist, jsubsetlist, list); // 
      countrun += currentrun;
      phaseblock=-1;
      excitblock += list.size();
      orbmapTime += psociTime() - temp2;

      // new countrun could be zero so we need to account for it.
      
      if ( currentrun != 0 ) {

	int nsefi = fetchi.nsefi;
	if ( cimat.size() != nsefi ) {
	  cimat.resize( nsefi );
	  icicol.resize( nsefi );
	  number.resize( nsefi );
	}

        temp3=psociTime();

	vector<int>::iterator jit;
       // cout << "jit list size is " << list.size() << endl;

	for( jit=list.begin(); jit!=list.end();++jit) { 
          ++didblock;
          int jconf = (*jit);
	  
          temp2 = psociTime();
          l_deters->fetchAndUnpackDeterminantData( jconf+1, info1, fetchj );
          fetchTime += psociTime() - temp2;
          ++totalFetches; // comment me out eventually
	  
	  //SOCI always clears sefsef
	  
	  ++allcount;
	  temp2 = psociTime();
	  
	  block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj, sefsef );
	  sociTime += psociTime() - temp2;

	  if ( block_status ) {
	    int numTest = packSefsef( fetchi, fetchj , sefsef, cimat, icicol, number );
	    numElems += numTest;
	  }
        } // jit

        totalCompute += psociTime() - temp3;
	
	/* Update CIMAT. For a given I we need ALL contributing J values */
	if ( numElems > 0 ) { // For this to not be true would be unlikely (but possible with disjoint Jblocks
	  cout << "Iter total numE " << totalElems<<" "<<numElems<<endl;
	  totalElems += numElems;
	  int nsefi = number.size();
	    for(int ii=0; ii< nsefi; ++ii ) {
            cout << "Riow sparsity " << number[ii] << endl;
	    }
	  /* since we repeat loop over iconf we must be able to accumulate into H and not just put
	   */
	  //cant use this            uploadCIMATtoGAAcc( l_nsef[ fetchi.index - 1 ], cimat, icicol, number);
	}
      } // currentrun != 0
     //} // newjs > 0

#ifndef STATIC
#ifdef REVERSE
      nexttask = jobs.nxtaskRev( g_size, l_maxspatials );
#else
      nexttask = jobs.nxtask( g_size );
#endif
    }
#endif
  } // iconf


#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif

//ifndef REPLICATEDETVALUE
  } // jchunk
//#endif
    
// accumulate blocks

  if ( GA::nodeid() == 0 ) {
  cout << GA::nodeid() << " trianglular hamiltonian times " << endl;
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "fetchConfigsTime is " << fetchTime<<endl;
  cout << GA::nodeid() << "jfetchTime is " << jfetchTime<<endl;
  cout << "allcount is " << allcount<<endl; //this is correct
  cout << "total computes is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
  cout << "orbmap time is " << orbmapTime << endl;
  cout << "totalCompute exclude upload " << totalCompute << endl;
 cout << "upload time is " << uploadTime << endl;
  cout << "summed didblock phaseblock excitblock totalblock " <<didblock<<" "<< phaseblock<<" "<<excitblock<<" "<<totalblock<<endl;
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  //cout << "fetchT is " << fetchT << " compreT is " << compareT << endl;
  }
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif   

  char nop[] = "+";
  long gone=1;
//  GA::igop( &phaseblock, gone, nop);
//  if ( GA::nodeid() == 0 ) cout << "SUMMED phaseblock " << phaseblock<< endl;
/*
  GA::igop( &ldidblock, gone, nop);
  if ( GA::nodeid() == 0 ) cout << "SUMMED didblock " <<didblock<< endl;
  GA::igop( &ltotalFetches, gone, nop);
  if ( GA::nodeid() == 0 ) cout << "SUMMED fetches " << totalFetches<< endl;
*/
 // return( numberNonzeroElements( totalElems ) );
}
 
long PsociGAhamiltonian::constructGlobalHamiltonianAltOrbMapDirectChunkJmapInner(pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalHamiltonianAltOrbMapDirectChunkJmapInner();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( status );
}


/* newer method with lots of orbmap comms. In this case we inset the jchunk loop INSIDE the iconf loop.
   Thus we have lots calls to process Jsets which ( as expected, is slow 
   
   This method was intended as an experimental method not to be deployed for calculations
*/
long PsociGAhamiltonian::constructGlobalHamiltonianAltOrbMapDirectChunkJmapInner()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  if ( g_rank == 0 ) cout << "I am " << GA::nodeid() << " at the INNER JMAP CHUNK DIRECT SPLIT ORBMAP hamiltonian " << endl;
  
  int counter=0;
  long allcount=0;
  
  double temp1, temp2;
  double jTime=0.0, orbmapTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;

#ifdef SUPPLEMENTAL
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif
#ifdef SQUAREHAMILTONIAN
  if ( g_rank == 0 ) cout << "I am " << GA::nodeid() << " INNER JCHUNK SQUARE DIRECT HAMILTONIAN " << endl;
#endif
#ifdef NEWORBMAP
  GA::error(" cannot use constructGlobalHamiltonianAltOrbMapDirectChunkJmap method with NEWORBMAP",-1);
#endif
#ifdef STATIC
  GA::error(" Cannot do DIRECT and STATIC since g_det is not created",1);
#endif
  
  if ( GA::nodeid() ==0 ) cout << "FETCH jmaps as a subset " << endl;
  
  const int numjchunks=l_numJchunksHbuild; // indicator of chunkiness of J space 
  
  vector<int> list; // list of identically non-zero js per i
  
  int chunk_width;
  int l_maxspatials = l_deters->fetchMaxSpatials();
  
  (l_maxspatials % numjchunks==0)? chunk_width = l_maxspatials / numjchunks: chunk_width = 1 + (l_maxspatials / numjchunks);
  if ( GA::nodeid() ==0 ) cout << "INNER Jchunk stats num: "<<numjchunks<<" width " << chunk_width << endl;
  int jchunk_index=0; // these track the actual real index numbers (0,maxefs-1)
  
  double jfetchTime = 0.0; 
  
  pair<int,double> info1;
  
  long totalElems=0;
  int block_status;
  int num_jchunks=0;
  
  double compareT=0.0, fetchT=0.0;
  
  int totalFetches=0; // Used for Tuning studies
  
  if ( g_rank == 0 ) cout << "INNER JSUBSET CHUNK: NXTVAL based load balancing scheme being used " << endl;

  // Build Hamiltonian from scratch
  // Just in case afterall we do an accumulate 

  flushH(); // this simply zeros out cimat,icicol, num and the supp spaces
  
  PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
  jobs.initNxtask();
  long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
  int ilo = 0;
  int ihi = l_maxspatials;
  
  JOUTFG fetchi; 
  JOUTFG fetchj;
  
  // Build Hamiltonian from scratch
  // Just in case afterall we do an accumulate 
  // flushH(); // this simply zeros out cimat,icicol, num and the supp spaces

  int countme=0;
  int countrun=0;
  
  vector<COMPRESS_SIZE> imap( l_nbf+2  ); // a spatial configuration
  vector<COMPRESS_SIZE> jmap( l_nbf+2  );
  
  vector<COMPRESS_SIZE> jsubsetlist; // jmap subset buffer gets sized in fetchJsubsetlist 
  vector<int> jindexlist;
  
  for(int iconf = ilo; iconf < ihi; ++iconf ) {
    
    int numElems = 0;
    
    if ( nexttask == iconf ) { // in the loop stride is implemented 
      // do this above
      
      vector<vector<double> > cimat(1);
      vector<vector<int> > icicol(1);
      vector<int> number(1);
      
      temp2 = psociTime();
      list.clear();
      
      double temp3 = psociTime();
      l_deters->fetchOrbMap( iconf+1, imap);
      fetchT += psociTime() - temp3;
      temp3 = psociTime();
      
      // Find which local Js have nonzero overlap with the current imap
      /* Here we may only have part of the Js. so we need to ensure we filter correctly
       */
      // new countrun could be zero so we need to account for it.
      
      // move this up in the code
      temp2 = psociTime();
      l_deters->computeConfigsGA( iconf+1, imap, fetchi ); // also does a fetchOrb underneath
      fetchTime += psociTime() - temp2;
      ++totalFetches;
      
      int nsefi = fetchi.nsefi;
      if ( cimat.size() != nsefi ) {
	cimat.resize( nsefi );
	icicol.resize( nsefi );
	number.resize( nsefi );
      }
      
      for(int jchunk=0; jchunk < numjchunks; ++jchunk ) {
	
	++num_jchunks;
	
	double tempj = psociTime();
	int jnumlist = fetchJsubsetlist( jchunk, jindexlist, jsubsetlist, chunk_width, l_maxspatials ); 
	jfetchTime += psociTime() - tempj;
	
	temp3 = psociTime();
	int currentrun =  l_deters->compareLocalJSubsetOrbMap(jnumlist, iconf+1, imap, jindexlist, jsubsetlist, list); // 
	countrun += currentrun;
	//cout << GA::nodeid() << " " << iconf << " " << "subset returned " << countrun << endl;
	orbmapTime += psociTime() - temp2;
	
	vector<int>::iterator jit;
	for( jit=list.begin(); jit!=list.end();++jit) { 
	  int jconf = (*jit);
	  
	  temp2 = psociTime();
	  l_deters->computeConfigsGA( jconf+1, jmap, fetchj );
          fetchTime += psociTime() - temp2;
          ++totalFetches; // comment me out eventually
	  
	  //SOCI always clears sefsef
	  
	  ++allcount;
	  temp2 = psociTime();
	  
	  block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj, sefsef );
	  sociTime += psociTime() - temp2;
	  
	  if ( block_status ) {
	    int numTest = packSefsef( fetchi, fetchj , sefsef, cimat, icicol, number );
	    numElems += numTest;
	  }
        } // jit
      } // jchunk
      
      /* Update CIMAT. For a given I we need ALL contributing J values */
      if ( numElems > 0 ) { // For this to not be true would be somewhat rare
	totalElems += numElems;
	
	/* since we repeat loop over iconf we must be able to accumulate into H and not just put
	 */
#ifdef SUPPLEMENTAL
        temp2 = psociTime();
        uploadCIMATtoGASupplemental( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
        uploadTime += psociTime() - temp2;
#else
        uploadCIMATtoGA( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
#endif
      }
      nexttask = jobs.nxtask( g_size );
    } // inexttask  
  } // iconf
  
  jobs.destroyNxtask(); //Close'er up
  
  if ( GA::nodeid() == 0 ) {
#ifdef SQUAREHAMILTONIAN
  cout << GA::nodeid() << " square hamiltonian times " << endl;
#else 
  cout << GA::nodeid() << " trianglular hamiltonian times " << endl;
#endif
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "computeConfigsTime is " << fetchTime<<endl;
  cout << "allcount is " << allcount<<endl; //this is correct
  cout << "total computes is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
  cout << "orbmap time is " << orbmapTime << endl;
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  cout << "fetcT is " << fetchT << " compreT is " << compareT << endl;
  cout << "num jchunk is " << num_jchunks << endl;
  }
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif   
  
  return( numberNonzeroElements( totalElems ) );
}


// This version of Exp has an alternative model for fetching vectors
/*
long PsociGAhamiltonian::constructGlobalHamiltonianExpAlternative()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  
  cout << "I am " << GA::nodeid() << " at the exp-ALTERNATIVE  SPLIT hamiltonian " << endl;
  GA::SERVICES.sync();
  
  int counter=0;
  long allcount=0;

  double temp, temp2;
  double jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;
  
#ifdef PREALLOCATE
  int l_maxlength = l_deters->fetchMaxLength();
  vector<COMPRESS_SIZE> buffer( l_granularity * l_maxlength );
#endif

#ifdef SUPPLEMENTAL
  // Initialize suppval;
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif
  
#ifdef STATIC
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposiiton is contiguous. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used " << endl;
  jobs.initNxtask();// Initialize jobs.nxtask
  long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
  long inexttask=-1;
  int ilo = 0;
  int ihi = l_maxspatials;
#endif

#ifdef PREALLOCATE
  if ( g_rank == 0 ) cout << " USING PREALLOCATED FETCH METHODS " << endl;
#endif
  
  pair<int,double> info1;

// Preallocate space for fetching determinant information from GA

  JOUTFG fetchi;

#ifdef PREALLOCATE
  l_deters->preallocateTempJoutfgObject( fetchi );
  vector<JOUTFG> fetchj( l_granularity, fetchi );
#else
  vector<JOUTFG> fetchj;
#endif

  //   Put the fetch configs process here.......
  //   for(int iconf=0; iconf< l_maxspatials; ++iconf ) {
  //   for(int iconf=g_rank; iconf< l_maxspatials; iconf += g_size ) { 
  
  int numElems=0;
  long totalElems=0;
  int block_status;
  
  // Build Hamiltonian from scratch

  int istart, jfinal;

  for(int iconf = ilo; iconf < ihi; ++iconf ) {  
    
#ifndef STATIC
    ++inexttask; 
    if ( inexttask == nexttask ) {
//    cout << GA::nodeid() << " Doing a nexttask I J are " << iconf <<  " " << inexttask << " " << nexttask << endl;
#endif
      
#ifdef PREALLOCATE
      l_deters->fetchAndUnpackDeterminantDataPreallocate( iconf+1, info1, buffer, fetchi );
#else
      l_deters->fetchAndUnpackDeterminantData( iconf+1, info1, fetchi );
#endif
 
      // One configuration at a time for now
      vector<vector<double> > cimat( fetchi.nsefi );
      vector<vector<int> > icicol( fetchi.nsefi );
      vector<int> number( fetchi.nsefi );
      
      //Need to block these gets up into larger groups l_granularity in width

        numElems = 0;
        istart = -1;
        for(int jconf=0; jconf <= iconf; ++jconf )  
	{ 
             temp2 = psociTime();
             if ( istart == -1 ) {
                jfinal = min( jconf + l_granularity, l_maxspatials );
#ifdef PREALLOCATE
	        l_deters->fetchAndUnpackDeterminantData2DPreallocate( jconf+1, jfinal, buffer, fetchj );
#else
                fetchj.clear();
                l_deters->fetchAndUnpackDeterminantData2D( jconf+1, jfinal, fetchj );
#endif
                istart = 0;
             }
             fetchTime += psociTime() - temp2;

             ++allcount;

             temp2 = psociTime();
	     block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj[istart], sefsef );
              sociTime += psociTime() - temp2;
          
             if ( block_status ) {
                      numElems += packSefsef( fetchi, fetchj[istart] , sefsef, cimat, icicol, number );
             }

        ++istart;
        if ( istart >= l_granularity ) istart = -1;
	}
      
      if ( numElems > 0 ) { // For this to not be true would be somewhat rare
        totalElems += numElems;
#ifdef SUPPLEMENTAL
	uploadCIMATtoGASupplemental( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
#else
	uploadCIMATtoGA( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
#endif
      }
#ifndef STATIC
      nexttask = jobs.nxtask( g_size );
    }
#endif
    
  }
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "fetchTime is " << fetchTime<<endl;

  // Print all of them for now since they indicate load balance
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  GA::SERVICES.sync();
  
#ifdef SUPPLEMENTAL
  //Destroy suppval;
  destroySuppVal();
#endif
  
#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif
  
#ifdef RCRMONITOR
  stats =  RCRpopLocation();
  stats =  RCRstopLogging();
#endif
  
  //return( totalElems ); // switched for performance measurment
  return( numberNonzeroElements( totalElems ) );
}
*/

//Timer wrapped version 
/*
long PsociGAhamiltonian::constructGlobalHamiltonianExpAll( pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalHamiltonianExpAll();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid(); 
  return( status );
} 
*/

/* This version fully prefetches ALL joutfg data
   It is fast but jhow big can it go?
*/
/*
long PsociGAhamiltonian::constructGlobalHamiltonianExpAll()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  
  cout << "I am " << GA::nodeid() << " at the GOT ALL SPLIT hamiltonian " << endl;
  GA::SERVICES.sync();
  
  int counter=0;
  long allcount=0;

  double temp, temp2;
  double jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;

#ifdef SUPPLEMENTAL
  // Initialize suppval;
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model " << endl;
  initSuppVal();
#endif
  
#ifdef STATIC
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposiiton is contiguous. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used " << endl;
  jobs.initNxtask();// Initialize jobs.nxtask
  long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
  long inexttask=-1;
  int ilo = 0;
  int ihi = l_maxspatials;
#endif
  
  pair<int,double> info1;

  JOUTFG fetchi, fetchj;
  vector<JOUTFG> fetchtemp;
  l_deters->fetchAndUnpackDeterminantData2D( 1, l_maxspatials, fetchtemp );
  

  //   Put the fetch configs process here.......
  //   for(int iconf=0; iconf< l_maxspatials; ++iconf ) {
  //   for(int iconf=g_rank; iconf< l_maxspatials; iconf += g_size ) { 
  
  int numElems=0;
  long totalElems=0;
  int block_status;
  
  // Build Hamiltonian from scratch

  for(int iconf = ilo; iconf < ihi; ++iconf ) {  
    
#ifndef STATIC
    ++inexttask; 
    if ( inexttask == nexttask ) {
//    cout << GA::nodeid() << " Doing a nexttask I J are " << iconf <<  " " << inexttask << " " << nexttask << endl;
#endif

      // One configuration at a time for now
      vector<vector<double> > cimat(fetchtemp[iconf].nsefi );
      vector<vector<int> > icicol(fetchtemp[iconf].nsefi );
      vector<int> number(fetchtemp[iconf].nsefi );
      
      //Need to block these gets up into larger groups l_granularity in width

        numElems = 0;
        for(int jconf=0; jconf <= iconf; ++jconf )  
	{ 

	     ++allcount;
             temp2 = psociTime();
	     block_status = l_gaHamiltonian->generateSOCIblocks( fetchtemp[iconf], fetchtemp[jconf], sefsef );
             sociTime += psociTime() - temp2;
          
	     if ( block_status ) {
                      numElems += packSefsef( fetchtemp[iconf], fetchtemp[jconf] , sefsef, cimat, icicol, number );
             }
	}
      
      if ( numElems > 0 ) { // For this to not be true would be somewhat rare
        totalElems += numElems;
#ifdef SUPPLEMENTAL
	uploadCIMATtoGASupplemental( l_nsef[ fetchtemp[iconf].index - 1 ] , cimat, icicol, number );
#else
	uploadCIMATtoGA( l_nsef[ fetchtemp[iconf].index - 1 ] , cimat, icicol, number );
#endif
      }
#ifndef STATIC
      nexttask = jobs.nxtask( g_size );
    }
#endif
    
  }
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "fetchTime is " << fetchTime<<endl;

  // Print all of them for now since they indicate load balance
  cout << GA::nodeid() << " Done after loops " <<  totalElems << endl;
  GA::SERVICES.sync();
  
#ifdef SUPPLEMENTAL
  //Destroy suppval;
  destroySuppVal();
#endif
  
#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif
  
#ifdef RCRMONITOR
  stats =  RCRpopLocation();
  stats =  RCRstopLogging();
#endif
  
  //return( totalElems ); // switched for performance measurment
  return( numberNonzeroElements( totalElems ) );
}
*/

//Local call
int PsociGAhamiltonian::fetchNumberSefs()
{
  return( l_maxsef );
}

int PsociGAhamiltonian::fetchNBF()
{
  return( l_nbf );
}

// Local  Specifies what confs are LOCAL not neede anymore
// Not really used anymore

int PsociGAhamiltonian::setLocalConfigs( vector<JOUTFG> & confs )
{
  l_confs = confs;
  l_numconfs = l_confs.size();
  return( l_numconfs );
}

/*
  Local: On return, all computed sefsef arrays are shoved into 
  GA space (cimat, icicol, number) and can be processed by the
  H*V methods
*/

// Note: The local confs must interact with ALL j confs (in GA)
// Only for developmental purposes


// This may NOT EVEN work anymore ands was used for early validation
/*
int PsociGAhamiltonian::constructLocalHamiltonian()
{
  int totalCount=0;
  
  // Compute offdiagonal terms
  if ( l_numconfs > 1 ) {
    for(int i=0; i< l_numconfs; ++i ){
      for(int j=0; j<i; ++j ) {
	++totalCount;
	cout << "XXXXXXXXXXXXXXXXXX NEXT PAIR XXXX " << endl;
	cout << " Pair are i and j = " << i+1 << " " << j+1 << endl;
	l_gaHamiltonian->generateSOCIblocks( l_confs[i], l_confs[j], sefsef );
      }
    }
  }
  // Compute diagonal blocks
  if ( l_numconfs > 0 ) {
    for(int i=0; i< l_numconfs; ++i ){
      ++totalCount;
      cout << "XXXXXXXXXXXXXXXXXX NEXT DIAGONAL XXXX " << endl;
      cout << " Pair are i and j = " << i+1 << " " << i+1 << endl;
      l_gaHamiltonian->generateSOCIblocks( l_confs[i], l_confs[i], sefsef );
    }
  }
  return( totalCount);
}
*/

/* Every node needs a global map to where in the hamiltonian(sef,sef) a particular i and j config starts.
   We take the ordered data in l_nsef, generate cummulative sums and reuse l_nsef
   
   I think we could scramble order to better balance H
   This is a Collective call
*/

int PsociGAhamiltonian::computeNsefEntryPoints( vector<int> & nsef )
{
#if 0
  cout << "Node " << GA::nodeid() << " at the NsefEntryPoints " << endl;
#endif
  
  int num_nsefi = 0;
  if ( nsef.size() != l_maxspatials ) {
    cerr << "nsef.size() != l_maxspatials " << nsef.size() << " " << l_maxspatials  << endl;
    GA::error( "nsef.size() != l_maxspatials ", nsef.size() );
  }
  
  //Compute entry points for the set of sefs for configuration i. (absolute index)
  if ( GA::nodeid() == 0 ) {
    int nsefi;
    for(int i=0; i< l_maxspatials; ++i) {
      nsefi = nsef[i];
      nsef[i] = num_nsefi;
      num_nsefi += nsefi;
    }
  }
  //brdcst results out to all.
  GA::brdcst( &nsef[0], nsef.size()*sizeof(int), 0) ;
  
  l_nsef = nsef;  
  return( l_nsef.size() );
}

//Local call: generally for debugging purposes
void PsociGAhamiltonian::printNsefEntryPoints()
{
  vector<int>::iterator it;
  cout << "Current set of Entry Points for l_nsef " << endl;
  int index=0;
  for(it=l_nsef.begin(); it!=l_nsef.end(); ++it) {
    cout << " " << index+1 << " " << (*it) << endl;
  }
}

/* Take the current sefsef and pack it into a dense array.  The (ii) terms will be reused. We can use local indexed
   values of ii (1,nsefi). Only just before pushing to GA do we need to offset this to an absolute value.
   
   Many jj terms will get pushed into the dense array by repeatedly calling this routine.
   
   l_nsef[i] begins the set of sefs for the configuration i.  
*/

/* cimat:  Local array containing H elements (ii,jj)
   icicol: Local array that maps jj to jjabsolute
   number: Local array that contains total number of non-zero elements in cimat[ii][*]
   
   maxElems reports the current (local) maximum number of JJ per II. Presumably for the
   given config; I
*/

int PsociGAhamiltonian::packSefsef( JOUTFG & iconf, JOUTFG & jconf , vector<double> & sefsef,
				    vector<vector<double> > & cimat, vector<vector<int> > & icicol, vector<int> & number )
{
  int idex = iconf.index;
  int jdex = jconf.index;
  
  int nsefi = iconf.nsefi;
  int nsefj = jconf.nsefi;
  
  int isef = l_nsef[ idex - 1 ]; // These are the absolute entry points for the start of confgs I
  int jsef = l_nsef[ jdex - 1 ];
  int jmax;
  
//#ifdef DETAILEDHAMCHECK
  if ( sefsef.size() != nsefi * nsefj ) {
    cerr << "Erroneous sizes in packSefsef: Probable mismatch: size is " << sefsef.size() << "nsefi and nsefj are " << nsefi << " " << nsefj <<  endl;
    GA::error(  "Erroneous sizes in packSefsef: Probable mismatch:", -1);
  }
  
  if ( cimat.size() != nsefi || icicol.size() != nsefi || number.size() != nsefi ) {
    cerr << "pack matrices not of the correct leading dimension " << endl;
    GA::error( "pack matrices not of the correct leading dimension ", cimat.size() );
  }
//#endif
  
  // cimat[nsefi][*], icicol[nsefi][*], number[nsefi];
  
  int ncol;
  int indSefsef;
  
  // Perhaps we should estimate the length of row to minimize the push_backs 
  
  int numElems = 0;

  for(int i=0; i< nsefi; ++i ) {
#ifdef SQUAREHAMILTONIAN
    jmax=nsefj-1;
#else
    (isef == jsef) ? jmax=i : jmax=nsefj-1;
#endif

/* get the diagonal elements here 
*/
    ncol=0;
    indSefsef = i * nsefj;
    for (int j=0; j<= jmax; ++j ) {
      if ( abs( sefsef[ indSefsef + j ] ) >= MIN_DOUBLE_VALUE || (isef == jsef && i == j) ) {
	cimat[i].push_back( sefsef[ indSefsef + j ]);
	icicol[i].push_back( ncol+jsef ); 
	++number[i]; 
        ++numElems;
	//cout << "CIMAT indexe sefsef  value is " << indSefsef<<" "<<sefsef[ indSefsef + j ]<<" "<<number[i] << " i and j are " << i << " " << j << endl;
	//cout << " absolute J index is " << ncol+jsef; 
	//cout << " individual: i and number[i] are " << i << " " << number[i] << endl;
      }
      ++ncol;       
    }
  }
  return( numElems );
}

// Collective call: Single GA space - should merge this CODE with the split version
void PsociGAhamiltonian::createGAhamiltonian()
{
  const int DATA_NDIM = MAX_CI_VECTOR_DIMS;
  int dims[DATA_NDIM];
  int chunk[DATA_NDIM];
  
  const int DISTRIBUTE_EVENLY = -1; //TODO need to optimize this at some point
  
  dims[0] = l_maxsef;
  dims[1] = l_maxsparse;
  
  chunk[0] = DISTRIBUTE_EVENLY;
  //  chunk[1] = DISTRIBUTE_EVENLY;
  chunk[1] = l_maxsparse; //Just like the original - horizontal strips
  

  // No guarantee that C_DBL/C_INT equal double/int GA should provide a sizeof method....
  long totalWords = l_maxsparse * (long) l_maxsef;
  long totalMemory = totalWords * sizeof(double) + totalWords * sizeof(int) + l_maxsef * sizeof(int);
  
  long size = GA::nodes();
  long divisor = 1000L*1000L;
  if ( GA::nodeid() == 0 ) {
    cout << "l_maxsef for the GA is " << l_maxsef << endl;
    cout << "The estimated amount of global Hamiltonian memory (B) to allocate is " << totalMemory << endl;
    cout << "The estimated amount of per-core Hamiltonian memory (MB) to allocate is " << totalMemory / (divisor * size ) << endl;
  }
  
  char local_title[] = "g_cimat";
  //cout << "Creating g_cimat of size " << totalWords * sizeof(double) << endl; 
  l_g_cimat = GA::createGA( C_DBL, DATA_NDIM, dims, (char *) local_title, chunk);
  l_g_cimat->zero();
  
  char local_title2[] = "g_icicol";
  l_g_icicol = GA::createGA( C_INT, DATA_NDIM, dims, (char *) local_title2, chunk);
  l_g_icicol->zero();
  
  dims[1] = 1;
  chunk[1] = -1;
  char local_title3[] = "g_num";
  l_g_number = GA::createGA( C_INT, DATA_NDIM, dims, (char *) local_title3, chunk);
  l_g_number->zero();
  
  char local_title4[] = "g_diag_sef";
  l_g_diag_sef = GA::createGA( C_DBL, DATA_NDIM, dims, (char *) local_title4, chunk);
  l_g_diag_sef->zero();
}

/* this is used oNLY with the direct matrix-vector product method
*/
void PsociGAhamiltonian::createGAhamiltonianDiagOnly()
{
  const int DATA_NDIM = MAX_CI_VECTOR_DIMS;
  int dims[DATA_NDIM];
  int chunk[DATA_NDIM];
  
  const int DISTRIBUTE_EVENLY = -1; //TODO need to optimize this at some point
  
  dims[0] = l_maxsef;
  dims[1] = l_maxsparse;
  
  chunk[0] = DISTRIBUTE_EVENLY;
  //  chunk[1] = DISTRIBUTE_EVENLY;
  chunk[1] = l_maxsparse; //Just like the original - horizontal strips
  

  // No guarantee that C_DBL/C_INT equal double/int GA should provide a sizeof method....
  long totalWords = l_maxsparse * (long) l_maxsef;
  long totalMemory = totalWords * sizeof(double) + totalWords * sizeof(int) + l_maxsef * sizeof(int);
  
  long size = GA::nodes();
  long divisor = 1000L*1000L;
  if ( GA::nodeid() == 0 ) {
    cout << "l_maxsef for the GA is " << l_maxsef << endl;
    cout << "The estimated amount of global Diag Hamiltonian memory (B) to allocate is " << totalMemory << endl;
    cout << "The estimated amount of per-core Diag Hamiltonian memory (MB) to allocate is " << totalMemory / (divisor * size ) << endl;
  }

  char local_title4[] = "g_diag_sef";
  l_g_diag_sef = GA::createGA( C_DBL, DATA_NDIM, dims, (char *) local_title4, chunk);
  l_g_diag_sef->zero();
}

// Collective call : open up multiple GA spaces: Split data model
void PsociGAhamiltonian::createGAhamiltonian2Array()
{
  const int DATA_NDIM = MAX_CI_VECTOR_DIMS;
  int dims[DATA_NDIM];
  int chunk[DATA_NDIM];
  
  const int DISTRIBUTE_EVENLY = -1; //TODO need to optimize this at some point
  
  // Open up MAIN array 
  
  dims[0] = l_maxsef;
  dims[1] = l_maxsparse;
  chunk[0] = DISTRIBUTE_EVENLY; // Need to do performance tuning on this

  chunk[1] = l_maxsparse; //Just like the original
  
  long totalWords = l_maxsparse * (long) l_maxsef;
  long totalMemory = totalWords * sizeof(double) + totalWords * sizeof(int) + l_maxsef * sizeof(int);
  
  long size = GA::nodes();
  long divisor = 1000L*1000L;
  

  if (GA::nodeid() == 0 ) cout << "Entered create H"<< endl;

  if ( GA::nodeid() == 0 ) {
    cout << GA::nodeid() <<  " l_maxsef for the MAIN GA is " << l_maxsef << endl;
    cout << GA::nodeid() <<  " The estimated amount of main global Hamiltonian memory (B) to allocate is " << totalMemory << endl;
    cout << GA::nodeid() <<   " The estimated amount of main per-core Hamiltonian memory (MB) to allocate is " << totalMemory / (divisor * size ) << endl;
  }
  
if (GA::nodeid() == 0 ) cout << " g_cimat"<<endl;
  char local_title[] = "g_cimat";
  l_g_cimat = GA::createGA( C_DBL, DATA_NDIM, dims, (char *) local_title, chunk);
  l_g_cimat->zero();

if (GA::nodeid() == 0 ) cout << " g_icicol"<<endl;
  char local_title2[] = "g_icicol";
  l_g_icicol = GA::createGA( C_INT, DATA_NDIM, dims, (char *) local_title2, chunk);
  l_g_icicol->zero();
  
  dims[1] = 1;
  chunk[1] = 1;
if (GA::nodeid() == 0 ) cout << " g_num"<<endl;
  char local_title3[] = "g_num";
  l_g_number = GA::createGA( C_INT, DATA_NDIM, dims, (char *) local_title3, chunk);
  l_g_number->zero();
  
if (GA::nodeid() == 0 ) cout << " g_diag_sef"<<endl;
  char local_title4[] = "g_diag_sef";
  l_g_diag_sef = GA::createGA( C_DBL, DATA_NDIM, dims, (char *) local_title4, chunk);
  l_g_diag_sef->zero();
  
  // Open up Supplementary array 
  
  dims[0] = l_maxwidth;
  dims[1] = l_maxsparse_supp;
  chunk[0] = DISTRIBUTE_EVENLY;
  chunk[1] = l_maxsparse_supp;
  
  totalWords += l_maxsparse_supp * (long) l_maxwidth;
  totalMemory += totalWords * sizeof(double) + totalWords * sizeof(int) + l_maxwidth * sizeof(int);
  
//  cout << "l_maxwidth for the SUPPLEMENTAL GA is " << l_maxwidth << endl;
  if ( GA::nodeid() == 0 ) {
    cout << "The estimated amount of TOTAL (main+supplemental)  global (MAIN + SUPP) Hamiltonian memory (B) to allocate is " << totalMemory << endl;
    cout << "The estimated amount of TOTAL (main+supplemental)  per-core (MAIN + SUPP) Hamiltonian memory (MB) to allocate is " << totalMemory / (divisor*size) << endl;
  }
  
if (GA::nodeid() == 0 ) cout << " g_cimat_sup"<<endl;

  char local_title5[] = "g_cimat_sup";
//  cout << "Creating g_cimat_sup of size " << totalWords * sizeof(double) << endl;
  l_g_cimat_supp = GA::createGA( C_DBL, DATA_NDIM, dims, (char *) local_title5, chunk);
  l_g_cimat_supp->zero();
  
if (GA::nodeid() == 0 ) cout << " g_icicol_sup"<<endl;

  char local_title6[] = "g_icicol_sup";
  l_g_icicol_supp = GA::createGA( C_INT, DATA_NDIM, dims, (char *) local_title6, chunk);
  l_g_icicol_supp->zero();
  
  dims[1] = 1;
  chunk[1] = 1;
if (GA::nodeid() == 0 ) cout << " g_num_sup"<<endl;

  char local_title7[] = "g_num_sup";
  l_g_number_supp = GA::createGA( C_INT, DATA_NDIM, dims, (char *) local_title7, chunk);
  l_g_number_supp->zero();
}


//Local call
void PsociGAhamiltonian::hamiltonianDistribution()
{
  if (GA::nodeid() == 0 ) cout << "Print out GA Hamiltonian objects distribution " << endl;
  l_g_cimat->printDistribution();
  l_g_icicol->printDistribution();
  l_g_number->printDistribution();
#ifdef SUPPLEMENTAL
  l_g_cimat_supp->printDistribution();
  l_g_icicol_supp->printDistribution();
  l_g_number_supp->printDistribution();
#endif
}

// Local: single GA array version
// Local call - based a specific distribution whereby we have all jj terms for the current configuration; I (ii=1:nsefi)
// currently expect NO spaces in the data set

// Also populate the g_diag_sef terms for use by the preconditioner Diag[isef][sef]= value

// Local 
int PsociGAhamiltonian::uploadCIMATtoGA( int isef, vector<vector<double> > &cimat, vector<vector<int> > &icicol, vector<int> & number )
{
  // I think we should COMBINE all [ii] terms into a single row. Then upload the whole thing.
  // We can keep [jj] the same. We will always know how wide the [ii] is.
  // TODO: Need a statistical analysis of the likely number of zeros in such a data object
  
  /* For all current load balancing algorithms, i=1,n, j =1,i. So that LAST entry in every buffer is jsef==isef.
   */
  
#ifdef DETAILEDHAMCHECK
  if ( number.size() != cimat.size() || number.size() != icicol.size() ) {
    cerr << " Incorrect sizes " << number.size() << endl;
    GA::error( " Incorrect sizes ", number.size() );
  }
#endif
  
  int nsefi = number.size();
  int irow;
  int lo[MAX_CI_VECTOR_DIMS];
  int hi[MAX_CI_VECTOR_DIMS];
  int n,n1;
  
  // For now push out a single I columns worth 
  
  double dValue=0.0;
  int jsef;
  
  for(int ii=0; ii< nsefi; ++ii ) {
    irow  = isef + ii;
    lo[0] = irow;
    lo[1] = 0; // 
    hi[0] = irow; //
    int numii = number[ii];
    hi[1] = numii-1;
    
    if ( numii > 0 ) {

    if ( numii-1 >= l_maxsparse ) {
      cerr << "number[ii]-1 exceeds specified maxsparse: Will not be able to store to GA." << endl;
      cerr << " number[ii]-1 = " << number[ii]-1 << " maxsparsity is " << l_maxsparse << endl;
      GA::error( "number[ii]-1 exceeds specified maxsparse", numii-1);
    }
    
    n = numii; // Not sure of this !
    n1 = 1;
    
    if ( cimat[ii].size() > 0 ) {

#ifdef GATRACE
      trace_stime_();
#endif
      l_g_cimat->put( lo, hi, &cimat[ii][0], &n);
      l_g_icicol->put( lo, hi, &icicol[ii][0], &n);
      
      hi[1]=0;
      
//      n = numii; // Not sure of this either !!
      l_g_number->put( lo, hi, &numii, &n1);
      
      /* DO some checks for j in case someone changes the load balancing scheme....
	 If this become a bottleneck then simply shove cimat value to GA
      */
#ifdef SQUAREHAMILTONIAN
      for(int jtest=0; jtest< numii-1; ++jtest ) {
         jsef = icicol[ii][jtest];
          if (jsef == irow ) {
          dValue = cimat[ii][ jtest ];
          break;
      }
#else
      jsef = icicol[ii][ numii-1 ]; //This is already an ABSOLUTE location
      if (jsef == irow ) {
	//cout << GA::nodeid() << " PUTTING a diag " <<  cimat[ii][ number[ii]-1 ] << endl;
	dValue = cimat[ii][ numii-1 ]; 
#endif

     cout << "DIAG " <<  dValue << endl;
	l_g_diag_sef->put( lo, hi, &dValue, &n1 );
      }
#ifdef GATRACE
      trace_etime_();
      int op = 1;
      int * pop = & op;
      trace_genrec_( (int *) l_g_cimat, &lo[0], &hi[0], &lo[1],&hi[1], pop);
      trace_genrec_( (int *) l_g_icicol, &lo[0], &hi[0], &lo[1],&hi[1], pop);
#endif
    }
   } // numii > 0 
  }
  return(0);
}
/* Method to ONLY store H diag elementsa to GA
*/
int PsociGAhamiltonian::uploadCIMATdiagElemsOnlytoGA( int isef, vector<vector<double> > &cimat, vector<vector<int> > &icicol, vector<int> & number )
{
#ifndef DIRECTMATRIXPRODUCTS
  GA::error(" sociGAhamiltonian::uploadCIMATdiagElemsOnlytoGA only used with -DDIRECTMATRIXPRODUCTS",1);
#endif
  
  /* For all current load balancing algorithms, i=1,n, j =1,i. So that LAST entry in every buffer is jsef==isef.
   */
  
  int nsefi = number.size();
  int irow;
  int lo[MAX_CI_VECTOR_DIMS];
  int hi[MAX_CI_VECTOR_DIMS];
  int n,n1;
  
  // For now push out a single I columns worth 
  
  double dValue=0.0;
  int jsef;
  
  for(int ii=0; ii< nsefi; ++ii ) {
    irow  = isef + ii;
    
    if ( cimat[ii].size() <= 0 ) GA::error("cimat[ii].size() <= 0  should never happen with a diag",1);
      
    lo[0] = irow;
    lo[1] = 0; // 
    hi[0] = irow; //
    hi[1] = 0;
    int numii = number[ii];

    jsef = icicol[ii][ numii-1 ]; //This is already an ABSOLUTE location

    if (jsef == irow ) {
     dValue = cimat[ii][ numii-1 ]; 
     //cout << "ONLY DIAG " <<  dValue << endl;
	l_g_diag_sef->put( lo, hi, &dValue, &n1 );
      }
  }
  return(0);
}

/* A new version of upload for the distributed orbmap methods
   We use accumulate here since we may now be looping over jchunks
*/
int PsociGAhamiltonian::uploadCIMATtoGAAcc( int isef, vector<vector<double> > &cimat, vector<vector<int> > &icicol, vector<int> & number )
{
  // I think we should COMBINE all [ii] terms into a single row. Then upload the whole thing.
  // We can keep [jj] the same. We will always know how wide the [ii] is.
  // TODO: Need a statistical analysis of the likely number of zeros in such a data object
  
  /* For all current load balancing algorithms, i=1,n, j =1,i. So that LAST entry in every buffer is jsef==isef.
   */
  
  double dOne=1.0;
  int iOne=1;

#ifdef DETAILEDHAMCHECK
  if ( number.size() != cimat.size() || number.size() != icicol.size() ) {
    cerr << " Incorrect sizes " << number.size() << endl;
    GA::error( " Incorrect sizes ", number.size() );
  }
#endif
  
  int nsefi = number.size();
  int irow;
  int lo[MAX_CI_VECTOR_DIMS];
  int hi[MAX_CI_VECTOR_DIMS];
  int n,n1;
  
  // For now push out a single I columns worth 
  
  double dValue=0.0;
  int jsef;
  
  for(int ii=0; ii< nsefi; ++ii ) {
    irow  = isef + ii;
    lo[0] = irow;
    lo[1] = 0; // 
    hi[0] = irow; //
    int numii = number[ii];
    hi[1] = numii-1;
    
    if ( numii > 0 ) {


    if ( numii-1 >= l_maxsparse ) {
      cerr << "number[ii]-1 exceeds specified maxsparse: Will not be able to store to GA." << endl;
      cerr << " number[ii]-1 = " << number[ii]-1 << " maxsparsity is " << l_maxsparse << endl;
      GA::error( "number[ii]-1 exceeds specified maxsparse", numii-1);
    }
    
    n = numii; // Not sure of this !
    n1 = 1;
    
    if ( cimat[ii].size() > 0 ) {

#ifdef GATRACE
      trace_stime_();
#endif
      l_g_cimat->acc( lo, hi, &cimat[ii][0], &n, &dOne);
      l_g_icicol->acc( lo, hi, &icicol[ii][0], &n, &iOne);
      
      hi[1]=0;
      
//      n = numii; // Not sure of this either !!
      l_g_number->acc( lo, hi, &numii, &n1, &iOne);
      
      /* DO some checks for j in case someone changes the load balancing scheme....
	 If this become a bottleneck then simply shove cimat value to GA
      */
#ifdef SQUAREHAMILTONIAN
      for(int jtest=0; jtest< numii-1; ++jtest ) {
         jsef = icicol[ii][jtest];
          if (jsef == irow ) {
          dValue = cimat[ii][ jtest ];
          break;
      }
#else
      jsef = icicol[ii][ numii-1 ]; //This is already an ABSOLUTE location
      if (jsef == irow ) {
	//cout << GA::nodeid() << " PUTTING a diag " <<  cimat[ii][ number[ii]-1 ] << endl;
	dValue = cimat[ii][ numii-1 ]; 
#endif

	l_g_diag_sef->put( lo, hi, &dValue, &n1 );
      }
#ifdef GATRACE
      trace_etime_();
      int op = 1;
      int * pop = & op;
      trace_genrec_( (int *) l_g_cimat, &lo[0], &hi[0], &lo[1],&hi[1], pop);
      trace_genrec_( (int *) l_g_icicol, &lo[0], &hi[0], &lo[1],&hi[1], pop);
#endif
    }
  }

  } // numii > 0
  return(0);
}

//A Collective call useful only for debugging
void PsociGAhamiltonian::dumpCimat()
{
  if ( GA::nodeid() == 0 ) {
    cout << "Dump out cimat,icicol,number matrices: I am " << GA::nodeid() << endl;
    cout << "CIMAT" << endl;
    l_g_cimat->print();
    cout << "ICICOL" << endl;
    l_g_icicol->print();
    cout << "NUMBER" << endl;
    l_g_number->print();
#ifdef SUPPLEMENTAL
    l_g_cimat_supp->print();
    l_g_icicol_supp->print();
    l_g_number_supp->print();
#endif
  } else {
    cout << "dumpCimat: only node==0 actually prints out the matrix " << endl;
  }
  GA::SERVICES.sync();
}

//Desired statistics for the global run: total number of processed H elements above threshold
long PsociGAhamiltonian::numberNonzeroElements( long l_number )
{
  char op[] = "+";
  long global_numElems  = l_number;
  GA::gop( &global_numElems, 1, op );
  return( global_numElems );
}

// Collective
// Not to be used for regular computations. This is simply acquire data on the number of non-zeros 

// A specialized method that simply prints out sparsity data: Not for production runs
// Specifically REMOVED the full matrix representation as memory needs got too big too fast

/* Not really intended to be used and in fact may not be broken. The internal data structures are
   not scalable and get very bery large on a single core.

   NO NEED TO USE THIS ANYMORE
*/
/*
long PsociGAhamiltonian::constructGlobalHamiltonianSparseStats()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  
  cout << "I am " << GA::nodeid() << " at the hamiltonian " << endl;
  GA::SERVICES.sync();
  
  // Each node constructs their lists of configurations to process.
  
  int nnodes = GA::nodes();
  
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  
  pair<int,double> info1;
  vector<JOUTFG> fetchi, fetchj;
  vector<int> grab1; // The set of configurations. All diagional and offdiagonal terms are computed.

// Too big
  vector<vector<int> > matrix( l_maxsef,vector<int>(l_maxsef,0) ); // A array to fill with 0s and 1s
  
  int numElems=0;
  long totalElems=0;
  
  long sumElems=0, sumNonzeros=0;
  
 //   Put the fetch configs process here.......
 //   for(int iconf=0; iconf< l_maxspatials; ++iconf ) {
 //   for(int iconf=g_rank; iconf< l_maxspatials; iconf += g_size ) { 
  
  for(int iconf=ihi; iconf>= ilo; --iconf ) {
    fetchi.resize(1);
    l_deters->fetchAndUnpackDeterminantData( iconf+1, info1, fetchi[0] );
    
    // One at a time for now
    vector<vector<double> > cimat(fetchi[0].nsefi );
    vector<vector<int> > icicol(fetchi[0].nsefi );
    vector<int> number(fetchi[0].nsefi );

    // Statistics data
    // pair<int,int>( totalnum, num nonZero) info

    vector<pair<long,long> > stats( fetchi[0].nsefi ); 
    
    //  We need to change this. Do not upload data to GA until ALL jconf have been processed.

    // Also we want to gang up these into a single call. We probably should simply grab whatever exists on the next node.......
    
    vector<JOUTFG>::iterator jit;
    for(int jconf=0; jconf <= iconf; ++jconf ) 
      {
	fetchj.resize(1);
	l_deters->fetchAndUnpackDeterminantData( jconf+1, info1, fetchj[0] );
	l_gaHamiltonian->generateSOCIblocks( fetchi[0], fetchj[0], sefsef );
	numElems = packSefsef( fetchi[0], fetchj[0] , sefsef, cimat, icicol, number, stats, matrix );
      }
    
    if ( totalElems > 0 ) { // For this to not ber true would be very rare indeed.
      uploadCIMATtoGA( l_nsef[ fetchi[0].index - 1 ] , cimat, icicol, number );
    }
    for(int k=0; k< fetchi[0].nsefi; ++k) {
      cout << "stats " << stats[k].first << " " << stats[k].second << endl;
      sumElems+= stats[k].first;
      sumNonzeros +=  stats[k].second;
    }
  }
  aggregateStatsData(l_maxsef, matrix );
  
  cout << "Print out MATRIX ii,jj,value " << matrix.size() << endl;
  if ( g_rank == 0 ) {
    for(int i=0; i< l_maxsef; ++i ) {
      for(int j=0; j< l_maxsef; ++j ) {
	cout << i << " " << j << " " << matrix[i][j] << endl;
      }
      cout << endl;
    }
  }
  GA::SERVICES.sync();
  
  return( totalElems ); // return LOCAL totalElems rather than global total elems produced by numberNonzeroElements method
  //  return( numberNonzeroElements( totalElems ) );
}
*/


// Only useful in collecting statistics for understanding nonzero behavior

int PsociGAhamiltonian::packSefsef( JOUTFG & iconf, JOUTFG & jconf , vector<double> & sefsef,
				    vector<vector<double> > & cimat, vector<vector<int> > & icicol, vector<int> & number, 
				    vector<pair<long,long> > & stats, vector<vector<int> > & matrix )
{
  int idex = iconf.index;
  int jdex = jconf.index;
  
  int nsefi = iconf.nsefi;
  int nsefj = jconf.nsefi;
  
  int isef = l_nsef[ idex - 1 ]; // These are the absolute entry points for the start of confgs I
  int jsef = l_nsef[ jdex - 1 ];
  int jmax;
  
#ifdef DETAILEDHAMCHECK
  if ( sefsef.size() != nsefi * nsefj ) {
    cerr << "Erroneous sizes in packSefsef: Probable mismatch: size is " << sefsef.size() << "nsefi and nsefj are " << nsefi << " " << nsefj <<  endl;
    GA::error( "Erroneous sizes in packSefsef:", sefsef.size() );
  }
  
  if ( cimat.size() != nsefi || icicol.size() != nsefi || number.size() != nsefi ) {
    cerr << "pack matrices not of the correct leading dimension " << endl;
    GA::error( "pack matrices not of the correct leading dimension ",cimat.size() );
  }
  
  if (stats.size() != nsefi ) {
    cerr << "Inconsistent size of stats = " << stats.size() << endl;
    GA::error( "Inconsistent size of stats ",stats.size() );
  }
#endif
  
  // cimat[nsefi][*], icicol[nsefi][*], number[nsefi];
  // int numb; //Total number of non-zero elements for each ii,jj

  int ncol;  
  int indSefsef;
  int numElems = 0;
  
  for(int i=0; i< nsefi; ++i ) {
#ifdef SQUAREHAMILTONIAN
    jmax=nsefj-1;
#else
    (isef == jsef) ? jmax=i : jmax=nsefj-1;
#endif

    ncol=0;
    for (int j=0; j<= jmax; ++j ) {
      indSefsef = i * nsefj;
      if ( abs( sefsef[ indSefsef + j ] ) >= MIN_DOUBLE_VALUE || (isef == jsef && i == j) ) {
	cimat[i].push_back( sefsef[ indSefsef + j ]);
	icicol[i].push_back( ncol+jsef ); 
	++number[i]; 
        ++numElems;
        matrix[isef+i][jsef+ncol] += number[i];
	//cout << "CIMAT value is " << sefsef[ indSefsef + j ] << " i and j are " << i << " " << j;
	//cout << " absolute J index is " << ncol+jsef; 
	//cout << " individual: i and number[i] are " << i << " " << number[i] << endl;
      }
      ++ncol;       
    }
    stats[i].first += ncol;
    stats[i].second = number[i];
  }
  return( numElems );
}
// Not for production use
// gather distributed sparsity stats back top node 0
// Eliminated the full matrix sized buffer as we need to look at larger problems

// Not intended toi be used. THis was part of the stats collection suite of methods and is not scalable
void  PsociGAhamiltonian::aggregateStatsData(int l_maxsef, vector<vector<int> > & matrix )
{
  vector<int> longArray( l_maxsef * l_maxsef );
  
  int index;
  for(int i=0; i< l_maxsef; ++i ) {
    index = i * l_maxsef;
    for(int j=0; j< l_maxsef; ++j ) {
      longArray[index + j ] = matrix[i][j];
    }
  }
  
  char op[] = "+";
  long l_max =  l_maxsef * l_maxsef;
  
  GA::igop( (int *) &longArray[0], l_max, op );

  for(int i=0; i< l_maxsef; ++i ) {
    index = i * l_maxsef;
    for(int j=0; j< l_maxsef; ++j ) {
      matrix[i][j] = longArray[index + j ];
    }
  }
}

// Collective:  Timer wrapped version
int PsociGAhamiltonian::matrixVectorProductsIncore(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int hvcount = matrixVectorProductsIncore(start, end, g_v, v_dims, g_Hv, hv_dims, v_lo, v_hi );
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return(hvcount);
}

/* This matrix vector product attempts an update in a highly localized data parallel manner. As a result
   both CIMAT and v (Hv) must have the first index  distributed in GA in the same way
   This constraint will be relaxed at a later date
   
   CIMAT[maxsef][maxsparse], V[maxsef], Hv[maxsef]
*/

//TODO I am pretty sure we can now relax the conforming maxsef distribution constraint
// Collective
// distributions of v and Hv must conform

int PsociGAhamiltonian::matrixVectorProductsIncore( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
						    int * hv_dims, int * v_lo, int * v_hi )
{
#ifdef DETAILEDHAMCHECK
  char * vectorMessage = "PsociGAhamiltonian:matrixVectorProductsIncore:Vector";
  g_v->checkHandle( vectorMessage );
  
  char * HvectorMessage = "PsociGAhamiltonian:matrixVectorProductsIncore:H-Vector";
  g_Hv->checkHandle( HvectorMessage );
  
  char * cimatMessage = "PsociGAhamiltonian:matrixVectorProductsIncore:cimat";
  l_g_cimat->checkHandle( cimatMessage );
#endif
  
  int g_rank = GA::nodeid();
  
  // Check basic array dimensions
  
  int type;
  int dims[ 2 ]; //used by all
  int idim;
  
  //  l_g_cimat->inquire( &type, &idim , dims );
  int l_maxsefs = local_dims[0];
  
  /*
    g_v->inquire( &type, &idim , dims );
    int v_sefs = dims[0];
    int maxroots = dims[1];
  */
  int v_sefs = *v_dims;
  int maxroots = *v_dims+1;
  
  if ( l_maxsefs != v_sefs ) {
    cerr << "Nonconforming cimat and vector array " << l_maxsefs << " " << v_sefs << endl;
    GA::error( "Nonconforming cimat and vector array ", -1);
  } 
  
  /*
    g_Hv->inquire( &type, &idim , dims );
    int Hv_sefs = dims[0];
    int Hvmaxroots = dims[1];
  */

  int Hv_sefs = *hv_dims;
  int Hvmaxroots = *hv_dims+1;
  
  if ( l_maxsefs != Hv_sefs ) {
    cerr << "Nonconforming cimat and Hvector array " << l_maxsefs << " " << Hv_sefs << endl;
    GA::error( "Nonconforming cimat and Hvector array ", -1);
  } 
  
  /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
     within GA.  In the future this will be relaxed.
     Check them out here
  */
  
  int lo[2], hi[2];
  int cilo[2], cihi[2];
  
//  l_g_cimat->distribution( GA::nodeid(), cilo, cihi );
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  cilo[0] = local_cilo[0];
  cilo[1] = local_cilo[1];
  cihi[0] = local_cihi[0];
  cihi[1] = local_cihi[1];
  
#if 0
  cout << "CIMAT new ISEF distribution is " << clow << " " << chi << endl;
#endif
  
//  g_v->distribution( GA::nodeid(), lo, hi );
//  int vlow = lo[0];
//  int vhi  = hi[0];

  int vlow = *v_lo;
  int vhi  = *v_hi;
/*
  lo[0] = *v_lo;
  lo[1] = *v_lo+1;
  hi[0] = *v_hi;
  hi[1] = *v_hi+1;
*/

#if 0
  cout << "VECTOR maxsef distribution is " << vlow << " " << vhi << endl;
#endif
  
  /* Can eliminate this check since PsociGAbasis already does it 
   */
  /*
    g_Hv->distribution( GA::nodeid(), lo, hi );
    int hvlow = lo[0];
    int hvhi  = hi[0];
  */
  int hvlow = vlow;
  int hvhi  = vhi;

#if 0
  cout << "H-VECTOR maxsef distribution is " << hvlow << " " << hvhi << endl;
  g_Hv->printDistribution();
#endif
  
  //TODO check if this causes performance issues
  
#ifdef DETAILEDHAMCHECK
  if ( clow != vlow || clow != hvlow || vlow != hvlow ) {
    cout << "WARNING: LOs are not equal" << endl;
    //    GA::error("LOs are wrong ",-1);
  }
  if ( chi != vhi || chi != hvhi || vhi != hvhi ) {
    cout << "WARNING: HIs are not equal " << endl;
    //GA::error("HIs are wrong ",-1);
  }
#endif
  
  // Okay we are equivelently distributed accross the cores.
  // Now we can use a data-parallel update scheme which is very fast
  // Start the processing Only process one H*v at a time for now
  
  if ( start < 0 || end <= start || end > maxroots ) {
    cerr << "absurd value for start or end " << start << " " << end << endl;
    GA::error("absurd value for start or end ", end-start );
  }
  
  // Loop over all sefs, fetching LOCAL cimat for processing.
  
  vector<double> bufV( l_maxsefs, 0);
  
  lo[0] = 0; // v and Hv related terms that basically never change
  hi[0] = l_maxsefs - 1;
  lo[1] = 0; 
  hi[1] = 0;
  
  vector<int> icicol( l_maxsparse, 0 );
  vector<double> cimat( l_maxsparse, 0.0 );
  vector<double> testVector( l_maxsefs, 0.0 );// Vector buffer for individual H*v
  
  // Start the work 
  
  int hvcount=0;
  int n;
  
  for(int ivec=start; ivec< end; ++ivec ) 
    {
#if 0
      cout << "IVEC HAM " << ivec << endl;
#endif
      lo[1]=ivec;
      hi[1]=ivec;
      n=1;
      
      vector<double> bufV( l_maxsefs, 0); // Is this faster or slower than simply zeroing out?
      
// Eventually you want to chunkup these vector reads

      if (g_rank == 0 ) g_v->get( lo, hi, &testVector[0], &n ); // Everybody needs one
      GA::brdcst( &testVector[0], testVector.size()*sizeof(double), 0);
      
      // Begin processing compute the H*v term and shove result into Hv
      
      int number; // how long was the relevant CI row ?
      int n;     // used for ld
      double dOne = 1.0;
      
      //cout << "ilo ihi are " << clow << " " << chi << endl;
      
      /*
	vector<double> rtemp( local_nbasis, 0.0 );
	vector<double> rhv( local_nbasis, 0.0 );
	int ld=1;
	int lo[2],hi[2];
	lo[0]=0;
	lo[1]=0;
	hi[0]=local_nbasis-1;
	hi[1]=0;
	g_rnorms->get( lo, hi, &rtemp[0], &ld);
	double etemp;;
	local_g_Hv->get( lo, hi, &rhv[0], &ld);
	hi[0]=0;
	g_evals->get(lo,hi,&etemp, &ld);
	hi[0]=local_nbasis-1;
      */

      for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
	cilo[0] = i;
	cihi[0] = i; 
	cilo[1] = 0;
	cihi[1] = 0;
	l_g_number->get( cilo, cihi, &number, &n );
	
	if ( number > l_maxsparse ) {
	  cerr << "number exceeds specified maxsparse: Will not be able to fetch the GA." << endl;
	  cerr << " number = " << number << " maxsparsity is " << l_maxsparse << endl;
	  GA::error("number exceeds specified maxsparse:", l_maxsparse);
	}
	
	cihi[1] = number - 1;  
	n = number;
	l_g_icicol->get( cilo, cihi, &icicol[ 0 ], &n);
	l_g_cimat->get( cilo, cihi, &cimat[ 0 ], &n);
	
// Way too much effort being expended here
	
	int jhi;
	for(int j=0; j< number; ++j ) {
	  jhi = icicol[ j ]; 
	  bufV[i] += cimat[ j ] * testVector[ jhi ]; 
	  //cout << "bufV is " << i << " " << j  << " " << jhi << " " <<  cimat[ j ] << " " << testVector[ jhi ] <<  " " << bufV[ i ] << endl;
	}
	
	for(int j=0; j< number-1; ++j ) {
	  jhi = icicol[j];
	  bufV[ jhi ] += cimat[ j ] * testVector[ i ];
          //cout << "bufV-2 is " << j << " " << jhi << " " << testVector[ i ] << " " << cimat[ j ] << " " << bufV[ jhi ] << endl;
	}
	
      }
      // Push back to ivec in H*v 
      /* At some point we may want to use deltaHs instead of Hs but for now this should be okay.
	 Using the acc method, however, forces you to specify any preexisting H*v terms using
	 PsociGAbasis::setpresent_matrixProducts( numPrecomputedRoots ); 
	 
	 Be careful about preexisting H*v will be accumulated causing errors.
	 We MUST accumulate here  
      */
      
      n=1;
      g_Hv->acc( lo, hi, &bufV[ 0 ], &n, &dOne );
      
      ++hvcount;
    }  
  return(hvcount);
}

// Alternative version of the matrixVector product 
int PsociGAhamiltonian::matrixVectorProductsIncoreGOP(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int hvcount = matrixVectorProductsIncoreGOP(start, end, g_v, v_dims, g_Hv, hv_dims, v_lo, v_hi );
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return(hvcount);
}

int PsociGAhamiltonian::matrixVectorProductsIncoreGOPSplit(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int hvcount = matrixVectorProductsIncoreGOPSplit(start, end, g_v, v_dims, g_Hv, hv_dims, v_lo, v_hi );
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return(hvcount);
}



/* This matrix vector product attempts an update in a highly localized data parallel manner. As a result
   both CIMAT and v (Hv) must have the first index  distributed in GA in the same way
   This constraint will be relaxed at a later date
   
   CIMAT[maxsef][maxsparse], V[maxsef], Hv[maxsef]
*/

int PsociGAhamiltonian::matrixVectorProductsIncoreGOP( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
						    int * hv_dims, int * v_lo, int * v_hi )
{
#ifdef DETAILEDHAMCHECK
  char * vectorMessage = "PsociGAhamiltonian:matrixVectorProductsIncore:Vector";
  g_v->checkHandle( vectorMessage );
  
  char * HvectorMessage = "PsociGAhamiltonian:matrixVectorProductsIncore:H-Vector";
  g_Hv->checkHandle( HvectorMessage );
  
  char * cimatMessage = "PsociGAhamiltonian:matrixVectorProductsIncore:cimat";
  l_g_cimat->checkHandle( cimatMessage );
#endif
  
  // Check basic array dimensions
  
  int type;
  int dims[ 2 ]; //used by all
  int idim;
  
  int l_maxsefs = local_dims[0];

  char op[] = "+";
  long numElems = l_maxsefs;
  
  int v_sefs = *v_dims;
  int maxroots = *v_dims+1;
  
#ifdef DETAILEDHAMCHECK
  if ( l_maxsefs != v_sefs ) {
    cerr << "Nonconforming cimat and vector array " << l_maxsefs << " " << v_sefs << endl;
    GA::error( "Nonconforming cimat and vector array ", -1);
  } 
#endif

  int Hv_sefs = *hv_dims;
  int Hvmaxroots = *hv_dims+1;
  
#ifdef DETAILEDHAMCHECK
  if ( l_maxsefs != Hv_sefs ) {
    cerr << "Nonconforming cimat and Hvector array " << l_maxsefs << " " << Hv_sefs << endl;
    GA::error( "Nonconforming cimat and Hvector array ", -1);
  } 
#endif
  
  /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
     within GA.  In the future this will be relaxed.
     Check them out here
  */
  
  int lo[2], hi[2];
  int cilo[2], cihi[2];
  
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  cilo[0] = local_cilo[0];
  cilo[1] = local_cilo[1];
  cihi[0] = local_cihi[0];
  cihi[1] = local_cihi[1];
  
  int vlow = *v_lo;
  int vhi  = *v_hi;
  
  int hvlow = vlow;
  int hvhi  = vhi;

  double tempTime= psociTime();
  double tempTime2, getTime = 0.0, accTime=0.0, brdTime=0.0;
  
#ifdef DETAILEDHAMCHECK
  if ( clow != vlow || clow != hvlow || vlow != hvlow ) {
    cout << "WARNING: LOs are not equal" << endl;
    //    GA::error("LOs are wrong ",-1);
  }
  if ( chi != vhi || chi != hvhi || vhi != hvhi ) {
    cout << "WARNING: HIs are not equal " << endl;
    //GA::error("HIs are wrong ",-1);
  }
#endif
  
  // Okay we are equivelently distributed accross the cores.
  // Now we can use a data-parallel update scheme which is very fast
  // Start the processing Only process one H*v at a time for now
  
#ifdef DETAILEDHAMCHECK
  if ( start < 0 || end <= start || end > maxroots ) {
    cerr << "absurd value for start or end " << start << " " << end << endl;
    GA::error("absurd value for start or end ", end-start );
  }
#endif
  
  // Loop over all sefs, fetching LOCAL cimat for processing.
  
  lo[0] = v_lo[0];
  hi[0] = v_hi[0];
  lo[1] = 0; 
  hi[1] = 0;
  
  int index = v_lo[0]; // For bufV and testvector
  int width = v_hi[0] - v_lo[0] + 1; // width of local vector read/put

  vector<int> icicol( l_maxsparse, 0 );
  vector<double> cimat( l_maxsparse, 0.0 );
  
  int hvcount=0;
  int n;
  double vecti;
  int number; // how long was the relevant CI row ?
  double dOne = 1.0;

// TOO BIG
  vector<double> bufV( l_maxsefs); // Is this faster or slower than simply zeroing out?
  vector<double> testVector( l_maxsefs );// Vector buffer for individual H*v

  //cout << " INDEX " << lo[0]<<" "<<lo[1]<<" "<<hi[0]<<" "<<hi[1]<<endl;

  for(int ivec=start; ivec< end; ++ivec ) 
    {
      lo[1]=ivec;
      hi[1]=ivec;
      n=1;
      
      //tempTime2 = psociTime();
      zeroVector( bufV );
      zeroVector( testVector );

//      vector<double> bufV( l_maxsefs,0.0); // Is this faster or slower than simply zeroing out?
//      vector<double> testVector( l_maxsefs,0.0 );// Vector buffer for individual H*v

      g_v->get( lo, hi, &testVector[index], &n ); // Everybody needs one
      GA::gop( &testVector[0], numElems, op );
      //brdTime += psociTime() - tempTime2;
      
      // Begin processing compute the H*v term and shove result into Hv
      
      int jhi;
      for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
	cilo[0] = i;
	cihi[0] = i; 
	cilo[1] = 0;
	cihi[1] = 0;

        //tempTime2 = psociTime();
	l_g_number->get( cilo, cihi, &number, &n );
	
#ifdef DETAILEDHAMCHECK
	if ( number > l_maxsparse ) {
	  cerr << "number exceeds specified maxsparse: Will not be able to fetch the GA." << endl;
	  cerr << " number = " << number << " maxsparsity is " << l_maxsparse << endl;
	  GA::error("number exceeds specified maxsparse:", l_maxsparse);
	}
#endif
	
        vecti = testVector[ i ];

	cihi[1] = number - 1;  
	n = number;
	l_g_icicol->get( cilo, cihi, &icicol[ 0 ], &n);
	l_g_cimat->get( cilo, cihi, &cimat[ 0 ], &n);
        //getTime += psociTime() - tempTime2;
	
        int jhi;
        double testVec = testVector[ i ];
        for(int j=0; j< number; ++j ) {
          jhi = icicol[ j ];
          bufV[i] += cimat[ j ] * testVector[ jhi ];
        //  cout << "Loop One " << bufV[ i ]<<" "<<j<<" "<<cimat[ j ]<<" "<<testVector[ jhi ]<<endl;
        }

        for(int j=0; j< number-1; ++j ) {
          jhi = icicol[j];
          bufV[ jhi ] += cimat[ j ] * testVec;
         // cout << "Loop Two " << bufV[ jhi ]<<" "<<j<<" "<<cimat[ j ]<<" "<<testVec<<endl;
        }
      }

      // Push back to ivec in H*v 
      /* At some point we may want to use deltaHs instead of Hs but for now this should be okay.
	 Using the acc method, however, forces you to specify any preexisting H*v terms using
	 PsociGAbasis::setpresent_matrixProducts( numPrecomputedRoots ); 
	 
	 Be careful about preexisting H*v will be accumulated causing errors.
	 We MUST accumulate here  
      */
      n=1;

      //tempTime2 = psociTime();
      GA::gop( &bufV[0], numElems, op);
      g_Hv->put( lo, hi, &bufV[ index ], &n);
      //accTime += psociTime() - tempTime2;
      
      ++hvcount;

      //double outerTime = psociTime() - tempTime;
      //if (GA::nodeid() == 0 ) cout << "GOP-OUTER No split" << outerTime<<" ACCUM "<< accTime<<" GET " << getTime<< " BRD-now-GOP "<<brdTime<<endl;
    }  
  return(hvcount);
}

int PsociGAhamiltonian::matrixVectorProductsIncoreGatherScatter(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int hvcount = matrixVectorProductsIncoreGatherScatter(start, end, g_v, v_dims, g_Hv, hv_dims, v_lo, v_hi );
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return(hvcount);
}

/* NON SPLIT versionof the GATHER SCATTER methods
   compile code as -DGATHERHV
*/
int PsociGAhamiltonian::matrixVectorProductsIncoreGatherScatter( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
								 int * hv_dims, int * v_lo, int * v_hi )
{
#ifndef SUPPLEMENTAL
  GA::error(" Must compile code with -DSUPPLEMENTAL ",-1);
#endif
#ifdef SQUAREHAMILTONIAN
  GA::error(" Must NOT compile with -DSQUAREHAMILTONIAN",1);
#endif

  // Check basic array dimensions
  
  int type;
  int dims[ 2 ]; //used by all
  int idim;
  
  int l_maxsefs = local_dims[0];
  
  char op[] = "+";
  long numElems = l_maxsefs;
  
  /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
     within GA.  In the future this will be relaxed.
     Check them out here
  */
  
  /* cimat and friends are distributed over isefs but no breaking along sparity axis */
  
  
  int lo[2], hi[2];
  int cilo[2], cihi[2];
  
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  cilo[0] = clow;
  cilo[1] = local_cilo[1];
  cihi[0] = chi;
  cihi[1] = local_cihi[1];
  
  // these are the local vector piece
  
  double temp;
  double tempTime= psociTime();
  double tempTime2, getTime = 0.0, accTime=0.0, brdTime=0.0;
  double temp2, bufTime=0.0;
  double mapTime=0.0;
  double pairTime=0.0;
  double minimapTime=0.0;
  
  double vecTime=0.0;
  double vecSort=0.0;
  
  // Loop over all sefs, fetching LOCAL cimat for processing.
  // Gather Scatter objects can never be bigger than this
  
 // vector<double> bufV( max(l_maxsparse,l_maxsparse_supp), 0.0  );
  
  const int nchunk = l_nchunk;
  int chunk_width;
  int ch_start=0;
  int ch_end=0;
  (l_maxsefs % nchunk==0)? chunk_width = l_maxsefs / nchunk: chunk_width = 1 + (l_maxsefs / nchunk);
  if ( GA::nodeid() == 0 ) cout << "VECTOR get CHUNK is " << nchunk << " " << chunk_width << " "<<MAX_AGGREGATE<<" "<<MAX_DEFLATE_ELEMENTS<<" "<<MINIMAP_MAX<<endl;
  
  vector<double> testVector( chunk_width , 0.0 );
  
  int hvcount=0;
  int n=1;
  int number; // how long was the relevant CI row ?
  int fullnumber=0;
  double dOne = 1.0;
  double alpha = 1.0;
  double  listTime=0.0,compTime1=0.0, compTime2=0.0, scatTime=0.0, getVTime=0.0, deflateTime=0.0;
  double deflateListTime=0.0, deflateVecTime=0.0;
  
  double temptime, cimatTime=0.0;
  int chunk_lo[2];
  int chunk_hi[2];
  
  int newnumber=0;
  
  int numdeflates;
  
  int width = MAX_AGGREGATE; 
 
  //if ( my_chunklist.size() != 1 ) GA::error("chunklist must be 1 " ,1);
  
  int deflate_index;
  int partialnumber;
  int nummapcontracts=0; // number of time minimap was copied to buffer_map
  
  double scatbuf;

// multimaps MUST NOT SPAN multiple ivecs
  for(int ivec=start; ivec< end; ++ivec ) 
    {
      multimap<int, double > buffer_map;
      multimap<int, double > mini_map;
      
      // current ly only works for 1 core
      chunk_lo[1] = ivec;
      chunk_hi[1] = ivec;
      
      //pretend loop for now
      
      numdeflates=0;
      vector<int>::iterator it;
      
      for(it=my_chunklist.begin(); it!=my_chunklist.end(); ++it) { // this potentially staggers accross cores 
	
	ch_start = chunk_width*(*it); 
	ch_end = min( ch_start + chunk_width, l_maxsefs );
	
	//    ch_start=0;
	//    ch_end = l_maxsefs -1 ;
	chunk_lo[0] = ch_start;
	chunk_hi[0] = ch_end-1;
	
	// Step One: Grab the current GA-vector chunk for processing
	
	temp = psociTime();
#ifdef STAGGERCHUNKLIST
	g_v->get( chunk_lo, chunk_hi, &testVector[0], &n );
#else
	if ( GA::nodeid() == 0 ) g_v->get( chunk_lo, chunk_hi, &testVector[0], &n );
	GA::brdcst( &testVector[0], testVector.size()*sizeof(double), 0);
#endif
	getVTime += psociTime() - temp;

// alternative

//      tempTime2=psociTime();
//      zeroVector( testVector );

//      g_v->get( chunk_lo, chunk_hi, &testVector[0], &n ); // we are only fetching the local piece into a full size space
//      GA::gop( &testVector[0], chunkwidth, op );
//      getVTime += psociTime() - tempTime2;

// end alternative
	
	// Hamiltonian reduction
	
	int jhi;
	
/* Changed the meaning to the maximum number of elements before a scatter
*/
	deflate_index = 0; // how many deflations/scatters occured
	
	/* 
	   how many to minimap insertions ( per ivec )
	   Resets to zero at minimap >= MINIMAP_MAX. THis minimizes time
	   in a hot loop a sorted insertion ( nlogn ). Independant 
	   on the value of i
	*/
	
	
	int buf_index=-1;

// Treat the map elements counter.
        int max_deflate_elem=0;
        int minimap = 0;
	
#ifdef MATRIXBALANCE
	int fetch_cimat_count=-1;
	
	// for(int it=0; it < cimat_list.size(); ++it ) {
	for(int it=0; it < cimat_list.size(); it+= MAX_AGGREGATE ){ // a factor of 5 in time variation
	  int mxnum = cimat_list.size();
	  int max_num = MAX_AGGREGATE;
	  if ( ( mxnum - it ) < MAX_AGGREGATE ) max_num = mxnum - it;
	  
	  vector<int> newlist;
          temptime=psociTime();
	  for(int j=0; j< max_num; ++j) {
	    newlist.push_back( cimat_list[ it+j ]);
	  }
          fetchsubListHamiltonianData( newlist ); // cimat still populate starting at 0 
          buf_index=-1;
	  minimap = 0; 
          cimatTime += psociTime() - temptime;
	  
	  for(int it2=0; it2< max_num; ++it2) { // extra loop
	    int i = newlist[ it2 ]; 
#else
	    for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
#endif
	      ++buf_index;
	      
#ifndef MATRIXBALANCE
              // local data is grabbed in this case
	      buf_index=0; //override and simply grab local i
	      fetchLocalHblock( i, local_cimat[buf_index], local_icicol[buf_index], local_number[buf_index] ); // It actually does not need to be local but beware performance
#endif
	      
	      int start_rng=local_icicol[buf_index][0];
	      int end_rng=local_icicol[buf_index][local_number[buf_index]-1];
	      number = local_number[buf_index];
	      zeroVector( bufV, number );
	      
	      double dtemp=0.0;
	      
	      if (!ch_end-1 < start_rng  || !ch_start > end_rng ) { //none within the range
		
               // ++minimap; // grabs full size later
		temp = psociTime();
		
		for(int j=0; j< number; ++j ) {
		 // temp2=psociTime();
		  jhi = local_icicol[buf_index][j];
		  if ( jhi >= ch_start && jhi < ch_end ) {
		    dtemp += local_cimat[buf_index][ j ] * testVector[ jhi - ch_start]; //globally accesses testVector
		  }
		}
                
		mini_map.insert( pair<int, double>( local_icicol[buf_index][number-1], dtemp ) );
	//	minimapTime+= psociTime() - temp2;
	      }
	      
	      // Need to minimize the length of the multimap otherwise NlogN bites you
	      
//	      compTime1 += psociTime() - temp;
//	      temp = psociTime();
	      
	      if ( i >= ch_start && i < ch_end ) {
		double testVec = testVector[ i - ch_start ]; 
		
		for(int j=0; j< number-1; ++j ) {
		  int jhi = local_icicol[buf_index][j];
//		  temp2=psociTime();
		  bufV[ j ] += local_cimat[buf_index][ j ] * testVec;
//		  bufTime += psociTime()-temp2;
		  
		  if ( abs(bufV[j]) > MIN_HVEC_TOKEEP ) { 
//        	    temp2=psociTime();
		    mini_map.insert(pair<int, double>( jhi, bufV[ j ] ) );
//		    minimapTime += psociTime()-temp2;
		  }
		}
	      } // ch_start screen
              minimap = mini_map.size();
	      
              compTime2 += psociTime() - temp;

              temp2=psociTime();
	      if ( minimap >= MINIMAP_MAX ) { //other cases get handled below
		++nummapcontracts;
//		temp = psociTime();
		buffer_map.insert( mini_map.begin(), mini_map.end() );
                max_deflate_elem += mini_map.size();
              //  dump_count += mini_map.size();
		mini_map.clear();
//		mapTime+= psociTime() - temp2;
		minimap = 0;
	      }
              mapTime+= psociTime() - temp2;
	      
//	      compTime2 += psociTime() - temp;
	      
#ifdef MATRIXBALANCE
	      if ( max_deflate_elem >= MAX_DEFLATE_ELEMENTS || i >= cimat_list[ cimat_list.size() - 1 ]  ) {
#else
		if ( max_deflate_elem >= MAX_DEFLATE_ELEMENTS || i >= chi ) {
#endif
		  temp = psociTime();
		  if ( mini_map.size() >= 1 ) {
		    buffer_map.insert( mini_map.begin(), mini_map.end() ); // ensure no strays
		    ++nummapcontracts;
                    minimap=0;
		    mini_map.clear();
		  }
		  mapTime+= psociTime() - temp;
		  
		  temp=psociTime();
		  newnumber = deflateMultimapToVector( buffer_map, dataarray, indexarray, ivec ); 
                  buffer_map.clear();
		  deflateTime +=psociTime()-temp;

		  if ( newnumber > 0 ) {
		    temp=psociTime();
		    NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&dataarray[0], (int*)&indexarray[0], newnumber, &alpha);
		    deflate_index=0;
                    scatbuf = psociTime() - temp;
                    scatTime += scatbuf;
                    //if (GA::nodeid() == 0 ) cout << "do scattr with " << newnumber << " " << scatbuf << endl;
		    //scatTime += psociTime() - temp;
		    ++numdeflates;
		  } 
		  //dump_count = 0;
                  max_deflate_elem=0;
                  
		  
		} // dump count
#ifdef MATRIXBALANCE
	      } // extra loop: MAX_AGGREGATE blocks	
#endif
	    } // clow
	  } // chunk_list 
	  ++hvcount;
	}  // ivec 
	
//	cout << GA::nodeid() << "AGGREGATE ciTime getV comp1 comp2 bufV mapT minimapT deflate SUM  mapcontracts numdeflates "<< MAX_AGGREGATE <<" "<< cimatTime<<" "<<getVTime<<" "<<compTime1 <<" "<<compTime2 << " "<<bufTime<<" " <<mapTime <<" "<< minimapTime << " " <<deflateTime<<" "<<compTime1+compTime2+cimatTime+getVTime+mapTime+minimapTime+pairTime+deflateTime<<" "<<nummapcontracts<<" " << numdeflates<< endl;

       cout << GA::nodeid() << "cimatTime getV comp2 bufV mapT minimapT deflate scattTime SUM scatters"<<cimatTime<<" "<<  getVTime<<" "<<compTime2 << " "<<mapTime <<" "<< " " <<deflateTime<<" "<<scatTime<<" "<<compTime2+cimatTime+getVTime+mapTime+deflateTime+scatTime<<" "<<numdeflates<<endl;

	
	GA::SERVICES.sync();
	return(hvcount);
      }

/* Distributed input vector but replicated bufV

   The current DEFAULT for non-square calculations
*/
int PsociGAhamiltonian::matrixVectorProductsSplitChunkRepAccumulate(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int hvcount = matrixVectorProductsSplitChunkRepAccumulate(start, end, g_v, v_dims, g_Hv, hv_dims, v_lo, v_hi );
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return(hvcount);
}
 
int PsociGAhamiltonian::matrixVectorProductsSplitChunkScatterAccumulate(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int hvcount = matrixVectorProductsSplitChunkScatterAccumulate(start, end, g_v, v_dims, g_Hv, hv_dims, v_lo, v_hi );
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return(hvcount);
}

/* return approximate number of bytes per-core required for the matrixVectorProductsSplitChunkRepAccumulate method
   TODO this has become incorrect over time
*/
 long PsociGAhamiltonian::matrixVectorProductsMemoryPerCore()
{
   int l_maxsefs = local_dims[0];
   long word=sizeof( double );
   long intword=sizeof( int );

// Note these values are somewhat accurate for JCHUNK method; not the scatter methods
#ifdef NEWACCUMULATEDSCATTER
   long sum = MAX_AGGREGATE*intword*max(l_maxsparse,l_maxsparse_supp)+MAX_AGGREGATE*word*max(l_maxsparse,l_maxsparse_supp)
              + word*(2*(l_maxsefs / l_nchunk));
#else
   long sum = MAX_AGGREGATE*intword*max(l_maxsparse,l_maxsparse_supp)+2*MAX_AGGREGATE*word*max(l_maxsparse,l_maxsparse_supp)
              + word*(2*(l_maxsefs / l_nchunk));
#endif
   return( sum );
}

int PsociGAhamiltonian::matrixVectorProductsSplitChunkRepAccumulate( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
								      int * hv_dims, int * v_lo, int * v_hi )
 {
#ifndef SUPPLEMENTAL
   GA::error(" Must compile code with -DSUPPLEMENTAL ",-1);
#endif
#ifdef SQUAREHAMILTONIAN
   GA::error(" Must NOT compile with -DSQUAREHAMILTONIAN",1);
#endif
   
   // Check basic array dimensions
   
   int type;
   int dims[ 2 ]; //used by all
   int idim;
   int l_maxsefs = local_dims[0];
   
   char op[] = "+";
   
   /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
      within GA.  In the future this will be relaxed.
      Check them out here
   */
   
   /* cimat and friends are distributed over isefs but no breaking along sparity axis */
  
  
  int lo[2], hi[2];
  
  int cilo[2], cihi[2];
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  
  cilo[0] = clow;
  cilo[1] = local_cilo[1];
  cihi[0] = chi;
  cihi[1] = local_cihi[1];
  
  // these are the local vector piece
  
  double temp;
  double tempTime= psociTime();
  double tempTime2, getTime = 0.0, accTime=0.0, brdTime=0.0;
  double temp2;
  
  
  // Loop over all sefs, fetching LOCAL cimat for processing.
  // Gather Scatter objects can never be bigger than this
  const int nchunk = l_nchunk;
  
  int chunk_width;
  int ch_start=0;
  int ch_end=0;
  (l_maxsefs % nchunk==0)? chunk_width = l_maxsefs / nchunk: chunk_width = 1 + (l_maxsefs / nchunk);
  if ( GA::nodeid() == 0 ) cout << "exp VECTOR get CHUNK is " << nchunk << " " << chunk_width << " "<<MAX_AGGREGATE<<endl;
#ifdef USEMPI 
  if ( GA::nodeid() == 0 ) cout << "MPI collectives for H*v " << endl;
#endif
  
  vector<double> testVector( chunk_width , 0.0 );
  
  int hvcount=0;
  int n=1;
  int number; // how long was the relevant CI row ?
  double dOne = 1.0;
  double alpha = 1.0;
  double  listTime=0.0,compTime1=0.0, compTime2=0.0, scatTime=0.0, getVTime=0.0;
  
  double temptime, cimatTime=0.0;
  int chunk_lo[2];
  int chunk_hi[2];
 
  int width = MAX_AGGREGATE;  // pertain to H prefeting
  
  //if ( my_chunklist.size() != 1 ) GA::error("chunklist must be 1 " ,1);
  
  vector<double> bufV(max(l_maxsefs,l_maxwidth) ); // Is this faster or slower than simply zeroing out?
  
  /* For bufV upload after GOP*/
  lo[0] = v_lo[0];
  hi[0] = v_hi[0];
  lo[1] = 0;
  hi[1] = 0;
  int index = v_lo[0]; // For bufV and testvector

  for(int ivec=start; ivec< end; ++ivec ) 
    {
      chunk_lo[1] = ivec;
      chunk_hi[1] = ivec;
      
      zeroVector( bufV );
      
      vector<int>::iterator it;
      for(it=my_chunklist.begin(); it!=my_chunklist.end(); ++it) { // this potentially staggers accross cores 
	
        ch_start = chunk_width*(*it);
        ch_end = min( ch_start + chunk_width, l_maxsefs );

        //    ch_start=0;
        //    ch_end = l_maxsefs -1 ;
        chunk_lo[0] = ch_start;
        chunk_hi[0] = ch_end-1;

        // Step One: Grab the current GA-vector chunk for processing

	temp = psociTime();
        fetchAndReplicateVectorChunk((*it), ivec, chunk_width, v_lo, v_hi, testVector, g_v );
	getVTime += psociTime() - temp;
	
	// Hamiltonian reduction
	
	int jhi;
	int buf_index=-1;
	
#ifdef MATRIXBALANCE
	int fetch_cimat_count=-1;
	
	// for(int it=0; it < cimat_list.size(); ++it ) {
	
	for(int it=0; it < cimat_list.size(); it+= MAX_AGGREGATE ){ // a factor of 5 in time variation
	  int mxnum = cimat_list.size();
	  int max_num = MAX_AGGREGATE;
	  if ( ( mxnum - it ) < MAX_AGGREGATE ) max_num = mxnum - it;
	  
	  vector<int> newlist;
          temptime=psociTime();
	  for(int j=0; j< max_num; ++j) {
	    newlist.push_back( cimat_list[ it+j ]);
	  }
          fetchsubListHamiltonianData( newlist ); // cimat still populate starting at 0 
          buf_index=-1;
          cimatTime += psociTime() - temptime;
	  
	  for(int it2=0; it2< max_num; ++it2) { // extra loop
	    int i = newlist[ it2 ]; 
#else
	    for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
#endif
	      ++buf_index;
	      
#ifndef MATRIXBALANCE
              // local data is grabbed in this case

// WE see significant load imbal;ance probably caused by this

	      buf_index=0; //override and simply grab local i
	      fetchLocalHblock( i, local_cimat[buf_index], local_icicol[buf_index], local_number[buf_index] ); // It actually does not need to be local but beware performance
#endif
	      int start_rng=local_icicol[buf_index][0];
	      int end_rng=local_icicol[buf_index][local_number[buf_index]-1];
	      number = local_number[buf_index];
	      
	      double dtemp=0.0;
	      if (!ch_end-1 < start_rng  || !ch_start > end_rng ) { //none within the range
		temp = psociTime();
		for(int j=0; j< number; ++j ) {
		  jhi = local_icicol[buf_index][j];
		  if ( jhi >= ch_start && jhi < ch_end ) {
		    dtemp += local_cimat[buf_index][ j ] * testVector[ jhi - ch_start]; //globally accesses testVector
		  }
		}
                bufV[i] += dtemp; // not 100% we should be summing here
		compTime1+=psociTime() - temp;
	      }
	      if ( i >= ch_start && i < ch_end ) {
                temp = psociTime();
		double testVec = testVector[ i - ch_start ]; 
		for(int j=0; j< number-1; ++j ) {
		  int jhi = local_icicol[buf_index][j];
		  bufV[ jhi ] += local_cimat[buf_index][ j ] * testVec;
		}
		compTime2+=psociTime() - temp;
	      } // ch_start screen
#ifdef MATRIXBALANCE
	    } // extra loop: MAX_AGGREGATE blocks	
#endif
	  } // clow
	} // chunk_list 

	tempTime2 = psociTime();

/*
cout << v_lo[0]<<" "<<v_lo[1]<<" "<<v_hi[0]<<" "<<v_hi[1]<<" " << ivec<<" "<<l_maxsefs<<endl;
        replicateAndPutHV( ivec, l_maxsefs, v_lo, v_hi, bufV, g_Hv );
*/

	n=1;
	lo[1] = ivec;
	hi[1] = ivec;
	GA::gop( &bufV[0], l_maxsefs, op);
	g_Hv->put( lo, hi, &bufV[index], &n);
	accTime += psociTime() - tempTime2;
	
	++hvcount;
      }  // ivec 
      
      GA::SERVICES.sync();
      if ( GA::nodeid() == 0 ) cout << GA::nodeid() << " REP: cimatTime getV accTime compTime1 comptTime2 "<<cimatTime<<" "<<  getVTime<<" "<<accTime << " "<<compTime1<<" "<<compTime2<<endl;
      
      return(hvcount);
    }

int PsociGAhamiltonian::matrixVectorProductsSplitChunkScatterAccumulate( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
								      int * hv_dims, int * v_lo, int * v_hi )
 {
#ifndef SUPPLEMENTAL
   GA::error(" Must compile code with -DSUPPLEMENTAL ",-1);
#endif
#ifdef SQUAREHAMILTONIAN
   GA::error(" Must NOT compile with -DSQUAREHAMILTONIAN",1);
#endif
   
   // Check basic array dimensions
   
   int type;
   int dims[ 2 ]; //used by all
   int idim;
   int l_maxsefs = local_dims[0];
   
   char op[] = "+";
   
   /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
      within GA.  In the future this will be relaxed.
      Check them out here
   */
   /* cimat and friends are distributed over isefs but no breaking along sparity axis */
  
  int lo[2], hi[2];
  int cilo[2], cihi[2];
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  
  cilo[0] = clow;
  cilo[1] = local_cilo[1];
  cihi[0] = chi;
  cihi[1] = local_cihi[1];
  
  // these are the local vector piece
  
  double temp;
  double tempTime= psociTime();
  double tempTime2, getTime = 0.0, accTime=0.0, brdTime=0.0;
  double temp2;
  
  // Loop over all sefs, fetching LOCAL cimat for processing.
  // Gather Scatter objects can never be bigger than this
  const int nchunk = l_nchunk;
  
  int chunk_width;
  int ch_start=0;
  int ch_end=0;
  (l_maxsefs % nchunk==0)? chunk_width = l_maxsefs / nchunk: chunk_width = 1 + (l_maxsefs / nchunk);
  if ( GA::nodeid() == 0 ) cout << "exp VECTOR get CHUNK is " << nchunk << " " << chunk_width << " "<<MAX_AGGREGATE<<endl;
#ifdef USEMPI 
  if ( GA::nodeid() == 0 ) cout << "MPI collectives for H*v " << endl;
#endif
  
  double tempfull = psociTime(); 
  vector<double> testVector( chunk_width , 0.0 );
  
  int hvcount=0;
  int n=1;
  int number; // how long was the relevant CI row ?
  double dOne = 1.0;
  double alpha = 1.0;
  double packTime=0.0, scatbuf=0.0, deflateTime=0.0,  mapTime=0.0, listTime=0.0,compTime1=0.0, compTime2=0.0, scatTime=0.0, getVTime=0.0;
  
  double temptime, cimatTime=0.0;
  int chunk_lo[2];
  int chunk_hi[2];
 
  int width = MAX_AGGREGATE;  // pertain to H prefeting
  
  /* For bufV upload after GOP*/
  lo[0] = v_lo[0];
  hi[0] = v_hi[0];
  lo[1] = 0;
  hi[1] = 0;
  int index = v_lo[0]; // For bufV and testvector

// New scatter based stuff

  int numdeflates;
  int deflate_index;
  int newnumber;
  int nummapcontracts=0; // number of time minimap was copied to buffer_map

  double zeroTime = 0.0;

  int numscatters = 0;

  int zerowidth = chi-clow+1;
  zeroMatrix( bufV2, zerowidth  );
  //zeroMatrix( bufVIndex );


  
  for(int ivec=start; ivec< end; ++ivec ) 
    {
      chunk_lo[1] = ivec;
      chunk_hi[1] = ivec;
/*
  vector<vector<double> > bufV2(chi-clow+1,vector<double>( max(l_maxsparse,l_maxsparse_supp), 0.0  ) ); // catch a local I's worth of data then scatter
  vector<vector<int> > bufVIndex(chi-clow+1, vector<int>(  max(l_maxsparse,l_maxsparse_supp), 0  ) );
  vector<int> bufVnumber(chi-clow+1,0 );
*/

// tes tthis

#ifdef MATRIXBALANCE
//cout << "MATBal " << chi-clow+1 << endl; 
      int zerowidth = chi-clow+1;;
      zeroMatrix( bufV2, zerowidth  );
      //zeroMatrix( bufVIndex );
#else
//cout << "NO MATBAL " << endl;
      int zerowidth = chi-clow+1;
//      zeroMatrix( bufV2, zerowidth  );
      //zeroMatrix( bufVIndex );
#endif

int indexnewarray = 0;


#ifdef MATRIXBALANCE
//prefetch CIMAT data based on the some kind of chunking

// Moving this inside of the chunk loop increases the chance of better load balancing

       int newistart = -1;

       int max_num;
       int cimat_size = 0;
       vector<int> newlist;
       int indexsum = 0; // new

        int fetch_cimat_count=-1;

        for(int itci=0; itci < cimat_list.size(); itci+= MAX_AGGREGATE ){ // a factor of 5 in time variation
          int mxnum = cimat_list.size();
          max_num = MAX_AGGREGATE;
          if ( ( mxnum - itci ) < MAX_AGGREGATE ) max_num = mxnum - itci;

          temptime=psociTime();
          newlist.clear();
          for(int j=0; j< max_num; ++j) {
            newlist.push_back( cimat_list[ itci+j ]);
          }
          fetchsubListHamiltonianData( newlist ); // cimat still populate starting at 0 

          cimat_size += newlist.size(); // sums for the start of the clo iterator
          cimatTime += psociTime() - temptime;

// zero Out array only as needed.
// at this point we've prefetch as many as MAX_AGGREGATE H rows into cimat,icicol,number

#endif

      vector<int>::iterator it;
      for(it=my_chunklist.begin(); it!=my_chunklist.end(); ++it) { // this potentially staggers accross cores 
	
        temp = psociTime();
	zeroTime += psociTime() - temp;

        ch_start = chunk_width*(*it);
        ch_end = min( ch_start + chunk_width, l_maxsefs );

        //    ch_start=0;
        //    ch_end = l_maxsefs -1 ;
        chunk_lo[0] = ch_start;
        chunk_hi[0] = ch_end-1;

        // Step One: Grab the current GA-vector chunk for processing

	temp = psociTime();
        fetchAndReplicateVectorChunk((*it), ivec, chunk_width, v_lo, v_hi, testVector, g_v );
	getVTime += psociTime() - temp;
	
// Hamiltonian reduction
	
	int jhi;
	int buf_index=-1; // always restart
        // indexnewarray = -1

#ifdef MATRIXBALANCE
        indexnewarray = indexsum-1;
#else
        indexnewarray = -1; // because we always do the entire local list
#endif
          vector<int> tempindex;
#ifdef MATRIXBALANCE
	  for(int it2=0; it2< max_num; ++it2) { // extra loop
	    int i = newlist[ it2 ]; 
#else
	  for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
#endif
	      ++buf_index;
              ++indexnewarray;
#ifndef MATRIXBALANCE
              // local data is grabbed in this case
              // if not here then we already prefetched at the head of the method

	      buf_index=0; //override and simply grab local i
	      fetchLocalHblock( i, local_cimat[buf_index], local_icicol[buf_index], local_number[buf_index] ); // It actually does not need to be local but beware performance
#endif
	      int start_rng=local_icicol[buf_index][0];
	      int end_rng=local_icicol[buf_index][local_number[buf_index]-1];
	      number = local_number[buf_index];
	      
	      double dtemp=0.0;
              bufVnumber[ indexnewarray ] = number; // how many were found
	      if (!ch_end-1 < start_rng  || !ch_start > end_rng ) { //none within the range
		temp = psociTime();
		for(int j=0; j< number; ++j ) {
		  jhi = local_icicol[buf_index][j];
		  if ( jhi >= ch_start && jhi < ch_end ) {
		    dtemp += local_cimat[buf_index][ j ] * testVector[ jhi - ch_start]; 
		  }
		}
                bufV2[ indexnewarray ][ number-1 ] += dtemp;
                bufVIndex[ indexnewarray ][ number-1 ] = local_icicol[buf_index][number-1]; // no sum it is always the same
		compTime1+=psociTime() - temp;
	      }

	      if ( i >= ch_start && i < ch_end ) {
                temp = psociTime();
		double testVec = testVector[ i - ch_start ]; 
		for(int j=0; j< number-1; ++j ) {
	     	    int jhi = local_icicol[buf_index][j];
		    bufV2[ indexnewarray ][ j ] += local_cimat[buf_index][ j ] * testVec;
                    bufVIndex[ indexnewarray ][ j ] = jhi; // no need to sum
		}
		compTime2+=psociTime() - temp;
	      } // ch_start screen

                if ( indexnewarray >= MAX_AGGREGATE-1 ) {
                 ++numscatters;
                temp = psociTime();
                int realnumber = packDataNoZeros( ivec, bufVnumber, bufVIndex, bufV2, indexarray, dataarray );
                //cout << "EXPERIMANTEL test 2 realnumbe ris " << realnumber << endl;
                packTime = psociTime() - temp;
                temp = psociTime();
                if ( realnumber > 0 ) {
                     NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&dataarray[0], (int*)&indexarray[0], realnumber, &alpha);
                ++numscatters;
                }
                 scatTime += psociTime() - temp;
                 indexnewarray = -1;
             }
          } // clow

// We could conceivably do the scatter here as well

#ifdef INNERSCATTER
      temp = psociTime();
      int realnumber = packDataNoZeros( ivec, bufVnumber, bufVIndex, bufV2, indexarray, dataarray );
      //cout << GA::nodeid() << " test 2 realnumbe ris " << realnumber << endl;
      packTime = psociTime() - temp;
      temp = psociTime();
      if ( realnumber > 0 ) {
        NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&dataarray[0], (int*)&indexarray[0], realnumber, &alpha);
      ++numscatters;
      }
      scatTime += psociTime() - temp;
#endif

        } // chunk_list 

#ifdef MATRIXBALANCE
         indexsum += newlist.size(); // set to beginniing of next chunk
         }  
#endif 
          indexnewarray = -1;

/* Now we want to take the entire set of results, pack them into dataarray an dindex array 
   and scatter back to the nodes
   for a GIVEN ivec
*/
       
/*
#ifdef MATRIXBALANCE
cout << "AT THE DUMP " << endl;
          //for(int it2=0; it2< max_num; ++it2) { // extra loop
            for(int it2=0; it2 < cimat_list.size(); ++it2 ) {
            int i = cimat_list[ it2];
          
           // int i = newlist[ it2 ];
#else
            for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
#endif
cout << "dump next clow " << i << " max " << max_num << endl;
              double dtest = 0.0;
              ++indexnewarray;
                  int number = bufVnumber[ indexnewarray ] ;
                  for(int j=0; j< number; ++j ) {
                  double  dvalue = bufV2[ indexnewarray ][ j ];
                  dataarray[ j ] = dvalue;
                  dtest += abs(dvalue);
                  int offset = 2*j;
                  indexarray[ offset ] = bufVIndex[ indexnewarray ][ j ];
                  indexarray[ offset+1 ] = ivec;
 cout << "dump dataarray" << dataarray[ j ] <<" "<<indexarray[ offset ]<<" "<<indexarray[ offset+1 ]<<endl;
                  }
                  if ( dtest > 0.0 ) {
               //      ++doneloops;
                     NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&dataarray[0], (int*)&indexarray[0], number, &alpha);
                  }    
cout << "end "<<endl << endl;
             }
*/

#ifndef INNERSCATTER
      temp = psociTime();
      int realnumber = packDataNoZeros( ivec, bufVnumber, bufVIndex, bufV2, indexarray, dataarray );
     // cout << "OUTER realnumbe ris " << realnumber << endl;
      packTime = psociTime() - temp;
      temp = psociTime();
      if ( realnumber > 0 ) {
         NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&dataarray[0], (int*)&indexarray[0], realnumber, &alpha);
      }
      scatTime += psociTime() - temp;
#endif
      ++hvcount;

      }  // ivec 

      GA::SERVICES.sync();
      tempTime2 = psociTime() - tempfull;
      if ( GA::nodeid() == 0 ) cout << GA::nodeid() << " SCAT: numscatters cimatTime getV accTime compTime1 comptTime2 tempTime2 scatTime packTime "<< numscatters << " " << cimatTime<<" "<<  getVTime<<" "<<accTime << " "<<compTime1<<" "<<compTime2<<" "<< tempTime2 << " " <<  scatTime<<" "<<packTime <<" " << endl;
      return(hvcount);
    }

/* This version attmeps to call ga_scatter at a much lower level inside of GA
   to get around excessive memory allocations
*/
int PsociGAhamiltonian::matrixVectorProductsSplitChunkScatterAccumulateLowLevel( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
								      int * hv_dims, int * v_lo, int * v_hi )
 {
#ifndef SUPPLEMENTAL
   GA::error(" Must compile code with -DSUPPLEMENTAL ",-1);
#endif
#ifdef SQUAREHAMILTONIAN
   GA::error(" Must NOT compile with -DSQUAREHAMILTONIAN",1);
#endif
   
   GANbhdl * g_handle; // Try loop over tiny nbAcc calls instead of scatter

   // Check basic array dimensions
   
   int type;
   int dims[ 2 ]; //used by all
   int idim;
   int l_maxsefs = local_dims[0];
   
   char op[] = "+";
   
   /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
      within GA.  In the future this will be relaxed.
      Check them out here
   */
   /* cimat and friends are distributed over isefs but no breaking along sparity axis */
  
  int lo[2], hi[2];
  int cilo[2], cihi[2];
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  
  cilo[0] = clow;
  cilo[1] = local_cilo[1];
  cihi[0] = chi;
  cihi[1] = local_cihi[1];
  
  // these are the local vector piece
  
  double temp;
  double tempTime= psociTime();
  double tempTime2, getTime = 0.0, accTime=0.0, brdTime=0.0;
  double temp2;
  
  // Loop over all sefs, fetching LOCAL cimat for processing.
  // Gather Scatter objects can never be bigger than this
  const int nchunk = l_nchunk;
  
  int chunk_width;
  int ch_start=0;
  int ch_end=0;
  (l_maxsefs % nchunk==0)? chunk_width = l_maxsefs / nchunk: chunk_width = 1 + (l_maxsefs / nchunk);
  if ( GA::nodeid() == 0 ) cout << "exp VECTOR get CHUNK is " << nchunk << " " << chunk_width << " "<<MAX_AGGREGATE<<endl;
#ifdef USEMPI 
  if ( GA::nodeid() == 0 ) cout << "MPI collectives for H*v " << endl;
#endif
  
  double tempfull = psociTime(); 
  vector<double> testVector( chunk_width , 0.0 );
  
  int hvcount=0;
  int n=1;
  int number; // how long was the relevant CI row ?
  double dOne = 1.0;
  double alpha = 1.0;
  double packTime=0.0, scatbuf=0.0, deflateTime=0.0,  mapTime=0.0, listTime=0.0,compTime1=0.0, compTime2=0.0, scatTime=0.0, getVTime=0.0;
  
  double temptime, cimatTime=0.0;
  int chunk_lo[2];
  int chunk_hi[2];
 
  int width = MAX_AGGREGATE;  // pertain to H prefeting
  
  /* For bufV upload after GOP*/
  lo[0] = v_lo[0];
  hi[0] = v_hi[0];
  lo[1] = 0;
  hi[1] = 0;
  int index = v_lo[0]; // For bufV and testvector

// New scatter based stuff

  int numdeflates;
  int deflate_index;
  int newnumber;
  int nummapcontracts=0; // number of time minimap was copied to buffer_map

  double zeroTime = 0.0;

  int numscatters = 0;
  int newrealnumber = 0;

  int zerowidth = MAX_AGGREGATE;
  zeroMatrix( bufV2, zerowidth  );
  
  for(int ivec=start; ivec< end; ++ivec ) 
    {
      chunk_lo[1] = ivec;
      chunk_hi[1] = ivec;

#ifdef MATRIXBALANCE
      zeroMatrix( bufV2, zerowidth  );
#endif

int indexnewarray = 0;

#ifdef MATRIXBALANCE
//prefetch CIMAT data based on the some kind of chunking

// Moving this inside of the chunk loop increases the chance of better load balancing

       int newistart = -1;

       int max_num;
       int cimat_size = 0;
       vector<int> newlist;
       int indexsum = 0; // new

        int fetch_cimat_count=-1;

        for(int itci=0; itci < cimat_list.size(); itci+= MAX_AGGREGATE ){ // a factor of 5 in time variation
          int mxnum = cimat_list.size();
          max_num = MAX_AGGREGATE;
          if ( ( mxnum - itci ) < MAX_AGGREGATE ) max_num = mxnum - itci;

          temptime=psociTime();
          newlist.clear();
          for(int j=0; j< max_num; ++j) {
            newlist.push_back( cimat_list[ itci+j ]);
          }
          fetchsubListHamiltonianData( newlist ); // cimat still populate starting at 0 

          cimat_size += newlist.size(); // sums for the start of the clo iterator
          cimatTime += psociTime() - temptime;

// zero Out array only as needed.
// at this point we've prefetch as many as MAX_AGGREGATE H rows into cimat,icicol,number

#endif

      vector<int>::iterator it;
      for(it=my_chunklist.begin(); it!=my_chunklist.end(); ++it) { // this potentially staggers accross cores 
	
        temp = psociTime();
	zeroTime += psociTime() - temp;

        ch_start = chunk_width*(*it);
        ch_end = min( ch_start + chunk_width, l_maxsefs );

        //    ch_start=0;
        //    ch_end = l_maxsefs -1 ;
        chunk_lo[0] = ch_start;
        chunk_hi[0] = ch_end-1;

        // Step One: Grab the current GA-vector chunk for processing

	temp = psociTime();
        fetchAndReplicateVectorChunk((*it), ivec, chunk_width, v_lo, v_hi, testVector, g_v );
	getVTime += psociTime() - temp;
	
// Hamiltonian reduction
	
	int jhi;
	int buf_index=-1; // always restart
        // indexnewarray = -1

#ifdef MATRIXBALANCE
        indexnewarray = indexsum-1;
#else
        indexnewarray = -1; // because we always do the entire local list
#endif
          vector<int> tempindex;
#ifdef MATRIXBALANCE
	  for(int it2=0; it2< max_num; ++it2) { // extra loop
	    int i = newlist[ it2 ]; 
#else
	  for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
#endif
	      ++buf_index;
              ++indexnewarray;
#ifndef MATRIXBALANCE
              // local data is grabbed in this case
              // if not here then we already prefetched at the head of the method

	      buf_index=0; //override and simply grab local i
	      fetchLocalHblock( i, local_cimat[buf_index], local_icicol[buf_index], local_number[buf_index] ); // It actually does not need to be local but beware performance
#endif
	      int start_rng=local_icicol[buf_index][0];
	      int end_rng=local_icicol[buf_index][local_number[buf_index]-1];
	      number = local_number[buf_index];
	      
#ifdef LOOPNPUT
              if ( g_handle != NULL ) GA::nbWait( g_handle );
#endif

	      double dtemp=0.0;
              bufVnumber[ indexnewarray ] = number; // how many were found
	      if (!ch_end-1 < start_rng  || !ch_start > end_rng ) { //none within the range
		temp = psociTime();
		for(int j=0; j< number; ++j ) {
		  jhi = local_icicol[buf_index][j];
		  if ( jhi >= ch_start && jhi < ch_end ) {
		    dtemp += local_cimat[buf_index][ j ] * testVector[ jhi - ch_start]; 
		  }
		}
                bufV2[ indexnewarray ][ number-1 ] += dtemp;
                bufVIndex[ indexnewarray ][ number-1 ] = local_icicol[buf_index][number-1]; // no sum it is always the same
		compTime1+=psociTime() - temp;
	      }

	      if ( i >= ch_start && i < ch_end ) {
                temp = psociTime();
		double testVec = testVector[ i - ch_start ]; 
		for(int j=0; j< number-1; ++j ) {
	     	    int jhi = local_icicol[buf_index][j];
		    bufV2[ indexnewarray ][ j ] += local_cimat[buf_index][ j ] * testVec;
                    bufVIndex[ indexnewarray ][ j ] = jhi; // no need to sum
		}
		compTime2+=psociTime() - temp;
	      } // ch_start screen

                if ( indexnewarray >= MAX_AGGREGATE-1 ) {
                 ++numscatters;
                temp = psociTime();

#ifdef TACC
//Best for TACC
      newrealnumber = packDataNoZerosTransposed(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray ); // pnga_* directly
#endif
#ifdef LOOPNPUT
      newrealnumber = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray );
#endif
#ifdef SCATFLAT
      newrealnumber = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray );
#endif
#ifdef SCATNOFLAT
      indexarray2d = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, newrealnumber, dataarray ); // nga_scatter_acc()
#endif

     //          indexarray2d = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, newrealnumber, dataarray ); // nga_scatter_acc()
     //          newrealnumber = packDataNoZerosTransposed(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray ); // pnga_* directly
     //          newrealnumber = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray );
                
                packTime += psociTime() - temp;
                temp = psociTime();

                if ( newrealnumber > 0 ) {
#ifdef TACC
#ifdef GAEXP 
//not so much hopper as it is using ga-exp1 woth a different interface

      indexToInteger( indexarray, indexarrayInteger );
      pnga_scatter_acc( g_Hv->handle(), (void *)&dataarray[0], (Integer*)&indexarrayInteger[0], newrealnumber, &alpha);
#else
      pnga_scatter_acc( g_Hv->handle(), (void *)&dataarray[0], &indexarray[0], newrealnumber, &alpha); // low level call
#endif
#endif
#ifdef LOOPNPUT
      LoopNputMethod(g_Hv,indexarray, dataarray, g_handle );
#endif
#ifdef SCATFLAT
      NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&dataarray[0], (int*)&indexarray[0], newrealnumber, &alpha);
#endif
#ifdef SCATNOFLAT
      NGA_Scatter_acc(g_Hv->handle(), (void *)&dataarray[0], indexarray2d, newrealnumber, &alpha);
#endif

     //           NGA_Scatter_acc(g_Hv->handle(), (void *)&dataarray[0], indexarray2d, newrealnumber, &alpha);
     //           pnga_scatter_acc( g_Hv->handle(), (void *)&dataarray[0], &indexarray[0], newrealnumber, &alpha); 
     //           NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&dataarray[0], (int*)&indexarray[0], newnumber, &alpha);
     //           LoopNputMethod(g_Hv,indexarray, dataarray, g_handle );
                ++numscatters;
                }
#ifdef SCATNOFLAT
     emptyIndexarray( newrealnumber, indexarray2d ); //  nga_scatter_acc()
#endif
                 scatTime += psociTime() - temp;
                 indexnewarray = -1;
             }
          } // clow

// We could conceivably do the scatter here as well

#ifdef INNERSCATTER
      temp = psociTime();
#ifdef TACC
      newrealnumber = packDataNoZerosTransposed(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray ); // pnga_* directly
#endif
#ifdef LOOPNPUT
      newrealnumber = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray );
#endif
#ifdef SCATFLAT
      newrealnumber = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray );
#endif
#ifdef SCATNOFLAT
      indexarray2d = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, newrealnumber, dataarray ); // nga_scatter_acc()
#endif

//      indexarray2d = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, newrealnumber, dataarray );
//      newrealnumber = packDataNoZerosTransposed(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray ); // pnga_* directly
//      newrealnumber = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray );

      packTime += psociTime() - temp;
      temp = psociTime();
      if ( newrealnumber > 0 ) {
#ifdef TACC
#ifdef GAEXP 
      indexToInteger( indexarray, indexarrayInteger );
      pnga_scatter_acc( g_Hv->handle(), (void *)&dataarray[0], (Integer*)&indexarrayInteger[0], newrealnumber, &alpha);
#else
      pnga_scatter_acc( g_Hv->handle(), (void *)&dataarray[0], &indexarray[0], newrealnumber, &alpha); // low level call
#endif
#endif
#ifdef LOOPNPUT
      LoopNputMethod(g_Hv,indexarray, dataarray, g_handle );
#endif
#ifdef SCATFLAT
      NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&dataarray[0], (int*)&indexarray[0], newrealnumber, &alpha);
#endif
#ifdef SCATNOFLAT
      NGA_Scatter_acc(g_Hv->handle(), (void *)&dataarray[0], indexarray2d, newrealnumber, &alpha);
#endif
//    NGA_Scatter_acc(g_Hv->handle(), (void *)&dataarray[0], indexarray2d, newrealnumber, &alpha);
//    pnga_scatter_acc( g_Hv->handle(), (void *)&dataarray[0], &indexarray[0], newrealnumber, &alpha); 
//    LoopNputMethod(g_Hv,indexarray, dataarray, g_handle );

      ++numscatters;
      }

#ifdef SCATNOFLAT
      emptyIndexarray( newrealnumber, indexarray2d );
#endif
      scatTime += psociTime() - temp;
#endif
        } // chunk_list 

#ifdef MATRIXBALANCE
         indexsum += newlist.size(); // set to beginniing of next chunk
         }  
#endif 
          indexnewarray = -1;

/* Now we want to take the entire set of results, pack them into dataarray an dindex array 
   and scatter back to the nodes
   for a GIVEN ivec
*/

#ifndef INNERSCATTER
      temp = psociTime();
#ifdef TACC
      newrealnumber = packDataNoZerosTransposed(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray ); // pnga_* directly
#endif
#ifdef LOOPNPUT
      newrealnumber = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray );
#endif
#ifdef SCATFLAT
      newrealnumber = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray );
#endif
#ifdef SCATNOFLAT
      indexarray2d = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, newrealnumber, dataarray ); // nga_scatter_acc()
#endif
//    indexarray2d = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, newrealnumber, dataarray );
//    newrealnumber = packDataNoZerosTransposed(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray ); // pnga_* directly
//    newrealnumber = packDataNoZeros(ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray );

      packTime += psociTime() - temp;
      temp = psociTime();
      if ( newrealnumber > 0 ) {
#ifdef TACC
#ifdef GAEXP 
      indexToInteger( indexarray, indexarrayInteger );
      pnga_scatter_acc( g_Hv->handle(), (void *)&dataarray[0], (Integer*)&indexarrayInteger[0], newrealnumber, &alpha);
#else
      pnga_scatter_acc( g_Hv->handle(), (void *)&dataarray[0], &indexarray[0], newrealnumber, &alpha); // low level call
#endif
#endif
#ifdef LOOPNPUT
      LoopNputMethod(g_Hv,indexarray, dataarray, g_handle );
#endif
#ifdef SCATFLAT
      NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&dataarray[0], (int*)&indexarray[0], newrealnumber, &alpha);
#endif
#ifdef SCATNOFLAT
      NGA_Scatter_acc(g_Hv->handle(), (void *)&dataarray[0], indexarray2d, newrealnumber, &alpha);
#endif
//      NGA_Scatter_acc(g_Hv->handle(), (void *)&dataarray[0], indexarray2d, newrealnumber, &alpha);
//      pnga_scatter_acc( g_Hv->handle(), (void *)&dataarray[0], &indexarray[0], newrealnumber, &alpha); // low level call
//      LoopNputMethod(g_Hv,indexarray, dataarray, g_handle );

      ++numscatters;
      }
#ifdef SCATNOFLAT
    emptyIndexarray( newrealnumber, indexarray2d );
#endif
      scatTime += psociTime() - temp;

#endif
      ++hvcount;

      }  // ivec 

      GA::SERVICES.sync();
      tempTime2 = psociTime() - tempfull;
      if ( GA::nodeid() == 0 ) cout << GA::nodeid() << " SCAT: numscatters cimatTime getV accTime compTime1 comptTime2 tempTime2 scatTime packTime "<< numscatters << " " << cimatTime<<" "<<  getVTime<<" "<<accTime << " "<<compTime1<<" "<<compTime2<<" "<< tempTime2 << " " <<  scatTime<<" "<<packTime <<" " << endl;


      indexarray.clear();
      dataarray.clear();
      indexarrayInteger.clear();

      return(hvcount);
    }

/* HIGHLY experimental DIRECT vector method
*/

/* New LOwLevel method that does a better job of managing merging and accumuation of dataarray elements
 * The -DSCATNOFLAT flag has been removed since it never really was shown to be faster
 * -DLOOPNPUT has also been remove since it was never shown to be faster
 */
int PsociGAhamiltonian::matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
								      int * hv_dims, int * v_lo, int * v_hi )
 {
#ifndef SUPPLEMENTAL
   GA::error(" Must compile code with -DSUPPLEMENTAL ",-1);
#endif
#ifdef SQUAREHAMILTONIAN
   GA::error(" Must NOT compile with -DSQUAREHAMILTONIAN",1);
#endif
#ifdef SCATNOFLAT
GA::error(" SCATNOFLAT not supported in matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate",1);
#endif
#ifdef LOOPNPUT
GA::error(" LOOPNPUT not supported in matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate",1);
#endif
#ifdef USEMPI
GA::error(" USEMPI not supported in matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate",1);
#endif

   // Check basic array dimensions
   
   int type;
   int dims[ 2 ]; //used by all
   int idim;
   int l_maxsefs = local_dims[0];
   
   char op[] = "+";
   
   /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
      within GA.  In the future this will be relaxed.
      Check them out here
   */
   /* cimat and friends are distributed over isefs but no breaking along sparity axis */
  
  int lo[2], hi[2];
  int cilo[2], cihi[2];
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  
  cilo[0] = clow;
  cilo[1] = local_cilo[1];
  cihi[0] = chi;
  cihi[1] = local_cihi[1];
  
  // these are the local vector piece
  
  double temp;
  double tempTime= psociTime();
  double fetchLocal = 0.0, tempTime2, getTime = 0.0, accTime=0.0, brdTime=0.0;
  double temp2;
  
  // Loop over all sefs, fetching LOCAL cimat for processing.
  // Gather Scatter objects can never be bigger than this
  const int nchunk = l_nchunk;
  
  int chunk_width;
  int ch_start=0;
  int ch_end=0;
  (l_maxsefs % nchunk==0)? chunk_width = l_maxsefs / nchunk: chunk_width = 1 + (l_maxsefs / nchunk);
  if ( GA::nodeid() == 0 ) cout << "exp VECTOR get CHUNK is " << nchunk << " " << chunk_width << " "<<MAX_AGGREGATE<<endl;
  
  double tempfull = psociTime(); 

  vector<double> testVector( chunk_width , 0.0 );
  
  int hvcount=0;
  int n=1;
  int number; // how long was the relevant CI row ?
  double dOne = 1.0;
  double alpha = 1.0;
  double packTime=0.0, scatbuf=0.0, deflateTime=0.0,  mapTime=0.0, listTime=0.0,compTime1=0.0, compTime2=0.0, scatTime=0.0, getVTime=0.0;
  
  double temptime, cimatTime=0.0;
  int chunk_lo[2];
  int chunk_hi[2];
 
  int width = MAX_AGGREGATE;  // pertain to H prefeting
  
  /* For bufV upload after GOP*/
  lo[0] = v_lo[0];
  hi[0] = v_hi[0];
  lo[1] = 0;
  hi[1] = 0;
  int index = v_lo[0]; // For bufV and testvector

// New scatter based stuff

  int numdeflates;
  int deflate_index;
  int newnumber;
  int nummapcontracts=0; // number of time minimap was copied to buffer_map

  double zeroTime = 0.0;

  int numscatters = 0;
  int newrealnumber = 0;

  int zerowidth = MAX_AGGREGATE;
  zeroMatrix( bufV2, zerowidth  ); // suspicious why not AGGREGATE size?

  int me = GA::nodeid();
  int size = GA::nodes();
  if ( me == 0 ) cout << "LL STRIDED CIMAT TEST" << endl;

#ifdef TACC
      bool transpose=true;
#else
      bool transpose=false;
#endif
  
int estimate_elements=0;

  for(int ivec=start; ivec< end; ++ivec ) 
    {
      chunk_lo[1] = ivec;
      chunk_hi[1] = ivec;

//cout << "NEW VECTOR " << ivec << endl;

//  these are constrained to MAX_AGGREGATE when building for -DNEWACCUMULATEDSCATTER
//
#ifdef MATRIXBALANCE
      zeroMatrix( bufV2, zerowidth  );
#endif

  int indexnewarray = 0;

#ifdef MATRIXBALANCE
// prefetch CIMAT data based on the some kind of chunking
// Moving this inside of the chunk loop increases the chance of better load balancing

       int newistart = -1;

       int max_num;
       int cimat_size = 0;
       vector<int> newlist;
       int indexsum = 0; // new

        int fetch_cimat_count=-1;

        for(int itci=0; itci < cimat_list.size(); itci+= MAX_AGGREGATE ){ // a factor of 5 in time variation
          int mxnum = cimat_list.size();
          max_num = MAX_AGGREGATE;
          if ( ( mxnum - itci ) < MAX_AGGREGATE ) max_num = mxnum - itci;

          temptime=psociTime();
          newlist.clear();
          for(int j=0; j< max_num; ++j) {
            newlist.push_back( cimat_list[ itci+j ]);
          }
          fetchsubListHamiltonianData( newlist ); // cimat still populate starting at 0 

          cimat_size += newlist.size(); // sums for the start of the clo iterator
          cimatTime += psociTime() - temptime;

// zero Out array only as needed.
// at this point we've prefetch as many as MAX_AGGREGATE H rows into cimat,icicol,number
#endif

      vector<int>::iterator it;
      for(it=my_chunklist.begin(); it!=my_chunklist.end(); ++it) { // this potentially staggers accross cores 
	
        temp = psociTime();
	zeroTime += psociTime() - temp;

        ch_start = chunk_width*(*it);
        ch_end = min( ch_start + chunk_width, l_maxsefs );
        chunk_lo[0] = ch_start;
        chunk_hi[0] = ch_end-1;

        // Step One: Grab the current GA-vector chunk for processing

//maybe scatter chiunk reads but that wasnot too helpful in the past

	temp = psociTime();
        fetchAndReplicateVectorChunk((*it), ivec, chunk_width, v_lo, v_hi, testVector, g_v );
	getVTime += psociTime() - temp;
	
// Hamiltonian reduction
	
	int jhi;
	int buf_index=-1; // always restart
        // indexnewarray = -1

#ifdef MATRIXBALANCE
        indexnewarray = indexsum-1;
#else
        indexnewarray = -1; // because we always do the entire local list
#endif
          vector<int> tempindex;
#ifdef MATRIXBALANCE
	  for(int it2=0; it2< max_num; ++it2) { // extra loop
	    int i = newlist[ it2 ]; 
#else
          estimate_elements=0; // to trigger dump phase

	  for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
       //   chi = l_maxsefs-1;
       //   for(int i=me; i< l_maxsefs; i+= size ) { 
#endif
	      ++buf_index;
              ++indexnewarray;
#ifndef MATRIXBALANCE
              // local data is grabbed in this case
              // if not here then we already prefetched at the head of the method

	      buf_index=0; //override and simply grab local i
              temp = psociTime();
	      fetchLocalHblock( i, local_cimat[buf_index], local_icicol[buf_index], local_number[buf_index] ); // It actually does not need to be local but beware performance
              fetchLocal += psociTime() - temp;
#endif
	      int start_rng=local_icicol[buf_index][0];
	      int end_rng=local_icicol[buf_index][local_number[buf_index]-1];
	      number = local_number[buf_index];

	      double dtemp=0.0;
              bufVnumber[ indexnewarray ] = number; // how many were found
              temp = psociTime();
	      if (!ch_end-1 < start_rng  || !ch_start > end_rng ) { //none within the range
		for(int j=0; j< number; ++j ) {
		  jhi = local_icicol[buf_index][j];
		  if ( jhi >= ch_start && jhi < ch_end ) {
		    dtemp += local_cimat[buf_index][ j ] * testVector[ jhi - ch_start]; 
		  }
		}
                bufV2[ indexnewarray ][ number-1 ] += dtemp;
                bufVIndex[ indexnewarray ][ number-1 ] = local_icicol[buf_index][number-1]; // no sum it is always the same

                ++estimate_elements;
	      }
              compTime1+=psociTime() - temp;


	      if ( i >= ch_start && i < ch_end ) {
                temp = psociTime();
                ++estimate_elements;
		double testVec = testVector[ i - ch_start ]; 
		for(int j=0; j< number-1; ++j ) {
	     	    int jhi = local_icicol[buf_index][j];
		    bufV2[ indexnewarray ][ j ] += local_cimat[buf_index][ j ] * testVec;
                    bufVIndex[ indexnewarray ][ j ] = jhi; // no need to sum
		}
		compTime2+=psociTime() - temp;
	      } // ch_start screen
                   
// see if we can eliminate the multidimentionsl bufV2 and friends.
//
//
/*
      temp = psociTime();
      newrealnumber = packDataNoZerosAccumulateSingle(transpose, indexnewarray, ivec, bufVnumber[indexnewarray], bufVIndex[indexnewarray], bufV2[indexnewarray], indexarray, dataarray );
      packTime += psociTime() - temp;
*/


      const int maxelements = 1*1024*1044; // limits scatter stack to something like 16 MB
         if ( (indexnewarray >= MAX_AGGREGATE-1 ) || estimate_elements > maxelements || i >= chi) { // 100000 value should tie to GA stack

                 temp = psociTime();
                 newrealnumber = packDataNoZerosAccumulate(transpose, indexnewarray, ivec, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray );
                 packTime += psociTime() - temp;

                 temp = psociTime();
                 //cout << GA::nodeid() << " trigger a dump at " << indexnewarray << " with estimate_elements " << estimate_elements << endl;

                 if ( newrealnumber > 0 ) { // double check
                     ++numscatters;
#ifdef TACC
#ifdef GAEXP 
      indexToInteger( indexarray, indexarrayInteger );
      pnga_scatter_acc( g_Hv->handle(), (void *)&dataarray[0], (Integer*)&indexarrayInteger[0], newrealnumber, &alpha);
#else
      pnga_scatter_acc( g_Hv->handle(), (void *)&dataarray[0], &indexarray[0], newrealnumber, &alpha); // low level call
#endif
#endif
#ifdef SCATFLAT
      NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&dataarray[0], (int*)&indexarray[0], newrealnumber, &alpha);
#endif
                indexarray.clear(); // last step cleanup
                dataarray.clear();
                estimate_elements = 0;

                }
                 scatTime += psociTime() - temp;
                 indexnewarray = -1;
             }
          } // clow
        } // chunk_list 

#ifdef MATRIXBALANCE
         indexsum += newlist.size(); // set to beginniing of next chunk
         }  
#endif 
          indexnewarray = -1;

      ++hvcount;
      if ( dataarray.size() > 0 || indexarray.size() > 0  ) GA::error("?Index/Dataarray should never be zero at the end " ,dataarray.size());
      }  // ivec 

    //  GA::SERVICES.sync();
      tempTime2 = psociTime() - tempfull;
      if( GA::nodeid() == 0 ) cout << GA::nodeid() << " NEW ACCUMULATE SCAT: numscatters fetchLocal cimatTime getV accTime compTime1 comptTime2 tempTime2 scatTime packTime "<< numscatters << " "<<fetchLocal<<" " << cimatTime<<" "<<  getVTime<<" "<<accTime << " "<<compTime1<<" "<<compTime2<<" "<< tempTime2 << " " <<  scatTime<<" "<<packTime <<" " << endl;

      GA::SERVICES.sync();

      indexarray.clear(); // doesn;t hurt to do it again
      dataarray.clear();
      indexarrayInteger.clear();

      return(hvcount);
    }

/* In this method we break up the inut VECTOR into several chunks. Also, we break up bufV into several chunks
   and then do a simple ACC of bufV at the end of the I loop.
   THough we repeat the ILOOP jchunk times, it will still be much faster than doing a scatter
*/
int PsociGAhamiltonian::matrixVectorProductsSplitChunkMultipleJchunks( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
								      int * hv_dims, int * v_lo, int * v_hi )
 {
#ifndef SUPPLEMENTAL
   GA::error(" Must compile code with -DSUPPLEMENTAL ",-1);
#endif
#ifdef SQUAREHAMILTONIAN
   GA::error(" Must NOT compile with -DSQUAREHAMILTONIAN",1);
#endif
#ifdef SCATNOFLAT
GA::error(" SCATNOFLAT not supported in matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate",1);
#endif
#ifdef LOOPNPUT
GA::error(" LOOPNPUT not supported in matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate",1);
#endif
#ifdef MATRIXBALANCE
GA::error(" MATRIXBALANCE not supported",1);
#endif

   // Check basic array dimensions
   
   int my_rank = GA::nodeid();
   int my_size = GA::nodes();
#ifdef NXTASKVECTORS
   if (my_rank == 0 ) cout << "NXTASK for M*V" << endl;
#endif

#ifdef REVERSE
    if (my_rank == 0 ) cout << "REVERSE loop order" << endl;
#endif

   int type;
   int dims[ 2 ]; //used by all
   int idim;
   int l_maxsefs = local_dims[0];
   
   char op[] = "+";
   
   /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
      within GA.  In the future this will be relaxed.
      Check them out here
   */
   /* cimat and friends are distributed over isefs but no breaking along sparity axis */
  
  int lo[2], hi[2];
  int cilo[2], cihi[2];
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  
  cilo[0] = clow;
  cilo[1] = local_cilo[1];
  cihi[0] = chi;
  cihi[1] = local_cihi[1];
  
  // these are the local vector piece
  
  double temp;
  double tempTime= psociTime();
  double fetchLocal = 0.0, tempTime2, getTime = 0.0, accTime=0.0, brdTime=0.0;
  double gopTime=0.0;
  double temp2;
  
  // Loop over all sefs, fetching LOCAL cimat for processing.
  // Gather Scatter objects can never be bigger than this
  const int nchunk = l_nchunk;
  
  int chunk_width;
  int ch_start=0;
  int ch_end=0;

  int jch_start=0;
  int jch_end=0;

  (l_maxsefs % nchunk==0)? chunk_width = l_maxsefs / nchunk: chunk_width = 1 + (l_maxsefs / nchunk);
  if ( GA::nodeid() == 0 ) cout << " exp VECTOR get CHUNK is " << nchunk << " " << chunk_width << " "<<MAX_AGGREGATE<<endl;
  if ( GA::nodeid() == 0 ) cout << " Will break up Jchunk loop into chunks " << nchunk << endl;
  if ( GA::nodeid() == 0 ) cout << " NEW OCT FETCH LOOPS ARE DEEPER" << endl;
  //if ( GA::nodeid() == 0 ) cout << "NON_LOCAL H lops FETCH LOOPS ARE DEEPER" << endl;
  
// THis applies to both vector and bufV
  int my_vector_low = v_lo[0];
  int my_vector_hi = v_hi[0]; // this does not need subtraction by 1
  int newlow[2], newhi[2];

  double tempfull = psociTime(); 

  vector<double> testVector( chunk_width , 0.0 );

//  vector<double> bufV( chunk_width , 0.0 );
  
  int hvcount=0;
  int n=1;
  int number; // how long was the relevant CI row ?
  double dOne = 1.0;
  double alpha = 1.0;
  double packTime=0.0, scatbuf=0.0, deflateTime=0.0,  mapTime=0.0, listTime=0.0,compTime1=0.0, compTime2=0.0, scatTime=0.0, getVTime=0.0;
  
  double temptime, cimatTime=0.0;
  int chunk_lo[2];
  int chunk_hi[2];
  int width = MAX_AGGREGATE;  // pertain to H prefeting
  
  int jchunk_lo[2];
  int jchunk_hi[2];

  /* For bufV upload after GOP*/
  lo[0] = v_lo[0];
  hi[0] = v_hi[0];
  lo[1] = 0;
  hi[1] = 0;
  int index = v_lo[0]; // For bufV and testvector

// New scatter based stuff

  int numdeflates;
  int deflate_index;
  int newnumber;
  int nummapcontracts=0; // number of time minimap was copied to buffer_map

  double zeroTime = 0.0;

  int numscatters = 0;
  int newrealnumber = 0;
  int zerowidth = MAX_AGGREGATE;

  int me = GA::nodeid();
  int size = GA::nodes();
  if ( me == 0 ) cout << "NEW STRIDED CIMAT TEST: matrixVectorProductsSplitChunkMultipleJchunk" << endl;

 // start NXTVAL driver
#ifdef NXTASKVECTORS
  PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
#ifdef REVERSE
    jobs.initNxtask( l_maxsefs );
#else
    jobs.initNxtask();
#endif
#endif

  vector<double> bufV( chunk_width , 0.0 );

      int numaccs=0, numgets=0;


  for(int ivec=start; ivec< end; ++ivec ) 
    {
      chunk_lo[1] = ivec;
      chunk_hi[1] = ivec;
   
      newlow[1]=ivec;
      newhi[1]=ivec;

//
      int indexnewarray = 0;

// add outer loop over JCHUNK/bufV loops
// Step Zero: Specify the current GA-bufV chunk for processing

      vector<int>::iterator jit;
      //vector<double> bufV( chunk_width , 0.0 ); 

      //if ( me == 0 ) cout << "TEST chunkwidth " <<GA::nodeid() << " " << chunk_width<<" "<<v_lo[0]<<" "<<v_hi[0]<<" "<<v_lo[1]<<" "<<v_hi[1]<< endl;

int testcount=0;

      if ( nchunk != my_chunklist.size() ) GA::error("nchunk != my_chunklist.size() ",my_chunklist.size());


      for(jit=my_chunklist.begin(); jit!=my_chunklist.end(); ++jit) { // this potentially staggers accross cores 

      int jchunk_width = chunk_width;

      zeroVector( bufV ); //  do we need this ?

      jch_start = jchunk_width*(*jit);
      jch_end = min( jch_start + jchunk_width, l_maxsefs );
      jchunk_lo[0] = jch_start;
      jchunk_hi[0] = jch_end-1;

// Inner loop

      vector<int>::iterator it;
      for(it=my_chunklist.begin(); it!=my_chunklist.end(); ++it) { // this potentially staggers accross cores 

        ch_start = chunk_width*(*it);
        ch_end = min( ch_start + chunk_width, l_maxsefs );
        chunk_lo[0] = ch_start;
        chunk_hi[0] = ch_end-1;

        // Step One: Grab the current GA-vector chunk for processing

	temp = psociTime();
        if ( chunk_width > testVector.size() ) GA::error("  chunk_width > testVector.size()",chunk_width);

        //cout << GA::nodeid() << "CHECK jchunk_width chunk_width and  testVector.size() " << jchunk_width << " " << chunk_width<<" "<< testVector.size() <<" "<<ivec<<" "<<endl;


// What none are actually required ?
// Perhaps we should wait ?

#ifdef USEMPI
        fetchAndReplicateVectorChunkMPI((*it), ivec, chunk_width, v_lo, v_hi, testVector, g_v );
#else
        fetchAndReplicateVectorChunk((*it), ivec, chunk_width, v_lo, v_hi, testVector, g_v );
#endif
        numgets++;

        ++testcount; // delete me

	getVTime += psociTime() - temp;
	
// Hamiltonian reduction
	
	int jhi;
	//int buf_index=-1; // always restart
        // indexnewarray = -1

        indexnewarray = -1; // because we always do the entire local list
        //  vector<int> tempindex;

//cout << GA::nodeid() << " I loop " << endl; 
//Poor load balancxe

//        for (int i=my_rank; i< l_maxsefs; i+=my_size ) { 
//	  for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef

#ifdef NXTASKVECTORS
#ifdef REVERSE
          jobs.nxtaskReset( l_maxsefs );
          long nexttask = jobs.nxtaskRev( my_size, l_maxsef );
#else
          jobs.nxtaskReset();
          long nexttask = jobs.nxtask( my_size ); // 
#endif
#endif

         int start_rng=0;
         int end_rng=0;

#ifdef NXTASKVECTORS
#ifdef REVERSE
          for (int i=l_maxsefs-1; i>0; --i ) {
#else
          for (int i=0; i< l_maxsefs; ++i ) { 
#endif
             if ( nexttask == i ) { // in the loop stride is implemented 
             //cout << my_rank << " nexttask " << i <<" "<<nexttask << endl;

#else
//              for (int i=my_rank; i< l_maxsefs; i+=my_size ) { 
             for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
#endif

/*
              temp = psociTime();
	      fetchLocalHblock( i, local_cimat[buf_index], local_icicol[buf_index], local_number[buf_index] ); 
              fetchLocal += psociTime() - temp;
	      start_rng=local_icicol[buf_index][0];
	      end_rng=local_icicol[buf_index][local_number[buf_index]-1];
	      number = local_number[buf_index];
*/

// what J range are we looking at?

	      double dtemp=0.0;
              temp = psociTime();
           
              const int buf_index=0; //override and simply grab local i
              if ( i >= jch_start && i < jch_end || i >= ch_start && i < ch_end ) {
                temp = psociTime();
                fetchLocalHblock( i, local_cimat[buf_index], local_icicol[buf_index], local_number[buf_index] );
                fetchLocal += psociTime() - temp;
                start_rng=local_icicol[buf_index][0];
                end_rng=local_icicol[buf_index][local_number[buf_index]-1];
                number = local_number[buf_index];
              } 


              if ( i >= jch_start && i < jch_end ) {
	      if (!ch_end-1 < start_rng  || !ch_start > end_rng ) { //none within the range
		for(int j=0; j< number; ++j ) {
		  jhi = local_icicol[buf_index][j];
		  if ( jhi >= ch_start && jhi < ch_end ) {
		    dtemp += local_cimat[buf_index][ j ] * testVector[ jhi - ch_start]; 
                  }
		}
                bufV[ i-jch_start ] += dtemp; // 
	      }
              }// i > jch

              compTime1+=psociTime() - temp;

	      if ( i >= ch_start && i < ch_end ) {
                temp = psociTime();
		double testVec = testVector[ i - ch_start ]; 

		for(int j=0; j< number-1; ++j ) {
	     	    int jhi = local_icicol[buf_index][j];
                    if ( jhi >= jch_start && jhi < jch_end ) {
                    bufV[ jhi-jch_start ] += local_cimat[buf_index][ j ] * testVec;
                    //cout << "bufV " << bufV[ jhi-jch_start ] << endl;
                    }
		}
		compTime2+=psociTime() - temp;
	      } // ch_start screen

#ifdef NXTASKVECTORS

#ifdef REVERSE
              nexttask = jobs.nxtaskRev( my_size, l_maxsefs );
#else
              nexttask = jobs.nxtask( my_size );
#endif
             } // nexttask==i
#endif

           } // clow -- poor load balance still

        } // vector chunk_list 

        double dAlpha=1.0;

       tempTime2 = psociTime();

      if ( jchunk_width > bufV.size() ) GA::error("jchunk_width > bufV.size()",bufV.size() );
      temp = psociTime();
      accumulateAndPutProductChunk((*jit), ivec, jchunk_width, v_lo, v_hi, bufV, g_Hv );
      numaccs++;
      accTime += psociTime() - tempTime2; 

       } // Jchunk_list

       ++hvcount;
//if ( GA::nodeid() == 0 ) cout << "TEMPCOUNT FOR ivec " << ivec<<" is" << testcount<< endl;

      }  // ivec 

#ifdef NXTASKVECTORS
      jobs.destroyNxtask(); //Close'er up
#endif

      tempTime2 = psociTime() - tempfull;
      if ( GA::nodeid() == 0 ) cout << GA::nodeid() << "CHUNKED JCHUNKS ACC: fetchLocal cimatTime getV accTime compTime1 comptTime2 gopTime tempTime2 "<< " "<<fetchLocal<<" " << cimatTime<<" "<<  getVTime<<" "<<accTime << " "<<compTime1<<" "<<compTime2<<" "<< gopTime << " " << tempTime2 << " " <<  endl;
      // if ( GA::nodeid() <= 300 ) cout << " num calls " << numgets<<" "<<numaccs << endl;

      // not needed GA::SERVICES.sync();

      return(hvcount);
    }

/* EXPERIMENTAL DIRECT METHOD
*/
int PsociGAhamiltonian::matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulateDirect( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
								      int * hv_dims, int * v_lo, int * v_hi )
 {
#ifndef SUPPLEMENTAL
   GA::error(" Must compile code with -DSUPPLEMENTAL ",-1);
#endif
#ifdef SQUAREHAMILTONIAN
   GA::error(" Must NOT compile with -DSQUAREHAMILTONIAN",1);
#endif
#ifdef SCATNOFLAT
GA::error(" SCATNOFLAT not supported in matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate",1);
#endif
#ifdef LOOPNPUT
GA::error(" LOOPNPUT not supported in matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate",1);
#endif
#ifdef USEMPI
GA::error(" USEMPI not supported in matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate",1);
#endif

  int g_rank = GA::nodeid();
  int g_size = GA::nodes();

  if ( g_rank == 0 ) cout << "construct DIRECT matrix vector products" << endl;

  char handleMessage[] = "matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulateDirect: current CI vectors: g_v";
  g_v->checkHandle( handleMessage );
  
  double temp, temp2, temp3;
  double packTime=0.0, scatterTime=0.0, sociTime=0.0, updateTime =0.0, compareTime=0.0, computeConfsTime=0.0 ;
  double assemTime=0.0;
  int totalFetches=0;
  long allcount=0;
  long block_status;

  int hvcount=0;
  
/* bufVnumber, bufV2, bufVIndex are all preallocated to be MAX_AGGREGATE int 
   indexnewarray is used to index it
*/

  pair<int,double> info1;    // Need to loop over the H, fetch it into sefsef blocks and pass that to the below method
  
  //PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
  const int ilo = 0;
  const int ihi = l_maxspatials;
  
  if ( g_rank == 0 ) cout << " RECOMPUTE DIRECT SOCI instead of simply fetch them " << endl;

#ifdef NEWORBMAP
  GA::error(" PsociGAhamiltonian::matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulateDirect must not use -DNEWORBMAP",1);
#endif

  int countrun;

  int l_maxsefs = local_dims[0];
  int maxsefperconf = l_deters->fetchGlobalMaxSefPerSpatial();

  int numscatters=0;
  double alpha=1.0;

#ifdef TACC
    bool transpose=true;
#else
    bool transpose=false;
#endif

/*
    vector<vector<double> > cimat(1);
    vector<vector<int> > icicol(1);
    vector<int> number(1);
*/

// could also check on even odd for screening
      
    int didblock=0;
    int totalblock=0;
    int excitblock=0;
    int phaseblock=0;
    int possibleblock=0;

    vector<double> testVector; // hold maxsefperdet worht of vector data 

//must check parallel runs
// doesn't do anything
    vector<char> * rep_phases = l_deters->fetchReplicatedPhasesHandle();
    cout << "rep phase is size " << rep_phases->size()<<endl; // returns iopen

    int emptyfetches=0;
      
    for(int iroot=start; iroot < end; ++iroot ) {
      countrun=0;
      cout << "IN roots" << endl; 

      PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
      jobs.initNxtask();
      long nexttask = jobs.nxtask( g_size );
      
      int indexnewarray = -1; // indexarray preallocate to be MAX_AGGREGATE in length
      int estimate_elements=0;
      
          vector<double> testVector( l_maxsefs, 0.0);
//        testVector.resize( l_maxsefs );
          fetchAllVector( iroot, testVector, g_v ); // Note the value of vend

//      cout << "use j[hase directly " << endl;

      for(int iconf = ilo; iconf < ihi; ++iconf ) {
        int numElems = 0;

        totalblock += iconf; 
	
        if ( nexttask == iconf ) { // in the loop stride is implemented  

         vector<COMPRESS_SIZE> imap( l_nbf+2  );
         JOUTFG fetchi;

         temp2 = psociTime();
         l_deters->computeConfigsGA( iconf+1, imap, fetchi ); // resiuze imap then computeConfigs
         computeConfsTime += psociTime() - temp2;
         ++totalFetches;

// Now first identify all  easily removed abs(iopen-jopen) > 4 entries
// drops list down to 12% for a simple case

         temp2 = psociTime();
         int iphase = rep_phases->at(fetchi.index-1);
         vector<int> jphases;
         int numJs = assemblePhasePossibles(iconf, iphase, jphases, rep_phases );
         phaseblock += jphases.size();
         //cout << "assemblePhasePossibles sinple > 4 test returned " << numJs << " out of a total of " << iconf << endl;
         assemTime += psociTime() - temp2;
        
// Now check two things. 1 for a given comfigure are we more than 2 exciations different AND
// Do we have specific projections that match beyond simply even or odd?
// drops list down to < 1% for a simple case
  
          vector<int> list;

          temp2 = psociTime();
//already fetch above l_deters->fetchOrbMap(iconf+1, imap); // resize imap
//          countrun +=  l_deters->compareAllOrbMapMatrixVectorProducts(iconf+1, imap, list ); // no change to imap : fetch jlist for each i usually 100-300 entries
// this really should only be done ONCE

          countrun +=  l_deters->compareAllOrbMapMatrixVectorProducts(iconf+1, imap, jphases, list );
          jphases.clear();
          excitblock += list.size();
          compareTime += psociTime() - temp2;

          int width = fetchi.nsefi;
          vector<vector<double> > cimat( width );
          vector<vector<int> > icicol( width );
          vector<int> number( width );
/*
          set<int> iprojections;
          for(int k=0; k< fetchi.ndeti; ++k ) {
            iprojections.insert( fetchi.spin_projection[k] ) ;
          }
*/

          int totalElems=0;
          vector<int>::iterator jit;
          for(jit = list.begin(); jit != list.end(); ++jit )
            {
              int jconf = (*jit);
              vector<COMPRESS_SIZE> jmap( l_nbf+2  );
              JOUTFG fetchj;

              temp2 = psociTime();
              l_deters->computeConfigsGA( jconf+1, jmap, fetchj ); // jmap, fetchj allocated internally
              ++totalFetches;
              computeConfsTime += psociTime() - temp2;

/*
              set<int> jprojections;
              for(int k=0; k< fetchj.ndeti; ++k ) {
                jprojections.insert( fetchj.spin_projection[k] );
              }
*/

// Does the new J have ANY specific projections in common with the input iconf?

/*
           //bool process=false;
           for (set<int>::iterator sit=iprojections.begin(); sit != iprojections.end(); ++sit) {
              int temp=(*sit);
              for (set<int>::iterator sjt=jprojections.begin(); sjt != jprojections.end(); ++sjt) {
               cout << "projs "<<temp<<" "<<(*sjt)<<endl;
                  if ( temp == (*sjt) ) {
                  process=true;
                   cout <<"projs inner at " << temp<<" "<<(*sjt)<<endl;
                  goto processBlocks;
                  }
                }
            }
*/
           // cout << "if here no projection match " << iconf<<" "<<jconf<<endl;

           processBlocks: 
           bool process=true;

           if ( process ) {

//cout << "Compute a bocks" << endl;

              ++didblock;
              temp2 = psociTime();
              block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj, sefsef ); // sefsefallocated internally
              sociTime += psociTime() - temp2;

// Need to work through the zweros here

              int numTest=0;

if ( block_status == 0 ) ++emptyfetches;

              if ( block_status )
                { // return 0 or 1 where 1 is GOOD: found something
                  temp3 = psociTime();
                  numTest = packSefsef( fetchi, fetchj , sefsef, cimat, icicol, number ); // push_back onto cimat
                  totalElems += numTest;
                  packTime += psociTime() - temp3;
		} // block status
               //cout << "numTest at jconf " << block_status <<" "<< jconf<<" "<<numTest<<endl;
             } // process blocks
	    } // jconf do all js for a given iconf

            if ( totalElems > 0 ) { // this should be uncommon as it is for the entire jconf
                  //cout << "MVP " << cimat.size()<<" "<<icicol.size()<<" "<<number.size() << endl;
                  temp2 = psociTime();
                  estimate_elements += matrixVectorUpdate(iconf, indexnewarray, fetchi.nsefi, bufVnumber, bufVIndex, bufV2, icicol, cimat, number, testVector );
                  updateTime += psociTime() - temp2;
            }

                  //const int maxelements = 1024*1024; // limits scatter stack to something like 16 MB
             //if ( (indexnewarray >= 5) || estimate_elements > maxelements || indexnewarray >= l_maxsefs-1) {
             if ( totalElems > 0  ) {
             //cout << "scatter " << totalElems<<" "<<indexnewarray << endl;
             temp2 = psociTime();
             hvcount += matrixProductScatter( iroot, indexnewarray, estimate_elements, bufVnumber, bufVIndex, bufV2, indexarray, dataarray, g_Hv );
             scatterTime += psociTime() - temp2;
             }

	  nexttask = jobs.nxtask( g_size );

	} // nexttask
      } // iconf
      
      jobs.destroyNxtask(); //Close'er up

      cout << "end of root didblock phaseblock excitblock totalblock is " << didblock<<" phaseblock " << phaseblock<< " excit block " << excitblock<< " totalblock is " << totalblock<<endl;
      
     } // iroot
      
      if ( GA::nodeid() == 0 ) {
	cout << GA::nodeid() << " sociTime is " << sociTime<<endl;
	cout << GA::nodeid() << " assemTime is " << assemTime<<endl;
	cout << GA::nodeid() << " compare is " << compareTime<<endl;
	cout << GA::nodeid() << " computeConfsTime is " << computeConfsTime << endl;
	cout << GA::nodeid() << " updateTime is " << updateTime << endl;
        cout << GA::nodeid() << " packTime is " << packTime << endl;
	cout << GA::nodeid() << " scatterTime is " << scatterTime << endl;
	cout << "total matrix Fetches is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
        cout << "Empty fetches was " << emptyfetches << endl;
      }
      
#ifdef SUPPLEMENTAL
      //      destroySuppVal();
#endif
      
      GA::SERVICES.sync();
      
      return(0);
 }

/* BUild the M*V direct BUT fetch dets instead of compute them
*/
int PsociGAhamiltonian::matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulateDirectFetch( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
								      int * hv_dims, int * v_lo, int * v_hi )
 {
#ifndef SUPPLEMENTAL
   GA::error(" Must compile code with -DSUPPLEMENTAL ",-1);
#endif
#ifdef SQUAREHAMILTONIAN
   GA::error(" Must NOT compile with -DSQUAREHAMILTONIAN",1);
#endif
#ifdef SCATNOFLAT
GA::error(" SCATNOFLAT not supported in matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate",1);
#endif
#ifdef LOOPNPUT
GA::error(" LOOPNPUT not supported in matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate",1);
#endif
#ifdef USEMPI
GA::error(" USEMPI not supported in matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate",1);
#endif

  int g_rank = GA::nodeid();
  int g_size = GA::nodes();

  if ( g_rank == 0 ) cout << "construct DIRECT matrix vector products" << endl;

  char handleMessage[] = "matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulateDirect: current CI vectors: g_v";
  g_v->checkHandle( handleMessage );
  
  double temp, temp2, temp3;
  double packTime=0.0, scatterTime=0.0, sociTime=0.0, updateTime =0.0, compareTime=0.0, computeConfsTime=0.0 ;
  double assemTime=0.0;
  int totalFetches=0;
  long allcount=0;
  long block_status;

  int hvcount=0;
  
/* bufVnumber, bufV2, bufVIndex are all preallocated to be MAX_AGGREGATE int 
   indexnewarray is used to index it
*/

  pair<int,double> info1;    // Need to loop over the H, fetch it into sefsef blocks and pass that to the below method
  
  //PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
  const int ilo = 0;
  const int ihi = l_maxspatials;
  
  if ( g_rank == 0 ) cout << " RECOMPUTE DIRECT SOCI instead of simply fetch them " << endl;

#ifdef NEWORBMAP
  GA::error(" PsociGAhamiltonian::matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulateDirect must not use -DNEWORBMAP",1);
#endif

  int countrun;

  int l_maxsefs = local_dims[0];
  int maxsefperconf = l_deters->fetchGlobalMaxSefPerSpatial();

  int numscatters=0;
  double alpha=1.0;

#ifdef TACC
    bool transpose=true;
#else
    bool transpose=false;
#endif

/*
    vector<vector<double> > cimat(1);
    vector<vector<int> > icicol(1);
    vector<int> number(1);
*/

// could also check on even odd for screening
      
    int didblock=0;
    int totalblock=0;
    int excitblock=0;
    int phaseblock=0;
    int possibleblock=0;

    vector<double> testVector; // hold maxsefperdet worht of vector data 

//must check parallel runs
// doesn't do anything
    vector<char> * rep_phases = l_deters->fetchReplicatedPhasesHandle();
    cout << "rep phase is size " << rep_phases->size()<<endl; // returns iopen

    int emptyfetches=0;
      
    for(int iroot=start; iroot < end; ++iroot ) {
      countrun=0;
      cout << "IN roots" << endl; 

      PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
      jobs.initNxtask();
      long nexttask = jobs.nxtask( g_size );
      
      int indexnewarray = -1; // indexarray preallocate to be MAX_AGGREGATE in length
      int estimate_elements=0;
      
          vector<double> testVector( l_maxsefs, 0.0);
//        testVector.resize( l_maxsefs );
          fetchAllVector( iroot, testVector, g_v ); // Note the value of vend

//      cout << "use j[hase directly " << endl;

      for(int iconf = ilo; iconf < ihi; ++iconf ) {
        int numElems = 0;

        totalblock += (iconf+1); 
	
        if ( nexttask == iconf ) { // in the loop stride is implemented  

         vector<COMPRESS_SIZE> imap( l_nbf+2  );
         JOUTFG fetchi;

         temp2 = psociTime();
         l_deters->fetchAndUnpackDeterminantData( iconf+1, fetchi );
         l_deters->fetchOrbMap( iconf+1, imap);
         computeConfsTime += psociTime() - temp2;
         ++totalFetches;

// Now first identify all  easily removed abs(iopen-jopen) > 4 entries
// drops list down to 12% for a simple case

         temp2 = psociTime();
         int iphase = rep_phases->at(fetchi.index-1);
         vector<int> jphases;
         int numJs = assemblePhasePossibles(iconf, iphase, jphases, rep_phases );
         phaseblock += jphases.size();
         //cout << "assemblePhasePossibles sinple > 4 test returned " << numJs << " out of a total of " << iconf << endl;
         assemTime += psociTime() - temp2;
        
// Now check two things. 1 for a given comfigure are we more than 2 exciations different AND
// Do we have specific projections that match beyond simply even or odd?
// drops list down to < 1% for a simple case
  
          vector<int> list;

          temp2 = psociTime();
//already fetch above l_deters->fetchOrbMap(iconf+1, imap); // resize imap
//          countrun +=  l_deters->compareAllOrbMapMatrixVectorProducts(iconf+1, imap, list ); // no change to imap : fetch jlist for each i usually 100-300 entries
// this really should only be done ONCE

          countrun +=  l_deters->compareAllOrbMapMatrixVectorProducts(iconf+1, imap, jphases, list );
          jphases.clear();
          excitblock += list.size();
          compareTime += psociTime() - temp2;

          int width = fetchi.nsefi;
          vector<vector<double> > cimat( width );
          vector<vector<int> > icicol( width );
          vector<int> number( width );
/*
          set<int> iprojections;
          for(int k=0; k< fetchi.ndeti; ++k ) {
            iprojections.insert( fetchi.spin_projection[k] ) ;
          }
*/

          int totalElems=0;
          vector<int>::iterator jit;
          for(jit = list.begin(); jit != list.end(); ++jit )
            {
              int jconf = (*jit);
              vector<COMPRESS_SIZE> jmap( l_nbf+2  );
              JOUTFG fetchj;

              temp2 = psociTime();
              l_deters->fetchAndUnpackDeterminantData( jconf+1, fetchj );
              //l_deters->fetchOrbMap( jconf+1, jmap);
              ++totalFetches;
              computeConfsTime += psociTime() - temp2;

/*
              set<int> jprojections;
              for(int k=0; k< fetchj.ndeti; ++k ) {
                jprojections.insert( fetchj.spin_projection[k] );
              }
*/

// Does the new J have ANY specific projections in common with the input iconf?

/*
           //bool process=false;
           for (set<int>::iterator sit=iprojections.begin(); sit != iprojections.end(); ++sit) {
              int temp=(*sit);
              for (set<int>::iterator sjt=jprojections.begin(); sjt != jprojections.end(); ++sjt) {
               cout << "projs "<<temp<<" "<<(*sjt)<<endl;
                  if ( temp == (*sjt) ) {
                  process=true;
                   cout <<"projs inner at " << temp<<" "<<(*sjt)<<endl;
                  goto processBlocks;
                  }
                }
            }
*/
           // cout << "if here no projection match " << iconf<<" "<<jconf<<endl;

           processBlocks: 
           bool process=true;

           if ( process ) {

//cout << "Compute a bocks" << endl;

              ++didblock;
              temp2 = psociTime();
              block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj, sefsef ); // sefsefallocated internally
              sociTime += psociTime() - temp2;

// Need to work through the zweros here

              int numTest=0;

if ( block_status == 0 ) ++emptyfetches;

              if ( block_status )
                { // return 0 or 1 where 1 is GOOD: found something
                  temp3 = psociTime();
                  numTest = packSefsef( fetchi, fetchj , sefsef, cimat, icicol, number ); // push_back onto cimat
                  totalElems += numTest;
                  packTime += psociTime() - temp3;
		} // block status
               //cout << "numTest at jconf " << block_status <<" "<< jconf<<" "<<numTest<<endl;
             } // process blocks
	    } // jconf do all js for a given iconf

            if ( totalElems > 0 ) { // this should be uncommon as it is for the entire jconf
                  //cout << "MVP " << cimat.size()<<" "<<icicol.size()<<" "<<number.size() << endl;
                  temp2 = psociTime();
                  estimate_elements += matrixVectorUpdate(iconf, indexnewarray, fetchi.nsefi, bufVnumber, bufVIndex, bufV2, icicol, cimat, number, testVector );
                  updateTime += psociTime() - temp2;
            }

                  //const int maxelements = 1024*1024; // limits scatter stack to something like 16 MB
             //if ( (indexnewarray >= 5) || estimate_elements > maxelements || indexnewarray >= l_maxsefs-1) {
             if ( totalElems > 0  ) {
             //cout << "scatter " << totalElems<<" "<<indexnewarray << endl;
             temp2 = psociTime();
             hvcount += matrixProductScatter( iroot, indexnewarray, estimate_elements, bufVnumber, bufVIndex, bufV2, indexarray, dataarray, g_Hv );
             scatterTime += psociTime() - temp2;
             }

	  nexttask = jobs.nxtask( g_size );

	} // nexttask
      } // iconf
      
      jobs.destroyNxtask(); //Close'er up

      cout << "end of root didblock phaseblock excitblock totalblock is " << didblock<<" phaseblock " << phaseblock<< " excit block " << excitblock<< " totalblock is " << totalblock<<endl;
      
     } // iroot
      
      if ( GA::nodeid() == 0 ) {
	cout << GA::nodeid() << " sociTime is " << sociTime<<endl;
	cout << GA::nodeid() << " assemTime is " << assemTime<<endl;
	cout << GA::nodeid() << " compare is " << compareTime<<endl;
	cout << GA::nodeid() << " computeConfsTime is " << computeConfsTime << endl;
	cout << GA::nodeid() << " updateTime is " << updateTime << endl;
        cout << GA::nodeid() << " packTime is " << packTime << endl;
	cout << GA::nodeid() << " scatterTime is " << scatterTime << endl;
	cout << "total matrix Fetches is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
        cout << "Empty fetches was " << emptyfetches << endl;
      }
      
#ifdef SUPPLEMENTAL
      //      destroySuppVal();
#endif
      
      GA::SERVICES.sync();
      
      return(0);
 }
 

//int PsociGAhamiltonian::packDataNoZeros(int ivec, vector<int> & number,  vector<vector<int> >& index, vector<vector<double> >& buf, vector<int> & indexarray, vector<double> & dataarray )

/* this seems to have some mistake */
int PsociGAhamiltonian::matrixVectorProductsSplitLocalVectorRepAccumulate( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
								      int * hv_dims, int * v_lo, int * v_hi )
 {
#ifndef SUPPLEMENTAL
   GA::error(" Must compile code with -DSUPPLEMENTAL ",-1);
#endif
#ifdef SQUAREHAMILTONIAN
   GA::error(" Must NOT compile with -DSQUAREHAMILTONIAN",1);
#endif
   
   // Check basic array dimensions
   
   int type;
   int dims[ 2 ]; //used by all
   int idim;
   int numcores = GA::nodes();
   
   int l_maxsefs = local_dims[0];
   
   char op[] = "+";
   
   /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
      within GA.  In the future this will be relaxed.
      Check them out here
   */
   
   /* cimat and friends are distributed over isefs but no breaking along sparity axis */
  
  
  int lo[2], hi[2];
  
  int cilo[2], cihi[2];
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  
  cilo[0] = clow;
  cilo[1] = local_cilo[1];
  cihi[0] = chi;
  cihi[1] = local_cihi[1];
  
  // these are the local vector piece
  
  double temp;
  double tempTime= psociTime();
  double tempTime2, getTime = 0.0, accTime=0.0, brdTime=0.0;
  double temp2;
  
  
  // Loop over all sefs, fetching LOCAL cimat for processing.
  // Gather Scatter objects can never be bigger than this
  

  int chunk_width = v_hi[0] - v_lo[0] + 1;
  if ( chunk_width <= 0 ) {
     GA::error(" chunk_width <= 0", chunk_width ); 
  }

  
  vector<double> testVector( chunk_width , 0.0 ); // need to worry about this it will eventually segfault 
  
  int hvcount=0;
  int n=1;
  int number; // how long was the relevant CI row ?
  double dOne = 1.0;
  double alpha = 1.0;
  double  listTime=0.0,compTime1=0.0, compTime2=0.0, scatTime=0.0, getVTime=0.0;
  
  double temptime, cimatTime=0.0;
 
  int width = MAX_AGGREGATE;  // pertain to H prefeting - here we only need a value of 1
  
  //if ( my_chunklist.size() != 1 ) GA::error("chunklist must be 1 " ,1);
  
  vector<double> bufV(max(l_maxsefs,l_maxwidth) ); // Is this faster or slower than simply zeroing out?
  
  /* For bufV upload after GOP*/
  lo[0] = v_lo[0];
  hi[0] = v_hi[0];
  lo[1] = 0;
  hi[1] = 0;
  int index = v_lo[0]; // For bufV 
  int vec_index = 0;

  int ch_end = v_hi[0]; // Note do not subtract one from here
  int ch_start = v_lo[0];

  for(int ivec=start; ivec< end; ++ivec ) 
    {
      lo[1] = ivec;
      hi[1] = ivec;
      int icountloops = 0;
      
      zeroVector( bufV );
//      zeroVector( testVector );
	
	temp = psociTime();

// Everybody gets their local piece of V

        for(int icore=0; icore < numcores; ++icore ) {

        ++icountloops;

        zeroVector( testVector );
        ch_start = -1;
        ch_end = -1;

	if ( GA::nodeid() == icore ) {
            g_v->get( lo, hi, &testVector[0], &n );
            ch_end = v_hi[0]; // Note do not subtract one from here
            ch_start = v_lo[0];
        }
//maybe put a resize function here

    //    cout << GA::nodeid() <<" before the brcsts" << endl;
        GA::brdcst( &ch_start, sizeof(int), icore);
        GA::brdcst( &ch_end, sizeof(int), icore );
        testVector.resize( ch_end - ch_start + 1 ) ;
        GA::brdcst( &testVector[0], testVector.size()*sizeof(double), icore);
     ///   cout << GA::nodeid() << " vector get " << lo[0]<<" "<<hi[0]<<" "<<ch_start<<" "<<ch_end<<endl;

	getVTime += psociTime() - temp;
	
	int jhi;
	int buf_index=-1;
	
#ifdef MATRIXBALANCE
	int fetch_cimat_count=-1;
	
	// for(int it=0; it < cimat_list.size(); ++it ) {
	
	for(int it=0; it < cimat_list.size(); it+= MAX_AGGREGATE ){ // a factor of 5 in time variation
	  int mxnum = cimat_list.size();
	  int max_num = MAX_AGGREGATE;
	  if ( ( mxnum - it ) < MAX_AGGREGATE ) max_num = mxnum - it;
	  
	  vector<int> newlist;
          temptime=psociTime();
	  for(int j=0; j< max_num; ++j) {
	    newlist.push_back( cimat_list[ it+j ]);
	  }
          fetchsubListHamiltonianData( newlist ); // cimat still populate starting at 0 
          buf_index=-1;
          cimatTime += psociTime() - temptime;
	  
	  for(int it2=0; it2< max_num; ++it2) { // extra loop
	    int i = newlist[ it2 ]; 
#else
	    for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
#endif
	      ++buf_index;
	      
#ifndef MATRIXBALANCE
              // local data is grabbed in this case
	      buf_index=0; //override and simply grab local i
	      fetchLocalHblock( i, local_cimat[buf_index], local_icicol[buf_index], local_number[buf_index] ); // It actually does not need to be local but beware performance
#endif
	      int start_rng=local_icicol[buf_index][0];
	      int end_rng=local_icicol[buf_index][local_number[buf_index]-1];
	      number = local_number[buf_index];

	      double dtemp=0.0;
	      if (!ch_end < start_rng  || !ch_start > end_rng ) { //none within the range
                //cout << GA::nodeid() << " MAPS " << jhi <<" " <<  ch_start<<" "<<ch_end<< endl;
		temp = psociTime();
		for(int j=0; j< number; ++j ) {
		  jhi = local_icicol[buf_index][j];
		  if ( jhi >= ch_start && jhi <= ch_end ) {
               //        cout << GA::nodeid() << " MAPS " << jhi <<" " <<  ch_start<<" "<<ch_end<< endl;
		  dtemp += local_cimat[buf_index][ j ] * testVector[ jhi - ch_start]; //globally accesses testVector
		  }
		}
                bufV[i] += dtemp;
		compTime1+=psociTime() - temp;
	      }
	      if ( i >= ch_start && i <= ch_end ) {
                temp = psociTime();
		double testVec = testVector[ i - ch_start ]; 
		for(int j=0; j< number-1; ++j ) {
		  int jhi = local_icicol[buf_index][j];
		  bufV[ jhi ] += local_cimat[buf_index][ j ] * testVec;
		}
		compTime2+=psociTime() - temp;
	      } // ch_start screen
#ifdef MATRIXBALANCE
	    } // extra loop: MAX_AGGREGATE blocks	
#endif
	  } // clow
        } // numcore experiment
	
	tempTime2 = psociTime();
	n=1;
	lo[1] = ivec;
	hi[1] = ivec;
	GA::gop( &bufV[0], l_maxsefs, op);
	g_Hv->put( lo, hi, &bufV[index], &n);
	accTime += psociTime() - tempTime2;
	
	++hvcount;

  //    cout << "count loops " << icountloops << " " << endl;

      }  // ivec 
      
      if ( GA::nodeid() == 0 ) cout << GA::nodeid() << "cimatTime getV accTime compTime1 comptTime2"<<cimatTime<<" "<<  getVTime<<" "<<accTime << " "<<compTime1<<" "<<compTime2<<endl;
      
      GA::SERVICES.sync();
      return(hvcount);
    }
  
/* experimental new version using chunked gathers instead of gets 
*/

int PsociGAhamiltonian::matrixVectorProductsIncoreGatherScatterNoVectorGet(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int hvcount = matrixVectorProductsIncoreGatherScatterNoVectorGet(start, end, g_v, v_dims, g_Hv, hv_dims, v_lo, v_hi );
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return(hvcount);
}

/* NON SPLIT versionof the GATHER SCATTER methods
   compile code as -DGATHERHV
*/
int PsociGAhamiltonian::matrixVectorProductsIncoreGatherScatterNoVectorGet( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
								 int * hv_dims, int * v_lo, int * v_hi )
{
#ifndef SUPPLEMENTAL
  GA::error(" Must compile code with -DSUPPLEMENTAL ",-1);
#endif
#ifdef SQUAREHAMILTONIAN
  GA::error(" Must NOT compile with -DSQUAREHAMILTONIAN",1);
#endif

  // Check basic array dimensions
  
  int type;
  int dims[ 2 ]; //used by all
  int idim;
  
  int l_maxsefs = local_dims[0];
  
  char op[] = "+";
  long numElems = l_maxsefs;
  
  /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
     within GA.  In the future this will be relaxed.
     Check them out here
  */
  
  /* cimat and friends are distributed over isefs but no breaking along sparity axis */
  
  
  int lo[2], hi[2];

/*
  int cilo[2], cihi[2];
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  cilo[0] = clow;
  cilo[1] = local_cilo[1];
  cihi[0] = chi;
  cihi[1] = local_cihi[1];
*/
  
  // these are the local vector piece
  
  double temp;
  double tempTime= psociTime();
  double tempTime2, getTime = 0.0, accTime=0.0, brdTime=0.0;
  double temp2, bufTime=0.0;
  double mapTime=0.0;
  double pairTime=0.0;
  double minimapTime=0.0;
  
  double vectorTime=0.0;
  double vecSort=0.0;

  double temp4,allTime=0.0;
  
  // Loop over all sefs, fetching LOCAL cimat for processing.
  // Gather Scatter objects can never be bigger than this
  
 // vector<double> bufV( max(l_maxsparse,l_maxsparse_supp), 0.0  );
  
/*
  const int nchunk = l_nchunk;
  int chunk_width;
  int ch_start=0;
  int ch_end=0;
  (l_maxsefs % nchunk==0)? chunk_width = l_maxsefs / nchunk: chunk_width = 1 + (l_maxsefs / nchunk);
  if ( GA::nodeid() == 0 ) cout << "chunkliness is  " <<MAX_AGGREGATE<<" "<<MAX_DEFLATE_ELEMENTS<<" "<<MINIMAP_MAX<<endl;
*/

//  vector<double> testVector( chunk_width , 0.0 );
  
  int hvcount=0;
  int n=1;
  int number; // how long was the relevant CI row ?
  int fullnumber=0;
  double dOne = 1.0;
  
  double alpha = 1.0;
  double  listTime=0.0,compTime1=0.0, compTime2=0.0, scatTime=0.0, getVTime=0.0, deflateTime=0.0;
  double deflateListTime=0.0, deflateVecTime=0.0;
  
  double temptime, cimatTime=0.0;
 // int chunk_lo[2];
 // int chunk_hi[2];
  
  int newnumber=0;
  int numdeflates;
  
// Results data
  int width = MAX_AGGREGATE; 


//  vector<double> dataarray( width * max(l_maxsparse,l_maxsparse_supp), 0.0  ); // memory registration too expensive to update in deflation
//  vector<int> indexarray( width * 2*max(l_maxsparse,l_maxsparse_supp),0 ); // maximum possible 
//  vector<vector<double> > vectors( width, vector<double>( max(l_maxsparse,l_maxsparse_supp) ) );
  
  
  int deflate_index;
  int partialnumber;
  
  int nummapcontracts=0; // number of time minimap was copied to buffer_map

  //local_icicol.resize( width, vector<int>( max(l_maxsparse,l_maxsparse_supp), 0 ) );
  //local_cimat.resize( width, vector<double>( max(l_maxsparse,l_maxsparse_supp), 0.0 ));
  //local_number.resize( width, 0 );
  
#ifndef MATRIXBALANCE
  GA::error(" Must define -DMATRIXBALANCE to use :matrixVectorProductsIncoreGatherScatterNoVectorGet",1); 
#endif

#ifdef SQUAREHAMILTONIAN
  GA::error(" Must NOT define -DSQUAREHAMILTONIAN to use :matrixVectorProductsIncoreGatherScatterNoVectorGet",1);
#endif

// multimaps MUST NOT SPAN multiple ivecs

  temp4=psociTime();
  for(int ivec=start; ivec< end; ++ivec ) 
    {
      multimap<int, double > buffer_map;
      multimap<int, double > mini_map;
      
      numdeflates=0;
      vector<int>::iterator it;
	
	// Hamiltonian reduction
	
	int jhi;
	
	deflate_index = 0; // how many deflations/scatters occured
	
	int buf_index=-1;
        int max_deflate_elem=0;
        int minimap = 0;
	
	int fetch_cimat_count=-1;
	
// remember MAX_AGGREGATE is basically rows at a time and blocks over gather calls and hamiltonian fetch calls

	for(int it=0; it < cimat_list.size(); it+= MAX_AGGREGATE ){ 

	  int mxnum = cimat_list.size();
	  int max_num = MAX_AGGREGATE;
	  if ( ( mxnum - it ) < MAX_AGGREGATE ) max_num = mxnum - it;
	  
	  vector<int> newlist;
          temptime=psociTime();

	  for(int j=0; j< max_num; ++j) {
	    newlist.push_back( cimat_list[ it+j ]);
	  }
//-> H
          fetchsubListHamiltonianData( newlist ); // cimat still populate starting at 0 
          buf_index=-1;
	  minimap = 0; 
          cimatTime += psociTime() - temptime;
//-> V
          temptime=psociTime();
          int num = gatherVectorSet(g_v, ivec, newlist, vectors, dataarray, indexarray);
          vectorTime += psociTime() - temptime;

/* fetch set of vectors as well */

          int vecindex = -1;

	  for(int it2=0; it2< max_num; ++it2) { // extra loop
	    int i = newlist[ it2 ]; 
	      ++buf_index; // niot needed use it2
	      
	      number = local_number[buf_index];
	      zeroVector( bufV, number );
	      
	      double dtemp=0.0;
		
               // ++minimap; // grabs full size later
		temp = psociTime();
		for(int j=0; j< number; ++j ) {
                  ++vecindex;
		  jhi = local_icicol[buf_index][j];
		    //dtemp += local_cimat[buf_index][ j ] * vectors[buf_index][j]; //globally accesses testVector
                    dtemp += local_cimat[buf_index][ j ] * vectors[vecindex];
		}

                temp2=psociTime();
		mini_map.insert( pair<int, double>( local_icicol[buf_index][number-1], dtemp ) );
		minimapTime+= psociTime() - temp2;
	      
	      // Need to minimize the length of the multimap otherwise NlogN bites you
	      
	      compTime1 += psociTime() - temp;
	      temp = psociTime();
	      
		//double testVec = vectors[buf_index][ number-1 ]; // only works for triangular H 
                double testVec = vectors[vecindex]; // same as ending of last loop 
		
		for(int j=0; j< number-1; ++j ) {
		  int jhi = local_icicol[buf_index][j];
		  temp2=psociTime();
		  bufV[ j ] += local_cimat[buf_index][ j ] * testVec;
		  bufTime += psociTime()-temp2;
		  
		  if ( abs(bufV[j]) > MIN_HVEC_TOKEEP ) { 
		    temp2=psociTime();
		    mini_map.insert(pair<int, double>( jhi, bufV[ j ] ) );
                    minimap = mini_map.size(); 
		    minimapTime += psociTime()-temp2;
		  }
		}
	      
	      if ( minimap >= MINIMAP_MAX ) { //other cases get handled below
		++nummapcontracts;
		temp2 = psociTime();
		buffer_map.insert( mini_map.begin(), mini_map.end() );
                max_deflate_elem += mini_map.size();
              //  dump_count += mini_map.size();
		mini_map.clear();
		mapTime+= psociTime() - temp2;
		minimap = 0;
	      }
	      
	      compTime2 += psociTime() - temp;
	      
	      if ( max_deflate_elem >= MAX_DEFLATE_ELEMENTS || i >= cimat_list[ cimat_list.size() - 1 ]  ) {

		  temp = psociTime();
		  if ( mini_map.size() >= 1 ) {
		    buffer_map.insert( mini_map.begin(), mini_map.end() ); // ensure no strays
		    ++nummapcontracts;
                    minimap=0;
		    mini_map.clear();
		  }
		  mapTime+= psociTime() - temp;
		  
// data and indexarray get reused

		  temp=psociTime();
		  newnumber = deflateMultimapToVector( buffer_map, dataarray, indexarray, ivec ); 
                  buffer_map.clear();
		  deflateTime +=psociTime()-temp;

		  if ( newnumber > 0 ) {
//cout << "do scattr with " << newnumber << "size of index is " << indexarray.size() << endl;
		    temp2=psociTime();
		    NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&dataarray[0], (int*)&indexarray[0], newnumber, &alpha);
		    deflate_index=0;
		    scatTime += psociTime() - temp2;
		    ++numdeflates;
		  } 
		  //dump_count = 0;
                  max_deflate_elem=0;
                  
		  
		} // dump count
	      } // extra loop: MAX_AGGREGATE blocks	
	    } // clow
	  ++hvcount;
	}  // ivec 
	
        allTime = psociTime()- temp4;
	cout << GA::nodeid() << "alltim  AGGREGATE ciTime vectorT comp1 comp2 bufV mapT minimapT deflate SUM  mapcontracts num deflates"<< allTime<<" "<<MAX_AGGREGATE <<" "<< cimatTime<<" "<<vectorTime<<" "<<compTime1 <<" "<<compTime2 << " "<<bufTime<<" " <<mapTime <<" "<< minimapTime << " " <<deflateTime<<" "<<cimatTime+getVTime+mapTime+minimapTime+pairTime+deflateTime<<" "<<nummapcontracts<<" " << numdeflates<< endl;
	
	GA::SERVICES.sync();
	return(hvcount);
      }
      
      
/* The SPLIT version of this method
*/
int PsociGAhamiltonian::matrixVectorProductsIncoreGOPSplit( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
						    int * hv_dims, int * v_lo, int * v_hi )
{
#ifdef DETAILEDHAMCHECK
  char * vectorMessage = "PsociGAhamiltonian:matrixVectorProductsIncore:Vector";
  g_v->checkHandle( vectorMessage );
  
  char * HvectorMessage = "PsociGAhamiltonian:matrixVectorProductsIncore:H-Vector";
  g_Hv->checkHandle( HvectorMessage );
  
  char * cimatMessage = "PsociGAhamiltonian:matrixVectorProductsIncore:cimat";
  l_g_cimat->checkHandle( cimatMessage );
#endif
  
#ifndef SUPPLEMENTAL
  GA::error(" Must compile code with -DSUPPLEMENTAL ",-1);
#endif

  // Check basic array dimensions
  
  int type;
  int dims[ 2 ]; //used by all
  int idim;
  
  int l_maxsefs = local_dims[0];

  char op[] = "+";
  long numElems = l_maxsefs;
  
  int v_sefs = *v_dims;
  int maxroots = *v_dims+1;
  
#ifdef DETAILEDHAMCHECK
  if ( l_maxsefs != v_sefs ) {
    cerr << "Nonconforming cimat and vector array " << l_maxsefs << " " << v_sefs << endl;
    GA::error( "Nonconforming cimat and vector array ", -1);
  } 
#endif

  int Hv_sefs = *hv_dims;
  int Hvmaxroots = *hv_dims+1;
  
#ifdef DETAILEDHAMCHECK
  if ( l_maxsefs != Hv_sefs ) {
    cerr << "Nonconforming cimat and Hvector array " << l_maxsefs << " " << Hv_sefs << endl;
    GA::error( "Nonconforming cimat and Hvector array ", -1);
  } 
#endif
  
  /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
     within GA.  In the future this will be relaxed.
     Check them out here
  */
  
/* cimat and friends are distributed over isefs but no breaking along sparity axis */

  int lo[2], hi[2];
  int cilo[2], cihi[2];
  
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  cilo[0] = local_cilo[0];
  cilo[1] = local_cilo[1];
  cihi[0] = local_cihi[0];
  cihi[1] = local_cihi[1];
  
// these are the local vector piece

  int vlow = *v_lo;
  int vhi  = *v_hi;
  int hvlow = vlow;
  int hvhi  = vhi;

  double tempTime= psociTime();
  double tempTime2, bufTime=0.0, getTime = 0.0, accTime=0.0, brdTime=0.0;
  
#ifdef DETAILEDHAMCHECK
  if ( clow != vlow || clow != hvlow || vlow != hvlow ) {
    cout << "WARNING: LOs are not equal" << endl;
    //    GA::error("LOs are wrong ",-1);
  }
  if ( chi != vhi || chi != hvhi || vhi != hvhi ) {
    cout << "WARNING: HIs are not equal " << endl;
    //GA::error("HIs are wrong ",-1);
  }
#endif
  
  // Okay we are equivelently distributed accross the cores.
  // Now we can use a data-parallel update scheme which is very fast
  // Start the processing Only process one H*v at a time for now
  
#ifdef DETAILEDHAMCHECK
  if ( start < 0 || end <= start || end > maxroots ) {
    cerr << "absurd value for start or end " << start << " " << end << endl;
    GA::error("absurd value for start or end ", end-start );
  }
#endif
  
  // Loop over all sefs, fetching LOCAL cimat for processing.
  
  lo[0] = v_lo[0];
  hi[0] = v_hi[0];
  lo[1] = 0; 
  hi[1] = 0;
  
  int index = v_lo[0]; // For bufV and testvector
  int width = v_hi[0] - v_lo[0] + 1; // width of local vector read/put

 // cout << GA::nodeid() << " l_maxsparse " << l_maxsparse<<" l_maxsparse_supp " << l_maxsparse_supp<<endl;
  vector<int> icicol( max(l_maxsparse,l_maxsparse_supp), 0 );
  vector<double> cimat( max(l_maxsparse,l_maxsparse_supp), 0.0 );
  
  int hvcount=0;
  int n;
  int number; // how long was the relevant CI row ?

  // cout << GA::nodeid() << "MATRIX VECTOR INDEX " << lo[0]<<" "<<lo[1]<<" "<<hi[0]<<" "<<hi[1]<<endl;

  vector<double> bufV(max(l_maxsefs,l_maxwidth) ); // Is this faster or slower than simply zeroing out?
  vector<double> testVector(max(l_maxsefs,l_maxwidth) );// Vector buffer for individual H*v

// Data parallel approach with CIMAT blocked

  double tempTime3, cimatTime=0.0;
  for(int ivec=start; ivec< end; ++ivec ) 
    {
      lo[1]=ivec;
      hi[1]=ivec;
      
      
      tempTime2=psociTime();
      zeroVector( bufV );
      zeroVector( testVector );

      n=1;
      g_v->get( lo, hi, &testVector[index], &n ); // we are only fetching the local piece into a full size space
      GA::gop( &testVector[0], numElems, op );
      getTime += psociTime() - tempTime2;

      // Begin processing compute the H*v term and shove result into Hv
      /* In this scenario it is possible that the first level storage points
         to the split storage. We will need to re-factor this but as a test 
         this should still work
      */

      int jhi;
      for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
	cilo[0] = i;
	cihi[0] = i; // This may not apply to the supp space: need a check at creation
	cilo[1] = 0;
	cihi[1] = 0;

       tempTime3=psociTime();
       fetchLocalHblock( i, cimat, icicol, number ); // It actually does not need to be local but beware performance
       cimatTime+= psociTime() - tempTime3;

        int jhi;
        tempTime2=psociTime();
        double testVec = testVector[ i ];
        for(int j=0; j< number; ++j ) {
          jhi = icicol[ j ];
           bufV[i] += cimat[ j ] * testVector[ jhi ];
        }
        for(int j=0; j< number-1; ++j ) {
          jhi = icicol[j];
          bufV[ jhi ] += cimat[ j ] * testVec; // globally updates bufV
        }
        bufTime += psociTime()-tempTime2;
       }

/* New method to test */
 
      // Push back to ivec in H*v 
      /* At some point we may want to use deltaHs instead of Hs but for now this should be okay.
	 Using the acc method, however, forces you to specify any preexisting H*v terms using
	 PsociGAbasis::setpresent_matrixProducts( numPrecomputedRoots ); 
	 
	 DISABLED: Might need to revisit this for -DCONFIGSCRAMBLE: Be careful about preexisting H*v will be accumulated causing errors.
	 We MUST accumulate here  
      */

      tempTime2 = psociTime();
      n=1;
      GA::gop( &bufV[0], numElems, op);
      g_Hv->put( lo, hi, &bufV[ index ], &n);
      accTime += psociTime() - tempTime2;
      
      ++hvcount;

      double outerTime = psociTime() - tempTime;
      cout << "GOP-OUTER split" << outerTime<<" cimat " << cimatTime << " ACCUM "<< accTime<<" GET " << getTime<< " buftime "<<bufTime<<endl;

    } // ivec 
  return(hvcount);
}

int PsociGAhamiltonian::matrixVectorProductsIncoreGOPSplitChunk(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int hvcount = matrixVectorProductsIncoreGOPSplitChunk(start, end, g_v, v_dims, g_Hv, hv_dims, v_lo, v_hi );
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return(hvcount);
}

/* The SPLIT version of this method
   NOTE is selecting SQUAREHAMILTONIAN then always do chunking
   nchunk is a caller specified chunking ofr vec ands bufV

   This version may only be accessed via -DSQUAREHAMILTONIAN
   Performance note: The time for a H*v grows linearly ( or slightly more) with nchunk

   MUST use square hamiltonian

*/
int PsociGAhamiltonian::matrixVectorProductsIncoreGOPSplitChunk( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
						    int * hv_dims, int * v_lo, int * v_hi )
{
#ifdef DETAILEDHAMCHECK
  char * vectorMessage = "PsociGAhamiltonian:matrixVectorProductsIncore:Vector";
  g_v->checkHandle( vectorMessage );
  
  char * HvectorMessage = "PsociGAhamiltonian:matrixVectorProductsIncore:H-Vector";
  g_Hv->checkHandle( HvectorMessage );
  
  char * cimatMessage = "PsociGAhamiltonian:matrixVectorProductsIncore:cimat";
  l_g_cimat->checkHandle( cimatMessage );
#endif

#ifndef SUPPLEMENTAL
  GA::error(" Must compile code with -DSUPPLEMENTAL ",1);
#endif
  
#ifndef SQUAREHAMILTONIAN
  GA::error(" Must compile code with -DSQUAREHAMILTONIAN ",1);
#endif

  // Check basic array dimensions
  
  int l_maxsefs = local_dims[0];
  
  // long numElems = l_maxsefs;
  
  
#ifdef DETAILEDHAMCHECK
  if ( l_maxsefs != v_sefs ) {
    cerr << "Nonconforming cimat and vector array " << l_maxsefs << " " << v_sefs << endl;
    GA::error( "Nonconforming cimat and vector array ", -1);
  } 
#endif
  
  //int Hv_sefs = *hv_dims;
  //int Hvmaxroots = *hv_dims+1;
  
#ifdef DETAILEDHAMCHECK
  if ( l_maxsefs != Hv_sefs ) {
    cerr << "Nonconforming cimat and Hvector array " << l_maxsefs << " " << Hv_sefs << endl;
    GA::error( "Nonconforming cimat and Hvector array ", -1);
  } 
#endif
  
  /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
     within GA.  In the future this will be relaxed.
     Check them out here
  */

  int buf_lo[2], buf_hi[2];
  //int cilo[2], cihi[2];
  
  int clow = local_cilo[0];
  int chi  = local_cihi[0];

/*
  cilo[0] = local_cilo[0];
  cilo[1] = local_cilo[1];
  cihi[0] = local_cihi[0];
  cihi[1] = local_cihi[1];
*/

  double tempTime= psociTime();
  double getTime = 0.0, accTime=0.0;
  double compTime = 0.0, outerTime = 0.0;
  
#ifdef DETAILEDHAMCHECK
  if ( clow != vlow || clow != hvlow || vlow != hvlow ) {
    cout << "WARNING: LOs are not equal" << endl;
    //    GA::error("LOs are wrong ",-1);
  }
  if ( chi != vhi || chi != hvhi || vhi != hvhi ) {
    cout << "WARNING: HIs are not equal " << endl;
    //GA::error("HIs are wrong ",-1);
  }
#endif
  
  // Okay we are equivelently distributed accross the cores.
  // Now we can use a data-parallel update scheme which is very fast
  // Start the processing Only process one H*v at a time for now
  
#ifdef DETAILEDHAMCHECK
  if ( start < 0 || end <= start || end > maxroots ) {
    cerr << "absurd value for start or end " << start << " " << end << endl;
    GA::error("absurd value for start or end ", end-start );
  }
#endif
  
// cimat, icicol, a number are now prefetched


  int width = chi - clow + 1;  

// This really shoukld only work well at-scale when width is small

// Maybe check for data?
/*
  vector<vector<int> > local_icicol( width, vector<int>( max(l_maxsparse,l_maxsparse_supp), 0 ) );
  vector<vector<double> > local_cimat( width, vector<double>( max(l_maxsparse,l_maxsparse_supp), 0.0 ));
  /vector<int> local_number( width, 0 );

  vector<int> icicol( max(l_maxsparse,l_maxsparse_supp), 0 );
  vector<double> cimat( max(l_maxsparse,l_maxsparse_supp), 0.0 );
  int number; // how long was the relevant CI row ?
*/
  
  int hvcount=0;
  int n;

  const int nchunk = l_nchunk;
  int chunk_width;
  int ch_start=0;
  int ch_end=0;

  (l_maxsefs % nchunk==0)? chunk_width = l_maxsefs / nchunk: chunk_width = 1 + (l_maxsefs / nchunk);

#ifdef DETAILEDCHECK
  cout << GA::nodeid() << " chunk_width " << chunk_width << "nchunk is " << nchunk << endl;
  cout << GA::nodeid() << " nchunk is " << nchunk << " maxsef is " << l_maxsefs <<endl;
#endif
  
  // Do this once from the claling program
  // vector<int> my_chunklist;
  // generateStaggeredChunkList( GA::nodeid(), nchunk, my_chunklist );

   vector<double> bufV( width ); // as wide as the H loop - varies with parallelism
   vector<double> testVector( chunk_width );// Vector buffer for individual H*v

// Data parallel approach with CIMAT blocked

  int chunkfetch=0;
  
  int chunk_lo[2];
  int chunk_hi[2];
  double temp;

  double fetchTime=0.0;

// We could simply download the entire local H matrix stuff here to save time.
/*
  temp2=psociTime();
  fetchAllLocalHblock( clow, chi, local_cimat, local_icicol, local_number ); 
  fetchTime += psociTime() - temp2;
*/
  int ham_index=0;

  for(int ivec=start; ivec< end; ++ivec ) 
    {
      compTime=0.0;
      accTime=0.0;
      getTime=0.0;
      outerTime=0.0;
      fetchTime=0.0;

      tempTime = psociTime();

      buf_lo[1]=ivec;
      buf_hi[1]=ivec;
      
      zeroVector( bufV );
      
      chunk_lo[1] = ivec;
      chunk_hi[1] = ivec;
      
      ch_start=0;
      ch_end=0;

//scramble order in which nodes request data else we simply all call the same GA nodes at about the same time.

      vector<int>::iterator it;
      chunkfetch = 0;
      for(it=my_chunklist.begin(); it!=my_chunklist.end(); ++it) { // this staggers the way get communications occur.
 
        ch_start = chunk_width*(*it); 
	ch_end = min( ch_start + chunk_width, l_maxsefs );

	chunk_lo[0] = ch_start;
	chunk_hi[0] = ch_end-1;
        //int numelems = ch_end - ch_start;
	n=1;

        temp = psociTime();
#ifdef STAGGERCHUNKLIST
        g_v->get( chunk_lo, chunk_hi, &testVector[0], &n );
#else
        if ( GA::nodeid() == 0 ) g_v->get( chunk_lo, chunk_hi, &testVector[0], &n );
        GA::brdcst( &testVector[0], testVector.size()*sizeof(double), 0);
#endif
        getTime += psociTime() - temp;

	++chunkfetch;

      // Begin processing compute the H*v term and shove result into Hv
      /* In this scenario it is possible that the first level storage points
         to the split storage. We will need to re-factor this but as a test 
         this should still work
      */

//TODO reverse index for cimat, et al.

	int jhi;
	int buf_index=-1;

        temp= psociTime();

	for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
	  ++buf_index;

	  int start_rng=local_icicol[buf_index][0];
	  int end_rng=local_icicol[buf_index][local_number[buf_index]-1];
	  if (!ch_end-1 < start_rng  || !ch_start > end_rng ) {

	    int jhi;
	    double testVec = testVector[ i-ch_start ];
	    
	    for(int j=0; j< local_number[buf_index]; ++j ) {
	      jhi = local_icicol[buf_index][ j ];
	      if ( jhi >= ch_start && jhi < ch_end ) {
		bufV[ buf_index ] += local_cimat[buf_index][ j ] * testVector[ jhi - ch_start ]; //globally accesses testVector
	      }
	    }
	  }
          ++ham_index;
	} // clo chi
        compTime += psociTime() - temp;
      }

      buf_lo[0] = clow;
      buf_hi[0] = chi;
      n=1;
      temp = psociTime();
// not needed     g_Hv->nbPut( buf_lo, buf_hi, &bufV[ 0 ], &n, bufhandle);
      g_Hv->put( buf_lo, buf_hi, &bufV[ 0 ], &n);
      accTime += psociTime() - temp;
      
      ++hvcount;

      outerTime += psociTime() - tempTime;
      if (GA::nodeid() == 0 ) cout << "square chunkfetch is " << chunkfetch << "COMP TIME is " << compTime << " GOP-OUTER " << outerTime<<" ACCUM "<< accTime<<" GET " << getTime<< " FETCHN time is " << fetchTime << endl;
      
    } // ivec 
  
  return(hvcount);
}


int PsociGAhamiltonian::matrixVectorProductsIncoreGOPSplitGatherScatter(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int hvcount = matrixVectorProductsIncoreGOPSplitGatherScatter(start, end, g_v, v_dims, g_Hv, hv_dims, v_lo, v_hi );
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return(hvcount);
}

/* The SPLIT version of this method
   Uses Gather scatter instead of GOPs


   COMPILE CODE with -DGATHERHV -DSUPPLEMENTAL to use this method

*/
int PsociGAhamiltonian::matrixVectorProductsIncoreGOPSplitGatherScatter( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, 
									 int * hv_dims, int * v_lo, int * v_hi )
{
#ifdef DETAILEDHAMCHECK
  char * vectorMessage = "PsociGAhamiltonian:matrixVectorProductsIncoreGatherScatter:Vector";
  g_v->checkHandle( vectorMessage );
  
  char * HvectorMessage = "PsociGAhamiltonian:matrixVectorProductsIncoreGatherScatter:H-Vector";
  g_Hv->checkHandle( HvectorMessage );
  
  char * cimatMessage = "PsociGAhamiltonian:matrixVectorProductsIncoreGatherScatter:cimat";
  l_g_cimat->checkHandle( cimatMessage );
#endif
  
#ifndef SUPPLEMENTAL
  GA::error(" Must compile code with -DSUPPLEMENTAL ",-1);
#endif
  
  // Check basic array dimensions
  
  int l_maxsefs = local_dims[0];
  
#ifdef DETAILEDHAMCHECK
  if ( l_maxsefs != v_sefs ) {
    cerr << "Nonconforming cimat and vector array " << l_maxsefs << " " << v_sefs << endl;
    GA::error( "Nonconforming cimat and vector array ", -1);
  } 
#endif

  
#ifdef DETAILEDHAMCHECK
  if ( l_maxsefs != Hv_sefs ) {
    cerr << "Nonconforming cimat and Hvector array " << l_maxsefs << " " << Hv_sefs << endl;
    GA::error( "Nonconforming cimat and Hvector array ", -1);
  } 
#endif
  
  /* NOTE sefs dimension SHOULD be identically mapped to cores. Doing so results in lots of locality
     within GA.  In the future this will be relaxed.
     Check them out here
  */
  
/* cimat and friends are distributed over isefs but no breaking along sparity axis */

  int cilo[2], cihi[2];
  
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  cilo[0] = clow;
  cilo[1] = local_cilo[1];
  cihi[0] = chi;
  cihi[1] = local_cihi[1];

  // these are the local vector piece
  
  int vlow = *v_lo;
  int vhi  = *v_hi;
  
//  cout << GA::nodeid() << " vlow and vhi are " << v_lo[0]<<" "<<v_hi[0]<<endl;

  double temp;
  double tempTime= psociTime();
  
  // Loop over all sefs, fetching LOCAL cimat for processing.
  
/*
  lo[0] = v_lo[0];
  hi[0] = v_hi[0];
  lo[1] = 0; 
  hi[1] = 0;
*/
  
  vector<int> icicol( max(l_maxsparse,l_maxsparse_supp), 0 );
  vector<double> cimat( max(l_maxsparse,l_maxsparse_supp), 0.0 );
  
  // Gather Scatter objects can never be bigger than this
 
  vector<double> bufV( max(l_maxsparse,l_maxsparse_supp), 0.0  );
  vector<double> testVector( max(l_maxsparse,l_maxsparse_supp), 0.0 );
  vector<int> subsArray(2*max(l_maxsparse,l_maxsparse_supp),0 );
  
  int hvcount=0;
  int n;
  int number; // how long was the relevant CI row ?
  int suppRow;

  double alpha = 1.0;

// Data parallel approach with CIMAT blocked
  
  double scatTime=0.0, getVTime=0.0;
  for(int ivec=start; ivec< end; ++ivec ) 
    {
/*
      lo[1]=ivec;
      hi[1]=ivec;
*/
      n=1;
      
      //int jhi;
      
      for(int i=clow; i<=chi; ++i ) { //data parallel across maxsef
	
	cilo[0] = i;
	cihi[0] = i; // This may not apply to the supp space: need a check at creation
	cilo[1] = 0;
	cihi[1] = 0;
	
	l_g_number->get( cilo, cihi, &number, &n );
        if ( number == 0 ) { //Then this points to an absolute area in supplemental
	  cout << GA::nodeid() << i<< " ZERO ICICOL never supposed to happen " << cilo[0]<<" "<<cilo[1]<<" "<<cihi[0]<<" "<<cihi[1]<<endl;
        }
	
        if ( number >= 0 ) {
#ifdef DETAILEDHAMCHECK
	  if ( number > l_maxsparse ) {
	    cout << "number exceeds specified maxsparse: Will not be able to fetch the main GA. Try supplemental ." << endl;
	    cout << " number = " << number << " maxsparsity is " << l_maxsparse << endl;
	  }    
#endif
	  cihi[1] = number - 1;   
	  n = number;
	  
	  l_g_icicol->get( cilo, cihi, &icicol[ 0 ], &n); 
	  l_g_cimat->get( cilo, cihi, &cimat[ 0 ], &n); 
        
	} else if ( number < 0 ) {
	  suppRow = -number - 1;
	  
	  cilo[0] = suppRow;
	  cihi[0] = suppRow; // This may not apply to the supp space: need a check at creation
	  l_g_number_supp->get( cilo, cihi, &number, &n );
	  
#ifdef DETAILEDHAMCHECK
	  if ( number > l_maxsparse_supp ) {
	    cerr << "number exceeds specified maxsparse_supp: Will not be able to fetch the GA_supp." << endl;
	    cerr << " number = " << number << " maxsparsity_supp is " << l_maxsparse_supp << endl;
	    GA::error("number exceeds specified maxsparse_supp:", l_maxsparse_supp);
	  }
#endif
	  cihi[1] = number - 1;
	  n = number;
	  l_g_icicol_supp->get( cilo, cihi, &icicol[ 0 ], &n);
	  l_g_cimat_supp->get( cilo, cihi, &cimat[ 0 ], &n);
        }
	
	// Gather array indexing one could do simple chunking as well
	
	for (int j=0; j<number; ++j) {
	  int offset = j*2;
	  subsArray[offset] = icicol[j];
	  subsArray[offset+1] = ivec;
	}
        temp=psociTime();
        NGA_Gather_flat(g_v->handle(), (void *)&testVector[0], (int *)&subsArray[0], number);
       getVTime += psociTime() - temp;
        zeroVector( bufV, number ); // need to zero using this method
	
	//NGA_Gather_flat(g_Hv->handle(), (void *)&bufV[0], (int *)&subsArray[0], number);
	
        for(int j=0; j< number; ++j ) {
          bufV[ number-1 ] += cimat[ j ] * testVector[ j ]; //globally accesses testVector

        }
        double testVec = testVector[ number-1 ];
        for(int j=0; j< number-1; ++j ) {
           bufV[ j ] += cimat[ j ] * testVec;
        }
	
        temp=psociTime();
	NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&bufV[0], (int*)&subsArray[0], number, &alpha);
        scatTime += psociTime() - temp;
	//NGA_Scatter_flat(g_Hv->handle(), (void *)&bufV[0], (int*)&subsArray[0], number );
	
      } // clow
      ++hvcount;
      GA::SERVICES.sync();
    }  // ivec 
    cout << "get-2 other scat times " << getVTime<<" "<<scatTime<<" "<<endl;

//  GA::SERVICES.sync(); // TODO CHeck if this is needed
  return(hvcount);
}

// Local call lo to hi INCLUSIVE output in diags
// Local
int PsociGAhamiltonian::fetchListOfDiags( int lo, int hi, vector<double> & diags )
{
  if ( hi < lo || lo < 1 || hi > l_maxsef ) {
    cerr << "Error fetching diags lo and hi are " << lo << " " << hi << endl;
    GA::error("Error fetching diags lo and hi are", -1);
  }
  
  int n=1;
  int dlo[2], dhi[2];
  int num = hi-lo+1;
  if ( diags.size() != num) diags.resize( num );
  
  dlo[0] = lo - 1;
  dhi[0] = hi - 1;
  dlo[1]=0;
  dhi[1]=0;
  
  l_g_diag_sef->get( dlo, dhi, &diags[0], &n );
  
  return( diags.size() );
}

// Collective call  to get ALL diags
int PsociGAhamiltonian::fetchDiags( vector<double> & diags )
{
  int n=1;
  int dlo[2], dhi[2];
  if ( diags.size() != l_maxsef ) diags.resize( l_maxsef );
  
  // Probably faster to have root gather and brdcst
  
  if ( GA::nodeid() == 0 ) {
    dlo[0] = 0;
    dhi[0] = l_maxsef - 1;
    dlo[1] = 0;
    dhi[1] = 0;
    
    l_g_diag_sef->get( dlo, dhi, &diags[0], &n );
  }
  
  GA::brdcst( &diags[0], diags.size()*sizeof(double), 0);
  return( diags.size() );
}


// Local
void PsociGAhamiltonian::printDiags() 
{
  int n=1;
  int dlo[2], dhi[2];
  int num = l_maxsef; 
  vector<double> diags( num );
  
  double temp=0.0;
  
  dlo[0] = 0;
  dhi[0] = num - 1;
  dlo[1]=0;
  dhi[1]=0;
  
  l_g_diag_sef->get( dlo, dhi, &diags[0], &n );
  
  vector<double>::iterator it;
  int i=0;
  for(it = diags.begin(); it!=diags.end(); ++it) {
    cout << "I and diag are " << i << " " << (*it) << endl;;
    temp += (*it);
    ++i;
  }
//#if 0
  cout << " SUM of diags is " << temp << endl;
//#endif
}

// Apply current preconditioner to the current vector, g_v. Transform in-place
// Collective 

//Deprecated - Code no resides in PsociGABasis::

int PsociGAhamiltonian::precondition( double shift, GA::GlobalArray * g_v, vector<double> & diags )
{
  cerr << "This method has been disabled: Refer to new version i PsociGAbasis class " << endl;
  GA::error( "EXIT disabled precond ", 0 );
  return(-1);
  
  /*
    
  char * handleMessage = "PsociGAhamiltonian:precondition";
  g_v->checkHandle( handleMessage );
  
  int lo[2], hi[2];
  g_v->distribution( GA::nodeid(), lo, hi );
  cout << "VECTOR: precondition: distribution is " << lo[0] << " " << hi[0] << endl;
  
  if ( hi[1] != 0 ) {
  cerr << "precondition: g_v has more than 1 existing vector: Must be only 1 " << GA::nodeid() << " " << hi[1] << endl;
  GA::Terminate();
  }
  if ( diags.size() > hi[0] ) {
  cerr << "precondition: diags.size > hi[0]  " << endl;
  GA::error("precondition: diag size:", diags.size());
  }
  
  cerr <<"precondition not used anynmore " << endl;
  GA::Terminate();
  
  //cout << "At the precond " << endl;
  int n = 1; // Need to double check this
  int ilo = lo[0];
  // int ihi = hi[0]; // should be numSefs
  
  if ( hi[0] == -1 ) return( 0 ); // I don't own any of g_v
  num = hi[0] - lo[0] + 1;
  
  vector<double> vec( num, 0 ); 
  g_v->get( lo, hi, &vec[0], &n );
  
  vector<double> vdiags( num, 1); 
  
  int index = ilo;
  double temp;
  for(int i=0; i< num; ++i) {
  temp = diags[index+i]-shift;
  if( temp >= MIN_DENOMINATOR ) vdiags[i] = temp;
  }
  
  for(int i=0; i<=num; ++i) {
  vec[i] /= vdiags[i];
  }
  
  g_v->put( lo, hi, &vec[0], &n ); 
  */

  return( -1 );
  }

/* Begin addition of potentially useful method for managing a supplemental GA array 
   These codes will be used to track the upload of hamiltonian rows that are too long
   for insertion into g_cimat. This supplemental matrix keeps track of the current number of 
   supplemental storage action that were taken.

   Pending statistical analysis of H matrix structure, this could be more or less deep hierarchically
*/

void PsociGAhamiltonian::initSuppVal()
{
  // GA::GlobalArray * g_suppval
  
  int DATA_TYPE = C_LONG; // 
  int DATA_NDIM = 1;
  int dims[1];
  dims[0] = 1;
  
  int chunk[1];
  chunk[0] = -1;
  
  char local_title[] = "supplement array ";
  g_suppval = GA::createGA( DATA_TYPE, DATA_NDIM, dims, (char *) local_title, chunk);
  
  int lo[1];
  lo[0] = 0;
  int hi[1];
  hi[0] = 0;
  int n=1;
  
  int iZero=0;
  g_suppval->put( lo, hi, &iZero, &n);
}

void PsociGAhamiltonian::destroySuppVal()
{
  if ( GA::nodeid() == 0 ) {
    cout << "Total number of Hamiltonian rows was " << l_maxsef << endl; 
    cout << "Sparsity settings were " << l_maxsparse << " and " << l_maxsparse_supp << endl;
    cout << "Total number of supplemental rows required was " << nextSupplimentIrow() << endl;
  }
  g_suppval->destroy();
}


long PsociGAhamiltonian::nextSupplimentIrow()
{
  static int icount = 0;
  int sub = 0;
  const int inc = 1;
  icount = g_suppval->readInc( &sub, inc );
  return( icount );
}



// Also populate the g_diag_sef terms for use by the preconditioner Diag[isef][sef]= value

// Local
int PsociGAhamiltonian::uploadCIMATtoGASupplemental( int isef , vector<vector<double> > &cimat, vector<vector<int> > &icicol, vector<int> & number )
{
  // I think we should COMBINE all [ii] terms into a single row. Then upload the whole thing.
  // We can keep [jj] the same. We will always know how wide the [ii] is.
  // TODO: Need a statistical analysis of the likely number of zeros in such a data object
  
  /* For all current load balancing algorithms, i=1,n, j =1,i. So that LAST entry in every buffer is jsef==isef.
   */

//TODO move the ii data to registers.....

  if ( number.size() != cimat.size() || number.size() != icicol.size() ) {
    cerr << " Incorrect sizes " << number.size() << endl;
    GA::error(" Incorrect sizes ", number.size());
  }
  
  int nsefi = number.size();
  int irow;
  int lo[MAX_CI_VECTOR_DIMS];
  int hi[MAX_CI_VECTOR_DIMS];
  int n,n1;
  
  int irowBig;            // Extra stuff for managing the Supplemental storage
  int losupp[MAX_CI_VECTOR_DIMS];
  int hisupp[MAX_CI_VECTOR_DIMS];
  
  // For now push out a single I columns worth 
  
  double dValue=0.0;
  int jsef;
  
  for(int ii=0; ii< nsefi; ++ii ) {
    irow  = isef + ii;
    lo[0] = irow;
    lo[1] = 0; // 
    hi[0] = irow; //
    hi[1] = number[ii]-1;

    //dValue = cimat[ii][ number[ii]-1 ]; // diagnal term for this ii
    
    // We assume maxsparse <= maxsparseBig
    
    if ( number[ii] > l_maxsparse && number[ii] > l_maxsparse_supp ) {
      cout << "Current number of " << number[ii] << " is greater than maxsparse values " << l_maxsparse << " " << l_maxsparse_supp << endl;
      GA::error("greater than maxsparse values ", number[ii]);
    }
    
    if ( l_maxsparse_supp > l_maxsparse && number[ii] > l_maxsparse ) { 
      // Data will not fit in the default array -> push to supplemental  
#ifdef DETAILEDCHECK
      cout << "number[ii] exceeds specified MAIN GA maxsparse: Will not be able to store to MAIN GA: Attempting supplemental store." << endl;
#endif
      irowBig = nextSupplimentIrow(); // Do we have any space left in the supplental area?
      
      if ( irowBig < l_maxwidth ) {
	
	losupp[0]=irowBig;
	hisupp[0]=irowBig;
	losupp[1]=0;
	hisupp[1]=number[ii]-1;
	n = number[ii];
	n1 = 1;
	
	//cout << "Attempt adding data to the SUPP array at ii, location " << irow << " " << irowBig << endl;
	
	if ( cimat[ii].size() > 0 ) {
	  //cout << losupp[0] << " " <<  losupp[1] << " " <<  hisupp[0] << " " <<  hisupp[1] << endl;
	  l_g_cimat_supp->put( losupp, hisupp, &cimat[ii][0], &n);
	  l_g_icicol_supp->put( losupp, hisupp, &icicol[ii][0], &n);
	  hisupp[1]=0;
	  l_g_number_supp->put( losupp, hisupp, &number[ii], &n1);
	  
	  int origii = -1*(irowBig+1); // Add one to handle: using a negative value as a trigger
	  //cout << " SIZE is " << origii << " " << number[ii] << endl;
	  //cout << "Storing to SUPP array " << irow << " " << irowBig << " " << number[ii] << endl;
	  hi[1]=0;
	  l_g_number->put( lo, hi, &origii, &n1); // 
	}
      } else {
        cerr << " irowBig >= l_maxwidth " << irowBig << endl;
        GA::error(" irowBig >= l_maxwidth ",irowBig);
      }
    } else {
      if ( cimat[ii].size() > 0 ) {

//cout << lo[0] << " " <<  lo[1] << " " <<  hi[0] << " " <<  hi[1] << endl;
//cout <<"UPLOAD junk " << cimat[ii].size()<<" to irow " << irow  << endl;

for(int k=0; k< cimat[ii].size(); ++k ) {
cout << " " << cimat[ii][k]<<" "<< icicol[ii][k]<<" "<<number[ii] << endl;
}
cout << endl;

        l_g_cimat->put( lo, hi, &cimat[ii][0], &n);
        l_g_icicol->put( lo, hi, &icicol[ii][0], &n);
        hi[1]=0;
        l_g_number->put( lo, hi, &number[ii], &n1);
      }
    }
    /* DO some checks for j in case someone changes the load balancing scheme....
       If this become a bottleneck then simply shove cimat value to GA
    */

    int numii = number[ii]; // total number in this row

#ifdef SQUAREHAMILTONIAN
    dValue = -1.0;
      for(int jtest=0; jtest< numii; ++jtest ) { // now find the diagonal term
         jsef = icicol[ii][jtest];
          if (jsef == irow )  {
          dValue = cimat[ii][ jtest ];
      l_g_diag_sef->put( lo, hi, &dValue, &n1 );
        break;
        }
      }
#else
      jsef = icicol[ii][ numii-1 ]; //This is already an ABSOLUTE location
      if (jsef == irow ) {
        //cout << GA::nodeid() << " PUTTING a diag " <<  cimat[ii][ number[ii]-1 ] << endl;
        dValue = cimat[ii][ numii-1 ];
      l_g_diag_sef->put( lo, hi, &dValue, &n1 );
    }
#endif

  }
  return(0);
}

/* New accumulate based upload. In this method we need to do additional checking to see if the data must actually
   reside in supplemental. Since the calling is now looping over iconf repeatedly two things must occur here.
   1) use Accumulate GA to update rows.
   2) Check the current size of data in non-supplemental storage to see if data needs to be moved to supp.

   number is the number of entries for this pass of iconf/jchunk

  Some constraints apply:
  1) the current H construction algorithms sum all J contributions for the given i. Adding the jchunk
  loop breaks up these contributions BUT no J column is ever revisited (e.g., no duplicates exist). Thus
  we do NOT need to use ga->acc for cimat nor icicol. Only for num.
*/
int PsociGAhamiltonian::uploadCIMATtoGASupplementalAcc( int isef , vector<vector<double> > &cimat, vector<vector<int> > &icicol, vector<int> & number )
{
  // I think we should COMBINE all [ii] terms into a single row. Then upload the whole thing.
  // We can keep [jj] the same. We will always know how wide the [ii] is.
  // TODO: Need a statistical analysis of the likely number of zeros in such a data object
  
  /* For all current load balancing algorithms, i=1,n, j =1,i. So that LAST entry in every buffer is jsef==isef.
   */

//TODO move the ii data to registers.....

  if ( number.size() != cimat.size() || number.size() != icicol.size() ) {
    cerr << " Incorrect sizes " << number.size() << endl;
    GA::error(" Incorrect sizes ", number.size());
  }
  
  //double dZero=0.0;
  //int iZero=0;

  int nsefi = number.size();
  int irow;
  int lo[MAX_CI_VECTOR_DIMS];
  int hi[MAX_CI_VECTOR_DIMS];
  int n,n1;
  
  int irowBig;            // Extra stuff for managing the Supplemental storage
  int losupp[MAX_CI_VECTOR_DIMS];
  int hisupp[MAX_CI_VECTOR_DIMS];
  
  // For now push out a single I columns worth 
  
  double dValue=0.0;
  double dOne=1.0;
  double dZero=0.0;
  int iZero=0;

  int iOne=1;
  int jsef;
  int jsef_shift=0; // could be zero on first assignment indicating current supp usage
  
  //dValue = cimat[ii][ number[ii]-1 ]; // diagnal term for this ii
  
  for(int ii=0; ii< nsefi; ++ii ) {
    irow  = isef + ii;
    lo[0] = irow;
    lo[1] = 0; // these may require possible shifting
    hi[0] = irow; //
    int numii = number[ii];
    hi[1] = numii-1; // these may require possible shifting

    if ( numii > 0 ) { // to the end

/* with jchunks turned on it is possible to get an empty column here. 
   if true, then we must skip it
*/

   //cout << endl << "XXXXX I loop " << ii << " isef is " << isef << endl;

   int existing_num = 0;
   int supp_existing_num = -1; // only set to valid IF real supp data exists
   int supp_existing_row = 0; // this needs lots of testing !

   if ( cimat[ii].size() > 0 ) {
      int nex=1;
      hi[1] = 0; // refer to unshifted original location
      l_g_number->get( lo, hi, &existing_num, &nex); // guaranteed to be zero if none exist
      if ( existing_num < 0 ) { // guarentees a valid supp row exists
          supp_existing_row = (-1*existing_num) - 1;
          losupp[0]= supp_existing_row;
          hisupp[0]= supp_existing_row;
          losupp[1]=0;
          hisupp[1]=0;
          l_g_number_supp->get( losupp, hisupp, &supp_existing_num, &nex );
       }
    }

   //cout << "init nums " << existing_num <<" "<<supp_existing_num<<" "<<numii<<endl;

// hi[1],lo[1] will require resetting below
// Need to make checks below based on current Ga numvalue and  new values
// We are not yet checking to determine memory pool to use - just the final size

    int num_numii = 0;
    if ( existing_num >= 0 ) {
        num_numii = numii + existing_num; // what if it is negative?
    }

    if (  existing_num < 0 && supp_existing_num >= 0 ) {
       num_numii =  numii + supp_existing_num;
    }

// Choose memory location: regular, supp or transition from reg to supp
    
    if ( num_numii > l_maxsparse && num_numii > l_maxsparse_supp ) {
      cout << "Current number (sum) of " << num_numii << " is greater than maxsparse values " << l_maxsparse << " " << l_maxsparse_supp << endl;
      GA::error("Acc: (sum) greater than maxsparse values ", num_numii);
    }
    
// three cases to be treated,

    //cout << "Before choice num_numii " << num_numii << " " << l_maxsparse << " " << l_maxsparse_supp << endl;

    if ( l_maxsparse_supp > l_maxsparse && num_numii > l_maxsparse ) { // regardless of current situation: final data must be in supp
     if ( existing_num < 0 ) { // if empty and too big the scenrio is the same

//case one existing data already exists in supp
        if ( cimat[ii].size() > 0 ) {

          n = numii;
          int newnumii = numii + supp_existing_num;
          losupp[0]= supp_existing_row; // shift up 
          hisupp[0]= supp_existing_row;
          hisupp[1]= newnumii -1;
          losupp[1]= supp_existing_num;
          //cout << "Basic supp upload " << losupp[0] << " " <<  losupp[1] << " " <<  hisupp[0] << " " <<  hisupp[1] << endl;
          l_g_cimat_supp->put( losupp, hisupp, &cimat[ii][0], &n );
          l_g_icicol_supp->put( losupp, hisupp, &icicol[ii][0], &n );

          losupp[1]=0; // shifted back for number
          hisupp[1]=0;
          l_g_number_supp->put( losupp, hisupp, &newnumii, &n1);
        }

     } else { // existing >= 0  

//case two: existing data in regular GA, but new data moves to supp
//Need to open up a new supp row for this AND preferably zero out old data row

      irowBig = nextSupplimentIrow(); // Do we have any space left in the supplental area?

      if ( irowBig < l_maxwidth ) {
// fetch regular data into a local array

        vector<double> ext_cimat;
        vector<int> ext_icicol;

        if ( existing_num > 0 ) {
          ext_cimat.resize(existing_num, 0.0 );
          ext_icicol.resize(existing_num, 0 );
          lo[0] = irow;
          lo[1] = 0; //  get all from the beginning
          hi[0] = irow; //
          int numii = number[ii];
          n1 = 1;
          hi[1] = existing_num-1; // these may require possible shifting
          //cout << "GET existing row " << lo[0]<<" "<<hi[0]<<" "<<lo[1]<<" "<<hi[1]<<endl;
          l_g_cimat->get( lo, hi, &ext_cimat[0], &n1 );
          l_g_icicol->get( lo, hi, &ext_icicol[0], &n1 );
         }
// Now take current icicol, cimat data and merge them

// everyone needs this
        for (int i=0; i< numii; ++i ) {
            ext_cimat.push_back( cimat[ii][i] );
            ext_icicol.push_back( icicol[ii][i] );
        }

        if ( ext_cimat.size() != numii+existing_num ) {
              cout << "numii existing num are " << numii<<" "<<existing_num<<endl;
              GA::error(" ext_cimat.size() != numii+existing_num",ext_cimat.size());        
         }
// push to new supp space

        //cout <<" a total upload size of " << ext_cimat.size() <<" and "<< ext_icicol.size() << endl;
        int newnumii = numii + existing_num;
	losupp[0]=irowBig;
	hisupp[0]=irowBig;
	losupp[1]=0; // no shifting req: first time here
	hisupp[1]=newnumii-1;
	n1 = 1;
	//cout << "Attempt adding data to the SUPP array at ii, location " << irow << " " << irowBig << endl;
	if ( cimat[ii].size() > 0 ) {
	  //cout << "Put index " << losupp[0] << " " <<  losupp[1] << " " <<  hisupp[0] << " " <<  hisupp[1] << endl;
	  l_g_cimat_supp->put( losupp, hisupp, &ext_cimat[0], &n1 );
	  l_g_icicol_supp->put( losupp, hisupp, &ext_icicol[0], &n1 );
        
	  hisupp[1]=0;
	  l_g_number_supp->put( losupp, hisupp, &newnumii, &n1); // first time no acc needed
	  
	  int origii = -1*(irowBig+1); // Add one to handle: using a negative value as a trigger
	  //cout << "flush  SIZE is " << origii << " " << number[ii] << endl;
	  //cout << "Storing to SUPP array " << irow << " " << irowBig << " " << number[ii] << endl;
	  hi[1]=0;
	  l_g_number->put( lo, hi, &origii, &n1); // 

// zero out prior row

        n1 = 1;

        if ( existing_num > 0 ) {
          hi[1] = existing_num-1; //shift back
          vector<double> emptyd(existing_num, dZero);
          vector<int> emptyi(existing_num, iZero);
          l_g_cimat->put(lo, hi, &emptyd[0], &n1);
          l_g_icicol->put(lo, hi, &emptyi[0], &n1 );
	}

       }
      } else {
        cerr << " irowBig >= l_maxwidth " << irowBig << endl;
        GA::error("Acc: irowBig >= l_maxwidth ",irowBig);
      }
     } // new or old supp usage

//last case
   } else { // was using regular memory and will contonue to do so
      
      n = numii;
      if ( cimat[ii].size() > 0 ) {
	hi[1] = numii-1 + existing_num; // shift up: could use numii + existing_num;
	lo[1] = existing_num;
	
        //cout << "exis " << existing_num << " " << numii << " " << cimat[ii].size() << endl;
        //cout <<" lo hi " <<  lo[0] << " " <<  lo[1] << " " <<  hi[0] << " " <<  hi[1] << endl;
        l_g_cimat->put( lo, hi, &cimat[ii][0], &n );
        l_g_icicol->put( lo, hi, &icicol[ii][0], &n );
	
        lo[1]=0; // shifted back for number
        hi[1]=0;
        //cout << "update reg num " << endl;
        int newnum = numii + existing_num;
        l_g_number->put( lo, hi, &newnum, &n1 ); // correct
      }
    }
    /* DO some checks for j in case someone changes the load balancing scheme....
       If this become a bottleneck then simply shove cimat value to GA
    */
    
// required by all cases
    lo[1] =0; // shifted back  for diag
    hi[1] =0;

#ifdef SQUAREHAMILTONIAN
    dValue = -1.0;
    for(int jtest=0; jtest< numii; ++jtest ) { // now find the diagonal term
      jsef = icicol[ii][jtest];
      if (jsef == irow )  {
	dValue = cimat[ii][ jtest ];
	l_g_diag_sef->put( lo, hi, &dValue, &n1 );
        break;
      }
    }
#else
    jsef = icicol[ii][ numii-1 ]; //This is already an ABSOLUTE location
    if (jsef == irow ) {
       //cout << GA::nodeid() << " PUTTING a diag " <<  cimat[ii][ number[ii]-1 ] << endl;
      dValue = cimat[ii][ numii-1 ];
      l_g_diag_sef->put( lo, hi, &dValue, &n1 );
    }
#endif
    
   } // check for numii
  } // 
  return(0);
}

void PsociGAhamiltonian::zeroMatrix( vector<vector<double> >& matrix ) 
{
     int width = matrix.size();
     for(int i=0; i< width; ++i ) {
        zeroVector( matrix[i] );
     }
}

void PsociGAhamiltonian::zeroMatrix( vector<vector<int> >& matrix )
{
     int width = matrix.size();
     for(int i=0; i< width; ++i ) {
        zeroVector( matrix[i] );
     }
}

void PsociGAhamiltonian::zeroMatrix( vector<vector<double> >& matrix, int number )
{
     int width = matrix.size();
     for(int i=0; i< width; ++i ) {
        zeroVector( matrix[i] );
     }
}

void PsociGAhamiltonian::zeroMatrix( vector<vector<int> >& matrix, int number )
{
     int width = matrix.size();
     for(int i=0; i< width; ++i ) {
        zeroVector( matrix[i] );
     }
}

// Simple methoid to zero out an existig vector of doubles
// Is this faster than a reallocation?

void PsociGAhamiltonian::zeroVector( vector<double> & vector )
{
  const double dZero = 0.0;
  int num = vector.size();
  for(int i=0; i< num; ++i ) 
    {
     vector[i] = dZero;
    }
//doesnt speedup    memset(&vector[0], dZero, num*sizeof(double));
}

// Simple method to zero out subset of entries
// Is this faster than a reallocation?

void PsociGAhamiltonian::zeroVector( vector<double> & vector, int number )
{
  const double dZero = 0.0;
  int num = number;
  for(int i=0; i< num; ++i )
    {
      vector[i] = dZero;
    }
}

void PsociGAhamiltonian::zeroVector( vector<int> & vector )
{
  const int iZero = 0;
  int num = vector.size();
  for(int i=0; i< num; ++i ) 
    {
      vector[i] = iZero;
    }
}

// Simple method to zero out subset of entries
// Is this faster than a reallocation?

void PsociGAhamiltonian::zeroVector( vector<int> & vector, int number )
{
  const double iZero = 0;
  int num = number;
  for(int i=0; i< num; ++i )
    {
      vector[i] = iZero;
    }
}


/* This method destroys all GA space associated with the current Hamiltonian
   THis is not a waste since the sam edata has already been put to DRA
*/

void PsociGAhamiltonian::freeHamiltonian()
{
  if ( GA::nodeid() == 0 ) 
    { 
      cout << "Destroying Hamiltonian GA space " << endl;
      GA::printStats();
    }
  
  l_g_cimat->destroy();
  l_g_icicol->destroy();
  l_g_number->destroy();
  l_g_diag_sef->destroy();
  
#ifdef SUPPLEMENTAL
  l_g_cimat_supp->destroy();
  l_g_icicol_supp->destroy();
  l_g_number_supp->destroy();
#endif
}


///// BEGIN methods for computing the CI DENSITY

//Timer wrapped version
//long PsociGAhamiltonian::constructGlobalHamiltonianForCIDEN2WAY(  int numRoots,  GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & info )

/* Not using this anymore for the population work but it is still valid */

/*
long PsociGAhamiltonian::constructGlobalHamiltonianForCIDEN2WAY( pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long status = constructGlobalHamiltonianForCIDEN2WAY();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( status );
}
*/

/* A 2D H update scheme. This constructs a partial H suitable for building a 1e CIDENSITY 
   Several more layers of indexing are required BUT at the benefit of many fewer calls to fetch
   
   BEST compiler options are : -DTWOWAY -DPREALLOCATE -DFETCH2D -DSUPPLEMENTAL - more local memory 
                2nd BEST are : -DTWOWAY -DFETCH2D -DSUPPLEMENTAL          
                     min are : -DTWOWAY -DFETCH2D - looses SPLIT-GA storage

  Store the MUCH smaller H-matrix to the same GA as the real one. The data transformations can then be 
  performed later.

  This should be a private method eventually 

  This is pretty heavy machinery for this type of H. Nevertheless it should work fine.

  Moving forward the only thing different is the call the the new sociblocks code. flush could be used for all

  TODO combine with the regular CI code

*/
//long PsociGAhamiltonian::constructGlobalHamiltonianForCIDEN2WAY(int numRoots,  GA::GlobalArray * g_v, vector<vector<double> > & density )

/*
long PsociGAhamiltonian::constructGlobalHamiltonianForCIDEN2WAY()
{
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  if ( g_rank == 0 ) {
    cout << "I am " << GA::nodeid() << " at the fetch2d ForCIDEN 2WAY-BLOCKED SPLIT hamiltonian " << endl;
  }
  
  int counter=0;
  long allcount=0;
  
  double temp1, temp2;
  double jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0;
  
#ifdef PREALLOCATE
  int l_maxlength = l_deters->fetchMaxLength();
  vector<COMPRESS_SIZE> buffer( max(l_granularity,l_nxtval_chunk) * l_maxlength );
#endif
  
#ifdef SUPPLEMENTAL
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model ForCIDEN" << endl;
  initSuppVal();
#endif
  
#ifdef STATIC
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposition is contiguous and inclusive. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used " << endl;
  PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
  jobs.initNxtask();
  long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
  //long inexttask=-1;
  int ilo = 0;
  int ihi = l_maxspatials;
#endif
  
#ifdef PREALLOCATE
  if ( g_rank == 0 ) cout << " USING PREALLOCATED FETCH METHODS " << endl;
#endif
  
#ifdef PREALLOCATE
  JOUTFG temp;
  l_deters->preallocateTempJoutfgObject( temp );
  vector<JOUTFG> fetchi( l_nxtval_chunk, temp );
  vector<JOUTFG> fetchj( l_granularity, temp );
#else
  vector<JOUTFG> fetchi;
  vector<JOUTFG> fetchj;
#endif


  if ( GA::nodeid() == 0 ) cout << "Flush H spaces " << endl;
  flushH(); // this simply zeros out cimat,icicol, num and the supplimental spaces

  pair<int,double> info1;
  
  int numElems=0;
  long totalElems=0;
  int block_status;
  
  int totalFetches=0; // Used for Tuning studies
  
  // Build 1e appropriate Hamiltonian from scratch
  
  int countme=0;
  for(int iconf = ilo; iconf < ihi; ++iconf ) { 
    
#ifndef STATIC
    if ( nexttask == iconf ) { // in the loop stride is implemented 
#endif
      
      int ifinal = min( l_nxtval_chunk, ihi - iconf ); // ihi already adds a ONE
      
      vector<vector<vector<double> > > cimat(ifinal);
      vector<vector<vector<int> > > icicol( ifinal);
      vector<vector<int> > number(ifinal);
      
      temp2 = psociTime();
#ifdef PREALLOCATE
      l_deters->fetchAndUnpackDeterminantData2DPreallocate( iconf+1, iconf+ifinal, buffer, fetchi );
#else
      l_deters->fetchAndUnpackDeterminantData2D( iconf+1, iconf+ifinal, fetchi );
#endif
      fetchTime += psociTime() - temp2;
      ++totalFetches;

      vector<int> iindex( ifinal); // populate absolute I index array - added code but of fixed work;
      for(int ii=0; ii< ifinal; ++ii ) {
	iindex[ii]=iconf+ii;
      }
      
      vector<int> numElems(ifinal, 0); 
      
      for(int iconf_block=0; iconf_block < ifinal; ++iconf_block )
	{
	  int nsefi = fetchi[iconf_block].nsefi;
	  if ( cimat[iconf_block].size() != nsefi ) {
	    cimat[iconf_block].resize(nsefi );
	    icicol[iconf_block].resize(nsefi );
	    number[iconf_block].resize(nsefi );
	  }
	}
      
      for(int jconf=0; jconf < iconf+ifinal ; jconf += l_granularity )
	{
	  int kfinal = min( l_granularity, (iconf+ifinal-1) - jconf + 1  );
	  temp2 = psociTime();
#ifdef PREALLOCATE
	  l_deters->fetchAndUnpackDeterminantData2DPreallocate( jconf+1, jconf + kfinal, buffer, fetchj );
#else
	  fetchj.clear();
	  l_deters->fetchAndUnpackDeterminantData2D( jconf+1, jconf + kfinal, fetchj );
#endif

	  fetchTime += psociTime() - temp2;
	  ++totalFetches; // comment me out eventually
	  
	  vector<int> jindex( kfinal); // populate absolute J index array
	  for(int jj=0; jj< kfinal; ++jj ) {
	    jindex[jj]=jconf+jj;
	  }
	  
	  //SOCI always clears sefsef
	  
	  for(int iconf_block=0; iconf_block < ifinal; ++iconf_block ) //always ONE for now
	    {
	      for(int jconf_block=0; jconf_block < kfinal; ++jconf_block )
		{
		  if ( jindex[jconf_block] <= iindex[iconf_block] ) {
		    ++allcount;
		    temp2 = psociTime();
		    block_status = l_gaHamiltonian->generateSOCIblocksForCIDEN( fetchi[iconf_block], fetchj[jconf_block], sefsef ); 

		    sociTime += psociTime() - temp2;
		    
		    if ( block_status ) {
		    int numTest = packSefsef( fetchi[iconf_block], fetchj[jconf_block] , sefsef, cimat[iconf_block], icicol[iconf_block], number[iconf_block] );
		      numElems[iconf_block] += numTest;
		    }
		  } 
		}
	    }
	} // jconf
      
      for(int iconf_block=0; iconf_block < ifinal; ++iconf_block )
	{
	  if ( numElems[iconf_block] > 0 ) { // For this to not be true would be somewhat rare
	    totalElems += numElems[iconf_block];
#ifdef SUPPLEMENTAL
	    uploadCIMATtoGASupplemental( l_nsef[ fetchi[iconf_block].index - 1 ], cimat[iconf_block], icicol[iconf_block], number[iconf_block]);
#else
	    uploadCIMATtoGA( l_nsef[ fetchi[iconf_block].index - 1 ], cimat[iconf_block], icicol[iconf_block], number[iconf_block]);
#endif
	  }
	} 
#ifndef STATIC
      nexttask = jobs.nxtask( g_size );
   //   iconf = iconf + (l_nxtval_chunk - 1); // skip to next possible value from jobs.nxtask
    }
#endif
  } // iconf
  
  cout << GA::nodeid() << "CIDEN sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "CIDEN fetchTime is " << fetchTime<<endl;
  cout << "CIDEN allcount is " << allcount<<endl; //this is correct
  cout << "total CIDEN Fetches is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
  cout << GA::nodeid() << " Done after CIDEN loops " <<  totalElems << endl;
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif
  
#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif
  
  return( numberNonzeroElements( totalElems ) );
}
*/





// Additional methods required by constructGlobalHamiltonianForCIDEN2WAY
/* Trace the SEFSEF into an MO density matrix
   We need to recompute the iex/jex because the original architecture had those terms buried deep
   This is only required if iconf <> jconf
   
   density has already been sized and zero'd do not do so again.
*/
int PsociGAhamiltonian::transformSEFSEFtoMO( int numRoots, JOUTFG & iconf, JOUTFG & jconf , vector<double> & sefsef, GA::GlobalArray * g_v, vector<vector<double> > & density )
{
  int status = 0; // meaningless for now
  double dZero = 0.0;

  int indexi = iconf.index;
  int indexj = jconf.index;
  int nsefi = iconf.nsefi; // Size of SEFSEF
  int nsefj = jconf.nsefi;
  
  vector<int> iex, jex; 
  
  int istart = l_nsef[ indexi - 1 ]; // These posiiton to the start of the current sef-block
  int jstart = l_nsef[ indexj - 1 ];
  
  vector<double> vec_i, vec_j; //Need to address the need for Fdouble
  
  long iex_status;
  
  double value;
  int indSefsef;

  if ( indexi == indexj ) {
    
    vector<pair<int,int> > molist;

// Loop over all basis MOs - better to simply construct a new iterator for this.

    pair<int,int> pvalue;
    vector<int>::iterator it;
    int occupation=2;
    for(it=iconf.doubly.begin(); it!=iconf.doubly.end(); ++it) {
      pvalue.first=(*it);
      pvalue.second=occupation;
      molist.push_back(pvalue);
    }
    occupation=1;
    for(it=iconf.singly.begin(); it!=iconf.singly.end(); ++it) {
      pvalue.first=(*it);
      pvalue.second=occupation;
      molist.push_back(pvalue);
    }
    
    vector<pair<int,int> >::iterator pit;
    
    for(pit=molist.begin(); pit!=molist.end(); ++pit) {
      
      value = dZero;
      
      int orb=(*pit).first;
      int occupation=(*pit).second;
      
      int dens_index = twoDindex(orb,orb); // Get absolute offset into density MO array
      
      for(int iroot=0; iroot< numRoots; ++iroot ) {
        
	value = dZero;
	int num_i = fetchFullCIvector(iroot+1, g_v, vec_i );
	
	for(int isef=0; isef<nsefi;++isef) {
	  double  temp = vec_i[istart + isef];
	  indSefsef  = isef * nsefj;
          for(int jsef=0; jsef<nsefj;++jsef) {
	    value += sefsef[ indSefsef + jsef] * temp * vec_i[jstart + jsef] * occupation;
	  }
	} 
	density[iroot][dens_index] += value;
      }
    }
    
  } else { // indexi != indexj 
    
    iex_status = generateIEXJEXforCIDEN( iconf, jconf, iex, jex );
    if ( iex.size() > 1 || iex.size() > 1 ) return(iex_status); //GA::error("Error in IEX/JEX",-1);
    if (iex_status) {
      
      int dens_index = twoDindex(iex[0],jex[0]);
      for(int iroot=0; iroot < numRoots; ++iroot ) {
	value = dZero;
	int num_i = fetchFullCIvector(iroot+1, g_v, vec_i );
	
	for(int isef=0; isef<nsefi;++isef) {
	  double  temp = vec_i[istart + isef];
	  indSefsef  = isef * nsefj;
          for(int jsef=0; jsef<nsefj;++jsef) {
	    value += sefsef[ indSefsef + jsef] * temp * vec_i[jstart + jsef];
	  }
	}
	density[iroot][dens_index] += value;
      }
    }
  }
  return( status );
}

/* Trace the SEFSEF into an MO density matrix

   An alternative method. In this case we fetch components of vec_i
   but only Jsef worth of data
   
   density has already been sized and zero'd do not do so again.
*/
int PsociGAhamiltonian::transformSEFSEFtoMOGET(int iroot, JOUTFG & iconf, JOUTFG & jconf , vector<double> & sefsef, GA::GlobalArray * g_v , vector<double> & density )
{
  int status = 0; // meaningless for now
  double dZero = 0.0;

  int indexi = iconf.index;
  int indexj = jconf.index;
  int nsefi = iconf.nsefi; // Size of SEFSEF
  int nsefj = jconf.nsefi;
  
  vector<int> iex, jex; 
  
  int istart = l_nsef[ indexi - 1 ]; // position to the start of the current sef-block
  int jstart = l_nsef[ indexj - 1 ];
  
  long iex_status;
  
  double value = dZero;
  int indSefsef;

  if ( indexi == indexj ) {
    
    vector<double> vec_j( nsefj ); //Need to addess the GA space size of Fdouble

    int num_j = fetchPartialCIvector(iroot+1, jstart, nsefj, g_v, vec_j );

    vector<pair<int,int> > molist;

// Loop over all basis MOs - better to simply construct a new iterator for this.

    pair<int,int> pvalue;
    vector<int>::iterator it;
    int occupation=2;
    for(it=iconf.doubly.begin(); it!=iconf.doubly.end(); ++it) {
      pvalue.first=(*it);
      pvalue.second=occupation;
      molist.push_back(pvalue);
    }
    occupation=1;
    for(it=iconf.singly.begin(); it!=iconf.singly.end(); ++it) {
      pvalue.first=(*it);
      pvalue.second=occupation;
      molist.push_back(pvalue);
    }
    
    vector<pair<int,int> >::iterator pit;
    
    for(pit=molist.begin(); pit!=molist.end(); ++pit) {
      
      value = dZero;
      
      int orb=(*pit).first;
      int occupation=(*pit).second;
      
      int dens_index = twoDindex(orb,orb); // Get absolute offset into density MO array
      
	value = dZero;
	for(int isef=0; isef<nsefi;++isef) {
	  double  temp = vec_j[isef];
	  indSefsef  = isef * nsefj;
          for(int jsef=0; jsef<nsefj;++jsef) {
	    value += sefsef[ indSefsef + jsef] * temp * vec_j[jsef] * occupation;
	  }
	} 
	density[dens_index] += value;
    }
    
  } else { // indexi != indexj 
    
    vector<double> vec_j( nsefj ); //Need to addess the GA space tye of Fdouble
    int num_j = fetchPartialCIvector(iroot+1, jstart, nsefj, g_v, vec_j );
    vector<double> vec_i( nsefi ); //Need to addess the GA space tye of Fdouble
    int num_i = fetchPartialCIvector(iroot+1, istart, nsefi, g_v, vec_i );

    iex_status = generateIEXJEXforCIDEN( iconf, jconf, iex, jex );
    if ( iex.size() > 1 || iex.size() > 1 ) return(iex_status); //GA::error("Error in IEX/JEX",-1);
    if (iex_status) {
      
      int dens_index = twoDindex(iex[0],jex[0]);
	value = dZero;
	
	for(int isef=0; isef<nsefi;++isef) {
	  double  temp = vec_i[isef];
	  indSefsef  = isef * nsefj;
          for(int jsef=0; jsef<nsefj;++jsef) {
	    value += sefsef[ indSefsef + jsef] * temp * vec_j[jsef];
	  }
	}
	density[dens_index] += value;
    }
  }
  return( status );
}

/* Trace the SEFSEF into an MO density matrix

   An alternative method. In this case we only process a single-CI vector worth of density
   Further, the CI vector is prefetched into vec_i. HPCToolKit identified the fetchFullCIvector
   as a HUGE bottleneck

   We need to recompute the iex/jex because the original architecture had those terms buried deep
   This is only required if iconf <> jconf
   
   density has already been sized and zero'd do not do so again.
*/
int PsociGAhamiltonian::transformSEFSEFtoMOnoGET(JOUTFG & iconf, JOUTFG & jconf , vector<double> & sefsef, vector<double> & vec_i, vector<double> & density )
{
  int status = 0; // meaningless for now
  double dZero = 0.0;

  int indexi = iconf.index;
  int indexj = jconf.index;
  int nsefi = iconf.nsefi; // Size of SEFSEF
  int nsefj = jconf.nsefi;
  
  vector<int> iex, jex; 
  
  int istart = l_nsef[ indexi - 1 ]; // These posiiton to the start of the current sef-block
  int jstart = l_nsef[ indexj - 1 ];
  
  long iex_status;
  
  double value;
  int indSefsef;

  if ( indexi == indexj ) {
    
    vector<pair<int,int> > molist;

// Loop over all basis MOs - better to simply construct a new iterator for this.

    pair<int,int> pvalue;
    vector<int>::iterator it;
    int occupation=2;
    for(it=iconf.doubly.begin(); it!=iconf.doubly.end(); ++it) {
      pvalue.first=(*it);
      pvalue.second=occupation;
      molist.push_back(pvalue);
    }
    occupation=1;
    for(it=iconf.singly.begin(); it!=iconf.singly.end(); ++it) {
      pvalue.first=(*it);
      pvalue.second=occupation;
      molist.push_back(pvalue);
    }
    
    vector<pair<int,int> >::iterator pit;
    
    for(pit=molist.begin(); pit!=molist.end(); ++pit) {
      
      value = dZero;
      
      int orb=(*pit).first;
      int occupation=(*pit).second;
      
      int dens_index = twoDindex(orb,orb); // Get absolute offset into density MO array
      
	value = dZero;
	for(int isef=0; isef<nsefi;++isef) {
	  double  temp = vec_i[istart + isef];
	  indSefsef  = isef * nsefj;
          for(int jsef=0; jsef<nsefj;++jsef) {
	    value += sefsef[ indSefsef + jsef] * temp * vec_i[jstart + jsef] * occupation;
	  }
	} 
	density[dens_index] += value;
    }
    
  } else { // indexi != indexj 
    
    iex_status = generateIEXJEXforCIDEN( iconf, jconf, iex, jex );
    if ( iex.size() > 1 || iex.size() > 1 ) return(iex_status); //GA::error("Error in IEX/JEX",-1);
    if (iex_status) {
      
      int dens_index = twoDindex(iex[0],jex[0]);
	value = dZero;
	
	for(int isef=0; isef<nsefi;++isef) {
	  double  temp = vec_i[istart + isef];
	  indSefsef  = isef * nsefj;
          for(int jsef=0; jsef<nsefj;++jsef) {
	    value += sefsef[ indSefsef + jsef] * temp * vec_i[jstart + jsef];
	  }
	}
	density[dens_index] += value;
    }
  }
  return( status );
}

/* A truncated process ONLY useful to generating a CIDEN: Excitations are truncated here
   DO NOT USE for a generalized H matrix
*/

long PsociGAhamiltonian:: generateIEXJEXforCIDEN(JOUTFG & iconf, JOUTFG & jconf, vector<int> & iex, vector<int> & jex )
{
  long have_processed = 0; //None worth using
   iex.clear();
   jex.clear();
  
  vector<int> iocc( l_nbf, 0 );
  int sumelec = 0;
  
  for( int ia=0; ia < iconf.num_doubly; ++ia ) {
    sumelec += 2;
    iocc[ iconf.doubly[ia] ] = 2;
  }
  for( int ia=0; ia < iconf.num_singly; ++ia ) {
    ++sumelec;
    iocc[ iconf.singly[ia] ] = 1;
  }
  
  if ( sumelec != l_nelec ) {
    cerr << "GA: Probable corrupted ioutfg object. Sum iconf electrons != nelec " << sumelec << endl;
    GA::Terminate();
  }
  
  sumelec = 0;
  vector<int> jocc( l_nbf, 0 );
  
  for( int ja=0; ja < jconf.num_doubly; ++ja ) {
    sumelec += 2;
    jocc[ jconf.doubly[ja] ] = 2;
  }
  for( int ja=0; ja < jconf.num_singly; ++ja ) {
    ++sumelec;
    jocc[ jconf.singly[ja] ] = 1;
  }
  
  if ( sumelec != l_nelec ) {
    cerr << "Probable corrupted joutfg object. Sum jconf electrons != nelec " << sumelec << endl;
    GA::Terminate();
  }
  int jmarr;
  int sumiex=0;
  int sumjex=0;
  
  for(int i=0; i< l_nbf; ++i) {
    jmarr = iocc[i] - jocc[i];

    if ( jmarr < 0 ) {
      if ( (sumjex-jmarr) > 1 ) return( have_processed );
      for(int j=-1; j>=jmarr; --j) {
	jex.push_back( i );
	++sumjex;
      }
    } else if ( jmarr > 0 ) {
      if ( (sumiex+jmarr) > 1 ) return( have_processed );
      
      for(int j=1; j<=jmarr; ++j) {
	iex.push_back( i );
	++sumiex;
      }
    }
  }

// TODO is this next case really possible. It is not a performance problem so more checks is okay with me
// THis was needed when the incorrectly constructed 11 elec dets were constructed

    if (iex.size()==0 || jex.size()==0 ) {
      GA::error(" iex.size()==0 || jex.size()==0", -1 );
    }

  return( 1 ) ;//TODO render me useless or foresake me
}

/* Values are C-style and begin at 0, ID computes in a Fortran style way
   Canonical indexing of a 2D array
*/

inline int PsociGAhamiltonian::twoDindex( int mC, int nC ) {
  int m = mC + 1;
  int n = nC + 1;
  int value = ((max(m,n))*((max(m,n))-1))/2+min(m,n);
  return( value - 1 ); // convert back to C-style
}

/* Simply grab the full CI vector for index. Return the size of index as a modest check

   We anticipate all nodes needing this
*/

int PsociGAhamiltonian::fetchFullCIvector(int index, GA::GlobalArray * g_v, vector<double> & vect )
{
  int lo[2],hi[2];
  
  int n = 1;
  lo[0]=0;
  hi[0]=l_maxsef-1;

  lo[1]=index-1;
  hi[1]=index-1;
  
  vect.resize( l_maxsef );
  g_v->get( lo, hi, &vect[0], &n );
  
  return( vect.size() );
}

/* In this version we only fetch a portion of the CI vector. For the input value of iroot (1,,n)
   we simply grab the contribution for the current jsef

   vec_j must be sized as vec_j( nsefj ) or more;

*/
int PsociGAhamiltonian::fetchPartialCIvector(const int iroot, const int jstart, const int nsefj, GA::GlobalArray * g_v, vector<double> & vec_j )
{
  int lo[2],hi[2];

  int n = 1;
  lo[0]=jstart;
  hi[0]=jstart+nsefj-1;

  lo[1]=iroot-1;
  hi[1]=iroot-1;

  g_v->get( lo, hi, &vec_j[0], &n );

  return( vec_j.size() );
}


/* Simply grab the full CI vector for index. Return the size of index as a modest check
   Do this as a BRDCST

   We anticipate all nodes needing this
*/
int PsociGAhamiltonian::fetchFullCIvectorBRD(int index, GA::GlobalArray * g_v, vector<double> & vect )
{
  int lo[2],hi[2];
  
  int n = 1;
  lo[0]=0;
  hi[0]=l_maxsef-1;

  lo[1]=index-1;
  hi[1]=index-1;
  
// Better still we can do a PARALLEL READ AND A gop

  vect.resize( l_maxsef );

  if ( GA::nodeid() == 0 ) {
  g_v->get( lo, hi, &vect[0], &n );
  }
  GA::brdcst( &vect[0], vect.size()*sizeof(double), 0) ;
  
  return( vect.size() );
}

/* Timer wrapper version 
*/
/*
long PsociGAhamiltonian::constructCIdensity( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & info  )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long num = constructCIdensity(numRoots, g_v, density);
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( num );
}
*/


/* New test method to fetch and properly assemble SEFSEF data for constructing a
   distributed density array. The data are GOP'd later
   
   the dimensions and size of density are already established do not change them here.
   creates a partially filled density for all roots. A subsequent GOP is required
   
   The 1e-H is already computed. We just need to grab them from GA space
   
   Data parallel code
   
   Serious performance problems here
   
   NOTE: this needs to be done for ALL the requested roots (numroots)


   DO NOT USE THIS ONE.
*/
/*
long PsociGAhamiltonian::constructCIdensity( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density )
{
  int g_rank = GA::nodeid();   
  int g_size = GA::nodes();
  
  if ( g_rank == 0 ) cout << "constructCIdensity using numroots = " << numRoots<<endl;
  if ( g_rank == 0 ) cout << "ALTERNATIVE SEFSEF noGET " << endl;
  
  char handleMessage[] = "constructCIdensity: current CI vectors: g_v";
  g_v->checkHandle( handleMessage );
  
  double temp1, temp2;
  double jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0, fetchCITime=0.0, sefsefTime=0.0;
  int totalFetches=0;
  int allcount=0;
  long block_status;
  
  // Need to loop over the H, fetch it into sefsef blocks and pass that to the below method
  // Much like was done in actually constructing H
  
#ifdef PREALLOCATE
  int l_maxlength = l_deters->fetchMaxLength();
  vector<COMPRESS_SIZE> buffer( max(l_granularity,l_nxtval_chunk) * l_maxlength );
#endif
  
#ifdef SUPPLEMENTAL
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model For CI density" << endl;
  initSuppVal();
#endif
  
#ifdef STATIC 
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used for CI density " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposition is contiguous and inclusive. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used for CI density" << endl;
  PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
  jobs.initNxtask();
  long nexttask = jobs.nxtaskAlt( g_size ); // Keep J's agglomerated and dynamically balance I for now
  long inexttask=-1;
  int ilo = 0;
  int ihi = l_maxspatials;
#endif
  
#ifdef PREALLOCATE
  if ( g_rank == 0 ) cout << " USING PREALLOCATED FETCH METHODS for CI dendity " << endl;
#endif

  if ( g_rank == 0 ) cout << " RECOMPUTE SOCI instead of simply fetch them " << endl;

#ifdef PREALLOCATE
  JOUTFG temp;
  l_deters->preallocateTempJoutfgObject( temp );
  vector<JOUTFG> fetchi( l_nxtval_chunk, temp );
  vector<JOUTFG> fetchj( l_granularity, temp );
#else
  vector<JOUTFG> fetchi;
  vector<JOUTFG> fetchj;
#endif
  

  flushH(); // this simply zeros out cimat,icicol, num and the supplimental spaces
  
  int fetchCIvec = 0;
  vector<double> vec_i; // Gets resized in fechFullCIvector
  
  for(int iconf = ilo; iconf < ihi; ++iconf ) {
    
#ifndef STATIC
    if ( nexttask == iconf ) { // in the loop stride is implemented  
#endif
      
      int ifinal = min( l_nxtval_chunk, ihi - iconf ); // ihi already adds a ONE
      
      temp2 = psociTime();
#ifdef PREALLOCATE
      l_deters->fetchAndUnpackDeterminantData2DPreallocate( iconf+1, iconf+ifinal, buffer, fetchi );
#else
      fetchi.clear();
      l_deters->fetchAndUnpackDeterminantData2D( iconf+1, iconf+ifinal, fetchi );
#endif
      
      fetchTime += psociTime() - temp2;
      ++totalFetches;
      vector<int> iindex( ifinal); // populate absolute I index array - added code but of fixed work;
      
      for(int ii=0; ii< ifinal; ++ii ) {
	iindex[ii]=iconf+ii;
      }
      

      for(int jconf=0; jconf < iconf+ifinal ; jconf += l_granularity )
        {
	  int got_read=0; //     do/dont reread fetch's
	  int got_sefsef=0; //   do/dont recalculate sefsef
	  
	  int kfinal = min( l_granularity, (iconf+ifinal-1) - jconf + 1  );
	  vector<int> jindex( kfinal); // populate absolute J index array
	  vector<int> block_status_list (ifinal * kfinal, 0.0);
	  
	  for(int iroot=0; iroot < numRoots; ++iroot ) 
	    {
	      double timein = psociTime();
	      int num_i = fetchFullCIvector(iroot+1, g_v, vec_i );
	      ++fetchCIvec;
	      fetchCITime += psociTime() - timein;
	      
	      if ( !got_read ) { // only fetch on first root
		
		temp2 = psociTime();
#ifdef PREALLOCATE
		l_deters->fetchAndUnpackDeterminantData2DPreallocate( jconf+1, jconf + kfinal, buffer, fetchj );
#else
		fetchj.clear();
		l_deters->fetchAndUnpackDeterminantData2D( jconf+1, jconf + kfinal, fetchj );
#endif
		fetchTime += psociTime() - temp2;
		++totalFetches; // comment me out eventually
		
		for(int jj=0; jj< kfinal; ++jj ) {
		  jindex[jj]=jconf+jj;
		}
	      } // got_read
	      got_read = 1;
	      
	      // need an array to hold a group of block statuses for checking
	      // block_status[ifinal][kfinal] in size
	      // first block of i/j get all
	      
	      int block_index=0;
	      
	      for(int iconf_block=0; iconf_block < ifinal; ++iconf_block ) //always ONE for now
		{
		  for(int jconf_block=0; jconf_block < kfinal; ++jconf_block )
		    {
		      if ( jindex[jconf_block] <= iindex[iconf_block] ) {
			++allcount;
			temp2 = psociTime();
			if ( !got_sefsef ) {
			  block_status_list[ block_index ]  = l_gaHamiltonian->generateSOCIblocksForCIDEN( fetchi[iconf_block], fetchj[jconf_block], sefsef );
			  sociTime += psociTime() - temp2;
			}
			if ( block_status_list[block_index] ) { // return 0 or 1 where 1 is GOOD
			  int new_stat =  transformSEFSEFtoMOnoGET(numRoots,fetchi[iconf_block],fetchj[jconf_block],sefsef, vec_i, density[iroot]);
			}
			++block_index;
		      }
		    }
		}
	      got_sefsef=1;
	    } // new iroot
        } // jconf
      
#ifndef STATIC
      nexttask = jobs.nxtaskAlt( g_size );
      iconf = iconf + (l_nxtval_chunk - 1); // skip to next possible value from jobs.nxtask
    } // inext
#endif
  } // iconf
  
  cout << GA::nodeid() << "CIDEN sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "CIDEN fetchTime is " << fetchTime<<endl;
  cout << "CIDEN allcount is " << allcount<<endl; //this is correct
  cout << "total CIDEN Fetches is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
  //cout << "fetch CI " << fetchCIvec<< endl;
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif

#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif
  
  GA::SERVICES.sync();
  return(0);
}
*/

/* Timer wrapper version  
*/
/*
long PsociGAhamiltonian::constructCIdensityAlt( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & info  )
{
  int g_rank = GA::nodeid(); 
  double timein = psociTime();
  long num = constructCIdensityAlt(numRoots, g_v, density);
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( num );
}
*/

/* New test method to fetch and properly assemble SEFSEF data for constructing a
   distributed density array. The data are GOP'd later

   the dimensions and size of density are already established do not change them here.
   creates a partially filled density for all roots. A subsequent GOP is required

   The 1e-H is already computed. We just need to grab them from GA space

   Data parallel code

   Serious performance problems here

   NOTE: this needs to be done for ALL the requested roots (numroots)

   -DALTERNATIVEDENSITY 
   Alternative version: should be faster for small numbers of roots and maybe more so generally
   This is not a bad choice BUT it doesn't orbmap prescreewning. YOu should use that one instead.

   -DALTERNATIVEDENSITY
   DO NOT USE anymore 
*/
/*
long PsociGAhamiltonian::constructCIdensityAlt( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density )
{
  int g_rank = GA::nodeid();   
  int g_size = GA::nodes();
  
  if ( g_rank == 0 ) cout << "constructCIdensity ALT using numroots = " << numRoots<<endl;
  if ( g_rank == 0 ) cout << "ALTERNATIVE SEFSEF noGET " << endl;
  
  char handleMessage[] = "constructCIdensity: current CI vectors: g_v";
  g_v->checkHandle( handleMessage );
  
  double temp1, temp2;
  double jTime=0.0, sociTime=0.0, packTime=0.0, fetchTime =0.0, uploadTime=0.0, fetchCITime=0.0, sefsefTime=0.0;
  int totalFetches=0;
  long allcount=0;
  long block_status;
  
  // Need to loop over the H, fetch it into sefsef blocks and pass that to the below method
  // Much like was done in actually constructing H
  
#ifdef PREALLOCATE
  int l_maxlength = l_deters->fetchMaxLength();
  vector<COMPRESS_SIZE> buffer( max(l_granularity,l_nxtval_chunk) * l_maxlength );
#endif
  
#ifdef SUPPLEMENTAL
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model For CI density" << endl;
  initSuppVal();
#endif
  
#ifdef STATIC 
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used for CI density " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposition is contiguous and inclusive. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used for CI density" << endl;
  PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
  jobs.initNxtask();
  long nexttask = jobs.nxtaskAlt( g_size ); // Keep J's agglomerated and dynamically balance I for now
  long inexttask=-1;
  int ilo = 0;
  int ihi = l_maxspatials;
#endif
  
#ifdef PREALLOCATE
  if ( g_rank == 0 ) cout << " USING PREALLOCATED FETCH METHODS for CI density " << endl;
#endif

  if ( g_rank == 0 ) cout << " RECOMPUTE SOCI instead of simply fetch them " << endl;

#ifdef PREALLOCATE
  JOUTFG temp;
  l_deters->preallocateTempJoutfgObject( temp );
  vector<JOUTFG> fetchi( l_nxtval_chunk, temp );
  vector<JOUTFG> fetchj( l_granularity, temp );
#else
  vector<JOUTFG> fetchi;
  vector<JOUTFG> fetchj;
#endif
  

  flushH(); // this simply zeros out cimat,icicol, num and the supplimental spaces
  
  int fetchCIvec = 0;
  vector<double> vec_i; // Gets resized in fechFullCIvector
  
  for(int iconf = ilo; iconf < ihi; ++iconf ) {
    
#ifndef STATIC
    if ( nexttask == iconf ) { // in the loop stride is implemented  
#endif

      int ifinal = min( l_nxtval_chunk, ihi - iconf ); // ihi already adds a ONE
      
      temp2 = psociTime();
#ifdef PREALLOCATE
      l_deters->fetchAndUnpackDeterminantData2DPreallocate( iconf+1, iconf+ifinal, buffer, fetchi );
#else
      fetchi.clear();
      l_deters->fetchAndUnpackDeterminantData2D( iconf+1, iconf+ifinal, fetchi );
#endif
      
      fetchTime += psociTime() - temp2;
      totalFetches += ifinal;
      vector<int> iindex( ifinal); // populate absolute I index array - added code but of fixed work;
      
      for(int ii=0; ii< ifinal; ++ii ) {
	iindex[ii]=iconf+ii;
      }
      

      for(int iroot=0; iroot < numRoots; ++iroot )
	{
	  double timein = psociTime();
	  int num_i = fetchFullCIvector(iroot+1, g_v, vec_i );
	  ++fetchCIvec;
	  fetchCITime += psociTime() - timein;
	  
	  for(int jconf=0; jconf < iconf+ifinal ; jconf += l_granularity )
	    {
	      int got_read=0; //     do/dont reread fetch's
	      int got_sefsef=0; //   do/dont recalculate sefsef
	      
	      int kfinal = min( l_granularity, (iconf+ifinal-1) - jconf + 1  );
	      vector<int> jindex( kfinal); // populate absolute J index array
	      vector<int> block_status_list (ifinal * kfinal, 0.0);
	  
	      temp2 = psociTime();
#ifdef PREALLOCATE
	      l_deters->fetchAndUnpackDeterminantData2DPreallocate( jconf+1, jconf + kfinal, buffer, fetchj );
#else
	      fetchj.clear();
	      l_deters->fetchAndUnpackDeterminantData2D( jconf+1, jconf + kfinal, fetchj );
#endif
	      fetchTime += psociTime() - temp2;
	      ++totalFetches; // comment me out eventually
	      
	      for(int jj=0; jj< kfinal; ++jj ) {
		jindex[jj]=jconf+jj;
	      }
	      
	      // need an array to hold a group of block statuses for checking
	      // block_status[ifinal][kfinal] in size
	      // first block of i/j get all
	      
	      int block_index=0;
	      
	      for(int iconf_block=0; iconf_block < ifinal; ++iconf_block ) //always ONE for now
		{
		  for(int jconf_block=0; jconf_block < kfinal; ++jconf_block )
		    {
		      if ( jindex[jconf_block] <= iindex[iconf_block] ) {
			++allcount;
			temp2 = psociTime();
			block_status_list[ block_index ]  = l_gaHamiltonian->generateSOCIblocksForCIDEN( fetchi[iconf_block], fetchj[jconf_block], sefsef );
			sociTime += psociTime() - temp2;
		      }
		      if ( block_status_list[block_index] ) { // return 0 or 1 where 1 is GOOD
                        int new_stat =  transformSEFSEFtoMOnoGET(numRoots,fetchi[iconf_block],fetchj[jconf_block],sefsef, vec_i, density[iroot]);
		      }
		      ++block_index;
		    }
		}
	    } // jconf
	} // new iroot
      
#ifndef STATIC
      nexttask = jobs.nxtaskAlt( g_size );
      iconf = iconf + (l_nxtval_chunk - 1); // skip to next possible value from jobs.nxtask
    } // inext
#endif
  } // iconf
  
  cout << GA::nodeid() << "CIDEN sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "CIDEN fetchTime is " << fetchTime<<endl;
  cout << "CIDEN allcount is " << allcount<<endl; //this is correct
  cout << "total CIDEN Fetches is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif

#ifndef STATIC
  jobs.destroyNxtask(); //Close'er up
#endif
  
  GA::SERVICES.sync();
  return(0);
}
*/

/* Timer wrapper version  
*/
long PsociGAhamiltonian::constructCIdensityAltOrbMap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & info  ) 
{
  int g_rank = GA::nodeid(); 
  double timein = psociTime();
  long num = constructCIdensityAltOrbMap(numRoots, g_v, density);
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( num );
}

/* new version: removed chunking and includes checks on nxtval_chunk. This version includes ORBMAP style prescreening.
   The orbmap construction difference form the standard H in that only 1e substitutions are relevant.  This should be the default.

   Moved the iroot loop to the outermost iterator. We rarely statistices on lots of roots and without caching
   we were doing many rereads as a function of iconf

   THIS IS THE CORRECT ONE TO USE
*/
long PsociGAhamiltonian::constructCIdensityAltOrbMap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density )
{
  int g_rank = GA::nodeid();   
  int g_size = GA::nodes();
  
  if ( g_rank == 0 ) cout << "FIXED Xform constructCIdensity ALT ORBMAP version using numroots = " << numRoots<<endl;
  if ( g_rank == 0 ) cout << "ALTERNATIVE SEFSEF noGET " << endl;
  
  char handleMessage[] = "constructCIdensity: current CI vectors: g_v";
  g_v->checkHandle( handleMessage );
  
  double temp2, temp3;
  double transTime=0.0, mapTime=0.0,sociTime=0.0, fetchTime =0.0, fetchCITime=0.0;
  int totalFetches=0;
  long allcount=0;
  long block_status;
  
  pair<int,double> info1;
  
  // Need to loop over the H, fetch it into sefsef blocks and pass that to the below method
  // Much like was done in actually constructing H
  
#ifdef XPREALLOCATE
  int l_maxlength = l_deters->fetchMaxLength();
  vector<COMPRESS_SIZE> buffer( l_maxlength );
#endif
  
#ifdef SUPPLEMENTAL
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model For CI density" << endl;
  initSuppVal();
#endif
  
#ifdef STATIC 
  if ( g_rank == 0 ) cout << " Simple static load balancing scheme being used for CI density " << endl;
  int ilo, ihi;
  l_deters->localConfRange( ilo, ihi );
  cout << "Local Determinant range (inclusive) is " << ilo << " " << ihi << endl;
  ++ihi; // This decomposition is contiguous and inclusive. We add ONE so that subsequent loop argument can be '<' instead of '<='
#else
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used for CI density" << endl;
  PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
  //jobs.initNxtask();
  int ilo = 0;
  int ihi = l_maxspatials;
#endif
  
#ifdef XPREALLOCATE
  if ( g_rank == 0 ) cout << " USING PREALLOCATED FETCH METHODS for CI dendity " << endl;
#endif
  if ( g_rank == 0 ) cout << " RECOMPUTE SOCI instead of simply fetch them " << endl;

  if ( g_rank == 0 ) cout << " MODIFIED get version of Xform sefseftoMO " << endl;
  
#ifdef NEWORBMAP
  if ( GA::nodeid() ==0 ) cout << "USING NEW GOP based ORBMAP approach " <<  l_maxspatials<< endl;
  vector<vector<pair<short,short> > > ijpairlist( l_maxspatials );
  int listsize = l_deters->fetchAllCompressedOrbMap( ijpairlist ); // Collective: performs a BRDCST underneath
  if ( GA::nodeid() ==0 ) cout << "replicated orbmap size is " << listsize << endl;
#endif

#ifndef NEWORBMAP
//    vector<COMPRESS_SIZE> imap( l_nbf  );
#endif
    vector<int> list;
    
    JOUTFG fetchi;
    JOUTFG fetchj;
    
/* NOTE: THis approach is much sparser than the original H matrix.
   so we may want to considfer modification of chunk sizes
*/

#ifndef DIRECTMATRIXPRODUCTS
// probably do not need this at all
    flushH(); // this simply zeros out cimat,icicol, num and the supplimental spaces
#endif

    
    vector<double> vec_i; // Gets resized in fechFullCIvector
    
    int countrun;

/* Read the vector First. In the future we could consider doing this by part s as well
*/
    for(int iroot=0; iroot < numRoots; ++iroot )
      {
      vector<COMPRESS_SIZE> imap( l_nbf  ); // gets hidden below
      //vector<COMPRESS_SIZE> jmap( l_nbf+2  );

      JOUTFG fetchi;
      JOUTFG fetchj;
      vector<double> vec_i; // Gets resized in fechFullCIvector

      countrun=0;

#ifndef STATIC
        jobs.initNxtask();
//        jobs.nxtaskReset();
        int nexttask = jobs.nxtask( g_size ); 
#endif
	for(int iconf = ilo; iconf < ihi; ++iconf ) {

#ifndef STATIC
	  if ( nexttask == iconf ) { // in the loop stride is implemented  
#endif
	    temp2 = psociTime();
#ifdef XPREALLOCATE
	    l_deters->fetchAndUnpackDeterminantDataPreallocate( iconf+1, buffer, fetchi );
#else
	    l_deters->fetchAndUnpackDeterminantData( iconf+1, fetchi );
#endif
	    fetchTime += psociTime() - temp2;
	    ++totalFetches;
	    
	    // Compute fullon map
	    temp3 = psociTime();
#ifndef NEWORBMAP
	    vector<COMPRESS_SIZE> imap( l_nbf+2  );
	    l_deters->fetchOrbMap( iconf+1, imap);
	    countrun +=  l_deters->compareAllOrbMapCIDEN(iconf+1, imap, list); // fetch jlist for each i usually 100-300 entries
#else
	    countrun +=  l_deters->compareAllOrbMapRemoveZerosCIDEN(iconf+1, ijpairlist, list );
#endif
	    mapTime += psociTime() - temp3;
	    allcount += iconf+1;
	    
	    vector<int>::iterator jit;
	    
	    /* eventually we will need to rewrite this to eliminate the fetch and store of an 
	       entire CI vector. It will be too long
	    */
	    for(jit = list.begin(); jit != list.end(); ++jit ) 
	      {
		int jconf = (*jit);
		
		temp2 = psociTime();
#ifdef XPREALLOCATE
		l_deters->fetchAndUnpackDeterminantDataPreallocate( jconf+1, buffer, fetchj ); 
#else
		l_deters->fetchAndUnpackDeterminantData( jconf+1, fetchj );
#endif
		fetchTime += psociTime() - temp2;
		++totalFetches; // comment me out eventually
		
		temp2 = psociTime();
		block_status = l_gaHamiltonian->generateSOCIblocksForCIDEN( fetchi, fetchj, sefsef );
		sociTime += psociTime() - temp2;
	      //}
	    if ( block_status ) 
	      { // return 0 or 1 where 1 is GOOD
                temp3 = psociTime();
//		int new_stat =  transformSEFSEFtoMOnoGET(fetchi,fetchj,sefsef, vec_i, density[iroot]);
                int new_stat =  transformSEFSEFtoMOGET(iroot, fetchi, fetchj,sefsef, g_v, density[iroot]);
                transTime += psociTime() - temp3;

	      } 
              }
	    
#ifndef STATIC
	    nexttask = jobs.nxtask( g_size );
	  } // inext
#endif
	} // iconf
#ifndef STATIC
      jobs.destroyNxtask(); //Close'er up
#endif
      } // new iroot
    
    if ( GA::nodeid() == 0 ) {
    cout << GA::nodeid() << " new CIDEN sociTime is " << sociTime<<endl;
    cout << GA::nodeid() << " new CIDEN fetchTime is " << fetchTime<<endl;
    cout << GA::nodeid() << " new CIDEN mapTime is " << mapTime<<endl;
    cout << GA::nodeid() << " countrun is " << countrun << endl;
    cout << GA::nodeid() << " fetchCIvectime is " << fetchCITime << endl;
    cout << GA::nodeid() << " transTime is " << transTime << endl;
    cout << "CIDEN allcount is " << allcount<<endl; //this is correct
    cout << "total CIDEN Fetches is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
    }
    
#ifdef SUPPLEMENTAL
    destroySuppVal();
#endif
    
#ifndef STATIC
    jobs.destroyNxtask(); //Close'er up
#endif
    
    GA::SERVICES.sync();
    return(0);
}


/* Need a method that will grab all elements of H(iconf,jconf) and pack them into 
   sefsef[nsefi][nsefj]. sefsef must be perfectly sized. sefsef is subsequently used to 
   construct CI density over MOs

   NOTE we muct account for a possibly SPLIT data layout.

   Performance has not been a fctor in the design of this routine - but it should be...

   return the number of elements found
*/
int PsociGAhamiltonian::fetchSEFSEF( JOUTFG & iconf, JOUTFG & jconf, vector<double>  & sefsef )
{
#ifndef SUPPLEMENTAL
  GA::error("fetchSEFSEF:: Must compile code with -DSUPPLEMENTAL ",-1);
#endif
  
  // TODO if we read iconf try to reuse this with a static declarator.
  // Then simply re-mine the icicol array for the proper values

  int count = 0; // empty I/J interaction (only comapres against input jconf not the remaining Js
  int have_processed = 0; // None processed
  
  sefsef.clear();
  
  int indexi = iconf.index;
  int indexj = jconf.index;
  int nsefi = iconf.nsefi; // Size of SEFSEF
  int nsefj = jconf.nsefi;
    
  // determine range of indexing that applis to this request
  
  int istart = l_nsef[ indexi - 1 ]; // absolute index for starting point of Iconf
  
  int jstart = l_nsef[ indexj - 1 ]; // absolute index for starting point of Jconf
  
  sefsef.resize ( nsefi * nsefj );

  // Buffer for holding H related stuff: Note for a rerad I, we get ALL J's not just those request - need to leverage that fact
  
  vector<int> icicol( l_maxsparse, 0 );
  vector<double> cimat( l_maxsparse, 0.0 );
  int cilo[2], cihi[2];
  int n=0;
  
  int number;
  int suppRow;

  for(int i=0;i< nsefi; ++i ) { // Loop over all sefs for the given iconf  sefsef[iconf][*]
  
    int iabs = i + istart; // current absolute isef value

    cilo[0] = iabs;
    cihi[0] = iabs;
    cilo[1] = 0; // always true
    cihi[1] = 0;

    l_g_number->get( cilo, cihi, &number, &n ); // How big or is it in supplemental?

    if ( number >= 0 ) { //Then this points to an absolute area in supplemental
#ifdef DETAILEDHAMCHECK
      if ( number > l_maxsparse ) {
	cout << "CI density: number exceeds specified maxsparse: Will not be able to fetch the main GA. Try supplemental ." << endl;
	cout << "CI density: number = " << number << " maxsparsity is " << l_maxsparse << endl;
        }
#endif
      cihi[1] = number - 1;
      n = number;
      l_g_icicol->get( cilo, cihi, &icicol[ 0 ], &n);
      l_g_cimat->get( cilo, cihi, &cimat[ 0 ], &n);
      
    } else {
      
      suppRow = -number - 1;
      cilo[0] = suppRow;
      cihi[0] = suppRow; 
      
      l_g_number_supp->get( cilo, cihi, &number, &n );
#ifdef DETAILEDHAMCHECK
      if ( number > l_maxsparse_supp ) {
	cerr << "CI density: number exceeds specified maxsparse_supp: Will not be able to fetch the GA_supp." << endl;
	cerr << "CI density:  number = " << number << " maxsparsity_supp is " << l_maxsparse_supp << endl;
	GA::error("CI density: number exceeds specified maxsparse_supp:", l_maxsparse_supp);
      }
#endif
      cihi[1] = number - 1;
      n = number;
      l_g_icicol_supp->get( cilo, cihi, &icicol[ 0 ], &n);
      l_g_cimat_supp->get( cilo, cihi, &cimat[ 0 ], &n);
    }

   
    
    /* At this point we have a singler ISEF(and all non-zero JSEFS associated wiht it)
       Now we need to disentagle the values into dense groups of nsefi and nsefj
       Push current cimat into a sefsef style array 
       current absolute I is  iabs. So we want absolute J values (0,iabs) only;
       Also for a given Iconf, we only have JCONF <= ICONF at the most The very last value is ICONF
    */

// Two cases: jconf < iconf and jconf == iconf.

    int jmax;
    ( indexi != indexj )? jmax = jstart + nsefj: jmax = iabs+1; //a "+1" for a consistent "<" check  
    
    int jindex=0;
    int itemp;
    for(int ib=0; ib<number; ++ib ) {
        itemp =  icicol[ib];
        if ( itemp >= jmax ) {
	  break;
        } else if (  itemp >= jstart && itemp < jmax ) {
	  ++count;
	  sefsef[ i*nsefj + (itemp - jstart) ] = cimat[ ib ];
        }
        ++jindex;
    }
  } // end i
  
#if 0
  int index=0;
  for(int i=0; i<nsefi; ++i ) {
    for(int j=0; j<nsefj; ++j ) {
      
      index = i * nsefj + j;
      cout << "SEFSEF is " <<indexi<<" "<<indexj<<" "<<i<<" "<<j<<" "<<index<<" "<< sefsef[ index] <<endl;
    }
  }
#endif
  
  if ( count > 0 ) have_processed = 1;
  return ( have_processed );
}

/* Timer wrapper version  
   for the new DIRECT-based version
*/

long PsociGAhamiltonian::constructCIdensityAltOrbMapDirect( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & info  )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long num = constructCIdensityAltOrbMapDirect(numRoots, g_v, density);
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( num );
}

/* new version: removed chunking and includes checks on nxtval_chunk. This version includes ORBMAP style prescreening.
   The orbmap construction difference form the standard H in that only 1e substitutions are relevant.  This should be the default.

   or compile with the following. -DDIRECT -DSINGLENODECONFIGS


   Moved the iroot loop to the outermost iterator. We rarely want statistics on lots of roots and without caching
   we were doing many rereads as a function of iconf

*/

long PsociGAhamiltonian::constructCIdensityAltOrbMapDirect( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density )
{
  int g_rank = GA::nodeid();   
  int g_size = GA::nodes();
  
  if ( g_rank == 0 ) cout << "constructCIdensity ALT ORBMAP DIRECT dets version using numroots = " << numRoots<<endl;
  if ( g_rank == 0 ) cout << "ALTERNATIVE SEFSEF noGET " << endl;
  
  char handleMessage[] = "constructCIdensity: current CI vectors: g_v";
  g_v->checkHandle( handleMessage );
  
#ifdef STATIC
  if ( g_rank == 0 ) cout << "CIdensity: static load balance" << endl;
#endif
  double temp2, temp3;
  double transTime=0.0, mapTime=0.0, sociTime=0.0, fetchTime =0.0, fetchCITime=0.0;
  int totalFetches=0;
  long allcount=0;
  long block_status;
  
  pair<int,double> info1;
  
  // Need to loop over the H, fetch it into sefsef blocks and pass that to the below method
  // Much like was done in actually constructing H
  
#ifdef SUPPLEMENTAL
  if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model For CI density" << endl;
  initSuppVal();
#endif
  
#ifndef STATIC
  if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used for CI density" << endl;
  PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
//  int ilo = 0;
//  int ihi = l_maxspatials;
#endif
  int ilo = 0;
  int ihi = l_maxspatials;

  
  if ( g_rank == 0 ) cout << " RECOMPUTE DIRECT SOCI instead of simply fetch them " << endl;
  
#ifdef NEWORBMAP
  if ( GA::nodeid() ==0 ) cout << "USING NEW GOP based ORBMAP approach " << endl;
  vector<vector<pair<short,short> > > ijpairlist( l_maxspatials );
  int listsize = l_deters->fetchAllCompressedOrbMap( ijpairlist ); // Collective: performs a BRDCST underneath
  if ( GA::nodeid() ==0 ) cout << "replicated orbmap size is " << listsize << endl;
#endif

    
/* NOTE: THis approach is much sparser than the original H matrix.
   so we may want to consider modification of chunk sizes
*/

#ifndef DIRECTMATRIXPRODUCTS
// Not sure why this is here: we do not use it
  flushH(); // this simply zeros out cimat,icicol, num and the supplimental spaces
#endif
  
  vector<int> list;
  int countrun;
  
  for(int iroot=0; iroot < numRoots; ++iroot )
    {
      vector<COMPRESS_SIZE> imap( l_nbf+2  );
      vector<COMPRESS_SIZE> jmap( l_nbf+2  );
      
      JOUTFG fetchi;
      JOUTFG fetchj;
     // vector<double> vec_i; // Gets resized in fechFullCIvector
      
      countrun=0;
      
#ifndef STATIC
      jobs.initNxtask();
      long nexttask = jobs.nxtask( g_size );     
#endif

#ifdef STATIC
      for(int iconf = g_rank; iconf < ihi; iconf += g_size ) {
#endif

#ifndef STATIC
        for(int iconf = ilo; iconf < ihi; ++iconf ) {
	if ( nexttask == iconf ) { // in the loop stride is implemented  
#endif

	  
	  temp3 = psociTime();
#ifndef NEWORBMAP
//	  vector<COMPRESS_SIZE> imap( l_nbf+2  );
	  l_deters->fetchOrbMap( iconf+1, imap);
	  countrun +=  l_deters->compareAllOrbMapCIDEN(iconf+1, imap, list); // fetch jlist for each i usually 100-300 entries
#else
	  countrun +=  l_deters->compareAllOrbMapRemoveZerosCIDEN(iconf+1, ijpairlist, list );
#endif
	  temp2 = psociTime();
	  l_deters->computeConfigsGA( iconf+1, imap, fetchi );
	  fetchTime += psociTime() - temp2;
	  ++totalFetches;
	  
	  mapTime += psociTime() - temp3;
	  allcount += iconf+1;
	  
	  vector<int>::iterator jit;
	  
	  for(jit = list.begin(); jit != list.end(); ++jit ) 
	    {
	      int jconf = (*jit);
	      
	      temp2 = psociTime();
	      
	      l_deters->computeConfigsGA( jconf+1, jmap, fetchj );
	      
	      fetchTime += psociTime() - temp2;
	      ++totalFetches; // comment me out eventually
	      
	      temp2 = psociTime();
	      block_status = l_gaHamiltonian->generateSOCIblocksForCIDEN( fetchi, fetchj, sefsef );
	      sociTime += psociTime() - temp2;
	      if ( block_status ) 
		{ // return 0 or 1 where 1 is GOOD
		  temp3 = psociTime();
		//  int new_stat =  transformSEFSEFtoMOnoGET(fetchi,fetchj,sefsef, vec_i, density[iroot]);
		  int new_stat =  transformSEFSEFtoMOGET(iroot, fetchi, fetchj,sefsef, g_v, density[iroot]);
		  transTime += psociTime() - temp3;
		  
		} 
	    }
	  
#ifndef STATIC
	  nexttask = jobs.nxtask( g_size );
	} // inext
#endif
      } // iconf
#ifndef STATIC
      jobs.destroyNxtask(); //Close'er up
#endif
    } // new iroot
  
  if ( GA::nodeid() == 0 ) {
  cout << GA::nodeid() << " new direct CIDEN sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << " new direct CIDEN fetchTime is " << fetchTime<<endl;
  cout << GA::nodeid() << " new direct CIDEN mapTime is " << mapTime<<endl;
  cout << GA::nodeid() << " countrun is " << countrun << endl;
  cout << GA::nodeid() << " fetchCIvectime is " << fetchCITime << endl;
  cout << GA::nodeid() << " transTime is " << transTime << endl;
  cout << "CIDEN allcount is " << allcount<<endl; //this is correct
  cout << "total CIDEN Fetches is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
  }
  
#ifdef SUPPLEMENTAL
  destroySuppVal();
#endif
  
  GA::SERVICES.sync();
  return(0);
}

/* new versions of the Direct CI density using a distributed chunking scheme
   Not yet fully tested
*/
long PsociGAhamiltonian::constructCIdensityAltOrbMapDirectChunkJmap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & info  )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long num = constructCIdensityAltOrbMapDirectChunkJmap(numRoots, g_v, density);
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( num );
}

long PsociGAhamiltonian::constructCIdensityAltOrbMapChunkJmap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & info  )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long num = constructCIdensityAltOrbMapChunkJmap(numRoots, g_v, density);
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( num );
}
 
/* Added basic screening stuff to the CI step - direct only 
*/
 long PsociGAhamiltonian::constructCIdensityAltOrbMapDirectChunkJmap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density )
 {
   int g_rank = GA::nodeid();   
   int g_size = GA::nodes();
   
   if ( g_rank == 0 ) cout << "constructCIdensity new JCHUNK ORBMAP DIRECT dets version using numroots = " << numRoots<<endl;
   if ( g_rank == 0 ) cout << "ALTERNATIVE SEFSEF noGET " << endl;
   
   char handleMessage[] = "constructCIdensity: current CI vectors: g_v";
   g_v->checkHandle( handleMessage );
  
   double temp2, temp3;
   double orbmapTime=0.0, transTime=0.0, mapTime=0.0, sociTime=0.0, fetchTime =0.0, fetchCITime=0.0;
   int totalFetches=0;
   long allcount=0;
   long block_status;

  int didblock=0;
  long totalblock=0;
  int excitblock=0;
  int phaseblock=0;
  int possibleblock=0;

   
   pair<int,double> info1;
   
   // Need to loop over the H, fetch it into sefsef blocks and pass that to the below method
   // Much like was done in actually constructing H
   
#ifdef NEWORBMAP
   GA::error(" do not use constructGlobalHamiltonianAltOrbMapDirectChunkJmap method with NEWORBMAP",-1);
#endif
   
#ifdef SUPPLEMENTAL
   if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model For CI density" << endl;
   initSuppVal();
#endif
   
#ifdef STATIC 
    if ( g_rank == 0 ) cout << "CIdensity: STATIC "<< endl; 
#endif

#ifdef REVERSE
  if ( g_rank == 0 ) cout << "DIRECT: CIdensity Using reverse iconf indexing " << endl;
#endif

#ifndef STATIC
   if ( g_rank == 0 ) cout << " NXTVAL based load balancing scheme being used for CI density" << endl;
#endif

   int ilo;
   int ihi;
   
   if ( g_rank == 0 ) cout << " RECOMPUTE DIRECT SOCI instead of simply fetch them " << endl;
    
   /* NOTE: THis approach is much sparser than the original H matrix.
      so we may want to consider modification of chunk sizes
   */
   
#ifdef REPLICATEDETVALUE
  if ( g_rank == 0 ) cout << "USING NEW JPHASE code for CIDEN build"<<endl; ;
   vector<char> * rep_phases = l_deters->fetchReplicatedPhasesHandle();
 //  cout << "H build rep phase is size " << rep_phases->size()<<endl; // returns iopen
#endif

   //CHECK ME flushH(); // this simply zeros out cimat,icicol, num and the supplimental spaces
   
   vector<int> list;
   int currentrun=0, countrun=0;
   
   double jfetchTime = 0.0;
   const int numjchunks=l_numJchunksHbuild; // carve up jmap space into 10 chunks
   
   int chunk_width;
   int l_maxspatials = l_deters->fetchMaxSpatials();
   
   (l_maxspatials % numjchunks==0)? chunk_width = l_maxspatials / numjchunks: chunk_width = 1 + (l_maxspatials / numjchunks);
   if ( GA::nodeid() ==0 ) cout << "Jchunk stats num: "<<numjchunks<<" width " << chunk_width << endl;
   int jchunk_index=0; // these track the actual real index numbers (0,maxefs-1)
   
   for(int iroot=0; iroot < numRoots; ++iroot )
     {
       vector<COMPRESS_SIZE> imap( l_nbf+2  );
       vector<COMPRESS_SIZE> jmap( l_nbf+2  );
       
       JOUTFG fetchi;
       JOUTFG fetchj;
       // vector<double> vec_i; // Gets resized in fechFullCIvector
       
       countrun=0;
       for(int jchunk=0; jchunk < numjchunks; ++jchunk ) { // data parallel
	 
	 double tempj = psociTime();
	 vector<COMPRESS_SIZE> jsubsetlist; // jmap subset buffer gets sized in fetchJsubsetlist 
	 vector<int> jindexlist;
	 int jnumlist = fetchJsubsetlist( jchunk, jindexlist, jsubsetlist, chunk_width, l_maxspatials ); //GA on the inside. resize width inside to be perfect
	 jfetchTime += psociTime() - tempj;
	 
	 // start original code
	 
#ifndef STATIC
	 PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
#endif
	 ilo = 0;
	 ihi = l_maxspatials;


#ifndef STATIC

#ifdef REVERSE
        jobs.initNxtask( l_maxspatials );
#else
        jobs.initNxtask();
#endif

#ifdef REVERSE
         long nexttask = jobs.nxtaskRev( g_size, l_maxspatials );
#else
         long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
#endif
#endif

#ifndef STATIC
#ifdef REVERSE
         for( int iconf = ihi-1; iconf >= ilo; --iconf ) {
#else
         for( int iconf = ilo; iconf < ihi; ++iconf ) {
#endif
#endif
#ifdef STATIC
#ifdef REVERSE
         for( int iconf = ihi-1-g_rank; iconf >= ilo; iconf-=g_size ) {
#else
         for( int iconf = g_rank; iconf < ihi; iconf +=g_size) {
#endif
#endif

           totalblock += (iconf+1);

#ifndef STATIC
	   if ( nexttask == iconf ) { // in the loop stride is implemented  
#endif
	     
	     //	  vector<COMPRESS_SIZE> imap( l_nbf+2  );
	     l_deters->fetchOrbMap( iconf+1, imap);
	     
	     temp2 = psociTime();
	     int currentrun =  l_deters->compareLocalJSubsetOrbMapCIDEN(jnumlist, iconf+1, rep_phases, imap, jindexlist, jsubsetlist, list); // 
	     countrun += currentrun;
	     orbmapTime += psociTime() - temp2;
	     
	     if ( currentrun != 0 ) {
	       
	       l_deters->computeConfigsGA( iconf+1, imap, fetchi );
	       fetchTime += psociTime() - temp2;
	       ++totalFetches;
	       allcount += iconf+1;
	       
	       vector<int>::iterator jit;
	       for(jit = list.begin(); jit != list.end(); ++jit ) 
		 {
                   ++didblock;
		   int jconf = (*jit);
		   temp2 = psociTime();
		   l_deters->computeConfigsGA( jconf+1, jmap, fetchj );
		   fetchTime += psociTime() - temp2;
		   ++totalFetches; // comment me out eventually
		   
		   temp2 = psociTime();
		   block_status = l_gaHamiltonian->generateSOCIblocksForCIDEN( fetchi, fetchj, sefsef );
		   sociTime += psociTime() - temp2;
		   if ( block_status ) 
		     { // return 0 or 1 where 1 is GOOD
		       temp3 = psociTime();
		       //  int new_stat =  transformSEFSEFtoMOnoGET(fetchi,fetchj,sefsef, vec_i, density[iroot]);
		       int new_stat =  transformSEFSEFtoMOGET(iroot, fetchi, fetchj,sefsef, g_v, density[iroot]);
		       transTime += psociTime() - temp3;
		     } 
		 }
	     } // currentrun 
#ifndef STATIC
#ifdef REVERSE
             nexttask = jobs.nxtaskRev( g_size, l_maxspatials );
#else
             nexttask = jobs.nxtask( g_size );
#endif
	   } // inext
#endif
	 } // iconf
#ifndef STATIC
	 jobs.destroyNxtask(); //Close'er up
#endif
       } // jchunk
     } // new iroot
   
   if ( GA::nodeid() == 0 ) {
     cout << GA::nodeid() << " new direct CIDEN sociTime is " << sociTime<<endl;
     cout << GA::nodeid() << " new direct CIDEN fetchTime is " << fetchTime<<endl;
     cout << GA::nodeid() << " countrun is " << countrun << endl;
     cout << GA::nodeid() << " fetchCIvectime is " << fetchCITime << endl;
     cout << GA::nodeid() << " transTime is " << transTime << endl;
     cout << GA::nodeid() << " jfetchTime is " << jfetchTime << endl;
     cout << GA::nodeid() << " orbmapTime is " << orbmapTime << endl;
     cout << GA::nodeid() << " didblock is " << didblock<<" "<<totalblock << endl;
     cout << "CIDEN allcount is " << allcount<<endl; //this is correct
     cout << "total CIDEN Fetches is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
   }
   
#ifdef SUPPLEMENTAL
   destroySuppVal();
#endif
   
   GA::SERVICES.sync();

   return(0);
 }
 long PsociGAhamiltonian::constructCIdensityAltOrbMapChunkJmap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density )
 {
   int g_rank = GA::nodeid();   
   int g_size = GA::nodes();
   
   if ( g_rank == 0 ) cout << "constructCIdensity new JCHUNK ORBMAP NONDIRECT dets version using numroots = " << numRoots<<endl;
   if ( g_rank == 0 ) cout << "ALTERNATIVE SEFSEF noGET " << endl;
   
   char handleMessage[] = "constructCIdensity: current CI vectors: g_v";
   g_v->checkHandle( handleMessage );
  
   double temp2, temp3;
   double orbmapTime=0.0, transTime=0.0, mapTime=0.0, sociTime=0.0, fetchTime =0.0, fetchCITime=0.0;
   int totalFetches=0;
   long allcount=0;
   long block_status;

  int didblock=0;
  long totalblock=0;
  int excitblock=0;
  int phaseblock=0;
  int possibleblock=0;

   
   pair<int,double> info1;
   
   // Need to loop over the H, fetch it into sefsef blocks and pass that to the below method
   // Much like was done in actually constructing H
   
#ifdef NEWORBMAP
   GA::error(" do not use constructGlobalHamiltonianAltOrbMapChunkJmap method with NEWORBMAP",-1);
#endif
   
#ifdef SUPPLEMENTAL
   if ( g_rank == 0 ) cout << " Using SUPPLEMENTAL storage model For CI density" << endl;
   initSuppVal();
#endif
   
#ifdef STATIC 
   if ( g_rank == 0 ) cout <<" CIdensity: with NONDIRECT and STATIC"<< endl;
#endif

#ifdef REVERSE
  if ( g_rank == 0 ) cout << "NONDIRECT: CIdensity Using reverse iconf indexing " << endl;
#endif

#ifndef STATIC
   if ( g_rank == 0 ) cout <<" CIdensity: NXTVAL" << endl;
#endif

   int ilo;
   int ihi;
   
   if ( g_rank == 0 ) cout << " simply fetch CI dets " << endl;
    
   /* NOTE: THis approach is much sparser than the original H matrix.
      so we may want to consider modification of chunk sizes
   */
   
#ifdef REPLICATEDETVALUE
  if ( g_rank == 0 ) cout << "USING NEW JPHASE code for CIDEN build"<<endl; ;
   vector<char> * rep_phases = l_deters->fetchReplicatedPhasesHandle();
 //  cout << "H build rep phase is size " << rep_phases->size()<<endl; // returns iopen
#endif

   //CHECK ME flushH(); // this simply zeros out cimat,icicol, num and the supplimental spaces
   
  //pair<int,double> info1;

   vector<int> list;
   int currentrun=0, countrun=0;
   
   double jfetchTime = 0.0;
   const int numjchunks=l_numJchunksHbuild; // carve up jmap space into 10 chunks
   
   int chunk_width;
   int l_maxspatials = l_deters->fetchMaxSpatials();
   
   (l_maxspatials % numjchunks==0)? chunk_width = l_maxspatials / numjchunks: chunk_width = 1 + (l_maxspatials / numjchunks);
   if ( GA::nodeid() ==0 ) cout << "Jchunk stats num: "<<numjchunks<<" width " << chunk_width << endl;
   int jchunk_index=0; // these track the actual real index numbers (0,maxefs-1)
   
   for(int iroot=0; iroot < numRoots; ++iroot )
     {
       vector<COMPRESS_SIZE> imap( l_nbf+2  );
       vector<COMPRESS_SIZE> jmap( l_nbf+2  );
       
       JOUTFG fetchi;
       JOUTFG fetchj;
       // vector<double> vec_i; // Gets resized in fechFullCIvector
       
       countrun=0;
       for(int jchunk=0; jchunk < numjchunks; ++jchunk ) { // data parallel
	 
	 double tempj = psociTime();
	 vector<COMPRESS_SIZE> jsubsetlist; // jmap subset buffer gets sized in fetchJsubsetlist 
	 vector<int> jindexlist;
	 int jnumlist = fetchJsubsetlist( jchunk, jindexlist, jsubsetlist, chunk_width, l_maxspatials ); //GA on the inside. resize width inside to be perfect
	 jfetchTime += psociTime() - tempj;
	 
	 // start original code
	 
	 PsociTaskManager jobs(l_nxtval_chunk, g_nxtval );
	 ilo = 0;
	 ihi = l_maxspatials;
#ifndef STATIC
#ifdef REVERSE
         jobs.initNxtask( l_maxspatials );
#else
         jobs.initNxtask();
#endif

#ifdef REVERSE
          long nexttask = jobs.nxtaskRev( g_size, l_maxspatials );
#else
          long nexttask = jobs.nxtask( g_size ); // Keep J's agglomerated and dynamically balance I for now
#endif
#endif

#ifndef STATIC
#ifdef REVERSE
         for( int iconf = ihi-1; iconf >= ilo; --iconf ) {
#else
         for( int iconf = ilo; iconf < ihi; ++iconf ) {
#endif
#endif
#ifdef STATIC
#ifdef REVERSE
         for( int iconf = ihi-1-g_rank; iconf >= ilo; iconf-=g_size ) {
#else
         for( int iconf = ilo+g_rank; iconf < ihi; iconf+=g_size ) {
#endif
#endif
           totalblock += (iconf+1);

#ifndef STATIC
	   if ( nexttask == iconf ) { // in the loop stride is implemented  
#endif
	     
	     //	  vector<COMPRESS_SIZE> imap( l_nbf+2  );
	     l_deters->fetchOrbMap( iconf+1, imap);
	     
	     temp2 = psociTime();
             l_deters->fetchOrbMap( iconf+1, imap);
	     int currentrun =  l_deters->compareLocalJSubsetOrbMapCIDEN(jnumlist, iconf+1, rep_phases, imap, jindexlist, jsubsetlist, list); // 
	     countrun += currentrun;
	     orbmapTime += psociTime() - temp2;
	     
	     if ( currentrun != 0 ) {
	       
               l_deters->fetchAndUnpackDeterminantData( iconf+1, info1, fetchi );
	       fetchTime += psociTime() - temp2;
	       ++totalFetches;
	       allcount += iconf+1;
	       
	       vector<int>::iterator jit;
	       for(jit = list.begin(); jit != list.end(); ++jit ) 
		 {
                   ++didblock;
		   int jconf = (*jit);
		   temp2 = psociTime();
                   l_deters->fetchAndUnpackDeterminantData( jconf+1, info1, fetchj );
		   fetchTime += psociTime() - temp2;
		   ++totalFetches; // comment me out eventually
		   
		   temp2 = psociTime();
		   block_status = l_gaHamiltonian->generateSOCIblocksForCIDEN( fetchi, fetchj, sefsef );
		   sociTime += psociTime() - temp2;
		   if ( block_status ) 
		     { // return 0 or 1 where 1 is GOOD
		       temp3 = psociTime();
		       //  int new_stat =  transformSEFSEFtoMOnoGET(fetchi,fetchj,sefsef, vec_i, density[iroot]);
		       int new_stat =  transformSEFSEFtoMOGET(iroot, fetchi, fetchj,sefsef, g_v, density[iroot]);
		       transTime += psociTime() - temp3;
		     } 
		 }
	     } // currentrun 
#ifndef STATIC
#ifdef REVERSE
             nexttask = jobs.nxtaskRev( g_size, l_maxspatials );
#else
             nexttask = jobs.nxtask( g_size );
#endif
	   } // inext
#endif
	 } // iconf
#ifndef STATIC
	 jobs.destroyNxtask(); //Close'er up
#endif
       } // jchunk
     } // new iroot
   
   if ( GA::nodeid() == 0 ) {
     cout << GA::nodeid() << " nondirect CIDEN sociTime is " << sociTime<<endl;
     cout << GA::nodeid() << " nondirect CIDEN fetchTime is " << fetchTime<<endl;
     cout << GA::nodeid() << " countrun is " << countrun << endl;
     cout << GA::nodeid() << " fetchCIvectime is " << fetchCITime << endl;
     cout << GA::nodeid() << " transTime is " << transTime << endl;
     cout << GA::nodeid() << " jfetchTime is " << jfetchTime << endl;
     cout << GA::nodeid() << " orbmapTime is " << orbmapTime << endl;
     cout << GA::nodeid() << " didblock is " << didblock<<" "<<totalblock << endl;
     cout << "CIDEN allcount is " << allcount<<endl; //this is correct
     cout << "total CIDEN Fetches is " << totalFetches<<endl; // Print all of them for now since they indicate load balance
   }
   
#ifdef SUPPLEMENTAL
   destroySuppVal();
#endif
   
   GA::SERVICES.sync();

   return(0);
 }
 
void PsociGAhamiltonian::setNxtvalChunk( int nchunk )
{
    l_nxtval_chunk = nchunk;
}

/* NOT really local I can be any valid row
*/
void PsociGAhamiltonian::fetchLocalHblock( int i, vector<double> & cimat, vector<int> & icicol, int & number )
{
        int cilo[2];
        int cihi[2];
        cilo[0] = i;
        cihi[0] = i; // This may not apply to the supp space: need a check at creation
        cilo[1] = 0; 
        cihi[1] = 0;
        int n;

        n=1;
//cout << GA::nodeid() << " call number " << i << endl;
        l_g_number->get( cilo, cihi, &number, &n );
      
        if ( number == 0 ) { //Then this points to an absolute area in supplemental
        GA::error("fetchLocalHblock: ZERO ICICOL never supposed to happen ",i);
        }

        if ( number >= 0 ) {
// what happens to ListH is it is zero ?

#ifdef DETAILEDHAMCHECK
        if ( number > l_maxsparse ) {
          cout << "number exceeds specified maxsparse: Will not be able to fetch the main GA. Try supplemental ." << endl;
          cout << " number = " << number << " maxsparsity is " << l_maxsparse << endl;
        }
#endif
        cihi[1] = number - 1;
        n = number;

        l_g_icicol->get( cilo, cihi, &icicol[ 0 ], &n);
        l_g_cimat->get( cilo, cihi, &cimat[ 0 ], &n);

// else

        } else if ( number < 0 ) {
        int suppRow = -number - 1;

        cilo[0] = suppRow;
        cihi[0] = suppRow; // This may not apply to the supp space: need a check at creation

        l_g_number_supp->get( cilo, cihi, &number, &n );

#ifdef DETAILEDHAMCHECK
        if ( number > l_maxsparse_supp ) {
          cerr << "number exceeds specified maxsparse_supp: Will not be able to fetch the GA_supp." << endl;
          cerr << " number = " << number << " maxsparsity_supp is " << l_maxsparse_supp << endl;
          GA::error("number exceeds specified maxsparse_supp:", l_maxsparse_supp);
        }
#endif
        cihi[1] = number - 1;
        n = number;
        l_g_icicol_supp->get( cilo, cihi, &icicol[ 0 ], &n);
        l_g_cimat_supp->get( cilo, cihi, &cimat[ 0 ], &n);
        }
}

/* Wrapper around fetcn to get ALL local H related objects. THis is really only meant to be helpful
   at-scale

   fetch H elements clow,chi inclusive
   For now The buffers have been PREALLOCATED in the calling method
*/

void PsociGAhamiltonian::fetchAllLocalHblock( int clow, int chi, vector<vector<double> > & local_cimat, vector<vector<int> > & local_icicol, vector<int> & local_number )
{
//     int width = chi - clow + 1;
     int index=0;
     for(int i=clow; i<=chi; ++i )  {
        fetchLocalHblock( i, local_cimat[index], local_icicol[index], local_number[index] );
        ++index;
     }
}

/* Wrapper around fetcn to get ALL (not necc local)  H related objects. 

   Grab all guys in the preestablished list
   For now The buffers have been PREALLOCATED in the calling method
*/

void PsociGAhamiltonian::fetchListHblock( vector<int> & list, vector<vector<double> > & local_cimat, vector<vector<int> > & local_icicol, vector<int> & local_number )
{
     int index=0;
     vector<int>::iterator it;
#pragma ivdep
     for(it=list.begin(); it!=list.end(); ++it  )  {
        fetchLocalHblock( (*it), local_cimat[index], local_icicol[index], local_number[index] );
        ++index;
     }
}


/* Determine the set of current usabvle icicol between start and end-1
   Exploit the orde rof icicol.
*/
int PsociGAhamiltonian::fetchListWithinRange( int start, int end, vector<int> & icicol, int icinum, vector<int> & list )
{
    int nsize = icinum;
    list.clear();
    for(int i=0; i< nsize; ++i ) {
       int jhi = icicol[i];
       if ( jhi >= start && jhi < end ) list.push_back(jhi);
    }
    return ( list.size() );
}
// Everyone needs to do this
void PsociGAhamiltonian::setVectorChunk( int nchunk )
{
     if ( GA::nodeid() == 0 ) l_nchunk = nchunk;
     GA::brdcst( &l_nchunk, sizeof(int), 0) ;
}

void PsociGAhamiltonian::printVectorChunk()
{
     if ( GA::nodeid() == 0 ) cout << "Current Vector Chunk is " << l_nchunk << endl;
}

/* Since the M*V is rather SPSD, there is a great change that with eveyine fetching g_vec, cores will all
   be hitting the same GA core at the same time. Since the operationin M*V commute w.r.t this we can simply 
   generate staggered lists which will help eleviate the queueing.

   for(int chunk=0; chunk<nchunk; ++chunk )
   list == (0,nchunk-1)

*/
void PsociGAhamiltonian::generateStaggeredChunkList(int nodeid, int nchunk )
{
     my_chunklist.clear(); 

#ifdef RANDOM 
     srand( nodeid ); // not really needed
     int index = rand() % nchunk; // between 0 and ihi
#else
     int index = nodeid % nchunk;
#endif

#ifndef STAGGERCHUNKLIST
     index = 0; // override and make the list not-staggered; we then use brdcsts in the calling app
#endif

     for (int i=0; i< nchunk; ++i ) {
         my_chunklist.push_back( index );
         ++index;
         if ( index >= nchunk ) index = 0;
     } 

#ifdef DETAILEDCHECK
for(int i=0; i< nchunk; ++i ) {
  cout << GA::nodeid() << " list " << GA::nodeid() << " index " << my_chunklist[i] << endl;
}
#endif
}

void PsociGAhamiltonian::destroyStaggeredChunkList()
{
  my_chunklist.clear();
}

       
/* It is possible some core could have no local data 
   Also nothing precludes you from grabbing H data that is not local
   -- it is just much more expensive to do so --
*/
void PsociGAhamiltonian::fetchLocalHamiltonianData()
{
  int clow = local_cilo[0];
  int chi  = local_cihi[0];
  int width = chi - clow + 1;

  local_icicol.resize( width, vector<int>( max(l_maxsparse,l_maxsparse_supp), 0 ) );
  local_cimat.resize( width, vector<double>( max(l_maxsparse,l_maxsparse_supp), 0.0 ));
  local_number.resize( width, 0 );

  fetchAllLocalHblock( clow, chi, local_cimat, local_icicol, local_number );
}

/* construct a list of indexes each core will use to process H*v elements
*/
void PsociGAhamiltonian::generateListHamiltonianData()
{
  GA::SERVICES.sync(); // must carve up configuration space pretty evenly

  cimat_list.clear();
  int me = GA::nodeid();
  int size = GA::nodes();

// for(int i=local_cilo[0]; i<=local_cihi[0]; ++i ) {
// for(int i=local_cihi[0]; i>=local_cilo[0]; --i ) {

//    for(int i=me; i< l_maxsef; i+=size) {
//     for(int i=l_maxsef-1-me; i >= 0; i-=size) {

  for(int i=me; i< l_maxsef; i+=size) {
       cimat_list.push_back( i ) ;
  }
}


/* 
  Populate local data structures with the lisat of H data
  Collective 

  Ultimately construct in reverse
*/
void PsociGAhamiltonian::fetchListHamiltonianData()
{

  int width = cimat_list.size();
  if ( width <= 0 ) GA::error("fetchListHamiltonianData unexpected list size",width );

//  if ( width <= 0 ) GA::error(" width is zero " ,size );

//  local_icicol.resize( width, vector<int>( max(l_maxsparse,l_maxsparse_supp), 0 ) );
//  local_cimat.resize( width, vector<double>( max(l_maxsparse,l_maxsparse_supp), 0.0 ));
//  local_number.resize( width, 0 );

  fetchListHblock( cimat_list, local_cimat, local_icicol, local_number );

// NEED THIS to process H*v  list.clear();
}

/* 
  Populate local data structures with the lisat of H data
  Collective 

  Ultimately construct in reverse
*/
void PsociGAhamiltonian::fetchsubListHamiltonianData(vector<int> & list )
{  

  int width = list.size();;
  if ( width <= 0 ) GA::error("fetchsubListHamiltonianData unexpected list size",width );
   
//  if ( width <= 0 ) GA::error(" width is zero " ,size );

//  local_icicol.resize( width, vector<int>( max(l_maxsparse,l_maxsparse_supp), 0 ) );
//  local_cimat.resize( width, vector<double>( max(l_maxsparse,l_maxsparse_supp), 0.0 ));
//  local_number.resize( width, 0 );

  fetchListHblock( list, local_cimat, local_icicol, local_number );
}

void PsociGAhamiltonian::destroyLocalHamiltonianData()
{
    local_icicol.clear();
    local_cimat.clear();
    local_number.clear();
    cimat_list.clear();
}

/* new method for deflating multimap data assembled in the Gather H*v method

   memory operations are too expensive
   
*/
int PsociGAhamiltonian::deflateMultimapToVector( multimap< int, double > & map, vector<double> & data, vector<int> & indexarray, const int & ivec )
{

//    if ( map.size () == 0 ) return (0);

    int number = 0;

// Loop overmap similar index pairs get explicitely summed

   multimap< int, double >::iterator in = map.begin();
   multimap< int, double >::iterator out;
   multimap< int, double >::iterator sumit;

   double dvalue;

   int offset;
   int index = 0;

   while ( in != map.end() ) {
     out = map.upper_bound( (*in).first ); // explicit sums coming 
      dvalue=0.0;
      for( sumit=in; sumit !=out; ++sumit) {
         dvalue += (*sumit).second;
      }

/* check for zero values... actually not needed anymore since the
   calling app also dies this
*/
   if ( abs(dvalue) > MIN_HVEC_TOKEEP ) { 
     offset = 2*index;
     data[index] = dvalue;
     indexarray[offset] = (*in).first ;
     indexarray[offset+1] = ivec; 
     ++index;
     ++number;
   }
   in = out;
   }
//cout << "deflate % kept is " << (100*number)/map.size()  << endl;
   return( number );
} 


/*
   Here we need to sort manually and then accumulate

   NOT USED ANYMORE

*/
int PsociGAhamiltonian::deflateListToVector( list<VECTRIPLE> & list_buf, vector<double> & data, vector<int> & indexarray)
{
    list_buf.sort( compare_vectortriplets );
    int initialsize = list_buf.size();
    return ( initialsize );
}

/*

   NOT USED ANYMORE

*/
int PsociGAhamiltonian::deflateVectorToVector( vector<VECTRIPLE> & vector_buf, vector<double> & data, vector<int> & indexarray)
{
    sort( vector_buf.begin(), vector_buf.end(), vector_vectortriplets() );
    int initialsize = vector_buf.size();
    return ( initialsize );
}


/* VECTRIPLE
   row == jhi
   col == ivec // by construction this is always the same
   value == data
*/
bool compare_vectortriplets( const VECTRIPLE & left, const VECTRIPLE & right )
{
      return left.row < right.row; 
}

bool compare_vectorpairs( const pair<int,double> & left, const pair<int,double> & right )
{
      return left.first < right.first;
}


/* 
   New method to gather vector elements for all rows specified in list
   for vector ivec.

   This requires that the cimat/icicol/number block has already been prefetched
   must have preallocated vectorts to be
        vectors[MAX_AGGREGATE][ max(l_maxsparse,l_maxsparse_supp)]
        data[ max(l_maxsparse,l_maxsparse_supp) ]
        scratch [ 2* max(l_maxsparse,l_maxsparse_supp) ]


*/
int PsociGAhamiltonian::gatherVectorSet(GA::GlobalArray * g_v, const int & ivec, vector<int> & list, vector<double> & vectors,  vector<double> & scratch, vector<int> & indexarray )
{
  int num = 0;
  int width = list.size();
  int index=0;

//  vector<double> testVector( max(l_maxsparse,l_maxsparse_supp) ,0.0 );
/*
  if ( width <= 0 ) return(0);
  if ( width > vectors.size() ) GA::error("width > vectors.size()",width);
*/

  for(int i=0; i< width; ++i ) {
   int number = local_number[ i ];
   for (int j=0; j<number; ++j) {
     int offset = j*2;
      indexarray[offset] = local_icicol[i][j];
      indexarray[offset+1] = ivec; 
   }
   //double temp=psociTime();
   // NGA_Gather_flat(g_v->handle(), (void *)&scratch[0], (int *)&indexarray[0], number);
   NGA_Gather_flat(g_v->handle(), (void *)&vectors[index], (int *)&indexarray[0], number);
   index += number;  
   //getVTime += psociTime() - temp;
   num +=  number;
   }
   return( num);
}

/* pull memory registration out of the matrix vector products

   CIMAT and friends nver get wider than MAX_AGGREGATE
*/
void PsociGAhamiltonian::allocateLocalScratchMatrixVector() 
{
    int width = MAX_AGGREGATE;
#ifndef NEWACCUMULATEDSCATTER
    bufV.resize( max(l_maxsparse,l_maxsparse_supp), 0.0  );
#endif
// dataarray and indexarray are potentially used for scratch in two methods. Need to ensure size works for both

   int basicsize = max(l_maxsparse,l_maxsparse_supp);
   int dsize = max( width * basicsize, MAX_DEFLATE_ELEMENTS+basicsize );
   int isize = max( width * 2 * basicsize, MAX_DEFLATE_ELEMENTS+basicsize );

// MAX_DEFLATE_TRIGGER can be low by as much as a row
/* MAY cause failure with older scatter-based routines

   dataarray.resize( dsize,0.0  ); // memory registration too expensive to update in deflation
   indexarray.resize( isize,0 ); // maximum poss
*/

// THe following is only needed for the newest methods

    vectors.resize( width * max(l_maxsparse,l_maxsparse_supp) );
}

void PsociGAhamiltonian::destroyLocalScratchMatrixVector()
{
   dataarray.clear();
   indexarray.clear();
   vectors.clear();
   bufV.clear();
   bufV2.clear();
   bufVIndex.clear();
   bufVnumber.clear();
}

/* pre allocate one or more rows worth of cimat,icicol, num
*/

void PsociGAhamiltonian::allocateLocalCIMATarray()
{
#ifdef MATRIXBALANCE
  int width = MAX_AGGREGATE;
#else
  int width = 1; // if noit maxtrixbalance then we simply do local reads
#endif
  local_icicol.resize( width, vector<int>( max(l_maxsparse,l_maxsparse_supp), 0 ) );
  local_cimat.resize( width, vector<double>( max(l_maxsparse,l_maxsparse_supp), 0.0 ));
  local_number.resize( width, 0 );
}

/* experimental accumu;ate method
accumulateAndPutProductChunk((*jit), ivec, jchunk_width, v_lo, v_hi, bufV, g_Hv );
By construction bufV and testVector are distributed the same
*/
void PsociGAhamiltonian::accumulateAndPutProductChunk( int jchunk, int ivec, int width, int * v_lo, int * v_hi, vector<double> & bufV, GA::GlobalArray * g_Hv )
{ 
  int l_maxsefs = local_dims[0];
  int my_vector_low = v_lo[0];
  int my_vector_hi = v_hi[0]; // this does not need subtraction by 1

   if ( width > bufV.size() ) GA::error(" width > bufV.size()",width);
#ifdef USEMPI
  MPI_Comm mpi_comm; // for comparing GAdgop to MPI dgop
  mpi_comm = GA_MPI_Comm();
#endif

#ifdef STAGGERCHUNKLIST
  GA::error(" Cannot use accumulateAndPutProductChunk method with staggered lists",1);
#endif

  int jchunk_lo[2];
  int jchunk_hi[2];
  
  int newlo[2];
  int newhi[2];
  double dAlpha=1.0;

  newlo[1] = ivec;
  newhi[1] = ivec;
  jchunk_lo[1] = ivec;
  jchunk_hi[1] = ivec;

  int ch_start = width*(jchunk);
  int ch_end = min( ch_start + width, l_maxsefs );
  int n = 1;
  char op='+';

  jchunk_lo[0] = ch_start;
  jchunk_hi[0] = ch_end-1; // this does need subtraction by 1

#ifdef GETANDBRDCSTHVEC
// tested and worls well on Lonestar
#ifdef USEMPI
  MPI_Allreduce( MPI_IN_PLACE, &bufV[0], width, MPI_DOUBLE, MPI_SUM, mpi_comm);
#else
  GA::dgop( &bufV[0], width, &op);
#endif
 if (GA::nodeid() == 0 ) g_Hv->acc( jchunk_lo, jchunk_hi, &bufV[0], &n, &dAlpha);

#endif

#ifndef GETANDBRDCSTHVEC
// Required for Hopper
  vector<int> l_list;
  for(int i=ch_start; i< ch_end; ++i ) { // loop entire actual chunk range
     if ( i >= my_vector_low && i <= my_vector_hi ) l_list.push_back( i ); // do I have any of it.
  }
// assemble others
#ifdef USEMPI
  MPI_Allreduce( MPI_IN_PLACE, &bufV[0], bufV.size(), MPI_DOUBLE, MPI_SUM, mpi_comm);
#else
  GA::dgop( &bufV[0], bufV.size(), &op);
#endif

//careful here: we do not want to doulbe acc on values
  if ( l_list.size() > 0 ) {
     newlo[0] = l_list[0];
     newhi[0] = l_list[ l_list.size() - 1 ];
     int index = newlo[0] - ch_start;
     g_Hv->acc( &newlo[0], &newhi[0], &bufV[index], &n, &dAlpha );
  }
#endif

  return;
}
/* encapsulate the process to find a chunk or the vector ivec. have each core fetch its local contribution
   followed by a GOP

   not sure yet what to do about zeroing local space for now let's do it.
   For now this simply replicates what was known to work.
*/
void PsociGAhamiltonian::fetchAndReplicateVectorChunk( int ichunk, int ivec, int width, int * v_lo, int * v_hi, vector<double> & chunk, GA::GlobalArray * g_v )
{ 
  int l_maxsefs = local_dims[0];
  int my_vector_low = v_lo[0];
  int my_vector_hi = v_hi[0]; // this does not need subtraction by 1

   if ( width > chunk.size() ) GA::error(" width > chunk.size()",width);
#ifdef USEMPI
  GA::error("fetchAndReplicateVectorChunk: cannot be used with USEMPI",1);
#endif

  int chunk_lo[2];
  int chunk_hi[2];
  
  int newlow[2];
  int newhi[2];

  newlow[1] = ivec;
  newhi[1] = ivec;
  chunk_lo[1] = ivec;
  chunk_hi[1] = ivec;

  int ch_start = width*(ichunk);
  int ch_end = min( ch_start + width, l_maxsefs );
  int n = 1;
  char op='+';

  chunk_lo[0] = ch_start;
  chunk_hi[0] = ch_end-1; // this does need subtraction by 1

#ifdef GETANDBRDCST

#ifdef STAGGERCHUNKLIST
        g_v->get( chunk_lo, chunk_hi, &chunk[0], &n );
#else
        if ( GA::nodeid() == 0 ) { 
               g_v->get( chunk_lo, chunk_hi, &chunk[0], &n );
            }
        GA::brdcst( &chunk[0], chunk.size()*sizeof(double), 0);
#endif

#endif


/* We need an alternative approach. Each core should compare the requested chunk data
   to its own local data. THen fetch its local contribution into the local buffer. 
   THis is followed by a GOP which should be pretty fast IF the data object is not too big.

   Note we must worry about inserting local data into the right absolute location of the chunk.
*/
   
#ifndef GETANDBRDCST

#ifdef STAGGERCHUNKLIST
  GA::error(" Cannot use new gop method with staggered lists",1);
#endif
  zeroVector( chunk ); // find a way around this

  vector<int> l_list;
  for(int i=ch_start; i< ch_end; ++i ) { // loop entire actual chunk range
     if ( i >= my_vector_low && i <= my_vector_hi ) l_list.push_back( i ); // do I have any of it.
  }

//cout << "NO GETANDFBRD " <<  my_vector_low << " " << my_vector_hi << " " << l_list.size() << endl;
  if ( l_list.size() > 0 ) {
     newlow[0] = l_list[0];
     newhi[0] = l_list[ l_list.size() - 1 ];
     int index = newlow[0] - ch_start;
//cout << "NO GETANDFBRD " <<  my_vector_low << " " << my_vector_hi << " " << l_list.size() << " chunk size " << chunk.size() << " " << index <<endl;
     g_v->get( &newlow[0], &newhi[0], &chunk[index], &n );
  }
  GA::dgop( &chunk[0], chunk.size(), &op);

#endif

  return;
}

void PsociGAhamiltonian::fetchAndReplicateVectorChunkMPI(int ichunk, int ivec, int width, int * v_lo, int * v_hi, vector<double> & chunk, GA::GlobalArray * g_v )
{ 
#ifndef USEMPI
   GA::error("fetchAndReplicateVectorChunkMPI: MUST compile with USEMPI ",1);
#endif

#ifdef USEMPI
  MPI_Comm mpi_comm; // for comparing GAdgop to MPI dgop
  mpi_comm = GA_MPI_Comm();
#endif

  int l_maxsefs = local_dims[0];
  int my_vector_low = v_lo[0];
  int my_vector_hi = v_hi[0]; // this does not need subtraction by 1

  if ( width > chunk.size() ) GA::error("MPI: width > chunk.size()",width);

  int chunk_lo[2];
  int chunk_hi[2];
  int newlow[2];
  int newhi[2];

  newlow[1] = ivec;
  newhi[1] = ivec;
  chunk_lo[1] = ivec;
  chunk_hi[1] = ivec;

  int ch_start = width*(ichunk);
  int ch_end = min( ch_start + width, l_maxsefs );
  int n = 1;
  char op='+';

  chunk_lo[0] = ch_start;
  chunk_hi[0] = ch_end-1; // this does need subtraction by 1

  //int cwidth = ch_end-ch_start;
   int cwidth=chunk.size();// in case not everyone has the same size width

#ifdef USEMPI

#ifdef GETANDBRDCST
  if ( GA::nodeid() == 0 ) { 
         g_v->get( chunk_lo, chunk_hi, &chunk[0], &n );
  }
  int iroot=0;
  MPI_Bcast( &chunk[0], cwidth, MPI_DOUBLE, iroot, mpi_comm );
#endif


#ifndef GETANDBRDCST
  zeroVector( chunk ); // find a way around this

  vector<int> l_list;
  for(int i=ch_start; i< ch_end; ++i ) { // loop entire actual chunk range
     if ( i >= my_vector_low && i <= my_vector_hi ) l_list.push_back( i ); // do I have any of it.
  }
  if ( l_list.size() > 0 ) {
     newlow[0] = l_list[0];
     newhi[0] = l_list[ l_list.size() - 1 ];
     int index = newlow[0] - ch_start;
     g_v->get( &newlow[0], &newhi[0], &chunk[index], &n );
  }
#ifdef USEMPI
  MPI_Allreduce( MPI_IN_PLACE, &chunk[0], cwidth, MPI_DOUBLE, MPI_SUM, mpi_comm);
#endif
#endif

#endif // USEMPI
  return;
}

void PsociGAhamiltonian::replicateAndPutHV( int ivec, int maxsefs, int * v_lo, int * v_hi, vector<double> & bufv, GA::GlobalArray * g_Hv ) 
{
  int lo[2], hi[2];
  int n = 1;
  char op ='+';
  int index = v_lo[0];

  lo[0] = v_lo[0];
  hi[0] = v_hi[0];
  lo[1] = ivec;
  hi[1] = ivec;

//cout <<" at dump "<<lo[0]<<" "<<lo[1]<<" "<<hi[0]<<" "<<hi[1]<<" "<<ivec<<" "<<maxsefs<<endl;
  GA::dgop( &bufV[0], maxsefs, &op);
  g_Hv->put( lo, hi, &bufV[index], &n);
  GA::SERVICES.sync(); // really needed ?

}


/* THis method requires CLOW AND CHI to have been defined
*/
void PsociGAhamiltonian::resizeBufferScatterMethod()
{
#ifndef NEWACCUMULATEDSCATTER
  GA::error(" resizeBufferScatterMethod only usable with -DNEWACCUMULATEDSCATTER",-1);
#endif

  //cout << "ALLOCATE bufV2 with " << l_maxsparse<<" "<<l_maxsparse_supp <<endl;
  int clow = local_cilo[0];
  int chi  = local_cihi[0];

// chi-clow is way too big

/*
#ifdef MATRIXBALANCE
cout << "TEST size for MATRIXBALANCE " << endl;
   // int bufWidth = MAX_AGGREGATE;
    int bufWidth = chi - clow +1;
#else
    int bufWidth = chi - clow +1;
#endif
*/

  if ( GA::nodeid() == 0 ) cout << "OVERRIDE bufV2 size to " << MAX_AGGREGATE;
  int bufWidth = MAX_AGGREGATE;

  bufV2.resize( bufWidth );
  bufVIndex.resize( bufWidth );
  bufVnumber.resize( bufWidth );

  
  for(int i=0; i< bufWidth; ++i ) {
    bufV2[i].resize( max(l_maxsparse,l_maxsparse_supp ) );
    bufVIndex[i].resize( max(l_maxsparse,l_maxsparse_supp) );
  }
/*
  vector<vector<double> > bufV2(chi-clow+1,vector<double>( max(l_maxsparse,l_maxsparse_supp), 0.0  ) ); // catch a local I's worth of data then scatter
  vector<vector<int> > bufVIndex(chi-clow+1, vector<int>(  max(l_maxsparse,l_maxsparse_supp), 0  ) );
  vector<int> bufVnumber(chi-clow+1,0 );
*/
}

/* Populate the dataarray and indexarray only with non-zero controbutions. THis is ultimately for feeding back to ga scatter
   For now open up a new vector array type 
   We zero out buf here since the data are alrerady in cache 
*/
int PsociGAhamiltonian::packDataNoZeros(int ivec, vector<int> & number,  vector<vector<int> >& index, vector<vector<double> >& buf, vector<int> & indexarray, vector<double> & dataarray ) 
{
    indexarray.clear();
    dataarray.clear();
    const double dZero = 0.0;
    vector<pair<int,double> > combined;
    int itotal = 0;

    for(int i=0; i< buf.size(); ++i ) {
       int num = number[i];
       for(int j=0; j< num; ++j ) {
          ++itotal;
          if ( abs(buf[i][j]) > dZero ) {
            combined.push_back( pair<int,double>( index[i][j], buf[i][j] ) );
            buf[i][j] = dZero; // need to do this anyway
          }
       }
    } 

//cout << "compressed ration " << combined.size() << " " << itotal << endl;

    sort( combined.begin(), combined.end(), compare_vectorpairs );


//cout << "Compare sizes " << buf.size() <<" "<<combined.size() << endl;
/* populate data array and index array - and remove zeros */

   int pushed_index=-1;
   double data, olddata ;
   int row, oldrow=-1 ;

   int iactual = 0;
   for(int i=0; i< combined.size(); ++i ) {
      row = combined[i].first;
      data = combined[i].second;

      if ( row == oldrow ) {
         dataarray[ pushed_index ] += data;
      } else {
      if ( abs(data) > 0.0 ) {
         dataarray.push_back( data );
         indexarray.push_back( row );
         indexarray.push_back( ivec );
         ++pushed_index;
         ++iactual;
      }
      }

      oldrow = row;
      olddata = data;
   }

// cout << "iactual" << iactual << " itotasl " << itotal <<  " % reduction is " << 100.0 * iactual / (double) itotal << endl;
  return ( dataarray.size() );
}

/* Populate the dataarray and indexarray only with non-zero controbutions. THis is ultimately for feeding back to ga scatter
   For now open up a new vector array type 
   We zero out buf here since the data are alrerady in cache 

   This build a TRANSPOSED and (+1) version of the array suitable for direct calls to pnga_scatter_acc
*/
int PsociGAhamiltonian::packDataNoZerosTransposed(int ivec, vector<int> & number,  vector<vector<int> >& index, vector<vector<double> >& buf, vector<int> & indexarray, vector<double> & dataarray ) 
{
    indexarray.clear();
    dataarray.clear();
    const double dZero = 0.0;
    vector<pair<int,double> > combined;
    int itotal = 0;

    for(int i=0; i< buf.size(); ++i ) {
       int num = number[i];
       for(int j=0; j< num; ++j ) {
          ++itotal;
          if ( abs(buf[i][j]) > dZero ) {
            combined.push_back( pair<int,double>( index[i][j], buf[i][j] ) );
            buf[i][j] = dZero; // need to do this anyway
          }
       }
    } 

//cout << "compressed ration " << combined.size() << " " << itotal << endl;

    sort( combined.begin(), combined.end(), compare_vectorpairs );

//cout << "Compare sizes " << buf.size() <<" "<<combined.size() << endl;
/* populate data array and index array - and remove zeros */

   int pushed_index=-1;
   double data, olddata ;
   int row, oldrow=-1 ;
   const int ndim = 2;

   int iactual = 0;
   for(int i=0; i< combined.size(); ++i ) {
      row = combined[i].first;
      data = combined[i].second;

      if ( row == oldrow ) {
         dataarray[ pushed_index ] += data;
      } else { dataarray.push_back( data );
       indexarray.push_back( ivec+1 ); // flip order
       indexarray.push_back( row+1 );
       ++pushed_index;
       ++iactual;
      }
      oldrow = row;
      olddata = data;
   }

// Should be fortran transposed for here

/*
   cout << "SIZE " << indexarray.size() << " data is " << dataarray.size() << endl;
   int ind=0;
   for(int j=0; j< 2*dataarray.size(); ++j ) {
   cout << "trans j " <<  indexarray[j] << endl;
   }
*/

// cout << "iactual" << iactual << " itotasl " << itotal <<  " % left is " << 100.0 * iactual / (double) itotal << endl;
  return ( dataarray.size() );
}


/* WARNING indexarray memory can easily be leaking here
   NOTE API has changed alot here we should rename it
*/
int** PsociGAhamiltonian::packDataNoZeros(int ivec, vector<int> & number,  vector<vector<int> >& index, vector<vector<double> >& buf, int & actualsize, vector<double> & dataarray )
{
// Double check that indexarray is empty ?
    dataarray.clear();

    const double dZero = 0.0;
    vector<pair<int,double> > combined;
    int itotal = 0;
    for(int i=0; i< buf.size(); ++i ) {
       int num = number[i];
       for(int j=0; j< num; ++j ) {
          ++itotal;
          if ( abs(buf[i][j]) > dZero ) {
            combined.push_back( pair<int,double>( index[i][j], buf[i][j] ) );
            buf[i][j] = dZero; // need to do this anyway
          }
       }
    }
   sort( combined.begin(), combined.end(), compare_vectorpairs );

   int pushed_index=-1;
   double data, olddata ;
   int row, oldrow=-1 ;

   int iactual = -1;
   int ** indexarray = new int*[ combined.size() ]; // Max number of entries

/* CONVERT TO FORTRAN indexing */
   
   int FORTSHIFT=0;
#ifdef WNGA
   FORTSHIFT=1;
#endif

   for(int i=0; i< combined.size(); ++i ) {

      row = combined[i].first;
      data = combined[i].second;

      if ( row == oldrow ) {
         dataarray[ pushed_index ] += data;
      } else {
      if ( abs(data) > 0.0 ) {
         ++iactual;
         indexarray[iactual] = new int[2];
         dataarray.push_back( data );
         indexarray[iactual][0] = row;
         indexarray[iactual][1] = ivec;
         ++pushed_index;
       }
      }
      oldrow = row;
      olddata = data;
   }
  actualsize = dataarray.size();
  return( indexarray );
}

/* going to try and eliminate the SCATTER by doing a non-blocking loop of single element puts

*/
void PsociGAhamiltonian::LoopNputMethod(GA::GlobalArray * g_a, vector<int> & indexarray, vector<double> & values, GANbhdl * g_handle ) 
{
     int num = values.size();
     int numindex = indexarray.size();
     int n=1;
     double dalpha = 1.0;
     int lo[2],hi[2];

     if (numindex != 2*num ) GA::error(" array mismatch: LoopNputMethod",numindex);     

     int index=-1;
     for(int i=0; i<num; ++i ) {
         ++index;
         lo[0] = indexarray[index]; 
         hi[0] = lo[0];
         ++index;
         lo[1] = indexarray[index];
         hi[1] = lo[1];
         double dvalue = values[i];
//         dataarray.push_back( data );
//         indexarray.push_back( row );
//         indexarray.push_back( ivec );

         g_a->nbAcc( lo, hi, (void *)&dvalue, &n, (void *)&dalpha, g_handle );
//cout << "PUT "<<lo[0]<<" "<<hi[0]<<"  "<<lo[1]<<" "<<hi[1]<<" "<<dvalue<<endl;
      }
}


/* Temporary method for use with experimental matrixVector product scheme
   Need to include some better checking
*/
void PsociGAhamiltonian::emptyIndexarray( int num, int**data )
{
    for( int i=0; i< num; ++i ) {
       delete [] data[i];
    }
    delete [] data;
}

void PsociGAhamiltonian::indexToInteger( vector<int> & indexarray, vector<Integer> & indexarrayInteger )
{
     indexarrayInteger.clear();
     vector<int>::iterator it;
     for(it=indexarray.begin(); it != indexarray.end(); ++it ) {
        indexarrayInteger.push_back( *it );
     }
     indexarray.clear();
     return;
}

// This method simply  finds the specific set of J's belonging to the current set.
/* This list will be passed to the fetchGA methjods to actually get jmaps for comparison. 
   NOTE: using actual values here NOT Fortranized (+1) values 

   jsubsetlist are the actual ormaps for the current subset ( trangular packed )

   jsubsetlist has not been sized
   jlist contains the actual J index ( GA format )
*/

int PsociGAhamiltonian::fetchJsubsetlist( int jchunk, vector<int> & jlist, vector<COMPRESS_SIZE> & jsubsetlist, int chunk_width, int l_maxspatials )
{
//Determine list of actual Js to fetch

    jlist.clear();
    int jstart = jchunk * chunk_width;
    int jend = jstart + chunk_width - 1;
    if ( jend > l_maxspatials-1 ) jend = l_maxspatials - 1; // jend is INCLUSIVE in this case


    for(int j=jstart; j<=jend; ++j ) {
       jlist.push_back( j );
    }


    COMPRESS_SIZE czero=0;
    jsubsetlist.resize( jlist.size() * (l_nbf+2 ),czero ); 

// now fetch actual orbmaps it's okay to do a simple loop approach

    double tempt=psociTime();
    double comptime=0.0;

    int jnumlist = jlist.size();

    int lo[2], hi[2];
    vector<int>::iterator jit;
    int index=0;
    for(jit=jlist.begin(); jit!= jlist.end(); ++jit ) {
       int jindex = *jit;
       l_deters->fetchOrbMap( jindex+1, &jsubsetlist[index]  );
     //cout <<"partial subs " << jindex+1 << " " <<index<<" " <<jsubsetlist[index] <<" "<<jsubsetlist[index+1] << endl;
    index += l_nbf+2;
    }

    comptime = psociTime() - tempt;
  return( jlist.size() );
}

/* Combines calls to pack zero and permits accumu.ation of indexarray and dataarray
 *  * note if transpose then this also IMPLIES fortran indexing (1,2,n) and not C (0,1,n-1)
 *   *
 *    * the size of bufV2 is set in resizeBufferScatterMethod. THis is set to be MAX_AGGREGATE
 *     * that sould be okay since it is possible that many rows got inserted
 *     */
int PsociGAhamiltonian::packDataNoZerosAccumulate(bool transpose, int indexnewarray, int ivec, vector<int> & number,  vector<vector<int> >& index, vector<vector<double> >& buf, vector<int> & indexarray, vector<double> & dataarray )
{
    const double dZero = 0.0;
    vector<pair<int,double> > combined;
    int itotal = 0;
    int offs = 0;
    int jseftemp, ivectemp;
    int realbufsize = indexnewarray+1; // calling app starts at zero

    double temp, first=0.0, second=0.0, third=0.0;

    if ( dataarray.size() > 0 ) { // read pre-existing and cache to that we can sort later and merge values for the same destination 

       int length=dataarray.size(); // indexarray should be 2* this

       temp = psociTime();
       for(int k=0; k< length; ++k ) {

          offs = 2*k;

          double dval = dataarray[k];
          if ( transpose==true ) {
            ivectemp = indexarray[ offs ]-1;
            jseftemp = indexarray[ offs+1 ]-1;

          } else {
            jseftemp = indexarray[ offs ];
            ivectemp = indexarray[ offs+1 ];
          }
          combined.push_back( pair<int,double>( jseftemp, dval ) );
        }
        indexarray.clear();
        dataarray.clear();
        first = psociTime() - temp;

   }
   itotal = combined.size();
    if ( buf.size() < realbufsize ) GA::error("buf.size() < realbufsize",realbufsize);

    temp = psociTime();

    for(int i=0; i< realbufsize; ++i ) {
       int num = number[i];
       for(int j=0; j< num; ++j ) {
          ++itotal;
          if ( abs(buf[i][j]) > dZero ) {
            combined.push_back( pair<int,double>( index[i][j], buf[i][j] ) );
            buf[i][j] = dZero; // need to do this anyway
          }
       }
    }
    sort( combined.begin(), combined.end(), compare_vectorpairs );

   second = psociTime() - temp;

   int pushed_index=-1;
   double data, olddata ;
   int row, oldrow=-1 ;
   const int ndim = 2;

  temp = psociTime();

   int iactual = 0;
   for(int i=0; i< combined.size(); ++i ) {
      row = combined[i].first;
      data = combined[i].second;

      if ( row == oldrow ) {
         dataarray[ pushed_index ] += data;
      } else { dataarray.push_back( data );
       if ( transpose ) {
        indexarray.push_back( ivec+1 ); // flip order
        indexarray.push_back( row+1 );
       } else {
        indexarray.push_back( row );
        indexarray.push_back( ivec );
       }
       ++pushed_index;
       ++iactual;
      }
      oldrow = row;
      olddata = data;
   }

   third = psociTime() - temp;
//   cout << "NEW TIMERS " << first<<" "<<second<<" "<<third<<endl;
    //cout << "iactual" << iactual << " itotal " << itotal <<  " % left is " << 100.0 * iactual / (double) itotal << endl;
    return ( dataarray.size() );
    }


int PsociGAhamiltonian::packDataNoZerosAccumulateSingle(bool transpose, int indexnewarray, int ivec, int  number,  vector<int> & index, vector<double> & buf, vector<int> & indexarray, vector<double> & dataarray )
{
    const double dZero = 0.0;
    vector<pair<int,double> > combined;
    int itotal = 0;
    int offs = 0;
    int jseftemp, ivectemp;
    int realbufsize = indexnewarray+1; // calling app starts at zero

    double temp, first=0.0, second=0.0, third=0.0;

    if ( dataarray.size() > 0 ) { // read pre-existing and cache to that we can sort later and merge values for the same destination 
       int length=dataarray.size(); // indexarray should be 2* this
       temp = psociTime();
       for(int k=0; k< length; ++k ) {
          offs = 2*k;
          double dval = dataarray[k];
          if ( transpose==true ) {
            ivectemp = indexarray[ offs ]-1;
            jseftemp = indexarray[ offs+1 ]-1;

          } else {
            jseftemp = indexarray[ offs ];
            ivectemp = indexarray[ offs+1 ];
          }
          combined.push_back( pair<int,double>( jseftemp, dval ) );
        }
        indexarray.clear();
        dataarray.clear();
        first += psociTime() - temp;
   }
   itotal = combined.size();

   temp = psociTime();

   if ( buf.size() < realbufsize ) GA::error("buf.size() < realbufsize",realbufsize);

   for(int i=0; i< realbufsize; ++i ) {
       int num = number;
       for(int j=0; j< num; ++j ) {
          ++itotal;
          if ( abs(buf[j]) > dZero ) {
            combined.push_back( pair<int,double>( index[j], buf[j] ) );
            buf[j] = dZero; // need to do this anyway
          }
       }
    }

    sort( combined.begin(), combined.end(), compare_vectorpairs );

    second += psociTime() - temp;
    temp = psociTime();

   int pushed_index=-1;
   double data, olddata ;
   int row, oldrow=-1 ;
   const int ndim = 2;

   int iactual = 0;
   for(int i=0; i< combined.size(); ++i ) {
      row = combined[i].first;
      data = combined[i].second;

      if ( row == oldrow ) {
         dataarray[ pushed_index ] += data;
      } else { dataarray.push_back( data );
       if ( transpose ) {
        indexarray.push_back( ivec+1 ); // flip order
        indexarray.push_back( row+1 );
       } else {
        indexarray.push_back( row );
        indexarray.push_back( ivec );
       }
       ++pushed_index;
       ++iactual;
      }
      oldrow = row;
      olddata = data;
   }
   third += psociTime() - temp;

//    cout << "TIMERS " << first<<" "<<second<<" "<<third<<endl;
    //cout << "iactual" << iactual << " itotal " << itotal <<  " % left is " << 100.0 * iactual / (double) itotal << endl;
    return ( dataarray.size() );
    }

/* moving forward we want to account for boundary crossing and lcoality. but for now simply fetch the data
   and see if the direct matrix-vector concept works at all
*/
void PsociGAhamiltonian::fetchPatchVector( int ivec, int ch_start, int ch_end, vector<double> & testvector, GA::GlobalArray * g_v )
{
     int lo[2], hi[2];
     int istart = ch_start;
     int iend = ch_end;
     int n=1;

     lo[1]=ivec;
     hi[1]=ivec;
     lo[0]=ch_start;
     hi[0]= ch_end-1;
     
     g_v->get( lo, hi, &testvector[0], &n);
}

/* moving forward we want to account for boundary crossing and lcoality. but for now simply fetch the data
   and see if the direct matrix-vector concept works at all
*/
void PsociGAhamiltonian::fetchAllVector( int ivec, vector<double> & testvector, GA::GlobalArray * g_v )
{
     int lo[2], hi[2];
     int n=1;

     int l_maxsefs = local_dims[0];
     lo[1]=ivec;
     hi[1]=ivec;
     lo[0]=0;
     hi[0]=l_maxsefs-1;

     g_v->get( lo, hi, &testvector[0], &n);
}

/* estimate_elements is wrong
   bring in FULL testVector 

   return the estimate_elements  for this call
*/

 int PsociGAhamiltonian::matrixVectorUpdate(int iconf,  int & indexnewarray, int nsefi,  vector<int> & bufVnumber, vector<vector<int> >& bufVIndex, vector<vector<double> >& bufV2, 
					    vector<vector<int> > & icicol, vector<vector<double> > & cimat, vector<int> & number, vector<double> & testVector )
 {
   
   int estimate_elements=0;

/* ensure that bufV2 ( and friends) are big enough to complete
   this call. We want to defer scatter to a higher level scope
   the scatter routine will prune back bufV2 (etc) to MAX_AGGREGATE
   as needed
*/

/*
cout << "In matrixVectorUpdate" << endl;
cout << "Sizes " << cimat.size()<<" "<<icicol.size() << " " << number.size() << endl;
cout << "Bufsizes " << bufV2.size() << " " << bufVIndex.size() << " "<< bufVnumber.size() << endl;
cout << "testVec " << testVector.size()<<endl;
cout << "Ints " << iconf<<" "<<indexnewarray<<" "<<nsefi<<endl;
cout << "ON entry cimat [0][0] is " << cimat[0][0] << endl;
*/

//NOTE: a resize causes a copy to be put on the stack


   if ( indexnewarray+nsefi >= MAX_AGGREGATE-1 ) {
     cout << "Boost size of bufV2 " << indexnewarray+nsefi << endl;
     GA::error("bufV2 size boost disabled ",1);
     int newsize = indexnewarray+nsefi + 1; //need +1 since indexnewarray is trailing
     bufVIndex.resize( newsize );
     bufV2.resize( newsize );
     bufVnumber.resize( newsize );
     for(int i=indexnewarray+1; i< newsize ; ++i ) { // need to allocated both dimensions
       bufV2[i].resize( max(l_maxsparse,l_maxsparse_supp ) );
       bufVIndex[i].resize( max(l_maxsparse,l_maxsparse_supp) );
     }
   }

   int isefstart = l_nsef[ iconf]; // vector or H indexing ? 
   int isefend = l_nsef[iconf] + nsefi;

   for(int in=0; in< nsefi; ++in ) {
     ++indexnewarray; // increment for the next entry
     ++estimate_elements;

     int istart = in + l_nsef[ iconf ]; // real address for I
     int numb = number[in];

     //cout << "MVP update numb is " << numb << " istart is " << istart <<" iconf is " << iconf << endl; 
     //if ( numb == 0 ) cout << " numb == 0 " << endl;

    // if ( numb > 0 ) {
     double dtemp=0.0;

     for(int j=0; j< numb; ++j ) {
       int jhi = icicol[in][j];
       dtemp += cimat[in][ j ] * testVector[ jhi ]; //istart offsets over sefs for the given iconf
       //cout << "dtemp " << dtemp <<" "<<in<<" "<< j<<" "<<jhi<<" "<<cimat[in][ j ]<<" "<<testVector[ jhi] << endl;

     }

     bufVnumber[ indexnewarray ] = numb; // how many were found
     bufV2[ indexnewarray ][ numb-1 ] += dtemp;
     bufVIndex[ indexnewarray ][ numb-1 ] = icicol[in][numb-1]; // no sum it is always the same
     //cout << "update First " << numb<<" "<< in<<" "<< indexnewarray<<" "<<bufV2[ indexnewarray ][ numb-1 ]<<" "<<bufVIndex[ indexnewarray ][ numb-1 ] << endl;


//bufVIndex is not always correct 
     
     double testVec = testVector[ istart ]; //need ith coef for the current reduction

     for(int j=0; j< numb-1; ++j ) {
       int jhi = icicol[in][j];
       bufV2[ indexnewarray ][ j ] += cimat[in][ j ] * testVec;
       bufVIndex[ indexnewarray ][ j ] = jhi; // no need to sum
       //cout << "update Second " << in<<" "<<j<<" "<<indexnewarray<<" "<<bufV2[ indexnewarray ][ numb-1 ]<<" "<<bufVIndex[ indexnewarray ][ j ] << endl;
     }
    // } // numb
   } // in loop
   return( estimate_elements );
 }
 
/* dump data to scatter
*/

 int PsociGAhamiltonian::matrixProductScatter( int iroot, int & indexnewarray, int & estimate_elements, vector<int> & bufVnumber, vector<vector<int> > & bufVIndex, vector<vector<double> > & bufV2, vector<int> & indexarray, vector<double> & dataarray, GA::GlobalArray * g_Hv )
 {

#ifdef TACC
   bool transpose=true;
#else
   bool transpose=false;
#endif
   double alpha=1.0;
   
// do I need this>?
   dataarray.clear(); // double check

// packDataNoZerosAccumulate empties dataarray and indexarray internaly
   int newrealnumber = packDataNoZerosAccumulate(transpose,  indexnewarray, iroot, bufVnumber,  bufVIndex, bufV2, indexarray, dataarray );

//what to do about bufV2 not needed anymore

   if ( newrealnumber > 0 ) { // double check
  char * HvectorMessage = "PsociGAhamiltonian:matrixVectorProductsScatter";
  g_Hv->checkHandle( HvectorMessage );

#ifdef TACC
     pnga_scatter_acc( g_Hv->handle(), (void *)&dataarray[0], &indexarray[0], newrealnumber, &alpha); // low level call
#endif
#ifdef SCATFLAT
     NGA_Scatter_acc_flat(g_Hv->handle(), (void *)&dataarray[0], (int*)&indexarray[0], newrealnumber, &alpha);
#endif
   }
   
   if ( bufV2.size() != MAX_AGGREGATE ) {
     GA::error("bufV2 size down prohibited",1);
     cout << "reset bufV2 sizes " << MAX_AGGREGATE << endl;
     bufVIndex.resize( MAX_AGGREGATE );
     bufV2.resize( MAX_AGGREGATE );
     bufVnumber.resize( MAX_AGGREGATE );
   }

   estimate_elements=0;
   indexnewarray=-1;

   indexarray.clear(); // last step cleanup
   dataarray.clear();

   int hvcount = 1;
   return( hvcount );
 }
 
/* Use a simple STATIC-ONLY approach to compute
   H diag elements for the direct matrix product method
*/
long PsociGAhamiltonian::computeOnlyDiagElements()
{
#ifndef DIRECTMATRIXPRODUCTS
  GA::error(" PsociGAhamiltonian::computeOnlyDiagElements only for -DDIRECTMATRIXPRODUCTS",1);
#endif

  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  
  //cout << "I am " << GA::nodeid() << " at the hamiltonian " << endl;
  //GA::SERVICES.sync();
  
  int counter=0;
  long allcount=0;
  
  if ( g_rank == 0 ) cout << " Only Diag terms: Simple static load balancing scheme being used " << endl;

  int ilo=0;
  int ihi=l_maxspatials;

  
  long totalElems=0;
  
  // Build Hamiltonian Diagonals from scratch
  
  double temp, temp2;
  double sociTime=0.0, packTime=0.0, fetchTime =0.0, diagTime=0.0;
  int totalFetches=0;
  int block_status;
  
  JOUTFG fetchi, fetchj;
  vector<COMPRESS_SIZE> imap( l_nbf+2  ); // gets hidden below
  vector<COMPRESS_SIZE> jmap( l_nbf+2  );
  
  int me = GA::nodeid();
  int size = GA::nodes();
  for(int iconf = ilo+me; iconf < ihi; iconf= iconf+size ) {  
      
    temp2 = psociTime();
    l_deters->computeConfigsGA( iconf+1, imap, fetchi ); // also does a fetchOrb underneath
    fetchTime += psociTime() - temp2;
     ++totalFetches;
 
      // One configuration at a time for now
      vector<vector<double> > cimat(fetchi.nsefi );
      vector<vector<int> > icicol(fetchi.nsefi );
      vector<int> number(fetchi.nsefi );
      
      temp = psociTime();

      int jconf = iconf; // only want the diagonl terms

      l_deters->computeConfigsGA( jconf+1, jmap, fetchj ); // also does a fetchOrb underneath
      fetchTime += psociTime() - temp2;
      ++totalFetches;
	  
      temp2 = psociTime() ;
      block_status = l_gaHamiltonian->generateSOCIblocks( fetchi, fetchj, sefsef );
      sociTime += psociTime() - temp2;

      temp2 = psociTime();
      totalElems += packSefsef( fetchi, fetchj , sefsef, cimat, icicol, number );
      packTime += psociTime() - temp2;

      temp2 = psociTime();
      uploadCIMATdiagElemsOnlytoGA( l_nsef[ fetchi.index - 1 ] , cimat, icicol, number );
      diagTime += psociTime() - temp2;
  }

  // This stuff still has meaning to me in the stdout.
  cout << GA::nodeid() << "packTime is " << packTime<<endl;
  cout << GA::nodeid() << "sociTime is " << sociTime<<endl;
  cout << GA::nodeid() << "diagTime is " << diagTime<<endl;
  cout << GA::nodeid() << "fetchTime is " << fetchTime<<endl;
  cout << GA::nodeid() << " Done after diag loops " <<  totalElems << endl;
  GA::SERVICES.sync();
  
  //return( totalElems ); // switched for performance measurment
  return( numberNonzeroElements( totalElems ) );
}

/* Use iopen now fopr a check
*/
//Note also: No need to return Js that are > i sincew we only do triangular matrices
int PsociGAhamiltonian::assemblePhasePossibles( int iconf, int iphase, vector<int> & jphases, vector<char> * replicated_phases )
{
#ifdef SQUAREHAMILTONIAN
    GA::error(" PsociGAhamiltonian::assemblePhasePossibles not for SQUARE H",1);
#endif

    int iopen = iphase;
    int jopen;

    int numphase=replicated_phases->size();

    for( int i=0; i< numphase; ++i ) { // Loops over maxspatials
         jopen = replicated_phases->at(i);
//cout << "RAW PHASES " << iconf << " "<<iphase<<" "<<i<<" "<< replicated_phases->at(i)<<endl;
// if ( replicated_phases->at(i) == iphase && (i <= iconf) ) jphases.push_back( i ); 

         if ( (abs(iopen - jopen) <= 4) && (i <= iconf) ) jphases.push_back( i );    
    }

    return( jphases.size() );
}
    
//Note also: No need to return Js that are > i sincew we only do triangular matrices
int PsociGAhamiltonian::assemblePhasePossiblesCIDEN( int iconf, int iphase, vector<int> & jphases, vector<char> * replicated_phases )
{
#ifdef SQUAREHAMILTONIAN
    GA::error(" PsociGAhamiltonian::assemblePhasePossibles not for SQUARE H",1);
#endif

    int iopen = iphase;
    int jopen;
    int numphase=replicated_phases->size();

    for( int i=0; i< numphase; ++i ) { // Loops over maxspatials
         jopen = replicated_phases->at(i);
         if ( (abs(iopen - jopen) <= 2) && (i <= iconf) ) jphases.push_back( i );
    }
    return( jphases.size() );
}

// Fetch the ICICOL data row by row and determine the number of neighbors: Look only
// to the right ( or down )
// This could be as big as a few GB of data
// It is possible that we can not worry about direction but we simply want to know who was adjacent to me

// OPen a text file if size nsef by nsef (symmetric) with non zeros indicasted by 1s


   int PsociGAhamiltonian::generateNodeRanking() 
{
  string tempfile = "nodeRankSquare.txt";
  int g_rank = GA::nodeid();

  cout << g_rank <<  " Start generateNodeRanking: using bitset" << endl;
// NOTE the ORDER is reversed in the bit so "1" is last

// only node zero actually writes
  if ( g_rank == 0 ) {

    fstream olocalfile;
    const ios_base::openmode VECTOR_READ  = std::ios::in;
    const ios_base::openmode VECTOR_WRITE = std::ios::out ;
    olocalfile.open( tempfile.c_str(), VECTOR_WRITE );
    if (!olocalfile.is_open()) {
      cerr << "Write: Cannot open nodeRankSquare.txt write " << tempfile << " Aborting write " << endl;
      exit(1);
    } 

    
// Loop over number then fetch values. If number is <0 it reports the row in supp matrices
// to find relevant data

   int biglo[2], bighi[2];
   int cilo[2], cihi[2];
   int lo[2], hi[2];
   int stride = 1;

   const char zero='0';
   const char one='1';
   int width=max(l_maxsparse_supp, l_maxsparse);

   for(int i=0; i< l_maxsef; ++i ) {
      vector<int> icicol( width, 0 ); 
      lo[0] = i;
      hi[0] = i;
      lo[1]= 0;
      hi[1]= 0;
      int n=-1;

      l_g_number->get( lo, hi, &n, &stride );
      if ( n >= 0 ) {
         cilo[0] = i;
         cihi[0] = i;
         cilo[1] = 0;
         cihi[1] = n-1;
         l_g_icicol->get( cilo,cihi,&icicol[0],&stride);

      } else {
         int suppRow = -n - 1;
         cilo[0] = suppRow;
         cihi[0] = suppRow;
         cilo[1]=0;
         cihi[1]=0; 
         int nbig=-1;
         l_g_number_supp->get( biglo, bighi, &n, &stride );
         cilo[1] = n-1;
         cihi[1] = n-1;
         l_g_icicol->get( cilo,cihi,&icicol,&stride);
      }
      //vector<char> temp( l_maxsef,(char) zero );
      
      boost::dynamic_bitset<unsigned long long> temp(l_maxsef); // zeros by default
         for (int ibyte=0; ibyte<n; ibyte++ ) { 
	     int icol = icicol.at(ibyte);
             temp[icol]=1; // sets value to 1
         }
      int count=0;
      olocalfile << temp << endl; // note order will be flipped 

     } // loop over GA fetches
     olocalfile.close();


//  Read the data and build a full symmetric matrix - rewrite symmetrized matrix to disk to disk
//  We can never fit all data in memory: So go ahead and sweep through as required
//  Star on the bottom as that row is full by default: then update newrow cummulatively

     boost::dynamic_bitset<unsigned long long> row1(l_maxsef); // zeros by default
     boost::dynamic_bitset<unsigned long long> row2(l_maxsef); // zeros by default
     boost::dynamic_bitset<unsigned long long> newrow(l_maxsef); // zeros by default
  
// open the nodeRankSquare.txt and begin symmetrization
     ifstream ilocalfile1;
     ilocalfile1.open( tempfile.c_str(), VECTOR_READ );
     ifstream ilocalfile2;
     ilocalfile2.open( tempfile.c_str(), VECTOR_READ );
     string tempfileout="symmetrizedNodeRank.txt";
     olocalfile.open( tempfileout.c_str(), VECTOR_WRITE );

    if (!ilocalfile1.is_open() || !ilocalfile2.is_open()) {
      cerr << "Write: Cannot open nodeRankSquare.txt 1 or 2 for read " << tempfile << " Aborting write " << endl;
      exit(1);
    } else {
      cout << "tempfile is open " << endl;
    }
    if (!olocalfile.is_open()) {
      cerr << "Write: Cannot open nodeRankSquare.txt write " << tempfile << " Aborting write " << endl;
      exit(1);
    }

// Read file from the FRONT and then add all
// ROW 1 is the row to be corrected and symmetrized
// NOTE: because of the use of bitset, the entries are "flipped" left-right
// so entry #1 is the rightmost entry in the array. 

  for ( int i=0; i< l_maxsef; ++i) {
    ilocalfile1 >> row1;
    ilocalfile2.clear();
    ilocalfile2.seekg(0,ios::beg);
    for (int iskip = 0; iskip <= i; iskip++) {
      ilocalfile2.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    for ( int j=i+1; j< l_maxsef; ++j ) {
        ilocalfile2 >> row2;
        row1[j]=row2[i];
    }
    //cout << "SYM ROW 1: "<<row1<<endl;
    olocalfile << row1 << endl;
}
    olocalfile.close();
    ilocalfile1.close();
    ilocalfile2.close();
     
/*
 tempfileout now contains the symmetrized icimat matrix (but in flipped order)
 Now reprocess that data to generate the node order information
*/

   ilocalfile1.open( tempfileout.c_str(), VECTOR_READ );
   if (!ilocalfile1.is_open() ) {
     cerr << "Write: Cannot open symmetrized.txt 1 or 2 for read " << tempfileout << " Aborting write " << endl;
     exit(1);
   } else {
     cout << "tempfile is open " << endl;
   }

   vector<int> nodeRank( l_maxsef,0 );
   for( int inode=0; inode < l_maxsef; inode++ ) { // search rows
      boost::dynamic_bitset<unsigned long long> row(l_maxsef);
      ilocalfile1 >> row;
      set<char> neighbors;
      for (int j=0; j< l_maxsef; ++j ) {
          if ( row[j]==1 ) neighbors.insert(j);
      }
      nodeRank[inode]=neighbors.size() - 1; // remove the self-terms
   }
   ilocalfile1.close();

   string tempfileorder="nodeAdjacencyRank.txt";
   olocalfile.open( tempfileorder.c_str(), VECTOR_WRITE );
   fstream olocalfile4;
   string tempfile4="nodeAdjacencyRank.txt";
   olocalfile4.open( tempfile4.c_str(), VECTOR_WRITE );
     for (int i=0; i< l_maxsef; i++ ) {
         olocalfile4 << "inode="<<i<<" size="<<nodeRank[i]<<endl;
     }
   olocalfile4.close();


  } //me
  GA::SERVICES.sync();
  cout << "Finished with RANK node " << endl;
  return(0);
}
       
// A new version that is HEAVY on reading the GA row but doesn't store data to disk
// other than a temp row and the ranking matrix
// IN PARALLEL

// different value swhen don ein PARALLEL
//FIX ME

   int PsociGAhamiltonian::generateNodeRankingInMemory() 
{
  int g_rank = GA::nodeid();

  cout << "HERE at InMemory" << endl;
  const ios_base::openmode VECTOR_READ  = std::ios::in;
  const ios_base::openmode VECTOR_WRITE = std::ios::out ;

  vector<int> nodeRank( l_maxsef,0 );

  if ( g_rank == 0 ) {
	cout << g_rank <<  "Start generateNodeRankingInMemory" << endl;
  }

  int cilo = local_cilo[0];
  int chi = local_cihi[0]; // get local range, carve up work 

  cout << " I am " << g_rank<<" "<<cilo<<" "<<chi<<endl;
  vector<set<int> > icimat;
  set<int> totalinodes;
  vector<int> loopindex;

  double timein = psociTime();
     fetchLocalIcicolGArows(loopindex, totalinodes, icimat ); 
     cout << g_rank << " size of totalinodes " << totalinodes.size() << endl;
     cout << g_rank <<"size of loopindex size " << loopindex.size() << endl;
  double fetchTime = psociTime() - timein;

  cout << GA::nodeid() << "fetch time is " << fetchTime << endl;
  GA::SERVICES.sync();
  cout << g_rank << " Dne with fetch time " << endl;

  timein = psociTime();

 // for(int inode=0;inode<l_maxsef; ++inode) {
    
    cout << g_rank <<"total node list is " << totalinodes.size() << endl;
    cout << g_rank <<" loopindex size " << loopindex.size() << endl;

    set<int>::iterator iter;

    for (iter = totalinodes.begin(); iter != totalinodes.end(); ++iter) {
    int inode = *iter;
    int localcount=0;
    double timeouter=0;
    double timeinner=0;

    int index=0;

// this can be an arbitrary list
//    for(int i=cilo; i<=chi; ++i ) {
    for(int ilist=0; ilist< loopindex.size(); ++ilist ) {

      int i = loopindex[ilist];
      if ( i >= inode ) {

     int nsize = icimat[index].size(); // hopefull prefetches 
     if ( i != inode ) {
      double temp=psociTime();
      if ( icimat[index].find(inode) != icimat[index].end()) localcount++;
      timeouter += psociTime()-temp;
     } else { 
        double temp=psociTime();
       localcount += nsize;
       timeinner += psociTime()-temp;
     } //if else inode

    } // ony look right (or self)
      ++index;
    }
    //cout << g_rank <<" times inner outer "<< timeinner<<" "<<timeouter<<endl;
    nodeRank[inode]=localcount-1;
  } // inode

  double processTime = psociTime() - timein;
  cout << "process time " << processTime << endl;


  char op[] = "+";
  GA::igop( (int *) &nodeRank[0], l_maxsef, op );
  
  /*
  if ( GA::nodeid() == 0 ) {
     for(int k=0; k< l_maxsef; ++k ) {
      cout << k<<" node rank - 1 " << nodeRank[k]-1<<endl;
     }
  }
  */

//  cout << g_rank << " nodes " << nodeRank.size() << " "<< nodeRank[0]  << endl;
//  cout << g_rank << "twice  nodes " << nodeRank.size() << endl;

// Now output ranks to a file

   if ( g_rank == 0 ) {
     fstream olocalfile;
     string tempfile="nodeAdjacencyRank.txt";
     olocalfile.open( tempfile.c_str(), VECTOR_WRITE );
     if (!olocalfile.is_open()) {
      cerr << "Write: Cannot open nodeRankSquare.txt write " << tempfile << " Aborting write " << endl;
      exit(1);
    } else {
      cout << " file is open " << tempfile << endl;
    }
    
     for (int i=0; i< l_maxsef; i++ ) {
         olocalfile << "inode="<<i<<" size="<<nodeRank[i]<<endl;
     }
     olocalfile.close();
   }

  cout << "Finished with parallel RANK node " << endl;
  return(0);
}
       
/* New methods fore the INCORE symmetrization routine: Fetch a single row from GA
 * and return it as a dynamic_bitset operation
 * Requires file to have been previously opened
 *
 * Return results in a SET so we can easily do a FIND for a specific value
 * totalinodes is the maxmimum possible list that the local node must loop over for the given local data
 */

    void PsociGAhamiltonian::fetchLocalIcicolGArows( vector<int> & listinodes, set<int> & totalinodes, vector<set<int> > & icimat )
{
// Loop over number then fetch values. If number is <0 it reports the row in supp matrices
// to find relevant data

   int clow = local_cilo[0];
   int chi  = local_cihi[0];

   //int numrows = chi-clow+1;
   //icimat.resize(numrows); 

   totalinodes.clear();

   int biglo[2], bighi[2];
   int cilo[2], cihi[2];
   int lo[2], hi[2];
   int stride = 1;

   int n=0;
   int number=0;

   int index=0;
   int resize=1;

   int g_rank = GA::nodeid();
   int g_size = GA::nodes();

   cout << "junk " << g_rank<<" "<<g_size<<endl;
   //for (int i=clow; i<=chi; ++i) {
   
//   for(int i=g_rank; i<l_maxsef; i=i+g_size) {
//   How to bias the distribution of rows without using nxtval which performs badly at scale

    int bunch = 1; // used to construct a bunching parameter

//    for(int i=l_maxsef-1; i>=0; i=i-g_size) {
      for(int i=g_rank; i< l_maxsef; i+=g_size ) {

      icimat.resize(resize++);
      listinodes.push_back(i);
      cout << g_rank <<" Inserting " << i << endl;
      vector<int> temp(max(l_maxsparse,l_maxsparse_supp),0);

      //icimat[index].resize(max(l_maxsparse,l_maxsparse_supp));
      //cout << GA::nodeid() <<" size " << icimat[index].size() << endl;

      cilo[0] = i;
      cihi[0] = i; // This may not apply to the supp space: need a check at creation
      cilo[1] = 0;
      cihi[1] = 0;
      n=1;
      l_g_number->get( cilo, cihi, &number, &n );
      if ( number >= 0 ) {
         cihi[1] = number - 1;
         n = number;
         l_g_icicol->get( cilo, cihi, &temp[0], &n);
         temp.resize(number);
      } else {
         int suppRow = -number - 1;
         cilo[0] = suppRow;
         cihi[0] = suppRow; 
         int n=1;
         l_g_number_supp->get( cilo, cihi, &number, &n );
         cihi[1] = number - 1;
         n = number;
         l_g_icicol_supp->get( cilo, cihi, &temp[0], &n);
         temp.resize(number);
       }

       //cout << "size of vector TO set is " << temp.size() << endl;
       set<int> newset(temp.begin(), temp.end());
       //cout << "size of set is " << newset.size()<< endl;
       icimat[index]=newset;

       set<int>::iterator iter;
       for (iter = newset.begin(); iter != newset.end(); ++iter) {
           totalinodes.insert(*iter);
       //    cout << "iter " << *iter << endl;
       }
/*
      cout << GA::nodeid() << " new row, " << " i is " << i << " n is " << n << endl;
      for(int r=0; r<n; ++r ) {
        cout << icimat[index][r] << " ";
      }
      cout << endl;
*/
      index++;
   }
   cout << " size of icimat " << icimat.size() << endl;
   return;
}


/*
 * To not be used
 */
    void PsociGAhamiltonian::fetchGArowFulllength(vector<int> & temp, int irow )
{
// Loop over number then fetch values. If number is <0 it reports the row in supp matrices
// to find relevant data

   int biglo[2], bighi[2];
   int cilo[2], cihi[2];
   int lo[2], hi[2];
   int stride = 1;

   const char zero='0';
   const char one='1';

   int width=max(l_maxsparse_supp, l_maxsparse);

//cout << "irow is: How big is width " << irow <<" "<<width << endl;
//cout << l_g_icicol_supp->ndim() << endl;

     int i = irow;
      vector<int> icicol( width, 0 ); 
      lo[0] = i;
      hi[0] = i;
      lo[1]= 0;
      hi[1]= 0;
      int n=-1;

      l_g_number->get( lo, hi, &n, &stride );

      if ( n >= 0 ) {
         cilo[0] = i;
         cihi[0] = i;
         cilo[1] = 0;
         cihi[1] = n-1;
         stride=n;
         l_g_icicol->get( cilo,cihi,&icicol[0],&stride);

      } else {
         int suppRow = (-1*n) - 1;
         cilo[0] = suppRow;
         cihi[0] = suppRow;
         cilo[1]=0;
         cihi[1]=0; 
         int nbig=-1;
         l_g_number_supp->get( cilo, cihi, &nbig, &stride );
         cilo[1] = 0;
         cihi[1] = nbig-1;
         stride = nbig;

//cout << GA::nodeid() << " fetch ICI " << suppRow<<" "<<nbig-1<< endl;
//cout << "TEST " << icicol[0]<<" "<<icicol[1]<<endl;

         l_g_icicol_supp->get( cilo,cihi,&icicol[0],&stride);
         n = nbig;
      }
         for (int ibyte=0; ibyte<n; ibyte++ ) { 
	     int icol = icicol.at(ibyte);
             //temp[icol]='1'; // sets value to !0 
             temp[icol]= icol; // sets value to !0 
         }

      int count=0;
      return;
}



