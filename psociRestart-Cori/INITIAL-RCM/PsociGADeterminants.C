//TODO remove the superfluous GA::sync() calls used for debugging
//TODO clean up the sillyness in the initialization code
 
/**************************************************************************************

* Copyright (c) 2010,2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: 
 
 Methods to manage the put/get of determinant data to/fro local memory to GA distributed memeory
 Mostly this means put and get operations and a data assembly disassembly step ( formally called pck and unpck )
             
 History:

 Added new 2D methods to fetch blocks of contious data from GA to decrease fetch latency. The calling app
 is responsible for a rational choice of subblock sizes.

 Removed default use of the 2D code becasue it was slower

 Added new methods to support ORBMAP processing for H and CIDENSITY objects

**************************************************************************************/
/**
 *   @file PsociGADeterminants.C
 *
 */

#include <cmath> 
#include <set>

#include "PsociTimer.hpp"
#include "PsociGADeterminants.hpp"

using namespace std;

//Preexisting GA handle required - sizes are allocated here.

/* This version creates a persistant PsociDetermanants object that may be used
   for the pseudo-direct CI methods. Remember to FREE this on the destructor

   ALTERNATIVE constructor required for computeConfig method for getting determinants.
*/
PsociGADeterminants::PsociGADeterminants( PsociConfigs * configs, GA::GlobalArray * global_det, GA::GlobalArray * global_det_num , GA::GlobalArray * global_nsef, GA::GlobalArray * global_orbmap )
{
  g_det = global_det;
  g_det_num = global_det_num;
  g_nsef = global_nsef;
  g_orbmap = global_orbmap;
  
  local_config = configs;

  l_maxSpatials = 0;
  
  //if ( GA::nodeid() == 0 ) cout << "Defauting to DIRECT and SmallMemoryConfigs model " << endl;
  local_det = new PsociDeterminants( local_config );

//   PsociDeterminants deters( &configs ); 
//  local_det = new PsociDeterminants( &configs );

// TODO add destroy on deters

}

/* Needed with the ALTERNATIVE constructor */
void PsociGADeterminants::destroyLocalDeterminantsData()
{
     cout <<" Detroying local determinant data " << endl;
     delete local_det;
}

/* Preexisting GA handle required - sizes are allocated here.

   STANDARD constructor for using fetch methods
*/
PsociGADeterminants::PsociGADeterminants( GA::GlobalArray * global_det, GA::GlobalArray * global_det_num , GA::GlobalArray * global_nsef, GA::GlobalArray * global_orbmap )
{
  g_det = global_det;
  g_det_num = global_det_num;
  g_nsef = global_nsef;
  g_orbmap = global_orbmap;
  
  l_maxSpatials = 0;
}

/* for use with ALTERNATIVE methods */

void PsociGADeterminants::destroyJoutfgGA()
{
  g_det->destroy();
  g_det_num->destroy();
  g_nsef->destroy();
  g_orbmap->destroy();
#ifdef REPLICATEDETVALUE
  replicated_phases.clear();
#endif
}

/* Intended for use with ALTERNATIVE methods */
PsociDeterminants * PsociGADeterminants::fetchInternalDeters()
{
   return( local_det );
}

vector<char> * PsociGADeterminants::fetchReplicatedPhasesHandle()
{
   return( &replicated_phases );
}


int PsociGADeterminants::fetchMaxSpatials()
{
  return( l_maxSpatials );
}

GA::GlobalArray * PsociGADeterminants::fetchDeterminantHandle()
{
  return( g_det );
}

GA::GlobalArray * PsociGADeterminants::fetchDeterminantNumberHandle()
{
  return( g_det_num );
}


// Put some checks into this
void PsociGADeterminants::printSpatialsData( JOUTFG & ioutfg )
{
  local_det->printSpatialsData(ioutfg );
  return;
}

int PsociGADeterminants::fetchNumElectrons()
{
  return( l_nelec ) ;
}

int PsociGADeterminants::fetchNumBasisFunctions()
{
  return( l_nbfn );
}

int PsociGADeterminants::fetchGlobalMaxSefPerSpatial()
{
  return(  l_maxsefperconf  );
}

int PsociGADeterminants::fetchGlobalMaxSefs()
{
  return(  l_maxsef );
}

int PsociGADeterminants::fetchGlobalMaxDetPerSpatial()
{
  return(  l_maxdetperconf );
}

int PsociGADeterminants::fetchGlobalSymmetry()
{
  return( l_ksym );
}


//Collective
//Because we changed the manner in which parameters are fetched we no longer can use this directly

//Timed wrappers
void PsociGADeterminants::createDistributedDeterminants(  PsociDeterminants * deters, pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  createDistributedDeterminants( deters ); 
  info.second = psociTime() - timein;
  info.first = g_rank;
}

void PsociGADeterminants::createDistributedDeterminants( PsociDeterminants * deters )  
{
  local_det = deters;
  
  local_det->assembleGlobalJobParameters(); //Ensure a consistent set of job parameters
  
  setParameters( deters ); // Everyone comes through here - set the relevent parameters
  
  // Keep these assignments for now
  
  int maxSpatials = l_maxSpatials;
  int mxseti = l_maxsefperconf;
  int mxdeti = l_maxdetperconf;
  int nelec = l_nelec;
  int nbf = l_nbfn;
  int mxopn = l_mxopn;
  
  int mxspatials = maxSpatials;
  
  if ( GA::nodeid() == 0 ) {
    cout << "Generating new GA for determinants using: mxopn = " << mxopn;
    cout << " mxdeti = " << mxdeti << " nelec = " << nelec << " nbf = " << nbf << endl;
    cout << " mxseti = " << mxseti << " nelec = " << nelec << " nbf = " << nbf << endl;
    cout << "Maximum number of spatials is " << mxspatials << endl;
  }
  
  int mxocc = ( nelec + mxopn ) /2;
  int mxclosed = nelec / 2; 
  
  l_maxlength = computeLength( mxdeti, mxocc, mxclosed, mxopn ); //maximum possible buffer length in ints
  
  if ( GA::nodeid() == 0 ) {
    cout << "The estimated amount of total GA memory (B) to allocate for Dets is " << sizeof(COMPRESS_SIZE) * maxSpatials * l_maxlength << endl; 
    cout << "The estimated amount of per-core GA memory (MB) to allocate for Dets is " << ( sizeof(COMPRESS_SIZE) * maxSpatials * l_maxlength ) / (1000*1000*GA::nodes()) << endl;
  }
  
#if 0
  cout << "New flipped orientation code " << endl;
#endif

  int DATA_TYPE = C_INT; 
  int DATA_NDIM = 2;
  int dims[2];
  dims[1] = l_maxlength;
  dims[0] = mxspatials;

  const int DISTRIBUTE_EVENLY = -1; //TODO need to optimize this at some point
  int chunk[2];

  
#ifndef DIRECT
  
  chunk[1] = dims[1]; //We want distribution by rows (length)
  chunk[0] = DISTRIBUTE_EVENLY;
  
  char local_title[] = "determinants per conf";
  g_det = GA::createGA( DATA_TYPE, DATA_NDIM, dims, (char *) local_title, chunk); 
  
  char handleGAMessage[] = "PsociGADeterminants::createDistributedDeterminants";
  g_det->checkHandle( handleGAMessage); 
  g_det->zero();
  
  dims[1] = 1;
  dims[0] = maxSpatials;
  chunk[1] = DISTRIBUTE_EVENLY;
  chunk[0] = DISTRIBUTE_EVENLY;
  
  char local_title2[] = "g_num_det ";
  g_det_num = GA::createGA( DATA_TYPE, DATA_NDIM, dims, (char *) local_title2, chunk);
  
  char handleGAMessage2[] = "PsociGADeterminants::createDistributedDeterminants:g_number";
  g_det_num->checkHandle( handleGAMessage2 );
  g_det_num->zero();
#endif
  
  dims[1] = 1;
  dims[0] = maxSpatials;
  chunk[1] = DISTRIBUTE_EVENLY;
  chunk[0] = DISTRIBUTE_EVENLY;

  char local_title3[] = "g_nsef ";
  g_nsef = GA::createGA( DATA_TYPE, DATA_NDIM, dims, (char *) local_title3, chunk);
  
  char handleGAMessage3[] = "PsociGADeterminants::createDistributedDeterminants:g_nsef";
  g_nsef->checkHandle( handleGAMessage3 );
  g_nsef->zero();

  DATA_TYPE = C_INT;
  DATA_NDIM = 2; 
  dims[1] = nbf+2; // append jsym and iconf to the end of list
  dims[0] = mxspatials;
  if ( GA::nodeid() == 0 ) cout << "Create a GA-ORBMAP of dims " << dims[0] << " by " << dims[1] << endl;
  chunk[1] = dims[1]; //We want distribution by rows (length)
  chunk[0] = DISTRIBUTE_EVENLY;
    
  char local_orb_title[] = "orbital maps + jsym and iconf";
  g_orbmap = GA::createGA( DATA_TYPE, DATA_NDIM, dims, (char *) local_orb_title, chunk);

  char handleGAMessage4[] = "PsociGADeterminants::createDistributedDeterminants:g_orbmap";
  g_orbmap->checkHandle( handleGAMessage4 );
  g_orbmap->zero();

  int ndet_loaded = packAndLoadDeterminantData();
  
#ifndef DIRECT
  assembleDimensions(); // Everyone should have this
#endif

#ifdef DIRECT
  cout << "DIRECT: skip GA space for Determinants " << endl;
#endif
  
  /* Now allocate relatively persistant buffers for use by download method
*/
  
  l_doubly.resize( mxocc );
  l_singly.resize( mxopn );
  l_spin_projection.resize( mxdeti );
  l_phases.resize( mxdeti );
  l_sumsa.resize( mxdeti );
  l_sumsb.resize( mxdeti );
  
  alpha_spins.resize(  mxdeti  );
  beta_spins.resize(  mxdeti  );
  aspins.resize( mxocc );
  bspins.resize( mxocc );
}

int PsociGADeterminants::fetchMaxLength()
{
  return( l_maxlength );
}

void PsociGADeterminants::printDistribution()
{
  int g_rank = GA::nodeid();
#ifndef DIRECT
  if ( g_rank == 0 ) cout << " Psoci: Distributed Determinant array " << endl;
  g_det->printDistribution();
  if ( g_rank == 0 ) cout << " Psoci: Distributed Determinant number array " << endl;
  g_det_num->printDistribution();
#endif
  if ( g_rank == 0 ) cout << " Psoci: Distributed nsef array " << endl;
  g_nsef->printDistribution();
  if ( g_rank == 0 ) cout << " Psoci: Distributed orbmap array " << endl;
  g_orbmap->printDistribution();
}

void PsociGADeterminants::fetchNumDimensions( int & rows, int & cols )
{
  int DATA_NDIM = 2;
  int type, ndim;
  
  int dims[ DATA_NDIM ];
  g_det->inquire( &type, &ndim, dims );  
  
  if ( ndim == 2 ) {
    rows = dims[0];
    cols = dims[1];
  } else {
    cerr << "fetchDimensions: Erroneous dims: " << ndim << " Aborting " << endl;
    GA::Terminate();
  }
}

void PsociGADeterminants::assembleDimensions()
{
  int DATA_NDIM = 2;
  int lo[ DATA_NDIM ];
  int hi[ DATA_NDIM ];
  
  g_det->distribution( GA::nodeid(), lo, hi );
  
  l_lrow=lo[0];
  l_hrow=hi[0];
  
  l_lcol=lo[1];
  l_hcol=hi[1];
}

void PsociGADeterminants::printDimensions()
{
  cout << "Local determinant coordinates for ID = " << GA::nodeid() << endl;
  cout << "l_lrow is " << l_lrow << " l_hrow is " << l_hrow << endl;
  cout << "l_lcol is " << l_lcol << " l_hcol is " << l_hcol << endl;
}

/* assemble, pack, and dump my local set of spatials to the GA-DET array
   Need to perform a contiguous packing into an array of COMPRESS_SIZE elements
*/

// GA space is always INTS while buffer may be compressed_ints

/* A constraint exists here. You MUST call this with member arrays that are of the correct length.
   This is not a problem since one never calls this after a has begun. Note, the associated FETCH method
   called fetchAndUnpackPreallocate() is opposite. That is sizes are NOT checked and must be preallocated
*/

int PsociGADeterminants::packAndLoadDeterminantData() 
{
  int num_configs_loaded=0;
  int ksym = fetchGlobalSymmetry();
  
  int nelec = l_nelec;
  int mxndeti = l_maxdetperconf;
  int mxopn = l_mxopn;
  
  int mxocc = ( nelec + mxopn ) /2;
  int mxclosed = nelec / 2;
  
  int maxsize = computeLength( mxndeti, mxocc, mxclosed, mxopn );
  

  vector<JOUTFG> * joutfg = local_det->fetchJoutfg();  // vector of locally stored spatials+determinants
  vector<JOUTFG>::iterator jit;

  for(jit = joutfg->begin(); jit != joutfg->end(); ++jit ) { // Loop over local spatials
    
    int ndeti = (*jit).ndeti;
    
/* The 13_1.1.117 compiler doens correctly read the bookl used for this test. 
   It also doesn't it. So for now skip it.
*/
#ifndef STAMPEDE
    if ( (*jit).completed_state == 0 ) {
      GA::error("Determinents data are incomplete at index = ",(*jit).index);
    }
#endif
    // Load the parameters into buffer
    
    int nopn = (*jit).num_singly;
    int nocc = ( nelec + nopn ) /2;
    int nclosed = nocc - nopn;
    int tstsize = computeLength( ndeti, nocc, nclosed, nopn ); 
    
//cout << "PACK " << nocc<<" "<<ndeti<<" "<<nclosed<<" "<<nopn<<endl;
//cout << "tst mxsize is " << tstsize<<" "<<maxsize<<endl;

    int index = (*jit).index; //For picking GA space
    
    //cout << "tstsize and mxsize are " << tstsize << " " << maxsize << endl;
    if ( tstsize > maxsize ) {
      cerr << "tstsize > maxsize for spatial = " << (*jit).index << "Aborting " << endl;
      GA::Terminate();
    }
    
/* At the time of this change, the Stampede icc composer_xe_2013.1.117 compiler was generating
   bad code when using the vector stl. THis was confirmed with allan porterfield. The case was
   too difficult to generate a tight little test case. So for now assume stampede is c99 compliant
*/
    vector<COMPRESS_SIZE> orbmap( (*jit).nbf+2, 0 ); // New code to generate orbital maps for screening
    vector<COMPRESS_SIZE> l_buffer(tstsize); // don;t use class member: this guarantees contiguous storage

    l_buffer[0] = (*jit).mx_openshells;
    l_buffer[1] = (*jit).nbf;
    l_buffer[2] = (*jit).nelec;
    l_buffer[3] = (*jit).ksym;
    l_buffer[4] = (*jit).jsym;
    l_buffer[5] = (*jit).nsefi;
    l_buffer[6] = nopn;
    l_buffer[7] = (*jit).num_doubly;
    l_buffer[8] = (*jit).ndeti;
    l_buffer[9] = index;
    int i = 10;
    
    // :doubly
    vector<int>::iterator vit;
    for(vit=(*jit).doubly.begin(); vit != (*jit).doubly.end(); ++vit) {
      l_buffer[i] = (*vit);
      orbmap[ (*vit) ] = 2; // this is redundant to the the check donw in SOCIgenerate
      ++i; 
    }
    // :singly
    for(vit=(*jit).singly.begin(); vit != (*jit).singly.end(); ++vit) {
      l_buffer[i] = (*vit);
      orbmap[ (*vit) ] = 1;
      ++i;
    }
    // :spin_projection ( formally called ms )
    for(vit=(*jit).spin_projection.begin(); vit != (*jit).spin_projection.end(); ++vit) {
      l_buffer[i] = (*vit);
      ++i;
    }
    // :phase 
    for(vit=(*jit).phases.begin(); vit != (*jit).phases.end(); ++vit) {
      l_buffer[i] = (*vit);
      ++i;
    }
    // :sumsa 
    for(vit=(*jit).sumsa.begin(); vit != (*jit).sumsa.end(); ++vit) {
      l_buffer[i] = (*vit);
      ++i;
    }
    // :sumsb
    for(vit=(*jit).sumsb.begin(); vit != (*jit).sumsb.end(); ++vit) {
      l_buffer[i] = (*vit);
      ++i;
    }
    // :alphaorbs
    vector<vector<int> >::iterator gavec;
    for(gavec = (*jit).alpha_spins.begin(); gavec != (*jit).alpha_spins.end(); ++gavec ) {
      l_buffer[i] =(*gavec).size(); // We could use the new nalpha array but this should be fine
      ++i;
      vector<int>::iterator it;
      for(it=(*gavec).begin(); it != (*gavec).end(); ++it) {
	l_buffer[i] = (*it);
	++i;
      }  
    }
    // :betaorbs
    for(gavec = (*jit).beta_spins.begin(); gavec != (*jit).beta_spins.end(); ++gavec ) {
      l_buffer[i] =(*gavec).size(); // We could use the new nbeta array but this should be fine
      ++i;
      vector<int>::iterator it;
      for(it=(*gavec).begin(); it != (*gavec).end(); ++it) {
	l_buffer[i] = (*it);
	++i;
      }  
    }
    
    orbmap[ (*jit).nbf ] =  (*jit).jsym;
    orbmap[ (*jit).nbf + 1 ] =  index; // training info helpful for the DIRECT version of the code

    int obslength = i; // in int words
    
    if ( obslength > tstsize ) {
      cout << "Warning: inappropriate size obslength != tstlength. tstlngth " << tstsize << " obs " << obslength << endl;
      cout << "This should never happen and probably indicates an error: Aborting " << endl;
      GA::Terminate();
    }
    
#ifndef DIRECT
    if ( uploadBufferToGA( index, tstsize, l_buffer ) != 0 ) {
      cerr << "Failed to upload configuration number " << index << endl;
      GA::Terminate();
    }
#endif
    
    uploadNsefToGA( index, (*jit).nsefi ); // Needed for hamiltonian construction steps

/* Now push the orbmap to GA for index ). This will be used later for quick and dirty 
   screening of the iconf/jconf loops
*/
   if ( uploadOrbMaptoGA( index,  orbmap ) != 0 ) {
      GA::error(" Failed to upload orbmap at", index);
   }

    ++num_configs_loaded;
  } //configs
  
  return( num_configs_loaded );
}

/* Warning: ints versus compressed_ints may be a problem
   Unfortunately, for the time being GA has no C_SHORT data type
*/

int PsociGADeterminants::uploadBufferToGA( int index, int tstsize, vector<COMPRESS_SIZE> & buffer) 
{
  char handleMessage[] = "PsociGADeterminants::uploadBufferToGA";
  g_det->checkHandle( handleMessage);
  
  int lo[2];
  lo[1] = 0;
  lo[0] = index-1; // 
  int hi[2];
  hi[1] = tstsize-1;
  hi[0] = index-1;
  int n;
  
/*
    cout << " UPLOAD " << tstsize << " " << index << endl;
    for(int i=0; i< tstsize; ++i ) {
    cout << buffer[i] << " " ;
    }
    cout << endl;
*/

  n = tstsize;
  g_det->put( lo, hi, &buffer[0], &n);
  
  hi[1] = 0;
  n = 1;
  g_det_num->put( lo, hi, &tstsize, &n);
  
  return(0);
}

/* The application must determine the proper length of the column to fetch. Note: we do not want to fetch more 
   than required as this stresses the network. Will probably need another data object that maps index to length.
*/
// GA itself will trap in appropriate indices ...

int PsociGADeterminants::downloadBufferFromGA( int index, int maxsize, int & obs_size, vector<COMPRESS_SIZE> & buffer)
{

#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGADeterminants::downBufferFromGA";
  g_det->checkHandle( handleMessage);
  
  char * handleMessage2 = "PsociGADeterminants::downBufferFromGA: number array ";
  g_det_num->checkHandle( handleMessage2);
#endif
  
  int lo[2];
  int hi[2];
  int n;
  lo[1] = 0;
  lo[0] = index-1;  
  hi[1] = 0;
  hi[0] = index-1;
  int tstsize;
  n = 1;
  g_det_num->get( lo, hi, &tstsize, &n);
  
#ifdef DETAILEDCHECK
  if ( tstsize == 0 || tstsize > maxsize ) {
    cerr << "At fetch tstsize > maxsize or it is  0 " << tstsize << " " << maxsize << endl;
    GA::Terminate();
  }
#endif
  
  if ( buffer.size() != tstsize ) buffer.resize( tstsize );
  
  lo[1] = 0;
  lo[0] = index-1; // 
  hi[1] = tstsize-1;
  hi[0] = index-1;
  n = tstsize;
  
  g_det->get( lo, hi, &buffer[0], &n);
  
  obs_size = tstsize;
  return(0);
}

/* A newer attempt to speed things up. In this method we PRESUME the buffer is preallocated and thus do no checks 
   nor resizing. It speeds things up abit. You need to compile test codes with -DPREALLOCATE
*/ 
int PsociGADeterminants::downloadBufferFromGAPreallocate( int index, int maxsize, int & obs_size, vector<COMPRESS_SIZE> & buffer)
{

#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGADeterminants::downBufferFromGA";
  g_det->checkHandle( handleMessage);
  
  char * handleMessage2 = "PsociGADeterminants::downBufferFromGA: number array ";
  g_det_num->checkHandle( handleMessage2);
#endif
  
  int lo[2];
  int hi[2];
  int n;
  lo[1] = 0;
  lo[0] = index-1;  
  hi[1] = 0;
  hi[0] = index-1;
  int tstsize;
  n = 1;
  g_det_num->get( lo, hi, &tstsize, &n);
  
#ifdef DETAILEDCHECK
  if ( tstsize == 0 || tstsize > maxsize ) {
    cerr << "At fetch tstsize > maxsize or it is  0 " << tstsize << " " << maxsize << endl;
    GA::Terminate();
  }
#endif
  
  lo[1] = 0;
  lo[0] = index-1; // 
  hi[1] = tstsize-1;
  hi[0] = index-1;
  n = tstsize;
  
  g_det->get( lo, hi, &buffer[0], &n);
  
  obs_size = tstsize;
  return(0);
}

/* This fetch method downloads a 2D block of determinant data.It is marginally helpful. To download a block, GA requires
   a rectangular array. As a result ( for examplek if one really large determinant is in that block) A large ( heavily padded) 
   block gets read. In real tests, the amounts of data read can increase significaNTLY.
*/ 
int PsociGADeterminants::downloadBufferFromGA2D( int index1, int index2, int maxsize, int & obs_size, vector<COMPRESS_SIZE> & buffer)
{

#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGADeterminants::downBufferFromGA2D";
  g_det->checkHandle( handleMessage);
  
  char * handleMessage2 = "PsociGADeterminants::downBufferFromGA2D: number array ";
  g_det_num->checkHandle( handleMessage2);
#endif

  int tstsize=0;

  int lo[2];
  int hi[2];
  int n;
  lo[1] = 0;
  hi[1] = 0;

/* This seems to run out of allocation on Hopper - MPIU UNSAFE ALLOC EFFORT
   is the error thrown by MPI2. The reason is as yet undetermined. I have a suspicion
   that this loop may be causing too many small messages to be in-flight.

  int width = index2 - index1  + 1; // We grab index1-index2 inclusive
  for(int i=index1; i<=index2; ++i)
  {
  index = i;
  lo[0] = index-1;
  hi[0] = index-1;
  n = 1;
  g_det_num->get( lo, hi, &tempsz, &n);
  tstsize = max( tstsize, tempsz ); // Who is the longest member of this block ?
  cout << "Original; version " << tempsz<<" "<<tstsize<<endl;
  }
*/

// Grab a bunch of values at one time
  int width = index2 - index1  + 1; // We grab index1-index2 inclusive
  vector<int> templist( width, 0);
  lo[0] = index1-1;
  hi[0] = index2-1;
  n = 1;
  g_det_num->get( lo, hi, &templist[0], &n);

  for(int i=0; i<width; ++i)
  {
     tstsize = max( tstsize, templist[i] ); // Who is the longest member of this block ?
  }

#ifdef DETAILEDCHECK
  if ( tstsize == 0 || tstsize > maxsize ) {
    cerr << "At fetch tstsize > maxsize or it is  0 " << tstsize << " " << maxsize << endl;
    GA::Terminate();
  }
#endif
  
  if ( buffer.size() != tstsize * width ) buffer.resize( tstsize * width );
  
  lo[1] = 0;
  lo[0] = index1-1; // 
  hi[1] = tstsize-1; //
  hi[0] = index2-1;
  
  n = tstsize;
  
  g_det->get( lo, hi, &buffer[0], &n);

//cout << GA::nodeid() << " downloadBufferFromGA2D index is " << index1 << " "<<index2<<" value was " << buffer[9] << endl;

  obs_size = tstsize;
  
  return(0);
}

/* The PREALLOCATE analogue of the GA2D block fetch. IN this case all arrays must be PREALLOCATED. It is a little faster.
 */
int PsociGADeterminants::downloadBufferFromGA2DPreallocate( int index1, int index2, int maxsize, int & obs_size, vector<COMPRESS_SIZE> & buffer)
{

#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGADeterminants::downBufferFromGA2DPreallocate";
  g_det->checkHandle( handleMessage);
  
  char * handleMessage2 = "PsociGADeterminants::downBufferFromGA2DPreallocate: number array ";
  g_det_num->checkHandle( handleMessage2);
#endif

  int tstsize=0, tempsz;
  int index;

  int lo[2];
  int hi[2];
  int n;
  lo[1] = 0;
  hi[1] = 0;

//TODO Replace thie with a GA MAX call

  for(int i=index1; i<=index2; ++i)
  {
  index = i;
  lo[0] = index-1;
  hi[0] = index-1;
  n = 1;
  g_det_num->get( lo, hi, &tempsz, &n);
  tstsize = max( tstsize, tempsz );
  }

#ifdef DETAILEDCHECK
  if ( tstsize == 0 || tstsize > maxsize ) {
    cerr << "At fetch tstsize > maxsize or it is  0 " << tstsize << " " << maxsize << endl;
    GA::Terminate();
  }
#endif
  
  
  lo[1] = 0;
  lo[0] = index1-1; // 
  hi[1] = tstsize-1; //
  hi[0] = index2-1;

  n = tstsize;
  
  g_det->get( lo, hi, &buffer[0], &n);
  obs_size = tstsize;

//cout << "downloadBufferFromGA2DPreallocate " <<index<<" value was " << buffer[9]<<endl;
  
  return(0);
}

//Timed wrappers
int PsociGADeterminants::fetchAndUnpackDeterminantData( int list, pair<int,double> & info, JOUTFG & joutfg  )
{
  int status;
  int g_rank = GA::nodeid();
  double timein = psociTime();
  status = fetchAndUnpackDeterminantData(list, joutfg );
  info.second = psociTime() - timein;
  info.first = g_rank;
  return( status );
}

int PsociGADeterminants::fetchAndUnpackDeterminantDataPreallocate( int list, pair<int,double> & info, vector<COMPRESS_SIZE> & buffer, JOUTFG & joutfg  )
{
  int status;
  int g_rank = GA::nodeid();
  double timein = psociTime();
  status = fetchAndUnpackDeterminantDataPreallocate(list, buffer, joutfg );
  info.second = psociTime() - timein;
  info.first = g_rank;
  return( status );
}

/* Unpack a fetched column form GA into a local joutfg object.
   We do not expect to need the JOUTFG objects anymore. So simply
   unpack into a generic JOUFG object
*/
//TODO correct useless return value

int PsociGADeterminants::fetchAndUnpackDeterminantData( int index, JOUTFG & joutfg  ) 
{
  int mxlength = l_maxlength; // inefficient for now
  int obslength;
  
  //  vector<COMPRESS_SIZE> buffer; // buffer gets properly sized in downloadBufferFromGA
  //  Ensure the size of the read data and the parmaters are consistent;
  
  if ( downloadBufferFromGA( index, mxlength, obslength, buffer ) != 0 ) { 
    cerr << "Failed to download configuration number " << index << endl;
    GA::Terminate();
  }  
  
  //  JOUTFG ioutfg; 
  
  int mxopn = buffer[0];
  int nbf   = buffer[1];
  int nelec = buffer[2];
  int ksym  = buffer[3];
  int jsym  =  buffer[4];
  int nseti = buffer[5]; 
  int nopn = buffer[6];
  int num_doubly = buffer[7];
  int ndeti = buffer[8];
  int config_index = buffer[9]; //absolute index for this configuration
  int i = 10;
  
  //cout << nelec << " " << mxopn << " " << nbf << " " << jsym << " " << num_doubly << endl;
  //int nocc = ( nelec + nopn ) /2;
  //int nclosed = nelec / 2;
  
  // Read and format data
  // :l_doubly;

  if ( joutfg.doubly.size() != num_doubly ) joutfg.doubly.resize( num_doubly );
  for(int ii=0; ii < num_doubly; ++ii ) {
    joutfg.doubly[ii] = buffer[i+ii];
  }
  i += num_doubly;

  // :singly
  if ( joutfg.singly.size() != nopn ) joutfg.singly.resize( nopn );
  for(int ii=0; ii < nopn; ++ii ) {
    joutfg.singly[ii] = buffer[i+ii];
  }
  i += nopn;
  
  if ( joutfg.spin_projection.size() != ndeti ) {
    joutfg.spin_projection.resize( ndeti );
    joutfg.phases.resize( ndeti );
    joutfg.sumsa.resize( ndeti );
    joutfg.sumsb.resize( ndeti );
  }
  
  for(int ii=0; ii < ndeti; ++ii ) {
    joutfg.spin_projection[ii] = buffer[i+ii];
    joutfg.phases[ii] = buffer[ndeti + i+ii];
    joutfg.sumsa[ii] = buffer[2*ndeti + i+ii];
    joutfg.sumsb[ii] = buffer[3*ndeti +i+ii];
  }
  i += 4 * ndeti;

/* These are the original versions which are suspiciously slow
*/

/*
// :doubly
vector<int> doubly( &buffer[i], &buffer[i+num_doubly] );
i += num_doubly;

// :singly
vector<int> singly( &buffer[i], &buffer[i+nopn] );
i += nopn;

// :spin_projection 
vector<int> spin_projection( &buffer[i], &buffer[i+ndeti] );
i += ndeti;

// :phase
vector<int> phases( &buffer[i], &buffer[i+ndeti] );
i += ndeti;

// :sumsa
vector<int> sumsa( &buffer[i], &buffer[i+ndeti] );
i += ndeti;

// :sumsb
vector<int> sumsb( &buffer[i], &buffer[i+ndeti] );
  i += ndeti;
*/

/*
// These are the original implementations BUT they are consistently 
// tagged as bottlenecks by HpcToolKit. I don't yet kow why

// :alphaorbs
vector<vector<int> > alpha_spins;
for(int spins=0; spins < ndeti; ++spins ) {
int nalpha = buffer[i] ;
++i;
vector<int> aspins( &buffer[i], &buffer[i+nalpha] );
i += nalpha;
alpha_spins.push_back( aspins );
}

// :betaorbs
vector<vector<int> > beta_spins;
for(int spins=0; spins < ndeti; ++spins ) {
int nbeta = buffer[i] ;
++i;
vector<int> bspins( &buffer[i], &buffer[i+nbeta] );
i += nbeta;
beta_spins.push_back( bspins );
}
*/
  
  
  if ( joutfg.alpha_spins.size() != ndeti ) joutfg.alpha_spins.resize( ndeti );
  if ( joutfg.beta_spins.size() != ndeti ) joutfg.beta_spins.resize( ndeti );
  
  //Continue to process separately
  
  for(int spins=0; spins < ndeti; ++spins ) {
    int nalpha = buffer[i] ;
    if ( joutfg.nalpha.size () != nalpha ) joutfg.nalpha.resize( ndeti );
    if ( aspins.size() != nalpha ) aspins.resize( nalpha );
    ++i;
    for(int ii=0; ii< nalpha; ++ii) {
      aspins[ii] = buffer[i+ii];
    }
    i += nalpha;
    joutfg.nalpha[spins] = nalpha;
    joutfg.alpha_spins[spins] = aspins;
  }
  
  for(int spins=0; spins < ndeti; ++spins ) {
    int nbeta = buffer[i] ;
    if ( joutfg.nbeta.size () != nbeta ) joutfg.nbeta.resize( ndeti );
    if ( bspins.size() != nbeta ) bspins.resize( nbeta );
    ++i;
    for(int ii=0; ii< nbeta; ++ii) {
      bspins[ii] = buffer[i+ii];
    }
    i += nbeta;
    joutfg.nbeta[spins] = nbeta;
    joutfg.beta_spins[spins] = bspins;
  }
  
  //unpackTime2 += psociTime() - temp;
  // Add elements to local joutfg object and push into the stack
  // Not all elements need to be populated
  
  joutfg.index = config_index; //Want the absolute index
  joutfg.jsym = jsym;
  joutfg.ksym = ksym;
  joutfg.mx_openshells = mxopn;
  joutfg.num_doubly = num_doubly;
  joutfg.nelec = nelec;
  joutfg.nbf = nbf;
  joutfg.num_singly = nopn;
  joutfg.nsefi = nseti;
  joutfg.ndeti = ndeti;
  joutfg.completed_state=1; // Keep as true since we have the phases,sums,projections,etc.
  
  return( 0 );
}

/* Use a name of toutfg to indicate a temporary-max size joutfg
   These fetch methods are the PREALLOCATE API to the PREALLOCATE download methods
*/

int PsociGADeterminants::fetchAndUnpackDeterminantDataPreallocate( int index, vector<COMPRESS_SIZE> & buffer, JOUTFG & toutfg  ) 
{
  int mxlength = l_maxlength; // inefficient for now
  int obslength;
  
  int ndld = downloadBufferFromGAPreallocate( index, mxlength, obslength, buffer );
  
// Regardless of the size of member arrays downloaded from GA, stuff into this fixed size joutfg

  int mxopn = buffer[0];
  int nbf   = buffer[1];
  int nelec = buffer[2];
  int ksym  = buffer[3];
  int jsym  =  buffer[4];
  int nseti = buffer[5]; 
  int nopn = buffer[6];
  int num_doubly = buffer[7];
  int ndeti = buffer[8];
  int config_index = buffer[9]; //absolute index for this configuration
  int i = 10;
  
  //cout << nelec << " " << mxopn << " " << nbf << " " << jsym << " " << num_doubly << endl;
  //int nocc = ( nelec + nopn ) /2;
  //int nclosed = nelec / 2;
  
  // Read and format data
  // :l_doubly;
  
  for(int ii=0; ii < num_doubly; ++ii ) {
    toutfg.doubly[ii] = buffer[i+ii];
  }
  i += num_doubly;
  
  // :singly
  for(int ii=0; ii < nopn; ++ii ) {
    toutfg.singly[ii] = buffer[i+ii];
  }
  i += nopn;
  
  for(int ii=0; ii < ndeti; ++ii ) {
    toutfg.spin_projection[ii] = buffer[i+ii];
    toutfg.phases[ii] = buffer[ndeti + i+ii];
    toutfg.sumsa[ii] = buffer[2*ndeti + i+ii];
    toutfg.sumsb[ii] = buffer[3*ndeti +i+ii];
  }
  i += 4 * ndeti;
  

//Continue to process separately

  for(int spins = 0; spins < ndeti; ++spins ) {
    int nalpha = buffer[i] ;  
    toutfg.nalpha[spins] = nalpha;
    ++i;
    for(int ii=0; ii< nalpha; ++ii) {
    toutfg.alpha_spins[spins][ii] = buffer[i+ii];
    }
    i += nalpha;
  }

  for(int spins = 0; spins < ndeti; ++spins ) {
    int nbeta = buffer[i] ;
    toutfg.nbeta[spins] = nbeta;
    ++i;
    for(int ii=0; ii< nbeta; ++ii) {
      toutfg.beta_spins[spins][ii] = buffer[i+ii];
    }
    i += nbeta;
  }
  
  toutfg.index = config_index; //Want the absolute index
  toutfg.jsym = jsym;
  toutfg.ksym = ksym;
  toutfg.mx_openshells = mxopn;
  toutfg.num_doubly = num_doubly;
  toutfg.nelec = nelec;
  toutfg.nbf = nbf;
  toutfg.num_singly = nopn;
  toutfg.nsefi = nseti;
  toutfg.ndeti = ndeti;
  toutfg.completed_state=1; // Keep as true since we have the phases,sums,projections,etc.
  return( 0 );
}

// Check into status of BUFFER
// 't' nomenclature continues to indicate preallocation is performed

int PsociGADeterminants::fetchAndUnpackDeterminantData2DPreallocate( int index1, int index2, vector<COMPRESS_SIZE> & buffer, vector<JOUTFG> & toutfg  ) 
{
  int mxlength = l_maxlength; 
  int obslength;
  
  //vector<COMPRESS_SIZE> buffer; // buffer gets properly sized in downloadBufferFromGA2D(max(tstsize) * width)

  int num = downloadBufferFromGA2DPreallocate( index1, index2, mxlength, obslength, buffer );
  
  int width = index2 - index1 + 1; // Mostly validated by the calling app
  
  for(int iconf=0; iconf < width; ++iconf ) {
    int shift = iconf * obslength;
    
    int mxopn = buffer[shift+0];
    int nbf   = buffer[shift+1];
    int nelec = buffer[shift+2];
    int ksym  = buffer[shift+3];
    int jsym  =  buffer[shift+4];
    int nseti = buffer[shift+5]; 
    int nopn = buffer[shift+6];
    int num_doubly = buffer[shift+7];
    int ndeti = buffer[shift+8];
    int config_index = buffer[shift+9]; //absolute index for this configuration
    
    int i = shift + 10;
    
    // Read and format data
    
    // :l_doubly;
  for(int ii=0; ii < num_doubly; ++ii ) {
    toutfg[iconf].doubly[ii] = buffer[i+ii];
  }   
  i += num_doubly;
  
  // :singly
  for(int ii=0; ii < nopn; ++ii ) {
    toutfg[iconf].singly[ii] = buffer[i+ii];  
  }
  i += nopn;
  
  for(int ii=0; ii < ndeti; ++ii ) {
    toutfg[iconf].spin_projection[ii] = buffer[i+ii];
    toutfg[iconf].phases[ii] = buffer[ndeti + i+ii];
    toutfg[iconf].sumsa[ii] = buffer[2*ndeti + i+ii];
    toutfg[iconf].sumsb[ii] = buffer[3*ndeti +i+ii];
  }
  i += 4 * ndeti;
  
  //unpackTime1 += psociTime() - temp;
  //temp = psociTime();
  
  
  //Continue to process separately
  
  for(int spins = 0; spins < ndeti; ++spins ) {
    int nalpha = buffer[i] ;  
    toutfg[iconf].nalpha[spins] = nalpha;
    ++i;
    for(int ii=0; ii< nalpha; ++ii) {
      toutfg[iconf].alpha_spins[spins][ii] = buffer[i+ii];
      
    }
    i += nalpha;
  } 
  
  for(int spins = 0; spins < ndeti; ++spins ) {
    int nbeta = buffer[i] ;
    toutfg[iconf].nbeta[spins] = nbeta;
    ++i;
    for(int ii=0; ii< nbeta; ++ii) {
      toutfg[iconf].beta_spins[spins][ii] = buffer[i+ii];
    }
    i += nbeta;
  }
  
  toutfg[iconf].index = config_index; //Want the absolute index
  toutfg[iconf].jsym = jsym;
  toutfg[iconf].ksym = ksym;
  toutfg[iconf].mx_openshells = mxopn;
  toutfg[iconf].num_doubly = num_doubly;
  toutfg[iconf].nelec = nelec;
  toutfg[iconf].nbf = nbf;
  toutfg[iconf].num_singly = nopn;
  toutfg[iconf].nsefi = nseti;
  toutfg[iconf].ndeti = ndeti;
  toutfg[iconf].completed_state=1; // Keep as true since we have the phases,sums,projections,etc.

  }
  return( 0 );
}


/* All joutfg members are properly resized as needed */

int PsociGADeterminants::fetchAndUnpackDeterminantData2D( int index1, int index2, vector<JOUTFG> & joutfg  ) 
{
  int mxlength = l_maxlength; 
  int obslength;
  
//vector<COMPRESS_SIZE> buffer; // buffer gets properly sized in downloadBufferFromGA2D(max(tstsize) * width)
  
  int num = downloadBufferFromGA2D( index1, index2, mxlength, obslength, buffer );
  int width = index2 - index1 + 1; // Mostly validated by the calling app

  for(int iconf=0; iconf < width; ++iconf ) {
    int shift = iconf * obslength;
    
    JOUTFG ioutfg; 
    
    int mxopn = buffer[shift+0];
    int nbf   = buffer[shift+1];
    int nelec = buffer[shift+2];
    int ksym  = buffer[shift+3];
    int jsym  =  buffer[shift+4];
    int nseti = buffer[shift+5]; 
    int nopn = buffer[shift+6];
    int num_doubly = buffer[shift+7];
    int ndeti = buffer[shift+8];
    int config_index = buffer[shift+9]; //absolute index for this configuration
    int i = shift + 10;
    
    // Read and format data
    
  // :l_doubly;
  if ( ioutfg.doubly.size() != num_doubly ) ioutfg.doubly.resize( num_doubly );
  for(int ii=0; ii < num_doubly; ++ii ) {
     ioutfg.doubly[ii] = buffer[i+ii];
  }
  i += num_doubly;

  // :singly
  if ( ioutfg.singly.size() != nopn ) ioutfg.singly.resize( nopn );
  for(int ii=0; ii < nopn; ++ii ) {
     ioutfg.singly[ii] = buffer[i+ii];
  }
  i += nopn;

  if ( ioutfg.spin_projection.size() != ndeti ) {
        ioutfg.spin_projection.resize( ndeti );
        ioutfg.phases.resize( ndeti );
        ioutfg.sumsa.resize( ndeti );
        ioutfg.sumsb.resize( ndeti );
  }

  for(int ii=0; ii < ndeti; ++ii ) {
     ioutfg.spin_projection[ii] = buffer[i+ii];
     ioutfg.phases[ii] = buffer[ndeti + i+ii];
     ioutfg.sumsa[ii] = buffer[2*ndeti + i+ii];
     ioutfg.sumsb[ii] = buffer[3*ndeti +i+ii];
  }
  i += 4 * ndeti;


  if ( ioutfg.alpha_spins.size() != ndeti ) ioutfg.alpha_spins.resize( ndeti );
  if ( ioutfg.beta_spins.size() != ndeti ) ioutfg.beta_spins.resize( ndeti );

  if ( ioutfg.nalpha.size() != ndeti ) ioutfg.nalpha.resize( ndeti );
  if ( ioutfg.nbeta.size() != ndeti ) ioutfg.nbeta.resize( ndeti );

//Continue to process separately

  for(int spins=0; spins < ndeti; ++spins ) {
    int nalpha = buffer[i] ;
    if ( aspins.size() != nalpha ) aspins.resize( nalpha );
    ioutfg.nalpha[spins] = nalpha;
    ++i;
    for(int ii=0; ii< nalpha; ++ii) {
       aspins[ii] = buffer[i+ii];
    }
    i += nalpha;
    ioutfg.alpha_spins[spins] = aspins;
  }
  
  for(int spins=0; spins < ndeti; ++spins ) {
    int nbeta = buffer[i] ;
    if ( bspins.size() != nbeta ) bspins.resize( nbeta );
    ioutfg.nbeta[spins] = nbeta;
    ++i;
    for(int ii=0; ii< nbeta; ++ii) {
      bspins[ii] = buffer[i+ii];
    }
    i += nbeta;
    ioutfg.beta_spins[spins] = bspins;
  }
  
  ioutfg.index = config_index; //Want the absolute index
  ioutfg.jsym = jsym;
  ioutfg.ksym = ksym;
  ioutfg.mx_openshells = mxopn;
  ioutfg.num_doubly = num_doubly;
  ioutfg.nelec = nelec;
  ioutfg.nbf = nbf;
  ioutfg.num_singly = nopn;   
  ioutfg.nsefi = nseti;
  ioutfg.ndeti = ndeti;
  ioutfg.completed_state=1; // Keep as true since we have the phases,sums,projections,etc.
  
  joutfg.push_back( ioutfg );
  }
  
  return( 0 );
}

/* encapsulate this calculation because it is important and subject to the specifics of 
   the buffer construction which will change over time.

   returns length in INTEGER words (elements in buffer)
*/
inline int PsociGADeterminants::computeLength( int deti, int occ, int closed, int opn )
{
  //  return( ( 2*deti*occ + 6*deti + closed + opn + 10 ) * sizeof( COMPRESS_SIZE ) );
  return( ( 2*deti*occ + 6*deti + closed + opn + 10 ) );
}

//Local call Not sure if "me" can be anyone ( see ga++ docs )
void PsociGADeterminants::localConfRange( int & ilo, int & ihi )
{
  int DATA_NDIM = 2;
  int lo[ DATA_NDIM ];
  int hi[ DATA_NDIM ];
  g_det->printDistribution(); 
  g_det->distribution( GA::nodeid(), lo, hi );
  
  ilo=lo[0];
  ihi=hi[0];
}

//Local call
void PsociGADeterminants::uploadNsefToGA( int index, int nsefi )  
{
  int DATA_NDIM = 2;
  int lo[ DATA_NDIM ];
  int hi[ DATA_NDIM ];
  
  char * handleMessage = "PsociGADeterminants::uploadNsefToGA";
  g_nsef->checkHandle( handleMessage);
  
  int indexG = index-1;
  lo[1] = 0;
  lo[0] = indexG; // 
  hi[1] = 0;
  hi[0] = indexG;
  int n=1;
  
  g_nsef->put( lo, hi, &nsefi, &n);

}

// TODO integrate the several assembly methods into a single method
// assembles and replicates the g_nsef array to all nodes.
/* I choose not to use "gather' for now because data may be shuffled, defeating the purpose
   of this array which maintains order in the Hamiltonian matrix 
*/
// This array was previously called g_nsef giving rise to isef, jsef

void PsociGADeterminants::assembleAndBrdcstNsef( vector<int> & local_nsef )
{
//This data structure is TOO BIG

  //  cout << " maxspatials is " << l_maxSpatials << endl;
  local_nsef.resize(l_maxSpatials,0);
  
  if ( GA::nodeid() == 0 ) {

/*
    vector<int> workList;
    for( int i=0; i< l_maxSpatials; ++i ) { 
       workList.push_back( i );
    }
   cout << "Perform scrambling step " << endl;
   random_shuffle(workList.begin(), workList.end());
   random_shuffle(workList.begin(), workList.end());
   for(int i=0; i< 10; ++i ) {
   cout << "List "<< workList[i] << endl;
   }
   cout << "Big change: assembleAndBrdcstNsef: build a scrambled l_nsef object for now just reverse it"<<endl;
*/

    int DATA_NDIM = 2;
    int lo[ DATA_NDIM ];
    int hi[ DATA_NDIM ];
    
    lo[1] = 0;
    hi[1] = 0;
    int n;
    int temp;
    
/* old version
*/
// for(int i=0; i< l_maxSpatials; ++i ) {
//       for(int i=l_maxSpatials-1; i>= 0; --i ) {

/*
    vector<int>::iterator it;
    for(it =workList.begin(); it!=workList.end(); ++it) {
      int i = (*it);
       lo[0] = i; //  
       hi[0] = i;
       n = 1;
       
       g_nsef->get( lo, hi, &temp, &n);
       local_nsef[i] = temp;
       }       
*/
    
    // New version
    lo[0] = 0; 
    hi[0] = l_maxSpatials - 1 ;
    n = 1;
    g_nsef->get( lo, hi, &local_nsef[0], &n);
  }

/* broken
   if ( GA::nodeid() == 0 ) {
   cout << "CHECK : Perform scrambling step on lnsef" << endl;
   random_shuffle(local_nsef.begin(), local_nsef.end());
   }
*/

  GA::brdcst( &local_nsef[0], sizeof(int), 0) ;
}

// Return local per-core memory 
// Temporarily only havine groot==0 do any processing 

long PsociGADeterminants::generateDistributedDeterminants(int crap)
{
  long localsize = 0;
  int g_rank = GA::nodeid();
  int g_size= GA::nodes();
  
//First step create distributed (local) determinants

  PsociDeterminants deters( local_config ); // Everyone has a full copy of configs
  if ( g_rank == 0 ) deters.printParams();
  int maxspatials = deters.maxSpatials();
  
  deters.printFinalDeterminants(); // just a switch for subsequent print calls    
  
// This next step processes in parallel the configs-> determinants
// using a naive distribution - parallelism, however, is optional
  
  pair<int,double> time;
  int mylo =1;
  int myhi =1;
  int numb = ( maxspatials / g_size ) + 1;
  
  mylo = (g_rank * numb ) + 1; // Can never be ZERO
  myhi = min( (g_rank + 1) * numb, maxspatials );
  
  // Fetch : data parallel
#if 0
  cout << "generateDistributedDeterminants: mylo myhi max are " << GA::nodeid() <<" "<<  mylo << " " << myhi << " " << maxspatials << endl;
#endif

  if ( mylo <= maxspatials ) {
    deters.fetchConfigs( mylo, myhi, time );
  } 

  if (GA::nodeid() == 0 ) cout << "I am " << time.first << " Time to process deters is " << time.second << endl;
  
  deters.distributeParameters();
  
//TODO looks good to here
  if ( GA::nodeid() == 0 ) {
  cout << "Print spin parity is " << deters.fetchSpinParity() << endl;
  cout << "num spatials readin by PsociConfigs method is " << deters.maxSpatials() << endl;
  cout << "Max num actually processed globally are " << deters.fetchGlobalMaxSpatials() << endl;
  }
  
  if (deters.maxSpatials() != deters.fetchGlobalMaxSpatials()) {
    cout << "WARNING: it seems not all spatials were processed into determinants " << endl;
  }
  
  if ( maxspatials < MAX_LINES_PRINT ) {
      for(int i=0; i<= myhi-mylo; ++i ) {
	deters.printSpatialsData( i+1 ); // NOTE: a relative index ( per core )
    }
  }

  if ( GA::nodeid() == 0 ) cout << endl;
  localsize += ( myhi - mylo + 1 ) * sizeof( int ); 
  
  // Create GA space and push data to it.
  
  pair<int,double> uptime;
  createDistributedDeterminants( &deters, uptime ); // GA created and all local dets are pushed to GA
  if (GA::nodeid() == 0 ) cout << "I am " << time.first << " Time to process and UPLOAD deters is " << uptime.second+time.second << endl;

  //Basic output for tests

  printDistribution();

#if 0
  printDimensions();
  cout << " Me is " << GA::nodeid() << " FETCHED global sym is " << l_ksym << endl;
  cout << " Me is " << GA::nodeid() << " FETCHED maxsef is " << l_maxsef << endl;
  cout << " Me is " << GA::nodeid() << " FETCHED testconf " << l_maxSpatials << endl;
#endif

  deters.tearDownDeterminants();// We can delete for NON-direct methods
  
  return( localsize );
}

/* 
   Small memory version to fetch configs and build determinants and push them to GA
   This is a transition method that will bew refactored after the next paper ( yea, right!)

   In this case, we presume ONLY rank = 0 has the populated configs array.
*/ 
// THIS is not too useful it ends up being VERY slow

long PsociGADeterminants::generateDistributedDeterminantsSmallMemory( PsociConfigs & configs )
{
  long localsize = 0;
  int g_rank = GA::nodeid();
  int g_size= GA::nodes();
  
  //First step create distributed (local) determinants
  
  GA::error("generateDistributedDeterminantsSmallMemory: not recommended: Aborting", -1);
  
#ifndef SINGLENODECONFIGS
  GA::error("This method should only be used if code is compiled with -DSINGLENODECONFIGS",-1);
#endif
  if ( GA::nodeid() == 0 ) {
    cout << "generateDistributedDeterminants: single node configs processing " << endl;
    cout << "MAXWIDTH_CONFIGS is " << MAXWIDTH_CONFIGS << endl;
  }
  
  PsociDeterminants deters( &configs ); // Not everyone has the configs BUT we all need objects
  
  if ( GA::nodeid() == 0 ) {
    deters.findParams(); // This reprocesses the entire list 
  }
  
  deters.distributeParameters();
  l_maxSpatials =  deters.maxSpatials();
  
  deters.overrideGlobalMaxSpatials( l_maxSpatials );
  if ( GA::nodeid() == 0 ) {
    cout << "Print spin parity is " << deters.fetchSpinParity() << endl;
    cout << "num spatials readin by PsociConfigs method is " << deters.maxSpatials() << endl;
    cout << "Max num actually processed globally are " << deters.fetchGlobalMaxSpatials() << endl;
  }
  
  if ( g_rank == 0 ) deters.printParams();
  
  int maxspatials = l_maxSpatials; // okay to here
  deters.printFinalDeterminants(); // just a switch for subsequent print calls    
  
  pair<int,double> time;
  
  vector<JOUTFG> * poutfg;
  initDistributedDeterminants( &deters ); 
  
  double timein;
  double det_time=0.0, uptime=0.0;
  
  int width = 0;
  
  if ( GA::nodeid() == 0 ) {
    
    poutfg = deters.fetchJoutfg();
    for(int i=1; i<=maxspatials; i+= MAXWIDTH_CONFIGS ) // NOTE inclusive
      {
	int istart = i;
	width += MAXWIDTH_CONFIGS;
	if ( width >  maxspatials) width = maxspatials;
	int istop = width;
	timein = psociTime();
	poutfg->clear(); //Not the best way but pretty easy to reuse existing ciode
	deters.specialFetchConfigs( istart, istop );
	det_time += psociTime() - timein;
	timein = psociTime();
	uploadSubsetDeterminants(); // DIRECT ifdefs inside
	uptime += psociTime() - timein;
      }
  } // rank == 0
  
  if (GA::nodeid() == 0 ) cout << "I am " << GA::nodeid() << " SmallMemory: Time to process and UPLOAD deters is " << det_time<<" "<<uptime<< endl;
  
  if (deters.maxSpatials() != deters.fetchGlobalMaxSpatials()) {
    cout << "WARNING: it seems not all spatials were processed into determinants " << endl;
    cout << "data are " << deters.maxSpatials() <<" "<<deters.fetchGlobalMaxSpatials() << endl;
  }
  
  GA::sync();
  deters.tearDownDeterminants();// We can delete for all methods 
  return( localsize );
}

//Timed wrappers
long PsociGADeterminants::generateDistributedDeterminantsSmallMemoryGA(pair<int,double> & info )
{
  long status;
  int g_rank = GA::nodeid();
  double timein = psociTime();
  status = generateDistributedDeterminantsSmallMemoryGA(); 
  info.second = psociTime() - timein;
  info.first = g_rank;
  return( status );
}

/* 
   Small memory version to fetch configs and build determinants and push them to GA
   This uses new methods for reading directly into GA. Everyone potentially has 
   a small piece of configs ( or is empty)
   
   To use this compile the code with 
   
   -DSINGLENODECONFIGS 
*/ 
// Remember at this point everyone has only a subet of the configs data 
long PsociGADeterminants::generateDistributedDeterminantsSmallMemoryGA()
{
  long localsize = 0;
  int g_rank = GA::nodeid();
  int g_size= GA::nodes();
  
//First step create distributed (local) determinants

  if ( GA::nodeid() == 0 ) {
    cout << "generateDistributedDeterminantsSmallMemoryGA: single node configs processing " << endl;
  }
  
  //PsociDeterminants local_det-> local_config );  // does scope get lost here for DIRECT?

  local_det->findParams(); 
  local_det->distributeParameters(); // internal GOPs

  l_maxSpatials =  local_det->maxSpatials();
  local_det->overrideGlobalMaxSpatials( l_maxSpatials );
  
  local_det->printFinalDeterminants(); // just a switch for subsequent print calls
  
#if 0
  cout << "l_maxSpatials is " <<  local_det->maxSpatials() << endl;
  cout << "l_maxdetperconf " << l_maxdetperconf << endl;
#endif
  
  if ( g_rank  == 0 ) {
    cout << "SmallMemoryGA: Print spin parity is " << local_det->fetchSpinParity() << endl;
    cout << "SmallMemoryGA: num spatials readin by PsociConfigs method is " << local_det->maxSpatials() << endl;
    cout << "SmallMemoryGA: Max num actually processed globally are " << local_det->fetchGlobalMaxSpatials() << endl;
  }
  
  if ( g_rank == 0 ) local_det->printParams();
  
#ifdef REPLICATEDETVALUE
// Build replicated space for holding phase into
   if ( GA::nodeid() == 0 ) {
      cout << "Using replicated PHASE vector: Resize to " << l_maxSpatials << endl;
   }
   replicated_phases.resize( l_maxSpatials, 0 );
#endif
  
  pair<int,double> time;
  
  vector<JOUTFG> * poutfg;
  initDistributedDeterminants( local_det ); // everyone is neede for this
  
/* at this point everyone should read and process their set of configs
   also everyone can update their local values of replicated_phases: it will be
   GOPed later
*/

  double timein;
  double det_time=0.0, uptime=0.0;
  
  poutfg = local_det->fetchJoutfg();
  poutfg->clear(); //not the best way but pretty easy to reuse existing code
  
  timein = psociTime();
#ifdef REPLICATEDETVALUE
   local_det->specialDistributedFetchConfigs( replicated_phases );
#else
  local_det->specialDistributedFetchConfigs();
#endif
  det_time += psociTime() - timein;
  timein = psociTime();
  uploadSubsetDeterminants();
  uptime += psociTime() - timein;
  
  if (local_det->maxSpatials() != local_det->fetchGlobalMaxSpatials()) {
    cout << "WARNING: it seems not all spatials were processed into determinants " << endl;
    cout << "data are " << local_det->maxSpatials() <<" "<<local_det->fetchGlobalMaxSpatials() << endl;
  }

  /* On lots of nodes this output will need some kind of per-core organization 
   */
  
  if ( l_maxSpatials < MAX_LINES_PRINT ) {
    local_det->printLocalSpatialsData(); // 
  }
  GA::sync(); // we may be having a race condition here

/* NOW we need to GOP the replicated space*/
#ifdef REPLICATEDETVALUE
  char * op = "+";
  long Elems  = l_maxSpatials;

//#ifdef USEMPI
    MPI_Allreduce( MPI_IN_PLACE, &replicated_phases[0] , Elems, MPI_CHAR, MPI_SUM, MPI_COMM_WORLD );
//#else
//    GA::error("Cannot not use igop since replicate structure is now of size char",1);
//    GA::igop( &replicated_phases[0], Elems, op );
//#endif

    if ( GA::nodeid() == 0 ) cout << "Finished GOP on replicated_phase structure " << endl;

/*
    if ( GA::nodeid() == 0 ) {
cout << "XX replicate array is size " << replicated_phases.size() << endl;
for(int k=0; k< replicated_phases.size(); ++k ) {
cout << k<<" "<<replicated_phases[k]<<endl;
}
cout << endl;
}
*/


#endif


#ifndef DIRECT
/* But also no need to actually create GA spaces for dets
*/
  local_det->tearDownDeterminants();// We can delete for NON-direct methods
#endif
  
  return( localsize );
}

/* ON entry we expact only a few configs to be provided ( or none ) 
   Careful: list = (1,n) inclusive)
*/
long PsociGADeterminants::generateDistributedDeterminantsSmallMemoryDistributed( PsociConfigs & configs, vector<int> & list )
{

  long localsize = 0;
  int g_rank = GA::nodeid();
  int g_size= GA::nodes();
  
  //First step create distributed (local) determinant
  
  GA::error(" generateDistributedDeterminantsSmallMemoryDistributed disabled ",1);
#ifndef SINGLENODECONFIGS
  GA::error("This method should only be used if code is compiled with -DSINGLENODECONFIGS",-1);
#endif
  if ( GA::nodeid() == 0 ) {
    cout << "generateDistributedDeterminantsSmallMemoryDistributed: single node configs processing " << endl;
    cout << "MAXWIDTH_CONFIGS is " << MAXWIDTH_CONFIGS << endl;
  }
  
  PsociDeterminants deters( &configs ); // Not everyone has the configs BUT we all need objects
  
  if ( GA::nodeid() == 0 ) {
    deters.findParams(); // This reprocesses the entire list 
  }

  deters.distributeParameters();
  l_maxSpatials =  deters.maxSpatials();
  deters.overrideGlobalMaxSpatials( l_maxSpatials );

  //GA::sync(); // TODO REMOVE ME
  //cout << GA::nodeid() << " l_maxdetperconf " << l_maxdetperconf << " l_maxsefperconf " << l_maxsefperconf << endl;
  
  if ( GA::nodeid() == 0 ) {
    cout << "SmallMemoryDistributed:  Print spin parity is " << deters.fetchSpinParity() << endl;
    cout << "num spatials readin by PsociConfigs method is " << deters.maxSpatials() << endl;
    cout << "Max num actually processed globally are " << deters.fetchGlobalMaxSpatials() << endl;
  }
  
  if ( g_rank == 0 ) deters.printParams();
  
  deters.printFinalDeterminants(); // just a switch for subsequent print calls    
  
  pair<int,double> time;
  
  vector<JOUTFG> * poutfg;
  initDistributedDeterminants( &deters ); 
  
  double timein;
  double det_time=0.0, uptime=0.0;
  
  int lsize = list.size();
  
  if (lsize <= 0 ) { 
    cout <<GA::nodeid() << " return with empty list " << endl;
    return(0);
  }
  int width = list[ lsize - 1 ];
  
  poutfg = deters.fetchJoutfg();
  
  int istart = list[0];
  int istop = list[ width ];
  
  timein = psociTime();
  poutfg->clear(); //Not the best way but pretty easy to reuse existing ciode
  deters.specialFetchConfigs( istart, istop );
  det_time += psociTime() - timein;
  timein = psociTime();
#ifndef DIRECT
  uploadSubsetDeterminants();
#endif
  uptime += psociTime() - timein;
  
  //cout << "I am " << GA::nodeid() << " SmallMemoryGA: Time to process and UPLOAD deters is " << istop << " "<< det_time<<" "<<uptime<< endl;
  
  if (deters.maxSpatials() != deters.fetchGlobalMaxSpatials()) {
    cout << "WARNING: it seems not all spatials were processed into determinants " << endl;
    cout << "data are " << deters.maxSpatials() <<" "<<deters.fetchGlobalMaxSpatials() << endl;
  }
  
  GA::sync();
  
  deters.tearDownDeterminants();// We can delete for NON-direct methods
  return( localsize );
}

// Return local per-core memory 
/* A new version that presumes the deters data already exists in local_det
   This requires the alternative constructor
   
   New distributed approach
*/
long PsociGADeterminants::generateDistributedDeterminants()
{
  long localsize = 0;
  int g_rank = GA::nodeid();
  int g_size= GA::nodes();
  
  if ( g_rank == 0 ) local_det->printParams();
  int maxspatials = local_det->maxSpatials();
  
  local_det->printFinalDeterminants(); // just a switch for subsequent print calls    
  
  pair<int,double> time;
  int mylo;
  int myhi;
  int numb = ( maxspatials / g_size ) + 1;
  
  mylo = (g_rank * numb ) + 1; // Can never be ZERO
  myhi = min( (g_rank + 1) * numb, maxspatials );
  
  // Fetch : data parallel
#if 0
  cout << "mylo myhi max are " << mylo << " " << myhi << " " << maxspatials << endl;
#endif

  if ( mylo <= maxspatials ) {
    local_det->fetchConfigs( mylo, myhi, time );
  } 

  if (GA::nodeid() == 0 ) cout << "I am " << time.first << " Time to process Alternative local_det->is " << time.second << endl;
  
  local_det->distributeParameters();
  
  if ( GA::nodeid() == 0 ) {
  cout << "Print spin parity is " << local_det->fetchSpinParity() << endl;
  cout << "num spatials readin by PsociConfigs method is " << local_det->maxSpatials() << endl;
  cout << "Max num actually processed globally are " << local_det->fetchGlobalMaxSpatials() << endl;
  }
  
  if (local_det->maxSpatials() != local_det->fetchGlobalMaxSpatials()) {
    cout << "WARNING: Alternative: it seems not all spatials were processed into determinants " << endl;
  }
  
  if ( maxspatials < MAX_LINES_PRINT ) {
      for(int i=0; i<= myhi-mylo; ++i ) {
	local_det->printSpatialsData( i+1 ); // NOTE: a relative index ( per core )
    }
  }
  if ( GA::nodeid() == 0 ) cout << endl;
  localsize += ( myhi - mylo + 1 ) * sizeof( int ); 
  
  // Create GA space and push data to it.
  
  createDistributedDeterminants( local_det); // GA created and all local dets are pushed to GA
  
  printDistribution();
#if 0
  printDimensions();
  cout << " Me is " << GA::nodeid() << " ALT: FETCHED global sym is " << l_ksym << endl;
  cout << " Me is " << GA::nodeid() << " ALT: FETCHED maxsef is " << l_maxsef << endl;
  cout << " Me is " << GA::nodeid() << " ALT: FETCHED testconf " << l_maxSpatials << endl;
#endif
  
  return( localsize );
}

//Local
void PsociGADeterminants::setParameters( PsociDeterminants * deters )
{
  //These all need to be set prior to tearing down deters->
  
  l_ksym = deters->fetchGlobalSymmetry();
  l_maxSpatials = deters->fetchGlobalMaxSpatials();
  l_maxsef = deters->fetchGlobalMaxSef();
  l_maxdet =  deters->fetchGlobalMaxDet();
  l_maxsefperconf = deters->fetchGlobalMaxSefPerSpatial();
  l_maxdetperconf =  deters->fetchGlobalMaxDetPerSpatial();
  l_nelec = deters->fetchNumElectrons();
  l_nbfn = deters->fetchNumBasisFunctions();
  l_mxopn = deters->fetchMaxOpenShells();
}

//Note index is the Fortran style indexing ( 1<=index<=nconfigs) 
//We want to return the result as a "string" object
int PsociGADeterminants::constructConfigurationData( JOUTFG & joutfg, string & configs )
{
  int nbf = joutfg.nbf;

// Empty and clear the configs list.

  vector<int> doubly = joutfg.doubly;
  vector<int> singly = joutfg.singly;
  vector<int>::iterator it;
  
  configs.clear(); // Empty 
  configs.resize(nbf,'0'); // Append nbf worth of 0 
  char num2='2';
  char num1='1';
  
  for(it=doubly.begin(); it!=doubly.end(); ++it) {
    configs[ (*it) ] = num2;
  }
  for(it=singly.begin(); it!=singly.end(); ++it) {
    configs[ (*it) ] = num1;
  }

#if 0
  cout << "construct configs: INDEX " << index << endl;
  cout << "CONFIGS is " << configs;
#endif
  return(0);
}

/* PREALLOCATE method: A new routine to take an existing JOUTFG and resize it. The plan being to 
   make it just big enough to handle ANY fetched GA-> joutfg object. 
   This way we do not need to keep checking and resizing which is becoming a bottleneck
*/
void PsociGADeterminants::preallocateTempJoutfgObject( JOUTFG & toutfg )
{
  int nelec = l_nelec;
  int mxndeti = l_maxdetperconf;
  int mxopn = l_mxopn;
  int mxocc = ( nelec + mxopn ) /2;
  int mxclosed = nelec / 2;
  
  // first 10 ints are allocated at constructor;

  toutfg.doubly.resize( mxclosed, 0 );
  toutfg.singly.resize( mxopn, 0 );
  toutfg.spin_projection.resize( mxndeti, 0 );
  toutfg.phases.resize( mxndeti, 0 );
  toutfg.sumsa.resize( mxndeti, 0 );
  toutfg.sumsb.resize( mxndeti, 0 );

  // open up 2D arrays
  toutfg.alpha_spins.resize( mxndeti, vector<int>(mxocc,0));
  toutfg.beta_spins.resize( mxndeti, vector<int>(mxocc,0));

  //open up associated arrays
  toutfg.nalpha.resize( mxndeti, 0);
  toutfg.nbeta.resize( mxndeti, 0); 
}

/* Independant way to openup GA space for determinants
   Only needed if using SINGLENODE configs style
*/
void PsociGADeterminants::initDistributedDeterminants( PsociDeterminants * deters )
{
  local_det = deters;
  local_det->assembleGlobalJobParameters(); //Ensure a consistent set of job parameters
  setParameters( deters ); // Everyone comes through here - set the relevent parameters
  
  int l_maxSpatials = deters->maxSpatials(); 
  int maxSpatials = l_maxSpatials;
  int mxseti = l_maxsefperconf;
  int mxdeti = l_maxdetperconf;
  int nelec = l_nelec;
  int nbf = l_nbfn;
  int mxopn = l_mxopn;
  
  int mxspatials = l_maxSpatials;
  
  if ( GA::nodeid() == 0 ) {
    cout << "initDistributedDeterminants " << endl;
    cout << "Generating new GA for determinants using: mxopn = " << mxopn;
    cout << "mxdeti = " << mxdeti << " nelec = " << nelec << " nbf = " << nbf << endl;
    cout << "mxseti = " << mxseti << " nelec = " << nelec << " nbf = " << nbf << endl;
    cout << "Maximum number of spatials is " << mxspatials << endl;
  }
  
  int mxocc = ( nelec + mxopn ) /2;
  int mxclosed = nelec / 2;
  
  l_maxlength = computeLength( mxdeti, mxocc, mxclosed, mxopn ); //maximum possible buffer length in ints
  
  if ( GA::nodeid() == 0 ) {
    cout << "INIT: The estimated amount of total GA memory (B) to allocate for Dets is " << sizeof(COMPRESS_SIZE) * maxSpatials * l_maxlength << endl;
    cout << "INIT: The estimated amount of per-core GA memory (MB) to allocate for Dets is " << ( sizeof(COMPRESS_SIZE) * maxSpatials * l_maxlength ) / (1000*1000*GA::nodes()) << endl;
  }
  
  int DATA_TYPE = C_INT;
  int DATA_NDIM = 2;
  int dims[2];
  dims[1] = l_maxlength;
  dims[0] = mxspatials;

  
  const int DISTRIBUTE_EVENLY = -1; //TODO need to optimize this at some point
  int chunk[2];

#ifndef DIRECT
 
  chunk[1] = dims[1]; //We want distribution by rows (length)
  chunk[0] = DISTRIBUTE_EVENLY;
  
  char local_title[] = "determinants per conf";
  
  g_det = GA::createGA( DATA_TYPE, DATA_NDIM, dims, (char *) local_title, chunk);
  
  char handleGAMessage[] = "PsociGADeterminants::createDistributedDeterminants";
  g_det->checkHandle( handleGAMessage);
  g_det->zero();
  
  dims[1] = 1;
  dims[0] = maxSpatials;
  chunk[1] = DISTRIBUTE_EVENLY;
  chunk[0] = DISTRIBUTE_EVENLY;
  
  char local_title2[] = "g_num_det ";
  g_det_num = GA::createGA( DATA_TYPE, DATA_NDIM, dims, (char *) local_title2, chunk);
  
  char handleGAMessage2[] = "PsociGADeterminants::createDistributedDeterminants:g_number";
  g_det_num->checkHandle( handleGAMessage2 );
  g_det_num->zero();
#endif

  dims[1] = 1;
  dims[0] = maxSpatials;
  chunk[1] = DISTRIBUTE_EVENLY;
  chunk[0] = DISTRIBUTE_EVENLY;
  
  char local_title3[] = "g_nsef ";
  g_nsef = GA::createGA( DATA_TYPE, DATA_NDIM, dims, (char *) local_title3, chunk);
  
  char handleGAMessage3[] = "PsociGADeterminants::createDistributedDeterminants:g_nsef";
  g_nsef->checkHandle( handleGAMessage3 );
  g_nsef->zero();
  
  DATA_TYPE = C_INT;
  DATA_NDIM = 2; 
  dims[1] = nbf+2; // added jsym and index
  dims[0] = mxspatials;
  if ( GA::nodeid() == 0 ) cout << "Create a GA-ORBMAP of dims " << dims[0] << " by " << dims[1] << endl;
  chunk[1] = dims[1]; //We want distribution by rows (length)
  chunk[0] = DISTRIBUTE_EVENLY;
    
  char local_orb_title[] = "orbital maps";
  g_orbmap = GA::createGA( DATA_TYPE, DATA_NDIM, dims, (char *) local_orb_title, chunk);

#ifndef DIRECT
  assembleDimensions(); // Everyone should have this
//  printDistribution();
#endif

  printDistribution();

#ifdef DIRECT
  if ( GA::nodeid() == 0 ) cout << "DIRECT methods:  NO local determinant space created " << endl;
#endif

}

void PsociGADeterminants::uploadSubsetDeterminants()
{
    int ndet_loaded = packAndLoadDeterminantData();
}
// We presume the calling app calls index startig at (1,2,n)
// we presume orbmap is already sized to nbf in length.
// orbmap2d[index][nbf]. Starts at 0
void PsociGADeterminants::fetchOrbMap( int index, vector<COMPRESS_SIZE> & orbmap )
{
    int lo[ 2 ];
    int hi[ 2 ];
    if ( orbmap.size() < l_nbfn + 2 ) orbmap.resize( l_nbfn + 2 ); // it gets sized back to pure nbf in computeConfigs

    lo[0] = index-1;
    hi[0] = index-1;
    lo[1] = 0;
    hi[1] = l_nbfn +1; // jsym and index
    int n = 1;
       g_orbmap->get( lo, hi, &orbmap[0], &n);
}

void PsociGADeterminants::fetchOrbMap( int index, COMPRESS_SIZE * orbmap )
{
    int lo[ 2 ];
    int hi[ 2 ];


    lo[0] = index-1;
    hi[0] = index-1;
    lo[1] = 0;
    hi[1] = l_nbfn +1; // jsym and index
    int n = 1;
       g_orbmap->get( lo, hi, orbmap, &n);
}


void PsociGADeterminants::fetchOrbMap( int index, vector<vector<COMPRESS_SIZE> > & orbmap2d )
{
    int lo[ 2 ];
    int hi[ 2 ];
    vector<COMPRESS_SIZE> orbmap(l_nbfn+2);

    lo[0] = index-1;
    hi[0] = index-1;
    lo[1] = 0;
    hi[1] = l_nbfn + 1; // jsym and index
    int n = 1;
    g_orbmap->get( lo, hi, &orbmap[0], &n);
    orbmap2d.push_back( orbmap );
}


// Local call
int PsociGADeterminants::uploadOrbMaptoGA( int index,  vector<COMPRESS_SIZE> orbmap )
{
  int length = orbmap.size();
  const int DATA_NDIM = 2;
  int lo[ DATA_NDIM ];
  int hi[ DATA_NDIM ];
  
  int indexG = index-1;
  lo[1] = 0;
  lo[0] = indexG; // 
  hi[1] = length-1;
  hi[0] = indexG;
  int n=1;
  g_orbmap->put( lo, hi, &orbmap[0], &n);
  return(0);
}

//compare orbmaps: same as code in the Hamiltonian methods.

/* Perhaps we can add additional conditions: symmetry; type T; whatever
*/
int PsociGADeterminants::compareOrbMap( vector<COMPRESS_SIZE> & imap, vector<COMPRESS_SIZE> &jmap )
{
  int jmarr, sumjex=0, sumiex=0;
  
cout << "IN coimpareOrbMap " << endl;
  for(int i=0; i< l_nbfn; ++i) {
    jmarr = imap[i] - jmap[i];
    if ( jmarr < 0 ) {
      if ( (sumjex-jmarr) > 2 ) {
	return( 0 );
      }
      
      for(int j=-1; j>=jmarr; --j) {
	++sumjex;
      }
    } else if ( jmarr > 0 ) {
      if ( (sumiex+jmarr) > 2 ) {
	return( 0 );
      }
      
      for(int j=1; j<=jmarr; ++j) {
	++sumiex;
      }
    }
  }
  return( 1 );
}

/* Only used for direct matrix-vector product experiments
*/
int PsociGADeterminants::compareOrbMapMVP( vector<COMPRESS_SIZE> & imap, vector<COMPRESS_SIZE> &jmap )
{
//compare move than just tyopes
  int jmarr, sumjex=0, sumiex=0;

  for(int i=0; i< l_nbfn; ++i) {
    jmarr = imap[i] - jmap[i];
    if ( jmarr < 0 ) {
      if ( (sumjex-jmarr) > 2 ) {
        return( 0 );
      }

      for(int j=-1; j>=jmarr; --j) {
        ++sumjex;
      }
    } else if ( jmarr > 0 ) {
      if ( (sumiex+jmarr) > 2 ) {
        return( 0 );
      }

      for(int j=1; j<=jmarr; ++j) {
        ++sumiex;
      }
    }
  }
  return( 1 );
}


//compare orbmaps: same as code in the Hamiltonian methods.
//USED only for CIDEN processing
int PsociGADeterminants::compareOrbMapCIDEN( COMPRESS_SIZE * imap, COMPRESS_SIZE * jmap )
{
  int jmarr, sumjex=0, sumiex=0;
  
  for(int i=0; i< l_nbfn; ++i) {
    jmarr = imap[i] - jmap[i];
    if ( jmarr < 0 ) {
      if ( (sumjex-jmarr) > 1 ) {
	return( 0 );
      }
      for(int j=-1; j>=jmarr; --j) {
	++sumjex;
      }
    } else if ( jmarr > 0 ) {
      if ( (sumiex+jmarr) > 1 ) {
	return( 0 );
      }
      
      for(int j=1; j<=jmarr; ++j) {
	++sumiex;
      }
    }
  }
  return( 1 );
}

/* compare orbmaps: same as code in the Hamiltonian methods.
   This is redundant to what occurs within sociBlocks. But it doesn't
   seem to impact speed and redunadnt checks are not a horrible thing
*/
int PsociGADeterminants::compareOrbMap( COMPRESS_SIZE * imap, COMPRESS_SIZE * jmap )
{
  int jmarr, sumjex=0, sumiex=0;
  
  for(int i=0; i< l_nbfn; ++i) {
    jmarr = imap[i] - jmap[i];
    
    if ( jmarr < 0 ) {
      if ( (sumjex-jmarr) > 2 ) {
	return( 0 );
      }
      
      for(int j=-1; j>=jmarr; --j) {
	++sumjex;
      }
    } else if ( jmarr > 0 ) {
      if ( (sumiex+jmarr) > 2 ) {
	return( 0 );
      }
      
      for(int j=1; j<=jmarr; ++j) {
	++sumiex;
      }
    }
  }
  return( 1 );
}

/* MODIFY here to winnow lists ofr direct matrix vector productsa
*/
int PsociGADeterminants::compareOrbMapMVP( COMPRESS_SIZE * imap, COMPRESS_SIZE * jmap )
{
  int jmarr, sumjex=0, sumiex=0;

  for(int i=0; i< l_nbfn; ++i) {
    jmarr = imap[i] - jmap[i];

    if ( jmarr < 0 ) {
      if ( (sumjex-jmarr) > 2 ) {
        return( 0 );
      }

      for(int j=-1; j>=jmarr; --j) {
        ++sumjex;
      }
    } else if ( jmarr > 0 ) {
      if ( (sumiex+jmarr) > 2 ) {
        return( 0 );
      }

      for(int j=1; j<=jmarr; ++j) {
        ++sumiex;
      }
    }
  }
  return( 1 );
}


/* Experimental approach: for the input index; i , return the list of 
   all j's that are potential interaction candidates.
   currently we always start at "1"
   indexing is 'C-style'
*/ 

int PsociGADeterminants::compareAllOrbMap( int index, vector<COMPRESS_SIZE> & imap, vector<int> & list ) 
{
  list.clear();
  int lo[ 2 ];
  int hi[ 2 ];
  
#ifdef SQUAREHAMILTONIAN
  int last = l_maxSpatials;
#else
  int last = index;
#endif
  
/* Overly complicated and not any faster: Moreover, 
   this doesn't always work without a segfaulrt
*/
  double gettime=0.0, compttime=0.0;

  double tempt = psociTime();

  vector<COMPRESS_SIZE> jmap( last * (l_nbfn+2) );

  lo[0] = 0;
  hi[0] = last - 1;
  lo[1] = 0;
  hi[1] = l_nbfn + 1 ;

  int n = l_nbfn + 2 ;
  //int n = 1;

  g_orbmap->get( lo, hi, &jmap[0], &n);

  gettime += psociTime() - tempt;
  tempt = psociTime();

  int linear_index=0;

  for(int jindex=0; jindex< last;++jindex) {
    if ( compareOrbMap( &imap[0], &jmap[linear_index]) != 0 ) { // 
      list.push_back( jindex );
    }
    linear_index += (l_nbfn+2); 
  }
  compttime += psociTime() - tempt;

  return( list.size()  );
}

/* attempt to support direct matrix-vector products
   on entry testvector is the portion of the total vector that woiuld be indexes 
   by the current set of H[i][*] terms. If all zero ( or near) then there is no 
   need to bother with this.

*/

int PsociGADeterminants::compareAllOrbMapMatrixVectorProducts( int index, vector<COMPRESS_SIZE> & imap, vector<int> & list ) 
{
  list.clear();
  int lo[ 2 ];
  int hi[ 2 ];
  
#ifdef SQUAREHAMILTONIAN
  GA::error("compareOrbMapMatrixVectorProducts not usable with SQUAREHAMITONOAHN",1);
#endif

  int last = index;
  double gettime=0.0, compttime=0.0;
  double tempt = psociTime();

  vector<COMPRESS_SIZE> jmap( last * (l_nbfn+2) );
  lo[0] = 0;
  hi[0] = last - 1;
  lo[1] = 0;
  hi[1] = l_nbfn + 1 ;

  int n = l_nbfn + 2 ;
  //int n = 1;

  g_orbmap->get( lo, hi, &jmap[0], &n);

  gettime += psociTime() - tempt;
  tempt = psociTime();

  int linear_index=0;

  for(int jindex=0; jindex< last;++jindex) {
    if ( compareOrbMapMVP( &imap[0], &jmap[linear_index]) != 0 ) { // 
      list.push_back( jindex );
    }
    linear_index += (l_nbfn+2); 
  }
  compttime += psociTime() - tempt;
 // compttime += psociTime() - tempt;

  return( list.size()  );
}

/* prune the list to contain values between ilow and ihi-1 */
/* Recall we must be careful regarduing conf versus sefs. */

// ilow ihi are in CONFIGURATION space
int PsociGADeterminants::pruneListVectorPatch( int ilow, int ihi, vector<int> & list ) 
{
#ifdef SQUAREHAMILTONIAN
  GA::error("pruneListVectorPatch not usable with SQUAREHAMITONOAHN",1);
#endif

  int size = list.size();

  vector<int> tempList;

//can probably use a break here
// ilow ihi are conf lists
//  cout<< GA::nodeid() << " list size on entry is " << list.size()<<" ranges " << ilow<<" "<<ihi  << endl;

//  for( vector<int>::iterator it = list.begin(); it != list.end(); ++it ){
//     int ival = *it;
//        cout << GA::nodeid() << " BEFORE kept list " << ival<<" ilof ihi are " << ilow<<" "<<ihi  << endl;
//    }

  for( vector<int>::iterator it = list.begin(); it != list.end(); ++it ){
     int ival = *it; 

/*
     if ( ival < ilow || ival >= ihi ) {
        list.erase(it);
     } else {
        cout << GA::nodeid() << " kept a J value " << ival << endl;
        ++it;
     }
*/
     if ( ival >= ilow &&  ival < ihi ) {
     //cout << GA::nodeid() << "KEEP this list value " << ival << endl;
        tempList.push_back( ival );
       }
     }

//cout <<" Try to clear the list " << list.size() << endl;
     list.clear();
//cout << GA::nodeid() <<" tempList size is " << tempList.size() << endl;
     for(int i=0; i< tempList.size(); i++) {
        list.push_back( tempList[i] );
     }
     tempList.clear();

 //cout << GA::nodeid() << " My new list is of size " << list.size() << endl;

  return( list.size()  );
}

/* attempt to support direct matrix-vector products
   on entry testvector is the portion of the total vector that woiuld be indexes 
   by the current set of H[i][*] terms. If all zero ( or near) then there is no 
   need to bother with this.

   Pre limit search space tothos in jphase
   jphase has already been limited to values <= iconf
*/

int PsociGADeterminants::compareAllOrbMapMatrixVectorProducts( int index, vector<COMPRESS_SIZE> & imap, vector<int> & jphase, vector<int> & list ) 
{
  list.clear();
  int lo[ 2 ];
  int hi[ 2 ];
 
#ifdef SQUAREHAMILTONIAN
  GA::error("compareOrbMapMatrixVectorProducts not usable with SQUAREHAMITONOAHN",1);
#endif

  int last = index;
  double gettime=0.0, compttime=0.0;
  double tempt = psociTime();

  vector<COMPRESS_SIZE> jmap( l_nbfn+2 );


  for(vector<int>::iterator jit=jphase.begin(); jit!=jphase.end(); ++jit) {
    lo[0] = (*jit);
    hi[0] = (*jit);
    lo[1] = 0;
    hi[1] = l_nbfn + 1 ;
    int n = 1;

    g_orbmap->get( lo, hi, &jmap[0], &n);
    if ( compareOrbMapMVP( &imap[0], &jmap[0]) != 0 ) { // 
      list.push_back( *jit );
    }    
  }

  compttime += psociTime() - tempt;

  return( list.size()  );
}

/* J subset procesing method. Actual jmaps have been prefetched and jsubsetlist is
   perfectly wide

   we simply return the total possible number of matches. we do not account for potential triangular sparsity.
*/
int PsociGADeterminants::compareLocalJSubsetOrbMap(int numJs, int index, vector<COMPRESS_SIZE> & imap, vector<int> &  jindexlist, vector<COMPRESS_SIZE> & jsubsetlist, vector<int> & list ) 
{
  list.clear();
  int lo[ 2 ];
  int hi[ 2 ];
  
#ifdef SQUAREHAMILTONIAN
  int last = l_maxSpatials;
#else
  int last = index;
#endif

  double gettime=0.0, compttime=0.0;
  double tempt = psociTime();
  int linear_index=0;


// correct to here

 // if ( jindexlist[0]+1 > index ) return(0); // broadbased check current lowest J index

  for(int jindex=0; jindex< numJs;++jindex) { // compare to the rest
    if ( jindexlist[jindex]+1 <= last ) {
           if ( compareOrbMap( &imap[0], &jsubsetlist[linear_index]) != 0 ) { // 
              //cout << "SubsetOr numJs " << numJs<<" "<<index <<" " << jindex << " " << jindexlist[jindex]+1 << endl;
              list.push_back( jindexlist[jindex] );
           }
      linear_index += (l_nbfn+2); 
    }
  }
  compttime += psociTime() - tempt;

  //cout << "subset compareLocalJSubsetOrbMapcomptiume " << compttime << endl;
  return( list.size()  );
}


/* In this version prior to fetching a jmap do a preliminary check using rep_phase
*/



int PsociGADeterminants::compareLocalJSubsetOrbMap(int numJs, int index, vector<char> * rep_phase, vector<COMPRESS_SIZE> & imap, vector<int> &  jindexlist, vector<COMPRESS_SIZE> & jsubsetlist, vector<int> & list ) 
{
  list.clear();
  int lo[ 2 ];
  int hi[ 2 ];
#ifdef SQUAREHAMILTONIAN
  GA::error("PsociGADeterminants::compareLocalJSubsetOrbMap: No use of -DSQUAREHAMILTONIAN",1);
#endif

  int last = index;
  double gettime=0.0, compttime=0.0;
  double tempt = psociTime();
  int linear_index=0;

  int iopen = rep_phase->at(index-1);
  int count_phase=0;

 // correct to here

 // if ( jindexlist[0]+1 > index ) return(0); // broadbased check current lowest J index

  for(int jindex=0; jindex< numJs;++jindex) { // compare to the rest
    int jopen = rep_phase->at(jindexlist[jindex]);
//    cout << "iindex jindex are " << index-1 << " " << jindexlist[jindex] << " " << iopen<<" "<<jopen<<endl;

        if ( (jindexlist[jindex]+1 <= last) ) {
           if ( abs(iopen-jopen) <= 4 ) {
           if ( compareOrbMap( &imap[0], &jsubsetlist[linear_index]) != 0 ) { // 
              list.push_back( jindexlist[jindex] );
           }
          } 
           //} else {
           //  ++count_phase;
           // }
        linear_index += (l_nbfn+2); 
       }
  }
/*
  if ( count_phase != 0 ) {
     cout << "compareLocalJSubsetOrbMap: assemb screened outcnfigs with " << count_phase << endl;
  } 
*/

  compttime += psociTime() - tempt;
  return( list.size()  );
}
/* 1e CIdensity version  that looks for possible interacting orbmaps that are only relevant for a 1e H*/
int PsociGADeterminants::compareLocalJSubsetOrbMapCIDEN(int numJs, int index, vector<COMPRESS_SIZE> & imap, vector<int> &  jindexlist, vector<COMPRESS_SIZE> & jsubsetlist, vector<int> & list ) 
{
  list.clear();
  int lo[ 2 ];
  int hi[ 2 ];
  
#ifdef SQUAREHAMILTONIAN
  int last = l_maxSpatials;
#else
  int last = index;
#endif

  double gettime=0.0, compttime=0.0;
  double tempt = psociTime();
  int linear_index=0;

// correct to here
//linear index must be wrong

 // if ( jindexlist[0]+1 > index ) return(0); // broadbased check current lowest J index

  for(int jindex=0; jindex< numJs;++jindex) { // compare to the rest

    if ( jindexlist[jindex]+1 <= last ) {
           if ( compareOrbMapCIDEN( &imap[0], &jsubsetlist[linear_index]) != 0 ) { // 
              //cout << "CIDEN SubsetOr numJs " << numJs<<" "<<index <<" " << jindex << " " << jindexlist[jindex]+1 << endl;
              list.push_back( jindexlist[jindex] );
           }
      linear_index += (l_nbfn+2); 
    }
  }
  compttime += psociTime() - tempt;

  //cout << "subset compareLocalJSubsetOrbMapCIDENcomptiume " << compttime << endl;
  return( list.size()  );
}

int PsociGADeterminants::compareLocalJSubsetOrbMapCIDEN(int numJs, int index, vector<char> * rep_phase, vector<COMPRESS_SIZE> & imap, vector<int> &  jindexlist, vector<COMPRESS_SIZE> & jsubsetlist, vector<int> & list ) 
{
  list.clear();
  int last = index;
  double gettime=0.0, compttime=0.0;
  double tempt = psociTime();
  int linear_index=0;

  int iopen = rep_phase->at(index-1);
  int count_phase=0;

 // correct to here

 // if ( jindexlist[0]+1 > index ) return(0); // broadbased check current lowest J index

  for(int jindex=0; jindex< numJs;++jindex) { // compare to the rest
        if ( (jindexlist[jindex]+1 <= last) ) {
           int jopen = rep_phase->at(jindexlist[jindex]);

           if ( abs(iopen-jopen) <= 2 ) {
           if ( compareOrbMapCIDEN( &imap[0], &jsubsetlist[linear_index]) != 0 ) { // 
              list.push_back( jindexlist[jindex] );
           }
          }
           //} else {
           //  ++count_phase;
           // }
        linear_index += (l_nbfn+2); 
       }
  }
  compttime += psociTime() - tempt;
  return( list.size()  );
}


/* Experimental approach: for the input index; i , return the list of 
   all j's that are potential interaction candidates.
   currently we always start at "1"
   indexing is 'C-style'

   USED for CIDEN related work
*/ 
int PsociGADeterminants::compareAllOrbMapCIDEN( int index, vector<COMPRESS_SIZE> & imap, vector<int> & list ) 
{
  list.clear();
  int lo[ 2 ];
  int hi[ 2 ];
  
  int last = index;
  vector<COMPRESS_SIZE> jmap( last * (l_nbfn+2) );
  lo[0] = 0;
  hi[0] = last - 1;
  lo[1] = 0;
  hi[1] = l_nbfn +  1 ;
  int n = l_nbfn + 2;
  //int n = 1;
  g_orbmap->get( lo, hi, &jmap[0], &n);
  
  int linear_index=0;
  for(int jindex=0; jindex< last;++jindex) {
    if ( compareOrbMapCIDEN( &imap[0], &jmap[linear_index]) != 0 ) {
      list.push_back( jindex );
    }
    linear_index += (l_nbfn+2); 
  } 
  return( list.size()  );
}

/* Now we want nodes to download ALL orbital info but only store the non-zero terms
   Then we need an alternative method for comparing lists.

   Loop over all l_maxSpatials terms. We MAY need to have each core do this independantly 
   because we would need to FIX the size of the object before sending.

   Instead of having everyone do a GA get, have one core do the get and then BRDCST it
   to the rest
*/
// Collective operation

int PsociGADeterminants::fetchAllCompressedOrbMap( vector<vector<pair<short,short> > > & ijpairlist )
{
//    vector<vector<pair<short,short> > > ijpairlist( l_maxSpatials, vector<pair<short,short> > );
    if ( l_nbfn > 512 ) GA::error("fetchAllCOmpressedOrbMap only works for nbf <= 512 ",1);

#ifdef DETAILEDCHECK
    cout<< GA::nodeid() << " Entry to fetchAllCOmpressedOrbMap with ijpairlist size of " << ijpairlist.size() << endl;
#endif

    int g_rank = GA::nodeid();
    int lo[ 2 ];
    int hi[ 2 ];
    
    lo[0] = 0;
    hi[0] = 0;
    lo[1] = 0;
    hi[1] = l_nbfn + 1 ;
    int n = 1;

    vector<COMPRESS_SIZE> orbmap( l_nbfn+2 );
    for(int iconf=0; iconf < l_maxSpatials; ++iconf ) {
    lo[0] = iconf;
    hi[0] = iconf;
    if ( g_rank == 0 ) g_orbmap->get( lo, hi, &orbmap[0], &n);
    GA::brdcst( &orbmap[0], l_nbfn*sizeof(COMPRESS_SIZE), 0); 

    orbmapRemoveZeros( orbmap, ijpairlist[iconf] ); // appends to the list 
    }
    return( ijpairlist.size() );
}

/* compare orbmaps: same as code in the Hamiltonian methods.
   This is redundant to what occurs within sociBlocks. But it doesn't
   seem to impact speed and redunadnt checks are not a horrible thing

   if NBF > 255 we have a problem
   ALso we strip off jsym and iconf ( elements nbf, nbf+1)
*/
int PsociGADeterminants::orbmapRemoveZeros( vector<COMPRESS_SIZE> & orbmap, vector<pair<short,short> > & densepairlist )
{ 
  pair<short,short> temp;
  for(int i=0; i< l_nbfn; ++i) {
     if ( orbmap[i] != 0 ) {
            temp.first=(short) i;
            temp.second=(short) orbmap[i];
            densepairlist.push_back( temp );
     }
  }
  return( densepairlist.size() );
}

/* compare method using removedZeros orbmap
   This takes the configuration index ( index; 1,2,n ) and determines all Js (0,1,2,n-1) not identically zero
   The Js are returned in list and the number of Js is returned as the return int.

   Sorry about the confusion in indexing. But in building H all indexes are refered to as index+1, so we must
   return Js as C-style

*/
int PsociGADeterminants::compareAllOrbMapRemoveZeros( int index, vector<vector<pair<short,short> > > & ijpairlist, vector<int> & list )
{
  list.clear();

#ifdef SQUAREHAMILTONIAN
  int last = l_maxSpatials;
#else
  int last = index;
#endif

  for(int jindex=0; jindex< last;++jindex) {
    if ( compareOrbMapRemoveZeros( index-1, jindex, ijpairlist ) != 0 ) {
      list.push_back( jindex );
    }
  }
  return( list.size() ) ;
}

int PsociGADeterminants::compareAllOrbMapRemoveZerosCIDEN( int index, vector<vector<pair<short,short> > > & ijpairlist, vector<int> & list )
{
  list.clear();

  int last = index;

  for(int jindex=0; jindex< last;++jindex) {
    if ( compareOrbMapRemoveZerosCIDEN( index-1, jindex, ijpairlist ) != 0 ) {
      list.push_back( jindex );
    }
  }
  return( list.size() ) ;
}



/* compare orbmaps: same as code in the Hamiltonian methods.
   but using compressed ( no zero) data objects
   This is basically redundant to what occurs within sociBlocks. 

   on entry iindex, jindex  are C_style (0,n-1)
*/
int PsociGADeterminants::compareOrbMapRemoveZeros( int iindex, int jindex, vector<vector<pair<short,short> > > & ijpairs )
{
  vector<COMPRESS_SIZE> imap(l_nbfn+2,0);
  vector<COMPRESS_SIZE> jmap(l_nbfn+2,0);

  vector<pair<short,short> >::iterator it;
 
//explode into a array of size nbf+2

   for(it=ijpairs[iindex].begin(); it != ijpairs[iindex].end(); ++it ) {
     imap[ (*it).first ] = (COMPRESS_SIZE) (*it).second;
   }

   for(it=ijpairs[jindex].begin(); it != ijpairs[jindex].end(); ++it ) {
      jmap[ (*it).first ] = (COMPRESS_SIZE)(*it).second;
   }
   return( compareOrbMap( &imap[0], &jmap[0] ) );
}


int PsociGADeterminants::compareOrbMapRemoveZerosCIDEN( int iindex, int jindex, vector<vector<pair<short,short> > > & ijpairs )
{
  vector<COMPRESS_SIZE> imap(l_nbfn+2,0);
  vector<COMPRESS_SIZE> jmap(l_nbfn+2,0);

  vector<pair<short,short> >::iterator it;

//explode into a array of size nbf

   for(it=ijpairs[iindex].begin(); it != ijpairs[iindex].end(); ++it ) {
      imap[ (*it).first ] = (*it).second;
   }

   for(it=ijpairs[jindex].begin(); it != ijpairs[jindex].end(); ++it ) {
      jmap[ (*it).first ] = (*it).second;
   }

   return( compareOrbMapCIDEN( &imap[0], &jmap[0] ) );
}

/* Not continued confusion with index: here we presume we';ve already started at '1'
*/
void PsociGADeterminants::computeConfigsGA( int index, vector<COMPRESS_SIZE> & orbmap, JOUTFG & toutfg )
{
     fetchOrbMap( index, orbmap );
     int nsefi = local_det->computeConfigs( orbmap, toutfg );
}




