/**************************************************************************************

* Copyright (c) 2010,2012 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: 

 History:



**************************************************************************************/
/**
 *   @file PsociTaskManager.C
 *
 */

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>

#include <ga++.h>

#ifdef TCGMSG 
#include <tcgmsg.h>
#endif

#include "PsociTaskManager.hpp"
GA::GlobalArray * g_nxtval;

PsociTaskManager::PsociTaskManager( int chunk, GA::GlobalArray * main_g_nxtval)
{
   g_nxtval = main_g_nxtval;
   l_nxtval_chunk = chunk ;

   nleft = 0;
   icount = 0;

}

// XXXXXXXXXXXXXXXXXXXXXXXX Dynamic processing methods

#ifndef TCGMSG 
void PsociTaskManager::initNxtask()
{
  // GA::GlobalArray * g_nxtval;
    
  //if ( GA::nodeid() == 0 ) cout << "initNxtask: Using GA-based readinc style  NXTVAL" << endl;
  int DATA_TYPE = C_LONG;
  int DATA_NDIM = 1;
  int dims[1];
  dims[0] = 1;
    
  int chunk[1];
  chunk[0] = -1;

  char local_title[] = "nxtval array ";
  g_nxtval = GA::createGA( DATA_TYPE, DATA_NDIM, dims, (char *) local_title, chunk);

  if ( GA::nodeid() == 0 ) {
    int lo[1];
    lo[0] = 0; 
    int hi[1];
    hi[0] = 0; 
    int n=1; 
    long iZero=0;
    g_nxtval->put( lo, hi, &iZero, &n);
  }   
  GA::sync(); 

  nleft = 0;
  icount = 0;
}       

// Initialize for a reverse looping of iconf
void PsociTaskManager::initNxtask( int topvalue )
{
  // GA::GlobalArray * g_nxtval;

  //if ( GA::nodeid() == 0 ) cout << "initNxtask: Using GA-based readinc style  NXTVAL" << endl;
  int DATA_TYPE = C_LONG;
  int DATA_NDIM = 1;
  int dims[1];
  dims[0] = 1;

  int chunk[1];
  chunk[0] = -1;

  char local_title[] = "nxtval array ";
  g_nxtval = GA::createGA( DATA_TYPE, DATA_NDIM, dims, (char *) local_title, chunk);

  if ( GA::nodeid() == 0 ) {
    int lo[1];
    lo[0] = 0;
    int hi[1];
    hi[0] = 0;
    int n=1;
    //int iVal= 1 + ( topvalue / l_nxtval_chunk ); // Top block not necc top values
    long iVal=0;
    iVal = ( topvalue % l_nxtval_chunk == 0 ) ? topvalue / l_nxtval_chunk  : 1 + ( topvalue / l_nxtval_chunk ); 
   // int iVal=topvalue;
    g_nxtval->put( lo, hi, &iVal, &n);
  }
  GA::sync();

  nleft = 0;
  icount = 0;
}

void PsociTaskManager::destroyNxtask()
{       
  g_nxtval->destroy();
}       
        
int PsociTaskManager::fetchNxtvalChunk() 
{     
    return( l_nxtval_chunk );
}   
    
void PsociTaskManager::setNxtvalChunk( int nchunk )
{ 
    if ( GA::nodeid() == 0 ) cout << "initNxtask: setting CHUNK value " << nchunk << endl;
    l_nxtval_chunk = nchunk;
}

long PsociTaskManager::nxtask( int nproc )
{
  const long ichunk = l_nxtval_chunk;
  //TODO this should be a LONG
  long sub=0;
  const long inc = 1;

  if ( nproc > 0 ) {
    if ( nleft == 0 ) {
      icount = g_nxtval->readInc( &sub, inc ) * ichunk;
      nleft = ichunk;
      //cout << "forward task " << GA::nodeid()<<" "<<nleft<<" "<<ichunk<<" "<<sub<<" "<<icount << endl;
    }
    ++icount;
    --nleft;
    return( icount -1 );
  } else {
    nleft = 0;
    return( 0 );
  }
}

// This method CAN return
// negative values. We don;t want to clamp it at zero
// since zero is a valid result. The calling app
// will need to deal with the negative case
long PsociTaskManager::nxtaskRev( int nproc, int topval )
{
  const long ichunk = l_nxtval_chunk;
  //TODO this should be a LONG
  long sub=0;
  long ltopval = topval;
  //int sub = topval;
  const long inc = -1; // decrementing blocks not necc values
  long itop=0;

  if ( nproc > 0 ) {
    if ( nleft == 0 ) {
      itop = g_nxtval->readInc( &sub, inc ) * ichunk;
      //cout << "revtask inner " << GA::nodeid()<<" "<<" "<<sub<<" "<<itop << endl;
      icount = min( itop, ltopval );
      nleft = min( ichunk, icount-( itop-ichunk ) );
      //cout << "task " << GA::nodeid()<<" "<<nleft<<" "<<itop<<" "<<topval<<" "<<ichunk<<" "<<sub<<" "<<icount << endl;
      icount++;
    }
    --icount;
    --nleft;
    return( icount-1 );
  } else {
    nleft = 0;
    return( 0 );
  }
}


/* An alternative version. The original would use chunk to specify
   what set of indexes should be next provided to the caller.
   Now, I simply want to return the next BEGINNING value. The intervening
   terms are handled implicitly in the new 2WAYExp methods
*/
/* readInc as set up here ALWAYS begins at the lowest post incremental value
   so if chunk=100, the first nxtval entry is 100. But, we want to compare
   in the NXTVAL loop to inexttask=(0,1,2,3..). We can accomplish this by
   simply subtracting out l_nxtval_chunk fromk the return value
*/
long PsociTaskManager::nxtaskAlt( int nproc )
{
  int sub=0;
  const long inc = l_nxtval_chunk;
  int icount = g_nxtval->readInc( &sub, inc );
  cout << GA::nodeid() << "ALT NEW icount is " << icount << endl; // correct
  return( icount );
}
#endif

//a collective operation
void PsociTaskManager::nxtaskReset()
{
   icount = 0;
   nleft = 0;
  if ( GA::nodeid() == 0 ) {
    int lo[1];
    lo[0] = 0; 
    int hi[1]; 
    hi[0] = 0; 
    int n=1; 
    int iZero=0;
    g_nxtval->put( lo, hi, &iZero, &n);
  }   
  GA::sync();
}

//a collective operation
void PsociTaskManager::nxtaskReset( int topval )
{
   icount = 0;
   nleft = 0;
  if ( GA::nodeid() == 0 ) {
    int lo[1];
    lo[0] = 0;
    int hi[1];
    hi[0] = 0;
    int n=1;
    int iVal=topval;
    g_nxtval->put( lo, hi, &iVal, &n);
  }
  GA::sync();
}

#ifdef TCGMSG
// Alternative NXTVAL routines that require TCGMSG to be build with GA
// If GA is not build correctly to do parallel nxtval, you'll see a message in stdout about a "silly distribution" 
// Initialize the routine

// You probably do not want this at all.
void PsociTaskManager::initNxtask()
{
  if ( GA::nodeid() == 0 ) cout << "Using real TCGMSG nxtval " << endl;
//TODO should we really cal nxtval here?
  nxtask( GA::nodeid() );

  nleft = 0;
  icount = 0;
}

//Wait until everyone is done.
void PsociTaskManager::destroyNxtask()
{
  nxtask( -(GA::nodeid()) );
}

// The NXTVAL routine - requires building GA5 --with-tcgmsg - and even that still doesn;t seem to do it
/* 
   Wrapper around nxtval to increase granularity
   Based on RJ Harrison tiny subroutine - cidbg4
*/
long PsociTaskManager::nxtask( int nproc )
{
  const int ichunk = l_nxtval_chunk;
  static int icount = 0, nleft = 0;
  if ( nproc > 0 ) {
    if ( nleft == 0 ) {
      icount = tcg_nxtval( nproc ) * ichunk;
      nleft = ichunk;
    }
    ++icount;
    --nleft;
    return( icount -1 );
  } else {
    nleft = 0;
    int junk = tcg_nxtval( -nproc );  //Dummy call
    return( 0 );
  }
}
#endif







