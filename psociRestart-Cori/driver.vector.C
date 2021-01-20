/**************************************************************************************
 * Copyright (c) 2010 RENCI.
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
// Now start testing for GA
#include <ga++.h>
//#include <GAServices.h>

#include <dra.h>
#define GA_DATA_TYPE C_DBL

//End ga++

#include "PsociTimer.hpp"
#include "PsociVector.hpp"
#include "PsociGArestart.hpp"
#include "PsociDRAservices.hpp"
#include "PsociDRArestart.hpp"
#include "PsociGAbasis.hpp"

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//
//  Start driver program 
//
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#define GCOUT if (g_rank==0) cout

void dumpVectors( vector<vector<double> > vecs ) {
  vector<vector<double> >::iterator outer, current;
  int64_t root_num=0;
  
  for( outer=vecs.begin(); outer!=vecs.end(); ++outer)
    {
      ++root_num;
      cout <<  root_num << endl;
      
      vector<double>::iterator it;
      cout << "Number of basis is " << (*outer).size() << endl;
      if ( (*outer).size() < MAXPRINTCOEFS ) {
	for (it=(*outer).begin(); it != (*outer).end(); ++it)
	  {
	    cout << (*it) << " "  ;
	  }
      }
      cout << endl;
    }
}

int main(int argc, char **argv)
{
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
  
#if Ustyle
  GCOUT << "PsociVector built to use the POSIX-style I/O" << endl;
#elif Cstyle
  GCOUT << "PsociVector built to use the C-style I/O" << endl;
#elif MPIstyle
  GCOUT << "PsociVector built to use the MPI-style I/O " << endl;
#else
  GCOUT << "PsociVector built to use the C++-style I/O (not recommended) " << endl;
#endif
  
  
  // Create a dummy array of eigenvectors for storage
  
  double base = 3.12345678;
  vector<vector<double> > testData;
  //  Simple test
  //  int64_t nroots = 3;
  //  int64_t nbasis = 10;
  
  // Try a realistic version 50 roots 1M basis
  
  int64_t nroots = 5;
  //int64_t nbasis = 1000000;
  int64_t nbasis = 10;
  
  double shift = 0.0;
  double fudge=2.0;

  int64_t index = -1;
  for( int64_t i=0; i< nroots; ++i) {
    shift += 1.0;
    double value;
    vector<double> root( nbasis, 0.0 );

    ++index;

    for ( int64_t nb=0; nb < nbasis; ++nb )
      {
        if ( nb == index ) {
	value = base*shift+fudge;;
        root[nb] = value/100.0;
        }
      }
    testData.push_back( root );
  }
  
  pair<int64_t,int64_t> info;
  info.first = nbasis;
  info.second = nroots; 
  
  cout << "Size of testData is " << testData.size() << endl;
  cout << "Start PsociVector methods " << endl;
  
  cout << endl << endl;
  
// Create a 2dim Global array of testData3: (nbfs x nroots+5)

  cout << endl << endl << "XXXXXXXXXXXXXXXXXXX test the GA restart routines " << endl;
  
  int VECTOR_TYPE = C_DBL;
  int VECTOR_NDIM = 2;
  int dims[2];
  cout << "info.first is " << info.first << " .second is " << info.second << endl;
  dims[0] = info.first;
  dims[1] = info.second+5; //Careful info is <int64_t,int64_t>
  cout << "Create a GA of dims " << dims[0] << " by " << dims[1] << endl;

  const int DISTRIBUTE_EVENLY = -1;

  int chunk[2];
  chunk[0] = DISTRIBUTE_EVENLY; //We want distribution across chunk0
  chunk[1] = DISTRIBUTE_EVENLY;
  
  int adims[2]; //A test extra vector set for appending....
  cout << "info.first is " << info.first << " .second is " << info.second << endl;
  adims[0] = info.first;
  adims[1] = 2; //Careful info is <int64_t,int64_t>

  string arrayName = "vector_array";
  //GA::GlobalArray *g_vector = GA::createGA( VECTOR_TYPE, VECTOR_NDIM, dims, (char *)"test", chunk); 
  GA::GlobalArray *g_vector = GA::createGA( VECTOR_TYPE, VECTOR_NDIM, dims, (char *)"test", NULL); 
  GA::GlobalArray *g_append = GA::createGA( VECTOR_TYPE, VECTOR_NDIM, adims, (char *)"test", NULL);
  g_vector->zero();

  int type_test;
  int ndims_test;
  int dims_test[2];
  
//Should be 20x10 at this point

  //int nbasis=info.first;
    vector<vector<double> >::iterator git;
    int lo[2];
    lo[0] = 0;
    lo[1] = 0;
    int hi[2];
    hi[0] = info.first-1;
    hi[1] = 0;
    
    int n = 1;
    
// Only grab the Last two vectors. Otherwise you'll get hit with orthonormality issues.
    /* Have only root zero do this part */
    if ( g_rank == 0 ) {
      git = testData.end();
      --git;
         g_append->put(lo, hi, (void *)&(*git)[0], &n);
      --git;
      ++lo[1];
      ++hi[1];
      g_append->put(lo, hi, (void *)&(*git)[0], &n);
    }
  //GA::sync();

  g_vector->printDistribution();

    lo[0] = 0;
    lo[1] = 0;
    hi[1] = 0;
    hi[0] = info.first-1;

//Fill in array for subsequent testing of appendVectors();

    if ( g_rank == 0 ) {
      for( git = testData.begin(); git != testData.end(); ++git)
        {
          cout << "XXXXX Put next vector " << endl;
          cout << "ROW lo vals " << lo[0] << " " << hi[0] << "COL lo vals " << lo[1] << " " << hi[1] << (*git)[0] << " " << (*git)[1] << " " << (*git)[19] << endl;
          g_vector->put(lo, hi, (void *)&(*git)[0], &n);
          ++lo[1];
          ++hi[1];
        }
    }
  //GA::sync();
  
  if ( nbasis <= 100 ) { 
    cout << "print array " << endl;
    g_vector->print();
  }
  
  // Test dumping GLobal Array test data to PsociVector disk.
  
  int alo[2];
  int ahi[2];

  /* works fine and as expected */
  /* every node does this */
  vector<double>  buf;
  buf.resize(nbasis);
  int ldv=1;
  alo[0]=0;
  alo[1]=0;
  ahi[1]=0;
  ahi[0]=nbasis-1;
  g_vector->get( alo, ahi, &buf[0], &ldv);
  double sum=0;
  for (int i=0; i<nbasis; ++i) {
    cout << "gotten data are " << buf[i] << endl;
    sum += buf[i] * buf[i];
  }
  cout << "Input vector overlaps as a check:  sum**2/sqrt(sum)1/sqrt(sum) are " << sum << " " << sqrt(sum) << " " << 1.0/sqrt(sum) << endl;
  //GA::sync();
  
  /* this one used to works */
    alo[0]=0;
    alo[1]=0;
    ahi[1]=0;
    ahi[0]=nbasis-1;
    
    double ddot1;
    ddot1 = g_vector->ddotPatch('N',alo,ahi, g_vector,'N', alo,ahi );
    cout << "rank is " << g_rank << "driver-ddot-test Value fo Vecvtor=1 ddot N*N is " << ddot1 <<  endl;
  
  // Try out new VECTOR routines
  
  PsociGAbasis vector_set( nroots, g_vector );
  vector<double> vecnorm;
  
  vector_set.getAllNorms( vecnorm );
  //GA::sync();
  for (int i=0; i<nroots; ++i) {
    cout << "rank is " << GA::nodeid() << " VECNORM data are " << vecnorm[i] << endl;
  }
  
  cout << "Now normalize vectors and print norms again " << endl;
  pair<int,double> time;
  vector<double> vecnorm2;
  vector_set.getAllNorms( vecnorm2, time);
  cout << "Times rank is " << time.first << " Norm time is " << time.second << endl;
  
  vector_set.normAllVectors();
  vector_set.getAllNorms( vecnorm );
  
    cout << "SECOND Size of vecnorm is " << vecnorm.size() << endl;
    for (int i=0; i<nroots; ++i) {
    cout << "SECOND VECNORM data are " << vecnorm[i] << endl;
    }
    g_vector->print();
  
  int or_status = vector_set.orthonormalizeAllVectors(); 
  if ( or_status != 0 ) {
    cout << "A linear dependancy was found at i " << or_status << endl;
    cout << "This currently is a hard failure: Aborting " << endl;
    GA::Terminate();
  }
  if ( nbasis < 100 ) {
    g_vector->print();
     }
  
  cout << " Check if all vectors are orthogonal " << endl;

  vector<pair<int,int> > overlap;
  if ( vector_set.checkAllOverlap( overlap ) != 0 )
    {
      cout << "Non-overlaps found " << endl;
    }
  
  // Now can we ADD a new vector to the g_vector set and orthonormalize it ?

  if ( nroots > 2 ) {
    int re_status = vector_set.resizeVectors( nroots - 2 );
  }
  cout << "vectors should now be smaller by TWO " << endl;
  vector<double> vecnorm3;
  vector_set.getAllNorms( vecnorm3);
  cout << "AFTER contraction: " << endl;
  cout << "Third Size of vecnorm is " << vecnorm3.size() << endl;
  for (int i=0; i<nroots-2; ++i) {
    cout << "Third VECNORM data are " << vecnorm3[i] << endl;
  }
  
  if ( nbasis < 100 ) { 
    g_vector->print();
  } 
  
  pair<int,double>  rtime;
  int re_status2 = vector_set.resizeVectors( nroots - 2 + 5, rtime );
  if ( re_status2 != 0 ) {
    cout << "A linear dependancy was found at i " << re_status2 << endl;
    cout << "This currently is a hard failure: Aborting " << endl;
    GA::Terminate();
  }
  
  cout << "Rank is " << g_rank << " resize time is " << rtime.second << endl;
  vector<pair<int,int> > overlap2;
  if ( vector_set.checkAllOverlap( overlap2 ) != 0 )
    {
      cout << "Non-overlaps found " << endl;
    }

//Renormalize a second time....
  or_status = vector_set.orthonormalizeAllVectors();
  if ( or_status != 0 ) {
    cout << "A linear dependancy was found at i " << or_status << endl;
    cout << "This currently is a hard failure: Aborting " << endl;
    GA::Terminate();
  }


  cout << "Contains augmented vectors " << endl;
  if ( nbasis < 100 ) {
    g_vector->print();
  }
  

  cout << " Test appendVectors() method " << endl;
  int ap_status = vector_set.resizeVectors( nroots - 2 );
  cout << "vectors should now be smaller by TWO ready for the append" << endl;
// Turn me on to show append FAIL
  int re_status3 = vector_set.resizeVectors( nroots - 2 + 5 ); // add one more to cause append to fail
  if ( re_status3 != 0 ) {
    cout << "A linear dependancy was found at i " << re_status3 << endl;
    cout << "This currently is a hard failure: Aborting " << endl;
    GA::Terminate();
  }
// End FAILURE
  ap_status = vector_set.appendVectors( g_append);
 
   if ( nbasis < 100 ) {
    g_vector->print();
  }



  GA::sync();
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX finished XXXXXXXXXXXXXXXXXXX " << endl;
  
  GA::Terminate();
  
}



