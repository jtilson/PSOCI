/**************************************************************************************
 * Copyright (c) 2010,2011 RENCI.
 * All rights reserved. This program and the accompanying materials
 * MAY BE available under the terms of the RENCI Open Source License
 * UNC at Chapel Hill which accompanies this distribution, and is available at
 * http://www.renci.org/resources/open-source-software-license

 * New implementation of PSOCI:

 Classes: 

 Description: 

  Test higher order aspects of the distributed iterative diagionalization routines.
 
 History:

**************************************************************************************/
/**
 *   @file 
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
//deprecated #include <GAServices.h>

#include <dra.h>
#define GA_DATA_TYPE C_DBL

//End ga++

#include "PsociTimer.hpp"
#include "PsociVector.hpp"
#include "PsociGArestart.hpp"
#include "PsociDRAservices.hpp"
#include "PsociDRArestart.hpp"
#include "PsociGAbasis.hpp"
#include "PsociGAhamiltonian.hpp"
#include "PsociGAsubspace.hpp"

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
  
  // Create a dummy array of eigenvectors for storage
  
  double base = 3.12345678;
  vector<vector<double> > testData;
  
  // Try a realistic version 50 roots 1M basis
  
  int64_t nroots = 3;
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
        if ( nb == index-1 ) {
        value = base*shift+fudge;
        cout << "2nd value " << value << endl;
        root[nb] = 2.0;
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
  
// Create a 2dim Global array of testData3: (nbfs x nroots)

  cout << endl << endl << "XXXXXXXXXXXXXXXXXXX test the GA restart routines " << endl;
  
  int VECTOR_TYPE = C_DBL;
  int VECTOR_NDIM = 2;
  int dims[2];
  cout << "info.first is " << info.first << " .second is " << info.second << endl;
  dims[0] = info.first;
  dims[1] = info.second+10; //Careful info is <int64_t,int64_t>
  cout << "Create a GA of dims " << dims[0] << " by " << dims[1] << endl;

  const int DISTRIBUTE_EVENLY = -1;

  int chunk[2];
  chunk[0] = DISTRIBUTE_EVENLY; //We want distribution across chunk0
  chunk[1] = DISTRIBUTE_EVENLY;
  
  int adims[2];
  cout << "info.first is " << info.first << " .second is " << info.second << endl;
  adims[0] = info.first;
  adims[1] = 2; //Careful info is <int64_t,int64_t>

  string arrayName = "vector_array";
  GA::GlobalArray *g_vector = GA::createGA( VECTOR_TYPE, VECTOR_NDIM, dims, (char *)"test", NULL); 
  GA::GlobalArray *g_Hvector = GA::createGA( VECTOR_TYPE, VECTOR_NDIM, dims, (char *)"testHv", NULL); //Could use duplicate....
  g_vector->zero();
  g_Hvector->zero(); //Not really needed actually

  int type_test;
  int ndims_test;
  int dims_test[2];
  
//Should be 20x10 at this point

  int nroot=info.first;
    vector<vector<double> >::iterator git;
    int lo[2];
    lo[0] = 0;
    lo[1] = 0;
    int hi[2];
    hi[1] = 0;
    hi[0] = info.first-1;
    int n = 1;
    
  g_vector->printDistribution();
  g_Hvector->printDistribution();

    lo[0] = 0;
    lo[1] = 0;
    hi[1] = 0;
    hi[0] = info.first-1;

//Fill in Vector array for testing

    if ( g_rank == 0 ) {
      for( git = testData.begin(); git != testData.end(); ++git)
        {
          cout << "XXXXX Put next vector " << endl;
          cout << "ROW lo vals " << lo[0] << " " << hi[0] << "COL lo vals " << lo[1] << " " << hi[1] << (*git)[0] << " " << (*git)[1] << " " << (*git)[5] << endl;
          g_vector->put(lo, hi, (void *)&(*git)[0], &n);
          ++lo[1];
          ++hi[1];
        }
    }
  GA::sync();

  if ( nbasis <= 100 ) {
    cout << "print array " << endl;
    g_vector->print();
  }

  cout << "Now fill fake Hv array " << endl;
  PsociGAhamiltonian g_h( g_Hvector );

  cout << "test build " << endl;
  g_h.testBuild();

  cout << "New H method Rows = " << g_h.CurrentNumRows() << endl;
  cout << "New H method Cols = " << g_h.CurrentNumCols() << endl;
  
  if ( nbasis <= 100 ) { 
    cout << "print Hv array " << endl;
    (g_h.fetchHandle())->print();
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
  GA::sync();
  
//Now orthoginalize the vector array only 

//Add GArestart to GAbasis with a .init and .read and .write method
  PsociGAbasis vector_set( nroots, g_vector );
  vector<double> vecnorm;
  
  cout << "EST fetch handle " << endl;
  vector_set.fetchHandle()->print();

  vector_set.getAllNorms( vecnorm );
  GA::sync();
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
  

// Prepare to compute eigenvectors
  int64_t nsub = vector_set.CurrentNumSoughtRoots();
  cout << " Size of current subspace is " << nsub << endl;

//Add vector space to the current set of vectors

 //int re_status = vector_set.resizeVectors( nroots + 1 );


// Construct subspace objects.

   PsociGAsubspace subspace( &vector_set, &g_h );
   cout << " maxsubspace is " << subspace.fetchMaxSubspace() << endl;
   cout << " subspace CI rows is " << subspace.numCIBasis() << endl;
   cout << " subspace Vec rows is " << subspace.numVectorBasis() << endl;
   cout << " subspace Vec cols is " << subspace.numRoots() << endl;
   cout << " subspace MaxRoots is " << subspace.fetchMaxSubspace() << endl;

   subspace.generateSubspace();
   subspace.diagSubspace();
   subspace.computeRnorms();
   subspace.destroySubspace();

  
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX finished XXXXXXXXXXXXXXXXXXX " << endl;
  GA::sync();
  GA::Terminate();
  
}



