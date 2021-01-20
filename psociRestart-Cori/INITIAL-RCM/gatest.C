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
 */
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>


#include <cmath>
#include <ga++.h>
//#include <mpi.h>

#define GA_DATA_TYPE C_DBL

#define GCOUT if (g_rank==0) cout

#include "PsociTimer.hpp"

extern "C" {
#include "armci.h"
#include "message.h"
#include "mpi.h"
#ifdef USEMPI
#include "ga-mpi.h"
#endif
}

//End ga++


// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//
//  Start driver program 
//
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


const unsigned long heap = 4L * 1024L * 1024L, stack = 9L * 1024L * 1024L ;
const unsigned long maxGAMemPerCore = 1.8L * 1024L * 1024L * 1024L;


using namespace std;

/* Basic functions */
void generateStaggeredChunkList(vector<int> & my_chunklist, int nchunk )
{
     my_chunklist.clear();

     int index = GA::nodeid() % nchunk;
#ifndef STAGGERCHUNKLIST
     index = 0; // override and make the list not-staggered; we then use brdcsts in the calling app
#endif
     for (int i=0; i< nchunk; ++i ) {
         my_chunklist.push_back( index );
         ++index;
         if ( index >= nchunk ) index = 0;
     }
}

int main(int argc, char **argv)
{
  // Set up global MPI fabric
  
  int g_size;
  int g_rank;
  

  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, maxGAMemPerCore );
  
  g_size = GA::nodes();
  g_rank = GA::nodeid();
  if ( GA::usesMA() ) cout << "GA memory is coming from MA " << endl;
  if ( GA::usesFAPI() ) cout << "GA is using Fortran indexing " << endl;

#ifdef USEMPI
  MPI_Comm mpi_comm; // for comparing GAdgop to MPI dgop
  mpi_comm = GA_MPI_Comm();
#endif
  
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  GCOUT << "RENCI Parallel GA test program " << endl;
  GCOUT << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
  GCOUT << endl;
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  
  GCOUT << endl;
  GCOUT << " -------------- setup ---------------------" << endl << endl;
  GCOUT << endl;

// Start the calculation
  
  GCOUT << "GA usesMA ? " << GA::usesMA() << endl;
  GCOUT << "GA memory limited ? " << GA::memoryLimited() << endl;
  GCOUT << "GA inquireMemory is " << GA::inquireMemory() << endl;
  GCOUT << "May not be GA array limit: GA memory avail " << GA::memoryAvailable() << endl;
  
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  GCOUT << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  GCOUT << "Parallel GA Run with " << g_size << " Total processors " << endl;
  
#ifdef USEMPI
  GCOUT << "REPLACE GA GOP and BRDCST with direct MPI calls " << endl;
#endif

// Create a test data object of size maxsefs.
 
 const int maxsefs = 35000000;
// const int maxsefs = 1000000;
 const int numrepetitions = 20;
// const int nchunk=1;
 const int nchunk=g_size;
//const int nchunk=3;


 double time, temp, getTime=0.0, totalTime=0.0, brdTime=0.0;

 int DATA_TYPE = C_DBL; // 
 int DATA_NDIM = 2;
     
 int vdim=2;
 int dims[2],chunk[2];
 dims[0] = maxsefs;
 dims[1] = 1;

 // Want distribution to be close to that of cimat and friends
 chunk[0] = -1;
 chunk[1] = 1; // horizontal strips

 time = psociTime();

 char * vector_title = "vector";
 GA::GlobalArray * g_v;
 g_v  = GA::createGA( C_DBL, vdim, dims, (char *) vector_title, chunk);
 g_v->zero();
 int dconst = 2.0;
 g_v->fill( &dconst );

 if (GA::nodeid() == 0 ) cout << "Print out GA Hamiltonian objects distribution " << endl;
 g_v->printDistribution();

// XXXXX start process 

  int ch_end=0;
  int chunk_width;
  (maxsefs % nchunk==0)? chunk_width = maxsefs / nchunk: chunk_width = 1 + (maxsefs / nchunk);
  GCOUT << "exp VECTOR get CHUNK is " << nchunk << " " << chunk_width <<endl;

  vector<int> my_chunklist;
  generateStaggeredChunkList( my_chunklist, nchunk);

//  vector<double> testVector( chunk_width , 0.0 );
  long totalbytes=0;

  int chunk_lo[2];
  int chunk_hi[2];
  vector<int>::iterator it;
  GCOUT << " start repetitions for GOP-base approach " << endl;

//    g_v->get( chunk_lo, chunk_hi, &testVector[0], &n );

    int n = 1;
    double dOne=1.0;

// alternative 1

#ifdef STAGGERCHUNKLIST
  GCOUT << "All get using staggered lists " << endl;
#endif

    brdTime = 0.0;
    getTime = 0.0;
    totalbytes = 0;

  vector<double> testVector( chunk_width , 0.0 );
  for(int irep=0; irep< numrepetitions; ++irep) {
  char op='+';
  for(it=my_chunklist.begin(); it!=my_chunklist.end(); ++it) { // this potentially staggers accross cores 

    int ch_start = chunk_width*(*it);
    int ch_end = min( ch_start + chunk_width, maxsefs );

    chunk_lo[0] = ch_start;
    chunk_hi[0] = ch_end-1;
    chunk_lo[1] = 0;
    chunk_hi[1] = 0;

    temp = psociTime();
    g_v->get( chunk_lo, chunk_hi, &testVector[0], &n );
    totalbytes += testVector.size()*sizeof(double);
    getTime += psociTime() - temp;
    }
   }

    GA::sync();
    totalTime = psociTime() - time;
    GCOUT << endl << "ALL get -----------------------------------------------------" << endl;
    GCOUT << "Total time + sync " << totalTime << endl;
    GCOUT << "Total GET time is " << getTime << endl;
    GCOUT << "numreps is " << numrepetitions<<endl;
    GCOUT << "total Bytes sent is " << totalbytes << endl;
    GCOUT << "GET effective MB/sec " << (double)totalbytes / (1048577.0 * (getTime) ) << endl;
    GCOUT << endl << endl;

    brdTime = 0.0;
    getTime = 0.0;
    totalbytes = 0;

// end alternative 1 


// alternative two
  for(int irep=0; irep< numrepetitions; ++irep) {
  char op='+';
  for(it=my_chunklist.begin(); it!=my_chunklist.end(); ++it) { // this potentially staggers accross cores 

    int ch_start = chunk_width*(*it);
    int ch_end = min( ch_start + chunk_width, maxsefs );

    chunk_lo[0] = ch_start;
    chunk_hi[0] = ch_end-1;
    chunk_lo[1] = 0;
    chunk_hi[1] = 0;

    vector<double> testVector( chunk_width , 0.0 );
    if ( GA::nodeid() == 0 ) g_v->get( chunk_lo, chunk_hi, &testVector[0], &n );
    getTime += psociTime() - temp;
    temp = psociTime();
    GA::brdcst( &testVector[0], testVector.size()*sizeof(double), 0);
    totalbytes += testVector.size()*sizeof(double);
    brdTime += psociTime() - temp;

    }
   }

    GA::sync();
    totalTime = psociTime() - time;
    GCOUT << endl << "RANK get + BRDCST-----------------------------------------------------" << endl;
    GCOUT << "Total time + sync " << totalTime << " Sum of collective is " << brdTime << endl;
    GCOUT << "Total GET time is " << getTime << endl;
    GCOUT << "numreps is " << numrepetitions<<" "<<"time/rep is " << brdTime/numrepetitions << endl;
    GCOUT << "total Bytes sent is " << totalbytes << endl;
    GCOUT << "BRD effective MB/sec " << (double)totalbytes / (1048577.0 * brdTime) << endl;
    GCOUT << "GET+BRD effective MB/sec " << (double)totalbytes / (1048577.0 * (brdTime+getTime) ) << endl;
    GCOUT << endl << endl;

    brdTime = 0.0;
    getTime = 0.0;
    totalbytes = 0;

//end alternative 2

// alternative-3

  for(int irep=0; irep< numrepetitions; ++irep) {
  char op='+';
  for(it=my_chunklist.begin(); it!=my_chunklist.end(); ++it) { // this potentially staggers accross cores 

    int ch_start = chunk_width*(*it);
    int ch_end = min( ch_start + chunk_width, maxsefs );

    chunk_lo[0] = ch_start;
    chunk_hi[0] = ch_end-1;
    chunk_lo[1] = 0;
    chunk_hi[1] = 0;

  int v_lo[2], v_hi[2];

  g_v->distribution( GA::nodeid(), v_lo, v_hi );
  //cout << "width and distribution is " << chunk_width << " " << v_lo[0]<<" "<<v_hi[0] << endl;

  int newlow[2];
  int newhi[2];
  int my_vector_low = v_lo[0];
  int my_vector_hi = v_hi[0]; // this does not need subtraction by 1

  //cout << "size of myvector " << my_vector_low<<" "<<my_vector_hi <<endl;
  int ivec = 0;
  newlow[1] = ivec;
  newhi[1] = ivec;

  int ichunk = (*it);

  int index = 0;
  vector<int> l_list;
  int indexstart = chunk_width * ichunk;

  for(int i=ch_start; i< ch_end; ++i ) { // loop entire actual chunk range
     if ( i >= my_vector_low && i <= my_vector_hi ) l_list.push_back( i ); // do I have any of it.
  }
//  vector<double> testVector( chunk_width , 0.0 );

  temp = psociTime();
  if ( l_list.size() > 0 ) {
     newlow[0] = l_list[0];
     newhi[0] = l_list[ l_list.size() - 1 ];
     index = newlow[0] - ch_start;
     g_v->get( &newlow[0], &newhi[0], &testVector[index], &n );
  }
  getTime += psociTime() - temp;
  temp = psociTime();

#ifdef USEMPI
//  MPI_Allreduce( &testVector[0], &testVector[0], testVector.size(), MPI_DOUBLE, MPI_SUM, mpi_comm);
  MPI_Allreduce( MPI_IN_PLACE, &testVector[0], testVector.size(), MPI_DOUBLE, MPI_SUM, mpi_comm);

#else
  GA::dgop( &testVector[0], testVector.size(), &op);
#endif

  totalbytes += testVector.size()*sizeof(double);
  brdTime += psociTime() - temp;
   } // my chunk list
   } // repetition list

// XXXXX clean up

 if ( GA::nodeid() == 0 ) GA::printStats();

 GA::sync();
 totalTime = psociTime() - time;
 GCOUT << endl << "----------------------------------------------------------" << endl;
 GCOUT << "Total time + sync " << totalTime << " Sum of collective is " << brdTime << endl;
 GCOUT << "Total GET time is " << getTime << endl;
 GCOUT << "numreps is " << numrepetitions<<" "<<"time/rep is " << brdTime/numrepetitions << endl;
 GCOUT << "total Bytes sent is " << totalbytes << endl;
 GCOUT << "GOP effective MB/sec " << (double)totalbytes / (1048577.0 * brdTime) << endl;
 GCOUT << "GET+GOP effective MB/sec " << (double)totalbytes / (1048577.0 * (brdTime+getTime) ) << endl;
 GCOUT << endl << endl;


// Compare to a BRDCST based approach

/*
  GCOUT << " start repetitions for BRDCST-base approach " << endl;

  brdTime=0.0;
  getTime=0.0;
  totalbytes=0;
  totalTime=0.0;
  
  for(int irep=0; irep< numrepetitions; ++irep) {
  char op='+'; 
  for(it=my_chunklist.begin(); it!=my_chunklist.end(); ++it) { // this potentially staggers accross cores 
        
    int ch_start = chunk_width*(*it);
    int ch_end = min( ch_start + chunk_width, maxsefs );

    chunk_lo[0] = ch_start;
    chunk_hi[0] = ch_end-1;
    chunk_lo[1] = 0;
    chunk_hi[1] = 0;

//    g_v->get( chunk_lo, chunk_hi, &testVector[0], &n );

    int n = 1;
    double dOne=1.0;

//    vector<double> testVector( chunk_width , 0.0 );

  temp = psociTime();
  int v_lo[2], v_hi[2];

  g_v->distribution( GA::nodeid(), v_lo, v_hi );
  //cout << "next width and distribution is " << chunk_width << " " << v_lo[0]<<" "<<v_hi[0] << endl;

  int newlow[2];
  int newhi[2];
  int my_vector_low = v_lo[0];
  int my_vector_hi = v_hi[0]; // this does not need subtraction by 1

  int ivec = 0;
  newlow[1] = ivec;
  newhi[1] = ivec;

  int ichunk = (*it);

  int index = 0;
  vector<int> l_list;
  int indexstart = chunk_width * ichunk;

  for(int i=ch_start; i< ch_end; ++i ) { // loop entire actual chunk range
     if ( i >= my_vector_low && i <= my_vector_hi ) l_list.push_back( i ); // do I have any of it.
  }
    vector<double> testVector( chunk_width , 0.0 );

  temp = psociTime();
  if ( l_list.size() > 0 ) {
     newlow[0] = l_list[0];
     newhi[0] = l_list[ l_list.size() - 1 ];
     index = newlow[0] - ch_start;
     g_v->get( &newlow[0], &newhi[0], &testVector[index], &n );
  }
  getTime += psociTime() - temp;

  temp = psociTime();
#ifdef USEMPI
  int num = MPI_Bcast( &testVector[0], testVector.size(), MPI_DOUBLE, 0, mpi_comm );
#else
  GA::brdcst( &testVector[0], testVector.size()*sizeof(double), 0);
#endif

  totalbytes += testVector.size()*sizeof(double);
  brdTime += psociTime() - temp;

   } // my chunk list
   } // repetition list

// XXXXX clean up

 if ( GA::nodeid() == 0 ) GA::printStats();
 GA::sync();
 totalTime = psociTime() - time;
 GCOUT << endl << "----------------------------------------------------------" << endl;
 GCOUT << "Total time + sync " << totalTime << " Sum of collective is " << brdTime << endl;
 GCOUT << "Total GET time is " << getTime << endl;
 GCOUT << "numreps is " << numrepetitions<<" "<<"time/rep is " << brdTime/numrepetitions << endl;
 GCOUT << "total Bytes sent is " << totalbytes << endl;
 GCOUT << "BRDCST effective MB/sec " << (double)totalbytes / (1048577.0 * brdTime) << endl;
 GCOUT << "GET+BRDCST effective MB/sec " << (double)totalbytes / (1048577.0 * (brdTime+getTime) ) << endl;
 GCOUT << endl << endl;

*/


 }







