/**************************************************************************************
 * Copyright (c) 2013 RENCI.
 * All rights reserved. This program and the accompanying materials
 * MAY BE available under the terms of the RENCI Open Source License
 * UNC at Chapel Hill which accompanies this distribution, and is available at
 * http://www.renci.org/resources/open-source-software-license
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

//left overs from using MPI _Bcast instead of GA brdcst.
extern "C"
{
#include <mpi.h>
#include <ga.h>
}

#define GCOUT if (g_rank==0) cout

using namespace std;

  const int numSoughtRoots=40;
  const int maxsubspace=200;
  const int maxsef=5000;


void fetchAndReplicateVectorChunk( int ichunk, int ivec, int width, int * v_lo, int * v_hi, vector<double> & chunk, GA::GlobalArray * g_v )
{
  int l_maxsefs = maxsef;

  int chunk_lo[2];
  int chunk_hi[2];

  chunk_lo[1] = ivec;
  chunk_hi[1] = ivec;

  int ch_start = width*(ichunk);
  int ch_end = min( ch_start + width, l_maxsefs );
  int n = 1;
  char op='+';

  chunk_lo[0] = ch_start;
  chunk_hi[0] = ch_end-1; // this does need subtraction by 1

  if ( GA::nodeid() == 0 ) cout << GA::nodeid() <<" perform the GET " << n<<" "<<chunk_lo[0]<<" "<<chunk_lo[1]<<" "<<chunk_hi[0]<<" "<<chunk_hi[1]<<" "<<chunk.size() << endl;
  if ( GA::nodeid() == 0 ) g_v->get( chunk_lo, chunk_hi, &chunk[0], &n );

  GA::sync();
  if ( GA::nodeid() == 0 ) cout << GA::nodeid() <<" perform the BRDCST: extra gsyn" << endl;

  GA::brdcst( &chunk[0], chunk.size()*sizeof(double), 0);
}

void fetchAndReplicateVectorChunk2( int ichunk, int ivec, int width, int * v_lo, int * v_hi,double *chunk, GA::GlobalArray * g_v )
{
  int l_maxsefs = maxsef;

  int chunk_lo[2];
  int chunk_hi[2];

  chunk_lo[1] = ivec;
  chunk_hi[1] = ivec;

  int ch_start = width*(ichunk);
  int ch_end = min( ch_start + width, l_maxsefs );
  int n = 1;
  char op='+';

  chunk_lo[0] = ch_start;
  chunk_hi[0] = ch_end-1; // this does need subtraction by 1

  if ( GA::nodeid() == 0 ) cout << GA::nodeid() <<" perform the GET " << n<<" "<<chunk_lo[0]<<" "<<chunk_lo[1]<<" "<<chunk_hi[0]<<" "<<chunk_hi[1]<<" "<<endl;
  if ( GA::nodeid() == 0 ) g_v->get( chunk_lo, chunk_hi, chunk, &n );

  if ( GA::nodeid() == 0 ) cout << GA::nodeid() <<" perform the BRDCST" << endl;

  GA::brdcst( &chunk[0], width*sizeof(double), 0);
}
// Start the work

int main(int argc, char **argv)
{

  const unsigned long heap = 128L * 1024L * 1024L, stack = 128L * 1024L * 1024L ;
  const unsigned long maxGAMemPerCore = 1.0L * 1024L * 1024L * 1024L;
//  int memoryStackWords = stack/8;
//  int memoryHeapWords  = heap/8;
//  int maxGAMemPerCoreWords = maxGAMemPerCore/8;

  int memoryStackWords = stack;
  int memoryHeapWords  = heap;
  int maxGAMemPerCoreWords = maxGAMemPerCore;


  MPI_Init(&argc, &argv);
  GA_Initialize();
  if ( GA_Uses_ma() ) {
    MA_init(C_DBL, memoryStackWords, memoryHeapWords+maxGAMemPerCoreWords );
  } else {
    MA_init(C_DBL, memoryStackWords, memoryHeapWords);
  }

  int g_rank = GA::nodeid();

  GCOUT << "Test of reading GA 10 times " << endl;

  char * vector_title = "vector";
  int vdim=2;
  int dims[2],chunk[2];
  dims[0] = maxsef;
  dims[1] = numSoughtRoots+maxsubspace;
  chunk[0] = -1;
  chunk[1] = numSoughtRoots+maxsubspace; // horizontal strips

  GA::GlobalArray * g_v;
  g_v  = GA::createGA( C_DBL, vdim, dims, (char *) vector_title, chunk);
  g_v->zero();

  GCOUT << " VECTOR distribution " << endl;
  g_v->printDistribution();

  int vtype,vidim;
  int vlocal_dims[2];
  int v_lo[2], v_hi[2];
  g_v->inquire( &vtype, &vidim , vlocal_dims );
  g_v->distribution( GA::nodeid(), v_lo, v_hi ); // not used here but mimics the real program

	// Start test

   int chunk_width = maxsef;
   int it=0;
   int chunk_lo[2], chunk_hi[2];

   int ch_start = 0;
   int ch_end = min( ch_start + chunk_width, maxsef );
   chunk_lo[0] = ch_start;
   chunk_hi[0] = ch_end-1;


        // Step One: Grab the current GA-vector chunk for processing

   int ivec=0;
   for(int kk=0; kk< 10; ++kk) {
    if ( GA::nodeid() == 0 ) cout << " fake read " << kk<< endl;
       // vector<double> testVector( chunk_width );  
       // fetchAndReplicateVectorChunk(it, ivec, chunk_width, v_lo, v_hi, testVector, g_v );

        double pTestVector[chunk_width];
        fetchAndReplicateVectorChunk2(it, ivec, chunk_width, v_lo, v_hi, pTestVector, g_v );

    }
    GA::sync();
    GCOUT << " done with fake read " << endl;
}
