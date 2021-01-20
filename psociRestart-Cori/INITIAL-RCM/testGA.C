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

// Start the work

int main(int argc, char **argv)
{

  const unsigned long heap = 128L * 1024L * 1024L, stack = 128L * 1024L * 1024L ;
  const unsigned long maxGAMemPerCore = 1.0L * 1024L * 1024L * 1024L;
  int memoryStackWords = stack/8;
  int memoryHeapWords  = heap/8;
  int maxGAMemPerCoreWords = maxGAMemPerCore/8;

  MPI_Init(&argc, &argv);
  GA_Initialize();
  if ( GA_Uses_ma() ) {
    MA_init(C_DBL, memoryStackWords, memoryHeapWords+maxGAMemPerCoreWords );
  } else {
    MA_init(C_DBL, memoryStackWords, memoryHeapWords);
  }

  int g_rank = GA::nodeid();

  GCOUT << "YAFT of GA: Test of put versus acc " << endl;

  char * vector_title = "vector";
  int vdim=2;
  int dims[2],chunk[2];
  dims[0] = 100;
  dims[1] = 10;
  chunk[0] = -1;
  chunk[1] = 100; // horizontal strips

  GA::GlobalArray * g_v;
  g_v  = GA::createGA( C_DBL, vdim, dims, (char *) vector_title, chunk);

	// Start test
        // generate some data and put or acc to g_v

 g_v->printDistribution();

// Now try again using an irreg distribution 


  int block[2];
  int g_size=GA::nnodes();

  block[0]=g_size;
  block[1]=1;

// int map={0,2,6,0,5};

  int map[5];
  map[0]=0;
  map[1]=2;
  map[2]=6;
  map[3]=0;
//  map[4]=5;


  GA::GlobalArray * g_v2;
  g_v2  = GA::createGA( C_DBL, vdim, dims, (char *) vector_title, block, map);

 g_v2->printDistribution();











/*
  int lo[2], hi[2];

   lo[0] = 0;
   hi[0] = 99;
   lo[1] = 1;
   hi[1] = 1;
   int n = 1;

// simple ACC yields empty matrix

   vector<double> data2(100);
   for( int kk=0; kk< 100; ++kk) {
      data2[kk] = 2.0;
   }

   double dAlpha=1.0;
   if ( GA::nodeid() == 0 ) {
    g_v->acc( lo, hi, &data2[0], &n, &dAlpha );
   }

   if ( GA::nodeid() == 0 ) cout << "ACC test" << endl;
   g_v->print();
*/

}

