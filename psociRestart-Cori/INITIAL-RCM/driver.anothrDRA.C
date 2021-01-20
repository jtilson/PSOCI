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

// Now start testing for GA
#include <ga++.h>
#include <GAServices.h>

#include <dra.h>
#define GA_DATA_TYPE MT_DBL

//End ga++

#include "PsociTimer.hpp"
#include "PsociVector.hpp"
#include "PsociGArestart.hpp"
#include "PsociDRAservices.hpp"
#include "PsociDRArestart.hpp"

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
  
  unsigned int heap=9000000, stack=9000000;
  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, 0);
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
  
  int64_t nroots = 10;
  //int64_t nbasis = 10000;
  int64_t nbasis = 20;
  
  double shift = 0.0;
  
  double fudge=2.0;

  for( int64_t i=0; i< nroots; ++i) {
    shift += 1.0;
    double value;
    vector<double> root;
    for ( int64_t nb=0; nb < nbasis; ++nb )
      {
        ++fudge;
	value = base*shift+fudge;;
	root.push_back( value/10000.0 );
      }
    testData.push_back( root );
  }
  
  pair<int64_t,int64_t> info;
  info.first = nroots;
  info.second = nbasis; 
  
  cout << "Size of testData is " << testData.size() << endl;
  cout << "Start PsociVector methods " << endl;
  
  // Directly  Test the formatted write/read routines of PSOCI
  
  PsociVector vec1 = PsociVector();
  cout << "Vector library version is " << vec1.getVersion() << endl;
  vec1.setMode( VECTOR_WRITE );
  vec1.OpenFile();
  cout << "Fetch filename to write is " << vec1.FetchFilename() << endl;
  cout << "MODE is " << vec1.fetchMode() << endl;
  vec1.printWriteType();
  cout << " Time for Print output of vector is " << vec1.PrintBurst( testData, info ) << endl;
  cout << " Time for WRITE output of vector is " << vec1.WriteBurst( testData, info ) << endl;
  vec1.CloseFile();
  cout << endl << endl;
  
  vector<vector<double> > testData2; 
  pair<int64_t,int64_t> info2;
  
  PsociVector vec2 = PsociVector();
  vec2.setMode( VECTOR_READ );
  vec2.OpenFile();
  cout << endl << "(ascii) Fetch filename to read is  " << vec2.FetchFilename() << endl;
  vec2.printWriteType();
  cout << "Time for ASCII INPUT of vector is " << vec2.ReadBurst( testData2, info2 ) << endl;
  vec2.CloseFile();
  cout << endl << endl;
  cout << "Dump ASCII vectors to stdout " << endl;
  vector<double> norm1;
  if ( !vec2.checkNormalization( testData2 , norm1) ) {
     cout << "Warning: at least one ascii unnormalized vector " << endl;
  }

// Create a 2dim Global array of testData3: (nbfs x nroots)

  cout << endl << endl << "XXXXXXXXXXXXXXXXXXX test the GA restart routines " << endl;
  
  int VECTOR_TYPE = MT_DBL;
  int VECTOR_NDIM = 2;
  int dims[2];
  dims[0] = info2.first;
  dims[1] = info2.second; //Careful info is <int64_t,int64_t>
  cout << "Create a GA of dims " << dims[0] << " by " << dims[1] << endl;

  const int DISTRIBUTE_EVENLY = -1;

  int chunk[2];
  chunk[0] = DISTRIBUTE_EVENLY; //We want distribution across chunk0
  chunk[1] = DISTRIBUTE_EVENLY;
  
  string arrayName = "vector_array";
  //GA::GlobalArray *g_vector = GA::SERVICES.createGA( VECTOR_TYPE, VECTOR_NDIM, dims, (char *)"test", chunk); 
  GA::GlobalArray *g_vector = GA::SERVICES.createGA( VECTOR_TYPE, VECTOR_NDIM, dims, (char *)"test", NULL); 
  g_vector->zero();
  
//Should be 10X20 at this point

  int nroot=info2.first;
  if ( g_rank == 0 ) {
    vector<vector<double> >::iterator git;
/*
    int lo[2];
    lo[0] = 0;
    lo[1] = 0;
    int hi[2];
    hi[0] = info3.first-1;
    hi[1] = 0; 
*/

    int lo[2];
    lo[0] = 0;
    lo[1] = 0;
    int hi[2];
    hi[0] = 0;
    hi[1] = info2.second-1;

    int n = 1;
    
    for(git = testData2.begin(); git != testData2.end(); ++git)
      { 
	cout << "XXXXX Put next vector " << endl;
	cout << "ROW lo vals " << lo[0] << " " << hi[0] << "COL lo vals " << lo[1] << " " << hi[1] << (*git)[0] << " " << (*git)[1] << " " << (*git)[19] << endl;
	g_vector->put(lo, hi, (void *)&(*git)[0], &n); 
	++lo[0];
	++hi[0];
      }
  }
  GA::sync();
  g_vector->printDistribution();
  
  if ( nbasis <= 1000 ) { 
    cout << "print array " << endl;
    g_vector->print();
  }
  


// Test dumping GLobal Array test data to PsociVector disk.

  PsociGArestart resVec( g_rank, g_vector );

  string filename2 = "testGAvector.txt";
  if ( g_rank == 0 ) {
    resVec.openVectorWrite( g_rank, filename2 );
    double time;
    resVec.dumpVector( g_rank, g_vector, time );
    cout << "Time to dump ga to disk (s)" << time << endl;
    resVec.closeVectorFile();
  }
  GA::sync();

// Test reading GA from PsociVector info g_read_vector

  GA::GlobalArray *g_read_vector = GA::SERVICES.createGA( *g_vector );
  PsociGArestart resVec2( g_rank, g_read_vector );
  string filename3 = "testGAvector.txt";
    
  if (g_rank == 0 ) {
    resVec2.openVectorRead( g_rank, filename3 );
    double timeR;
    resVec2.readVector( g_rank, g_read_vector, timeR );
    resVec2.closeVectorFile();
    cout << "timeR value is " << timeR << endl;
    cout << "print READ array " << endl;
  }
  GA::sync();





  GA::summarize( 1 ); 
  GA::printStats();
  char * checker = "testChecker";
  g_vector->checkHandle( checker );  //If check fails this method fills the job
  
  GA::sync();
  cout << "Begin DRA work: I am " << g_rank << endl;

  PsociDRAservices DRAservice;
  DRAservice.initDRA( g_rank );
  DRAservice.reportDRAIOservers();
  DRAservice.reportDRAParams();


  string draarrayname = "TESTarrayname";
  string drafilename = "filenameDRA.txt";

  cout << "Got the DRA initialized " << g_rank << endl;
  GA::sync();

  PsociDRArestart dravector( g_rank, draarrayname, drafilename );
//  PsociDRArestart dravector;

  dravector.reportName();
  dravector.reportFilename();
  dravector.setType( MT_DBL );
  dravector.printType();
  dravector.setMode( DRA_RW );

  GA::sync();
  cout << "On to the creation of ddims " << g_rank << endl;

  
  long ndims = 2;
  long maxsparse = 50;
  pair<dra_size_t,dra_size_t> ddims;
  ddims.first=nroots;
  ddims.second=nbasis;

  pair<dra_size_t,dra_size_t> rdims;
  rdims.first=-1;
  rdims.second=-1;

  dravector.setDimensions( ndims, ddims, rdims ); 

  cout << "Manual creation of DRA at driver " << g_rank << endl;

  cout << "Onto the CREATE Ranks are " << g_rank << endl;
  if ( dravector.createDRA() != 0 ) 
  {
     cout << "failed to create the DRA: Aborting app " << endl;
     GA::Terminate();
  }

  dravector.printInternals();

//For now dump the entire GA into the first nroots oclumns of DRA

  vector< pair<int,int> > ga_patch;
  ga_patch.push_back( pair<int,int>( 1,nroots));
  ga_patch.push_back( pair<int,int>( 1,nbasis ));

  vector< pair<dra_size_t,dra_size_t> > dra_patch;
  dra_patch.push_back( pair<dra_size_t,dra_size_t>( 1,nroots ));
  dra_patch.push_back( pair<dra_size_t,dra_size_t>( 1,nbasis ));

  if ( dravector.dumpGASectionDRA( ga_patch, dra_patch, g_vector) != 0 ) {
     cout << "Dump failed " << endl;
  }

//Now dump same patch into the next set or root area
  
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX finished XXXXXXXXXXXXXXXXXXX " << endl;
  GA::sync();
  //DRAservice.exitDRA();
  GA::Terminate();


  
  

}



