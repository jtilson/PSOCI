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

int main(int argc, char **argv)
{
  
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

  int VECTOR_TYPE = MT_DBL;
  int VECTOR_NDIM = 2;
  int dims[2];
  dims[0] = 10;
  dims[1] = 20;
  cout << "Create a GA of dims " << dims[0] << " by " << dims[1] << endl;

  const int DISTRIBUTE_EVENLY = -1;

  int chunk[2];
  chunk[0] = DISTRIBUTE_EVENLY; //We want distribution across chunk0
  chunk[1] = DISTRIBUTE_EVENLY;
  
  string arrayName = "vector_array";
  GA::GlobalArray *g_vector = GA::SERVICES.createGA( VECTOR_TYPE, VECTOR_NDIM, dims, (char *)"test", NULL); 
  g_vector->zero();
  
  GA::summarize( 0 ); 
   
  char * checker = "check dump GA vector";
  g_vector->checkHandle( checker );  //If check fails this method fills the job
  
  cout << "Begin DRA work: I am " << g_rank << endl;

  PsociDRAservices DRAservice;
  DRAservice.initDRA( g_rank );
  DRAservice.reportDRAIOservers();
  DRAservice.reportDRAParams();

  cout << "Got the DRA initialized " << g_rank << endl;
  cout << "Manual creation of DRA at driver " << g_rank << endl;

  string draarrayname = "TESTarrayname";
  string drafilename = "filenameDRA.txt";
  PsociDRArestart dravector( g_rank, draarrayname, drafilename );
  dravector.reportName();
  dravector.reportFilename();
  dravector.setType( MT_DBL );
  dravector.printType();
  dravector.setMode( DRA_RW );

  long ndims = 2;
  long maxsparse = 50;
  pair<dra_size_t,dra_size_t> ddims;
  ddims.first=10;
  ddims.second=20;
    
  pair<dra_size_t,dra_size_t> rdims;
  rdims.first=-1;
  rdims.second=-1;
  dravector.setDimensions( ndims, ddims, rdims );

  cout << "Onto the CREATE Ranks are " << g_rank << endl;
  if ( dravector.createDRA() != 0 )
  {
     cout << "failed to create the DRA: Aborting app " << endl;
     GA::Terminate();
  }

  dravector.printInternals();

//For now dump the entire GA into the first nroots oclumns of DRA

  vector< pair<int,int> > ga_patch;
  ga_patch.push_back( pair<int,int>( 1,10));
  ga_patch.push_back( pair<int,int>( 1,20 ));

  vector< pair<dra_size_t,dra_size_t> > dra_patch;
  dra_patch.push_back( pair<dra_size_t,dra_size_t>( 1,10 ));
  dra_patch.push_back( pair<dra_size_t,dra_size_t>( 1,20 ));

  if ( dravector.dumpGASectionDRA( ga_patch, dra_patch, g_vector) != 0 ) {
     cout << "Dump failed " << endl;
  }


  GA::sync();  
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX finished XXXXXXXXXXXXXXXXXXX " << endl;
  GA::Terminate();

}



