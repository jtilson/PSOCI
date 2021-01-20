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
  
  int64_t nroots = 10;
  //int64_t nbasis = 1000000;
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
  info.first = nbasis;
  info.second = nroots;
  
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
#ifdef BINARY
  cout << " Time for BINARY WRITE output of vector is " << vec1.WriteBurstBinary( testData, info ) << endl;
#else
  cout << " Time for WRITE output of vector is " << vec1.WriteBurst( testData, info ) << endl;
#endif
  if ( nbasis < 1000 ) {
    vector<double> norm1;
    if ( !vec1.checkNormalization( testData , norm1) ) {
      cout << "Warning: at least one ascii unnormalized vector " << endl;
    }
    for (int i=0;i<nroots;++i) {
      cout << "ASCII " << norm1[i] << endl;
    }
  }
  
  vec1.CloseFile();
  cout << endl << endl;
  
  // Create a 2dim Global array of testData3: (nbfs x nroots)
  
  cout << endl << endl << "XXXXXXXXXXXXXXXXXXX test the GA restart routines " << endl;
  
  int VECTOR_TYPE = C_DBL;
  int VECTOR_NDIM = 2;
  int dims[2];
  dims[0] = info.first;
  dims[1] = info.second; //Careful info is <int64_t,int64_t>
  cout << "Create a GA of dims " << dims[0] << " by " << dims[1] << endl;
  
  const int DISTRIBUTE_EVENLY = -1;
  
  int chunk[2];
  chunk[0] = DISTRIBUTE_EVENLY; //We want distribution across chunk0
  chunk[1] = DISTRIBUTE_EVENLY;
  
  string arrayName = "vector_array";
  //GA::GlobalArray *g_vector = GA::createGA( VECTOR_TYPE, VECTOR_NDIM, dims, (char *)"test", chunk); 
  GA::GlobalArray *g_vector = GA::createGA( VECTOR_TYPE, VECTOR_NDIM, dims, (char *)"test", NULL); 
  g_vector->zero();
  
  int type_test;
  int ndims_test;
  int dims_test[2];
  
  
  //Should be 20x10 at this point
  
  //  int nbasis=info.first;
  
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
    hi[0] = info.first-1;
    hi[1] = 0;
    
    int n = 1;
    
    for(git = testData.begin(); git != testData.end(); ++git)
      { 
	cout << "XXXXX Put next vector " << endl;
	cout << "Basis lo vals " << lo[0] << " " << hi[0] << "Row lo vals " << lo[1] << " " << hi[1] << (*git)[0] << " " << (*git)[1] << " " << (*git)[19] << endl;
	g_vector->put(lo, hi, (void *)&(*git)[0], &n); 
	++lo[1];
	++hi[1];
      }
  }
  GA::sync();
  g_vector->printDistribution();
  
  if ( nbasis <= 1000 ) { 
    cout << "print array " << endl;
    g_vector->print();
  }
  
  // Test dumping GLobal Array test data to PsociVector disk.
  
  cout << "Do the next restart " << endl;
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
  
  int alo[2];
  int ahi[2];
  
  /* works fine and as expected */
  
  vector<double>  buf;
  buf.resize(nbasis);
  int ldv=1;
  alo[0]=0;
  ahi[0]=nbasis-1;
  alo[1]=0;
  ahi[1]=0;
  g_vector->get( alo, ahi, &buf[0], &ldv);
  
/*
  double sum=0;
  for (int i=0; i<nbasis; ++i) {
    cout << "gotten data are " << buf[i] << endl;
    sum += buf[i] * buf[i];
  }
  cout << " sum/norm are " << sum << " " << sqrt(sum) << endl;
*/
  
  
  /* this one also now works */
  
  alo[0]=0;
  ahi[0]=nbasis-1;
  alo[1]=0;
  ahi[1]=0;
  
  double ddot1;
  ddot1 = g_vector->ddotPatch('N',alo,ahi, g_vector,'N', alo,ahi );
  
  cout << "rank is " << g_rank << " Value fo Vecvtor=1 ddot N*N is " << ddot1 <<  endl;
  
  
  // Try out new VECTOR routines
  
  PsociGAbasis vector_set( g_vector );
  vector<double> vecnorm;
  
  cout << " on the the vecnorm " << endl;
  vector_set.getAllNorms( vecnorm );
  cout << "Size of vecnorm is " << vecnorm.size() << endl;
  
  for (int i=0; i<nroots; ++i) {
    cout << "VECNORM data are " << vecnorm[i] << endl;
  }
  
  
  // Test reading GA from PsociVector info g_read_vector
  
  GA::GlobalArray *g_read_vector = GA::createGA( *g_vector );
  
  PsociGArestart resVec2( g_rank, g_read_vector );
  string filename3 = "testGAvector.txt";
  
  //Now the openVectorRead is done by a single core BUT, we may need to 
  //have every core  now how many cores were read. 
  int numRootsRead;
  if (g_rank == 0 ) {
    resVec2.openVectorRead( g_rank, filename3 );
    double timeR;
    resVec2.readVector( g_rank, g_read_vector, timeR );
    resVec2.closeVectorFile();
    cout << "timeR value is " << timeR << endl;
    cout << "print READ array " << endl;
    numRootsRead = resVec2.getNumberVectors();
  }

  GA::brdcst( &numRootsRead, sizeof(numRootsRead) , 0 );
  cout << "I am rank " << g_rank << " and was told g_vector had this many roots " << numRootsRead << endl;
  GA::sync();
  
  //Dump array to DRA as a test 
  
  PsociDRAservices DRAservice;
  DRAservice.initDRA( g_rank );
  DRAservice.reportDRAIOservers();
  DRAservice.reportDRAParams();
  
  string draarrayname = "TESTarrayname";
  string drafilename = "filenameDRA.txt";
  
  cout << "Got the DRA initialized " << g_rank << endl;
  GA::sync();
  
  PsociDRArestart dravector( g_rank, draarrayname, drafilename );
  
  dravector.reportName();
  dravector.reportFilename();
  dravector.setType( C_DBL );
  dravector.printType();
  dravector.setMode( DRA_RW );
  
  long ndims = 2;
  long maxsparse = nroots;
  
  pair<dra_size_t,dra_size_t> ddims;
  ddims.first=nbasis;
  ddims.second=nroots;

  pair<dra_size_t,dra_size_t> rdims;
  rdims.first=-1;
  rdims.second=-1;
  dravector.setDimensions( ndims, ddims, rdims );

  if ( dravector.createDRA() != 0 )
  {
    cout << "failed to create the DRA: Aborting app " << endl;
    GA::Terminate();
  }
  
  dravector.printInternals();

  vector< pair<int,int> > ga_patch;
  ga_patch.push_back( pair<int,int>( 1,nroots));
  ga_patch.push_back( pair<int,int>( 1,nbasis ));

  vector< pair<dra_size_t,dra_size_t> > dra_patch;
  dra_patch.push_back( pair<dra_size_t,dra_size_t>( 1,nroots ));
  dra_patch.push_back( pair<dra_size_t,dra_size_t>( 1,nbasis ));
  
  if ( dravector.dumpGAtoDRA( g_vector ) != 0 ) {
    cout << "Whole GA Dump to DRA failed " << endl;
  }
  
  /*
    if ( dravector.dumpGASectionDRA( ga_patch, dra_patch, g_vector) != 0 ) {
    cout << "Dump to DRA failed " << endl;
    }
  */
  
  // Now try to READ a DRA into a fresh GA.
  
  GA::GlobalArray *g_read_dra = GA::createGA( *g_vector );
  
  string readdraarrayname = "TESTarrayname";
  string readdrafilename = "filenameDRA.txt";

  PsociDRArestart readdravector( g_rank, readdraarrayname, readdrafilename );
  readdravector.openDRAforREAD();
  readdravector.inquireDRAsetREAD();

/*
  ndims = 2;
  ddims.first=nroots;
  ddims.second=nbasis;
  rdims.first=-1;
  rdims.second=-1;
  readdravector.setDimensions( ndims, ddims, rdims );
*/

/* optional for a type of type C_DBL */
/*
  readdravector.reportName();
  readdravector.reportFilename();
  readdravector.setType( C_DBL );
  readdravector.printType();
  readdravector.setMode( DRA_R );
*/
  
  cout << "XXXXX on to the READ " << endl;
  if ( readdravector.readGAfromDRA( g_read_dra ) != 0 ) {
    cout << "Whole GA Read to DRA failed " << endl;
  } 
  
  readdravector.flickDRA();
  int read_status = readdravector.waitDRA();
  readdravector.closeDRA();
  cout << "print READ GA from DRA array " << endl;
  if ( nbasis < 1000 ) {
  g_read_dra->print();
  }
  
// Test writing and reading a patch
/* Use the GA (roots, basis), but create an array DRA(roots, 3*basis), then 
   the write out three copies of the GA into the DRA.
*/

// Get data to store
  GA::GlobalArray *g_patch_vector = GA::createGA( *g_vector );
  g_patch_vector->copy( g_vector ); 
  
  // Create a DRA (nroots, by n * bbasis );
  cout << "Dump to a DRA section " << endl;
  int ncopy = 3;
  
  string patchdraarrayname = "TESTpatcharrayname";
  string patchdrafilename = "filenamepatchDRA.txt";
  
  PsociDRArestart patchdravector( g_rank, patchdraarrayname, patchdrafilename );
  ndims = 2;
  ddims.first = nbasis;
  ddims.second = ncopy * nroots;
  rdims.first=-1;
  rdims.second=-1;
  patchdravector.setDimensions( ndims, ddims, rdims );

  patchdravector.reportName();
  patchdravector.reportFilename();
  patchdravector.setType( C_DBL );
  patchdravector.printType();
  patchdravector.setMode( DRA_RW );

  if ( patchdravector.createDRA() != 0 )
    {
      cout << "failed to create the patch DRA: Aborting app " << endl;
      GA::Terminate();
    }
  
  // Nopw we need to set the size of the GA to write and the placement patch into the DRA. (
  
  // USE FORTRAN style indexing. Conversion occurs underneath
  
  cout << "Prepare to DUMP patch to DRA section " << endl;
  
  vector< pair<int,int> > ga_section_patch;
  ga_section_patch.push_back( pair<int,int>( 1,nbasis));
  ga_section_patch.push_back( pair<int,int>( 1,nroots )); //specify GA size is the patch 

//Cause an intentional failure:
//ncopy = 4;

  double dscale=1.0;
  int patch_shift = 0;
  for (int i=0; i< ncopy; ++i) {
    
    cout << " Dump test for ncopy = " << i+1 << endl;
    
    g_patch_vector->scale( &dscale );
    
    vector< pair<dra_size_t,dra_size_t> > dra_section_patch;
    
    dra_section_patch.push_back( pair<dra_size_t,dra_size_t>( 1, nbasis));
    dra_section_patch.push_back( pair<dra_size_t,dra_size_t>( 1+patch_shift, nroots+patch_shift ));

    patchdravector.dumpGASectionDRA(ga_section_patch, dra_section_patch, g_patch_vector );
    patchdravector.flickDRA();
    patchdravector.waitDRA();
    patch_shift += nroots;
    ++dscale;
  }
  patchdravector.closeDRA();
  cout << "Done with dump section test " << endl;

  // Now we need to test reading the data back in sections.
  // Open one last GA and zero it out.
  
  // Get data to store
  GA::GlobalArray *g_patch_read = GA::createGA( *g_patch_vector );
  g_patch_read->zero();
  
  // Open existing DRA
  
  //  string patchdraarrayname = "TESTpatcharrayname";
  //  string patchdrafilename = "filenamepatchDRA.txt";
  
  PsociDRArestart readdrapatch( g_rank, patchdraarrayname, patchdrafilename );
  readdrapatch.openDRAforREAD();
  readdrapatch.inquireDRAsetREAD();
  
  // Get dimensions, checksize of columns and divide by ncopy.
  
  int patch_ndim;
  pair<dra_size_t,dra_size_t> patch_dims;
  pair<dra_size_t,dra_size_t> patch_rdims; 
  readdrapatch.getDimensions( patch_ndim, patch_dims, patch_rdims );
  
  cout << "On disk patch DRA has dimension " << patch_ndim << endl;
  cout << "On disk patch DRA num rows = " << patch_dims.first << " and num columns = " << patch_dims.second << endl;
  
  // Start reading the on-disk DRA in sections 
  // Always drop contents into the same GA size - so intentioanlly zero before each read
  
  vector< pair<int,int> > ga_section_patch2;
  ga_section_patch2.push_back( pair<int,int>( 1,nbasis));
  ga_section_patch2.push_back( pair<int,int>( 1,nroots ));
  
  patch_shift = 0;
  for (int i=0; i< ncopy; ++i) {
    g_patch_read->zero();
    
    cout << " Dump Read test for ncopy = " << i+1 << endl;
    
    vector< pair<dra_size_t,dra_size_t> > dra_section_patch;
    dra_section_patch.push_back( pair<dra_size_t,dra_size_t>( 1,nbasis ));
    dra_section_patch.push_back( pair<dra_size_t,dra_size_t>( 1+patch_shift,nroots+patch_shift ));
    
    readdrapatch.readGASectionDRA(ga_section_patch2, dra_section_patch, g_patch_read );
    readdrapatch.flickDRA();
    readdrapatch.waitDRA();
    patch_shift += nroots;
    cout << "XXXXXXXXXXXXXXXXXX" << endl;
    if ( nbasis < 1000 ) {
    g_patch_read-> print();
    }
  }
  readdrapatch.closeDRA();
  cout << "Done with read section test " << endl;
  
  
  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX finished XXXXXXXXXXXXXXXXXXX " << endl;
  
  DRAservice.exitDRA();

  cout << "Do the INITIAL guess approach " << endl;
  PsociGArestart resVecInit( g_rank, g_vector );
  string filename4 = "testGAInitGuessvector.txt";

 cout << "Try a parallel initial guess vectors method" << endl;
// resVecInit.specifyDiagonalGuessVectors( -1, g_vector );
  resVecInit.specifyDiagonalGuessVectors( 8, g_vector );

  if ( g_rank == 0 ) {
    resVec.openVectorWrite( g_rank, filename4 );
    double time;
    resVec.dumpVector( g_rank, g_vector, time );
    cout << "Time to dump GUESS ga to disk (s)" << time << endl;
    resVec.closeVectorFile();
  }

 g_vector->print();

  GA::Terminate();

}



