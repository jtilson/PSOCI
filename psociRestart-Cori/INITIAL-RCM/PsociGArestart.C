/**************************************************************************************

* Copyright (c) 2010, 2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: A class for managing restarts for PSOCI: Basically interfaces to PsociVectors
              GA to/from through usual Posix (or equiv) IO

              Generally speaking this GA -> disk could be anything but in practice it is
              intended as a CI vector of a specific format. We expect it as nbasis X nroots
              in rank and shape. Data are transposed for write to disk

              These are ALL LOCAL calls though you can use rank to modify the resultant filenames.
 History:

 Confirmed correct values for ld in GA get/put july 2011


**************************************************************************************/
/**
 *   @file PsociGArestart.C
 *
 */

#include "PsociTimer.hpp"
#include "PsociVector.hpp"
#include "PsociGArestart.hpp"

// Class definitions

using namespace std;

// We require GA is started already ( but not DRA )
// No good way to check ...

PsociGArestart::PsociGArestart()
{
}

PsociGArestart::PsociGArestart( int num, GA::GlobalArray * g_array )
{
#if 0
  cout << " Total number of GA procs is " << GA::nodes();
  cout << " I am process " << GA::nodeid() << endl;
#endif
  char * handleMessage = "PsociGArestart: constructor failed with g_array";
  g_array->checkHandle( handleMessage ); 
}

void PsociGArestart::openVectorWrite( int rank, string filename )
{    
    vectorRestart.setFilename( rank, filename );
    vectorRestart.setMode( VECTOR_WRITE );
    vectorRestart.OpenFile();
}

void PsociGArestart::openVectorRead( int rank, string filename )
{    
    vectorRestart.setFilename( rank, filename );
    vectorRestart.setMode( VECTOR_READ );
    vectorRestart.OpenFile();
}

// Close the file. You want to do this as late as possible to permit buffering to occur
void PsociGArestart::closeVectorFile()
{
    vectorRestart.CloseFile();
}

/** Take the input global array vector, assemble into a single local memory area and pass
    to the vector restart routine.  Someday we want to do this in a distributed fashion 
    but for now keep it simple -- do not have all nodes do this.....
*/
// Only need to open this file on write. You do not want to close immediately because we gain some
// computation/IO overlap. So defer it as long as possible.
// g_vector --> <double>g(nroots,nbasis) is the preferred orientation but no constraints are imposed 
// Store data to file as rows of nbasis separated by nroot

// Keep this a local operation. Ultimately some cores may have different restart needs

void PsociGArestart::dumpVector( int rank, GA::GlobalArray * g_vector, double &time )
{
    double timein = psociTime();
    dumpVector( rank, g_vector );
    time = psociTime() - timein ;
}

void PsociGArestart::dumpVector( int rank, GA::GlobalArray * g_vector )
{
    int type;
    int ndims;
    int dims[MAX_VECTOR_DIMS];
    char * handleMessage = "PsociGArestart: Write: dumpVector: ga handle invalid ";
    g_vector->checkHandle( handleMessage );
    g_vector->inquire( &type, &ndims, dims);

//#if 0
    cout << "DUMP ndims is " << ndims << " type is " << type << endl;
//#endif

    if ( type != C_DBL ) {
       cerr << "A GA of the wrong type was passed to dumpVector. Cannot save: skipping" << endl;
       return;
    }
    if ( ndims != 2 ) {
       cerr << "A GA of the wrong dimensions was passed to dumpVector. Must be = 2: skipping " << ndims << endl;
       return;
    }

// Grab the data columnss and shove into an object for restart.
// for now root == rows, basis == columns

   int ibasislo = 1; 
   int ibasishi = dims[0];

   // int irootlo = 1;
   int iroothi = dims[1];

   if ( iroothi > ibasishi ) {
      cout << "dumpVector: Unexpected matrix dims. iroothi > ibasishi: continuing " << iroothi << " " << ibasishi << endl;
   }
   // int numColumns = iroothi;
//#if 0
   cout << "dumpVector: iroothi is " << iroothi << " ibasishi is " << ibasishi << endl;
//#endif

// Grab the data one column at a time
// C indexing style

   vector<vector<double> > localVectors;
   int lo[MAX_VECTOR_DIMS];
   int hi[MAX_VECTOR_DIMS];
   int ld = 1;

   lo[0] = ibasislo-1; // I do this to keep clear the indexing issues
   hi[0] = ibasishi-1;

   for (int j = 0; j < iroothi; ++j ) {
       lo[1] = j;
       hi[1] = j;
       vector<double> row( ibasishi );
       g_vector->get( lo, hi, (void *)&row[0], &ld );
       localVectors.push_back( row );
   }
//#if 0
   cout << "Size of localVectors is " << localVectors.size() << endl;
   cout << "Current filename to dump " << vectorRestart.FetchFilename() << endl;
   vectorRestart.printWriteType();
//#endif

   pair<int64_t,int64_t> info;
   info.first = ibasishi;
   info.second = iroothi;
#ifdef BINARY
   cout << " Time for BINARY WRITE output of vector is " << vectorRestart.WriteBurstBinary( localVectors, info ) << endl;
#else
   cout << " Time for WRITE output of vector is " << vectorRestart.WriteBurst( localVectors, info ) << endl;
#endif
   vector<double> norm1;
   if ( !vectorRestart.checkNormalization( localVectors , norm1) ) {
      cout << "dumpVector: Warning: at least one unnormalized vector " << endl;
   }
#if 0
   cout << " Time for PRINT output of vector is " << vectorRestart.PrintBurst( localVectors, info ) << endl;
#endif
}

void PsociGArestart::readVector( int rank, GA::GlobalArray * g_vector, double &time )
{
    double timein = psociTime();
    readVector( rank, g_vector );
    time = psociTime() - timein; 
}

/** Take the input global array vector and populate it with the data in the vector restart file 
*/
// Only need to open this file on read. You do not want to close immediately because we may gain some
// computation/IO overlap. SO defer it as long as possible.
// g_vector --> <double>g(nbasis,nroots) is the preferred approach but no constraints are imposed 

// Keep this a local operation. Ultimately some cores may have different restart needs

// A problem: this requires having created an array of the proper dimensions already......
// TODO put in a restartVector augmentation ( gram-schmidt + normalization )

void PsociGArestart::readVector( int rank, GA::GlobalArray * g_vector )
{
    int type;
    int ndims;
    int dims[MAX_VECTOR_DIMS];

    char * handleMessage = "PsociGArestart: Read Vector: ga handle invalid ";
    g_vector->checkHandle( handleMessage );
    g_vector->inquire( &type, &ndims, dims);
#if 0
    cout << "ndims is " << ndims << " type is " << type << endl;
#endif
    if ( type != C_DBL ) {
       cerr << "A GA of the wrong type was passed to readVector. Cannot save: skipping" << endl;
       return;
    }
    if ( ndims != 2 ) {
       cerr << "A GA of the wrong dimensions was passed to readVector. Must be = 2: skipping " << ndims << endl;
       return;
    }

// Read disk data into a localarray

   vector<vector<double> > localVector;
   pair<int64_t,int64_t> info;

#if 0
  cout << endl << "Fetch filename to read is  " << vectorRestart.FetchFilename() << endl;
  vectorRestart.printWriteType();
#endif

#ifdef BINARY
  cout << "Time for BINARY INPUT of vector is " << vectorRestart.ReadBurstBinary( localVector, info ) << endl;;
#else
  cout << "Time for ASCII INPUT of vector is " << vectorRestart.ReadBurst( localVector, info ) << endl;;
#endif
  vector<double> norm1;
  if ( !vectorRestart.checkNormalization( localVector , norm1) ) {
     cout << "Warning: On Read: at least one unnormalized vector " << endl;
  }
//Compare dimensions
   if  ( (info.first != dims[0]) || (info.second != dims[1]) )
   {
      cerr << "readVector error: readVector and g_vector not conforming: aborting " << endl;
      cerr << "restart dims are " << info.first << " " << info.second << " ";
      cerr << "g_vector dims are " << dims[0] << " " << dims[1] << endl;
      GA::Terminate();
   }

// Push the local data ito the preallocated g_array
   int ibasislo = 0; 
   int ibasishi = dims[0]-1;

   // int irootlo = 0;
   int iroothi = dims[1]-1;

// Compare dimensions...
   if ( iroothi > ibasishi ) {
      cout << "dumpVector: Unexpected matrix dims. iroothi > ibasishi " << iroothi << " " << ibasishi << endl;
   }

   int numColumns = iroothi;
#if 0
   cout << "dumpVector: iroothi is " << iroothi << " ibasishi is " << ibasishi << endl;
#endif
// Grab the data one column at a time
// C indexing style

   vector<vector<double> >::iterator it;
   int lo[MAX_VECTOR_DIMS];
   int hi[MAX_VECTOR_DIMS];
   int ld = 1;

   lo[0] = ibasislo; //  C indexing......
   hi[0] = ibasishi;

//Simply place vectors in the same order.... might want to perform a reshuffle at some point but not here
   int root=-1;
   for( it = localVector.begin(); it != localVector.end(); ++it ) {
       ++root;
       lo[1] = root; 
       hi[1] = root; 
       g_vector->put( lo, hi, (void *)&(*it)[0], &ld );
    }
   numReadRoots = root + 1;
   localVector.clear();
}

int PsociGArestart::getNumberVectors()
{
    return( numReadRoots );
}

void PsociGArestart::setNumberVectors( int resetVectorNum)
{
    numReadRoots = resetVectorNum;
}

/* Need a method for simply creating an initial guess. Do not dump the results to disk
   because we do not want to presume the various dI/O related settings have been made

   We still need to ensure the g_vector is of the right size 

   The size of the GA[basis][nroots]
*/

// nroots can be less than the total number
void PsociGArestart::specifyDiagonalGuessVectors( int numroots, GA::GlobalArray * g_array )
{
  char * handleMessage = "specifyDiagonalGuessVectors: non-existant g_array";
  g_array->checkHandle( handleMessage );
  
  if (GA::nodeid() != 0 ) return;
  
  int type;
  int ndims;
  int dims[MAX_VECTOR_DIMS];
  
  g_array->inquire( &type, &ndims, dims);
  
  if ( ndims != 2 ) {
    GA::error( "DiagGuessVec: ndims not 2", -1);
  }
  
  int iroothi = dims[1]-1;
  if ( numroots > dims[1] || numroots > dims[0] ) {
    GA::error("Requested Guess root > GA size", numroots);
  }  
  
  if ( numroots < 0 ) numroots = dims[1]; 
  
  g_array->zero();
  
  int lo[MAX_VECTOR_DIMS];
  int hi[MAX_VECTOR_DIMS];
  int ld = 1;
  
  double dOne = 1.0;
  
  for(int i=0; i< numroots; ++i ) {
    lo[0] = i; 
    hi[0] = i;
    lo[1] = i; 
    hi[1] = i; 
    g_array->put( lo, hi, &dOne, &ld );
  }
  numReadRoots = numroots;
}

