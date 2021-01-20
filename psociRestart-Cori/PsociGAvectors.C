/**************************************************************************************

* Copyright (c) 2010, 2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: A class for managing vectors for PSOCI: Basically interfaces to PsociVectors
              GA to/from through usual Posix (or equiv) IO

              Generally speaking this GA -> disk could be anything but in practice it is
              intended as a CI vector of a specific format. We expect it as nbasis X nroots
              in rank and shape. Data are transposed for write to disk

              These are ALL LOCAL calls though you can use rank to modify the resultant filenames.
 History:

 Confirmed correct values for ld in GA get/put july 2011


**************************************************************************************/
/**
 *   @file PsociGAvectors.C
 *
 */

#include "PsociTimer.hpp"
#include "PsociVector.hpp"
#include "PsociGAvectors.hpp"

// Class definitions

using namespace std;

// We require GA is started already ( but not DRA )
// No good way to check ...

PsociGAvectors::PsociGAvectors( int nrootsSought, GA::GlobalArray * g_array )
{
  char handleMessage[] = "PsociGAvectors: constructor failed with g_array";
  g_array->checkHandle( handleMessage ); 
  
  if ( nrootsSought < 1 ) {
    cerr << "Nonsensical numrootsSought " << endl;
    GA::error( "Nonsensical numrootsSought", nrootsSought);
  }
  numSoughtRoots = nrootsSought;
  numReadRoots = 0;
  
  int type;
  int ndims;
  int dims[MAX_VECTOR_DIMS];
  g_array->inquire( &type, &ndims, dims);
  if ( ndims != 2 ) { 
    GA::error( "g_array: dims not 2", -1);
  }
  
  if ( dims[1] < numSoughtRoots || 10*dims[1] < numSoughtRoots ) {
    cout << "Advisory: dims[1] is " << dims[1] << " and SoughtRoots is " << numSoughtRoots << endl;
    cout << "This means the maximum subspace is fairly small and may inhibit/preclude convergence " << endl;
  }  
}

void PsociGAvectors::openVectorWrite( int rank, string filename )
{    
  if ( rank != 0 ) return;
  vectorRestart.setFilename( rank, filename );
  vectorRestart.setMode( VECTOR_WRITE );
  vectorRestart.OpenFile();
}

void PsociGAvectors::openVectorRead( int rank, string filename )
{    
  if ( rank != 0 ) return;
  vectorRestart.setFilename( rank, filename );
  vectorRestart.setMode( VECTOR_READ );
  vectorRestart.OpenFile();
}

// Close the file. You want to do this as late as possible to permit buffering to occur
void PsociGAvectors::closeVectorFile()
{
  if ( GA::nodeid() != 0 ) return;
  vectorRestart.CloseFile();
}

/** Take the input global array vector, assemble into a single local memory area and pass
    to the vector vectors routine.  Someday we want to do this in a distributed fashion 
    but for now keep it simple -- do not have all nodes do this.....
*/
// Only need to open this file on write. You do not want to close immediately because we gain some
// computation/IO overlap. So defer it as long as possible.
// g_vector --> <double>g(nroots,nbasis) is the preferred orientation but no constraints are imposed 
// Store data to file as rows of nbasis separated by nroot

// Keep this a local operation. Ultimately some cores may have different vectors needs

void PsociGAvectors::dumpVector( int rank, GA::GlobalArray * g_vector, double &time )
{
  double timein = psociTime();
  dumpVector( rank, g_vector );
  time = psociTime() - timein ;
}

void PsociGAvectors::dumpVector( int rank, int num, GA::GlobalArray * g_vector, double &time )
{   
  double timein = psociTime();
  dumpVector( rank, num, g_vector );
  time = psociTime() - timein ;
} 
//Local BUT as a convenience only rank 0 does anything
//Only print lowest NUM vectors to disk

void PsociGAvectors::dumpVector( int rank, int num, GA::GlobalArray * g_vector )
{
  if ( rank != 0 ) return;

  int type;
  int ndims;
  int dims[MAX_VECTOR_DIMS];
  char handleMessage[] = "PsociGAvectors: Write: dumpVector: ga handle invalid ";
  g_vector->checkHandle( handleMessage );
  g_vector->inquire( &type, &ndims, dims);
  
#if 0
  cout << "DUMP ndims is " << ndims << " type is " << type << endl;
#endif
  
  if ( type != C_DBL ) {
    cerr << "A GA of the wrong type was passed to dumpVector. Cannot save: skipping" << endl;
    return;
  }
  if ( ndims != 2 ) {
    cerr << "A GA of the wrong dimensions was passed to dumpVector. Must be = 2: skipping " << ndims << endl;
    return;
  }
  
  // Grab the data columnss and shove into an object for vectors.
  // for now root == rows, basis == columns
  
  int ibasislo = 1; 
  int ibasishi = dims[0];
  
  // int irootlo = 1;
  int iroothi = dims[1]; 
  if ( num > iroothi ) {
     GA::error(" num > iroothi ", num);
  }
  int soughtHi = num;
  
#if 0
  cout << "dumpVector-2: num is " << soughtHi << " ibasishi is " << ibasishi << endl;
#endif
  
  // Grab the data one column at a time
  // C indexing style
  
  vector<vector<double> > localVectors;
  int lo[MAX_VECTOR_DIMS];
  int hi[MAX_VECTOR_DIMS];
  int ld = 1;
  
  lo[0] = ibasislo-1; // I do this to keep clear the indexing issues
  hi[0] = ibasishi-1;
  
  for (int j = 0; j < soughtHi; ++j ) {
    lo[1] = j;
    hi[1] = j;
    vector<double> row( ibasishi );
    g_vector->get( lo, hi, &row[0], &ld );
    localVectors.push_back( row );
  }

#if 0
  cout << "Size of localVectors is " << localVectors.size() << endl;
  cout << "Current filename to dump " << vectorRestart.FetchFilename() << endl;
  vectorRestart.printWriteType();
#endif
  
  pair<int64_t,int64_t> info;
  info.first = ibasishi;
  info.second = soughtHi;
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

/* New scalable write version 
   More Scalable in a memory sense. Root 0 still does all the fetching and assmebly work

   We can dump a partial vector depending on the value of width

*/
void PsociGAvectors::dumpVectorScalable( int rank, int num, GA::GlobalArray * g_vector )
{
  if ( rank != 0 ) return;

  int type;
  int ndims;
  int dims[MAX_VECTOR_DIMS];
  char handleMessage[] = "PsociGAvectors: Write: dumpVectorScalable: ga handle invalid ";
  g_vector->checkHandle( handleMessage );
  g_vector->inquire( &type, &ndims, dims);
  
#if 0
  cout << "DUMP/Scalable ndims is " << ndims << " type is " << type << endl;
#endif
  
  if ( type != C_DBL ) {
    cerr << "A GA of the wrong type was passed to dumpVectorScalable. Cannot save: skipping" << endl;
    return;
  }
  if ( ndims != 2 ) {
    cerr << "A GA of the wrong dimensions was passed to dumpVectorScalable. Must be = 2: skipping " << ndims << endl;
    return;
  }
  
  // Grab the data columnss and shove into an object for vectors.
  // for now root == rows, basis == columns
  
  int ibasislo = 1; 
  int ibasishi = dims[0];
  
  // int irootlo = 1;
  int iroothi = dims[1]; 
  if ( num > iroothi ) {
     GA::error(" num > iroothi:scalable ", num);
  }
  int soughtHi = num;
  
  int lo[MAX_VECTOR_DIMS];
  int hi[MAX_VECTOR_DIMS];
  int ld = 1;
  
  pair<int64_t,int64_t> info;
  info.first = ibasishi; 
  info.second = soughtHi;

#if 0
  cout << "dumpVector-2 Scalable: num is " << soughtHi << " ibasishi is " << ibasishi << endl;
#endif
  
  // Grab the data one column at a time
  // C indexing style
  
// Find a "typical" width per node
// g_vector->distribution( GA::nodeid(), lo, hi );

#ifdef BINARY
 cout << " SCALABLE BINARY WRITE HEADER output of vector is " << endl;
 int wnum = vectorRestart.scalable_WriteHeaderBurstBinary( info );
#else
 cout << " SCALABLE WRITE HEADER output of vector is " << endl;
 int wnum = vectorRestart.scalable_WriteHeaderBurst( info ) ;
#endif

//  int width = ibasishi;
//  vector<double> localVectors(width); // no more than one vector at a time; generally much less than that
  
// Arbitrarily choose a maximunm width of a single CI vector to fetch and dump

  int width = MAX_VECTOR_RESTART_WIDTH; // CHANGE ME

  vector<double> localVectors(width); // no more than one vector at a time; generally much less than that

 cout << "dump sought " << soughtHi << endl;
  for (int j = 0; j < soughtHi; ++j ) {
    lo[1] = j;
    hi[1] = j;

    bool start = true;
    int rootnum = j;

    for( int ibas=0; ibas< ibasishi; ibas+=width ) 
   { 
    lo[0] = ibas; // I do this to keep clear the indexing issues
    int hirange = min( ibas+width, ibasishi );
    hi[0] = hirange-1;
    localVectors.resize( hirange - ibas ); // Burst needs it properly sized
    g_vector->get( lo, hi, &localVectors[0], &ld );

#if 0
  cout << "Size of localVectors is " << localVectors.size() << endl;
  cout << "Current filename to dump " << vectorRestart.FetchFilename() << endl;
  vectorRestart.printWriteType();
#endif
  
#ifdef BINARY
    wnum = vectorRestart.WriteBurstBinaryScalable(start, rootnum+1, localVectors, info );
#else
    wnum = vectorRestart.WriteBurstScalable(start, rootnum+1, localVectors, info );
#endif

    start = false;

    } // hirange loop
#ifndef BINARY
    vectorRestart.WriteBurstScalableEOL();
#endif
   } //soughtHi

// Dump trailing metadata same code BINARY or ASCII
// title defined using setTitle() call

#ifdef BINARY
    wnum = vectorRestart.WriteBurstScalableBinaryTrailer();
#else
    wnum = vectorRestart.WriteBurstScalableTrailer();
#endif

    cout << "Time for Scalable CI vector Write is " << wnum << endl;
}

//Local BUT as a convenience only rank 0 does anything
void PsociGAvectors::dumpVector( int rank, GA::GlobalArray * g_vector )
{
  if ( rank != 0 ) return;

  int type;
  int ndims;
  int dims[MAX_VECTOR_DIMS];
  char handleMessage[] = "PsociGAvectors: Write: dumpVector: ga handle invalid ";
  g_vector->checkHandle( handleMessage );
  g_vector->inquire( &type, &ndims, dims);
  
#if 0
  cout << "DUMP ndims is " << ndims << " type is " << type << endl;
#endif
  
  if ( type != C_DBL ) {
    cerr << "A GA of the wrong type was passed to dumpVector. Cannot save: skipping" << endl;
    return;
  }
  if ( ndims != 2 ) {
    cerr << "A GA of the wrong dimensions was passed to dumpVector. Must be = 2: skipping " << ndims << endl;
    return;
  }
  
  // Grab the data columnss and shove into an object for vectors.
  // for now root == rows, basis == columns
  
  int ibasislo = 1; 
  int ibasishi = dims[0];
  
  // int irootlo = 1;
  int iroothi = dims[1]; 
  
  int soughtHi = iroothi;
  if ( iroothi != numSoughtRoots ) {
    cout << "Attn: iroothi != numSoughtRoot: Dumping numSoughtRoots to disk " << numSoughtRoots << endl;
    soughtHi = numSoughtRoots;
  }
  
  if ( iroothi > soughtHi ) {
    cout << "dumpVector: Unexpected matrix dims. iroothi > isoughtHi: continuing " << iroothi << " " << soughtHi << endl;
  }
  // int numColumns = iroothi;
#if 0
  cout << "dumpVector: iroothi is " << soughtHi << " ibasishi is " << ibasishi << endl;
#endif
  
  // Grab the data one column at a time
  // C indexing style
  
  vector<vector<double> > localVectors;
  int lo[MAX_VECTOR_DIMS];
  int hi[MAX_VECTOR_DIMS];
  int ld = 1;
  
  lo[0] = ibasislo-1; // I do this to keep clear the indexing issues
  hi[0] = ibasishi-1;
  
  for (int j = 0; j < soughtHi; ++j ) {
    lo[1] = j;
    hi[1] = j;
    vector<double> row( ibasishi );
    g_vector->get( lo, hi, &row[0], &ld );
    localVectors.push_back( row );
  }
#if 0
  cout << "Size of localVectors is " << localVectors.size() << endl;
  cout << "Current filename to dump " << vectorRestart.FetchFilename() << endl;
  vectorRestart.printWriteType();
#endif
  
  pair<int64_t,int64_t> info;
  info.first = ibasishi;
  info.second = soughtHi;
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

//Local

void PsociGAvectors::readVector( int rank, GA::GlobalArray * g_vector, double &time )
{
  double timein = psociTime();
  readVector( rank, g_vector );
  time = psociTime() - timein; 
}

/** Take the input global array vector and populate it with the data in the vector vectors file 
*/
// Only need to open this file on read. You do not want to close immediately because we may gain some
// computation/IO overlap. SO defer it as long as possible.
// g_vector --> <double>g(nbasis,nroots) is the preferred approach but no constraints are imposed 

// Keep this a local operation. Ultimately some cores may have different vectors needs

// A problem: this requires having created an array of the proper dimensions already......
// TODO put in a vectorsVector augmentation ( gram-schmidt + normalization )

/* Read all vectors but PUT only numSoughtRoots from disk: g_array was only created for that many
   This way we can read the trailing TITLE file in the future
*/
void PsociGAvectors::readVector( int rank, GA::GlobalArray * g_vector )
{
  if ( rank != 0 ) return;

  int type;
  int ndims;
  int dims[MAX_VECTOR_DIMS];
  
  char handleMessage[] = "PsociGAvectors: Read Vector: ga handle invalid ";
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
  if  ( (info.first != dims[0])  )
    {
      cerr << "readVector error: readVector and g_vector not conforming: aborting " << endl;
      cerr << "vectors dims are " << info.first << " " << info.second << " ";
      cerr << "g_vector dims are " << dims[0] << " " << dims[1] << endl;
      GA::Terminate();
    }
  if  ( (info.second > dims[1]) )
    {
      cout << "Will only KEEP this many = " << dims[1] << endl;;
    }
  
  // Push the local data into the preallocated g_array
  int ibasislo = 0; 
  int ibasishi = dims[0]-1;
  
  // int irootlo = 0;
  int iroothi = dims[1]-1;
  
  // Compare dimensions...
  if ( iroothi > ibasishi ) {
    cout << "dumpVector: Unexpected matrix dims. iroothi > ibasishi " << iroothi << " " << ibasishi << endl;
  }
  
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
//truncate order
  int root=-1;
  for( it = localVector.begin(); it != localVector.end(); ++it ) {
    if ( root < dims[1]-1 ) {
    ++root;
    lo[1] = root; 
    hi[1] = root; 
    g_vector->put( lo, hi, &(*it)[0], &ld );
    } // else do not upload any more vectors to GA
  }
  numReadRoots = root + 1;
  localVector.clear();
}

/* A new method to read CI vectors from disk to GA in a way to NOT
   require reading it all into memory

*/
void PsociGAvectors::readVectorScalable( int rank, GA::GlobalArray * g_vector )
{
  if ( rank != 0 ) return;

  int type;
  int ndims;
  int dims[MAX_VECTOR_DIMS];
  
  char handleMessage[] = "PsociGAvectorsScalable: Read Vector: ga handle invalid ";
  g_vector->checkHandle( handleMessage );
  g_vector->inquire( &type, &ndims, dims);

#if 0
  cout << "ndims is " << ndims << " type is " << type << endl;
#endif

  if ( type != C_DBL ) {
    cerr << "Scalable: A GA of the wrong type was passed to readVector. Cannot save: skipping" << endl;
    return;
  }
  if ( ndims != 2 ) {
    cerr << "Scalable: A GA of the wrong dimensions was passed to readVector. Must be = 2: skipping " << ndims << endl;
    return;
  }
  
  // Read disk header into local objects
  
  pair<int64_t,int64_t> info;

#ifdef BINARY
 int wnum = vectorRestart.scalable_ReadHeaderBurstBinary( info ) ;
#else
 int wnum = vectorRestart.scalable_ReadHeaderBurst( info ) ;
#endif

  // Grab the data one column at a time
  // C indexing style

#if 0
  cout << endl << "Fetch filename to read is  " << vectorRestart.FetchFilename() << endl;
  vectorRestart.printWriteType();
#endif
  
  //deleted normalization check

  int ntotalRootsOnDisk = info.second;

  //Compare dimensions
  if  ( (info.first != dims[0])  )
    {
      cerr << "vectors dims are " << info.first << " " << info.second << " ";
      cerr << "g_vector dims are " << dims[0] << " " << dims[1] << endl;
      GA::error( " Scalable readVector error: readVector and g_vector not conforming: aborting ",1);
    }
  if  ( (info.second > dims[1]) )
    {
      cout << "Scalable: Will only KEEP this many vectors on read = " << dims[1] << endl;;
    }

  int nroots = min( ntotalRootsOnDisk,  dims[1] );
  
  // Push the local data into the preallocated g_array
  int ibasislo = 0; 
  int ibasishi = dims[0]-1;
  int irootlo = 0;
  int iroothi = dims[1]-1;
  int nbasis = ibasishi + 1;
  
#if 0
  cout << "readVector: iroothi is " << iroothi << " ibasishi is " << ibasishi << endl;
#endif
  
  vector<vector<double> >::iterator it;
  int lo[MAX_VECTOR_DIMS];
  int hi[MAX_VECTOR_DIMS];
  int ld = 1;

  int width = MAX_VECTOR_RESTART_WIDTH; // arbitrary width for testing

  vector<double> localVectors(width); // no more than one vector at a time; generally much less than that
  
  for (int j = 0; j < nroots; ++j ) {

    lo[1] = j;
    hi[1] = j;

    bool start = true;
    int rootnum = j;

    for( int ibas=0; ibas< nbasis; ibas+=width )
   {
    lo[0] = ibas; // I do this to keep clear the indexing issues
    int hirange = min( ibas+width, nbasis );
    hi[0] = hirange-1;
    localVectors.resize( hirange - ibas ); // Burst needs it properly sized

    int numelems = hirange - ibas;
#ifdef BINARY
    wnum = vectorRestart.ReadBurstScalableBinary(start, numelems, localVectors, info );
#else
    wnum = vectorRestart.ReadBurstScalable(start, numelems, localVectors, info );
#endif

    start = false; //trigger for reading root number
    g_vector->put( lo, hi, &localVectors[0], &ld );

    } // hirange loop

    ++numReadRoots;

   } //soughtHi

    cout << "Time for Scalable CI vector READ is " << wnum << endl;


/* At this point we MAY want to read the final meta data 
*/

}


int PsociGAvectors::getNumberVectors()
{
  return( numReadRoots );
}

int PsociGAvectors::getNumberVectorsSought()
{
  return( numSoughtRoots );
}

void PsociGAvectors::setNumberVectors( int resetVectorNum)
  {
  numSoughtRoots = resetVectorNum;
  }

/* Need a method for simply creating an initial guess. Do not dump the results to disk
   because we do not want to presume the various dI/O related settings have been made
   
   We still need to ensure the g_vector is of the right size 
   The size of the GA[basis][nroots]
*/

// nroots can be less than the total number AND we can append guess on top of read vectors
// We do not orthogonalize them here, though.

//Local - BUT only rank == 0 actually does anything here.

void PsociGAvectors::specifyDiagonalGuessVectors(GA::GlobalArray * g_array )
{
  char handleMessage[] = "specifyDiagonalGuessVectors: non-existant g_array";
  g_array->checkHandle( handleMessage );
  
  if (GA::nodeid() != 0 ) return;
  
  /* If numSoughtRoot <=  numReadRoot ignore.
     Else APPEND diagonal guess vectors (non-orthogonal) to the end up to numSoughtRoots.
  */

  int type;
  int ndims;
  int dims[MAX_VECTOR_DIMS];
  g_array->inquire( &type, &ndims, dims);
  if ( ndims != 2 ) {
    GA::error( "DiagGuessVec: ndims not 2", -1);
  }
  
  if ( numSoughtRoots > numReadRoots ) {
    
#if 0
    cout << "ADDING TO RESTART " << endl;
#endif
    
    int lo[MAX_VECTOR_DIMS];
    int hi[MAX_VECTOR_DIMS];
    int ld = 1;
    
    double dOne = 1.0;
    
    for(int i=numReadRoots; i< numSoughtRoots; ++i ) {
      lo[0] = i; 
      hi[0] = i;
      lo[1] = i; 
      hi[1] = i; 
      g_array->put( lo, hi, &dOne, &ld );
    }
  }
}

// Use this to only get the numSoughtRoots, else simply do a ga_print.

//Local BUT only root 0 does anything.

void PsociGAvectors::printCurrentSoughtRoots( GA::GlobalArray * g_array )  
{
  if ( GA::nodeid() != 0 ) return;
  
  if ( numSoughtRoots > 0 ) {
    
    int type;
    int ndims;
    int dims[MAX_VECTOR_DIMS];
    g_array->inquire( &type, &ndims, dims);
    if ( ndims != 2 ) {
      GA::error( "g_array: dims not 2", -1);
    }
    
    int lo[MAX_VECTOR_DIMS];
    int hi[MAX_VECTOR_DIMS];
    
    int nbasis = dims[0];
    
    if ( dims[1] < numSoughtRoots || 10*dims[1] < numSoughtRoots ) {
      cout << "Advisory: dims[1] is " << dims[1] << " and SoughtRoots is " << numSoughtRoots << endl;
      cout << "This means the maximum subspace is fairly small and may inhibit/preclude convergence " << endl;
    }
    
    cout << "Printing out the Sought Roots " << dims[0] << " " << numSoughtRoots << endl;
    
    vector<double> vector( dims[0],0 );
    for(int i=0; i< numSoughtRoots; ++i ) {
      //cout << "Data for diagonal root number " << i << endl;
      
      lo[0]=0;
      hi[0]=dims[0]-1;
      
      lo[1]=i;
      hi[1]=i;
      int ldv=1;
      
      g_array->get( lo, hi, &vector[0], &ldv);
      
      for (int j=0; j<nbasis; ++j) {
	cout << vector[j] << " ";
      }
      cout << endl;
    }
  }  
}



