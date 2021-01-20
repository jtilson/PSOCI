/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: A simple class for managing restarts for PSOCI: Basically this writes and reads
              CI vectors 

 History:



**************************************************************************************/
/**
 *   @file PsociVector.C
 *
 */

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

/* For future considerations
   #include <cstdio> // Captures the standard C file IO
   #include <fcntl.h> // Captures the standard C file IO
   #include <unistd.h> // Captures Unix IO on the BG/L
*/

#include "PsociTimer.hpp"
#include "PsociVector.hpp"

//fstream olocalfile; 

// Class definitions

using namespace std;

PsociVector::PsociVector()
{
  blk_size = 8192;  
  set_number = 0;
  string filename = "psoci_vector.txt";
  precision.first = MAXWIDTH;
  precision.second = MAXDECIMAL;
  current_title="default vector test format ";
  Filename( set_number, filename );
}

PsociVector::PsociVector( const string filename )
{
  blk_size = 8192;
  set_number = -1;
  precision.first = MAXWIDTH;
  precision.second = MAXDECIMAL;
  current_title="default vector test format ";
  Filename( set_number, filename );
}

PsociVector::PsociVector( const int num, const string filename )
{
  blk_size = 8192;
  set_number = num;
  precision.first = MAXWIDTH;
  precision.second = MAXDECIMAL;
  current_title="default";
  Filename( set_number, filename );
}

PsociVector::PsociVector(const int num, const string filename, const pair<int,int> prec )
{
  blk_size = 8192;
  set_number = num;
  precision.first = prec.first;
  precision.second = prec.second;
  current_title="default";
  Filename( set_number, filename );
}

void PsociVector::setFilename( const int num, const string filename )
{
  set_number = num;
  current_title="default";
  Filename( set_number, filename );
}

void PsociVector::setFilename( const int num, const string filename, const pair<int,int> prec )
{
  precision.first = prec.first;
  precision.second = prec.second;
  set_number = num;
  current_title="default";
  Filename( set_number, filename );
}


void PsociVector::setTitle( const string title )
{ 
  current_title = title;
}

void PsociVector::printTitle()
{
  cout << "Current title is " << current_title << endl;
}

void PsociVector::Params() 
{
  cout << "---------------------------------------" << endl;
  cout << "Block size (B)   = " << blk_size << endl;
  cout << "set_number       = " << set_number << endl;
  cout << "Current filename = " << current_filename << endl;
  cout << "Current data title = " << current_title << endl;
  cout << "---------------------------------------" << endl << endl;
}

inline int PsociVector::OpenFileSeq()
{
#if 0
  cout << "C++ style file IO " << endl;
#endif
  
  if ( current_mode != VECTOR_WRITE && current_mode != VECTOR_READ ) {
    cerr << "Error: You MUST specify PsociVector::setMode(VECTOR_WRITE || VECTOR_READ) before attempting an open: I am going to abort to protect any pre-existing files from damage " << endl;
    return ( PV_OPEN_FAILED );
  }
  olocalfile.open( current_filename.c_str(), current_mode );
  if (!olocalfile.is_open()) {
    cerr << "OpenFile: Cannot open PsociVector output file " << current_filename << " Aborting write " << endl;
    return ( PV_OPEN_FAILED );
  }
  return ( PV_SUCCESS );
}


void PsociVector::CloseFileSeq()
{
  if ( current_mode == VECTOR_WRITE ) {
    olocalfile.flush();
  }
  //sync();
  olocalfile.close();
#ifdef REMOVEFILES
  // remove( current_filename.c_str()  );
#endif
}

void PsociVector::setMode( ios_base::openmode MODE )
{
  current_mode = MODE;
}

ios_base::openmode PsociVector::fetchMode()
{
  return (current_mode);
}

int PsociVector::OpenFile()
{
  return ( OpenFileSeq() );
}

void PsociVector::CloseFile()
{
  CloseFileSeq();
}

/** Construct a reasonable filename for the text: Assume we need to 
    prevent name clashing in a parallel environment
    One can explicitly incorporate the path as well. Insert the
    quantity: -(num) before the last . in .ext. If no extension
    was specified, then -num is simply appended

    if num == -1, then no processing is done at all.
*/

void PsociVector::Filename( const int num, const string filename )
{
  //Convert num to a string with a prepended hyphen .....are you kidding me 
  if ( num == -1 ) {
    current_filename = filename;
    return;
  }
  
  string tempString = filename;
  string scount;
  stringstream sout;
  sout << num;
  scount = sout.str();
  scount.insert(0,"-");
  //yuk. Now find location to prepend 
  
  int locate = tempString.find_last_of(".");
  if ( locate != -1 ) {
    tempString.insert( locate, scount );
  } else {
    tempString.append( scount );
  }
  current_filename = tempString;
}

string PsociVector::FetchFilename()
{
  return ( current_filename );
}

// TODO error checking

/* vector<vector<double> > is basically a 2D array of CI vector coefficients
   They are ordered as V(root,basis). So each outer iterator points to a CI set
   of coefficients.  The basis are currently on a unit metric and precision needs to
   be enforced.
   
   For now we desire roots to be ordered from lowest E to highest.
   
   Check to ensure consistent values for info
*/

//C++ I/O (streams): Print out the coefficients:  ASCII approach Not useful for large arrays

void PsociVector::PrintBurstCPPversion( vector<vector<double> > & vecs, pair<int64_t,int64_t> info )
{
  int64_t nbasis = info.first;
  int width = precision.first;
  int prec  = precision.second;
  
#if 0
  cout << "CurrentTitle is " << current_title << endl;
  cout << "Current precision info (width/decimals) " << width << " " << prec << endl;
  cout << "Number of roots is " << vecs.size() << endl;
#endif

  if ( nbasis < MAXPRINTCOEFS ) {
    int root_num=0;
    cout << setprecision(prec) << setiosflags(ios::fixed);
    
    vector<vector<double> >::const_iterator outer;
    for( outer=vecs.begin(); outer!=vecs.end(); ++outer) 
      {
	++root_num;
	cout <<  root_num << endl;
	
	vector<double>::const_iterator it;
	if ( (*outer).size() != nbasis ) {
	  cerr << "PrintBurstCPP: Error writing vector: nbasis and actual are different " <<  (*outer).size()  << endl;
	  exit(1);
	}
#if 0
	cout << "Number of basis is " << (*outer).size() << endl;
#endif
	for (it=(*outer).begin(); it != (*outer).end(); ++it)
	  {
	    cout << setw(width) << (*it) ;
	  }
	cout << endl;
      }
  } else
    {
      cout << "nbasis is larger than MAXPRINTCOEFS: " << nbasis << " " << MAXPRINTCOEFS << " Printing coefs to stdout is cancelled " << endl;
    }
}

//C++ approach I/O (streams) Write coefs to file: ASCII approach
int PsociVector::WriteBurstCPPversion( vector<vector<double> > & vecs, const pair<int64_t,int64_t> info )
{
  if (!olocalfile.is_open()) {
    cerr << "Write: Cannot open PsociVector output file " << current_filename << " Aborting write " << endl;    
    exit(1);
  } 
  cout << "WriteBurstCPPversion" << endl;
  int64_t nbasis = info.first;
  int64_t nroots = info.second;
  int width = precision.first;
  int prec  = precision.second;
  
  olocalfile << nbasis << " " << nroots << endl;
  olocalfile << width << " " << prec << endl;
  
  int root_num=0;
  olocalfile << setprecision(prec) << setiosflags(ios::fixed);
  
  vector<vector<double> >::const_iterator outer;
  for( outer=vecs.begin(); outer!=vecs.end(); ++outer) 
    {
      ++root_num;
      olocalfile <<  root_num << endl;
      
      vector<double>::const_iterator it;
      if ( (*outer).size() != nbasis ) {
	cerr << "WriteBurstCPP: Error writing vector: nbasis and actual are different " <<  (*outer).size() << endl;
	return ( PV_WRITE_FAILED );
      }
      for (it=(*outer).begin(); it != (*outer).end(); ++it)
	{
	 olocalfile << setw(width) << (*it) ;
	}
      olocalfile << endl;
    }
  olocalfile << current_title << endl;
  return(PV_SUCCESS);
}

int PsociVector::WriteBurstScalableEOL()
{
  if (!olocalfile.is_open()) {
    cerr << "Cannot open PsociVector Scalable output file " << current_filename << " Aborting write " << endl;
    return( PV_WRITE_FAILED );
  }
  olocalfile << endl;
  return(PV_SUCCESS);
}

int PsociVector::WriteBurstCPPversionScalable(bool start, int rootnum,  vector<double> & vecs, const pair<int64_t,int64_t> info )
{
  if (!olocalfile.is_open()) {
    cerr << "Cannot open PsociVector Scalable output file " << current_filename << " Aborting write " << endl;
    return( PV_WRITE_FAILED );
  }
  int64_t nbasis = info.first;
  int width = precision.first;
  int prec  = precision.second;

  vector<double>::iterator it;
  int root_num=rootnum;
  if ( start == 1 ) {
  olocalfile << setprecision(prec) << setiosflags(ios::fixed);
  olocalfile <<  root_num << endl;
  }
  for(it=vecs.begin(); it!=vecs.end(); ++it) {
      olocalfile << setw(width) << (*it) ;
  }
//  olocalfile << endl; //We've carved up each CI vector and now must do this separately.
  return( PV_SUCCESS );
}

//C++ Read vector coefs from disk:  I/O (streams) ASCII
int PsociVector::ReadBurstCPPversion( vector<vector<double> > & vecs, pair<int64_t,int64_t> & info )
{
  if (!olocalfile.is_open()) {
    cerr << "Cannot open PsociVector file for READ" << current_filename << " Aborting write " << endl;
    return( PV_READ_FAILED );
  }
  
  // int width = precision.first;
  // int prec  = precision.second;
  int64_t nroots, nbasis;
  int readwidth, readprec;
  
  olocalfile >> nbasis >> nroots;
  olocalfile >> readwidth >>  readprec ;
  
  info.first = nbasis;
  info.second = nroots;
  
#if 0
  cout << "ASCII Read nroots and nbasis as " << nroots << " " << nbasis << endl;
  cout << "ASCII Read data width and prec are " << readwidth << " " << readprec << endl;
#endif

  int64_t current_root;
  double value;
  
  for ( int64_t i=0; i < nroots; ++i) {
    olocalfile >> current_root;
    vector<double> coefs;
    for ( int64_t j=0; j< nbasis; ++j ) {
      olocalfile >> value;
      coefs.push_back(value);
    }
    vecs.push_back( coefs );
  }
  getline(olocalfile, current_title);
  return( PV_SUCCESS );
}

//C++ approach I/O BINARY approach
int PsociVector::ReadBurstCPPversionBinary( vector<vector<double> > & vecs, pair<int64_t,int64_t> & info )
{
  if (!olocalfile.is_open()) {
    cerr << "Cannot open PsociVector file for READ" << current_filename << " Aborting read " << endl;
    return( PV_READ_FAILED );
  }
  pair<int,int> precision;
  
  // Start reading the data file
  
  int64_t infodata[2];
  olocalfile.read((char*)infodata, 2 * sizeof(int64_t) );
  info.first=infodata[0];
  info.second=infodata[1];
  int64_t nbasis = info.first;
  int64_t nroots = info.second;
  
  int precdata[2];
  olocalfile.read((char*)&precdata, 2 * sizeof(int) );
  // int width = precdata[0];
  // int prec = precdata[1];
  
#if 0
  cout << "Read nroots and nbasis as " << nroots << " " << nbasis << endl;
  cout << "Read data width and prec are " << width << " " << prec << endl;
#endif

  int64_t current_root;
  
  for (int i=0;i<nroots; ++i) {
    vector<double> tempvec;
    tempvec.resize(nbasis);
    
    olocalfile.read( (char *) &current_root, 1 * sizeof(int64_t) );
    olocalfile.read( (char *)(&tempvec[0]), nbasis * sizeof(double) );
    vecs.push_back( tempvec );
  }
  getline(olocalfile, current_title);
  return( PV_SUCCESS );
}

/* The header must have been read already and the calling app is taking care of when to treat EOL, etc.
*/
int PsociVector::ReadBurstCPPversionScalable(bool start, int range, vector<double> & vecs, pair<int64_t,int64_t> & info )
{
  if (!olocalfile.is_open()) {
    cerr << "Cannot open Scalable PsociVector file for READ" << current_filename << " Aborting write " << endl;
    return( PV_READ_FAILED );
  }
  
  vecs.clear();
  int64_t nroots, nbasis;
  int readwidth, readprec;
  
  int hi = range;
  int64_t current_root;
  double value;
  
  if ( start == 1 ) {
    olocalfile >> current_root;
  }

  int count=0;
  for ( int j=0; j< hi; ++j ) {
  ++count;
    double value;
    olocalfile >> value;
    vecs.push_back(value);
  }

  return( PV_SUCCESS );
}

/* BINARY version 
   On entry vecs has already been properly sized for the indicated range of basis functions
*/
int PsociVector::ReadBurstCPPversionScalableBinary(bool start, int range, vector<double> & vecs, pair<int64_t,int64_t> & info )
{
  if (!olocalfile.is_open()) {
    cerr << "Cannot open Scalable BINARY PsociVector file for READ" << current_filename << " Aborting write " << endl;
    return( PV_READ_FAILED );
  }
  
  vecs.clear();
  int64_t nroots, nbasis;
  int readwidth, readprec;
  
  int hi = range;
  int64_t current_root;
  double value;
  
  if ( start == 1 ) {
    olocalfile.read( (char *) &current_root, 1 * sizeof(int64_t) );
  }

  olocalfile.read( (char *)( &vecs[0]), range * sizeof(double) );

  return( PV_SUCCESS );
}

    
//C++ approach I/O (streams) BINARY approach
int PsociVector::WriteBurstCPPversionBinary( vector<vector<double> > & vecs, const pair<int64_t,int64_t> info )
{
  if (!olocalfile.is_open()) {
    cerr << "Cannot open PsociVector output file " << current_filename << " Aborting write " << endl;
    return( PV_WRITE_FAILED );
  }
  
  // Lay down the values for nroots, nbasis, width, and prec
  int64_t infodata[2];
  int64_t nbasis = info.first;
  int64_t nroots = info.second;
  infodata[0] = nbasis;
  infodata[1] = nroots;
  olocalfile.write((char*)infodata, 2 * sizeof(int64_t) );
  
  int precdata[2];
  int width = precision.first;
  int prec  = precision.second;
  precdata[0] = width;
  precdata[1] = prec;
  olocalfile.write((char*)&precdata, 2 * sizeof(int) );
  
  int root_num=0;
  vector<vector<double> >::const_iterator outer;
  for( outer=vecs.begin(); outer!=vecs.end(); ++outer)
    {
      ++root_num;
      olocalfile.write( (char*)&root_num, sizeof(int64_t)*1);
      
      if ( (*outer).size() != nbasis ) {
	cerr << "Error writing vector: nbasis and actual are different " << endl;
	return( PV_WRITE_FAILED );
      }
      olocalfile.write( (char *)&(*outer)[0], (*outer).size()*sizeof(double) );
    }
  // Put titles at the end
  olocalfile.write(current_title.c_str(),current_title.size() );
  return( PV_SUCCESS );
}

//C++ approach I/O (streams) BINARY approach
/* New piece for a scalable Vector write
*/
int PsociVector::scalable_WriteHeaderBurstBinary(const pair<int64_t,int64_t> info )
{
  if (!olocalfile.is_open()) {
    cerr << "Header: Cannot open PsociVector output file " << current_filename << " Aborting write " << endl;
    return( PV_WRITE_FAILED );
  }
  
  // Lay down the values for nroots, nbasis, width, and prec
  int64_t infodata[2];
  int64_t nbasis = info.first;
  int64_t nroots = info.second;
  infodata[0] = nbasis;
  infodata[1] = nroots;
  olocalfile.write((char*)infodata, 2 * sizeof(int64_t) );
  
  int precdata[2];
  int width = precision.first;
  int prec  = precision.second;
  precdata[0] = width;
  precdata[1] = prec;
  olocalfile.write((char*)&precdata, 2 * sizeof(int) );
  
  return( PV_SUCCESS );
}

int PsociVector::scalable_WriteHeaderBurst(const pair<int64_t,int64_t> info )
{
  if (!olocalfile.is_open()) {
    cerr << "Cannot open PsociVector Scalable output file " << current_filename << " Aborting write " << endl;
    return( PV_WRITE_FAILED );
  }
  int64_t nbasis = info.first;
  int64_t nroots = info.second;
  int width = precision.first;
  int prec  = precision.second;

  olocalfile << nbasis << " " << nroots << endl;
  olocalfile << width << " " << prec << endl;

  return( PV_SUCCESS );
}

int PsociVector::scalable_ReadHeaderBurst(pair<int64_t,int64_t> & info )
{
  if (!olocalfile.is_open()) {
    cerr << "Cannot read PsociVector Scalable output file " << current_filename << " Aborting write " << endl;
    return( PV_WRITE_FAILED );
  }
  int64_t nroots, nbasis;
  int readwidth, readprec;

  olocalfile >> nbasis >> nroots;
  olocalfile >> readwidth >>  readprec ;

  precision.first = readwidth;
  precision.second = readprec;

  info.first = nbasis;
  info.second = nroots;

  return( PV_SUCCESS );
}

int PsociVector::scalable_ReadHeaderBurstBinary(pair<int64_t,int64_t> & info )
{
  if (!olocalfile.is_open()) {
    cerr << "Cannot read PsociVector Scalable output file " << current_filename << " Aborting write " << endl;
    return( PV_WRITE_FAILED );
  }
  int64_t infodata[2];
  olocalfile.read((char*)infodata, 2 * sizeof(int64_t) );
  info.first=infodata[0];
  info.second=infodata[1];

  int precdata[2];
  olocalfile.read((char*)&precdata, 2 * sizeof(int) );
  precision.first = precdata[0];
  precision.second = precdata[1]; 

  return( PV_SUCCESS );
}


/*
  C++ approach I/O (streams) BINARY approach
  In this case we need to remove the req of having an entire vector held in memory
  performance is not an issue here

  We use the enumeration state to decide which bits and pieces must get written to the file
  One constraint though, a ROW to be written must not span multiple invocatin as here we simply
  apend rows to the output file

  NOTE: we mostly likely arrive here with partial rows that need to be appended.

  Header and footer information get written by the calling app
*/

int PsociVector::WriteBurstCPPversionBinaryScalable(bool start, int root_num, vector<double> & vecs, const pair<int64_t,int64_t> info )
{
  if (!olocalfile.is_open()) {
    cerr << "Cannot open PsociVector Scalable output file " << current_filename << " Aborting write " << endl;
    return( PV_WRITE_FAILED );
  }
  if ( start == 1 ) {
     olocalfile.write( (char*)&root_num, sizeof(int64_t)*1);
  }
  olocalfile.write((char *)&(vecs)[0], (vecs).size()*sizeof(double) );
  return( PV_SUCCESS );
}

/* New piece for Scalable CI vevtor processing 
   Only writes one vector row at a time
*/

/*
int PsociVector::WriteBurstCPPversionScalable( int root_num, vector<double> & vecs, const pair<int64_t,int64_t> info )
{
  if (!olocalfile.is_open()) {
    cerr << "Write: Cannot open PsociVector output file " << current_filename << " Aborting write " << endl;    
    exit(1);
  } 
  cout << "WriteBurstCPPversionScalable" << endl;
  int64_t nbasis = info.first;
  int64_t nroots = info.second;
  int width = precision.first;
  int prec  = precision.second;
  
  olocalfile << setprecision(prec) << setiosflags(ios::fixed);

  if ( vecs.size() != nbasis ) { 
    cerr << "WriteBurstCPPScalable: Error writing vector: nbasis and actual are different " <<  vecs.size() << endl;
    return ( PV_WRITE_FAILED );
  }

  vector<double>::const_iterator it;
  olocalfile <<  root_num << endl;
    for (it=vecs.begin(); it != vecs.end(); ++it)
    {   
       olocalfile << setw(width) << (*it) ;
    } 

  return(PV_SUCCESS);
}
*/

int PsociVector::WriteBurstScalableTrailer()
{
  if (!olocalfile.is_open()) {
    cerr << "Write: Cannot open PsociVector output file  for trailer write " << current_filename << " Aborting write " << endl;
    exit(1);
  } 
  olocalfile << current_title << endl;
  return( PV_SUCCESS );
}

int PsociVector::WriteBurstScalableBinaryTrailer()
{
  if (!olocalfile.is_open()) {
    cerr << "Write: Cannot open PsociVector output file  for trailer write " << current_filename << " Aborting write " << endl;
    exit(1);
  }
  olocalfile.write(current_title.c_str(),current_title.size() );
  return( PV_SUCCESS );
}


/* 
  C++ approach I/O (streams) BINARY approach
  region is START, ONGOING, END
*/

/** A timer wrapped version of write burst: Use of a vector is 
    overdoing it for now, but I have plans.....For example an outer
    loop that drive multiple WriteBursts
*/

// No exit because a write may recoverable
double PsociVector::WriteBurst( vector<vector<double> > & vecs, const pair<int64_t,int64_t> info  )
{
  double timein = psociTime();
  if ( WriteBurstCPPversion( vecs, info ) != PV_SUCCESS ){
    cerr << "Failure to writeBurst " << endl;
  }
  return( psociTime() - timein );
}

double PsociVector::ReadBurst( vector<vector<double> > & vecs, pair<int64_t,int64_t> & info  )
{
  double timein = psociTime();
  if ( ReadBurstCPPversion( vecs, info ) != PV_SUCCESS ) {
    cerr << "Failure to readBurst " << endl;
    exit(1);
  }
  return( psociTime() - timein );
}

double PsociVector::ReadBurstScalable(bool start, int range, vector<double> & vecs, pair<int64_t,int64_t> & info  )
{
  double timein = psociTime();
  if ( ReadBurstCPPversionScalable( start, range, vecs, info ) != PV_SUCCESS ) {
    cerr << "Failure to readBurst " << endl;
    exit(1);
  }
  return( psociTime() - timein );
}

double PsociVector::ReadBurstScalableBinary(bool start, int range, vector<double> & vecs, pair<int64_t,int64_t> & info  )
{
  double timein = psociTime();
  if ( ReadBurstCPPversionScalableBinary( start, range, vecs, info ) != PV_SUCCESS ) {
    cerr << "Failure to readBurst " << endl;
    exit(1);
  }
  return( psociTime() - timein );
}

double PsociVector::WriteBurstBinary( vector<vector<double> > & vecs, const pair<int64_t,int64_t> info  )
{
  double timein = psociTime();
  if ( WriteBurstCPPversionBinary( vecs, info ) != PV_SUCCESS ){
    cerr << "Failure to writeBurst Binary " << endl;
  }
  return( psociTime() - timein );
}

double PsociVector::WriteBurstBinaryScalable(bool start, int rootnum,  vector<double> & vecs, const pair<int64_t,int64_t> info  )
{
  double timein = psociTime();
  if ( WriteBurstCPPversionBinaryScalable(start, rootnum, vecs, info ) != PV_SUCCESS ){
    cerr << "Failure to writeBurst Binary Scalable " << endl;
  }
  return( psociTime() - timein );
}

double PsociVector::WriteBurstScalable(bool start, int rootnum,  vector<double> & vecs, const pair<int64_t,int64_t> info  )
{
  double timein = psociTime();
  cout << " going to the Scalable ASCII  write burst " << endl;
  if ( WriteBurstCPPversionScalable(start, rootnum, vecs, info ) != PV_SUCCESS ){
    cerr << "Failure to writeBurst Binary Scalable " << endl;
  }
  return( psociTime() - timein );
}

double PsociVector::ReadBurstBinary( vector<vector<double> > & vecs, pair<int64_t,int64_t> & info  )
{
  double timein = psociTime();
  cout << " going to the binary read burst " << endl;
  if ( ReadBurstCPPversionBinary( vecs, info ) != PV_SUCCESS ) {
    cerr << "Failure to readBurst " << endl;
    exit(1);
  }
  return( psociTime() - timein );
}


// Used for double checking values: Not useful in production runs

double PsociVector::PrintBurst( vector<vector<double> > & vecs, const pair<int64_t,int64_t> info  )
{
  double timein = psociTime();
  PrintBurstCPPversion( vecs, info );
  return( psociTime() - timein );
}

void PsociVector::printWriteType()
{
  cout << "C++ IO style being used " << endl;
}
string PsociVector::getVersion() 
{
   return(PSOCIVECTOR_VERSION); 
}

/** Check the normalization ( unit metric ) for the subset of vectors
    read on this process. Use the read precision value as a guide
    for normality precision
*/

/* I didn't want to overcomplicate things so the DDOT is manual 
   and sqrt comes from cmath. THe actual application does a better
   level of checking anyway. Return all norms in norms 

   norm = euclidean norm values
*/

//TODO a decent Euclidean norm method
int PsociVector::checkNormalization( vector<vector<double> > & vecs, vector<double> & norms )
{
  double dONE = 1.0;
  int ALLNORMALIZED=1;
  vector<vector<double> >::const_iterator outer;
  int64_t root_num=0;
  for( outer = vecs.begin(); outer !=vecs.end(); ++outer)
  {
    ++root_num;
    double sum2 = 0.0;
    vector<double>::const_iterator it;
    for (it=(*outer).begin(); it != (*outer).end(); ++it)
      {
	sum2 += (*it)*(*it);
      }
    double norm = sqrt( sum2 );
    norms.push_back( norm );
    if ( (norm - dONE) > precision.second ) {
       ALLNORMALIZED = 0;
       }
  }
  return(ALLNORMALIZED);
}






