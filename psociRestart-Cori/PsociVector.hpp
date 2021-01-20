/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociVector:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociVector.hpp
 *
 */

#ifndef PSOCIVECTOR_H
#define PSOCIVECTOR_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include <sstream>
#include <fstream>

using namespace std;

const int MAXWIDTH = 18;
const int MAXDECIMAL = 14; // The more we have the better the restart

const int MAXPRINTCOEFS = 100;

const int PV_SUCCESS=101;
const int PV_OPEN_FAILED=102;
const int PV_WRITE_FAILED=103;
const int PV_READ_FAILED=104;

const ios_base::openmode VECTOR_WRITE = std::ios::out | std::ios::binary;
const ios_base::openmode VECTOR_READ  = std::ios::in | std::ios::binary;

const string PSOCIVECTOR_VERSION = "1.0";

class PsociVector {
  
private:
  int64_t blk_size;
  string current_filename;
  string current_title;
  int set_number;
  pair<int,int> precision;
  ios_base::openmode current_mode;
 
  int WriteBurstCPPversion( vector<vector<double> > & vecs, const pair<int64_t,int64_t> info );
  void PrintBurstCPPversion( vector<vector<double> > & vecs, const pair<int64_t,int64_t> info );
  int ReadBurstCPPversion( vector<vector<double> > & vecs, pair<int64_t,int64_t> & info );
  int WriteBurstCPPversionBinary( vector<vector<double> > & vecs, const pair<int64_t,int64_t> info );
  int ReadBurstCPPversionBinary( vector<vector<double> > & vecs, pair<int64_t,int64_t> & info );
  //int WriteBurstCPPversionBinaryScalable(vector<double> & vecs, const pair<int64_t,int64_t> info );
  //int WriteBurstCPPversionScalable(bool start, int rootnum,  vector<double> & vecs, const pair<int64_t,int64_t> info );
  
  int OpenFileSeq( void );
  void CloseFileSeq( void );
  
  void Filename( const int num, const string filname );
  void setTitle( const string title );
  void printTitle();

  fstream olocalfile;
  
public:
  PsociVector();
  PsociVector( const int num, const string filename );
  PsociVector( const string filename );
  PsociVector( const int num, const string filename, const pair<int,int> info );
  void Params( void );
  void setMode( ios_base::openmode MODE );
  ios_base::openmode fetchMode();
  string FetchFilename();
  double WriteBurst( vector<vector<double> > & vecs, const pair<int64_t,int64_t> info );
  double WriteBurstBinary( vector<vector<double> > & vecs, const pair<int64_t,int64_t> info );
  void printWriteType();
  double PrintBurst( vector<vector<double> > & vecs, const pair<int64_t,int64_t> info );
  double ReadBurst( vector<vector<double> > & vecs, pair<int64_t,int64_t> & info );
  double ReadBurstBinary( vector<vector<double> > & vecs, pair<int64_t,int64_t> & info );
  double WriteBurstBinaryScalable(bool start, int rootnum,  vector<double> & vecs, const pair<int64_t,int64_t> info  );
  double WriteBurstScalable(bool start, int rootnum,  vector<double> & vecs, const pair<int64_t,int64_t> info  );
  int OpenFile();
  void CloseFile();
  string getVersion( void );
  int checkNormalization( vector<vector<double> > & vecs, vector<double> & norms );
  void setFilename( const int num, const string filename );
  void setFilename( const int num, const string filename, const pair<int,int> prec );
/*
  int scalable_WriteHeaderBurstBinary(const pair<int64_t,int64_t> info );
  int scalable_WriteHeaderBurst(const pair<int64_t,int64_t> info );
  int WriteBurstCPPversionBinaryScalable(bool start, int root_num, vector<double> & vecs, const pair<int64_t,int64_t> info );
  int WriteBurstCPPversionScalable(bool start, int root_num, vector<double> & vecs, const pair<int64_t,int64_t> info );
  int WriteBurstScalableTrailer();
  int WriteBurstScalableBinaryTrailer();
  int WriteBurstScalableEOL();
  double ReadBurstScalable(bool start, int range, vector<double> & vecs, pair<int64_t,int64_t> & info  );
  double ReadBurstScalableBinary(bool start, int range, vector<double> & vecs, pair<int64_t,int64_t> & info  );
  int ReadBurstCPPversionScalable(bool start, int range, vector<double> & vecs, pair<int64_t,int64_t> & info );
  int ReadBurstCPPversionScalableBinary(bool start, int range, vector<double> & vecs, pair<int64_t,int64_t> & info );
*/
  int ReadBurstCPPversionScalable(bool start, int range, vector<double> & vecs, pair<int64_t,int64_t> & info );
  int ReadBurstCPPversionScalableBinary(bool start, int range, vector<double> & vecs, pair<int64_t,int64_t> & info );
  int WriteBurstCPPversionBinaryScalable(bool start, int root_num, vector<double> & vecs, const pair<int64_t,int64_t> info );
  int WriteBurstCPPversionScalable(bool start, int root_num, vector<double> & vecs, const pair<int64_t,int64_t> info );
  int WriteBurstScalableTrailer();
  int WriteBurstScalableBinaryTrailer();
  int WriteBurstScalableEOL();
  int scalable_WriteHeaderBurstBinary(const pair<int64_t,int64_t> info );
  int scalable_WriteHeaderBurst(const pair<int64_t,int64_t> info );
  int scalable_ReadHeaderBurstBinary(pair<int64_t,int64_t> & info );
  int scalable_ReadHeaderBurst(pair<int64_t,int64_t> & info );
  double ReadBurstScalable(bool start, int range, vector<double> & vecs, pair<int64_t,int64_t> & info  );
  double ReadBurstScalableBinary(bool start, int range, vector<double> & vecs, pair<int64_t,int64_t> & info  );

};

#endif //PSOCIVECTOR_H
