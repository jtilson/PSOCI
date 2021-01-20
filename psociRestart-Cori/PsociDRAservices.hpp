/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociGAservices:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociDRAservices.hpp
 *
 */

#ifndef PSOCIDRASERVICE_H
#define PSOCIDRASERVICE_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <ga++.h>

extern "C" {
#include <dra.h>
#include <macdecls.h>
}

using namespace std;

#define DRA_MAX_ARRAY 10

//const int DEF_PSOCI_MAX_ARRAYS = 4;
const int DEF_PSOCI_MAX_ARRAYS = 7; // need seven if using SUPPLEMENTAL
const double DEF_PSOCI_MAX_ARRAYSIZE = -1;
const double DEF_PSOCI_TOTALDISKSPACE = -1;
const double DEF_PSOCI_MAX_MEMORY = -1;

//const int MAX_DRA_DIMS = 2;

class PsociDRAservices {
  
private:
  int64_t test;
  int PSOCI_MAX_ARRAYS;
  double PSOCI_MAX_ARRAYSIZE;
  double PSOCI_TOTALDISKSPACE;
  double PSOCI_MAX_MEMORY;
  int PSOCI_NUM_DRA_IO_FILES;
  int PSOCI_NUM_DRA_IO_SERVERS;

  int DRA_status;
  

public:
  PsociDRAservices();
  void reportDRAParams( void );
  void setDRAParams( int max_arrays );
  void setDRAParams( int max_arrays, int max_array_size, int total_disk, int max_memory );
  void numIOservers( int numIOfiles, int numIOservers );
  void reportDRAIOservers( void );
  void initDRA( int rank );
  void initDRA( int rank, int servers );
  void exitDRA();
};

#endif //PSOCIDRASERVICE_H
