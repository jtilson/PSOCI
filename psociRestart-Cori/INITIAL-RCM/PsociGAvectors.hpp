/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociGAvectors:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociGAvectors.hpp
 *
 */

#ifndef PSOCIGAVECTORS_H
#define PSOCIGAVECTORS_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <ga++.h>
#include <GAServices.h>

using namespace std;

const int MAX_VECTOR_DIMS = 2;

const int MAX_VECTOR_RESTART_WIDTH = 100000; // Preclude reading/writing too much data to/from CI restart vec file

class PsociGAvectors {
  
private:
  int64_t test;
  int numReadRoots;
  int numSoughtRoots;

  PsociVector vectorRestart; 

public:
  explicit PsociGAvectors( int nrootsSought, GA::GlobalArray * g_array );
  void closeVectorFile(void);
  void openVectorWrite( int rank, string filename );
  void openVectorRead( int rank, string filename );
  void dumpVector(int rank, GA::GlobalArray * g_vector );
  void dumpVector(int rank, int num, GA::GlobalArray * g_vector );
  void readVector(int rank, GA::GlobalArray * g_vector );
  void dumpVector(int rank, GA::GlobalArray * g_vector, double &time );
  void dumpVector(int rank, int num, GA::GlobalArray * g_vector, double &time );
  void dumpVectorScalable( int rank, int num, GA::GlobalArray * g_vector );
  void readVector(int rank, GA::GlobalArray * g_vector, double &time );
  int getNumberVectors( void );
  int getNumberVectorsSought( void );
  void setNumberVectors( int resetVectorNum );
  void specifyDiagonalGuessVectors( GA::GlobalArray * g_array );
  void printCurrentSoughtRoots( GA::GlobalArray * g_array ) ;
  void readVectorScalable( int rank, GA::GlobalArray * g_vector );

};

#endif //PSOCIGAVECTORS_H
