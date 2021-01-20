/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociGArestart:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociGArestart.hpp
 *
 */

#ifndef PSOCIGARESTART_H
#define PSOCIGARESTART_H

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

class PsociGArestart {
  
private:
  int64_t test;
  int numReadRoots;

  PsociVector vectorRestart; 

public:
  PsociGArestart();
  PsociGArestart( int num, GA::GlobalArray * g_array );
  void closeVectorFile(void);
  void openVectorWrite( int rank, string filename );
  void openVectorRead( int rank, string filename );
  void dumpVector(int rank, GA::GlobalArray * g_vector );
  void readVector(int rank, GA::GlobalArray * g_vector );
  void dumpVector(int rank, GA::GlobalArray * g_vector, double &time );
  void readVector(int rank, GA::GlobalArray * g_vector, double &time );
  int getNumberVectors( void );
  void setNumberVectors( int resetVectorNum );
  void specifyDiagonalGuessVectors( int numroots, GA::GlobalArray * g_array );
  
};

#endif //PSOCIGARESTART_H
