/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociHamiltonian.hpp
 *
 */

#ifndef PSOCIHAMILTONIAN_H
#define PSOCIHAMILTONIAN_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "ga++.h"

#include "PsociDeterminants.hpp"
#include "PsociIntegrals.hpp"

using namespace std;

const double MIN_DOUBLE_VALUE = 1.0e-8;  //packSefsef

class PsociHamiltonian {
  
private:
  int state;
  int l_nbf;
  int l_nelec;
  int l_maxdet; // per spatial
  int l_maxsef; // per spatial
  int l_ksym;
  
  PsociIntegrals * l_integrals;


  int processSingleExcitation( JOUTFG & iconf, JOUTFG & jconf , vector<int> & iocc, vector<int> & jocc, vector<int> & iex, vector<int> & jex, vector<double> & vnt, vector<double> & detdet  );
  int processDoubleExcitation( vector<int> & iex, vector<int> & jex, vector<double> & vnt );
  int processDiagExcitation( JOUTFG & iconf, vector<int> & iocc, vector<double> & vnt );
  int socig(JOUTFG & iconf, JOUTFG & jconf, vector<int> & iex, vector<int> & jex, vector<double> & vnt, vector<double> & detdet );
  void print2Darray( int nrows, int ncols, vector<double> & detdet );
  void transformCDCt( int ncpi, int ncpj, int nopeni, int nopenj, int nsefi, int ndeti, int nsefj, int ndetj, int & rows, int & cols, vector<double>  & sefsef, vector<pair<int,int> > & coef, vector<double>  & detdet);
  void transformCDCt( int ncpi, int ncpj, int nopeni, int nopenj, int nsefi, int ndeti, int nsefj, int ndetj, int & rows, int & cols, vector<double>  & sefsef, vector<pair<int,int> > & coef, vector<double>  & detdet, pair<int,double> & info);
  int generateSEFCIG(JOUTFG & iconf, JOUTFG & jconf, vector<int> & iex, vector<int> & jex, vector<double> & vnt, vector<double> & detdet );
  int generateSEFCIGCIDEN(JOUTFG & iconf, JOUTFG & jconf, vector<int> & iex, vector<int> & jex, vector<double> & detdet );


public:
  PsociHamiltonian( int nbf, int nelec );
  void specifyIntegrals( PsociIntegrals * ints );
  int generateOffdiagBlocks( JOUTFG & iconf, JOUTFG & jconf, vector<int> & iocc, vector<int> & jocc, vector<int> & iex, vector<int> & jex, vector<double> & vnt, vector<double> & detdet  );
  int generateOffdiagBlocks( JOUTFG & iconf, JOUTFG & jconf, vector<int> & iocc, vector<int> & jocc, vector<int> & iex, vector<int> & jex, vector<double> & vnt, vector<double> & detdet , pair<int,double> & info );
  int generateDiagonalBlocks( JOUTFG & iconf, vector<int> & iocc, vector<double> & vnt );
  int generateDiagonalBlocks(JOUTFG & iconf, JOUTFG & jconf, vector<int> & iocc, vector<double> & vnt );
  int generateDiagonalBlocks( JOUTFG & iconf, vector<int> & iocc, vector<double> & vnt, pair<int,double> & info );
  int generateSOCIblocks(JOUTFG & iconf, JOUTFG & jconf, vector<double> & sefsef );
  int generateSOCIblocksForCIDEN(JOUTFG & iconf, JOUTFG & jconf, vector<double> & sefsef );
  int determinantToDoublegroupMaps( vector<pair<int,int> > & sumcoef );

  void setPerSpatialValues( int maxdet, int maxsef );
  void printPerSpatialValues();
  void printKsym();
  void setKsym( int ksym );
  void ciout( JOUTFG & iconf, JOUTFG & jconf, vector<double> & detdet, vector<double> & sefsef );

  
};

#endif //PSOCIHAMILTONIAN_H
