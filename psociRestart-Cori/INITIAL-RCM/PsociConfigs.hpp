/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociConfigs:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociConfigs.hpp
 *
 */

#ifndef PSOCICONFIGS_H
#define PSOCICONFIGS_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>

#include "ga++.h"

using namespace std;

const string DEFAULT_FILENAME="configs";
const int MAX_CONFIG_PRINT = 1000; 

class PsociConfigs {
  
private:

  GA::GlobalArray * g_configs; // New stuff for readConfigsGA and related methods

  int n_openshells; //supposed to be the maximum number of open shells in the data - we should check this however.
  int n_basisftns;
  int n_electrons;
  int st_symmetry;
  int n_configurations;
  void OpenFile();
  void CloseFile();
  int spin_parity; // SPIN is even or odd: applies to all configs. = nelec % 2 (ODD implies E symmetry )
  
  ifstream ilocalfile;
  
  string local_title;
  string local_filename;
  
  multimap<int,string> configs;

//  vector<int> configs_index; // only for new GA-based methods
//  multimap<int,int> configs_index;

/* Special object for distributed configs
   since duplicate (isym) order is not guarenteed we cannot simply shove indexes into another
   object. They must be tethered to the configuration
    
*/

  vector<pair<int,string> > vec_configs;// Yikes we want this to more easily do random acces in fetchConfigsAlt
  vector<int> vec_configs_index;
 
  
public:
  PsociConfigs();
  PsociConfigs( string filename );
  void printFilename();
  int readConfigs();
  int readConfigs( pair<int,double> & info );
  int readConfigsToGA();
  int readConfigsToGA( pair<int,double> & info );
  void printParams();
  void fetchParams( vector<int> & params );
  void setParams( vector<int> & params );
  void printConfigurations();
  void printDistributedConfigurations();
  void brdcstConfigs( int root );
  void brdcstConfigs( int root, pair<int,double> & info );
  void brdcstParams( int root );
  void convertMultimapToVector( vector<int> & isym, vector<short> & confs, vector<long> & sizes );
  void convertVectorToMultimap( vector<int> & isym, vector<short> & confs );
  void convertVectorToVectormap( vector<int> & isym, vector<short> & confs );
  void validateConfigurationSymmetry();
  int fetchConfig( int start_iconf, int stop_iconf , vector<pair<int,string> > & selcfgs );
  int fetchConfigAlt( int start_iconf, int stop_iconf, vector<pair<int,string> > & selcfgs );
  int numTotalSpatials();
  int numBasisFunctions();
  void tearDownConfigs();
  int fetchSpinParity();
  int fetchMaxOpenShells();
  int fetchStateSymmetry();
  int fetchNumBasisFtns();
  int fetchNumElectrons();
  void copyInternalMapToVector();
  void convertStringToShortArray( string word, vector<short> & vector );
  void convertShortArrayToString( vector<short> & vector, string & word );
  void convertStringToIntArray( string word, vector<int> & vector );
  void convertStringToIntArray( vector<char> & word, vector<int> & vector );
  void convertIntArrayToString( vector<int> & vector, string & word );
  void shuffleConfigurations(int root, vector<int> & isym, vector<short> & configs );

/* new GA methods */
  void createGAconfigs( int num_basis, int num_configs );
  void destroyGAconfigs(); 
  void pushConfigsToGA( int indexGA, int nbf, int isym, vector<int> & configs );
  int getLocalConfigsFromGA( );
  int fetchConfigGA( vector<pair<int,string> > & selcfgs, vector<int> & index );
  int localConfigsSize();
  int localVecConfigsSize();

};

#endif //PSOCICONFIGS_H
