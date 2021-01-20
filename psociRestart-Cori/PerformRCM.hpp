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
 *   @file PerformRCM.hpp 
 *
 */

#ifndef PSOCIPERFORMRCM_H
#define PSOCIPERFORMRCM_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <ga++.h>

extern "C" {
#include <dra.h>
#include <macdecls.h>
}

using namespace std;

typedef int DRA;

const int RCM_MAX_DRA_DIMS = 5;
const int RCM_DEF_DRA_MODE = DRA_RW;
const int RCM_DEF_DRA_TYPE = MT_DBL;

#define DRA_TRUE (logical)1
#define DRA_FALSE (logical)0

class PerformRCM {

private:
  int64_t test;
  int DRA_status;
  DRA local_dra;
  int l_ndim;
  int l_mode; 
  int l_dim0,l_dim1,l_rdim0,l_rdim1;
  string l_arrayname, l_filename;
  int l_type;

  int l_req; // Request handle 

  GA::GlobalArray * g_icicol;
  GA::GlobalArray * l_g_icicol_supp;


public:
  PerformRCM();
  PerformRCM( int rank, string arrayname, string filename );
  void reportName();
  void reportFilename();
  void setType( int type );
  void printType();
  void setMode( int mode );
  void setDimensions( int ndim, pair<dra_size_t,dra_size_t> dims, pair<dra_size_t,dra_size_t> rdims );
  void getDimensions( int &ndim, pair<dra_size_t,dra_size_t> &dims, pair<dra_size_t,dra_size_t> &rdims );
  int createDRA();
  void printInternals();
  int dumpGASectionDRA( vector<pair<int,int> > ga_patch, vector<pair<dra_size_t,dra_size_t> > dra_patch, GA::GlobalArray * g_vector );
  int readGASectionDRA( vector<pair<int,int> > ga_patch, vector< pair<dra_size_t,dra_size_t> > dra_patch,  GA::GlobalArray * g_matrix);
  int dumpGAtoDRA( GA::GlobalArray * g_matrix );
  int readGAfromDRA( GA::GlobalArray * g_matrix );
  int waitDRA();
  int openDRAforREAD();
  int closeDRA();
  int deleteDRA();
  void flickDRA();
  int inquireDRAsetREAD()
;

};

#endif //PSOCIPERFORMRCM_H
