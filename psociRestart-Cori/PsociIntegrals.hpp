/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociIntegrals:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociIntegrals.hpp
 *
 */

#ifndef PSOCIINTEGRALS_H
#define PSOCIINTEGRALS_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>

#include "ga++.h"

//Fortran conversions

#ifdef FORT64
typedef int64_t Fint;
#else
typedef int Fint; //simply grab local size
#endif

typedef double Fdouble;

const int MAX_TITLE_LENGTH = 80;
const int MAX_LABEL_LENGTH = 4; 
const int MAX_BFNLAB_LENGTH = 8;
const int MAX_TYPE_LENGTH = 8;
const int MAX_ONEINT_HEADING = 32;
const int TENMILLION = 10 * 1000 * 1000;
const int MAX_TWOE_PRINTSIZE = TENMILLION;;

const double MIN_DOUBLE_PRINTED = 1.0e-8;

//Columbus SIFS calls need to open/close the integral file.
//void ftocxxstringlist_(string & list,const int & size, const int & length);

using namespace std;

extern "C" {
  void ftocxxopen_(Fint &unit, const char * filename, Fint &namesize,Fint &istatus );
  void ftocxxclose_(Fint &unit);
  void sifrh1_(Fint &unit, Fint &ntitle, Fint &nsym, Fint &nbft,
               Fint &ninfo, Fint &nenrgy, Fint &nmap, Fint &ierr);
  
  void sifrh2_( Fint &unit, Fint &ntitle, Fint &nsym, Fint &nbft, 
		Fint &ninfo, Fint &nenrgy, Fint &nmap, char * title,  
		Fint * nbpsy, char * slabels, Fint * info, char * bfnlabs, 
		Fint * ietype, Fdouble * energy, Fint * imtype, Fint * map, //Verify settings for MAP
		Fint &ierr );

  void sifrd1_( Fint &unit, Fint * info, Fint & nipv, Fint & iretbv, 
                Fdouble * buffer, Fint & num, Fint &  last, Fint &  itypea,
                Fint &  itypeb, Fint & ifmt, Fint & ibvtyp, 
                Fdouble * values, Fint * labels, Fdouble & fcore, Fint * ibitv, Fint & ierr ); 

  void sifrd2_( Fint &unit, Fint * info, Fint & nipv, Fint & iretbv, 
                Fdouble * buffer, Fint & num, Fint &  last, Fint &  itypea,
                Fint &  itypeb, Fint & ifmt, Fint & ibvtyp, 
                Fdouble * values, Fint * labels, Fint * ibitv, Fint & ierr ); 

  void siftyp_( Fint & itypea, Fint & itypeb, char chrtyp[MAX_TYPE_LENGTH] );

/*
  void sifrsh_( Fint & aoints, Fint * info, Fdouble * buffer, Fdouble * values,
                Fint * labels, Fint & nsym, Fint * nbpsy, Fint * map,
                Fint & nnbft, Fdouble * s1h1, Fdouble & score, Fdouble & hcore, 
                Fint * symb, Fint &ierr );
*/


/*
        call sifrd2(moint2,info,nipv,iretbv,buff2,num,last,itypea,
     &    itypeb,ifmt,ibvtyp,value2,label2,idummy,ierr)
*/

}

/*
   char title[20][80];
   Fint nbpsy[4];
   char slabel[4][4];
   int info[5];
   char bfnlab[86][8];
   int ietype[1];
   double energy[1];
   int imtype[86];
   int map[86][0];
*/


const int DEFAULT_MOINTS_UNIT = 2;  
const string DEFAULT_MOINTS_FILENAME = "moints";

class PsociIntegrals {
  
private:

  int ROOT_READER; // process handling the read of the integrals
  string l_moints; // Name of the current MO integral file
  Fint l_unit;      // Numerical unit number for MOINTS used by Fortran routines.
  Fint l_ntitle, l_nsym, l_nbft, l_ninfo, l_nenrgy, l_nmap; //header 1 info

  vector<string> l_title; //Therse get passed back/forth to sifs routines and so have 'F' size data types
  vector<string> l_slabel;
  vector<string> l_bfnlab;
  vector<Fint> l_nbpsy;
  vector<Fint> l_info;
  vector<Fint> l_ietype; 
  vector<vector<Fint> > l_map;
  vector<Fint> l_imtype; 
  vector<Fdouble> l_energy;

  vector<Fdouble> l_valone; // Filled out 1e integral array 
  vector<Fdouble> l_valtwo;
  vector<Fdouble> l_s1h1; // For new overlap methods


  bool got_1e;
  bool got_2e;
 
  bool perform_2e_print; // These are stopgaps that preclude dumping lots of data even when calling print...
  bool perform_1e_print;

  int TwoEIndex(int i, int j); // Call using Cstyle indexing
  int canonicalIndexF77(int i,int j,int k,int el); // Call; for indexes/labels starting at 1
  int canonicalIndex(int i,int j,int k,int el); // Call using indexes that start at 0


public:
  PsociIntegrals( int reader, Fint unit, string filename );
  PsociIntegrals();
  void set1ePrint();
  void set2ePrint();
  int fetchUnit();
  void printFilename();
  void OpenFile();
  void printSifrh2();
  void CloseFile();
  int sifrh1();
  void printSifrh1();
  int sifrh2();
  int fetchOneElectronInts();
  int fetchAndProcessIntegrals();
  int fetchTwoElectronInts();
  int fetchOverlapInts();
  int fetchLabels( vector<string> & labels );
  int fetchBasisLabels( vector<string> & labels );
  void printOneElectronInts();
  void printTwoElectronInts();
  void brdcstParams();
  void brdcstSifrh2();
  int brdcstTwoElectronInts();
  int brdcstOneElectronInts();
  double fetch2eIntegralValue(int i,int j,int k,int el);
  double fetch1eIntegralValue(int i,int j);
  int fetchAndProcessOverlapIntegrals();
  double fetchCoreEnergy();

};

#endif //PSOCIINTEGRALS_H
