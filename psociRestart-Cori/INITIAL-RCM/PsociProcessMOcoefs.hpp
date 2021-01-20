/**************************************************************************************

* Copyright (c) 2010,2011,2012 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociProcessMOcoefs:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociProcessMOcoefs.hpp
 *
 */

#ifndef PSOCIPROCESSMOCOEFS_H
#define PSOCIPROCESSMOCOEFS_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <vector>

#include "ga++.h"

//Fortran conversions

#ifdef FORT64
typedef int64_t Fint;
#else
typedef int Fint; //simply grab local size
#endif
typedef double Fdouble;

using namespace std;

extern "C" {
  void ftocxxopen_(Fint &unit, const char * filename, Fint &namesize, Fint &istatus );
  void ftocxxformattedopen_(Fint &unit, const char * filename, Fint &namesize, Fint &istatus );
  void ftocxxformattednewopen_(Fint &unit, const char * filename, Fint &namesize, Fint &istatus );
  void ftocxxclose_(Fint &unit);
  void moread_( Fint &unit, Fint & flcod, Fint & fierr, Fint & syserr, 
                Fint & maxtitle, Fint & ntitle, char * titles, char * afmt, 
                Fint & nsym, Fint * nbfpsy, Fint * nmopsy, char * labels, Fint & clen, 
                Fdouble * coefs ); 

/* We resort to the C-Fortran added information; array sizes, here. (stitle, saftm)
   primarily becasue in mowrit aftm is a character constant. Thus we MUST pass in 
   the size of the array. The stitle is their simply to keep order.

   All character I/O (except afmt) are assumable arrays such as character*80 (*)
   For these we can simply assume their size and processs accordingly. This works fine.

   afmr, however is a character constant character *(*). THis is NOT an assumable
   array and ordinarily the Fortran compiler would write code to properly assign this memory
   base don the calling app. Calling, from C++, however, precludes that. SO we MUST provide
   size information in an alternative way.
*/

  void mowrit_( Fint & unit, Fint & flcod, Fint & fierr, Fint & syserr,
                Fint & ntitle, char * titles, char * afmt,
                Fint & nsym, Fint * nbfpsy, Fint * nmopsy, char * labels,
                Fdouble * coefs, int stitle, int saftm );
}

const string MO_DEFAULT_FILENAME="mocoef";
const Fint MO_DEFAULT_MOCOEF_UNIT=5;
const int MO_DEFAULT_ROOT_READER=0;

const int MO_MAX_TITLES=10;
const int MO_MAX_TITLE_LENGTH = 80;
const int MO_MAX_FORMAT_LENGTH = 80;
const int MO_MAX_LABEL_LENGTH = 4; 
const int MO_MAX_LABELS = 8;

class PsociProcessMOcoefs {

private:

 vector<string> l_title;
 vector<string> l_labels;
 string l_afmt;
 vector<Fint> l_nbfpsy;
 vector<Fint> l_nmopsy;
 int READER;
 
 vector<Fdouble> l_coefs;

 string l_filename; // Intended for MOs
 vector<string> l_list_filenames; // Used for NOs
 int l_unit, l_nsym, l_ntitle;
 int l_base_unit;
 vector<int> l_list_units;

 int  l_nbf;
 int  l_nbfsq; // coeff array can be no bigger than nbf*nbf.
 int l_numroots;


public:

   PsociProcessMOcoefs(int nbf) ;
   PsociProcessMOcoefs( int nbf,int unit );
   PsociProcessMOcoefs( int nbf, int unit, string filename );
   PsociProcessMOcoefs( int nbf, int unit, int numroots, string filename ); // Intended for using multiple files for NOs


 void OpenMOfile( string filename );
 void OpenMOfile( );
 void OpenNOfiles();
 void OpenNOfiles( int iroot );
 void CloseMOfile( );
 void CloseNOfiles();
 void CloseNOfiles( int iroot );
 void WriteNOheaders( string title );
 int WriteNOcoefs(string title, int iroot, vector<Fdouble> & evals, vector<Fdouble> & coefs );
 int MOcoefficientRead( void );
 void printMOparameters( void);
 void printMOcoefficients( void );
 vector<Fdouble> * fetchCoefsHandle();
 int fetchNsym();
 int fetch_nmopsy(int isym);
 int fetch_nbfpsy(int isym);
 double fetchMocoefs( int index );

 int copyOutHeaderInfo( int & nsym, vector<string> & title, vector<string> & labels, string & afmt, vector<Fint> & nbfpsy, vector<Fint> & nmopsy );
 int copyInHeaderInfo( int & nsym, int & ntitle, vector<string> & title, vector<string> & labels, string & afmt, vector<Fint> & nbfpsy, vector<Fint> & nmopsy );

};

#endif //PSOCIPROCESSMOCOEFS_H
