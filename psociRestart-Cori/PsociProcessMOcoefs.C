//TODO Use of Fint and Fdouble is clumsy and needs to be cleaned up
/**************************************************************************************

* Copyright (c) 2010,2011,2012 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: A class to process MO coefficients for optimized MOs/NOs
              Currently uses MOREAD methods from Columbus
 History:

 

**************************************************************************************/
/**
 *   @file PsociProcessMOcoefs.C
 *
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include "PsociProcessMOcoefs.hpp"

using namespace std;

PsociProcessMOcoefs::PsociProcessMOcoefs( int nbf )
{
   READER = MO_DEFAULT_ROOT_READER;
   l_unit = MO_DEFAULT_MOCOEF_UNIT;
   l_filename =  MO_DEFAULT_FILENAME;
   l_nbf = nbf;
   l_nbfsq = nbf * nbf;
   l_numroots = 0;
}

PsociProcessMOcoefs::PsociProcessMOcoefs( int nbf, int unit )
{
   READER = MO_DEFAULT_ROOT_READER;
   l_filename =  MO_DEFAULT_FILENAME;
   l_unit = unit;
   l_nbf = nbf;
   l_nbfsq = nbf * nbf;
   l_numroots = 0;
}

PsociProcessMOcoefs::PsociProcessMOcoefs( int nbf, int unit, string filename )
{
   READER = MO_DEFAULT_ROOT_READER;
   l_filename =  filename;
   l_unit = unit;
   l_nbf = nbf;
   l_nbfsq = nbf * nbf;
   l_numroots = 0;
}

PsociProcessMOcoefs::PsociProcessMOcoefs( int nbf, int unit, int numroots, string filename )
{
   READER = MO_DEFAULT_ROOT_READER;
   l_filename =  filename;
   l_base_unit = unit; // Special case: we will construct a RANGE of units (base, base+numroots);
   l_nbf = nbf;
   l_nbfsq = nbf * nbf;
   l_numroots = numroots;
}


void PsociProcessMOcoefs::OpenMOfile()
{
  if ( GA::nodeid() == READER ) {
    Fint size = l_filename.size();
    Fint unit = l_unit;
    Fint istatus;
    ftocxxformattedopen_( unit, l_filename.c_str(), size, istatus);
    if ( istatus != 0 ) GA::error(" Failed to open MO file ",unit );
  }
}

void PsociProcessMOcoefs::OpenNOfiles() // Take the base filename and process for one file per root
{
  if ( GA::nodeid() == READER ) {
    cout << "Open one NO file for each CI root " << l_numroots << endl;

// Build an array of NO filenames (l_list_filenames )
// Need an associated list of UNITs as well: start with l_base_unit

  int temp_unit = l_base_unit;
  for(int i=0; i< l_numroots; ++i ) {
    string tempString = l_filename;
    string scount;
    stringstream sout;
  sout << i;
  scount = sout.str();
  scount.insert(0,"-");
  int locate = tempString.find_last_of(".");
  if ( locate != -1 ) {
        tempString.insert( locate, scount );
  } else {
        tempString.append( scount );
  }
  string amended_filename = tempString;
  l_list_units.push_back( temp_unit );
  l_list_filenames.push_back( amended_filename );
  ++temp_unit;
  }

#ifdef DETAILEDCHECK
  cout << "Size of NO list of fienames is " << l_list_filenames.size();
  cout << "Size of UNITS list is " << l_list_units.size() << endl;
#endif

  for(int i=0; i< l_numroots; ++i ) {
    Fint size = l_list_filenames[i].size();
    Fint unit = l_list_units[i];
    Fint istatus;
    ftocxxformattednewopen_( unit, l_list_filenames[i].c_str(), size, istatus);
    if (istatus != 0 ) GA::error(" Failed to open NO file",unit);
  }
}
}

/* Alternative: open one file ata time 
*/
void PsociProcessMOcoefs::OpenNOfiles(int iroot ) // Take the base filename and process for one file per root
{
  if ( GA::nodeid() == READER ) {
    cout << "Open one NO file for CI root " << iroot << endl;

// Build an array of NO filenames (l_list_filenames )
// Need an associated list of UNITs as well: start with l_base_unit

    string tempString = l_filename;
    string scount;
    stringstream sout;
  sout << iroot;
  scount = sout.str();
  scount.insert(0,"-");
  int locate = tempString.find_last_of(".");
  if ( locate != -1 ) {
        tempString.insert( locate, scount );
  } else {
        tempString.append( scount );
  }
  string amended_filename = tempString;

#ifdef DETAILEDCHECK
  cout << "Size of NO list of fienames is " << l_list_filenames.size();
  cout << "Size of UNITS list is " << l_list_units.size() << endl;
#endif

    Fint size = amended_filename.size();
    Fint unit = l_base_unit;
    Fint istatus;
    ftocxxformattednewopen_( unit, amended_filename.c_str(), size, istatus);
    if (istatus != 0 ) GA::error(" Failed to open alt-NO file",unit);
}
}

/* This method basically replicates the header information from the already-opened
   MOCOEF file and replicates that to all the NO files.  These data are all the same (nbf, nsym, etc.)
   The only different will be the actual coefficients and the inclusion of the occupations

   Add a customed additional title to the list of titles
*/
void PsociProcessMOcoefs::WriteNOheaders( string s_title ) 
{
  if ( GA::nodeid() != READER ) return;

  Fint nsym = 0;

// Scratch Columbus requires strickly conforming array space

  nsym = l_nsym;

  vector<char> title( MO_MAX_TITLES * MO_MAX_TITLE_LENGTH,' ');  // 10 max titles at max 80 chars each
  vector<char> afmt( MO_MAX_FORMAT_LENGTH,' ' );                  // max 80 chars
  vector<char> labels( MO_MAX_LABEL_LENGTH * MO_MAX_LABELS,' ' ); // 8 max labels at max 4 chars e

  vector<Fint> nbfpsy( MO_MAX_LABELS,0 ); // 8 maximum labels
  vector<Fint> nmopsy( MO_MAX_LABELS,0 );

  for(int i=0; i< nsym; ++i) {
    nbfpsy[i] = l_nbfpsy[i];
  }

  for(int i=0; i< nsym; ++i) {
    nmopsy[i] = l_nmopsy[i];
  }

  for(int j=0; j< MO_MAX_FORMAT_LENGTH; ++j ) {
     afmt[ j ] = l_afmt[j];
  }

  for(int i=0; i< nsym; ++i) {
    string word = l_labels[i];
    for(int j=0; j< word.size(); ++j ) {
      labels[ MO_MAX_LABEL_LENGTH * i + j ] = word[ j ]; 
    }
  }

// Treat the titles INSIDE the loop over roots for simplicity

 for( int iroot=0; iroot< l_numroots; ++iroot ) { // process headers for all files.
   int ntitle = l_ntitle;

   vector<string> temp_title = l_title;
   string tempString = s_title;
   string scount;
   stringstream sout;
   sout << iroot;
   scount = sout.str();
   scount.insert(0,"-"); // preceding dash
   tempString.append( scount );
   cout << " NEW TITLE " << tempString <<  endl;

   if ( ntitle < MO_MAX_TITLES ) {
           temp_title.push_back( tempString );
           ++ntitle;
   } else { //overwrite last entry
           temp_title[ MO_MAX_TITLES -1 ] =  tempString;
   }

  for(int i=0; i< ntitle; ++i) {
    string word = temp_title[i];
    for(int j=0; j< word.size(); ++j ) { 
       title[ MO_MAX_TITLE_LENGTH * i +  j ] = word[ j ];
    }   
  }

// Dump to nocoef-<n>

   vector<Fdouble> l_coefs; // Not used here but a pointer is still required

   Fint unit = l_list_units[iroot];

/*
  cout << "WRITING NO file header for root " << iroot<<" " << ntitle<<endl;
            mowrit_( unit, FLCOD, fierr, syserr, ntitle,
            &title[0], &afmt[0], nsym, &nbfpsy[0], &nmopsy[0], &labels[0],
            &l_coefs[0] );
*/
}
}


/* Dump the current set of coefficients and occupations to the irroot-th file.
   The must must have its header information prewritten and be open
*/
int PsociProcessMOcoefs::WriteNOcoefs( string s_title, int iroot, vector<Fdouble> & evals, vector<Fdouble> & coefs )
{
  if ( GA::nodeid() != READER ) return(0);

  Fint fierr = 0, syserr = 0;
  Fint ntitle = l_ntitle;
  Fint nsym = 0;

  vector<char> title( MO_MAX_TITLES * MO_MAX_TITLE_LENGTH,' ');  // 10 max titles at max 80 chars each

   vector<string> temp_title = l_title;

   string tempString = s_title;
   string scount;
   stringstream sout;
   sout << iroot;
   scount = sout.str();
   scount.insert(0,"-"); // preceding dash
   tempString.append( scount );

   tempString.resize( MO_MAX_TITLE_LENGTH,' '); 

   if ( ntitle < MO_MAX_TITLES ) {
           temp_title.push_back( tempString );
           ++ntitle;
   } else { //overwrite last entry
           temp_title[ MO_MAX_TITLES -1 ] =  tempString;
   }

// Push into an array

  for(int i=0; i< ntitle; ++i) {
    string word = temp_title[i];
    for(int j=0; j< word.size(); ++j ) {
       title[ MO_MAX_TITLE_LENGTH * i +  j ] = word[ j ];
    }
  }

// Prune to perfect size

  title.resize( ntitle * MO_MAX_TITLE_LENGTH );

  nsym = l_nsym;

  vector<char> afmt( MO_MAX_FORMAT_LENGTH,' ' );                   // max 80 chars
  vector<char> labels( MO_MAX_LABEL_LENGTH * nsym ,' ' ); // 8 max labels at max 4 chars e
  vector<Fint> nbfpsy( nsym,0 ); // 8 maximum labels
  vector<Fint> nmopsy( nsym,0 );

  for(int i=0; i< nsym; ++i) {
    nbfpsy[i] = l_nbfpsy[i];
  }

  for(int i=0; i< nsym; ++i) {
    nmopsy[i] = l_nmopsy[i];
  }

  for(int j=0; j< MO_MAX_FORMAT_LENGTH; ++j ) {
     afmt[ j ] = l_afmt[j];
  //   cout << "AFMT " << afmt[ j ] << endl;
  }

  for(int i=0; i< nsym; ++i) {
    string word = l_labels[i];
    for(int j=0; j< word.size(); ++j ) {
      labels[ MO_MAX_LABEL_LENGTH * i + j ] = word[ j ]; 
    }
  }

// Dump to nocoef-<n>

   vector<Fdouble> l_coefs; // Not used here but a pointer is still required

   Fint unit = l_base_unit;

// Ugh, push format data into a new string.

  string str_afmt;
  for(int i=0; i< MO_MAX_FORMAT_LENGTH; ++i ) {
     str_afmt.append(1,afmt[i] );
     if ( afmt[i] == ' ' ) // stop at the first space
        break;
     }
  int afmt_size = str_afmt.size();

  Fint FLCOD=10;
            mowrit_( unit, FLCOD, fierr, syserr, ntitle,
            &title[0], &str_afmt[0], nsym, &nbfpsy[0], &nmopsy[0], &labels[0],
            &coefs[0],  ntitle*MO_MAX_TITLE_LENGTH, afmt_size );


  if ( fierr != 0 || syserr != 0 ) {
     cout << "fierr and syserr should be zero. One or both are not " << endl;
     cout << "fierr is " << fierr<<" syserr is " << syserr << endl;
     return(-1);
  }

  
  FLCOD=20;
            mowrit_( unit, FLCOD, fierr, syserr, ntitle,
            &title[0],  &afmt[0], nsym, &nbfpsy[0], &nmopsy[0], &labels[0],
            &coefs[0], ntitle*MO_MAX_TITLE_LENGTH, afmt_size  );


  if ( fierr != 0 || syserr != 0 ) {
     cout << "fierr and syserr should be zero. One or both are not " << endl;
     cout << "fierr is " << fierr<<" syserr is " << syserr << endl;
     return(-1);
  }

  FLCOD=40;
            mowrit_( unit, FLCOD, fierr, syserr, ntitle,
            &title[0], &afmt[0], nsym, &nbfpsy[0], &nmopsy[0], &labels[0],
            &evals[0], ntitle*MO_MAX_TITLE_LENGTH, afmt_size  );


  if ( fierr != 0 || syserr != 0 ) {
     cout << "fierr and syserr should be zero. One or both are not " << endl;
     cout << "fierr is " << fierr<<" syserr is " << syserr << endl;
     return(-1);
  }

   return(0);
}



vector<Fdouble> * PsociProcessMOcoefs::fetchCoefsHandle() 
{
    return ( &l_coefs );
}

double PsociProcessMOcoefs::fetchMocoefs( int index )
{
    return ( l_coefs[index] );
}

int PsociProcessMOcoefs::fetchNsym()
{
     return( l_nsym );
}

int PsociProcessMOcoefs::fetch_nmopsy(int isym) 
{
     return( l_nmopsy[isym] );
}

int PsociProcessMOcoefs::fetch_nbfpsy(int isym)
{
     return( l_nbfpsy[isym] );
}


void PsociProcessMOcoefs::CloseMOfile()
{
  if ( GA::nodeid() == READER ) {
    Fint unit = l_unit;
    ftocxxclose_( unit );
    cout << "Closed MO coefficient file " << endl;
  }
}

void PsociProcessMOcoefs::CloseNOfiles()
{
    if ( GA::nodeid() == READER ) {
    vector<int>::iterator it;
    for(it=l_list_units.begin(); it !=l_list_units.end(); ++it ) {
       Fint unit = (*it);
       ftocxxclose_( unit );
     } 
    }
}

void PsociProcessMOcoefs::CloseNOfiles(int iroot )
{
    if ( GA::nodeid() == READER ) {
       Fint unit = l_base_unit;
       ftocxxclose_( unit );
    }
}



int PsociProcessMOcoefs::MOcoefficientRead()
{
  Fint FLCOD=20;
  Fint MAXTITLE = MO_MAX_TITLES; // cast change

  Fint fierr = 0, syserr = 0;
  Fint ntitle = 0, nsym = 0;

// These never get too big so simply allocate based on the maximum possible for Columbus

// Scratch
  vector<char> title( MO_MAX_TITLES * MO_MAX_TITLE_LENGTH ); // 10 max titles at max 80 chars each
  vector<char> afmt( MO_MAX_FORMAT_LENGTH ); // max 80 chars
  vector<char> labels( MO_MAX_LABEL_LENGTH * MO_MAX_LABELS ); // 8 max labels at max 4 chars e

  vector<Fint> nbfpsy(MO_MAX_LABELS,0 ); // 8 maximum labels
  vector<Fint> nmopsy(MO_MAX_LABELS,0 );

//TODO perhaps we can resize the coes to clen at the end ?

// Begin work
  l_coefs.resize( l_nbfsq ); // DUMMY for now
  Fint clen = l_coefs.size(); 
  
  Fint unit = l_unit;

  if ( GA::nodeid() == READER ) {  
            moread_( unit, FLCOD, fierr, syserr, MAXTITLE, ntitle,  
            &title[0], &afmt[0], nsym, &nbfpsy[0], &nmopsy[0], &labels[0], 
            clen, &l_coefs[0] );
  }

  l_nsym = nsym;
  l_ntitle = ntitle;

  if ( fierr != 0 || syserr != 0 ) {
     cout << "fierr and syserr should be zero. One or both are not " << endl;
     return(-1);
  }

// Process fortran created objects into what we want
  
  for(int i=0; i< ntitle; ++i) {
    string word;
    for(int j=0; j< MO_MAX_TITLE_LENGTH; ++j ) {
      word.push_back( title[ MO_MAX_TITLE_LENGTH * i +  j ] );
    }
    l_title.push_back( word );
  }

  //Symmetry distribution
  for(int i=0; i< nsym; ++i) {
    l_nbfpsy.push_back( nbfpsy[i] );
  }

  for(int i=0; i< nsym; ++i) {
    l_nmopsy.push_back( nmopsy[i] );
  }

  for(int i=0; i< nsym; ++i) {
    string word;
    for(int j=0; j< MO_MAX_LABEL_LENGTH; ++j ) {
      word.push_back( labels[ MO_MAX_LABEL_LENGTH * i +  j ] );
    }
    l_labels.push_back( word );
  }

  for(int j=0; j< MO_MAX_FORMAT_LENGTH; ++j ) {
    l_afmt.push_back( afmt[ j ] );
  }

  return( fierr );
}

void PsociProcessMOcoefs::printMOparameters() 
{
  if (GA::nodeid() != READER ) return;

  cout << "------------------------------------------------------" << endl << endl;
  cout << "MOCOEF parameters read from filename = " << l_filename << endl;
  cout << "  unit: " << l_unit<<endl;
  cout << "ntitle: " << l_ntitle<<endl;
  cout << "  nsym: " << l_nsym<<endl;
  cout << "mx cln: " << l_nbfsq << endl;

  vector<string>::iterator sit; 
  cout << endl << "Set of read titles " << endl;
  for(sit=l_title.begin(); sit!=l_title.end(); ++sit) {
     cout << (*sit) << endl;
  }

  cout << endl << "Set of symmetry labels " << endl;
  for(sit=l_labels.begin(); sit!=l_labels.end(); ++sit) {
     cout << (*sit) << " ";
  }
  cout << endl;

  vector<Fint>::iterator it;
  cout << endl << "Set nbfpsy " << endl;
  for(it=l_nbfpsy.begin(); it!=l_nbfpsy.end(); ++it) {
     cout << (*it) << " ";
  }
  cout << endl;

  cout << endl << "Set nmopsy " << endl;
  for(it=l_nmopsy.begin(); it!=l_nmopsy.end(); ++it) {
     cout << (*it) << " ";
  }
  cout << endl;


 cout << endl << "coefficients were read with FORMAT of ";
 cout << l_afmt << endl;
 

 cout << endl << "Current size of coefs array ";
 cout << l_coefs.size() << endl;

   
 return;
}

/* Coeffs are packed into a linear array 
*/

// COEFS[basis][orbs]

void PsociProcessMOcoefs::printMOcoefficients()
{
     if (GA::nodeid() != READER ) return;

     int index = 0;

     for(int isym=0; isym < l_nsym ; ++isym) {
     cout << endl<< " SYMMETRY BLOCK " << isym+1 << " " << l_labels[isym] << endl << endl;
     int numorbs = l_nmopsy[isym];
     int numbasis = l_nbfpsy[isym];
        for(int iorb = 0; iorb < numorbs; ++iorb ) {
           for(int i=0; i< numbasis; ++i ) {
           cout << l_coefs[ index ]<< " ";
           ++index; // they are arranged in canonical order
           }
           cout << endl;
         }
      }
}

/* Copy out HEADER information. This way we can grab the data from INPUT orbitals
   and write out to compute (natural) orbitals
*/

int PsociProcessMOcoefs::copyOutHeaderInfo( int & nsym, vector<string> & title, vector<string> & labels, string & afmt, vector<Fint> & nbfpsy, vector<Fint> & nmopsy )
{
  if (GA::nodeid() != READER ) return(nsym);

  nsym = l_nsym;

  title.clear();
  title = l_title;

  labels.clear();
  labels = l_labels;

  afmt = l_afmt;

  nbfpsy.clear();
  nbfpsy = l_nbfpsy;

  nmopsy.clear();
  nmopsy = l_nmopsy;

  return( l_ntitle );
}

int PsociProcessMOcoefs::copyInHeaderInfo( int & nsym, int & ntitle, vector<string> & title, vector<string> & labels, string & afmt, vector<Fint> & nbfpsy, vector<Fint> & nmopsy )
{
  if (GA::nodeid() != READER ) return(nsym);

  l_nsym = nsym;
  l_ntitle = ntitle;

  l_title.clear();
  l_title = title;

  l_labels.clear();
  l_labels = labels;

  l_afmt = afmt;

  l_nbfpsy.clear();
  l_nbfpsy = nbfpsy;

  l_nmopsy.clear();
  l_nmopsy = nmopsy;

  return( nsym );
}

