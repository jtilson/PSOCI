//TODO remove the superfluous GA::sync() calls used for debugging
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
 *   @file PsociIntegrals.C
 *
 */

#include "PsociTimer.hpp"
#include "PsociIntegrals.hpp"

#define NNDXF(i) (i * (i - 1)) / 2
#define NNDAF(i) (i * (i + 1)) / 2 


// Class definitions - nearly all are collectives

using namespace std;

//Not everyone needs to call but it is advisable.
PsociIntegrals::PsociIntegrals( int reader, Fint unit, string filename )
{
  l_moints = filename;
  l_unit = unit;
  perform_1e_print=false;
  perform_2e_print=false;
  
  ROOT_READER = reader;
}

//Local

PsociIntegrals::PsociIntegrals()
{
  l_unit = DEFAULT_MOINTS_UNIT;
  l_moints = DEFAULT_MOINTS_FILENAME;
  perform_1e_print=false;
  perform_2e_print=false;
  ROOT_READER = 0;
}

int PsociIntegrals::fetchUnit(){
  return( l_unit );
}

void PsociIntegrals::set1ePrint()
{
  perform_1e_print = true;
}


void PsociIntegrals::set2ePrint()
{
  perform_2e_print = true;
}

//Local
void PsociIntegrals::printFilename()
{
  if (GA::nodeid() == ROOT_READER ) cout << "Current Integral file is named " << l_moints << " with the unit number = " <<  l_unit << endl;
}


// Only get the first core energy
double PsociIntegrals::fetchCoreEnergy()
{
    return( (double) l_energy[0] );
} 

int PsociIntegrals::fetchLabels( vector<string> & labels ) 
{
     labels = l_slabel;
     return ( labels.size() );
}

int PsociIntegrals::fetchBasisLabels( vector<string> & labels )
{
     labels = l_bfnlab;
     return ( labels.size() );
}

void PsociIntegrals::OpenFile()
{
  if ( GA::nodeid() == ROOT_READER ) {
    Fint size = l_moints.size();
    Fint istatus;
    ftocxxopen_( l_unit, l_moints.c_str(), size, istatus);
    if ( istatus != 0 ) GA::error(" Failed to open unformatted file ",l_unit );
  }
}

//Local
void PsociIntegrals::CloseFile()
{
  if ( GA::nodeid() == ROOT_READER ) {
    ftocxxclose_( l_unit );
  }
}

/*
  sifrh1_(l_unit, ntitle, nsym, nbft,
  ninfo, nenrgy, nmap, ierr);
*/

/* relevant text snipped from colib8.f
   c  input:
   c  aoints = input file unit number.
   c
   c  output:
   c  ntitle = number of titles.
   c  nsym = number of symmetry blocks.
   c  nbft = total number of basis functions.
   c  ninfo = number of record-definition parameters.
   c  nenrgy = number of core energies.  the first element must be the
   c           nuclear repulsion energy.
   c  nmap = number of optional map vectors.
   c  ierr = error return code.  0 for normal return.
   c
   c  26-jun-89 written by ron shepard.
*/

int PsociIntegrals::sifrh1()
{
  Fint unit = l_unit;
  Fint fierr = 0; 
  
  if (GA::nodeid() == ROOT_READER ) {
    sifrh1_( unit, l_ntitle, l_nsym, l_nbft,
	     l_ninfo, l_nenrgy, l_nmap, fierr); // SIfs Fortran77 call uses Integer*4
  }
  return( fierr ); 
}

void PsociIntegrals::brdcstParams() 
{
  Fint l_buffer[5];
  if ( GA::nodeid() == ROOT_READER ) {
    l_buffer[0]=l_nsym;
    l_buffer[1]=l_nbft;
    l_buffer[2]=l_ninfo;
    l_buffer[3]=l_nenrgy;
    l_buffer[4]=l_nmap;
  }
  GA::brdcst( l_buffer, sizeof(Fint)*5, ROOT_READER );
  if ( GA::nodeid() != ROOT_READER ) {
    l_nsym = l_buffer[0];
    l_nbft = l_buffer[1];
    l_ninfo = l_buffer[2];
    l_nenrgy = l_buffer[3];
    l_nmap = l_buffer[4];
  }
}

// Dump header 1 information
void PsociIntegrals::printSifrh1()
{
  if ( GA::nodeid() == ROOT_READER ) {
    cout << " Parameters from the SIFRH1 method reading filename = " << l_moints << endl;
    cout << " Number of titles: ntitle is " << l_ntitle << endl;
    cout << " Number of symmetry blocks: nsym is " << l_nsym << endl;
    cout << " Number of basis functions: nbft is " << l_nbft << endl;
    cout << " Number of record-definition blocks: ninfo is " << l_ninfo << endl;
    cout << " Number of core energies (nuc repulsion is first): nenrgy is " << l_nenrgy << endl;
    cout << " Number of (optional) map vectors: nmap is " << l_nmap << endl;
  }
}

// Wrapper around the Fortran sifrh2 method

int PsociIntegrals::sifrh2()
{
  
  /* More relevent text snipped from Columbus@colib8.f
     c  read header_2 from the standard integral file.
     c
     c  input:
     c  aoints = input file unit number.
     c  ntitle = number of titles.
     c  nsym = number of symmetry blocks.
     c  nbft = total number of basis functions.
     c  ninfo = number of record-definition parameters.
     c  nenrgy = number of core energies.  the first element must be the
     c           nuclear repulsion energy.
     c  nmap = number of optional map vectors.
     c
     c  output:
     c  title*80(1:ntitle) = identifying titles.
     c  nbpsy(1:nsym) = number of basis functions per symmetry block.
     c  slabel*4(1:nsym) = symmetry labels.
     c  info(1:ninfo) = record-definition parameters.
     c  bfnlab*8(1:nbft) = basis function labels.
     c  ietype(1:nenrgy) = core energy types.
     c  energy(1:nenrgy) = core energies.
     c  imtype(1:nmap) = map vector types.
     c  map(1:nbft,1:nmap) = basis function map vectors.
     c  ierr = error return code.  0 for normal return.
     c
     c  26-jun-89 written by ron shepard.
     c
  */

  if ( GA::nodeid() != ROOT_READER ) { // all others simply return 
    return( 0 ); //Not an error to be the wrong guy.
  }
  
  Fint unit = l_unit;
  Fint fierr;
  
  vector<char> title( l_ntitle * MAX_TITLE_LENGTH );
  vector<char> bfnlab( l_nbft * MAX_BFNLAB_LENGTH );
  vector<char> slabel( l_nsym * MAX_LABEL_LENGTH );
  vector<Fint> map( l_nbft * l_nmap );
  
  vector<Fint> nbpsy(l_nsym,0 );
  vector<Fint> info( l_ninfo,0 );
  vector<Fint> ietype(l_nenrgy,0 );
  vector<Fdouble> energy( l_nenrgy ,0 );
  vector<Fint> imtype( l_nmap,0 );
  
#if 0
  cout << "Unit is " << unit << " onto the sifrh2: l_ntitle, l_nsym, l_nbft, l_ninfo, l_nenrgy, l_nmap are : " ;
  cout << " " << l_ntitle << " " <<  l_nsym  << " " <<  l_nbft  << " " <<  l_ninfo << " " <<  l_nenrgy << " " <<   l_nmap;
  cout << endl;
#endif
  
  sifrh2_( unit, l_ntitle, l_nsym, l_nbft,
	   l_ninfo, l_nenrgy, l_nmap, &title[0], 
	   &nbpsy[0], &slabel[0], &info[0], &bfnlab[0], 
	   &ietype[0], &energy[0], &imtype[0], &map[0], 
	   fierr ); 
  
  //Current set of titles
  for(int i=0; i< l_ntitle; ++i) {
    string word;
    for(int j=0; j< MAX_TITLE_LENGTH; ++j ) {
      word.push_back( title[ MAX_TITLE_LENGTH * i +  j ] ); 
    }
    l_title.push_back( word );
  }
  
  //Symmetry distribution
  for(int i=0; i< l_nsym; ++i) {
    l_nbpsy.push_back( nbpsy[i] );
  }
  
  // User specified sym labels
  for(int i=0; i< l_nsym; ++i) {
    string word;
    for(int j=0; j< MAX_LABEL_LENGTH; ++j ) {
      word.push_back( slabel[ MAX_LABEL_LENGTH * i +  j ] );
    }
    l_slabel.push_back( word );
  }
  // File record parameters required for fetching integrals
  for(int i=0; i< l_ninfo; ++i) {
    l_info.push_back( info[i] );
  }
  
  //Basis function labels 
  for(int i=0; i< l_nbft; ++i) {
    string word;
    for(int j=0; j< MAX_BFNLAB_LENGTH; ++j ) {
      word.push_back( bfnlab[ MAX_BFNLAB_LENGTH * i + j] );
    }
    l_bfnlab.push_back( word );
  }
  
  // Available core energies (-1 is the nuclear repulsion energy )
  for(int i=0; i< l_nenrgy; ++i) {
    l_ietype.push_back( ietype[i] );
  }
  
  // core energies ( atomic units )
  for(int i=0; i< l_nenrgy; ++i) {
    l_energy.push_back( energy[i] );
  }
  
  // Available map vector types
  for(int i=0; i< l_nmap; ++i) {
    l_imtype.push_back( imtype[i] );
    //cout << "IMTYPE is " <<  imtype[i] << endl;
  }
  
// Need to transpose data around for proper layout ( F77 expect this as (nbft)(nmap)
// Vector maps (nbft,nmap): We should never have a null set
// package up all basis ftns for a particular nsym


  int index=0;
  for(int i=0; i< l_nmap; ++i) {
//     vector<Fint> temp(1,0); //Ensures at least one zero is added
       vector<Fint> temp;
         for(int j=0; j< l_nbft; ++j) {
          index = i * l_nbft + j;
          temp.push_back( map[ index ] );
    }
    l_map.push_back( temp ); // one temp per nmap
  }


/*
  for(int i=0; i< l_nbft; ++i) {
    vector<Fint> temp(1,0); //Ensures at least one zero is added
    for(int j=0; j< l_nmap; ++j) {
      
      cout << "MAP i j v "<<i<<" "<<j<<" "<<  map[ l_nbft * i + j ] << endl;
      temp.push_back( map[ l_nbft * i + j ] );
    }
    l_map.push_back( temp );
  }
*/


  return( 0 );
}

/* BRDCST (some) sifrh2 information
   Namely the ONE core energy
*/
//Collective
void PsociIntegrals::brdcstSifrh2()
{
   int num = 1;
   if (GA::nodeid() != ROOT_READER ) l_energy.resize(1);
   int size = num * sizeof(Fdouble);
   GA::brdcst( &l_energy[0], size, ROOT_READER );
}


// Dump header 2 information
void PsociIntegrals::printSifrh2()
{
  if ( GA::nodeid() != ROOT_READER ) { // all others simply return 
    return;
  } 

  cout << " Parameters from the SIFRH2 method reading filename = " << l_moints << endl;
  vector<string>::iterator sit;
  vector<Fint>::iterator vit;
  vector<vector<Fint> >::iterator v2it;
  vector<Fdouble>::iterator dit;
  
  //
  cout << "Set of " << l_title.size() << " current TITLES. " << endl;
  for(sit=l_title.begin(); sit != l_title.end(); ++sit) {
    cout << (*sit) << endl; 
  }
  cout << endl;
  //
  cout << "Set of " << l_nbpsy.size() << " current SYMMETRIES. " << endl;
  for(vit=l_nbpsy.begin(); vit != l_nbpsy.end(); ++vit) {
    cout << (*vit) << endl; 
  }
  cout << endl;
  //
  cout << "Set of " << l_slabel.size() << " current symmetry LABELS  " << endl;
  for(sit=l_slabel.begin(); sit != l_slabel.end(); ++sit) {
    cout << (*sit) << endl; 
  }
  cout << endl;
  //
  cout << "Set of " << l_info.size() << " current FILE parameters. " << endl;
  for(vit=l_info.begin(); vit != l_info.end(); ++vit) {
    cout << (*vit) << endl; 
  }
  cout << endl;
  //Basis function labels 
  
  cout << "Set of " << l_bfnlab.size() << " basis ftns LABELS  " << endl;
  for(sit=l_bfnlab.begin(); sit != l_bfnlab.end(); ++sit) {
    cout << (*sit) << " ";
  }
  cout << endl;
  
  // Available core energies (-1 is the nuclear repulsion energy )
  cout << "Set of " << l_ietype.size() << " available core energies. " << endl;
  for(vit=l_ietype.begin(); vit != l_ietype.end(); ++vit) {
    cout << (*vit) << endl; 
  }
  cout << endl;
  // core energies ( atomic units )
  cout << "Set of " << l_energy.size() << " core energies (au) " << endl;
  for(dit=l_energy.begin(); dit != l_energy.end(); ++dit) {
    cout << (*dit) << endl;
  }
  cout << endl;
  // Available map vector types
  cout << "Set of " << l_imtype.size() << " available imtype vectors. " << endl;
  for(vit=l_imtype.begin(); vit != l_imtype.end(); ++vit) {
    cout << (*vit) << endl;
  }
  // Vector maps (nbft,nmap)
  cout << "Set of " << l_imtype.size() << " available mappints per bftn. " << endl;
  if (l_imtype.size() > 0 ) {
    for(v2it=l_map.begin(); v2it != l_map.end(); ++v2it) {
      for( vit = (*v2it).begin(); vit != (*v2it).end(); ++vit ) {
	cout << (*vit) << " ";
      }
      cout << endl;
    }
  }
}  

// Fetch the one-electron integrals ( nuclear replsion and Im SO(x,y,z). Skip the others
// return number of read integrals ( or <= 0 on error )


/* ROOT_READER reads and process the integrals from moints.  The raw data are
   unpacked and processed and summed into the object l_valone. Data are NOT brdcst to the 
   other cores
*/
int PsociIntegrals::fetchOneElectronInts()
{
  Fint ierr;
  const int FETCH_FAILED=-1;
  
  if ( GA::nodeid() != ROOT_READER ) {
    return( GA::nodeid() );
  }
  
  // This is the only scenario needed by the SOCI code
  
  Fint nipv = 2; //return one orbital label in each label entry
  Fint iretbv = 0; //no bitvectors are returned
  
  if ( iretbv != 0) { 
    cerr << "fetchOneElectronInts only support an iretbv type of 0: specified was " << iretbv << endl;
    return( FETCH_FAILED );
  }
  
  //input
  Fint l1rec = l_info[1]; 
  Fint n1max = l_info[2];
  
  Fint num ;
  Fint last ;
  Fint itypea, itypeb;
  Fint ifmt ;
  Fint ibvtyp;
  Fint ibitv[2]; //ditto
  
  Fdouble fcore = 0.0;
  
  const int nomore = 2;
  
  l_valone.resize( NNDAF( l_nbft ), 0 ); // maximum size 1e of integral array 
  
  //Local scratch used by the Fortran routines
  
  vector<Fdouble> buffer(l1rec,0);
  vector<Fdouble> values(n1max,0);
  vector<Fint> labels( n1max*nipv, 0 ); //Treat as pooled storage for use by F77
  
  //Loop over records
  
  last = 0;
  while( last != nomore  ) {
    
    sifrd1_( l_unit, &l_info[0], nipv, iretbv, &buffer[0],
	     num, last, itypea, itypeb, ifmt, ibvtyp,
	     &values[0], &labels[0], fcore, ibitv, ierr );
    
    //This does the right thing easy conversion from Fortran to C++
    
    int index=0;
    
    // Transpose data back to C++ orientation
    vector<vector<Fint> >relabels(nipv, vector<Fint>( num, 0 ));
    for(int i=0; i < 2*num; i += 2) {
      relabels[0][index]=labels[i];
      relabels[1][index]=labels[i+1];
      //cout << i << " INT " << values[i] << " " << labels[i] << " " << labels[1+i] << " " << relabels[0][index] << " " << relabels[1][index] << endl;
      ++index;
    }
    
    //Process this batch of integrals
    
    int ijndx;
    int ilab, jlab;
    double ecore = 0.0;
    
    //
    //   char chrtyp[MAX_TYPE_LENGTH];
    //   siftyp_( itypea, itypeb, chrtyp );
    //   cout << chrtyp[0] << " " << chrtyp[1] << " " << chrtyp[2] << endl;
    //   string wchrtyp = chrtyp;
    //   cout << "chrtyp value is " << wchrtyp << endl; // doesn;t work correctly
    //   cout << "Start the loops " << endl;
    //   cout << "itypea ityperb " << itypea << " " << itypeb << endl;
    
    if ( itypea == 0 ) {
      if ( itypeb != 0 ) {
        for(int i = 0; i < num; ++i ) {
          //cout << "num is " << i << " Labels are " << relabels[0][i] << " " << relabels[1][i] << endl;
          ilab = max( relabels[0][i], relabels[1][i] );
          jlab = min( relabels[0][i], relabels[1][i] );
          ijndx = NNDXF(ilab) + jlab;
          // cout << "ijndx-1 " << values[i]<<" " <<l_valone[ ijndx-1 ]<< endl;;
          l_valone[ ijndx-1 ] += values[i];
        }
	
	if( last != 0 ) {
	  ecore += fcore;
	}
      }
    } else if (itypea == 2 ) {
      if ( itypeb <= 2 ) {
        for(int i = 0; i < num; ++i ) {
          ilab = relabels[0][i];
          jlab = relabels[1][i];
          // cout << "itypea==2: Labels are " << relabels[0][i] << " " << relabels[1][i] << endl;
          ijndx = NNDXF(ilab) + jlab;
          if ( ilab < jlab ) {
	    l_valone[ NNDXF(jlab) + ilab-1 ] = -values[i] ;
          } else {
	    l_valone[ NNDXF(ilab) + jlab-1 ] = values[i] ;
          }
	}
      } else {
	cout << "Skipping unknown itypea condition " << itypea << endl;
      }
    }
  } //last flag
  got_1e = true; //Need this to inform fetchTwoE that the file position is set
  return (ierr );
}

/* BRDCST all 1e integrals from ROOT_READER to all nodes
 */
int PsociIntegrals::brdcstOneElectronInts()
{
  if ( GA::nodeid() != ROOT_READER ) {
     l_valone.resize( NNDAF( l_nbft ) );
  }

  int estimate = NNDAF( l_nbft ) * sizeof(Fdouble);
  int size = l_valone.size() * sizeof(Fdouble);

//Test ensure ROOTS data are the correct size
  if ( size != estimate ) {
    cerr << "Error: estimate and size are not equal estimate = " << estimate << " size = " << size << endl;
    GA::Terminate();
  }
  
  GA::brdcst( &l_valone[0],  size,  ROOT_READER );
  got_1e = true; // Now everyone gets this set 
  return( l_valone.size() );
}

/* Local call: print out values of the 1e integrals
   Generally this is not too useful. Ints MUST exist. So if node<>RootREADER attempts to print
   You must have brdcst already 
 */
void PsociIntegrals::printOneElectronInts()
{
  if ( GA::nodeid() != ROOT_READER ) return;

  if ( got_1e != true ) {
    cerr << "printOneElectronInts(): " << GA::nodeid() << " was requested to print a 1e array but got_1e != true: Skipping " << endl;
    return;
  }
  
  //char heading[MAX_ONEINT_HEADING] = "h(core)+ hso(x)+ hso(y)+ hso(z)";
  
  int size = l_valone.size();
  if ( size != NNDAF( l_nbft ) ) {
    cerr << " Wrong number of 1e integrals: was " << size << " should be: " << NNDAF(l_nbft) << endl;
    exit(1);
  }
  cout << "Current number of 1e integrals is " << size << endl;;
  cout << "Integrals below the value " << MIN_DOUBLE_PRINTED << " are not displayed " << endl;
  
  cout << "Current headers for the printed integral files " << endl;
  cout << "---------------------------------------------- " << endl;
  vector<string>::iterator sit;
  for(sit=l_title.begin(); sit != l_title.end(); ++sit) {
    cout << (*sit) << endl;
  }   
  cout << endl;
  if ( perform_1e_print ) 
    {
      cout << " i j integrals " << endl;
      cout << " -------------" << endl;
      
      //Build i,j,k,el indexing based on canonical order.
      
      int ij = 0;

//INDEXING is still not perfect

      Fdouble vnt; 
      for(int i=0; i< l_nbft; ++i) {
	for(int j=0; j<= i; ++j) {
          vnt = fetch1eIntegralValue(i,j);
	  if ( abs(vnt) >= MIN_DOUBLE_PRINTED ) { 
	    cout << setw(4) << i << setw(4) << j << setw(15) << setiosflags(ios::fixed) << setprecision(9) << vnt << endl;
	  }
	  ++ij;
	}
      }
    } else {
    cout << "PrintOneElectronInts was requested but perform_1e_print not set to true: skipping " << endl;
  }
}

/* Fetch and construct 2e integrals from MOINTS file
 */

int PsociIntegrals::fetchTwoElectronInts()
{
  /*
    Reproduced from columbus:colib8.f
    c
    c  read and decode a 2-e integral record.
    c
    c  input:
    c  aoint2 = input file unit number.
    c  info(*) = info array for this file.
    c  nipv = number of integers per value to be returned.
    c       = 0 only unpack dword.  values(*), labels(*), and ibitv(*)
    c           are not referenced.
    c       = 1 return four orbital labels packed in each labels(*) entry.
    c       = 2 return two orbital labels packed in each labels(*) entry.
    c       = 4 return one orbital label in each labels(*) entry.
    c  iretbv = bit vector request type.
    c     if ( iretbv=0 ) then
    c         null request, don't return ibitv(*).
    c     elseif ( iretbv=ibvtyp ) then
    c         request return of the bit-vector of type iretbv.
    c     elseif ( iretbv=-1 .and. ibvtyp<>0 ) then
    c         return any type of bit-vector that is on the record.
    c     else
    c        error. requested bit-vector is not available in buffer(*).
    c     endif
    c  buffer(1:l2rec) = packed  buffer with l2rec=info(4).
    c
    c  output:
    c  num = actual number of values in the packed buffer.
    c  last = integral continuation parameter.
    c  itypea,itypeb = generic and specific integral types.
    c  ifmt = format of the packed buffer.
    c  ibvtyp = type of packed bit-vector.
    c  values(1:num) = values (referenced only if nipv.ne.0).
    c  labels(1:nipv,1:num) = integral labels
    c           (referenced only if nipv.ne.0).
    c           note: if ifmt=0, then as many as ((nipv*n2max+7)/8)*8
    c                 elements of labels(*) are referenced.
    c  ibitv(*) = unpacked bit vector (referenced only if iretbv.ne.0).
    c             note: as many as ((n2max+63)/64)*64 elements of this
    c                   array are referenced.
    c  ierr = error return code.  0 for normal return.
    c
    c  26-jun-89 written by ron shepard.
    c
  */
  Fint ierr;
  const int FETCH_FAILED=-1;
  
  if ( GA::nodeid() != ROOT_READER ) {
    return( 0 );
  }
  
  // This are the only scenario currently needed by the SOCI code
  
  Fint nipv = 4; //return one orbital label in each label entry
  Fint iretbv = 0; //no bitvectors are returned
  
  if ( iretbv != 0) { 
    cerr << "fetchTwoElectronInts only support an iretbv type of 0: specified was " << iretbv << endl;
    //Else you will need to open up memory for ibitv = ((n2max+63)/64)*64
    return( FETCH_FAILED );
  }
  
  //input
  Fint l2rec = l_info[3]; 
  Fint n2max = l_info[4];
  
  Fint num ;
  Fint last ;
  Fint itypea, itypeb;
  Fint ifmt ;
  Fint ibvtyp;
  Fint ibitv[1]; //ditto
  
  const int nomore = 2;
  
  l_valtwo.resize( NNDAF( NNDAF( l_nbft ) ) ); //Canonical array of 2E ints 
  
  cout << "Maximum size of valtwo is set to be " << NNDAF( NNDAF( l_nbft ) ) << endl;
  
  //Local scratch used by the Fortran routines
  
  vector<Fdouble> buffer( l2rec,0 );
  vector<Fdouble> values( n2max,0 );
  vector<Fint> labels( n2max*nipv, 0 ); //Treat as pooled storage for use by F77
  
  // Ron Treats 2e integrals potentially differently..... Prefers to use SIFS to open close files
  // reproduced from the colib8:sif2cf method:
  /*
    c  the correct operation order in the calling program is:
    c     open(unit=aoints,...)        # standard open for the 1-e file.
    c     call sifo2f(..aoint2.)       # open the 2-e file.
    c     call sifc2f(aoint2...)       # close the 2-e file.
    c     close(unit=aoints...)        # close the 1-e file.
  */
  // We can skip this if we ensure got_1e == true and FSPLIT==1
  
  if ( got_1e != true ) {
    cerr << "fetchTwoElectronInts only uses the FSPLIT=1 (combined file) integral processing ";
    cerr << "This requires the file stream to be positioned just after the 1E read. It is currently not: Aborting " << endl;
    GA::Terminate();
  }
  
  last = 0;
  while( last != nomore  ) {
    
    sifrd2_( l_unit, &l_info[0], nipv, iretbv, &buffer[0],
	     num, last, itypea, itypeb, ifmt, ibvtyp,
	     &values[0], &labels[0], ibitv, ierr );
    
    //This does the right thing easy conversion from Fortran to C++
    int index=0;
    
    //Transpose data back to C++ order;
    //NOT FULLY TESTED - depends on being in canonical order
    vector<vector<Fint> >relabels(nipv, vector<Fint>( num, 0 ));
    for(int i=0; i < 4*num; i += 4) {
      relabels[0][index]=labels[i];
      relabels[1][index]=labels[i+1];
      relabels[2][index]=labels[i+2];
      relabels[3][index]=labels[i+3];
      ++index;
    }
    
// THis is probablky the ONLY time labels begining at one wil be used
    int ilab, jlab, klab, ellab;
    for(int i = 0; i < num; ++i ) { // These labels are fetched as Fortran style. Subtract 1 and never look back
      ilab = relabels[0][i];
      jlab = relabels[1][i]; 
      klab = relabels[2][i];
      ellab = relabels[3][i];
      //cout << "labels " << ilab << " " << jlab << endl;
      l_valtwo[ canonicalIndexF77(ilab,jlab,klab,ellab)-1 ] = values[i];
    }
  }
  got_2e = true;
  return( 0 );
}


/* BRDCST all 2e integrals from ROOT_READER to all nodes
 */
int PsociIntegrals::brdcstTwoElectronInts()
{
  if ( GA::nodeid() != ROOT_READER ) { 
     l_valtwo.resize( NNDAF( NNDAF( l_nbft ) ) );
  }
  int size = l_valtwo.size() * sizeof(Fdouble);
  int estimate = NNDAF( NNDAF( l_nbft ) ) * sizeof(Fdouble);
  
  if ( size != estimate ) {
    cerr << "Error: 2electron estimate and size are not equal estimate = " << estimate << " size = " << size << endl;
    GA::Terminate();
  }
  
  GA::brdcst( &l_valtwo[0],  size,  ROOT_READER );
  got_2e = true; // Now everyone gets this set 
  return( l_valtwo.size() );
}


/* Local call: Dump the current set of two-electron integrals to stdout. This can be a fairly sizable file
   Useful for ensuring (debugging) a parallel app but not intended for general usage
*/
void PsociIntegrals::printTwoElectronInts()
{
  if ( GA::nodeid() != ROOT_READER ) return;
  if ( got_2e != true ) {
    cerr << "printTwoElectronInts(): " << GA::nodeid() << " was requested ot print a 2e array but got_2e != true: Skipping " << endl;
    return;
  }
  
  if ( l_valtwo.size() > MAX_TWOE_PRINTSIZE ) {
    cout << "l_valtwo.size() > MAX_TWOE_PRINTSIZE and is rediculous to print: Cancelling. If you really need to print this ";
    cout << " Adjust the parameters MAX_TWOE_PRINTSIZE to be greater than " << l_valtwo.size() << " and recompile " << endl;
    return;
  }
  
  Fdouble vnt;

  if ( perform_2e_print ) {
    
    int size = l_valtwo.size();
    if ( size != NNDAF( NNDAF( l_nbft ) ) ) {
      cerr << " Wrong number of 2e integrals: was " << size << " should be: " << NNDAF( NNDAF(l_nbft) ) << endl;
      exit(1);
    }
    cout << endl << "Current total number of 2e integrals is " << size << endl;; 
    cout << "Integrals below the value " << MIN_DOUBLE_PRINTED << " are not displayed " << endl;
    
    cout << "Current headers for the printed integral files " << endl;
    cout << "---------------------------------------------- " << endl;
    vector<string>::iterator sit;
    for(sit=l_title.begin(); sit != l_title.end(); ++sit) {
      cout << (*sit) << endl;
    }
    cout << endl;
    cout << " i j k l integrals " << endl;
    cout << " -----------------" << endl;
    
    //Build i,j,k,el indexing based on canonical order.
    
    int elmax;
    
    for(int i=0; i< l_nbft; ++i) {
      for(int j=0; j<= i; ++j) {
	for(int k=0; k<= i; ++k) {
	  (k == i) ? elmax = j: elmax = k ;
	  for(int el=0; el <= elmax; ++el) {
            vnt = fetch2eIntegralValue(i, j, k, el);
	    if ( abs(vnt)  >= MIN_DOUBLE_PRINTED ) {
	      cout << setw(4) << i << setw(4) << j << setw(4) << k << setw(4) << el << setw(15) << setiosflags(ios::fixed) << setprecision(8) <<  vnt << endl;
	    }
	  }
	}
      }
    }
  } else {
    cout << "PrintTwoElectronInts was requested but perform_2e_print not set to true: skipping " << endl;
  }
}

/* TODO check for performance of accesing intergral in the calling program.
   I think the callig program amortizes this call over lots of flops.
   Nevertheless, performance should be measured.
*/

//Index using C-style indexing
double PsociIntegrals::fetch2eIntegralValue(int i, int j, int k, int el)
{
   int index = canonicalIndex(i, j, k, el ); // C-indexing style
   if (index > l_valtwo.size() ) {
      cerr << "Erroneous valtwo index: Aborting " << index<<" " <<l_valtwo.size() << " "<<i<<" " << j << " " << k << " " << el << endl;
      GA::Terminate();
   }
   double value = l_valtwo[ index - 1 ];
   return( value );
}

double PsociIntegrals::fetch1eIntegralValue(int i, int j)
{
   int index = TwoEIndex( i, j );
   if (index > l_valone.size() ) {
      cerr << "Erroneous valone index: Aborting " << i << " " << j << endl;
      GA::Terminate();
   }
   double value = l_valone[ index - 1 ];
   return( value );
}


/* Compute a linearized index from a 4-index canonical order
   intend to use C-style numbering Here we enter with F77 style indexing (starting at 1 )
*/
int PsociIntegrals::canonicalIndexF77(int i, int j, int k, int el)
{
  int ij = (((max(i,j)) - 1)*(max(i,j)))/2 + min(i,j);
  int kl = (((max(k,el)) - 1)*(max(k,el)))/2 + min(k,el);
  return( (((max(ij,kl)) - 1)*(max(ij,kl)))/2 + min(ij,kl) );
}

/* Compute a linearized index from a 4-index canonical order
   intend to use C-style numbering */
// Here we enter using C indexing style: added 
int PsociIntegrals::canonicalIndex(int iC, int jC, int kC, int elC)
{
  int i = iC + 1;
  int j = jC + 1;
  int k = kC + 1;
  int el = elC + 1;
  int ij = (((max(i,j)) - 1)*(max(i,j)))/2 + min(i,j);
  int kl = (((max(k,el)) - 1)*(max(k,el)))/2 + min(k,el);
  return( (((max(ij,kl)) - 1)*(max(ij,kl)))/2 + min(ij,kl) );
}


/* Compute a linearized index from a 2-index canonical order
   Requires using indexing that begins at ZERO */
int PsociIntegrals::TwoEIndex(int iC, int jC)
{
  int i = iC + 1;
  int j = jC + 1;
  int ij = ( (max(i,j) - 1)*max(i,j) )/2 + min(i,j);
  return( ij );
}

// Collective. THis method basically grabs and brdcst's all the integrals.
// TODO add some checks to ensure file was opened 

int  PsociIntegrals::fetchAndProcessIntegrals()
{
   int status = 1;
   int PRINT1eINTS=0;
   int PRINT2eINTS=0;

   printFilename();
   Fint unit = fetchUnit();
   OpenFile();

   if ( sifrh1() != 0 ) {
     cerr << "sifrh1 not zero value aborting: " << endl;
     GA::error("sifrh1:probably incompatible 32bit/64bit builds",-1);
   }
   if ( sifrh2() != 0 ) {
     cerr << "sifrh2 not zero value aborting: " << endl;
     cerr << "probably incompatible 32bit/64bit builds " << endl;
     GA::error("sifrh2: probably incompatible 32bit/64bit builds",-1); 
   }

// Ensure all nodes now have the parameters

   brdcstParams();
   brdcstSifrh2(); //For now only core energies

   if ( GA::nodeid() == ROOT_READER ) cout << endl << endl;

   printSifrh1();
   if ( GA::nodeid() == ROOT_READER ) cout << endl << endl;

   printSifrh2();

   status = fetchOneElectronInts();
   status = brdcstOneElectronInts();

   if( GA::nodeid() == ROOT_READER ) PRINT1eINTS = 1 ;
   if( GA::nodeid() == ROOT_READER ) PRINT2eINTS = 0 ;


   if( PRINT1eINTS == 1 ) set1ePrint();
   printOneElectronInts();

   fetchTwoElectronInts();
   int statb = brdcstTwoElectronInts();
   if( PRINT2eINTS == 1 ) set2ePrint();
   printTwoElectronInts();
   return( status );
}

/* NOTE: in constrast to the regular method for fetching 1e integrals
   this method 1) is NOT parallel and thus does not BRDCST out data and
   2) only fetches the overlap ints.
*/
int  PsociIntegrals::fetchAndProcessOverlapIntegrals()
{
   if ( GA::nodeid() != ROOT_READER ) return( 0 );
   int status = 1;
   int PRINT1eINTS=0;

   printFilename();
   Fint unit = fetchUnit();
   OpenFile();

   if ( sifrh1() != 0 ) {
     cerr << "sifrh1 not zero value aborting: " << endl;
     GA::error("sifrh1:probably incompatible 32bit/64bit builds",-1);
   }
   if ( sifrh2() != 0 ) {
     cerr << "sifrh2 not zero value aborting: " << endl;
     cerr << "perhaps incompatible 32bit/64bit builds " << endl;
     GA::error("sifrh2: possibly file was alreayd open or incompatible 32bit/64bit builds",-1); 
   }

#ifdef DETAILEDCHECK
   cout << endl << endl;
   printSifrh1();
   cout << endl << endl;
   printSifrh2();
#endif
   status = fetchOverlapInts();

   if( GA::nodeid() == ROOT_READER ) PRINT1eINTS = 0 ;

   if( PRINT1eINTS == 1 ) set1ePrint();
   
   if( PRINT1eINTS == 1 )  printOneElectronInts(); // prints out some info

   return( status );
}



// Fetch the overlap integrals for use by the Mullikan population analysis


/* ROOT_READER reads and process the integrals from aoints.  The raw data are
   unpacked and processed and summed into the object l_valone. Data are NOT brdcst to the 
   other cores

   For now the filename MUST be aoints. I don't think aoints2 formatting works
*/

int PsociIntegrals::fetchOverlapInts()
{
  Fint ierr;
  const int FETCH_FAILED=-1;
  
  if ( GA::nodeid() != ROOT_READER ) GA::nodeid();
  
  // This is the only scenario needed by the SOCI code
  

  Fint nipv = 2; //return one orbital label in each label entry
  Fint iretbv = 0; //no bitvectors are returned
  
  if ( iretbv != 0) { 
    cerr << "fetchOverkapInts only support an iretbv type of 0: specified was " << iretbv << endl;
    return( FETCH_FAILED );
  }
  
  
  // This is the only scenario needed by the SOCI code
  if ( iretbv != 0) { 
    cerr << "fetchOverlapInts only support an iretbv type of 0: specified was " << iretbv << endl;
    return( FETCH_FAILED );
  }
  
  //input
  Fint l1rec = l_info[1]; 
  Fint n1max = l_info[2];
  
  Fint num ;
  Fint last ;
  Fint itypea, itypeb;
  Fint ifmt ;
  Fint ibvtyp;
  Fint ibitv[2]; //ditto
  
  Fdouble fcore = 0.0;
  
  const int nomore = 2;
  
  l_valone.resize( NNDAF( l_nbft ), 0 ); // Reuse for overlap integrals maximum size 1e of integral array 
  
  //Local scratch used by the Fortran routines
  
  vector<Fdouble> buffer(l1rec,0);
  vector<Fdouble> values(n1max,0);
  vector<Fint> labels( n1max*nipv, 0 ); //Treat as pooled storage for use by F77
  
  //Loop over records
  
  last = 0;
  while( last != nomore  ) {
    
    sifrd1_( l_unit, &l_info[0], nipv, iretbv, &buffer[0],
	     num, last, itypea, itypeb, ifmt, ibvtyp,
	     &values[0], &labels[0], fcore, ibitv, ierr );
    
    //This does the right thing easy conversion from Fortran to C++
    
    int index=0;
    
    // Transpose data back to C++ orientation
    vector<vector<Fint> >relabels(nipv, vector<Fint>( num, 0 ));
    for(int i=0; i < 2*num; i += 2) {
      relabels[0][index]=labels[i];
      relabels[1][index]=labels[i+1];
      //cout << i << " INT " << values[i] << " " << labels[i] << " " << labels[1+i] << " " << relabels[0][index] << " " << relabels[1][index] << endl;
      ++index;
    }
    
    //Process the overlap integrals: Array is nbf*(nbf+1)/2: So it is packed but not symmnetry blocked
    
    int ijndx;
    int ilab, jlab;
    
    if ( itypea == 0 ) { // THis is the only scenrio of interest
      if ( itypeb == 0 ) {
        for(int i = 0; i < num; ++i ) {
          //cout << "num is " << i << " Labels are " << relabels[0][i] << " " << relabels[1][i] << endl;
          ilab = max( relabels[0][i], relabels[1][i] );
          jlab = min( relabels[0][i], relabels[1][i] );
          ijndx = NNDXF(ilab) + jlab;
          l_valone[ ijndx-1 ] += values[i];
        }
      }
    }

  } //nomore

  got_1e = true; //ALso need this to print out integrals if desired. Need this to inform fetchTwoE that the file position is set
  return (ierr );
}

