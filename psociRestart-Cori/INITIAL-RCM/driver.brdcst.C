/**************************************************************************************
 * Copyright (c) 2010 RENCI.
 * All rights reserved. This program and the accompanying materials
 * MAY BE available under the terms of the RENCI Open Source License
 * UNC at Chapel Hill which accompanies this distribution, and is available at
 * http://www.renci.org/resources/open-source-software-license

 * New implementation of PSOCI:

 Classes: 

 Description: 
 
 History:

**************************************************************************************/
/**
 *   @file driverTestVector.C
 *
 */
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>


#include <cmath>
// Now start testing for GA
#include <ga++.h>
//#include <GAServices.h>

#include <dra.h>
#define GA_DATA_TYPE C_DBL

using namespace std;


// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//
//  Start driver program 
//
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#define GCOUT if (g_rank==0) cout

int main(int argc, char **argv)
{
  // Set up global MPI fabric
  
  int g_size;
  int g_rank;
  
  // Set up Global Array fabric
  // heap and stack are per-core quantities - for collective operations
  // ALlocate in terms of DOUBLES for now

  const long wordSize = sizeof(double);
  const long OneGigaByte = 128*1024*1024 * wordSize; // A total of 1 Gig but Generally thinking interms of doubles (words)... 
  long maxMemPerCore = OneGigaByte;

  unsigned int heap=9000000, stack=9000000;
  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, 0);
  GA::setMemoryLimit( maxMemPerCore );
  g_size = GA::nodes();
  g_rank = GA::nodeid();
  if ( GA::usesMA() ) cout << "GA memory is coming from MA " << endl;
  if ( GA::usesFAPI() ) cout << "GA is using Fortran indexing " << endl;
  
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  GCOUT << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
  GCOUT << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl << endl;
  
  GCOUT << "Parallel GA Run with " << g_size << " Total processors " << endl;
  

//chars form forks fine

/*
  vector<char> listchars;
  listchars.reserve(3);

  if ( g_rank == 0 ) {
    listchars.push_back('1');
    listchars.push_back('3');
    listchars.push_back('2');
  }
  cout << g_rank << " BF brd name 0 " <<  listchars[0] << endl;
  cout << g_rank << " BF brd name 1 " <<  listchars[1] << endl;
  cout << "Size of arrays is " << listchars.size() << endl;

  GA::brdcst( &listchars[0], 3*sizeof(char), 0 );

  cout << g_rank << "after brd name 0 " <<  listchars[0] << " end" <<  endl;
  cout << g_rank << "after brd name 1 " <<  listchars[1] << " end" <<  endl;
*/

// Try a simple string send
/* THis doesn;t work as desired */

/*
  string test;
  if ( g_rank == 0 ) {
           test.append("johns test");
  } else {
           test.reserve(10);
  }
  cout << g_rank << "before john test = " << test << endl;
  cout << g_rank << "before john test size = " << test.length() << endl;
  
  GA::brdcst( &test, test.size(), 0 );

  cout << g_rank << "after john test =  " << test << endl;

//Another test 
  cout << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV " << endl;

  char ctest[20]="0000000000000000";

  if ( g_rank == 0 ) {
     ctest[0]='j';
     ctest[1]='o';
     ctest[2]='h';
     ctest[3]='n';
     ctest[4]='s';
     ctest[5]=' ';
     ctest[6]='t';
     ctest[7]='e';
     ctest[8]='s';
     ctest[9]='t';
   }

  cout << "Current size of the char array is " << strlen( ctest ) << endl;
  cout << g_rank <<  " contents are " << ctest << endl;
  
  GA::brdcst( &ctest[0], 20*sizeof(char), 0 );

  cout << "After brdcst: Current size of the char array is " << strlen( ctest ) << endl;
  cout << g_rank << " contents are " << ctest << endl;
*/

//Try with ints this fails as well

  cout << "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" << endl;

  vector<short> intchars;
  intchars.reserve(3);
    if ( g_rank == 0 ) {
    intchars.push_back(1);
    intchars.push_back(3);
    intchars.push_back(2);
  }


/*
  int intchars[3];
  if ( g_rank == 0 ) {
    intchars[0]=1;
    intchars[1]=3;
    intchars[2]=2;
  }
*/

  cout << g_rank << " INT brd name 0 " <<  intchars[0] << endl;
  cout << g_rank << " INT brd name 1 " <<  intchars[1] << endl;
//  cout << "Size of INT arrays is " << sizeof(intchars) << endl;
  cout << "Size of INT arrays is " << intchars.size() << endl;

  GA::brdcst( &intchars[0], 3*sizeof(short), 0 );

  cout << g_rank << "after INTbrd name 0 " <<  intchars[0] << " end" <<  endl;
  cout << g_rank << "after INTbrd name 1 " <<  intchars[1] << " end" <<  endl;



//Try using strings for now
/* These fail
  cout << " NEXT TEST " << endl;

  string name;
  name.resize(2);
  int nbf=2;
  int ncfgs=3;
  //vector<string> names(5); //subsequent push_backs ADD to this  
  vector<string> names; 
  names.reserve(3);

  if ( g_rank == 0 ) {
    names.push_back("13");
    ames.push_back("17");
    names.push_back("23");
  }

  cout << g_rank << "2nd BF brd name 1 " << names[0] << " end " << endl;
  cout << g_rank << "2nd BF brd name 2 " << names[1] << " end " << endl;

  vector<string>::iterator vit;
  for (vit=names.begin(); vit != names.end(); ++vit) {
      cout << g_rank << " value " << (*vit) << endl; 
  }
  cout << g_rank << "sizeof 2nd BF brd name 1 " << sizeof(names[0]) << endl;
  cout << g_rank << "size of 2nd BF brd name 2 " << sizeof(names[1]) << endl;

  cout << "Size of arrays is " << names.size() << endl;

  GA::brdcst( &names[0], 5*sizeof(char), 0 );


  cout << "AFTER BRDCST " << endl;
  for (vit=names.begin(); vit != names.end(); ++vit) {
      cout << g_rank << " value " << (*vit) << endl;
  }
  cout << g_rank << "2nd name 1 " << names[0] << endl;
  cout << g_rank << "2nd name 2 " << names[1] << endl;
*/


//New string attempts

/*
  vector<char> conf;
  cout << "Both getting the reserve " << endl;
  conf.reserve(20);
if ( g_rank == 0 ) {
         conf.push_back('2');
         conf.push_back('2');
         conf.push_back('2');
         conf.push_back('2');
         conf.push_back('0');
         conf.push_back('0');
         conf.push_back('1');
         conf.push_back('1');
   }

   cout << g_rank << " On to the BRDCST " << endl;
   GA::brdcst( &conf[0], 16*sizeof(char), 0 );
   cout << g_rank << " size of cit is " << conf.size() << endl;
   vector<char>::iterator cit;
   for (cit = conf.begin(); cit != conf.end(); ++cit) {
       cout << g_rank << " value " << (*cit) << endl;
   }
*/

/* also fails

  char conf[8];
  if ( g_rank == 0 ) {
    conf[0]='2';
    conf[1]='0';
    conf[2]='1';
    conf[3]='1';
    conf[4]='1';
    conf[5]='1';
    conf[6]='1';
    conf[7]='1';
   }

   for (int i=0; i<8; ++i) {
    cout << g_rank << " BEFORE the  BRD " << conf[i] << endl;
   }

   GA::brdcst( &conf[0], 8*sizeof(char), 0 );

   for (int i=0; i<8; ++i) {
    cout << g_rank << " after BRD " << conf[i] << endl;
   }


// more complicated exaple
   cout << " Go for a big example " << endl;

   vector<char> vconf;
   vconf.reserve(1);

   if ( g_rank == 0 ) {
     vconf.push_back(conf);
   }

   vector<char>::iterator vvit;
   for (vvit=vconf.begin(); vvit != vconf.end(); ++vvit ) {
         for (int i=0; i<8; ++i) {
             cout << g_rank << " before vconf " << (*vvit)[i] << endl;
         }
   }

   GA::brdcst( &vconf[0], 8*sizeof(char), 0 );

   for (vvit=vconf.begin(); vvit != vconf.end(); ++vvit ) {
         for (int i=0; i<8; ++i) {
             cout << g_rank << " after vconf " << (*vvit)[i] << endl;
         }
   }

*/

  cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX finished XXXXXXXXXXXXXXXXXXX " << endl;
  
  GA::Terminate();
  
}



