#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <iostream>
#include <cstdio>
#include <cmath>
using namespace std;
#include "ga++.h"


extern "C" {
#include "armci.h"
#include "message.h"
#include "mpi.h"
}

#define GCOUT if (me==0) cout
#define GA_DATA_TYPE MT_F_REAL


int
main(int argc, char *argv[]) {
  
  int ONE=1;   /* useful constants */
  int me, nproc;
  int i, row;
  
  int heap  = 200000;
  int stack = 200000;
  
  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, 0);
  me = GA_Nodeid();
  nproc = GA_Nnodes();

  GCOUT << "Parallel GA Run with " << nproc << " Total processors " << endl;


  int a = 0;
  for(row=me; row<me*10000; row+= 1){
      a=1;
  }


  cout << "NEW TESTING m" << me << endl;
  GA::SERVICES.sync();

  for( int i=0; i< nproc; ++i ) {
    if ( i == me ) {
        cout << "     my rank is " << me << endl;
    }
    GA::SERVICES.sync();
  }
  GA::SERVICES.sync();
   
  if(me==0) cout << "Terminating...\n";
  GA::Terminate();
  
}
