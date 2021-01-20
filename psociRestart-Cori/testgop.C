#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "mpi.h"


/* Sinple code to test the time for an basic MPI brdcst and gop.
   GA++ was about 0.8sec for a 10 MB object (768 cores) which seems 
   high 
*/

#define CHUNKSIZE 10 
#define NUMTEST 10
#define TESTSIZE 100*1024*1024
#define TESTSIZE2 10*1024*1024


#define GCOUT if (rank==0) cout 
using namespace std;

int main(int argc, char ** argv) {

  MPI_Init(&argc, &argv);
  int rank,size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &size);

    GCOUT << "================================================="<< endl;
    GCOUT << "=  Running with MPI                             ="<< endl;
    GCOUT << "================================================="<< endl;
    GCOUT << " rank is " << rank << " size is " << size << endl;

//Test one BRDCST out a vector of data 

    double tempTime = MPI_Wtime();
    vector<double> brdvector( TESTSIZE, 0.0 );
    vector<double> rcvvector( TESTSIZE, 0.0 );

    if (rank == 0 ) {
       for(int i=0; i< TESTSIZE; ++i ) {
       brdvector[i] = -2;
       }
    }
    MPI_Barrier( MPI_COMM_WORLD );

// DO the brdcst a number of times and average results
    
   double dcount =0;

   long runsize = TESTSIZE;
   int MAXTEST = 10;

   int testrun = 1;
   double tempLoop, tempLoop1;
   while ( runsize > 1 && testrun < MAXTEST ) {
   tempLoop =  MPI_Wtime();
   for (int itest=0; itest<NUMTEST; ++itest) {
       ++dcount;
       MPI_Bcast( &brdvector[0], runsize, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   }
   ++testrun;
   double loopTime = (  MPI_Wtime() - tempLoop ) / dcount;
   GCOUT << "LOOP is " << testrun << " Me is " << rank << " SIZE (W) " << runsize << " Bcast time (s)(ave) " << loopTime << endl;

//   runsize /= CHUNKSIZE;
   }

   runsize = TESTSIZE2;
   MPI_Barrier( MPI_COMM_WORLD );
   if ( rank == 0 ) cout << " XXXXX " << endl << endl;
   MPI_Barrier( MPI_COMM_WORLD );

/*
   while ( runsize > 1 ) {
   tempLoop =  MPI_Wtime();
   for (int itest=0; itest<NUMTEST; ++itest) {
       ++dcount;
       MPI_Bcast( &brdvector[0], runsize, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   }
   double loopTime = (  MPI_Wtime() - tempLoop ) / dcount;
   GCOUT << "TEST 2-SIZE (W) " << runsize << " Bcast time (s)(ave) " << loopTime << endl;
   runsize /= CHUNKSIZE;
   }

*/
   
// Do the GOP/Reduce test 

   double rcount =0;
   double tempReduce = MPI_Wtime();
/*
   for (int itest=0; itest<NUMTEST; ++itest) {
       ++rcount;
       MPI_Allreduce( &brdvector[0], &rcvvector[0], TESTSIZE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   }
*/
   double reduceTime = (  MPI_Wtime() - tempReduce ) / rcount;

    double wholeTime = MPI_Wtime() - tempTime;

    GCOUT << "Whole time " << wholeTime << " Reduce time " << reduceTime <<  endl;

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}


