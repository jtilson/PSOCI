/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: A rudementary class for managing servicess for DRA: Basically interfaces 
              to the DRA library 

 History:

 

**************************************************************************************/
/**
 *   @file PsociDRArestart.C
 *
 */

#include "PsociTimer.hpp"
#include "PsociDRArestart.hpp"

// Class definitions

using namespace std;

// Mostly these are collective calls.....

//One instance PER DRA array

PsociDRArestart::PsociDRArestart( int rank, string arrayname, string filename )
{
    //GA::sync();
    l_arrayname = arrayname;
    l_filename = filename;
    l_mode = DEF_DRA_MODE;
    l_type = DEF_DRA_TYPE;
}

//Local call
void PsociDRArestart::reportName()
{
    if (GA::nodeid() == 0 ) cout << "Current Array name is " << l_arrayname;
}

//Local call
void PsociDRArestart::reportFilename()
{
   if (GA::nodeid() == 0 ) cout << "Base filename to be used is " << l_filename;
}
//Collective call
void PsociDRArestart::setType( int type )
{
    //GA::sync();
    if ( (type != MT_INT) && (type != C_DBL) && (type != C_LONG)) {
       cerr << "DRA type must be either C_LONG, C_INT or C_DBL " << endl;
       return;
    }
    l_type = type;
}
//Local call
void PsociDRArestart::printType()
{
    if ( GA::nodeid() != 0 ) return;

    switch (l_type)
    {
    case C_DBL:
      cout << "Current array is of type C_DBL " << endl;
      break;
    case C_INT:
      cout << "Current array is of type C_INT " << endl;
      break;
    case C_LONG:
      cout << "Current array is of type C_LONG " << endl;
      break;
    default:
      cout << "Wrong type for DRAs " << l_type;
      break;
    }
}
//Collective call: must be DRA_R, DRA_W, or DRA_RW.
void PsociDRArestart::setMode( int mode )
{
   //GA::sync();
   l_mode = mode;
}

//Collective
void PsociDRArestart::setDimensions( int ndim, pair<dra_size_t,dra_size_t> dims, pair<dra_size_t,dra_size_t> rdims ) 
{
    //GA::sync();
    if ( ndim > MAX_DRA_DIMS ) {
       cerr << "Requested ndim > MAX_DRA_DIMS " << ndim << " MAX is " << MAX_DRA_DIMS << endl;
       return;
    } 
    l_ndim = ndim;
    l_dim0 = dims.first;
    l_dim1 = dims.second;
    l_rdim0 = rdims.first;
    l_rdim1 = rdims.second;
}

//Local call
void PsociDRArestart::getDimensions( int &ndim, pair<dra_size_t,dra_size_t> &dims, pair<dra_size_t,dra_size_t> &rdims )
{
    ndim = l_ndim;
    dims.first = l_dim0;
    dims.second = l_dim1;
    rdims.first = l_rdim0;
    rdims.second = l_rdim1;
}

//Collective 
int PsociDRArestart::createDRA()
{
    //GA::sync();

#if 0
    cout << "ndims is " << l_ndim << endl;
    cout << "Dim0 is " << l_dim0 << " dim1 is " << l_dim1 << endl;
    cout << "rim0 is " << l_rdim0 << " rdims1 is " << l_rdim1 << endl;
#endif

    dra_size_t ddims[MAX_DRA_DIMS];
    ddims[0] = l_dim0;
    ddims[1] = l_dim1;

    dra_size_t  rdims[MAX_DRA_DIMS];
    rdims[0] = l_rdim0;
    rdims[1] = l_rdim1;

   char * local_filename = new char[ l_filename.size()+1];
   strcpy (local_filename, l_filename.c_str() );

   char * local_name = new char[ l_arrayname.size()+1];
   strcpy (local_name, l_arrayname.c_str() );

#if 0
   cout << "l_type " << l_type << " " << ddims[0] << " " << ddims[1];
   cout << " " << local_name << " " << local_filename  << " " << l_mode;
   cout << " " << rdims[0] << " " << rdims[1];
   cout << "rank is " << GA::nodeid() << "DRA handle is " << local_dra << "END " << endl;
#endif

   NDRA_Create( l_type, l_ndim, ddims, local_name, local_filename, l_mode, rdims, &local_dra);
  // NDRA_Create( C_DBL, 2, ddims, local_name, local_filename, DRA_RW, rdims, &local_dra);
   return(0);
}

//Collective only node zero does any work
void PsociDRArestart::printInternals()
{
    //GA::sync();
    DRA_Print_internals( local_dra );
}

//Collective Dump a GA of size (n,m)  to a DRA of size (n,m);
int PsociDRArestart::dumpGAtoDRA( GA::GlobalArray * g_matrix )
{
    //GA::sync();

// Fetch and Check GA handle
   char handleMessage[] = "dumpGAtoDRA: ga handle invalid ";
   g_matrix->checkHandle( handleMessage );

//DRA patches  - should have been preset at DRA construction
/*
   dra_size_t dlo[l_ndim], dhi[l_ndim];
   dlo[0] = 0;
   dlo[1] = 0;
   dhi[0] = l_dim0-1;
   dhi[1] = l_dim1-1;
*/
#if 0
   cout << "dilo and dihi are " << dlo[0] << " " << dhi[0] << endl;
   cout << "djlo and djhi are " << dlo[1] << " " << dhi[1] << endl;
#endif

//GA patches more uglyness.. the DRA argument list expects INTS for glo and ghi  

   int ga_ndim;
   int dims[2];
   int type;
   g_matrix->inquire( &type, &ga_ndim, dims);
   if ( ga_ndim != 2 || ga_ndim != l_ndim ) {
      cerr << "dumpGAtoDRA: only a 2D GA can be dumped. Inquery reports this " << ga_ndim << endl;
      return(-1);
   }
   if ( type != C_DBL && type != C_INT && type != C_LONG ) {
      cerr << "Error: dumpGAtoDRA only handles types C_DBL, C_INT, or C_LONG " << endl;
      exit(1);
      return(-1);
   }
 
/*
   int glo[ga_ndim], ghi[ga_ndim];
   glo[0] = 0;
   ghi[0] = dims[0]-1;
   glo[1] = 0;
   ghi[1] = dims[1]-1;
*/

#if 0
   cout << "gilo and gihi are " << glo[0] << " " << ghi[0] << endl;
   cout << "gjlo and gjhi are " << glo[1] << " " << ghi[1] << endl;
#endif
      
   return( NDRA_Write(g_matrix->handle(),local_dra, &l_req) );
}
 
/*Collective its up to the appl to ensure the GA and DRA conform. I don't have time to deal with this just now.
  The vectors are of the form I(lo,hi), J(lo,hi);
  This is complicated by the fact that only part of the GA may be stored

*/
int PsociDRArestart::dumpGASectionDRA( vector<pair<int,int> > ga_patch, vector< pair<dra_size_t,dra_size_t> > dra_patch,  GA::GlobalArray * g_matrix )
{
    int status = -1;
    //GA::sync();

// Fetch and Check dra patches
   if ( dra_patch.size() != 2 || ga_patch.size() != 2 ) {
      cerr << "Unexpected size of the DRA or GA patches: Going to abort " << dra_patch.size();
      cerr << " GA patch is " << ga_patch.size() << endl;
      return(status);
   } 
   char handleMessage[] = "dumpGASectionDRA: ga handle invalid ";
   g_matrix->checkHandle( handleMessage );

//DRA patches 
   dra_size_t dlo[MAX_DRA_DIMS], dhi[MAX_DRA_DIMS];
   vector<pair<dra_size_t,dra_size_t> >::iterator dit;
   dit = dra_patch.begin();
   dlo[0] = (*dit).first-1;
   dhi[0] = (*dit).second-1;
   ++dit;
   dlo[1] = (*dit).first-1;
   dhi[1] = (*dit).second-1;
#if 0
   cout << "dilo and dihi are " << dlo[0] << " " << dhi[0] << endl;
   cout << "djlo and djhi are " << dlo[1] << " " << dhi[1] << endl;
#endif

//GA patches ... more uglyness.. the DRA argument list expects INTS for glo and ghi  

   int glo[MAX_DRA_DIMS], ghi[MAX_DRA_DIMS];
   vector<pair<int,int> >::iterator git;
   git = ga_patch.begin();
   glo[0] = (*git).first-1;
   ghi[0] = (*git).second-1;
   ++git;
   glo[1] = (*git).first-1;
   ghi[1] = (*git).second-1;

#if 0
   cout << "gilo and gihi are " << glo[0] << " " << ghi[0] << endl;
   cout << "gjlo and gjhi are " << glo[1] << " " << ghi[1] << endl;
#endif

// The first indeces should have the same size

   if ( ghi[0] > dhi[0] ) {
      cerr << "GA rows longer than DRA rows: Not necessarily an error " << endl;
   }

   if ( ghi[1] > dhi[1] ) {
      cerr << "GA columns longer than DRA columns: Not necessarily an error " << endl;
   }

 return( NDRA_Write_section(DRA_FALSE, g_matrix->handle(), glo, ghi, local_dra, dlo, dhi, &l_req) );
}

// Read a patch
int PsociDRArestart::readGASectionDRA( vector<pair<int,int> > ga_patch, vector< pair<dra_size_t,dra_size_t> > dra_patch,  GA::GlobalArray * g_matrix )
{
    int status=-1;
    //GA::sync();

// Fetch and Check dra patches
   if ( dra_patch.size() != 2 || ga_patch.size() != 2 ) {
      cerr << "Read: Unexpected size of the DRA or GA patches: Going to abort " << dra_patch.size();
      cerr << " Read: GA patch is " << ga_patch.size() << endl;
      return(status);
   } 
   char handleMessage[] = "readGASectionDRA: ga handle invalid ";
   g_matrix->checkHandle( handleMessage );

//DRA patches 
   dra_size_t dlo[MAX_DRA_DIMS], dhi[MAX_DRA_DIMS];
   vector<pair<dra_size_t,dra_size_t> >::iterator dit;
   dit = dra_patch.begin();
   dlo[0] = (*dit).first-1;
   dhi[0] = (*dit).second-1;
   ++dit;
   dlo[1] = (*dit).first-1;
   dhi[1] = (*dit).second-1;
#if 0
   cout << "dilo and dihi are " << dlo[0] << " " << dhi[0] << endl;
   cout << "djlo and djhi are " << dlo[1] << " " << dhi[1] << endl;
#endif

//GA patches ... more uglyness.. the DRA argument list expects INTS for glo and ghi  

   int glo[MAX_DRA_DIMS], ghi[MAX_DRA_DIMS];
   vector<pair<int,int> >::iterator git;
   git = ga_patch.begin();
   glo[0] = (*git).first-1;
   ghi[0] = (*git).second-1;
   ++git;
   glo[1] = (*git).first-1;
   ghi[1] = (*git).second-1;

#if 0
   cout << "gilo and gihi are " << glo[0] << " " << ghi[0] << endl;
   cout << "gjlo and gjhi are " << glo[1] << " " << ghi[1] << endl;
#endif

// The first indecies should have the same size

   if ( ghi[0] > dhi[0] ) {
      cerr << "GA rows longer than DRA rows: Not necessarily a read  error " << endl;
   }

   if ( ghi[1] > dhi[1] ) {
      cerr << "GA columns longer than DRA columns: Not necessarily a read error " << endl;
   }
 return( NDRA_Read_section(DRA_FALSE, g_matrix->handle(), glo, ghi, local_dra, dlo, dhi, &l_req) );
}

// Read a DRA of size D(n,m) into a GA of size (n,m)
// Collective Read a GA of size (n,m)  to a DRA of size (n,m);
int PsociDRArestart::readGAfromDRA( GA::GlobalArray * g_matrix )
{
    int status=-1;
    //GA::sync();

// Fetch and Check GA handle
   char handleMessage[] = "readGAfromDRA: ga handle invalid ";
   g_matrix->checkHandle( handleMessage );

   int ga_ndim;
   int dims[2];
   int type;
   g_matrix->inquire( &type, &ga_ndim, dims);
   if ( ga_ndim != 2 || ga_ndim != l_ndim ) {
      cerr << "readGAfromDRA: only a 2D GA can be read. Inquery reports GA dim = " << ga_ndim ;
      cerr << " and DRA dim is " << l_ndim << endl;
      cerr << " Did you remember to setdimension(new) or inquireDRAsetRead(restart)  ? "<< endl;
      cerr << " This could also mean the file doesn't exist/corrupted especially if DRA dim is large neg number " << endl;
      return(-1);
   }
   if ( type != C_DBL && type != C_INT && type != C_LONG ) {
      cerr << "Error: readGAfromDRA only handles types C_DBL or C_INT or C_LONG " << endl;
      exit(1);
      return(-1);
   }

   status = NDRA_Read(g_matrix->handle(), local_dra, &l_req);
   DRA_Flick();
   return ( status );
}
int PsociDRArestart::waitDRA()
{
    DRA_Flick();
    return( DRA_Wait( l_req ) );
}
int PsociDRArestart::openDRAforREAD()
{
    int status;
    //cout << "on the read filename is " << l_filename << endl;
    char * local_filename = new char[ l_filename.size()+1];
    strcpy (local_filename, l_filename.c_str() );
    //cout << "Calling the read " << " mode is " << l_mode << endl;
    status = DRA_Open(local_filename, l_mode, &local_dra );
    if ( status != 0 ) {
       cout << " DRA_Open is not zero: this is probably a non-existant file: " << status << endl;
    }
    return( status ); 
}
    
int PsociDRArestart::closeDRA()
{
    return( DRA_Close( local_dra ) );
}

int PsociDRArestart::deleteDRA()
{
    return( DRA_Delete( local_dra ) ) ;
}
void PsociDRArestart::flickDRA()
{
    DRA_Flick();
    return;
}

//Useful only on a restart when we need to set the dimensions
int PsociDRArestart::inquireDRAsetREAD()
{
    int status = -1;
    char filename[200];
    char name[80];
    int type, ndim;
    dra_size_t dims[2];

    status = NDRA_Inquire(local_dra, &type, &ndim, dims, name, filename);
    
    if ( ndim > MAX_DRA_DIMS ) {
       cerr << "Requested ndim > MAX_DRA_DIMS " << ndim << " MAX is " << MAX_DRA_DIMS << endl;
       return(status);
    }
    
    l_ndim = ndim;
    l_dim0 = dims[0];
    l_dim1 = dims[1];
    l_rdim0 = -1;
    l_rdim1 = -1;
#if 0
    cout << "DRA inquire parameters are "<< l_ndim << " " << l_dim0 << " " << l_dim1 << endl;
#endif
    status = 0;
    return( status );
}

   
     
