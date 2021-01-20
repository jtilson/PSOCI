/**************************************************************************************

* Copyright (c) 2010,2011,2012 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: A class to process the current set of results into populations analysis
              One population analysis is performed for each provided CI vector
              As part of this method we include
		1) Reading MO coefs
		2) writing NO coefs in the MO basis(the CI density)
		3) Writing NO coefs in the AO basis
		4) Mulliken population analysis for each root
    		5) Optionally (TODO) a density computation over R^3 suitable for visualization



 History:

**************************************************************************************/
/**
 *   @file PsociNOpopulations.C
 *
 */

#include <list>
#include "PsociBlasLapack.hpp"
#include "PsociNOpopulations.hpp"

#define NNDAF(i) (i * (i + 1)) / 2 

using namespace std;

PsociNOpopulations::PsociNOpopulations(int numRoots, GA::GlobalArray * g_v, PsociGAhamiltonian * hamilton )
{
  l_hamilton = hamilton;
  l_nbf = l_hamilton->fetchNBF();
  l_numroots = numRoots;
  l_maxsef = hamilton->fetchNumberSefs();
  l_g_v = g_v;
  l_nelec = hamilton->fetchNELEC();
}

void PsociNOpopulations::selectMOs( string title, int mo_unit, string mo_filename )
{
  if ( GA::nodeid() != 0 ) return;
  cout << "Allocate Molecular Orbital object for "<< title<<endl;
  cout << mo_unit<<" "<<mo_filename<<endl;
  l_mocoef = new PsociProcessMOcoefs( l_nbf, mo_unit, mo_filename );
}

void PsociNOpopulations::selectNOs( string title, int no_unit, int numroots, string no_filename )
{
  if ( GA::nodeid() != 0 ) return;
  cout << "Allocate Natural Orbital object for "<< title<<endl;
  cout << no_unit<<" "<<"Base filename is "<<no_filename<<endl;
  l_nocoef = new PsociProcessMOcoefs( l_nbf, no_unit, numroots, no_filename );
}

void PsociNOpopulations::removeMOs()
{
  if ( GA::nodeid() != 0 ) return;
  delete l_mocoef;
}
void PsociNOpopulations::removeNOs()
{
  if ( GA::nodeid() != 0 ) return;
  delete l_nocoef;
}


PsociProcessMOcoefs * PsociNOpopulations::fetchMOcoefs()
{
  return( l_mocoef );
}

PsociProcessMOcoefs * PsociNOpopulations::fetchNOcoefs()
{
  return( l_nocoef );
}

/* READ coefficients (etc) from a COlumbus formatted coef file
*/
int PsociNOpopulations::readMOs()
{
  if ( GA::nodeid() != 0 ) return(0);
  
  l_mocoef->OpenMOfile();
  int numrd = l_mocoef->MOcoefficientRead();
  
  l_mocoef->printMOparameters();
  l_mocoef->CloseMOfile();
  
#ifdef DETAILEDCHECK
  l_mocoef->printMOcoefficients();
#endif
    return( numrd );
}

/* Take header information existing in the "l_" arrays. Unpack them
   into a format suitable for pasing to F77 routines. Then use the
   Columbus mowrit method to write the headers. Writing the coefs will come later

   Space permitting add a new header to the title list ELSE overwrite the
   last (10th) entry with the new one
*/

void PsociNOpopulations::writeNOheaders( string title )
{
  if ( GA::nodeid() != 0 ) return;

  l_nocoef->printMOparameters();
  l_nocoef->WriteNOheaders( title );
}


/* Open a set of files for output
*/
void PsociNOpopulations::openNOs()
{
  if ( GA::nodeid() != 0 ) return;
  l_nocoef->OpenNOfiles();
}

void PsociNOpopulations::closeNOs()
{
  if ( GA::nodeid() != 0 ) return;
  l_nocoef->CloseNOfiles();
}


/* Timed entry point */

long PsociNOpopulations::populationDensityMatrix( pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  long num = populationDensityMatrix();
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return( num );
}


/* entry point to the main CI density code
*/

long PsociNOpopulations::populationDensityMatrix()
{
  int num = constructDensityMatrix( l_numroots, l_g_v, *l_mocoef );

   if ( GA::nodeid() == 0 ){
      GA::printStats();
   }
  return( num ); 
}


/* Take the existing set of NOs in the MO basis ( for the given symetry ) and 
   transform to NOs in the AO basis
   
   dens_vectorAO_sym gets set to numorbs*numorbs: [AO][NO]
   input: dens_vector_sym : [MO][NO] these are already diagonalized and normalized on I.
*/

// I bet this is wrong even though it gives the same answers as Columbus::density.x !

void PsociNOpopulations::transformNOstoAObasis( int isym, vector<double> & no_occupations, vector<double> & dens_vector_sym, vector<double> & dens_vectorAO_sym, PsociProcessMOcoefs & coefs )
{
  if ( GA::nodeid() != 0 ) return;
  
  int nsym = coefs.fetchNsym();
  
  vector<int> nstart(nsym,0);
  vector<int> nsize(nsym,0);
  /* 
     Compute entry points for density(ij) style array
     nstart[isym] return the absolute i to index on.
  */
  //Build offsets etc
  int noffset = 0;
  for(int jsym=0; jsym< nsym; ++jsym) {
    nstart[jsym] = noffset;
    nsize[jsym] = coefs.fetch_nmopsy(jsym);
    noffset += coefs.fetch_nmopsy(jsym) ;
  }
  int noffset_mos=0;
  int noffset_density=0;
  
  for(int i=0; i< isym; ++i ) {
    int ntemp = nsize[i];
    noffset_mos += ntemp * ntemp;
    noffset_density += NNDAF( ntemp );
  }
  int numbasis = nsize[isym];
  
  if ( dens_vectorAO_sym.size() != numbasis * numbasis ) dens_vectorAO_sym.resize( numbasis * numbasis ); // Will overwrite anything 
  
  /* In the below we perform excessive numbers of memory copies. Moving forward we may go ahead
     and combine work. But for now, because the data objects are packed and blocked and generally
     confusing ( F77 style ordering, etc.) this version is easy to read
  */
   vector<double> row_scratch( numbasis, 0.0 );
   
   const int ld = 1;
   
   int index_mos;
   int index = 0;

   for (int i=0; i< numbasis; ++i ) {
     
     vector<double> vec_scratch( nsize[isym], 0.0 );
     for(int k=0; k< numbasis; ++k ) { // indexes MO
       vec_scratch[k] = dens_vector_sym[ i*numbasis + k]; // ith vector over all curent MOs
     //  cout <<"vec_scratch "<<  k << " " << vec_scratch[k]<< endl;
     }
     
     for( int j=0; j< numbasis; ++j ) {
       vector<double> mo_scratch( nsize[isym], 0.0 ); // densMO is packed and triangular
       index_mos = j;
       for(int k=0; k< numbasis; ++k ) {              // indexes MO
         mo_scratch[k] = coefs.fetchMocoefs( index_mos + noffset_mos  );
         index_mos += numbasis;
       }
       row_scratch[j]  = ddotPsoci(numbasis, mo_scratch, ld, vec_scratch, ld );
     }
     for(int k =0; k< numbasis; ++k ) {
       dens_vectorAO_sym[index + k]  = row_scratch[k];
     }
     
#ifdef DETAILEDCHECK
    for(int k =0; k< numbasis; ++k ) {
      cout << "specifics " << dens_vectorAO_sym[index + k] << endl;
    }
#endif
    
    //Push back into NO array
    index+= numbasis;
   }

/* Print out NOs in the AO basis as a check
   For the test case these are identical those from Columbus::density.x. 
*/
   string title="NO in the AO basis ";
   outputOccupationsSetvectors( title, isym, no_occupations, dens_vectorAO_sym );
   
   return;
}
/* The CI eigenvectors need to be sorted for subsequent methods
   such as the NOtoAO transformation.  The eigenvectors and values
   must conform in the order. Zero valued vectors must be retained
   
   Degenerate values for now are randomly ordered. In the future we may need 
   to subsort based on coefficients
   
   Use a simple multimap to sort on keys of eigenvalues. Then push data
   back into the evals annd evec arrays
   
   input data are manipulated on return
*/

//TODO double check that the sort operations are okay, e.g., do we need sort by abs(eval)?

void PsociNOpopulations::sortNOinMObasis(int isym, vector<double> & evals, vector<double> & vectors )
{

  if ( GA::nodeid() != 0 ) return;
  cout << "Sort NOs in the MObasis by Evalue. Largest Magnitude first" << endl;
  int numevals = evals.size();
  if ( numevals <= 0 ) GA::error("sortNOinMObasis:bad evals size ",numevals );
  if ( vectors.size() != numevals*numevals ) GA::error("sortNOinMObasis: incompatible vector size ",numevals );

  //cout << "numevals is " << numevals << endl;

  
  list< pair<double, vector<double> > > sort_buffer;

  
  // Not the best way but very easy to follow and not a big deal overal (w.r.t memory or flops )
  
  for(int i=0; i< numevals; ++i ) {
     double dkey = evals[i];
     vector<double> buffer( numevals,0.0);
     for(int k=0; k< numevals; ++k ) {
       buffer[k] = vectors[i*numevals+k];
     }  
     sort_buffer.push_back( pair<double,vector<double> >(dkey, buffer) );
  }  
  
// Now unpack back into two vector arrays

  evals.clear();
  vectors.clear();
  
  sort_buffer.sort( compare_evalues );
  
  list< pair<double,vector<double> > >::iterator it;
  vector<vector<double> > vec_scratch;
  
  for( it = sort_buffer.begin(); it != sort_buffer.end(); ++it ) {
    double dkey = (*it).first;
    evals.push_back( dkey );
    vec_scratch.push_back( (*it).second );
  }
#ifdef DETAILEDCHECK
  cout << "SIZE of pushed vector " << vec_scratch.size() << endl;
  cout << "SIZE of pushed evals is " << evals.size() << endl;
#endif
  
  // Now take set of scratch vectors and condense into a single vector 
  // evals is already sorted
  
  vector<vector<double> >::iterator vit;
  vectors.resize(numevals*numevals,0.0); // for output
  
  int index = 0;
  for(vit = vec_scratch.begin(); vit != vec_scratch.end(); ++vit) {
    vector<double> scratch = (*vit);
    for(int k=0; k< numevals; ++k ) {
      vectors[index] = scratch[k];
      ++index;
    }
  }
  
#ifdef DETAILEDCHECK
  for(int i=0; i< numevals; ++i ) {
    cout << "sort: eigenvalues after sort " << evals[i]<<endl;
  }
#endif
}

// Sort from highest to lowest 
bool compare_evalues( pair< double,vector<double> > left, pair< double,vector<double> > right )
{
  double dleft = left.first;
  double dright = right.first;
  if (dleft > dright ) return( true );
  return ( false );
}

/* Construct the CI density matrix in the MO basis.  
   nsize() is used to calculate offsets that allow mergin of symetry blocked (packed) data sets (uhg)
  
   on exit the array densityMO is populated (packed format) for the CI elements: NOTE
   do not resize this array at this time
*/

void PsociNOpopulations::computeMOdensity( int isym, vector<double> & densityMO, vector<double> densityCI, PsociProcessMOcoefs & coefs )
{
  if ( GA::nodeid() != 0 ) return;
  double dZero = 0.0;
   double dtemp;
   
   // Compute our current offset
   
   int noffset_mos = 0;
   int noffset_density = 0;
   
   int ntemp;
   
   int nsym = coefs.fetchNsym();
   
   vector<int> nstart(nsym,0);
   vector<int> nsize(nsym,0);
   
   /* 
      Compute entry points for density(ij) style array
      nstart[isym] return the absolute i to index on.
   */
   
   int noffset = 0;
   for(int jsym=0; jsym< nsym; ++jsym) {
     nstart[jsym] = noffset;
     nsize[jsym] = coefs.fetch_nmopsy(jsym);
     noffset += coefs.fetch_nmopsy(jsym) ;
   }
   
   for(int i=0; i< isym; ++i ) {
      ntemp = nsize[i];
      noffset_mos += ntemp * ntemp;
      noffset_density += NNDAF( ntemp );
   }
   
   int numbasis = nsize[isym];
   
   vector<double> dens_buffer_sym( numbasis*numbasis, dZero );
   vector<double> scratch( numbasis*numbasis, dZero );
   int index = 0;
   for(int i=nstart[isym]; i< nstart[isym]+numbasis; ++i ) {
     for(int j=nstart[isym]; j< nstart[isym]+numbasis; ++j) {
//       for(int j=nstart[isym]; j<=i; ++j) {
       int ij = twoDindex(i,j);
       dens_buffer_sym[index] = densityCI[ij];
       ++index;
     }
   }
   
   // index_mos = noffset_mos; // start MOs at the next symmetry block - normal square packing
   // Perform first half transform for the current root and symmetry
   // D(i,k)*coef(k,j) = buf(mo,basis)
   
   int index_scratch = 0;
   int index_mos;
   int index_densao;
   int index_density;
   
// convert this to something more clear

   for (int i=0; i< numbasis; ++i ) {
     for(int j=0; j< numbasis; ++j ) {
       dtemp = dZero;
       index_mos = j;
       for(int k=0; k< numbasis; ++k ) { // loop over orbitals jth entry in each k orb
	 index_density = i*numbasis+k;
	 dtemp += dens_buffer_sym[index_density] * coefs.fetchMocoefs( index_mos +noffset_mos  );
        // cout << "xx COEF for this set " << coefs.fetchMocoefs( index_mos +noffset_mos  ) << endl;
	 index_mos += numbasis; // selects jth basis function from orbital k
       }
       //cout << "partial dens " << j<<" "<<i<<" "<<dtemp << endl;
       scratch[index_scratch] = dtemp;
       ++index_scratch;
     }
   }
   
   // Perform second half transform for the current root and symmetry
   // coef(i,k) * buf(j,k) = dens(basis,basis)
   
   index_densao=0;
   for (int i=0; i< numbasis; ++i ) {
     for(int j=0; j<=i; ++j ) {
       dtemp = dZero;
       index_mos = j;
       index_scratch = i;
       index_densao = twoDindex(i,j);
       
       for(int k=0; k< numbasis; ++k ) {
	 dtemp += scratch[index_scratch] * coefs.fetchMocoefs( index_mos +noffset_mos  );
         index_mos += numbasis; // selects jth basis function from orbital k
         index_scratch += numbasis;
       }
       densityMO[index_densao+noffset_density] = dtemp;
     }
   }
   return;
}

/* Diagonalize the input CI density (array) and return NOs in the MO basis and associated occupations 
   So far it is clear the NOs in the MO basis are orthogonal consistent with dsygv documentation
   The strange population/norm problem smust come after the NOtoAO transformation
*/
int PsociNOpopulations::diagonalizeNOinMObasis( int isym, int nstart, int numorbs, vector<double> & no_occupations, vector<double> & dens_vector_sym, vector<double> & array ) {
  if ( GA::nodeid() != 0 ) return(0);
  cout << "Diagonalize the CI density matrix " << endl;
  
  int index = 0; 
  for(int i=nstart; i< nstart+numorbs; ++i ) {
    for(int j=nstart; j< nstart+numorbs; ++j) {
      int ij = twoDindex(i,j);
      dens_vector_sym[index] = array[ij];
      ++index;
    }
  }
  
  int itype = 1; // A*x = B*x*lambda
  char jobz = 'V'; 
  char uplo = 'L'; //Lower triangle
  int lwork = 50 * numorbs; //TODO Need to optimize this
  int info=-1;
  
  // Create a Unit metric  
  
  vector<double> unit_matrix_t( numorbs * numorbs, 0.0 );
  index=0;
  for (int i=0; i< numorbs; ++i ) {
    for (int j=0; j< numorbs; ++j ) {
      if ( i == j ) unit_matrix_t[index] = 1.0; 
      ++index;
    }
  }
  
/*
   1:	A*x = (lambda)*B*x
   The	eigenvectors are normalized as follows:	if ITYPE = 1, Z**T*B*Z = I
*/
  vector<double> work(numorbs*numorbs, 0.0);

  dsygvPsoci( itype, jobz, uplo, numorbs, dens_vector_sym , numorbs, unit_matrix_t, numorbs, no_occupations, work, lwork, info );

  if ( info != 0 ) GA::error("PsociNOpopulations::dsygv info not zero ",info);

  return( info );
}


/* OUTPUT the NOs in the MO basis including evals    On entry vectors is SQUARE with an LDA of numevals   dsygv return the VECTORS as columns.*/ 

void PsociNOpopulations::outputOccupationsSetvectors( string title, int isym, vector<double> evals, vector<double> vectors )
{
  if ( GA::nodeid() != 0 ) return;
  
  int numevals = evals.size();
  if ( numevals <= 0 ) GA::error("bad evals size ",numevals );
  if ( vectors.size() != numevals * numevals ) GA::error("incompatible vector size ",numevals );
  
  /* Since we now have an explicit sort method DO NOT try to reorder here
   */
  
  cout << "eigenvectors : Symmetry block " << isym << endl;
  cout << "display occupations >= "<< setprecision(10)<<MIN_DOUBLE_EVAL_PRINT<<endl<<endl;
  cout << title << endl;
  
#ifdef DETAILEDCHECK
  for(vit=evals.begin(); vit!=evals.end(); ++vit) {
    cout << "NO Occupation " << (*vit) << endl;
  }
#endif
  
  int index=0;
  
  /* Chemists like to read important NOs first in the list
     Also we need to explicitely order the roots (by energy) for subsequent methods
     such as the NOtoAO transformations
  */
  
  for(int i=0; i< numevals; ++i ) {  
    if ( evals[i] >= MIN_DOUBLE_EVAL_PRINT ) {  
      cout << "Eigenroot: occupation is " << evals[i] << endl;
      for(int j=0; j< numevals; ++j ) { // indexes MO
	cout << vectors[i*numevals+j] << " "; 
	++index;
      }    
      cout << endl << endl;
    }
  }
}

void PsociNOpopulations::constructNaturalOrbitals( bool square, int nbf, int nroots, vector<vector<double> > & array, PsociProcessMOcoefs & coefs,pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  constructNaturalOrbitals( square, nbf, nroots, array, coefs );
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
}

/* 
   Take the input ARRAY and construct NOS for each CI root. We require eigenvectors and eigenvalue (occupations ) for all
   We stipulate that all roots were constructed using the same basis of MOs
   
   The MOs are assembled into symmetry blocks so this complicates the assembly a little bit
   For now we want to fold the AO coefficients into the full bfn array. Since the wavefunction
   
   the full H method has such a high complexity it is hard to imagine the nbf even getting too big.
   
   on entry is the density over MOs (array) stored as packed: array( nroots, nbf*(nbf+1)/2 )
   Only need to diagonalize symmetry blocks to get NOs and occupations
   
   One node computation: it is pretty fast 
   
   array == dij
*/

void PsociNOpopulations::constructNaturalOrbitals( bool square, int nbf, int nroots, vector<vector<double> > & array, PsociProcessMOcoefs & coefs )
{
  if ( GA::nodeid() != 0 ) return;
  
  cout << "Sequential processing: Constructing Natural Orbitals for " << nroots<<" eigenroots " << endl;
  
  if ( nroots != array.size() ) cout << "Requested number of roots < array.size() " << array.size() << endl;
  
  const double dZero = 0.0;
  int nsym = coefs.fetchNsym();
  
  // The subsequent Mulliken analysis requires the full (but packed) densAO array. Not a symmetry blocked version
  // Construct matrices over roots then symmetry blocks. 
  
  // mocoef stored as (numbasis)(nummos) 
  
  vector<vector<double> > densAO(nroots);
  
  for(int iroot=0; iroot<nroots; ++iroot) {
    cout << "------------------------------------------" << endl;
    cout << "Processing NOs for Root                   " << iroot << endl;
    
// is failed already
// array is empty here at iroot = 1
#ifdef DETAILEDCHECK 
    string title1="CI density title";
    outputPacked2Dmatrix( false, nbf, array[iroot], title1 ); // print for all roots
#endif
    
    /* 
       Compute entry points for density(ij) style array
       nstart[isym] return the absolute i to index on.
    */

    vector<int> nstart(nsym,0);
    vector<int> nsize(nsym,0);
    int noffset = 0;
    for(int isym=0; isym< nsym; ++isym) {
      nstart[isym] = noffset;
      nsize[isym] = coefs.fetch_nmopsy(isym);
      noffset += coefs.fetch_nmopsy(isym) ;
    }
    
    vector<double> scratch_aos( NNDAF(nbf) , dZero );
    
    int real_index = 0; // Used towards the end for packing up coefs for printout.
    int eval_index = 0; // ditto for the evals

    vector<Fdouble> packed_coefs( l_nbf*l_nbf, dZero ); // These get written to disk via F77 routines
    vector<Fdouble> packed_evals( l_nbf, dZero );

    for(int isym=0; isym< nsym; ++isym) {
      
      int numorbs = coefs.fetch_nmopsy(isym);
      int numbasis = coefs.fetch_nbfpsy(isym);
      
      // Not sure if we can treat this case so bail if uncertain
      if ( numorbs != numbasis ) {
	GA::error("constructNaturalOrbitals: numorbs != numbasis ",-1);
      }
      /* 
	 Compute NOs/Occupations in the MO basis for the current symmetry block
	 Push CI density info square form then diagonalize. Follow this with
	 a reverse orbering of the results to output
      */
      vector<double> dens_vector_sym( numorbs * numorbs, dZero );
      vector<double> no_occupations( numorbs, dZero );

      int info = diagonalizeNOinMObasis( isym, nstart[isym], numorbs, no_occupations, dens_vector_sym, array[iroot] );
      
      sortNOinMObasis(isym, no_occupations, dens_vector_sym );
      
      string title="NOs in the MO basis";
      outputOccupationsSetvectors( title,  isym, no_occupations, dens_vector_sym );
      /* 
	 Compute the density matrix in the AO basis for the current symmetry. Ultimately we need all symmetries in the same object
	 Only compute the triangle elements and pack into the densAO memory space
	 This object is desired for generation of the Mulliken population analysis
      */

      computeMOdensity( isym, scratch_aos, array[iroot], coefs );
      vector<double> dens_vectorAO_sym( numorbs*numorbs, dZero );
      transformNOstoAObasis( isym, no_occupations, dens_vector_sym, dens_vectorAO_sym,coefs );

/* Grab current set of occupations(isym) and vectors(isym) and pack them into a new object that
   will get printed to an NO file after the last symmetry is processed for root = iroot
   Eigenvectors/eigenvalues are already sorted 
*/

  int local_index=0;
  for(int i=0; i< numorbs; ++i ) {
         packed_evals[ eval_index ] = no_occupations[ i ];
      for(int j=0; j< numorbs; ++j ) {
         packed_coefs[ real_index ] = dens_vectorAO_sym[ local_index ]; 
        ++local_index;
        ++real_index;
      }
    ++eval_index;
    }

    } // end symmetry

    int dens_packed_size=0;
    for(int isym=0; isym< nsym; ++isym) {
      dens_packed_size += NNDAF( coefs.fetch_nmopsy(isym) );
    }

    densAO[iroot] = scratch_aos; // so far so good.  it should be triangle and symmetry blocked
    densAO[iroot].resize( dens_packed_size ); // Will remove zeros from the back

#ifdef DETAILEDCHECK
    cout << "densAO total size for root is " << iroot<<" "<<densAO[iroot].size()<<endl;
    string title3="densAO scratch title";
    outputPacked1Darray(dens_packed_size, densAO[iroot], title3 ); // print for all roots
#endif

/* dump NOs for the current root to the nocoef files set. By construction we have a sufficient number of files
   prepared for writting. Simply pass the IROOT index to the method to write NOs for the IROOTth CI vector.
*/

   string stitle="PSOCI CI Natural Orbitals and occupations";
   l_nocoef->OpenNOfiles(iroot ); 
   int numcoefs = l_nocoef->WriteNOcoefs( stitle, iroot, packed_evals,  packed_coefs );
   l_nocoef->CloseNOfiles(iroot );

/* Perform population analyusis for each root: densAO[iroot] is trangle AND symmetry blocked
*/

#ifdef DETAILEDCHECK
 cout << "check MO normalization " << endl;
   checkOrthonormality(l_nbf, coefs   );
 cout << "end check MO normalization " << endl;
#endif

 cout << "Call Mulliken Population Analysis on the 1e Density Matrix" << endl;
 cout << "----------------------------------------------------------" << " CI root is " << iroot << endl;
 cout << endl;
   constructMullikenPopulationDensAO(l_nbf, iroot, densAO[iroot], coefs);
 cout << endl;
 cout << "Finished with ---------" << " CI root is " << iroot << endl << endl;


/*
 cout << "Call Mulliken Population Analysis for Occupied Orbitals" << endl;
 cout << "-------------------------------------------------------" << " CI root is " << iroot << endl;
   constructMullikenPopulationNOs(l_nbf, packed_evals,  packed_coefs , coefs);
 cout << endl;
 cout << "Finished with ---------" << " CI root is " << iroot << endl << endl;
*/

  } // next root
  
  return;
}


/* Only rank==0 can print data
   square==false: array should be of the form (n*(n+1)/2)
   square==true:  array should be of the form (n*n)

*/

void PsociNOpopulations::outputPacked2Dmatrix(bool square, int nbf, vector<double> & array, string title )
{
  if ( GA::nodeid() != 0 ) return;
  
  cout << "Printing array "<< title<<endl;
  int length = array.size();
  if ( !square )  {
    if ( NNDAF(nbf)  != length ) {
      GA::error("outputPacked2Dmatrix: unexpected packed array length ",length );
    }  
  } else {
    if ( nbf*nbf != length ) {
      GA::error("outputPacked2Dmatrix: unexpected square array length ",length );
    }  
  }

  /*   
       double sum = accumulate( array.begin(), array.end(), 0.0 );
       cout << "accumulate is " << sum <<endl; 
  */
  
  if ( square ) {
    
    cout << title<<endl;
    cout << "Square format " << endl;
    
    for(int i=0; i< nbf; ++i ) {
      cout << "Row "<<i<<endl; 
      for(int j=0; j< nbf; ++j) {
        int index  = i * nbf + j; 
        
        cout << array[ index ]<<" ";
      }
      cout << endl;
    }
    
  } else {
    
    cout << title<<endl;
    cout << nbf<< " packed format " << title.size()<< endl;
    for(int i=0; i< nbf; ++i ) {
      cout << "Row "<<i<<endl;
      for(int j=0; j<=i; ++j) {
        cout << array[ twoDindex(i,j) ]<<" ";
      }
      cout << endl;
    }
  }
}

/* 
   A method to simply dump the total contents of a vector to stdout
   This is used mostly for debugging for now. As an example densAO[nroot][]
   arranges values in memory as a succession of packed triangles.
*/

void PsociNOpopulations::outputPacked1Darray(int nelem, vector<double> & array, string title )
{
  if ( GA::nodeid() != 0 ) return;
  if ( nelem > array.size() ) cout << "outputPacked1Darray: warning. Array not that big " << array.size() <<endl; 
  double sum=0.0;
  for(int i=0; i< nelem; ++i ) {
    cout << i << " "<<array[i] << endl;
    sum += array[i];
  }  
}   

/* 
   A method to simply dump the AO density to the output
   This method can resized and unresized densAO arrays.
   nstart[nsym] report the absolute value of I the begins the next symmetry block
   nsize is an array that reports numbasis per symmetry
      by construction numorbs = numbasis else we would have stopped by now
*/

void PsociNOpopulations::outputPrintDensityAObasis( int nsym, vector<int> & nsize, vector<int> & nstart, vector<double> & density, string title ){
  
  if ( GA::nodeid() != 0 ) return;
  cout << "Print out density matrix for "<<title<<endl;
  
  int index = 0;
  for(int isym = 0; isym < nsym; ++isym ) {
    int istart = nstart[isym];
    int numbasis = nsize[isym];
    
    for(int i = istart; i < istart+numbasis; ++i ) {
      for(int j = istart; j <= i; ++ j ) {
	cout << i<<" "<<j<<" "<<density[index] << endl;
	++index;
      }
    }  
  }  
  
  return;
}

/* Values are C-style and begin at 0, ID computes in a Fortran style way
   Canonical indexing of a 2D array
*/

inline int PsociNOpopulations::twoDindex( int mC, int nC ) {
  int m = mC + 1;
  int n = nC + 1;
  int value = ((max(m,n))*((max(m,n))-1))/2+min(m,n);
  return( value - 1 ); // convert back to C-style
}

/* Main entry to the DENSITY code. THis is a parallel code and all must enter.
   The primary results are:
       1) Density over the AOs for each root.
       2) Natural orbitals
       3) Mulliken population of various degrees. ( Gross, net, atomic, orbital, fragment, etc.)

      On entry provide the number of Roots. This correspnds to the number of roots in the provided
      GA that are of actualy desire. Note: GA is likely much bigger as it could be allocated to include
      the full subspace size.

      THe calling app MUST ensure the numRoots are the converged roots of interest. One could use
      PsociGAvectors method to read off disk into GA space
  
      DENSITY (the CI density) can be replicated since it is only of size ( roots x ( nbf*nbf) ). This will
      not get very big for now. We can always do a distributed version later. 
*/ 

/* This is a PARALLEL method in this class. It requires coefs to be replicated
*/

long PsociNOpopulations::constructDensityMatrix( int numRoots, GA::GlobalArray * g_v, PsociProcessMOcoefs & coefs )
{  
  double dZero = 0.0;
  
  int g_rank = GA::nodeid();
  
  if ( g_rank == 0 ) {
    cout << endl;
    cout << "------------------------------------------------------------------"<<endl;
    cout << "I am " << GA::nodeid() << " Computing the Density Matrix " << endl;
    cout << "Number of roots to process is "<<numRoots<<" number of AO BasisFtns is "<<l_nbf<<endl;        
    cout << "------------------------------------------------------------------"<<endl;
    cout << endl;
  }
  
  // Check state of the vector 
  char * vectorMessage = "PsociNOpopulations:constructDensityMatrix:Vector";
  g_v->checkHandle( vectorMessage );
  
  int type, ndim, dims[2];
  g_v->inquire( &type, &ndim , dims );
  
  //   g_v->print();
  
  int maxsefs = dims[0];
  int maxroots = dims[1];
  
  if ( numRoots > maxroots ) {
    GA::error("Density code: g_v smaller than specified roots: maxroots is ", maxroots );
  } 
  if ( maxsefs != l_maxsef ) {
     GA::error("Density code: g_v reports in appropriate maxsefs", maxsefs );
  }

/* Open up DENSITY array. Each root is of size nbf*(nbf+1)/2. This will be processed in parallel. So a 
   GOP operation will be forthcoming.
*/
// Compute and store into GA, the CI Hamiltonian appropriate for a 1e Density. That is no excitations > 2, etc.

 pair<int,double> CIDinfo;

/* Not needed since CIdensity is also computing them
 long numHelems = l_hamilton->constructGlobalHamiltonianForCIDEN2WAY(CIDinfo );

 if ( GA::nodeid() == 0 ) {
         cout << endl <<"1e H: Total number of non-zero elements is " << GA::nodeid()<<" "<<numHelems << endl;
         cout << "Timing:1e H: Me is " << CIDinfo.first << " time is "<< CIDinfo.second << endl << endl;
   }
*/
//   Now construct the Density Matrix
//   This is the 1d DENSITY over MOs : sized correctly

 vector<vector<double> > density( numRoots, vector<double>( NNDAF(l_nbf), dZero ) );

#ifdef ALTERNATIVEDENSITY
//long numHelems = l_hamilton->constructCIdensityAlt( numRoots, g_v, density, CIDinfo  );
  GA::error(" DO NOT compile with ALTERNATIVEDENSITY anymore " ,1);
#endif

// both give the same error
#ifdef DIRECT 
// Need to build a distributed Jchunk version of the CI density

#ifdef NEWORBMAP
if ( GA::nodeid() == 0 ) cout << "Replicated Orbmap CI density " << endl;
long numHelems = l_hamilton->constructCIdensityAltOrbMapDirect(numRoots, g_v, density, CIDinfo  );
#else
if ( GA::nodeid() == 0 ) cout << "New CI ChunkJmap density " << endl;
long numHelems = l_hamilton->constructCIdensityAltOrbMapDirectChunkJmap(numRoots, g_v, density, CIDinfo  );
#endif

#else
//long numHelems = l_hamilton->constructCIdensityAltOrbMap(numRoots, g_v, density, CIDinfo  );
if ( GA::nodeid() == 0 ) cout << "New CI ChunkJmap NONDIRECT density " << endl;
long numHelems = l_hamilton->constructCIdensityAltOrbMapChunkJmap(numRoots, g_v, density, CIDinfo  );
#endif

// empty to here
// uncomment if you want to print the CI density matrix
/*
 if ( g_rank == 0 ) {
   cout << "BEFORE GOP SIZE of density in MO basis is " << density.size() << endl;
   bool square=false;
   for(int i2=0; i2< numRoots; ++i2 ) {
     cout << "PREGOP DUMP ROOT is " << i2<<endl;
     string title="ciden MO test title";
     outputPacked2Dmatrix( square, l_nbf, density[i2], title ); // print for all roots
   }
 }
*/

/*
#ifdef TWOWAY
#ifdef ALTERNATIVEDENSITY
// long numHelems = l_hamilton->constructCIdensityAlt( numRoots, g_v, density, CIDinfo  );
long numHelems = l_hamilton->constructCIdensityAltOrbMap(numRoots, g_v, density, CIDinfo  );
#else
 long numHelems = l_hamilton->constructCIdensity( numRoots, g_v, density, CIDinfo  );
#endif
#else
cout << "constructCIdensityAltOrbMap " << endl;
long numHelems = l_hamilton->constructCIdensityAltOrbMap(numRoots, g_v, density, CIDinfo  );
 //long numHelems = l_hamilton->constructCIdensityAlt( numRoots, g_v, density, CIDinfo  );
#endif
*/

 if ( GA::nodeid() == 0 ) { 
        cout << "Timing:CIdensity: Me is " << CIDinfo.first << " time is "<< CIDinfo.second << endl << endl;   }

 bool square=false;
 pair<int,double> Repinfo;
 replicateCIdensity( square, l_nbf, numRoots, density, Repinfo ); // Yes because parallel CI constructs in blocks 
 if ( GA::nodeid() == 0 ) cout << endl << "Time to replicate CI density is " << Repinfo.second<<endl ;

 
/* uncomment if you want to print the CI density matrix
*/
// already broken at here

 // We might want to consider replicating coefs since only node zero has them  
 
 pair<int,double> Popinfo;
 if ( g_rank == 0 ) {
cout << "before constructNOs "<< l_nbf<<" "<<numRoots<<" "<<density.size()<<endl;

         constructNaturalOrbitals(square, l_nbf, numRoots, density, coefs, Popinfo ); 
         cout << endl << "Time to construct NOs and populations is " << Popinfo.second<<endl;

 } 
 
 return( numHelems );
}

void PsociNOpopulations::replicateCIdensity(bool square, int nbf, int nroots, vector<vector<double> > & array, pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  replicateCIdensity(square, nbf, nroots, array );
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
}


/* In the construction of density, we use a distributed SEFSEWFtoMO loop
   
   we loop over all vectors and perform a GOP +

   Note also that the density COULD be square OR packed.

   input ARRAY is modified on output
*/

/* This is a PARALLEL method: it replicats the current CI density in the MO basis
   Not sure anymore why this is needed since all the other NO methods are single processor.
*/
void PsociNOpopulations::replicateCIdensity(bool square, int nbf, int nroots, vector<vector<double> > & array )
{

#ifdef DETAILEDCHECK
  if ( GA::nodeid() == 0 ) {
    cout << "Assembling/scatter array to all nodes " << endl;
  }       
#endif

  if ( nroots > array.size() )  {
    GA::error(" replicateCIdensity: unexpected array size ",array.size() );
  }       
  
  int expected_length; // in doubles
  ( square ) ? expected_length = nbf*nbf : expected_length = NNDAF(nbf) ;
  
  // We loop over vectors because  only the inner loop is contiguous
  
  char * op = "+";
  long global_Elems  = expected_length;
  
  //cout << "expected_length " << global_Elems << endl;
  for(int i=0; i< nroots; ++i ) { // Not the most efficient way but this wil never be a big deal either.


    GA::gop( &array[i][0], global_Elems, op );

/* irroot > 0 empty to here
    if ( GA::nodeid() == 0 ) {
    string title1="CI density title";
    cout << "CI density root is " << i << endl;
    outputPacked2Dmatrix( false, nbf, array[i], title1 ); // print for all roots
    }
    GA::sync();
*/

  }
}

/* The intention here is to copy headers from the INPUT orbital file(mocoef) and
   reproduce them in the output (set) of NO orbitals file(s).

   Title (metadata) will be updated at the writwe stage
*/

void PsociNOpopulations::fetchOrbitalFileHeader()
{
 if ( GA::nodeid() != 0 ) return;

   int nsym;
   vector<string> title;
   vector<string> labels;
   string afmt;
   vector<Fint> nbfpsy;
   vector<Fint> nmopsy;

   int ntitle = l_mocoef->copyOutHeaderInfo( nsym, title, labels, afmt, nbfpsy, nmopsy );
   int num = l_nocoef->copyInHeaderInfo( nsym, ntitle, title, labels, afmt, nbfpsy, nmopsy );
}

/* Construct Mulliken populations for all roots. Prepare to dump all results to stdout.
   coefs is needed mostly to provide symmetry blocking information

   densAO is symetry blocked, triangular packed.
   overlap comes from aoints as triangular packed but not symmetry blocked.

  DO only one densAO root at a time
  densAO comes in as symmetry packed and triangular

*/

void PsociNOpopulations::constructMullikenPopulationDensAO(int nbf, int iroot, vector<double> & densAO, PsociProcessMOcoefs & coefs  )
{
  if (GA::nodeid() != 0 ) return;

  cout <<"Mulliken Population analysis: Atomic resolution" << endl;
  
  const int ROOT_READER=0;
  Fint l_unit = 4;
  
  Fdouble charge;
  Fdouble totalcharge = 0.0;
  
  /* Build an array of symmtry blocks: used by all roots.
   */
   int nsym = coefs.fetchNsym();
   vector<int> nstart(nsym,0);
   vector<int> nsize(nsym,0);
   
   int noffset = 0; 
   int nnbft = 0; // total length of a sym-blocked triangulat matrix
   for(int jsym=0; jsym< nsym; ++jsym) {
     nstart[jsym] = noffset;
     nsize[jsym] = coefs.fetch_nmopsy(jsym);
     noffset += coefs.fetch_nmopsy(jsym) ; 
     nnbft += NNDAF( nsize[jsym] );
   }
   
   string aoints_filename="aoints";
   PsociIntegrals aos( ROOT_READER, l_unit, aoints_filename );
   aos.fetchAndProcessOverlapIntegrals();
   
   /* Build a label to basis number array
      These are already in proper order so we can use this to carve up population space
   */
   vector<string> labels; 
   aos.fetchBasisLabels( labels );
   
   vector<Fdouble> overlap( nnbft, 0.0 );
   
   int idens_shift=0;
   for(int isym=0; isym< nsym; ++isym ) {
     
     int numbasis = nsize[isym];
     vector<vector<Fdouble> >overlap( numbasis, vector<Fdouble>( numbasis, 0.0 ) );
   
     int index = 0;
     int i1=0;
     int j1;
     for(int i=nstart[isym]; i< nstart[isym]+numbasis; ++i ) {
       j1=0;
       for(int j=nstart[isym]; j< nstart[isym]+numbasis; ++j ) {
	 overlap[i1][j1] = aos.fetch1eIntegralValue(i,j);
	 ++index;//TODO remove me 
	 ++j1;
       }
       ++i1;
     }

     // Read densAO in a full matrix of size numbasis*numbasis
     
     vector<vector<double > > density( numbasis, vector<Fdouble>(numbasis, 0.0) );
     vector<vector<double > > population( numbasis, vector<Fdouble>(numbasis, 0.0) );
     
     index = 0;
     for( int mu=0; mu< numbasis; ++mu ) { // Put a threshhold on the eval
       for(int nu=0; nu< numbasis; ++nu ) {
	 density[nu][mu] = densAO[ twoDindex(mu,nu) + idens_shift ]; // Local coef in square array (0,numbasis)
	 ++index; // TODO remove me
       }
     }

#ifdef DETAILEDCHECK
     for( int mu=0; mu< numbasis; ++mu ) { // Put a threshhold on the eval
       for(int nu=0; nu< numbasis; ++nu ) {
	 cout << "squ density densAO" << mu<<" "<<nu<<" "<< density[mu][nu] << endl;
       }
     }
#endif
     idens_shift += NNDAF( numbasis);
     
     charge = 0.0;
     for(int mu=0; mu< numbasis; ++mu ) {
       for(int nu=0; nu< numbasis; ++nu ) {
	 population[nu][mu] += density[nu][mu] * overlap[nu][mu];
	 // temp = population[nu][mu];
	 //++ov_index;
       }
     }
     
     for(int mu=0; mu< numbasis; ++mu ) {
       for(int nu=0; nu< numbasis; ++nu ) {
	 charge += population[nu][mu];
       }
     }
     
     totalcharge += charge;

     //cout << "densAO charge densAO " << charge << endl;

   cout << endl << "Current Symmetry block is " << isym << endl << endl;
     
   double at_temp = atomicPopulations(nstart[isym], nstart[isym]+numbasis, labels, population );
   double sb_temp = subshellPopulations(nstart[isym], nstart[isym]+numbasis, labels, population );
   double bd_temp = atomicOverlapPopulations(nstart[isym], nstart[isym]+numbasis, labels, population );

   //int num = clusterBasisLabels( condense,  labels , outlabels);

   } // isym

   cout << "Total electron charge for iroot = " << iroot << " is " << totalcharge << " nelec is " << l_nelec << endl;

  aos.CloseFile();

}

/* Alternative population approach. We know the NOs are correct. 
   coefs are symmetry blocked but square
   
   This is NOT intended for natural orbitals as they have already been derived from the density matrix

   NEED this method to compute properties on a per-orbital basis. 

*/

void PsociNOpopulations::constructMullikenPopulationNOs(int nbf, vector<double> & pk_evals, vector<double> & pk_coefs, PsociProcessMOcoefs & coefs  )
{
  if (GA::nodeid() != 0 ) return;
  cout <<"Mulliken Population analysis: Orbital resolution: Only orbitals with occupation >= " << setprecision(6) << MIN_EVAL_BOND_POPULATION <<" are resolved " << endl;
  
  const int ROOT_READER=0;
  Fint l_unit = 2;
  
  /* Build an array of symmtry blocks: used by all roots.
   */
  int noffset = 0;
  int nsym = coefs.fetchNsym();
  
  vector<int> nstart(nsym,0);
  vector<int> nsize(nsym,0);
  
  int nnbft = 0; // total length of a sym-blocked triangular matrix
  for(int jsym=0; jsym< nsym; ++jsym) {
    nstart[jsym] = noffset;
    nsize[jsym] = coefs.fetch_nmopsy(jsym);
    noffset += coefs.fetch_nmopsy(jsym) ; 
    nnbft += NNDAF( nsize[jsym] );
  }

   string aoints_filename="aoints";
   PsociIntegrals aos( ROOT_READER, l_unit, aoints_filename );
   aos.fetchAndProcessOverlapIntegrals();

   /* Build a label to basis number array
      These are already in proper order so we can use this to carve up population space
   */
   vector<string> labels;
   aos.fetchBasisLabels( labels );

   vector<Fdouble> overlap( nnbft, 0.0 );

/* useful in fetching MOcoefs */
   int icharge_shift = 0;
   int ev_index=0;
   int evalue_shift = 0;

   for(int isym=0; isym< nsym; ++isym ) {

// Fetch ints
     int numbasis = nsize[isym];
     vector<vector<Fdouble> >overlap( numbasis, vector<Fdouble>( numbasis, 0.0 ) );
     int index = 0;
     int i1=0;
     int j1;
     for(int i=nstart[isym]; i< nstart[isym]+numbasis; ++i ) {
       j1=0;
       for(int j=nstart[isym]; j< nstart[isym]+numbasis; ++j ) {
         overlap[i1][j1] = aos.fetch1eIntegralValue(i,j);
         ++index;
         ++j1;
       }
       ++i1;
     }
    
    vector<Fdouble> vec_coefs( numbasis, 0.0); // do these need to be Fdoubles?
    Fdouble allevals = 0.0;
    
// Build a coulson-like matrix for the specific orbital. Then call subshellPopulations for a resolution of the charge into shells

    index = 0;
    ev_index=0;

    int numNOs = numbasis;

    for( int iorb=0; iorb < numNOs; ++iorb ) { // Put a threshhold on the eval
      allevals = pk_evals[ev_index + evalue_shift ];
      if ( allevals >= MIN_EVAL_BOND_POPULATION ) { 

      cout <<"Process orbital " << iorb << " occ is " << allevals << endl;

      for(int mu=0; mu< numbasis; ++mu ) {
	vec_coefs[mu] = pk_coefs[ index + icharge_shift ]; // Local coef in square array (0,numbasis)
	++index;
      }

/* build a poplation matrix for this orbital : still need all basis fnts represented */
     vector<vector<double > > population( numbasis, vector<Fdouble>(numbasis, 0.0) );

     for(int mu=0; mu< numbasis; ++mu ) {
       for(int nu=0; nu< numbasis; ++nu ) {
         population[nu][mu] += vec_coefs[mu] * vec_coefs[nu] * allevals * overlap[nu][mu];
         //temp = population[nu][mu];
       }
     }

/* double check this should equal allevals */
     double charge = 0.0;
     for(int mu=0; mu< numbasis; ++mu ) { 
       for(int nu=0; nu< numbasis; ++nu ) {
         charge += population[nu][mu];
       }
     }
    // if ( charge != allevals ) GA::error(" Bond population: charge != allevals", charge );
    cout << "charge " << charge <<" allevals " << allevals << endl;

/* resolve into basis ftns */
    double sb_temp = subshellPopulations(nstart[isym], nstart[isym]+numbasis, labels, population );

/* resolve into Atomic centers ftns */
//    double sb_temp = subshellPopulations(nstart[isym], nstart[isym]+numbasis, labels, population );

      } // eval condition
      ++ev_index;
    } // Next orbital
    
    icharge_shift += numbasis*numbasis;
    evalue_shift += numbasis;
  } // isym
  
    aos.CloseFile();

}// end

/* Method to check orthonormality of the orbitals in "coefs" over the AO overlap metrix
   the coef are square but symmetry packed.
   
   If the orbs are NOT normalized this could be a symptom of combining the wriong aoints, moints
   and mocoef input files to the problem. Maybe this should be flagged and indicated to the user
*/

void PsociNOpopulations::checkOrthonormality(int nbf, PsociProcessMOcoefs & coefs  ) 
{
  if ( GA::nodeid() != 0 ) return;
  
  Fint l_unit = 3;
  
  /* Build an array of symmtry blocks: used by all roots.
   */
  int noffset = 0;
  int nsym = coefs.fetchNsym();
  int realorb_no=0;
  
  vector<int> nstart(nsym,0); 
  vector<int> nsize(nsym,0);
  
  int nnbft = 0; // total length of a sym-blocked triangular matrix
  for(int jsym=0; jsym< nsym; ++jsym) {
    nstart[jsym] = noffset;
    nsize[jsym] = coefs.fetch_nmopsy(jsym);
    noffset += coefs.fetch_nmopsy(jsym) ; 
    nnbft += NNDAF( nsize[jsym] );
  }
  
  // Get the AOs: put them into a convenient and square format
  
  string aoints_filename="aoints";
  PsociIntegrals aos( 0, l_unit, aoints_filename );
  aos.fetchAndProcessOverlapIntegrals();
  
  int ishift = 0;
  for(int isym=0; isym< nsym; ++isym ) { // check all orbs  
    
    int numbasis = nsize[isym];
    vector<Fdouble> overlap( numbasis * numbasis, 0.0 );
    
    int index = 0;
    for(int i=nstart[isym]; i< nstart[isym]+numbasis; ++i ) {
      for(int j=nstart[isym]; j< nstart[isym]+numbasis; ++j ) {
	overlap[ index ] = aos.fetch1eIntegralValue(i,j);
	++index;
      }
    }
    /*
      cout << "square overlap for sym block fortran index " << endl;
      index = 0;
      for(int i=0; i< numbasis; ++i ){
      for(int j=0; j< numbasis; ++j ) {
      cout << i<<" "<<j<<" "<<overlap[ index ] << endl;
      ++index;
      }
      }
    */
    
    // Get the coef handle, coefs and put them into a convenient handle one NO at a time.
    
    vector<Fdouble> * mocoef = coefs.fetchCoefsHandle();
    
    index = 0; // start with first orb
    for(int iorb = 0; iorb < nsize[isym]; ++iorb ) { // check norm of al orbitals
      
      // int ishift = nstart[isym];
      vector<Fdouble> vector( numbasis, 0.0 );
      
      for (int i=0; i< numbasis; ++i ) {
	vector[i] = coefs.fetchMocoefs( index + ishift );
	++index;
      }
      
      //Compute the normalization C**t*S*C - skip orthogonality check for now
      
      int o_index = 0;
      Fdouble norm = 0.0;
      for (int i=0; i< numbasis; ++i ) {
	Fdouble Cmu = vector[i];
	for (int j=0; j< numbasis; ++j ) { 
	  Fdouble Cnu = vector[j];   
	  Fdouble ovrlp = overlap[o_index];
	  norm += Cmu * Cnu * ovrlp;
	  ++o_index;
	}
      }
      //cout <<" computed MO norm for MO = " << realorb_no << " is " << norm << endl;
      ++realorb_no;
    } // iorb
    ishift += numbasis*numbasis; // For next go'round. next sym of coefs
  } // isym
  cout << " Finished with MO normalization check" << endl;
  aos.CloseFile();
  
}

/* We need a method to take the list of basis functions and performing some combining
   
Given a basis set labels list as follows

58Ru_4f
59Ru_4f
60Ru_4f
61Ru_4f


We can do various levels of condensation via the enum cond_level
choices include:

NONE ( does nothing ) 

SUBSHELL ( strips off preceding nums ) This permits the calling app to  identify common population labels
for summing. For example all Ru_1s will be summed

BASIS ( replaces preceding nums  ) 

ATOMIC ( strips off preceding nums AND trailing subshell indicators ) E.g., 60Ru_4f -> Ru

Currently this code is very ARGOS specific
As a result EACH input label is a ful 8Bytes long with padding blanks possible preceding the numerical values

*/

int PsociNOpopulations::clusterBasisLabels( CONDENSE condense,  vector<string> & inlabels , vector<string> & outlabels)
{
  if ( GA::nodeid() != 0 ) return( 0 );
  
  int startsize = inlabels.size();
  outlabels.clear();
  
  vector<string>::iterator it;
  if ( condense == SUBSHELL ) {
    int istart;
    for(it = inlabels.begin(); it != inlabels.end(); ++it ) {
      string word = (*it);
      int wordsize = word.size();
      istart = 0;
      for(int i=0; i< wordsize; ++i ) {
	if (isdigit(word[i]) || !isalnum(word[i]) ) { // include numbers AND blanks
	  ++istart;
	} else {
	  break;
	}
      }
      string subword = word.substr( istart, word.size()); // stops at the end
      outlabels.push_back( subword );
    }
    if ( outlabels.size() != inlabels.size() ) GA::error(" outlabels.size != inlabels.size ", outlabels.size() );
    
  } else if ( condense == ATOMIC ) {
    
    int istart;
    for(it = inlabels.begin(); it != inlabels.end(); ++it ) {
      string word = (*it);
      int wordsize = word.size();
      
      istart = 0;
      for(int i=0; i< wordsize; ++i ) {
	if (isdigit(word[i]) || !isalnum(word[i]) ) {
	  ++istart;
	} else {
	  break;
	}
      }
      int length = word.find_first_of('_') - istart;
      string subword = word.substr( istart, length); // stops at the end
      outlabels.push_back( subword );     
    }
    
    if ( outlabels.size() != inlabels.size() ) GA::error("ATOMIC outlabels.size != inlabels.size ", outlabels.size() );
  } else { 
    
    // NONE is implied
    outlabels = inlabels;
  }
  
  return (  outlabels.size() );
}

/* Assemble the ATOMIC charge distributions based on a mulliken population 

   nstart/nend refer to the basis fnt offsets that are required here
   population is [numbasis][numbasis[ is size
*/

double PsociNOpopulations::atomicPopulations(int nstart, int nend, vector<string> & labels, vector<vector<double> > population )
{
  if ( GA::nodeid() != 0 ) return( 0.0 );

#ifdef DETAILEDCHECK
  cout << "Compute atomic populations " << nstart<<" "<<nend<< endl;
#endif
  
  CONDENSE condense = ATOMIC;
  vector<string> outlabels; //should be the same length as labels
  
  int num = clusterBasisLabels( condense,  labels , outlabels ); // grab list of ATOMIC oriented labels
  
#ifdef DETAILEDCHECK
  vector<string>::iterator it;
  for(it=outlabels.begin(); it != outlabels.end(); ++it ) {
    cout << " "<<(*it)<<endl;
  }
#endif
  
  /* Loop over all basis or all labels ?
     How to best carve up the domain  ?
     Note labels can be MIXED and currently we sum all terms for all common atoms NO distinction between C1 and C2
     But if the labels are mixed so are the basis functions. So what is the best way to access them? 
  */
  
  /* 
     Build an array of indexes that are used to update a particular bfn_label
     Then we can simply loop through the set of indexes
  */
  
  // First find the unique labels and associated with a set of indices
  // also only grab indeces within the range nstart, nend-1
  // and offset for the fact that population is a relative indexing...
  
  int nsize = outlabels.size(); // bfns are in-order !
  
  // TODO This is fully generated for each symmetry block ( fix this )
  
  multimap<string, int> labelIndex;
  for(int i=0; i< nsize; ++i ) {
    string lab = outlabels[i];
    if ( i >= nstart && i < nend ) labelIndex.insert( pair<string,int>(lab, i-nstart) ); // note the offset
  }
  
  // Now build a loop based on the keys. followed by processing the index loops
  
  multimap<string, int>::iterator it2, ot2;
  multimap<string, int>::iterator inneri, innerj;
  
  double Fsum;
  double Ftotal=0.0;
  
  it2 = labelIndex.begin();
  
  cout << "Atomic label: " << "Mulliken population analysis (partial) " << endl;
  cout << "-----------------------------------------------------------" << endl;

  while ( it2 != labelIndex.end() ) // Select unique atoms
    {
      Fsum = 0.0;
      string bfn_label = (*it2).first;
      ot2 = labelIndex.upper_bound(bfn_label);
      for( inneri=it2; inneri != ot2; ++inneri ) // choose an atom center
	{
          int i_index = (*inneri).second;
          string i_label = (*inneri).first;
          for( innerj=labelIndex.begin(); innerj != labelIndex.end(); ++innerj )
	    {
	      int j_index = (*innerj).second;
              string j_label = (*innerj).first;
              Fsum += population[j_index][i_index];
	    }
	}
      Ftotal += Fsum;
      it2 = ot2;
      cout << setw(10) << bfn_label<< setw(10)<< setprecision(3)<< Fsum << endl;
      //cout << "Fsum label is " << bfn_label<<" sum is " << Fsum << endl;
    }
  cout << setw(10) << "sum" << setw(10)<< setprecision(3)<< Ftotal << endl;
  
  return( Fsum );
}
/* Assemble the ATOMIC bond distributions based on a mulliken population 
   
nstart/nend refer to the basis fnt offsets that are required here
population is [numbasis][numbasis[ is size
*/
double PsociNOpopulations::atomicOverlapPopulations(int nstart, int nend, vector<string> & labels, vector<vector<double> > population )
{
  if ( GA::nodeid() != 0 ) return( 0.0 );
  
#ifdef DETAILEDCHECK
  cout << "Compute atomic overlap populations " << nstart<<" "<<nend<< endl;
#endif
  
    CONDENSE condense = ATOMIC;
    vector<string> outlabels; //should be the same length as labels
    int num = clusterBasisLabels( condense,  labels , outlabels ); // grab list of ATOMIC oriented labels
    
#ifdef DETAILEDCHECK
    vector<string>::iterator it;
    for(it=outlabels.begin(); it != outlabels.end(); ++it ) {
      cout << " "<<(*it)<<endl;
    }
#endif
    
    /* 
       Loop over all unique PAIRS of labels ( atomic-centric are the only ones that make sense here )
       Then for each compute the overlap - Q? do we still do mulliken overlap partitioning withing this group       
       -- I think NO. But note the values of Left: right and right:left are identical and SHOULD NOT be summed.
    */
    
    int nsize = outlabels.size(); // bfns are in-order !
    
    // This is fully generated for each symmetry block ( fix this )
    
    multimap<string, int> labelIndex;
    for(int i=0; i< nsize; ++i ) {
      string lab = outlabels[i];
      if ( i >= nstart && i < nend ) labelIndex.insert( pair<string,int>(lab, i-nstart) ); // note the offset
    }
    
    // Now build a loop based on the keys. followed by processing the index loops
    
    multimap<string, int>::iterator atL,atR;
    multimap<string, int>::iterator atm_left_out, atm_right_out;
    multimap<string, int>::iterator index_left, index_right;
    
    double Fsum = 0.0;
    double Ftotal=0.0;
    
    cout << "overlap pair: " << "Partial Mulliken overlap population " << endl;
    cout << "--------------------------------------------------------" << endl;
    
    string atm_left, atm_right;
    
    atL = labelIndex.begin();
    Ftotal = 0.0;

/* We do not want to printout both L-R and R-L pairs as this can be misleading
   So only do the triangle of pairs and multiply by 2
*/
    double factor=0.0;
    while ( atL != labelIndex.end() ) // select all atoms
      {
	atm_left = (*atL).first;
	atR = labelIndex.begin();
	atm_left_out = labelIndex.upper_bound(atm_left);

        atR = atL; // triangle

	while ( atR != labelIndex.end() )
	  {
            atm_right = (*atR).first;
            atm_right_out = labelIndex.upper_bound(atm_right);
	    
            Fsum = 0.0;
           // if ( atm_left != atm_right ) {
	      
	      for( index_left=atL; index_left != atm_left_out; ++index_left )
		{
		  for( index_right=atR; index_right != atm_right_out; ++index_right )
		    {
		      ( atm_left == atm_right ) ? factor = 1.0: factor = 2.0;
                      Fsum += factor * population[ (*index_left).second][(*index_right).second];
		    }
		}
            cout << setw(4) << atm_left<< " : "<< setw(4) << atm_right<< setw(10)<< setprecision(3)<< Fsum << endl;
            Ftotal += Fsum;
	    //}
	    atR = atm_right_out;
	  }
	atL = atm_left_out;
      }
    // cout << setw(4) << "sum"<< " : "<< setw(4) << "all"<< setw(10)<< setprecision(3)<< Ftotal << endl;
    
    return( Fsum );
}

/* Assemble the SUBSHELL charge distributions based on a mulliken population 
   nstart/nend refer to the basis fnt offsets that are required here
   population is [numbasis][numbasis[ is size
*/
double PsociNOpopulations::subshellPopulations(int nstart, int nend, vector<string> & labels, vector<vector<double> > population )
{
  if ( GA::nodeid() != 0 ) return( 0.0 );

#ifdef DETAILEDCHECK
  cout << "Compute subshell populations " << nstart<<" "<<nend<< endl;
#endif
  
  CONDENSE condense = SUBSHELL;
  vector<string> outlabels; //should be the same length as labels
  
  int num = clusterBasisLabels( condense,  labels , outlabels ); // grab list of SUBSHELL oriented labels
  
#ifdef DETAILEDCHECK
  vector<string>::iterator it;
  for(it=outlabels.begin(); it != outlabels.end(); ++it ) {
    cout << " "<<(*it)<<endl;
  }
#endif
  
  /* Loop over all basis or all labels ?
     How to best carve up the domain  ?
     Note labels can be MIXED and currently we sum all terms for all common atoms NO distinction between C1 and C2
     But if the labels are mixed so are the basis functions. So what is the best way to access them? 
  */
  
  // First find the unique labels and associated with a set of indices
  // also only grab indeces within the range nstart, nend-1
  // and offset for the fact that population is a relative indexing...
  
  int nsize = outlabels.size(); // bfns are in-order !
  
  // This is fully generated for each symmetry block ( fix this )
  
  multimap<string, int> labelIndex;
  for(int i=0; i< nsize; ++i ) {
    string lab = outlabels[i];
    if ( i >= nstart && i < nend ) labelIndex.insert( pair<string,int>(lab, i-nstart) ); // note the offset
  }
  
  // Now build a loop based on the keys. followed by processing the index loops
  
  multimap<string, int>::iterator it2, ot2;
  multimap<string, int>::iterator inneri, innerj;
  
  double Fsum;
  double Ftotal=0.0;
  
  it2 = labelIndex.begin();
  
  cout << "subshell label: " << "Mulliken population analysis (partial)" << endl;
  cout << "------------------------------------------------------------" << endl;
  while ( it2 != labelIndex.end() )
    {
      Fsum = 0.0;
      string bfn_label = (*it2).first;
      ot2 = labelIndex.upper_bound(bfn_label);
      for( inneri=it2; inneri != ot2; ++inneri )
	{
          int i_index = (*inneri).second;
          string i_label = (*inneri).first;

          for( innerj=labelIndex.begin(); innerj != labelIndex.end(); ++innerj )
            {
              int j_index = (*innerj).second;
              string j_label = (*innerj).first;
              Fsum += population[j_index][i_index];
            }
	}
      Ftotal += Fsum;
      it2 = ot2;
      cout << setw(10) << bfn_label<< setw(10)<< setprecision(3)<< Fsum << endl;
    }
  cout << setw(10) << "sum" << setw(10)<< setprecision(3)<< Ftotal << endl;  
  return( Fsum );
}
