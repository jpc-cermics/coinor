/* Nsp
 * Copyright (C) 2014-2014 Jean-Philippe Chancelier Enpc/Cermics
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/*  use COIN-OR File Reading Routines */

#include <limits>
#include "CoinMpsIO.hpp"
#include "CoinLpIO.hpp"
#include "CoinModel.hpp"
#include "CoinMessageHandler.hpp"

extern "C" {
#include <nsp/interf.h>
#include <nsp/matrix.h>
#include <nsp/smatrix.h>
#include <nsp/hash.h>
#include <nsp/spcolmatrix.h>
#include <nsp/cells.h>
#include <nsp/sciio.h>
#include <mex/mex.h>
#include <nsp/system.h> /* FSIZE */
#include <nsp/nsptcl.h> /* get_extension */

  int  int_coinmp_readlp(Stack stack, int rhs, int opt, int lhs);
  static NspSpColMatrix *spcolmatrix_from_triplet(const char *name,const CoinBigIndex *jc,
						  const int *ir,const double *pr,int nelem,int m, int n);
  double nsp_coin_dbl_max(void);
}

using namespace std;

double nsp_coin_dbl_max(void)
{
  return COIN_DBL_MAX;
}

/* provide our own print handlers */

class NReadDerivedHandler : public CoinMessageHandler {
public:
  virtual int print() ;
};

int NReadDerivedHandler::print()
{
  Sciprintf(messageBuffer());
  Sciprintf("\n");
  return 0;
}

class NReadDummyHandler : public CoinMessageHandler {
public:
  virtual int print() ;
};

int NReadDummyHandler::print()
{
  return 0;
}

/* utility */

static NspSpColMatrix *spcolmatrix_from_triplet(const char *name,const CoinBigIndex *jc,const int *ir,const double *pr,int nelem,int m, int n)
  {
    NspSpColMatrix *A;
    int i;
    if ((A =nsp_spcolmatrix_create(name,'r',m,n) ) == NULLSPCOL) return NULL;
    if ( nsp_spcol_alloc_col_triplet(A,nelem)== FAIL) return NULL;
    for(i = 0; i <= n; i++)   A->triplet.Jc[i] = (mwIndex)jc[i]; 
    for(i = 0; i < nelem; i++) A->triplet.Ir[i] = ir[i];
    memcpy(A->triplet.Pr,pr, nelem*sizeof(double));
    if ( nsp_spcol_update_from_triplet(A) == FAIL) return NULL;
    return A;
  }

/* Reader using coin functions can read a quadratic pb 
 */

int int_coinmp_readlp(Stack stack, int rhs, int opt, int lhs)
{
  const CoinPackedMatrix *pm;
  const CoinBigIndex *jc;
  const int *ir;
  int err = -1, *ijc, *iir;
  double *elements;
  NReadDerivedHandler *mexprinter = NULL;  
  NReadDummyHandler *dumprinter = NULL;  
  int nbr_sos_sets = 0;
  CoinSet **sets;

  NspHash *D;
  NspMatrix *Obj_offset;
  NspCells *Sos = NULLCELLS;
  NspSMatrix *Sense = NULLSMAT;
  const char *row_sense;
  NspSMatrix *var_type = NULLSMAT;
  NspMatrix *C=NULLMAT, *Lb=NULLMAT, *Ub=NULLMAT;
  NspMatrix *Rl=NULLMAT, *Ru=NULLMAT;
  NspSpColMatrix *A=NULLSPCOLMAT, *H=NULLSPCOLMAT;

  const char *type_str=NULL, *type_str_def = "mps";
  size_t i;
  int_types T[] = {string_c, new_opts, t_end} ;
  nsp_option opts[] ={{"type",string_c,NULLOBJ,-1},
		      {"printlevel",s_int,NULLOBJ,-1},
 		      { NULL,t_end,NULLOBJ,-1}};
  int printLevel = 0;
  char Fname_expanded[FSIZE+1];
  char *Fname;
  const char *extension;

  if ( GetArgs(stack,rhs,opt,T,&Fname, &opts, &type_str, &printLevel) == FAIL) 
    return RET_BUG;
  nsp_expand_file_with_exec_dir(&stack,Fname,Fname_expanded);
    
  extension = nsp_get_extension(Fname);
  if ( type_str == NULL && extension != NULL ) 
    {
      type_str = extension +1;
      /* Sciprintf("using %s for type\n",type_str); */
    }

  if ( type_str == NULL) type_str = type_str_def;

  if( !strcmp(type_str,"mps") || !strcmp(type_str,"qps") || !strcmp(type_str,"mod") || !strcmp(type_str,"gms")) 
    {
      CoinMpsIO m;
      size_t ncol, nrow;
      double infinity = numeric_limits<double>::infinity();

      if(printLevel) 
	{
	  mexprinter = new NReadDerivedHandler();
	  mexprinter->setLogLevel(printLevel);
	  m.passInMessageHandler(mexprinter);
	}
      else 
	{
	  dumprinter = new NReadDummyHandler();
	  dumprinter->setLogLevel(0);
	  m.passInMessageHandler(dumprinter);
	}
      /* Setup Options */
      m.setInfinity(infinity);   
      
      /* read pb */
      try 
	{
	  if(!strcmp(type_str,"mps") || !strcmp(type_str,"qps"))
	    {
	      err = m.readMps(Fname_expanded,"mps",nbr_sos_sets,sets);
	    }
	  else if(!strcmp(type_str,"mod"))
	    {
#ifdef COIN_HAS_GLPK
	      char *data = NULL;
	      err = m.readGMPL(Fname_expanded,data,false);
#else 
	      Scierror("Error: failed to read file %s, coin library was not compiled with glpk \n",Fname);
	      return RET_BUG;
#endif 
	    }
	  else if(!strcmp(type_str,"gms"))
	    {
	      err = m.readGms(Fname_expanded,"gms",false);
	    }
	  if ( err ) 
	    {
	      Scierror("Error: failed to read file %s, Error Code: %d\n",Fname,err);
	      return RET_BUG;
	    }
	}
      catch (CoinError e) {
	Scierror("Error: while reading Lp file \n"); // e.message());
	return RET_BUG;
      }
      catch (...) {
	Scierror("Error: while reading Lp file \n");
	return RET_BUG;
      }
      
      if ( ( D = nsp_hash_create(NVOID, 20) ) == NULLHASH) return RET_BUG;
      
      ncol = m.getNumCols();
      nrow = m.getNumRows();
      
      if ( ncol <= 0 ) 
	{
	  Scierror("Error: number of columns is 0\n");
	  return RET_BUG;
	}
      
      /* Get Objective Coefficients */
      if ((C = nsp_matrix_create("c",'r',ncol,1) ) == NULLMAT ) return FAIL;
      memcpy(C->R,m.getObjCoefficients(),ncol*sizeof(double));

      if (nsp_hash_enter(D,NSP_OBJECT(C))== FAIL) return RET_BUG;
      
      /* Get Sparse Matrix A a la matlab */
      pm = m.getMatrixByCol();
      size_t nelem = pm->getNumElements();
      jc = pm->getVectorStarts();
      ir = pm->getIndices();
      A = spcolmatrix_from_triplet("A",jc,ir,pm->getElements(),nelem,nrow,ncol);
      if ( A == NULL) return FAIL;

      if (nsp_hash_enter(D,NSP_OBJECT(A))== FAIL) return RET_BUG;
      
      /* we return the RowLower, the RowUpper */
      /* 
       * (lower,upper) -> (sense, right,range) 
       * (finite, finite) -> ('E' | 'R', upper, upper-lower)
       * (finite, infinite)-> ('G' , lower, .)
       * (-infinity,finite)-> ('L',  upper, .)
       * (-infinity,+infinity)-> ('N',0,0);
       *
       * (sense, right,range) -> (lower,upper) 
       * 'E'                  -> (right,right)
       * 'L'                  -> (-infinity,right)
       * 'G'                  -> (right, infinity)
       * 'R'                  -> (right-range, right)
       * 'N':                 -> (-infinity,+infinity)
      */
      
      if ((Sense = nsp_smatrix_create("Sense", nrow, 1, "L", 1))  == NULLSMAT ) return FAIL;
      row_sense = m.getRowSense();
      for ( i = 0 ; i < nrow ; i++)
	{
	  Sense->S[i][0]=row_sense[i];
	}
      
      if (nsp_hash_enter(D,NSP_OBJECT(Sense))== FAIL) return RET_BUG;


      /* Row lower bounds */
      if ((Rl = nsp_matrix_create("rl",'r',nrow,1) ) == NULLMAT ) return FAIL;
      if( m.getRowLower() != NULL)
	{
	  memcpy(Rl->R,m.getRowLower(),nrow*sizeof(double));
	}
      else
	{
	  for (i=0 ; i < nrow ; i++) Rl->R[i] = - infinity;
	}
      
      if (nsp_hash_enter(D,NSP_OBJECT(Rl))== FAIL) return RET_BUG;

      /* Row upper bounds */
      if ((Ru = nsp_matrix_create("ru",'r',nrow,1) ) == NULLMAT ) return FAIL;
      if ( m.getRowUpper() != NULL)
	{
	  memcpy(Ru->R,m.getRowUpper(),nrow*sizeof(double));
	}
      else
	{
	  for (i=0 ; i < nrow ; i++) Ru->R[i] = infinity;
	}
      
      if (nsp_hash_enter(D,NSP_OBJECT(Ru))== FAIL) return RET_BUG;
      
      /* lb : column lower bounds */
      if ((Lb = nsp_matrix_create("cl",'r',ncol,1) ) == NULLMAT ) return FAIL;
      if(m.getColLower() != NULL) memcpy(Lb->R,m.getColLower(),ncol*sizeof(double));
      /* ub : column upper bounds */
      if ((Ub = nsp_matrix_create("cu",'r',ncol,1) ) == NULLMAT ) return FAIL;
      if(m.getColUpper() != NULL) memcpy(Ub->R,m.getColUpper(),ncol*sizeof(double)); 

      /* proper infinity for mod */
      if(!strcmp(type_str,"mod")) 
	{
	  for(i=0;i<ncol;i++)
	    {
	      if(Lb->R[i] < -1e300) Lb->R[i] = - infinity;
	      if(Ub->R[i] > 1e300)  Ub->R[i] = infinity;
	    }
	  for(i=0;i<nrow;i++)
	    {
	      if(Rl->R[i] < -1e300) Rl->R[i] = - infinity;
	      if(Ru->R[i] > 1e300)  Ru->R[i] = infinity;
	    }
	}

      if (nsp_hash_enter(D,NSP_OBJECT(Lb))== FAIL) return RET_BUG;
      if (nsp_hash_enter(D,NSP_OBJECT(Ub))== FAIL) return RET_BUG;
      
      /* integer/continuous variables  */
      if ((var_type = nsp_smatrix_create("ivar", ncol, 1, "C", 1))  == NULLSMAT ) return FAIL;
      if(m.integerColumns() != NULL)
	{
	  int i;
	  for (i=0; i < var_type->mn ; i++) 
	    var_type->S[i][0]= (m.integerColumns()[i] == 0) ? 'C' : 'I';
	}

      if (nsp_hash_enter(D,NSP_OBJECT(var_type))== FAIL) return RET_BUG;

      /* quadratic objective if present */
      err = m.readQuadraticMps(NULL, ijc, iir, elements, 0);
      if(!err) 
	{ 
	  /* read a quadratic objective */
	  nelem = (size_t)ijc[ncol];
	  H = spcolmatrix_from_triplet("H",ijc,iir,elements,nelem, ncol,ncol);
        }
      else //No Quad Section
	{
	  if ((H =nsp_spcolmatrix_create("H",'r',ncol,ncol) ) == NULLSPCOL) return FAIL;
	}

      if (nsp_hash_enter(D,NSP_OBJECT(H))== FAIL) return RET_BUG;
      
      
      if ((Sos = nsp_cells_create("Sos",nbr_sos_sets,1)) == NULLCELLS ) return FAIL;
      for ( i=0 ; i <  (unsigned int) Sos->mn ; i++)
	{
	  int j;
	  int sos_size;
	  const double *sos_weights;
	  const int *sos_index;
	  NspMatrix *Sos_index, *Sos_type, *Sos_weights;
	  NspCells *Elt; 
	  if ((Elt = nsp_cells_create("elt",3,1)) == NULLCELLS ) return FAIL;
	  if ((Sos_type = nsp_matrix_create("type",'r',1,1)) == NULLMAT ) return FAIL; 
	  Sos_type->R[0] = (double)sets[i]->setType();
	  Elt->objs[0] = NSP_OBJECT(Sos_type);
	  /* Sos_index matrix */
	  sos_size = sets[i]->numberEntries();
	  sos_index  = sets[i]->which();
	  if ((Sos_index = nsp_matrix_create("index",'r',sos_size,1)) == NULLMAT ) return FAIL;
	  for ( j=0 ; j < sos_size;j++) Sos_index->R[j]= (double)sos_index[j] + 1;
	  Elt->objs[1] = NSP_OBJECT(Sos_index);
	  /* weights */
	  if ((Sos_weights = nsp_matrix_create("weights",'r',sos_size,1)) == NULLMAT ) return FAIL;
	  sos_weights= sets[i]->weights();
	  memcpy(Sos_weights->R,sos_weights,sos_size *sizeof(double));
	  Elt->objs[1] = NSP_OBJECT(Sos_weights);
	  Sos->objs[i]= NSP_OBJECT(Elt);
	}
      
      if (nsp_hash_enter(D,NSP_OBJECT(Sos))== FAIL) return RET_BUG;

      if ((Obj_offset = nsp_matrix_create("offset",'r',1,1) ) == NULLMAT ) return FAIL;
      Obj_offset->R[0]= - m.objectiveOffset();

      if (nsp_hash_enter(D,NSP_OBJECT(Obj_offset))== FAIL) return RET_BUG;

      MoveObj(stack,1,NSP_OBJECT(D));
      return Max(lhs,1);      
    }
  else if(!strcmp(type_str,"lp"))
    {
      CoinLpIO l;
      size_t ncol, nrow;
      double infinity = numeric_limits<double>::infinity();
      
      if(printLevel) 
	{
	  mexprinter = new NReadDerivedHandler();
	  mexprinter->setLogLevel(printLevel);
	  l.passInMessageHandler(mexprinter);
	}  
      else
	{
	  dumprinter = new NReadDummyHandler();
	  dumprinter->setLogLevel(0);
	  l.passInMessageHandler(dumprinter);
	}
      
      /* Setup Options */
      l.setInfinity(numeric_limits<double>::infinity());
      /* Read Problem */
      try 
	{
	  l.readLp(Fname_expanded);
	}
      catch (CoinError e) {
	Scierror("Error: while reading Lp file \n"); // e.message());
	return RET_BUG;
      }
      catch (...) {
	Scierror("Error: while reading Lp file \n");
	return RET_BUG;
      }


      l.readLp(Fname_expanded);

      if ( ( D = nsp_hash_create(NVOID, 20) ) == NULLHASH) return RET_BUG;
      
      ncol = l.getNumCols();
      nrow = l.getNumRows();

      if(ncol <= 0) 
	{
	  Scierror("Error: number of columns is 0\n");
	  return RET_BUG;
	}

      /* Get Objective Coefficients */
      if ((C = nsp_matrix_create("c",'r',ncol,1) ) == NULLMAT ) return FAIL;
      memcpy(C->R,l.getObjCoefficients(),ncol*sizeof(double));
      
      if (nsp_hash_enter(D,NSP_OBJECT(C))== FAIL) return RET_BUG;
      
      /* Get Sparse Matrix A a la matlab */
      pm = l.getMatrixByCol();
      size_t nelem = pm->getNumElements();
      jc = pm->getVectorStarts();
      ir = pm->getIndices();
      A = spcolmatrix_from_triplet("A",jc,ir,pm->getElements(),nelem,nrow,ncol);
      if ( A == NULL) return FAIL;

      if (nsp_hash_enter(D,NSP_OBJECT(A))== FAIL) return RET_BUG;
      
      /* Row lower bounds */
      if ((Rl = nsp_matrix_create("rl",'r',nrow,1) ) == NULLMAT ) return FAIL;
      if( l.getRowLower() != NULL)
	{
	  memcpy(Rl->R,l.getRowLower(),nrow*sizeof(double));
	}
      else
	{
	  for (i=0 ; i < nrow ; i++) Rl->R[i] = - infinity;
	}
      
      if (nsp_hash_enter(D,NSP_OBJECT(Rl))== FAIL) return RET_BUG;

      /* Row upper bounds */
      if ((Ru = nsp_matrix_create("ru",'r',nrow,1) ) == NULLMAT ) return FAIL;
      if ( l.getRowUpper() != NULL)
	{
	  memcpy(Ru->R,l.getRowUpper(),nrow*sizeof(double));
	}
      else
	{
	  for (i=0 ; i < nrow ; i++) Ru->R[i] = infinity;
	}
      
      if (nsp_hash_enter(D,NSP_OBJECT(Ru))== FAIL) return RET_BUG;
      
      /* lb : column lower bounds */
      if ((Lb = nsp_matrix_create("cl",'r',ncol,1) ) == NULLMAT ) return FAIL;
      if(l.getColLower() != NULL) memcpy(Lb->R,l.getColLower(),ncol*sizeof(double));
      /* ub : column upper bounds */
      if ((Ub = nsp_matrix_create("cu",'r',ncol,1) ) == NULLMAT ) return FAIL;
      if(l.getColUpper() != NULL) memcpy(Ub->R,l.getColUpper(),ncol*sizeof(double)); 

      if (nsp_hash_enter(D,NSP_OBJECT(Lb))== FAIL) return RET_BUG;
      if (nsp_hash_enter(D,NSP_OBJECT(Ub))== FAIL) return RET_BUG;

      /* integer/continuous variables  */
      if ((var_type = nsp_smatrix_create("ivar", ncol, 1, "C", 1))  == NULLSMAT ) return FAIL;
      if(l.integerColumns() != NULL)
	{
	  int i;
	  for (i=0; i < var_type->mn ; i++) 
	    var_type->S[i][0]= (l.integerColumns()[i] == 0) ? 'C' : 'I';
	}

      /* no quadratic objective  */
      if ((H =nsp_spcolmatrix_create("H",'r',ncol,ncol) ) == NULLSPCOL) return FAIL;
      
      if (nsp_hash_enter(D,NSP_OBJECT(H))== FAIL) return RET_BUG;

      if ((Sos = nsp_cells_create("Sos",0,1)) == NULLCELLS ) return FAIL;
      
      if (nsp_hash_enter(D,NSP_OBJECT(Sos))== FAIL) return RET_BUG;

      if ((Obj_offset = nsp_matrix_create("offset",'r',1,1) ) == NULLMAT ) return FAIL;
      Obj_offset->R[0]= - l.objectiveOffset();

      if (nsp_hash_enter(D,NSP_OBJECT(Obj_offset))== FAIL) return RET_BUG;

      MoveObj(stack,1,NSP_OBJECT(D));
      return Max(lhs,1);      
    }
  else
    {
      Scierror("Error: failed to read file %s with unrecognized type %s\n",Fname,type_str);
      return RET_BUG;
    }

  return RET_BUG;
}

double nsp_coimp_void1()
{
  return COIN_DBL_MIN + COIN_INT_MAX + COIN_INT_MAX_AS_DOUBLE;
}





