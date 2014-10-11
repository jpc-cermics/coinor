/* Nsp
 * Copyright (C) 2014-2014 J.-Ph. Chancelier Cermics/ENPC 
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
 *
 */

#include "nspclp.h"
#include "Coin_C_defines.h"
#include "CoinMessageHandler.hpp"
#include "ClpSimplex.hpp"
#include "ClpInterior.hpp"
#include "ClpCholeskyBase.hpp"

/* use Sciprintf to print messages */

class DerivedHandler :
  public CoinMessageHandler {
public:
  virtual int print() ;
};

int DerivedHandler::print()
{
  Sciprintf(messageBuffer());
  Sciprintf("\n");
  return 0;
}


int nsp_clp_solve(nsp_clp_params *options,int sense, int ncols, int nrows, int neq, 
		  NspIMatrix*Cmatbeg, NspIMatrix *Cmatind, NspMatrix *Cmatval, 
		  NspMatrix *lower, NspMatrix *upper, NspMatrix *Objective,
		  NspIMatrix*Qmatbeg, NspIMatrix *Qmatind, NspMatrix *Qmatval, 
		  NspMatrix *Rhs,  NspMatrix *Lhs,char *var_type[],  NspMatrix *X,NspMatrix *Lambda,
		  NspMatrix *RetCost,NspMatrix *Retcode,
		  const char *filename,int save_only)
{
  double *primal, *dual;

  ClpSimplex *modelByColumn = NULL;
  ClpInterior *modelPrimalDual = NULL;
  DerivedHandler * mexprinter = NULL;

  if ( Qmatbeg != NULL ) options->solverchoice=3;
  
  switch (options->solverchoice)
    {
    default:
      {				
	modelByColumn = new ClpSimplex();
	modelByColumn->loadProblem(ncols,nrows,(const int*)Cmatbeg->Iv,(const int*)Cmatind->Iv,
				   Cmatval->R,lower->R,upper->R,Objective->R,Lhs->R,Rhs->R);
	if ( var_type != NULL) 
	  {
	    int i;
	    for ( i=0; i < ncols ; i++)
	      {
		if ( var_type[i] != NULL &&  var_type[i][0]== 'I') 
		  modelByColumn->setInteger(i);
	      }
	  }
	modelByColumn->setOptimizationDirection((sense==0) ? 1: -1);
	if ( Qmatbeg != NULL ) 
	  {	
	    modelByColumn->loadQuadraticObjective(ncols,(const int*) Qmatbeg->Iv,(const int*) Qmatind->Iv,Qmatval->R);
	  }	
	break;
      }
    case 3:
      {						
	modelPrimalDual = new ClpInterior();
	modelPrimalDual->loadProblem(ncols,nrows,(const int*)Cmatbeg->Iv,(const int*)Cmatind->Iv,
				     Cmatval->R,lower->R,upper->R,Objective->R,Lhs->R,Rhs->R);
	if ( var_type != NULL) 
	  {
	    int i;
	    for ( i=0; i < ncols ; i++)
	      {
		if ( var_type[i] != NULL &&  var_type[i][0]== 'I') 
		  modelByColumn->setInteger(i);
	      }
	  }
	modelPrimalDual->setOptimizationDirection((sense==0) ? 1: -1);
	if ( Qmatbeg != NULL ) 
	  {	
	    modelPrimalDual->loadQuadraticObjective(ncols,(const int*) Qmatbeg->Iv,(const int*) Qmatind->Iv,Qmatval->R);
	  }	
	break;
      }
    }

  /* change handler for printing */
  mexprinter = new DerivedHandler(); 
  mexprinter->setLogLevel(options->loglevel);		 

  if ( filename != NULL )
    {
      if (options->solverchoice==3)
	{	
	  modelPrimalDual->writeMps(filename);
	}
      else
	{	
	  modelByColumn->writeMps(filename);
	}
      
    }

  if ( save_only == TRUE ) 
    {
      /* do not solve the pb just return */
      goto ok;
    }

  
  if (options->solverchoice == 3)
    {		
      modelPrimalDual->setMaximumIterations(options->maxnumiterations);
      modelPrimalDual->setMaximumSeconds(options->maxnumseconds);
      modelPrimalDual->setPrimalTolerance(options->primaltolerance);
      modelPrimalDual->setDualTolerance(options->dualtolerance);	
      modelPrimalDual->passInMessageHandler(mexprinter);		    
    }
  else
    {
      modelByColumn->setMaximumIterations(options->maxnumiterations);
      modelByColumn->setMaximumSeconds(options->maxnumseconds);
      modelByColumn->setPrimalTolerance(options->primaltolerance);
      modelByColumn->setDualTolerance(options->dualtolerance);				
      modelByColumn->passInMessageHandler(mexprinter);	
    }
  
  switch (options->solverchoice)
    {
    default:
    case 1:			
      {	
	/* Sciprintf("using primal\n");*/
	modelByColumn->primal();		
	break;
      }
    case 2:			
      {				
	/* Sciprintf("using dual\n"); */
	modelByColumn->dual();
	break;
      }
    case 3:
      {			
	/* Sciprintf("using primal-dual\n"); */
	ClpCholeskyBase * cholesky = new ClpCholeskyBase();			
	cholesky->setKKT(true);		
	modelPrimalDual->setCholesky(cholesky);	
	if (modelPrimalDual->primalDual())
	  {
	    Sciprintf("Failed\n");
	  }
	break;
      }
    }
	
  if (options->solverchoice==3)
    {	
      primal = modelPrimalDual->primalColumnSolution();
      dual = modelPrimalDual->dualRowSolution();
      RetCost->R[0] = modelPrimalDual->getObjValue(); 
      Retcode->R[0] = modelPrimalDual->status();			
      /* modelPrimalDual->writeMps("poo.mps"); */
    }
  else
    {	
      primal = modelByColumn->primalColumnSolution();
      dual = modelByColumn->dualRowSolution();
      RetCost->R[0] =  modelByColumn->getObjValue(); 
      Retcode->R[0] = modelByColumn->status();		
      /* modelByColumn->writeMps("poo.mps"); */
    }
  
  /* 
  switch ((int) Retcode->R[0]) 
  {
  case 0:  OK 
  case 1:  infeasible 
  case 2:  unbounded 
  }
  */

  /* variables */

  if (primal != NULL) 
    {
      memcpy(X->R,primal,   ncols*sizeof(double));
    }
  if (dual   != NULL) 
    { 
      /* change the order since we have here 
       * equality constraints and the inequalities 
       * and we want the same order as in glpk interface 
       */
      /* first the inequality constraints */
      memcpy(Lambda->R,dual+neq,(nrows-neq)*sizeof(double));
      /* then the equality constraints */
      memcpy(Lambda->R+(nrows-neq),dual,(neq)*sizeof(double));
    }

 ok:

  /* Delete C++ allocated objects */
  if (modelByColumn   != NULL){delete(modelByColumn);}	
  if (modelPrimalDual != NULL){delete(modelPrimalDual);}		
  if (mexprinter      != NULL){delete(mexprinter);}		
  return OK;
}

