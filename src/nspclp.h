#ifndef NSP_CLP_CPP_H 
#define NSP_CLP_CPP_H 

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

#ifdef __cplusplus
extern "C" {
#endif 

#include <nsp/nsp.h>
#include <nsp/objects.h>
#include <nsp/imatrix.h>

  typedef struct _nsp_clp_params nsp_clp_params ;
  
  struct _nsp_clp_params
  {
    int solverchoice ;
    int maxnumiterations ;
    int loglevel ;
    int primalpivot ;
    int dualpivot ;
    double maxnumseconds ;
    double primaltolerance ;
    double dualtolerance ;
  };

  /* exported functions from C++ which are used in C */
  
  extern double nsp_coin_dbl_max(void);

  extern int nsp_clp_solve(nsp_clp_params *options,int sense, int ncols, int nrows, int neq, 
			   NspIMatrix*Cmatbeg, NspIMatrix *Cmatind, NspMatrix *Cmatval, 
			   NspMatrix *lower, NspMatrix *upper, NspMatrix *Objective,
			   NspIMatrix*Qmatbeg, NspIMatrix *Qmatind, NspMatrix *Qmatval, 
			   NspMatrix *Rhs,  NspMatrix *Lhs, char *var_type[], NspMatrix *X,NspMatrix *Lambda,
			   NspMatrix *RetCost,NspMatrix *Retcode,
			   const char *filename,int save_only);

#ifdef __cplusplus
}
#endif 

#endif /* CPL_CPP_H  */


