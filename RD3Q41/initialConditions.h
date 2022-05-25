#ifndef _INITIALCONDITIONS_H_
#define _INITIALCONDITIONS_H_

/** @file */

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>
#include"gridBCC3DBasic.h"
#include"paniniAliases.h"
#include"D3RQ41.h"
#include"procInfo3D.h"
#include"simdHelper.h"
#include"equilibrium41.h"


  template <int N,int numblock, typename dataType1>
  void initialConditions(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,dataType1 u_inlet,int FLUID,int nP2, int *coord,dataType1 rhoDSMC);
  
  template <int N,int numblock, typename dataType1>
  void initializeKidaParallel(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 u_inlet, int *coord,int nP1,int nP2,int nP3);

#endif
