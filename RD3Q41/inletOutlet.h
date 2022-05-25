#ifndef _INLETOUTLET_H_
#define _INLETOUTLET_H_

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>
#include"gridBCC3DBasic.h"
#include"paniniAliases.h"
#include"D3RQ41.h"
#include"simdHelper.h"
#include"equilibrium41.h"
#include"boundaryCondition.h"

  template <int N,int numblock, typename dataType1>
  void X1BeginInlet(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 u_inlet, dataType1 cs,int step,int VECT_LENGTH);

  template <int N,int numblock, typename dataType1>
  void X1EndOutlet(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH);

#endif
