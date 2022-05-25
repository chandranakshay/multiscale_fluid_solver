#ifndef _PVTS_H
#define _PVTS_H

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
#include"procInfo3D.h"



template <int N,int numblock, typename dataType1>
void dumpResultsNodeToPVTI(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,procInfo3D  CartersianGrid3D,int myRank,int x1,int x2,int x3,int nP1, int nP2, int nP3,int counter,dataType1 F1, dataType1 F2, dataType1 F3, dataType1 dt);

void writeResultsNodeMasterPVTI(int x1,int x2,int x3,int nP1, int nP2, int nP3,int counter);

template <typename dataType1>
void dumpMarkerNodeToPVTI(lbmRD3Q41<dataType1> &lbModel,gridBCC3D<1, 1, int> &marker,procInfo3D  CartersianGrid3D,int myRank,int x1,int x2,int x3,int nP1, int nP2, int nP3,int counter);

void writeResultsNodeMarkerPVTI(int x1,int x2,int x3,int nP1, int nP2, int nP3,int counter);

template <int N,int numblock, typename dataType1>
void averageVelocity(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step,int size,int myRank,dataType1 F1,dataType1 F2,dataType1 F3,dataType1 uRef,dataType1 &averageVelocity1,dataType1 dt);

#endif
