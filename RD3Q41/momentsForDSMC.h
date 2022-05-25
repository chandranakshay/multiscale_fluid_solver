#ifndef _MOMENTS_DSMC_H
#define _MOMENTS_DSMC_H

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
void printMoments(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 u_inlet,dataType1 dt, dataType1 beta,int step,dataType1 convectionTime,int size,int myRank, std::ofstream &file, int *coord,dataType1 F1, dataType1 F2, dataType1 F3,int centerX,int centerY, int centerZ, double f2g_factor,double length, double totLength, int countTimo, int, int, int, int , int, int);
template <int N,int numblock, typename dataType1>
void printMomentsCell(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 u_inlet,dataType1 dt, dataType1 beta,int step,dataType1 convectionTime,int size,int myRank, int *coord,dataType1 F1, dataType1 F2, dataType1 F3,int centerX,int centerY, int centerZ);
template <int N,int numblock, typename dataType1>
void printMomentsLocal(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 u_inlet,dataType1 dt, dataType1 beta,int step,dataType1 convectionTime,int size,int myRankZ, std::ofstream &file, int *coord,dataType1 F1, dataType1 F2, dataType1 F3,int centerX,int centerY, int centerZ, double f2g_factor, double length, double totLength, int countTimo);

#endif
