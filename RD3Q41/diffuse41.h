#ifndef _DIFFUSE_H_
#define _DIFFUSE_H_

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
void prepareDiffuseTopWall(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);

template <int N,int numblock, typename dataType1>
void correctDiffuseF0massManipTopWall(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int VECT_LENGTH, dataType1 tempTop, double uRef);

template <int N,int numblock, typename dataType1>
void prepareDiffuseBottomWall(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)   ;

template <int N,int numblock, typename dataType1>
void correctDiffuseF0massManipBottomWall(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int VECT_LENGTH, dataType1 tempBottom, double uRef);


template<typename dataType1>
void getHydroMomentSinglePoint(lbmRD3Q41<dataType1> &lbModel, dataType1 &rho, dataType1 &uX, dataType1 &uY, dataType1 &uZ, dataType1 &theta);
 
 
template <int N,int numblock, typename dataType1>
void getTempProfile(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int step, double visc, double delT,  int myRank,  double time);

template <int N,int numblock, typename dataType1>
void globalMass(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int step, int myRank);



#endif
