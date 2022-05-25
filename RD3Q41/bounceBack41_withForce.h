#ifndef _BOUNCEBACK41_WITHFORCE_H_
#define _BOUNCEBACK41_WITHFORCE_H_

/** @file */

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>
#include"gridBCC3DBasic.h"
#include"paniniAliases.h"
#include"D3RQ41.h"
#include"simdHelper.h"

  template <typename dataType1> 
  inline void SWAP(dataType1 &a,dataType1 &b);
  
  template <int N,int numblock, typename dataType1>
  void bounceBackG12(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);
  
  template <int N,int numblock, typename dataType1>
  void bounceBackG34(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);

  template <int N,int numblock, typename dataType1>
  void bounceBackG56(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);

  template <int N,int numblock, typename dataType1>
  void bounceBackG78(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);

  template <int N,int numblock, typename dataType1>
  void bounceBackG910(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);

  template <int N,int numblock, typename dataType1>
  void bounceBackG12Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);        

  template <int N,int numblock, typename dataType1>
  void bounceBackG34Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);        

  template <int N,int numblock, typename dataType1>
  void bounceBackG56Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);        

  template <int N,int numblock, typename dataType1>
  void bounceBackG78Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);        
  
  template <int N,int numblock, typename dataType1>
  void bounceBackG910Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);        
  
template <int N,int numblock, typename dataType1>                                                                                                          
void bounceBackWithForces(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ,int step);        
  
  template <int N,int numblock, typename dataType1>
  void bounceBackG12Solid_beforeAdvection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);        

  template <int N,int numblock, typename dataType1>
  void bounceBackG34Solid_beforeAdvection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);        

  template <int N,int numblock, typename dataType1>
  void bounceBackG56Solid_beforeAdvection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);        

  template <int N,int numblock, typename dataType1>
  void bounceBackG78Solid_beforeAdvection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);        
  
  template <int N,int numblock, typename dataType1>
  void bounceBackG910Solid_beforeAdvection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);        
  
template <int N,int numblock, typename dataType1>                                                                                                          
void bounceBackWithForces_beforeAdvection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ,int step);      


 #endif
