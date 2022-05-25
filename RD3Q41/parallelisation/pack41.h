#include<stdio.h>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>
#include"gridBCC3DBasic.h"
#include"paniniAliases.h"
#include"D3RQ41.h"
#include"simdHelper.h"  
#include"ioStreamGrid.h"  


#ifndef _PACK41_H_
#define _PACK41_H_


PANINI_REAL* allocateForPacking(int m1, int m2);

template <int N,int numblock, typename dataType1>
void packTo1P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void unpackFrom1P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void packTo1M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void unpackFrom1M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void packMSgTo1Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp1Psend, dataType1 *temp1Msend);

template <int N,int numblock, typename dataType1> 
void recvMSgFrom1Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp1Preceive, dataType1 *temp1Mreceive);

template <int N,int numblock, typename dataType1>
void packTo2P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void unpackFrom2P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void packTo2M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void unpackFrom2M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void packMSgTo2Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp2Psend, dataType1 *temp2Msend);

template <int N,int numblock, typename dataType1> 
void recvMSgFrom2Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp2Preceive, dataType1 *temp2Mreceive);

template <int N,int numblock, typename dataType1>
void packTo3P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void unpackFrom3P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void packTo3M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void unpackFrom3M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void packMSgTo3Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp3Psend, dataType1 *temp3Msend);

template <int N,int numblock, typename dataType1> 
void recvMSgFrom3Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp3Preceive, dataType1 *temp3Mreceive);

#endif
