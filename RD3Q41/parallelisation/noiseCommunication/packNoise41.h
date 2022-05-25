#ifndef _PACK41_NOISE_H_
#define _PACK41_NOISE_H_

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
#include"ioStreamGrid_noise.h"  
#include"pack41.h"  
#include"communicationNoise.h"  



PANINI_REAL* allocateForPacking_Noise(int m1, int m2);

template <int N,int numblock, typename dataType1>
void packTo1P_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void unpackFrom1P_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void packTo1M_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void unpackFrom1M_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void packMSgTo1Neb_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp1Psend, dataType1 *temp1Msend);

template <int N,int numblock, typename dataType1> 
void recvMSgFrom1Neb_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp1Preceive, dataType1 *temp1Mreceive);

template <int N,int numblock, typename dataType1>
void packTo2P_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void unpackFrom2P_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void packTo2M_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void unpackFrom2M_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void packMSgTo2Neb_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp2Psend, dataType1 *temp2Msend);

template <int N,int numblock, typename dataType1> 
void recvMSgFrom2Neb_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp2Preceive, dataType1 *temp2Mreceive);

template <int N,int numblock, typename dataType1>
void packTo3P_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void unpackFrom3P_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void packTo3M_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void unpackFrom3M_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void packMSgTo3Neb_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp3Psend, dataType1 *temp3Msend);

template <int N,int numblock, typename dataType1> 
void recvMSgFrom3Neb_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp3Preceive, dataType1 *temp3Mreceive);


#endif
