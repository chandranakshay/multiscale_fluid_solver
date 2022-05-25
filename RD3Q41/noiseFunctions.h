#ifndef _NOISEFUNCTIONS_H_
#define _NOISEFUNCTIONS_H_

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>
#include"gridBCC3DBasic.h"
#include"paniniAliases.h"
#include"D3RQ41.h"
#include"equilibrium41.h"
#include"procInfo3D.h"
#include"perlinNoise.h"

// template <int N,int numblock, typename dataType1>
// void curl_noise_communicateAfterThis(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &testGrid, int nP1, int nP2, int nP3, dataType1 &normg,int *coordinates,int VECT_LENGTH,double Lx, double Ly, double Lz);
//
// template <int N,int numblock, int numblock1, typename dataType1>
// void curl_noise_communiateBeforeThis(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, gridBCC3D<N, numblock1, dataType1> &testGrid, int nP1, int nP2, int nP3, dataType1 &normg,int *coordinates,int VECT_LENGTH,double Lx, double Ly, double Lz);
//
// template <int N,int numblock, typename dataType1>
// void averageVelocityAndTKE(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int nP1, int nP2, int nP3,int *coordinates,int VECT_LENGTH,dataType1 F1,dataType1 F2,dataType1 F3,int step,int myRank,int size,dataType1 U0);
//
template <int N,int numblock, typename dataType1>
void perturb_parab(int *coordinates, lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int nP1, int nP2, int nP3, dataType1 normg,dataType1 U0,int VECT_LENGTH, dataType1 density, dataType1 uRef);


template <int N,int numblock, typename dataType1>
void curl_noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int nP1, int nP2, int nP3, dataType1 &normg,int *coordinates,int VECT_LENGTH, dataType1 density, dataType1 uRef);


#endif
