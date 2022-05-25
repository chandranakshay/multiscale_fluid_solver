#ifndef _ADVECTION41_H_
#define _ADVECTION41_H_

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>
#include"gridBCC3DBasic.h"
#include"paniniAliases.h"
#include"D3RQ41.h"
#include"simdHelper.h"

template <int N,int numblock, typename dataType1>
void advectionG1(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);            // ZZP,ZZP2,ZP2Z,P2ZZ
template <int N,int numblock, typename dataType1>
void advectionG2(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);            // ZZM,ZZM2,ZM2Z,M2ZZ
template <int N,int numblock, typename dataType1>                                                
void advectionG3(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);          // ZPP,ZMP,PZP,MZP 
template <int N,int numblock, typename dataType1>
void advectionG4(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);           // ZPM,ZMM,PZM,MZM
template <int N,int numblock, typename dataType1>
void advectionG5(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);          // ZPM,ZMM,PZM,MZM
template <int N,int numblock, typename dataType1>
void advectionG6(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);            // ZZP,ZZT,ZTZ,TZZ
template <int N,int numblock, typename dataType1>
void advectionG7(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);          // ZPM,ZMM,PZM,MZM
template <int N,int numblock, typename dataType1>
void advectionG8(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);          // ZPM,ZMM,PZM,MZM
template <int N,int numblock, typename dataType1>
void advectionG9(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);           //PPP,MPP,MMP,PMP
template <int N,int numblock, typename dataType1>
void advectionG10(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);           //PPP,MPP,MMP,PMP
template <int N,int numblock, typename dataType1>
void advection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);

#endif
