#ifndef _COLLIDE41_H_
#define _COLLIDE41_H_

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


void timestamp ( void );

template <int N,int numblock, typename dataType1>
void collide(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 twoBeta,dataType1 oneMinustwoBeta);

template <typename dataType1>
void getMomentsWithForce(lbmRD3Q41<dataType1> &lbModel,dataType1 *f, double F1, double F2, double F3);

template <int N,int numblock, typename dataType1>
void collideWithForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH, dataType1 twoBeta, dataType1 tau , dataType1 F1, dataType1 F2, dataType1 F3);

template <int N,int numblock, typename dataType1>
void copyFromNodeTo4Pt41Array(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3, dataType1 feq[][4]);

template <int N,int numblock, typename dataType1>
void copyFromCellTo4Pt41Array(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3, dataType1 feq[][4]);

template <typename dataType1>
void getfTempFSIMD(lbmRD3Q41<dataType1> &lbModel, dataType1 *rho, dataType1 F1, dataType1 F2, dataType1 F3);

template <int N,int numblock, typename dataType1>
void collideWithForceSIMD(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH, dataType1 twoBeta, dataType1 tau , dataType1 F1, dataType1 F2, dataType1 F3,dataType1 dt);

template<typename dataType1>
void calculateAlpha(lbmRD3Q41<dataType1> &lbModel,dataType1* x_i, double f_i[][41],dataType1 &alpha,int index);

template<typename dataType1>
void calculateAlphaSIMD(lbmRD3Q41<dataType1> &lbModel,
                        dataType1 x_i[][4],
                        dataType1 f_i[][4],
                        dataType1 beta,
                        dataType1* alpha);

template<int N, int numblock, typename dataType1>
void get_xi_SIMD_Node(lbmRD3Q41<dataType1> &lbModel,
                      gridBCC3D<N, numblock, dataType1> &myGrid,
                      int i1, int i2, int i3,
                      int VECT_LENGTH,
                      dataType1 x_SIMD_i[][4]);

template<int N, int numblock, typename dataType1>
void get_xi_SIMD_Cell(lbmRD3Q41<dataType1> &lbModel,
                      gridBCC3D<N, numblock, dataType1> &myGrid,
                      int VECT_LENGTH,
                      int i1, int i2, int i3,
                      dataType1 x_SIMD_i[][4]);

template <int N,int numblock, typename dataType1>
void collideWithForceEntropicSIMD(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH, dataType1 twoBeta, dataType1 tau , dataType1 F1, dataType1 F2, dataType1 F3,dataType1 dt,int step,int procRank);

#endif
