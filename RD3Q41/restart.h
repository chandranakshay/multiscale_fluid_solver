#ifndef _RESTART_H_
#define _RESTART_H_

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
void dumpRestart(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,char* restartFile);

template <int N,int numblock, typename dataType1>
void restartIC(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,char* restartFile);

#endif
