#ifndef _IO_STREAM_GRID_H_
#define _IO_STREAM_GRID_H_

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
void probeComputationalGridPlane12Node(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i3, int k);

template <int N,int numblock, typename dataType1>
void probeComputationalGridPlane12Cell(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i3, int k);

template <int N,int numblock, typename dataType1>
void probeComputationalGridPlane23Node(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i1, int k);

template <int N,int numblock, typename dataType1>
void probeComputationalGridPlane23Cell(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i1, int k);

template <int N,int numblock, typename dataType1>
void probeComputationalGridPlane31Node(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i2, int k);

template <int N,int numblock, typename dataType1>
void probeComputationalGridPlane31Cell(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i2, int k);

template <int N,int numblock, typename dataType1>
void modifyComputationalGridPlane12Node(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i3, int k);

template <int N,int numblock, typename dataType1>
void modifyComputationalGridPlane12Cell(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i3, int k);

template <int N,int numblock, typename dataType1>
void modifyComputationalGridPlane23Node(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i1, int k);

template <int N,int numblock, typename dataType1>
void modifyComputationalGridPlane23Cell(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i1, int k);

template <int N,int numblock, typename dataType1>
void modifyComputationalGridPlane31Node(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i2, int k);

template <int N,int numblock, typename dataType1>
void modifyComputationalGridPlane31Cell(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i2, int k);

#endif
