/*********************************************************************************************
 * Copyright (c) <2015>, <Santosh Ansumali@JNCASR>                                           *
 *  All rights reserved.                                                                     *
 *   Redistribution and use in source and binary forms, with or without modification, are    *
 *   permitted provided that the following conditions are met:                               *
 *                                                                                           *
 *    1. Redistributions of source code must retain the above copyright notice, this list of *
 *       conditions and the following disclaimer.                                            *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list *
 *       of conditions and the following disclaimer in the documentation and/or other        *
 *       materials provided with the distribution.                                           *
 *    3. Neither the name of the <JNCASR> nor the names of its contributors may be used to   *
 *       endorse or promote products derived from this software without specific prior       *
 *       written permission.                                                                 *
 *                                                                                           *
 *       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND     *
 *       ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED       *
 *       WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  *
 *       IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,    *
 *       INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,      *
 *       BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,       *
 *       DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF     *
 *       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE     *
 *       OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED   *
 *       OF THE POSSIBILITY OF SUCH DAMAGE.                                                  *
 *                                                                                           *
 *       Suggestions:          ansumali@jncasr.ac.in                                         *
 *       Bugs:                 ansumali@jncasr.ac.in                                         *
 *                                                                                           *
 *********************************************************************************************/

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


#ifndef _PARTIAL_PACK41_H_
#define _PARTIAL_PACK41_H_


PANINI_REAL* partial_allocateForPacking(int m1, int m2);

template <int N,int numblock, typename dataType1>
void partial_packTo1P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void partial_unpackFrom1P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void partial_packTo1M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void partial_unpackFrom1M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void partial_packMSgTo1Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp1Psend, dataType1 *temp1Msend);

template <int N,int numblock, typename dataType1> 
void partial_recvMSgFrom1Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp1Preceive, dataType1 *temp1Mreceive);

template <int N,int numblock, typename dataType1>
void partial_packTo2P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void partial_unpackFrom2P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void partial_packTo2M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void partial_unpackFrom2M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void partial_packMSgTo2Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp2Psend, dataType1 *temp2Msend);

template <int N,int numblock, typename dataType1> 
void partial_recvMSgFrom2Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp2Preceive, dataType1 *temp2Mreceive);

template <int N,int numblock, typename dataType1>
void partial_packTo3P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void partial_unpackFrom3P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1>
void partial_packTo3M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void partial_unpackFrom3M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index);

template <int N,int numblock, typename dataType1> 
void partial_packMSgTo3Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp3Psend, dataType1 *temp3Msend);

template <int N,int numblock, typename dataType1> 
void partial_recvMSgFrom3Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp3Preceive, dataType1 *temp3Mreceive);


#endif
