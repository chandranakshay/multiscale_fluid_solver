 /*********************************************************************************************
 *   Copyright (c) <2015>, <Santosh Ansumali@JNCASR>                                         *
 *   All rights reserved.                                                                    *
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

#ifndef _PARTIAL_COMMUNICATION41_H_
#define _PARTIAL_COMMUNICATION41_H_

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>
#include"gridBCC3DBasic.h"
#include"paniniAliases.h"
#include"D3RQ41.h"
#include"procInfo3D.h"
#include"simdHelper.h"


void partial_recvCommunicate3DFace1(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataReceive1M, PANINI_REAL *tmpDataReceive1P, int m1, int m2, int m3,procInfo3D  CartersianGrid3D);

void partial_sendCommunicateFace1(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataSend1P, PANINI_REAL *tmpDataSend1M, int m1, int m2, int m3,procInfo3D  CartersianGrid3D);

void partial_recvCommunicate3DFace2(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataReceive2M, PANINI_REAL *tmpDataReceive2P, int m1, int m2, int m3,procInfo3D  CartersianGrid3D);

void partial_sendCommunicateFace2(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataSend2P, PANINI_REAL *tmpDataSend2M, int m1, int m2, int m3,procInfo3D  CartersianGrid3D);

void partial_recvCommunicate3DFace3(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataReceive3M, PANINI_REAL *tmpDataReceive3P, int m1, int m2, int m3,procInfo3D  CartersianGrid3D);

void partial_sendCommunicateFace3(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataSend3P, PANINI_REAL *tmpDataSend3M, int m1, int m2, int m3,procInfo3D  CartersianGrid3D);

#endif
