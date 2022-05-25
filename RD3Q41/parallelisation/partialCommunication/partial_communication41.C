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
#include<partial_communication41.h>


void partial_recvCommunicate3DFace1(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataReceive1M, PANINI_REAL *tmpDataReceive1P, int m1, int m2, int m3,procInfo3D  CartersianGrid3D)
{
 int dataSizePlane23 = 30*m2*m3;
//std::cout<<"dataSize: "<<m1<<"\t"<<m2<<"\t"<<m3<<std::endl;
 int tag = 0;
 MPI_Irecv(tmpDataReceive1P, dataSizePlane23, MPI_DOUBLE, myneighbours[NEB_M1_ZERO_ZERO], tag,MPI_COMM_WORLD,&request[0]);
 tag = 1;
 MPI_Irecv(tmpDataReceive1M, dataSizePlane23, MPI_DOUBLE, myneighbours[NEB_P1_ZERO_ZERO], tag,MPI_COMM_WORLD,&request[1]);
}

void partial_sendCommunicateFace1(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataSend1P, PANINI_REAL *tmpDataSend1M, int m1, int m2, int m3,procInfo3D  CartersianGrid3D)
{
 int dataSizePlane23 = 30*m2*m3;
 int tag = 0;
  MPI_Isend(tmpDataSend1P,dataSizePlane23 , MPI_DOUBLE, myneighbours[NEB_P1_ZERO_ZERO], tag,MPI_COMM_WORLD,&request[0]);
 tag  = 1;
  MPI_Isend(tmpDataSend1M,dataSizePlane23 , MPI_DOUBLE, myneighbours[NEB_M1_ZERO_ZERO], tag,MPI_COMM_WORLD,&request[1]);
}

void partial_recvCommunicate3DFace2(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataReceive2M, PANINI_REAL *tmpDataReceive2P, int m1, int m2, int m3,procInfo3D  CartersianGrid3D)
{
 int dataSizePlane13 = 30*m1*m3;
 int tag = 0;
 MPI_Irecv(tmpDataReceive2P, dataSizePlane13, MPI_DOUBLE, myneighbours[NEB_ZERO_M1_ZERO], tag,MPI_COMM_WORLD,&request[0]);
 tag = 1;
 MPI_Irecv(tmpDataReceive2M, dataSizePlane13, MPI_DOUBLE, myneighbours[NEB_ZERO_P1_ZERO], tag,MPI_COMM_WORLD,&request[1]);
}

void partial_sendCommunicateFace2(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataSend2P, PANINI_REAL *tmpDataSend2M, int m1, int m2, int m3,procInfo3D  CartersianGrid3D)
{
 int dataSizePlane13 = 30*m1*m3;
 int tag = 0;
  MPI_Isend(tmpDataSend2P,dataSizePlane13 , MPI_DOUBLE, myneighbours[NEB_ZERO_P1_ZERO], tag,MPI_COMM_WORLD,&request[0]);
 tag  = 1;
  MPI_Isend(tmpDataSend2M,dataSizePlane13 , MPI_DOUBLE, myneighbours[NEB_ZERO_M1_ZERO], tag,MPI_COMM_WORLD,&request[1]);
}

void partial_recvCommunicate3DFace3(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataReceive3M, PANINI_REAL *tmpDataReceive3P, int m1, int m2, int m3,procInfo3D  CartersianGrid3D)
{
 int dataSizePlane12 = 30*m1*m2;
 int tag = 0;
 MPI_Irecv(tmpDataReceive3P, dataSizePlane12, MPI_DOUBLE, myneighbours[NEB_ZERO_ZERO_M1], tag,MPI_COMM_WORLD,&request[0]);
 tag = 1;
 MPI_Irecv(tmpDataReceive3M, dataSizePlane12, MPI_DOUBLE, myneighbours[NEB_ZERO_ZERO_P1], tag,MPI_COMM_WORLD,&request[1]);
}

void partial_sendCommunicateFace3(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataSend3P, PANINI_REAL *tmpDataSend3M, int m1, int m2, int m3,procInfo3D  CartersianGrid3D)
{
 int dataSizePlane12 = 30*m1*m2;
 int tag = 0;
  MPI_Isend(tmpDataSend3P,dataSizePlane12 , MPI_DOUBLE, myneighbours[NEB_ZERO_ZERO_P1], tag,MPI_COMM_WORLD,&request[0]);
 tag  = 1;
  MPI_Isend(tmpDataSend3M,dataSizePlane12 , MPI_DOUBLE, myneighbours[NEB_ZERO_ZERO_M1], tag,MPI_COMM_WORLD,&request[1]);
}
