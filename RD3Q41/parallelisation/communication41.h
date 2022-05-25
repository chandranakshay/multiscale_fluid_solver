#ifndef _COMMUNICATION41_H_
#define _COMMUNICATION41_H_

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

void recvCommunicate3DFace1(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataReceive1M, PANINI_REAL *tmpDataReceive1P, int m1, int m2, int m3,procInfo3D  CartersianGrid3D);

void sendCommunicateFace1(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataSend1P, PANINI_REAL *tmpDataSend1M, int m1, int m2, int m3,procInfo3D  CartersianGrid3D);

void recvCommunicate3DFace2(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataReceive2M, PANINI_REAL *tmpDataReceive2P, int m1, int m2, int m3,procInfo3D  CartersianGrid3D);

void sendCommunicateFace2(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataSend2P, PANINI_REAL *tmpDataSend2M, int m1, int m2, int m3,procInfo3D  CartersianGrid3D);

void recvCommunicate3DFace3(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataReceive3M, PANINI_REAL *tmpDataReceive3P, int m1, int m2, int m3,procInfo3D  CartersianGrid3D);

void sendCommunicateFace3(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataSend3P, PANINI_REAL *tmpDataSend3M, int m1, int m2, int m3,procInfo3D  CartersianGrid3D);

#endif
