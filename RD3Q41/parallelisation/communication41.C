#include<communication41.h>


void recvCommunicate3DFace1(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataReceive1M, PANINI_REAL *tmpDataReceive1P, int m1, int m2, int m3,procInfo3D  CartersianGrid3D)
{
 int dataSizePlane23 = 84*m2*m3;
 int tag = 0;
 MPI_Irecv(tmpDataReceive1P, dataSizePlane23, MPI_DOUBLE, myneighbours[NEB_M1_ZERO_ZERO], tag,MPI_COMM_WORLD,request);
 tag = 1;
 MPI_Irecv(tmpDataReceive1M, dataSizePlane23, MPI_DOUBLE, myneighbours[NEB_P1_ZERO_ZERO], tag,MPI_COMM_WORLD,request);
}

void sendCommunicateFace1(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataSend1P, PANINI_REAL *tmpDataSend1M, int m1, int m2, int m3,procInfo3D  CartersianGrid3D)
{
 int dataSizePlane23 = 84*m2*m3;
 int tag = 0;
  MPI_Isend(tmpDataSend1P,dataSizePlane23 , MPI_DOUBLE, myneighbours[NEB_P1_ZERO_ZERO], tag,MPI_COMM_WORLD,request);
 tag  = 1;
  MPI_Isend(tmpDataSend1M,dataSizePlane23 , MPI_DOUBLE, myneighbours[NEB_M1_ZERO_ZERO], tag,MPI_COMM_WORLD,request);
}

void recvCommunicate3DFace2(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataReceive2M, PANINI_REAL *tmpDataReceive2P, int m1, int m2, int m3,procInfo3D  CartersianGrid3D)
{
 int dataSizePlane13 = 84*m1*m3;
 int tag = 0;
 MPI_Irecv(tmpDataReceive2P, dataSizePlane13, MPI_DOUBLE, myneighbours[NEB_ZERO_M1_ZERO], tag,MPI_COMM_WORLD,request);
 tag = 1;
 MPI_Irecv(tmpDataReceive2M, dataSizePlane13, MPI_DOUBLE, myneighbours[NEB_ZERO_P1_ZERO], tag,MPI_COMM_WORLD,request);
}

void sendCommunicateFace2(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataSend2P, PANINI_REAL *tmpDataSend2M, int m1, int m2, int m3,procInfo3D  CartersianGrid3D)
{
 int dataSizePlane13 = 84*m1*m3;
 int tag = 0;
  MPI_Isend(tmpDataSend2P,dataSizePlane13 , MPI_DOUBLE, myneighbours[NEB_ZERO_P1_ZERO], tag,MPI_COMM_WORLD,request);
 tag  = 1;
  MPI_Isend(tmpDataSend2M,dataSizePlane13 , MPI_DOUBLE, myneighbours[NEB_ZERO_M1_ZERO], tag,MPI_COMM_WORLD,request);
}

void recvCommunicate3DFace3(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataReceive3M, PANINI_REAL *tmpDataReceive3P, int m1, int m2, int m3,procInfo3D  CartersianGrid3D)
{
 int dataSizePlane12 = 84*m1*m2;
 int tag = 0;
 MPI_Irecv(tmpDataReceive3P, dataSizePlane12, MPI_DOUBLE, myneighbours[NEB_ZERO_ZERO_M1], tag,MPI_COMM_WORLD,request);
 tag = 1;
 MPI_Irecv(tmpDataReceive3M, dataSizePlane12, MPI_DOUBLE, myneighbours[NEB_ZERO_ZERO_P1], tag,MPI_COMM_WORLD,request);
}

void sendCommunicateFace3(int myRank,  int *myneighbours, MPI_Request *request, PANINI_REAL *tmpDataSend3P, PANINI_REAL *tmpDataSend3M, int m1, int m2, int m3,procInfo3D  CartersianGrid3D)
{
 int dataSizePlane12 = 84*m1*m2;
 int tag = 0;
  MPI_Isend(tmpDataSend3P,dataSizePlane12 , MPI_DOUBLE, myneighbours[NEB_ZERO_ZERO_P1], tag,MPI_COMM_WORLD,request);
 tag  = 1;
  MPI_Isend(tmpDataSend3M,dataSizePlane12 , MPI_DOUBLE, myneighbours[NEB_ZERO_ZERO_M1], tag,MPI_COMM_WORLD,request);
}
