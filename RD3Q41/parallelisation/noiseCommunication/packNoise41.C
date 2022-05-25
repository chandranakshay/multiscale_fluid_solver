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


#include<packNoise41.h>


PANINI_REAL* allocateForPacking_Noise(int m1, int m2)
{
 PANINI_REAL *temp;  
 temp = new PANINI_REAL[m1*m2*4];
 return temp;
} 


template <int N,int numblock, typename dataType1>
void packTo1P_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n2*myGrid.n3;  
 int stride = 1; 

 /////////////
 //   Node_Noise  //
 ///////////// 
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
}

template <int N,int numblock, typename dataType1> 
void unpackFrom1P_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n2*myGrid.n3;  
 int stride = 1; 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
}

template <int N,int numblock, typename dataType1>
void packTo1M_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n2*myGrid.n3;  
 int stride = 1; 

 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
}

template <int N,int numblock, typename dataType1> 
void unpackFrom1M_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n2*myGrid.n3;  
 int stride = 1; 
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane23Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
}

template <int N,int numblock, typename dataType1> 
void packMSgTo1Neb_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp1Psend, dataType1 *temp1Msend)
{
 int index = 0;
 packTo1P_Noise(lbModel,myGrid,temp1Psend,index);  
 index = 0;
 packTo1M_Noise(lbModel,myGrid,temp1Msend,index); 
}

template <int N,int numblock, typename dataType1> 
void recvMSgFrom1Neb_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp1Preceive, dataType1 *temp1Mreceive)
{
 int index = 0;
 unpackFrom1P_Noise(lbModel,myGrid,temp1Preceive,index);  
 index = 0;
 unpackFrom1M_Noise(lbModel,myGrid,temp1Mreceive,index); 
}




template <int N,int numblock, typename dataType1>
void packTo2P_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n3;  
 int stride = 1; 

 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
}

template <int N,int numblock, typename dataType1> 
void unpackFrom2P_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n3;  
 int stride = 1; 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
}

template <int N,int numblock, typename dataType1>
void packTo2M_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n3;  
 int stride = 1; 
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB2,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB2,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB2,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB2,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
}

template <int N,int numblock, typename dataType1> 
void unpackFrom2M_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n3;  
 int stride = 1; 
  ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
}



template <int N,int numblock, typename dataType1> 
void packMSgTo2Neb_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp2Psend, dataType1 *temp2Msend)
{
 int index = 0;
 packTo2P_Noise(lbModel,myGrid,temp2Psend,index);  
 index = 0;
 packTo2M_Noise(lbModel,myGrid,temp2Msend,index); 
}

template <int N,int numblock, typename dataType1> 
void recvMSgFrom2Neb_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp2Preceive, dataType1 *temp2Mreceive)
{
 int index = 0;
 unpackFrom2P_Noise(lbModel,myGrid,temp2Preceive,index);  
 index = 0;
 unpackFrom2M_Noise(lbModel,myGrid,temp2Mreceive,index); 
}



template <int N,int numblock, typename dataType1>
void packTo3P_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n2;  
 int stride = 1; 
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
}

template <int N,int numblock, typename dataType1> 
void unpackFrom3P_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n2;  
 int stride = 1; 
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 }

template <int N,int numblock, typename dataType1>
void packTo3M_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n2;  
 int stride = 1; 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB3,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size;                                                                        
 probeComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB3,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;                                                                        
 probeComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB3,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB3,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
}

template <int N,int numblock, typename dataType1> 
void unpackFrom3M_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n2;  
 int stride = 1; 

 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size;                                                                         
 modifyComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;                                                                         
 modifyComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane12Node_Noise(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
}

template <int N,int numblock, typename dataType1> 
void packMSgTo3Neb_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp3Psend, dataType1 *temp3Msend)
{
 int index = 0;
 packTo3P_Noise(lbModel,myGrid,temp3Psend,index);  
 index = 0;
 packTo3M_Noise(lbModel,myGrid,temp3Msend,index); 
}

template <int N,int numblock, typename dataType1> 
void recvMSgFrom3Neb_Noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp3Preceive, dataType1 *temp3Mreceive)
{
 int index = 0;
 unpackFrom3P_Noise(lbModel,myGrid,temp3Preceive,index);  
 index = 0;
 unpackFrom3M_Noise(lbModel,myGrid,temp3Mreceive,index); 
}


// Explicit declaration
template void packTo1P_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *,int &);
template void unpackFrom1P_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *,int &);
template void packTo1M_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *,int &);
template void unpackFrom1M_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *,int &);
template void packMSgTo1Neb_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *, double *);
template void recvMSgFrom1Neb_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *, double *);

template void packTo2P_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *,int &);
template void unpackFrom2P_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *,int &);
template void packTo2M_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *,int &);
template void unpackFrom2M_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *,int &);
template void packMSgTo2Neb_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *, double *);
template void recvMSgFrom2Neb_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *, double *);

template void packTo3P_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *,int &);
template void unpackFrom3P_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *,int &);
template void packTo3M_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *,int &);
template void unpackFrom3M_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *,int &);
template void packMSgTo3Neb_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *, double *);
template void recvMSgFrom3Neb_Noise<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, double *, double *);

