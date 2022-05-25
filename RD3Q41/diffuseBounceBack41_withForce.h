#ifndef _DIFFUSEBOUNCEBACK_WITHFORCE_H_
#define _DIFFUSEBOUNCEBACK_WITHFORCE_H_

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
void getRhoNodeSinglePoint(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,dataType1 &rho, int i1,int i2,int i3);
template <int N,int numblock, typename dataType1>
void getRhoCellSinglePoint(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,dataType1 &rho, int i1,int i2,int i3);
template <typename dataType1>
void getFEqG12(lbmRD3Q41<dataType1> &lbModel,dataType1 rho,dataType1 uX,dataType1 uY,dataType1 uZ,dataType1 theta);
template <typename dataType1>
void getFEqG34(lbmRD3Q41<dataType1> &lbModel,dataType1 rho,dataType1 uX,dataType1 uY,dataType1 uZ,dataType1 theta);
template <typename dataType1>
void getFEqG56(lbmRD3Q41<dataType1> &lbModel,dataType1 rho,dataType1 uX,dataType1 uY,dataType1 uZ,dataType1 theta);
template <typename dataType1>
void getFEqG78910(lbmRD3Q41<dataType1> &lbModel,dataType1 rho,dataType1 uX,dataType1 uY,dataType1 uZ,dataType1 theta);
template <int N,int numblock, typename dataType1>
void getRhoForDiffuseBounceBackNodeG12(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID);            
template <int N,int numblock, typename dataType1>
void getRhoForDiffuseBounceBackNodeG34(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID);            
template <int N,int numblock, typename dataType1>
void getRhoForDiffuseBounceBackNodeG56(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID);            
template <int N,int numblock, typename dataType1>
void getRhoForDiffuseBounceBackNodeG910(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID);            
template <int N,int numblock, typename dataType1>
void getRhoForDiffuseBounceBackNodeCellG78(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID);            
template <int N,int numblock, typename dataType1>
void getRhoForDiffuseBounceBackCellG12(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID);           
template <int N,int numblock, typename dataType1>
void getRhoForDiffuseBounceBackCellG34(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID);            
template <int N,int numblock, typename dataType1>
void getRhoForDiffuseBounceBackCellG56(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID);            
template <int N,int numblock, typename dataType1>
void getRhoForDiffuseBounceBackCellG910(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID);            
template <int N,int numblock, typename dataType1>
void diffusebounceBackNodeG12withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                
template <int N,int numblock, typename dataType1>
void diffusebounceBackNodeG34withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                
template <int N,int numblock, typename dataType1>
void diffusebounceBackNodeG56withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                
template <int N,int numblock, typename dataType1>
void diffusebounceBackNodeG910withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);              
template <int N,int numblock, typename dataType1>
void diffusebounceBackCellG12withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                
template <int N,int numblock, typename dataType1>
void diffusebounceBackCellG34withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);               
template <int N,int numblock, typename dataType1>
void diffusebounceBackCellG56withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                 
template <int N,int numblock, typename dataType1>
void diffusebounceBackCellG910withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);               
template <int N,int numblock, typename dataType1>
void diffusebounceBackNodeCellG78withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                
template <int N,int numblock, typename dataType1> 
void diffusebounceBackwithForces(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);            
template <int N,int numblock, typename dataType1>
void getRhoAll(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, dataType1> &rhoGrid);            
template <int N,int numblock, typename dataType1>
void diffusebounceBackNodeG12Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                
template <int N,int numblock, typename dataType1>
void diffusebounceBackCellG12Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                
template <int N,int numblock, typename dataType1>
void diffusebounceBackNodeG34Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                
template <int N,int numblock, typename dataType1>
void diffusebounceBackCellG34Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                
template <int N,int numblock, typename dataType1>
void diffusebounceBackNodeG56Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                 
template <int N,int numblock, typename dataType1>
void diffusebounceBackCellG56Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                 
template <int N,int numblock, typename dataType1>
void diffusebounceBackNodeG910Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);              
template <int N,int numblock, typename dataType1>
void diffusebounceBackCellG910Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);              
template <int N,int numblock, typename dataType1>
void diffusebounceBackNodeCellG78Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);                
template <int N,int numblock, typename dataType1> 
void diffusebounceBackSolid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ);            
#endif
