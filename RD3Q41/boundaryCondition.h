#ifndef _BOUNDARYCONDITION_H_
#define _BOUNDARYCONDITION_H_

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
    void topWall(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);

    template <int N,int numblock, typename dataType1>
    void bottomWall(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid);

    template<typename dataType1>
    void getGradThermalSinglePoint(lbmRD3Q41<dataType1> &lbModel , dataType1 rho, dataType1 uX, dataType1 uY, dataType1 uZ, dataType1 theta, dataType1 pXX, dataType1 pYY, dataType1 pZZ, dataType1 pXY, dataType1 pYZ, dataType1 pZX ,dataType1 *grad);

    template <int N,int numblock, typename dataType1>
    void copyFromArrayToNode(lbmRD3Q41<dataType1> &lbModel ,dataType1 *grad, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3);

    template <int N,int numblock, typename dataType1>
    void copyFromArrayToCell(lbmRD3Q41<dataType1> &lbModel ,dataType1 *grad, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3);

    template <int N,int numblock, typename dataType1>
    void copyFromArrayToNode2(lbmRD3Q41<dataType1> &lbModel ,dataType1 *grad, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3,int myRank);

    template <int N,int numblock, typename dataType1>
    void copyFromArrayToCell2(lbmRD3Q41<dataType1> &lbModel ,dataType1 *grad, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3,int myRank);

    template<typename dataType1>
    void getGradThermalSinglePoint_heatAlpha(lbmRD3Q41<dataType1> &lbModel , dataType1 rho, dataType1 uX, dataType1 uY, dataType1 uZ, dataType1 theta, dataType1 pXX, dataType1 pYY, dataType1 pZZ, dataType1 pXY, dataType1 pYZ, dataType1 pZX, dataType1 q1, dataType1 q2, dataType1 q3,dataType1 *grad);



#endif
