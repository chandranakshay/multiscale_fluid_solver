#ifndef _EQUILIBRIUM41_H_
#define _EQUILIBRIUM41_H_

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>
#include"gridBCC3DBasic.h"
#include"paniniAliases.h"
#include"D3RQ41.h"
#include"simdHelper.h"

template<typename dataType1>
void getFEq(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta);

template<typename dataType1>
void getFEqSIMD(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, int);

    template<typename dataType1>
    void getIterativeFEq(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta);
    
    template<typename dType>
void getFEqNewIterative(lbmRD3Q41<dType> &lbModel,int VECT_LENGTH, dType *rho, dType *uX, dType *uY, dType *uZ, dType *theta);

template<typename dataType1>
void getFEqSinglePointIntoArray(lbmRD3Q41<dataType1> &lbModel , dataType1 rho, dataType1 uX, dataType1 uY, dataType1 uZ, dataType1 theta,dataType1*  array);

template<typename dataType1>
void getGrad(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx);

    template<typename dataType1>
    void getGradSinglePointIntoArray(lbmRD3Q41<dataType1> &lbModel , dataType1 rho, dataType1 uX, dataType1 uY, dataType1 uZ, dataType1 Pxx, dataType1 Pyy, dataType1 Pzz, dataType1 Pxy, dataType1 Pyz, dataType1 Pzx ,dataType1*  array);
    
#endif
