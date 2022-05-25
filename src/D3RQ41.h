#ifndef _D3Q41_H_
#define _D3Q41_H_

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>
#include"gridBCC3DBasic.h"
#include"paniniAliases.h"
#include"simdHelper.h"

/** @file */


/*! \brief Defines the 41 velocity Lattice Boltzmann model.
 *  It defines the following parameters used by the model.
 *  - weights, velocities and temperature of the model
 *  - number of groups and velocities.
 *     We have 11 groups with 4 velocities in each group. The fisrt one being the Zero velocity shell.
 *  - enumerators for velocities (Global and group wise)
 *  - other variables which are very oftem computed in caculating the FEquilibrium
 *  - temperory 2D arrays to store FEq values at 4 points for each group.(4x4 2D array)
 *  - temporary variables to store density,velocity and temperature
 */

template<typename dataType>
struct lbmRD3Q41
{
    /********************************************************************************
     * Initial LBM with unit cell size of BCC                                        *
     *********************************************************************************/
    typedef dataType lbDataType;
    lbmRD3Q41(dataType _cellSize)
    {
        cellSize= _cellSize;
        //    getLatticeParameter();
    }

    dataType cellSize;
    dataType T0;
    dataType theta0;
    dataType oneByTheta0;
    dataType w0;
    dataType wSC1;
    dataType wFCC;
    dataType wSC2;
    dataType wBCC1;
    dataType wBCC;
    dataType fTemp0[4][4]  __attribute__ ((aligned(8192)));
    dataType fTemp1[4][4]  __attribute__ ((aligned(8192)));
    dataType fTemp2[4][4]  __attribute__ ((aligned(8192)));
    dataType fTemp3[4][4]  __attribute__ ((aligned(8192)));
    dataType fTemp4[4][4]  __attribute__ ((aligned(8192)));
    dataType fTemp5[4][4]  __attribute__ ((aligned(8192)));
    dataType fTemp6[4][4]  __attribute__ ((aligned(8192)));
    dataType fTemp7[4][4]  __attribute__ ((aligned(8192)));
    dataType fTemp8[4][4]  __attribute__ ((aligned(8192)));
    dataType fTemp9[4][4]  __attribute__ ((aligned(8192)));
    dataType fTemp10[4][4] __attribute__ ((aligned(8192)));

    dataType forceTemp0[4][4]  __attribute__ ((aligned(8192)));
    dataType forceTemp1[4][4]  __attribute__ ((aligned(8192)));
    dataType forceTemp2[4][4]  __attribute__ ((aligned(8192)));
    dataType forceTemp3[4][4]  __attribute__ ((aligned(8192)));
    dataType forceTemp4[4][4]  __attribute__ ((aligned(8192)));
    dataType forceTemp5[4][4]  __attribute__ ((aligned(8192)));
    dataType forceTemp6[4][4]  __attribute__ ((aligned(8192)));
    dataType forceTemp7[4][4]  __attribute__ ((aligned(8192)));
    dataType forceTemp8[4][4]  __attribute__ ((aligned(8192)));
    dataType forceTemp9[4][4]  __attribute__ ((aligned(8192)));
    dataType forceTemp10[4][4] __attribute__ ((aligned(8192)));

    dataType rho;
    dataType uX;
    dataType uY;
    dataType uZ;
    dataType theta;

    // Define enumeration for velocity set
    const static int  DV_ZERO_ZERO_ZERO    =  0;        const static int CENTER_DV_ZERO_ZERO_ZERO  = 0;

    const static int  DV_ZERO_ZERO_P1      =  1;        const static int G1_DV_ZERO_ZERO_P1        = 0;
    const static int  DV_ZERO_ZERO_P2      =  2;        const static int G1_DV_ZERO_ZERO_P2        = 1;
    const static int  DV_ZERO_P2_ZERO      =  3;        const static int G1_DV_ZERO_P2_ZERO        = 2;
    const static int  DV_P2_ZERO_ZERO      =  4;        const static int G1_DV_P2_ZERO_ZERO        = 3;

    const static int  DV_ZERO_ZERO_M1      =  5;        const static int G2_DV_ZERO_ZERO_M1        = 0;
    const static int  DV_ZERO_ZERO_M2      =  6;        const static int G2_DV_ZERO_ZERO_M2        = 1;
    const static int  DV_ZERO_M2_ZERO      =  7;        const static int G2_DV_ZERO_M2_ZERO        = 2;
    const static int  DV_M2_ZERO_ZERO      =  8;        const static int G2_DV_M2_ZERO_ZERO        = 3;

    const static int  DV_ZERO_P1_P1        =  9;        const static int G3_DV_ZERO_P1_P1          = 0;
    const static int  DV_ZERO_M1_P1        = 10;        const static int G3_DV_ZERO_M1_P1          = 1;
    const static int  DV_P1_ZERO_P1        = 11;        const static int G3_DV_P1_ZERO_P1          = 2;
    const static int  DV_M1_ZERO_P1        = 12;        const static int G3_DV_M1_ZERO_P1          = 3;

    const static int  DV_ZERO_M1_M1        = 13;        const static int G4_DV_ZERO_M1_M1          = 0;
    const static int  DV_ZERO_P1_M1        = 14;        const static int G4_DV_ZERO_P1_M1          = 1;
    const static int  DV_M1_ZERO_M1        = 15;        const static int G4_DV_M1_ZERO_M1          = 2;
    const static int  DV_P1_ZERO_M1        = 16;        const static int G4_DV_P1_ZERO_M1          = 3;

    const static int  DV_P1_P1_ZERO        = 17;        const static int G5_DV_P1_P1_ZERO          = 0;
    const static int  DV_M1_P1_ZERO        = 18;        const static int G5_DV_M1_P1_ZERO          = 1;
    const static int  DV_ZERO_P1_ZERO      = 19;        const static int G5_DV_ZERO_P1_ZERO        = 2;
    const static int  DV_P1_ZERO_ZERO      = 20;        const static int G5_DV_P1_ZERO_ZERO        = 3;

    const static int  DV_M1_M1_ZERO        = 21;        const static int G6_DV_M1_M1_ZERO          = 0;
    const static int  DV_P1_M1_ZERO        = 22;        const static int G6_DV_P1_M1_ZERO          = 1;
    const static int  DV_ZERO_M1_ZERO      = 23;        const static int G6_DV_ZERO_M1_ZERO        = 2;
    const static int  DV_M1_ZERO_ZERO      = 24;        const static int G6_DV_M1_ZERO_ZERO        = 3;

    const static int  DV_P_P_P             = 25;        const static int G7_DV_P_P_P               = 0;
    const static int  DV_M_P_P             = 26;        const static int G7_DV_M_P_P               = 1;
    const static int  DV_M_M_P             = 27;        const static int G7_DV_M_M_P               = 2;
    const static int  DV_P_M_P             = 28;        const static int G7_DV_P_M_P               = 3;

    const static int  DV_M_M_M             = 29;        const static int G8_DV_M_M_M               = 0;
    const static int  DV_P_M_M             = 30;        const static int G8_DV_P_M_M               = 1;
    const static int  DV_P_P_M             = 31;        const static int G8_DV_P_P_M               = 2;
    const static int  DV_M_P_M             = 32;        const static int G8_DV_M_P_M               = 3;

    const static int  DV_P1_P1_P1          = 33;        const static int G9_DV_P1_P1_P1            = 0;
    const static int  DV_M1_P1_P1          = 34;        const static int G9_DV_M1_P1_P1            = 1;
    const static int  DV_M1_M1_P1          = 35;        const static int G9_DV_M1_M1_P1            = 2;
    const static int  DV_P1_M1_P1          = 36;        const static int G9_DV_P1_M1_P1            = 3;

    const static int  DV_M1_M1_M1          = 37;        const static int G10_DV_M1_M1_M1           = 0;
    const static int  DV_P1_M1_M1          = 38;        const static int G10_DV_P1_M1_M1           = 1;
    const static int  DV_P1_P1_M1          = 39;        const static int G10_DV_P1_P1_M1           = 2;
    const static int  DV_M1_P1_M1          = 40;        const static int G10_DV_M1_P1_M1           = 3;

    const static int dvN = 41;
    const static int G0  =  0;
    const static int G1  =  1;
    const static int G2  =  2;
    const static int G3  =  3;
    const static int G4  =  4;
    const static int G5  =  5;
    const static int G6  =  6;
    const static int G7  =  7;
    const static int G8  =  8;
    const static int G9  =  9;
    const static int G10  =  10;

    dataType wt[dvN] __attribute__ ((aligned(4096)));
    dataType cx[dvN] __attribute__ ((aligned(4096)));
    dataType cy[dvN] __attribute__ ((aligned(4096)));
    dataType cz[dvN] __attribute__ ((aligned(4096)));
    dataType cc[dvN] __attribute__ ((aligned(4096)));
    dataType cxcy[dvN] __attribute__ ((aligned(4096)));
    dataType cycz[dvN] __attribute__ ((aligned(4096)));
    dataType czcx[dvN] __attribute__ ((aligned(4096)));
    dataType cx2[dvN] __attribute__ ((aligned(4096)));
    dataType cy2[dvN] __attribute__ ((aligned(4096)));
    dataType cz2[dvN] __attribute__ ((aligned(4096)));
    dataType cxcsq[dvN] __attribute__ ((aligned(4096)));
    dataType cycsq[dvN] __attribute__ ((aligned(4096)));
    dataType czcsq[dvN] __attribute__ ((aligned(4096)));
    dataType csq2[dvN] __attribute__ ((aligned(4096)));
    dataType fTemp[dvN] __attribute__ ((aligned(4096)));
    dataType yMinus3_CENTER;
    dataType yTenMinusYSqrPlus15_CENTER;
    dataType yMinus3_SC2;
    dataType yTenMinusYSqrPlus15_SC2;
    dataType yMinus3_SC1;
    dataType yTenMinusYSqrPlus15_SC1;
    dataType yMinus3_FCC;
    dataType yTenMinusYSqrPlus15_FCC;
    dataType yMinus3_BCC;
    dataType yTenMinusYSqrPlus15_BCC;
    dataType yMinus3_BCC1;
    dataType yTenMinusYSqrPlus15_BCC1;


};


template <typename dataType>
void getLatticeParameter(lbmRD3Q41<dataType>& lbModel);

template <typename dataType>
void setModelParameters(lbmRD3Q41<dataType>& lbModel);

template <int N,int numblock, typename dataType1>
void copyFromNode(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3);

template <int N,int numblock, typename dataType1>
void copyToNode(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3);

template <int N,int numblock, typename dataType1>
void copyFromCell(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3);

template <int N,int numblock, typename dataType1>
void copyToCell(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3);

template <int N,int numblock, typename dataType1>
void copyToNodeSinglePoint(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3,int point);

template <int N,int numblock, typename dataType1>
void copyToCellSinglePoint(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3,int point);

template <int N,int numblock, typename dataType1>
void copyFromNodeSinglePoint(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3,int point);

template <int N,int numblock, typename dataType1>
void copyFromCellSinglePoint(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3,int point);


template <int N,int numblock, typename dataType1>
void getHydroMomentsFromNode(lbmRD3Q41<dataType1> &lbModel,gridBCC3D<N, numblock, dataType1> &myGrid,int i1, int i2, int i3, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta);

template <int N,int numblock, typename dataType1>
void getHydroMomentsFromNode(lbmRD3Q41<dataType1> &lbModel,gridBCC3D<N, numblock, dataType1> &myGrid,int i1, int i2, int i3, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta);

template <int N,int numblock, typename dataType1>
void getHydroMomentsFromCell(lbmRD3Q41<dataType1> &lbModel,gridBCC3D<N, numblock, dataType1> &myGrid,int i1, int i2, int i3, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta);

template <int N,int numblock, typename dataType1>
void getHydroMomentsFromNodeWithForce(lbmRD3Q41<dataType1> &lbModel,gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta,dataType1 F1, dataType1 F2, dataType1 F3,dataType1 dt);

template <int N,int numblock, typename dataType1>
void getHydroMomentsFromCellWithForce(lbmRD3Q41<dataType1> &lbModel,gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta,dataType1 F1, dataType1 F2, dataType1 F3,dataType1 dt);

template <int N,int numblock, typename dataType1>
void globalMassTotalFluid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step, std::ofstream &file,int size,int myRank);

template <int N,int numblock, typename dataType1>
void globalMassTotalSolid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step, std::ofstream &file,int size,int myRank);

template <int N,int numblock, typename dataType1>
void globalMassInternal(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step, std::ofstream &file,int size,int myRank, int *coord);

template <int N,int numblock, typename dataType1>
void globalMassExternal(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step, std::ofstream &file,int size,int myRank, int *coord);

template<typename dataType1>
void getGrad11Moment(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx);

template<typename dataType1>
void getEnstrophyMoments(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx, dataType1 *Sxx, dataType1 *Syy, dataType1 *Szz, dataType1 *Sxy, dataType1 *Syz, dataType1 *Szx,dataType1 dt, dataType1 beta);

template<typename dataType1>
void getEnstrophyMomentsSIMD(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx, dataType1 *Sxx, dataType1 *Syy, dataType1 *Szz, dataType1 *Sxy, dataType1 *Syz, dataType1 *Szx,dataType1 dt, dataType1 beta);

template <int N,int numblock, typename dataType1>
void printEnstrophy(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 u_inlet,dataType1 dt, dataType1 beta,int step,dataType1 convectionTime,int size,int myRank,dataType1 &rhoSum,dataType1 &thetasum,dataType1 &kineticEnergy,dataType1 &internalEnergy,dataType1 &totalEnergy,dataType1 &enstrophySum,dataType1 &enstrophyMax);

template<typename dataType1>
void getGrad11MomentSIMD(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx);

template <int N,int numblock, typename dataType1>
void copyFromNodeToArray(lbmRD3Q41<dataType1> &lbModel ,dataType1 *grad, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3);

template <int N,int numblock, typename dataType1>
void copyFromCellToArray(lbmRD3Q41<dataType1> &lbModel ,dataType1 *grad, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3);

template<typename dataType1>
void getFEqSinglePointIntoArray(lbmRD3Q41<dataType1> &lbModel , dataType1 rho, dataType1 uX, dataType1 uY, dataType1 uZ, dataType1 theta,dataType1*  array);


template <int N,int numblock, typename dataType1>
void printTransverseKineticEnergy(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 u_inlet,dataType1 dt, dataType1 beta,int step,dataType1 convectionTime,int size,int myRank, std::ofstream &file,int yLocation, int *coord,dataType1 F1, dataType1 F2, dataType1 F3, int coresYtot);

template<typename dataType1>
void getHydroMomentSinglePoint(lbmRD3Q41<dataType1> &lbModel, dataType1 &rho, dataType1 &uX, dataType1 &uY, dataType1 &uZ, dataType1 &theta);

template<typename dataType1>
  void getHydroMomentSinglePointWithForce(lbmRD3Q41<dataType1> &lbModel, dataType1 &rho, dataType1 &uX, dataType1 &uY, dataType1 &uZ, dataType1 &theta,dataType1 F1, dataType1 F2, dataType1 F3,dataType1 dt );

template<int N,int numblock, typename dataType1>
void printVelocityAtTheCenterOfDomain(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int *coord, dataType1 F1, dataType1 F2, dataType1 F3,int centerX, int centerY, int centerZ, std::ofstream &file,int step,dataType1 convectionTime,dataType1 dt) ;

template <int N,int numblock, typename dataType1>
void globalMassTotalFluid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step,int size,int myRank);

template <int N,int numblock, typename dataType1>
void globalMassTotalSolid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step,int size,int myRank);

template<typename dataType1>
void getGrad11Moment_heatAlpha(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx, dataType1 *q1, dataType1 *q2, dataType1 *q3);

template<typename dataType1>
void getGrad11Moment_heatAlpha_YminusFlux(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx, dataType1 *q1, dataType1 *q2, dataType1 *q3);

template<typename dataType1>
void getGrad11Moment_heatAlpha_YplusFlux(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx, dataType1 *q1, dataType1 *q2, dataType1 *q3);

template <int N,int numblock, typename dataType1>
void calcLBTopFlux(lbmRD3Q41<dataType1> &lbModel41,int VECT_LENGTH, gridBCC3D<N, numblock, dataType1> &gridLBM41, int LBpointsX, int LBpointsY, int numMoments, int numSpaceAvgCells_LB_X, int numSpaceAvgCells_LB_Y, int yIndex, dataType1* momentsLBFlux, int myRank);

template <int N,int numblock, typename dataType1>
void calcLBBotFlux(lbmRD3Q41<dataType1> &lbModel41,int VECT_LENGTH, gridBCC3D<N, numblock, dataType1> &gridLBM41, int LBpointsX, int LBpointsY, int numMoments, int numSpaceAvgCells_LB_X, int numSpaceAvgCells_LB_Y, int yIndex, dataType1* momentsLBFlux);

#endif
