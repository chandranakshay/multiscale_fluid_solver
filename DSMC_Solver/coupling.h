/*********************************************************************************************
 *   Copyright (c) <2018>, <Santosh Ansumali@JNCASR>                                         *
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


/*********************************************************************************************

 *  @Author: Akshay Chandran and Santosh Ansumali                                            *

 *********************************************************************************************/

#ifndef COUPLING_H_INCLUDED
#define COUPLING_H_INCLUDED

#include "hardSphere.h"
#include "gridBCC3DBasic.h"
#include "D3RQ41.h"
#include <mpi.h>

template<int DIM,typename real,typename dof>
void calcSendLBmomentsLowerBufferZone(hardSphere<DIM,real,dof>& p, lbmRD3Q41<real> &lbModel41, gridBCC3D<4, 11, real> &gridLBM41, real* moments_LB_Send,
		real weightLB_P1, real weightLB_P2, int LBbyDSMCpointsX, int LBbyDSMCpointsY, int numSpaceAvgCells_LB_X, int numSpaceAvgCells_LB_Y, int numSpaceAvgCells_LB, 
		int numSpaceAvgCells_DSMC_X, int numSpaceAvgCells_DSMC_Y, int numSpaceAvgCells_DSMC, int myRank, int myRankX, int myRankY, int myRankZ, int myRankX_LB, int myRankY_LB, 
		int myRankZ_LB, int count, int lowerCoupleDSMC)
{
  int VECT_LENGTH = 4;
  double rho_LB_1Temp2[4],uX_LB_1Temp2[4],uY_LB_1Temp2[4],uZ_LB_1Temp2[4],theta_LB_1Temp2[4],Pxx_LB_1Temp2[4],Pyy_LB_1Temp2[4],Pzz_LB_1Temp2[4],Pxy_LB_1Temp2[4],Pyz_LB_1Temp2[4],Pzx_LB_1Temp2[4];
  double q1_LB_1Temp2[4],q2_LB_1Temp2[4],q3_LB_1Temp2[4];
  for(int cellY = 0; cellY < (p.YcellSizePerCore/LBbyDSMCpointsY); cellY = cellY + numSpaceAvgCells_LB_Y)
  {
    for(int cellX = 0; cellX < (p.XcellSizePerCore/LBbyDSMCpointsX); cellX = cellX + numSpaceAvgCells_LB_X)
    {
      int cellNumberLocal = ((int)(p.XcellSizePerCore/(numSpaceAvgCells_DSMC_X)))*((int)(cellY/numSpaceAvgCells_LB_Y)) + (int)(cellX/numSpaceAvgCells_LB_X);

      // std::cout<<"dsmcParticleFlux1: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<p.dsmcParticleFlux[cellNumberLocal]<<"\t"<<((double)p.dsmcParticleFlux[cellNumberLocal]/(p.volume*numSpaceAvgCells_DSMC_X*numSpaceAvgCells_DSMC_Y))<<std::endl;

      moments_LB_Send[p.numMoments*cellNumberLocal     ] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal +  1] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal +  2] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal +  3] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal +  4] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal +  5] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal +  6] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal +  7] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal +  8] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal +  9] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal + 10] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal + 11] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal + 12] = 0.;
      moments_LB_Send[p.numMoments*cellNumberLocal + 13] = 0.;

      // Difference in notations of LB and DSMC
      // LB(x,y,z) = DSMC(y,z,x)

      for(int i = 0; i < numSpaceAvgCells_LB_X; i++)
      {
        for(int j = 0; j < numSpaceAvgCells_LB_Y; j++)
        {
          int xIndex = cellY + gridLBM41.nB1 + j;
          int zIndex = cellX + gridLBM41.nB3 + i;
          int yIndex =         gridLBM41.nB2 + 3;//(0.5*p.nCoresZ - 1);

          copyFromNode(lbModel41,gridLBM41,VECT_LENGTH,xIndex,yIndex,zIndex);
          getGrad11Moment_heatAlpha(lbModel41,VECT_LENGTH,rho_LB_1Temp2,uX_LB_1Temp2,uY_LB_1Temp2,uZ_LB_1Temp2,theta_LB_1Temp2,Pxx_LB_1Temp2,Pyy_LB_1Temp2,Pzz_LB_1Temp2,Pxy_LB_1Temp2,Pyz_LB_1Temp2,Pzx_LB_1Temp2,q1_LB_1Temp2,q2_LB_1Temp2,q3_LB_1Temp2);

          moments_LB_Send[p.numMoments*cellNumberLocal     ] += weightLB_P2 * rho_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  1] += weightLB_P2 * uX_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  2] += weightLB_P2 * uY_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  3] += weightLB_P2 * uZ_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  4] += weightLB_P2 * theta_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  5] += weightLB_P2 * Pxx_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  6] += weightLB_P2 * Pyy_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  7] += weightLB_P2 * Pzz_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  8] += weightLB_P2 * Pxy_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  9] += weightLB_P2 * Pyz_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal + 10] += weightLB_P2 * Pzx_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal + 11] += weightLB_P2 * q1_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal + 12] += weightLB_P2 * q2_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal + 13] += weightLB_P2 * q3_LB_1Temp2[0];

          if((cellNumberLocal == 0) && (myRankY_LB == 5)  && (myRankX_LB == 5))
          {
              std::cout<<"recvdFromLB1: "<<myRankZ_LB<<"\t"<<rho_LB_1Temp2[0]<<"\t"<<weightLB_P2 * rho_LB_1Temp2[0]<<"\t"<<weightLB_P2<<std::endl;
          }

          copyFromNode(lbModel41,gridLBM41,VECT_LENGTH,xIndex,yIndex - 1,zIndex);
          getGrad11Moment_heatAlpha(lbModel41,VECT_LENGTH,rho_LB_1Temp2,uX_LB_1Temp2,uY_LB_1Temp2,uZ_LB_1Temp2,theta_LB_1Temp2,Pxx_LB_1Temp2,Pyy_LB_1Temp2,Pzz_LB_1Temp2,Pxy_LB_1Temp2,Pyz_LB_1Temp2,Pzx_LB_1Temp2,q1_LB_1Temp2,q2_LB_1Temp2,q3_LB_1Temp2);

          moments_LB_Send[p.numMoments*cellNumberLocal     ] += weightLB_P1 * rho_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  1] += weightLB_P1 * uX_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  2] += weightLB_P1 * uY_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  3] += weightLB_P1 * uZ_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  4] += weightLB_P1 * theta_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  5] += weightLB_P1 * Pxx_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  6] += weightLB_P1 * Pyy_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  7] += weightLB_P1 * Pzz_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  8] += weightLB_P1 * Pxy_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal +  9] += weightLB_P1 * Pyz_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal + 10] += weightLB_P1 * Pzx_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal + 11] += weightLB_P1 * q1_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal + 12] += weightLB_P1 * q2_LB_1Temp2[0];
          moments_LB_Send[p.numMoments*cellNumberLocal + 13] += weightLB_P1 * q3_LB_1Temp2[0];

          if((cellNumberLocal == 0) && (myRankY_LB == 5)  && (myRankX_LB == 5))
          {
              std::cout<<"recvdFromLB2: "<<myRankZ_LB<<"\t"<<rho_LB_1Temp2[0]<<"\t"<<weightLB_P1 * rho_LB_1Temp2[0]<<"\t"<<weightLB_P1<<std::endl;
          }
        }
      }

      moments_LB_Send[p.numMoments*cellNumberLocal     ] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal +  1] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal +  2] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal +  3] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal +  4] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal +  5] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal +  6] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal +  7] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal +  8] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal +  9] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal + 10] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal + 11] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal + 12] /= ((double)(numSpaceAvgCells_LB));
      moments_LB_Send[p.numMoments*cellNumberLocal + 13] /= ((double)(numSpaceAvgCells_LB));

      double rho_LB_1TempSpaceAvgYminus   = 0.;
      // double uX_LB_1TempSpaceAvgYminus    = 0.;
      // double uY_LB_1TempSpaceAvgYminus    = 0.;
      // double uZ_LB_1TempSpaceAvgYminus    = 0.;
      // double theta_LB_1TempSpaceAvgYminus = 0.;
      // double Pxx_LB_1TempSpaceAvgYminus   = 0.;
      // double Pyy_LB_1TempSpaceAvgYminus   = 0.;
      // double Pzz_LB_1TempSpaceAvgYminus   = 0.;
      // double Pxy_LB_1TempSpaceAvgYminus   = 0.;
      // double Pyz_LB_1TempSpaceAvgYminus   = 0.;
      // double Pzx_LB_1TempSpaceAvgYminus   = 0.;
      // double q1_LB_1TempSpaceAvgYminus    = 0.;
      // double q2_LB_1TempSpaceAvgYminus    = 0.;
      // double q3_LB_1TempSpaceAvgYminus    = 0.;

      for(int i = 0; i < numSpaceAvgCells_LB_X; i++)
      {
        for(int j = 0; j < numSpaceAvgCells_LB_Y; j++)
        {
          int xIndex = cellY + gridLBM41.nB1 + j;
          int zIndex = cellX + gridLBM41.nB3 + i;
          int yIndex =         gridLBM41.nB2 + 2;//(0.5*p.nCoresZ - 1);

          copyFromNode(lbModel41,gridLBM41,VECT_LENGTH,xIndex,yIndex,zIndex);
          getGrad11Moment_heatAlpha_YminusFlux(lbModel41,VECT_LENGTH,rho_LB_1Temp2,uX_LB_1Temp2,uY_LB_1Temp2,uZ_LB_1Temp2,theta_LB_1Temp2,Pxx_LB_1Temp2,Pyy_LB_1Temp2,Pzz_LB_1Temp2,Pxy_LB_1Temp2,Pyz_LB_1Temp2,Pzx_LB_1Temp2,q1_LB_1Temp2,q2_LB_1Temp2,q3_LB_1Temp2);

          rho_LB_1TempSpaceAvgYminus    += weightLB_P2 * rho_LB_1Temp2[0];
          // uX_LB_1TempSpaceAvgYminus     += weightLB_P2 * uX_LB_1Temp2[0];
          // uY_LB_1TempSpaceAvgYminus     += weightLB_P2 * uY_LB_1Temp2[0];
          // uZ_LB_1TempSpaceAvgYminus     += weightLB_P2 * uZ_LB_1Temp2[0];
          // theta_LB_1TempSpaceAvgYminus  += weightLB_P2 * theta_LB_1Temp2[0];
          // Pxx_LB_1TempSpaceAvgYminus    += weightLB_P2 * Pxx_LB_1Temp2[0];
          // Pyy_LB_1TempSpaceAvgYminus    += weightLB_P2 * Pyy_LB_1Temp2[0];
          // Pzz_LB_1TempSpaceAvgYminus    += weightLB_P2 * Pzz_LB_1Temp2[0];
          // Pxy_LB_1TempSpaceAvgYminus    += weightLB_P2 * Pxy_LB_1Temp2[0];
          // Pyz_LB_1TempSpaceAvgYminus    += weightLB_P2 * Pyz_LB_1Temp2[0];
          // Pzx_LB_1TempSpaceAvgYminus    += weightLB_P2 * Pzx_LB_1Temp2[0];
          // q1_LB_1TempSpaceAvgYminus     += weightLB_P2 * q1_LB_1Temp2[0];
          // q2_LB_1TempSpaceAvgYminus     += weightLB_P2 * q2_LB_1Temp2[0];
          // q3_LB_1TempSpaceAvgYminus     += weightLB_P2 * q3_LB_1Temp2[0];

          copyFromNode(lbModel41,gridLBM41,VECT_LENGTH,xIndex,yIndex - 1,zIndex);
          getGrad11Moment_heatAlpha_YminusFlux(lbModel41,VECT_LENGTH,rho_LB_1Temp2,uX_LB_1Temp2,uY_LB_1Temp2,uZ_LB_1Temp2,theta_LB_1Temp2,Pxx_LB_1Temp2,Pyy_LB_1Temp2,Pzz_LB_1Temp2,Pxy_LB_1Temp2,Pyz_LB_1Temp2,Pzx_LB_1Temp2,q1_LB_1Temp2,q2_LB_1Temp2,q3_LB_1Temp2);

          rho_LB_1TempSpaceAvgYminus    += weightLB_P1 * rho_LB_1Temp2[0];
          // uX_LB_1TempSpaceAvgYminus     += weightLB_P1 * uX_LB_1Temp2[0];
          // uY_LB_1TempSpaceAvgYminus     += weightLB_P1 * uY_LB_1Temp2[0];
          // uZ_LB_1TempSpaceAvgYminus     += weightLB_P1 * uZ_LB_1Temp2[0];
          // theta_LB_1TempSpaceAvgYminus  += weightLB_P1 * theta_LB_1Temp2[0];
          // Pxx_LB_1TempSpaceAvgYminus    += weightLB_P1 * Pxx_LB_1Temp2[0];
          // Pyy_LB_1TempSpaceAvgYminus    += weightLB_P1 * Pyy_LB_1Temp2[0];
          // Pzz_LB_1TempSpaceAvgYminus    += weightLB_P1 * Pzz_LB_1Temp2[0];
          // Pxy_LB_1TempSpaceAvgYminus    += weightLB_P1 * Pxy_LB_1Temp2[0];
          // Pyz_LB_1TempSpaceAvgYminus    += weightLB_P1 * Pyz_LB_1Temp2[0];
          // Pzx_LB_1TempSpaceAvgYminus    += weightLB_P1 * Pzx_LB_1Temp2[0];
          // q1_LB_1TempSpaceAvgYminus     += weightLB_P1 * q1_LB_1Temp2[0];
          // q2_LB_1TempSpaceAvgYminus     += weightLB_P1 * q2_LB_1Temp2[0];
          // q3_LB_1TempSpaceAvgYminus     += weightLB_P1 * q3_LB_1Temp2[0];
        }
      }

      rho_LB_1TempSpaceAvgYminus    /= ((double)(numSpaceAvgCells_LB));
      // uX_LB_1TempSpaceAvgYminus     /= ((double)(numSpaceAvgCells_LB));
      // uY_LB_1TempSpaceAvgYminus     /= ((double)(numSpaceAvgCells_LB));
      // uZ_LB_1TempSpaceAvgYminus     /= ((double)(numSpaceAvgCells_LB));
      // theta_LB_1TempSpaceAvgYminus  /= ((double)(numSpaceAvgCells_LB));
      // Pxx_LB_1TempSpaceAvgYminus    /= ((double)(numSpaceAvgCells_LB));
      // Pyy_LB_1TempSpaceAvgYminus    /= ((double)(numSpaceAvgCells_LB));
      // Pzz_LB_1TempSpaceAvgYminus    /= ((double)(numSpaceAvgCells_LB));
      // Pxy_LB_1TempSpaceAvgYminus    /= ((double)(numSpaceAvgCells_LB));
      // Pyz_LB_1TempSpaceAvgYminus    /= ((double)(numSpaceAvgCells_LB));
      // Pzx_LB_1TempSpaceAvgYminus    /= ((double)(numSpaceAvgCells_LB));
      // q1_LB_1TempSpaceAvgYminus     /= ((double)(numSpaceAvgCells_LB));
      // q2_LB_1TempSpaceAvgYminus     /= ((double)(numSpaceAvgCells_LB));
      // q3_LB_1TempSpaceAvgYminus     /= ((double)(numSpaceAvgCells_LB));

      PANINI_REAL dsmcDensityFlux = ((double)p.dsmcParticleFlux[cellNumberLocal]/(p.volume*numSpaceAvgCells_DSMC_X*numSpaceAvgCells_DSMC_Y));
      // std::cout<<"dsmcDensityFlux: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<p.dsmcParticleFlux[cellNumberLocal]<<"\t"<<dsmcDensityFlux<<std::endl;

      // moments_LB_Send[p.numMoments*cellNumberLocal     ] += (rho_LB_1TempSpaceAvgYminus  *DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal +  1] += (uX_LB_1TempSpaceAvgYminus   *DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal +  2] += (uY_LB_1TempSpaceAvgYminus   *DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal +  3] += (uZ_LB_1TempSpaceAvgYminus   *DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal +  4] += (theta_LB_1TempSpaceAvgYminus*DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal +  5] += (Pxx_LB_1TempSpaceAvgYminus  *DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal +  6] += (Pyy_LB_1TempSpaceAvgYminus  *DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal +  7] += (Pzz_LB_1TempSpaceAvgYminus  *DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal +  8] += (Pxy_LB_1TempSpaceAvgYminus  *DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal +  9] += (Pyz_LB_1TempSpaceAvgYminus  *DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal + 10] += (Pzx_LB_1TempSpaceAvgYminus  *DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal + 11] += (q1_LB_1TempSpaceAvgYminus   *DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal + 12] += (q2_LB_1TempSpaceAvgYminus   *DSMC_CountInv);
      // moments_LB_Send[p.numMoments*cellNumberLocal + 13] += (q3_LB_1TempSpaceAvgYminus   *DSMC_CountInv);



      if((cellNumberLocal == 0) && (myRankY_LB == ((int)(0.5*p.nCoresY))) && (myRankX_LB == ((int)(0.5*p.nCoresX))))
      {
        std::cout<<"====================================================================================================================="<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 0<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal     ]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 1<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  1]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 2<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  2]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 3<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  3]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 4<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  4]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 5<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  5]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 6<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  6]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 7<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  7]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 8<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  8]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 9<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  9]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<<10<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal + 10]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<<11<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal + 11]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<<12<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal + 12]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"LB2DSMC_moments_lower: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<<13<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal + 13]<<"\t"<<dsmcDensityFlux<<std::endl;
        std::cout<<"====================================================================================================================="<<std::endl;
      }

    }
  }

}

template<int DIM,typename real,typename dof>
void calcRecvLBmomentsLowerBufferZone(hardSphere<DIM,real,dof>& p, lbmRD3Q41<real> &lbModel41, gridBCC3D<4, 11, real> &gridLBM41, real* moments_LB_Recv, real* momentsSave,
		real c, real f2g_factor, real weightLB_P1, real weightLB_P2, int LBbyDSMCpointsX, int LBbyDSMCpointsY, int numSpaceAvgCells_LB_X, int numSpaceAvgCells_LB_Y, int numSpaceAvgCells_LB, 
		int numSpaceAvgCells_DSMC_X, int numSpaceAvgCells_DSMC_Y, int numSpaceAvgCells_DSMC, int myRank, int myRankX, int myRankY, int myRankZ, int myRankX_LB, int myRankY_LB, 
		int myRankZ_LB, int count, int lowerCoupleDSMC)
{
   for(int cellY = 0; cellY < (p.YcellSizePerCore/LBbyDSMCpointsY); cellY = cellY + numSpaceAvgCells_LB_Y)
   {
     for(int cellX = 0; cellX < (p.XcellSizePerCore/LBbyDSMCpointsX); cellX = cellX + numSpaceAvgCells_LB_X)
     {
       int cellNumberLocal = ((int)(p.XcellSizePerCore/(numSpaceAvgCells_DSMC_X)))*((int)(cellY/numSpaceAvgCells_LB_Y)) + (int)(cellX/numSpaceAvgCells_LB_X);

       double rho_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal     ];
       double uX_LB_1TempSpaceAvg    = moments_LB_Recv[p.numMoments*cellNumberLocal +  1];
       double uY_LB_1TempSpaceAvg    = moments_LB_Recv[p.numMoments*cellNumberLocal +  2];
       double uZ_LB_1TempSpaceAvg    = moments_LB_Recv[p.numMoments*cellNumberLocal +  3];
       double theta_LB_1TempSpaceAvg = moments_LB_Recv[p.numMoments*cellNumberLocal +  4];
       double Pxx_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal +  5];
       double Pyy_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal +  6];
       double Pzz_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal +  7];
       double Pxy_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal +  8];
       double Pyz_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal +  9];
       double Pzx_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal + 10];
       double q1_LB_1TempSpaceAvg    = moments_LB_Recv[p.numMoments*cellNumberLocal + 11];
       double q2_LB_1TempSpaceAvg    = moments_LB_Recv[p.numMoments*cellNumberLocal + 12];
       double q3_LB_1TempSpaceAvg    = moments_LB_Recv[p.numMoments*cellNumberLocal + 13];

       PANINI_REAL rho_DSMC,uX_DSMC,uY_DSMC,uZ_DSMC,sigmaXX_DSMC,sigmaYY_DSMC,sigmaZZ_DSMC,sigmaXY_DSMC,sigmaYZ_DSMC,sigmaZX_DSMC,p_DSMC,qX_DSMC,qY_DSMC,qZ_DSMC;

       rho_DSMC       = rho_LB_1TempSpaceAvg;

       p_DSMC         = rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg*c;

       uX_DSMC        = uZ_LB_1TempSpaceAvg*c;
       uY_DSMC        = uX_LB_1TempSpaceAvg*c;
       uZ_DSMC        = uY_LB_1TempSpaceAvg*c;

       PANINI_REAL u2 = uX_LB_1TempSpaceAvg*uX_LB_1TempSpaceAvg + uY_LB_1TempSpaceAvg*uY_LB_1TempSpaceAvg + uZ_LB_1TempSpaceAvg*uZ_LB_1TempSpaceAvg;

       sigmaXX_DSMC   = (Pzz_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*uZ_LB_1TempSpaceAvg*uZ_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg)*f2g_factor;
       sigmaYY_DSMC   = (Pxx_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*uX_LB_1TempSpaceAvg*uX_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg)*f2g_factor;
       sigmaYZ_DSMC   = (Pxy_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*uX_LB_1TempSpaceAvg*uY_LB_1TempSpaceAvg)*f2g_factor;
       sigmaZZ_DSMC   = (Pyy_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*uY_LB_1TempSpaceAvg*uY_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg)*f2g_factor;
       sigmaXY_DSMC   = (Pzx_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*uZ_LB_1TempSpaceAvg*uX_LB_1TempSpaceAvg)*f2g_factor;
       sigmaZX_DSMC   = (Pyz_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*uY_LB_1TempSpaceAvg*uZ_LB_1TempSpaceAvg)*f2g_factor;

       qX_DSMC = ((q3_LB_1TempSpaceAvg) - ((rho_LB_1TempSpaceAvg*u2*uZ_LB_1TempSpaceAvg) + 5.0*rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg*uZ_LB_1TempSpaceAvg + 2.0*(uX_LB_1TempSpaceAvg*sigmaXY_DSMC + uY_LB_1TempSpaceAvg*sigmaZX_DSMC + uZ_LB_1TempSpaceAvg*sigmaXX_DSMC)/f2g_factor))*(c*c*c)*f2g_factor;
       qY_DSMC = ((q1_LB_1TempSpaceAvg) - ((rho_LB_1TempSpaceAvg*u2*uX_LB_1TempSpaceAvg) + 5.0*rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg*uX_LB_1TempSpaceAvg + 2.0*(uX_LB_1TempSpaceAvg*sigmaYY_DSMC + uY_LB_1TempSpaceAvg*sigmaYZ_DSMC + uZ_LB_1TempSpaceAvg*sigmaXY_DSMC)/f2g_factor))*(c*c*c)*f2g_factor;
       qZ_DSMC = ((q2_LB_1TempSpaceAvg) - ((rho_LB_1TempSpaceAvg*u2*uY_LB_1TempSpaceAvg) + 5.0*rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg*uY_LB_1TempSpaceAvg + 2.0*(uX_LB_1TempSpaceAvg*sigmaYZ_DSMC + uY_LB_1TempSpaceAvg*sigmaZZ_DSMC + uZ_LB_1TempSpaceAvg*sigmaZX_DSMC)/f2g_factor))*(c*c*c)*f2g_factor;

       sigmaXX_DSMC   *= (c*c);
       sigmaXY_DSMC   *= (c*c);
       sigmaZX_DSMC   *= (c*c);
       sigmaYY_DSMC   *= (c*c);
       sigmaYZ_DSMC   *= (c*c);
       sigmaZZ_DSMC   *= (c*c);


       if((cellNumberLocal == 0) && (myRankY == ((int)(0.5*p.nCoresY))) && (myRankX == ((int)(0.5*p.nCoresX))))
       {
         std::cout<<"====================================================================================================================="<<std::endl;
         cout<<"rho_DSMC_lower@: "<<count<<" "<<myRank<<":"<<rho_DSMC<<std::endl;
         cout<<"uX_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<uX_DSMC<<std::endl;
         cout<<"uY_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<uY_DSMC<<std::endl;
         cout<<"uZ_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<uZ_DSMC<<std::endl;
         cout<<"t0_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<theta_LB_1TempSpaceAvg<<std::endl;
         cout<<"sXX_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<sigmaXX_DSMC<<std::endl;
         cout<<"SXY_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<sigmaXY_DSMC<<std::endl;
         cout<<"sYZ_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<sigmaYZ_DSMC<<std::endl;
         cout<<"sYY_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<sigmaYY_DSMC<<std::endl;
         cout<<"sZZ_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<sigmaZZ_DSMC<<std::endl;
         cout<<"sZX_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<sigmaZX_DSMC<<std::endl;
         cout<<"qX_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<qX_DSMC<<std::endl;
         cout<<"qY_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<qY_DSMC<<std::endl;
         cout<<"qZ_DSMC_lower@ : "<<count<<" "<<myRank<<":"<<qZ_DSMC<<std::endl;
         std::cout<<"====================================================================================================================="<<std::endl;
       }

       momentsSave[p.numMoments*cellNumberLocal     ] = rho_DSMC    ;
       momentsSave[p.numMoments*cellNumberLocal +  1] = p_DSMC      ;
       momentsSave[p.numMoments*cellNumberLocal +  2] = uX_DSMC     ;
       momentsSave[p.numMoments*cellNumberLocal +  3] = uY_DSMC     ;
       momentsSave[p.numMoments*cellNumberLocal +  4] = uZ_DSMC     ;
       momentsSave[p.numMoments*cellNumberLocal +  5] = sigmaXX_DSMC;
       momentsSave[p.numMoments*cellNumberLocal +  6] = sigmaXY_DSMC;
       momentsSave[p.numMoments*cellNumberLocal +  7] = sigmaZX_DSMC;
       momentsSave[p.numMoments*cellNumberLocal +  8] = sigmaYY_DSMC;
       momentsSave[p.numMoments*cellNumberLocal +  9] = sigmaYZ_DSMC;
       momentsSave[p.numMoments*cellNumberLocal + 10] = sigmaZZ_DSMC;
       momentsSave[p.numMoments*cellNumberLocal + 11] = qX_DSMC     ;
       momentsSave[p.numMoments*cellNumberLocal + 12] = qY_DSMC     ;
       momentsSave[p.numMoments*cellNumberLocal + 13] = qZ_DSMC     ;

       // p.totalMomentsWithFlux(momentsSave,cellNumberLocal,myRank,count);
       // p.generateDistributionFromMoments(&momentsSave[p.numMoments*cellNumberLocal],cellX/numSpaceAvgCells_LB_X,cellY/numSpaceAvgCells_LB_Y,myRank,count);

     }
   }

   // for(int cellY = 0; cellY < (YcellSizePerCore/LBbyDSMCpointsY); cellY = cellY + numSpaceAvgCells_LB_Y)
   // {
   //   for(int cellX = 0; cellX < (XcellSizePerCore/LBbyDSMCpointsX); cellX = cellX + numSpaceAvgCells_LB_X)
   //   {
   //     int cellNumberLocal = ((int)(XcellSizePerCore/(numSpaceAvgCells_DSMC_X)))*((int)(cellY/numSpaceAvgCells_LB_Y)) + (int)(cellX/numSpaceAvgCells_LB_X);
   //
   //     p.poissonParticleGen(&momentsSave[p.numMoments*cellNumberLocal],cellX/numSpaceAvgCells_LB_X,cellY/numSpaceAvgCells_LB_Y,myRank,tempCount);
   //
   //   }
   // }
}

template<int DIM,typename real,typename dof>
void calcSendLBmomentsUpperBufferZone(hardSphere<DIM,real,dof>& p, lbmRD3Q41<real> &lbModel41, gridBCC3D<4, 11, real> &gridLBM41, real* moments_LB_Send,
		real weightLB_P1, real weightLB_P2, int LBbyDSMCpointsX, int LBbyDSMCpointsY, int numSpaceAvgCells_LB_X, int numSpaceAvgCells_LB_Y, int numSpaceAvgCells_LB, 
		int numSpaceAvgCells_DSMC_X, int numSpaceAvgCells_DSMC_Y, int numSpaceAvgCells_DSMC, int myRank, int myRankX, int myRankY, int myRankZ, int myRankX_LB, int myRankY_LB, 
		int myRankZ_LB, int count, int upperCoupleDSMC)
{
   int VECT_LENGTH = 4;
   double rho_LB_1Temp2[4],uX_LB_1Temp2[4],uY_LB_1Temp2[4],uZ_LB_1Temp2[4],theta_LB_1Temp2[4],Pxx_LB_1Temp2[4],Pyy_LB_1Temp2[4],Pzz_LB_1Temp2[4],Pxy_LB_1Temp2[4],Pyz_LB_1Temp2[4],Pzx_LB_1Temp2[4];
   double q1_LB_1Temp2[4],q2_LB_1Temp2[4],q3_LB_1Temp2[4];
   for(int cellY = 0; cellY < (p.YcellSizePerCore/LBbyDSMCpointsY); cellY = cellY + numSpaceAvgCells_LB_Y)
   {
     for(int cellX = 0; cellX < (p.XcellSizePerCore/LBbyDSMCpointsX); cellX = cellX + numSpaceAvgCells_LB_X)
     {
       int cellNumberLocal = ((int)(p.XcellSizePerCore/(numSpaceAvgCells_DSMC_X)))*((int)(cellY/numSpaceAvgCells_LB_Y)) + (int)(cellX/numSpaceAvgCells_LB_X);

       // std::cout<<"dsmcParticleFlux1: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<p.dsmcParticleFlux[cellNumberLocal]<<"\t"<<((double)p.dsmcParticleFlux[cellNumberLocal]/(p.volume*numSpaceAvgCells_DSMC_X*numSpaceAvgCells_DSMC_Y))<<std::endl;
       moments_LB_Send[p.numMoments*cellNumberLocal     ] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal +  1] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal +  2] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal +  3] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal +  4] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal +  5] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal +  6] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal +  7] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal +  8] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal +  9] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal + 10] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal + 11] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal + 12] = 0.;
       moments_LB_Send[p.numMoments*cellNumberLocal + 13] = 0.;

       // Difference in notations of LB and DSMC
       // LB(x,y,z) = DSMC(y,z,x)

       for(int i = 0; i < numSpaceAvgCells_LB_X; i++)
       {
         for(int j = 0; j < numSpaceAvgCells_LB_Y; j++)
         {
           int xIndex = cellY + gridLBM41.nB1 + j;
           int zIndex = cellX + gridLBM41.nB3 + i;
           int yIndex =         gridLBM41.nE2 - 3;//(0.5*p.nCoresZ - 1);

           copyFromNode(lbModel41,gridLBM41,VECT_LENGTH,xIndex,yIndex,zIndex);
           getGrad11Moment_heatAlpha(lbModel41,VECT_LENGTH,rho_LB_1Temp2,uX_LB_1Temp2,uY_LB_1Temp2,uZ_LB_1Temp2,theta_LB_1Temp2,Pxx_LB_1Temp2,Pyy_LB_1Temp2,Pzz_LB_1Temp2,Pxy_LB_1Temp2,Pyz_LB_1Temp2,Pzx_LB_1Temp2,q1_LB_1Temp2,q2_LB_1Temp2,q3_LB_1Temp2);

           moments_LB_Send[p.numMoments*cellNumberLocal     ] += weightLB_P2 * rho_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  1] += weightLB_P2 * uX_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  2] += weightLB_P2 * uY_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  3] += weightLB_P2 * uZ_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  4] += weightLB_P2 * theta_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  5] += weightLB_P2 * Pxx_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  6] += weightLB_P2 * Pyy_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  7] += weightLB_P2 * Pzz_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  8] += weightLB_P2 * Pxy_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  9] += weightLB_P2 * Pyz_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal + 10] += weightLB_P2 * Pzx_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal + 11] += weightLB_P2 * q1_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal + 12] += weightLB_P2 * q2_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal + 13] += weightLB_P2 * q3_LB_1Temp2[0];

           copyFromNode(lbModel41,gridLBM41,VECT_LENGTH,xIndex,yIndex + 1,zIndex);
           getGrad11Moment_heatAlpha(lbModel41,VECT_LENGTH,rho_LB_1Temp2,uX_LB_1Temp2,uY_LB_1Temp2,uZ_LB_1Temp2,theta_LB_1Temp2,Pxx_LB_1Temp2,Pyy_LB_1Temp2,Pzz_LB_1Temp2,Pxy_LB_1Temp2,Pyz_LB_1Temp2,Pzx_LB_1Temp2,q1_LB_1Temp2,q2_LB_1Temp2,q3_LB_1Temp2);

           moments_LB_Send[p.numMoments*cellNumberLocal     ] += weightLB_P1 * rho_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  1] += weightLB_P1 * uX_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  2] += weightLB_P1 * uY_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  3] += weightLB_P1 * uZ_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  4] += weightLB_P1 * theta_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  5] += weightLB_P1 * Pxx_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  6] += weightLB_P1 * Pyy_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  7] += weightLB_P1 * Pzz_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  8] += weightLB_P1 * Pxy_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal +  9] += weightLB_P1 * Pyz_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal + 10] += weightLB_P1 * Pzx_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal + 11] += weightLB_P1 * q1_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal + 12] += weightLB_P1 * q2_LB_1Temp2[0];
           moments_LB_Send[p.numMoments*cellNumberLocal + 13] += weightLB_P1 * q3_LB_1Temp2[0];
         }
       }

       moments_LB_Send[p.numMoments*cellNumberLocal     ] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal +  1] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal +  2] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal +  3] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal +  4] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal +  5] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal +  6] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal +  7] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal +  8] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal +  9] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal + 10] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal + 11] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal + 12] /= ((double)(numSpaceAvgCells_LB));
       moments_LB_Send[p.numMoments*cellNumberLocal + 13] /= ((double)(numSpaceAvgCells_LB));

       double rho_LB_1TempSpaceAvgYplus   = 0.;
       // double uX_LB_1TempSpaceAvgYplus    = 0.;
       // double uY_LB_1TempSpaceAvgYplus    = 0.;
       // double uZ_LB_1TempSpaceAvgYplus    = 0.;
       // double theta_LB_1TempSpaceAvgYplus = 0.;
       // double Pxx_LB_1TempSpaceAvgYplus   = 0.;
       // double Pyy_LB_1TempSpaceAvgYplus   = 0.;
       // double Pzz_LB_1TempSpaceAvgYplus   = 0.;
       // double Pxy_LB_1TempSpaceAvgYplus   = 0.;
       // double Pyz_LB_1TempSpaceAvgYplus   = 0.;
       // double Pzx_LB_1TempSpaceAvgYplus   = 0.;
       // double q1_LB_1TempSpaceAvgYplus    = 0.;
       // double q2_LB_1TempSpaceAvgYplus    = 0.;
       // double q3_LB_1TempSpaceAvgYplus    = 0.;

       for(int i = 0; i < numSpaceAvgCells_LB_X; i++)
       {
         for(int j = 0; j < numSpaceAvgCells_LB_Y; j++)
         {
           int xIndex = cellY + gridLBM41.nB1 + j;
           int zIndex = cellX + gridLBM41.nB3 + i;
           int yIndex =         gridLBM41.nE2 - 2;//(0.5*p.nCoresZ - 1);

           copyFromNode(lbModel41,gridLBM41,VECT_LENGTH,xIndex,yIndex,zIndex);
           getGrad11Moment_heatAlpha_YplusFlux(lbModel41,VECT_LENGTH,rho_LB_1Temp2,uX_LB_1Temp2,uY_LB_1Temp2,uZ_LB_1Temp2,theta_LB_1Temp2,Pxx_LB_1Temp2,Pyy_LB_1Temp2,Pzz_LB_1Temp2,Pxy_LB_1Temp2,Pyz_LB_1Temp2,Pzx_LB_1Temp2,q1_LB_1Temp2,q2_LB_1Temp2,q3_LB_1Temp2);

           rho_LB_1TempSpaceAvgYplus    += weightLB_P2 * rho_LB_1Temp2[0];
           // uX_LB_1TempSpaceAvgYplus     += weightLB_P2 * uX_LB_1Temp2[0];
           // uY_LB_1TempSpaceAvgYplus     += weightLB_P2 * uY_LB_1Temp2[0];
           // uZ_LB_1TempSpaceAvgYplus     += weightLB_P2 * uZ_LB_1Temp2[0];
           // theta_LB_1TempSpaceAvgYplus  += weightLB_P2 * theta_LB_1Temp2[0];
           // Pxx_LB_1TempSpaceAvgYplus    += weightLB_P2 * Pxx_LB_1Temp2[0];
           // Pyy_LB_1TempSpaceAvgYplus    += weightLB_P2 * Pyy_LB_1Temp2[0];
           // Pzz_LB_1TempSpaceAvgYplus    += weightLB_P2 * Pzz_LB_1Temp2[0];
           // Pxy_LB_1TempSpaceAvgYplus    += weightLB_P2 * Pxy_LB_1Temp2[0];
           // Pyz_LB_1TempSpaceAvgYplus    += weightLB_P2 * Pyz_LB_1Temp2[0];
           // Pzx_LB_1TempSpaceAvgYplus    += weightLB_P2 * Pzx_LB_1Temp2[0];
           // q1_LB_1TempSpaceAvgYplus     += weightLB_P2 * q1_LB_1Temp2[0];
           // q2_LB_1TempSpaceAvgYplus     += weightLB_P2 * q2_LB_1Temp2[0];
           // q3_LB_1TempSpaceAvgYplus     += weightLB_P2 * q3_LB_1Temp2[0];

           copyFromNode(lbModel41,gridLBM41,VECT_LENGTH,xIndex,yIndex + 1,zIndex);
           getGrad11Moment_heatAlpha_YplusFlux(lbModel41,VECT_LENGTH,rho_LB_1Temp2,uX_LB_1Temp2,uY_LB_1Temp2,uZ_LB_1Temp2,theta_LB_1Temp2,Pxx_LB_1Temp2,Pyy_LB_1Temp2,Pzz_LB_1Temp2,Pxy_LB_1Temp2,Pyz_LB_1Temp2,Pzx_LB_1Temp2,q1_LB_1Temp2,q2_LB_1Temp2,q3_LB_1Temp2);

           rho_LB_1TempSpaceAvgYplus    += weightLB_P1 * rho_LB_1Temp2[0];
           // uX_LB_1TempSpaceAvgYplus     += weightLB_P1 * uX_LB_1Temp2[0];
           // uY_LB_1TempSpaceAvgYplus     += weightLB_P1 * uY_LB_1Temp2[0];
           // uZ_LB_1TempSpaceAvgYplus     += weightLB_P1 * uZ_LB_1Temp2[0];
           // theta_LB_1TempSpaceAvgYplus  += weightLB_P1 * theta_LB_1Temp2[0];
           // Pxx_LB_1TempSpaceAvgYplus    += weightLB_P1 * Pxx_LB_1Temp2[0];
           // Pyy_LB_1TempSpaceAvgYplus    += weightLB_P1 * Pyy_LB_1Temp2[0];
           // Pzz_LB_1TempSpaceAvgYplus    += weightLB_P1 * Pzz_LB_1Temp2[0];
           // Pxy_LB_1TempSpaceAvgYplus    += weightLB_P1 * Pxy_LB_1Temp2[0];
           // Pyz_LB_1TempSpaceAvgYplus    += weightLB_P1 * Pyz_LB_1Temp2[0];
           // Pzx_LB_1TempSpaceAvgYplus    += weightLB_P1 * Pzx_LB_1Temp2[0];
           // q1_LB_1TempSpaceAvgYplus     += weightLB_P1 * q1_LB_1Temp2[0];
           // q2_LB_1TempSpaceAvgYplus     += weightLB_P1 * q2_LB_1Temp2[0];
           // q3_LB_1TempSpaceAvgYplus     += weightLB_P1 * q3_LB_1Temp2[0];
         }
       }

       rho_LB_1TempSpaceAvgYplus    /= ((double)(numSpaceAvgCells_LB));
       // uX_LB_1TempSpaceAvgYplus     /= ((double)(numSpaceAvgCells_LB));
       // uY_LB_1TempSpaceAvgYplus     /= ((double)(numSpaceAvgCells_LB));
       // uZ_LB_1TempSpaceAvgYplus     /= ((double)(numSpaceAvgCells_LB));
       // theta_LB_1TempSpaceAvgYplus  /= ((double)(numSpaceAvgCells_LB));
       // Pxx_LB_1TempSpaceAvgYplus    /= ((double)(numSpaceAvgCells_LB));
       // Pyy_LB_1TempSpaceAvgYplus    /= ((double)(numSpaceAvgCells_LB));
       // Pzz_LB_1TempSpaceAvgYplus    /= ((double)(numSpaceAvgCells_LB));
       // Pxy_LB_1TempSpaceAvgYplus    /= ((double)(numSpaceAvgCells_LB));
       // Pyz_LB_1TempSpaceAvgYplus    /= ((double)(numSpaceAvgCells_LB));
       // Pzx_LB_1TempSpaceAvgYplus    /= ((double)(numSpaceAvgCells_LB));
       // q1_LB_1TempSpaceAvgYplus     /= ((double)(numSpaceAvgCells_LB));
       // q2_LB_1TempSpaceAvgYplus     /= ((double)(numSpaceAvgCells_LB));
       // q3_LB_1TempSpaceAvgYplus     /= ((double)(numSpaceAvgCells_LB));

       PANINI_REAL dsmcDensityFlux = ((double)p.dsmcParticleFlux[cellNumberLocal]/(p.volume*numSpaceAvgCells_DSMC_X*numSpaceAvgCells_DSMC_Y));
       // std::cout<<"dsmcDensityFlux: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<p.dsmcParticleFlux[cellNumberLocal]<<"\t"<<dsmcDensityFlux<<std::endl;

       // moments_LB_Send[p.numMoments*cellNumberLocal     ] += (rho_LB_1TempSpaceAvgYplus  *DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal +  1] += (uX_LB_1TempSpaceAvgYplus   *DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal +  2] += (uY_LB_1TempSpaceAvgYplus   *DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal +  3] += (uZ_LB_1TempSpaceAvgYplus   *DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal +  4] += (theta_LB_1TempSpaceAvgYplus*DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal +  5] += (Pxx_LB_1TempSpaceAvgYplus  *DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal +  6] += (Pyy_LB_1TempSpaceAvgYplus  *DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal +  7] += (Pzz_LB_1TempSpaceAvgYplus  *DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal +  8] += (Pxy_LB_1TempSpaceAvgYplus  *DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal +  9] += (Pyz_LB_1TempSpaceAvgYplus  *DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal + 10] += (Pzx_LB_1TempSpaceAvgYplus  *DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal + 11] += (q1_LB_1TempSpaceAvgYplus   *DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal + 12] += (q2_LB_1TempSpaceAvgYplus   *DSMC_CountInv);
       // moments_LB_Send[p.numMoments*cellNumberLocal + 13] += (q3_LB_1TempSpaceAvgYplus   *DSMC_CountInv);

       if((cellNumberLocal == 0) && (myRankY_LB == ((int)(0.5*p.nCoresY))) && (myRankX_LB == ((int)(0.5*p.nCoresX))))
       {
         std::cout<<"====================================================================================================================="<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 0<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal     ]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 1<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  1]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 2<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  2]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 3<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  3]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 4<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  4]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 5<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  5]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 6<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  6]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 7<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  7]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 8<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  8]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<< 9<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal +  9]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<<10<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal + 10]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<<11<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal + 11]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<<12<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal + 12]<<"\t"<<dsmcDensityFlux<<std::endl;
         std::cout<<"LB2DSMC_moments_upper: "<<"\t"<<count<<"\t"<<myRank<<"\t"<<myRankX<<"\t"<<myRankY<<"\t"<<cellNumberLocal<<"\t"<<13<<"\t"<<moments_LB_Send[p.numMoments*cellNumberLocal + 13]<<"\t"<<dsmcDensityFlux<<std::endl;
           std::cout<<"====================================================================================================================="<<std::endl;
       }
     }
   }
}

template<int DIM,typename real,typename dof>
void calcRecvLBmomentsUpperBufferZone(hardSphere<DIM,real,dof>& p, lbmRD3Q41<real> &lbModel41, gridBCC3D<4, 11, real> &gridLBM41, real* moments_LB_Recv, real* momentsSave,
		real c, real f2g_factor, real weightLB_P1, real weightLB_P2, int LBbyDSMCpointsX, int LBbyDSMCpointsY, int numSpaceAvgCells_LB_X, int numSpaceAvgCells_LB_Y, int numSpaceAvgCells_LB, 
		int numSpaceAvgCells_DSMC_X, int numSpaceAvgCells_DSMC_Y, int numSpaceAvgCells_DSMC, int myRank, int myRankX, int myRankY, int myRankZ, int myRankX_LB, int myRankY_LB, 
		int myRankZ_LB, int count, int upperCoupleDSMC)
{
   for(int cellY = 0; cellY < (p.YcellSizePerCore/LBbyDSMCpointsY); cellY = cellY + numSpaceAvgCells_LB_Y)
   {
     for(int cellX = 0; cellX < (p.XcellSizePerCore/LBbyDSMCpointsX); cellX = cellX + numSpaceAvgCells_LB_X)
     {
       int cellNumberLocal = ((int)(p.XcellSizePerCore/(numSpaceAvgCells_DSMC_X)))*((int)(cellY/numSpaceAvgCells_LB_Y)) + (int)(cellX/numSpaceAvgCells_LB_X);

       PANINI_REAL rho_DSMC,uX_DSMC,uY_DSMC,uZ_DSMC,sigmaXX_DSMC,sigmaYY_DSMC,sigmaZZ_DSMC,sigmaXY_DSMC,sigmaYZ_DSMC,sigmaZX_DSMC,p_DSMC,qX_DSMC,qY_DSMC,qZ_DSMC;

       double rho_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal     ];
       double uX_LB_1TempSpaceAvg    = moments_LB_Recv[p.numMoments*cellNumberLocal +  1];
       double uY_LB_1TempSpaceAvg    = moments_LB_Recv[p.numMoments*cellNumberLocal +  2];
       double uZ_LB_1TempSpaceAvg    = moments_LB_Recv[p.numMoments*cellNumberLocal +  3];
       double theta_LB_1TempSpaceAvg = moments_LB_Recv[p.numMoments*cellNumberLocal +  4];
       double Pxx_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal +  5];
       double Pyy_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal +  6];
       double Pzz_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal +  7];
       double Pxy_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal +  8];
       double Pyz_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal +  9];
       double Pzx_LB_1TempSpaceAvg   = moments_LB_Recv[p.numMoments*cellNumberLocal + 10];
       double q1_LB_1TempSpaceAvg    = moments_LB_Recv[p.numMoments*cellNumberLocal + 11];
       double q2_LB_1TempSpaceAvg    = moments_LB_Recv[p.numMoments*cellNumberLocal + 12];
       double q3_LB_1TempSpaceAvg    = moments_LB_Recv[p.numMoments*cellNumberLocal + 13];

       rho_DSMC       = rho_LB_1TempSpaceAvg;

       p_DSMC         = rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg*c;

       uX_DSMC        = uZ_LB_1TempSpaceAvg*c;
       uY_DSMC        = uX_LB_1TempSpaceAvg*c;
       uZ_DSMC        = uY_LB_1TempSpaceAvg*c;

       PANINI_REAL u2 = uX_LB_1TempSpaceAvg*uX_LB_1TempSpaceAvg + uY_LB_1TempSpaceAvg*uY_LB_1TempSpaceAvg + uZ_LB_1TempSpaceAvg*uZ_LB_1TempSpaceAvg;

       sigmaXX_DSMC   = (Pzz_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*uZ_LB_1TempSpaceAvg*uZ_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg)*f2g_factor;
       sigmaYY_DSMC   = (Pxx_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*uX_LB_1TempSpaceAvg*uX_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg)*f2g_factor;
       sigmaYZ_DSMC   = (Pxy_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*uX_LB_1TempSpaceAvg*uY_LB_1TempSpaceAvg)*f2g_factor;
       sigmaZZ_DSMC   = (Pyy_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*uY_LB_1TempSpaceAvg*uY_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg)*f2g_factor;
       sigmaXY_DSMC   = (Pzx_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*uZ_LB_1TempSpaceAvg*uX_LB_1TempSpaceAvg)*f2g_factor;
       sigmaZX_DSMC   = (Pyz_LB_1TempSpaceAvg - rho_LB_1TempSpaceAvg*uY_LB_1TempSpaceAvg*uZ_LB_1TempSpaceAvg)*f2g_factor;

       qX_DSMC = ((q3_LB_1TempSpaceAvg) - ((rho_LB_1TempSpaceAvg*u2*uZ_LB_1TempSpaceAvg) + 5.0*rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg*uZ_LB_1TempSpaceAvg + 2.0*(uX_LB_1TempSpaceAvg*sigmaXY_DSMC + uY_LB_1TempSpaceAvg*sigmaZX_DSMC + uZ_LB_1TempSpaceAvg*sigmaXX_DSMC)/f2g_factor))*(c*c*c)*f2g_factor;
       qY_DSMC = ((q1_LB_1TempSpaceAvg) - ((rho_LB_1TempSpaceAvg*u2*uX_LB_1TempSpaceAvg) + 5.0*rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg*uX_LB_1TempSpaceAvg + 2.0*(uX_LB_1TempSpaceAvg*sigmaYY_DSMC + uY_LB_1TempSpaceAvg*sigmaYZ_DSMC + uZ_LB_1TempSpaceAvg*sigmaXY_DSMC)/f2g_factor))*(c*c*c)*f2g_factor;
       qZ_DSMC = ((q2_LB_1TempSpaceAvg) - ((rho_LB_1TempSpaceAvg*u2*uY_LB_1TempSpaceAvg) + 5.0*rho_LB_1TempSpaceAvg*theta_LB_1TempSpaceAvg*uY_LB_1TempSpaceAvg + 2.0*(uX_LB_1TempSpaceAvg*sigmaYZ_DSMC + uY_LB_1TempSpaceAvg*sigmaZZ_DSMC + uZ_LB_1TempSpaceAvg*sigmaZX_DSMC)/f2g_factor))*(c*c*c)*f2g_factor;

       sigmaXX_DSMC   *= (c*c);
       sigmaXY_DSMC   *= (c*c);
       sigmaZX_DSMC   *= (c*c);
       sigmaYY_DSMC   *= (c*c);
       sigmaYZ_DSMC   *= (c*c);
       sigmaZZ_DSMC   *= (c*c);


       if((cellNumberLocal == 0) && (myRankY == ((int)(0.5*p.nCoresY))) && (myRankX == ((int)(0.5*p.nCoresX))))
       {
         std::cout<<"====================================================================================================================="<<std::endl;
         cout<<"rho_DSMC_upper@: "<<count<<" "<<myRank<<":"<<rho_DSMC<<std::endl;
         cout<<"uX_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<uX_DSMC<<std::endl;
         cout<<"uY_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<uY_DSMC<<std::endl;
         cout<<"uZ_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<uZ_DSMC<<std::endl;
         cout<<"t0_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<theta_LB_1TempSpaceAvg<<std::endl;
         cout<<"sXX_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<sigmaXX_DSMC<<std::endl;
         cout<<"SXY_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<sigmaXY_DSMC<<std::endl;
         cout<<"sYZ_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<sigmaYZ_DSMC<<std::endl;
         cout<<"sYY_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<sigmaYY_DSMC<<std::endl;
         cout<<"sZZ_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<sigmaZZ_DSMC<<std::endl;
         cout<<"sZX_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<sigmaZX_DSMC<<std::endl;
         cout<<"qX_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<qX_DSMC<<std::endl;
         cout<<"qY_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<qY_DSMC<<std::endl;
         cout<<"qZ_DSMC_upper@ : "<<count<<" "<<myRank<<":"<<qZ_DSMC<<std::endl;
           std::cout<<"====================================================================================================================="<<std::endl;
       }

       momentsSave[p.numMoments*cellNumberLocal     ] = rho_DSMC    ;
       momentsSave[p.numMoments*cellNumberLocal +  1] = p_DSMC      ;
       momentsSave[p.numMoments*cellNumberLocal +  2] = uX_DSMC     ;
       momentsSave[p.numMoments*cellNumberLocal +  3] = uY_DSMC     ;
       momentsSave[p.numMoments*cellNumberLocal +  4] = uZ_DSMC     ;
       momentsSave[p.numMoments*cellNumberLocal +  5] = sigmaXX_DSMC;
       momentsSave[p.numMoments*cellNumberLocal +  6] = sigmaXY_DSMC;
       momentsSave[p.numMoments*cellNumberLocal +  7] = sigmaZX_DSMC;
       momentsSave[p.numMoments*cellNumberLocal +  8] = sigmaYY_DSMC;
       momentsSave[p.numMoments*cellNumberLocal +  9] = sigmaYZ_DSMC;
       momentsSave[p.numMoments*cellNumberLocal + 10] = sigmaZZ_DSMC;
       momentsSave[p.numMoments*cellNumberLocal + 11] = qX_DSMC     ;
       momentsSave[p.numMoments*cellNumberLocal + 12] = qY_DSMC     ;
       momentsSave[p.numMoments*cellNumberLocal + 13] = qZ_DSMC     ;

       // p.totalMomentsWithFlux(&momentsSave[p.numMoments*cellNumberLocal],cellNumberLocal,myRank,count);
       // p.generateDistributionFromMoments(&momentsSave[p.numMoments*cellNumberLocal],cellX/numSpaceAvgCells_LB_X,cellY/numSpaceAvgCells_LB_Y,myRank,count);

     }
   }
}

template<int DIM,typename real,typename dof>
void calcSendDSMCmomentsLowerBufferZone(hardSphere<DIM,real,dof>& p, real* moments_DSMC_Send, int DSMC2LBCouplingCellZ, int LBbyDSMCpointsX, int LBbyDSMCpointsY, int numSpaceAvgCells_LB_X, 
		int numSpaceAvgCells_LB_Y, int numSpaceAvgCells_LB, int numSpaceAvgCells_DSMC_X, int numSpaceAvgCells_DSMC_Y, int numSpaceAvgCells_DSMC, int myRank, int count)
{
   	for(int cellX = 0; cellX < (p.XcellSizePerCore/LBbyDSMCpointsX); cellX++)
	{
	  for(int cellY = 0; cellY < (p.YcellSizePerCore/LBbyDSMCpointsY); cellY++)
	  {
	      int cellNumberLocal           = ((int)(p.XcellSizePerCore/numSpaceAvgCells_DSMC_X))*((int)(cellY/numSpaceAvgCells_LB_Y)) + ((int)(cellX/numSpaceAvgCells_LB_X));
	      int cellNumberLocalDSMCcouple = DSMC2LBCouplingCellZ*(p.XcellSizePerCore/numSpaceAvgCells_DSMC_X)*(p.YcellSizePerCore/numSpaceAvgCells_DSMC_Y) + cellNumberLocal;

	      moments_DSMC_Send[p.numMoments*cellNumberLocal     ] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple     ];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal     ]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal +  1] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  1];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  1]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal +  2] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  2];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  2]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal +  3] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  3];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  3]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal +  4] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  4];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  4]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal +  5] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  5];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  5]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal +  6] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  6];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  6]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal +  7] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  7];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  7]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal +  8] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  8];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  8]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal +  9] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  9];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  9]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal + 10] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple + 10];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal + 10]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal + 11] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple + 11];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal + 11]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal + 12] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple + 12];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal + 12]*LB_CountInv;
	      moments_DSMC_Send[p.numMoments*cellNumberLocal + 13] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple + 13];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal + 13]*LB_CountInv;

	      if((myRank == 0) && (cellNumberLocal == 0))
	      {
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<0 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal     ]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<1 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  1]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<2 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  2]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<3 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  3]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<4 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  4]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<5 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  5]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<6 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  6]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<7 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  7]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<8 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  8]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<9 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  9]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<10<<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal + 10]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<11<<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal + 11]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<12<<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal + 12]<<std::endl;;
		std::cout<<"momentsLBRecvLower: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<13<<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal + 13]<<std::endl;;
	      }
	  }
	}
}

template<int DIM,typename real,typename dof>
void calcRecvDSMCmomentsLowerBufferZone(hardSphere<DIM,real,dof>& p, lbmRD3Q41<real> &lbModel41, gridBCC3D<4, 11, real> &gridLBM41, real* moments_DSMC_Recv, real c, real f2g_factor, real dt, 
	real tau, int DSMC2LBCouplingCellZ, int LBbyDSMCpointsX, int LBbyDSMCpointsY, int numSpaceAvgCells_LB_X, int numSpaceAvgCells_LB_Y, int numSpaceAvgCells_LB, int numSpaceAvgCells_DSMC_X, 
	int numSpaceAvgCells_DSMC_Y, int numSpaceAvgCells_DSMC, int myRank, int myRankX_LB, int myRankY_LB, int myRankZ_LB, int step)
{
      int VECT_LENGTH = 4;
      double rho_LB_1Temp2[4],uX_LB_1Temp2[4],uY_LB_1Temp2[4],uZ_LB_1Temp2[4],theta_LB_1Temp2[4],Pxx_LB_1Temp2[4],Pyy_LB_1Temp2[4],Pzz_LB_1Temp2[4],Pxy_LB_1Temp2[4],Pyz_LB_1Temp2[4],Pzx_LB_1Temp2[4];
      double q1_LB_1Temp2[4],q2_LB_1Temp2[4],q3_LB_1Temp2[4];
      for(int cellX = 0; cellX < (p.XcellSizePerCore/LBbyDSMCpointsX); cellX++)
      {
          for(int cellY = 0; cellY < (p.YcellSizePerCore/LBbyDSMCpointsY); cellY++)
          {
              int cellNumberLocal           = ((int)(p.XcellSizePerCore/numSpaceAvgCells_DSMC_X))*((int)(cellY/numSpaceAvgCells_LB_Y)) + ((int)(cellX/numSpaceAvgCells_LB_X));
              int cellNumberLocalDSMCcouple = DSMC2LBCouplingCellZ*(p.XcellSizePerCore/numSpaceAvgCells_DSMC_X)*(p.YcellSizePerCore/numSpaceAvgCells_DSMC_Y) + cellNumberLocal;

              PANINI_REAL rho_DSMC,uX_DSMC,uY_DSMC,uZ_DSMC,sigmaXX_DSMC,sigmaYY_DSMC,sigmaZZ_DSMC,sigmaXY_DSMC,sigmaYZ_DSMC,sigmaZX_DSMC,p_DSMC,qX_DSMC,qY_DSMC,qZ_DSMC;

              rho_DSMC      = moments_DSMC_Recv[p.numMoments*cellNumberLocal     ];
              p_DSMC        = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  1];
              uX_DSMC       = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  2];
              uY_DSMC       = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  3];
              uZ_DSMC       = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  4];
              sigmaXX_DSMC  = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  5];
              sigmaXY_DSMC  = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  6];
              sigmaZX_DSMC  = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  7];
              sigmaYY_DSMC  = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  8];
              sigmaYZ_DSMC  = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  9];
              sigmaZZ_DSMC  = moments_DSMC_Recv[p.numMoments*cellNumberLocal + 10];
              qX_DSMC       = moments_DSMC_Recv[p.numMoments*cellNumberLocal + 11];
              qY_DSMC       = moments_DSMC_Recv[p.numMoments*cellNumberLocal + 12];
              qZ_DSMC       = moments_DSMC_Recv[p.numMoments*cellNumberLocal + 13];


              PANINI_REAL grad[lbModel41.dvN];
              PANINI_REAL fEq[lbModel41.dvN];

              // Difference in notations pf LB and DSMC
              // LB(x,y,z) = DSMC(y,z,x)

              PANINI_REAL rho_LB   = rho_DSMC;

              PANINI_REAL uX_LB   = uY_DSMC/c;
              PANINI_REAL uY_LB   = uZ_DSMC/c;
              PANINI_REAL uZ_LB   = uX_DSMC/c;
              PANINI_REAL u2      = uX_LB*uX_LB + uY_LB*uY_LB + uZ_LB*uZ_LB;

              PANINI_REAL theta_LB = (p_DSMC/rho_DSMC)/(c*c);

              PANINI_REAL q1_LB_g = (qY_DSMC + (rho_LB*u2*uX_LB) + 2.0*(uX_LB*sigmaYY_DSMC + uY_LB*sigmaYZ_DSMC + uZ_LB*sigmaXY_DSMC))/(c*c*c*f2g_factor) + ((3.0 - (dt/tau))*rho_LB*theta_LB*uX_LB)/(c*c*c);
              PANINI_REAL q2_LB_g = (qZ_DSMC + (rho_LB*u2*uY_LB) + 2.0*(uX_LB*sigmaYZ_DSMC + uY_LB*sigmaZZ_DSMC + uZ_LB*sigmaZX_DSMC))/(c*c*c*f2g_factor) + ((3.0 - (dt/tau))*rho_LB*theta_LB*uY_LB)/(c*c*c);
              PANINI_REAL q3_LB_g = (qX_DSMC + (rho_LB*u2*uZ_LB) + 2.0*(uX_LB*sigmaXY_DSMC + uY_LB*sigmaZX_DSMC + uZ_LB*sigmaXX_DSMC))/(c*c*c*f2g_factor) + ((3.0 - (dt/tau))*rho_LB*theta_LB*uZ_LB)/(c*c*c);

              PANINI_REAL Pxx_LB_g  = ((sigmaYY_DSMC/f2g_factor + rho_DSMC*uY_DSMC*uY_DSMC - (dt/(2.*tau))*(rho_DSMC*theta_LB))/(c*c));
              PANINI_REAL Pyy_LB_g  = ((sigmaZZ_DSMC/f2g_factor + rho_DSMC*uZ_DSMC*uZ_DSMC - (dt/(2.*tau))*(rho_DSMC*theta_LB))/(c*c));
              PANINI_REAL Pzz_LB_g  = ((sigmaXX_DSMC/f2g_factor + rho_DSMC*uX_DSMC*uX_DSMC - (dt/(2.*tau))*(rho_DSMC*theta_LB))/(c*c));
              PANINI_REAL Pxy_LB_g  = ((sigmaYZ_DSMC/f2g_factor + rho_DSMC*uY_DSMC*uZ_DSMC)/(c*c));
              PANINI_REAL Pyz_LB_g  = ((sigmaZX_DSMC/f2g_factor + rho_DSMC*uZ_DSMC*uX_DSMC)/(c*c));
              PANINI_REAL Pzx_LB_g  = ((sigmaXY_DSMC/f2g_factor + rho_DSMC*uX_DSMC*uY_DSMC)/(c*c));

              int xIndex = cellY + gridLBM41.nB1;
              int zIndex = cellX + gridLBM41.nB3;
              int yIndex =         gridLBM41.nB2 + 2;// + (0.25*p.nCoresZ - 1);

              copyFromNode(lbModel41,gridLBM41,VECT_LENGTH,xIndex,yIndex,zIndex);
              getGrad11Moment_heatAlpha_YplusFlux(lbModel41,VECT_LENGTH,rho_LB_1Temp2,uX_LB_1Temp2,uY_LB_1Temp2,uZ_LB_1Temp2,theta_LB_1Temp2,Pxx_LB_1Temp2,Pyy_LB_1Temp2,Pzz_LB_1Temp2,Pxy_LB_1Temp2,Pyz_LB_1Temp2,Pzx_LB_1Temp2,q1_LB_1Temp2,q2_LB_1Temp2,q3_LB_1Temp2);

              if((cellNumberLocal == 0) && (myRankY_LB == ((int)(0.5*p.nCoresY))) && (myRankX_LB == ((int)(0.5*p.nCoresX))))
              {
                  cout<<"rho_plusFlux_Lower@: "<<step<<" "<<myRank<<":"<<rho_LB<<"\t"<<rho_LB_1Temp2[0]<<std::endl;
                  cout<<"uX_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<uX_LB<<"\t"<<uX_LB_1Temp2[0]<<std::endl;
                  cout<<"uY_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<uY_LB<<"\t"<<uY_LB_1Temp2[0]<<std::endl;
                  cout<<"uZ_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<uZ_LB<<"\t"<<uZ_LB_1Temp2[0]<<std::endl;
                  cout<<"t0_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<theta_LB<<"\t"<<theta_LB_1Temp2[0]<<std::endl;
                  cout<<"sXX_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<Pxx_LB_g<<"\t"<<Pxx_LB_1Temp2[0]<<std::endl;
                  cout<<"SXY_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<Pxy_LB_g<<"\t"<<Pxy_LB_1Temp2[0]<<std::endl;
                  cout<<"sYZ_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<Pyz_LB_g<<"\t"<<Pyz_LB_1Temp2[0]<<std::endl;
                  cout<<"sYY_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<Pyy_LB_g<<"\t"<<Pyy_LB_1Temp2[0]<<std::endl;
                  cout<<"sZZ_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<Pzz_LB_g<<"\t"<<Pzz_LB_1Temp2[0]<<std::endl;
                  cout<<"sZX_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<Pzx_LB_g<<"\t"<<Pzx_LB_1Temp2[0]<<std::endl;
                  cout<<"qX_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<q1_LB_g<<"\t"<<q1_LB_1Temp2[0]<<std::endl;
                  cout<<"qY_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<q2_LB_g<<"\t"<<q2_LB_1Temp2[0]<<std::endl;
                  cout<<"qZ_plusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<q3_LB_g<<"\t"<<q3_LB_1Temp2[0]<<std::endl;
              }

              copyFromNode(lbModel41,gridLBM41,VECT_LENGTH,xIndex,yIndex,zIndex);
              getGrad11Moment_heatAlpha_YminusFlux(lbModel41,VECT_LENGTH,rho_LB_1Temp2,uX_LB_1Temp2,uY_LB_1Temp2,uZ_LB_1Temp2,theta_LB_1Temp2,Pxx_LB_1Temp2,Pyy_LB_1Temp2,Pzz_LB_1Temp2,Pxy_LB_1Temp2,Pyz_LB_1Temp2,Pzx_LB_1Temp2,q1_LB_1Temp2,q2_LB_1Temp2,q3_LB_1Temp2);

              if((cellNumberLocal == 0) && (myRankY_LB == ((int)(0.5*p.nCoresY))) && (myRankX_LB == ((int)(0.5*p.nCoresX))))
              {
                  cout<<"rho_minusFlux_Lower@: "<<step<<" "<<myRank<<":"<<rho_LB<<"\t"<<rho_LB_1Temp2[0]<<std::endl;
                  cout<<"uX_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<uX_LB<<"\t"<<uX_LB_1Temp2[0]<<std::endl;
                  cout<<"uY_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<uY_LB<<"\t"<<uY_LB_1Temp2[0]<<std::endl;
                  cout<<"uZ_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<uZ_LB<<"\t"<<uZ_LB_1Temp2[0]<<std::endl;
                  cout<<"t0_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<theta_LB<<"\t"<<theta_LB_1Temp2[0]<<std::endl;
                  cout<<"sXX_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<Pxx_LB_g<<"\t"<<Pxx_LB_1Temp2[0]<<std::endl;
                  cout<<"SXY_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<Pxy_LB_g<<"\t"<<Pxy_LB_1Temp2[0]<<std::endl;
                  cout<<"sYZ_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<Pyz_LB_g<<"\t"<<Pyz_LB_1Temp2[0]<<std::endl;
                  cout<<"sYY_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<Pyy_LB_g<<"\t"<<Pyy_LB_1Temp2[0]<<std::endl;
                  cout<<"sZZ_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<Pzz_LB_g<<"\t"<<Pzz_LB_1Temp2[0]<<std::endl;
                  cout<<"sZX_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<Pzx_LB_g<<"\t"<<Pzx_LB_1Temp2[0]<<std::endl;
                  cout<<"qX_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<q1_LB_g<<"\t"<<q1_LB_1Temp2[0]<<std::endl;
                  cout<<"qY_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<q2_LB_g<<"\t"<<q2_LB_1Temp2[0]<<std::endl;
                  cout<<"qZ_minusFlux_Lower@ : "<<step<<" "<<myRank<<":"<<q3_LB_g<<"\t"<<q3_LB_1Temp2[0]<<std::endl;
              }

              if((cellNumberLocal == 0) && (myRankY_LB == ((int)(0.5*p.nCoresY))) && (myRankX_LB == ((int)(0.5*p.nCoresX))))
              {
                  cout<<"rho_LB_Lower@: "<<step<<" "<<myRank<<":"<<rho_LB<<"\t"<<rho_LB_1Temp2[0]<<std::endl;
                  cout<<"uX_LB_Lower@ : "<<step<<" "<<myRank<<":"<<uX_LB<<"\t"<<uX_LB_1Temp2[0]<<std::endl;
                  cout<<"uY_LB_Lower@ : "<<step<<" "<<myRank<<":"<<uY_LB<<"\t"<<uY_LB_1Temp2[0]<<std::endl;
                  cout<<"uZ_LB_Lower@ : "<<step<<" "<<myRank<<":"<<uZ_LB<<"\t"<<uZ_LB_1Temp2[0]<<std::endl;
                  cout<<"t0_LB_Lower@ : "<<step<<" "<<myRank<<":"<<theta_LB<<"\t"<<theta_LB_1Temp2[0]<<std::endl;
                  cout<<"sXX_LB_Lower@ : "<<step<<" "<<myRank<<":"<<Pxx_LB_g<<"\t"<<Pxx_LB_1Temp2[0]<<std::endl;
                  cout<<"SXY_LB_Lower@ : "<<step<<" "<<myRank<<":"<<Pxy_LB_g<<"\t"<<Pxy_LB_1Temp2[0]<<std::endl;
                  cout<<"sYZ_LB_Lower@ : "<<step<<" "<<myRank<<":"<<Pyz_LB_g<<"\t"<<Pyz_LB_1Temp2[0]<<std::endl;
                  cout<<"sYY_LB_Lower@ : "<<step<<" "<<myRank<<":"<<Pyy_LB_g<<"\t"<<Pyy_LB_1Temp2[0]<<std::endl;
                  cout<<"sZZ_LB_Lower@ : "<<step<<" "<<myRank<<":"<<Pzz_LB_g<<"\t"<<Pzz_LB_1Temp2[0]<<std::endl;
                  cout<<"sZX_LB_Lower@ : "<<step<<" "<<myRank<<":"<<Pzx_LB_g<<"\t"<<Pzx_LB_1Temp2[0]<<std::endl;
                  cout<<"qX_LB_Lower@ : "<<step<<" "<<myRank<<":"<<q1_LB_g<<"\t"<<q1_LB_1Temp2[0]<<std::endl;
                  cout<<"qY_LB_Lower@ : "<<step<<" "<<myRank<<":"<<q2_LB_g<<"\t"<<q2_LB_1Temp2[0]<<std::endl;
                  cout<<"qZ_LB_Lower@ : "<<step<<" "<<myRank<<":"<<q3_LB_g<<"\t"<<q3_LB_1Temp2[0]<<std::endl;
              }

              getGradThermalSinglePoint_heatAlpha (lbModel41,rho_LB,uX_LB,uY_LB,uZ_LB,theta_LB,Pxx_LB_g,Pyy_LB_g,Pzz_LB_g,Pxy_LB_g,Pyz_LB_g,Pzx_LB_g,q1_LB_g,q2_LB_g,q3_LB_g,grad);

              // for(int i = 0; i < numSpaceAvgCells_LB_X; i++)
              {
                  // for(int j = 0; j < numSpaceAvgCells_LB_Y; j++)
                  {

                      copyFromArrayToNode(lbModel41,grad, gridLBM41,xIndex,yIndex,zIndex);
                      copyFromArrayToCell(lbModel41,grad, gridLBM41,xIndex,yIndex,zIndex);

                      copyFromArrayToNode(lbModel41,grad, gridLBM41,xIndex,yIndex - 1,zIndex);
                      copyFromArrayToCell(lbModel41,grad, gridLBM41,xIndex,yIndex - 1,zIndex);


                      copyFromArrayToNode(lbModel41,grad, gridLBM41,xIndex,yIndex - 2,zIndex);
                      copyFromArrayToCell(lbModel41,grad, gridLBM41,xIndex,yIndex - 2,zIndex);

                      copyFromArrayToNode(lbModel41,grad, gridLBM41,xIndex,yIndex - 3,zIndex);
                      copyFromArrayToCell(lbModel41,grad, gridLBM41,xIndex,yIndex - 3,zIndex);


                      copyFromArrayToNode(lbModel41,grad, gridLBM41,xIndex,yIndex - 4,zIndex);
                      copyFromArrayToCell(lbModel41,grad, gridLBM41,xIndex,yIndex - 4,zIndex);
                  }
              }
          }
      }
}

template<int DIM,typename real,typename dof>
void calcSendDSMCmomentsUpperBufferZone(hardSphere<DIM,real,dof>& p, real* moments_DSMC_Send, int DSMC2LBCouplingCellZ, int LBbyDSMCpointsX, int LBbyDSMCpointsY, int numSpaceAvgCells_LB_X, 
		int numSpaceAvgCells_LB_Y, int numSpaceAvgCells_LB, int numSpaceAvgCells_DSMC_X, int numSpaceAvgCells_DSMC_Y, int numSpaceAvgCells_DSMC, int myRank, int count)
{
      for(int cellX = 0; cellX < (p.XcellSizePerCore/LBbyDSMCpointsX); cellX++)
      {
          for(int cellY = 0; cellY < (p.YcellSizePerCore/LBbyDSMCpointsY); cellY++)
          {
              int cellNumberLocal           = ((int)(p.XcellSizePerCore/numSpaceAvgCells_DSMC_X))*((int)(cellY/numSpaceAvgCells_LB_Y)) + ((int)(cellX/numSpaceAvgCells_LB_X));
              int cellNumberLocalDSMCcouple = DSMC2LBCouplingCellZ*(p.XcellSizePerCore/numSpaceAvgCells_DSMC_X)*(p.YcellSizePerCore/numSpaceAvgCells_DSMC_Y) + cellNumberLocal;

              moments_DSMC_Send[p.numMoments*cellNumberLocal     ] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple     ];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal     ]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal +  1] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  1];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  1]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal +  2] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  2];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  2]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal +  3] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  3];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  3]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal +  4] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  4];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  4]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal +  5] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  5];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  5]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal +  6] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  6];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  6]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal +  7] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  7];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  7]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal +  8] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  8];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  8]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal +  9] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple +  9];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal +  9]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal + 10] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple + 10];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal + 10]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal + 11] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple + 11];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal + 11]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal + 12] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple + 12];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal + 12]*LB_CountInv;
              moments_DSMC_Send[p.numMoments*cellNumberLocal + 13] = p.moments_Average[p.numMoments*cellNumberLocalDSMCcouple + 13];// + p.moments_Average_Flux[p.numMoments*cellNumberLocal + 13]*LB_CountInv;

              if((myRank == 0) && (cellNumberLocal == 0))
              {
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<0 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal     ]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<1 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  1]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<2 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  2]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<3 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  3]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<4 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  4]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<5 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  5]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<6 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  6]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<7 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  7]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<8 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  8]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<9 <<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal +  9]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<10<<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal + 10]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<11<<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal + 11]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<12<<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal + 12]<<std::endl;;
                std::cout<<"momentsLBRecvUpper: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<13<<"\t"<<p.moments_Average_Flux[p.numMoments*cellNumberLocal + 13]<<std::endl;;
              }
          }
      }
	
}

template<int DIM,typename real,typename dof>
void calcRecvDSMCmomentsUpperBufferZone(hardSphere<DIM,real,dof>& p, lbmRD3Q41<real> &lbModel41, gridBCC3D<4, 11, real> &gridLBM41, real* moments_DSMC_Recv, real c, real f2g_factor, real dt, 
	real tau, int DSMC2LBCouplingCellZ, int LBbyDSMCpointsX, int LBbyDSMCpointsY, int numSpaceAvgCells_LB_X, int numSpaceAvgCells_LB_Y, int numSpaceAvgCells_LB, int numSpaceAvgCells_DSMC_X, 
	int numSpaceAvgCells_DSMC_Y, int numSpaceAvgCells_DSMC, int myRank, int myRankX_LB, int myRankY_LB, int myRankZ_LB, int step)
{
      for(int cellX = 0; cellX < (p.XcellSizePerCore/LBbyDSMCpointsX); cellX++)
      {
          for(int cellY = 0; cellY < (p.YcellSizePerCore/LBbyDSMCpointsY); cellY++)
          {
              int cellNumberLocal           = ((int)(p.XcellSizePerCore/numSpaceAvgCells_DSMC_X))*((int)(cellY/numSpaceAvgCells_LB_Y)) + ((int)(cellX/numSpaceAvgCells_LB_X));
              int cellNumberLocalDSMCcouple = DSMC2LBCouplingCellZ*(p.XcellSizePerCore/numSpaceAvgCells_DSMC_X)*(p.YcellSizePerCore/numSpaceAvgCells_DSMC_Y) + cellNumberLocal;

              PANINI_REAL rho_DSMC,uX_DSMC,uY_DSMC,uZ_DSMC,sigmaXX_DSMC,sigmaYY_DSMC,sigmaZZ_DSMC,sigmaXY_DSMC,sigmaYZ_DSMC,sigmaZX_DSMC,p_DSMC,qX_DSMC,qY_DSMC,qZ_DSMC;

              rho_DSMC      = moments_DSMC_Recv[p.numMoments*cellNumberLocal     ];
              p_DSMC        = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  1];
              uX_DSMC       = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  2];
              uY_DSMC       = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  3];
              uZ_DSMC       = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  4];
              sigmaXX_DSMC  = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  5];
              sigmaXY_DSMC  = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  6];
              sigmaZX_DSMC  = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  7];
              sigmaYY_DSMC  = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  8];
              sigmaYZ_DSMC  = moments_DSMC_Recv[p.numMoments*cellNumberLocal +  9];
              sigmaZZ_DSMC  = moments_DSMC_Recv[p.numMoments*cellNumberLocal + 10];
              qX_DSMC       = moments_DSMC_Recv[p.numMoments*cellNumberLocal + 11];
              qY_DSMC       = moments_DSMC_Recv[p.numMoments*cellNumberLocal + 12];
              qZ_DSMC       = moments_DSMC_Recv[p.numMoments*cellNumberLocal + 13];

              PANINI_REAL grad[lbModel41.dvN];
              PANINI_REAL fEq[lbModel41.dvN];

              // Difference in notations pf LB and DSMC
              // LB(x,y,z) = DSMC(y,z,x)

              PANINI_REAL rho_LB   = rho_DSMC;

              PANINI_REAL uX_LB   = uY_DSMC/c;
              PANINI_REAL uY_LB   = uZ_DSMC/c;
              PANINI_REAL uZ_LB   = uX_DSMC/c;
              PANINI_REAL u2      = uX_LB*uX_LB + uY_LB*uY_LB + uZ_LB*uZ_LB;

              PANINI_REAL theta_LB = (p_DSMC/rho_DSMC)/(c*c);

              PANINI_REAL q1_LB_g = (qY_DSMC + (rho_LB*u2*uX_LB) + 2.0*(uX_LB*sigmaYY_DSMC + uY_LB*sigmaYZ_DSMC + uZ_LB*sigmaXY_DSMC))/(c*c*c*f2g_factor) + ((3.0 - (dt/tau))*rho_LB*theta_LB*uX_LB)/(c*c*c);
              PANINI_REAL q2_LB_g = (qZ_DSMC + (rho_LB*u2*uY_LB) + 2.0*(uX_LB*sigmaYZ_DSMC + uY_LB*sigmaZZ_DSMC + uZ_LB*sigmaZX_DSMC))/(c*c*c*f2g_factor) + ((3.0 - (dt/tau))*rho_LB*theta_LB*uY_LB)/(c*c*c);
              PANINI_REAL q3_LB_g = (qX_DSMC + (rho_LB*u2*uZ_LB) + 2.0*(uX_LB*sigmaXY_DSMC + uY_LB*sigmaZX_DSMC + uZ_LB*sigmaXX_DSMC))/(c*c*c*f2g_factor) + ((3.0 - (dt/tau))*rho_LB*theta_LB*uZ_LB)/(c*c*c);

              PANINI_REAL Pxx_LB_g  = ((sigmaYY_DSMC/f2g_factor + rho_DSMC*uY_DSMC*uY_DSMC - (dt/(2.*tau))*(rho_DSMC*theta_LB))/(c*c));
              PANINI_REAL Pyy_LB_g  = ((sigmaZZ_DSMC/f2g_factor + rho_DSMC*uZ_DSMC*uZ_DSMC - (dt/(2.*tau))*(rho_DSMC*theta_LB))/(c*c));
              PANINI_REAL Pzz_LB_g  = ((sigmaXX_DSMC/f2g_factor + rho_DSMC*uX_DSMC*uX_DSMC - (dt/(2.*tau))*(rho_DSMC*theta_LB))/(c*c));
              PANINI_REAL Pxy_LB_g  = ((sigmaYZ_DSMC/f2g_factor + rho_DSMC*uY_DSMC*uZ_DSMC)/(c*c));
              PANINI_REAL Pyz_LB_g  = ((sigmaZX_DSMC/f2g_factor + rho_DSMC*uZ_DSMC*uX_DSMC)/(c*c));
              PANINI_REAL Pzx_LB_g  = ((sigmaXY_DSMC/f2g_factor + rho_DSMC*uX_DSMC*uY_DSMC)/(c*c));

              if((cellNumberLocal == 0) && (myRankY_LB == ((int)(0.5*p.nCoresY))) && (myRankX_LB == ((int)(0.5*p.nCoresX))))
              {
                  cout<<"rho_LB_Upper@: "<<step<<" "<<myRank<<":"<<rho_LB<<std::endl;
                  cout<<"uX_LB_Upper@ : "<<step<<" "<<myRank<<":"<<uX_LB<<std::endl;
                  cout<<"uY_LB_Upper@ : "<<step<<" "<<myRank<<":"<<uY_LB<<std::endl;
                  cout<<"uZ_LB_Upper@ : "<<step<<" "<<myRank<<":"<<uZ_LB<<std::endl;
                  cout<<"t0_LB_Upper@ : "<<step<<" "<<myRank<<":"<<theta_LB<<std::endl;
                  cout<<"sXX_LB_Upper@ : "<<step<<" "<<myRank<<":"<<Pxx_LB_g<<std::endl;
                  cout<<"SXY_LB_Upper@ : "<<step<<" "<<myRank<<":"<<Pxy_LB_g<<std::endl;
                  cout<<"sYZ_LB_Upper@ : "<<step<<" "<<myRank<<":"<<Pyz_LB_g<<std::endl;
                  cout<<"sYY_LB_Upper@ : "<<step<<" "<<myRank<<":"<<Pyy_LB_g<<std::endl;
                  cout<<"sZZ_LB_Upper@ : "<<step<<" "<<myRank<<":"<<Pzz_LB_g<<std::endl;
                  cout<<"sZX_LB_Upper@ : "<<step<<" "<<myRank<<":"<<Pzx_LB_g<<std::endl;
                  cout<<"qX_LB_Upper@ : "<<step<<" "<<myRank<<":"<<q1_LB_g<<std::endl;
                  cout<<"qY_LB_Upper@ : "<<step<<" "<<myRank<<":"<<q2_LB_g<<std::endl;
                  cout<<"qZ_LB_Upper@ : "<<step<<" "<<myRank<<":"<<q3_LB_g<<std::endl;
              }

              getGradThermalSinglePoint_heatAlpha (lbModel41,rho_LB,uX_LB,uY_LB,uZ_LB,theta_LB,Pxx_LB_g,Pyy_LB_g,Pzz_LB_g,Pxy_LB_g,Pyz_LB_g,Pzx_LB_g,q1_LB_g,q2_LB_g,q3_LB_g,grad);

              int xIndex = cellY + gridLBM41.nB1;
              int zIndex = cellX + gridLBM41.nB3;
              int yIndex =         gridLBM41.nE2 - 2;// - (0.25*p.nCoresZ - 1);

              copyFromArrayToNode(lbModel41,grad, gridLBM41,xIndex,yIndex,zIndex);
              copyFromArrayToCell(lbModel41,grad, gridLBM41,xIndex,yIndex,zIndex);

              copyFromArrayToNode(lbModel41,grad, gridLBM41,xIndex,yIndex + 1,zIndex);
              copyFromArrayToCell(lbModel41,grad, gridLBM41,xIndex,yIndex + 1,zIndex);

              copyFromArrayToNode(lbModel41,grad, gridLBM41,xIndex,yIndex + 2,zIndex);
              copyFromArrayToCell(lbModel41,grad, gridLBM41,xIndex,yIndex + 2,zIndex);


              copyFromArrayToNode(lbModel41,grad, gridLBM41,xIndex,yIndex + 3,zIndex);
              copyFromArrayToCell(lbModel41,grad, gridLBM41,xIndex,yIndex + 3,zIndex);

              copyFromArrayToNode(lbModel41,grad, gridLBM41,xIndex,yIndex + 4,zIndex);
              copyFromArrayToCell(lbModel41,grad, gridLBM41,xIndex,yIndex + 4,zIndex);

          }
      }
}

#endif // COUPLING_H_INCLUDED
