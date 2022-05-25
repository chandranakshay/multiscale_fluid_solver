/*********************************************************************************************
 *   Copyright (c) <2022>, <Santosh Ansumali@JNCASR>                                         *
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


#ifndef HARDSPHERE_H_INCLUDED
#define HARDSPHERE_H_INCLUDED

#include "cellList.h"

#include <time.h>
#include <mpi.h>
#include <chrono>
#include <random>

#define pi  3.141592653589793238462643383279502884197169
#define dia 1.
#define M   1.
#define kT  0.2948964908710636

std::mt19937 mrand(10);

template<int DIM,typename real,typename dof>
class hardSphere
{
    public:
        hardSphere(int XSize, int YSize, int ZSize, int nCx, int nCy, int nCz, int nPCell, int actualZcell)
        {
          XcellSize                        = XSize;
          YcellSize                        = YSize;
          ZcellSize                        = ZSize;
          nCoresX                          = nCx;
          nCoresY                          = nCy;
          nCoresZ                          = nCz;
          cellSizeZ                        = actualZcell*nCoresZ;
          numCells                         = XcellSize*YcellSize*cellSizeZ;
          nCores                           = (nCoresX*nCoresY*nCoresZ);
          numCellCore                      = (XcellSize/nCoresX)*(YcellSize/nCoresY)*(cellSizeZ/nCoresZ);
          numParticlesPerCell              = nPCell;
          nParticlesPerCore                = numCellCore*numParticlesPerCell;
          nParticles                       = cellSizeZ*XcellSize*YcellSize*numParticlesPerCell;
          XcellSizePerCore                 = XcellSize/nCoresX;
          YcellSizePerCore                 = YcellSize/nCoresY;
          ZcellSizePerCore                 = cellSizeZ/nCoresZ;
          xsubi[0]                         = 10;
          numMoments                       = 14;

          /************************PARTICLE DATA INITIALIZATION*********************/
          Oparticle.initializeParticles(XSize, YSize, actualZcell, nCx, nCy, nCz, nPCell);
          /*************************************************************************/

          cList                            = (cellList<DIM,real,dof>*) _mm_malloc(sizeof(cellList<DIM,real,dof>)*numCellCore, 4096);
          outParticlePlusY                 = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          outParticleMinsY                 = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          outParticlePlusPeriodicY         = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          outParticleMinsPeriodicY         = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          outParticlePlusX                 = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          outParticleMinsX                 = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          outParticlePlusPeriodicX         = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          outParticleMinsPeriodicX         = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);

          partCurrRankRecvPlusOneY         = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          partCurrRankRecvMinsOneY         = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          partCurrRankRecvPlusOnePeriodicY = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          partCurrRankRecvMinsOnePeriodicY = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          partCurrRankRecvPlusOneX         = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          partCurrRankRecvMinsOneX         = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          partCurrRankRecvPlusOnePeriodicX = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);
          partCurrRankRecvMinsOnePeriodicX = (real*) _mm_malloc(DIM*2*nParticlesPerCore*sizeof(real)/4, 4096);

          reGenerateData                   = (real*) _mm_malloc(DIM*2*2*nParticlesPerCore*sizeof(real)/ZcellSizePerCore, 4096);
          flagBoundary                     = (int*)  _mm_malloc(nParticlesPerCore*sizeof(int), 4096);
          boundaryTime                     = (real*) _mm_malloc(nParticlesPerCore*sizeof(int), 4096);
        }

        ~hardSphere()
        {
            _mm_free(cList);
            _mm_free(spaceList);
            _mm_free(voidList);
            _mm_free(pressureDSMC);
            _mm_free(cellStress);
            _mm_free(cellheatFl);
            _mm_free(pressureDSMCFlux);
            _mm_free(cellStressFlux);
            _mm_free(cellheatFlFlux);
            _mm_free(moments_Average);
            _mm_free(moments_Average_Flux);
            _mm_free(moments_Average_Flux_Temp);
            _mm_free(spaceAvgFlux);
            _mm_free(velDSMC);
            _mm_free(velDSMCFlux);

            _mm_free(outParticlePlusY);
            _mm_free(outParticleMinsY);
            _mm_free(outParticlePlusPeriodicY);
            _mm_free(outParticleMinsPeriodicY);
            _mm_free(outParticlePlusX);
            _mm_free(outParticleMinsX);
            _mm_free(outParticlePlusPeriodicX);
            _mm_free(outParticleMinsPeriodicX);

            _mm_free(partCurrRankRecvPlusOneY);
            _mm_free(partCurrRankRecvMinsOneY);
            _mm_free(partCurrRankRecvPlusOnePeriodicY);
            _mm_free(partCurrRankRecvMinsOnePeriodicY);
            _mm_free(partCurrRankRecvPlusOneX);
            _mm_free(partCurrRankRecvMinsOneX);
            _mm_free(partCurrRankRecvPlusOnePeriodicX);
            _mm_free(partCurrRankRecvMinsOnePeriodicX);

            _mm_free(reGenerateData);
            _mm_free(flagBoundary);
            _mm_free(boundaryTime);
            _mm_free(dsmcParticleFlux);
	          _mm_free(tempDataStorage);
        }

        void initHardSphere()
        {
	          volume = deltaX*deltaY*deltaZ;
            diameter = dia;
            extent[0] = Xbound;
            extent[1] = Ybound;
            extent[2] = height;

            mass = M;
        }

        void spaceMemoryAlloc()
        {
            spaceList            = (cellList<DIM,real,dof>*) _mm_malloc(sizeof(cellList<DIM,real,dof>)*numCellCore/(spaceAvgCells_X*spaceAvgCells_Y), 4096);
            moments_Average      = (real*) _mm_malloc(sizeof(real)*numCellCore*numMoments/(spaceAvgCells_X*spaceAvgCells_Y), 4096);
            moments_Average_Flux = (real*) _mm_malloc(sizeof(real)*numCellCore*numMoments/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y), 4096);
            moments_Average_Flux_Temp = (real*) _mm_malloc(sizeof(real)*numCellCore*numMoments/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y), 4096);
            velDSMC              = (real*) _mm_malloc(sizeof(real)*DIM*numCellCore/(spaceAvgCells_X*spaceAvgCells_Y), 4096);
            pressureDSMC         = (real*) _mm_malloc(sizeof(real)*numCellCore/(spaceAvgCells_X*spaceAvgCells_Y), 4096);
            cellStress           = (real*) _mm_malloc(sizeof(real)*DIM*2*numCellCore/(spaceAvgCells_X*spaceAvgCells_Y), 4096);
            cellheatFl           = (real*) _mm_malloc(sizeof(real)*DIM*numCellCore/(spaceAvgCells_X*spaceAvgCells_Y), 4096);

            velDSMCFlux          = (real*) _mm_malloc(sizeof(real)*DIM*numCellCore/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y), 4096);
            pressureDSMCFlux     = (real*) _mm_malloc(sizeof(real)*numCellCore/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y), 4096);
            cellStressFlux       = (real*) _mm_malloc(sizeof(real)*DIM*2*numCellCore/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y), 4096);
            cellheatFlFlux       = (real*) _mm_malloc(sizeof(real)*DIM*numCellCore/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y), 4096);

            spaceAvgFlux         = (real*) _mm_malloc(sizeof(real)*2*DIM*nParticlesPerCore/(ZcellSizePerCore), 4096);
            voidList             = (int*)  _mm_malloc(sizeof(int)*nParticlesPerCore*2./ZcellSizePerCore, 4096);

            dsmcParticleFlux     = (int*)  _mm_malloc((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y)*sizeof(int), 4096);
	          tempDataStorage		    = (real*) _mm_malloc(DIM*2*2*nParticlesPerCore*sizeof(real), 4096);
        }

        void extraVarInit(int myRank)
        {
            delXInv        = 1./deltaX;
            delYInv        = 1./deltaY;
            delZInv        = 1./deltaZ;

            nCoresXInv     = 1./((real)nCoresX);
            nCoresYInv     = 1./((real)nCoresY);
            nCoresZInv     = 1./((real)nCoresZ);
            nCoresInv      = 1./((real)nCores);

            coreNumberZ    = myRank/(nCoresX*nCoresY);
            coreNumberSlab = myRank%(nCoresX*nCoresY);
            coreNumberY    = coreNumberSlab/nCoresX;
            coreNumberX    = coreNumberSlab%nCoresX;

            XboundInv      = 1./Xbound;
            YboundInv      = 1./Ybound;
            heightInv      = 1./height;
        }

        real totKinEnergy();

        void MomentChainSpaceAvg(int,int,int&,int);

        void fgradPassSpaceAvg(int, real, real*, real*, real*, real*, int , int&, int, int);

        void MomentChainSpaceAvgFlux(int,int,int&,int,int);

        void fgradPassSpaceAvgFlux(int, real, real*, real*, real*, real*, int , int&, int, int);

        void evolveSystemDiffuse(int,int,int);

        void updateVelocityBottomWallX(int particleIndex);

        void updateVelocityBottomWallY(int particleIndex);

        void updateVelocityBottomWallZ(int particleIndex);

        void updateVelocityTopWall(int particleIndex);

        void updateVelocityTopWallX(int);

        void updateVelocityTopWallY(int);

        void updateVelocityTopWallZ(int);

        void XperiodicBoundary(int,int);

        void YperiodicBoundary(int,int);

        void clearCellParticles(int);

        void initPositions(int);

        void initVelocity();

        void updateCellCoord(int, int);

        real gaussianRandom();

        real relSpeed(int particleOne, int particleTwo);

        void Collision(int, int);

        void processCollision(int particleOne, int particleTwo);

        void initVelocitySteadyState(int coreNumber, real velScale);

        void DSMCclearAvgQts();

        void generateDistributionFromMoments(real* ,int cellX,int cellY, int myRank,int);

        real gaussianRandomLB2DSMC();

        void initGenerateDistributionFromMoments(real ,real ,int ,real);

        void projDistributionFromMoments(real* , int , int);

        void saveGeneratedData();

        void reGenerateParticles(int, int, int, real);

        void poissonParticleGen(real*, int, int, int, int);

	      void storeParticleDataToFile(int);

        void retrieveParticleDataFromFile(int);

        void pushParticleIndices(int, int, int, int);

        void sortVoidList();

        void totalMomentsWithFlux(real*, int, int, int);

        void clearMomentsFlux(int);

    //private:
        real Mach;
        real Re;
        real Kn;
        real gamma;
        real vSound;
        real Uc;
        real viscosity;
        real density;
        real mfp;
        real gy;
        int numMoments;
        int numCells;
        int numCellCore;
        int nParticles;
        int nParticlesPerCore;
        int XcellSize;
        int YcellSize;
        int ZcellSize;
        int cellSizeZ;
        real __n0;
        cellList<DIM,real,dof>* cList;
        cellList<DIM,real,dof>* spaceList;
        int* voidList;
        particle<DIM,real,dof> Oparticle;

        real uWall;
        int* flagBoundary;
        real* boundaryTime;

        int nE;

        real mass;
        real volume;
        real deltat;
        real height;
        real Xbound;
        real Ybound;
        real deltaX;
        real deltaY;
        real deltaZ;
        real diameter;
        real extent[DIM];
        real* vCells;
        real time_r = 0.;

        real* outParticlePlusY;
        real* outParticleMinsY;
        real* outParticlePlusPeriodicY;
        real* outParticleMinsPeriodicY;
        real* outParticlePlusX;
        real* outParticleMinsX;
        real* outParticlePlusPeriodicX;
        real* outParticleMinsPeriodicX;


        real* reGenerateData;
        real* tempReGenerateData;
        int*  particlesInSpaceAvgCells;
        int*  tempParticlesInSpaceAvgCells;
        int spaceAvgCells_X;
        int spaceAvgCells_Y;

        int outPlusCountY = 0;
        int outMinsCountY = 0;
        int outPlusCountPeriodicY = 0;
        int outMinsCountPeriodicY = 0;
        int outPlusCountX = 0;
        int outMinsCountX = 0;
        int outPlusCountPeriodicX = 0;
        int outMinsCountPeriodicX = 0;
        int countParticleCurrCore = 0;
        int totParticleCurrCore = 0;

        real* partCurrRankRecvPlusOneY;
        real* partCurrRankRecvMinsOneY;
        real* partCurrRankRecvPlusOnePeriodicY;
        real* partCurrRankRecvMinsOnePeriodicY;
        real* partCurrRankRecvPlusOneX;
        real* partCurrRankRecvMinsOneX;
        real* partCurrRankRecvPlusOnePeriodicX;
        real* partCurrRankRecvMinsOnePeriodicX;

        real* moments_Average;
        real* moments_Average_Flux;
        real* moments_Average_Flux_Temp;
	      unsigned short xsubi[3];
        int sumParticles;

        real* rhoAvg;
        real* vCellsAvg;
        real* stress;
        real* stressAvg;
        real* heatFlux;
        real* pressureDSMC;
        real* velDSMC;
        real* cellStress;
        real* cellheatFl;
        real* pressureDSMCFlux;
        real* velDSMCFlux;
        real* cellStressFlux;
        real* cellheatFlFlux;
        real* spaceAvgFlux;
        real* pressureAvg;
        real* heatFluxAvg;
        real fourthOrderMoment;
        real rho;

        int totParticleCurrCoreSave;
        real* tempDataStorage;

        int XcellSizePerCore;
        int YcellSizePerCore;
        int ZcellSizePerCore;

        int nCoresX;
        int nCoresY;
        int nCoresZ;
        int nCores;

        int numParticlesPerCell;

        real delXInv      ;
        real delYInv      ;
        real delZInv      ;
        real nCoresXInv   ;
        real nCoresYInv   ;
        real nCoresZInv   ;
        real nCoresInv    ;
        int coreNumberZ   ;
        int coreNumberSlab;
        int coreNumberY   ;
        int coreNumberX   ;

        real XboundInv    ;
        real YboundInv    ;
        real heightInv    ;

        int reGenIndex    ;
        int spaceAvgFluxIndex;
        int voidListIndex;
        int voidListSize;

        int* dsmcParticleFlux;

        md::rng molrng;
};

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::storeParticleDataToFile(int myRank)
{
    char restartFile[150];
    sprintf(restartFile,"./../data/particleData/particleData_%d.dat",myRank);

    std::ofstream outData;
    outData.open(restartFile,std::ios::out|std::ios::binary);
    outData.write((char *)&totParticleCurrCore,   sizeof(int));
    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
        {
            outData.write((char *)&(Oparticle(particleIndex,dim).position()),   sizeof(real));
            outData.write((char *)&(Oparticle(particleIndex,dim).velocity()),   sizeof(real));
            //outData<<Oparticle(particleIndex,dim).position()<<endl;
            //outData<<Oparticle(particleIndex,dim).velocity()<<endl;
        }
    }

    outData.close();
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::retrieveParticleDataFromFile(int myRank)
{
    char restartFile[150];
    sprintf(restartFile,"./../data/particleData/particleData_%d.dat",myRank);
    //std::string fileNameData = "./../particleData/particleData" + to_string(myRank) + ".dat";
    std::ifstream inData;
    inData.open(restartFile,std::ios::in|std::ios::binary);
    inData.read((char *)&totParticleCurrCore,   sizeof(int));
    //std::string tempData;
    std::string line;
    cout<<"dataCount: "<<totParticleCurrCore<<endl;

    //int dataCount = 0;

    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
        {
            inData.read((char *)&(Oparticle(particleIndex,dim).position()),   sizeof(real));
            inData.read((char *)&(Oparticle(particleIndex,dim).velocity()),   sizeof(real));
            //outData<<Oparticle(particleIndex,dim).position()<<endl;
            //outData<<Oparticle(particleIndex,dim).velocity()<<endl;
        }
    }
    //inData>>tempData;
    //tempDataStorage[dataCount] = atof(tempData.c_str());
    /*if(inData.is_open())
    {
    	while(getline(inData,line))
    	{
        	//inData.read((char*)&tempDataStorage[dataCount], sizeof(real));
        	tempDataStorage[dataCount] = atof(line.c_str());
        	dataCount++;
    	}
    }*/


    inData.close();

    //totParticleCurrCore = dataCount/6;

    /*for(int dim = 1; dim <= DIM; dim++)
    {
        for(int i = 0; i < totParticleCurrCore; i++)
        {
            Oparticle(i,dim).position() = tempDataStorage[totParticleCurrCore*2*(dim - 1) + 2*i    ];
            Oparticle(i,dim).velocity() = tempDataStorage[totParticleCurrCore*2*(dim - 1) + 2*i + 1];
        }
    }*/



    for(int i = 0; i < totParticleCurrCore; i++)
    {
            if((i == 2) && (myRank == 0))
            {
                cout<<"retrievecheck: "<<myRank<<"\t"<<i<<"\t"<<Oparticle(i,1).position()<<"\t"<<Oparticle(i,2).velocity()<<"\t"<<2.*Oparticle(i,3).position()<<"\t"<<Oparticle(i,1).velocity()*2.<<endl;
            }

	    int cellCoordX         = (int)((Oparticle(i,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
            int cellCoordY         = (int)((Oparticle(i,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
	    real heightNew;
            if(coreNumberZ < ((int)(0.5*nCoresZ)))
            {
                heightNew = Oparticle(i,3).position();
            }
            else
            {
                heightNew = height - Oparticle(i,3).position();
            }
            int cellCoordZ      = (int)(heightNew*delZInv);
            //int numSpaceAvgCells   = numCellCore/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y);

            //int cellNumberSpaceAvg = numSpaceAvgCells*cellCoordZ + (XcellSizePerCore/spaceAvgCells_X)*(cellCoordY/spaceAvgCells_Y) + (cellCoordX/spaceAvgCells_X);
            int cellNumberSpaceAvg = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);
            //if(myRank == 0)
	    //{
		//	std::cout<<"cellextra: "<<myRank<<"\t"<<cellNumberSpaceAvg<<"\t"<<numCellCore/(spaceAvgCells_X*spaceAvgCells_Y)<<"\t"<<cellCoordZ<<"\t"<<cellCoordX<<"\t"<<cellCoordY<<std::endl;
	    //}

            spaceList[cellNumberSpaceAvg].cellParticleList[spaceList[cellNumberSpaceAvg].cellParticleSize] = i;
            spaceList[cellNumberSpaceAvg].cellParticleSize++;
            /*if(myRank == 0)
            {
            	std::cout<<"spaceSize: "<<myRank<<"\t"<<cellNumberSpaceAvg<<"\t"<<spaceList[cellNumberSpaceAvg].cellParticleList.size()<<std::endl;
            }*/
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::saveGeneratedData()
{
    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
        {
            reGenerateData[totParticleCurrCore*2*(dim - 1) + 2*particleIndex    ] = Oparticle(particleIndex,dim).position();
            reGenerateData[totParticleCurrCore*2*(dim - 1) + 2*particleIndex + 1] = Oparticle(particleIndex,dim).velocity();
        }
    }

    totParticleCurrCoreSave = totParticleCurrCore;
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::pushParticleIndices(int myRank, int count, int spaceCellNumberNew, int particlesToBeRemoved)
{
    int countParticleCurrCoreSpace = 0;
    int spaceParticleIndexSize     = spaceList[spaceCellNumberNew].cellParticleSize;
    int spaceParticleEndIndex      = spaceParticleIndexSize - particlesToBeRemoved;

    for(int particleIndexIterator = (spaceParticleIndexSize - 1); particleIndexIterator >= spaceParticleEndIndex; particleIndexIterator--)
    {
        int particleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndexIterator];
        for(int currParticle = (totParticleCurrCore - 1); currParticle >= particleIndex; currParticle--)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                Oparticle(currParticle,dim).position() = Oparticle(currParticle + 1,dim).position();
                Oparticle(currParticle,dim).velocity() = Oparticle(currParticle + 1,dim).velocity();
            }
        }
    }

    countParticleCurrCoreSpace += spaceParticleIndexSize;

    totParticleCurrCore -= countParticleCurrCoreSpace;

    // // if(myRank == 0)
    // {
    //     std::cout<<"pushParticleIndices: "<<count<<"\t"<<myRank<<"\t"<<countParticleCurrCoreSpace<<"\t"<<totParticleCurrCore<<std::endl;
    // }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::reGenerateParticles(int myRank, int cellX, int cellY, real rho_DSMC)
{
    int cellZ              = ZcellSizePerCore - 1;
    int spaceCellNumber    = (XcellSizePerCore/spaceAvgCells_X)*cellY + cellX;
    int numSpaceAvgCells   = numCellCore/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y);
    int spaceCellNumberNew = numSpaceAvgCells*cellZ + spaceCellNumber;

    int totP    = ((int)(rho_DSMC*volume))*spaceAvgCells_X*spaceAvgCells_Y + 1;

    int spaceCellSizeOld   = spaceList[spaceCellNumberNew].cellParticleSize;

    for(int numP = 0; numP < totP; numP++)
    {
	      int particleIndex;
        if(numP < spaceCellSizeOld)
        {
            particleIndex = spaceList[spaceCellNumberNew].cellParticleList[numP];
        }
        else
        {
            particleIndex = totParticleCurrCore;
            spaceList[spaceCellNumberNew].cellParticleList[spaceList[spaceCellNumberNew].cellParticleSize] = particleIndex;
            spaceList[spaceCellNumberNew].cellParticleSize++;
            totParticleCurrCore++;
        }

        Oparticle(particleIndex,1).position() = reGenerateData[(reGenIndex + numP)*DIM*2    ];
        Oparticle(particleIndex,1).velocity() = reGenerateData[(reGenIndex + numP)*DIM*2 + 1];

        Oparticle(particleIndex,2).position() = reGenerateData[(reGenIndex + numP)*DIM*2 + 2];
        Oparticle(particleIndex,2).velocity() = reGenerateData[(reGenIndex + numP)*DIM*2 + 3];

        Oparticle(particleIndex,3).position() = reGenerateData[(reGenIndex + numP)*DIM*2 + 4];
        Oparticle(particleIndex,3).velocity() = reGenerateData[(reGenIndex + numP)*DIM*2 + 5];
    }

    reGenIndex += totP;
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::generateDistributionFromMoments(real* momentsRecv,int cellX2, int cellY2, int myRank, int countTime)
{
    real rho_DSMC          = momentsRecv[ 0];
    real p_DSMC            = momentsRecv[ 1];
    real uX_DSMC           = momentsRecv[ 2];
    real uY_DSMC           = momentsRecv[ 3];
    real uZ_DSMC           = momentsRecv[ 4];
    real sigmaXX_DSMC      = momentsRecv[ 5];
    real sigmaXY_DSMC      = momentsRecv[ 6];
    real sigmaZX_DSMC      = momentsRecv[ 7];
    real sigmaYY_DSMC      = momentsRecv[ 8];
    real sigmaYZ_DSMC      = momentsRecv[ 9];
    real sigmaZZ_DSMC      = momentsRecv[10];
    real qX_DSMC           = momentsRecv[11];
    real qY_DSMC           = momentsRecv[12];
    real qZ_DSMC           = momentsRecv[13];

    int cellZ2             = ZcellSizePerCore - 1;
    int spaceCellNumber    = (XcellSizePerCore/spaceAvgCells_X)*cellY2 + cellX2;
    int numSpaceAvgCells   = numCellCore/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y);
    int spaceCellNumberNew = numSpaceAvgCells*cellZ2 + spaceCellNumber;

    int spaceCellSizeOld   = spaceList[spaceCellNumberNew].cellParticleSize;

    // if(myRank == 0)
    // {
    //     std::cout<<"spaceCellCountBeforeGen: "<<spaceCellNumberNew<<"\t"<<spaceCellSizeOld<<std::endl;
    // }

    real heightRestricted  = cellSizeZ*height/((real)ZcellSize);      // For Top and Bottom slabs combined

    int rhoNumParticles    = ((int)(rho_DSMC*volume))*spaceAvgCells_X*spaceAvgCells_Y + 2;
    std::default_random_engine generator;
    std::poisson_distribution<int> distribution(rhoNumParticles);

    int rhoPoissonNumParticles = distribution(mrand);
    real rhoPoisson = (real)rhoPoissonNumParticles/((double)(volume*spaceAvgCells_X*spaceAvgCells_Y));

    int numP = 0;
    int totP = rhoPoissonNumParticles;

    // if((myRank == 0) && (countTime >= 240))
    // {
    //     std::cout<<"genParticle: "<<myRank<<"\t"<<cellX2<<"\t"<<cellY2<<"\t"<<rho_DSMC<<"\t"<<volume<<"\t"<<spaceAvgCells_X<<"\t"<<spaceAvgCells_Y<<"\t"<<rhoNumParticles<<"\t"<<rhoPoissonNumParticles<<"\t"<<spaceCellSizeOld<<"\t"<<totParticleCurrCore<<"\t"<<nParticlesPerCore<<std::endl;
    // }

    real fneq;
    real theta_DSMC = p_DSMC/rhoPoisson;
    real bigM = 0.;

    {
        real a0 = rhoPoisson;
        if(fabs(a0) > bigM)
        {
            bigM = a0;
        }
        real a1 = rhoPoisson*uX_DSMC;
        if(fabs(a1) > bigM)
        {
            bigM = a1;
        }
        real a2 = rhoPoisson*uY_DSMC;
        if(fabs(a2) > bigM)
        {
            bigM = a2;
        }
        real a3 = rhoPoisson*uZ_DSMC;
        if(fabs(a3) > bigM)
        {
            bigM = a3;
        }
        real a11 = sigmaXX_DSMC - rhoPoisson*theta_DSMC;
        if(fabs(a11) > bigM)
        {
            bigM = a11;
        }
        real a22 = sigmaYY_DSMC - rhoPoisson*theta_DSMC;
        if(fabs(a22) > bigM)
        {
            bigM = a22;
        }
        real a33 = sigmaZZ_DSMC - rhoPoisson*theta_DSMC;
        if(fabs(a33) > bigM)
        {
            bigM = a33;
        }
        real a12 = sigmaXY_DSMC;
        if(fabs(a12) > bigM)
        {
            bigM = a12;
        }
        real a23 = sigmaYZ_DSMC;
        if(fabs(a23) > bigM)
        {
            bigM = a23;
        }
        real a31 = sigmaZX_DSMC;
        if(fabs(a31) > bigM)
        {
            bigM = a31;
        }
    }

    real c = 2.*bigM;

    real theta_DSMCInvSq   = 1./(((real)theta_DSMC*theta_DSMC));
    real theta_DSMCInvCube = theta_DSMCInvSq/((real)theta_DSMC);

    real uSq = uX_DSMC*uX_DSMC + uY_DSMC*uY_DSMC + uZ_DSMC*uZ_DSMC;

    while(numP < totP)
    {
    		real zetaX = molrng.normal()*sqrt(theta_DSMC);
        real zetaY = molrng.normal()*sqrt(theta_DSMC);
        real zetaZ;

        if(coreNumberZ < ((int)(0.5*nCoresZ)))
        {
            zetaZ = -molrng.normal()*sqrt(theta_DSMC);
        }
        else
        {
            zetaZ = molrng.normal()*sqrt(theta_DSMC);
        }

        real zetaSq = zetaX*zetaX + zetaY*zetaY + zetaZ*zetaZ;

        real term2 = .5*((sigmaXX_DSMC)*(zetaX*zetaX - theta_DSMC)
                   + 2.*sigmaXY_DSMC*zetaX*zetaY + (sigmaYY_DSMC)*(zetaY*zetaY - theta_DSMC)
                   + 2.*sigmaYZ_DSMC*zetaY*zetaZ + 2.*sigmaZX_DSMC*zetaZ*zetaX
                   + (sigmaZZ_DSMC)*(zetaZ*zetaZ - theta_DSMC))*theta_DSMCInvSq;

        real term3 = .1*(qX_DSMC*zetaX + qY_DSMC*zetaY + qZ_DSMC*zetaZ)*(zetaSq - 5.0*theta_DSMC)*theta_DSMCInvCube;

        fneq = (rhoPoisson + term2 + term3);

    		real selectCriteriaRand = c*molrng.uniform();    //rho_DSMC is the amplitude of the new distribution

    		if(selectCriteriaRand <= fneq)
    		{
            int particleIndex;
            if(numP < spaceCellSizeOld)
            {
                particleIndex = spaceList[spaceCellNumberNew].cellParticleList[numP];
            }
            else
            {
                particleIndex = totParticleCurrCore;
                spaceList[spaceCellNumberNew].cellParticleList[spaceList[spaceCellNumberNew].cellParticleSize] = particleIndex;
                spaceList[spaceCellNumberNew].cellParticleSize++;
                totParticleCurrCore++;
            }

            Oparticle(particleIndex,1).velocity() = zetaX + uX_DSMC;
            Oparticle(particleIndex,2).velocity() = zetaY + uY_DSMC;
            Oparticle(particleIndex,3).velocity() = zetaZ + uZ_DSMC;

            Oparticle(particleIndex,1).position() = Xbound*(real)coreNumberX*nCoresXInv + ((real)cellX2 + molrng.uniform())*deltaX*spaceAvgCells_X;
            Oparticle(particleIndex,2).position() = Ybound*(real)coreNumberY*nCoresYInv + ((real)cellY2 + molrng.uniform())*deltaY*spaceAvgCells_Y;

            if(coreNumberZ < ((int)(0.5*nCoresZ)))
            {
                Oparticle(particleIndex,3).position() = heightRestricted*(real)(coreNumberZ*nCoresZInv) + (cellZ2 + molrng.uniform())*deltaZ;
            }
            else
            {
                int coreNumberNewZ = (nCoresZ - 1 - coreNumberZ);
                Oparticle(particleIndex,3).position() = height - (heightRestricted*(real)(coreNumberNewZ*nCoresZInv) + (cellZ2 + molrng.uniform())*deltaZ);
            }

            reGenerateData[(reGenIndex + numP)*DIM*2    ] = Oparticle(particleIndex,1).position();
            reGenerateData[(reGenIndex + numP)*DIM*2 + 1] = Oparticle(particleIndex,1).velocity();

            reGenerateData[(reGenIndex + numP)*DIM*2 + 2] = Oparticle(particleIndex,2).position();
            reGenerateData[(reGenIndex + numP)*DIM*2 + 3] = Oparticle(particleIndex,2).velocity();

            reGenerateData[(reGenIndex + numP)*DIM*2 + 4] = Oparticle(particleIndex,3).position();
            reGenerateData[(reGenIndex + numP)*DIM*2 + 5] = Oparticle(particleIndex,3).velocity();

      			numP++;
    		}
    }

    if(totP < spaceCellSizeOld)
    {
        for(int indexBeyondtotP = 0; indexBeyondtotP < (spaceCellSizeOld - totP); indexBeyondtotP ++)
        {
            voidList[voidListSize] = spaceList[spaceCellNumberNew].cellParticleList[spaceCellSizeOld - indexBeyondtotP - 1];
            voidListSize++;
            spaceList[spaceCellNumberNew].cellParticleSize--;
        }
    }

    // if(myRank == 0)
    // {
    //     std::cout<<"spaceCellCountAfterGen: "<<spaceCellNumberNew<<"\t"<<spaceList[spaceCellNumberNew].cellParticleSize<<std::endl;
    // }

    reGenIndex += totP;

    real uSum[DIM] = {};

    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = 0; particleIndex < totP; particleIndex++)
        {
            int spaceParticleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndex];
            uSum[dim - 1] += Oparticle(spaceParticleIndex,dim).velocity();
        }
        uSum[dim - 1] = (uSum[dim - 1])/((real)totP);
    }

    uSum[0] = (uX_DSMC - uSum[0]);
    uSum[1] = (uY_DSMC - uSum[1]);
    uSum[2] = (uZ_DSMC - uSum[2]);

    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = 0; particleIndex < totP; particleIndex++)
        {
            int spaceParticleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndex];
            Oparticle(spaceParticleIndex,dim).velocity() += uSum[dim - 1];
        }
    }

    real uxAnalytical = uX_DSMC;
    real uyAnalytical = uY_DSMC;
    real uzAnalytical = uZ_DSMC;
    real KE = 0.;

    for(int particleIndex = 0; particleIndex < totP; particleIndex++)
    {
        int spaceParticleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndex];

        KE += (Oparticle(spaceParticleIndex,1).velocity() - uxAnalytical)*(Oparticle(spaceParticleIndex,1).velocity() - uxAnalytical) + (Oparticle(spaceParticleIndex,2).velocity() - uyAnalytical)*(Oparticle(spaceParticleIndex,2).velocity() - uyAnalytical) + (Oparticle(spaceParticleIndex,3).velocity() - uzAnalytical)*(Oparticle(spaceParticleIndex,3).velocity() - uzAnalytical);
    }

    real velTotSqScale = sqrt(KE/(totP*3.*theta_DSMC));
    real velTotSqScaleInv = 1./velTotSqScale;


    for(int particleIndex = 0; particleIndex < totP; particleIndex++)
    {
        int spaceParticleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndex];
        Oparticle(spaceParticleIndex,1).velocity() = (Oparticle(spaceParticleIndex,1).velocity() - uxAnalytical)*velTotSqScaleInv + uxAnalytical;
    }
    for(int particleIndex = 0; particleIndex < totP; particleIndex++)
    {
        int spaceParticleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndex];
        Oparticle(spaceParticleIndex,2).velocity() = (Oparticle(spaceParticleIndex,2).velocity() - uyAnalytical)*velTotSqScaleInv + uyAnalytical;
    }
    for(int particleIndex = 0; particleIndex < totP; particleIndex++)
    {
        int spaceParticleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndex];
        Oparticle(spaceParticleIndex,3).velocity() = (Oparticle(spaceParticleIndex,3).velocity() - uzAnalytical)*velTotSqScaleInv + uzAnalytical;
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::poissonParticleGen(real* momentsRecv, int cellX, int cellY, int myRank, int tempCount)
{
    real rho_DSMC          = momentsRecv[ 0];
    real p_DSMC            = momentsRecv[ 1];
    real uX_DSMC           = momentsRecv[ 2];
    real uY_DSMC           = momentsRecv[ 3];
    real uZ_DSMC           = momentsRecv[ 4];
    real sigmaXX_DSMC      = momentsRecv[ 5];
    real sigmaXY_DSMC      = momentsRecv[ 6];
    real sigmaZX_DSMC      = momentsRecv[ 7];
    real sigmaYY_DSMC      = momentsRecv[ 8];
    real sigmaYZ_DSMC      = momentsRecv[ 9];
    real sigmaZZ_DSMC      = momentsRecv[10];
    real qX_DSMC           = momentsRecv[11];
    real qY_DSMC           = momentsRecv[12];
    real qZ_DSMC           = momentsRecv[13];
    real heightRestricted  = cellSizeZ*height/((real)ZcellSize);      // For Top and Bottom slabs combined

    int rhoNumParticles    = ((int)(rho_DSMC*volume))*spaceAvgCells_X*spaceAvgCells_Y;
    std::default_random_engine generator;
    std::poisson_distribution<int> distribution(rhoNumParticles);

    int rhoPoissonNumParticles = distribution(mrand) + 10;

    int cellZ              = ZcellSizePerCore - 1;
    int spaceCellNumber    = (XcellSizePerCore/spaceAvgCells_X)*cellY + cellX;
    int numSpaceAvgCells   = numCellCore/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y);
    int spaceCellNumberNew = numSpaceAvgCells*cellZ + spaceCellNumber;

    real rhoPoisson        = (real)rhoPoissonNumParticles/((double)(volume*spaceAvgCells_X*spaceAvgCells_Y));
    real theta_DSMC        = p_DSMC/rhoPoisson;

    int spaceCellSizeOld = spaceList[spaceCellNumberNew].cellParticleSize;
    int totP = rhoPoissonNumParticles - spaceCellSizeOld;
    int countVoidList = 0;
    int voidListSizeOld = voidListSize;

    if(totP < 0)
    {
        for(int currVoid = 0; currVoid < (-totP); currVoid++)
        {
            voidList[voidListSize] = spaceList[spaceCellNumberNew].cellParticleList[spaceCellSizeOld - currVoid - 1];
            voidListSize++;
            spaceList[spaceCellNumberNew].cellParticleSize--;
        }
    }
    else
    {
        real fneq;
        real bigM = 0.;

        {
            real a0 = rhoPoisson;
            if(fabs(a0) > bigM)
            {
                bigM = a0;
            }
            real a1 = rhoPoisson*uX_DSMC;
            if(fabs(a1) > bigM)
            {
                bigM = a1;
            }
            real a2 = rhoPoisson*uY_DSMC;
            if(fabs(a2) > bigM)
            {
                bigM = a2;
            }
            real a3 = rhoPoisson*uZ_DSMC;
            if(fabs(a3) > bigM)
            {
                bigM = a3;
            }
            real a11 = sigmaXX_DSMC - rhoPoisson*theta_DSMC;
            if(fabs(a11) > bigM)
            {
                bigM = a11;
            }
            real a22 = sigmaYY_DSMC - rhoPoisson*theta_DSMC;
            if(fabs(a22) > bigM)
            {
                bigM = a22;
            }
            real a33 = sigmaZZ_DSMC - rhoPoisson*theta_DSMC;
            if(fabs(a33) > bigM)
            {
                bigM = a33;
            }
            real a12 = sigmaXY_DSMC;
            if(fabs(a12) > bigM)
            {
                bigM = a12;
            }
            real a23 = sigmaYZ_DSMC;
            if(fabs(a23) > bigM)
            {
                bigM = a23;
            }
            real a31 = sigmaZX_DSMC;
            if(fabs(a31) > bigM)
            {
                bigM = a31;
            }
        }

        real c = 2.*bigM;

        real theta_DSMCInvSq   = 1./(((real)theta_DSMC*theta_DSMC));
        real theta_DSMCInvCube = theta_DSMCInvSq/((real)theta_DSMC);

        real uSq = uX_DSMC*uX_DSMC + uY_DSMC*uY_DSMC + uZ_DSMC*uZ_DSMC;

        int numP = 0;

        while(numP < totP)
        {
            real zetaX = molrng.normal()*sqrt(theta_DSMC);
            real zetaY = molrng.normal()*sqrt(theta_DSMC);
            real zetaZ = molrng.normal()*sqrt(theta_DSMC);

            real zetaSq = zetaX*zetaX + zetaY*zetaY + zetaZ*zetaZ;

            real term2 = .5*((sigmaXX_DSMC)*(zetaX*zetaX - theta_DSMC)
                       + 2.*sigmaXY_DSMC*zetaX*zetaY + (sigmaYY_DSMC)*(zetaY*zetaY - theta_DSMC)
                       + 2.*sigmaYZ_DSMC*zetaY*zetaZ + 2.*sigmaZX_DSMC*zetaZ*zetaX
                       + (sigmaZZ_DSMC)*(zetaZ*zetaZ - theta_DSMC))*theta_DSMCInvSq;

            real term3 = .1*(qX_DSMC*zetaX + qY_DSMC*zetaY + qZ_DSMC*zetaZ)*(zetaSq - 5.0*theta_DSMC)*theta_DSMCInvCube;

            fneq = (rhoPoisson + term2 + term3);

            real selectCriteriaRand = c*molrng.uniform();    //rho_DSMC is the amplitude of the new distribution

            if(selectCriteriaRand <= fneq)
            {
                int particleIndex;
                if(numP < voidListSizeOld)
                {
                    particleIndex = voidList[voidListSize - 1];
                    voidListSize--;
                    countVoidList++;
                }
                else
                {
                    particleIndex = totParticleCurrCore;
                    totParticleCurrCore++;
                }

                Oparticle(particleIndex,1).velocity() = zetaX + uX_DSMC;
                Oparticle(particleIndex,2).velocity() = zetaY + uY_DSMC;
                Oparticle(particleIndex,3).velocity() = zetaZ + uZ_DSMC;

                Oparticle(particleIndex,1).position() = Xbound*(real)coreNumberX*nCoresXInv + ((real)cellX + molrng.uniform())*deltaX*spaceAvgCells_X;
                Oparticle(particleIndex,2).position() = Ybound*(real)coreNumberY*nCoresYInv + ((real)cellY + molrng.uniform())*deltaY*spaceAvgCells_Y;

                if(coreNumberZ < (0.5*nCoresZ))
                {
                    Oparticle(particleIndex,3).position() = heightRestricted*(real)(coreNumberZ*nCoresZInv) + (cellZ + molrng.uniform())*deltaZ;
                }
                else
                {
                    int coreNumberNewZ = (nCoresZ - 1 - coreNumberZ);
                    Oparticle(particleIndex,3).position() = height - (heightRestricted*(real)(coreNumberNewZ*nCoresZInv) + (cellZ + molrng.uniform())*deltaZ);
                }

                spaceList[spaceCellNumberNew].cellParticleList[spaceList[spaceCellNumberNew].cellParticleSize] = particleIndex;
                spaceList[spaceCellNumberNew].cellParticleSize++;

                numP++;
            }
        }
    }

    real uSum[DIM] = {};

    totP = spaceList[spaceCellNumberNew].cellParticleSize;

    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = 0; particleIndex < totP; particleIndex++)
        {
            int spaceParticleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndex];
            uSum[dim - 1] += Oparticle(spaceParticleIndex,dim).velocity();
        }
        uSum[dim - 1] = (uSum[dim - 1])/((real)totP);
    }

    uSum[0] = (uX_DSMC - uSum[0]);
    uSum[1] = (uY_DSMC - uSum[1]);
    uSum[2] = (uZ_DSMC - uSum[2]);

    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = 0; particleIndex < totP; particleIndex++)
        {
            int spaceParticleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndex];
            Oparticle(spaceParticleIndex,dim).velocity() += uSum[dim - 1];
        }
    }

    real uxAnalytical = uX_DSMC;
    real uyAnalytical = uY_DSMC;
    real uzAnalytical = uZ_DSMC;
    real KE = 0.;

    for(int particleIndex = 0; particleIndex < totP; particleIndex++)
    {
        int spaceParticleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndex];

        KE += (Oparticle(spaceParticleIndex,1).velocity() - uxAnalytical)*(Oparticle(spaceParticleIndex,1).velocity() - uxAnalytical) + (Oparticle(spaceParticleIndex,2).velocity() - uyAnalytical)*(Oparticle(spaceParticleIndex,2).velocity() - uyAnalytical) + (Oparticle(spaceParticleIndex,3).velocity() - uzAnalytical)*(Oparticle(spaceParticleIndex,3).velocity() - uzAnalytical);
    }

    real velTotSqScale = sqrt(KE/(totP*3.*theta_DSMC));
    real velTotSqScaleInv = 1./velTotSqScale;


    for(int particleIndex = 0; particleIndex < totP; particleIndex++)
    {
        int spaceParticleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndex];
        Oparticle(spaceParticleIndex,1).velocity() = (Oparticle(spaceParticleIndex,1).velocity() - uxAnalytical)*velTotSqScaleInv + uxAnalytical;
    }
    for(int particleIndex = 0; particleIndex < totP; particleIndex++)
    {
        int spaceParticleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndex];
        Oparticle(spaceParticleIndex,2).velocity() = (Oparticle(spaceParticleIndex,2).velocity() - uyAnalytical)*velTotSqScaleInv + uyAnalytical;
    }
    for(int particleIndex = 0; particleIndex < totP; particleIndex++)
    {
        int spaceParticleIndex = spaceList[spaceCellNumberNew].cellParticleList[particleIndex];
        Oparticle(spaceParticleIndex,3).velocity() = (Oparticle(spaceParticleIndex,3).velocity() - uzAnalytical)*velTotSqScaleInv + uzAnalytical;
    }

    // if(myRank == 0)
    // {
    //     std::cout<<"poissonParticleGen: "<<spaceCellNumberNew<<"\t"<<rhoNumParticles<<"\t"<<rhoPoissonNumParticles<<"\t"<<spaceCellSizeOld<<"\t"<<totP<<"\t"<<totParticleCurrCore<<"\t"<<voidListSize<<"\t"<<countVoidList<<"\t"<<rho_DSMC<<"\t"<<rho_DSMC*volume<<"\t"<<((int)(rho_DSMC*volume))*XcellSizePerCore*YcellSizePerCore*ZcellSizePerCore<<std::endl;
    // }

    if(totParticleCurrCore >= 2*nParticlesPerCore)
    {
        std::cout<<"totParticleCurrCoreExceed: "<<myRank<<"\t"<<totParticleCurrCore<<"\t"<<2.*nParticlesPerCore<<std::endl;
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::sortVoidList()
{
    // for(int i = 0; i < voidListSize; i++)
    // {
    //     std::cout<<"beforesorting: "<<voidList[0].cellParticleList[i]<<std::endl;
    // }

    for(int currElementGlobal = 0; currElementGlobal < (voidListSize - 1); currElementGlobal++)
    {
      for(int currElement = 0; currElement < (voidListSize - currElementGlobal - 1); currElement++)
      {
          int elementOne = voidList[currElement    ];
          int elementTwo = voidList[currElement + 1];

          if(elementOne > elementTwo)
          {
              voidList[currElement + 1] = elementOne;
              voidList[currElement    ] = elementTwo;
          }
      }
    }

    // for(int i = 0; i < voidListSize; i++)
    // {
    //     std::cout<<"aftersorting : "<<voidList[0].cellParticleList[i]<<std::endl;
    // }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::initGenerateDistributionFromMoments(real rho_DSMC, real p_DSMC, int myRank, real velScale)
{
    real uX_DSMC = 0.;
    real uY_DSMC = 0.;
    real uZ_DSMC = 0.;

    real heightByCores = height/((real)nCores);

    real heightRestricted = cellSizeZ*height/((real)ZcellSize);      // For Top and Bottom slabs combined

    for(int cellNumberLocal = 0; cellNumberLocal < numCellCore; cellNumberLocal++)
    {
        int numP = 0;
        int totP;

        int cellNumberPlane = cellNumberLocal%(XcellSizePerCore*YcellSizePerCore);
        int cellX           = cellNumberPlane%XcellSizePerCore;
        int cellY           = cellNumberPlane/XcellSizePerCore;
        int cellZ           = cellNumberLocal/(XcellSizePerCore*YcellSizePerCore);

        if(coreNumberZ < ((int)(0.5*nCoresZ)))
        {
            uY_DSMC = -velScale*(((double)cellZ*deltaZ/height)*((double)cellZ*deltaZ/height) - 0.5)*2.*Uc;

            if((myRank == 0) && (cellNumberLocal == 0))
            {
            	std::cout<<"uY_DSMCGen: "<<uY_DSMC<<"\t"<<velScale<<"\t"<<uWall<<std::endl;
            }
        }
        else
        {
            uY_DSMC = velScale*(((double)cellZ*deltaZ/height)*((double)cellZ*deltaZ/height) - 0.5)*2.*Uc;
        }

        int rhoNumParticles = (int)(rho_DSMC*volume) + 1;           // Per  DSMC cell number of particles

        if((myRank == 0) && (cellNumberLocal == 0))
        {
            std::cout<<"Initial generalised particles per cell: "<<rhoNumParticles<<std::endl;
        }

        totP = rhoNumParticles;

        //if(myRank == 0)
        //{
        //     std::cout<<"initGenParticle: "<<myRank<<"\t"<<cellX<<"\t"<<cellY<<"\t"<<rho_DSMC<<"\t"<<volume<<"\t"<<spaceAvgCells_X<<"\t"<<spaceAvgCells_Y<<"\t"<<rhoNumParticles<<std::endl;
        //}

        real fneq;
        real theta_DSMC = p_DSMC/rho_DSMC;

        real c = 2.*rho_DSMC;

        real theta_DSMCInvSq   = 1./(((real)theta_DSMC*theta_DSMC));
        real theta_DSMCInvCube = theta_DSMCInvSq/((real)theta_DSMC);

        real uSq = uX_DSMC*uX_DSMC + uY_DSMC*uY_DSMC + uZ_DSMC*uZ_DSMC;

        while(numP < totP)
        {
        		real zetaX = molrng.normal()*sqrt(theta_DSMC);
            real zetaY = molrng.normal()*sqrt(theta_DSMC);
            real zetaZ = molrng.normal()*sqrt(theta_DSMC);

            real zetaSq = zetaX*zetaX + zetaY*zetaY + zetaZ*zetaZ;

            real term2 = 0.;

            real term3 = 0.;

            fneq = (rho_DSMC + term2 + term3);

        		real selectCriteriaRand = c*molrng.uniform();    //rho_DSMC is the amplitude of the new distribution

        		if(selectCriteriaRand <= fneq)
        		{
                int particleIndex = totParticleCurrCore;
                totParticleCurrCore++;

		//std::cout<<"beforegen: "<<myRank<<std::endl;

                Oparticle(particleIndex,1).position() = Xbound*(real)(coreNumberX*nCoresXInv) + (cellX + molrng.uniform())*deltaX;
                Oparticle(particleIndex,1).velocity() = zetaX + uX_DSMC;
                Oparticle(particleIndex,2).position() = Ybound*(real)(coreNumberY*nCoresYInv) + (cellY + molrng.uniform())*deltaY;
                Oparticle(particleIndex,2).velocity() = zetaY + uY_DSMC;
		//std::cout<<"beforegen2: "<<myRank<<std::endl;

                if(coreNumberZ < ((int)(0.5*nCoresZ)))
                {
                    Oparticle(particleIndex,3).position() = heightRestricted*(real)(coreNumberZ*nCoresZInv) + (cellZ + molrng.uniform())*deltaZ;
                }
                else
                {
                    int coreNumberNewZ = nCoresZ - 1 - coreNumberZ;
                    Oparticle(particleIndex,3).position() = height - (heightRestricted*(real)(coreNumberNewZ*nCoresZInv) + (cellZ + molrng.uniform())*deltaZ);
                }
                Oparticle(particleIndex,3).velocity() = zetaZ + uZ_DSMC;
		//std::cout<<"beforegen3: "<<myRank<<std::endl;

                //int cellCoordX         = (int)((Oparticle(particleIndex,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
                //int cellCoordY         = (int)((Oparticle(particleIndex,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);

                int numSpaceAvgCells   = numCellCore/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y);

                int cellNumberSpaceAvg = numSpaceAvgCells*cellZ + (XcellSizePerCore/spaceAvgCells_X)*(cellY/spaceAvgCells_Y) + (cellX/spaceAvgCells_X);

		if(cellNumberSpaceAvg >= numCellCore/(spaceAvgCells_X*spaceAvgCells_Y))
		{
			std::cout<<"cellextra: "<<myRank<<"\t"<<cellNumberSpaceAvg<<"\t"<<numCellCore/(spaceAvgCells_X*spaceAvgCells_Y)<<"\t"<<cellZ<<"\t"<<cellX<<"\t"<<cellY<<std::endl;
		}

                spaceList[cellNumberSpaceAvg].cellParticleList[spaceList[cellNumberSpaceAvg].cellParticleSize] = particleIndex;
                spaceList[cellNumberSpaceAvg].cellParticleSize++;

          			numP++;
        		}
        }

    	real uxAnalytical = uX_DSMC;
    	real uyAnalytical = uY_DSMC;
    	real uzAnalytical = uZ_DSMC;
    	real KE = 0.;

    	for(int particleIndex = (totParticleCurrCore - totP); particleIndex < totParticleCurrCore; particleIndex++)
    	{
        	KE += (Oparticle(particleIndex,1).velocity() - uxAnalytical)*(Oparticle(particleIndex,1).velocity() - uxAnalytical) + (Oparticle(particleIndex,2).velocity() - uyAnalytical)*(Oparticle(particleIndex,2).velocity() - uyAnalytical) + (Oparticle(particleIndex,3).velocity() - uzAnalytical)*(Oparticle(particleIndex,3).velocity() - uzAnalytical);
    	}

    	real velTotSqScale = sqrt(KE/(totP*3.*theta_DSMC));
    	real velTotSqScaleInv = 1./velTotSqScale;


    	for(int particleIndex = (totParticleCurrCore - totP); particleIndex < totParticleCurrCore; particleIndex++)
    	{
        	Oparticle(particleIndex,1).velocity() = (Oparticle(particleIndex,1).velocity() - uxAnalytical)*velTotSqScaleInv + uxAnalytical;
    	}
    	for(int particleIndex = (totParticleCurrCore - totP); particleIndex < totParticleCurrCore; particleIndex++)
    	{
        	Oparticle(particleIndex,2).velocity() = (Oparticle(particleIndex,2).velocity() - uyAnalytical)*velTotSqScaleInv + uyAnalytical;
    	}
    	for(int particleIndex = (totParticleCurrCore - totP); particleIndex < totParticleCurrCore; particleIndex++)
    	{
        	Oparticle(particleIndex,3).velocity() = (Oparticle(particleIndex,3).velocity() - uzAnalytical)*velTotSqScaleInv + uzAnalytical;
    	}

        real uSum[DIM] = {};

    	for(int dim = 1; dim <= DIM; dim++)
    	{
        	for(int particleIndex = (totParticleCurrCore - totP); particleIndex < totParticleCurrCore; particleIndex++)
        	{
            		uSum[dim - 1] += Oparticle(particleIndex,dim).velocity();
        	}
        	uSum[dim - 1] = (uSum[dim - 1])/((real)totP);
    	}

    	uSum[0] = (uX_DSMC - uSum[0]);
    	uSum[1] = (uY_DSMC - uSum[1]);
    	uSum[2] = (uZ_DSMC - uSum[2]);

    	for(int dim = 1; dim <= DIM; dim++)
    	{
        	for(int particleIndex = (totParticleCurrCore - totP); particleIndex < totParticleCurrCore; particleIndex++)
        	{
            		Oparticle(particleIndex,dim).velocity() += uSum[dim - 1];
        	}
    	}
    }



    //std::cout<<"outsideInitGen: "<<myRank<<"\t"<<totParticleCurrCore<<"\t"<<0<<"\t"<<spaceList[numCellCore*(ZcellSizePerCore - 1)/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y)].cellParticleList.size()<<std::endl;
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::projDistributionFromMoments(real* momentsDSMCProj, int currCell, int myRank)
{
    // uX_DSMC = 0.;
    // uZ_DSMC = 0.;

    int numSpaceAvgCells_DSMC_X = spaceAvgCells_X;
    int numSpaceAvgCells_DSMC_Y = spaceAvgCells_Y;

    int numCellSpace_X = XcellSize/spaceAvgCells_X;
    int numCellSpace_Y = (YcellSize/nCoresY)/spaceAvgCells_Y;
    int numTotSpace    = numCellSpace_X*numCellSpace_Y;

    int cellX2 = currCell%numCellSpace_X;
    int cellY2 = currCell/numCellSpace_X;

    real rho_DSMC      = momentsDSMCProj[numMoments*currCell +  0];
    real p_DSMC        = momentsDSMCProj[numMoments*currCell +  1];
    real uX_DSMC       = momentsDSMCProj[numMoments*currCell +  2];
    real uY_DSMC       = momentsDSMCProj[numMoments*currCell +  3];
    real uZ_DSMC       = momentsDSMCProj[numMoments*currCell +  4];
    real sigmaXX_DSMC  = momentsDSMCProj[numMoments*currCell +  5];
    real sigmaXY_DSMC  = momentsDSMCProj[numMoments*currCell +  6];
    real sigmaZX_DSMC  = momentsDSMCProj[numMoments*currCell +  7];
    real sigmaYY_DSMC  = momentsDSMCProj[numMoments*currCell +  8];
    real sigmaYZ_DSMC  = momentsDSMCProj[numMoments*currCell +  9];
    real sigmaZZ_DSMC  = momentsDSMCProj[numMoments*currCell + 10];
    real qX_DSMC       = momentsDSMCProj[numMoments*currCell + 11];
    real qY_DSMC       = momentsDSMCProj[numMoments*currCell + 12];
    real qZ_DSMC       = momentsDSMCProj[numMoments*currCell + 13];

    int numP = 0;
    int totP;

    real delXInv = 1./deltaX;
    real delYInv = 1./deltaY;
    real delZInv = 1./deltaZ;

    int numCellsZ = (int)(height*delZInv);
    real nCoresYInv = 1./((real)nCoresY);
    real nCoresZInv = 1./((real)nCoresZ);
    int coreNumberY = myRank%nCoresY;
    int coreNumberZ = myRank/nCoresY;


    real heightRestricted = height*.5*nCoresZ/nCores; //only one side

    int numCellCoreZ = cellSizeZ/nCores;
    real nCoresInv = 1./((real)nCores);


    int rhoNumParticles = (int)(rho_DSMC*volume*numSpaceAvgCells_DSMC_X*numSpaceAvgCells_DSMC_Y);
    std::default_random_engine generator;

    real rhoPoisson = (real)rhoNumParticles/((double)(volume*numSpaceAvgCells_DSMC_X*numSpaceAvgCells_DSMC_Y));
    real preFactor = sqrt(2.*kT/M);

    {
        // cout<<"flag1"<<endl;
        totP = rhoNumParticles;
    }
    int spaceCellNumber = (XcellSize/spaceAvgCells_X)*cellY2 + cellX2;
    // if(spaceCellNumber == 0)
    // {
    //     particlesInSpaceAvgCells[spaceCellNumber] = totP;   //Cumulative values. Will be useful when regenerating particles
    // }
    // else
    // {
    //     particlesInSpaceAvgCells[spaceCellNumber] = particlesInSpaceAvgCells[spaceCellNumber - 1] + totP;   //Cumulative values. Will be useful when regenerating particles
    // }
    real heightByCores = height/((real)nCores);		//Don't forget to change to nCoresZ when converting to 2D parallel
    real scaleFactor = 1.5;
    real fneq;
    real theta_DSMC = p_DSMC/rho_DSMC;
    real bigM = 0.;


    {
        real a0 = rhoPoisson;
        if(fabs(a0) > bigM)
        {
            bigM = a0;
        }
        real a1 = rhoPoisson*uX_DSMC;
        if(fabs(a1) > bigM)
        {
            bigM = a1;
        }
        real a2 = rhoPoisson*uY_DSMC;
        if(fabs(a2) > bigM)
        {
            bigM = a2;
        }
        real a3 = rhoPoisson*uZ_DSMC;
        if(fabs(a3) > bigM)
        {
            bigM = a3;
        }
        real a11 = sigmaXX_DSMC - rhoPoisson*theta_DSMC;
        if(fabs(a11) > bigM)
        {
            bigM = a11;
        }
        real a22 = sigmaYY_DSMC - rhoPoisson*theta_DSMC;
        if(fabs(a22) > bigM)
        {
            bigM = a22;
        }
        real a33 = sigmaZZ_DSMC - rhoPoisson*theta_DSMC;
        if(fabs(a33) > bigM)
        {
            bigM = a33;
        }
        real a12 = sigmaXY_DSMC;
        if(fabs(a12) > bigM)
        {
            bigM = a12;
        }
        real a23 = sigmaYZ_DSMC;
        if(fabs(a23) > bigM)
        {
            bigM = a23;
        }
        real a31 = sigmaZX_DSMC;
        if(fabs(a31) > bigM)
        {
            bigM = a31;
        }
    }

    real c = 2.*bigM;

    // int cellNumberPlane = cellNumberLocal%(XcellSize*(YcellSize/nCoresY));
    // int cellX = cellNumberPlane%XcellSize;
    // int cellY = cellNumberPlane/XcellSize;
    // int cellZ = cellNumberLocal/(XcellSize*(YcellSize/nCoresY));
    int cellNumberLocal = XcellSize*cellY2*numSpaceAvgCells_DSMC_Y + cellX2*numSpaceAvgCells_DSMC_X;
    int cellZ2 = 0;
    // cout<<"C FACTOR: "<<c<<" "<<bigM<<endl;

    if(cellNumberLocal == (int)(0.5*numCellCore))
    {
        // cout<<"Print Input Moments: "<<myRank<<"\t"<<totP<<"\t"<<rho_DSMC<<" "<<p_DSMC<<" "<<uX_DSMC<<" "<<uY_DSMC<<" "<<uZ_DSMC<<" "<<sigmaXX_DSMC<<" "<<sigmaXY_DSMC<<" "<<sigmaZX_DSMC<<" "<<sigmaYY_DSMC<<" "<<sigmaYZ_DSMC<<" "<<sigmaZZ_DSMC<<" "<<qX_DSMC<<" "<<qY_DSMC<<" "<<qZ_DSMC<<" "<<rhoPoisson<<" "<<theta_DSMC<<endl;
    }

    real theta_DSMCInvSq   = 1./(((real)theta_DSMC*theta_DSMC));
    real theta_DSMCInvCube = theta_DSMCInvSq/((real)theta_DSMC);

    real uSq = uX_DSMC*uX_DSMC + uY_DSMC*uY_DSMC + uZ_DSMC*uZ_DSMC;


	  real KE = 0.;
    real uStream[3] = {};
    real stressGen[6] = {};
    real heatGen[3] = {};

    while(numP < totP)
    {
        // if(((numP%100) == 0) && (myRank == 105))
        // {
        //     cout<<"Generating..."<<coreNumberZ<<"\t"<<coreNumberY<<"\t"<<myRank<<"\t"<<cellX2<<"\t"<<cellY2<<"\t"<<cellNumberLocal<<"\t"<<numP<<"\t"<<totP<<"\t"<<c<<"\t"<<rhoPoisson<<endl;
        // }
    		real zetaX = molrng.normal()*sqrt(theta_DSMC);
        real zetaY = molrng.normal()*sqrt(theta_DSMC);
        real zetaZ = molrng.normal()*sqrt(theta_DSMC);

        real zetaSq = zetaX*zetaX + zetaY*zetaY + zetaZ*zetaZ;

         real term2 = .5*((sigmaXX_DSMC)*(zetaX*zetaX - theta_DSMC)
         + 2.*sigmaXY_DSMC*zetaX*zetaY + (sigmaYY_DSMC)*(zetaY*zetaY - theta_DSMC)
         + 2.*sigmaYZ_DSMC*zetaY*zetaZ + 2.*sigmaZX_DSMC*zetaZ*zetaX
         + (sigmaZZ_DSMC)*(zetaZ*zetaZ - theta_DSMC))*theta_DSMCInvSq;

        real term3 = .1*(qX_DSMC*zetaX + qY_DSMC*zetaY + qZ_DSMC*zetaZ)*(zetaSq - 5.0*theta_DSMC)*theta_DSMCInvCube;

        fneq = (rhoPoisson + term2 + term3);

    		real selectCriteriaRand = c*molrng.uniform();    //rho_DSMC is the amplitude of the new distribution

    		if(selectCriteriaRand <= fneq)
    		{
            int particleIndex = totParticleCurrCore;
            // cList[cellNumberLocal].cellParticleList.push_back(particleIndex);
            totParticleCurrCore++;

            Oparticle(particleIndex,1).velocity() = zetaX + uX_DSMC;
            Oparticle(particleIndex,2).velocity() = zetaY + uY_DSMC;
            Oparticle(particleIndex,3).velocity() = zetaZ + uZ_DSMC;

            Oparticle(particleIndex,1).position() = ((real)cellX2 + molrng.uniform())*deltaX*numSpaceAvgCells_DSMC_X;
            Oparticle(particleIndex,2).position() = Ybound*(real)coreNumberY*nCoresYInv + ((real)cellY2 + molrng.uniform())*deltaY*numSpaceAvgCells_DSMC_Y;

            if(coreNumberZ < (0.5*nCoresZ))
            {
                Oparticle(particleIndex,3).position() = (coreNumberZ + cellZ2 + molrng.uniform())*deltaZ;
                // cout<<"pos: "<<myRank<<"  :"<<Oparticle(particleIndex,3).position()<<endl;
            }
            else
            {
                int coreNumberNewZ = (nCoresZ - 1 - coreNumberZ);
                Oparticle(particleIndex,3).position() = height - ((coreNumberNewZ + cellZ2 + molrng.uniform())*deltaZ);
                //cout<<"pos: "<<myRank<<"  :"<<Oparticle(particleIndex,3).position()<<endl;
            }

            // int cellCoordXGuess      = (int)(Oparticle(particleIndex,1).position()*delXInv);
            // int cellCoordYGuess      = (int)((Oparticle(particleIndex,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
            //
            // int cellNumberSpaceAvg = (XcellSize/spaceAvgCells_X)*(cellCoordYGuess/spaceAvgCells_Y) + (cellCoordXGuess/spaceAvgCells_X);
            //
            // // std::cout<<"GuessValue: "<<myRank<<"\t"<<cellCoordXGuess<<"\t"<<cellX2<<"\t"<<cellCoordYGuess<<"\t"<<cellY2<<"\t"<<particleIndex<<std::endl;
            //
            // std::cout<<"GuessValue: "<<myRank<<"\t"<<cellCoordXGuess<<"\t"<<cellCoordYGuess<<"\t"<<cellNumberSpaceAvg<<"\t"<<currCell<<"\t"<<particleIndex<<std::endl;

			      KE += zetaSq;
            uStream[0] += Oparticle(particleIndex,1).velocity();
            uStream[1] += Oparticle(particleIndex,2).velocity();
            uStream[2] += Oparticle(particleIndex,3).velocity();

        		stressGen[0] += zetaX*zetaX;
        		stressGen[1] += zetaX*zetaY;
        		stressGen[2] += zetaY*zetaY;
        		stressGen[3] += zetaY*zetaZ;
        		stressGen[4] += zetaZ*zetaZ;
        		stressGen[5] += zetaZ*zetaX;

        		heatGen[0] += zetaSq*zetaX;
        	 	heatGen[1] += zetaSq*zetaY;
        		heatGen[2] += zetaSq*zetaZ;

      			numP++;
    		}
    }

	  KE /= (3.*totP);
    uStream[0] /= ((real)totP);
    uStream[1] /= ((real)totP);
    uStream[2] /= ((real)totP);

		for(int i = 0; i < 6; i++)
		{
			stressGen[i] /= ((real)volume*numSpaceAvgCells_DSMC_X*numSpaceAvgCells_DSMC_Y);
		}
		for(int i = 0; i < 3; i++)
		{
			heatGen[i] /= ((real)volume*numSpaceAvgCells_DSMC_X*numSpaceAvgCells_DSMC_Y);
		}

    // if(cellNumberLocal == (int)(0.5*numCellCore))
    {
      	//cout<<"thetaNew: "<<myRank<<"\t"<<cellNumberLocal<<" "<<KE<<" "<<totParticleCurrCore<<" "<<cList[cellNumberLocal].cellParticleList.size()<<" "<<uStream[0]<<" "<<uX_DSMC<<" "<<uStream[1]<<" "<<uY_DSMC<<" "<<uStream[2]<<" "<<uZ_DSMC<<" "<<stressGen[0]-rhoPoisson*theta_DSMC<<" "<<sigmaXX_DSMC<<" "<<stressGen[1]<<" "<<sigmaXY_DSMC<<" "<<stressGen[2]-rhoPoisson*theta_DSMC<<" "<<sigmaYY_DSMC<<" "<<stressGen[3]<<" "<<sigmaYZ_DSMC<<" "<<stressGen[4]-rhoPoisson*theta_DSMC<<" "<<sigmaZZ_DSMC<<" "<<stressGen[5]<<" "<<sigmaZX_DSMC<<" "<<heatGen[0]<<" "<<qX_DSMC<<" "<<heatGen[1]<<" "<<qY_DSMC<<" "<<heatGen[2]<<" "<<qZ_DSMC<<endl;
    }

    real uSum[DIM] = {};
    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = (totParticleCurrCore - totP); particleIndex < totParticleCurrCore; particleIndex++)
        {
            uSum[dim - 1] += Oparticle(particleIndex,dim).velocity();
        }
        uSum[dim - 1] = (uSum[dim - 1])/((real)totP);
    }

    // uSum[0] = (uSum[0])/((real)totParticleCurrCore);
    // uSum[1] = (uSum[1])/((real)totParticleCurrCore);
    // uSum[2] = (uSum[2])/((real)totParticleCurrCore);

    uSum[0] = (uX_DSMC - uSum[0]);
    uSum[1] = (uY_DSMC - uSum[1]);
    uSum[2] = (uZ_DSMC - uSum[2]);

    // cout<<"DifferenceVelocity: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<uSum[0]<<endl;
    // cout<<"DifferenceVelocity: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<uSum[1]<<endl;
    // cout<<"DifferenceVelocity: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<uSum[2]<<endl;
    /*real uSumNew[DIM] = {};
    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = (totParticleCurrCore - totP); particleIndex < totParticleCurrCore; particleIndex++)
        {
            uSumNew[dim - 1] += Oparticle(particleIndex,dim).velocity();
        }
        uSumNew[dim - 1] = (uSumNew[dim - 1])/((real)totP);
    }

    cout<<"beforevelocity: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<uSumNew[0]<<"\t"<<uX_DSMC<<"\t"<<uSumNew[1]<<"\t"<<uY_DSMC<<"\t"<<uSumNew[2]<<"\t"<<uZ_DSMC<<endl;
    */
    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = (totParticleCurrCore - totP); particleIndex < totParticleCurrCore; particleIndex++)
        {
            Oparticle(particleIndex,dim).velocity() += uSum[dim - 1];
        }
    }

    real uxAnalytical = uX_DSMC;
    real uyAnalytical = uY_DSMC;
    real uzAnalytical = uZ_DSMC;
    KE = 0.;
    for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    {
        KE += (Oparticle(particleIndex,1).velocity() - uxAnalytical)*(Oparticle(particleIndex,1).velocity() - uxAnalytical) + (Oparticle(particleIndex,2).velocity() - uyAnalytical)*(Oparticle(particleIndex,2).velocity() - uyAnalytical) + (Oparticle(particleIndex,3).velocity() - uzAnalytical)*(Oparticle(particleIndex,3).velocity() - uzAnalytical);
    }

    real velTotSqScale = sqrt(KE/(totParticleCurrCore*3.*theta_DSMC));
    real velTotSqScaleInv = 1./velTotSqScale;



    for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    {
        Oparticle(particleIndex,1).velocity() = (Oparticle(particleIndex,1).velocity() - uxAnalytical)*velTotSqScaleInv + uxAnalytical;
    }
    for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    {
        Oparticle(particleIndex,2).velocity() = (Oparticle(particleIndex,2).velocity() - uyAnalytical)*velTotSqScaleInv + uyAnalytical;
    }
    for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    {
        Oparticle(particleIndex,3).velocity() = (Oparticle(particleIndex,3).velocity() - uzAnalytical)*velTotSqScaleInv + uzAnalytical;
    }
    /*uSumNew[0] = 0.;
    uSumNew[1] = 0.;
    uSumNew[2] = 0.;
    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = (totParticleCurrCore - totP); particleIndex < totParticleCurrCore; particleIndex++)
        {
            uSumNew[dim - 1] += Oparticle(particleIndex,dim).velocity();
        }
        uSumNew[dim - 1] = (uSumNew[dim - 1])/((real)totP);
    }

    cout<<"aftervelocity: "<<myRank<<"\t"<<cellNumberLocal<<"\t"<<uSumNew[0]<<"\t"<<uX_DSMC<<"\t"<<uSumNew[1]<<"\t"<<uY_DSMC<<"\t"<<uSumNew[2]<<"\t"<<uZ_DSMC<<endl;

    if(cellNumberLocal == (0.5*numCellCore))
    {
        real tSum;
        real uDSMC[DIM] = {uX_DSMC,uY_DSMC,uZ_DSMC};
        for(int dim = 1; dim <= DIM; dim++)
        {
            tSum = 0.;
            for(int particleIndex = (totParticleCurrCore - totP); particleIndex < totParticleCurrCore; particleIndex++)
            {
                tSum += Oparticle(particleIndex,dim).velocity();
            }
            cout<<"Velocity200: "<<dim-1<<" "<<tSum/((real)totParticleCurrCore)<<" "<<uDSMC[dim - 1]<<" "<<uSum[dim - 1]<<endl;
        }
    }

    if(cellNumberLocal == 800)
        cout<<"cellNumberGen2: "<<cellNumberLocal<<" "<<cList[cellNumberLocal].cellParticleList.size()<<endl;

    real uxAnalytical = uX_DSMC;
    real uyAnalytical = uY_DSMC;
    real uzAnalytical = uZ_DSMC;
    KE = 0.;
    real uStreamNew[3] = {};
    real stressGenNew[6] = {};
    real heatGenNew[3] = {};

    for(int particleIndex = (totParticleCurrCore - totP); particleIndex < totParticleCurrCore; particleIndex++)
    {
        real zetaX = Oparticle(particleIndex,1).velocity() - uxAnalytical;
        real zetaY = Oparticle(particleIndex,2).velocity() - uyAnalytical;
        real zetaZ = Oparticle(particleIndex,3).velocity() - uzAnalytical;

        real zetaSq = (zetaX)*(zetaX) + (zetaY)*(zetaY) + (zetaZ)*(zetaZ);

        KE += zetaSq;

        uStreamNew[0] += Oparticle(particleIndex,1).velocity();
        uStreamNew[1] += Oparticle(particleIndex,2).velocity();
        uStreamNew[2] += Oparticle(particleIndex,3).velocity();

        stressGenNew[0] += zetaX*zetaX;
        stressGenNew[1] += zetaX*zetaY;
        stressGenNew[2] += zetaY*zetaY;
        stressGenNew[3] += zetaY*zetaZ;
        stressGenNew[4] += zetaZ*zetaZ;
        stressGenNew[5] += zetaZ*zetaX;

        heatGenNew[0] += zetaSq*zetaX;
        heatGenNew[1] += zetaSq*zetaY;
        heatGenNew[2] += zetaSq*zetaZ;
    }

    KE /= (3.*totP);
    uStreamNew[0] /= ((real)totP);
    uStreamNew[1] /= ((real)totP);
    uStreamNew[2] /= ((real)totP);

		for(int i = 0; i < 6; i++)
		{
			stressGenNew[i] /= ((real)volume*numSpaceAvgCells_DSMC_X*numSpaceAvgCells_DSMC_Y);
		}
		for(int i = 0; i < 3; i++)
		{
			heatGenNew[i] /= ((real)volume*numSpaceAvgCells_DSMC_X*numSpaceAvgCells_DSMC_Y);
		}

    cout<<"Distribution generated! "<<KE<<endl;

    // if(cellNumberLocal == (int)(0.5*numCellCore))
    {
      	cout<<"thetaNew2: "<<myRank<<"\t"<<cellNumberLocal<<" "<<KE<<" "<<totParticleCurrCore<<" "<<cList[cellNumberLocal].cellParticleList.size()<<" "<<uStreamNew[0]<<" "<<uX_DSMC<<" "<<uStreamNew[1]<<" "<<uY_DSMC<<" "<<uStreamNew[2]<<" "<<uZ_DSMC<<" "<<stressGenNew[0]-rhoPoisson*theta_DSMC<<" "<<sigmaXX_DSMC<<" "<<stressGenNew[1]<<" "<<sigmaXY_DSMC<<" "<<stressGenNew[2]-rhoPoisson*theta_DSMC<<" "<<sigmaYY_DSMC<<" "<<stressGenNew[3]<<" "<<sigmaYZ_DSMC<<" "<<stressGenNew[4]-rhoPoisson*theta_DSMC<<" "<<sigmaZZ_DSMC<<" "<<stressGenNew[5]<<" "<<sigmaZX_DSMC<<" "<<heatGenNew[0]<<" "<<qX_DSMC<<" "<<heatGenNew[1]<<" "<<qY_DSMC<<" "<<heatGenNew[2]<<" "<<qZ_DSMC<<endl;
    }*/

    // std::cout<<"CurrSpaceCell: "<<myRank<<"\t"<<currCell<<"\t"<<numCellCore/(spaceAvgCells_X*spaceAvgCells_Y)<<"\t"<<numCellCore<<"\t"<<spaceAvgCells_X<<"\t"<<spaceAvgCells_Y<<"\t"<<numCellSpace_X<<"\t"<<numCellSpace_Y<<std::endl;
}

// template<int DIM,typename real,typename dof>
// void hardSphere<DIM,real,dof>::projGenerateDistributionFromMoments(real* momentsDSMCProj, int myRank)
// {
//     // uX_DSMC = 0.;
//     // uZ_DSMC = 0.;
//
//     real nCoresYInv = 1./((real)nCoresY);
//     real nCoresZInv = 1./((real)nCoresZ);
//     int coreNumberY = myRank%nCoresY;
//     int coreNumberZ = myRank/nCoresY;
//
//     real heightByCores = height/((real)nCores);
//
//     real delXInv = 1./deltaX;
//     real delYInv = 1./deltaY;
//     real delZInv = 1./deltaZ;
//     int numCellsZ = (int)(height*delZInv);
//
//     real heightRestricted = nCoresZ*deltaZ; //only one side
//
//
//     for(int cellNumberLocal = 0;cellNumberLocal < (numCellCore/(spaceAvgCells_X*spaceAvgCells_Y));cellNumberLocal++)
//     {
//
//         real rho_DSMC      = momentsDSMCProj[numMoments*cellNumberLocal +  0];
//         real p_DSMC        = momentsDSMCProj[numMoments*cellNumberLocal +  1];
//         real uX_DSMC       = momentsDSMCProj[numMoments*cellNumberLocal +  2];
//         real uY_DSMC       = momentsDSMCProj[numMoments*cellNumberLocal +  3];
//         real uZ_DSMC       = momentsDSMCProj[numMoments*cellNumberLocal +  4];
//         real sigmaXX_DSMC  = momentsDSMCProj[numMoments*cellNumberLocal +  5];
//         real sigmaXY_DSMC  = momentsDSMCProj[numMoments*cellNumberLocal +  6];
//         real sigmaZX_DSMC  = momentsDSMCProj[numMoments*cellNumberLocal +  7];
//         real sigmaYY_DSMC  = momentsDSMCProj[numMoments*cellNumberLocal +  8];
//         real sigmaYZ_DSMC  = momentsDSMCProj[numMoments*cellNumberLocal +  9];
//         real sigmaZZ_DSMC  = momentsDSMCProj[numMoments*cellNumberLocal + 10];
//         real qX_DSMC       = momentsDSMCProj[numMoments*cellNumberLocal + 11];
//         real qY_DSMC       = momentsDSMCProj[numMoments*cellNumberLocal + 12];
//         real qZ_DSMC       = momentsDSMCProj[numMoments*cellNumberLocal + 13];
//
//         coreNumberZ = myRank/nCoresY;
//         int numP = 0;
//         int totP;
//         int cellNumberPlane = cellNumberLocal%(XcellSize*(YcellSize/nCoresY));
//         int cellX = cellNumberPlane%XcellSize;
//         int cellY = cellNumberPlane/XcellSize;
//         int cellZ = cellNumberLocal/(XcellSize*(YcellSize/nCoresY));
//         int numCellCoreZ = cellSizeZ/nCores;
//         real nCoresInv = 1./((real)nCores);
//
//         // real volumeCore = Xbound*Ybound*(height/(real)nCores);
//
//         int rhoNumParticles = (int)(rho_DSMC*volume);
//         // cout<<"GeneratedParticles: "<<myRank<<"\t"<<rhoNumParticles<<endl;
//         std::default_random_engine generator;
//         std::poisson_distribution<int> distribution(rhoNumParticles);
//
//         int rhoPoissonNumParticles = distribution(mrand);
//         real rhoPoisson = (real)rhoNumParticles/volume;
//         real preFactor = sqrt(2.*kT/M);
//         {
//             totP = rhoNumParticles;
//         }
//
//         real scaleFactor = 1.5;
//         real fneq;
//         real theta_DSMC = p_DSMC/rhoPoisson;
//         real bigM = 0.;
//
//         real delZInv = 1./deltaZ;
//
//
//         {
//             real a0 = rhoPoisson;
//             if(fabs(a0) > bigM)
//             {
//                 bigM = a0;
//             }
//             real a1 = rhoPoisson*uX_DSMC;
//             if(fabs(a1) > bigM)
//             {
//                 bigM = a1;
//             }
//             real a2 = rhoPoisson*uY_DSMC;
//             if(fabs(a2) > bigM)
//             {
//                 bigM = a2;
//             }
//             real a3 = rhoPoisson*uZ_DSMC;
//             if(fabs(a3) > bigM)
//             {
//                 bigM = a3;
//             }
//             real a11 = sigmaXX_DSMC - rhoPoisson*theta_DSMC;
//             if(fabs(a11) > bigM)
//             {
//                 bigM = a11;
//             }
//             real a22 = sigmaYY_DSMC - rhoPoisson*theta_DSMC;
//             if(fabs(a22) > bigM)
//             {
//                 bigM = a22;
//             }
//             real a33 = sigmaZZ_DSMC - rhoPoisson*theta_DSMC;
//             if(fabs(a33) > bigM)
//             {
//                 bigM = a33;
//             }
//             real a12 = sigmaXY_DSMC;
//             if(fabs(a12) > bigM)
//             {
//                 bigM = a12;
//             }
//             real a23 = sigmaYZ_DSMC;
//             if(fabs(a23) > bigM)
//             {
//                 bigM = a23;
//             }
//             real a31 = sigmaZX_DSMC;
//             if(fabs(a31) > bigM)
//             {
//                 bigM = a31;
//             }
//         }
//
//         real c = 2.*bigM;
//
//         // if(cellNumberLocal == 0.5*numCellCore)
//         {
//             // cout<<"Print Input Moments: "<<rho_DSMC<<" "<<p_DSMC<<" "<<uX_DSMC<<" "<<uY_DSMC<<" "<<uZ_DSMC<<" "<<sigmaXX_DSMC<<" "<<sigmaXY_DSMC<<" "<<sigmaZX_DSMC<<" "<<sigmaYY_DSMC<<" "<<sigmaYZ_DSMC<<" "<<sigmaZZ_DSMC<<" "<<qX_DSMC<<" "<<qY_DSMC<<" "<<qZ_DSMC<<" "<<rhoPoisson<<" "<<theta_DSMC<<endl;
//         }
//
//         real theta_DSMCInvSq   = 1./(((real)theta_DSMC*theta_DSMC));
//         real theta_DSMCInvCube = theta_DSMCInvSq/((real)theta_DSMC);
//
//         real uSq = uX_DSMC*uX_DSMC + uY_DSMC*uY_DSMC + uZ_DSMC*uZ_DSMC;
//
//
//     	  real KE = 0.;
//         real uStream[3] = {};
//
//         while(numP < totP)
//         {
//         		real zetaX = molrng.normal()*sqrt(kT);
//             real zetaY = molrng.normal()*sqrt(kT);
//             real zetaZ = molrng.normal()*sqrt(kT);
//
//             real zetaSq = zetaX*zetaX + zetaY*zetaY + zetaZ*zetaZ;
//
//             real term2 = .5*((sigmaXX_DSMC)*(zetaX*zetaX - theta_DSMC)
//                        + 2.*sigmaXY_DSMC*zetaX*zetaY + (sigmaYY_DSMC)*(zetaY*zetaY - theta_DSMC)
//                        + 2.*sigmaYZ_DSMC*zetaY*zetaZ + 2.*sigmaZX_DSMC*zetaZ*zetaX
//                        + (sigmaZZ_DSMC)*(zetaZ*zetaZ - theta_DSMC))*theta_DSMCInvSq;
//
//             real term3 = .1*(qX_DSMC*zetaX + qY_DSMC*zetaY + qZ_DSMC*zetaZ)*(zetaSq - 5.0*theta_DSMC)*theta_DSMCInvCube;
//
//             fneq = (rhoPoisson + term2 + term3);
//
//         		real selectCriteriaRand = c*molrng.uniform();    //rho_DSMC is the amplitude of the new distribution
//
//         		if(selectCriteriaRand <= fneq)
//         		{
//                 int particleIndex = totParticleCurrCore;
//                 totParticleCurrCore++;
//
//                 Oparticle(particleIndex,1).velocity() = zetaX + uX_DSMC;
//                 Oparticle(particleIndex,2).velocity() = zetaY + uY_DSMC;
//                 Oparticle(particleIndex,3).velocity() = zetaZ + uZ_DSMC;
//
//                 Oparticle(particleIndex,1).position() = ((real)cellX + molrng.uniform())*deltaX;
//                 Oparticle(particleIndex,2).position() = Ybound*(real)(coreNumberY*nCoresYInv) + (cellY + molrng.uniform())*deltaY;
//
//                 if(coreNumberZ < (0.5*nCoresZ))
//                 {
//                     Oparticle(particleIndex,3).position() = (coreNumberZ + cellZ + molrng.uniform())*deltaZ;
//                     //cout<<"pos: "<<myRank<<"  :"<<Oparticle(particleIndex,3).position()<<endl;
//                 }
//                 else
//                 {
//                     int coreNumberNewZ = nCoresZ - 1 - (myRank/nCoresY);
//                     Oparticle(particleIndex,3).position() = height - ((coreNumberNewZ + cellZ + molrng.uniform())*deltaZ);
//
//                     //cout<<"pos: "<<myRank<<"  :"<<Oparticle(particleIndex,3).position()<<endl;
//                 }
//
//     			      KE += zetaSq;
//                 uStream[0] += Oparticle(particleIndex,1).velocity();
//                 uStream[1] += Oparticle(particleIndex,2).velocity();
//                 uStream[2] += Oparticle(particleIndex,3).velocity();
//
//           			numP++;
//         		}
//         }
//
//     	  KE /= (3.*totP);
//         uStream[0] /= ((real)totP);
//         uStream[1] /= ((real)totP);
//         uStream[2] /= ((real)totP);
//
//     }
//
//     // cout<<"initialised particles: "<<coreNumberZ<<"\t"<<coreNumberY<<"\t"<<myRank<<"\t"<<totParticleCurrCore<<endl;
//
//
// }

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::DSMCclearAvgQts()
{
  for(int i = 0; i < numMoments*numCellCore/(spaceAvgCells_X*spaceAvgCells_Y); i++)
	{
		  moments_Average[i] = 0.;
	}

  // for(int i = 0; i < numMoments*numCellCore/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y); i++)
	// {
	// 	  moments_Average_Flux[i] = 0.;
	// }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::fgradPassSpaceAvg(int cellNumber, real rhoLocal, real* velDSMC2, real* cellStress2, real* pressure2, real* cellheatFl2, int DSMC_Count, int& momCount, int count, int coreNumber)
{
    //Do not forget to divide by the number of iterations
    //Create a local structure, not a global dynamic allocation. Helps in reducing
    //amount of memory used as well as for clearing out data at the start of every
    //new iteration

    real DSMC_CountInv = 1./((real)DSMC_Count);

    int numParticleCell = rhoLocal*volume*spaceAvgCells_X*spaceAvgCells_Y;

    moments_Average[numMoments*cellNumber    ] += rhoLocal*DSMC_CountInv;
    int cellNumberZ = cellNumber/(XcellSizePerCore*YcellSizePerCore);

    /*if((cellNumberZ == (ZcellSizePerCore - 1)) && (count == (DSMC_Count - 1)))
    {
        std::cout<<"rhoLocal: "<<moments_Average[numMoments*cellNumber    ]<<"\t"<<coreNumberZ<<"\t"<<cellNumberZ<<"\t"<<rhoLocal<<"\t"<<density<<std::endl;
    }*/

    moments_Average[numMoments*cellNumber + 1] += pressure2[cellNumber]*DSMC_CountInv;

    moments_Average[numMoments*cellNumber + 2] += velDSMC2[DIM*cellNumber + 0]*DSMC_CountInv;
    moments_Average[numMoments*cellNumber + 3] += velDSMC2[DIM*cellNumber + 1]*DSMC_CountInv;
    moments_Average[numMoments*cellNumber + 4] += velDSMC2[DIM*cellNumber + 2]*DSMC_CountInv;


    for(int stressIt = 0; stressIt < 2*DIM; stressIt++)
    {
         moments_Average[numMoments*cellNumber + 5 + stressIt] += cellStress2[6*cellNumber + stressIt]*DSMC_CountInv;
    }

    for(int heatIt = 0; heatIt < DIM; heatIt++)
    {
         moments_Average[numMoments*cellNumber + 11 + heatIt] += cellheatFl2[3*cellNumber + heatIt]*DSMC_CountInv;
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::MomentChainSpaceAvg(int coreNumber,int DSMC_Count, int& momCount, int count)
{
    fourthOrderMoment = 0.;
    real delXInv = 1./((real)deltaX);
    real delYInv = 1./((real)deltaY);

    //Averaging along the width(X) of the channel
    int numParticlesInPlane = 0;
    int numSpaceAvgCells = numCellCore/(spaceAvgCells_X*spaceAvgCells_Y);
    for(int cellNumber = 0; cellNumber < numSpaceAvgCells; cellNumber++)
    {
        velDSMC[3*cellNumber + 0] = 0.;
        velDSMC[3*cellNumber + 1] = 0.;
        velDSMC[3*cellNumber + 2] = 0.;

        pressureDSMC[cellNumber] = 0.;

        cellStress[6*cellNumber + 0] = 0.;
        cellStress[6*cellNumber + 1] = 0.;
        cellStress[6*cellNumber + 2] = 0.;
        cellStress[6*cellNumber + 3] = 0.;
        cellStress[6*cellNumber + 4] = 0.;
        cellStress[6*cellNumber + 5] = 0.;

        cellheatFl[3*cellNumber + 0] = 0.;
        cellheatFl[3*cellNumber + 1] = 0.;
        cellheatFl[3*cellNumber + 2] = 0.;
    }

    real heightNew;
    for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    {
        if(coreNumberZ < ((int)(0.5*nCoresZ)))
        {
            heightNew = Oparticle(particleIndex,3).position();
        }
        else
        {
            heightNew = height - Oparticle(particleIndex,3).position();
        }

        int cellCoordX      = (int)((Oparticle(particleIndex,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
        int cellCoordY      = (int)((Oparticle(particleIndex,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
        int cellCoordZ      = (int)(heightNew*delZInv);

        int cellNumber = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

        velDSMC[DIM*cellNumber + 0] += Oparticle(particleIndex,1).velocity();
        velDSMC[DIM*cellNumber + 1] += Oparticle(particleIndex,2).velocity();
        velDSMC[DIM*cellNumber + 2] += Oparticle(particleIndex,3).velocity();
    }

    for(int cellNumber = 0; cellNumber < numSpaceAvgCells; cellNumber++)
    {
        velDSMC[DIM*cellNumber + 0] /= (spaceList[cellNumber].cellParticleSize);
        velDSMC[DIM*cellNumber + 1] /= (spaceList[cellNumber].cellParticleSize);
        velDSMC[DIM*cellNumber + 2] /= (spaceList[cellNumber].cellParticleSize);
    }

    for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    {
        if(coreNumberZ < ((int)(0.5*nCoresZ)))
        {
            heightNew = Oparticle(particleIndex,3).position();
        }
        else
        {
            heightNew = height - Oparticle(particleIndex,3).position();
        }

        int cellCoordX      = (int)((Oparticle(particleIndex,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
        int cellCoordY      = (int)((Oparticle(particleIndex,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
        int cellCoordZ      = (int)(heightNew*delZInv);

        int cellNumber = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

        real cxRel = (Oparticle(particleIndex,1).velocity() - velDSMC[DIM*cellNumber + 0]);
        real cyRel = (Oparticle(particleIndex,2).velocity() - velDSMC[DIM*cellNumber + 1]);
        real czRel = (Oparticle(particleIndex,3).velocity() - velDSMC[DIM*cellNumber + 2]);
        real cxSq = cxRel*cxRel;
        real cySq = cyRel*cyRel;
        real czSq = czRel*czRel;

        // LB(x,y,z) = DSMC(y,z,x)
        cellStress[6*cellNumber + 0] += cxSq;                                                                                             //sigma11
        cellStress[6*cellNumber + 1] += (cxRel*cyRel);       //sigma12
        cellStress[6*cellNumber + 2] += (cxRel*czRel);       //sigma13
        cellStress[6*cellNumber + 3] += cySq;                                                                                             //sigma22
        cellStress[6*cellNumber + 4] += (cyRel*czRel);       //sigma23
        cellStress[6*cellNumber + 5] += czSq;                                                                                             //sigma33                                                                                                    //sigma33

        real vSpeedSq = cxSq + cySq + czSq;

        pressureDSMC[cellNumber] += vSpeedSq;

        cellheatFl[3*cellNumber + 0] += (vSpeedSq*cxRel);  //q1
        cellheatFl[3*cellNumber + 1] += (vSpeedSq*cyRel);  //q2
        cellheatFl[3*cellNumber + 2] += (vSpeedSq*czRel);  //q3
    }

    for(int cellNumber = 0; cellNumber < numSpaceAvgCells; cellNumber++)
    {
        real rhoLocal = (real)(spaceList[cellNumber].cellParticleSize/(volume*spaceAvgCells_X*spaceAvgCells_Y));

        for(int i = 0; i < 2*DIM; i++)
        {
            cellStress[6*cellNumber + i] /= (volume*spaceAvgCells_X*spaceAvgCells_Y);
        }

        for(int i = 0; i < DIM; i++)
        {
            cellheatFl[3*cellNumber + i] /= (volume*spaceAvgCells_X*spaceAvgCells_Y);
        }

        pressureDSMC[cellNumber] *= 1./(3.*volume*spaceAvgCells_X*spaceAvgCells_Y);

        fgradPassSpaceAvg(cellNumber,rhoLocal,velDSMC,cellStress,pressureDSMC,cellheatFl,DSMC_Count,momCount,count,coreNumber);
    }

    if((count == (DSMC_Count - 1)) && (coreNumberY == ((int)(0.5*nCoresY))) && (coreNumberX == ((int)(0.5*nCoresX))))
    {
        int cellX                     = (int)(0.5*XcellSizePerCore);
        int cellY                     = (int)(0.5*YcellSizePerCore);
        int cellNumberLocal           = ((int)(XcellSizePerCore/spaceAvgCells_X))*((int)(cellY/spaceAvgCells_Y)) + ((int)(cellX/spaceAvgCells_X));

        for(int Zcell = 0; Zcell < ZcellSizePerCore; Zcell++)
        {
          int cellNumberLocalDSMCcouple = Zcell*(XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y) + cellNumberLocal;

          std::cout<<"momentsAvg: "<<coreNumber<<"\t"<<Zcell<<"\t"<<cellNumberLocalDSMCcouple<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 0]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 1]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 2]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 3]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 4]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 5]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 6]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 7]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 8]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 9]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 10]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 11]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 12]<<"\t"<<moments_Average[numMoments*cellNumberLocalDSMCcouple + 13]<<std::endl;
        }
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::fgradPassSpaceAvgFlux(int cellNumber, real rhoLocal, real* velDSMC2, real* cellStress2, real* pressure2, real* cellheatFl2, int DSMC_Count, int& momCount, int count, int coreNumber)
{
    //Do not forget to divide by the number of iterations
    //Create a local structure, not a global dynamic allocation. Helps in reducing
    //amount of memory used as well as for clearing out data at the start of every
    //new iteration

    real DSMC_CountInv = 1./((real)DSMC_CountInv);

    int numParticleCell = rhoLocal*volume*spaceAvgCells_X*spaceAvgCells_Y;

    moments_Average_Flux[numMoments*cellNumber    ] += rhoLocal*DSMC_CountInv;

    moments_Average_Flux[numMoments*cellNumber + 1] += pressure2[cellNumber]*DSMC_CountInv;

    moments_Average_Flux[numMoments*cellNumber + 2] += velDSMC2[DIM*cellNumber + 0]*DSMC_CountInv;
    moments_Average_Flux[numMoments*cellNumber + 3] += velDSMC2[DIM*cellNumber + 1]*DSMC_CountInv;
    moments_Average_Flux[numMoments*cellNumber + 4] += velDSMC2[DIM*cellNumber + 2]*DSMC_CountInv;


    for(int stressIt = 0; stressIt < 2*DIM; stressIt++)
    {
         moments_Average_Flux[numMoments*cellNumber + 5 + stressIt] += cellStress2[6*cellNumber + stressIt]*DSMC_CountInv;
    }

    for(int heatIt = 0; heatIt < DIM; heatIt++)
    {
         moments_Average_Flux[numMoments*cellNumber + 11 + heatIt] += cellheatFl2[3*cellNumber + heatIt]*DSMC_CountInv;
    }

    if((coreNumberY == ((int)(0.5*nCoresY))) && (coreNumberX == ((int)(0.5*nCoresX))) && (cellNumber == 0))
    {
        std::cout<<"====================================================================================================================="<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<0 <<"\t"<<rhoLocal<<"\t"<<numParticleCell<<"\t"<<moments_Average_Flux[0]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<1 <<"\t"<<pressure2[cellNumber]<<"\t"<<moments_Average_Flux[1]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<2 <<"\t"<<velDSMC2[3*cellNumber + 0]<<"\t"<<moments_Average_Flux[2]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<3 <<"\t"<<velDSMC2[3*cellNumber + 1]<<"\t"<<moments_Average_Flux[3]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<4 <<"\t"<<velDSMC2[3*cellNumber + 2]<<"\t"<<moments_Average_Flux[4]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<5 <<"\t"<<cellStress2[6*cellNumber + 0]<<"\t"<<moments_Average_Flux[5]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<6 <<"\t"<<cellStress2[6*cellNumber + 1]<<"\t"<<moments_Average_Flux[6]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<7 <<"\t"<<cellStress2[6*cellNumber + 2]<<"\t"<<moments_Average_Flux[7]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<8 <<"\t"<<cellStress2[6*cellNumber + 3]<<"\t"<<moments_Average_Flux[8]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<9 <<"\t"<<cellStress2[6*cellNumber + 4]<<"\t"<<moments_Average_Flux[9]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<10<<"\t"<<cellStress2[6*cellNumber + 5]<<"\t"<<moments_Average_Flux[10]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<11<<"\t"<<cellheatFl2[3*cellNumber + 0]<<"\t"<<moments_Average_Flux[11]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<12<<"\t"<<cellheatFl2[3*cellNumber + 1]<<"\t"<<moments_Average_Flux[12]<<std::endl;
        std::cout<<"momentsNewFluxFunc: "<<count<<"\t"<<coreNumberZ<<"\t"<<13<<"\t"<<cellheatFl2[3*cellNumber + 2]<<"\t"<<moments_Average_Flux[13]<<std::endl;
        std::cout<<"====================================================================================================================="<<std::endl;
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::MomentChainSpaceAvgFlux(int coreNumber,int DSMC_Count, int& momCount, int count, int DSMC2LBCouplingCellZ)
{
    fourthOrderMoment = 0.;
    real delXInv = 1./((real)deltaX);
    real delYInv = 1./((real)deltaY);

    //Averaging along the width(X) of the channel
    int numParticlesInPlane = 0;
    const int numSpaceAvgCells = numCellCore/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y);
    for(int cellNumber = 0; cellNumber < numSpaceAvgCells; cellNumber++)
    {
        velDSMCFlux[3*cellNumber + 0] = 0.;
        velDSMCFlux[3*cellNumber + 1] = 0.;
        velDSMCFlux[3*cellNumber + 2] = 0.;

        pressureDSMCFlux[cellNumber] = 0.;

        cellStressFlux[6*cellNumber + 0] = 0.;
        cellStressFlux[6*cellNumber + 1] = 0.;
        cellStressFlux[6*cellNumber + 2] = 0.;
        cellStressFlux[6*cellNumber + 3] = 0.;
        cellStressFlux[6*cellNumber + 4] = 0.;
        cellStressFlux[6*cellNumber + 5] = 0.;

        cellheatFlFlux[3*cellNumber + 0] = 0.;
        cellheatFlFlux[3*cellNumber + 1] = 0.;
        cellheatFlFlux[3*cellNumber + 2] = 0.;
    }

    // spaceAvgFluxIndex = 10;
    real heightNew;
    int countParticleSpaceCells[numSpaceAvgCells] = {};
    for(int particleIndex = 0; particleIndex < spaceAvgFluxIndex; particleIndex++)
    {
        int cellCoordX      = (int)((spaceAvgFlux[2*DIM*particleIndex] -  Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
        int cellCoordY      = (int)((spaceAvgFlux[2*DIM*particleIndex + 1] - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);

        int cellNumber = (XcellSizePerCore/spaceAvgCells_X)*((cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

        velDSMCFlux[DIM*cellNumber + 0] += spaceAvgFlux[2*DIM*particleIndex + 3];
        velDSMCFlux[DIM*cellNumber + 1] += spaceAvgFlux[2*DIM*particleIndex + 4];
        velDSMCFlux[DIM*cellNumber + 2] += spaceAvgFlux[2*DIM*particleIndex + 5];

        countParticleSpaceCells[cellNumber]++;
    }

    for(int cellNumber = 0; cellNumber < numSpaceAvgCells; cellNumber++)
    {
        velDSMCFlux[DIM*cellNumber + 0] /= (countParticleSpaceCells[cellNumber]);
        velDSMCFlux[DIM*cellNumber + 1] /= (countParticleSpaceCells[cellNumber]);
        velDSMCFlux[DIM*cellNumber + 2] /= (countParticleSpaceCells[cellNumber]);
    }

    for(int particleIndex = 0; particleIndex < spaceAvgFluxIndex; particleIndex++)
    {
        int cellCoordX      = (int)((spaceAvgFlux[2*DIM*particleIndex] - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
        int cellCoordY      = (int)((spaceAvgFlux[2*DIM*particleIndex + 1] - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);

        int cellNumber = (XcellSizePerCore/spaceAvgCells_X)*((cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

        real cxRel = (spaceAvgFlux[2*DIM*particleIndex + 3] - velDSMCFlux[DIM*cellNumber + 0]);
        real cyRel = (spaceAvgFlux[2*DIM*particleIndex + 4] - velDSMCFlux[DIM*cellNumber + 1]);
        real czRel = (spaceAvgFlux[2*DIM*particleIndex + 5] - velDSMCFlux[DIM*cellNumber + 2]);
        real cxSq = cxRel*cxRel;
        real cySq = cyRel*cyRel;
        real czSq = czRel*czRel;

        // LB(x,y,z) = DSMC(y,z,x)
        cellStressFlux[6*cellNumber + 0] += cxSq;                                                                                             //sigma11
        cellStressFlux[6*cellNumber + 1] += (cxRel*cyRel);       //sigma12
        cellStressFlux[6*cellNumber + 2] += (cxRel*czRel);       //sigma13
        cellStressFlux[6*cellNumber + 3] += cySq;                                                                                             //sigma22
        cellStressFlux[6*cellNumber + 4] += (cyRel*czRel);       //sigma23
        cellStressFlux[6*cellNumber + 5] += czSq;                                                                                             //sigma33                                                                                                    //sigma33

        real vSpeedSq = cxSq + cySq + czSq;

        pressureDSMCFlux[cellNumber] += vSpeedSq;

        cellheatFlFlux[3*cellNumber + 0] += (vSpeedSq*cxRel);  //q1
        cellheatFlFlux[3*cellNumber + 1] += (vSpeedSq*cyRel);  //q2
        cellheatFlFlux[3*cellNumber + 2] += (vSpeedSq*czRel);  //q3
    }

    for(int cellNumber = 0; cellNumber < numSpaceAvgCells; cellNumber++)
    {
        real rhoLocal = (real)(countParticleSpaceCells[cellNumber])/(volume*spaceAvgCells_X*spaceAvgCells_Y);

        for(int i = 0; i < 2*DIM; i++)
        {
            cellStressFlux[6*cellNumber + i] /= (volume*spaceAvgCells_X*spaceAvgCells_Y);
        }

        for(int i = 0; i < DIM; i++)
        {
            cellheatFlFlux[3*cellNumber + i] /= (volume*spaceAvgCells_X*spaceAvgCells_Y);
        }

        pressureDSMCFlux[cellNumber] *= 1./(3.*volume*spaceAvgCells_X*spaceAvgCells_Y);

        fgradPassSpaceAvgFlux(cellNumber,rhoLocal,velDSMCFlux,cellStressFlux,pressureDSMCFlux,cellheatFlFlux,DSMC_Count,momCount,count,coreNumber);
    }

    // std::cout<<"finishedFluxMoments!: "<<count<<"\t"<<coreNumber<<std::endl;
    // std::cout<<"fluxquantity: "<<count<<"\t"<<coreNumber<<"\t"<<(real)(countParticleSpaceCells[0])/(volume*spaceAvgCells_X*spaceAvgCells_Y)<<std::endl;
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::totalMomentsWithFlux(real* momentsRecv, int cellNumber, int myRank, int count)
{
    int cellNumberGlobal = ((int)((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y)*(ZcellSizePerCore - 1))) + cellNumber;
    real rhoLocal = (real)(spaceList[cellNumberGlobal].cellParticleSize/(volume*spaceAvgCells_X*spaceAvgCells_Y));

    if((coreNumberY == ((int)(0.5*nCoresY))) && (coreNumberX == ((int)(0.5*nCoresX))) && (cellNumber == 0))
    {
        std::cout<<"====================================================================================================================="<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<0 <<"\t"<<momentsRecv[numMoments*cellNumber     ]<<"\t"<<rhoLocal<<"\t"<<moments_Average_Flux[0]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<1 <<"\t"<<momentsRecv[numMoments*cellNumber + 1 ]<<"\t"<<pressureDSMC[cellNumberGlobal]<<"\t"<<moments_Average_Flux[1]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<2 <<"\t"<<momentsRecv[numMoments*cellNumber + 2 ]<<"\t"<<velDSMC[3*cellNumberGlobal + 0]<<"\t"<<moments_Average_Flux[2]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<3 <<"\t"<<momentsRecv[numMoments*cellNumber + 3 ]<<"\t"<<velDSMC[3*cellNumberGlobal + 1]<<"\t"<<moments_Average_Flux[3]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<4 <<"\t"<<momentsRecv[numMoments*cellNumber + 4 ]<<"\t"<<velDSMC[3*cellNumberGlobal + 2]<<"\t"<<moments_Average_Flux[4]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<5 <<"\t"<<momentsRecv[numMoments*cellNumber + 5 ]<<"\t"<<cellStress[6*cellNumberGlobal + 0]<<"\t"<<moments_Average_Flux[5]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<6 <<"\t"<<momentsRecv[numMoments*cellNumber + 6 ]<<"\t"<<cellStress[6*cellNumberGlobal + 1]<<"\t"<<moments_Average_Flux[6]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<7 <<"\t"<<momentsRecv[numMoments*cellNumber + 7 ]<<"\t"<<cellStress[6*cellNumberGlobal + 2]<<"\t"<<moments_Average_Flux[7]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<8 <<"\t"<<momentsRecv[numMoments*cellNumber + 8 ]<<"\t"<<cellStress[6*cellNumberGlobal + 3]<<"\t"<<moments_Average_Flux[8]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<9 <<"\t"<<momentsRecv[numMoments*cellNumber + 9 ]<<"\t"<<cellStress[6*cellNumberGlobal + 4]<<"\t"<<moments_Average_Flux[9]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<10<<"\t"<<momentsRecv[numMoments*cellNumber + 10]<<"\t"<<cellStress[6*cellNumberGlobal + 5]<<"\t"<<moments_Average_Flux[10]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<11<<"\t"<<momentsRecv[numMoments*cellNumber + 11]<<"\t"<<cellheatFl[3*cellNumberGlobal + 0]<<"\t"<<moments_Average_Flux[11]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<12<<"\t"<<momentsRecv[numMoments*cellNumber + 12]<<"\t"<<cellheatFl[3*cellNumberGlobal + 1]<<"\t"<<moments_Average_Flux[12]<<std::endl;
        std::cout<<"momentsNewFlux: "<<count<<"\t"<<coreNumberZ<<"\t"<<13<<"\t"<<momentsRecv[numMoments*cellNumber + 13]<<"\t"<<cellheatFl[3*cellNumberGlobal + 2]<<"\t"<<moments_Average_Flux[13]<<std::endl;
        std::cout<<"====================================================================================================================="<<std::endl;
    }

    // moments_Average_Flux[numMoments*cellNumber     ] = 0.000001;

    moments_Average_Flux[numMoments*cellNumber     ] = momentsRecv[numMoments*cellNumber     ];//rhoLocal;
    moments_Average_Flux[numMoments*cellNumber + 1 ] = momentsRecv[numMoments*cellNumber + 1 ];//pressureDSMC[cellNumberGlobal];
    moments_Average_Flux[numMoments*cellNumber + 2 ] = momentsRecv[numMoments*cellNumber + 2 ];//velDSMC[3*cellNumberGlobal + 0];
    moments_Average_Flux[numMoments*cellNumber + 3 ] = momentsRecv[numMoments*cellNumber + 3 ];//velDSMC[3*cellNumberGlobal + 1];
    moments_Average_Flux[numMoments*cellNumber + 4 ] = momentsRecv[numMoments*cellNumber + 4 ];//velDSMC[3*cellNumberGlobal + 2];
    moments_Average_Flux[numMoments*cellNumber + 5 ] = momentsRecv[numMoments*cellNumber + 5 ];//cellStress[6*cellNumberGlobal + 0];
    moments_Average_Flux[numMoments*cellNumber + 6 ] = momentsRecv[numMoments*cellNumber + 6 ];//cellStress[6*cellNumberGlobal + 1];
    moments_Average_Flux[numMoments*cellNumber + 7 ] = momentsRecv[numMoments*cellNumber + 7 ];//cellStress[6*cellNumberGlobal + 2];
    moments_Average_Flux[numMoments*cellNumber + 8 ] = momentsRecv[numMoments*cellNumber + 8 ];//cellStress[6*cellNumberGlobal + 3];
    moments_Average_Flux[numMoments*cellNumber + 9 ] = momentsRecv[numMoments*cellNumber + 9 ];//cellStress[6*cellNumberGlobal + 4];
    moments_Average_Flux[numMoments*cellNumber + 10] = momentsRecv[numMoments*cellNumber + 10];//cellStress[6*cellNumberGlobal + 5];
    moments_Average_Flux[numMoments*cellNumber + 11] = momentsRecv[numMoments*cellNumber + 11];//cellheatFl[3*cellNumberGlobal + 0];
    moments_Average_Flux[numMoments*cellNumber + 12] = momentsRecv[numMoments*cellNumber + 12];//cellheatFl[3*cellNumberGlobal + 1];
    moments_Average_Flux[numMoments*cellNumber + 13] = momentsRecv[numMoments*cellNumber + 13];//cellheatFl[3*cellNumberGlobal + 2];
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::clearMomentsFlux(int coreNumber)
{
    int numSpaceAvgCells = numCellCore/(ZcellSizePerCore*spaceAvgCells_X*spaceAvgCells_Y);
    for(int cellNumber = 0; cellNumber < numSpaceAvgCells; cellNumber++)
    {
        for(int currMoment = 0; currMoment < numMoments; currMoment++)
        {
            moments_Average_Flux_Temp[numMoments*cellNumber + currMoment] = moments_Average_Flux[numMoments*cellNumber + currMoment];
            moments_Average_Flux[numMoments*cellNumber + currMoment] = 0.;
        }
    }
}


template<int DIM,typename real,typename dof>
real hardSphere<DIM,real,dof>::totKinEnergy()
{
    real KE = 0.;
    for(int dim = 1; dim <= DIM; dim++)
    {
        for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
        {
            KE += (Oparticle(particleIndex,dim).velocity())*(Oparticle(particleIndex,dim).velocity());
        }
    }

    return KE/(3.*totParticleCurrCore);
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::initVelocitySteadyState(int coreNumber, real velScale)
{
    //X direction
    int numP = 0;
    int totP = nParticles/nCores;
    real heightByCores = height/((real)nCores);		//Don't forget to change to nCoresZ when converting to 2D parallel
    real scaleFactor = 1.;
    real uxAnalytical = 0.;
    real rho = nParticles/((real)Xbound*Ybound*height);
    //real sqTwokbTbyMInv = 1./sqrt(3.*kT/(real)M);

    while(numP < totP)
    {
    		real vxTemp = molrng.normal()*sqrt(kT) + uxAnalytical;

    		//cout<<"numPX: "<<numP<<" "<<totP<<" "<<coreNumber<<endl;

    		real vRelSq = (vxTemp - uxAnalytical)*(vxTemp - uxAnalytical);
    		real selectCriteriaRand = rho*molrng.uniform();

    		real gaussianDist = 1.;

    		if(selectCriteriaRand <= gaussianDist)
    		{
    			Oparticle(numP,1).velocity() = vxTemp;

    			numP++;
    		}
    }

    //Y direction
    numP = 0;
    real uyAnalytical = -velScale*.5*gy*heightByCores*heightByCores*coreNumber*(coreNumber - nCores + 1)/viscosity;

    while(numP < totP)
    {
    		real vyTemp = molrng.normal()*sqrt(kT) + uyAnalytical;

    		//cout<<"numPY: "<<numP<<" "<<totP<<" "<<coreNumber<<endl;

    		real vRelSq = (vyTemp - uyAnalytical)*(vyTemp - uyAnalytical);
    		real selectCriteriaRand = rho*molrng.uniform();

    		real gaussianDist = 1.;

    		if(selectCriteriaRand <= gaussianDist)
    		{
    			Oparticle(numP,2).velocity() = vyTemp;

    			numP++;
    		}
    }

    //Z direction
    numP = 0;
    real uzAnalytical = 0.;
    // if(coreNumber == lowerCoupleDSMC)
    // {
    //     uzAnalytical = 0.1;
    // }
    // else if(coreNumber == upperCoupleDSMC)
    // {
    //     uzAnalytical = -0.1;
    // }

    while(numP < totP)
    {
    		real vzTemp = molrng.normal()*sqrt(kT) + uzAnalytical;
    		//cout<<"numPZ: "<<numP<<" "<<totP<<" "<<coreNumber<<endl;

    		real vRelSq = (vzTemp - uzAnalytical)*(vzTemp - uzAnalytical);
    		real selectCriteriaRand = rho*molrng.uniform();

    		real gaussianDist = 1.;

    		if(selectCriteriaRand <= gaussianDist)
    		{
    			Oparticle(numP,3).velocity() = vzTemp;

    			numP++;
    		}
    }

    // real KE = 0.;
    // for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    // {
    //     KE += (Oparticle(particleIndex,1).velocity() - uxAnalytical)*(Oparticle(particleIndex,1).velocity() - uxAnalytical) + (Oparticle(particleIndex,2).velocity() - uyAnalytical)*(Oparticle(particleIndex,2).velocity() - uyAnalytical) + (Oparticle(particleIndex,3).velocity() - uzAnalytical)*(Oparticle(particleIndex,3).velocity() - uzAnalytical);
    // }
    //
    // real nParticlesPerCore = nParticles/((real)nCores);
    // real velTotSqScale = sqrt(KE/(nParticlesPerCore*3.*kT));
    // real velTotSqScaleInv = 1./velTotSqScale;
    //
    //
    // for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    // {
    //     Oparticle(particleIndex,1).velocity() = (Oparticle(particleIndex,1).velocity() - uxAnalytical)*velTotSqScaleInv + uxAnalytical;
    // }
    // for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    // {
    //     Oparticle(particleIndex,2).velocity() = (Oparticle(particleIndex,2).velocity() - uyAnalytical)*velTotSqScaleInv + uyAnalytical;
    // }
    // for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    // {
    //     Oparticle(particleIndex,3).velocity() = (Oparticle(particleIndex,3).velocity() - uzAnalytical)*velTotSqScaleInv + uzAnalytical;
    // }


}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::updateVelocityBottomWallX(int particleIndex)
{
  	real kTbyM = sqrt(kT/mass);
  	Oparticle(particleIndex,1).updateVelWall(kTbyM*molrng.normal());
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::updateVelocityBottomWallY(int particleIndex)
{
  	real kTbyM = sqrt(kT/mass);
  	Oparticle(particleIndex,2).updateVelWall(kTbyM*molrng.normal() + uWall);
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::updateVelocityBottomWallZ(int particleIndex)
{
    real kTbyM = sqrt(kT/mass);
    Oparticle(particleIndex,3).updateVelWall(kTbyM*pow(-2.*log(molrng.uniform()),.5));
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::updateVelocityTopWallX(int particleIndex)
{
  	real kTbyM = sqrt(kT/mass);
  	Oparticle(particleIndex,1).updateVelWall(kTbyM*molrng.normal());
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::updateVelocityTopWallY(int particleIndex)
{
  	real kTbyM = sqrt(kT/mass);
  	Oparticle(particleIndex,2).updateVelWall(kTbyM*molrng.normal() - uWall);
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::updateVelocityTopWallZ(int particleIndex)
{
        real kTbyM = sqrt(kT/mass);
        Oparticle(particleIndex,3).updateVelWall(-kTbyM*pow(-2.*log(molrng.uniform()),.5));       //negative for top wall
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::evolveSystemDiffuse(int coreNumber, int count, int DSMC2LBCouplingCellZ)
{
    real deltatWall;

    real heightRestricted = cellSizeZ*height/((real)ZcellSize);      // For Top and Bottom slabs combined
    real heightByCoresInv = 1./(heightRestricted*nCoresZInv);

    int VECT_LENGTH   = 8;
    spaceAvgFluxIndex = 0;

    for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    {
        real posNewZ = Oparticle(particleIndex,3).position() + Oparticle(particleIndex,3).velocity()*deltat;

        if(posNewZ < 0.)
        {
            deltatWall = fabs(Oparticle(particleIndex,3).position()/Oparticle(particleIndex,3).velocity());
            Oparticle(particleIndex,1).updatePos(deltatWall);

            updateVelocityBottomWallX(particleIndex);

            Oparticle(particleIndex,1).updatePos(deltat - deltatWall);
            XperiodicBoundary(particleIndex,1);
        }
        else if(posNewZ > extent[2])
        {
            deltatWall = fabs((height - Oparticle(particleIndex,3).position())/Oparticle(particleIndex,3).velocity());
            Oparticle(particleIndex,1).updatePos(deltatWall);

            updateVelocityTopWallX(particleIndex);

            Oparticle(particleIndex,1).updatePos(deltat - deltatWall);
            XperiodicBoundary(particleIndex,1);
        }
        else
        {
            Oparticle(particleIndex,1).updatePos(deltat);
            XperiodicBoundary(particleIndex,1);
        }
    }

    real gydeltatSqByTwo = .5*gy*deltat*deltat;
    real gydeltat   = gy*deltat;

    for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    {
        real posNewZ = Oparticle(particleIndex,3).position() + Oparticle(particleIndex,3).velocity()*deltat;

        if(posNewZ < 0.)
        {
            deltatWall = fabs(Oparticle(particleIndex,3).position()/Oparticle(particleIndex,3).velocity());
            Oparticle(particleIndex,2).updatePosAcc(deltatWall,gy);

            updateVelocityBottomWallY(particleIndex);

            Oparticle(particleIndex,2).updatePosAcc(deltat - deltatWall,gy);
            Oparticle(particleIndex,2).updateVelAcc(gy*(deltat - deltatWall));
            YperiodicBoundary(particleIndex,2);
        }
        else if(posNewZ > extent[2])
        {
            deltatWall = fabs((height - Oparticle(particleIndex,3).position())/Oparticle(particleIndex,3).velocity());
            Oparticle(particleIndex,2).updatePosAcc(deltatWall,gy);

            updateVelocityTopWallY(particleIndex);

            Oparticle(particleIndex,2).updatePosAcc(deltat - deltatWall,gy);
            Oparticle(particleIndex,2).updateVelAcc(gy*(deltat - deltatWall));
            YperiodicBoundary(particleIndex,2);
        }
        else
        {
            Oparticle(particleIndex,2).updatePosAcc(deltat,gy);
            Oparticle(particleIndex,2).updateVelAcc(gy*deltat);
            YperiodicBoundary(particleIndex,2);
        }
    }

    //Z direction -- Wall Boundary
    outPlusCountY         = 0;
    outMinsCountY         = 0;
    outPlusCountPeriodicY = 0;
    outMinsCountPeriodicY = 0;
    outPlusCountX         = 0;
    outMinsCountX         = 0;
    outPlusCountPeriodicX = 0;
    outMinsCountPeriodicX = 0;

    int extraOutParticles = 0;
    countParticleCurrCore = 0;
    voidListIndex         = 0;

    for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    {
        real posNewZ = Oparticle(particleIndex,3).position() + Oparticle(particleIndex,3).velocity()*deltat;

        if(posNewZ < 0.)
        {
            deltatWall = fabs(Oparticle(particleIndex,3).position()/Oparticle(particleIndex,3).velocity());
            Oparticle(particleIndex,3).updatePos(deltatWall);

            updateVelocityBottomWallZ(particleIndex);

            Oparticle(particleIndex,3).updatePos(deltat - deltatWall);
        }
        else if(posNewZ > extent[2])
        {
            deltatWall = fabs((height - Oparticle(particleIndex,3).position())/Oparticle(particleIndex,3).velocity());
            Oparticle(particleIndex,3).updatePos(deltatWall);

            updateVelocityTopWallZ(particleIndex);

            Oparticle(particleIndex,3).updatePos(deltat - deltatWall);
        }
        else
        {
            Oparticle(particleIndex,3).position() = posNewZ;
        }

        // Packing for MPI Communication
        if(voidListSize == 0)
        {
            // std::cout<<"beforePacking: "<<coreNumber<<"\t"<<particleIndex<<"\t"<<voidListIndex<<"\t"<<voidListSize<<"\t"<<voidList[voidListIndex]<<std::endl;
            real heightNew;
            int cNumberZ;
            if(coreNumberZ < ((int)(0.5*nCoresZ)))
            {
                heightNew = Oparticle(particleIndex,3).position();
                cNumberZ  = (int)(heightNew*heightByCoresInv);
                if(cNumberZ >= ((int)(0.5*nCoresZ)))
                {
                    cNumberZ = nCoresZ;

                    // int cellCoordX      = (int)((Oparticle(particleIndex,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
                    // int cellCoordY      = (int)((Oparticle(particleIndex,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
                    // int cellNumber      = (XcellSizePerCore/spaceAvgCells_X)*(cellCoordY/spaceAvgCells_Y) + (cellCoordX/spaceAvgCells_X);
                    // // dsmcParticleFlux[cellNumber]++;
                    //
                    // if(cellNumber >= ((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y)))
                    // {
                    //     std::cout<<"cellOutofCore01: "<<coreNumber<<"\t"<<cellNumber<<"\t"<<((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y))<<std::endl;
                    // }
                }
            }
            else
            {
                heightNew = height - Oparticle(particleIndex,3).position();
                cNumberZ  = nCoresZ - 1 - (int)(heightNew*heightByCoresInv);
                if(cNumberZ < ((int)(0.5*nCoresZ)))
                {
                    cNumberZ = nCoresZ;

                    // int cellCoordX      = (int)((Oparticle(particleIndex,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
                    // int cellCoordY      = (int)((Oparticle(particleIndex,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
                    // int cellNumber      = (XcellSizePerCore/spaceAvgCells_X)*(cellCoordY/spaceAvgCells_Y) + (cellCoordX/spaceAvgCells_X);
                    // // dsmcParticleFlux[cellNumber]++;
                    //
                    // if(cellNumber >= ((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y)))
                    // {
                    //     std::cout<<"cellOutofCore02: "<<coreNumber<<"\t"<<cellNumber<<"\t"<<((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y))<<std::endl;
                    // }
                }
            }

            int cNumberX = (int)(Oparticle(particleIndex,1).position()*nCoresX*XboundInv);
            int cNumberY = (int)(Oparticle(particleIndex,2).position()*nCoresY*YboundInv);
            int cNumber  = nCoresX*(nCoresY*cNumberZ + cNumberY) + cNumberX;

            if(cNumber >= nCores)
            {
                // cout<<"not in core"<<coreNumber<<endl;
                extraOutParticles++;
            }
            else if(cNumber == coreNumber)
            {
                if((cNumber < 0) || (cNumber > nCores))
                {
                    cout<<"core beyond scope"<<endl;
                }

                real tempVel[2*DIM];

                for(int dim = 1; dim <= DIM; dim++)
                {
                    tempVel[dim - 1]                                = Oparticle(particleIndex,dim).position();
                    Oparticle(countParticleCurrCore,dim).position() = tempVel[dim - 1];
                    tempVel[DIM + dim - 1]                          = Oparticle(particleIndex,dim).velocity();
                    Oparticle(countParticleCurrCore,dim).velocity() = tempVel[DIM + dim - 1];
                }

                int cellCoordX      = (int)((Oparticle(countParticleCurrCore,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
                int cellCoordY      = (int)((Oparticle(countParticleCurrCore,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
                int cellCoordZ      = (int)(heightNew*delZInv);
                int cellNumber      = XcellSizePerCore*(YcellSizePerCore*cellCoordZ + cellCoordY) + cellCoordX;

                int cellNumberSpaceAvg = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

                if((cellNumber >= numCellCore) || (cellNumber < 0))
                {
                    cout<<"cellscream: "<<count<<"\t"<<cellNumber<<"\t"<<numCellCore<<"\t"<<coreNumber<<"\t"<<coreNumberX<<"\t"<<coreNumberY<<"\t"<<coreNumberZ<<"\t"<<particleIndex<<"\t"<<Oparticle(countParticleCurrCore,1).position()<<"\t"<<Oparticle(countParticleCurrCore,2).position()<<"\t"<<Oparticle(countParticleCurrCore,3).position()<<"\t"<<cellCoordX<<"\t"<<cellCoordY<<"\t"<<cellCoordZ<<"\t"<<heightNew<<"\t"<<deltaZ<<"\t"<<delZInv<<"\t"<<XcellSizePerCore<<"\t"<<YcellSizePerCore<<"\t"<<ZcellSizePerCore<<"\t"<<XcellSizePerCore*YcellSizePerCore*cellCoordZ<<"\t"<<XcellSizePerCore*cellCoordY<<"\t"<<endl;
                    cellNumber = numCellCore - 1;
                }

                // if((cellCoordZ == DSMC2LBCouplingCellZ) && (cellCoordZ > cellNumberZOld))
                // {
                //     spaceAvgFlux[2*DIM*spaceAvgFluxIndex    ] = tempVel[0];
                //     spaceAvgFlux[2*DIM*spaceAvgFluxIndex + 1] = tempVel[1];
                //     spaceAvgFlux[2*DIM*spaceAvgFluxIndex + 2] = tempVel[2];
                //
                //     spaceAvgFlux[2*DIM*spaceAvgFluxIndex + 3] = tempVel[3];
                //     spaceAvgFlux[2*DIM*spaceAvgFluxIndex + 4] = tempVel[4];
                //     spaceAvgFlux[2*DIM*spaceAvgFluxIndex + 5] = tempVel[5];
                //
                //     spaceAvgFluxIndex++;
                // }

                cList[cellNumber].cellParticleList[cList[cellNumber].cellParticleSize] = countParticleCurrCore;
                cList[cellNumber].cellParticleSize++;
                spaceList[cellNumberSpaceAvg].cellParticleList[spaceList[cellNumberSpaceAvg].cellParticleSize] = countParticleCurrCore;
                spaceList[cellNumberSpaceAvg].cellParticleSize++;
                Oparticle.getCellCoord(countParticleCurrCore) = cellNumber;
                countParticleCurrCore++;
            }
            else if((cNumberY - coreNumberY) == 1)
            {
                outParticlePlusY[outPlusCountY                    ] = Oparticle(particleIndex,1).position();
                outParticlePlusY[outPlusCountY + 1                ] = Oparticle(particleIndex,1).velocity();
                outParticlePlusY[outPlusCountY + 2                ] = Oparticle(particleIndex,2).position();
                outParticlePlusY[outPlusCountY + 3                ] = Oparticle(particleIndex,2).velocity();
                outParticlePlusY[outPlusCountY + 4                ] = Oparticle(particleIndex,3).position();
                outParticlePlusY[outPlusCountY + 5                ] = Oparticle(particleIndex,3).velocity();

                outPlusCountY += DIM*2;
            }
            else if((cNumberY - coreNumberY) == -1)
            {
                outParticleMinsY[outMinsCountY                    ] = Oparticle(particleIndex,1).position();
                outParticleMinsY[outMinsCountY + 1                ] = Oparticle(particleIndex,1).velocity();
                outParticleMinsY[outMinsCountY + 2                ] = Oparticle(particleIndex,2).position();
                outParticleMinsY[outMinsCountY + 3                ] = Oparticle(particleIndex,2).velocity();
                outParticleMinsY[outMinsCountY + 4                ] = Oparticle(particleIndex,3).position();
                outParticleMinsY[outMinsCountY + 5                ] = Oparticle(particleIndex,3).velocity();

                outMinsCountY += DIM*2;
            }
            else if((cNumberY == 0) && (coreNumberY == (nCoresY - 1)))
            {
                outParticlePlusPeriodicY[outPlusCountPeriodicY    ] = Oparticle(particleIndex,1).position();
                outParticlePlusPeriodicY[outPlusCountPeriodicY + 1] = Oparticle(particleIndex,1).velocity();
                outParticlePlusPeriodicY[outPlusCountPeriodicY + 2] = Oparticle(particleIndex,2).position();
                outParticlePlusPeriodicY[outPlusCountPeriodicY + 3] = Oparticle(particleIndex,2).velocity();
                outParticlePlusPeriodicY[outPlusCountPeriodicY + 4] = Oparticle(particleIndex,3).position();
                outParticlePlusPeriodicY[outPlusCountPeriodicY + 5] = Oparticle(particleIndex,3).velocity();

                outPlusCountPeriodicY += DIM*2;
            }
            else if((cNumberY == (nCoresY - 1)) && (coreNumberY == 0))
            {
                outParticleMinsPeriodicY[outMinsCountPeriodicY    ] = Oparticle(particleIndex,1).position();
                outParticleMinsPeriodicY[outMinsCountPeriodicY + 1] = Oparticle(particleIndex,1).velocity();
                outParticleMinsPeriodicY[outMinsCountPeriodicY + 2] = Oparticle(particleIndex,2).position();
                outParticleMinsPeriodicY[outMinsCountPeriodicY + 3] = Oparticle(particleIndex,2).velocity();
                outParticleMinsPeriodicY[outMinsCountPeriodicY + 4] = Oparticle(particleIndex,3).position();
                outParticleMinsPeriodicY[outMinsCountPeriodicY + 5] = Oparticle(particleIndex,3).velocity();

                outMinsCountPeriodicY += DIM*2;
            }
            else if((cNumberX - coreNumberX) == 1)
            {
                outParticlePlusX[outPlusCountX                    ] = Oparticle(particleIndex,1).position();
                outParticlePlusX[outPlusCountX + 1                ] = Oparticle(particleIndex,1).velocity();
                outParticlePlusX[outPlusCountX + 2                ] = Oparticle(particleIndex,2).position();
                outParticlePlusX[outPlusCountX + 3                ] = Oparticle(particleIndex,2).velocity();
                outParticlePlusX[outPlusCountX + 4                ] = Oparticle(particleIndex,3).position();
                outParticlePlusX[outPlusCountX + 5                ] = Oparticle(particleIndex,3).velocity();

                outPlusCountX += DIM*2;
            }
            else if((cNumberX - coreNumberX) == -1)
            {
                outParticleMinsX[outMinsCountX                    ] = Oparticle(particleIndex,1).position();
                outParticleMinsX[outMinsCountX + 1                ] = Oparticle(particleIndex,1).velocity();
                outParticleMinsX[outMinsCountX + 2                ] = Oparticle(particleIndex,2).position();
                outParticleMinsX[outMinsCountX + 3                ] = Oparticle(particleIndex,2).velocity();
                outParticleMinsX[outMinsCountX + 4                ] = Oparticle(particleIndex,3).position();
                outParticleMinsX[outMinsCountX + 5                ] = Oparticle(particleIndex,3).velocity();

                outMinsCountX += DIM*2;
            }
            else if((cNumberX == 0) && (coreNumberX == (nCoresX - 1)))
            {
                outParticlePlusPeriodicX[outPlusCountPeriodicX    ] = Oparticle(particleIndex,1).position();
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 1] = Oparticle(particleIndex,1).velocity();
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 2] = Oparticle(particleIndex,2).position();
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 3] = Oparticle(particleIndex,2).velocity();
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 4] = Oparticle(particleIndex,3).position();
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 5] = Oparticle(particleIndex,3).velocity();

                outPlusCountPeriodicX += DIM*2;
            }
            else if((cNumberX == (nCoresX - 1)) && (coreNumberX == 0))
            {
                outParticleMinsPeriodicX[outMinsCountPeriodicX    ] = Oparticle(particleIndex,1).position();
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 1] = Oparticle(particleIndex,1).velocity();
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 2] = Oparticle(particleIndex,2).position();
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 3] = Oparticle(particleIndex,2).velocity();
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 4] = Oparticle(particleIndex,3).position();
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 5] = Oparticle(particleIndex,3).velocity();

                outMinsCountPeriodicX += DIM*2;
            }
            else
            {
                std::cout<<"ExtraParticles: "<<coreNumber<<"\t"<<cNumber<<"\t"<<coreNumberX<<"\t"<<coreNumberY<<"\t"<<coreNumberZ<<"\t"<<particleIndex<<"\t"<<Oparticle(particleIndex,1).position()<<"\t"<<Oparticle(particleIndex,2).position()<<"\t"<<Oparticle(particleIndex,3).position()<<"\t"<<Oparticle(particleIndex,1).velocity()<<"\t"<<Oparticle(particleIndex,2).velocity()<<"\t"<<Oparticle(particleIndex,3).velocity()<<std::endl;
            }
        }
        else if(particleIndex != voidList[voidListIndex])
        {
            real heightNew;
            int cNumberZ;
            if(coreNumberZ < ((int)(0.5*nCoresZ)))
            {
                heightNew = Oparticle(particleIndex,3).position();
                cNumberZ  = (int)(heightNew*heightByCoresInv);
                if(cNumberZ >= ((int)(0.5*nCoresZ)))
                {
                    cNumberZ = nCoresZ;

                    // int cellCoordX      = (int)((Oparticle(particleIndex,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
                    // int cellCoordY      = (int)((Oparticle(particleIndex,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
                    // int cellNumber      = (XcellSizePerCore/spaceAvgCells_X)*(cellCoordY/spaceAvgCells_Y) + (cellCoordX/spaceAvgCells_X);
                    // // dsmcParticleFlux[cellNumber]++;
                    //
                    // if(cellNumber >= ((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y)))
                    // {
                    //     std::cout<<"cellOutofCore1: "<<coreNumber<<"\t"<<cellNumber<<"\t"<<((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y))<<"\t"<<Oparticle(particleIndex,1).position()<<"\t"<<Xbound*(real)(coreNumberX*nCoresXInv)<<"\t"<<Xbound<<std::endl;
                    // }
                }
            }
            else
            {
                heightNew = height - Oparticle(particleIndex,3).position();
                cNumberZ  = nCoresZ - 1 - (int)(heightNew*heightByCoresInv);
                if(cNumberZ < ((int)(0.5*nCoresZ)))
                {
                    cNumberZ = nCoresZ;

                    // int cellCoordX      = (int)((Oparticle(particleIndex,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
                    // int cellCoordY      = (int)((Oparticle(particleIndex,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
                    // int cellNumber      = (XcellSizePerCore/spaceAvgCells_X)*(cellCoordY/spaceAvgCells_Y) + (cellCoordX/spaceAvgCells_X);
                    // // dsmcParticleFlux[cellNumber]++;
                    //
                    // if(cellNumber >= ((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y)))
                    // {
                    //     std::cout<<"cellOutofCore2: "<<coreNumber<<"\t"<<cellNumber<<"\t"<<((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y))<<std::endl;
                    // }
                }
            }

            int cNumberX = (int)(Oparticle(particleIndex,1).position()*nCoresX*XboundInv);
            int cNumberY = (int)(Oparticle(particleIndex,2).position()*nCoresY*YboundInv);
            int cNumber  = nCoresX*(nCoresY*cNumberZ + cNumberY) + cNumberX;

            if(cNumber >= nCores)
            {
                // cout<<"not in core"<<coreNumber<<endl;
                extraOutParticles++;
            }
            else if(cNumber == coreNumber)
            {
                if((cNumber < 0) || (cNumber > nCores))
                {
                    cout<<"core beyond scope"<<endl;
                }

                real tempVel[2*DIM];

                for(int dim = 1; dim <= DIM; dim++)
                {
                    tempVel[dim - 1]                                = Oparticle(particleIndex,dim).position();
                    Oparticle(countParticleCurrCore,dim).position() = tempVel[dim - 1];
                    tempVel[DIM + dim - 1]                          = Oparticle(particleIndex,dim).velocity();
                    Oparticle(countParticleCurrCore,dim).velocity() = tempVel[DIM + dim - 1];
                }

                int cellCoordX      = (int)((Oparticle(countParticleCurrCore,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
                int cellCoordY      = (int)((Oparticle(countParticleCurrCore,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
                int cellCoordZ      = (int)(heightNew*delZInv);
                int cellNumber      = XcellSizePerCore*(YcellSizePerCore*cellCoordZ + cellCoordY) + cellCoordX;

                int cellNumberSpaceAvg = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

                if((cellNumber >= numCellCore) || (cellNumber < 0))
                {
                    cout<<"cellscream: "<<count<<"\t"<<cellNumber<<"\t"<<numCellCore<<"\t"<<coreNumber<<"\t"<<coreNumberX<<"\t"<<coreNumberY<<"\t"<<coreNumberZ<<"\t"<<particleIndex<<"\t"<<Oparticle(countParticleCurrCore,1).position()<<"\t"<<Oparticle(countParticleCurrCore,2).position()<<"\t"<<Oparticle(countParticleCurrCore,3).position()<<"\t"<<cellCoordX<<"\t"<<cellCoordY<<"\t"<<cellCoordZ<<"\t"<<heightNew<<"\t"<<deltaZ<<"\t"<<delZInv<<"\t"<<XcellSizePerCore<<"\t"<<YcellSizePerCore<<"\t"<<ZcellSizePerCore<<"\t"<<XcellSizePerCore*YcellSizePerCore*cellCoordZ<<"\t"<<XcellSizePerCore*cellCoordY<<"\t"<<endl;
                    cellNumber = numCellCore - 1;
                }

                // if((cellCoordZ == DSMC2LBCouplingCellZ) && (cellCoordZ > cellNumberZOld))
                // {
                //     spaceAvgFlux[2*DIM*spaceAvgFluxIndex    ] = tempVel[0];
                //     spaceAvgFlux[2*DIM*spaceAvgFluxIndex + 1] = tempVel[1];
                //     spaceAvgFlux[2*DIM*spaceAvgFluxIndex + 2] = tempVel[2];
                //
                //     spaceAvgFlux[2*DIM*spaceAvgFluxIndex + 3] = tempVel[3];
                //     spaceAvgFlux[2*DIM*spaceAvgFluxIndex + 4] = tempVel[4];
                //     spaceAvgFlux[2*DIM*spaceAvgFluxIndex + 5] = tempVel[5];
                //
                //     spaceAvgFluxIndex++;
                // }

                cList[cellNumber].cellParticleList[cList[cellNumber].cellParticleSize] = countParticleCurrCore;
                cList[cellNumber].cellParticleSize++;
                spaceList[cellNumberSpaceAvg].cellParticleList[spaceList[cellNumberSpaceAvg].cellParticleSize] = countParticleCurrCore;
                spaceList[cellNumberSpaceAvg].cellParticleSize++;
                Oparticle.getCellCoord(countParticleCurrCore) = cellNumber;
                countParticleCurrCore++;
            }
            else if((cNumberY - coreNumberY) == 1)
            {
                outParticlePlusY[outPlusCountY                    ] = Oparticle(particleIndex,1).position();
                outParticlePlusY[outPlusCountY + 1                ] = Oparticle(particleIndex,1).velocity();
                outParticlePlusY[outPlusCountY + 2                ] = Oparticle(particleIndex,2).position();
                outParticlePlusY[outPlusCountY + 3                ] = Oparticle(particleIndex,2).velocity();
                outParticlePlusY[outPlusCountY + 4                ] = Oparticle(particleIndex,3).position();
                outParticlePlusY[outPlusCountY + 5                ] = Oparticle(particleIndex,3).velocity();

                outPlusCountY += DIM*2;
            }
            else if((cNumberY - coreNumberY) == -1)
            {
                outParticleMinsY[outMinsCountY                    ] = Oparticle(particleIndex,1).position();
                outParticleMinsY[outMinsCountY + 1                ] = Oparticle(particleIndex,1).velocity();
                outParticleMinsY[outMinsCountY + 2                ] = Oparticle(particleIndex,2).position();
                outParticleMinsY[outMinsCountY + 3                ] = Oparticle(particleIndex,2).velocity();
                outParticleMinsY[outMinsCountY + 4                ] = Oparticle(particleIndex,3).position();
                outParticleMinsY[outMinsCountY + 5                ] = Oparticle(particleIndex,3).velocity();

                outMinsCountY += DIM*2;
            }
            else if((cNumberY == 0) && (coreNumberY == (nCoresY - 1)))
            {
                outParticlePlusPeriodicY[outPlusCountPeriodicY    ] = Oparticle(particleIndex,1).position();
                outParticlePlusPeriodicY[outPlusCountPeriodicY + 1] = Oparticle(particleIndex,1).velocity();
                outParticlePlusPeriodicY[outPlusCountPeriodicY + 2] = Oparticle(particleIndex,2).position();
                outParticlePlusPeriodicY[outPlusCountPeriodicY + 3] = Oparticle(particleIndex,2).velocity();
                outParticlePlusPeriodicY[outPlusCountPeriodicY + 4] = Oparticle(particleIndex,3).position();
                outParticlePlusPeriodicY[outPlusCountPeriodicY + 5] = Oparticle(particleIndex,3).velocity();

                outPlusCountPeriodicY += DIM*2;
            }
            else if((cNumberY == (nCoresY - 1)) && (coreNumberY == 0))
            {
                outParticleMinsPeriodicY[outMinsCountPeriodicY    ] = Oparticle(particleIndex,1).position();
                outParticleMinsPeriodicY[outMinsCountPeriodicY + 1] = Oparticle(particleIndex,1).velocity();
                outParticleMinsPeriodicY[outMinsCountPeriodicY + 2] = Oparticle(particleIndex,2).position();
                outParticleMinsPeriodicY[outMinsCountPeriodicY + 3] = Oparticle(particleIndex,2).velocity();
                outParticleMinsPeriodicY[outMinsCountPeriodicY + 4] = Oparticle(particleIndex,3).position();
                outParticleMinsPeriodicY[outMinsCountPeriodicY + 5] = Oparticle(particleIndex,3).velocity();

                outMinsCountPeriodicY += DIM*2;
            }
            else if((cNumberX - coreNumberX) == 1)
            {
                outParticlePlusX[outPlusCountX                    ] = Oparticle(particleIndex,1).position();
                outParticlePlusX[outPlusCountX + 1                ] = Oparticle(particleIndex,1).velocity();
                outParticlePlusX[outPlusCountX + 2                ] = Oparticle(particleIndex,2).position();
                outParticlePlusX[outPlusCountX + 3                ] = Oparticle(particleIndex,2).velocity();
                outParticlePlusX[outPlusCountX + 4                ] = Oparticle(particleIndex,3).position();
                outParticlePlusX[outPlusCountX + 5                ] = Oparticle(particleIndex,3).velocity();

                outPlusCountX += DIM*2;
            }
            else if((cNumberX - coreNumberX) == -1)
            {
                outParticleMinsX[outMinsCountX                    ] = Oparticle(particleIndex,1).position();
                outParticleMinsX[outMinsCountX + 1                ] = Oparticle(particleIndex,1).velocity();
                outParticleMinsX[outMinsCountX + 2                ] = Oparticle(particleIndex,2).position();
                outParticleMinsX[outMinsCountX + 3                ] = Oparticle(particleIndex,2).velocity();
                outParticleMinsX[outMinsCountX + 4                ] = Oparticle(particleIndex,3).position();
                outParticleMinsX[outMinsCountX + 5                ] = Oparticle(particleIndex,3).velocity();

                outMinsCountX += DIM*2;
            }
            else if((cNumberX == 0) && (coreNumberX == (nCoresX - 1)))
            {
                outParticlePlusPeriodicX[outPlusCountPeriodicX    ] = Oparticle(particleIndex,1).position();
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 1] = Oparticle(particleIndex,1).velocity();
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 2] = Oparticle(particleIndex,2).position();
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 3] = Oparticle(particleIndex,2).velocity();
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 4] = Oparticle(particleIndex,3).position();
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 5] = Oparticle(particleIndex,3).velocity();

                outPlusCountPeriodicX += DIM*2;
            }
            else if((cNumberX == (nCoresX - 1)) && (coreNumberX == 0))
            {
                outParticleMinsPeriodicX[outMinsCountPeriodicX    ] = Oparticle(particleIndex,1).position();
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 1] = Oparticle(particleIndex,1).velocity();
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 2] = Oparticle(particleIndex,2).position();
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 3] = Oparticle(particleIndex,2).velocity();
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 4] = Oparticle(particleIndex,3).position();
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 5] = Oparticle(particleIndex,3).velocity();

                outMinsCountPeriodicX += DIM*2;
            }
            else
            {
                std::cout<<"ExtraParticles: "<<coreNumber<<"\t"<<cNumber<<"\t"<<coreNumberX<<"\t"<<coreNumberY<<"\t"<<coreNumberZ<<"\t"<<particleIndex<<"\t"<<Oparticle(particleIndex,1).position()<<"\t"<<Oparticle(particleIndex,2).position()<<"\t"<<Oparticle(particleIndex,3).position()<<"\t"<<Oparticle(particleIndex,1).velocity()<<"\t"<<Oparticle(particleIndex,2).velocity()<<"\t"<<Oparticle(particleIndex,3).velocity()<<std::endl;
            }
        }
        else
        {
            voidListIndex++;
        }
    }

    voidListSize = 0;

    totParticleCurrCore = countParticleCurrCore;

    // density flux calculation
    // if((coreNumberX == ((int)(0.5*nCoresX))))
    // {
    //     int valueToSend;
    //     if(coreNumberY == 0)
    //     {
    //         valueToSend = outMinsCountPeriodicY/6;
    //     }
    //     else
    //     {
    //         valueToSend = outMinsCountY/6;
    //     }
    //
    //     if(coreNumber != 0)
    //     {
    //       MPI_Send(&(valueToSend),1,MPI_INT,0,6111,MPI_COMM_WORLD);
    //     }
    // }
    //
    // // if(count == 0)
    // {
    //     if(coreNumber == 0)
    //     {
    //       int recvP;
    //       int cX = ((int)(0.5*nCoresX));
    //       for(int i = 0; i <= (2*nCoresY - 1); i++)
    //       {
    //         int cNumberDensityFlux = nCoresX*(nCoresY*((int)(i/nCoresY)) + (i%nCoresY)) + cX;
    //         MPI_Recv(&recvP,1,MPI_INT,cNumberDensityFlux,6111,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //         densityFluxFileMins<<count<<"\t"<<cNumberDensityFlux<<"\t"<<recvP<<std::endl;
    //       }
    //     }
    // }
    //
    // if((coreNumberX == ((int)(0.5*nCoresX))))
    // {
    //     int valueToSend;
    //     if(coreNumberY == (nCoresY - 1))
    //     {
    //         valueToSend = outPlusCountPeriodicY/6;
    //     }
    //     else
    //     {
    //         valueToSend = outPlusCountY/6;
    //     }
    //
    //     if(coreNumber != 0)
    //     {
    //       MPI_Send(&(valueToSend),1,MPI_INT,0,6112,MPI_COMM_WORLD);
    //     }
    // }
    //
    // // if(count == 0)
    // {
    //     if(coreNumber == 0)
    //     {
    //       int recvP;
    //       int cX = ((int)(0.5*nCoresX));
    //       for(int i = 0; i <= (2*nCoresY - 1); i++)
    //       {
    //         int cNumberDensityFlux = nCoresX*(nCoresY*((int)(i/nCoresY)) + (i%nCoresY)) + cX;
    //         MPI_Recv(&recvP,1,MPI_INT,cNumberDensityFlux,6112,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //         densityFluxFilePlus<<count<<"\t"<<cNumberDensityFlux<<"\t"<<recvP<<std::endl;
    //       }
    //     }
    // }
    //
    // if((coreNumberX == ((int)(0.5*nCoresX))) && (coreNumberY == ((int)(0.5*nCoresY))))
    // {
    //     if(coreNumberZ != 0)
    //     {
    //         int cNumberDensityFlux = nCoresX*(coreNumberY) + coreNumberX;
    //         MPI_Send(densityFluxChannelHeight,ZcellSizePerCore,MPI_INT,cNumberDensityFlux,6113,MPI_COMM_WORLD);
    //     }
    //
    //     if(coreNumberZ == 0)
    //     {
    //         for(int i = 0; i < ZcellSizePerCore; i++)
    //         {
    //             densityFluxFileTops<<count<<"\t"<<i<<"\t"<<densityFluxChannelHeight[i]<<std::endl;
    //         }
    //         int recvP[ZcellSizePerCore];
    //         int cNumberDensityFlux = nCoresX*(nCoresY + coreNumberY) + coreNumberX;
    //         MPI_Recv(recvP,ZcellSizePerCore,MPI_INT,cNumberDensityFlux,6113,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //
    //         for(int i = 0; i < ZcellSizePerCore; i++)
    //         {
    //             densityFluxFileTops<<count<<"\t"<<i+ZcellSizePerCore<<"\t"<<recvP[i]<<std::endl;
    //         }
    //     }
    // }

    // for(int currCell = 0; currCell < ((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y)); currCell++)
    //     std::cout<<"dsmcParticleFlux: "<<coreNumber<<"\t"<<currCell<<"\t"<<dsmcParticleFlux[currCell]<<"\t"<<((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y))<<std::endl;

    // if((totParticleCurrCore - (countParticleCurrCore+(outPlusCountY/6)+(outMinsCountY/6)+(outPlusCountPeriodicY/6)+(outMinsCountPeriodicY/6)+(outPlusCountX/6)+(outMinsCountX/6)+(outPlusCountPeriodicX/6)+(outMinsCountPeriodicX/6)+extraOutParticles)) != 0)
    // {
    //     cout<<"Number of Particles: "<<count<<"\t"<<coreNumberZ<<"\t"<<coreNumberY<<"\t"<<coreNumber<<"\t"<<totParticleCurrCore<<"\t"<<(countParticleCurrCore+(outPlusCountY/6)+(outMinsCountY/6)+(outPlusCountPeriodicY/6)+(outMinsCountPeriodicY/6)+(outPlusCountX/6)+(outMinsCountX/6)+(outPlusCountPeriodicX/6)+(outMinsCountPeriodicX/6)+extraOutParticles)<<"\t"<<countParticleCurrCore<<"\t"<<outPlusCountY/6<<"\t"<<outMinsCountY/6<<"\t"<<outPlusCountPeriodicY/6<<"\t"<<outMinsCountPeriodicY/6<<"\t"<<outPlusCountX/6<<"\t"<<outMinsCountX/6<<"\t"<<outPlusCountPeriodicX/6<<"\t"<<outMinsCountPeriodicX/6<<"\t"<<extraOutParticles<<endl;
    // }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::updateVelocityTopWall(int particleIndex)
{
    //X and Y velocity components remain unchanged
    Oparticle(particleIndex,3).velocity() = -Oparticle(particleIndex,3).velocity();
    //reversal of normal component after specular collision with wall
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::XperiodicBoundary(int particleIndex,int dim)
{
    //periodic only in the X and Y directions
    {
        Oparticle(particleIndex,dim).setPositionPeriodic(Xbound);
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::YperiodicBoundary(int particleIndex,int dim)
{
    //periodic only in the X and Y directions
    {
        Oparticle(particleIndex,dim).setPositionPeriodic(Ybound);
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::clearCellParticles(int coreNumber)
{
    //int numCellPerCore = numCells/nCores;
    for(int currCell = 0; currCell < numCellCore; currCell++)
    {
        //cList[currCell].cellParticleList = {};
        // std::cout<<"cells: "<<coreNumber<<"\t"<<currCell<<"\t"<<numCellCore<<"\t"<<cList[currCell].cellParticleSize<<std::endl;
        cList[currCell].cellParticleSize = 0;
    }

    int numSpaceAvgCells = numCellCore/(spaceAvgCells_X*spaceAvgCells_Y);
    for(int currCell = 0; currCell < numSpaceAvgCells; currCell++)
    {
        //cList[currCell].cellParticleList = {};
        spaceList[currCell].cellParticleSize = 0;
    }

    // for(int currCell = 0; currCell < ((XcellSizePerCore/spaceAvgCells_X)*(YcellSizePerCore/spaceAvgCells_Y)); currCell++)
    // {
    //     dsmcParticleFlux[currCell] = 0;
    // }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::initPositions(int coreNumber)
{
    //numParticlePerCore[coreNumber] = nParticles/nCores;
    real nCoresInv = 1./((real)nCores);

    for(int cellNumber = 0; cellNumber < numCellCore; cellNumber++)
    {
        int cellNumberPlane = cellNumber%(XcellSize*YcellSize);
        int cellX = cellNumberPlane%XcellSize;
        int cellY = cellNumberPlane/XcellSize;
        int cellZ = cellNumber/(XcellSize*YcellSize);

        for(int particleIndex = 500*cellNumber; particleIndex < 500*(cellNumber + 1); particleIndex++)
        {
            Oparticle(particleIndex,1).position() = ((real)cellX + molrng.uniform())*deltaX;
        }
        for(int particleIndex = 500*cellNumber; particleIndex < 500*(cellNumber + 1); particleIndex++)
        {
            Oparticle(particleIndex,2).position() = ((real)cellY + molrng.uniform())*deltaY;
        }
        for(int particleIndex = 500*cellNumber; particleIndex < 500*(cellNumber + 1); particleIndex++)
        {
            Oparticle(particleIndex,3).position() = coreNumber*height*nCoresInv + ((real)cellZ + molrng.uniform())*deltaZ;
            //Height initialisation depends on rank of core
        }
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::initVelocity()
{
    real scaleFactor = 0.;  //shifting mean to scalefactor value
    real sqTwokbTbyMInv = 1./sqrt(3.*kT/(real)M);
    for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    {
        Oparticle(particleIndex,1).velocity() = molrng.normal()*sqTwokbTbyMInv + scaleFactor;
    }
    for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    {
        Oparticle(particleIndex,2).velocity() = molrng.normal()*sqTwokbTbyMInv + scaleFactor;
    }
    for(int particleIndex = 0; particleIndex < totParticleCurrCore; particleIndex++)
    {
        Oparticle(particleIndex,3).velocity() = molrng.normal()*sqTwokbTbyMInv + scaleFactor;
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::updateCellCoord(int coreNumber, int DSMC2LBCouplingCellZ)
{
    int countPlusOne              = 0;
    int countMinsOne              = 0;
    int countPlusOneRecv          = 0;
    int countMinsOneRecv          = 0;
    int countPlusOneRecv_MAX      = DIM*2*nParticlesPerCore;
    int countMinsOneRecv_MAX      = DIM*2*nParticlesPerCore;
    int partCountMinsOneY         = 0;
    int partCountPlusOneY         = 0;
    int partCountMinsOnePeriodicY = 0;
    int partCountPlusOnePeriodicY = 0;
    int partCountMinsOneX         = 0;
    int partCountPlusOneX         = 0;
    int partCountMinsOnePeriodicX = 0;
    int partCountPlusOnePeriodicX = 0;

    MPI_Request requestF1[4];
    MPI_Status statusF1[4];



    /*******************************************************************************************************/
    /***********************************SEND TO ONE Y CORE LOWER********************************************/
    /*******************************************************************************************************/
    if(coreNumberY != 0)
    {
        MPI_Send(outParticleMinsY,outMinsCountY,MPI_DOUBLE,coreNumber - nCoresX,1,MPI_COMM_WORLD);
        MPI_Send(&outMinsCountY,1,MPI_INT,coreNumber - nCoresX,3,MPI_COMM_WORLD);
    }

    if(coreNumberY != (nCoresY - 1))
    {
        MPI_Recv(partCurrRankRecvPlusOneY,countPlusOneRecv_MAX,MPI_DOUBLE,coreNumber + nCoresX,1,MPI_COMM_WORLD,statusF1);
        MPI_Recv(&partCountPlusOneY,1,MPI_INT,coreNumber + nCoresX,3,MPI_COMM_WORLD,statusF1);

        partCountPlusOneY /= 6;
    }
    /*******************************************************************************************************/
    /*******************************************************************************************************/
    /*******************************************************************************************************/



    /*******************************************************************************************************/
    /***********************************SEND TO ONE Y CORE HIGHER*******************************************/
    /*******************************************************************************************************/
    if(coreNumberY != (nCoresY - 1))
    {
        MPI_Send(outParticlePlusY,outPlusCountY,MPI_DOUBLE,coreNumber + nCoresX,0,MPI_COMM_WORLD);
        MPI_Send(&outPlusCountY,1,MPI_INT,coreNumber + nCoresX,2,MPI_COMM_WORLD);
    }

    if(coreNumberY != 0)
    {
        MPI_Recv(partCurrRankRecvMinsOneY,countMinsOneRecv_MAX,MPI_DOUBLE,coreNumber - nCoresX,0,MPI_COMM_WORLD,statusF1);
        MPI_Recv(&partCountMinsOneY,1,MPI_INT,coreNumber - nCoresX,2,MPI_COMM_WORLD,statusF1);

        partCountMinsOneY /= 6;
    }
    /*******************************************************************************************************/
    /*******************************************************************************************************/
    /*******************************************************************************************************/



    /*******************************************************************************************************/
    /**************************SEND TO ZEROETH Y CORE FROM LAST Y CORE PERIODIC*****************************/
    /*******************************************************************************************************/
    if(coreNumberY == (nCoresY - 1))
    {
        MPI_Send(outParticlePlusPeriodicY,outPlusCountPeriodicY,MPI_DOUBLE,coreNumber - (nCoresX*(nCoresY - 1)),50,MPI_COMM_WORLD);
        MPI_Send(&outPlusCountPeriodicY,1,MPI_INT,coreNumber - (nCoresX*(nCoresY - 1)),52,MPI_COMM_WORLD);
    }

    if(coreNumberY == 0)
    {
        MPI_Recv(partCurrRankRecvMinsOnePeriodicY,countMinsOneRecv_MAX,MPI_DOUBLE,coreNumber + (nCoresX*(nCoresY - 1)),50,MPI_COMM_WORLD,statusF1);
        MPI_Recv(&partCountMinsOnePeriodicY,1,MPI_INT,coreNumber + (nCoresX*(nCoresY - 1)),52,MPI_COMM_WORLD,statusF1);

        partCountMinsOnePeriodicY /= 6;
    }
    /*******************************************************************************************************/
    /*******************************************************************************************************/
    /*******************************************************************************************************/



    /*******************************************************************************************************/
    /**************************SEND TO LAST Y CORE FROM ZEROETH Y CORE PERIODIC*****************************/
    /*******************************************************************************************************/
    if(coreNumberY == 0)
    {
        MPI_Send(outParticleMinsPeriodicY,outMinsCountPeriodicY,MPI_DOUBLE,coreNumber + (nCoresX*(nCoresY - 1)),60,MPI_COMM_WORLD);
        MPI_Send(&outMinsCountPeriodicY,1,MPI_INT,coreNumber + (nCoresX*(nCoresY - 1)),63,MPI_COMM_WORLD);
    }

    if(coreNumberY == (nCoresY - 1))
    {
        MPI_Recv(partCurrRankRecvPlusOnePeriodicY,countMinsOneRecv_MAX,MPI_DOUBLE,coreNumber - (nCoresX*(nCoresY - 1)),60,MPI_COMM_WORLD,statusF1);
        MPI_Recv(&partCountPlusOnePeriodicY,1,MPI_INT,coreNumber - (nCoresX*(nCoresY - 1)),63,MPI_COMM_WORLD,statusF1);

        partCountPlusOnePeriodicY /= 6;
    }
    /*******************************************************************************************************/
    /*******************************************************************************************************/
    /*******************************************************************************************************/

    int cNumberX;

    /*******************************************************************************************************/
    /********************************UNPACK ROUTINE FOR THE Y DIRECTION*************************************/
    /*******************************************************************************************************/

    for(int recvPartPO = 0; recvPartPO < partCountPlusOneY; recvPartPO++)
    {
        real heightNew;
        if(coreNumberZ < ((int)(0.5*nCoresZ)))
        {
            heightNew = partCurrRankRecvPlusOneY[DIM*2*recvPartPO + 4];
        }
        else
        {
            heightNew = height - partCurrRankRecvPlusOneY[DIM*2*recvPartPO + 4];
        }

        cNumberX = (int)(partCurrRankRecvPlusOneY[DIM*2*recvPartPO]*nCoresX*XboundInv);

        if(cNumberX == coreNumberX)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                Oparticle(countParticleCurrCore,dim).position() = partCurrRankRecvPlusOneY[DIM*2*recvPartPO + 2*(dim - 1)];           //position at 2*(dim - 1)
                Oparticle(countParticleCurrCore,dim).velocity() = partCurrRankRecvPlusOneY[DIM*2*recvPartPO + 2*(dim - 1) + 1];       //velocity at 2*(dim - 1) + 1
            }

            int cellCoordX      = (int)((Oparticle(countParticleCurrCore,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
            int cellCoordY      = (int)((Oparticle(countParticleCurrCore,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
            int cellCoordZ      = (int)(heightNew*delZInv);
            int cellNumber      = XcellSizePerCore*(YcellSizePerCore*cellCoordZ + cellCoordY) + cellCoordX;

            int cellNumberSpaceAvg = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

            // if(cellCoordZ == DSMC2LBCouplingCellZ)
            // {
            //     spaceAvgFlux[DIM*]
            // }

          	if((cellNumber >= numCellCore) || (cellNumber < 0))
          	{
          		cout<<"cellscream2: "<<cellNumber<<"\t"<<numCellCore<<"\t"<<coreNumber<<endl;
          		cellNumber = numCellCore - 1;
          	}
            cList[cellNumber].cellParticleList[cList[cellNumber].cellParticleSize] = countParticleCurrCore;
            cList[cellNumber].cellParticleSize++;
            spaceList[cellNumberSpaceAvg].cellParticleList[spaceList[cellNumberSpaceAvg].cellParticleSize] = countParticleCurrCore;
            spaceList[cellNumberSpaceAvg].cellParticleSize++;
            Oparticle.getCellCoord(countParticleCurrCore) = cellNumber;
            countParticleCurrCore++;
        }
        else if((cNumberX - coreNumberX) == 1)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticlePlusX[outPlusCountX    ] = partCurrRankRecvPlusOneY[DIM*2*recvPartPO + 2*(dim - 1)];
                outParticlePlusX[outPlusCountX + 1] = partCurrRankRecvPlusOneY[DIM*2*recvPartPO + 2*(dim - 1) + 1];

                outPlusCountX += 2;
            }
        }
        else if((cNumberX - coreNumberX) == -1)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticleMinsX[outMinsCountX    ] = partCurrRankRecvPlusOneY[DIM*2*recvPartPO + 2*(dim - 1)];
                outParticleMinsX[outMinsCountX + 1] = partCurrRankRecvPlusOneY[DIM*2*recvPartPO + 2*(dim - 1) + 1];

                outMinsCountX += 2;
            }
        }
        else if((cNumberX == 0) && (coreNumberX == (nCoresX - 1)))
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticlePlusPeriodicX[outPlusCountPeriodicX    ] = partCurrRankRecvPlusOneY[DIM*2*recvPartPO + 2*(dim - 1)];
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 1] = partCurrRankRecvPlusOneY[DIM*2*recvPartPO + 2*(dim - 1) + 1];

                outPlusCountPeriodicX += 2;
            }
        }
        else if((cNumberX == (nCoresX - 1)) && (coreNumberX == 0))
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticleMinsPeriodicX[outMinsCountPeriodicX    ] = partCurrRankRecvPlusOneY[DIM*2*recvPartPO + 2*(dim - 1)];
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 1] = partCurrRankRecvPlusOneY[DIM*2*recvPartPO + 2*(dim - 1) + 1];

                outMinsCountPeriodicX += 2;
            }
        }
    }

    for(int recvPartMO = 0; recvPartMO < partCountMinsOneY; recvPartMO++)
    {
        real heightNew;
        if(coreNumberZ < ((int)(0.5*nCoresZ)))
        {
            heightNew = partCurrRankRecvMinsOneY[DIM*2*recvPartMO + 4];
        }
        else
        {
            heightNew = height - partCurrRankRecvMinsOneY[DIM*2*recvPartMO + 4];
        }

        cNumberX = (int)(partCurrRankRecvMinsOneY[DIM*2*recvPartMO]*nCoresX*XboundInv);

        if(cNumberX == coreNumberX)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                Oparticle(countParticleCurrCore,dim).position() = partCurrRankRecvMinsOneY[DIM*2*recvPartMO + 2*(dim - 1)];      //position at 2*(dim - 1)
                Oparticle(countParticleCurrCore,dim).velocity() = partCurrRankRecvMinsOneY[DIM*2*recvPartMO + 2*(dim - 1) + 1];      //velocity at 2*(dim - 1) + 1
            }

            int cellCoordX      = (int)((Oparticle(countParticleCurrCore,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
            int cellCoordY      = (int)((Oparticle(countParticleCurrCore,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
            int cellCoordZ      = (int)(heightNew*delZInv);
            int cellNumber      = XcellSizePerCore*(YcellSizePerCore*cellCoordZ + cellCoordY) + cellCoordX;

            int cellNumberSpaceAvg = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

          	if((cellNumber >= numCellCore) || (cellNumber < 0))
          	{
            		cout<<"cellscream3: "<<Oparticle(countParticleCurrCore,1).position()<<"\t"<<Oparticle(countParticleCurrCore,2).position()<<"\t"<<Oparticle(countParticleCurrCore,3).position()<<"\t"<<cellNumber<<"\t"<<numCellCore<<"\t"<<coreNumber<<endl;
            		cellNumber = numCellCore - 1;
          	}
          	cList[cellNumber].cellParticleList[cList[cellNumber].cellParticleSize] = countParticleCurrCore;
            cList[cellNumber].cellParticleSize++;
            spaceList[cellNumberSpaceAvg].cellParticleList[spaceList[cellNumberSpaceAvg].cellParticleSize] = countParticleCurrCore;
            spaceList[cellNumberSpaceAvg].cellParticleSize++;
            Oparticle.getCellCoord(countParticleCurrCore) = cellNumber;
            countParticleCurrCore++;
        }
        else if((cNumberX - coreNumberX) == 1)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticlePlusX[outPlusCountX    ] = partCurrRankRecvMinsOneY[DIM*2*recvPartMO + 2*(dim - 1)];
                outParticlePlusX[outPlusCountX + 1] = partCurrRankRecvMinsOneY[DIM*2*recvPartMO + 2*(dim - 1) + 1];

                outPlusCountX += 2;
            }
        }
        else if((cNumberX - coreNumberX) == -1)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticleMinsX[outMinsCountX    ] = partCurrRankRecvMinsOneY[DIM*2*recvPartMO + 2*(dim - 1)];
                outParticleMinsX[outMinsCountX + 1] = partCurrRankRecvMinsOneY[DIM*2*recvPartMO + 2*(dim - 1) + 1];

                outMinsCountX += 2;
            }
        }
        else if((cNumberX == 0) && (coreNumberX == (nCoresX - 1)))
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticlePlusPeriodicX[outPlusCountPeriodicX    ] = partCurrRankRecvMinsOneY[DIM*2*recvPartMO + 2*(dim - 1)];
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 1] = partCurrRankRecvMinsOneY[DIM*2*recvPartMO + 2*(dim - 1) + 1];

                outPlusCountPeriodicX += 2;
            }
        }
        else if((cNumberX == (nCoresX - 1)) && (coreNumberX == 0))
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticleMinsPeriodicX[outMinsCountPeriodicX    ] = partCurrRankRecvMinsOneY[DIM*2*recvPartMO + 2*(dim - 1)];
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 1] = partCurrRankRecvMinsOneY[DIM*2*recvPartMO + 2*(dim - 1) + 1];

                outMinsCountPeriodicX += 2;
            }
        }
    }

    for(int recvPartPO = 0; recvPartPO < partCountPlusOnePeriodicY; recvPartPO++)
    {
        real heightNew;
        if(coreNumberZ < ((int)(0.5*nCoresZ)))
        {
            heightNew = partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO + 4];
        }
        else
        {
            heightNew = height - partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO + 4];
        }

        cNumberX = (int)(partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO]*nCoresX*XboundInv);

        if(cNumberX == coreNumberX)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                Oparticle(countParticleCurrCore,dim).position() = partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO + 2*(dim - 1)];      //position at 2*(dim - 1)
                Oparticle(countParticleCurrCore,dim).velocity() = partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO + 2*(dim - 1) + 1];      //velocity at 2*(dim - 1) + 1
            }

            int cellCoordX      = (int)((Oparticle(countParticleCurrCore,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
            int cellCoordY      = (int)((Oparticle(countParticleCurrCore,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
            int cellCoordZ      = (int)(heightNew*delZInv);
            int cellNumber      = XcellSizePerCore*(YcellSizePerCore*cellCoordZ + cellCoordY) + cellCoordX;

            int cellNumberSpaceAvg = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

          	if((cellNumber >= numCellCore) || (cellNumber < 0))
          	{
          		cout<<"cellscream2: "<<cellNumber<<"\t"<<numCellCore<<"\t"<<coreNumber<<endl;
          		cellNumber = numCellCore - 1;
          	}
            cList[cellNumber].cellParticleList[cList[cellNumber].cellParticleSize] = countParticleCurrCore;
            cList[cellNumber].cellParticleSize++;
            spaceList[cellNumberSpaceAvg].cellParticleList[spaceList[cellNumberSpaceAvg].cellParticleSize] = countParticleCurrCore;
            spaceList[cellNumberSpaceAvg].cellParticleSize++;
            Oparticle.getCellCoord(countParticleCurrCore) = cellNumber;
            countParticleCurrCore++;
        }
        else if((cNumberX - coreNumberX) == 1)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticlePlusX[outPlusCountX    ] = partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO + 2*(dim - 1)];
                outParticlePlusX[outPlusCountX + 1] = partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO + 2*(dim - 1) + 1];

                outPlusCountX += 2;
            }
        }
        else if((cNumberX - coreNumberX) == -1)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticleMinsX[outMinsCountX    ] = partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO + 2*(dim - 1)];
                outParticleMinsX[outMinsCountX + 1] = partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO + 2*(dim - 1) + 1];

                outMinsCountX += 2;
            }
        }
        else if((cNumberX == 0) && (coreNumberX == (nCoresX - 1)))
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticlePlusPeriodicX[outPlusCountPeriodicX    ] = partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO + 2*(dim - 1)];
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 1] = partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO + 2*(dim - 1) + 1];

                outPlusCountPeriodicX += 2;
            }
        }
        else if((cNumberX == (nCoresX - 1)) && (coreNumberX == 0))
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticleMinsPeriodicX[outMinsCountPeriodicX    ] = partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO + 2*(dim - 1)];
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 1] = partCurrRankRecvPlusOnePeriodicY[DIM*2*recvPartPO + 2*(dim - 1) + 1];

                outMinsCountPeriodicX += 2;
            }
        }
    }

    for(int recvPartMO = 0; recvPartMO < partCountMinsOnePeriodicY; recvPartMO++)
    {
        real heightNew;
        if(coreNumberZ < ((int)(0.5*nCoresZ)))
        {
            heightNew = partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO + 4];
        }
        else
        {
            heightNew = height - partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO + 4];
        }

        cNumberX = (int)(partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO]*nCoresX*XboundInv);

        if(cNumberX == coreNumberX)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                Oparticle(countParticleCurrCore,dim).position() = partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO + 2*(dim - 1)];      //position at 2*(dim - 1)
                Oparticle(countParticleCurrCore,dim).velocity() = partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO + 2*(dim - 1) + 1];      //velocity at 2*(dim - 1) + 1
            }

            int cellCoordX      = (int)((Oparticle(countParticleCurrCore,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
            int cellCoordY      = (int)((Oparticle(countParticleCurrCore,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
            int cellCoordZ      = (int)(heightNew*delZInv);
            int cellNumber      = XcellSizePerCore*(YcellSizePerCore*cellCoordZ + cellCoordY) + cellCoordX;

            int cellNumberSpaceAvg = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

          	if((cellNumber >= numCellCore) || (cellNumber < 0))
          	{
            		cout<<"cellscream3: "<<Oparticle(countParticleCurrCore,1).position()<<"\t"<<Oparticle(countParticleCurrCore,2).position()<<"\t"<<Oparticle(countParticleCurrCore,3).position()<<"\t"<<cellNumber<<"\t"<<numCellCore<<"\t"<<coreNumber<<endl;
            		cellNumber = numCellCore - 1;
          	}
          	cList[cellNumber].cellParticleList[cList[cellNumber].cellParticleSize] = countParticleCurrCore;
            cList[cellNumber].cellParticleSize++;
            spaceList[cellNumberSpaceAvg].cellParticleList[spaceList[cellNumberSpaceAvg].cellParticleSize] = countParticleCurrCore;
            spaceList[cellNumberSpaceAvg].cellParticleSize++;
            Oparticle.getCellCoord(countParticleCurrCore) = cellNumber;
            countParticleCurrCore++;
        }
        else if((cNumberX - coreNumberX) == 1)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticlePlusX[outPlusCountX    ] = partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO + 2*(dim - 1)];
                outParticlePlusX[outPlusCountX + 1] = partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO + 2*(dim - 1) + 1];

                outPlusCountX += 2;
            }
        }
        else if((cNumberX - coreNumberX) == -1)
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticleMinsX[outMinsCountX    ] = partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO + 2*(dim - 1)];
                outParticleMinsX[outMinsCountX + 1] = partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO + 2*(dim - 1) + 1];

                outMinsCountX += 2;
            }
        }
        else if((cNumberX == 0) && (coreNumberX == (nCoresX - 1)))
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticlePlusPeriodicX[outPlusCountPeriodicX    ] = partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO + 2*(dim - 1)];
                outParticlePlusPeriodicX[outPlusCountPeriodicX + 1] = partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO + 2*(dim - 1) + 1];

                outPlusCountPeriodicX += 2;
            }
        }
        else if((cNumberX == (nCoresX - 1)) && (coreNumberX == 0))
        {
            for(int dim = 1; dim <= DIM; dim++)
            {
                outParticleMinsPeriodicX[outMinsCountPeriodicX    ] = partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO + 2*(dim - 1)];
                outParticleMinsPeriodicX[outMinsCountPeriodicX + 1] = partCurrRankRecvMinsOnePeriodicY[DIM*2*recvPartMO + 2*(dim - 1) + 1];

                outMinsCountPeriodicX += 2;
            }
        }
    }

    /*******************************************************************************************************/
    /*******************************************************************************************************/
    /*******************************************************************************************************/



    /*******************************************************************************************************/
    /***********************************SEND TO ONE X CORE LOWER********************************************/
    /*******************************************************************************************************/
    if(coreNumberX != 0)
    {
        MPI_Send(outParticleMinsX,outMinsCountX,MPI_DOUBLE,coreNumber - 1,1,MPI_COMM_WORLD);
        MPI_Send(&outMinsCountX,1,MPI_INT,coreNumber - 1,3,MPI_COMM_WORLD);
    }

    if(coreNumberX != (nCoresX - 1))
    {
        MPI_Recv(partCurrRankRecvPlusOneX,countPlusOneRecv_MAX,MPI_DOUBLE,coreNumber + 1,1,MPI_COMM_WORLD,statusF1);
        MPI_Recv(&partCountPlusOneX,1,MPI_INT,coreNumber + 1,3,MPI_COMM_WORLD,statusF1);

        partCountPlusOneX /= 6;
    }
    /*******************************************************************************************************/
    /*******************************************************************************************************/
    /*******************************************************************************************************/



    /*******************************************************************************************************/
    /***********************************SEND TO ONE X CORE HIGHER*******************************************/
    /*******************************************************************************************************/
    if(coreNumberX != (nCoresX - 1))
    {
        MPI_Send(outParticlePlusX,outPlusCountX,MPI_DOUBLE,coreNumber + 1,0,MPI_COMM_WORLD);
        MPI_Send(&outPlusCountX,1,MPI_INT,coreNumber + 1,2,MPI_COMM_WORLD);
    }

    if(coreNumberX != 0)
    {
        MPI_Recv(partCurrRankRecvMinsOneX,countMinsOneRecv_MAX,MPI_DOUBLE,coreNumber - 1,0,MPI_COMM_WORLD,statusF1);
        MPI_Recv(&partCountMinsOneX,1,MPI_INT,coreNumber - 1,2,MPI_COMM_WORLD,statusF1);

        partCountMinsOneX /= 6;
    }
    /*******************************************************************************************************/
    /*******************************************************************************************************/
    /*******************************************************************************************************/



    /*******************************************************************************************************/
    /**************************SEND TO ZEROETH X CORE FROM LAST X CORE PERIODIC*****************************/
    /*******************************************************************************************************/
    if(coreNumberX == (nCoresX - 1))
    {
        MPI_Send(outParticlePlusPeriodicX,outPlusCountPeriodicX,MPI_DOUBLE,coreNumber - (nCoresX - 1),50,MPI_COMM_WORLD);
        MPI_Send(&outPlusCountPeriodicX,1,MPI_INT,coreNumber - (nCoresX - 1),52,MPI_COMM_WORLD);
    }

    if(coreNumberX == 0)
    {
        MPI_Recv(partCurrRankRecvMinsOnePeriodicX,countMinsOneRecv_MAX,MPI_DOUBLE,coreNumber + (nCoresX - 1),50,MPI_COMM_WORLD,statusF1);
        MPI_Recv(&partCountMinsOnePeriodicX,1,MPI_INT,coreNumber + (nCoresX - 1),52,MPI_COMM_WORLD,statusF1);

        partCountMinsOnePeriodicX /= 6;
    }
    /*******************************************************************************************************/
    /*******************************************************************************************************/
    /*******************************************************************************************************/



    /*******************************************************************************************************/
    /**************************SEND TO LAST X CORE FROM ZEROETH X CORE PERIODIC*****************************/
    /*******************************************************************************************************/
    if(coreNumberX == 0)
    {
        MPI_Send(outParticleMinsPeriodicX,outMinsCountPeriodicX,MPI_DOUBLE,coreNumber + (nCoresX - 1),60,MPI_COMM_WORLD);
        MPI_Send(&outMinsCountPeriodicX,1,MPI_INT,coreNumber + (nCoresX - 1),63,MPI_COMM_WORLD);
    }

    if(coreNumberX == (nCoresX - 1))
    {
        MPI_Recv(partCurrRankRecvPlusOnePeriodicX,countMinsOneRecv_MAX,MPI_DOUBLE,coreNumber - (nCoresX - 1),60,MPI_COMM_WORLD,statusF1);
        MPI_Recv(&partCountPlusOnePeriodicX,1,MPI_INT,coreNumber - (nCoresX - 1),63,MPI_COMM_WORLD,statusF1);

        partCountPlusOnePeriodicX /= 6;
    }
    /*******************************************************************************************************/
    /*******************************************************************************************************/
    /*******************************************************************************************************/




    /*******************************************************************************************************/
    /********************************UNPACK ROUTINE FOR THE X DIRECTION*************************************/
    /*******************************************************************************************************/

    for(int recvPartPO = 0; recvPartPO < partCountPlusOneX; recvPartPO++)
    {
        real heightNew;
        if(coreNumberZ < ((int)(0.5*nCoresZ)))
        {
            heightNew = partCurrRankRecvPlusOneX[DIM*2*recvPartPO + 4];
        }
        else
        {
            heightNew = height - partCurrRankRecvPlusOneX[DIM*2*recvPartPO + 4];
        }

        for(int dim = 1; dim <= DIM; dim++)
        {
            Oparticle(countParticleCurrCore,dim).position() = partCurrRankRecvPlusOneX[DIM*2*recvPartPO + 2*(dim - 1)];      //position at 2*(dim - 1)
            Oparticle(countParticleCurrCore,dim).velocity() = partCurrRankRecvPlusOneX[DIM*2*recvPartPO + 2*(dim - 1) + 1];       //velocity at 2*(dim - 1) + 1
        }

        int cellCoordX      = (int)((Oparticle(countParticleCurrCore,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
        int cellCoordY      = (int)((Oparticle(countParticleCurrCore,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
        int cellCoordZ      = (int)(heightNew*delZInv);
        int cellNumber      = XcellSizePerCore*(YcellSizePerCore*cellCoordZ + cellCoordY) + cellCoordX;

        int cellNumberSpaceAvg = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

      	if((cellNumber >= numCellCore) || (cellNumber < 0))
      	{
      		cout<<"cellscream2: "<<cellNumber<<"\t"<<numCellCore<<"\t"<<coreNumber<<endl;
      		cellNumber = numCellCore - 1;
      	}
        cList[cellNumber].cellParticleList[cList[cellNumber].cellParticleSize] = countParticleCurrCore;
        cList[cellNumber].cellParticleSize++;
        spaceList[cellNumberSpaceAvg].cellParticleList[spaceList[cellNumberSpaceAvg].cellParticleSize] = countParticleCurrCore;
        spaceList[cellNumberSpaceAvg].cellParticleSize++;
        Oparticle.getCellCoord(countParticleCurrCore) = cellNumber;
        countParticleCurrCore++;
    }

    for(int recvPartMO = 0; recvPartMO < partCountMinsOneX; recvPartMO++)
    {
        real heightNew;
        if(coreNumberZ < ((int)(0.5*nCoresZ)))
        {
            heightNew = partCurrRankRecvMinsOneX[DIM*2*recvPartMO + 4];
        }
        else
        {
            heightNew = height - partCurrRankRecvMinsOneX[DIM*2*recvPartMO + 4];
        }

        for(int dim = 1; dim <= DIM; dim++)
        {
            Oparticle(countParticleCurrCore,dim).position() = partCurrRankRecvMinsOneX[DIM*2*recvPartMO + 2*(dim - 1)];      //position at 2*(dim - 1)
            Oparticle(countParticleCurrCore,dim).velocity() = partCurrRankRecvMinsOneX[DIM*2*recvPartMO + 2*(dim - 1) + 1];      //velocity at 2*(dim - 1) + 1
        }

        int cellCoordX      = (int)((Oparticle(countParticleCurrCore,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
        int cellCoordY      = (int)((Oparticle(countParticleCurrCore,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
        int cellCoordZ      = (int)(heightNew*delZInv);
        int cellNumber      = XcellSizePerCore*(YcellSizePerCore*cellCoordZ + cellCoordY) + cellCoordX;

        int cellNumberSpaceAvg = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

      	if((cellNumber >= numCellCore) || (cellNumber < 0))
      	{
        		cout<<"cellscream3: "<<Oparticle(countParticleCurrCore,1).position()<<"\t"<<Oparticle(countParticleCurrCore,2).position()<<"\t"<<Oparticle(countParticleCurrCore,3).position()<<"\t"<<cellNumber<<"\t"<<numCellCore<<"\t"<<coreNumber<<endl;
        		cellNumber = numCellCore - 1;
      	}
      	cList[cellNumber].cellParticleList[cList[cellNumber].cellParticleSize] = countParticleCurrCore;
        cList[cellNumber].cellParticleSize++;
        spaceList[cellNumberSpaceAvg].cellParticleList[spaceList[cellNumberSpaceAvg].cellParticleSize] = countParticleCurrCore;
        spaceList[cellNumberSpaceAvg].cellParticleSize++;
        Oparticle.getCellCoord(countParticleCurrCore) = cellNumber;
        countParticleCurrCore++;
    }

    for(int recvPartPO = 0; recvPartPO < partCountPlusOnePeriodicX; recvPartPO++)
    {
        real heightNew;
        if(coreNumberZ < ((int)(0.5*nCoresZ)))
        {
            heightNew = partCurrRankRecvPlusOnePeriodicX[DIM*2*recvPartPO + 4];
        }
        else
        {
            heightNew = height - partCurrRankRecvPlusOnePeriodicX[DIM*2*recvPartPO + 4];
        }

        for(int dim = 1; dim <= DIM; dim++)
        {
            Oparticle(countParticleCurrCore,dim).position() = partCurrRankRecvPlusOnePeriodicX[DIM*2*recvPartPO + 2*(dim - 1)];      //position at 2*(dim - 1)
            Oparticle(countParticleCurrCore,dim).velocity() = partCurrRankRecvPlusOnePeriodicX[DIM*2*recvPartPO + 2*(dim - 1) + 1];      //velocity at 2*(dim - 1) + 1
        }

        int cellCoordX      = (int)((Oparticle(countParticleCurrCore,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
        int cellCoordY      = (int)((Oparticle(countParticleCurrCore,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
        int cellCoordZ      = (int)(heightNew*delZInv);
        int cellNumber      = XcellSizePerCore*(YcellSizePerCore*cellCoordZ + cellCoordY) + cellCoordX;

        int cellNumberSpaceAvg = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

      	if((cellNumber >= numCellCore) || (cellNumber < 0))
      	{
      		cout<<"cellscream2: "<<cellNumber<<"\t"<<numCellCore<<"\t"<<coreNumber<<endl;
      		cellNumber = numCellCore - 1;
      	}
        cList[cellNumber].cellParticleList[cList[cellNumber].cellParticleSize] = countParticleCurrCore;
        cList[cellNumber].cellParticleSize++;
        spaceList[cellNumberSpaceAvg].cellParticleList[spaceList[cellNumberSpaceAvg].cellParticleSize] = countParticleCurrCore;
        spaceList[cellNumberSpaceAvg].cellParticleSize++;
        Oparticle.getCellCoord(countParticleCurrCore) = cellNumber;
        countParticleCurrCore++;
    }

    for(int recvPartMO = 0; recvPartMO < partCountMinsOnePeriodicX; recvPartMO++)
    {
        real heightNew;
        if(coreNumberZ < ((int)(0.5*nCoresZ)))
        {
            heightNew = partCurrRankRecvMinsOnePeriodicX[DIM*2*recvPartMO + 4];
        }
        else
        {
            heightNew = height - partCurrRankRecvMinsOnePeriodicX[DIM*2*recvPartMO + 4];
        }

        for(int dim = 1; dim <= DIM; dim++)
        {
            Oparticle(countParticleCurrCore,dim).position() = partCurrRankRecvMinsOnePeriodicX[DIM*2*recvPartMO + 2*(dim - 1)];      //position at 2*(dim - 1)
            Oparticle(countParticleCurrCore,dim).velocity() = partCurrRankRecvMinsOnePeriodicX[DIM*2*recvPartMO + 2*(dim - 1) + 1];      //velocity at 2*(dim - 1) + 1
        }

        int cellCoordX      = (int)((Oparticle(countParticleCurrCore,1).position() - Xbound*(real)(coreNumberX*nCoresXInv))*delXInv);
        int cellCoordY      = (int)((Oparticle(countParticleCurrCore,2).position() - Ybound*(real)(coreNumberY*nCoresYInv))*delYInv);
        int cellCoordZ      = (int)(heightNew*delZInv);
        int cellNumber      = XcellSizePerCore*(YcellSizePerCore*cellCoordZ + cellCoordY) + cellCoordX;

        int cellNumberSpaceAvg = (XcellSizePerCore/spaceAvgCells_X)*((YcellSizePerCore/spaceAvgCells_Y)*cellCoordZ + (cellCoordY/spaceAvgCells_Y)) + (cellCoordX/spaceAvgCells_X);

      	if((cellNumber >= numCellCore) || (cellNumber < 0))
      	{
        		cout<<"cellscream3: "<<Oparticle(countParticleCurrCore,1).position()<<"\t"<<Oparticle(countParticleCurrCore,2).position()<<"\t"<<Oparticle(countParticleCurrCore,3).position()<<"\t"<<cellNumber<<"\t"<<numCellCore<<"\t"<<coreNumber<<endl;
        		cellNumber = numCellCore - 1;
      	}
      	cList[cellNumber].cellParticleList[cList[cellNumber].cellParticleSize] = countParticleCurrCore;
        cList[cellNumber].cellParticleSize++;
        spaceList[cellNumberSpaceAvg].cellParticleList[spaceList[cellNumberSpaceAvg].cellParticleSize] = countParticleCurrCore;
        spaceList[cellNumberSpaceAvg].cellParticleSize++;
        Oparticle.getCellCoord(countParticleCurrCore) = cellNumber;
        countParticleCurrCore++;
    }

    /*******************************************************************************************************/
    /*******************************************************************************************************/
    /*******************************************************************************************************/


    totParticleCurrCore = countParticleCurrCore;
}

template<int DIM,typename real,typename dof>
real hardSphere<DIM,real,dof>::gaussianRandom() //Box-Mueller Transform
{
	real uniRand  = molrng.uniform();
	real uniRand2 = molrng.uniform();

	return pow(-2.*log(uniRand),.5)*sin(2*pi*uniRand2);
}

template<int DIM,typename real,typename dof>
real hardSphere<DIM,real,dof>::gaussianRandomLB2DSMC() //Box-Mueller Transform
{
	real uniRand  = molrng.uniform();
	real uniRand2 = molrng.uniform();
  real insideFactor = 2.*pi*kT*sqrt(2.*pi*kT);

	return pow(-2.*kT*log(uniRand*insideFactor),.5)*sin(2*pi*uniRand2);
}

template<int DIM,typename real,typename dof>
real hardSphere<DIM,real,dof>::relSpeed(int particleOne, int particleTwo)
{
    real relSp = 0.;
    real p1Vel;
    real p2Vel;
    for(int dim = 1; dim <= DIM; dim++)
    {
        p1Vel = Oparticle(particleOne,dim).velocity();
        p2Vel = Oparticle(particleTwo,dim).velocity();
        relSp += (p1Vel - p2Vel)*(p1Vel - p2Vel);
    }

    return pow(relSp,.5);
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::processCollision(int particleOne, int particleTwo)
{
    //real p1Vel;
    //real p2Vel;
    //real p1VelNew;
    //real p2VelNew;
    real vCOM[DIM];
    real vRel[DIM];
    real relSp = relSpeed(particleOne,particleTwo);

    //for(int dim = 1; dim <= DIM; dim++)
    {
        //p1Vel = Oparticle(particleOne,dim).velocity();
        //p2Vel = Oparticle(particleTwo,dim).velocity();
        vCOM[0] = .5*(Oparticle(particleOne,1).velocity() + Oparticle(particleTwo,1).velocity());
      	vCOM[1] = .5*(Oparticle(particleOne,2).velocity() + Oparticle(particleTwo,2).velocity());
      	vCOM[2] = .5*(Oparticle(particleOne,3).velocity() + Oparticle(particleTwo,3).velocity());
    }

    real phiRand = molrng.uniform();
    real phi = 2.*pi*phiRand;

    real thetaRand = 2.*molrng.uniform() - 1.;
    real cosTheta = thetaRand;
    real sinTheta = pow(1. - cosTheta*cosTheta,.5);

    vRel[0] = relSp*sinTheta*cos(phi);
    vRel[1] = relSp*sinTheta*sin(phi);
    vRel[2] = relSp*cosTheta;

    //for(int dim = 1; dim <= DIM; dim++)
    {
        Oparticle(particleOne,1).velocity() = vCOM[0] + .5*vRel[0];
        Oparticle(particleTwo,1).velocity() = vCOM[0] - .5*vRel[0];

	      Oparticle(particleOne,2).velocity() = vCOM[1] + .5*vRel[1];
        Oparticle(particleTwo,2).velocity() = vCOM[1] - .5*vRel[1];

	      Oparticle(particleOne,3).velocity() = vCOM[2] + .5*vRel[2];
        Oparticle(particleTwo,3).velocity() = vCOM[2] - .5*vRel[2];

        //Oparticle(particleOne,dim).velocity() = p1VelNew;
        //Oparticle(particleTwo,dim).velocity() = p2VelNew;
    }
}

template<int DIM,typename real,typename dof>
void hardSphere<DIM,real,dof>::Collision(int coreNumber, int tempCount)
{
    for(int currCell = 0; currCell < numCellCore; currCell++)
    {
        real vrMax = 6.;
        real selectCriteria;
        real relSp;
        int cellPartSize = cList[currCell].cellParticleSize;
        real numCollisions = cellPartSize*(cellPartSize)*deltat*pi*diameter*diameter*vrMax*.5*nE/volume;
        //100 effective particles per particle
        int uRandPairSelect[2];

        for(int collNo = 1; collNo <= (int)numCollisions; collNo++)
        {
            uRandPairSelect[0] = (int)(cellPartSize*molrng.uniform());
            uRandPairSelect[1] = (int)(cellPartSize*molrng.uniform());
            int particleOne = cList[currCell].cellParticleList[uRandPairSelect[0]];
            int particleTwo = cList[currCell].cellParticleList[uRandPairSelect[1]];
            relSp = relSpeed(particleOne,particleTwo);
            selectCriteria = molrng.uniform();

            if((relSp/vrMax) > selectCriteria)
            {
                processCollision(particleOne,particleTwo);
            }
        }
    }
}

#endif // HARDSPHERE_H_INCLUDED
