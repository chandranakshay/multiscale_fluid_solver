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

 *  @Author: Akshay Chandran, Praveen Kolluru and Santosh Ansumali*

 *********************************************************************************************/

//////////
//  Uncomment this if you want to run COUPLED LB and DSMC
//////////
#define COUPLE 1

//////////
//  Uncomment this if you want to run DECOUPLED - LB ands DSMC
//////////
// #define DECOUPLE 1

//////////
//  Uncomment this if you want to run ONLY LB solver
//////////
// #define LBM_ONLY 1


//////////
// Keep only one of the above uncommented to avoid confusion
//////////


#include "coupling.h"


#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>

#include"procInfo3D.h"
#include"gridBCC3DBasic.h"
#include"gridBCC3DPeriodic.h"
#include"latticeDescriptor.h"
#include"D3RQ41.h"
#include"pack41.h"
#include"partial_pack41.h"
#include"marker.h"
#include"equilibrium41.h"
#include"initialConditions.h"
#include"collide41.h"
#include"communication41.h"
#include"partial_communication41.h"
#include"boundaryCondition.h"
#include"inletOutlet.h"
#include"advection41.h"
#include"bounceBack41_withForce.h"
#include"diffuseBounceBack41_withForce.h"
#include"diffuse41.h"
#include"perlinNoise.h"
#include"noiseFunctions.h"
#include<packNoise41.h>
#include<communicationNoise.h>
#include"pvts.h"
#include"restart.h"
#include"momentsForDSMC.h"
#include"euler.h"

int m1, m2, m3;

PANINI_REAL *temp1Psend, *temp1Prec, *temp2Psend, *temp2Prec, *temp3Psend, *temp3Prec, *temp1Msend, *temp1Mrec, *temp2Msend, *temp2Mrec, *temp3Msend, *temp3Mrec;
PANINI_REAL *partial_temp1Psend, *partial_temp1Prec, *partial_temp2Psend, *partial_temp2Prec, *partial_temp3Psend, *partial_temp3Prec, *partial_temp1Msend, *partial_temp1Mrec, *partial_temp2Msend, *partial_temp2Mrec, *partial_temp3Msend, *partial_temp3Mrec;
PANINI_REAL *temp1Psend_Noise, *temp1Prec_Noise, *temp2Psend_Noise, *temp2Prec_Noise, *temp3Psend_Noise, *temp3Prec_Noise, *temp1Msend_Noise, *temp1Mrec_Noise, *temp2Msend_Noise, *temp2Mrec_Noise, *temp3Msend_Noise, *temp3Mrec_Noise;

void allocate(int m1, int m2, int m3)
{
  temp1Psend = allocateForPacking(m2,m3);
  temp1Prec  = allocateForPacking(m2,m3);
  temp2Psend = allocateForPacking(m1,m3);
  temp2Prec  = allocateForPacking(m1,m3);
  temp3Psend = allocateForPacking(m1,m2);
  temp3Prec  = allocateForPacking(m1,m2);
  temp1Msend = allocateForPacking(m2,m3);
  temp1Mrec  = allocateForPacking(m2,m3);
  temp2Msend = allocateForPacking(m1,m3);
  temp2Mrec  = allocateForPacking(m1,m3);
  temp3Msend = allocateForPacking(m1,m2);
  temp3Mrec  = allocateForPacking(m1,m2);
}

void partial_allocate(int m1, int m2, int m3)
{
  partial_temp1Psend = partial_allocateForPacking(m2,m3);
  partial_temp1Prec  = partial_allocateForPacking(m2,m3);
  partial_temp2Psend = partial_allocateForPacking(m1,m3);
  partial_temp2Prec  = partial_allocateForPacking(m1,m3);
  partial_temp3Psend = partial_allocateForPacking(m1,m2);
  partial_temp3Prec  = partial_allocateForPacking(m1,m2);
  partial_temp1Msend = partial_allocateForPacking(m2,m3);
  partial_temp1Mrec  = partial_allocateForPacking(m2,m3);
  partial_temp2Msend = partial_allocateForPacking(m1,m3);
  partial_temp2Mrec  = partial_allocateForPacking(m1,m3);
  partial_temp3Msend = partial_allocateForPacking(m1,m2);
  partial_temp3Mrec  = partial_allocateForPacking(m1,m2);
}

/* Nondimensional parameter setting for simulations and DSMC region */
void DSMCParameters(hardSphere<3,double,dofParticle<double> >& p,int myRank,double& time_MAX)
{
    p.Mach                   = 0.2;
    p.gamma                  = 5./3.;
    p.vSound                 = sqrt(p.gamma*kT/(double)M);
    p.Uc                     = p.Mach*p.vSound;
    p.uWall		               = 0.;//p.Uc;
    p.Kn                     = 1./2400.;
    p.Re                     = 16.*p.Mach*sqrt(p.gamma)/(5.*sqrt(2.*pi)*p.Kn);

    double deltaZBymfp       = 1./(p.Kn*(double)(p.ZcellSize));
    double deltaXBymfp       = 1.666667*deltaZBymfp;
    double deltaYBymfp       = 5.*deltaZBymfp;
    double nParticlesPerCell = p.numParticlesPerCell;
    p.nE	             = 1.;

    p.height                 = sqrt(nParticlesPerCell*dia*dia*sqrt(2)*pi/(deltaXBymfp*deltaYBymfp*deltaZBymfp))/p.Kn;
    p.mfp                    = p.height*p.Kn;
    p.density                = 1./(sqrt(2.)*pi*p.mfp*dia*dia);
    p.deltaZ                 = deltaZBymfp*p.mfp;
    p.deltaY                 = deltaYBymfp*p.mfp;
    p.deltaX                 = deltaXBymfp*p.mfp;
    p.Xbound                 = p.XcellSize*p.deltaX;
    p.Ybound                 = p.YcellSize*p.deltaY;

    p.initHardSphere();

    /*=================================VISCOSITY CORRECTION=============================*/
    double viscTemp          = p.Uc*(p.height)/(p.Re);
    double mft               = p.mfp/sqrt(3.*kT/M);
    double deltatBymft       = 0.2;
    p.deltat                 = deltatBymft*mft;
    double viscCorr          = viscTemp*(1. + (32./(150.*pi))*deltatBymft*deltatBymft)*(1. + (32./(135.*pi))*deltaZBymfp*deltaZBymfp);
    /*==================================================================================*/

    p.viscosity              = viscCorr;
    p.gy                     = 12.*p.Uc*p.viscosity/((p.height)*(p.height));
    time_MAX                 = 10.*(p.height)*(p.height)/p.viscosity;

    p.deltat                 = deltatBymft*mft;
}
/***************************************************************************************************/

int main(int argc, char *argv[])
{
    int rank, size;
    /*Initialize DSMC solver (spanwise total cells, streamwise tot cells, wall normal tot (virtual) cells, span procs, stream procs, wall procs, 
    num particles per cell, num cells from wall)*/
    hardSphere<3,double,dofParticle<double> > p(1200,720,1000,20,10,2,30,4);
    /*********************************/
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double time_MAX;
    int simulationTime;
    /*Set DSMC parameters*/
    DSMCParameters(p,rank,time_MAX);

    /*Coupling region (buffer zone) distance from either walls (in terms of number of DSMC cells)*/
    PANINI_INT  DSMC2LBCouplingCellZ      = 0.5*p.ZcellSizePerCore;
    PANINI_REAL distToCouplingZone        = (DSMC2LBCouplingCellZ + 1)*p.deltaZ;
    /*********************************************************************************************/

    PANINI_REAL distToDSMCRegionTopBottom = 2.*distToCouplingZone;
    PANINI_REAL LBbyDSMCmfp               = 20.;

    int myRank;
    myRank = rank;
    /*Initialize LB solver (streamwise cells per proc, wall cells per proc, spanwise cells per proc, stream procs, wall procs, span procs, 
    stream ghost, wall ghost, span ghost)*/
    latticeInfo simParam(72,100,20,10,2,20,2,2,2);
    /*******************************/

    /*LBM RD3Q41 grid initialization*/
    gridBCC3D<4,11,PANINI_REAL>  gridLBM41(simParam.n1,simParam.n2,simParam.n3, simParam.nGhost1,simParam.nGhost2,simParam.nGhost3);
    gridBCC3D<4,1 ,PANINI_REAL>   testGrid(simParam.n1,simParam.n2,simParam.n3, simParam.nGhost1,simParam.nGhost2,simParam.nGhost3);
    gridBCC3D<1,1 ,PANINI_INT >     marker(simParam.n1,simParam.n2,simParam.n3, simParam.nGhost1,simParam.nGhost2,simParam.nGhost3);
    /********************************/

    /*LBM solid fluid marker indicators*/
    PANINI_INT SOLID = 1;
    PANINI_INT FLUID = 2;
    /***********************************/


    m1 = gridLBM41.n1; m2 = gridLBM41.n2; m3 = gridLBM41.n3;
    allocate(m1,m2,m3);
    partial_allocate(m1,m2,m3);

    procInfo3D  CartersianGrid3D(simParam.nP1,simParam.nP2,simParam.nP3, true, true, true, true);
    CartersianGrid3D.getnList();

    MPI_Request requestF1[4];
    MPI_Request requestF2[4];
    MPI_Request requestF3[4];

    MPI_Status statusF1[4];
    MPI_Status statusF2[4];
    MPI_Status statusF3[4];

    int procRank = CartersianGrid3D.coordinates[2]*simParam.nP1*simParam.nP2 + CartersianGrid3D.coordinates[1]*simParam.nP1 + CartersianGrid3D.coordinates[0];

    std::ofstream file1;   file1.open("ZPZ2H.txt");
    std::ofstream file2;   file2.open("ZP25H.txt");
    std::ofstream file3;   file3.open("ZP5H.txt");
    std::ofstream file4;   file4.open("velocityAtCenter.txt");
    std::ofstream file5;   file5.open("momentsAtCenter.txt");

    int ZPZ2H = (int)(0.02*simParam.n2*simParam.nP2) + gridLBM41.nB2;
    int ZP25H = (int)(0.25*simParam.n2*simParam.nP2) + gridLBM41.nB2;
    int ZP5H  = (int)(0.50*simParam.n2*simParam.nP2) + gridLBM41.nB2;


    char restartFile[502];
    int  Rcount=6;

    lbmRD3Q41<PANINI_REAL> lbModel41(1.0);
    setModelParameters(lbModel41);
    getLatticeParameter(lbModel41);

    /*LBM parameter initialization (Nondimensional numbers same as DSMC)*/
    PANINI_REAL cs2                = p.vSound*p.vSound;//5.0/3.0*lbModel41.theta0;//
    PANINI_REAL cs                 = sqrt(cs2);
    PANINI_REAL c                  = sqrt(cs2/(p.gamma*lbModel41.theta0));
    PANINI_REAL Ma                 = p.Mach;
    PANINI_REAL refLen             = p.height;// - distToDSMCRegionTopBottom;//(simParam.n2*simParam.nP2);//
    PANINI_REAL refLen2            = p.height;// - (4.*p.deltaZ);// - distToDSMCRegionTopBottom;//(simParam.n2*simParam.nP2);//
    PANINI_REAL u_inlet            = Ma*cs;
    const PANINI_REAL Re           = p.Re;
    PANINI_REAL rhoLB              = p.density;
    PANINI_REAL uRef               = u_inlet;
    PANINI_REAL kinVisc            = p.viscosity;
    PANINI_REAL tau                = kinVisc/(lbModel41.theta0);
    PANINI_REAL dx                 = refLen2/(simParam.n2*simParam.nP2);
    PANINI_REAL dt                 = dx/c;
    PANINI_REAL tauN3              = tau/dt;
    PANINI_REAL beta               = 1.0/(2.0*tauN3 + 1.0);
    PANINI_REAL twoBeta            = 2.0*beta;
    PANINI_REAL betaTau            = beta*tau;
    PANINI_REAL oneMinustwoBeta    = 1.0-twoBeta;
    PANINI_REAL convectionTime     = (refLen*refLen)/kinVisc;
    PANINI_REAL bodyConvectionTime = refLen/u_inlet;
    PANINI_REAL f2g_factor         = 2.*tau/(2.*tau+dt);
    PANINI_REAL dtByTwoTau         = (dt/(2.0*tau));
    /*******************************************************************/


    int VECT_LENGTH = 4;
    PANINI_REAL forceNodeBefore[3],forceCellBefore[3];
    PANINI_REAL forceNodeAfter[3],forceCellAfter[3];
    PANINI_REAL NetForce[3];
    PANINI_REAL ForceX(0.0),ForceY(0.0),ForceZ(0.0),cDOLD(0.0),cD(0.0),cLOLD(0.0),cL(0.0),cZOLD(0.0),cZ(0.0);
    PANINI_REAL uXBody = 0.0;
    PANINI_REAL uYBody = 0.0;
    PANINI_REAL uZBody = 0.0;

    int centerX =  (int)(0.5*simParam.n1*simParam.nP1) + gridLBM41.nB1;
    int centerY =  (int)(0.5*simParam.n2*simParam.nP2) + gridLBM41.nB2;
    int centerZ =  (int)(0.5*simParam.n3*simParam.nP3) + gridLBM41.nB3;

    /*External body force (taken from DSMC structure)*/
    PANINI_REAL bodyForceX = p.gy;//12.0*uRef*kinVisc/(refLen*refLen); //1.12598889221483e-08;
    PANINI_REAL bodyForceY = 0.0;
    PANINI_REAL bodyForceZ = 0.0;
    PANINI_REAL idealBodyForce = bodyForceX;
    /*************************************************/

    PANINI_INT  numPointsLBX = p.Ybound*simParam.n2*simParam.nP2/(refLen);
    PANINI_INT  numPointsLBY = dx*simParam.n2*simParam.nP2*simParam.n2*simParam.nP2/(refLen);
    PANINI_INT  numPointsLBZ = p.Xbound*simParam.n2*simParam.nP2/(refLen);

    /*Coupling processors*/
    int lowerCoupleLB     = 0;
    int upperCoupleLB     = (p.nCoresZ - 1);
    int lowerCoupleDSMC   = 0;
    int upperCoupleDSMC   = (p.nCoresZ - 1);
    int lowerCoupleDSMCLB = 0;
    int upperCoupleDSMCLB = (p.nCoresZ - 1);
    /*********************/

    if(myRank == 0)
    {
        std::cout<<"distanceMeasure: "<<p.height<<"\t"<<refLen<<"\t"<<refLen2<<"\t"<<distToDSMCRegionTopBottom<<"\t"<<dx*simParam.n2*simParam.nP2 + distToDSMCRegionTopBottom<<std::endl<<std::endl;

        std::cout<<"numPointsLBX   : "<<numPointsLBX<<"\t"<<simParam.n1*simParam.nP1<<"\t"<<dx<<std::endl;
        std::cout<<"numPointsLBY   : "<<numPointsLBY<<"\t"<<simParam.n2*simParam.nP2<<"\t"<<dx<<std::endl;
        std::cout<<"numPointsLBZ   : "<<numPointsLBZ<<"\t"<<simParam.n3*simParam.nP3<<"\t"<<dx<<std::endl;
    }

    if(myRank==0)
    {
       std::cout << "---------------------------------------------------------------------------"        <<std::endl;
       std::cout << "Parameters : " << std::setw(16) << "LB"             << std::setw(16) << "DSMC"      << std::setw(16) << "Ratio"                <<std::endl;
       std::cout << "---------------------------------------------------------------------------"        <<std::endl;
       std::cout << "Re         : " << std::setw(16) << Re               << std::setw(16) << p.Re        << std::setw(16) << Re/p.Re                <<std::endl;
       std::cout << "Ma         : " << std::setw(16) << Ma               << std::setw(16) << p.Mach      << std::setw(16) << Ma/p.Mach              <<std::endl;
       std::cout << "Kn         : " << std::setw(16) << p.Kn             << std::setw(16) << p.Kn      << std::setw(16) << p.Kn/p.Kn              <<std::endl;
       std::cout << "cs         : " << std::setw(16) << cs               << std::setw(16) << p.vSound    << std::setw(16) << cs/p.vSound            <<std::endl;
       std::cout << "U0         : " << std::setw(16) << uRef             << std::setw(16) << p.Uc        << std::setw(16) << uRef      /p.Uc        <<std::endl;
       std::cout << "Height     : " << std::setw(16) << refLen           << std::setw(16) << p.height    << std::setw(16) << refLen    /p.height    <<std::endl;
       std::cout << "Visc       : " << std::setw(16) << kinVisc          << std::setw(16) << p.viscosity << std::setw(16) << kinVisc   /p.viscosity <<std::endl;
       std::cout << "density    : " << std::setw(16) << rhoLB            << std::setw(16) << p.density   << std::setw(16) << rhoLB/p.density        <<std::endl;
       std::cout << "bodyForceX : " << std::setw(16) << bodyForceX       << std::setw(16) << p.gy        << std::setw(16) << bodyForceX/p.gy        <<std::endl<<std::endl;

       // Physical parameters
       std::cout << "-----------------------" << std::endl;
       std::cout << "Channel dimensions : "   << std::endl;
       std::cout << "-----------------------" << std::endl;
       std::cout << "LB                 : " << dx*simParam.n1*simParam.nP1 << " x "
                                            << dx*simParam.n2*simParam.nP2 << " x "
                                            << dx*simParam.n3*simParam.nP3 << std::endl;

       std::cout << "DSMC               : " << p.deltaY * p.YcellSize << " x "
                                            << p.deltaZ * p.ZcellSize << " x "
                                            << p.deltaX * p.XcellSize << std::endl<<"\t"<<p.Xbound<< std::endl;

       // LB(x,y,z) = DSMC(y,z,x)
       std::cout << "-----------------------------------------------"      << std::endl;
       std::cout << "Spatial discretization(dx)  " <<"LB" << std::setw(16) <<"    DSMC"      << std::endl;
       std::cout << "-----------------------------------------------"      << std::endl;
       std::cout << "Streamwise direction      : " << dx  << std::setw(16) << p.deltaY   << std::endl;
       std::cout << "Wall direction            : " << dx  << std::setw(16) << p.deltaZ   << std::endl;
       std::cout << "Spanwise direction        : " << dx  << std::setw(16) << p.deltaX   << std::endl<< std::endl;


       std::cout << "-----------------------------------------------"      << std::endl;
       std::cout << "Time discretization(dt)     " <<"LB" << std::setw(16) <<"    DSMC"      << std::endl;
       std::cout << "-----------------------------------------------"      << std::endl;
       std::cout << "                            " << dt  << std::setw(16) << p.deltat   << std::endl << std::endl;


       std::cout << "Ratio of time steps    : " << dt/p.deltat                                  << std::endl << std::endl;
       std::cout << "Diffusion Time         : " << (p.height)*(p.height)/p.viscosity            << std::endl << std::endl;
       std::cout << "Steps to diffusion Time: " << (p.height)*(p.height)/(p.viscosity*p.deltat) << std::endl << std::endl;


       std::cout << "---------------------------------------------------------------" << std::endl;
       std::cout << "Number of                       " << "LB points "                << std::setw(16) <<  "DSMC Cells"  << std::endl;
       std::cout << "---------------------------------------------------------------" << std::endl;
       std::cout << "Streamwise direction          : " << simParam.n1*simParam.nP1    << std::setw(16) <<  p.YcellSize  << std::endl;
       std::cout << "Wall direction                : " << simParam.n2*simParam.nP2    << std::setw(16) <<  p.ZcellSize  << std::endl;
       std::cout << "Spanwise direction            : " << simParam.n3*simParam.nP3    << std::setw(16) <<  p.XcellSize  << std::endl<< std::endl;


       double pkk_DSMCToLB_stream = (PANINI_REAL)(simParam.n1*simParam.nP1)/(PANINI_REAL)p.YcellSize;
       double pkk_DSMCToLB_wall   = (PANINI_REAL)(simParam.n2*simParam.nP2)/(PANINI_REAL)p.ZcellSize;
       double pkk_DSMCToLB_span   = (PANINI_REAL)(simParam.n3*simParam.nP3)/(PANINI_REAL)p.XcellSize;

       std::cout << "-----------------------------------------------"   << std::endl;
       std::cout << "Each DSMC cell corresponds to : " << std::setw(16) << "    LB points "   << std::endl;
       std::cout << "-----------------------------------------------"   << std::endl;
       std::cout << "Streamwise direction          : " << std::setw(10) << pkk_DSMCToLB_stream  << std::endl;
       std::cout << "Wall direction                : " << std::setw(10) << pkk_DSMCToLB_wall    << std::endl;
       std::cout << "Spanwise direction            : " << std::setw(10) << pkk_DSMCToLB_span    << std::endl<< std::endl;

       double pkk_LBToDSMC_stream = (PANINI_REAL)p.YcellSize/(PANINI_REAL)(simParam.n1*simParam.nP1) ;
       double pkk_LBToDSMC_wall   = (PANINI_REAL)p.ZcellSize/(PANINI_REAL)(simParam.n2*simParam.nP2) ;
       double pkk_LBToDSMC_span   = (PANINI_REAL)p.XcellSize/(PANINI_REAL)(simParam.n3*simParam.nP3) ;

       std::cout << "-----------------------------------------------"   << std::endl;
       std::cout << "Each LB point corresponds to  : " << std::setw(10) << "    Cells "   << std::endl;
       std::cout << "-----------------------------------------------"   << std::endl;
       std::cout << "Streamwise direction          : " << std::setw(10) << pkk_LBToDSMC_stream << std::endl;
       std::cout << "Wall direction                : " << std::setw(10) << pkk_LBToDSMC_wall   << std::endl;
       std::cout << "Spanwise direction            : " << std::setw(10) << pkk_LBToDSMC_span   << std::endl<< std::endl;

       std::cout << "Mean free path                : " << std::setw(16) << p.mfp            << std::endl;
    }

    int restartTimeint = (int)(2500.0*convectionTime);

    if(myRank==0)
    {
     std::cout<<"Parameters : "<<std::setw(16)<<"LB"<<std::setw(16)<<"DSMC"<<std::endl;
     std::cout<<"Re         : "<<std::setw(16)<<Re<<std::setw(16)<<p.Re<<std::endl;
     std::cout<<"Ma         : "<<std::setw(16)<<Ma<<std::setw(16)<<p.Mach<<std::endl;
     std::cout<<"dt:        : "<<std::setw(16)<<dt<<std::setw(16)<<p.deltat<<std::endl;
     std::cout<<"dt LB/DSMC : "<<std::setw(16)<<dt/p.deltat<<std::endl;
     std::cout<<"dt DSMC/LB : "<<std::setw(16)<<p.deltat/dt<<std::endl;
     std::cout<<"Visc       : "<<std::setw(16)<<kinVisc<<std::setw(16)<<p.viscosity<<std::endl;
     std::cout<<"U0         : "<<std::setw(16)<<uRef<<std::setw(16)<<p.Uc<<std::endl;
     std::cout<<"Height     : "<<std::setw(16)<<refLen<<std::setw(16)<<p.height<<std::endl;
     std::cout<<"bodyForceX : "<<std::setw(16)<<bodyForceX<<std::setw(16)<<p.gy<<std::endl;
     std::cout<<"mfp        : "<<std::setw(16)<<p.mfp<<std::endl;
     std::cout<<"mft        : "<<std::setw(16)<<p.mfp/sqrt(2.*kT/M)<<std::endl;
     std::cout<<"dt/mft     : "<<std::setw(16)<<p.deltat/(p.mfp/sqrt(2.*kT/M))<<std::endl;
     std::cout<<"density    : "<<std::setw(16)<<p.density<<std::endl;
     std::cout<<"mfpbydeltaZ: "<<std::setw(16)<<p.Xbound<<std::endl;
     std::cout<<"c          : "<<std::setw(16)<<c<<std::endl;
     std::cout<<"c2         : "<<std::setw(16)<<c*c<<std::endl;
     std::cout<<"f2g_factor : "<<std::setw(16)<<f2g_factor<<std::endl;
     std::cout<<"cs2        : "<<std::setw(16)<<p.vSound*p.vSound<<"   "<<cs2/(5.0*lbModel41.theta0/3.0) <<std::endl;
     std::cout<<"cs         : "<<std::setw(16)<<cs<<std::endl<<std::endl<<std::endl;

    }

    /*Number of LB steps within inner loop*/
    int LB_Count              = 100;
    /*Number of DSMC steps within LB steps*/
    int DSMC_Count            = (int)(LB_Count*dt/p.deltat);
    int DSMC_Count_Actual     = (int)(LB_Count*dt/p.deltat);
    /**************************************/

    PANINI_REAL DSMC_CountInv = 1./((PANINI_REAL)DSMC_Count);
    PANINI_REAL LB_CountInv   = 1./((PANINI_REAL)LB_Count);

    /*weights for interpolation of non-aligned cells in the buffer zone*/
    double weightLB_P1 = 0.;//((LBbyDSMCmfp*p.mfp) - distToCouplingZone)/(LBbyDSMCmfp*p.mfp);
    double weightLB_P2 = 1.;//distToCouplingZone/(LBbyDSMCmfp*p.mfp);
    /*******************************************************************/

    if(myRank==0)
    {
     std::cout<<"ratio of time steps DSMC/LB: "<<dt/p.deltat<<"  "<<(int)(dt/p.deltat)<<std::endl;
     std::cout<<"LB steps:                  : "<<LB_Count<<std::endl;
     std::cout<<"DSMC steps:                : "<<DSMC_Count<<std::endl<<std::endl<<std::endl;;
     std::cout<<"weights Interpolation      : "<<weightLB_P1<<"\t"<<weightLB_P2<<std::endl<<std::endl<<std::endl;
    }



    char LBrestartFile[150];
    sprintf(LBrestartFile,"./../data/LBData/LBdata_%d.dat",myRank);
    //string LBrestartFile = "./../LBdata" + to_string(myRank) + ".dat";
    //char const* LBrestartFile = LBrestartFileString.c_str();

    /*LBM initial conditions*/
    initializeMarker(marker,FLUID);

    PANINI_REAL velScale = 0.*1.5;
    initialConditions(lbModel41,gridLBM41,marker,VECT_LENGTH,uRef*velScale,FLUID,simParam.nP2,CartersianGrid3D.coordinates,p.density);
    //restartIC(lbModel41, gridLBM41,LBrestartFile);
    MPI_Barrier(MPI_COMM_WORLD);
    /************************/


    /* Perlin based curl noise (uncomment only for initial noise in velocity field)*/
    PANINI_REAL normg;
    curl_noise(lbModel41, gridLBM41, simParam.nP1, simParam.nP2, simParam.nP3,
             normg, CartersianGrid3D.coordinates, VECT_LENGTH, p.density, uRef*1.5);
    perturb_parab(CartersianGrid3D.coordinates, lbModel41, gridLBM41,
                simParam.nP1, simParam.nP2, simParam.nP3, normg,
                uRef*1.5, VECT_LENGTH, p.density, uRef*1.5);
    MPI_Barrier(MPI_COMM_WORLD);
    /*******************************************************************************/

    simulationTime = int(1.1*convectionTime);

    /* Temporary variables */
    int count      = 0;
    int momCount   = 1;
    int step       = 0;
    double totTime = 0.;

    clock_t time1 = 0.,time2 = 0.,time3 = 0.,time4 = 0.,time5 = 0.,time6 = 0.,time7 = 0.,time8 = 0.,time9 = 0.,time10 = 0.,timeDSMC = 0.,timeLB = 0.;
  double rho_LB_1Temp2[4],uX_LB_1Temp2[4],uY_LB_1Temp2[4],uZ_LB_1Temp2[4],theta_LB_1Temp2[4],Pxx_LB_1Temp2[4],Pyy_LB_1Temp2[4],Pzz_LB_1Temp2[4],Pxy_LB_1Temp2[4],Pyz_LB_1Temp2[4],Pzx_LB_1Temp2[4];
  double q1_LB_1Temp2[4],q2_LB_1Temp2[4],q3_LB_1Temp2[4];

    int LBbyDSMCpointsY = (int)((PANINI_REAL)p.YcellSize/(PANINI_REAL)(simParam.n1*simParam.nP1));
    int LBbyDSMCpointsZ = (int)((PANINI_REAL)p.ZcellSize/(PANINI_REAL)(simParam.n2*simParam.nP2));
    int LBbyDSMCpointsX = (int)((PANINI_REAL)p.XcellSize/(PANINI_REAL)(simParam.n3*simParam.nP3));

    int numMoments = 14;

    int numSpaceAvgCells_LB_Y = 1;
    int numSpaceAvgCells_LB_X = 1;

    int XcellSizePerCore = p.XcellSize/p.nCoresX;
    int YcellSizePerCore = p.YcellSize/p.nCoresY;

    int numSpaceAvgCells_DSMC_Y = LBbyDSMCpointsY*numSpaceAvgCells_LB_Y;
    int numSpaceAvgCells_DSMC_X = LBbyDSMCpointsX*numSpaceAvgCells_LB_X;

    p.spaceAvgCells_X = numSpaceAvgCells_DSMC_X;
    p.spaceAvgCells_Y = numSpaceAvgCells_DSMC_Y;

    p.spaceMemoryAlloc();

    int numSpaceAvgCells_LB   = numSpaceAvgCells_LB_X*numSpaceAvgCells_LB_Y;
    int numSpaceAvgCells_DSMC = numSpaceAvgCells_DSMC_X*numSpaceAvgCells_DSMC_Y;

    double moments_LB_Send[numMoments*(XcellSizePerCore*YcellSizePerCore/numSpaceAvgCells_DSMC)];
    double moments_LB_Recv[numMoments*(XcellSizePerCore*YcellSizePerCore/numSpaceAvgCells_DSMC)];

    double moments_DSMC_Send[numMoments*(XcellSizePerCore*YcellSizePerCore/numSpaceAvgCells_DSMC)];
    double moments_DSMC_Recv[numMoments*(XcellSizePerCore*YcellSizePerCore/numSpaceAvgCells_DSMC)];

    double momentsSave[numMoments*(XcellSizePerCore*YcellSizePerCore/numSpaceAvgCells_DSMC)];
    /***********************/

    /*Processor coordinates for DSMC and LB regions*/
    int myRankZ          = myRank/(p.nCoresX*p.nCoresY);
    int myRankSlab       = myRank%(p.nCoresX*p.nCoresY);
    int myRankY          = myRankSlab/p.nCoresX;
    int myRankX          = myRankSlab%p.nCoresX;

    int myRankY_LB       = myRank/(p.nCoresX*p.nCoresZ);
    int myRankSlab_LB    = myRank%(p.nCoresX*p.nCoresZ);
    int myRankZ_LB       = myRankSlab_LB/p.nCoresX;
    int myRankX_LB       = myRankSlab_LB%p.nCoresX;
    /************************************************/

    /*Cartesian convention for the two solvers*/
    //DSMC (x,y,z) = LB (z,x,y)

    /*Initialize additional dependent DSMC variables*/
    p.extraVarInit(myRank);

    int centreX_DSMC = (int)(0.5*p.XcellSize/numSpaceAvgCells_DSMC_X);
    int centreY_DSMC = (int)(0.5*(p.YcellSize/p.nCoresY)/numSpaceAvgCells_DSMC_Y);


    double moments_To_Print[numMoments];


    //restartIC(lbModel41, gridLBM41,LBrestartFile);
    /*Initial conditions for DSMC region (Gaussian for equilibrium)*/
    p.initGenerateDistributionFromMoments(p.density,p.density*kT,myRank,velScale);
    //p.retrieveParticleDataFromFile(myRank);

    std::ofstream fileTotalDensity("./densityTot.dat");

    while(totTime < time_MAX)
    {

       MPI_Barrier(MPI_COMM_WORLD);

       /* Calculating moments from the LB nodes within the buffer region and 
       	sending to the corresponding processor in the DSMC region (One time before DSMC evolution and collision for n time steps)*/
       {
         if(myRankZ_LB == lowerCoupleLB)
         {
           /* Calculate and send averaged moments in the lower LB buffer zone*/
           calcSendLBmomentsLowerBufferZone(p,lbModel41,gridLBM41,moments_LB_Send,weightLB_P1,weightLB_P2,LBbyDSMCpointsX,LBbyDSMCpointsY,numSpaceAvgCells_LB_X,numSpaceAvgCells_LB_Y,           				numSpaceAvgCells_LB,numSpaceAvgCells_DSMC_X,numSpaceAvgCells_DSMC_Y,numSpaceAvgCells_DSMC,myRank,myRankX,myRankY,myRankZ,
		           myRankX_LB,myRankY_LB,myRankZ_LB,count,lowerCoupleDSMC);
		           
  	   int rankToSendTo = p.nCoresX*(p.nCoresY*lowerCoupleDSMC + myRankY_LB) + myRankX_LB;
     	   MPI_Isend(moments_LB_Send,p.numMoments*(p.XcellSizePerCore*p.YcellSizePerCore/numSpaceAvgCells_DSMC),MPI_DOUBLE,rankToSendTo,201,MPI_COMM_WORLD,&requestF1[0]);
         }

         if(myRankZ == lowerCoupleDSMC)
         {
           /* Receive moments in the lower DSMC buffer zone*/
           int rankRecvFrom = p.nCoresX*(p.nCoresZ*myRankY + lowerCoupleLB) + myRankX;
              MPI_Irecv(moments_LB_Recv,numMoments*(XcellSizePerCore*YcellSizePerCore/numSpaceAvgCells_DSMC),MPI_DOUBLE,rankRecvFrom,201,MPI_COMM_WORLD,&requestF1[0]);
              MPI_Wait(&requestF1[0], &statusF1[0]);
              
           calcRecvLBmomentsLowerBufferZone(p,lbModel41,gridLBM41,moments_LB_Recv,momentsSave,c,f2g_factor,weightLB_P1,weightLB_P2,LBbyDSMCpointsX,LBbyDSMCpointsY,numSpaceAvgCells_LB_X,
           		numSpaceAvgCells_LB_Y,numSpaceAvgCells_LB,numSpaceAvgCells_DSMC_X,numSpaceAvgCells_DSMC_Y,numSpaceAvgCells_DSMC,myRank,myRankX,myRankY,myRankZ,
		           myRankX_LB,myRankY_LB,myRankZ_LB,count,lowerCoupleDSMC);
         }

         if(myRankZ_LB == upperCoupleLB)
         {
           /* Calculate and send averaged moments in the upper LB buffer zone*/
           calcSendLBmomentsUpperBufferZone(p,lbModel41,gridLBM41,moments_LB_Send,weightLB_P1,weightLB_P2,LBbyDSMCpointsX,LBbyDSMCpointsY,numSpaceAvgCells_LB_X,numSpaceAvgCells_LB_Y,           				numSpaceAvgCells_LB,numSpaceAvgCells_DSMC_X,numSpaceAvgCells_DSMC_Y,numSpaceAvgCells_DSMC,myRank,myRankX,myRankY,myRankZ,
		           myRankX_LB,myRankY_LB,myRankZ_LB,count,upperCoupleDSMC);

           int rankToSendTo = p.nCoresX*(p.nCoresY*upperCoupleDSMC + myRankY_LB) + myRankX_LB;
              MPI_Isend(moments_LB_Send,numMoments*(XcellSizePerCore*YcellSizePerCore/numSpaceAvgCells_DSMC),MPI_DOUBLE,rankToSendTo,205,MPI_COMM_WORLD,&requestF1[0]);
         }

         if(myRankZ == upperCoupleDSMC)
         {
           /* Receive moments in the upper DSMC buffer zone*/
           int rankRecvFrom = p.nCoresX*(p.nCoresZ*myRankY + upperCoupleLB) + myRankX;
              MPI_Irecv(moments_LB_Recv,numMoments*(XcellSizePerCore*YcellSizePerCore/numSpaceAvgCells_DSMC),MPI_DOUBLE,rankRecvFrom,205,MPI_COMM_WORLD,&requestF1[0]);
              MPI_Wait(&requestF1[0], &statusF1[0]);
              
           calcRecvLBmomentsUpperBufferZone(p,lbModel41,gridLBM41,moments_LB_Recv,momentsSave,c,f2g_factor,weightLB_P1,weightLB_P2,LBbyDSMCpointsX,LBbyDSMCpointsY,numSpaceAvgCells_LB_X,
           		numSpaceAvgCells_LB_Y,numSpaceAvgCells_LB,numSpaceAvgCells_DSMC_X,numSpaceAvgCells_DSMC_Y,numSpaceAvgCells_DSMC,myRank,myRankX,myRankY,myRankZ,
		           myRankX_LB,myRankY_LB,myRankZ_LB,count,upperCoupleDSMC);
	 }
       }

       /*Clear moment fluxes in the buffer zone*/
       p.clearMomentsFlux(myRank);

       /* Start of DSMC region simulation with LB moments in the buffer layer*/
       for(int tempCount = 0; tempCount < DSMC_Count; tempCount++)
       {
          auto begin1 = std::chrono::high_resolution_clock::now();

          // LB to DSMC Coupling Routing Starts here

          p.reGenIndex = 0;

          // if(tempCount == 0)
          {
              for(int cellY = 0; cellY < (YcellSizePerCore/LBbyDSMCpointsY); cellY = cellY + numSpaceAvgCells_LB_Y)
              {
                for(int cellX = 0; cellX < (XcellSizePerCore/LBbyDSMCpointsX); cellX = cellX + numSpaceAvgCells_LB_X)
                {
                  int cellNumberLocal = ((int)(XcellSizePerCore/(numSpaceAvgCells_DSMC_X)))*((int)(cellY/numSpaceAvgCells_LB_Y)) + (int)(cellX/numSpaceAvgCells_LB_X);

                  if((cellNumberLocal == 0) && (myRankY == ((int)(0.5*p.nCoresY))) && (myRankX == ((int)(0.5*p.nCoresX))))
                  {
                    std::cout<<"====================================================================================================================="<<std::endl;
                    cout<<"rho_DSMC@: "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal]<<std::endl;
                    cout<<"uX_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 2]<<std::endl;
                    cout<<"uY_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 3]<<std::endl;
                    cout<<"uZ_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 4]<<std::endl;
                    cout<<"t0_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 1]/momentsSave[numMoments*cellNumberLocal]<<std::endl;
                    cout<<"sXX_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 5]<<std::endl;
                    cout<<"SXY_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 6]<<std::endl;
                    cout<<"sYZ_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 7]<<std::endl;
                    cout<<"sYY_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 8]<<std::endl;
                    cout<<"sZZ_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 9]<<std::endl;
                    cout<<"sZX_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 10]<<std::endl;
                    cout<<"qX_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 11]<<std::endl;
                    cout<<"qY_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 12]<<std::endl;
                    cout<<"qZ_DSMC@ : "<<count<<" "<<myRankZ<<" "<<myRank<<":"<<momentsSave[numMoments*cellNumberLocal + 13]<<std::endl;
                    std::cout<<"====================================================================================================================="<<std::endl;
                  }

		  /* Generate new particles in the DSMC buffer layer based on moments received from LB buffer zone*/
                  p.generateDistributionFromMoments(&(momentsSave[numMoments*cellNumberLocal]),cellX/numSpaceAvgCells_LB_X,cellY/numSpaceAvgCells_LB_Y,myRank,count);

                }
              }
          }
 
 	  /*Averaged moments in the DSMC region. Averaged moments in the DSMC buffer layer will be sent to 
 	  LB buffer zone within LB loop*/
          p.MomentChainSpaceAvg(myRank,DSMC_Count,momCount,count);

          // LB to DSMC Coupling Routing Ends here

          // auto begin2 = std::chrono::high_resolution_clock::now();
          // time1 += std::chrono::duration_cast<std::chrono::nanoseconds>(begin2-begin1).count();

	  /*Voidlist data structure for adding new particles to particle list in every processor*/
          p.sortVoidList();
          /*Clearing particle IDs in cells*/
          p.clearCellParticles(myRank);
          if(myRank == 0)
          {
            cout<<"Cleared!"<<std::endl;
          }

	  /*Calculate total number of particles in the DSMC region*/
          if(myRank != 0)
          {
            MPI_Send(&(p.totParticleCurrCore),1,MPI_INT,0,61,MPI_COMM_WORLD);
          }

          if(myRank == 0)
          {
            p.sumParticles = p.totParticleCurrCore;
            // cout<<"Particles Per Core: "<<count<<" "<<0<<" "<<(p.totParticleCurrCore)<<" "<<p.volume<<std::endl;
            int recvP;
            for(int i = 1; i < p.nCores; i++)
            {
              MPI_Recv(&recvP,1,MPI_INT,i,61,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
              //MPI_Wait(&requestF1[0], &statusF1[0]);
              if((((i/p.nCoresY) <= lowerCoupleDSMC) || ((i/p.nCoresY) >= upperCoupleDSMC)) && ((i%p.nCoresY) == (0.5*p.nCoresY)))
              {
                cout<<"Particles Per Core Lower: "<<count<<" "<<i/p.nCoresY<<" "<<i<<" "<<recvP<<" "<<p.volume<<std::endl;
              }
              p.sumParticles += recvP;
            }
            cout<<"Total Particle Count Lower: "<<count<<" "<<p.sumParticles<<std::endl<<std::endl;;

            fileTotalDensity<<count<<"\t"<<p.sumParticles<<std::endl;
          }

	  /*DSMC region particle advection routine*/
          p.evolveSystemDiffuse(myRank,count,DSMC2LBCouplingCellZ);
          // p.voidListSize = 0;

          if(myRank == 0)
          {
            cout<<"Evolved!"<<std::endl;
          }

          // auto begin3 = std::chrono::high_resolution_clock::now();
          // time2 += std::chrono::duration_cast<std::chrono::nanoseconds>(begin3-begin2).count();

	  /*MPI communication for sending particles across processors*/
          p.updateCellCoord(myRank,DSMC2LBCouplingCellZ);

          if(myRank == 0)
          {
            cout<<"Updated!"<<std::endl;
          }

          // auto begin4 = std::chrono::high_resolution_clock::now();
          // time3 += std::chrono::duration_cast<std::chrono::nanoseconds>(begin4-begin3).count();
          //
          // auto begin5 = std::chrono::high_resolution_clock::now();
          // time4 += std::chrono::duration_cast<std::chrono::nanoseconds>(begin5-begin4).count();

	  /*DSMC region collision routine*/
          p.Collision(myRank,tempCount);

          if(myRank == 0)
          {
            cout<<"Collided!"<<std::endl;
          }

          // auto begin6 = std::chrono::high_resolution_clock::now();
          // time5 += std::chrono::duration_cast<std::chrono::nanoseconds>(begin6-begin5).count();


          if(myRank == lowerCoupleLB)
          {
            cout<<"count  : "<<count<<std::endl;
            cout<<"time   : "<<totTime<<std::endl;
            cout<<"totTime: "<<time_MAX<<std::endl;
            cout<<"transKE: "<<p.totKinEnergy()<<std::endl<<std::endl;;
          }

          MPI_Barrier(MPI_COMM_WORLD);


          auto begin7 = std::chrono::high_resolution_clock::now();
          // time6 += std::chrono::duration_cast<std::chrono::nanoseconds>(begin7-begin6).count();
          //
          timeDSMC += std::chrono::duration_cast<std::chrono::nanoseconds>(begin7-begin1).count();


          count++;
          momCount++;
          totTime += p.deltat;

       }

      /*Start of LB routine*/

      MPI_Barrier(MPI_COMM_WORLD);

      for(int tempCount=0;tempCount<LB_Count;tempCount++)
      {
          auto begin8 = std::chrono::high_resolution_clock::now();
          // auto time2 += std::chrono::duration_cast<std::chrono::nanoseconds>(begin3-begin2).count();

          /* Calculating moments from the DSMC space cells within the buffer layer and 
       	  sending to the corresponding processor in the LB region*/

          if(myRankZ == lowerCoupleDSMC)
          {
              /*Pack and send averaged DSMC moments in the lower buffer layer to the corresponding processors in the LB region*/
              calcSendDSMCmomentsLowerBufferZone(p,moments_DSMC_Send,DSMC2LBCouplingCellZ,LBbyDSMCpointsX,LBbyDSMCpointsY,numSpaceAvgCells_LB_X,numSpaceAvgCells_LB_Y,numSpaceAvgCells_LB, 
		numSpaceAvgCells_DSMC_X,numSpaceAvgCells_DSMC_Y,numSpaceAvgCells_DSMC,myRank,count);

              int rankToSendTo = p.nCoresX*(p.nCoresZ*myRankY + lowerCoupleLB) + myRankX;
              MPI_Isend(moments_DSMC_Send,numMoments*(XcellSizePerCore*YcellSizePerCore/numSpaceAvgCells_DSMC),MPI_DOUBLE,rankToSendTo,301,MPI_COMM_WORLD,&requestF1[0]);
          }

          if(myRankZ_LB == lowerCoupleLB)
          {
              /*Receive moments and calculate LB populations in the LB processor in the lower buffer layer*/
              int rankRecvFrom = p.nCoresX*(p.nCoresY*lowerCoupleDSMC + myRankY_LB) + myRankX_LB;
              MPI_Irecv(moments_DSMC_Recv,numMoments*(XcellSizePerCore*YcellSizePerCore/numSpaceAvgCells_DSMC),MPI_DOUBLE,rankRecvFrom,301,MPI_COMM_WORLD,&requestF1[0]);
              MPI_Wait(&requestF1[0], &statusF1[0]);

              calcRecvDSMCmomentsLowerBufferZone(p,lbModel41,gridLBM41,moments_DSMC_Recv,c,f2g_factor,dt,tau,DSMC2LBCouplingCellZ,LBbyDSMCpointsX,LBbyDSMCpointsY,numSpaceAvgCells_LB_X, 
		numSpaceAvgCells_LB_Y,numSpaceAvgCells_LB,numSpaceAvgCells_DSMC_X,numSpaceAvgCells_DSMC_Y,numSpaceAvgCells_DSMC,myRank,myRankX_LB,myRankY_LB, 
		myRankZ_LB,step);
          }

          if(myRankZ == upperCoupleDSMC)
          {
              /*Pack and send averaged DSMC moments in the upper buffer layer to the corresponding processors in the LB region*/
              calcSendDSMCmomentsUpperBufferZone(p,moments_DSMC_Send,DSMC2LBCouplingCellZ,LBbyDSMCpointsX,LBbyDSMCpointsY,numSpaceAvgCells_LB_X,numSpaceAvgCells_LB_Y,numSpaceAvgCells_LB, 
		numSpaceAvgCells_DSMC_X,numSpaceAvgCells_DSMC_Y,numSpaceAvgCells_DSMC,myRank,count);

              int rankToSendTo = p.nCoresX*(p.nCoresZ*myRankY + upperCoupleLB) + myRankX;
              MPI_Isend(moments_DSMC_Send,numMoments*(XcellSizePerCore*YcellSizePerCore/numSpaceAvgCells_DSMC),MPI_DOUBLE,rankToSendTo,303,MPI_COMM_WORLD,&requestF1[0]);
          }

          if(myRankZ_LB == upperCoupleLB)
          {
              /*Receive moments and calculate LB populations in the LB processor in the lower buffer layer*/
              int rankRecvFrom = p.nCoresX*(p.nCoresY*upperCoupleDSMC + myRankY_LB) + myRankX_LB;
              MPI_Irecv(moments_DSMC_Recv,numMoments*(XcellSizePerCore*YcellSizePerCore/numSpaceAvgCells_DSMC),MPI_DOUBLE,rankRecvFrom,303,MPI_COMM_WORLD,&requestF1[0]);
              MPI_Wait(&requestF1[0], &statusF1[0]);

              calcRecvDSMCmomentsUpperBufferZone(p,lbModel41,gridLBM41,moments_DSMC_Recv,c,f2g_factor,dt,tau,DSMC2LBCouplingCellZ,LBbyDSMCpointsX,LBbyDSMCpointsY,numSpaceAvgCells_LB_X, 
		numSpaceAvgCells_LB_Y,numSpaceAvgCells_LB,numSpaceAvgCells_DSMC_X,numSpaceAvgCells_DSMC_Y,numSpaceAvgCells_DSMC,myRank,myRankX_LB,myRankY_LB, 
		myRankZ_LB,step);
          }

          //collideWithForceSIMD(lbModel41,gridLBM41,VECT_LENGTH,twoBeta,tau,bodyForceX,bodyForceY,bodyForceZ,dt);

	  /*Entropic collision LB*/
          collideWithForceEntropicSIMD(lbModel41,gridLBM41,VECT_LENGTH,beta,betaTau,bodyForceX,bodyForceY,bodyForceZ,dt,step,myRank);

	  /*Packing for communication LB*/
          partial_recvCommunicate3DFace1(myRank, CartersianGrid3D.nList, requestF1,partial_temp1Mrec,partial_temp1Prec,m1,m2,m3,CartersianGrid3D);
          partial_packMSgTo1Neb(lbModel41, gridLBM41,partial_temp1Psend, partial_temp1Msend);
          partial_sendCommunicateFace1(myRank, CartersianGrid3D.nList, requestF1+2,partial_temp1Psend, partial_temp1Msend,m1,m2,m3,CartersianGrid3D);
          MPI_Waitall(4,requestF1,statusF1);
          partial_recvMSgFrom1Neb(lbModel41, gridLBM41,partial_temp1Prec,partial_temp1Mrec);

          partial_recvCommunicate3DFace2(myRank, CartersianGrid3D.nList, requestF2,partial_temp2Mrec,partial_temp2Prec,m1,m2,m3,CartersianGrid3D);
          partial_packMSgTo2Neb(lbModel41, gridLBM41,partial_temp2Psend, partial_temp2Msend);
          partial_sendCommunicateFace2(myRank, CartersianGrid3D.nList, requestF2+2,partial_temp2Psend, partial_temp2Msend,m1,m2,m3,CartersianGrid3D);
          MPI_Waitall(4,requestF2,statusF2);
          partial_recvMSgFrom2Neb(lbModel41, gridLBM41,partial_temp2Prec,partial_temp2Mrec);

          partial_recvCommunicate3DFace3(myRank, CartersianGrid3D.nList, requestF3,partial_temp3Mrec,partial_temp3Prec,m1,m2,m3,CartersianGrid3D);
          partial_packMSgTo3Neb(lbModel41, gridLBM41,partial_temp3Psend, partial_temp3Msend);
          partial_sendCommunicateFace3(myRank, CartersianGrid3D.nList, requestF3+2,partial_temp3Psend, partial_temp3Msend,m1,m2,m3,CartersianGrid3D);
          MPI_Waitall(4,requestF3,statusF3);
          partial_recvMSgFrom3Neb(lbModel41, gridLBM41,partial_temp3Prec,partial_temp3Mrec);

          int printFreq;
          // if((myRankZ >= lowerCoupleLB) && (myRankZ <= upperCoupleLB))
          {

	      /*LB advection*/
              advection(lbModel41, gridLBM41);

              printFreq = 10;

                    // if(step%500==0)
                    {
                     if(myRank==0)
                      timestamp();

		     /* Print Transverse kinetic energy LB*/
                     printTransverseKineticEnergy(lbModel41,gridLBM41,VECT_LENGTH,uRef,dt,beta,step,convectionTime,size,myRank,file1,ZPZ2H,
                     CartersianGrid3D.coordinates,bodyForceX,bodyForceY,bodyForceZ,p.nCoresY*p.nCoresX);
                     
                     printTransverseKineticEnergy(lbModel41,gridLBM41,VECT_LENGTH,uRef,dt,beta,step,convectionTime,size,myRank,file2,ZP25H,
                     CartersianGrid3D.coordinates,bodyForceX,bodyForceY,bodyForceZ,p.nCoresY*p.nCoresX);
                     
                     printTransverseKineticEnergy(lbModel41,gridLBM41,VECT_LENGTH,uRef,dt,beta,step,convectionTime,size,myRank,file3,ZP5H,
                     CartersianGrid3D.coordinates,bodyForceX,bodyForceY,bodyForceZ,p.nCoresY*p.nCoresX);
                     
                     printVelocityAtTheCenterOfDomain(lbModel41,gridLBM41,CartersianGrid3D.coordinates,bodyForceX,bodyForceY,bodyForceZ,centerX,centerY,centerZ,file4,step,convectionTime,dt);
                     globalMassTotalFluid(lbModel41,gridLBM41,marker,VECT_LENGTH,SOLID,FLUID,step,size,myRank);
                     // globalMassTotalSolid(lbModel41,gridLBM41,marker,VECT_LENGTH,SOLID,FLUID,step,size,myRank);
                     // PANINI_REAL domainAverage = 0.0;
                     // averageVelocity(lbModel41,gridLBM41,marker,VECT_LENGTH,SOLID,FLUID,step,size,myRank,bodyForceX,bodyForceY,bodyForceZ,uRef,domainAverage,dt);
                    }
          }



          MPI_Barrier(MPI_COMM_WORLD);

          if((tempCount%50==0))
          {
	    /* Print moments in the LB region*/
	    std::ofstream file24;
            char fileName[150];
            sprintf(fileName,"./moments/moments_node_%d_%d.dat",step,myRank);
            file24.open(fileName,std::ios::app);

            for (int localZ = gridLBM41.nB3; localZ <= gridLBM41.nE3; localZ++){
                for (int localY = gridLBM41.nB2; localY <= gridLBM41.nE2; localY++){
                        for (int localX = gridLBM41.nB1; localX <= gridLBM41.nE1; localX++){

                                int globalX = localX + simParam.n1*CartersianGrid3D.coordinates[0];
                                int globalY = localY + simParam.n2*CartersianGrid3D.coordinates[1];
                                int globalZ = localZ + simParam.n3*CartersianGrid3D.coordinates[2];

                    printMoments(lbModel41,gridLBM41,VECT_LENGTH,uRef,dt,beta,step,convectionTime,size,myRankZ_LB,file24,CartersianGrid3D.coordinates,bodyForceX,bodyForceY,bodyForceZ,centerX,centerY,centerZ,f2g_factor,p.deltaZ*LBbyDSMCpointsZ,p.height,count,localX,localY,localZ,globalX,globalY,globalZ);

                }
              }
            }

            file24.close();
            /******************************/

	    /*Print moments in the DSMC region*/
            std::ofstream file241;
            char fileName1[150];
            sprintf(fileName1,"./momentsLocal/moments_node1_%d.dat",step);
            file241.open(fileName1,std::ios::app);

            std::ofstream file242;
            char fileName2[150];
            sprintf(fileName2,"./momentsLocal/moments_node2_%d.dat",step);
            file242.open(fileName2,std::ios::app);


            if((myRankY_LB == ((int)(0.5*p.nCoresY))) && (myRankX_LB == ((int)(0.5*p.nCoresX))))
            {
            	/*Print local LB moments (probe line)*/
                printMomentsLocal(lbModel41,gridLBM41,VECT_LENGTH,uRef,dt,beta,step,convectionTime,size,myRankZ_LB,file241,
                CartersianGrid3D.coordinates,bodyForceX,bodyForceY,bodyForceZ,centerX,centerY,centerZ,f2g_factor,p.deltaZ*LBbyDSMCpointsZ,p.height,count);
            }

            if((myRankY_LB == ((int)(0.5*p.nCoresY))) && (myRankX_LB == 0))
            {
            	/*Print local LB moments (probe line)*/
                printMomentsLocal(lbModel41,gridLBM41,VECT_LENGTH,uRef,dt,beta,step,convectionTime,size,myRankZ_LB,file242,
                CartersianGrid3D.coordinates,bodyForceX,bodyForceY,bodyForceZ,centerX,centerY,centerZ,f2g_factor,p.deltaZ*LBbyDSMCpointsZ,p.height,count);
            }
            file241.close();
            file242.close();

            std::ofstream file24_dsmc;
            char fileName_dsmc[150];
            sprintf(fileName_dsmc,"./momentsLocal/moments_node_dsmc_%d.dat",step);
            file24_dsmc.open(fileName_dsmc,std::ios::app);



            if((myRankZ == lowerCoupleDSMC) && (myRankY == ((int)(0.5*p.nCoresY))) && (myRankX == ((int)(0.5*p.nCoresX))))
            {
                 int cellX                     = (int)(0.5*XcellSizePerCore);
                 int cellY                     = (int)(0.5*YcellSizePerCore);
                 int cellNumberLocal           = ((int)(XcellSizePerCore/numSpaceAvgCells_DSMC_X))*((int)(cellY/numSpaceAvgCells_DSMC_Y)) + ((int)(cellX/numSpaceAvgCells_DSMC_X));

                 for(int Zcell = 0; Zcell < p.ZcellSizePerCore; Zcell++)
                 {
                   int cellNumberLocalDSMCcouple = Zcell*(XcellSizePerCore/numSpaceAvgCells_DSMC_X)*(YcellSizePerCore/numSpaceAvgCells_DSMC_Y) + cellNumberLocal;

                   double rho_DSMC      = p.moments_Average[numMoments*cellNumberLocalDSMCcouple     ];
                   double p_DSMC        = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  1];
                   double t0_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  1]/p.moments_Average[numMoments*cellNumberLocalDSMCcouple     ];
                   double uX_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  2];
                   double uY_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  3];
                   double uZ_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  4];
                   double sigmaXX_DSMC  = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  5];
                   double sigmaXY_DSMC  = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  6];
                   double sigmaZX_DSMC  = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  7];
                   double sigmaYY_DSMC  = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  8];
                   double sigmaYZ_DSMC  = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  9];
                   double sigmaZZ_DSMC  = p.moments_Average[numMoments*cellNumberLocalDSMCcouple + 10];
                   double qX_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple + 11];
                   double qY_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple + 12];
                   double qZ_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple + 13];

                   file24_dsmc<<"\n"<<(Zcell + 0.5)*p.deltaZ<<"\t"<<rho_DSMC<<"\t"<<t0_DSMC<<"\t"<<p_DSMC<<" "<<uY_DSMC<<"\t"<<uZ_DSMC<<"\t"<<uX_DSMC<<"\t"<<sigmaYY_DSMC<<"\t"<<sigmaZZ_DSMC<<"\t"<<sigmaXX_DSMC<<"\t"<<sigmaYZ_DSMC<<"\t"<<sigmaZX_DSMC<<"\t"<<sigmaXY_DSMC<<"\t"<<qY_DSMC<<"\t"<<qZ_DSMC<<"\t"<<qX_DSMC<<std::endl;
                   std::cout<<"Printed here1: "<<count<<" "<<myRankZ<<" "<<Zcell<<"  "<<Zcell*p.deltaZ<<"  "<<rho_DSMC<<"  "<<t0_DSMC<<"  "<<p_DSMC<<" "<<uY_DSMC<<"  "<<uZ_DSMC<<"  "<<uX_DSMC<<"  "<<sigmaYY_DSMC<<"  "<<sigmaZZ_DSMC<<"  "<<sigmaXX_DSMC<<"  "<<sigmaYZ_DSMC<<"  "<<sigmaZX_DSMC<<"  "<<sigmaXY_DSMC<<"  "<<std::endl;
                 }
             }


             if((myRankZ == upperCoupleDSMC) && (myRankY == ((int)(0.5*p.nCoresY))) && (myRankX == ((int)(0.5*p.nCoresX))))
             {
                 int cellX                     = (int)(0.5*XcellSizePerCore);
                 int cellY                     = (int)(0.5*YcellSizePerCore);
                 int cellNumberLocal           = ((int)(XcellSizePerCore/numSpaceAvgCells_DSMC_X))*((int)(cellY/numSpaceAvgCells_DSMC_Y)) + ((int)(cellX/numSpaceAvgCells_DSMC_X));

                 for(int Zcell = 0; Zcell < p.ZcellSizePerCore; Zcell++)
                 {
                   int cellNumberLocalDSMCcouple = Zcell*(XcellSizePerCore/numSpaceAvgCells_DSMC_X)*(YcellSizePerCore/numSpaceAvgCells_DSMC_Y) + cellNumberLocal;

                   double rho_DSMC      = p.moments_Average[numMoments*cellNumberLocalDSMCcouple     ];
                   double p_DSMC        = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  1];
                   double t0_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  1]/p.moments_Average[numMoments*cellNumberLocalDSMCcouple     ];
                   double uX_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  2];
                   double uY_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  3];
                   double uZ_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  4];
                   double sigmaXX_DSMC  = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  5];
                   double sigmaXY_DSMC  = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  6];
                   double sigmaZX_DSMC  = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  7];
                   double sigmaYY_DSMC  = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  8];
                   double sigmaYZ_DSMC  = p.moments_Average[numMoments*cellNumberLocalDSMCcouple +  9];
                   double sigmaZZ_DSMC  = p.moments_Average[numMoments*cellNumberLocalDSMCcouple + 10];
                   double qX_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple + 11];
                   double qY_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple + 12];
                   double qZ_DSMC       = p.moments_Average[numMoments*cellNumberLocalDSMCcouple + 13];

                   file24_dsmc<<"\n"<<(p.height - (Zcell + 0.5)*p.deltaZ)<<"\t"<<rho_DSMC<<"\t"<<t0_DSMC<<"\t"<<p_DSMC<<" "<<uY_DSMC<<"\t"<<uZ_DSMC<<"\t"<<uX_DSMC<<"\t"<<sigmaYY_DSMC<<"\t"<<sigmaZZ_DSMC<<"\t"<<sigmaXX_DSMC<<"\t"<<sigmaYZ_DSMC<<"\t"<<sigmaZX_DSMC<<"\t"<<sigmaXY_DSMC<<"\t"<<qY_DSMC<<"\t"<<qZ_DSMC<<"\t"<<qX_DSMC<<std::endl;
                   std::cout<<"Printed here2: "<<count<<" "<<myRankZ<<" "<<Zcell<<"  "<<(p.height - Zcell*p.deltaZ)<<"  "<<rho_DSMC<<"  "<<t0_DSMC<<"  "<<p_DSMC<<" "<<uY_DSMC<<"  "<<uZ_DSMC<<"  "<<uX_DSMC<<"  "<<sigmaYY_DSMC<<"  "<<sigmaZZ_DSMC<<"  "<<sigmaXX_DSMC<<"  "<<sigmaYZ_DSMC<<"  "<<sigmaZX_DSMC<<"  "<<sigmaXY_DSMC<<"  "<<std::endl;
               }
             }

            MPI_Barrier(MPI_COMM_WORLD);

            file24.close();
            file24_dsmc.close();
            /************************************/


	    /*Dump visualization files (VTK)*/
            if((step%(20*LB_Count)==0))
            {
                dumpResultsNodeToPVTI(lbModel41,gridLBM41,CartersianGrid3D,procRank,simParam.n1,simParam.n2,simParam.n3,simParam.nP1,simParam.nP2,simParam.nP3,step,bodyForceX,bodyForceY,bodyForceZ,dt);
                if(myRank==0)
                    writeResultsNodeMasterPVTI(simParam.n1,simParam.n2,simParam.n3,simParam.nP1,simParam.nP2,simParam.nP3,step);
            }
            /*******************************/
          }


          // auto begin10 = std::chrono::high_resolution_clock::now();
          // time8 += std::chrono::duration_cast<std::chrono::nanoseconds>(begin10-begin9).count();



          auto begin11 = std::chrono::high_resolution_clock::now();
          // time9 += std::chrono::duration_cast<std::chrono::nanoseconds>(begin11-begin10).count();
          //
          timeLB += std::chrono::duration_cast<std::chrono::nanoseconds>(begin11-begin8).count();


          step++;
          /*End of m LB steps and correspinding n DSMC steps*/

      }


	    /*Store DSMC particle data and LB population data for entire grid for restart*/
            p.storeParticleDataToFile(myRank);
            dumpRestart(lbModel41, gridLBM41,LBrestartFile);


      //Clearing DSMC Moments
      {
      	  /*Clear DSMC averaged moments*/
          p.DSMCclearAvgQts();
      }

    }
    /*End of Simulation*/

    std::cout<<"Total_time: "<<myRank<<"\t"<<timeDSMC<<"\t"<<timeLB<<"\t"<<(timeDSMC)+(timeLB)<<std::endl;

    fileTotalDensity.close();

    if(rank==0)
    {
        std::cout<<"Executed! ";
        timestamp();
    }

    MPI_Finalize();


}
