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
//
// #include"bounceBack41.h"
#include"bounceBack41_withForce.h"
#include"diffuseBounceBack41_withForce.h"
#include"diffuse41.h"

// #include"ioStreamGrid.h"
//
#include"perlinNoise.h"
#include"noiseFunctions.h"
#include<packNoise41.h>
#include<communicationNoise.h>
#include"pvts.h"
#include"pvts.h"
#include"momentsForDSMC.h"

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

// void allocate_Noise(int m1, int m2, int m3)
// {
//   temp1Psend_Noise = allocateForPacking_Noise(m2,m3);
//   temp1Prec_Noise  = allocateForPacking_Noise(m2,m3);
//   temp2Psend_Noise = allocateForPacking_Noise(m1,m3);
//   temp2Prec_Noise  = allocateForPacking_Noise(m1,m3);
//   temp3Psend_Noise = allocateForPacking_Noise(m1,m2);
//   temp3Prec_Noise  = allocateForPacking_Noise(m1,m2);
//   temp1Msend_Noise = allocateForPacking_Noise(m2,m3);
//   temp1Mrec_Noise  = allocateForPacking_Noise(m2,m3);
//   temp2Msend_Noise = allocateForPacking_Noise(m1,m3);
//   temp2Mrec_Noise  = allocateForPacking_Noise(m1,m3);
//   temp3Msend_Noise = allocateForPacking_Noise(m1,m2);
//   temp3Mrec_Noise  = allocateForPacking_Noise(m1,m2);
// }



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


/**
 * @mainpage
 *
 * This is the documentation of the <em> <b> Higher order crystallographic lattice Boltzmann model<em> <b> with 41 velocities.
 */
int main()
{

  int myRank, size;
  MPI::Init();
  myRank = MPI::COMM_WORLD.Get_rank() ;
  size   = MPI::COMM_WORLD.Get_size() ;

  latticeInfo simParam(48,32,4,1,6,1,2,2,2);
//   double Lx(7.0),Ly(2.0),Lz(3.5);

  lbmRD3Q41<PANINI_REAL> lbModel41(1.0);
  gridBCC3D<4,11,PANINI_REAL>  gridLBM41(simParam.n1,simParam.n2,simParam.n3, simParam.nGhost1,simParam.nGhost2,simParam.nGhost3);
  gridBCC3D<4,1 ,PANINI_REAL>   testGrid(simParam.n1,simParam.n2,simParam.n3, simParam.nGhost1,simParam.nGhost2,simParam.nGhost3);
  gridBCC3D<1,1 ,PANINI_INT >     marker(simParam.n1,simParam.n2,simParam.n3, simParam.nGhost1,simParam.nGhost2,simParam.nGhost3);

  PANINI_INT SOLID = 1;
  PANINI_INT FLUID = 2;

  setModelParameters(lbModel41);
  getLatticeParameter(lbModel41);

  m1 = gridLBM41.n1; m2 = gridLBM41.n2; m3 = gridLBM41.n3;
  allocate(m1,m2,m3);
  partial_allocate(m1,m2,m3);
//   allocate_Noise(m1,m2,m3);

  procInfo3D  CartersianGrid3D(simParam.nP1,simParam.nP2,simParam.nP3, true, true, true, true);
  CartersianGrid3D.getnList();
  MPI::Request requestF1[4];
  MPI::Request requestF2[4];
  MPI::Request requestF3[4];
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

                                                                                                            //
                                                                                                            //  temperature   : 0.5
                                                                                                            //  Re            : 100
                                                                                                            //  Ma            : 0.2
  PANINI_REAL cs2      = 5.0*0.5/3.0;//c*c*lbModel41.theta0;                                                //  L             : 2000
  PANINI_REAL c        = 1.68103;                                                        //  cs_DSMC       : 0.912871
  PANINI_REAL cs       = sqrt(cs2);                                                                         //  u             : 0.182574
  PANINI_REAL Ma       = 0.2;                                                                               //  nu            : 3.65148
  PANINI_REAL refLen   = 2000.0;//(simParam.nP2*simParam.n2)+1.0;                                           //  Kn            : 0.00255323
  PANINI_REAL u_inlet  = Ma*cs;                                                                             //  lambda        : 5.10646
  const PANINI_REAL Re = 50.0;                                                                             //  dt_DSMC       : 7.22163
  PANINI_REAL uRef     = u_inlet;                                                                           //
  PANINI_REAL kinVisc  = uRef*refLen/Re;
  PANINI_REAL tau      = kinVisc/cs2;
  PANINI_REAL dx       = refLen/(simParam.n2*simParam.nP2);
  PANINI_REAL dt       = dx/c;                                                                              //  N             : 200
  PANINI_REAL tauNdim  = tau/dt;                                                                            //  dx_LB         : 10
  PANINI_REAL beta     = 1.0/(2.0*tauNdim + 1.0);                                                           //  theta0        : 0.294896
  PANINI_REAL twoBeta  = 2.0*beta;                                                                          //  c             : 1.68103
  PANINI_REAL betaTau  = beta*tau;                                                                          //  dt_LB         : 5.94875
  PANINI_REAL oneMinustwoBeta = 1.0-twoBeta;                                                                //  tau           : 4.38178
  PANINI_REAL convectionTime = (refLen*refLen)/kinVisc;                                                     //  beta          : 0.404338
  PANINI_REAL bodyConvectionTime =  refLen/u_inlet;                                                         //
  int VECT_LENGTH = 4;                                                                                      //
  PANINI_REAL forceNodeBefore[3],forceCellBefore[3];                                                        //
  PANINI_REAL forceNodeAfter[3],forceCellAfter[3];                                                          //  Re            : 100
  PANINI_REAL NetForce[3];                                                                                  //  Ma            : 0.2
  PANINI_REAL ForceX(0.0),ForceY(0.0),ForceZ(0.0),cDOLD(0.0),cD(0.0),cLOLD(0.0),cL(0.0),cZOLD(0.0),cZ(0.0); //  Kn            : 0.00255323
  PANINI_REAL uXBody = 0.0;                                                                                 //  N             : 200
  PANINI_REAL uYBody = 0.0;                                                                                 //  L             : 2000
  PANINI_REAL uZBody = 0.0;                                                                                 //  dt_DSMC       : 7.22163
                                                                                                            //  dt_LB         : 5.94875
                                                                                                            //  dt_LB/dt_DSMC : 0.82374


/*
                                                                                                            //
                                                                                                            //  temperature   : 0.5
                                                                                                            //  Re            : 100
                                                                                                            //  Ma            : 0.2
  PANINI_REAL cs2      = lbModel41.theta0;                                                //  L             : 2000
  PANINI_REAL c        = 1.0;//sqrt(cs2/lbModel41.theta0);                                                        //  cs_DSMC       : 0.912871
  PANINI_REAL cs       = sqrt(cs2);                                                                         //  u             : 0.182574
  PANINI_REAL Ma       = 0.2;                                                                               //  nu            : 3.65148
  PANINI_REAL refLen   = (simParam.nP2*simParam.n2) ;                                           //  Kn            : 0.00255323
  PANINI_REAL u_inlet  = Ma*cs;                                                                             //  lambda        : 5.10646
  const PANINI_REAL Re = 100.0;                                                                             //  dt_DSMC       : 7.22163
  PANINI_REAL uRef     = u_inlet;                                                                           //
  PANINI_REAL kinVisc  = uRef*refLen/Re;
  PANINI_REAL tau      = kinVisc/cs2;
//   PANINI_REAL dx       = refLen/simParam.n2;
  PANINI_REAL dt       = 1.0;//dx/c;                                                                              //  N             : 200
  PANINI_REAL tauNdim  = tau/dt;                                                                            //  dx_LB         : 10
  PANINI_REAL beta     = 1.0/(2.0*tauNdim + 1.0);                                                           //  theta0        : 0.294896
  PANINI_REAL twoBeta  = 2.0*beta;                                                                          //  c             : 1.68103
  PANINI_REAL betaTau  = beta*tau;                                                                          //  dt_LB         : 5.94875
  PANINI_REAL oneMinustwoBeta = 1.0-twoBeta;                                                                //  tau           : 4.38178
  PANINI_REAL convectionTime = refLen*refLen/kinVisc;                                                     //  beta          : 0.404338
  PANINI_REAL bodyConvectionTime =  refLen/u_inlet;                                                         //
  int VECT_LENGTH = 4;                                                                                      //
  PANINI_REAL forceNodeBefore[3],forceCellBefore[3];                                                        //
  PANINI_REAL forceNodeAfter[3],forceCellAfter[3];                                                          //  Re            : 100
  PANINI_REAL NetForce[3];                                                                                  //  Ma            : 0.2
  PANINI_REAL ForceX(0.0),ForceY(0.0),ForceZ(0.0),cDOLD(0.0),cD(0.0),cLOLD(0.0),cL(0.0),cZOLD(0.0),cZ(0.0); //  Kn            : 0.00255323
  PANINI_REAL uXBody = 0.0;                                                                                 //  N             : 200
  PANINI_REAL uYBody = 0.0;                                                                                 //  L             : 2000
  PANINI_REAL uZBody = 0.0;                                                                                 //  dt_DSMC       : 7.22163
                                                                                                            //  dt_LB         : 5.94875
                                                                                                            //  dt_LB/dt_DSMC : 0.82374


            */



  int centerX =  (int)(0.5*simParam.n1*simParam.nP1) + gridLBM41.nB1;
  int centerY =  (int)(0.5*simParam.n2*simParam.nP2) + gridLBM41.nB2;
  int centerZ =  (int)(0.5*simParam.n3*simParam.nP3) + gridLBM41.nB3;

  PANINI_REAL ForceXOld(0.0),ForceX_global(0.0);
  PANINI_REAL ForceYOld(0.0),ForceY_global(0.0);
  PANINI_REAL ForceZOld(0.0),ForceZ_global(0.0);

  PANINI_REAL bodyForceX = 8.0*uRef*kinVisc/(refLen*refLen); //1.12598889221483e-08;
  PANINI_REAL bodyForceY = 0.0;
  PANINI_REAL bodyForceZ = 0.0;
  PANINI_REAL idealBodyForce = bodyForceX;// = 8.0*uRef*kinVisc/(refLen*refLen); //1.12598889221483e-08;

  int restartTimeint = (int)(2500.0*convectionTime);

  if(myRank==0)
  {
    std::cout<<"Dimensions          : "<<gridLBM41.m1*simParam.nP1<<"x"<<gridLBM41.m2*simParam.nP2<<"x"<<gridLBM41.m3*simParam.nP3<<"\n";
    std::cout<<"Dimensions          : "<<gridLBM41.m1<<"*"<<simParam.nP1<<"x"<<gridLBM41.m2<<"*"<<simParam.nP2<<"x"<<gridLBM41.m3<<"*"<<simParam.nP3<<"\n";
    std::cout<<"Re                  : "<<Re<<"\n";
    std::cout<<"U_O                 : "<<u_inlet<<"\n";
    std::cout<<"Ma                  : "<<Ma<<"\n";
    std::cout<<"beta                : "<<beta<<"\n";
    std::cout<<"dt                  : "<<dt<<"\n";
    std::cout<<"convectionTime      : "<<convectionTime<<"\n";
    std::cout<<"convectionTime Body : "<<bodyConvectionTime<<"\n";
    std::cout<<"body force          : "<<bodyForceX<<"\n";
    std::cout<<"c                   : "<<c<<"\n";
 }


  initializeMarker(marker,FLUID);

//   markTopAndBottomWalls      (marker,FLUID,SOLID,CartersianGrid3D.coordinates,simParam.nP1,simParam.nP2,simParam.nP3);
//   markOnePercentBumpOnTheWall(marker,FLUID,SOLID,CartersianGrid3D.coordinates,simParam.nP1,simParam.nP2,simParam.nP3,simParam.n1,simParam.n2,simParam.n3);
  initialConditions(lbModel41,gridLBM41,marker,VECT_LENGTH,u_inlet,FLUID,simParam.nP2,CartersianGrid3D.coordinates);

  PANINI_REAL normg;

//   curl_noise_communicateAfterThis(lbModel41,testGrid,simParam.nP1,simParam.nP2,simParam.nP3,normg,CartersianGrid3D.coordinates,VECT_LENGTH,Lx,Ly,Lz) ;
//
//
//   recvCommunicate3DFace1_Noise(myRank, CartersianGrid3D.nList, requestF1,temp1Mrec_Noise,temp1Prec_Noise,m1,m2,m3,CartersianGrid3D);
//   packMSgTo1Neb_Noise(lbModel41, testGrid,temp1Psend_Noise, temp1Msend_Noise);
//   sendCommunicateFace1_Noise(myRank, CartersianGrid3D.nList, requestF1+2,temp1Psend_Noise, temp1Msend_Noise,m1,m2,m3,CartersianGrid3D);
//   MPI::Request::Waitall(4,requestF1);
//   recvMSgFrom1Neb_Noise(lbModel41, testGrid,temp1Prec_Noise,temp1Mrec_Noise);
//
//   recvCommunicate3DFace2_Noise(myRank, CartersianGrid3D.nList, requestF2,temp2Mrec_Noise,temp2Prec_Noise,m1,m2,m3,CartersianGrid3D);
//   packMSgTo2Neb_Noise(lbModel41, testGrid,temp2Psend_Noise, temp2Msend_Noise);
//   sendCommunicateFace2_Noise(myRank, CartersianGrid3D.nList, requestF2+2,temp2Psend_Noise, temp2Msend_Noise,m1,m2,m3,CartersianGrid3D);
//   MPI::Request::Waitall(4,requestF2);
//   recvMSgFrom2Neb_Noise(lbModel41, testGrid,temp2Prec_Noise,temp2Mrec_Noise);
//
//   recvCommunicate3DFace3_Noise(myRank, CartersianGrid3D.nList, requestF3,temp3Mrec_Noise,temp3Prec_Noise,m1,m2,m3,CartersianGrid3D);
//   packMSgTo3Neb_Noise(lbModel41, testGrid,temp3Psend_Noise, temp3Msend_Noise);
//   sendCommunicateFace3_Noise(myRank, CartersianGrid3D.nList, requestF3+2,temp3Psend_Noise, temp3Msend_Noise,m1,m2,m3,CartersianGrid3D);
//   MPI::Request::Waitall(4,requestF3);
//   recvMSgFrom3Neb_Noise(lbModel41, testGrid,temp3Prec_Noise,temp3Mrec_Noise);
//
//   curl_noise_communiateBeforeThis(lbModel41,gridLBM41,testGrid,simParam.nP1,simParam.nP2,simParam.nP3,normg,CartersianGrid3D.coordinates,VECT_LENGTH,Lx,Ly,Lz) ;
//


//  curl_noise_communiateBeforeThis(lbModel41,gridLBM41,testGrid,simParam.nP1,simParam.nP2,simParam.nP3,normg,CartersianGrid3D.coordinates,VECT_LENGTH,Lx,Ly,Lz) ;
// curl_noise(lbModel41,gridLBM41,simParam.nP1,simParam.nP2,simParam.nP3,normg,CartersianGrid3D.coordinates,VECT_LENGTH);

//  perturb_parab(CartersianGrid3D.coordinates,lbModel41,gridLBM41,simParam.nP1,simParam.nP2,simParam.nP3,normg,u_inlet,VECT_LENGTH);

  if(myRank==0)
   std::cout<<"normg: "<<normg<<"\n";

  //   //////////////////////////
  //   // READ FROM RESTART FILE
  //   //////////////////////////
  //   Rcount = 6;
  //   sprintf(restartFile,"./results/restart_%d_%d.dat",Rcount,procRank);
  //   restartIC(lbModel41,gridLBM41,restartFile);
  //
  //   startTime = ((Rcount+1)*(int)) + 1;


  recvCommunicate3DFace1(myRank, CartersianGrid3D.nList, requestF1,temp1Mrec,temp1Prec,m1,m2,m3,CartersianGrid3D);
  packMSgTo1Neb(lbModel41, gridLBM41,temp1Psend, temp1Msend);
  sendCommunicateFace1(myRank, CartersianGrid3D.nList, requestF1+2,temp1Psend, temp1Msend,m1,m2,m3,CartersianGrid3D);
  MPI::Request::Waitall(4,requestF1);
  recvMSgFrom1Neb(lbModel41, gridLBM41,temp1Prec,temp1Mrec);

  recvCommunicate3DFace2(myRank, CartersianGrid3D.nList, requestF2,temp2Mrec,temp2Prec,m1,m2,m3,CartersianGrid3D);
  packMSgTo2Neb(lbModel41, gridLBM41,temp2Psend, temp2Msend);
  sendCommunicateFace2(myRank, CartersianGrid3D.nList, requestF2+2,temp2Psend, temp2Msend,m1,m2,m3,CartersianGrid3D);
  MPI::Request::Waitall(4,requestF2);
  recvMSgFrom2Neb(lbModel41, gridLBM41,temp2Prec,temp2Mrec);

  recvCommunicate3DFace3(myRank, CartersianGrid3D.nList, requestF3,temp3Mrec,temp3Prec,m1,m2,m3,CartersianGrid3D);
  packMSgTo3Neb(lbModel41, gridLBM41,temp3Psend, temp3Msend);
  sendCommunicateFace3(myRank, CartersianGrid3D.nList, requestF3+2,temp3Psend, temp3Msend,m1,m2,m3,CartersianGrid3D);
  MPI::Request::Waitall(4,requestF3);
  recvMSgFrom3Neb(lbModel41, gridLBM41,temp3Prec,temp3Mrec);

  dumpResultsNodeToPVTI(lbModel41,gridLBM41,CartersianGrid3D,procRank,simParam.n1,simParam.n2,simParam.n3,simParam.nP1,simParam.nP2,simParam.nP3,0,0.0,0.0,0.0,dt);
  if(myRank==0)
    writeResultsNodeMasterPVTI(simParam.n1,simParam.n2,simParam.n3,simParam.nP1,simParam.nP2,simParam.nP3,0);

  if(myRank==0)
    timestamp();
  globalMassTotalFluid(lbModel41,gridLBM41,marker,VECT_LENGTH,SOLID,FLUID,0,size,myRank);
  globalMassTotalSolid(lbModel41,gridLBM41,marker,VECT_LENGTH,SOLID,FLUID,0,size,myRank);

  int simulationTime = int(1.1*convectionTime);


  PANINI_REAL domainAverage = 0.0;


//       averageVelocity(lbModel41,gridLBM41,marker,VECT_LENGTH,SOLID,FLUID,0,size,myRank,bodyForceX,bodyForceY,bodyForceZ,uRef,domainAverage);


//     std::cout<<myRank<<"  "<<domainAverage<<std::endl;

//      if(domainAverage<(2.0/3.0*uRef))
//       bodyForceX += 0.01*bodyForceX;
//
//      if(std::fabs(domainAverage-(2.0/3.0*uRef))<1e-3)
//       bodyForceX = idealBodyForce;
//
//      if(myRank==0)
//          std::cout<<0<<std::setw(24)<<std::setprecision(14)<<bodyForceX<<std::setw(24)<<domainAverage<<" Force "<<std::endl;



  for(int step=1;step<=simulationTime;step++)
  {



//     collideWithForceEntropicSIMD(lbModel41,gridLBM41,VECT_LENGTH,beta,betaTau,bodyForceX,bodyForceY,bodyForceZ);
    collideWithForceSIMD(lbModel41,gridLBM41,VECT_LENGTH,twoBeta,tau,bodyForceX,bodyForceY,bodyForceZ,dt);

    partial_recvCommunicate3DFace1(myRank, CartersianGrid3D.nList, requestF1,partial_temp1Mrec,partial_temp1Prec,m1,m2,m3,CartersianGrid3D);
    partial_packMSgTo1Neb(lbModel41, gridLBM41,partial_temp1Psend, partial_temp1Msend);
    partial_sendCommunicateFace1(myRank, CartersianGrid3D.nList, requestF1+2,partial_temp1Psend, partial_temp1Msend,m1,m2,m3,CartersianGrid3D);
    MPI::Request::Waitall(4,requestF1);
    partial_recvMSgFrom1Neb(lbModel41, gridLBM41,partial_temp1Prec,partial_temp1Mrec);



    partial_recvCommunicate3DFace2(myRank, CartersianGrid3D.nList, requestF2,partial_temp2Mrec,partial_temp2Prec,m1,m2,m3,CartersianGrid3D);
    partial_packMSgTo2Neb(lbModel41, gridLBM41,partial_temp2Psend, partial_temp2Msend);
    partial_sendCommunicateFace2(myRank, CartersianGrid3D.nList, requestF2+2,partial_temp2Psend, partial_temp2Msend,m1,m2,m3,CartersianGrid3D);
    MPI::Request::Waitall(4,requestF2);
    partial_recvMSgFrom2Neb(lbModel41, gridLBM41,partial_temp2Prec,partial_temp2Mrec);

    partial_recvCommunicate3DFace3(myRank, CartersianGrid3D.nList, requestF3,partial_temp3Mrec,partial_temp3Prec,m1,m2,m3,CartersianGrid3D);
    partial_packMSgTo3Neb(lbModel41, gridLBM41,partial_temp3Psend, partial_temp3Msend);
    partial_sendCommunicateFace3(myRank, CartersianGrid3D.nList, requestF3+2,partial_temp3Psend, partial_temp3Msend,m1,m2,m3,CartersianGrid3D);
    MPI::Request::Waitall(4,requestF3);
    partial_recvMSgFrom3Neb(lbModel41, gridLBM41,partial_temp3Prec,partial_temp3Mrec);


    if(CartersianGrid3D.coordinates[1]==0)
     prepareDiffuseBottomWall(lbModel41, gridLBM41);

    if(CartersianGrid3D.coordinates[1]==simParam.nP2-1)
     prepareDiffuseTopWall   (lbModel41, gridLBM41);

    advection(lbModel41, gridLBM41);

    if(CartersianGrid3D.coordinates[1]==0)
      correctDiffuseF0massManipBottomWall(lbModel41, gridLBM41,VECT_LENGTH,lbModel41.theta0,0.0);
//
    if(CartersianGrid3D.coordinates[1]==simParam.nP2-1)
      correctDiffuseF0massManipTopWall   (lbModel41, gridLBM41,VECT_LENGTH,lbModel41.theta0,0.0);
//
//
// //     bounceBackWithForces(lbModel41, gridLBM41,marker,FLUID,ForceX,ForceY,ForceZ,step);
//
//
    averageVelocity(lbModel41,gridLBM41,marker,VECT_LENGTH,SOLID,FLUID,step,size,myRank,bodyForceX,bodyForceY,bodyForceZ,uRef,domainAverage,dt);


//     std::cout<<myRank<<"  "<<domainAverage<<std::endl;

    //     if((step>5.0*bodyConvectionTime) &&(step%1000==0))
    //     {
    //      if(domainAverage<(2.0/3.0*uRef))
    //       bodyForceX += 0.001*bodyForceX;
    //
    //      if(domainAverage>(2.0/3.0*uRef))
    //       bodyForceX -= 0.001*bodyForceX;
    //
    //      if(std::fabs(domainAverage-(2.0/3.0*uRef))<1e-3)
    //       bodyForceX = idealBodyForce;
    //
    //      if(myRank==0)
    //          std::cout<<step<<std::setw(24)<<std::setprecision(14)<<bodyForceX<<std::setw(24)<<domainAverage<<" Force "<<std::endl;
    //     }
    //
    if(step%500==0)
    {
     if(myRank==0)
      timestamp();

     printTransverseKineticEnergy(lbModel41,gridLBM41,VECT_LENGTH,uRef,dt,beta,step,convectionTime,size,myRank,file1,ZPZ2H,CartersianGrid3D.coordinates,bodyForceX,bodyForceY,bodyForceZ);
     printTransverseKineticEnergy(lbModel41,gridLBM41,VECT_LENGTH,uRef,dt,beta,step,convectionTime,size,myRank,file2,ZP25H,CartersianGrid3D.coordinates,bodyForceX,bodyForceY,bodyForceZ);
     printTransverseKineticEnergy(lbModel41,gridLBM41,VECT_LENGTH,uRef,dt,beta,step,convectionTime,size,myRank,file3,ZP5H ,CartersianGrid3D.coordinates,bodyForceX,bodyForceY,bodyForceZ);
     printVelocityAtTheCenterOfDomain(lbModel41,gridLBM41,CartersianGrid3D.coordinates,bodyForceX,bodyForceY,bodyForceZ,centerX,centerY,centerZ,file4,step,convectionTime,dt);
     globalMassTotalFluid(lbModel41,gridLBM41,marker,VECT_LENGTH,SOLID,FLUID,step,size,myRank);
     globalMassTotalSolid(lbModel41,gridLBM41,marker,VECT_LENGTH,SOLID,FLUID,step,size,myRank);
     PANINI_REAL domainAverage = 0.0;
     
     

    }
    
    
    
    if(step%500==0)
    {
     std::ofstream file;
     char fileName[150];
     sprintf(fileName,"./moments/moments_node_%d.dat",step);
     file.open(fileName,std::ios::app);
     
     printMoments(lbModel41,gridLBM41,VECT_LENGTH,uRef,dt,beta,step,convectionTime,size,myRank,file,CartersianGrid3D.coordinates,bodyForceX,bodyForceY,bodyForceZ,centerX,centerY,centerZ);

     MPI_Barrier(MPI::COMM_WORLD);
     file.close();
    }

    if( (step%(int)(0.1*bodyConvectionTime)==0) )
    {
      recvCommunicate3DFace1(myRank, CartersianGrid3D.nList, requestF1,temp1Mrec,temp1Prec,m1,m2,m3,CartersianGrid3D);
      packMSgTo1Neb(lbModel41, gridLBM41,temp1Psend, temp1Msend);
      sendCommunicateFace1(myRank, CartersianGrid3D.nList, requestF1+2,temp1Psend, temp1Msend,m1,m2,m3,CartersianGrid3D);
      MPI::Request::Waitall(4,requestF1);
      recvMSgFrom1Neb(lbModel41, gridLBM41,temp1Prec,temp1Mrec);

      recvCommunicate3DFace2(myRank, CartersianGrid3D.nList, requestF2,temp2Mrec,temp2Prec,m1,m2,m3,CartersianGrid3D);
      packMSgTo2Neb(lbModel41, gridLBM41,temp2Psend, temp2Msend);
      sendCommunicateFace2(myRank, CartersianGrid3D.nList, requestF2+2,temp2Psend, temp2Msend,m1,m2,m3,CartersianGrid3D);
      MPI::Request::Waitall(4,requestF2);
      recvMSgFrom2Neb(lbModel41, gridLBM41,temp2Prec,temp2Mrec);

      recvCommunicate3DFace3(myRank, CartersianGrid3D.nList, requestF3,temp3Mrec,  temp3Prec,m1,m2,m3,CartersianGrid3D);
      packMSgTo3Neb(lbModel41, gridLBM41,temp3Psend, temp3Msend);
      sendCommunicateFace3(myRank, CartersianGrid3D.nList, requestF3+2,temp3Psend, temp3Msend,m1,m2,m3,CartersianGrid3D);
      MPI::Request::Waitall(4,requestF3);
      recvMSgFrom3Neb(lbModel41, gridLBM41,temp3Prec,temp3Mrec);

      dumpResultsNodeToPVTI(lbModel41,gridLBM41,CartersianGrid3D,procRank,simParam.n1,simParam.n2,simParam.n3,simParam.nP1,simParam.nP2,simParam.nP3,step,bodyForceX,bodyForceY,bodyForceZ,dt);
      if(myRank==0)
       writeResultsNodeMasterPVTI(simParam.n1,simParam.n2,simParam.n3,simParam.nP1,simParam.nP2,simParam.nP3,step);
    }

//     if( step%restartTimeint==0 )
//     {
//       //////////////////////////
//       // WRITE A RESTART FILE
//       //////////////////////////
//       sprintf(restartFile,"./results/restart_%d_%d.dat",Rcount,procRank);
//       dumpRestart(lbModel41,gridLBM41,restartFile);
//       Rcount++;
//     }

  }
  MPI::Finalize();
}
