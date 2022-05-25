#include<momentsForDSMC.h>


template <int N,int numblock, typename dataType1>
void copyFromNodeToArray(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int i1, int i2, int i3, dataType1 *fTemp)
{
   fTemp[lbModel.DV_ZERO_ZERO_ZERO] = myGrid(i1,i2,i3,lbModel.G0 ,myGrid.node[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO);

   fTemp[lbModel.DV_ZERO_ZERO_P1  ] = myGrid(i1,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P1      );
   fTemp[lbModel.DV_ZERO_ZERO_P2  ] = myGrid(i1,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P2      );
   fTemp[lbModel.DV_ZERO_P2_ZERO  ] = myGrid(i1,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_ZERO_P2_ZERO      );
   fTemp[lbModel.DV_P2_ZERO_ZERO  ] = myGrid(i1,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_P2_ZERO_ZERO      );

   fTemp[lbModel.DV_ZERO_ZERO_M1  ] = myGrid(i1,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M1      );
   fTemp[lbModel.DV_ZERO_ZERO_M2  ] = myGrid(i1,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M2      );
   fTemp[lbModel.DV_ZERO_M2_ZERO  ] = myGrid(i1,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_ZERO_M2_ZERO      );
   fTemp[lbModel.DV_M2_ZERO_ZERO  ] = myGrid(i1,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_M2_ZERO_ZERO      );

   fTemp[lbModel.DV_ZERO_P1_P1    ] = myGrid(i1,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_ZERO_P1_P1        );
   fTemp[lbModel.DV_ZERO_M1_P1    ] = myGrid(i1,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_ZERO_M1_P1        );
   fTemp[lbModel.DV_P1_ZERO_P1    ] = myGrid(i1,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_P1_ZERO_P1        );
   fTemp[lbModel.DV_M1_ZERO_P1    ] = myGrid(i1,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_M1_ZERO_P1        );

   fTemp[lbModel.DV_ZERO_M1_M1    ] = myGrid(i1,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_ZERO_M1_M1        );
   fTemp[lbModel.DV_ZERO_P1_M1    ] = myGrid(i1,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_ZERO_P1_M1        );
   fTemp[lbModel.DV_M1_ZERO_M1    ] = myGrid(i1,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_M1_ZERO_M1        );
   fTemp[lbModel.DV_P1_ZERO_M1    ] = myGrid(i1,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_P1_ZERO_M1        );

   fTemp[lbModel.DV_P1_P1_ZERO    ] = myGrid(i1,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_P1_P1_ZERO        );
   fTemp[lbModel.DV_M1_P1_ZERO    ] = myGrid(i1,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_M1_P1_ZERO        );
   fTemp[lbModel.DV_ZERO_P1_ZERO  ] = myGrid(i1,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_ZERO_P1_ZERO      );
   fTemp[lbModel.DV_P1_ZERO_ZERO  ] = myGrid(i1,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_P1_ZERO_ZERO      );

   fTemp[lbModel.DV_M1_M1_ZERO    ] = myGrid(i1,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_M1_M1_ZERO        );
   fTemp[lbModel.DV_P1_M1_ZERO    ] = myGrid(i1,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_P1_M1_ZERO        );
   fTemp[lbModel.DV_ZERO_M1_ZERO  ] = myGrid(i1,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_ZERO_M1_ZERO      );
   fTemp[lbModel.DV_M1_ZERO_ZERO  ] = myGrid(i1,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_M1_ZERO_ZERO      );

   fTemp[lbModel.DV_P_P_P         ] = myGrid(i1,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_P_P_P             );
   fTemp[lbModel.DV_M_P_P         ] = myGrid(i1,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_M_P_P             );
   fTemp[lbModel.DV_M_M_P         ] = myGrid(i1,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_M_M_P             );
   fTemp[lbModel.DV_P_M_P         ] = myGrid(i1,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_P_M_P             );

   fTemp[lbModel.DV_M_M_M         ] = myGrid(i1,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_M_M_M             );
   fTemp[lbModel.DV_P_M_M         ] = myGrid(i1,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_P_M_M             );
   fTemp[lbModel.DV_P_P_M         ] = myGrid(i1,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_P_P_M             );
   fTemp[lbModel.DV_M_P_M         ] = myGrid(i1,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_M_P_M             );

   fTemp[lbModel.DV_P1_P1_P1      ] = myGrid(i1,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_P1_P1_P1          );
   fTemp[lbModel.DV_M1_P1_P1      ] = myGrid(i1,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_M1_P1_P1          );
   fTemp[lbModel.DV_M1_M1_P1      ] = myGrid(i1,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_M1_M1_P1          );
   fTemp[lbModel.DV_P1_M1_P1      ] = myGrid(i1,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_P1_M1_P1          );

   fTemp[lbModel.DV_M1_M1_M1      ] = myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1         );
   fTemp[lbModel.DV_P1_M1_M1      ] = myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1         );
   fTemp[lbModel.DV_P1_P1_M1      ] = myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1         );
   fTemp[lbModel.DV_M1_P1_M1      ] = myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1         );
}



template <int N,int numblock, typename dataType1>
void printMoments(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 u_inlet,dataType1 dt, dataType1 beta,int step,dataType1 convectionTime,int size,int myRankZ, std::ofstream &file, int *coord,dataType1 F1, dataType1 F2, dataType1 F3,int centerX,int centerY, int centerZ, double f2g_factor, double length, double totLength, int countTimo, int localX, int localY, int localZ, int globalX, int globalY, int globalZ)
{
 dataType1 rho     ;
 dataType1 uX      ;
 dataType1 uY      ;
 dataType1 uZ      ;
 dataType1 Pxx     ;
 dataType1 Pyy     ;
 dataType1 Pzz     ;
 dataType1 Pxy     ;
 dataType1 Pyz     ;
 dataType1 Pzx     ;
 dataType1 qx      ;
 dataType1 qy      ;
 dataType1 qz      ;
 dataType1 zeta4   ;
 // double Pxx_temp;
 // double Pyy_temp;
 // double Pzz_temp;
 // double Pxy_temp;
 // double Pyz_temp;
 // double Pzx_temp;

 dataType1 f[41]   ;

 int nCoresZ = 2;
 int nCoresY = 10;


int y0 = myGrid.nB2;
int yF = myGrid.nE2;


int CenterXlocal = (int)(0.5*(myGrid.nE1 + myGrid.nB1));
int CenterYlocal = (int)(0.5*(myGrid.nE2 + myGrid.nB2));
int CenterZlocal = (int)(0.5*(myGrid.nE3 + myGrid.nB3));
//
// for (int localZ = myGrid.nB3; localZ <= myGrid.nE3; localZ++){
// 	for (int localY = myGrid.nB2; localY <= myGrid.nE2; localY++){
// 		for (int localX = myGrid.nB1; localX <= myGrid.nE1; localX++){

			// int globalX = localX + simParam.n1*coord[0];
			// int globalY = localY + simParam.n2*coord[1];
			// int globalZ = localZ + simParam.n3*coord[2];

			//if (globalX == CenterX && globalZ == CenterZ){
			// if (localX == CenterXlocal && localZ == CenterZlocal){

/*if(myRankZ == lowerCoupleLB)
{
    y0 = myGrid.nB2 + (0.5*nCoresZ - 1);
    yF = myGrid.nE2;
}
else if(myRankZ == upperCoupleLB)
{
    y0 = myGrid.nB2;
    yF = myGrid.nE2 - (0.5*nCoresZ - 1);
}*/

 //for(int y=y0;y<=yF;y++)
 {
  rho  = 0.0;
  uX   = 0.0;   uY   = 0.0;   uZ   = 0.0;
  Pxx  = 0.0;   Pxy  = 0.0;   qx   = 0.0;
  Pyy  = 0.0;   Pyz  = 0.0;   qy   = 0.0;
  Pzz  = 0.0;   Pzx  = 0.0;   qz   = 0.0;
  zeta4= 0.0;
 //  Pxx_temp  = 0.;
 // Pyy_temp = 0.;
 // Pzz_temp = 0.;
 // Pxy_temp = 0.;
 // Pyz_temp = 0.;
 // Pzx_temp = 0.;


  copyFromNodeToArray(lbModel,myGrid,localX,localY,localZ,f) ;

  for(int dv=0;dv<lbModel.dvN;dv++)
  {
    rho    += f[dv];
    uX     += f[dv]*lbModel.cx[dv]   ;
    uY     += f[dv]*lbModel.cy[dv]   ;
    uZ     += f[dv]*lbModel.cz[dv]   ;
    Pxx    += f[dv]*lbModel.cx2[dv]  ;
    Pyy    += f[dv]*lbModel.cy2[dv]  ;
    Pzz    += f[dv]*lbModel.cz2[dv]  ;
    Pxy    += f[dv]*lbModel.cxcy[dv] ;
    Pyz    += f[dv]*lbModel.cycz[dv] ;
    Pzx    += f[dv]*lbModel.czcx[dv] ;
    // Pxx_temp    += (f[dv] - getFEqSIMD*lbModel.cx2[dv]  ;
    // Pyy_temp    += (f[dv] - getFEqSIMD*lbModel.cy2[dv]  ;
    // Pzz_temp    += (f[dv] - getFEqSIMD*lbModel.cz2[dv]  ;
    // Pxy_temp    += (f[dv] - getFEqSIMD*lbModel.cxcy[dv] ;
    // Pyz_temp    += (f[dv] - getFEqSIMD*lbModel.cycz[dv] ;
    // Pzx_temp    += (f[dv] - getFEqSIMD*lbModel.czcx[dv] ;
    qx     += f[dv]*lbModel.cxcsq[dv];
    qy     += f[dv]*lbModel.cycsq[dv];
    qz     += f[dv]*lbModel.czcsq[dv];
    zeta4  += f[dv]*lbModel.cc[dv]*lbModel.cc[dv];
  }

  // dataType1 theta = lbModel.theta0;

  uX = uX/rho;
  uY = uY/rho;
  uZ = uZ/rho;

  // double array[41];
  // getFEqSinglePointIntoArray(lbModel ,rho, uX, uY, uZ, theta, array);
  // getFEqSIMD(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);
  // copyEqToArray(lbModel,i);

  // for(int dv=0;dv<lbModel.dvN;dv++)
  // {
  //   Pxx_temp    += (array[dv])*lbModel.cx2[dv]  ;
  //   Pyy_temp    += (array[dv])*lbModel.cy2[dv]  ;
  //   Pzz_temp    += (array[dv])*lbModel.cz2[dv]  ;
  //   Pxy_temp    += (array[dv])*lbModel.cxcy[dv] ;
  //   Pyz_temp    += (array[dv])*lbModel.cycz[dv] ;
  //   Pzx_temp    += (array[dv])*lbModel.czcx[dv] ;
  // }



  double u2 =uX*uX + uY*uY + uZ*uZ;
dataType1 theta = ((Pxx+Pyy+Pzz)-rho*(uX*uX + uY*uY + uZ*uZ))/(3.0*rho);
//   qx = qx - 0.5*rho*(u2*uX) - 1.5*rho*theta*uX ;

// uX += 0.5*dt*F1 ;
// uY += 0.5*dt*F2 ;
// uZ += 0.5*dt*F3 ;


double  sigma_Pxx = Pxx - rho*uX*uX;
double  sigma_Pyy = Pyy - rho*uY*uY;
double  sigma_Pzz = Pzz - rho*uZ*uZ;
double  sigma_Pxy = Pxy - rho*uX*uY;
double  sigma_Pyz = Pyz - rho*uY*uZ;
double  sigma_Pzx = Pzx - rho*uZ*uX;

double  qX_LB = (qx - ((rho*u2*uZ) + 3.0*rho*theta*uZ + 2.0*(uX*sigma_Pxy   + uY*sigma_Pzx + uZ*sigma_Pxx)));
double  qY_LB = (qy - ((rho*u2*uX) + 3.0*rho*theta*uX + 2.0*(uX*sigma_Pyy   + uY*sigma_Pyz + uZ*sigma_Pxy)));
double  qZ_LB = (qz - ((rho*u2*uY) + 3.0*rho*theta*uY + 2.0*(uX*sigma_Pyz   + uY*sigma_Pzz + uZ*sigma_Pzx)));



//   file<<globalY<<"  "<<rho<<"  "<<uX<<"  "<<uY<<"  "<<uZ<<"  "<<theta<<"  "<<Pxx<<"  "<<Pyy<<"  "<<Pzz<<"  "<<Pxy<<"  "<<Pyz<<"  "<<Pzx<<"  "<<qx<<"  "<<qy<<"  "<<qz<<"  "<<zeta4<<"\n";
     //  1            2         3        4          5         6             7          8          9          10        11        12         13        14        15         16       17
file<<globalX<<"\t"<<globalY<<"\t"<<globalZ<<"\t"<<rho<<"\t"<<theta<<"\t"<<rho*theta<<"\t"<<uX<<"\t"<<uY<<"\t"<<uZ<<"\t"<<sigma_Pxx*f2g_factor<<"\t"<<sigma_Pyy*f2g_factor<<"\t"<<sigma_Pzz*f2g_factor<<"\t"<<sigma_Pxy*f2g_factor<<"\t"<<sigma_Pyz*f2g_factor<<"\t"<<sigma_Pzx*f2g_factor<<std::endl;
       //  1          2         3         4          5              6               7               8                9              10                11            12           13          14                 15               16               17              18               19           20        21      22
//std::cout<<"Printed hereLB: "<<countTimo<<" "<<myRankZ<<" "<<(myRankZ*totLength/((double)nCoresZ)) + y*length*5.<<"  "<<rho<<"  "<<theta<<"  "<<rho*theta<<"  "<<uX<<"  "<<uY<<"  "<<uZ<<"  "<<sigma_Pxx*f2g_factor<<"  "<<sigma_Pyy*f2g_factor<<"  "<<sigma_Pzz*f2g_factor<<"  "<<sigma_Pxy*f2g_factor<<"  "<<sigma_Pyz*f2g_factor<<"  "<<sigma_Pzx*f2g_factor<<" "<<qX_LB<<" "<<qY_LB<<" "<<qZ_LB<<std::endl;
 }

// 			}
//
// 		}
// 	}
// }


}


template <int N,int numblock, typename dataType1>
void printMomentsLocal(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 u_inlet,dataType1 dt, dataType1 beta,int step,dataType1 convectionTime,int size,int myRankZ, std::ofstream &file, int *coord,dataType1 F1, dataType1 F2, dataType1 F3,int centerX,int centerY, int centerZ, double f2g_factor, double length, double totLength, int countTimo)
{
 dataType1 rho     ;
 dataType1 uX      ;
 dataType1 uY      ;
 dataType1 uZ      ;
 dataType1 Pxx     ;
 dataType1 Pyy     ;
 dataType1 Pzz     ;
 dataType1 Pxy     ;
 dataType1 Pyz     ;
 dataType1 Pzx     ;
 dataType1 qx      ;
 dataType1 qy      ;
 dataType1 qz      ;
 dataType1 zeta4   ;
 // double Pxx_temp;
 // double Pyy_temp;
 // double Pzz_temp;
 // double Pxy_temp;
 // double Pyz_temp;
 // double Pzx_temp;

 dataType1 f[41]   ;

 int nCoresZ = 2;
 int nCoresY = 10;


int y0 = myGrid.nB2;
int yF = myGrid.nE2;


int CenterXlocal = (int)(0.5*(myGrid.nE1 + myGrid.nB1));
int CenterYlocal = (int)(0.5*(myGrid.nE2 + myGrid.nB2));
int CenterZlocal = (int)(0.5*(myGrid.nE3 + myGrid.nB3));

for (int localZ = myGrid.nB3; localZ <= myGrid.nE3; localZ++){
	for (int localY = myGrid.nB2; localY <= myGrid.nE2; localY++){
		for (int localX = myGrid.nB1; localX <= myGrid.nE1; localX++){

			//int globalX = localX + simParam.n1*coord[0];
			//int globalY = localY + simParam.n2*coord[1];
			//int globalZ = localZ + simParam.n3*coord[2];

			//if (globalX == CenterX && globalZ == CenterZ){
			if (localX == CenterXlocal && localZ == CenterZlocal){

/*if(myRankZ == lowerCoupleLB)
{
    y0 = myGrid.nB2 + (0.5*nCoresZ - 1);
    yF = myGrid.nE2;
}
else if(myRankZ == upperCoupleLB)
{
    y0 = myGrid.nB2;
    yF = myGrid.nE2 - (0.5*nCoresZ - 1);
}*/

 //for(int y=y0;y<=yF;y++)
 {
  rho  = 0.0;
  uX   = 0.0;   uY   = 0.0;   uZ   = 0.0;
  Pxx  = 0.0;   Pxy  = 0.0;   qx   = 0.0;
  Pyy  = 0.0;   Pyz  = 0.0;   qy   = 0.0;
  Pzz  = 0.0;   Pzx  = 0.0;   qz   = 0.0;
  zeta4= 0.0;
 //  Pxx_temp  = 0.;
 // Pyy_temp = 0.;
 // Pzz_temp = 0.;
 // Pxy_temp = 0.;
 // Pyz_temp = 0.;
 // Pzx_temp = 0.;


  copyFromNodeToArray(lbModel,myGrid,localX,localY,localZ,f) ;

  for(int dv=0;dv<lbModel.dvN;dv++)
  {
    rho    += f[dv];
    uX     += f[dv]*lbModel.cx[dv]   ;
    uY     += f[dv]*lbModel.cy[dv]   ;
    uZ     += f[dv]*lbModel.cz[dv]   ;
    Pxx    += f[dv]*lbModel.cx2[dv]  ;
    Pyy    += f[dv]*lbModel.cy2[dv]  ;
    Pzz    += f[dv]*lbModel.cz2[dv]  ;
    Pxy    += f[dv]*lbModel.cxcy[dv] ;
    Pyz    += f[dv]*lbModel.cycz[dv] ;
    Pzx    += f[dv]*lbModel.czcx[dv] ;
    // Pxx_temp    += (f[dv] - getFEqSIMD*lbModel.cx2[dv]  ;
    // Pyy_temp    += (f[dv] - getFEqSIMD*lbModel.cy2[dv]  ;
    // Pzz_temp    += (f[dv] - getFEqSIMD*lbModel.cz2[dv]  ;
    // Pxy_temp    += (f[dv] - getFEqSIMD*lbModel.cxcy[dv] ;
    // Pyz_temp    += (f[dv] - getFEqSIMD*lbModel.cycz[dv] ;
    // Pzx_temp    += (f[dv] - getFEqSIMD*lbModel.czcx[dv] ;
    qx     += f[dv]*lbModel.cxcsq[dv];
    qy     += f[dv]*lbModel.cycsq[dv];
    qz     += f[dv]*lbModel.czcsq[dv];
    zeta4  += f[dv]*lbModel.cc[dv]*lbModel.cc[dv];
  }

  // dataType1 theta = lbModel.theta0;

  uX = uX/rho;
  uY = uY/rho;
  uZ = uZ/rho;

  // double array[41];
  // getFEqSinglePointIntoArray(lbModel ,rho, uX, uY, uZ, theta, array);
  // getFEqSIMD(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);
  // copyEqToArray(lbModel,i);

  // for(int dv=0;dv<lbModel.dvN;dv++)
  // {
  //   Pxx_temp    += (array[dv])*lbModel.cx2[dv]  ;
  //   Pyy_temp    += (array[dv])*lbModel.cy2[dv]  ;
  //   Pzz_temp    += (array[dv])*lbModel.cz2[dv]  ;
  //   Pxy_temp    += (array[dv])*lbModel.cxcy[dv] ;
  //   Pyz_temp    += (array[dv])*lbModel.cycz[dv] ;
  //   Pzx_temp    += (array[dv])*lbModel.czcx[dv] ;
  // }



  double u2 =uX*uX + uY*uY + uZ*uZ;
dataType1 theta = ((Pxx+Pyy+Pzz)-rho*(uX*uX + uY*uY + uZ*uZ))/(3.0*rho);
//   qx = qx - 0.5*rho*(u2*uX) - 1.5*rho*theta*uX ;

// uX += 0.5*dt*F1 ;
// uY += 0.5*dt*F2 ;
// uZ += 0.5*dt*F3 ;


double  sigma_Pxx = Pxx - rho*uX*uX;
double  sigma_Pyy = Pyy - rho*uY*uY;
double  sigma_Pzz = Pzz - rho*uZ*uZ;
double  sigma_Pxy = Pxy - rho*uX*uY;
double  sigma_Pyz = Pyz - rho*uY*uZ;
double  sigma_Pzx = Pzx - rho*uZ*uX;

double  qX_LB = (qx - ((rho*u2*uZ) + 3.0*rho*theta*uZ + 2.0*(uX*sigma_Pxy   + uY*sigma_Pzx + uZ*sigma_Pxx)));
double  qY_LB = (qy - ((rho*u2*uX) + 3.0*rho*theta*uX + 2.0*(uX*sigma_Pyy   + uY*sigma_Pyz + uZ*sigma_Pxy)));
double  qZ_LB = (qz - ((rho*u2*uY) + 3.0*rho*theta*uY + 2.0*(uX*sigma_Pyz   + uY*sigma_Pzz + uZ*sigma_Pzx)));



double globalY = (double)localY + double(myGrid.m2*coord[1]) ;
if(localY == myGrid.nB2)
{
    std::cout<<"localY: "<<myRankZ<<"\t"<<totLength<<"\t"<<(double)localY<<"\t"<<localY<<"\t"<<globalY<<"\t"<<(myRankZ*totLength/((double)nCoresZ)) + (localY-4.)*length<<"\t"<<(myRankZ*totLength/((double)nCoresZ))<<"\t"<<(localY-4.)*length<<std::endl;
}
//   file<<globalY<<"  "<<rho<<"  "<<uX<<"  "<<uY<<"  "<<uZ<<"  "<<theta<<"  "<<Pxx<<"  "<<Pyy<<"  "<<Pzz<<"  "<<Pxy<<"  "<<Pyz<<"  "<<Pzx<<"  "<<qx<<"  "<<qy<<"  "<<qz<<"  "<<zeta4<<"\n";
     //  1            2         3        4          5         6             7          8          9          10        11        12         13        14        15         16       17
file<<length*(globalY - 3.5)<<"\t"<<rho<<"\t"<<theta<<"\t"<<rho*theta<<"\t"<<uX<<"\t"<<uY<<"\t"<<uZ<<"\t"<<sigma_Pxx*f2g_factor<<"\t"<<sigma_Pyy*f2g_factor<<"\t"<<sigma_Pzz*f2g_factor<<"\t"<<sigma_Pxy*f2g_factor<<"\t"<<sigma_Pyz*f2g_factor<<"\t"<<sigma_Pzx*f2g_factor<<std::endl<<std::endl<<std::endl<<std::endl;
       //  1          2         3         4          5              6               7               8                9              10                11            12           13          14                 15               16               17              18               19           20        21      22
//std::cout<<"Printed hereLB: "<<countTimo<<" "<<myRankZ<<" "<<(myRankZ*totLength/((double)nCoresZ)) + y*length*5.<<"  "<<rho<<"  "<<theta<<"  "<<rho*theta<<"  "<<uX<<"  "<<uY<<"  "<<uZ<<"  "<<sigma_Pxx*f2g_factor<<"  "<<sigma_Pyy*f2g_factor<<"  "<<sigma_Pzz*f2g_factor<<"  "<<sigma_Pxy*f2g_factor<<"  "<<sigma_Pyz*f2g_factor<<"  "<<sigma_Pzx*f2g_factor<<" "<<qX_LB<<" "<<qY_LB<<" "<<qZ_LB<<std::endl;
 }

			}

		}
	}
}


}

template <int N,int numblock, typename dataType1>
void printMomentsCell(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 u_inlet,dataType1 dt, dataType1 beta,int step,dataType1 convectionTime,int size,int myRank, int *coord,dataType1 F1, dataType1 F2, dataType1 F3,int centerX,int centerY, int centerZ)
{
 dataType1 rho     ;
 dataType1 uX      ;
 dataType1 uY      ;
 dataType1 uZ      ;
 dataType1 Pxx     ;
 dataType1 Pyy     ;
 dataType1 Pzz     ;
 dataType1 Pxy     ;
 dataType1 Pyz     ;
 dataType1 Pzx     ;
 dataType1 qx      ;
 dataType1 qy      ;
 dataType1 qz      ;
 dataType1 zeta4   ;
 // double Pxx_temp;
 // double Pyy_temp;
 // double Pzz_temp;
 // double Pxy_temp;
 // double Pyz_temp;
 // double Pzx_temp;

 dataType1 f[41]   ;






 // for(int y=myGrid.nB2;y<=myGrid.nE2;y++)
 {
  rho  = 0.0;
  uX   = 0.0;   uY   = 0.0;   uZ   = 0.0;
  Pxx  = 0.0;   Pxy  = 0.0;   qx   = 0.0;
  Pyy  = 0.0;   Pyz  = 0.0;   qy   = 0.0;
  Pzz  = 0.0;   Pzx  = 0.0;   qz   = 0.0;
  zeta4= 0.0;
 //  Pxx_temp  = 0.;
 // Pyy_temp = 0.;
 // Pzz_temp = 0.;
 // Pxy_temp = 0.;
 // Pyz_temp = 0.;
 // Pzx_temp = 0.;


  copyFromNodeToArray(lbModel,myGrid,centerX,centerY,centerZ,f) ;

  for(int dv=0;dv<lbModel.dvN;dv++)
  {
    rho    += f[dv];
    uX     += f[dv]*lbModel.cx[dv]   ;
    uY     += f[dv]*lbModel.cy[dv]   ;
    uZ     += f[dv]*lbModel.cz[dv]   ;
    Pxx    += f[dv]*lbModel.cx2[dv]  ;
    Pyy    += f[dv]*lbModel.cy2[dv]  ;
    Pzz    += f[dv]*lbModel.cz2[dv]  ;
    Pxy    += f[dv]*lbModel.cxcy[dv] ;
    Pyz    += f[dv]*lbModel.cycz[dv] ;
    Pzx    += f[dv]*lbModel.czcx[dv] ;
    // Pxx_temp    += (f[dv] - getFEqSIMD*lbModel.cx2[dv]  ;
    // Pyy_temp    += (f[dv] - getFEqSIMD*lbModel.cy2[dv]  ;
    // Pzz_temp    += (f[dv] - getFEqSIMD*lbModel.cz2[dv]  ;
    // Pxy_temp    += (f[dv] - getFEqSIMD*lbModel.cxcy[dv] ;
    // Pyz_temp    += (f[dv] - getFEqSIMD*lbModel.cycz[dv] ;
    // Pzx_temp    += (f[dv] - getFEqSIMD*lbModel.czcx[dv] ;
    qx     += f[dv]*lbModel.cxcsq[dv];
    qy     += f[dv]*lbModel.cycsq[dv];
    qz     += f[dv]*lbModel.czcsq[dv];
    zeta4  += f[dv]*lbModel.cc[dv]*lbModel.cc[dv];
  }

  dataType1 theta = lbModel.theta0;

  uX = uX/rho;
  uY = uY/rho;
  uZ = uZ/rho;

  // double array[41];
  // getFEqSinglePointIntoArray(lbModel ,rho, uX, uY, uZ, theta, array);
  // getFEqSIMD(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);
  // copyEqToArray(lbModel,i);

  // for(int dv=0;dv<lbModel.dvN;dv++)
  // {
  //   Pxx_temp    += (array[dv])*lbModel.cx2[dv]  ;
  //   Pyy_temp    += (array[dv])*lbModel.cy2[dv]  ;
  //   Pzz_temp    += (array[dv])*lbModel.cz2[dv]  ;
  //   Pxy_temp    += (array[dv])*lbModel.cxcy[dv] ;
  //   Pyz_temp    += (array[dv])*lbModel.cycz[dv] ;
  //   Pzx_temp    += (array[dv])*lbModel.czcx[dv] ;
  // }



  double u2 =uX*uX + uY*uY + uZ*uZ;
//   qx = qx - 0.5*rho*(u2*uX) - 1.5*rho*theta*uX ;

  /*uX += 0.5*dt*F1 ;
  uY += 0.5*dt*F2 ;
  uZ += 0.5*dt*F3 ;*/


  double  sigma_Pxx = Pxx - rho*uX*uX;
  double  sigma_Pyy = Pyy - rho*uY*uY;
  double  sigma_Pzz = Pzz - rho*uZ*uZ;
  double  sigma_Pxy = Pxy - rho*uX*uY;
  double  sigma_Pyz = Pyz - rho*uY*uZ;
  double  sigma_Pzx = Pzx - rho*uZ*uX;

  double  qX_LB = (qx - ((rho*u2*uZ) + 3.0*rho*theta*uZ + 2.0*(uX*sigma_Pxy   + uY*sigma_Pzx + uZ*sigma_Pxx)));
  double  qY_LB = (qy - ((rho*u2*uX) + 3.0*rho*theta*uX + 2.0*(uX*sigma_Pyy   + uY*sigma_Pyz + uZ*sigma_Pxy)));
  double  qZ_LB = (qz - ((rho*u2*uY) + 3.0*rho*theta*uY + 2.0*(uX*sigma_Pyz   + uY*sigma_Pzz + uZ*sigma_Pzx)));



  // double globalY = (double)y + double(myGrid.m2*coord[1]) ;
//   file<<globalY<<"  "<<rho<<"  "<<uX<<"  "<<uY<<"  "<<uZ<<"  "<<theta<<"  "<<Pxx<<"  "<<Pyy<<"  "<<Pzz<<"  "<<Pxy<<"  "<<Pyz<<"  "<<Pzx<<"  "<<qx<<"  "<<qy<<"  "<<qz<<"  "<<zeta4<<"\n";
       //  1            2         3        4          5         6             7          8          9          10        11        12         13        14        15         16       17
  std::cout<<"PRINT MOMENTS CELL200: "<<step<<"  "<<myRank<<"  "<<rho<<"  "<<uX<<"  "<<uY<<"  "<<uZ<<"  "<<sigma_Pxx<<"  "<<sigma_Pyy<<"  "<<sigma_Pzz<<"  "<<sigma_Pxy<<"  "<<sigma_Pyz<<"  "<<sigma_Pzx<<" "<<qX_LB<<" "<<qY_LB<<" "<<qZ_LB<<"\n";
       //  1          2         3         4          5              6               7               8                9              10                11            12           13          14                 15               16               17              18               19           20        21      22
 }



}


template void printMoments<4,11,double>(lbmRD3Q41<double> &lbModel, gridBCC3D<4,11,double> &myGrid,int VECT_LENGTH,double u_inlet,double dt, double beta,int step,double convectionTime,int size,int myRank, std::ofstream &file, int *coord,double F1, double F2, double F3,int centerX,int centerY, int centerZ, double f2g_factor, double length,double totLength,int countTimo, int , int, int, int , int, int);
template void printMomentsCell<4,11,double>(lbmRD3Q41<double> &lbModel, gridBCC3D<4,11,double> &myGrid,int VECT_LENGTH,double u_inlet,double dt, double beta,int step,double convectionTime,int size,int myRank, int *coord,double F1, double F2, double F3,int centerX,int centerY, int centerZ);
template void printMomentsLocal<4,11,double>(lbmRD3Q41<double> &lbModel, gridBCC3D<4, 11, double> &myGrid,int VECT_LENGTH,double u_inlet,double dt, double beta,int step,double convectionTime,int size,int myRankZ, std::ofstream &file, int *coord,double F1, double F2, double F3,int centerX,int centerY, int centerZ, double f2g_factor, double length, double totLength, int countTimo);
