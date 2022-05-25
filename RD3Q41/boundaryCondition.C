#include"boundaryCondition.h"


    template<typename dataType1>
    void getFEqSinglePoint(lbmRD3Q41<dataType1> &lbModel , dataType1 rho, dataType1 uX, dataType1 uY, dataType1 uZ, dataType1 theta,int point)
    {
        dataType1 delThetaBy2;
        dataType1 delSqTheta0p125;
        dataType1 twoMinusU2,order1,order2,f0V,f0V1,uDotC,uDotCSq;
        dataType1 oneByTheta;

            delThetaBy2 = 0.5*(theta - lbModel.theta0)*lbModel.oneByTheta0;
            delSqTheta0p125 =  0.5* delThetaBy2 *delThetaBy2;

            oneByTheta  = 1.0/(theta);
            uX         *= oneByTheta;
            uY         *= oneByTheta;
            uZ         *= oneByTheta;
            rho        *= 0.5;


            twoMinusU2 = 2.0-(uX*uX + uY*uY + uZ*uZ)*theta;


            // CENTER
            order1       =  delThetaBy2      *lbModel.yMinus3_CENTER;
            order2       =  delSqTheta0p125  *lbModel.yTenMinusYSqrPlus15_CENTER;
            f0V          =  rho*lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO]*(1.0 + order1+ order2);
            lbModel.fTemp0[point][lbModel.CENTER_DV_ZERO_ZERO_ZERO]=  f0V*twoMinusU2;

            // G1 and G2
            // ZZP1 and ZZM1
            order1       =  delThetaBy2     *lbModel.yMinus3_SC1;
            order2       =  delSqTheta0p125 *lbModel.yTenMinusYSqrPlus15_SC1;
            f0V          =  rho*lbModel.wt[lbModel.DV_ZERO_ZERO_P1]*(1.0 + order1+ order2);
            uDotC        =  uZ;
            lbModel.fTemp1[point][lbModel.G1_DV_ZERO_ZERO_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            lbModel.fTemp2[point][lbModel.G2_DV_ZERO_ZERO_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));

            order1       =  delThetaBy2     *lbModel.yMinus3_SC2;
            order2       =  delSqTheta0p125 *lbModel.yTenMinusYSqrPlus15_SC2;
            f0V          =  rho*lbModel.wt[lbModel.DV_ZERO_ZERO_P2]*(1.0 + order1+ order2);
            // ZZP2 and ZZM2
            uDotC        =  2.0*uZ;
            lbModel.fTemp1[point][lbModel.G1_DV_ZERO_ZERO_P2] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            lbModel.fTemp2[point][lbModel.G2_DV_ZERO_ZERO_M2] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
            // ZP2Z and ZM2Z
            uDotC        =  2.0*uY;
            lbModel.fTemp1[point][lbModel.G1_DV_ZERO_P2_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            lbModel.fTemp2[point][lbModel.G2_DV_ZERO_M2_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
            // P2ZZ and M2ZZ
            uDotC        =  2.0*uX;
            lbModel.fTemp1[point][lbModel.G1_DV_P2_ZERO_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            lbModel.fTemp2[point][lbModel.G2_DV_M2_ZERO_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));



            // G3 and G4
            order1    =  delThetaBy2    *lbModel.yMinus3_FCC;
            order2    =  delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_FCC;
            f0V       =  rho*lbModel.wt[lbModel.DV_ZERO_P1_P1]*(1.0 + order1+ order2);
            //ZPP and ZMM
            uDotC     =  (uY + uZ);
            lbModel.fTemp3[point][lbModel.G3_DV_ZERO_P1_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            lbModel.fTemp4[point][lbModel.G4_DV_ZERO_M1_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
            //ZMP and ZPM
            uDotC     =  (uZ - uY);
            lbModel.fTemp3[point][lbModel.G3_DV_ZERO_M1_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            lbModel.fTemp4[point][lbModel.G4_DV_ZERO_P1_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
            //PZP and MZM
            uDotC     =  (uX + uZ);
            lbModel.fTemp3[point][lbModel.G3_DV_P1_ZERO_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            lbModel.fTemp4[point][lbModel.G4_DV_M1_ZERO_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
            //MZP and PZM
            uDotC     =  (uZ - uX);
            lbModel.fTemp3[point][lbModel.G3_DV_M1_ZERO_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            lbModel.fTemp4[point][lbModel.G4_DV_P1_ZERO_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));


            order1    =  delThetaBy2    *lbModel.yMinus3_FCC;
            order2    =  delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_FCC;
            f0V       =  rho*lbModel.wt[lbModel.DV_P1_P1_ZERO]*(1.0 + order1+ order2);
            //PPZ and MMZ
            uDotC     =  uX + uY;
            lbModel.fTemp5[point][lbModel.G5_DV_P1_P1_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            lbModel.fTemp6[point][lbModel.G6_DV_M1_M1_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
            //MPZ and PMZ
            uDotC     =  (uY - uX);
            lbModel.fTemp5[point][lbModel.G5_DV_M1_P1_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            lbModel.fTemp6[point][lbModel.G6_DV_P1_M1_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
            //ZPZ and ZMZ
            order1    =  delThetaBy2    *lbModel.yMinus3_SC1;
            order2    =  delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_SC1;
            f0V       =  rho*lbModel.wt[lbModel.DV_ZERO_P1_ZERO]*(1.0 + order1+ order2);
            uDotC     =  uY;
            lbModel.fTemp5[point][lbModel.G5_DV_ZERO_P1_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            lbModel.fTemp6[point][lbModel.G6_DV_ZERO_M1_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
            //PZZ aqnd MZZ
            uDotC     =   uX;
            lbModel.fTemp5[point][lbModel.G5_DV_P1_ZERO_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            lbModel.fTemp6[point][lbModel.G6_DV_M1_ZERO_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));



            order1    =  1.0 + delThetaBy2*lbModel.yMinus3_BCC +  delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_BCC;
            f0V       =  rho*lbModel.wt[lbModel.DV_P_P_P]*(order1);
            order2    =  1.0 + delThetaBy2*lbModel.yMinus3_BCC1 + delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_BCC1;
            f0V1      =  rho*lbModel.wt[lbModel.DV_P1_P1_P1]*(order2);
            //PPP1 and MMM1
            uDotC     = (uX + uY + uZ);
            uDotCSq   = uDotC*uDotC;
            lbModel.fTemp9 [point][lbModel.G9_DV_P1_P1_P1]   = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);
            lbModel.fTemp10[point][lbModel.G10_DV_M1_M1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);
            //PPP and MMM
            //uDotC     =   uDotC*0.5 ;
            lbModel.fTemp7[point][lbModel.G7_DV_P_P_P]   = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);
            lbModel.fTemp8[point][lbModel.G8_DV_M_M_M]   = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);

            //MPP1 and PMM1
            uDotC     = (uY - uX + uZ);
            uDotCSq   = uDotC*uDotC;
            lbModel.fTemp9 [point][lbModel.G9_DV_M1_P1_P1]  = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);
            lbModel.fTemp10[point][lbModel.G10_DV_P1_M1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);
            //MPP and PMM
            // uDotC     =  uDotC*0.5 ;
            lbModel.fTemp7[point][lbModel.G7_DV_M_P_P] = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);
            lbModel.fTemp8[point][lbModel.G8_DV_P_M_M] = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);

            //MMP1 and PPM 1
            uDotC     = (uZ - uX - uY);
            uDotCSq   = uDotC*uDotC;
            lbModel.fTemp9 [point][lbModel.G9_DV_M1_M1_P1] = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);
            lbModel.fTemp10[point][lbModel.G10_DV_P1_P1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);
            //MMP and PPM
            // uDotC     =   uDotC*0.5 ;
            lbModel.fTemp7[point][lbModel.G7_DV_M_M_P] = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);
            lbModel.fTemp8[point][lbModel.G8_DV_P_P_M] = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);

            //PMP and MPM
            uDotC     = (uX - uY + uZ);
            uDotCSq   = uDotC*uDotC;
            lbModel.fTemp9 [point][lbModel.G9_DV_P1_M1_P1]  = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);
            lbModel.fTemp10[point][lbModel.G10_DV_M1_P1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);
            //PMP and MPM
            // uDotC     =  uDotC *0.5 ;
            lbModel.fTemp7[point][lbModel.G7_DV_P_M_P] = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);
            lbModel.fTemp8[point][lbModel.G8_DV_M_P_M] = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);


    }



 template<typename dataType1>
 void getGradThermalSinglePoint(lbmRD3Q41<dataType1> &lbModel , dataType1 rho, dataType1 uX, dataType1 uY, dataType1 uZ, dataType1 theta, dataType1 pXX, dataType1 pYY, dataType1 pZZ, dataType1 pXY, dataType1 pYZ, dataType1 pZX ,dataType1 *grad)
{
//  dataType1 cS2 = theta;
//  dataType1 oneBycS2 = 1.0/cS2;
 dataType1 delTheta = (theta - lbModel.theta0)*lbModel.oneByTheta0;
 dataType1 dot(0.0),c2(0.0),cX2(0.0),cY2(0.0),cZ2(0.0),y[lbModel.dvN],f0V[lbModel.dvN],order1(0.0),order2(0.0);
 dataType1 A(0.0),B(0.0),C(0.0),D(0.0);
 dataType1 term1(0.0),/*term2(0.0),*/term3(0.0),term4(0.0),coeff1(0.0),coeff2(0.0),coeff3(0.0),coeff4(0.0);
 A = (rho*lbModel.theta0*lbModel.theta0)*(1.0+2.0*delTheta+0.47954*delTheta*delTheta);
 A = 1.0/A;
 A=0;

 for(int dv=0;dv<lbModel.dvN;dv++)
 {
  y[dv]     =  lbModel.oneByTheta0*lbModel.cc[dv];
  order1    =  0.5*delTheta*(y[dv]-3.0);
  order2    =  0.125*delTheta*delTheta*(15.0 - 10.0*y[dv] + y[dv]*y[dv]);
  f0V[dv]   =  rho*lbModel.wt[dv]*(1.0 + order1 + order2);
  cX2       =  lbModel.cx2[dv];
  cY2       =  lbModel.cy2[dv];
  cZ2       =  lbModel.cz2[dv];
//   c2        =  cX2 + cY2 + cZ2;
  A        +=  f0V[dv]*cX2*cY2;
  B        +=  f0V[dv]*(cX2 - theta)*(cX2 - theta);
  C        +=  f0V[dv]*(cY2 - cZ2)*(cY2 - cZ2);
  D        +=  f0V[dv]*(lbModel.cc[dv] - 3.0*theta)*(lbModel.cc[dv] - 3.0*theta);
  coeff1   +=  f0V[dv]*(cX2-cY2)*(cY2-cZ2);
 }

 A = 1.0/A;
 B = 1.0/B;
 C = 1.0/C;
 D = 1.0/D;
 coeff1 = coeff1*C;

 for(int dv=0;dv<lbModel.dvN;dv++)
 {
  y[dv]     =  lbModel.oneByTheta0*lbModel.cc[dv];
  order1    =  0.5*delTheta*(y[dv]-3.0);
  order2    =  0.125*delTheta*delTheta*(15.0 - 10.0*y[dv] + y[dv]*y[dv]);
  f0V[dv]   =  rho*lbModel.wt[dv]*(1.0 + order1 + order2);
  cX2       =  lbModel.cx2[dv];
  cY2       =  lbModel.cy2[dv];
  cZ2       =  lbModel.cz2[dv];
  c2        =  cX2 + cY2 + cZ2;
  coeff2    =  (cX2-cY2) - ((coeff1)*(cY2-cZ2));
  coeff4   +=  f0V[dv]*coeff2*coeff2;
 }

 coeff4 = 1.0/coeff4;
 for(int dv=0;dv<lbModel.dvN;dv++)
 {
  dot      =  lbModel.cx[dv]*uX + lbModel.cy[dv]*uY + lbModel.cz[dv]*uZ;
  cX2      =  lbModel.cx[dv]*lbModel.cx[dv];
  cY2      =  lbModel.cy[dv]*lbModel.cy[dv];
  cZ2      =  lbModel.cz[dv]*lbModel.cz[dv];
  c2       =  cX2 + cY2 + cZ2;
  coeff3   =  ((cX2-cY2) - coeff1*(cY2-cZ2))*((pXX-pYY) - coeff1*(pYY-pZZ));
  term1    =  A*(pXY*lbModel.cx[dv]*lbModel.cy[dv] + pYZ*lbModel.cy[dv]*lbModel.cz[dv] + pZX*lbModel.cz[dv]*lbModel.cx[dv]);
//   term2    =  B*(cX2-theta)*(pXX-rho*theta);
  term3    =  C*(cY2 - cZ2)*(pYY - pZZ);
  term4    =  D*(c2 - 3.0*theta)*((pXX+pYY+pZZ)-3.0*rho*theta);
  grad[dv] =  f0V[dv]*(1.0 +  dot/theta+ term1 + term3 + term4 +coeff4*coeff3);
 }



}

 template<typename dataType1>
 void getGradThermalSinglePoint_heatAlpha(lbmRD3Q41<dataType1> &lbModel , dataType1 rho, dataType1 uX, dataType1 uY, dataType1 uZ, dataType1 theta, dataType1 pXX, dataType1 pYY, dataType1 pZZ, dataType1 pXY, dataType1 pYZ, dataType1 pZX, dataType1 q1, dataType1 q2, dataType1 q3,dataType1 *grad)
{
//  dataType1 cS2 = theta;
//  dataType1 oneBycS2 = 1.0/cS2;
 dataType1 delTheta = (theta - lbModel.theta0)*lbModel.oneByTheta0;
//  std::cout<<delTheta<<std::endl;
 dataType1 dot(0.0),c2(0.0),cX2(0.0),cY2(0.0),cZ2(0.0),y[lbModel.dvN],f0V[lbModel.dvN],order1(0.0),order2(0.0);
 dataType1 A(0.0),B(0.0),C(0.0),D(0.0);
 dataType1 term1(0.0),/*term2(0.0),*/term3(0.0),term4(0.0),coeff1(0.0),coeff2(0.0),coeff3(0.0),coeff4(0.0);
 A = (rho*lbModel.theta0*lbModel.theta0)*(1.0+2.0*delTheta+0.47954*delTheta*delTheta);
 A = 1.0/A;
 A=0;


//  dataType1 q_1 = 0.0;
//  dataType1 q_2 = 0.0;
//  dataType1 q_3 = 0.0;

 for(int dv=0;dv<lbModel.dvN;dv++)
 {
  y[dv]     =  lbModel.oneByTheta0*lbModel.cc[dv];
  order1    =  0.5*delTheta*(y[dv]-3.0);
  order2    =  0.125*delTheta*delTheta*(15.0 - 10.0*y[dv] + y[dv]*y[dv]);
  f0V[dv]   =  rho*lbModel.wt[dv]*(1.0 + order1 + order2);
  cX2       =  lbModel.cx2[dv];
  cY2       =  lbModel.cy2[dv];
  cZ2       =  lbModel.cz2[dv];
//   c2        =  cX2 + cY2 + cZ2;
  A        +=  f0V[dv]*cX2*cY2;
  B        +=  f0V[dv]*(cX2 - theta)*(cX2 - theta);
  C        +=  f0V[dv]*(cY2 - cZ2)*(cY2 - cZ2);
  D        +=  f0V[dv]*(lbModel.cc[dv] - 3.0*theta)*(lbModel.cc[dv] - 3.0*theta);
  coeff1   +=  f0V[dv]*(cX2-cY2)*(cY2-cZ2);
 }

 A = 1.0/A;
 B = 1.0/B;
 C = 1.0/C;
 D = 1.0/D;
 coeff1 = coeff1*C;

 for(int dv=0;dv<lbModel.dvN;dv++)
 {
  y[dv]     =  lbModel.oneByTheta0*lbModel.cc[dv];
  order1    =  0.5*delTheta*(y[dv]-3.0);
  order2    =  0.125*delTheta*delTheta*(15.0 - 10.0*y[dv] + y[dv]*y[dv]);
  f0V[dv]   =  rho*lbModel.wt[dv]*(1.0 + order1 + order2);
  cX2       =  lbModel.cx2[dv];
  cY2       =  lbModel.cy2[dv];
  cZ2       =  lbModel.cz2[dv];
  c2        =  cX2 + cY2 + cZ2;
  coeff2    =  (cX2-cY2) - ((coeff1)*(cY2-cZ2));
  coeff4   +=  f0V[dv]*coeff2*coeff2;
 }

 coeff4 = 1.0/coeff4;

 dataType1 oneByTenTheta0Cube = 1.0/(10.0*lbModel.theta0*lbModel.theta0*lbModel.theta0);

 for(int dv=0;dv<lbModel.dvN;dv++)
 {
  dot      =  lbModel.cx[dv]*uX + lbModel.cy[dv]*uY + lbModel.cz[dv]*uZ;
  cX2      =  lbModel.cx[dv]*lbModel.cx[dv];
  cY2      =  lbModel.cy[dv]*lbModel.cy[dv];
  cZ2      =  lbModel.cz[dv]*lbModel.cz[dv];
  c2       =  cX2 + cY2 + cZ2;
  coeff3   =  ((cX2-cY2) - coeff1*(cY2-cZ2))*((pXX-pYY) - coeff1*(pYY-pZZ));
  term1    =  A*(pXY*lbModel.cx[dv]*lbModel.cy[dv] + pYZ*lbModel.cy[dv]*lbModel.cz[dv] + pZX*lbModel.cz[dv]*lbModel.cx[dv]);
//   term2    =  B*(cX2-theta)*(pXX-rho*theta);
  term3    =  C*(cY2 - cZ2)*(pYY - pZZ);
  term4    =  D*(c2 - 3.0*theta)*((pXX+pYY+pZZ)-3.0*rho*theta);

  dataType1 term5 = lbModel.wt[dv]*(q1- 5.0*rho*uX*theta) * lbModel.cx[dv] * (lbModel.cc[dv] - 5.0*lbModel.theta0) * oneByTenTheta0Cube;
  dataType1 term6 = lbModel.wt[dv]*(q2- 5.0*rho*uY*theta) * lbModel.cy[dv] * (lbModel.cc[dv] - 5.0*lbModel.theta0) * oneByTenTheta0Cube;
  dataType1 term7 = lbModel.wt[dv]*(q3- 5.0*rho*uZ*theta) * lbModel.cz[dv] * (lbModel.cc[dv] - 5.0*lbModel.theta0) * oneByTenTheta0Cube;

  grad[dv] =  f0V[dv] * ( 1.0 +  dot/theta+ term1 + term3 + term4 + coeff4*coeff3 )  + term5 + term6 + term7 ;
 }

}

    template <int N,int numblock, typename dataType1>
    void copyFromArrayToNode(lbmRD3Q41<dataType1> &lbModel ,dataType1 *grad, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3)
    {
       myGrid(i1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) = grad[lbModel.CENTER_DV_ZERO_ZERO_ZERO];


        int index = myGrid.getIndex(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv)  =  grad[1+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) =  grad[5+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) =  grad[9+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) =  grad[13+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv)=  grad[17+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv)=  grad[21+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) =  grad[25+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) =  grad[29+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],  0);
        for(int dv=0;dv<4;dv++)
        {
          myGrid.value(index+dv) =  grad[33+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) =  grad[37+dv];
        }
      }


    template <int N,int numblock, typename dataType1>
    inline void copyFromArrayToCell(lbmRD3Q41<dataType1> &lbModel ,dataType1 *grad, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3)
    {
       myGrid(i1,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) = grad[lbModel.CENTER_DV_ZERO_ZERO_ZERO];

        int index = myGrid.getIndex(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv)  =  grad[1+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) =  grad[5+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) =  grad[9+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) =  grad[13+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv)=  grad[17+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv)=  grad[21+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) =  grad[25+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) =  grad[29+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],  0);
        for(int dv=0;dv<4;dv++)
        {
          myGrid.value(index+dv) =  grad[33+dv];
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) =  grad[37+dv];
    }}

    template <int N,int numblock, typename dataType1>
    void copyFromArrayToNode2(lbmRD3Q41<dataType1> &lbModel ,dataType1 *grad, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3,int myRank)
    {
       myGrid(i1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) += grad[lbModel.CENTER_DV_ZERO_ZERO_ZERO];
       // if((myRank == 2) && (grad[lbModel.CENTER_DV_ZERO_ZERO_ZERO] > 10.))
       // {
       //    std::cout<<"grad value Node: "<<0<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[lbModel.CENTER_DV_ZERO_ZERO_ZERO]<<std::endl;
       // }

        int index = myGrid.getIndex(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv)  +=  grad[1+dv];
           // if((myRank == 2) && (grad[dv+1] > 10.))
           // {
           //    std::cout<<"grad value Node: "<<dv+1<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[1+dv]<<std::endl;
           // }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) +=  grad[5+dv];
           // if((myRank == 2) && (grad[dv+5] > 10.))
           // {
           //    std::cout<<"grad value Node: "<<dv+5<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[5+dv]<<std::endl;
           // }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) +=  grad[9+dv];
           // if((myRank == 2) && (grad[dv+9] > 10.))
           // {
           //    std::cout<<"grad value Node: "<<dv+9<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[9+dv]<<std::endl;
           // }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) +=  grad[13+dv];
            // if((myRank == 2) && (grad[dv+13] > 10.))
            // {
            //     std::cout<<"grad value Node: "<<dv+13<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[13+dv]<<std::endl;
            //  }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) +=  grad[17+dv];
           // if((myRank == 2) && (grad[dv+17] > 10.))
           // {
           //    std::cout<<"grad value Node: "<<dv+17<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[17+dv]<<std::endl;
           // }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) +=  grad[21+dv];
            // if((myRank == 2) && (grad[dv+21] > 10.))
            // {
            //     std::cout<<"grad value Node: "<<dv+21<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[21+dv]<<std::endl;
            //  }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) +=  grad[25+dv];
            // if((myRank == 2) && (grad[dv+25] > 10.))
            // {
            //     std::cout<<"grad value Node: "<<dv+25<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[25+dv]<<std::endl;
            //  }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) +=  grad[29+dv];
            // if((myRank == 2) && (grad[dv+29] > 10.))
            // {
            //     std::cout<<"grad value Node: "<<dv+29<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[29+dv]<<std::endl;
            //  }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],  0);
        for(int dv=0;dv<4;dv++)
        {
          myGrid.value(index+dv) +=  grad[33+dv];
          // if((myRank == 2) && (grad[dv+33] > 10.))
          // {
          //     std::cout<<"grad value Node: "<<dv+33<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[33+dv]<<std::endl;
          //  }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) +=  grad[37+dv];
           // if((myRank == 2) && (grad[dv+37] > 10.))
           // {
           //    std::cout<<"grad value Node: "<<dv+37<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[37+dv]<<std::endl;
           // }
    }}


    template <int N,int numblock, typename dataType1>
    inline void copyFromArrayToCell2(lbmRD3Q41<dataType1> &lbModel ,dataType1 *grad, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3,int myRank)
    {
       myGrid(i1,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) += grad[lbModel.CENTER_DV_ZERO_ZERO_ZERO];
       // if((myRank == 2) && (grad[lbModel.CENTER_DV_ZERO_ZERO_ZERO] > 10.))
       // {
       //    std::cout<<"grad value Cell: "<<0<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[lbModel.CENTER_DV_ZERO_ZERO_ZERO]<<std::endl;
       // }

        int index = myGrid.getIndex(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv)  +=  grad[1+dv];
           // if((myRank == 2) && (grad[dv+1] > 10.))
           // {
           //    std::cout<<"grad value Cell: "<<dv+1<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[1+dv]<<std::endl;
           // }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) +=  grad[5+dv];
           // if((myRank == 2) && (grad[dv+5] > 10.))
           // {
           //    std::cout<<"grad value Cell: "<<dv+5<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[5+dv]<<std::endl;
           // }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) +=  grad[9+dv];
           // if((myRank == 2) && (grad[dv+9] > 10.))
           // {
           //    std::cout<<"grad value Cell: "<<dv+9<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[9+dv]<<std::endl;
           // }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) +=  grad[13+dv];
            // if((myRank == 2) && (grad[dv+13] > 10.))
            // {
            //     std::cout<<"grad value Cell: "<<dv+13<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[13+dv]<<std::endl;
            //  }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) +=  grad[17+dv];
           // if((myRank == 2) && (grad[dv+17] > 10.))
           // {
           //    std::cout<<"grad value Cell: "<<dv+17<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[17+dv]<<std::endl;
           // }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) +=  grad[21+dv];
            // if((myRank == 2) && (grad[dv+21] > 10.))
            // {
            //     std::cout<<"grad value Cell: "<<dv+21<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[21+dv]<<std::endl;
            //  }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) +=  grad[25+dv];
            // if((myRank == 2) && (grad[dv+25] > 10.))
            // {
            //     std::cout<<"grad value Cell: "<<dv+25<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[25+dv]<<std::endl;
            //  }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],  0);
        for(int dv=0;dv<4;dv++)
        {
            myGrid.value(index+dv) +=  grad[29+dv];
            // if((myRank == 2) && (grad[dv+29] > 10.))
            // {
            //     std::cout<<"grad value Cell: "<<dv+29<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[29+dv]<<std::endl;
            //  }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],  0);
        for(int dv=0;dv<4;dv++)
        {
          myGrid.value(index+dv) +=  grad[33+dv];
          // if((myRank == 2) && (grad[dv+33] > 10.))
          // {
          //     std::cout<<"grad value Cell: "<<dv+33<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[33+dv]<<std::endl;
          //  }
        }

        index = myGrid.getIndex(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],  0);
        for(int dv=0;dv<4;dv++)
        {
           myGrid.value(index+dv) +=  grad[37+dv];
           // if((myRank == 2) && (grad[dv+37] > 10.))
           // {
           //    std::cout<<"grad value Cell: "<<dv+37<<" "<<i1<<" "<<i2<<" "<<i3<<" "<<grad[37+dv]<<std::endl;
           // }
        }
      }



    template <int N,int numblock, typename dataType1>
    inline void leftInletGrad(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 u_inlet,int VECT_LENGTH)
    {
       dataType1 *rho,*uX,*uY,*uZ,*theta,*Pxx,*Pyy,*Pzz,*Pxy,*Pyz,*Pzx;
       rho   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       uX    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       uY    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       uZ    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       theta = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       Pxx   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       Pyy   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       Pzz   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       Pxy   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       Pyz   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       Pzx   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       dataType1 grad[lbModel.dvN];

        for(int k=myGrid.ndB3; k<=myGrid.ndE3;k++)
            for(int j=myGrid.ndB2; j<=myGrid.ndE2;j++)
          {
             copyFromNode(lbModel,myGrid,VECT_LENGTH,myGrid.nB1,j,k);
             getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);
//              Pxx[0] = Pxx[0] - rho[0]*(uX[0]*uX[0] + theta[0]) + (u_inlet*u_inlet + lbModel.theta0);
//              Pxx[1] = Pxx[1] - rho[1]*(uX[1]*uX[1] + theta[1]) + (u_inlet*u_inlet + lbModel.theta0);

             for(int i=0;i<VECT_LENGTH;i++)
             {
              uX   [i] = u_inlet;//0.0;//
              uY   [i] = 0.0;
              uZ   [i] = 0.0;
              rho  [i] = 1.0;
              theta[i] = lbModel.theta0;
             }

              getGradThermalSinglePoint(lbModel, rho[0], uX[0],uY[0],uZ[0],theta[0],Pxx[0],Pyy[0],Pzz[0],Pxy[0],Pyz[0],Pzx[0],grad);
              for(int count=0;count<41;count++)
              {
                if(grad[count]<0.0)
                  std::cout<<"Grad blew up: "<<j<<"   "<<k<<"   "<<count<<std::endl;
              }
              copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nB1-1,j,k);
              copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nB1-2,j,k);
              copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nB1-3,j,k);
              copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nB1-4,j,k);


              myGrid(myGrid.nB1-2 ,j  ,k  ,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO) =  grad[lbModel.DV_P2_ZERO_ZERO];

              getGradThermalSinglePoint(lbModel, rho[1], uX[1],uY[1],uZ[1],theta[1],Pxx[1],Pyy[1],Pzz[1],Pxy[1],Pyz[1],Pzx[1],grad);
              myGrid(myGrid.nB1-1 ,j  ,k  ,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO) =  grad[lbModel.DV_P2_ZERO_ZERO];
          }

        for(int k=myGrid.ndB3; k<=myGrid.ndE3;k++)
            for(int j=myGrid.ndB2; j<=myGrid.ndE2;j++)
          {
             copyFromCell(lbModel,myGrid,VECT_LENGTH,myGrid.nB1,j,k);
             getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);
//              Pxx[0] = Pxx[0] - rho[0]*(uX[0]*uX[0] + theta[0]) + (u_inlet*u_inlet + lbModel.theta0);
//              Pxx[1] = Pxx[1] - rho[1]*(uX[1]*uX[1] + theta[1]) + (u_inlet*u_inlet + lbModel.theta0);
             for(int i=0;i<VECT_LENGTH;i++)
             {
              uX   [0] = u_inlet;//0.0;//
              uY   [0] = 0.0;
              uZ   [0] = 0.0;
              rho  [0] = 1.0;
              theta[0] = lbModel.theta0;
             }
              getGradThermalSinglePoint(lbModel, rho[0], uX[0],uY[0],uZ[0],theta[0],Pxx[0],Pyy[0],Pzz[0],Pxy[0],Pyz[0],Pzx[0],grad);
              copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-1,j,k);
              copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-2,j,k);
              copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-3,j,k);
              copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-4,j,k);

              myGrid(myGrid.nB1-2 ,j  ,k  ,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO) =  grad[lbModel.DV_P2_ZERO_ZERO];

              getGradThermalSinglePoint(lbModel, rho[1], uX[1],uY[1],uZ[1],theta[1],Pxx[1],Pyy[1],Pzz[1],Pxy[1],Pyz[1],Pzx[1],grad);
              myGrid(myGrid.nB1-1 ,j  ,k  ,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.DV_P2_ZERO_ZERO) =  grad[lbModel.DV_P2_ZERO_ZERO];

          }
        }
//

    template <int N,int numblock, typename dataType1>
    inline void rightOutletGrad(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH)
    {
       dataType1 *rho,*uX,*uY,*uZ,*theta,*Pxx,*Pyy,*Pzz,*Pxy,*Pyz,*Pzx;
       rho   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       uX    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       uY    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       uZ    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       theta = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       Pxx   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       Pyy   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       Pzz   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       Pxy   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       Pyz   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       Pzx   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
       dataType1 grad[lbModel.dvN];

        for(int k=myGrid.ndB3; k<=myGrid.ndE3;k++)
            for(int j=myGrid.ndB2; j<=myGrid.ndE2;j++)
            {
             copyFromNode(lbModel,myGrid,VECT_LENGTH,myGrid.nE1-3,j,k);
             getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);
             getGradThermalSinglePoint(lbModel, rho[3], uX[3],uY[3],uZ[3],theta[3],Pxx[3],Pyy[3],Pzz[3],Pxy[3],Pyz[3],Pzx[3],grad);
             copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+1,j,k);
             copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+2,j,k);
             copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+3,j,k);
             copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+4,j,k);

             myGrid(myGrid.nE1+2 ,j  ,k  ,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO) =  grad[lbModel.DV_M2_ZERO_ZERO];
             getGradThermalSinglePoint(lbModel, rho[2], uX[2],uY[2],uZ[2],theta[2],Pxx[2],Pyy[2],Pzz[2],Pxy[2],Pyz[2],Pzx[2],grad);
             myGrid(myGrid.nE1+1 ,j  ,k  ,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO) =  grad[lbModel.DV_M2_ZERO_ZERO];


             copyFromCell(lbModel,myGrid,VECT_LENGTH,myGrid.nE1-3,j,k);
             getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);
             getGradThermalSinglePoint(lbModel, rho[3], uX[3],uY[3],uZ[3],theta[3],Pxx[3],Pyy[3],Pzz[3],Pxy[3],Pyz[3],Pzx[3],grad);
             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+1,j,k);
             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+2,j,k);
             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+3,j,k);
             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+4,j,k);
             myGrid(myGrid.nE1+2 ,j  ,k  ,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO) =  grad[lbModel.DV_M2_ZERO_ZERO];
             getGradThermalSinglePoint(lbModel, rho[2], uX[2],uY[2],uZ[2],theta[2],Pxx[2],Pyy[2],Pzz[2],Pxy[2],Pyz[2],Pzx[2],grad);
             myGrid(myGrid.nE1+1 ,j  ,k  ,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO) =  grad[lbModel.DV_M2_ZERO_ZERO];
            }
    }


    template <int N,int numblock, typename dataType1>
    inline void leftInlet(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,dataType1 u_inlet,int VECT_LENGTH)
    {
        dataType1 rho[VECT_LENGTH] __attribute__ ((aligned(32)));
        dataType1 uX[VECT_LENGTH]  __attribute__ ((aligned(32)));
        dataType1 uY[VECT_LENGTH]  __attribute__ ((aligned(32)));
        dataType1 uZ[VECT_LENGTH]  __attribute__ ((aligned(32)));
        dataType1 theta[VECT_LENGTH] __attribute__ ((aligned(32)));

        for(int k=0; k<myGrid.n3;k++)
         for(int j=0; j<myGrid.n2;j++)
            {
//              getHydroMomentsFromNode(lbModel,myGrid,myGrid.nB1,j,k,rho,uX,uY,uZ,theta);
             for(int i=0;i<VECT_LENGTH;i++)
             {
              uX   [i] = u_inlet;
              uY   [i] = 0.0;
              uZ   [i] = 0.0;
              theta[i] = lbModel.theta0;///rho[i];
              rho  [i] = 1.0;
             }
              getFEq(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);

              copyToNodeSinglePoint(lbModel,myGrid,myGrid.nB1-1 ,j,k,0);
              copyToNodeSinglePoint(lbModel,myGrid,myGrid.nB1-2 ,j,k,0);
              copyToNodeSinglePoint(lbModel,myGrid,myGrid.nB1-3 ,j,k,0);
              copyToNodeSinglePoint(lbModel,myGrid,myGrid.nB1-4 ,j,k,0);

              copyToCellSinglePoint(lbModel,myGrid,myGrid.nB1-1 ,j,k,0);
              copyToCellSinglePoint(lbModel,myGrid,myGrid.nB1-2 ,j,k,0);
              copyToCellSinglePoint(lbModel,myGrid,myGrid.nB1-3 ,j,k,0);
              copyToCellSinglePoint(lbModel,myGrid,myGrid.nB1-4 ,j,k,0);

            }

    }


    template <int N,int numblock, typename dataType1>
    inline void rightOutlet(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH)
    {
        dataType1 rho[VECT_LENGTH] __attribute__ ((aligned(32)));
        dataType1 uX[VECT_LENGTH]  __attribute__ ((aligned(32)));
        dataType1 uY[VECT_LENGTH]  __attribute__ ((aligned(32)));
        dataType1 uZ[VECT_LENGTH]  __attribute__ ((aligned(32)));
        dataType1 theta[VECT_LENGTH] __attribute__ ((aligned(32)));

        for(int k=0; k<myGrid.n3;k++)
         for(int j=0; j<myGrid.n2;j++)
            {
             getHydroMomentsFromCell(lbModel,myGrid,myGrid.nE1-3,j,k,rho,uX,uY,uZ,theta);
             getFEq(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);

             copyToNodeSinglePoint(lbModel,myGrid,myGrid.nE1+1,j,k,3);
             copyToNodeSinglePoint(lbModel,myGrid,myGrid.nE1+2,j,k,3);
             copyToNodeSinglePoint(lbModel,myGrid,myGrid.nE1+3,j,k,3);
             copyToNodeSinglePoint(lbModel,myGrid,myGrid.nE1+4,j,k,3);

             copyToCellSinglePoint(lbModel,myGrid,myGrid.nE1+1,j,k,3);
             copyToCellSinglePoint(lbModel,myGrid,myGrid.nE1+2,j,k,3);
             copyToCellSinglePoint(lbModel,myGrid,myGrid.nE1+3,j,k,3);
             copyToCellSinglePoint(lbModel,myGrid,myGrid.nE1+4,j,k,3);
            }
    }

    template <int N,int numblock, typename dataType1>
    inline void replaceInletPopulations(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)
    {
        for(int k=0; k<myGrid.n3;k++)
            for(int j=0; j<myGrid.n2;j++)
          {
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO)    =   myGrid(myGrid.nB1-2 ,j  ,k  ,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO);
           myGrid(myGrid.nB1+1,j  ,k  ,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO)    =   myGrid(myGrid.nB1-2 ,j  ,k  ,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO);

           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_P1_ZERO_P1)    =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_P1_ZERO_P1)  ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_P1_ZERO_M1)    =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_P1_ZERO_M1)  ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_P1_P1_ZERO)    =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_P1_P1_ZERO)  ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_P1_ZERO_ZERO)  =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_P1_ZERO_ZERO);
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_P1_M1_ZERO)    =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_P1_M1_ZERO)  ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_P_P_P)         =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_P_P_P)       ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_P_M_P)         =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_P_M_P)       ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_P_M_M)         =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_P_M_M)       ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_P_P_M)         =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_P_P_M)       ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_P1_P1_P1)      =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_P1_P1_P1)    ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_P1_M1_P1)      =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_P1_M1_P1)    ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1)     =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1)   ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1)     =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1)   ;

           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO)    =   myGrid(myGrid.nB1-2 ,j  ,k  ,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO);
           myGrid(myGrid.nB1+1,j  ,k  ,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO)    =   myGrid(myGrid.nB1-2 ,j  ,k  ,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO);

           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_P1_ZERO_P1)    =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_P1_ZERO_P1)  ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_P1_ZERO_M1)    =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_P1_ZERO_M1)  ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_P1_P1_ZERO)    =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_P1_P1_ZERO)  ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_P1_ZERO_ZERO)  =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_P1_ZERO_ZERO);
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_P1_M1_ZERO)    =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_P1_M1_ZERO)  ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_P1_P1_P1)      =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_P1_P1_P1)    ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_P1_M1_P1)      =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_P1_M1_P1)    ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1)     =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1)   ;
           myGrid(myGrid.nB1  ,j  ,k  ,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1)     =    myGrid(myGrid.nB1-2,j  ,k  ,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1)   ;
          }

    }

    template <int N,int numblock, typename dataType1>
    inline void replaceOutletPopulations(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)
    {
        for(int k=0; k<myGrid.n3;k++)
            for(int j=0; j<myGrid.n2;j++)
          {
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO)    =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO);
           myGrid(myGrid.nE1-1,j  ,k  ,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO)    =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO);

           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_M1_ZERO_P1)    =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_M1_ZERO_P1)  ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_M1_ZERO_M1)    =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_M1_ZERO_M1)  ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_M1_P1_ZERO)    =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_M1_P1_ZERO)  ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_M1_M1_ZERO)    =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_M1_M1_ZERO)  ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_M1_ZERO_ZERO)  =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_M1_ZERO_ZERO);
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_M1_P1_P1)      =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_M1_P1_P1)    ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_M1_M1_P1)      =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_M1_M1_P1)    ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1)     =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1)   ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1)     =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1)   ;

           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO)    =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO);
           myGrid(myGrid.nE1-1,j  ,k  ,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO)    =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO);

           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_M1_ZERO_P1)    =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_M1_ZERO_P1)   ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_M1_ZERO_M1)    =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_M1_ZERO_M1)   ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_M1_P1_ZERO)    =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_M1_P1_ZERO)   ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_M1_M1_ZERO)    =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_M1_M1_ZERO)   ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_M1_ZERO_ZERO)  =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_M1_ZERO_ZERO) ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_M_P_P)         =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_M_P_P)        ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_M_M_P)         =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_M_M_P)        ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_M_M_M)         =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_M_M_M)        ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_M_P_M)         =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_M_P_M)        ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_M1_P1_P1)      =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_M1_P1_P1)     ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_M1_M1_P1)      =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_M1_M1_P1)     ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1)     =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1)    ;
           myGrid(myGrid.nE1  ,j  ,k  ,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1)     =   myGrid(myGrid.nE1+2,j  ,k  ,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1)    ;
          }

    }

    template <int N,int numblock, typename dataType1>
    void topWall(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)
  {
    // Cell - All Plus
    for(int k=myGrid.nB3; k<=myGrid.nE3;k++)
      for(int i=myGrid.nB1; i<=myGrid.nE1;i++)
      {
        myGrid(i  ,myGrid.nE2+2,k  ,lbModel.G2 ,myGrid.cell[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO) =  myGrid(i,myGrid.nE2  ,k,lbModel.G1 ,myGrid.cell[lbModel.G1]  ,lbModel.G1_DV_ZERO_P2_ZERO);
        myGrid(i  ,myGrid.nE2+1,k  ,lbModel.G2 ,myGrid.cell[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO) =  myGrid(i,myGrid.nE2-1,k,lbModel.G1 ,myGrid.cell[lbModel.G1]  ,lbModel.G1_DV_ZERO_P2_ZERO);
        myGrid(i  ,myGrid.nE2+1,k+1,lbModel.G4 ,myGrid.cell[lbModel.G4]  ,lbModel.G4_DV_ZERO_M1_M1  ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G3 ,myGrid.cell[lbModel.G3]  ,lbModel.G3_DV_ZERO_P1_P1  );
        myGrid(i  ,myGrid.nE2+1,k-1,lbModel.G3 ,myGrid.cell[lbModel.G3]  ,lbModel.G3_DV_ZERO_M1_P1  ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G4 ,myGrid.cell[lbModel.G4]  ,lbModel.G4_DV_ZERO_P1_M1  );
        myGrid(i+1,myGrid.nE2+1,k  ,lbModel.G6 ,myGrid.cell[lbModel.G6]  ,lbModel.G6_DV_M1_M1_ZERO  ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G5 ,myGrid.cell[lbModel.G5]  ,lbModel.G5_DV_P1_P1_ZERO  );
        myGrid(i-1,myGrid.nE2+1,k  ,lbModel.G6 ,myGrid.cell[lbModel.G6]  ,lbModel.G6_DV_P1_M1_ZERO  ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G5 ,myGrid.cell[lbModel.G5]  ,lbModel.G5_DV_M1_P1_ZERO  );
        myGrid(i  ,myGrid.nE2+1,k  ,lbModel.G6 ,myGrid.cell[lbModel.G6]  ,lbModel.G6_DV_ZERO_M1_ZERO) =  myGrid(i,myGrid.nE2  ,k,lbModel.G5 ,myGrid.cell[lbModel.G5]  ,lbModel.G5_DV_ZERO_P1_ZERO);
        myGrid(i+1,myGrid.nE2+1,k+1,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_M1_M1_M1   ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G9 ,myGrid.cell[lbModel.G9]  ,lbModel.G9_DV_P1_P1_P1    );
        myGrid(i-1,myGrid.nE2+1,k+1,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_P1_M1_M1   ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G9 ,myGrid.cell[lbModel.G9]  ,lbModel.G9_DV_M1_P1_P1    );
        myGrid(i+1,myGrid.nE2+1,k-1,lbModel.G9 ,myGrid.cell[lbModel.G9]  ,lbModel.G9_DV_M1_M1_P1    ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_P1_P1_M1   );
        myGrid(i-1,myGrid.nE2+1,k-1,lbModel.G9 ,myGrid.cell[lbModel.G9]  ,lbModel.G9_DV_P1_M1_P1    ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_M1_P1_M1   );
        // BCC
        myGrid(i+1,myGrid.nE2+1,k+1,lbModel.G8,myGrid.node[lbModel.G8]   ,lbModel.G8_DV_M_M_M       ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G7,myGrid.cell[lbModel.G7]   ,lbModel.G7_DV_P_P_P       );
        myGrid(i  ,myGrid.nE2+1,k+1,lbModel.G8,myGrid.node[lbModel.G8]   ,lbModel.G8_DV_P_M_M       ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G7,myGrid.cell[lbModel.G7]   ,lbModel.G7_DV_M_P_P       );
        myGrid(i+1,myGrid.nE2+1,k  ,lbModel.G7,myGrid.node[lbModel.G7]   ,lbModel.G7_DV_M_M_P       ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G8,myGrid.cell[lbModel.G8]   ,lbModel.G8_DV_P_P_M       );
        myGrid(i  ,myGrid.nE2+1,k  ,lbModel.G7,myGrid.node[lbModel.G7]   ,lbModel.G7_DV_P_M_P       ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G8,myGrid.cell[lbModel.G8]   ,lbModel.G8_DV_M_P_M       );
      }

      // Node - All Plus except BCC
  for(int k=myGrid.nB3; k<=myGrid.nE3;k++)
      for(int i=myGrid.nB1; i<=myGrid.nE1;i++)
        {
          myGrid(i  ,myGrid.nE2+2,k  ,lbModel.G2 ,myGrid.node[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO) =  myGrid(i,myGrid.nE2  ,k,lbModel.G1 ,myGrid.node[lbModel.G1]  ,lbModel.G1_DV_ZERO_P2_ZERO);
          myGrid(i  ,myGrid.nE2+1,k  ,lbModel.G2 ,myGrid.node[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO) =  myGrid(i,myGrid.nE2-1,k,lbModel.G1 ,myGrid.node[lbModel.G1]  ,lbModel.G1_DV_ZERO_P2_ZERO);
          myGrid(i  ,myGrid.nE2+1,k+1,lbModel.G4 ,myGrid.node[lbModel.G4]  ,lbModel.G4_DV_ZERO_M1_M1  ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G3 ,myGrid.node[lbModel.G3]  ,lbModel.G3_DV_ZERO_P1_P1  );
          myGrid(i  ,myGrid.nE2+1,k-1,lbModel.G3 ,myGrid.node[lbModel.G3]  ,lbModel.G3_DV_ZERO_M1_P1  ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G4 ,myGrid.node[lbModel.G4]  ,lbModel.G4_DV_ZERO_P1_M1  );
          myGrid(i+1,myGrid.nE2+1,k  ,lbModel.G6 ,myGrid.node[lbModel.G6]  ,lbModel.G6_DV_M1_M1_ZERO  ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G5 ,myGrid.node[lbModel.G5]  ,lbModel.G5_DV_P1_P1_ZERO  );
          myGrid(i-1,myGrid.nE2+1,k  ,lbModel.G6 ,myGrid.node[lbModel.G6]  ,lbModel.G6_DV_P1_M1_ZERO  ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G5 ,myGrid.node[lbModel.G5]  ,lbModel.G5_DV_M1_P1_ZERO  );
          myGrid(i  ,myGrid.nE2+1,k  ,lbModel.G6 ,myGrid.node[lbModel.G6]  ,lbModel.G6_DV_ZERO_M1_ZERO) =  myGrid(i,myGrid.nE2  ,k,lbModel.G5 ,myGrid.node[lbModel.G5]  ,lbModel.G5_DV_ZERO_P1_ZERO);
          myGrid(i+1,myGrid.nE2+1,k+1,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_M1_M1_M1   ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G9 ,myGrid.node[lbModel.G9]  ,lbModel.G9_DV_P1_P1_P1    );
          myGrid(i-1,myGrid.nE2+1,k+1,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_P1_M1_M1   ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G9 ,myGrid.node[lbModel.G9]  ,lbModel.G9_DV_M1_P1_P1    );
          myGrid(i+1,myGrid.nE2+1,k-1,lbModel.G9 ,myGrid.node[lbModel.G9]  ,lbModel.G9_DV_M1_M1_P1    ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_P1_P1_M1   );
          myGrid(i-1,myGrid.nE2+1,k-1,lbModel.G9 ,myGrid.node[lbModel.G9]  ,lbModel.G9_DV_P1_M1_P1    ) =  myGrid(i,myGrid.nE2  ,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_M1_P1_M1   );
        }

  }


  template <int N,int numblock, typename dataType1>
  void bottomWall(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)
  {

    // Cell - All Plus
    for(int k=myGrid.nB3; k<=myGrid.nE3;k++)
      for(int i=myGrid.nB1; i<=myGrid.nE1;i++)
      {
        myGrid(i  ,myGrid.nB2-2,k  ,lbModel.G1 ,myGrid.node[lbModel.G1]  ,lbModel.G1_DV_ZERO_P2_ZERO) =  myGrid(i,myGrid.nB2  ,k,lbModel.G2 ,myGrid.node[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO);
        myGrid(i  ,myGrid.nB2-1,k  ,lbModel.G1 ,myGrid.node[lbModel.G1]  ,lbModel.G1_DV_ZERO_P2_ZERO) =  myGrid(i,myGrid.nB2+1,k,lbModel.G2 ,myGrid.node[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO);
        myGrid(i  ,myGrid.nB2-1,k-1,lbModel.G3 ,myGrid.node[lbModel.G3]  ,lbModel.G3_DV_ZERO_P1_P1  ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G4 ,myGrid.node[lbModel.G4]  ,lbModel.G4_DV_ZERO_M1_M1  );
        myGrid(i  ,myGrid.nB2-1,k+1,lbModel.G4 ,myGrid.node[lbModel.G4]  ,lbModel.G4_DV_ZERO_P1_M1  ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G3 ,myGrid.node[lbModel.G3]  ,lbModel.G3_DV_ZERO_M1_P1  );
        myGrid(i-1,myGrid.nB2-1,k  ,lbModel.G5 ,myGrid.node[lbModel.G5]  ,lbModel.G5_DV_P1_P1_ZERO  ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G6 ,myGrid.node[lbModel.G6]  ,lbModel.G6_DV_M1_M1_ZERO  );
        myGrid(i+1,myGrid.nB2-1,k  ,lbModel.G5 ,myGrid.node[lbModel.G5]  ,lbModel.G5_DV_M1_P1_ZERO  ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G6 ,myGrid.node[lbModel.G6]  ,lbModel.G6_DV_P1_M1_ZERO  );
        myGrid(i  ,myGrid.nB2-1,k  ,lbModel.G5 ,myGrid.node[lbModel.G5]  ,lbModel.G5_DV_ZERO_P1_ZERO) =  myGrid(i,myGrid.nB2  ,k,lbModel.G6 ,myGrid.node[lbModel.G6]  ,lbModel.G6_DV_ZERO_M1_ZERO);
        myGrid(i-1,myGrid.nB2-1,k-1,lbModel.G9 ,myGrid.node[lbModel.G9]  ,lbModel.G9_DV_P1_P1_P1    ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_M1_M1_M1   );
        myGrid(i+1,myGrid.nB2-1,k-1,lbModel.G9 ,myGrid.node[lbModel.G9]  ,lbModel.G9_DV_M1_P1_P1    ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_P1_M1_M1   );
        myGrid(i-1,myGrid.nB2-1,k+1,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_P1_P1_M1   ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G9 ,myGrid.node[lbModel.G9]  ,lbModel.G9_DV_M1_M1_P1    );
        myGrid(i+1,myGrid.nB2-1,k+1,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_M1_P1_M1   ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G9 ,myGrid.node[lbModel.G9]  ,lbModel.G9_DV_P1_M1_P1    );
        // BCC             nB2
        myGrid(i-1,myGrid.nB2-1,k-1,lbModel.G7,myGrid.cell[lbModel.G7]   ,lbModel.G7_DV_P_P_P       ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G8,myGrid.node[lbModel.G8]   ,lbModel.G8_DV_M_M_M       );
        myGrid(i  ,myGrid.nB2-1,k-1,lbModel.G7,myGrid.cell[lbModel.G7]   ,lbModel.G7_DV_M_P_P       ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G8,myGrid.node[lbModel.G8]   ,lbModel.G8_DV_P_M_M       );
        myGrid(i-1,myGrid.nB2-1,k  ,lbModel.G8,myGrid.cell[lbModel.G8]   ,lbModel.G8_DV_P_P_M       ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G7,myGrid.node[lbModel.G7]   ,lbModel.G7_DV_M_M_P       );
        myGrid(i  ,myGrid.nB2-1,k  ,lbModel.G8,myGrid.cell[lbModel.G8]   ,lbModel.G8_DV_M_P_M       ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G7,myGrid.node[lbModel.G7]   ,lbModel.G7_DV_P_M_P       );
      }

      // Node - All Plus except BCC
    for(int k=myGrid.nB3; k<=myGrid.nE3;k++)
      for(int i=myGrid.nB1; i<=myGrid.nE1;i++)
        {
          myGrid(i  ,myGrid.nB2-2,k  ,lbModel.G1 ,myGrid.cell[lbModel.G1]  ,lbModel.G1_DV_ZERO_P2_ZERO) =  myGrid(i,myGrid.nB2  ,k,lbModel.G2 ,myGrid.cell[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO);
          myGrid(i  ,myGrid.nB2-1,k  ,lbModel.G1 ,myGrid.cell[lbModel.G1]  ,lbModel.G1_DV_ZERO_P2_ZERO) =  myGrid(i,myGrid.nB2+1,k,lbModel.G2 ,myGrid.cell[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO);
          myGrid(i  ,myGrid.nB2-1,k-1,lbModel.G3 ,myGrid.cell[lbModel.G3]  ,lbModel.G3_DV_ZERO_P1_P1  ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G4 ,myGrid.cell[lbModel.G4]  ,lbModel.G4_DV_ZERO_M1_M1  );
          myGrid(i  ,myGrid.nB2-1,k+1,lbModel.G4 ,myGrid.cell[lbModel.G4]  ,lbModel.G4_DV_ZERO_P1_M1  ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G3 ,myGrid.cell[lbModel.G3]  ,lbModel.G3_DV_ZERO_M1_P1  );
          myGrid(i-1,myGrid.nB2-1,k  ,lbModel.G5 ,myGrid.cell[lbModel.G5]  ,lbModel.G5_DV_P1_P1_ZERO  ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G6 ,myGrid.cell[lbModel.G6]  ,lbModel.G6_DV_M1_M1_ZERO  );
          myGrid(i+1,myGrid.nB2-1,k  ,lbModel.G5 ,myGrid.cell[lbModel.G5]  ,lbModel.G5_DV_M1_P1_ZERO  ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G6 ,myGrid.cell[lbModel.G6]  ,lbModel.G6_DV_P1_M1_ZERO  );
          myGrid(i  ,myGrid.nB2-1,k  ,lbModel.G5 ,myGrid.cell[lbModel.G5]  ,lbModel.G5_DV_ZERO_P1_ZERO) =  myGrid(i,myGrid.nB2  ,k,lbModel.G6 ,myGrid.cell[lbModel.G6]  ,lbModel.G6_DV_ZERO_M1_ZERO);
          myGrid(i-1,myGrid.nB2-1,k-1,lbModel.G9 ,myGrid.cell[lbModel.G9]  ,lbModel.G9_DV_P1_P1_P1    ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_M1_M1_M1   );
          myGrid(i+1,myGrid.nB2-1,k-1,lbModel.G9 ,myGrid.cell[lbModel.G9]  ,lbModel.G9_DV_M1_P1_P1    ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_P1_M1_M1   );
          myGrid(i-1,myGrid.nB2-1,k+1,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_P1_P1_M1   ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G9 ,myGrid.cell[lbModel.G9]  ,lbModel.G9_DV_M1_M1_P1    );
          myGrid(i+1,myGrid.nB2-1,k+1,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_M1_P1_M1   ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G9 ,myGrid.cell[lbModel.G9]  ,lbModel.G9_DV_P1_M1_P1    );
        }

    }


// Explicit declarations

    template void topWall<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &);
    template void bottomWall<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &);
    template void getGradThermalSinglePoint<double>(lbmRD3Q41<double>&,double,double,double,double,double,double,double,double,double,double,double,double *);
    template void  getGradThermalSinglePoint_heatAlpha<double> (lbmRD3Q41<double>&,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double *);
    template void copyFromArrayToNode<4,11,double>(lbmRD3Q41<double> & ,double *, gridBCC3D<4,11,double> &,int ,int ,int );
    template void copyFromArrayToCell<4,11,double>(lbmRD3Q41<double> & ,double *, gridBCC3D<4,11,double> &,int ,int ,int );
    template void copyFromArrayToNode2<4,11,double>(lbmRD3Q41<double> & ,double *, gridBCC3D<4,11,double> &,int ,int ,int ,int );
    template void copyFromArrayToCell2<4,11,double>(lbmRD3Q41<double> & ,double *, gridBCC3D<4,11,double> &,int ,int ,int ,int );
