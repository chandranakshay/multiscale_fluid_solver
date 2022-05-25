#include"noiseFunctions.h"


  template <int N,int numblock, typename dataType1>
  void curl_noise(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int nP1, int nP2, int nP3, dataType1 &normg,int *coordinates,int VECT_LENGTH, dataType1 density, dataType1 uRef)
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


    int i1, i2, i3;
    int nX1 = myGrid.m1*nP1;
    int nX2 = myGrid.m2*nP2;
    int nX3 = myGrid.m3*nP3;

    double Lx, Ly, Lz;

    Lx = 2.0*(nX1/nX2);
    Ly = 2.0*(nX2/nX2);
    Lz = 2.0*(nX3/nX2);

    #define  iXMacro(i)  (((double)(i-0.5)/(nX1))*Lx)
    #define  iYMacro(j)  (((double)(j-0.5)/(nX2))*Ly)
    #define  iZMacro(k)  (((double)(k-0.5)/(nX3))*Lz)
    #define ramp(r)  ((15.0 - 10.0*r*r + 3.0*r*r*r*r)*r/8.0)

    double oneBydx(1.0/nX1);//(nX1/Lx);
    double oneBydy(1.0/nX2);//(nX2/Ly);
    double oneBydz(1.0/nX3);//(nX3/Lz);

    double norm = 0 ;

    for(int i3= myGrid.nB3;  i3<= myGrid.nE3;  i3++)
      for(int i2= myGrid.nB2;  i2<= myGrid.nE2;  i2++)
        for(int i1= myGrid.nB1;  i1<= myGrid.nE1;  i1=i1+4)
        {
          int kGlobal = i3 + myGrid.m3*coordinates[2] ;
          int jGlobal = i2 + myGrid.m2*coordinates[1] ;
          for(int index=0;index<VECT_LENGTH;index++)
          {
            int iGlobal = (i1 + index) + myGrid.m1*coordinates[0] ;

            copyFromNode(lbModel,myGrid,VECT_LENGTH,i1,i2,i3);
            //         getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);

            if(jGlobal == myGrid.nB3 || jGlobal == myGrid.m2*nP2+3)
            {
              uX[index] = 0.0;
              uY[index] = 0.0;
              uZ[index] = 0.0 ;
            }

            else if(jGlobal <= (int)(0.5*(nX2+3)) )
            {
              uX[index] = ramp(iYMacro(jGlobal))  *( oneBydy* ( noise3(iXMacro(iGlobal)    ,iYMacro(jGlobal+1), iZMacro(kGlobal)  ) - noise3(iXMacro(iGlobal)    , iYMacro(jGlobal-1), iZMacro(kGlobal)  ) )   -  oneBydz*( noise3(iXMacro(iGlobal)    ,iYMacro(jGlobal)  ,iZMacro(kGlobal+1)) - noise3(iXMacro(iGlobal)     ,iYMacro(jGlobal)  , iZMacro(kGlobal-1)) ) );
              uY[index] = ramp(iYMacro(jGlobal))  *( oneBydz* ( noise3(iXMacro(iGlobal)    ,iYMacro(jGlobal)  , iZMacro(kGlobal+1)) - noise3(iXMacro(iGlobal)    , iYMacro(jGlobal)  , iZMacro(kGlobal-1)) )   -  oneBydx*( noise3(iXMacro(iGlobal+1)  ,iYMacro(jGlobal)  ,iZMacro(kGlobal)  ) - noise3(iXMacro(iGlobal-1)   ,iYMacro(jGlobal)  , iZMacro(kGlobal-1)) ) );
              uZ[index] = ramp(iYMacro(jGlobal))  *( oneBydx* ( noise3(iXMacro(iGlobal+1)  ,iYMacro(jGlobal)  , iZMacro(kGlobal)  ) - noise3(iXMacro(iGlobal-1)  , iYMacro(jGlobal)  , iZMacro(kGlobal)  ) )   -  oneBydy*( noise3(iXMacro(iGlobal)    ,iYMacro(jGlobal+1),iZMacro(kGlobal)  ) - noise3(iXMacro(iGlobal)     ,iYMacro(jGlobal-1), iZMacro(kGlobal)  ) ) );
            }
            else
            {
              uX[index] = ramp(iYMacro(nX2-jGlobal))  *( oneBydy* ( noise3(iXMacro(iGlobal)    ,iYMacro(jGlobal+1), iZMacro(kGlobal)  ) - noise3(iXMacro(iGlobal)    , iYMacro(jGlobal-1), iZMacro(kGlobal)  ) )   -   oneBydz* ( noise3(iXMacro(iGlobal)    ,iYMacro(jGlobal)  ,iZMacro(kGlobal+1)) - noise3(iXMacro(iGlobal)     ,iYMacro(jGlobal)  , iZMacro(kGlobal-1)) ) );
              uY[index] = ramp(iYMacro(nX2-jGlobal))  *( oneBydz* ( noise3(iXMacro(iGlobal)    ,iYMacro(jGlobal)  , iZMacro(kGlobal+1)) - noise3(iXMacro(iGlobal)    , iYMacro(jGlobal)  , iZMacro(kGlobal-1)) )   -   oneBydx* ( noise3(iXMacro(iGlobal+1)  ,iYMacro(jGlobal)  ,iZMacro(kGlobal)  ) - noise3(iXMacro(iGlobal-1)   ,iYMacro(jGlobal)  , iZMacro(kGlobal-1)) ) );
              uZ[index] = ramp(iYMacro(nX2+3-jGlobal))*( oneBydx* ( noise3(iXMacro(iGlobal+1)  ,iYMacro(jGlobal)  , iZMacro(kGlobal)  ) - noise3(iXMacro(iGlobal-1)  , iYMacro(jGlobal)  , iZMacro(kGlobal)  ) )   -   oneBydy* ( noise3(iXMacro(iGlobal)    ,iYMacro(jGlobal+1),iZMacro(kGlobal)  ) - noise3(iXMacro(iGlobal)     ,iYMacro(jGlobal-1), iZMacro(kGlobal)  ) ) );
            }
            norm += (uX[index]*uX[index] + uY[index]*uY[index] + uZ[index]*uZ[index]) ;

            rho[index]   = density;
            theta[index] = lbModel.theta0;
          }
          getFEq(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);
          copyToNode(lbModel,myGrid,VECT_LENGTH,i1,i2,i3);
          copyToCell(lbModel,myGrid,VECT_LENGTH,i1,i2,i3);
        }
        MPI_Allreduce(&norm,&normg,1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
        normg /= (nX1*nX2*nX3);
        normg = sqrt(normg) ;
}










































































// template <int N,int numblock, typename dataType1>
// void curl_noise_communicateAfterThis(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &testGrid, int nP1, int nP2, int nP3, dataType1 &normg,int *coordinates,int VECT_LENGTH,double Lx, double Ly, double Lz)
// {
//   int nX1 = testGrid.m1*nP1 + 2;
//   int nX2 = testGrid.m2*nP2 + 2;
//   int nX3 = testGrid.m3*nP3 + 2;
//
//
//
//   #define  iXMacro(i)  (((double)(i)/nX1)*Lx)
//   #define  iYMacro(j)  (((double)(j)/nX2)*Ly)
//   #define  iZMacro(k)  (((double)(k)/nX3)*Lz)
//   #define ramp(r)  ((15.0 - 10.0*r*r + 3.0*r*r*r*r)*r/8.0)
//
//   // Calculate Phis
//       for(int i3= testGrid.ndB3;  i3<= testGrid.ndE3;  i3++)
//         for(int i2= testGrid.ndB2;  i2<= testGrid.ndE2;  i2++)
//           for(int i1= testGrid.ndB1;  i1<= testGrid.ndE1;  i1++)
//        {
//         int iGlobal = i1 + testGrid.m1*coordinates[0] - 1 ;
//         int jGlobal = i2 + testGrid.m2*coordinates[1] - 1 ;
//         int kGlobal = i3 + testGrid.m3*coordinates[2] - 1 ;
//
//           testGrid(i1  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) = noise3(iXMacro(iGlobal)+Lx , iYMacro(jGlobal)    , iZMacro(kGlobal)    );
//           testGrid(i1  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) = noise3(iXMacro(iGlobal    ), iYMacro(jGlobal)+Ly , iZMacro(kGlobal)    );
//           testGrid(i1  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) = noise3(iXMacro(iGlobal    ), iYMacro(jGlobal)    , iZMacro(kGlobal)+Lz );
//        }
// }
//
// template <int N,int numblock, int numblock1, typename dataType1>
// void curl_noise_communiateBeforeThis(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, gridBCC3D<N, numblock1, dataType1> &testGrid, int nP1, int nP2, int nP3, dataType1 &normg,int *coordinates,int VECT_LENGTH,double Lx, double Ly, double Lz)
// {
//   dataType1 *rho,*uX,*uY,*uZ,*theta,*Pxx,*Pyy,*Pzz,*Pxy,*Pyz,*Pzx;
//   rho   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
//   uX    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
//   uY    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
//   uZ    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
//   theta = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
//
//   int nX1 = testGrid.m1*nP1;
//   int nX2 = testGrid.m2*nP2;
//   int nX3 = testGrid.m3*nP3;
//
//       double norm = 0 ;
//
//       for(int i3= testGrid.nB3;  i3<= testGrid.nE3;  i3++)
//         for(int i2= testGrid.nB2;  i2<= testGrid.nE2;  i2++)
//           for(int i1= testGrid.nB1;  i1<= testGrid.nE1;  i1=i1+4)
//           {
//             int kGlobal = i3 + testGrid.m3*coordinates[2] - 1 ;
//             int jGlobal = i2 + testGrid.m2*coordinates[1] - 1 ;
// //             int iGlobal = i1-testGrid.nB1 + testGrid.m1*coordinates[0] ;
//
//             for(int index=0;index<VECT_LENGTH;index++)
//             {
//               int iGlobal = (i1 + index) + testGrid.m1*coordinates[0] - 1 ;
//
//
//
//
// //               if (jGlobal <= (int)( 0.5*(testGrid.m2*nP2)) )
// //               {
// //                 uX[index] =    ( ramp(iYMacro(jGlobal  )) * testGrid(i1+index  ,i2  ,i3+1,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) - ramp(iYMacro(jGlobal)) * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) )
// //                              - ( ramp(iYMacro(jGlobal+1)) * testGrid(i1+index  ,i2+1,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) - ramp(iYMacro(jGlobal)) * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) );
// //
// //                 uY[index] =    ( ramp(iYMacro(jGlobal  )) * testGrid(i1+index+1,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) - ramp(iYMacro(jGlobal)) * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) )
// //                              - ( ramp(iYMacro(jGlobal  )) * testGrid(i1+index  ,i2  ,i3+1,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) - ramp(iYMacro(jGlobal)) * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) );
// //
// //                 uZ[index] =    ( ramp(iYMacro(jGlobal+1)) * testGrid(i1+index  ,i2+1,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) - ramp(iYMacro(jGlobal)) * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) )
// //                              - ( ramp(iYMacro(jGlobal  )) * testGrid(i1+index+1,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) - ramp(iYMacro(jGlobal)) * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) );
// //                }
// //
// //               if (jGlobal > (int)( 0.5*(testGrid.m2*nP2)) )
// //               {
// //                  uX[index] =    ( ramp(iYMacro(testGrid.m2*nP2 - jGlobal  )) * testGrid(i1+index  ,i2  ,i3+1,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) -  ramp(iYMacro(testGrid.m2*nP2 - jGlobal))  * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) )
// //                               - ( ramp(iYMacro(testGrid.m2*nP2 - jGlobal+1)) * testGrid(i1+index  ,i2+1,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) -  ramp(iYMacro(testGrid.m2*nP2 - jGlobal))  * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) );
// //
// //                  uY[index] =    ( ramp(iYMacro(testGrid.m2*nP2 - jGlobal  )) * testGrid(i1+index+1,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) -  ramp(iYMacro(testGrid.m2*nP2 - jGlobal))  * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) )
// //                               - ( ramp(iYMacro(testGrid.m2*nP2 - jGlobal  )) * testGrid(i1+index  ,i2  ,i3+1,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) -  ramp(iYMacro(testGrid.m2*nP2 - jGlobal))  * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) ) ;
// //
// //                  uZ[index] =    ( ramp(iYMacro(testGrid.m2*nP2 - jGlobal+1)) * testGrid(i1+index  ,i2+1,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) -  ramp(iYMacro(testGrid.m2*nP2 - jGlobal))  * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) )
// //                               - ( ramp(iYMacro(testGrid.m2*nP2 - jGlobal  )) * testGrid(i1+index+1,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) -  ramp(iYMacro(testGrid.m2*nP2 - jGlobal))  * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) );
// //                }
//               if (jGlobal <= (int)( 0.5*(testGrid.m2*nP2)) )
//               {
//                 uX[index] =    ( ramp(iYMacro(jGlobal  )) * testGrid(i1+index  ,i2  ,i3+1,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) - ramp(iYMacro(jGlobal)) * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) )
//                              - ( ramp(iYMacro(jGlobal  )) * testGrid(i1+index  ,i2+1,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) - ramp(iYMacro(jGlobal)) * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) );
//
//                 uY[index] =    ( ramp(iYMacro(jGlobal  )) * testGrid(i1+index+1,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) - ramp(iYMacro(jGlobal)) * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) )
//                              - ( ramp(iYMacro(jGlobal  )) * testGrid(i1+index  ,i2  ,i3+1,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) - ramp(iYMacro(jGlobal)) * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) );
//
//                 uZ[index] =    ( ramp(iYMacro(jGlobal  )) * testGrid(i1+index  ,i2+1,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) - ramp(iYMacro(jGlobal)) * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) )
//                              - ( ramp(iYMacro(jGlobal  )) * testGrid(i1+index+1,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) - ramp(iYMacro(jGlobal)) * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) );
//                }
//
//               if (jGlobal > (int)( 0.5*(testGrid.m2*nP2)) )
//               {
//                  uX[index] =    ( ramp(iYMacro(testGrid.m2*nP2 - jGlobal  )) * testGrid(i1+index  ,i2  ,i3+1,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) -  ramp(iYMacro(testGrid.m2*nP2 - jGlobal))  * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) )
//                               - ( ramp(iYMacro(testGrid.m2*nP2 - jGlobal  )) * testGrid(i1+index  ,i2+1,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) -  ramp(iYMacro(testGrid.m2*nP2 - jGlobal))  * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) );
//
//                  uY[index] =    ( ramp(iYMacro(testGrid.m2*nP2 - jGlobal  )) * testGrid(i1+index+1,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) -  ramp(iYMacro(testGrid.m2*nP2 - jGlobal))  * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) )
//                               - ( ramp(iYMacro(testGrid.m2*nP2 - jGlobal  )) * testGrid(i1+index  ,i2  ,i3+1,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) -  ramp(iYMacro(testGrid.m2*nP2 - jGlobal))  * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) ) ;
//
//                  uZ[index] =    ( ramp(iYMacro(testGrid.m2*nP2 - jGlobal  )) * testGrid(i1+index  ,i2+1,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) -  ramp(iYMacro(testGrid.m2*nP2 - jGlobal))  * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) )
//                               - ( ramp(iYMacro(testGrid.m2*nP2 - jGlobal  )) * testGrid(i1+index+1,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) -  ramp(iYMacro(testGrid.m2*nP2 - jGlobal))  * testGrid(i1+index  ,i2  ,i3  ,lbModel.G0,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) );
//                }
//
//
//              if(jGlobal == testGrid.nB2 || jGlobal == testGrid.nB2+testGrid.m2*nP2-1)
//              {
//                uX[index] = 0.0;
//                uY[index] = 0.0;
//                uZ[index] = 0.0 ;
//              }
//
//               norm += (uX[index]*uX[index] + uY[index]*uY[index] + uZ[index]*uZ[index]) ;
//
//               rho[index]   = 1.0;
//               theta[index] = lbModel.theta0;
//             }
//             getFEq(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);
//             copyToNode(lbModel,myGrid,VECT_LENGTH,i1,i2,i3);
//             copyToCell(lbModel,myGrid,VECT_LENGTH,i1,i2,i3);
//           }
//           MPI::COMM_WORLD.Allreduce(&norm,&normg,1, MPI::DOUBLE, MPI::SUM);
//           normg /= (nX1*nX2*nX3);
//           normg = sqrt(normg) ;
// }
//
//   template <int N,int numblock, typename dataType1>
//   inline void averageVelocityAndTKE(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int nP1, int nP2, int nP3,int *coordinates,int VECT_LENGTH,dataType1 F1,dataType1 F2,dataType1 F3,int step,int myRank,int size,dataType1 U0)
//   {
//     dataType1 *rho,*uX,*uY,*uZ,*theta,*Pxx,*Pyy,*Pzz,*Pxy,*Pyz,*Pzx;
//     rho   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
//     uX    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
//     uY    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
//     uZ    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
//     theta = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
//
//     int nX1 = myGrid.m1*nP1;
//     int nX2 = myGrid.m2*nP2;
//     int nX3 = myGrid.m3*nP3;
//
//         double TKE  = 0.0 ;         double TKEGlobalNode = 0.0 ;  double TKEGlobalCell = 0.0 ;
//
//         double KE  = 0.0 ;         double KEGlobalNode = 0.0 ;  double KEGlobalCell = 0.0 ;
//
//         double uAvgNode = 0.0 ;     double uAvgGlobalNode = 0.0 ;         double uAvgCell = 0.0 ;     double uAvgGlobalCell = 0.0 ;
//         double vAvgNode = 0.0 ;     double vAvgGlobalNode = 0.0 ;         double vAvgCell = 0.0 ;     double vAvgGlobalCell = 0.0 ;
//         double wAvgNode = 0.0 ;     double wAvgGlobalNode = 0.0 ;         double wAvgCell = 0.0 ;     double wAvgGlobalCell = 0.0 ;
//
//         for(int i3= myGrid.nB3;  i3<= myGrid.nE3;  i3++)
//           for(int i2= myGrid.nB2;  i2<= myGrid.nE2;  i2++)
//             for(int i1= myGrid.nB1;  i1<= myGrid.nE1;  i1=i1+4)
//             {
//              getHydroMomentsFromNodeWithForce(lbModel,myGrid,VECT_LENGTH,i1,i2,i3,rho,uX,uY,uZ,theta,F1,F2,F3);
//              for(int index=0;index<VECT_LENGTH;index++)
//              {
//               uAvgNode += uX[index];
//               vAvgNode += uY[index];
//               wAvgNode += uZ[index];
//               KE += ( uX[index]*uX[index] + uY[index]*uY[index] + uZ[index]*uZ[index]);
//              }
//             }
//
//             MPI::COMM_WORLD.Allreduce(&uAvgNode,&uAvgGlobalNode,1, MPI::DOUBLE, MPI::SUM);
//             MPI::COMM_WORLD.Allreduce(&vAvgNode,&vAvgGlobalNode,1, MPI::DOUBLE, MPI::SUM);
//             MPI::COMM_WORLD.Allreduce(&wAvgNode,&wAvgGlobalNode,1, MPI::DOUBLE, MPI::SUM);
//
//             MPI::COMM_WORLD.Allreduce(&KE,&KEGlobalNode,1, MPI::DOUBLE, MPI::SUM);
//
//
//         KE = 0.0;
//         for(int i3= myGrid.nB3;  i3<= myGrid.nE3;  i3++)
//           for(int i2= myGrid.nB2;  i2<= myGrid.nE2;  i2++)
//             for(int i1= myGrid.nB1;  i1<= myGrid.nE1;  i1=i1+4)
//             {
//              getHydroMomentsFromCellWithForce(lbModel,myGrid,VECT_LENGTH,i1,i2,i3,rho,uX,uY,uZ,theta,F1,F2,F3);
//              for(int index=0;index<VECT_LENGTH;index++)
//              {
//               uAvgCell += uX[index];
//               vAvgCell += uY[index];
//               wAvgCell += uZ[index];
//               KE += ( uX[index]*uX[index] + uY[index]*uY[index] + uZ[index]*uZ[index]);
//              }
//             }
//
//
//             MPI::COMM_WORLD.Allreduce(&uAvgCell,&uAvgGlobalCell,1, MPI::DOUBLE, MPI::SUM);
//             MPI::COMM_WORLD.Allreduce(&vAvgCell,&vAvgGlobalCell,1, MPI::DOUBLE, MPI::SUM);
//             MPI::COMM_WORLD.Allreduce(&wAvgCell,&wAvgGlobalCell,1, MPI::DOUBLE, MPI::SUM);
//
//
//             MPI::COMM_WORLD.Allreduce(&KE,&KEGlobalCell,1, MPI::DOUBLE, MPI::SUM);
//
// 	   double uAvgGlobal =  (uAvgGlobalNode +  uAvgGlobalCell)/(2.0*nX1*nX2*nX3);
// 	   double vAvgGlobal =  (vAvgGlobalNode +  vAvgGlobalCell)/(2.0*nX1*nX2*nX3);
// 	   double wAvgGlobal =  (wAvgGlobalNode +  wAvgGlobalCell)/(2.0*nX1*nX2*nX3);
//
//
// 	   KE =  (KEGlobalNode +  KEGlobalCell)/(nX1*nX2*nX3*U0*U0);
//
// /*
//         for(int i3= myGrid.nB3;  i3<= myGrid.nE3;  i3++)
//           for(int i2= myGrid.nB2;  i2<= myGrid.nE2;  i2++)
//             for(int i1= myGrid.nB1;  i1<= myGrid.nE1;  i1=i1+4)
//             {
//              getHydroMomentsFromNodeWithForce(lbModel,myGrid,VECT_LENGTH,i1,i2,i3,rho,uX,uY,uZ,theta,F1,F2,F3);
//              for(int index=0;index<VECT_LENGTH;index++)
//                TKE += ( (uX[index]-uAvgGlobal)*(uX[index]-uAvgGlobal) + (uY[index]-vAvgGlobal)*(uY[index]-vAvgGlobal) + (uZ[index]-wAvgGlobal)*(uZ[index]-wAvgGlobal) );
//              }
//        MPI::COMM_WORLD.Allreduce(&TKE,&TKEGlobalNode,1, MPI::DOUBLE, MPI::SUM);
//
//
//         TKE  = 0.0 ;
//
// 	for(int i3= myGrid.nB3;  i3<= myGrid.nE3;  i3++)
//           for(int i2= myGrid.nB2;  i2<= myGrid.nE2;  i2++)
//             for(int i1= myGrid.nB1;  i1<= myGrid.nE1;  i1=i1+4)
//             {
//              getHydroMomentsFromCellWithForce(lbModel,myGrid,VECT_LENGTH,i1,i2,i3,rho,uX,uY,uZ,theta,F1,F2,F3);
//              for(int index=0;index<VECT_LENGTH;index++)
//                TKE += ( (uX[index]-uAvgGlobal)*(uX[index]-uAvgGlobal) + (uY[index]-vAvgGlobal)*(uY[index]-vAvgGlobal) + (uZ[index]-wAvgGlobal)*(uZ[index]-wAvgGlobal) );
//              }
//              MPI::COMM_WORLD.Allreduce(&TKE,&TKEGlobalCell,1, MPI::DOUBLE, MPI::SUM);
//
//
//         TKE =  (TKEGlobalNode +  TKEGlobalCell)/(2.0*nX1*nX2*nX3);
//   */
//     if(myRank==0)
//     std::cout<<step<<"   "<< uAvgGlobal/U0 <<"       "<<vAvgGlobal/U0<<"       "<<wAvgGlobal/U0<<"       "<<KE<<std::endl;
//   }
//

// Actual Perturbation
  template <int N,int numblock, typename dataType1>
  void perturb_parab(int *coordinates, lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int nP1, int nP2, int nP3, dataType1 normg,dataType1 U0,int VECT_LENGTH,dataType1 density, dataType1 uRef)
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

    int nX1 = myGrid.m1*nP1;
    int nX2 = myGrid.m2*nP2;
    int nX3 = myGrid.m3*nP3;
    char fileName1[150],fileName2[150],fileName3[150];

    double mag = 0.5 ;

 dataType1 p = myGrid.nB2 + 1.5;
 dataType1 q = nX2 + myGrid.nB2 -1 - 1.5 ;// myGrid.nE2 - 1.5;
 dataType1 pPlusq  = p+q;
 dataType1 qMinusp = q-p;
 dataType1 Y;



    int numPoints = nX1*nX2*nX3;


    for(int i3= myGrid.nB3;  i3<= myGrid.nE3;  i3++)
      for(int i2= myGrid.nB2;  i2<= myGrid.nE2;  i2++)
      {
//        int kGlobal = i3 + myGrid.m3*coordinates[2] ;
       int jGlobal = i2 + myGrid.m2*coordinates[1] ;
       for(int i1= myGrid.nB1;  i1<= myGrid.nE1;  i1=i1+4)
       {
//           int iGlobal = (i1 + index) + myGrid.m1*coordinates[0] ;

          copyFromNode(lbModel,myGrid,VECT_LENGTH,i1,i2,i3);
          getHydroMomentsFromNode(lbModel,myGrid,i1,i2,i3,rho,uX,uY,uZ,theta);

          Y = (dataType1)((i2-3.5 + (myGrid.m2*coordinates[1])) );
          Y = (dataType1)(Y - 0.5*pPlusq);
          Y = Y  / (dataType1)(0.5*qMinusp);
//           lbModel.uX    =1.5*U0*(1.0 - Y*Y); // U0;//


          for(int index=0;index<VECT_LENGTH;index++)
          {


            uX[index]    = U0*(1 - Y*Y) + (uX[index]/normg)*U0*mag ;
            uY[index]    = (uY[index]/normg)*U0*mag ;
            uZ[index]    = 1.1*(uZ[index]/normg)*U0*mag ;
            rho[index]   = density;
            theta[index] = lbModel.theta0;

               if(jGlobal == myGrid.nB2 || jGlobal == myGrid.m2*nP2+3)
               {
                 uX[index] = 0.0;
                 uY[index] = 0.0;
                 uZ[index] = 0.0 ;
               }

          }

          getFEq(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);
          copyToNode(lbModel,myGrid,VECT_LENGTH,i1,i2,i3);
          copyToCell(lbModel,myGrid,VECT_LENGTH,i1,i2,i3);
        }
      }
  }


// template void curl_noise_communicateAfterThis<4,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,1,double> &, int , int , int , double &,int *,int ,double , double , double );

// template void curl_noise_communiateBeforeThis<4,11,1,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, gridBCC3D<4,1,double> &, int , int , int , double &,int *,int ,double , double , double );

// template void averageVelocityAndTKE<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, int , int , int ,int *,int ,double ,double ,double ,int ,int ,int ,double);


template void curl_noise<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, int , int , int , double &,int *,int , double, double);

template void perturb_parab<4,11,double>(int *, lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, int , int , int , double ,double ,int, double, double);
