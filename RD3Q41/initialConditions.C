#include"initialConditions.h"



/*! \brief Initialise the Grid with Density 1 and  \f$ \theta_0 \f$ and Zero velocity
 */

  // Initialise the Grid with Density 1 and theta0 and Zero velocity
  template <int N,int numblock, typename dataType1>
  void initialConditions(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,dataType1 u_inlet,int FLUID,int nP2, int *coord,dataType1 rhoDSMC)
  {
    dataType1 *rho,*uX,*uY,*uZ,*theta;
    rho   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uX    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uY    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uZ    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    theta = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);

    unsigned short xsubi[3];
    xsubi[0] = 10;

    dataType1 halfY =  myGrid.m2*nP2*0.5;

    // std::cout<<"halfY: "<<"\t"<<halfY<<std::endl;

    for(int k=0; k<myGrid.n3;k++)
    {
      for(int j=0; j<myGrid.n2;j++)
      {
        for(int i=0; i<myGrid.n1;i=i+4)
        {
          for(int index=0;index<VECT_LENGTH;index++)
          {
            dataType1 yGlobal = (static_cast<dataType1>(j) + myGrid.m2*coord[1]) - halfY - myGrid.nB2; //  /static_cast<dataType1>(n2*nP2);

            if(marker(i+index,j,k,0,0,0) != FLUID)
            {
              rho[index]   = 3.0;
              uX[index]    = 0.0;
              uY[index]    = 0.0;
              uZ[index]    = 0.0;
              theta[index] = lbModel.theta0;
            }
            else
            {
              rho[index]   = rhoDSMC;
              uX[index]    = u_inlet*(1 -  (yGlobal)/(halfY)*(yGlobal)/(halfY));//0.0;//0.0;//0.0;//u_inlet;
              uY[index]    = 0.0;
              uZ[index]    = 0.0;
              theta[index] = lbModel.theta0;
            }
          }
          getFEq(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);
          copyToNode(lbModel,myGrid,VECT_LENGTH,i,j,k);
        }
      }
    }

    for(int k=0; k<myGrid.n3;k++)
    {
      for(int j=0; j<myGrid.n2;j++)
      {
        for(int i=0; i<myGrid.n1;i=i+4)
        {
          for(int index=0;index<VECT_LENGTH;index++)
          {
            dataType1 yGlobal = (static_cast<dataType1>(j) + myGrid.m2*coord[1]) + 0.5 - halfY - myGrid.nB2; //  /static_cast<dataType1>(n2*nP2);
            if(marker(i+index,j,k,0,1,0) != FLUID)
            {
              rho[index]   = 3.0;
              uX[index]    = 0.0;
              uY[index]    = 0.0;
              uZ[index]    = 0.0;
              theta[index] = lbModel.theta0;
            }
            else
            {
              rho[index]   = rhoDSMC;
              uX[index]    = u_inlet*(1 - (yGlobal)/(halfY)*(yGlobal)/(halfY)) ;//0.0;//u_inlet;//0.0;//u_inlet;//0.0;//u_inlet;
              uY[index]    = 0.0;
              uZ[index]    = 0.0;
              theta[index] = lbModel.theta0;
            }
          }
          getFEq(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);
          copyToCell(lbModel,myGrid,VECT_LENGTH,i,j,k);
        }
      }
    }


  }









  template <int N,int numblock, typename dataType1>
  void initializeKidaParallel(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 u_inlet, int *coord,int nP1,int nP2,int nP3)
  {
    dataType1 x,y,z;
    int n1 = myGrid.m1;
    int n2 = myGrid.m2;
    int n3 = myGrid.m3;
    dataType1 *rho,*uX,*uY,*uZ,*theta;
    rho   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uX    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uY    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uZ    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    theta = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);

    for(int k=myGrid.nB3; k<=myGrid.nE3;k++)
    {
      for(int j=myGrid.nB2; j<=myGrid.nE2;j++)
      {
        for(int i=myGrid.nB1; i<=myGrid.nE1;i=i+4)
        {
          for(int index=0;index<VECT_LENGTH;index++)
          {
            x = (static_cast<dataType1>(i +index)-0.5 + n1*coord[0] ) *2.0*M_PI/static_cast<dataType1>(n1*nP1);
            y = (static_cast<dataType1>(j)       -0.5 + n2*coord[1] ) *2.0*M_PI/static_cast<dataType1>(n2*nP2);
            z = (static_cast<dataType1>(k)       -0.5 + n3*coord[2] ) *2.0*M_PI/static_cast<dataType1>(n3*nP3);
            rho[index] = 1.0;
            uX[index] = u_inlet*sin(x)*(cos(3.0*y)*cos(z) - cos(y)*cos(3.0*z));
            uY[index] = u_inlet*sin(y)*(cos(3.0*z)*cos(x) - cos(z)*cos(3.0*x));
            uZ[index] = u_inlet*sin(z)*(cos(3.0*x)*cos(y) - cos(x)*cos(3.0*y));
            theta[index] = lbModel.theta0;
          }
          getFEqSIMD(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,10);
          copyToNode(lbModel,myGrid,VECT_LENGTH,i,j,k);
        }
      }
    }


    //Cells
    for(int k=myGrid.nB3; k<=myGrid.nE3;k++)
    {
      for(int j=myGrid.nB2; j<=myGrid.nE2;j++)
      {
        for(int i=myGrid.nB1; i<=myGrid.nE1;i=i+4)
        {
          for(int index=0;index<VECT_LENGTH;index++)
          {
            x = (static_cast<dataType1>(i+index)  + n1*coord[0] ) *2.0*M_PI/static_cast<dataType1>(n1*nP1);
            y = (static_cast<dataType1>(j)        + n2*coord[1] ) *2.0*M_PI/static_cast<dataType1>(n2*nP2);
            z = (static_cast<dataType1>(k)        + n3*coord[2] ) *2.0*M_PI/static_cast<dataType1>(n3*nP3);
            rho[index] = 1.0;
            uX[index] = u_inlet*sin(x)*(cos(3.0*y)*cos(z) - cos(y)*cos(3.0*z));
            uY[index] = u_inlet*sin(y)*(cos(3.0*z)*cos(x) - cos(z)*cos(3.0*x));
            uZ[index] = u_inlet*sin(z)*(cos(3.0*x)*cos(y) - cos(x)*cos(3.0*y));
            theta[index] = lbModel.theta0;
          }
          getFEqSIMD(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,10);
          copyToCell(lbModel,myGrid,VECT_LENGTH,i,j,k);
        }
      }
    }
  }



  // Explicit definitions
  template void initialConditions(lbmRD3Q41<double> &lbModel, gridBCC3D<4, 11, double> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,double u_inlet,int FLUID,int nP2, int *coord,double rhoDSMC);
  template void initializeKidaParallel(lbmRD3Q41<double> &lbModel, gridBCC3D<4, 11, double> &myGrid,int VECT_LENGTH,double u_inlet, int *coord,int nP1,int nP2,int nP3);
