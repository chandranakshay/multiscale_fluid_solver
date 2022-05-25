#include"diffuseBounceBack41_withForce.h"
  

  template <int N,int numblock, typename dataType1>
  void getRhoNodeSinglePoint(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,dataType1 &rho, int i1,int i2,int i3)
  {
    rho = 0.0;
    
    rho += myGrid(i1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    
    
    int index = myGrid.getIndex(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=      myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
  }   
  
  template <int N,int numblock, typename dataType1>
  void getRhoCellSinglePoint(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,dataType1 &rho, int i1,int i2,int i3)
  {
    rho = 0.0;
    
    rho += myGrid(i1,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    
    
    int index = myGrid.getIndex(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=      myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
    index = myGrid.getIndex(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],  0);   
    for(int dv=0;dv<4;dv++)
      rho +=     myGrid.value(index+dv) ;
    
  }   
  
  
  template <typename dataType1>
  void getFEqG12(lbmRD3Q41<dataType1> &lbModel,dataType1 rho,dataType1 uX,dataType1 uY,dataType1 uZ,dataType1 theta)
  {
    dataType1 delThetaBy2 = 0.5*(theta - lbModel.theta0)*lbModel.oneByTheta0;
    dataType1 delSqTheta0p125 =  0.5* delThetaBy2 *delThetaBy2;   
    
    dataType1 oneByTheta  = 1.0/(theta);
    uX      *= oneByTheta;
    uY      *= oneByTheta;
    uZ      *= oneByTheta;
    rho     *= 0.5;
    
    dataType1 twoMinusU2 = 2.0-(uX*uX + uY*uY + uZ*uZ)*theta;
    
    
    // CENTER
    dataType1  order1       =  delThetaBy2      *lbModel.yMinus3_CENTER;     
    dataType1  order2       =  delSqTheta0p125  *lbModel.yTenMinusYSqrPlus15_CENTER;      
    dataType1  f0V          =  rho*lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO]*(1.0 + order1+ order2); 
    lbModel.fTemp0[0][lbModel.CENTER_DV_ZERO_ZERO_ZERO]=  f0V*twoMinusU2;    
    
    // G1 and G2
    // ZZP1 and ZZM1
    order1       =  delThetaBy2     *lbModel.yMinus3_SC1;        
    order2       =  delSqTheta0p125 *lbModel.yTenMinusYSqrPlus15_SC1;         
    f0V          =  rho*lbModel.wt[lbModel.DV_ZERO_ZERO_P1]*(1.0 + order1+ order2);        
    dataType1  uDotC        =  uZ;
    lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC)); 
    lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
    
    order1       =  delThetaBy2     *lbModel.yMinus3_SC2;        
    order2       =  delSqTheta0p125 *lbModel.yTenMinusYSqrPlus15_SC2;         
    f0V          =  rho*lbModel.wt[lbModel.DV_ZERO_ZERO_P2]*(1.0 + order1+ order2);   
    // ZZP2 and ZZM2
    uDotC        =  2.0*uZ;
    lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC)); 
    lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC)); 
    // ZP2Z and ZM2Z
    uDotC        =  2.0*uY;
    lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC)); 
    lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));         
    // P2ZZ and M2ZZ
    uDotC        =  2.0*uX;
    lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC)); 
    lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC)); 
    
  }   
  
  template <typename dataType1>
  void getFEqG34(lbmRD3Q41<dataType1> &lbModel,dataType1 rho,dataType1 uX,dataType1 uY,dataType1 uZ,dataType1 theta)
  {
    dataType1 delThetaBy2 ;
    dataType1 f0V(0.0),uDotC(0.0),twoMinusU2(0.0);
    dataType1 order1(0.0),order2(0.0);
    dataType1 oneByTheta ,  delSqTheta0p125 ;                                        
    
    oneByTheta  = 1.0/theta ;                                                             
    delThetaBy2 = 0.5*(theta - lbModel.theta0)*lbModel.oneByTheta0;                       
    delSqTheta0p125 =  0.5* delThetaBy2*delThetaBy2;                                              
    
    rho  *= 0.5;                                                                          
    dataType1 uX1 = uX*oneByTheta;                                           
    dataType1 uY1 = uY*oneByTheta;
    dataType1 uZ1 = uZ*oneByTheta;                                                         
    
    //         u2  =  (uX1*uX1+uY1*uY1+uZ1*uZ1)*oneByTheta;                                               
    //         twoMinusU2 = 2.0-u2;
    
    twoMinusU2 = 2.0-(uX1*uX1 + uY1*uY1 + uZ1*uZ1)*theta;
    
    // G3 and G4
    order1    =  delThetaBy2    *lbModel.yMinus3_FCC;        
    order2    =  delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_FCC;         
    f0V       =  rho*lbModel.wt[lbModel.DV_ZERO_P1_P1]*(1.0 + order1+ order2);        
    //ZPP and ZMM         
    uDotC     =  (uY1 + uZ1);        
    lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));        
    lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
    //ZMP and ZPM                                                                                            
    uDotC     =  (uZ1 - uY1);                                                                       
    lbModel.fTemp3[0][lbModel.G3_DV_ZERO_M1_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));        
    lbModel.fTemp4[0][lbModel.G4_DV_ZERO_P1_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
    //PZP and MZM                                                                                            
    uDotC     =  (uX1 + uZ1);                                                                       
    lbModel.fTemp3[0][lbModel.G3_DV_P1_ZERO_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));        
    lbModel.fTemp4[0][lbModel.G4_DV_M1_ZERO_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
    //MZP and PZM                                                                                            
    uDotC     =  (uZ1 - uX1);                                                                   
    lbModel.fTemp3[0][lbModel.G3_DV_M1_ZERO_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));    
    lbModel.fTemp4[0][lbModel.G4_DV_P1_ZERO_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC)); 
  }   
  
  template <typename dataType1>
  void getFEqG56(lbmRD3Q41<dataType1> &lbModel,dataType1 rho,dataType1 uX,dataType1 uY,dataType1 uZ,dataType1 theta)
  {
    dataType1 delThetaBy2 ;
    dataType1 f0V(0.0),uDotC(0.0),twoMinusU2(0.0);
    dataType1 order1(0.0),order2(0.0);
    dataType1 oneByTheta ,  delSqTheta0p125 ;                                        
    
    oneByTheta  = 1.0/theta ;                                                             
    delThetaBy2 = 0.5*(theta - lbModel.theta0)*lbModel.oneByTheta0;                       
    delSqTheta0p125 =  0.5* delThetaBy2*delThetaBy2;                                              
    
    rho  *= 0.5;                                                                          
    dataType1 uX1 = uX*oneByTheta;                                           
    dataType1 uY1 = uY*oneByTheta;
    dataType1 uZ1 = uZ*oneByTheta;                                                         
    
    //         u2  =  (uX1*uX1+uY1*uY1+uZ1*uZ1)*oneByTheta;                                               
    //         twoMinusU2 = 2.0-u2;
    twoMinusU2 = 2.0-(uX1*uX1 + uY1*uY1 + uZ1*uZ1)*theta;
    
    order1    =  delThetaBy2    *lbModel.yMinus3_FCC;
    order2    =  delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_FCC; 
    f0V       =  rho*lbModel.wt[lbModel.DV_P1_P1_ZERO]*(1.0 + order1+ order2);
    //PPZ and MMZ
    uDotC     =  uX1 + uY1;
    lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
    lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
    //MPZ and PMZ                                                                               
    uDotC     =  (uY1 - uX1);                                                             
    lbModel.fTemp5[0][lbModel.G5_DV_M1_P1_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
    lbModel.fTemp6[0][lbModel.G6_DV_P1_M1_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
    //ZPZ and ZMZ
    order1    =  delThetaBy2    *lbModel.yMinus3_SC1;
    order2    =  delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_SC1; 
    f0V       =  rho*lbModel.wt[lbModel.DV_ZERO_P1_ZERO]*(1.0 + order1+ order2);    
    uDotC     =  uY1;
    lbModel.fTemp5[0][lbModel.G5_DV_ZERO_P1_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
    lbModel.fTemp6[0][lbModel.G6_DV_ZERO_M1_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
    //PZZ aqnd MZZ                                                                                  
    uDotC     =   uX1;                                                                          
    lbModel.fTemp5[0][lbModel.G5_DV_P1_ZERO_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
    lbModel.fTemp6[0][lbModel.G6_DV_M1_ZERO_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));                      
  }                                                                                                            
  
  template <typename dataType1>
  void getFEqG78910(lbmRD3Q41<dataType1> &lbModel,dataType1 rho,dataType1 uX,dataType1 uY,dataType1 uZ,dataType1 theta)
  {
    dataType1 delThetaBy2 ;
    dataType1 f0V(0.0),f0V1(0.0),uDotC(0.0),twoMinusU2(0.0),uDotCSq(0.0);
    dataType1 order1(0.0),order2(0.0);
    dataType1 oneByTheta ,  delSqTheta0p125 ;                                        
    
    delThetaBy2 = 0.5*(theta - lbModel.theta0)*lbModel.oneByTheta0;                       
    delSqTheta0p125 =  0.5* delThetaBy2*delThetaBy2;                                              
    oneByTheta  = 1.0/theta ;                                                             
    
    rho  *= 0.5;                                                                          
    dataType1 uX1 = uX*oneByTheta;                                           
    dataType1 uY1 = uY*oneByTheta;
    dataType1 uZ1 = uZ*oneByTheta;                                                         
    
    twoMinusU2 = 2.0-(uX1*uX1 + uY1*uY1 + uZ1*uZ1)*theta;
    
    order1    =  1.0 + delThetaBy2*lbModel.yMinus3_BCC +  delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_BCC; 
    f0V       =  rho*lbModel.wt[lbModel.DV_P_P_P]*(order1);
    order2    =  1.0 + delThetaBy2*lbModel.yMinus3_BCC1 + delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_BCC1;         
    f0V1      =  rho*lbModel.wt[lbModel.DV_P1_P1_P1]*(order2);        
    //PPP1 and MMM1        
    uDotC     = (uX1 + uY1 + uZ1);      
    uDotCSq   = uDotC*uDotC;        
    lbModel.fTemp9 [0][lbModel.G9_DV_P1_P1_P1]   = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);
    lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1]  = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);   
    //PPP and MMM
    lbModel.fTemp7[0][lbModel.G7_DV_P_P_P]   = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);
    lbModel.fTemp8[0][lbModel.G8_DV_M_M_M]   = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);
    
    //MPP1 and PMM1                                                                             
    uDotC     = (uY1 - uX1 + uZ1); 
    uDotCSq   = uDotC*uDotC;
    lbModel.fTemp9 [0][lbModel.G9_DV_M1_P1_P1]  = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);
    lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);  
    //MPP and PMM                                                                           
    lbModel.fTemp7[0][lbModel.G7_DV_M_P_P] = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);     
    lbModel.fTemp8[0][lbModel.G8_DV_P_M_M] = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);       
    
    //MMP1 and PPM 1                                                                                     
    uDotC     = (uZ1 - uX1 - uY1);           
    uDotCSq   = uDotC*uDotC;       
    lbModel.fTemp9 [0][lbModel.G9_DV_M1_M1_P1] = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);       
    lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);        
    //MMP and PPM                                                                                  
    lbModel.fTemp7[0][lbModel.G7_DV_M_M_P] = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);       
    lbModel.fTemp8[0][lbModel.G8_DV_P_P_M] = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);       
    
    //PMP and MPM                                                                                      
    uDotC     = (uX1 - uY1 + uZ1);        
    uDotCSq   = uDotC*uDotC;       
    lbModel.fTemp9 [0][lbModel.G9_DV_P1_M1_P1]  = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);       
    lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);         
    //PMP and MPM                                                                                  
    lbModel.fTemp7[0][lbModel.G7_DV_P_M_P] = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);       
    lbModel.fTemp8[0][lbModel.G8_DV_M_P_M] = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq); 
  }
  
  
  template <int N,int numblock, typename dataType1>
  void getRhoForDiffuseBounceBackNodeG12(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID)            
  {
    dataType1 rho(0.0);
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                  
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) != FLUID)
          {
            getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2,i3);
            rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0) = rho; 
            
            // ZZP1
            if(marker(i1,i2,i3-1,0,nodeTYPE::NODE,0)==FLUID) 
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2,i3-1);
              rhoGrid(i1,i2,i3-1,0,nodeTYPE::NODE,0) = rho; 
            } 
            //ZZP2
            if((marker(i1,i2,i3-1,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2,i3-2,0,nodeTYPE::NODE,0)==FLUID)) 
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2,i3-2);
              rhoGrid(i1,i2,i3-2,0,nodeTYPE::NODE,0) = rho; 
            }     
            //ZZP2 Correction
            if( ((marker(i1,i2,i3-1,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeTYPE::NODE,0)==FLUID)) || ( (marker(i1,i2,i3-1,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeTYPE::NODE,0)!=FLUID) ) )
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2,i3-1);
              rhoGrid(i1,i2,i3-1,0,nodeTYPE::NODE,0) = rho; 
            }    
            //ZP2Z
            if( (marker(i1,i2-1,i3,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2-2,i3,0,nodeTYPE::NODE,0)==FLUID) )
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2-2,i3);
              rhoGrid(i1,i2-2,i3,0,nodeTYPE::NODE,0) = rho; 
            }        
            //ZP2Z Corrrection
            if( ((marker(i1,i2-1,i3,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2+1,i3,0,nodeTYPE::NODE,0)==FLUID)) || ( (marker(i1,i2-1,i3,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2+1,i3,0,nodeTYPE::NODE,0)!=FLUID) ) )
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2-1,i3);
              rhoGrid(i1,i2-1,i3,0,nodeTYPE::NODE,0) = rho; 
            }    
            //P2ZZ
            if( (marker(i1-1,i2,i3,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1-2,i2,i3,0,nodeTYPE::NODE,0)==FLUID) )
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1-2,i2,i3);
              rhoGrid(i1-2,i2,i3,0,nodeTYPE::NODE,0) = rho; 
            }     
            //P2ZZ Correction
            if( ((marker(i1-1,i2,i3,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1+1,i2,i3,0,nodeTYPE::NODE,0)==FLUID)) || ( (marker(i1-1,i2,i3,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1+1,i2,i3,0,nodeTYPE::NODE,0)!=FLUID) ) )
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1-1,i2,i3);
              rhoGrid(i1-1,i2,i3,0,nodeTYPE::NODE,0) = rho; 
            }       
            
            
            // ZZM1
            if( marker(i1,i2,i3+1,0,nodeTYPE::NODE,0)==FLUID )
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2,i3+1);
              rhoGrid(i1,i2,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }      
            //ZZM2
            if( (marker(i1,i2,i3+1,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2,i3+2,0,nodeTYPE::NODE,0)==FLUID) )
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2,i3+2);
              rhoGrid(i1,i2,i3+2,0,nodeTYPE::NODE,0) = rho; 
            }   
            //ZZM2 correction
            if( (marker(i1,i2,i3-1,0,nodeTYPE::NODE,0)!=FLUID ) && (marker(i1,i2,i3+1,0,nodeTYPE::NODE,0)==FLUID) )
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2,i3+1);
              rhoGrid(i1,i2,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }
            
            //ZM2Z
            if( (marker(i1,i2+1,i3,0,nodeTYPE::NODE,0)==FLUID) && (marker(i1,i2+2,i3,0,nodeTYPE::NODE,0)==FLUID) )
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2+2,i3);
              rhoGrid(i1,i2+2,i3,0,nodeTYPE::NODE,0) = rho; 
            }
            //ZM2Z correction
            if( (marker(i1,i2-1,i3,0,nodeTYPE::NODE,0)!=FLUID) && (marker(i1,i2+1,i3,0,nodeTYPE::NODE,0)==FLUID) )
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2+1,i3);
              rhoGrid(i1,i2+1,i3,0,nodeTYPE::NODE,0) = rho; 
            }
            //M2ZZ 
            if( (marker(i1+1,i2,i3,0,nodeTYPE::NODE,0)==FLUID) && (marker(i1+2,i2,i3,0,nodeTYPE::NODE,0)==FLUID) )
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+2,i2,i3);
              rhoGrid(i1+2,i2,i3,0,nodeTYPE::NODE,0) = rho; 
            }
            //M2ZZ Correction
            if( (marker(i1-1,i2,i3,0,nodeTYPE::NODE,0)!=FLUID) && (marker(i1+1,i2,i3,0,nodeTYPE::NODE,0)==FLUID) )
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2,i3);
              rhoGrid(i1+1,i2,i3,0,nodeTYPE::NODE,0) = rho; 
            }
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void getRhoForDiffuseBounceBackNodeG34(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID)            
  {
    dataType1 rho(0.0);
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                  
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) != FLUID)
          {
            getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3  );
            rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0) = rho; 
            //ZPP
            if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2-1,i3-1);
              rhoGrid(i1  ,i2-1,i3-1,0,nodeTYPE::NODE,0) = rho; 
            }
            //ZMP 
            if(marker(i1  ,i2+1,i3-1,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2+1,i3-1);
              rhoGrid(i1  ,i2+1,i3-1,0,nodeTYPE::NODE,0) = rho; 
            }
            //PZP
            if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1-1,i2  ,i3-1);
              rhoGrid(i1-1,i2  ,i3-1,0,nodeTYPE::NODE,0) = rho; 
            }                           
            //MZP
            if(marker(i1+1,i2  ,i3-1,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2  ,i3-1);
              rhoGrid(i1+1,i2  ,i3-1,0,nodeTYPE::NODE,0) = rho; 
            }   
            //ZMM
            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)                                                                                    
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2+1,i3+1);
              rhoGrid(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }    
            //ZPM
            if(marker(i1  ,i2-1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2-1,i3+1);
              rhoGrid(i1  ,i2-1,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }
            //MZM
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2  ,i3+1);
              rhoGrid(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }    
            //PZM
            if(marker(i1-1,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1-1,i2  ,i3+1);
              rhoGrid(i1-1,i2  ,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }                                     
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void getRhoForDiffuseBounceBackNodeG56(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID)            
  {
    dataType1 rho(0.0);
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) != FLUID)
          {  
            getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3  );
            rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0) = rho; 
            //PPZ
            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1-1,i2-1,i3  );
              rhoGrid(i1-1,i2-1,i3  ,0,nodeTYPE::NODE,0) = rho; 
            }    
            //ZPZ
            if(marker(i1+1,i2-1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2-1,i3  );
              rhoGrid(i1+1,i2-1,i3  ,0,nodeTYPE::NODE,0) = rho; 
            }    
            //ZPZ
            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2-1,i3  );
              rhoGrid(i1  ,i2-1,i3  ,0,nodeTYPE::NODE,0) = rho; 
            }  
            //PZZ
            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1-1,i2  ,i3  );
              rhoGrid(i1-1,i2  ,i3  ,0,nodeTYPE::NODE,0) = rho; 
            }    
            //MMZ
            if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2+1,i3  );
              rhoGrid(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0) = rho; 
            }    
            //PZM
            if(marker(i1-1,i2+1,i3  ,0,nodeTYPE::NODE,0) ==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1-1,i2+1,i3  );
              rhoGrid(i1-1,i2+1,i3  ,0,nodeTYPE::NODE,0) = rho; 
            }    
            //ZMZ
            if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2+1,i3  );
              rhoGrid(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0) = rho; 
            }    
            //MZZ
            if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2  ,i3  );
              rhoGrid(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0) = rho; 
            }              
          }   
        }
      }
    }
  }         
  
  template <int N,int numblock, typename dataType1>
  void getRhoForDiffuseBounceBackNodeG910(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID)            
  {
    dataType1 rho(0.0);
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) != FLUID)
          {  
            getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3  );
            rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0) = rho; 
            //P1P1P1
            if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1-1,i2-1,i3-1);
              rhoGrid(i1-1,i2-1,i3-1,0,nodeTYPE::NODE,0) = rho; 
            }    
            //M1P1P1
            if(marker(i1+1,i2-1,i3-1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2-1,i3-1);
              rhoGrid(i1+1,i2-1,i3-1,0,nodeTYPE::NODE,0) = rho; 
            }    
            //M1M1P1
            if(marker(i1+1,i2+1,i3-1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2+1,i3-1);
              rhoGrid(i1+1,i2+1,i3-1,0,nodeTYPE::NODE,0) = rho; 
            }    
            //P1M1P1
            if(marker(i1-1,i2+1,i3-1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1-1,i2+1,i3-1);
              rhoGrid(i1-1,i2+1,i3-1,0,nodeTYPE::NODE,0) = rho; 
            }    
            //M1M1M1
            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2+1,i3+1);
              rhoGrid(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }    
            //P1M1M1
            if(marker(i1-1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1-1,i2+1,i3+1);
              rhoGrid(i1-1,i2+1,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }    
            //P1P1M1
            if(marker(i1-1,i2-1,i3+1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1-1,i2-1,i3+1);
              rhoGrid(i1-1,i2-1,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }    
            //M1P1M1
            if(marker(i1+1,i2-1,i3+1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2-1,i3+1);
              rhoGrid(i1+1,i2-1,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void getRhoForDiffuseBounceBackNodeCellG78(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID)            
  {
    dataType1 rho(0.0);
    ///////////////////
    // Node --> Cell //
    ///////////////////        
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) != FLUID)
          {  
            getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3  );
            rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0) = rho; 
            //PPP
            if(marker(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3  );
              rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0) = rho;                             
              
              if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
                getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2+1,i3+1);
                rhoGrid(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0) = rho;                              
              }
            }
            //MPP
            if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2  ,i3  );
              rhoGrid(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0) = rho; 
              
              if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {   
                getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2+1,i3+1);
                rhoGrid(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0) = rho; 
              }
            }
            
            //MMP
            if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2+1,i3  );
              rhoGrid(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0) = rho; 
              
              if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
                getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3+1);
                rhoGrid(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0) = rho; 
              }
            }            
            //PMP
            if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2+1,i3  );
              rhoGrid(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0) = rho; 
              
              if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
                getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2  ,i3+1);
                rhoGrid(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0) = rho; 
              }    
            }            
            //MMM
            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2+1,i3+1);
              rhoGrid(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }                    
            //PMM
            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2+1,i3+1);
              rhoGrid(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }    
            //PPM
            if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3+1);
              rhoGrid(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }    
            //MPM
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              getRhoNodeSinglePoint(lbModel,myGrid,rho,i1+1,i2  ,i3+1);
              rhoGrid(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0) = rho; 
            }    
          }
        }
      }
    }
    
    ////////////////////
    // Cell --> Node  //
    ////////////////////             
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0)!=FLUID)
          {    
            getRhoNodeSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3  );
            rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0) = rho; 
            
            //PPP
            if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2-1,i3-1);
              rhoGrid(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0) = rho; 
              
              if(marker(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
              {
                getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3  );
                rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0) = rho; 
              }
            }
            //MPP
            if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {              
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2-1,i3-1);
              rhoGrid(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0) = rho; 
              
              if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
              {             
                getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2  ,i3  );
                rhoGrid(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0) = rho; 
              }             
            }             
            //MMP
            if(marker(i1  ,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {              
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3-1);
              rhoGrid(i1  ,i2  ,i3-1,0,nodeTYPE::CELL,0) = rho; 
              
              if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
              {             
                getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2-1,i3  );
                rhoGrid(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0) = rho; 
              }                             
            }     
            //PMP
            if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)                                                                                    
            {      
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2  ,i3-1);
              rhoGrid(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0) = rho; 
              
              if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
              {    
                getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2-1,i3  );
                rhoGrid(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0) = rho; 
              }      
            }       
            //MMM
            if(marker(i1  ,i2   ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {    
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2   ,i3  );
              rhoGrid(i1  ,i2   ,i3  ,0,nodeTYPE::CELL,0) = rho; 
            }                 
            //PMM
            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {    
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2  ,i3  );
              rhoGrid(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0) = rho; 
            }        
            //PPM
            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {    
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2-1,i3  );
              rhoGrid(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0) = rho; 
            }    
            //MPM
            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {    
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2-1,i3  );
              rhoGrid(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0) = rho; 
            }                
          }  
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void getRhoForDiffuseBounceBackCellG12(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID)            
  {
    dataType1 rho(0.0);
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) != FLUID)
          {
            getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2,i3);
            rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0) = rho; 
            
            // ZZP1
            if(marker(i1,i2,i3-1,0,nodeTYPE::CELL,0)==FLUID) 
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2,i3-1);
              rhoGrid(i1,i2,i3-1,0,nodeTYPE::CELL,0) = rho; 
            } 
            //ZZP2
            if((marker(i1,i2,i3-1,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2,i3-2,0,nodeTYPE::CELL,0)==FLUID)) 
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2,i3-2);
              rhoGrid(i1,i2,i3-2,0,nodeTYPE::CELL,0) = rho; 
            }     
            //ZZP2 Correction
            if( ((marker(i1,i2,i3-1,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeTYPE::CELL,0)==FLUID)) || ( (marker(i1,i2,i3-1,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeTYPE::CELL,0)!=FLUID) ) )
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2,i3-1);
              rhoGrid(i1,i2,i3-1,0,nodeTYPE::CELL,0) = rho; 
            }    
            //ZP2Z
            if( (marker(i1,i2-1,i3,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2-2,i3,0,nodeTYPE::CELL,0)==FLUID) )
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2-2,i3);
              rhoGrid(i1,i2-2,i3,0,nodeTYPE::CELL,0) = rho; 
            }        
            //ZP2Z Corrrection
            if( ((marker(i1,i2-1,i3,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2+1,i3,0,nodeTYPE::CELL,0)==FLUID)) || ( (marker(i1,i2-1,i3,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2+1,i3,0,nodeTYPE::CELL,0)!=FLUID) ) )
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2-1,i3);
              rhoGrid(i1,i2-1,i3,0,nodeTYPE::CELL,0) = rho; 
            }    
            //P2ZZ
            if( (marker(i1-1,i2,i3,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1-2,i2,i3,0,nodeTYPE::CELL,0)==FLUID) )
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-2,i2,i3);
              rhoGrid(i1-2,i2,i3,0,nodeTYPE::CELL,0) = rho; 
            }     
            //P2ZZ Correction
            if( ((marker(i1-1,i2,i3,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1+1,i2,i3,0,nodeTYPE::CELL,0)==FLUID)) || ( (marker(i1-1,i2,i3,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1+1,i2,i3,0,nodeTYPE::CELL,0)!=FLUID) ) )
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2,i3);
              rhoGrid(i1-1,i2,i3,0,nodeTYPE::CELL,0) = rho; 
            }       
            
            
            // ZZM1
            if( marker(i1,i2,i3+1,0,nodeTYPE::CELL,0)==FLUID )
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2,i3+1);
              rhoGrid(i1,i2,i3+1,0,nodeTYPE::CELL,0) = rho; 
            }      
            //ZZM2
            if( (marker(i1,i2,i3+1,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2,i3+2,0,nodeTYPE::CELL,0)==FLUID) )
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2,i3+2);
              rhoGrid(i1,i2,i3+2,0,nodeTYPE::CELL,0) = rho; 
            }   
            //ZZM2 correction
            if( (marker(i1,i2,i3-1,0,nodeTYPE::CELL,0)!=FLUID ) && (marker(i1,i2,i3+1,0,nodeTYPE::CELL,0)==FLUID) )
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2,i3+1);
              rhoGrid(i1,i2,i3+1,0,nodeTYPE::CELL,0) = rho; 
            }
            
            //ZM2Z
            if( (marker(i1,i2+1,i3,0,nodeTYPE::CELL,0)==FLUID) && (marker(i1,i2+2,i3,0,nodeTYPE::CELL,0)==FLUID) )
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2+2,i3);
              rhoGrid(i1,i2+2,i3,0,nodeTYPE::CELL,0) = rho; 
            }
            //ZM2Z correction
            if( (marker(i1,i2-1,i3,0,nodeTYPE::CELL,0)!=FLUID) && (marker(i1,i2+1,i3,0,nodeTYPE::CELL,0)==FLUID) )
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2+1,i3);
              rhoGrid(i1,i2+1,i3,0,nodeTYPE::CELL,0) = rho; 
            }
            //M2ZZ 
            if( (marker(i1+1,i2,i3,0,nodeTYPE::CELL,0)==FLUID) && (marker(i1+2,i2,i3,0,nodeTYPE::CELL,0)==FLUID) )
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1+2,i2,i3);
              rhoGrid(i1+2,i2,i3,0,nodeTYPE::CELL,0) = rho; 
            }
            //M2ZZ Correction
            if( (marker(i1-1,i2,i3,0,nodeTYPE::CELL,0)!=FLUID) && (marker(i1+1,i2,i3,0,nodeTYPE::CELL,0)==FLUID) )
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1+1,i2,i3);
              rhoGrid(i1+1,i2,i3,0,nodeTYPE::CELL,0) = rho; 
            }
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void getRhoForDiffuseBounceBackCellG34(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID)            
  {
    dataType1 rho(0.0);
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) != FLUID)
          {
            getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3  );
            rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0) = rho; 
            
            // ZPP   
            if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2-1,i3-1);
              rhoGrid(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0) = rho; 
            }
            // ZMP   
            if(marker(i1  ,i2+1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2+1,i3-1);
              rhoGrid(i1  ,i2+1,i3-1,0,nodeTYPE::CELL,0) = rho; 
            }
            // PZP   
            if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2  ,i3-1);
              rhoGrid(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0) = rho; 
            }                           
            // MZP   
            if(marker(i1+1,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1+1,i2  ,i3-1);
              rhoGrid(i1+1,i2  ,i3-1,0,nodeTYPE::CELL,0) = rho; 
            }             
            // ZMM   
            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::CELL,0)==FLUID)                                                                                    
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2+1,i3+1);
              rhoGrid(i1  ,i2+1,i3+1,0,nodeTYPE::CELL,0) = rho; 
            }                           
            // ZPM   
            if(marker(i1  ,i2-1,i3+1,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2-1,i3+1);
              rhoGrid(i1  ,i2-1,i3+1,0,nodeTYPE::CELL,0) = rho; 
            }         
            // MZM   
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1+1,i2  ,i3+1);
              rhoGrid(i1+1,i2  ,i3+1,0,nodeTYPE::CELL,0) = rho; 
            }         
            // PZM   
            if(marker(i1-1,i2  ,i3+1,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2  ,i3+1);
              rhoGrid(i1-1,i2  ,i3+1,0,nodeTYPE::CELL,0) = rho; 
            }                                     
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void getRhoForDiffuseBounceBackCellG56(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID)            
  {
    dataType1 rho(0.0);
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) != FLUID)
          {  
            getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3  );
            rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0) = rho; 
            
            //PPZ
            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2-1,i3  );
              rhoGrid(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0) = rho; 
            }         
            //MPZ
            if(marker(i1+1,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1+1,i2-1,i3  );
              rhoGrid(i1+1,i2-1,i3  ,0,nodeTYPE::CELL,0) = rho; 
            }             
            //ZPZ
            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2-1,i3  );
              rhoGrid(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0) = rho; 
            }              
            //PZZ
            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2  ,i3  );
              rhoGrid(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0) = rho; 
            }              
            //MMZ
            if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1+1,i2+1,i3  );
              rhoGrid(i1+1,i2+1,i3  ,0,nodeTYPE::CELL,0) = rho; 
            }              
            //PMZ
            if(marker(i1-1,i2+1,i3  ,0,nodeTYPE::CELL,0) ==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2+1,i3  );
              rhoGrid(i1-1,i2+1,i3  ,0,nodeTYPE::CELL,0) = rho; 
            }              
            //ZMZ
            if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2+1,i3  );
              rhoGrid(i1  ,i2+1,i3  ,0,nodeTYPE::CELL,0) = rho; 
            }              
            //MZZ
            if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1+1,i2  ,i3  );
              rhoGrid(i1+1,i2  ,i3  ,0,nodeTYPE::CELL,0) = rho; 
            }              
          }
        }
      }
    }
  }         
  
  template <int N,int numblock, typename dataType1>
  void getRhoForDiffuseBounceBackCellG910(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID)            
  {  
    dataType1 rho(0.0);
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) != FLUID)
          {  
            getRhoCellSinglePoint(lbModel,myGrid,rho,i1  ,i2  ,i3  );
            rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0) = rho; 
            
            //PPP
            if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2-1,i3-1);
              rhoGrid(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0) = rho; 
            }                  
            //MPP
            if(marker(i1+1,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1+1,i2-1,i3-1);
              rhoGrid(i1+1,i2-1,i3-1,0,nodeTYPE::CELL,0) = rho; 
            }                             
            //MMP
            if(marker(i1+1,i2+1,i3-1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1+1,i2+1,i3-1);
              rhoGrid(i1+1,i2+1,i3-1,0,nodeTYPE::CELL,0) = rho; 
            }                            
            //PMP
            if(marker(i1-1,i2+1,i3-1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2+1,i3-1);
              rhoGrid(i1-1,i2+1,i3-1,0,nodeTYPE::CELL,0) = rho; 
            }                             
            //MMM
            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1+1,i2+1,i3+1);
              rhoGrid(i1+1,i2+1,i3+1,0,nodeTYPE::CELL,0) = rho; 
            } 
            //PMM
            if(marker(i1-1,i2+1,i3+1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2+1,i3+1);
              rhoGrid(i1-1,i2+1,i3+1,0,nodeTYPE::CELL,0) = rho; 
            }                             
            //PPM
            if(marker(i1-1,i2-1,i3+1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1-1,i2-1,i3+1);
              rhoGrid(i1-1,i2-1,i3+1,0,nodeTYPE::CELL,0) = rho; 
            }                             
            //MPM
            if(marker(i1+1,i2-1,i3+1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1+1,i2-1,i3+1);
              rhoGrid(i1+1,i2-1,i3+1,0,nodeTYPE::CELL,0) = rho; 
            }
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackNodeG12withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                          
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) != FLUID)
          {
            
            //ZZP1
            if(marker(i1,i2,i3-1,0,nodeTYPE::NODE,0)==FLUID) 
            {
              
              rho = rhoGrid(i1,i2,i3-1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M1) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];
              ForceZ += myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);                            
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P1) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1];
            } 
            //ZZP2
            if((marker(i1,i2,i3-1,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2,i3-2,0,nodeTYPE::NODE,0)==FLUID)) 
            {
              rho = rhoGrid(i1,i2,i3-2,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3-2,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
              ForceZ += 2.0*myGrid(i1  ,i2  ,i3-2,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2);
              
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
            }        
            //ZZP2 - Correction
            if(    ( (marker(i1,i2,i3-1,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeTYPE::NODE,0)==FLUID) ) 
              || ( (marker(i1,i2,i3-1,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeTYPE::NODE,0)!=FLUID) ) )
            {
              rho = rhoGrid(i1,i2,i3-1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
              ForceZ += 2.0*myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2);
              
              rho = rhoGrid(i1,i2,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);                               
              myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
              ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);                               
              myGrid(i1  ,i2  ,i3,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
              myGrid(i1  ,i2  ,i3,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
            }    
            //ZP2Z
            if( (marker(i1,i2-1,i3,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2-2,i3,0,nodeTYPE::NODE,0)==FLUID) )
            {
              rho = rhoGrid(i1,i2-2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-2,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
              ForceY += 2.0*myGrid(i1  ,i2-2,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);                               
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
            }        
            //ZP2Z - Correction
            if(    ((marker(i1,i2-1,i3,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2+1,i3,0,nodeTYPE::NODE,0)==FLUID)) 
              || ((marker(i1,i2-1,i3,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2+1,i3,0,nodeTYPE::NODE,0)!=FLUID)) )
            {
              rho = rhoGrid(i1,i2-1,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-1,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
              ForceY += 2.0*myGrid(i1  ,i2-1,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO);
              
              rho = rhoGrid(i1,i2+1,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+1,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
              ForceY -= 2.0*myGrid(i1  ,i2+1,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);                               
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
            }    
            //P2ZZ
            if( (marker(i1-1,i2,i3,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1-2,i2,i3,0,nodeTYPE::NODE,0)==FLUID) )
            {
              rho = rhoGrid(i1-2,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-2,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
              ForceX += 2.0*myGrid(i1-2,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
            }        
            //P2ZZ - Correction
            if( ((marker(i1-1,i2,i3,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1+1,i2,i3,0,nodeTYPE::NODE,0)==FLUID)) || ( (marker(i1-1,i2,i3,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1+1,i2,i3,0,nodeTYPE::NODE,0)!=FLUID) ) )
            {
              rho = rhoGrid(i1-1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
              ForceX += 2.0*myGrid(i1-1,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO);
              
              rho = rhoGrid(i1+1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
              ForceX -= 2.0*myGrid(i1+1,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);                               
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
            }       
            
            
            
            
            
            
            
            //ZZM1
            if( marker(i1,i2,i3+1,0,nodeTYPE::NODE,0)==FLUID )
            {
              rho = rhoGrid(i1,i2,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P1) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1];
              ForceZ -= myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M1) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];
            }       
            //ZZM2
            if( (marker(i1,i2,i3+1,0,nodeTYPE::NODE,0)==FLUID ) && (marker(i1,i2,i3+2,0,nodeTYPE::NODE,0)==FLUID) )
            {
              rho = rhoGrid(i1,i2,i3+2,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
              ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
            }   
            //                         ZZM2 - Correction
            if( (marker(i1,i2,i3-1,0,nodeTYPE::NODE,0)!=FLUID ) && (marker(i1,i2,i3+1,0,nodeTYPE::NODE,0)==FLUID) )
            {
              rho = rhoGrid(i1,i2,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
              ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2);
              
              rho = rhoGrid(i1,i2,i3-1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
              ForceZ += 2.0*myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2);                            
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
            }
            
            //ZM2Z
            if( (marker(i1,i2+1,i3,0,nodeTYPE::NODE,0)==FLUID) && (marker(i1,i2+2,i3,0,nodeTYPE::NODE,0)==FLUID) )
            {
              rho = rhoGrid(i1,i2+2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+2,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
              ForceY -= 2.0*myGrid(i1  ,i2+2,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
            }
            //ZM2Z - Correction
            if( (marker(i1,i2-1,i3,0,nodeTYPE::NODE,0)!=FLUID) && (marker(i1,i2+1,i3,0,nodeTYPE::NODE,0)==FLUID) )
            {
              rho = rhoGrid(i1,i2+1,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+1,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
              ForceY -= 2.0*myGrid(i1  ,i2+1,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO);
              
              rho = rhoGrid(i1,i2-1,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-1,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
              ForceY += 2.0*myGrid(i1  ,i2-1,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
            }
            //M2ZZ
            if( (marker(i1+1,i2,i3,0,nodeTYPE::NODE,0)==FLUID) && (marker(i1+2,i2,i3,0,nodeTYPE::NODE,0)==FLUID) )
            {
              rho = rhoGrid(i1+2,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+2,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
              ForceX -= 2.0*myGrid(i1+2,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
            }
            //M2ZZ - Correction
            if( (marker(i1-1,i2,i3,0,nodeTYPE::NODE,0)!=FLUID) && (marker(i1+1,i2,i3,0,nodeTYPE::NODE,0)==FLUID) )
            {
              rho = rhoGrid(i1+1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
              ForceX -= 2.0*myGrid(i1+1,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO);
              
              rho = rhoGrid(i1-1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
              ForceX += 2.0*myGrid(i1-1,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
            }
          }
        }
      }
    }        
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackNodeG34withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) != FLUID)
          {
            if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2-1,i3-1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-1,i3-1,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_M1_M1) = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1];
              ForceY += myGrid(i1  ,i2-1,i3-1,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_M1_M1);
              ForceZ += myGrid(i1  ,i2-1,i3-1,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_M1_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_P1_P1) = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1];
            }
            if(marker(i1  ,i2+1,i3-1,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2+1,i3-1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+1,i3-1,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_P1_M1) = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_P1_M1];
              ForceY -= myGrid(i1  ,i2+1,i3-1,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_P1_M1);
              ForceZ += myGrid(i1  ,i2+1,i3-1,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_P1_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_M1_P1) = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_M1_P1];
            }
            if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1-1,i2  ,i3-1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2  ,i3-1,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_M1_ZERO_M1) = lbModel.fTemp4[0][lbModel.G4_DV_M1_ZERO_M1];
              ForceX += myGrid(i1-1,i2  ,i3-1,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_M1_ZERO_M1);
              ForceZ += myGrid(i1-1,i2  ,i3-1,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_M1_ZERO_M1);
              
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_P1_ZERO_P1) = lbModel.fTemp3[0][lbModel.G3_DV_P1_ZERO_P1];
            }                           
            if(marker(i1+1,i2  ,i3-1,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2  ,i3-1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2  ,i3-1,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_P1_ZERO_M1) = lbModel.fTemp4[0][lbModel.G4_DV_P1_ZERO_M1];
              ForceX -= myGrid(i1+1,i2  ,i3-1,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_P1_ZERO_M1);
              ForceZ += myGrid(i1+1,i2  ,i3-1,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_P1_ZERO_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_M1_ZERO_P1) = lbModel.fTemp3[0][lbModel.G3_DV_M1_ZERO_P1];
            } 
            
            
            
            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)                                                                                    
            {
              rho = rhoGrid(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+1,i3+1,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_P1_P1) = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1];
              ForceY -= myGrid(i1  ,i2+1,i3+1,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_P1_P1);
              ForceZ -= myGrid(i1  ,i2+1,i3+1,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_P1_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_M1_M1) = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1];
            }                           
            if(marker(i1  ,i2-1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2-1,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-1,i3+1,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_M1_P1) = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_M1_P1];
              ForceY += myGrid(i1  ,i2-1,i3+1,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_M1_P1);
              ForceZ -= myGrid(i1  ,i2-1,i3+1,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_M1_P1);
              
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_P1_M1) = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_P1_M1];
            }         
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2  ,i3+1,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_P1_ZERO_P1) = lbModel.fTemp3[0][lbModel.G3_DV_P1_ZERO_P1];
              ForceX -= myGrid(i1+1,i2  ,i3+1,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_P1_ZERO_P1);
              ForceZ -= myGrid(i1+1,i2  ,i3+1,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_P1_ZERO_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_M1_ZERO_M1) = lbModel.fTemp4[0][lbModel.G4_DV_M1_ZERO_M1];
            }         
            if(marker(i1-1,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1-1,i2  ,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2  ,i3+1,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_M1_ZERO_P1) = lbModel.fTemp3[0][lbModel.G3_DV_M1_ZERO_P1];
              ForceX += myGrid(i1-1,i2  ,i3+1,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_M1_ZERO_P1);
              ForceZ -= myGrid(i1-1,i2  ,i3+1,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_M1_ZERO_P1);
              
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_P1_ZERO_M1) = lbModel.fTemp4[0][lbModel.G4_DV_P1_ZERO_M1];
            }                                     
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackNodeG56withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                 
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) != FLUID)
          {  
            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1-1,i2-1,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2-1,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_M1_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO];
              ForceX += myGrid(i1-1,i2-1,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_M1_M1_ZERO);
              ForceY += myGrid(i1-1,i2-1,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_M1_M1_ZERO);
              
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_P1_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO];
            }         
            if(marker(i1+1,i2-1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2-1,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2-1,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_P1_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_P1_M1_ZERO];
              ForceX -= myGrid(i1+1,i2-1,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_P1_M1_ZERO);
              ForceY += myGrid(i1+1,i2-1,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_P1_M1_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_M1_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_M1_P1_ZERO];
            }             
            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2-1,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-1,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_ZERO_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_ZERO_M1_ZERO];
              ForceY += myGrid(i1  ,i2-1,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_ZERO_M1_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_ZERO_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_ZERO_P1_ZERO];
            }              
            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1-1,i2  ,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2  ,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_M1_ZERO_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_M1_ZERO_ZERO];
              ForceX += myGrid(i1-1,i2  ,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_M1_ZERO_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_P1_ZERO_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_P1_ZERO_ZERO];
            }              
            if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2+1,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_P1_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO];
              ForceX -= myGrid(i1+1,i2+1,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_P1_P1_ZERO);
              ForceY -= myGrid(i1+1,i2+1,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_P1_P1_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_M1_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO];
            }              
            if(marker(i1-1,i2+1,i3  ,0,nodeTYPE::NODE,0) ==FLUID)
            {
              rho = rhoGrid(i1-1,i2+1,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2+1,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_M1_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_M1_P1_ZERO];
              ForceX += myGrid(i1-1,i2+1,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_M1_P1_ZERO);
              ForceY -= myGrid(i1-1,i2+1,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_M1_P1_ZERO);
              
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_P1_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_P1_M1_ZERO];
            }              
            if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+1,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_ZERO_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_ZERO_P1_ZERO];
              ForceY -= myGrid(i1  ,i2+1,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_ZERO_P1_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_ZERO_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_ZERO_M1_ZERO];
            }              
            if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2  ,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_P1_ZERO_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_P1_ZERO_ZERO];
              ForceX -= myGrid(i1+1,i2  ,i3  ,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_P1_ZERO_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_M1_ZERO_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_M1_ZERO_ZERO];
            }              
          }   
        }
      }
    }
  }         
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackNodeG910withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)              
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) != FLUID)
          {  
            if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              rho = rhoGrid(i1-1,i2-1,i3-1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
              ForceX += myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1);
              ForceY += myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1);
              ForceZ += myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  , lbModel.G9,nodeTYPE::NODE, lbModel.G9_DV_P1_P1_P1) = lbModel.fTemp9[0] [lbModel.G9_DV_P1_P1_P1];
            }                  
            if(marker(i1+1,i2-1,i3-1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              rho = rhoGrid(i1+1,i2-1,i3-1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1];
              ForceX -= myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1);
              ForceY += myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1);
              ForceZ += myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  , lbModel.G9,nodeTYPE::NODE, lbModel.G9_DV_M1_P1_P1) = lbModel.fTemp9[0] [lbModel. G9_DV_M1_P1_P1];
            }                             
            if(marker(i1+1,i2+1,i3-1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              rho = rhoGrid(i1+1,i2+1,i3-1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1];
              ForceX -= myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1);
              ForceY -= myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1);
              ForceZ += myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  , lbModel.G9,nodeTYPE::NODE, lbModel.G9_DV_M1_M1_P1) = lbModel.fTemp9[0] [lbModel. G9_DV_M1_M1_P1];
            }                            
            if(marker(i1-1,i2+1,i3-1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              rho = rhoGrid(i1-1,i2+1,i3-1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1];
              ForceX += myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1);
              ForceY -= myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1);
              ForceZ += myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  , lbModel.G9,nodeTYPE::NODE, lbModel.G9_DV_P1_M1_P1) = lbModel.fTemp9[0] [lbModel. G9_DV_P1_M1_P1];
            }                             
            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              rho = rhoGrid(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2+1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_P1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1];
              ForceX -= myGrid(i1+1,i2+1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_P1_P1);
              ForceY -= myGrid(i1+1,i2+1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_P1_P1);
              ForceZ -= myGrid(i1+1,i2+1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_P1_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
            }                             
            if(marker(i1-1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              rho = rhoGrid(i1-1,i2+1,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2+1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_P1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_M1_P1_P1];
              ForceX += myGrid(i1-1,i2+1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_P1_P1);
              ForceY -= myGrid(i1-1,i2+1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_P1_P1);
              ForceZ -= myGrid(i1-1,i2+1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_P1_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1) = lbModel.fTemp10[0] [lbModel. G10_DV_P1_M1_M1];
            }                             
            if(marker(i1-1,i2-1,i3+1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              rho = rhoGrid(i1-1,i2-1,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2-1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_M1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_M1_M1_P1];
              ForceX += myGrid(i1-1,i2-1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_M1_P1);
              ForceY += myGrid(i1-1,i2-1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_M1_P1);
              ForceZ -= myGrid(i1-1,i2-1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_M1_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1) = lbModel.fTemp10[0] [lbModel. G10_DV_P1_P1_M1];
            }                             
            if(marker(i1+1,i2-1,i3+1,0,nodeTYPE::NODE,0)==FLUID)  
            {
              rho = rhoGrid(i1+1,i2-1,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2-1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_M1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_P1_M1_P1];
              ForceX -= myGrid(i1+1,i2-1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_M1_P1);
              ForceY += myGrid(i1+1,i2-1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_M1_P1);
              ForceZ -= myGrid(i1+1,i2-1,i3+1,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_M1_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1) = lbModel.fTemp10[0] [lbModel. G10_DV_M1_P1_M1];
            }
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackCellG12withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                          
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) != FLUID)
          {
            
            //ZZP1
            if(marker(i1,i2,i3-1,0,nodeTYPE::CELL,0)==FLUID) 
            {
              
              rho = rhoGrid(i1,i2,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M1) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];
              ForceZ += myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);                            
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P1) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1];
            } 
            //ZZP2
            if((marker(i1,i2,i3-1,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2,i3-2,0,nodeTYPE::CELL,0)==FLUID)) 
            {
              rho = rhoGrid(i1,i2,i3-2,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3-2,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
              ForceZ += 2.0*myGrid(i1  ,i2  ,i3-2,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2);
              
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
            }        
            //ZZP2 - Correction
            if(    ( (marker(i1,i2,i3-1,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeTYPE::CELL,0)==FLUID) ) 
              || ( (marker(i1,i2,i3-1,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeTYPE::CELL,0)!=FLUID) ) )
            {
              rho = rhoGrid(i1,i2,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
              ForceZ += 2.0*myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2);
              
              rho = rhoGrid(i1,i2,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);                               
              myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
              ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);                               
              myGrid(i1  ,i2  ,i3,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
              myGrid(i1  ,i2  ,i3,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
            }    
            //ZP2Z
            if( (marker(i1,i2-1,i3,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2-2,i3,0,nodeTYPE::CELL,0)==FLUID) )
            {
              rho = rhoGrid(i1,i2-2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-2,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
              ForceY += 2.0*myGrid(i1  ,i2-2,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);                               
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
            }        
            //ZP2Z - Correction
            if(    ((marker(i1,i2-1,i3,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2+1,i3,0,nodeTYPE::CELL,0)==FLUID)) 
              || ((marker(i1,i2-1,i3,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2+1,i3,0,nodeTYPE::CELL,0)!=FLUID)) )
            {
              rho = rhoGrid(i1,i2-1,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-1,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
              ForceY += 2.0*myGrid(i1  ,i2-1,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO);
              
              rho = rhoGrid(i1,i2+1,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+1,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
              ForceY -= 2.0*myGrid(i1  ,i2+1,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);                               
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
            }    
            //P2ZZ
            if( (marker(i1-1,i2,i3,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1-2,i2,i3,0,nodeTYPE::CELL,0)==FLUID) )
            {
              rho = rhoGrid(i1-2,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-2,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
              ForceX += 2.0*myGrid(i1-2,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
            }        
            //P2ZZ - Correction
            if( ((marker(i1-1,i2,i3,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1+1,i2,i3,0,nodeTYPE::CELL,0)==FLUID)) || ( (marker(i1-1,i2,i3,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1+1,i2,i3,0,nodeTYPE::CELL,0)!=FLUID) ) )
            {
              rho = rhoGrid(i1-1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
              ForceX += 2.0*myGrid(i1-1,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO);
              
              rho = rhoGrid(i1+1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
              ForceX -= 2.0*myGrid(i1+1,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);                               
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
            }       
            
            
            
            
            
            
            
            //ZZM1
            if( marker(i1,i2,i3+1,0,nodeTYPE::CELL,0)==FLUID )
            {
              rho = rhoGrid(i1,i2,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P1) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1];
              ForceZ -= myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M1) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];
            }       
            //ZZM2
            if( (marker(i1,i2,i3+1,0,nodeTYPE::CELL,0)==FLUID ) && (marker(i1,i2,i3+2,0,nodeTYPE::CELL,0)==FLUID) )
            {
              rho = rhoGrid(i1,i2,i3+2,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
              ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
            }   
            //                         ZZM2 - Correction
            if( (marker(i1,i2,i3-1,0,nodeTYPE::CELL,0)!=FLUID ) && (marker(i1,i2,i3+1,0,nodeTYPE::CELL,0)==FLUID) )
            {
              rho = rhoGrid(i1,i2,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
              ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2);
              
              rho = rhoGrid(i1,i2,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
              ForceZ += 2.0*myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2);                            
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
            }
            
            //ZM2Z
            if( (marker(i1,i2+1,i3,0,nodeTYPE::CELL,0)==FLUID) && (marker(i1,i2+2,i3,0,nodeTYPE::CELL,0)==FLUID) )
            {
              rho = rhoGrid(i1,i2+2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+2,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
              ForceY -= 2.0*myGrid(i1  ,i2+2,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
            }
            //ZM2Z - Correction
            if( (marker(i1,i2-1,i3,0,nodeTYPE::CELL,0)!=FLUID) && (marker(i1,i2+1,i3,0,nodeTYPE::CELL,0)==FLUID) )
            {
              rho = rhoGrid(i1,i2+1,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+1,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
              ForceY -= 2.0*myGrid(i1  ,i2+1,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO);
              
              rho = rhoGrid(i1,i2-1,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-1,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
              ForceY += 2.0*myGrid(i1  ,i2-1,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
            }
            //M2ZZ
            if( (marker(i1+1,i2,i3,0,nodeTYPE::CELL,0)==FLUID) && (marker(i1+2,i2,i3,0,nodeTYPE::CELL,0)==FLUID) )
            {
              rho = rhoGrid(i1+2,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+2,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
              ForceX -= 2.0*myGrid(i1+2,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
            }
            //M2ZZ - Correction
            if( (marker(i1-1,i2,i3,0,nodeTYPE::CELL,0)!=FLUID) && (marker(i1+1,i2,i3,0,nodeTYPE::CELL,0)==FLUID) )
            {
              rho = rhoGrid(i1+1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
              ForceX -= 2.0*myGrid(i1+1,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO);
              
              rho = rhoGrid(i1-1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
              ForceX += 2.0*myGrid(i1-1,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);   
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
            }
          }
        }
      }
    }        
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackCellG34withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)               
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) != FLUID)
          {
            if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-1,i3-1,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_M1_M1) = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1];
              ForceY += myGrid(i1  ,i2-1,i3-1,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_M1_M1);
              ForceZ += myGrid(i1  ,i2-1,i3-1,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_M1_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_P1_P1) = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1];
            }
            if(marker(i1  ,i2+1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2+1,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+1,i3-1,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_P1_M1) = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_P1_M1];
              ForceY -= myGrid(i1  ,i2+1,i3-1,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_P1_M1);
              ForceZ += myGrid(i1  ,i2+1,i3-1,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_P1_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_M1_P1) = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_M1_P1];
            }
            if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2  ,i3-1,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_M1_ZERO_M1) = lbModel.fTemp4[0][lbModel.G4_DV_M1_ZERO_M1];
              ForceX += myGrid(i1-1,i2  ,i3-1,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_M1_ZERO_M1);
              ForceZ += myGrid(i1-1,i2  ,i3-1,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_M1_ZERO_M1);
              
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_P1_ZERO_P1) = lbModel.fTemp3[0][lbModel.G3_DV_P1_ZERO_P1];
            }                           
            if(marker(i1+1,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2  ,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2  ,i3-1,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_P1_ZERO_M1) = lbModel.fTemp4[0][lbModel.G4_DV_P1_ZERO_M1];
              ForceX -= myGrid(i1+1,i2  ,i3-1,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_P1_ZERO_M1);
              ForceZ += myGrid(i1+1,i2  ,i3-1,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_P1_ZERO_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_M1_ZERO_P1) = lbModel.fTemp3[0][lbModel.G3_DV_M1_ZERO_P1];
            } 
            
            
            
            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::CELL,0)==FLUID)                                                                                    
            {
              rho = rhoGrid(i1  ,i2+1,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+1,i3+1,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_P1_P1) = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1];
              ForceY -= myGrid(i1  ,i2+1,i3+1,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_P1_P1);
              ForceZ -= myGrid(i1  ,i2+1,i3+1,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_P1_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_M1_M1) = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1];
            }                           
            if(marker(i1  ,i2-1,i3+1,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2-1,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-1,i3+1,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_M1_P1) = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_M1_P1];
              ForceY += myGrid(i1  ,i2-1,i3+1,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_M1_P1);
              ForceZ -= myGrid(i1  ,i2-1,i3+1,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_M1_P1);
              
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_P1_M1) = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_P1_M1];
            }         
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2  ,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2  ,i3+1,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_P1_ZERO_P1) = lbModel.fTemp3[0][lbModel.G3_DV_P1_ZERO_P1];
              ForceX -= myGrid(i1+1,i2  ,i3+1,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_P1_ZERO_P1);
              ForceZ -= myGrid(i1+1,i2  ,i3+1,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_P1_ZERO_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_M1_ZERO_M1) = lbModel.fTemp4[0][lbModel.G4_DV_M1_ZERO_M1];
            }         
            if(marker(i1-1,i2  ,i3+1,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1-1,i2  ,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2  ,i3+1,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_M1_ZERO_P1) = lbModel.fTemp3[0][lbModel.G3_DV_M1_ZERO_P1];
              ForceX += myGrid(i1-1,i2  ,i3+1,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_M1_ZERO_P1);
              ForceZ -= myGrid(i1-1,i2  ,i3+1,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_M1_ZERO_P1);
              
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);    
              myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_P1_ZERO_M1) = lbModel.fTemp4[0][lbModel.G4_DV_P1_ZERO_M1];
            }                                     
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackCellG56withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                 
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) != FLUID)
          {  
            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2-1,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_M1_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO];
              ForceX += myGrid(i1-1,i2-1,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_M1_M1_ZERO);
              ForceY += myGrid(i1-1,i2-1,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_M1_M1_ZERO);
              
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_P1_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO];
            }         
            if(marker(i1+1,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2-1,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2-1,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_P1_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_P1_M1_ZERO];
              ForceX -= myGrid(i1+1,i2-1,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_P1_M1_ZERO);
              ForceY += myGrid(i1+1,i2-1,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_P1_M1_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_M1_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_M1_P1_ZERO];
            }             
            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-1,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_ZERO_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_ZERO_M1_ZERO];
              ForceY += myGrid(i1  ,i2-1,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_ZERO_M1_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_ZERO_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_ZERO_P1_ZERO];
            }              
            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2  ,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_M1_ZERO_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_M1_ZERO_ZERO];
              ForceX += myGrid(i1-1,i2  ,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_M1_ZERO_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_P1_ZERO_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_P1_ZERO_ZERO];
            }              
            if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2+1,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2+1,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_P1_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO];
              ForceX -= myGrid(i1+1,i2+1,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_P1_P1_ZERO);
              ForceY -= myGrid(i1+1,i2+1,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_P1_P1_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_M1_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO];
            }              
            if(marker(i1-1,i2+1,i3  ,0,nodeTYPE::CELL,0) ==FLUID)
            {
              rho = rhoGrid(i1-1,i2+1,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2+1,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_M1_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_M1_P1_ZERO];
              ForceX += myGrid(i1-1,i2+1,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_M1_P1_ZERO);
              ForceY -= myGrid(i1-1,i2+1,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_M1_P1_ZERO);
              
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_P1_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_P1_M1_ZERO];
            }              
            if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2+1,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+1,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_ZERO_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_ZERO_P1_ZERO];
              ForceY -= myGrid(i1  ,i2+1,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_ZERO_P1_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_ZERO_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_ZERO_M1_ZERO];
            }              
            if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2  ,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2  ,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_P1_ZERO_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_P1_ZERO_ZERO];
              ForceX -= myGrid(i1+1,i2  ,i3  ,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_P1_ZERO_ZERO);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);  
              myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_M1_ZERO_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_M1_ZERO_ZERO];
            }              
          }   
        }
      }
    }
  }         
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackCellG910withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)               
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) != FLUID)
          {  
            if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              rho = rhoGrid(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
              ForceX += myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1);
              ForceY += myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1);
              ForceZ += myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  , lbModel.G9,nodeTYPE::CELL, lbModel.G9_DV_P1_P1_P1) = lbModel.fTemp9[0] [lbModel.G9_DV_P1_P1_P1];
            }                  
            if(marker(i1+1,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              rho = rhoGrid(i1+1,i2-1,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1];
              ForceX -= myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1);
              ForceY += myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1);
              ForceZ += myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  , lbModel.G9,nodeTYPE::CELL, lbModel.G9_DV_M1_P1_P1) = lbModel.fTemp9[0] [lbModel. G9_DV_M1_P1_P1];
            }                             
            if(marker(i1+1,i2+1,i3-1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              rho = rhoGrid(i1+1,i2+1,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1];
              ForceX -= myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1);
              ForceY -= myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1);
              ForceZ += myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  , lbModel.G9,nodeTYPE::CELL, lbModel.G9_DV_M1_M1_P1) = lbModel.fTemp9[0] [lbModel. G9_DV_M1_M1_P1];
            }                            
            if(marker(i1-1,i2+1,i3-1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              rho = rhoGrid(i1-1,i2+1,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1];
              ForceX += myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1);
              ForceY -= myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1);
              ForceZ += myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  , lbModel.G9,nodeTYPE::CELL, lbModel.G9_DV_P1_M1_P1) = lbModel.fTemp9[0] [lbModel. G9_DV_P1_M1_P1];
            }                             
            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              rho = rhoGrid(i1+1,i2+1,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2+1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_P1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1];
              ForceX -= myGrid(i1+1,i2+1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_P1_P1);
              ForceY -= myGrid(i1+1,i2+1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_P1_P1);
              ForceZ -= myGrid(i1+1,i2+1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_P1_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
            }                             
            if(marker(i1-1,i2+1,i3+1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              rho = rhoGrid(i1-1,i2+1,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2+1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_P1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_M1_P1_P1];
              ForceX += myGrid(i1-1,i2+1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_P1_P1);
              ForceY -= myGrid(i1-1,i2+1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_P1_P1);
              ForceZ -= myGrid(i1-1,i2+1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_P1_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1) = lbModel.fTemp10[0] [lbModel. G10_DV_P1_M1_M1];
            }                             
            if(marker(i1-1,i2-1,i3+1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              rho = rhoGrid(i1-1,i2-1,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2-1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_M1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_M1_M1_P1];
              ForceX += myGrid(i1-1,i2-1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_M1_P1);
              ForceY += myGrid(i1-1,i2-1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_M1_P1);
              ForceZ -= myGrid(i1-1,i2-1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_M1_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1) = lbModel.fTemp10[0] [lbModel. G10_DV_P1_P1_M1];
            }                             
            if(marker(i1+1,i2-1,i3+1,0,nodeTYPE::CELL,0)==FLUID)  
            {
              rho = rhoGrid(i1+1,i2-1,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2-1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_M1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_P1_M1_P1];
              ForceX -= myGrid(i1+1,i2-1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_M1_P1);
              ForceY += myGrid(i1+1,i2-1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_M1_P1);
              ForceZ -= myGrid(i1+1,i2-1,i3+1,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_M1_P1);
              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                              
              myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1) = lbModel.fTemp10[0] [lbModel. G10_DV_M1_P1_M1];
            }
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackNodeCellG78withForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                 
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    ///////////////////
    // Node --> Cell //
    ///////////////////        
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) != FLUID)
          {  
            //PPP
            if(marker(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M) = lbModel.fTemp8[0][lbModel.G8_DV_M_M_M];
              ForceX += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);
              ForceY += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);
              ForceZ += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);
              
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);                            
              myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P) = lbModel.fTemp7[0][lbModel.G7_DV_P_P_P];
              
              //PPP-1
              if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
                rho = rhoGrid(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_P1_P1) = lbModel.fTemp9[0][  lbModel.G9_DV_P1_P1_P1];
                myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
                
                ForceX -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_P1_P1);
                ForceY -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_P1_P1);
                ForceZ -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_P1_P1);
                
                ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1);
                ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1);
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1);
                
              }
            }
            //MPP
            if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M) = lbModel.fTemp8[0][lbModel.G8_DV_P_M_M];
              ForceX -= 0.5*myGrid(i1+1,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);
              ForceY += 0.5*myGrid(i1+1,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);
              ForceZ += 0.5*myGrid(i1+1,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta); 
              myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P) = lbModel.fTemp7[0][lbModel.G7_DV_M_P_P];
              //MPP-1
              if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {   
                rho = rhoGrid(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1  ,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_M1_P1_P1) =  lbModel.fTemp9[0][ lbModel.G9_DV_M1_P1_P1];
                myGrid(i1+1,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1];
                
                ForceX += myGrid(i1  ,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_M1_P1_P1);
                ForceY -= myGrid(i1  ,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_M1_P1_P1);
                ForceZ -= myGrid(i1  ,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_M1_P1_P1);
                
                ForceX -= myGrid(i1+1,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1);
                ForceY += myGrid(i1+1,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1);
                ForceZ += myGrid(i1+1,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1);                            
              }
            }
            //MMP
            if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M) = lbModel.fTemp8[0][lbModel.G8_DV_P_P_M];
              ForceX -= 0.5*myGrid(i1+1,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);
              ForceY -= 0.5*myGrid(i1+1,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);
              ForceZ += 0.5*myGrid(i1+1,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);
              
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta); 
              myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P) = lbModel.fTemp7[0][lbModel.G7_DV_M_M_P];
              //MMP-1
              if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
                rho = rhoGrid(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1  ,i2  ,i3+1, lbModel.G9,nodeTYPE::NODE, lbModel.G9_DV_M1_M1_P1) =  lbModel.fTemp9[0][ lbModel.G9_DV_M1_M1_P1];
                myGrid(i1+1,i2+1,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1];
                
                ForceX += myGrid(i1  ,i2  ,i3+1, lbModel.G9,nodeTYPE::NODE, lbModel.G9_DV_M1_M1_P1);
                ForceY += myGrid(i1  ,i2  ,i3+1, lbModel.G9,nodeTYPE::NODE, lbModel.G9_DV_M1_M1_P1);
                ForceZ -= myGrid(i1  ,i2  ,i3+1, lbModel.G9,nodeTYPE::NODE, lbModel.G9_DV_M1_M1_P1);
                
                ForceX -= myGrid(i1+1,i2+1,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1);
                ForceY -= myGrid(i1+1,i2+1,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1);
                ForceZ += myGrid(i1+1,i2+1,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1);                            
              }
            }            
            //PMP
            if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M) = lbModel.fTemp8[0][lbModel.G8_DV_M_P_M];
              ForceX += 0.5*myGrid(i1  ,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);
              ForceY -= 0.5*myGrid(i1  ,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);
              ForceZ += 0.5*myGrid(i1  ,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);
              
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta); 
              myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P) = lbModel.fTemp7[0][lbModel.G7_DV_P_M_P];
              //PMP-1
              if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
                rho = rhoGrid(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1+1,i2  ,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_M1_P1) =  lbModel.fTemp9[0][ lbModel.G9_DV_P1_M1_P1];
                myGrid(i1  ,i2+1,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1];
                
                ForceX -= myGrid(i1+1,i2  ,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_M1_P1);
                ForceY += myGrid(i1+1,i2  ,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_M1_P1);
                ForceZ -= myGrid(i1+1,i2  ,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_M1_P1);
                
                ForceX += myGrid(i1  ,i2+1,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1);
                ForceY -= myGrid(i1  ,i2+1,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1);
                ForceZ += myGrid(i1  ,i2+1,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1);                            
                
              }    
            }            
            
            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P) = lbModel.fTemp7[0][lbModel.G7_DV_P_P_P];
              ForceX -= 0.5*myGrid(i1+1,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);
              ForceY -= 0.5*myGrid(i1+1,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);
              ForceZ -= 0.5*myGrid(i1+1,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);
              
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta); 
              myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M) = lbModel.fTemp8[0][lbModel.G8_DV_M_M_M];
            }                    
            
            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P) = lbModel.fTemp7[0][lbModel.G7_DV_M_P_P];
              ForceX += 0.5*myGrid(i1  ,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);
              ForceY -= 0.5*myGrid(i1  ,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);
              ForceZ -= 0.5*myGrid(i1  ,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);
              
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta); 
              myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M) = lbModel.fTemp8[0][lbModel.G8_DV_P_M_M];
            }                        
            
            if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P) = lbModel.fTemp7[0][lbModel.G7_DV_M_M_P];
              ForceX += 0.5*myGrid(i1  ,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);
              ForceY += 0.5*myGrid(i1  ,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);
              ForceZ -= 0.5*myGrid(i1  ,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);
              
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta); 
              myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M) = lbModel.fTemp8[0][lbModel.G8_DV_P_P_M];
            }    
            
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              rho = rhoGrid(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P) = lbModel.fTemp7[0][lbModel.G7_DV_P_M_P];
              ForceX -= 0.5*myGrid(i1+1,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);
              ForceY += 0.5*myGrid(i1+1,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);
              ForceZ -= 0.5*myGrid(i1+1,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);
              
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta); 
              myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M) = lbModel.fTemp8[0][lbModel.G8_DV_M_P_M];
            }    
          }
        }
      }
    }
    
    ////////////////////
    // Cell --> Node  //
    ////////////////////             
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0)!=FLUID)
          {               
            if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {
              //PPP
              rho = rhoGrid(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M) = lbModel.fTemp8[0][lbModel.G8_DV_M_M_M];
              ForceX += 0.5*myGrid(i1-1,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);
              ForceY += 0.5*myGrid(i1-1,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);
              ForceZ += 0.5*myGrid(i1-1,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P) = lbModel.fTemp7[0][lbModel.G7_DV_P_P_P];
              //PPP-1
              if(marker(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
              {
                rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_P1_P1_P1) = lbModel.fTemp9[0][  lbModel.G9_DV_P1_P1_P1];
                myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
                
                ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_P1_P1_P1);
                ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_P1_P1_P1);
                ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_P1_P1_P1);
                
                ForceX += myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1);
                ForceY += myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1);
                ForceZ += myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1);                            
              }
            }
            //MPP
            if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {              
              rho = rhoGrid(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M) = lbModel.fTemp8[0][lbModel.G8_DV_P_M_M];
              ForceX -= 0.5*myGrid(i1  ,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);
              ForceY += 0.5*myGrid(i1  ,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);
              ForceZ += 0.5*myGrid(i1  ,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);
              
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P) = lbModel.fTemp7[0][lbModel.G7_DV_M_P_P];
              //MPP-1
              if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
              {             
                rho = rhoGrid(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1-1,i2  ,i3  ,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_P1_P1) = lbModel.fTemp9 [0][ lbModel.G9_DV_M1_P1_P1];
                myGrid(i1  ,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1];
                
                ForceX += myGrid(i1-1,i2  ,i3  ,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_P1_P1);
                ForceY -= myGrid(i1-1,i2  ,i3  ,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_P1_P1);
                ForceZ -= myGrid(i1-1,i2  ,i3  ,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_P1_P1);
                
                ForceX -= myGrid(i1  ,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1);
                ForceY += myGrid(i1  ,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1);
                ForceZ += myGrid(i1  ,i2-1,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1);                            
              }             
            }             
            //MMP
            if(marker(i1  ,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            {              
              rho = rhoGrid(i1  ,i2  ,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M) = lbModel.fTemp8[0][lbModel.G8_DV_P_P_M];
              ForceX -= 0.5*myGrid(i1  ,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);
              ForceY -= 0.5*myGrid(i1  ,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);
              ForceZ += 0.5*myGrid(i1  ,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P) = lbModel.fTemp7[0][lbModel.G7_DV_M_M_P];
              //MMP-1
              if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
              {             
                rho = rhoGrid(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1-1,i2-1,i3  ,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_M1_P1) =  lbModel.fTemp9[0][ lbModel.G9_DV_M1_M1_P1];
                myGrid(i1  ,i2  ,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1];
                
                ForceX += myGrid(i1-1,i2-1,i3  ,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_M1_P1);
                ForceY += myGrid(i1-1,i2-1,i3  ,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_M1_P1);
                ForceZ -= myGrid(i1-1,i2-1,i3  ,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_M1_P1);
                
                ForceX -= myGrid(i1  ,i2  ,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1);
                ForceY -= myGrid(i1  ,i2  ,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1);
                ForceZ += myGrid(i1  ,i2  ,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1);                            
              }                             
            }     
            //PMP
            if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)                                                                                    
            {      
              rho = rhoGrid(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M) = lbModel.fTemp8[0][lbModel.G8_DV_M_P_M];
              ForceX += 0.5*myGrid(i1-1,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);
              ForceY -= 0.5*myGrid(i1-1,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);
              ForceZ += 0.5*myGrid(i1-1,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);
              
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P) = lbModel.fTemp7[0][lbModel.G7_DV_P_M_P];
              //PMP-1
              if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
              {    
                rho = rhoGrid(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1  ,i2-1,i3  ,lbModel.G9, nodeTYPE::CELL, lbModel.G9_DV_P1_M1_P1) =  lbModel.fTemp9[0][ lbModel.G9_DV_P1_M1_P1];
                myGrid(i1-1,i2  ,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1];
                
                ForceX -= myGrid(i1  ,i2-1,i3  ,lbModel.G9, nodeTYPE::CELL, lbModel.G9_DV_P1_M1_P1);
                ForceY += myGrid(i1  ,i2-1,i3  ,lbModel.G9, nodeTYPE::CELL, lbModel.G9_DV_P1_M1_P1);
                ForceZ -= myGrid(i1  ,i2-1,i3  ,lbModel.G9, nodeTYPE::CELL, lbModel.G9_DV_P1_M1_P1);
                
                ForceX += myGrid(i1-1,i2  ,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1);
                ForceY -= myGrid(i1-1,i2  ,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1);
                ForceZ += myGrid(i1-1,i2  ,i3-1,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1);
                
              }      
            }       
            
            if(marker(i1  ,i2   ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {    
              rho = rhoGrid(i1  ,i2   ,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2   ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P) = lbModel.fTemp7[0][lbModel.G7_DV_P_P_P];
              ForceX -= 0.5*myGrid(i1  ,i2   ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);
              ForceY -= 0.5*myGrid(i1  ,i2   ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);
              ForceZ -= 0.5*myGrid(i1  ,i2   ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M) = lbModel.fTemp8[0][lbModel.G8_DV_M_M_M];
            }                 
            
            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {    
              rho = rhoGrid(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P) = lbModel.fTemp7[0][lbModel.G7_DV_M_P_P];
              ForceX += 0.5*myGrid(i1-1,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);
              ForceY -= 0.5*myGrid(i1-1,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);
              ForceZ -= 0.5*myGrid(i1-1,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);
              
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M) = lbModel.fTemp8[0][lbModel.G8_DV_P_M_M];
            }        
            
            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {    
              rho = rhoGrid(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P) = lbModel.fTemp7[0][lbModel.G7_DV_M_M_P];
              ForceX += 0.5*myGrid(i1-1,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);
              ForceY += 0.5*myGrid(i1-1,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);
              ForceZ -= 0.5*myGrid(i1-1,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);
              
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M) = lbModel.fTemp8[0][lbModel.G8_DV_P_P_M];
            }    
            
            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            {    
              rho = rhoGrid(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P) = lbModel.fTemp7[0][lbModel.G7_DV_P_M_P];
              ForceX -= 0.5*myGrid(i1  ,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);
              ForceY += 0.5*myGrid(i1  ,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);
              ForceZ -= 0.5*myGrid(i1  ,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);
              
              
              rho = rhoGrid(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M) = lbModel.fTemp8[0][lbModel.G8_DV_M_P_M];
            }                
          }  
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1> 
  void diffusebounceBackwithForces(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)            
  {
    getRhoForDiffuseBounceBackNodeG12 (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID);
    getRhoForDiffuseBounceBackNodeG34 (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID);
    getRhoForDiffuseBounceBackNodeG56 (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID);
    getRhoForDiffuseBounceBackNodeG910(lbModel,myGrid,marker,rhoGrid,FLUID,SOLID);         
    getRhoForDiffuseBounceBackCellG12 (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID);
    getRhoForDiffuseBounceBackCellG34 (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID);
    getRhoForDiffuseBounceBackCellG56 (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID);
    getRhoForDiffuseBounceBackCellG910(lbModel,myGrid,marker,rhoGrid,FLUID,SOLID);         
    getRhoForDiffuseBounceBackNodeCellG78(lbModel,myGrid,marker,rhoGrid,FLUID,SOLID);     
    
    
    diffusebounceBackNodeG12withForce (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackNodeG34withForce (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackNodeG56withForce (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackNodeG910withForce(lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);       
    diffusebounceBackCellG12withForce (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackCellG34withForce (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackCellG56withForce (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackCellG910withForce(lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackNodeCellG78withForce(lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);    
  }
  
  
  
  
  
  template <int N,int numblock, typename dataType1>
  void getRhoAll(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, dataType1> &rhoGrid)            
  {
    dataType1 rho(0.0);
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                          
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                  
        {
          getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2,i3);
          rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0) = rho; 
        }
      }
    }
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                          
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)  
        {    
              getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2,i3);
              rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0) = rho; 
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackNodeG12Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    
//     for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
//     {                                                                
//       for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                          
//       {                    
//         for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
//         {
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                          
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
        {    
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) == FLUID)
          {
            //ZZP1
            if(marker(i1,i2,i3+1,0,nodeTYPE::NODE,0)!=FLUID) 
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3 ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M1) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];

            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M1);
            } 
            
            //ZZP2
            if((marker(i1,i2,i3+2,0,nodeTYPE::NODE,0)!=FLUID)) 
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3 ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];

             if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceZ += 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2);
            } 
            if( (marker(i1,i2,i3+1,0,nodeTYPE::NODE,0) != FLUID && marker(i1,i2,i3+2,0,nodeTYPE::NODE,0) == FLUID) )
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3 ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];

              rho = rhoGrid(i1,i2,i3+2,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3+2 ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];

            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {              
              ForceZ += 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_ZERO_M2);
              ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2);
            }
            }
            
            
            //ZP2Z
            if( (marker(i1,i2+2,i3,0,nodeTYPE::NODE,0)!=FLUID))
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];

            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
             ForceY += 2.0*myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO);
            }
              if( (marker(i1,i2+1,i3,0,nodeTYPE::NODE,0) != FLUID && marker(i1,i2+2,i3,0,nodeTYPE::NODE,0) == FLUID) )
              {
                rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG12(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1  ,i2  ,i3 ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
               
                rho = rhoGrid(i1,i2+2,i3,0,nodeTYPE::NODE,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG12(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1  ,i2+2,i3 ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
                
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceY += 2.0*myGrid(i1  ,i2    ,i3,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_ZERO_M2_ZERO);
                ForceY -= 2.0*myGrid(i1  ,i2+2  ,i3,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO);
            }
              }                    
            
            //P2ZZ
            if((marker(i1+2,i2,i3,0,nodeTYPE::NODE,0)!=FLUID))
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceX += 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO);
            }    
            
              if( (marker(i1+1,i2,i3,0,nodeTYPE::NODE,0) != FLUID && marker(i1+2,i2,i3,0,nodeTYPE::NODE,0) == FLUID) )
              {
                rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG12(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1  ,i2  ,i3 ,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
               
                rho = rhoGrid(i1+2,i2,i3,0,nodeTYPE::NODE,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG12(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1+2,i2,i3 ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
                
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += 2.0*myGrid(i1  ,i2,i3,lbModel.G2,nodeTYPE::NODE,lbModel.G2_DV_M2_ZERO_ZERO);
                ForceX -= 2.0*myGrid(i1+2,i2,i3,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO);
            }   
              }              
            
            //ZZM1
            if(marker(i1,i2,i3-1,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P1) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P1);
            }       
            //ZZM2
            if((marker(i1,i2,i3-2,0,nodeTYPE::NODE,0)!=FLUID))
            {
              rho = rhoGrid(i1,i2,i3 ,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3 ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceZ -= 2.0*myGrid(i1  ,i2  ,i3 ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_ZERO_P2);
              
            }   
            //ZM2Z
            if((marker(i1,i2-2,i3,0,nodeTYPE::NODE,0)!=FLUID))
            {
              rho = rhoGrid(i1,i2 ,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2 ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceY -= 2.0*myGrid(i1  ,i2 ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_ZERO_P2_ZERO);
            }
            //M2ZZ
            if((marker(i1-2,i2,i3,0,nodeTYPE::NODE,0)!=FLUID))
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1 ,i2  ,i3  ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceX -= 2.0*myGrid(i1 ,i2  ,i3 ,lbModel.G1,nodeTYPE::NODE,lbModel.G1_DV_P2_ZERO_ZERO);
            }
          }
        }
      }
    }        
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackCellG12Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                          
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
        {  
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) == FLUID)
          {
            //ZZP1
            if(marker(i1,i2,i3+1,0,nodeTYPE::CELL,0)!=FLUID) 
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3 ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M1) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M1);
            } 
            
            //ZZP2
            if((marker(i1,i2,i3+2,0,nodeTYPE::CELL,0)!=FLUID)) 
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3 ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceZ += 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2);
            } 
            
            if( (marker(i1,i2,i3+1,0,nodeTYPE::CELL,0) != FLUID && marker(i1,i2,i3+2,0,nodeTYPE::CELL,0) == FLUID) )
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              
              myGrid(i1  ,i2  ,i3 ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];

              rho = rhoGrid(i1,i2,i3+2,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3+2 ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
              
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceZ += 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_ZERO_M2);
              ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2);
            }              
            }
            
            //ZP2Z
            if( (marker(i1,i2+2,i3,0,nodeTYPE::CELL,0)!=FLUID))
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceY += 2.0*myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO);
            }
            
              if( (marker(i1,i2+1,i3,0,nodeTYPE::CELL,0) != FLUID && marker(i1,i2+2,i3,0,nodeTYPE::CELL,0) == FLUID) )
              {
                rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG12(lbModel,rho,uX,uY,uZ,theta);
                
                myGrid(i1  ,i2  ,i3 ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
               
                rho = rhoGrid(i1,i2+2,i3,0,nodeTYPE::CELL,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG12(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1  ,i2+2,i3 ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
                
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceY += 2.0*myGrid(i1  ,i2    ,i3,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_ZERO_M2_ZERO);
                ForceY -= 2.0*myGrid(i1  ,i2+2  ,i3,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO);
            }
              }                    
            
            
            
            
            //P2ZZ
            if((marker(i1+2,i2,i3,0,nodeTYPE::CELL,0)!=FLUID))
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceX += 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO);
            }    
            
              if( (marker(i1+1,i2,i3,0,nodeTYPE::CELL,0) != FLUID && marker(i1+2,i2,i3,0,nodeTYPE::CELL,0) == FLUID) )
              {
                rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG12(lbModel,rho,uX,uY,uZ,theta);
                
                myGrid(i1  ,i2  ,i3 ,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO) = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
               
                rho = rhoGrid(i1+2,i2,i3,0,nodeTYPE::CELL,0); 
                uX = uXBody;
                uY = uYBody;
                uZ = uZBody;
                theta = lbModel.theta0;
                getFEqG12(lbModel,rho,uX,uY,uZ,theta);
                myGrid(i1+2,i2,i3 ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
                
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += 2.0*myGrid(i1  ,i2,i3,lbModel.G2,nodeTYPE::CELL,lbModel.G2_DV_M2_ZERO_ZERO);
                ForceX -= 2.0*myGrid(i1+2,i2,i3,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO);
              }   
              }              
            
            //ZZM1
            if(marker(i1,i2,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P1) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P1);
            }       
            //ZZM2
            if((marker(i1,i2,i3-2,0,nodeTYPE::CELL,0)!=FLUID))
            {
              rho = rhoGrid(i1,i2,i3 ,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3 ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceZ -= 2.0*myGrid(i1  ,i2  ,i3 ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_ZERO_P2);
              
            }   
            //ZM2Z
            if((marker(i1,i2-2,i3,0,nodeTYPE::CELL,0)!=FLUID))
            {
              rho = rhoGrid(i1,i2 ,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2 ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceY -= 2.0*myGrid(i1  ,i2 ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_ZERO_P2_ZERO);
            }
            //M2ZZ
            if((marker(i1-2,i2,i3,0,nodeTYPE::CELL,0)!=FLUID))
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG12(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1 ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO) = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceX -= 2.0*myGrid(i1 ,i2  ,i3  ,lbModel.G1,nodeTYPE::CELL,lbModel.G1_DV_P2_ZERO_ZERO);
            }
          }
        }
      }
    }        
  }
    
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackNodeG34Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                          
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
        {  
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) == FLUID)
          {
            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_M1_M1) = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceY += myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_M1_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_M1_M1);
            }
                
            }
            if(marker(i1  ,i2-1,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_P1_M1) = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_P1_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
              ForceY -= myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_P1_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_ZERO_P1_M1);
            }
            }
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_M1_ZERO_M1) = lbModel.fTemp4[0][lbModel.G4_DV_M1_ZERO_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
              ForceX += myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_M1_ZERO_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_M1_ZERO_M1);
            }
                
            }                           
            if(marker(i1-1,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_P1_ZERO_M1) = lbModel.fTemp4[0][lbModel.G4_DV_P1_ZERO_M1];

            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
              ForceX -= myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_P1_ZERO_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::NODE,lbModel.G4_DV_P1_ZERO_M1);
            } 
            }
            if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::NODE,0)!=FLUID)                                                                                    
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_P1_P1) = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
              ForceY -= myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_P1_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_P1_P1);
            }
                
            }                           
            if(marker(i1  ,i2+1,i3-1,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_M1_P1) = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_M1_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
              ForceY += myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_M1_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_ZERO_M1_P1);
            }
                
            }         
            if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_P1_ZERO_P1) = lbModel.fTemp3[0][lbModel.G3_DV_P1_ZERO_P1];

            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_P1_ZERO_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_P1_ZERO_P1);
            }  
                
            }      
            if(marker(i1+1,i2  ,i3-1,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_M1_ZERO_P1) = lbModel.fTemp3[0][lbModel.G3_DV_M1_ZERO_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_M1_ZERO_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::NODE,lbModel.G3_DV_M1_ZERO_P1);
            }
            }
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackCellG34Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                          
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
        {  
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) == FLUID)
          {
            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_M1_M1) = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceY += myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_M1_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_M1_M1);
            }
                
            }
            if(marker(i1  ,i2-1,i3+1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_P1_M1) = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_P1_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceY -= myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_P1_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_ZERO_P1_M1);
            }
                
            }
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_M1_ZERO_M1) = lbModel.fTemp4[0][lbModel.G4_DV_M1_ZERO_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_M1_ZERO_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_M1_ZERO_M1);
            }
                
            }                          
            if(marker(i1-1,i2  ,i3+1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_P1_ZERO_M1) = lbModel.fTemp4[0][lbModel.G4_DV_P1_ZERO_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_P1_ZERO_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G4,nodeTYPE::CELL,lbModel.G4_DV_P1_ZERO_M1);
            }
                
            } 
            if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0)!=FLUID)                                                                                    
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_P1_P1) = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceY -= myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_P1_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_P1_P1);
            }
            }
            if(marker(i1  ,i2+1,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_M1_P1) = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_M1_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceY += myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_M1_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_ZERO_M1_P1);
            }
            }
            if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_P1_ZERO_P1) = lbModel.fTemp3[0][lbModel.G3_DV_P1_ZERO_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {             
              ForceX -= myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_P1_ZERO_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_P1_ZERO_P1);
            }
            }
            if(marker(i1+1,i2  ,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG34(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_M1_ZERO_P1) = lbModel.fTemp3[0][lbModel.G3_DV_M1_ZERO_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_M1_ZERO_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G3,nodeTYPE::CELL,lbModel.G3_DV_M1_ZERO_P1);
            }   
                
            }                                  
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackNodeG56Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                 
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                          
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
        {  
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) == FLUID)
          {  
            if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_M1_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_M1_M1_ZERO);
              ForceY += myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_M1_M1_ZERO);
            }   
            }
            if(marker(i1-1,i2+1,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_P1_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_P1_M1_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {             
              ForceX -= myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_P1_M1_ZERO);
              ForceY += myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_P1_M1_ZERO);
            }
            }
            if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_ZERO_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_ZERO_M1_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
             ForceY += myGrid(i1,i2,i3 ,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_ZERO_M1_ZERO);
            }              
            if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_M1_ZERO_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_M1_ZERO_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
             ForceX += myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::NODE,lbModel.G6_DV_M1_ZERO_ZERO);
            }              
            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_P1_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_P1_P1_ZERO);
              ForceY -= myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_P1_P1_ZERO);
            } 
            }
            if(marker(i1+1,i2-1,i3  ,0,nodeTYPE::NODE,0) !=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_M1_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_M1_P1_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {             
              ForceX += myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_M1_P1_ZERO);
              ForceY -= myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_M1_P1_ZERO);
            }
            }
            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_ZERO_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_ZERO_P1_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
             ForceY -= myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_ZERO_P1_ZERO);
            }              
            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_P1_ZERO_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_P1_ZERO_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
             ForceX -= myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::NODE,lbModel.G5_DV_P1_ZERO_ZERO);
            }              
          }   
        }
      }
    }
  }         
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackCellG56Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                 
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                          
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
        {  
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) == FLUID)
          {  
            if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_M1_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_M1_M1_ZERO);
              ForceY += myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_M1_M1_ZERO);
            }   
            }
            if(marker(i1-1,i2+1,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_P1_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_P1_M1_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_P1_M1_ZERO);
              ForceY += myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_P1_M1_ZERO);
            }  
                
            }           
            if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_ZERO_M1_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_ZERO_M1_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              ForceY += myGrid(i1,i2,i3 ,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_ZERO_M1_ZERO);
            }              
            if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_M1_ZERO_ZERO) = lbModel.fTemp6[0][lbModel.G6_DV_M1_ZERO_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
             ForceX += myGrid(i1,i2,i3,lbModel.G6,nodeTYPE::CELL,lbModel.G6_DV_M1_ZERO_ZERO);
            }              
            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_P1_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_P1_P1_ZERO);
              ForceY -= myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_P1_P1_ZERO);
            }
            }
            if(marker(i1+1,i2-1,i3  ,0,nodeTYPE::CELL,0) !=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_M1_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_M1_P1_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {             
              ForceX += myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_M1_P1_ZERO);
              ForceY -= myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_M1_P1_ZERO);
            }
            }
            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_ZERO_P1_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_ZERO_P1_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
             ForceY -= myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_ZERO_P1_ZERO);
            }              
            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG56(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_P1_ZERO_ZERO) = lbModel.fTemp5[0][lbModel.G5_DV_P1_ZERO_ZERO];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
             ForceX -= myGrid(i1,i2,i3,lbModel.G5,nodeTYPE::CELL,lbModel.G5_DV_P1_ZERO_ZERO);
            }              
          }   
        }
      }
    }
  }         
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackNodeG910Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)              
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                          
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
        {  
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) == FLUID)
          {  
            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1);
              ForceY += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1);
            }   
            }
            if(marker(i1-1,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1);
              ForceY += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1);
            }   
            }
            if(marker(i1-1,i2-1,i3+1,0,nodeTYPE::NODE,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1);
              ForceY -= myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1);
            }   
            }
            if(marker(i1+1,i2-1,i3+1,0,nodeTYPE::NODE,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1);
              ForceY -= myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1);
            }  
            }
            
            
              if (marker(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID && marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
                
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1)  = lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
              
              rho = rhoGrid(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_P1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1];      
              
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1);
                ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_M1_M1);
                ForceX -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_P1_P1);
                ForceY -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_P1_P1);                       
                ForceZ -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_P1_P1);                
              }
              } 
              if((marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID && marker(i1-1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)) 
              {
                
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1)  = lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1];
              
              rho = rhoGrid(i1-1,i2+1,i3+1,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_M1_P1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_M1_P1_P1];  
              
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {                ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1);
                ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_M1_M1);  
                ForceX += myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_M1_P1_P1);
                ForceY -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_M1_P1_P1);                       
                ForceZ -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_M1_P1_P1);  
              }      
              } 
              if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID && marker(i1-1,i2-1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
               rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
               uX = uXBody;
               uY = uYBody;
               uZ = uZBody;
               theta = lbModel.theta0;
               getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
               myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1)  = lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1];
               
               rho = rhoGrid(i1-1,i2-1,i3+1,0,nodeTYPE::NODE,0); 
               uX = uXBody;
               uY = uYBody;
               uZ = uZBody;
               theta = lbModel.theta0;
               getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
               myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_M1_M1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_M1_M1_P1]; 
               
      
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1);
               ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1);                       
               ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_P1_P1_M1);                 
               ForceX += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_M1_M1_P1);
               ForceY += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_M1_M1_P1);                       
               ForceZ -= myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_M1_M1_P1); 
              }
              }
              
              if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID && marker(i1+1,i2-1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
               rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
               uX = uXBody;
               uY = uYBody;
               uZ = uZBody;
               theta = lbModel.theta0;
               getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
               myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1)  = lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1];
               
               rho = rhoGrid(i1+1,i2-1,i3+1,0,nodeTYPE::NODE,0); 
               uX = uXBody;
               uY = uYBody;
               uZ = uZBody;
               theta = lbModel.theta0;
               getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
               myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_M1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_P1_M1_P1]; 
               
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1);
                ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::NODE,lbModel.G10_DV_M1_P1_M1); 
                ForceX -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_M1_P1);
                ForceY += myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_M1_P1);                       
                ForceZ -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::NODE, lbModel.G9_DV_P1_M1_P1); 
              }            
            
              }            
            
            if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::NODE,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_P1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
              ForceX -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_P1_P1);
              ForceY -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_P1_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_P1_P1);
            }
            }
            if(marker(i1+1,i2-1,i3-1,0,nodeTYPE::NODE,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_P1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_M1_P1_P1];
             
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_P1_P1);
              ForceY -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_P1_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_P1_P1);
            }   
            }
            if(marker(i1+1,i2+1,i3-1,0,nodeTYPE::NODE,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_M1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_M1_M1_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_M1_P1);
              ForceY += myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_M1_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_M1_M1_P1);
            }   
            }
            if(marker(i1-1,i2+1,i3-1,0,nodeTYPE::NODE,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_M1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_P1_M1_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_M1_P1);
              ForceY += myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_M1_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::NODE,lbModel.G9_DV_P1_M1_P1);
            }
            }
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackCellG910Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)              
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                          
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
        {  
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0) == FLUID)
          {  
            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::CELL,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];

            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1);
              ForceY += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1);
            }   
            }
            if(marker(i1-1,i2+1,i3+1,0,nodeTYPE::CELL,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1];
             
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1);
              ForceY += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1);
            }   
            }
            if(marker(i1-1,i2-1,i3+1,0,nodeTYPE::CELL,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1);
              ForceY -= myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1);
            }   
            }
            if(marker(i1+1,i2-1,i3+1,0,nodeTYPE::CELL,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1) = lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1];
             
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1);
              ForceY -= myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1);
              ForceZ += myGrid(i1,i2,i3,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1);
            }   
            }
            
              if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID && marker(i1+1,i2+1,i3+1,0,nodeTYPE::CELL,0)==FLUID)
              {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1)  = lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
              
              rho = rhoGrid(i1+1,i2+1,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_P1_P1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1];      
              
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1);
                ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_M1_M1);
                ForceX -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_P1_P1_P1);
                ForceY -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_P1_P1_P1);                       
                ForceZ -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_P1_P1_P1);                
            } 
              }
              if( (marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID && marker(i1-1,i2+1,i3+1,0,nodeTYPE::CELL,0)==FLUID))
              {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1)  = lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1];
              
              rho = rhoGrid(i1-1,i2+1,i3+1,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_P1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_M1_P1_P1];  
              
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1);
                ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_M1_M1);  
                ForceX += myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_P1_P1);
                ForceY -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_P1_P1);                       
                ForceZ -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_P1_P1);  
            } 
               }
              if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID && marker(i1-1,i2-1,i3+1,0,nodeTYPE::CELL,0)==FLUID )
              {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
               uX = uXBody;
               uY = uYBody;
               uZ = uZBody;
               theta = lbModel.theta0;
               getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
               myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1)  = lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1];
               
               rho = rhoGrid(i1-1,i2-1,i3+1,0,nodeTYPE::CELL,0); 
               uX = uXBody;
               uY = uYBody;
               uZ = uZBody;
               theta = lbModel.theta0;
               getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
               myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_M1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_M1_M1_P1]; 
               
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {      
               ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1);
               ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1);                       
               ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_P1_P1_M1);                 
               ForceX += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_M1_P1);
               ForceY += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_M1_P1);                       
               ForceZ -= myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_M1_M1_P1); 
            }
              }
              if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID && marker(i1+1,i2-1,i3+1,0,nodeTYPE::CELL,0)==FLUID)
              {
               rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
               uX = uXBody;
               uY = uYBody;
               uZ = uZBody;
               theta = lbModel.theta0;
               getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
               myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1)  = lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1];
               
               rho = rhoGrid(i1+1,i2-1,i3+1,0,nodeTYPE::CELL,0); 
               uX = uXBody;
               uY = uYBody;
               uZ = uZBody;
               theta = lbModel.theta0;
               getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
               myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_P1_M1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_P1_M1_P1]; 
               
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1);
                ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeTYPE::CELL,lbModel.G10_DV_M1_P1_M1); 
                ForceX -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_P1_M1_P1);
                ForceY += myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_P1_M1_P1);                       
                ForceZ -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeTYPE::CELL, lbModel.G9_DV_P1_M1_P1); 
            } 
               }

               
            
            if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_P1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_P1_P1);
              ForceY -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_P1_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_P1_P1);
            }   
            }
            if(marker(i1+1,i2-1,i3-1,0,nodeTYPE::CELL,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_P1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_M1_P1_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_P1_P1);
              ForceY -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_P1_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_P1_P1);
            }   
            }
            if(marker(i1+1,i2+1,i3-1,0,nodeTYPE::CELL,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_M1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_M1_M1_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_M1_P1);
              ForceY += myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_M1_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_M1_M1_P1);
            }   
            }
            if(marker(i1-1,i2+1,i3-1,0,nodeTYPE::CELL,0)!=FLUID)  
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_M1_P1) = lbModel.fTemp9[0][lbModel.G9_DV_P1_M1_P1];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_M1_P1);
              ForceY += myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_M1_P1);
              ForceZ -= myGrid(i1,i2,i3,lbModel.G9,nodeTYPE::CELL,lbModel.G9_DV_P1_M1_P1);
            }
            }
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void diffusebounceBackNodeCellG78Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                 
  {
    dataType1 rho(0.0),uX(0.0),uY(0.0),uZ(0.0),theta(0.0); 
    ///////////////////
    // Node --> Cell //
    ///////////////////        
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                          
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
        {  
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) == FLUID)
          {  
            //PPP
            if(marker(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3 ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M) = lbModel.fTemp8[0][lbModel.G8_DV_M_M_M];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
              ForceX += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);
              ForceY += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);
              ForceZ += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);
            }
            }
            //MPP
            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M) = lbModel.fTemp8[0][lbModel.G8_DV_P_M_M];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);
              ForceY += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);
              ForceZ += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);
            }
            }
            //MMP
            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M) = lbModel.fTemp8[0][lbModel.G8_DV_P_P_M];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {             
              ForceX -= 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);
              ForceY -= 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);
              ForceZ += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);
            }
            }
            //PMP
            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M) = lbModel.fTemp8[0][lbModel.G8_DV_M_P_M];
             
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);
              ForceY -= 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);
              ForceZ += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);
            }   
            }
            
            if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P) = lbModel.fTemp7[0][lbModel.G7_DV_P_P_P];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {             
              ForceX -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);
              ForceY -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);
              ForceZ -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);
            }
            }
            
            if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P) = lbModel.fTemp7[0][lbModel.G7_DV_M_P_P];

            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);
              ForceY -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);
              ForceZ -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);
            }                        
            }
            
            if(marker(i1  ,i2  ,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P) = lbModel.fTemp7[0][lbModel.G7_DV_M_M_P];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);
              ForceY += 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);
              ForceZ -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);
            }   
            }
            
            if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P) = lbModel.fTemp7[0][lbModel.G7_DV_P_M_P];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {             
              ForceX -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);
              ForceY += 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);
              ForceZ -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);
            }
            }
          }
        }
      }
    }
    
    ////////////////////
    // Cell --> Node  //
    ////////////////////             
    for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                          
      {                    
        for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
        {  
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0)==FLUID)
          {               
            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            {
              //PPP
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M) = lbModel.fTemp8[0][lbModel.G8_DV_M_M_M];
             
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);
              ForceY += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);
              ForceZ += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);
            }
            }
            //MPP
            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            {              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M) = lbModel.fTemp8[0][lbModel.G8_DV_P_M_M];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);
              ForceY += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);
              ForceZ += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);
            }   
            }
            //MMP
            if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            {              
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M) = lbModel.fTemp8[0][lbModel.G8_DV_P_P_M];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);
              ForceY -= 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);
              ForceZ += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);
            }   
            }
            //PMP
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID)                                                                                    
            {      
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M) = lbModel.fTemp8[0][lbModel.G8_DV_M_P_M];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);
              ForceY -= 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);
              ForceZ += 0.5*myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);
            }   
            }
            
            if(marker(i1  ,i2   ,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            {    
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1  ,i2   ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P) = lbModel.fTemp7[0][lbModel.G7_DV_P_P_P];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= 0.5*myGrid(i1  ,i2   ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);
              ForceY -= 0.5*myGrid(i1  ,i2   ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);
              ForceZ -= 0.5*myGrid(i1  ,i2   ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);
            }                 
            }
            if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            {    
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P) = lbModel.fTemp7[0][lbModel.G7_DV_M_P_P];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);
              ForceY -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);
              ForceZ -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);
            }   
            }
            
            if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            {    
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P) = lbModel.fTemp7[0][lbModel.G7_DV_M_M_P];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX += 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);
              ForceY += 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);
              ForceZ -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);
            }    
            }
            if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            {    
              rho = rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0); 
              uX = uXBody;
              uY = uYBody;
              uZ = uZBody;
              theta = lbModel.theta0;
              getFEqG78910(lbModel,rho,uX,uY,uZ,theta);
              myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P) = lbModel.fTemp7[0][lbModel.G7_DV_P_M_P];
            if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
            {
                ForceX -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);
              ForceY += 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);
              ForceZ -= 0.5*myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);
            }   
            }
          }  
        }
      }
    }
  }
  
  
  template <int N,int numblock, typename dataType1> 
  void diffusebounceBackSolid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID, dataType1 uXBody, dataType1 uYBody, dataType1 uZBody, dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)            
  {
    getRhoAll(lbModel,myGrid,rhoGrid);
//     
    diffusebounceBackNodeG12Solid (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackNodeG34Solid (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackNodeG56Solid (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackNodeG910Solid(lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);       
    diffusebounceBackCellG12Solid (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackCellG34Solid (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackCellG56Solid (lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackCellG910Solid(lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);
    diffusebounceBackNodeCellG78Solid(lbModel,myGrid,marker,rhoGrid,FLUID,SOLID,uXBody,uYBody,uZBody,ForceX,ForceY,ForceZ);    
  }
  
  
// Explicit Declarations
template void getRhoNodeSinglePoint<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,double &, int ,int ,int );
template void getRhoCellSinglePoint<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,double &, int ,int ,int );
template void getFEqG12<double>(lbmRD3Q41<double> &,double ,double ,double ,double ,double );
template void getFEqG34<double>(lbmRD3Q41<double> &,double ,double ,double ,double ,double );
template void getFEqG56<double>(lbmRD3Q41<double> &,double ,double ,double ,double ,double );
template void getFEqG78910<double>(lbmRD3Q41<double> &,double ,double ,double ,double ,double );

template void getRhoForDiffuseBounceBackNodeG12<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> &,int ,int );            
template void getRhoForDiffuseBounceBackNodeG34<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> &,int ,int );            
template void getRhoForDiffuseBounceBackNodeG56<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> &,int ,int );            
template void getRhoForDiffuseBounceBackNodeG910<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> &,int ,int );            
template void getRhoForDiffuseBounceBackCellG12<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> &,int ,int );            
template void getRhoForDiffuseBounceBackCellG34<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> &,int ,int );            
template void getRhoForDiffuseBounceBackCellG56<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> &,int ,int );            
template void getRhoForDiffuseBounceBackCellG910<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> &,int ,int );            
template void getRhoForDiffuseBounceBackNodeCellG78<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> &,int ,int );            

template void diffusebounceBackNodeG12withForce<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> & ,int , int , double , double , double , double &, double &, double &);                
template void diffusebounceBackNodeG34withForce<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> & ,int , int , double , double , double , double &, double &, double &);                
template void diffusebounceBackNodeG56withForce<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> & ,int , int , double , double , double , double &, double &, double &);                
template void diffusebounceBackNodeG910withForce<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> & ,int , int , double , double , double , double &, double &, double &);                
template void diffusebounceBackCellG12withForce<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> & ,int , int , double , double , double , double &, double &, double &);                
template void diffusebounceBackCellG34withForce<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> & ,int , int , double , double , double , double &, double &, double &);                
template void diffusebounceBackCellG56withForce<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> & ,int , int , double , double , double , double &, double &, double &);                
template void diffusebounceBackCellG910withForce<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> & ,int , int , double , double , double , double &, double &, double &);                
template void diffusebounceBackNodeCellG78withForce<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> & ,int , int , double , double , double , double &, double &, double &);                
template void diffusebounceBackwithForces<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,int> &,gridBCC3D<1,1,double> & ,int , int , double , double , double , double &, double &, double &);            

template void getRhoAll<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1,1,double> &);            
template void diffusebounceBackNodeG12Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,gridBCC3D<1, 1, double> & ,int ,int , double , double , double , double &, double &, double &);                
template void diffusebounceBackNodeG34Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,gridBCC3D<1, 1, double> & ,int ,int , double , double , double , double &, double &, double &);           
template void diffusebounceBackNodeG56Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,gridBCC3D<1, 1, double> & ,int ,int , double , double , double , double &, double &, double &);          
template void diffusebounceBackNodeG910Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,gridBCC3D<1, 1, double> & ,int ,int , double , double , double , double &, double &, double &);            
template void diffusebounceBackCellG12Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,gridBCC3D<1, 1, double> & ,int ,int , double , double , double , double &, double &, double &);            
template void diffusebounceBackCellG34Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,gridBCC3D<1, 1, double> & ,int ,int , double , double , double , double &, double &, double &);            
template void diffusebounceBackCellG56Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,gridBCC3D<1, 1, double> & ,int ,int , double , double , double , double &, double &, double &);            
template void diffusebounceBackCellG910Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,gridBCC3D<1, 1, double> & ,int ,int , double , double , double , double &, double &, double &);      
template void diffusebounceBackNodeCellG78Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,gridBCC3D<1, 1, double> & ,int ,int , double , double , double , double &, double &, double &); 
template void diffusebounceBackSolid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,gridBCC3D<1, 1, double> & ,int , int , double , double , double , double &, double &, double &) ;           


