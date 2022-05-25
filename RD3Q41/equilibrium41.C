#include"equilibrium41.h"
#include"library.C"


template<typename dataType1>
void getFEq(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta)
{
    dataType1 delThetaBy2[VECT_LENGTH]       __attribute__ ((aligned(32)));
    dataType1 delSqTheta0p125[VECT_LENGTH]   __attribute__ ((aligned(32)));
    dataType1 twoMinusU2,order1,order2,f0V,f0V1,uDotC,uDotCSq;
    dataType1 oneByTheta;

    for(int i=0;i<VECT_LENGTH;i++)
    {    
        delThetaBy2[i]     =  0.5*(theta[i] - lbModel.theta0)*lbModel.oneByTheta0;
        delSqTheta0p125[i] =  0.5* delThetaBy2[i] *delThetaBy2[i];   
        
        oneByTheta  = 1.0/(theta[i]);
        uX[i]      *= oneByTheta;
        uY[i]      *= oneByTheta;
        uZ[i]      *= oneByTheta;
        rho[i] *= 0.5;
    }
    
    for(int i=0;i<VECT_LENGTH;i++)
    { 
        twoMinusU2 = 2.0-(uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i])*theta[i];
        
        // CENTER
        order1       =  delThetaBy2[i]      *lbModel.yMinus3_CENTER;     
        order2       =  delSqTheta0p125[i]  *lbModel.yTenMinusYSqrPlus15_CENTER;      
        f0V          =  rho[i]*lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO]*(1.0 + order1+ order2); 
        lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO]=  f0V*twoMinusU2;    
        // G1 and G2
        // ZZP1 and ZZM1
        order1       =  delThetaBy2[i]     *lbModel.yMinus3_SC1;        
        order2       =  delSqTheta0p125[i] *lbModel.yTenMinusYSqrPlus15_SC1;         
        f0V          =  rho[i]*lbModel.wt[lbModel.DV_ZERO_ZERO_P1]*(1.0 + order1+ order2);        
        uDotC        =  uZ[i];
        lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC)); 
        lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
        
        order1       =  delThetaBy2[i]     *lbModel.yMinus3_SC2;        
        order2       =  delSqTheta0p125[i] *lbModel.yTenMinusYSqrPlus15_SC2;         
        f0V          =  rho[i]*lbModel.wt[lbModel.DV_ZERO_ZERO_P2]*(1.0 + order1+ order2);   
        // ZZP2 and ZZM2
        uDotC        =  2.0*uZ[i];
        lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P2] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC)); 
        lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M2] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC)); 
        // ZP2Z and ZM2Z
        uDotC        =  2.0*uY[i];
        lbModel.fTemp1[i][lbModel.G1_DV_ZERO_P2_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC)); 
        lbModel.fTemp2[i][lbModel.G2_DV_ZERO_M2_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));         
        // P2ZZ and M2ZZ
        uDotC        =  2.0*uX[i];
        lbModel.fTemp1[i][lbModel.G1_DV_P2_ZERO_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC)); 
        lbModel.fTemp2[i][lbModel.G2_DV_M2_ZERO_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));   
    }
    
    for(int i=0;i<VECT_LENGTH;i++)
    { 
        twoMinusU2 = 2.0-(uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i])*theta[i];
        
        // G3 and G4
        order1    =  delThetaBy2[i]    *lbModel.yMinus3_FCC;        
        order2    =  delSqTheta0p125[i]*lbModel.yTenMinusYSqrPlus15_FCC;         
        f0V       =  rho[i]*lbModel.wt[lbModel.DV_ZERO_P1_P1]*(1.0 + order1+ order2);        
        //ZPP and ZMM         
        uDotC     =  (uY[i] + uZ[i]);        
        lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));        
        lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
        //ZMP and ZPM                                                                                            
        uDotC     =  (uZ[i] - uY[i]);                                                                       
        lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));        
        lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
        //PZP and MZM                                                                                            
        uDotC     =  (uX[i] + uZ[i]);                                                                       
        lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));        
        lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
        //MZP and PZM                                                                                            
        uDotC     =  (uZ[i] - uX[i]);                                                                   
        lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));    
        lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
        
    }
    
    //G5 and G6
    for(int i=0;i<VECT_LENGTH;i++)
    { 
        twoMinusU2 = 2.0-(uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i])*theta[i];
        
        order1    =  delThetaBy2[i]    *lbModel.yMinus3_FCC;
        order2    =  delSqTheta0p125[i]*lbModel.yTenMinusYSqrPlus15_FCC; 
        f0V       =  rho[i]*lbModel.wt[lbModel.DV_P1_P1_ZERO]*(1.0 + order1+ order2);
        //PPZ and MMZ
        uDotC     =  uX[i] + uY[i];
        lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
        lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
        //MPZ and PMZ                                                                               
        uDotC     =  (uY[i] - uX[i]);                                                             
        lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
        lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
        //ZPZ and ZMZ
        order1    =  delThetaBy2[i]    *lbModel.yMinus3_SC1;
        order2    =  delSqTheta0p125[i]*lbModel.yTenMinusYSqrPlus15_SC1; 
        f0V       =  rho[i]*lbModel.wt[lbModel.DV_ZERO_P1_ZERO]*(1.0 + order1+ order2);    
        uDotC     =  uY[i];
        lbModel.fTemp5[i][lbModel.G5_DV_ZERO_P1_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
        lbModel.fTemp6[i][lbModel.G6_DV_ZERO_M1_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
        //PZZ aqnd MZZ                                                                                  
        uDotC     =   uX[i];                                                                          
        lbModel.fTemp5[i][lbModel.G5_DV_P1_ZERO_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
        lbModel.fTemp6[i][lbModel.G6_DV_M1_ZERO_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
    }      
    
    // G7,G8,G9,G10
    for(int i=0;i<VECT_LENGTH;i++)
    { 
        twoMinusU2 = 2.0-(uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i])*theta[i];
        
        order1    =  1.0 + delThetaBy2[i]*lbModel.yMinus3_BCC +  delSqTheta0p125[i]*lbModel.yTenMinusYSqrPlus15_BCC; 
        f0V       =  rho[i]*lbModel.wt[lbModel.DV_P_P_P]*(order1);
        order2    =  1.0 + delThetaBy2[i]*lbModel.yMinus3_BCC1 + delSqTheta0p125[i]*lbModel.yTenMinusYSqrPlus15_BCC1;         
        f0V1      =  rho[i]*lbModel.wt[lbModel.DV_P1_P1_P1]*(order2);        
        //PPP1 and MMM1        
        uDotC     = (uX[i] + uY[i] + uZ[i]);      
        uDotCSq   = uDotC*uDotC;        
        lbModel.fTemp9 [i][lbModel.G9_DV_P1_P1_P1]   = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);
        lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]  = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);   
        //PPP and MMM
        //uDotC     =   uDotC*0.5 ;
        lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]   = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);
        lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]   = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);
        
        //MPP1 and PMM1                                                                             
        uDotC     = (uY[i] - uX[i] + uZ[i]); 
        uDotCSq   = uDotC*uDotC;
        lbModel.fTemp9 [i][lbModel.G9_DV_M1_P1_P1]  = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);
        lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);  
        //MPP and PMM                                                                           
        // uDotC     =  uDotC*0.5 ;                                           
        lbModel.fTemp7[i][lbModel.G7_DV_M_P_P] = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);     
        lbModel.fTemp8[i][lbModel.G8_DV_P_M_M] = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);       
        
        //MMP1 and PPM 1                                                                                     
        uDotC     = (uZ[i] - uX[i] - uY[i]);           
        uDotCSq   = uDotC*uDotC;       
        lbModel.fTemp9 [i][lbModel.G9_DV_M1_M1_P1] = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);       
        lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);        
        //MMP and PPM                                                                                  
        // uDotC     =   uDotC*0.5 ;                                                  
        lbModel.fTemp7[i][lbModel.G7_DV_M_M_P] = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);       
        lbModel.fTemp8[i][lbModel.G8_DV_P_P_M] = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);       
        
        //PMP and MPM                                                                                      
        uDotC     = (uX[i] - uY[i] + uZ[i]);        
        uDotCSq   = uDotC*uDotC;       
        lbModel.fTemp9 [i][lbModel.G9_DV_P1_M1_P1]  = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);       
        lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);         
        //PMP and MPM                                                                                  
        // uDotC     =  uDotC *0.5 ;                                                  
        lbModel.fTemp7[i][lbModel.G7_DV_P_M_P] = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);       
        lbModel.fTemp8[i][lbModel.G8_DV_M_P_M] = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);  
        
    }
    
    for(int i=0;i<VECT_LENGTH;i++)
    {    
        uX[i]      *= theta[i];
        uY[i]      *= theta[i];
        uZ[i]      *= theta[i];
        rho[i]     *= 2.0;   
    }
    
    
}


template<typename dataType1>
void getFEqSIMD(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, int countFeq)
{
    SIMD_REG _RHO, _UX, _UY, _UZ, _THETA;              //4
    SIMD_REG _temp1,_temp2,_temp3,_temp4;              //4
    SIMD_REG _temp5,_temp6,_temp7,_temp8;              //4       
    SIMD_REG _UDOTC,_twoReg,_oneReg;                   //3   = 16 Registers in Total    
    
    dataType1 *pointer1;     
    dataType1 tempArray[4];

    _oneReg = SET1_PD(1.0);     
    _twoReg = SET1_PD(2.0);
    _temp8 = SET1_PD(0.5);
    _temp7 = SET1_PD(lbModel.theta0);
    _temp6 = SET1_PD(lbModel.oneByTheta0);

	if(countFeq == 0)
	{
		std::cout<<""<<std::endl;
	}

    
    pointer1 = &rho[0];
    _RHO     = LOAD_PD(pointer1);  //rho      
    pointer1 = &uX[0];
    _UX = LOAD_PD(pointer1);       //uX     
    pointer1 = &uY[0];
    _UY = LOAD_PD(pointer1);       //uY
    pointer1 = &uZ[0];
    _UZ = LOAD_PD(pointer1);       //uZ    
    pointer1 = &theta[0];
    _THETA = LOAD_PD(pointer1);     //theta 

    _temp1 = SUB_PD(_THETA,_temp7);
    _temp1 = MUL_PD(_temp1,_temp6);
    _temp1 = MUL_PD(_temp1,_temp8); // delThetaBy2
    _temp2 = MUL_PD(_temp1,_temp1);
    _temp2 = MUL_PD(_temp2,_temp8); // delSqTheta0p125
    _temp7 = DIV_PD(_oneReg,_THETA);
    
    _UX    = MUL_PD(_temp7,_UX);
    _UY    = MUL_PD(_temp7,_UY);
    _UZ    = MUL_PD(_temp7,_UZ);
    _RHO   = MUL_PD(_temp8,_RHO);
    
   // CENTER
   _temp3  = MUL_PD(_UX,_UX);                                //uX^2
   _temp4  = MUL_PD(_UY,_UY);                                //uY^2 
   _temp3  = ADD_PD(_temp3,_temp4);                          //uX^2+uY^2 
   _temp5  = MUL_PD(_UZ,_UZ);                                //uZ^2 
   _temp3  = ADD_PD(_temp3,_temp5);                          //uX^2+uY^2+uZ^2
   _temp3  = MUL_PD(_temp3,_THETA);                          //(uX^2+uY^2+uZ^2)*theta
   _temp3  = SUB_PD(_twoReg,_temp3);  //twoMinusU2           // twoMinusU2 =  2.0 - (uX^2+uY^2+uZ^2)*theta
 
   
    // _temp1,_temp2,_temp3 will be used till the end of the function 
   
   _temp4  = SET1_PD(lbModel.yMinus3_CENTER);                
   _temp4  = MUL_PD(_temp1,_temp4);                          // order1
   _temp5  = SET1_PD(lbModel.yTenMinusYSqrPlus15_CENTER);    
   _temp5  = MUL_PD(_temp2,_temp5);                          //order2
   _temp4  = ADD_PD(_temp4,_temp5);                          //order1+ order2
   _temp4  = ADD_PD(_oneReg,_temp4);                         // 1.0 + order1+ order2
   _temp5  = SET1_PD(lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO]); 
   _temp5  = MUL_PD(_RHO,_temp5);                            //rho[i]*lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO]
   _temp5  = MUL_PD(_temp4,_temp5);                          // f0V
   _temp5  = MUL_PD(_temp3,_temp5);                          // f0V*twoMinusU2

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp5); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO]=  tempArray[i];    

//    pointer1  = &lbModel.fTemp0;
//    STORE_PD(pointer1,_temp5); 
//    for(int i=0;i<VECT_LENGTH;i++)
//     lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO]=  tempArray[i];    

   _temp4  = SET1_PD(lbModel.yMinus3_SC1);
   _temp4  = MUL_PD(_temp1,_temp4);                          // order1 
   _temp5  = SET1_PD(lbModel.yTenMinusYSqrPlus15_SC1);
   _temp5  = MUL_PD(_temp2,_temp5);                          // order2 
   _temp4  = ADD_PD(_temp4,_temp5);                          //order1+ order2
   _temp4  = ADD_PD(_oneReg,_temp4);                         // 1.0 + order1+ order2
   _temp5  = SET1_PD(lbModel.wt[lbModel.DV_ZERO_ZERO_P1]); 
   _temp5  = MUL_PD(_RHO,_temp5);                            //rho[i]*lbModel.wt[lbModel.DV_ZERO_ZERO_P1]
   _temp5  = MUL_PD(_temp4,_temp5);                          // f0V
   _temp6  = ADD_PD(_twoReg,_UZ);                            // 2.0 + uDotC       
   _temp6  = MUL_PD(_temp6,_UZ);                             // uDotC (2.0 + uDotC)  
   _temp6  = ADD_PD(_temp3,_temp6);                          // twoMinusU2 + uDotC (2.0 + uDotC)  
   _temp6  = MUL_PD(_temp5,_temp6);                          // f0V * (twoMinusU2 + uDotC (2.0 + uDotC))
   _temp7  = SUB_PD(_twoReg,_UZ);                            // 2.0 - uDotC
   _temp7  = MUL_PD(_temp7,_UZ);                             // uDotC (2.0 - uDotC)
   _temp7  = SUB_PD(_temp3,_temp7);                          // twoMinusU2 - uDotC (2.0 + uDotC)
   _temp7  = MUL_PD(_temp5,_temp7);                          // f0V * (twoMinusU2 - uDotC (2.0 - uDotC))

	//std::cout<<"getfeqSIMD3!"<<std::endl;
   
   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp6); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P1]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp7); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M1]=  tempArray[i];    

   _temp6  = ADD_PD(_twoReg,_UY);                            // 2.0 + uDotC       
   _temp6  = MUL_PD(_temp6,_UY);                             // uDotC (2.0 + uDotC)  
   _temp6  = ADD_PD(_temp3,_temp6);                          // twoMinusU2 + uDotC (2.0 + uDotC)  
   _temp6  = MUL_PD(_temp5,_temp6);                          // f0V * (twoMinusU2 + uDotC (2.0 + uDotC))
   _temp7  = SUB_PD(_twoReg,_UY);                            // 2.0 - uDotC
   _temp7  = MUL_PD(_temp7,_UY);                             // uDotC (2.0 - uDotC)
   _temp7  = SUB_PD(_temp3,_temp7);                          // twoMinusU2 - uDotC (2.0 + uDotC)
   _temp7  = MUL_PD(_temp5,_temp7);                          // f0V * (twoMinusU2 - uDotC (2.0 - uDotC))

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp6); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp5[i][lbModel.G5_DV_ZERO_P1_ZERO]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp7); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp6[i][lbModel.G6_DV_ZERO_M1_ZERO]=  tempArray[i];     
   
   _temp6  = ADD_PD(_twoReg,_UX);                            // 2.0 + uDotC       
   _temp6  = MUL_PD(_temp6,_UX);                             // uDotC (2.0 + uDotC)  
   _temp6  = ADD_PD(_temp3,_temp6);                          // twoMinusU2 + uDotC (2.0 + uDotC)  
   _temp6  = MUL_PD(_temp5,_temp6);                          // f0V * (twoMinusU2 + uDotC (2.0 + uDotC))
   _temp7  = SUB_PD(_twoReg,_UX);                            // 2.0 - uDotC
   _temp7  = MUL_PD(_temp7,_UX);                             // uDotC (2.0 - uDotC)
   _temp7  = SUB_PD(_temp3,_temp7);                          // twoMinusU2 - uDotC (2.0 + uDotC)
   _temp7  = MUL_PD(_temp5,_temp7);                          // f0V * (twoMinusU2 - uDotC (2.0 - uDotC))

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp6); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp5[i][lbModel.G5_DV_P1_ZERO_ZERO]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp7); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp6[i][lbModel.G6_DV_M1_ZERO_ZERO]=  tempArray[i];     

   
   _temp4  = SET1_PD(lbModel.yMinus3_SC2);
   _temp4  = MUL_PD(_temp1,_temp4);                          // order1 
   _temp5  = SET1_PD(lbModel.yTenMinusYSqrPlus15_SC2);
   _temp5  = MUL_PD(_temp2,_temp5);                          // order2 
   _temp4  = ADD_PD(_temp4,_temp5);                          //order1+ order2
   _temp4  = ADD_PD(_oneReg,_temp4);                         // 1.0 + order1+ order2
   _temp5  = SET1_PD(lbModel.wt[lbModel.DV_ZERO_ZERO_P2]); 
   _temp5  = MUL_PD(_RHO,_temp5);                            //rho[i]*lbModel.wt[lbModel.DV_ZERO_ZERO_P1]   
   _temp5  = MUL_PD(_temp4,_temp5);                          // f0V
   _UDOTC  = MUL_PD(_twoReg,_UZ);   
   _temp6  = ADD_PD(_twoReg,_UDOTC);                            // 2.0 + uDotC       
   _temp6  = MUL_PD(_temp6,_UDOTC);                             // uDotC (2.0 + uDotC)  
   _temp6  = ADD_PD(_temp3,_temp6);                          // twoMinusU2 + uDotC (2.0 + uDotC)  
   _temp6  = MUL_PD(_temp5,_temp6);                          // f0V * (twoMinusU2 + uDotC (2.0 + uDotC))
   _temp7  = SUB_PD(_twoReg,_UDOTC);                            // 2.0 - uDotC
   _temp7  = MUL_PD(_temp7,_UDOTC);                             // uDotC (2.0 - uDotC)
   _temp7  = SUB_PD(_temp3,_temp7);                          // twoMinusU2 - uDotC (2.0 + uDotC)
   _temp7  = MUL_PD(_temp5,_temp7);                          // f0V * (twoMinusU2 - uDotC (2.0 - uDotC))

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp6); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P2]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp7); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M2]=  tempArray[i];    

   _UDOTC  = MUL_PD(_twoReg,_UY);   
   _temp6  = ADD_PD(_twoReg,_UDOTC);                            // 2.0 + uDotC       
   _temp6  = MUL_PD(_temp6,_UDOTC);                             // uDotC (2.0 + uDotC)  
   _temp6  = ADD_PD(_temp3,_temp6);                          // twoMinusU2 + uDotC (2.0 + uDotC)  
   _temp6  = MUL_PD(_temp5,_temp6);                          // f0V * (twoMinusU2 + uDotC (2.0 + uDotC))
   _temp7  = SUB_PD(_twoReg,_UDOTC);                            // 2.0 - uDotC
   _temp7  = MUL_PD(_temp7,_UDOTC);                             // uDotC (2.0 - uDotC)
   _temp7  = SUB_PD(_temp3,_temp7);                          // twoMinusU2 - uDotC (2.0 + uDotC)
   _temp7  = MUL_PD(_temp5,_temp7);                          // f0V * (twoMinusU2 - uDotC (2.0 - uDotC))

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp6); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp1[i][lbModel.G1_DV_ZERO_P2_ZERO]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp7); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp2[i][lbModel.G2_DV_ZERO_M2_ZERO]=  tempArray[i];    

   _UDOTC  = MUL_PD(_twoReg,_UX);   
   _temp6  = ADD_PD(_twoReg,_UDOTC);                            // 2.0 + uDotC       
   _temp6  = MUL_PD(_temp6,_UDOTC);                             // uDotC (2.0 + uDotC)  
   _temp6  = ADD_PD(_temp3,_temp6);                          // twoMinusU2 + uDotC (2.0 + uDotC)  
   _temp6  = MUL_PD(_temp5,_temp6);                          // f0V * (twoMinusU2 + uDotC (2.0 + uDotC))
   _temp7  = SUB_PD(_twoReg,_UDOTC);                            // 2.0 - uDotC
   _temp7  = MUL_PD(_temp7,_UDOTC);                             // uDotC (2.0 - uDotC)
   _temp7  = SUB_PD(_temp3,_temp7);                          // twoMinusU2 - uDotC (2.0 + uDotC)
   _temp7  = MUL_PD(_temp5,_temp7);                          // f0V * (twoMinusU2 - uDotC (2.0 - uDotC))

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp6); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp1[i][lbModel.G1_DV_P2_ZERO_ZERO]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp7); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp2[i][lbModel.G2_DV_M2_ZERO_ZERO]=  tempArray[i];    

   
   
   
   
   
   
   
   
   
   
   
   
   _temp4  = SET1_PD(lbModel.yMinus3_FCC);
   _temp4  = MUL_PD(_temp1,_temp4);                          // order1 
   _temp5  = SET1_PD(lbModel.yTenMinusYSqrPlus15_FCC);
   _temp5  = MUL_PD(_temp2,_temp5);                          // order2 
   _temp4  = ADD_PD(_temp4,_temp5);                          //order1+ order2
   _temp4  = ADD_PD(_oneReg,_temp4);                         // 1.0 + order1+ order2
   _temp5  = SET1_PD(lbModel.wt[lbModel.DV_ZERO_P1_P1]); 
   _temp5  = MUL_PD(_RHO,_temp5);                            //rho[i]*lbModel.wt[lbModel.DV_ZERO_ZERO_P1]   
   _temp5  = MUL_PD(_temp4,_temp5);                          // f0V

   _UDOTC  = ADD_PD(_UY,_UZ);   
   _temp6  = ADD_PD(_twoReg,_UDOTC);                            // 2.0 + uDotC       
   _temp6  = MUL_PD(_temp6,_UDOTC);                             // uDotC (2.0 + uDotC)  
   _temp6  = ADD_PD(_temp3,_temp6);                          // twoMinusU2 + uDotC (2.0 + uDotC)  
   _temp6  = MUL_PD(_temp5,_temp6);                          // f0V * (twoMinusU2 + uDotC (2.0 + uDotC))
   _temp7  = SUB_PD(_twoReg,_UDOTC);                            // 2.0 - uDotC
   _temp7  = MUL_PD(_temp7,_UDOTC);                             // uDotC (2.0 - uDotC)
   _temp7  = SUB_PD(_temp3,_temp7);                          // twoMinusU2 - uDotC (2.0 + uDotC)
   _temp7  = MUL_PD(_temp5,_temp7);                          // f0V * (twoMinusU2 - uDotC (2.0 - uDotC))

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp6); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp7); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1]=  tempArray[i];       

   _UDOTC  = SUB_PD(_UZ,_UY);   
   _temp6  = ADD_PD(_twoReg,_UDOTC);                            // 2.0 + uDotC       
   _temp6  = MUL_PD(_temp6,_UDOTC);                             // uDotC (2.0 + uDotC)  
   _temp6  = ADD_PD(_temp3,_temp6);                          // twoMinusU2 + uDotC (2.0 + uDotC)  
   _temp6  = MUL_PD(_temp5,_temp6);                          // f0V * (twoMinusU2 + uDotC (2.0 + uDotC))
   _temp7  = SUB_PD(_twoReg,_UDOTC);                            // 2.0 - uDotC
   _temp7  = MUL_PD(_temp7,_UDOTC);                             // uDotC (2.0 - uDotC)
   _temp7  = SUB_PD(_temp3,_temp7);                          // twoMinusU2 - uDotC (2.0 + uDotC)
   _temp7  = MUL_PD(_temp5,_temp7);                          // f0V * (twoMinusU2 - uDotC (2.0 - uDotC))

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp6); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp7); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1]=  tempArray[i];       

   _UDOTC  = ADD_PD(_UX,_UZ);   
   _temp6  = ADD_PD(_twoReg,_UDOTC);                            // 2.0 + uDotC       
   _temp6  = MUL_PD(_temp6,_UDOTC);                             // uDotC (2.0 + uDotC)  
   _temp6  = ADD_PD(_temp3,_temp6);                          // twoMinusU2 + uDotC (2.0 + uDotC)  
   _temp6  = MUL_PD(_temp5,_temp6);                          // f0V * (twoMinusU2 + uDotC (2.0 + uDotC))
   _temp7  = SUB_PD(_twoReg,_UDOTC);                            // 2.0 - uDotC
   _temp7  = MUL_PD(_temp7,_UDOTC);                             // uDotC (2.0 - uDotC)
   _temp7  = SUB_PD(_temp3,_temp7);                          // twoMinusU2 - uDotC (2.0 + uDotC)
   _temp7  = MUL_PD(_temp5,_temp7);                          // f0V * (twoMinusU2 - uDotC (2.0 - uDotC))

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp6); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp7); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1]=  tempArray[i];       

   _UDOTC  = SUB_PD(_UZ,_UX);   
   _temp6  = ADD_PD(_twoReg,_UDOTC);                            // 2.0 + uDotC       
   _temp6  = MUL_PD(_temp6,_UDOTC);                             // uDotC (2.0 + uDotC)  
   _temp6  = ADD_PD(_temp3,_temp6);                          // twoMinusU2 + uDotC (2.0 + uDotC)  
   _temp6  = MUL_PD(_temp5,_temp6);                          // f0V * (twoMinusU2 + uDotC (2.0 + uDotC))
   _temp7  = SUB_PD(_twoReg,_UDOTC);                            // 2.0 - uDotC
   _temp7  = MUL_PD(_temp7,_UDOTC);                             // uDotC (2.0 - uDotC)
   _temp7  = SUB_PD(_temp3,_temp7);                          // twoMinusU2 - uDotC (2.0 + uDotC)
   _temp7  = MUL_PD(_temp5,_temp7);                          // f0V * (twoMinusU2 - uDotC (2.0 - uDotC))

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp6); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp7); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1]=  tempArray[i];       

   _UDOTC  = ADD_PD(_UX,_UY);   
   _temp6  = ADD_PD(_twoReg,_UDOTC);                            // 2.0 + uDotC       
   _temp6  = MUL_PD(_temp6,_UDOTC);                             // uDotC (2.0 + uDotC)  
   _temp6  = ADD_PD(_temp3,_temp6);                          // twoMinusU2 + uDotC (2.0 + uDotC)  
   _temp6  = MUL_PD(_temp5,_temp6);                          // f0V * (twoMinusU2 + uDotC (2.0 + uDotC))
   _temp7  = SUB_PD(_twoReg,_UDOTC);                            // 2.0 - uDotC
   _temp7  = MUL_PD(_temp7,_UDOTC);                             // uDotC (2.0 - uDotC)
   _temp7  = SUB_PD(_temp3,_temp7);                          // twoMinusU2 - uDotC (2.0 + uDotC)
   _temp7  = MUL_PD(_temp5,_temp7);                          // f0V * (twoMinusU2 - uDotC (2.0 - uDotC))

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp6); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp7); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO]=  tempArray[i];       

   _UDOTC  = SUB_PD(_UY,_UX);   
   _temp6  = ADD_PD(_twoReg,_UDOTC);                            // 2.0 + uDotC       
   _temp6  = MUL_PD(_temp6,_UDOTC);                             // uDotC (2.0 + uDotC)  
   _temp6  = ADD_PD(_temp3,_temp6);                          // twoMinusU2 + uDotC (2.0 + uDotC)  
   _temp6  = MUL_PD(_temp5,_temp6);                          // f0V * (twoMinusU2 + uDotC (2.0 + uDotC))
   _temp7  = SUB_PD(_twoReg,_UDOTC);                            // 2.0 - uDotC
   _temp7  = MUL_PD(_temp7,_UDOTC);                             // uDotC (2.0 - uDotC)
   _temp7  = SUB_PD(_temp3,_temp7);                          // twoMinusU2 - uDotC (2.0 + uDotC)
   _temp7  = MUL_PD(_temp5,_temp7);                          // f0V * (twoMinusU2 - uDotC (2.0 - uDotC))

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp6); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp7); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO]=  tempArray[i];       
   
    
 
   
   
   
   
   
   
   _temp4  = SET1_PD(lbModel.yMinus3_BCC);
   _temp4  = MUL_PD(_temp1,_temp4);                          // order1 
   _temp5  = SET1_PD(lbModel.yTenMinusYSqrPlus15_BCC);
   _temp5  = MUL_PD(_temp2,_temp5);                          // order2 
   _temp4  = ADD_PD(_temp4,_temp5);                          //order1+ order2
   _temp4  = ADD_PD(_oneReg,_temp4);                         // 1.0 + order1+ order2
   _temp5  = SET1_PD(lbModel.wt[lbModel.DV_P_P_P]); 
   _temp5  = MUL_PD(_RHO,_temp5);                            //rho[i]*lbModel.wt[lbModel.DV_ZERO_ZERO_P1]   
   _temp5  = MUL_PD(_temp4,_temp5);                          // f0V   
   
 
   _temp6  = SET1_PD(lbModel.yMinus3_BCC1);
   _temp6  = MUL_PD(_temp1,_temp6);                          // order1 
   _temp7  = SET1_PD(lbModel.yTenMinusYSqrPlus15_BCC1);
   _temp7  = MUL_PD(_temp2,_temp7);                          // order2 
   _temp6  = ADD_PD(_temp6,_temp7);                          //order1+ order2
   _temp6  = ADD_PD(_oneReg,_temp6);                         // 1.0 + order1+ order2
   _temp7  = SET1_PD(lbModel.wt[lbModel.DV_P1_P1_P1]); 
   _temp7  = MUL_PD(_RHO,_temp7);                            //rho[i]*lbModel.wt[lbModel.DV_ZERO_ZERO_P1]   
   _temp7  = MUL_PD(_temp6,_temp7);                          // f0V1  
   
   _UDOTC  = ADD_PD(_UX,_UY);   
   _UDOTC  = ADD_PD(_UDOTC,_UZ);   
   _temp4  = MUL_PD(_UDOTC,_UDOTC);                          // uDotCSq
   
   _temp1 = MUL_PD(_twoReg,_UDOTC);
   _temp1 = ADD_PD(_temp3,_temp1);
   _temp1 = ADD_PD(_temp4,_temp1);
   _temp1 = MUL_PD(_temp7,_temp1);


   _temp2 = MUL_PD(_twoReg,_UDOTC);
   _temp2 = SUB_PD(_temp3,_temp2);
   _temp2 = ADD_PD(_temp4,_temp2);
   _temp2 = MUL_PD(_temp7,_temp2);
   
   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp1); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp2); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]=  tempArray[i];       

   _temp1 = ADD_PD(_temp3,_UDOTC);
   _temp2 = SUB_PD(_temp3,_UDOTC);
   _temp6 = SET1_PD(0.25);
   _temp6 = MUL_PD(_temp6,_temp4);
   _temp1 = ADD_PD(_temp6,_temp1);
   _temp2 = ADD_PD(_temp6,_temp2);
   _temp1 = MUL_PD(_temp5,_temp1);
   _temp2 = MUL_PD(_temp5,_temp2);   
   
   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp1); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp2); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]=  tempArray[i];      


   _UDOTC  = SUB_PD(_UY,_UX);   
   _UDOTC  = ADD_PD(_UDOTC,_UZ);   
   _temp4  = MUL_PD(_UDOTC,_UDOTC);                          // uDotCSq
   
   _temp1 = MUL_PD(_twoReg,_UDOTC);
   _temp1 = ADD_PD(_temp3,_temp1);
   _temp1 = ADD_PD(_temp4,_temp1);
   _temp1 = MUL_PD(_temp7,_temp1);


   _temp2 = MUL_PD(_twoReg,_UDOTC);
   _temp2 = SUB_PD(_temp3,_temp2);
   _temp2 = ADD_PD(_temp4,_temp2);
   _temp2 = MUL_PD(_temp7,_temp2);
   
   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp1); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp2); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]=  tempArray[i];       

   _temp1 = ADD_PD(_temp3,_UDOTC);
   _temp2 = SUB_PD(_temp3,_UDOTC);
   _temp6 = SET1_PD(0.25);
   _temp6 = MUL_PD(_temp6,_temp4);
   _temp1 = ADD_PD(_temp6,_temp1);
   _temp2 = ADD_PD(_temp6,_temp2);
   _temp1 = MUL_PD(_temp5,_temp1);
   _temp2 = MUL_PD(_temp5,_temp2);   
   
   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp1); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp2); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]=  tempArray[i];      
   

   _UDOTC  = SUB_PD(_UZ,_UX);   
   _UDOTC  = SUB_PD(_UDOTC,_UY);   
   _temp4  = MUL_PD(_UDOTC,_UDOTC);                          // uDotCSq
   
   _temp1 = MUL_PD(_twoReg,_UDOTC);
   _temp1 = ADD_PD(_temp3,_temp1);
   _temp1 = ADD_PD(_temp4,_temp1);
   _temp1 = MUL_PD(_temp7,_temp1);


   _temp2 = MUL_PD(_twoReg,_UDOTC);
   _temp2 = SUB_PD(_temp3,_temp2);
   _temp2 = ADD_PD(_temp4,_temp2);
   _temp2 = MUL_PD(_temp7,_temp2);
   
   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp1); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp2); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]=  tempArray[i];       

   _temp1 = ADD_PD(_temp3,_UDOTC);
   _temp2 = SUB_PD(_temp3,_UDOTC);
   _temp6 = SET1_PD(0.25);
   _temp6 = MUL_PD(_temp6,_temp4);
   _temp1 = ADD_PD(_temp6,_temp1);
   _temp2 = ADD_PD(_temp6,_temp2);
   _temp1 = MUL_PD(_temp5,_temp1);
   _temp2 = MUL_PD(_temp5,_temp2);   
   
   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp1); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp2); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]=  tempArray[i];      


   _UDOTC  = SUB_PD(_UX,_UY);   
   _UDOTC  = ADD_PD(_UDOTC,_UZ);   
   _temp4  = MUL_PD(_UDOTC,_UDOTC);                          // uDotCSq
   
   _temp1 = MUL_PD(_twoReg,_UDOTC);
   _temp1 = ADD_PD(_temp3,_temp1);
   _temp1 = ADD_PD(_temp4,_temp1);
   _temp1 = MUL_PD(_temp7,_temp1);


   _temp2 = MUL_PD(_twoReg,_UDOTC);
   _temp2 = SUB_PD(_temp3,_temp2);
   _temp2 = ADD_PD(_temp4,_temp2);
   _temp2 = MUL_PD(_temp7,_temp2);
   
   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp1); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp2); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]=  tempArray[i];       

   _temp1 = ADD_PD(_temp3,_UDOTC);
   _temp2 = SUB_PD(_temp3,_UDOTC);
   _temp6 = SET1_PD(0.25);
   _temp6 = MUL_PD(_temp6,_temp4);
   _temp1 = ADD_PD(_temp6,_temp1);
   _temp2 = ADD_PD(_temp6,_temp2);
   _temp1 = MUL_PD(_temp5,_temp1);
   _temp2 = MUL_PD(_temp5,_temp2);   
   
   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp1); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]=  tempArray[i];    

   pointer1  = &tempArray[0];
   STORE_PD(pointer1 ,_temp2); 
   for(int i=0;i<VECT_LENGTH;i++)
    lbModel.fTemp8[i][lbModel.G8_DV_M_P_M]=  tempArray[i];      

   _UX = MUL_PD(_UX,_THETA);
   _UY = MUL_PD(_UY,_THETA);
   _UZ = MUL_PD(_UZ,_THETA);   
   _RHO = MUL_PD(_twoReg,_RHO);

   STORE_PD(rho  ,_RHO  );  
   STORE_PD(uX   ,_UX   ); 
   STORE_PD(uY   ,_UY   ); 
   STORE_PD(uZ   ,_UZ   ); 
}




    template<typename dataType1>
    void getIterativeFEq(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta)
    {
        dataType1 uSq, E;

        dataType1 rho1 ,mxy ,mxz ,myz ,bx ,by ,bz ;
        dataType1 mxx ,myy ,mzz ,Vx ,Vy ,Vz ,qx ,qy ,qz ,qm ,q2m ;
        dataType1 dBx ,dBy ,dBz ,dga , del_theta ,ga ;

        dataType1 dot;
        
        for(int i=0;i<VECT_LENGTH;i++)
        {
            uSq = uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i];
            E   = uSq + 3.0*theta[i];
            bx  = uX[i]/theta[i];
            by  = uY[i]/theta[i];
            bz  = uZ[i]/theta[i];
            del_theta = (theta[i]-lbModel.theta0)*lbModel.oneByTheta0 ;
            ga  = 0.5*lbModel.oneByTheta0*del_theta*( 1.0-del_theta*del_theta ) / (1.0+del_theta); 
        
            lbModel.fTemp[lbModel.DV_ZERO_ZERO_ZERO]    = lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_ZERO] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_ZERO] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_ZERO] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_ZERO]) ;
            lbModel.fTemp[lbModel.DV_ZERO_ZERO_P1  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_P1  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_P1  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_P1  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_P1  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_P1  ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_ZERO_P2  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_P2  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_P2  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_P2  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_P2  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_P2  ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_P2_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_P2_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P2_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_P2_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_P2_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_P2_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_P2_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_P2_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_P2_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_P2_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_P2_ZERO_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_ZERO_M1  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_M1  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_M1  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_M1  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_M1  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_M1  ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_ZERO_M2  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_M2  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_M2  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_M2  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_M2  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_M2  ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_M2_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_M2_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M2_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_M2_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_M2_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_M2_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_M2_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_M2_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_M2_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_M2_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_M2_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_M2_ZERO_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_P1_P1    ]    = lbModel.wt[lbModel.DV_ZERO_P1_P1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P1_P1    ] + by*lbModel.cy[lbModel.DV_ZERO_P1_P1    ] + bz*lbModel.cz[lbModel.DV_ZERO_P1_P1    ] + ga*lbModel.cc[lbModel.DV_ZERO_P1_P1    ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_M1_P1    ]    = lbModel.wt[lbModel.DV_ZERO_M1_P1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M1_P1    ] + by*lbModel.cy[lbModel.DV_ZERO_M1_P1    ] + bz*lbModel.cz[lbModel.DV_ZERO_M1_P1    ] + ga*lbModel.cc[lbModel.DV_ZERO_M1_P1    ]) ;
            lbModel.fTemp[lbModel.DV_P1_ZERO_P1    ]    = lbModel.wt[lbModel.DV_P1_ZERO_P1    ]* exp( bx*lbModel.cx[lbModel.DV_P1_ZERO_P1    ] + by*lbModel.cy[lbModel.DV_P1_ZERO_P1    ] + bz*lbModel.cz[lbModel.DV_P1_ZERO_P1    ] + ga*lbModel.cc[lbModel.DV_P1_ZERO_P1    ]) ;
            lbModel.fTemp[lbModel.DV_M1_ZERO_P1    ]    = lbModel.wt[lbModel.DV_M1_ZERO_P1    ]* exp( bx*lbModel.cx[lbModel.DV_M1_ZERO_P1    ] + by*lbModel.cy[lbModel.DV_M1_ZERO_P1    ] + bz*lbModel.cz[lbModel.DV_M1_ZERO_P1    ] + ga*lbModel.cc[lbModel.DV_M1_ZERO_P1    ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_P1_M1    ]    = lbModel.wt[lbModel.DV_ZERO_P1_M1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P1_M1    ] + by*lbModel.cy[lbModel.DV_ZERO_P1_M1    ] + bz*lbModel.cz[lbModel.DV_ZERO_P1_M1    ] + ga*lbModel.cc[lbModel.DV_ZERO_P1_M1    ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_M1_M1    ]    = lbModel.wt[lbModel.DV_ZERO_M1_M1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M1_M1    ] + by*lbModel.cy[lbModel.DV_ZERO_M1_M1    ] + bz*lbModel.cz[lbModel.DV_ZERO_M1_M1    ] + ga*lbModel.cc[lbModel.DV_ZERO_M1_M1    ]) ;
            lbModel.fTemp[lbModel.DV_P1_ZERO_M1    ]    = lbModel.wt[lbModel.DV_P1_ZERO_M1    ]* exp( bx*lbModel.cx[lbModel.DV_P1_ZERO_M1    ] + by*lbModel.cy[lbModel.DV_P1_ZERO_M1    ] + bz*lbModel.cz[lbModel.DV_P1_ZERO_M1    ] + ga*lbModel.cc[lbModel.DV_P1_ZERO_M1    ]) ;
            lbModel.fTemp[lbModel.DV_M1_ZERO_M1    ]    = lbModel.wt[lbModel.DV_M1_ZERO_M1    ]* exp( bx*lbModel.cx[lbModel.DV_M1_ZERO_M1    ] + by*lbModel.cy[lbModel.DV_M1_ZERO_M1    ] + bz*lbModel.cz[lbModel.DV_M1_ZERO_M1    ] + ga*lbModel.cc[lbModel.DV_M1_ZERO_M1    ]) ;
            lbModel.fTemp[lbModel.DV_P1_P1_ZERO    ]    = lbModel.wt[lbModel.DV_P1_P1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_P1_P1_ZERO    ] + by*lbModel.cy[lbModel.DV_P1_P1_ZERO    ] + bz*lbModel.cz[lbModel.DV_P1_P1_ZERO    ] + ga*lbModel.cc[lbModel.DV_P1_P1_ZERO    ]) ;
            lbModel.fTemp[lbModel.DV_M1_P1_ZERO    ]    = lbModel.wt[lbModel.DV_M1_P1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_M1_P1_ZERO    ] + by*lbModel.cy[lbModel.DV_M1_P1_ZERO    ] + bz*lbModel.cz[lbModel.DV_M1_P1_ZERO    ] + ga*lbModel.cc[lbModel.DV_M1_P1_ZERO    ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_P1_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_P1_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P1_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_P1_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_P1_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_P1_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_P1_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_P1_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_P1_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_P1_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_P1_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_P1_ZERO_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_P1_M1_ZERO    ]    = lbModel.wt[lbModel.DV_P1_M1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_P1_M1_ZERO    ] + by*lbModel.cy[lbModel.DV_P1_M1_ZERO    ] + bz*lbModel.cz[lbModel.DV_P1_M1_ZERO    ] + ga*lbModel.cc[lbModel.DV_P1_M1_ZERO    ]) ;
            lbModel.fTemp[lbModel.DV_M1_M1_ZERO    ]    = lbModel.wt[lbModel.DV_M1_M1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_M1_M1_ZERO    ] + by*lbModel.cy[lbModel.DV_M1_M1_ZERO    ] + bz*lbModel.cz[lbModel.DV_M1_M1_ZERO    ] + ga*lbModel.cc[lbModel.DV_M1_M1_ZERO    ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_M1_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_M1_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M1_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_M1_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_M1_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_M1_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_M1_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_M1_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_M1_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_M1_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_M1_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_M1_ZERO_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_P_P_P         ]    = lbModel.wt[lbModel.DV_P_P_P         ]* exp( bx*lbModel.cx[lbModel.DV_P_P_P         ] + by*lbModel.cy[lbModel.DV_P_P_P         ] + bz*lbModel.cz[lbModel.DV_P_P_P         ] + ga*lbModel.cc[lbModel.DV_P_P_P         ]) ;
            lbModel.fTemp[lbModel.DV_M_P_P         ]    = lbModel.wt[lbModel.DV_M_P_P         ]* exp( bx*lbModel.cx[lbModel.DV_M_P_P         ] + by*lbModel.cy[lbModel.DV_M_P_P         ] + bz*lbModel.cz[lbModel.DV_M_P_P         ] + ga*lbModel.cc[lbModel.DV_M_P_P         ]) ;
            lbModel.fTemp[lbModel.DV_M_M_P         ]    = lbModel.wt[lbModel.DV_M_M_P         ]* exp( bx*lbModel.cx[lbModel.DV_M_M_P         ] + by*lbModel.cy[lbModel.DV_M_M_P         ] + bz*lbModel.cz[lbModel.DV_M_M_P         ] + ga*lbModel.cc[lbModel.DV_M_M_P         ]) ;
            lbModel.fTemp[lbModel.DV_P_M_P         ]    = lbModel.wt[lbModel.DV_P_M_P         ]* exp( bx*lbModel.cx[lbModel.DV_P_M_P         ] + by*lbModel.cy[lbModel.DV_P_M_P         ] + bz*lbModel.cz[lbModel.DV_P_M_P         ] + ga*lbModel.cc[lbModel.DV_P_M_P         ]) ;
            lbModel.fTemp[lbModel.DV_P_P_M         ]    = lbModel.wt[lbModel.DV_P_P_M         ]* exp( bx*lbModel.cx[lbModel.DV_P_P_M         ] + by*lbModel.cy[lbModel.DV_P_P_M         ] + bz*lbModel.cz[lbModel.DV_P_P_M         ] + ga*lbModel.cc[lbModel.DV_P_P_M         ]) ;
            lbModel.fTemp[lbModel.DV_M_P_M         ]    = lbModel.wt[lbModel.DV_M_P_M         ]* exp( bx*lbModel.cx[lbModel.DV_M_P_M         ] + by*lbModel.cy[lbModel.DV_M_P_M         ] + bz*lbModel.cz[lbModel.DV_M_P_M         ] + ga*lbModel.cc[lbModel.DV_M_P_M         ]) ;
            lbModel.fTemp[lbModel.DV_M_M_M         ]    = lbModel.wt[lbModel.DV_M_M_M         ]* exp( bx*lbModel.cx[lbModel.DV_M_M_M         ] + by*lbModel.cy[lbModel.DV_M_M_M         ] + bz*lbModel.cz[lbModel.DV_M_M_M         ] + ga*lbModel.cc[lbModel.DV_M_M_M         ]) ;
            lbModel.fTemp[lbModel.DV_P_M_M         ]    = lbModel.wt[lbModel.DV_P_M_M         ]* exp( bx*lbModel.cx[lbModel.DV_P_M_M         ] + by*lbModel.cy[lbModel.DV_P_M_M         ] + bz*lbModel.cz[lbModel.DV_P_M_M         ] + ga*lbModel.cc[lbModel.DV_P_M_M         ]) ;
            lbModel.fTemp[lbModel.DV_P1_P1_P1      ]    = lbModel.wt[lbModel.DV_P1_P1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_P1_P1      ] + by*lbModel.cy[lbModel.DV_P1_P1_P1      ] + bz*lbModel.cz[lbModel.DV_P1_P1_P1      ] + ga*lbModel.cc[lbModel.DV_P1_P1_P1      ]) ;
            lbModel.fTemp[lbModel.DV_M1_P1_P1      ]    = lbModel.wt[lbModel.DV_M1_P1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_P1_P1      ] + by*lbModel.cy[lbModel.DV_M1_P1_P1      ] + bz*lbModel.cz[lbModel.DV_M1_P1_P1      ] + ga*lbModel.cc[lbModel.DV_M1_P1_P1      ]) ;
            lbModel.fTemp[lbModel.DV_M1_M1_P1      ]    = lbModel.wt[lbModel.DV_M1_M1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_M1_P1      ] + by*lbModel.cy[lbModel.DV_M1_M1_P1      ] + bz*lbModel.cz[lbModel.DV_M1_M1_P1      ] + ga*lbModel.cc[lbModel.DV_M1_M1_P1      ]) ;
            lbModel.fTemp[lbModel.DV_P1_M1_P1      ]    = lbModel.wt[lbModel.DV_P1_M1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_M1_P1      ] + by*lbModel.cy[lbModel.DV_P1_M1_P1      ] + bz*lbModel.cz[lbModel.DV_P1_M1_P1      ] + ga*lbModel.cc[lbModel.DV_P1_M1_P1      ]) ;
            lbModel.fTemp[lbModel.DV_P1_P1_M1      ]    = lbModel.wt[lbModel.DV_P1_P1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_P1_M1      ] + by*lbModel.cy[lbModel.DV_P1_P1_M1      ] + bz*lbModel.cz[lbModel.DV_P1_P1_M1      ] + ga*lbModel.cc[lbModel.DV_P1_P1_M1      ]) ;
            lbModel.fTemp[lbModel.DV_M1_P1_M1      ]    = lbModel.wt[lbModel.DV_M1_P1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_P1_M1      ] + by*lbModel.cy[lbModel.DV_M1_P1_M1      ] + bz*lbModel.cz[lbModel.DV_M1_P1_M1      ] + ga*lbModel.cc[lbModel.DV_M1_P1_M1      ]) ;
            lbModel.fTemp[lbModel.DV_M1_M1_M1      ]    = lbModel.wt[lbModel.DV_M1_M1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_M1_M1      ] + by*lbModel.cy[lbModel.DV_M1_M1_M1      ] + bz*lbModel.cz[lbModel.DV_M1_M1_M1      ] + ga*lbModel.cc[lbModel.DV_M1_M1_M1      ]) ;
            lbModel.fTemp[lbModel.DV_P1_M1_M1      ]    = lbModel.wt[lbModel.DV_P1_M1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_M1_M1      ] + by*lbModel.cy[lbModel.DV_P1_M1_M1      ] + bz*lbModel.cz[lbModel.DV_P1_M1_M1      ] + ga*lbModel.cc[lbModel.DV_P1_M1_M1      ]) ;

            for (int iter=0;iter<4;iter++)
            {
                rho1= 0.0;
                Vx  = 0.0;
                Vy  = 0.0;
                Vz  = 0.0;
                mxx = 0.0;
                myy = 0.0;
                mzz = 0.0;
                mxy = 0.0;
                mxz = 0.0;
                myz = 0.0;
                qx  = 0.0;
                qy  = 0.0;
                qz  = 0.0;
                qm  = 0.0;
                q2m = 0.0;

                mat <double,4,4> matA;
                mat <double,4,1> matB;

                for (int dv=0;dv<41;dv++) 
                {
                    rho1+= lbModel.fTemp[dv];
                    Vx  += lbModel.fTemp[dv]*lbModel.cx[dv] ;
                    Vy  += lbModel.fTemp[dv]*lbModel.cy[dv] ;
                    Vz  += lbModel.fTemp[dv]*lbModel.cz[dv] ;
                    mxx += lbModel.fTemp[dv]*lbModel.cx[dv]*lbModel.cx[dv] ;
                    myy += lbModel.fTemp[dv]*lbModel.cy[dv]*lbModel.cy[dv] ;
                    mzz += lbModel.fTemp[dv]*lbModel.cz[dv]*lbModel.cz[dv] ;
                    mxy += lbModel.fTemp[dv]*lbModel.cx[dv]*lbModel.cy[dv] ;
                    mxz += lbModel.fTemp[dv]*lbModel.cx[dv]*lbModel.cz[dv] ;
                    myz += lbModel.fTemp[dv]*lbModel.cz[dv]*lbModel.cy[dv] ;
                    qx  += lbModel.fTemp[dv]*lbModel.cx[dv]*lbModel.cc[dv];
                    qy  += lbModel.fTemp[dv]*lbModel.cy[dv]*lbModel.cc[dv];
                    qz  += lbModel.fTemp[dv]*lbModel.cz[dv]*lbModel.cc[dv];
                    qm  += lbModel.fTemp[dv]*lbModel.cc[dv];
                    q2m += lbModel.fTemp[dv]*lbModel.cc[dv]*lbModel.cc[dv];
                }

                mat <double,4,1> matX;

                matA(0,0) = mxx - uX[i]*Vx;
                matA(0,1) = mxy - uX[i]*Vy;
                matA(0,2) = mxz - uX[i]*Vz;
                matA(0,3) = qx  - uX[i]*qm;

                matA(1,0) = mxy - uY[i]*Vx;
                matA(1,1) = myy - uY[i]*Vy;
                matA(1,2) = myz - uY[i]*Vz;
                matA(1,3) = qy  - uY[i]*qm;

                matA(2,0) = mxz - uZ[i]*Vx;
                matA(2,1) = myz - uZ[i]*Vy;
                matA(2,2) = mzz - uZ[i]*Vz;
                matA(2,3) = qz  - uZ[i]*qm;

                matA(3,0) = qx - E*Vx;
                matA(3,1) = qy - E*Vy;
                matA(3,2) = qz - E*Vz;
                matA(3,3) = q2m- E*qm;

                matB(0,0) = rho1*uX[i] - Vx;
                matB(1,0) = rho1*uY[i] - Vy;
                matB(2,0) = rho1*uZ[i] - Vz;
                matB(3,0) = rho1*E - qm;

                matX = gaussElimination(matA,matB);

                dBx = matX(0,0);
                dBy = matX(1,0);
                dBz = matX(2,0);
                dga = matX(3,0);

                bx += dBx;
                by += dBy;
                bz += dBz;
                ga += dga;
            
                lbModel.fTemp[lbModel.DV_ZERO_ZERO_ZERO]     = lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_ZERO] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_ZERO] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_ZERO] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_ZERO]) ;
                lbModel.fTemp[lbModel.DV_ZERO_ZERO_P1  ]     = lbModel.wt[lbModel.DV_ZERO_ZERO_P1  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_P1  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_P1  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_P1  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_P1  ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_ZERO_P2  ]     = lbModel.wt[lbModel.DV_ZERO_ZERO_P2  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_P2  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_P2  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_P2  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_P2  ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_P2_ZERO  ]     = lbModel.wt[lbModel.DV_ZERO_P2_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P2_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_P2_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_P2_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_P2_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO  ]     = lbModel.wt[lbModel.DV_P2_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_P2_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_P2_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_P2_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_P2_ZERO_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_ZERO_M1  ]     = lbModel.wt[lbModel.DV_ZERO_ZERO_M1  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_M1  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_M1  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_M1  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_M1  ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_ZERO_M2  ]     = lbModel.wt[lbModel.DV_ZERO_ZERO_M2  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_M2  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_M2  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_M2  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_M2  ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_M2_ZERO  ]     = lbModel.wt[lbModel.DV_ZERO_M2_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M2_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_M2_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_M2_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_M2_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_M2_ZERO_ZERO  ]     = lbModel.wt[lbModel.DV_M2_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_M2_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_M2_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_M2_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_M2_ZERO_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_P1_P1    ]     = lbModel.wt[lbModel.DV_ZERO_P1_P1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P1_P1    ] + by*lbModel.cy[lbModel.DV_ZERO_P1_P1    ] + bz*lbModel.cz[lbModel.DV_ZERO_P1_P1    ] + ga*lbModel.cc[lbModel.DV_ZERO_P1_P1    ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_M1_P1    ]     = lbModel.wt[lbModel.DV_ZERO_M1_P1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M1_P1    ] + by*lbModel.cy[lbModel.DV_ZERO_M1_P1    ] + bz*lbModel.cz[lbModel.DV_ZERO_M1_P1    ] + ga*lbModel.cc[lbModel.DV_ZERO_M1_P1    ]) ;
                lbModel.fTemp[lbModel.DV_P1_ZERO_P1    ]     = lbModel.wt[lbModel.DV_P1_ZERO_P1    ]* exp( bx*lbModel.cx[lbModel.DV_P1_ZERO_P1    ] + by*lbModel.cy[lbModel.DV_P1_ZERO_P1    ] + bz*lbModel.cz[lbModel.DV_P1_ZERO_P1    ] + ga*lbModel.cc[lbModel.DV_P1_ZERO_P1    ]) ;
                lbModel.fTemp[lbModel.DV_M1_ZERO_P1    ]     = lbModel.wt[lbModel.DV_M1_ZERO_P1    ]* exp( bx*lbModel.cx[lbModel.DV_M1_ZERO_P1    ] + by*lbModel.cy[lbModel.DV_M1_ZERO_P1    ] + bz*lbModel.cz[lbModel.DV_M1_ZERO_P1    ] + ga*lbModel.cc[lbModel.DV_M1_ZERO_P1    ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_P1_M1    ]     = lbModel.wt[lbModel.DV_ZERO_P1_M1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P1_M1    ] + by*lbModel.cy[lbModel.DV_ZERO_P1_M1    ] + bz*lbModel.cz[lbModel.DV_ZERO_P1_M1    ] + ga*lbModel.cc[lbModel.DV_ZERO_P1_M1    ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_M1_M1    ]     = lbModel.wt[lbModel.DV_ZERO_M1_M1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M1_M1    ] + by*lbModel.cy[lbModel.DV_ZERO_M1_M1    ] + bz*lbModel.cz[lbModel.DV_ZERO_M1_M1    ] + ga*lbModel.cc[lbModel.DV_ZERO_M1_M1    ]) ;
                lbModel.fTemp[lbModel.DV_P1_ZERO_M1    ]     = lbModel.wt[lbModel.DV_P1_ZERO_M1    ]* exp( bx*lbModel.cx[lbModel.DV_P1_ZERO_M1    ] + by*lbModel.cy[lbModel.DV_P1_ZERO_M1    ] + bz*lbModel.cz[lbModel.DV_P1_ZERO_M1    ] + ga*lbModel.cc[lbModel.DV_P1_ZERO_M1    ]) ;
                lbModel.fTemp[lbModel.DV_M1_ZERO_M1    ]     = lbModel.wt[lbModel.DV_M1_ZERO_M1    ]* exp( bx*lbModel.cx[lbModel.DV_M1_ZERO_M1    ] + by*lbModel.cy[lbModel.DV_M1_ZERO_M1    ] + bz*lbModel.cz[lbModel.DV_M1_ZERO_M1    ] + ga*lbModel.cc[lbModel.DV_M1_ZERO_M1    ]) ;
                lbModel.fTemp[lbModel.DV_P1_P1_ZERO    ]     = lbModel.wt[lbModel.DV_P1_P1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_P1_P1_ZERO    ] + by*lbModel.cy[lbModel.DV_P1_P1_ZERO    ] + bz*lbModel.cz[lbModel.DV_P1_P1_ZERO    ] + ga*lbModel.cc[lbModel.DV_P1_P1_ZERO    ]) ;
                lbModel.fTemp[lbModel.DV_M1_P1_ZERO    ]     = lbModel.wt[lbModel.DV_M1_P1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_M1_P1_ZERO    ] + by*lbModel.cy[lbModel.DV_M1_P1_ZERO    ] + bz*lbModel.cz[lbModel.DV_M1_P1_ZERO    ] + ga*lbModel.cc[lbModel.DV_M1_P1_ZERO    ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_P1_ZERO  ]     = lbModel.wt[lbModel.DV_ZERO_P1_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P1_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_P1_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_P1_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_P1_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_P1_ZERO_ZERO  ]     = lbModel.wt[lbModel.DV_P1_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_P1_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_P1_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_P1_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_P1_ZERO_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_P1_M1_ZERO    ]     = lbModel.wt[lbModel.DV_P1_M1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_P1_M1_ZERO    ] + by*lbModel.cy[lbModel.DV_P1_M1_ZERO    ] + bz*lbModel.cz[lbModel.DV_P1_M1_ZERO    ] + ga*lbModel.cc[lbModel.DV_P1_M1_ZERO    ]) ;
                lbModel.fTemp[lbModel.DV_M1_M1_ZERO    ]     = lbModel.wt[lbModel.DV_M1_M1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_M1_M1_ZERO    ] + by*lbModel.cy[lbModel.DV_M1_M1_ZERO    ] + bz*lbModel.cz[lbModel.DV_M1_M1_ZERO    ] + ga*lbModel.cc[lbModel.DV_M1_M1_ZERO    ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_M1_ZERO  ]     = lbModel.wt[lbModel.DV_ZERO_M1_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M1_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_M1_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_M1_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_M1_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_M1_ZERO_ZERO  ]     = lbModel.wt[lbModel.DV_M1_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_M1_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_M1_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_M1_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_M1_ZERO_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_P_P_P         ]     = lbModel.wt[lbModel.DV_P_P_P         ]* exp( bx*lbModel.cx[lbModel.DV_P_P_P         ] + by*lbModel.cy[lbModel.DV_P_P_P         ] + bz*lbModel.cz[lbModel.DV_P_P_P         ] + ga*lbModel.cc[lbModel.DV_P_P_P         ]) ;
                lbModel.fTemp[lbModel.DV_M_P_P         ]     = lbModel.wt[lbModel.DV_M_P_P         ]* exp( bx*lbModel.cx[lbModel.DV_M_P_P         ] + by*lbModel.cy[lbModel.DV_M_P_P         ] + bz*lbModel.cz[lbModel.DV_M_P_P         ] + ga*lbModel.cc[lbModel.DV_M_P_P         ]) ;
                lbModel.fTemp[lbModel.DV_M_M_P         ]     = lbModel.wt[lbModel.DV_M_M_P         ]* exp( bx*lbModel.cx[lbModel.DV_M_M_P         ] + by*lbModel.cy[lbModel.DV_M_M_P         ] + bz*lbModel.cz[lbModel.DV_M_M_P         ] + ga*lbModel.cc[lbModel.DV_M_M_P         ]) ;
                lbModel.fTemp[lbModel.DV_P_M_P         ]     = lbModel.wt[lbModel.DV_P_M_P         ]* exp( bx*lbModel.cx[lbModel.DV_P_M_P         ] + by*lbModel.cy[lbModel.DV_P_M_P         ] + bz*lbModel.cz[lbModel.DV_P_M_P         ] + ga*lbModel.cc[lbModel.DV_P_M_P         ]) ;
                lbModel.fTemp[lbModel.DV_P_P_M         ]     = lbModel.wt[lbModel.DV_P_P_M         ]* exp( bx*lbModel.cx[lbModel.DV_P_P_M         ] + by*lbModel.cy[lbModel.DV_P_P_M         ] + bz*lbModel.cz[lbModel.DV_P_P_M         ] + ga*lbModel.cc[lbModel.DV_P_P_M         ]) ;
                lbModel.fTemp[lbModel.DV_M_P_M         ]     = lbModel.wt[lbModel.DV_M_P_M         ]* exp( bx*lbModel.cx[lbModel.DV_M_P_M         ] + by*lbModel.cy[lbModel.DV_M_P_M         ] + bz*lbModel.cz[lbModel.DV_M_P_M         ] + ga*lbModel.cc[lbModel.DV_M_P_M         ]) ;
                lbModel.fTemp[lbModel.DV_M_M_M         ]     = lbModel.wt[lbModel.DV_M_M_M         ]* exp( bx*lbModel.cx[lbModel.DV_M_M_M         ] + by*lbModel.cy[lbModel.DV_M_M_M         ] + bz*lbModel.cz[lbModel.DV_M_M_M         ] + ga*lbModel.cc[lbModel.DV_M_M_M         ]) ;
                lbModel.fTemp[lbModel.DV_P_M_M         ]     = lbModel.wt[lbModel.DV_P_M_M         ]* exp( bx*lbModel.cx[lbModel.DV_P_M_M         ] + by*lbModel.cy[lbModel.DV_P_M_M         ] + bz*lbModel.cz[lbModel.DV_P_M_M         ] + ga*lbModel.cc[lbModel.DV_P_M_M         ]) ;
                lbModel.fTemp[lbModel.DV_P1_P1_P1      ]     = lbModel.wt[lbModel.DV_P1_P1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_P1_P1      ] + by*lbModel.cy[lbModel.DV_P1_P1_P1      ] + bz*lbModel.cz[lbModel.DV_P1_P1_P1      ] + ga*lbModel.cc[lbModel.DV_P1_P1_P1      ]) ;
                lbModel.fTemp[lbModel.DV_M1_P1_P1      ]     = lbModel.wt[lbModel.DV_M1_P1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_P1_P1      ] + by*lbModel.cy[lbModel.DV_M1_P1_P1      ] + bz*lbModel.cz[lbModel.DV_M1_P1_P1      ] + ga*lbModel.cc[lbModel.DV_M1_P1_P1      ]) ;
                lbModel.fTemp[lbModel.DV_M1_M1_P1      ]     = lbModel.wt[lbModel.DV_M1_M1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_M1_P1      ] + by*lbModel.cy[lbModel.DV_M1_M1_P1      ] + bz*lbModel.cz[lbModel.DV_M1_M1_P1      ] + ga*lbModel.cc[lbModel.DV_M1_M1_P1      ]) ;
                lbModel.fTemp[lbModel.DV_P1_M1_P1      ]     = lbModel.wt[lbModel.DV_P1_M1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_M1_P1      ] + by*lbModel.cy[lbModel.DV_P1_M1_P1      ] + bz*lbModel.cz[lbModel.DV_P1_M1_P1      ] + ga*lbModel.cc[lbModel.DV_P1_M1_P1      ]) ;
                lbModel.fTemp[lbModel.DV_P1_P1_M1      ]     = lbModel.wt[lbModel.DV_P1_P1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_P1_M1      ] + by*lbModel.cy[lbModel.DV_P1_P1_M1      ] + bz*lbModel.cz[lbModel.DV_P1_P1_M1      ] + ga*lbModel.cc[lbModel.DV_P1_P1_M1      ]) ;
                lbModel.fTemp[lbModel.DV_M1_P1_M1      ]     = lbModel.wt[lbModel.DV_M1_P1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_P1_M1      ] + by*lbModel.cy[lbModel.DV_M1_P1_M1      ] + bz*lbModel.cz[lbModel.DV_M1_P1_M1      ] + ga*lbModel.cc[lbModel.DV_M1_P1_M1      ]) ;
                lbModel.fTemp[lbModel.DV_M1_M1_M1      ]     = lbModel.wt[lbModel.DV_M1_M1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_M1_M1      ] + by*lbModel.cy[lbModel.DV_M1_M1_M1      ] + bz*lbModel.cz[lbModel.DV_M1_M1_M1      ] + ga*lbModel.cc[lbModel.DV_M1_M1_M1      ]) ;
                lbModel.fTemp[lbModel.DV_P1_M1_M1      ]     = lbModel.wt[lbModel.DV_P1_M1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_M1_M1      ] + by*lbModel.cy[lbModel.DV_P1_M1_M1      ] + bz*lbModel.cz[lbModel.DV_P1_M1_M1      ] + ga*lbModel.cc[lbModel.DV_P1_M1_M1      ]) ;

            }    

                rho1 = lbModel.fTemp[lbModel.DV_ZERO_ZERO_ZERO]
                     + lbModel.fTemp[lbModel.DV_ZERO_ZERO_P1  ]
                     + lbModel.fTemp[lbModel.DV_ZERO_ZERO_P2  ]
                     + lbModel.fTemp[lbModel.DV_ZERO_P2_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_ZERO_ZERO_M1  ]
                     + lbModel.fTemp[lbModel.DV_ZERO_ZERO_M2  ]
                     + lbModel.fTemp[lbModel.DV_ZERO_M2_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_M2_ZERO_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_ZERO_P1_P1    ]
                     + lbModel.fTemp[lbModel.DV_ZERO_M1_P1    ]
                     + lbModel.fTemp[lbModel.DV_P1_ZERO_P1    ]
                     + lbModel.fTemp[lbModel.DV_M1_ZERO_P1    ]
                     + lbModel.fTemp[lbModel.DV_ZERO_P1_M1    ]
                     + lbModel.fTemp[lbModel.DV_ZERO_M1_M1    ]
                     + lbModel.fTemp[lbModel.DV_P1_ZERO_M1    ]
                     + lbModel.fTemp[lbModel.DV_M1_ZERO_M1    ]
                     + lbModel.fTemp[lbModel.DV_P1_P1_ZERO    ]
                     + lbModel.fTemp[lbModel.DV_M1_P1_ZERO    ]
                     + lbModel.fTemp[lbModel.DV_ZERO_P1_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_P1_ZERO_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_P1_M1_ZERO    ]
                     + lbModel.fTemp[lbModel.DV_M1_M1_ZERO    ]
                     + lbModel.fTemp[lbModel.DV_ZERO_M1_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_M1_ZERO_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_P_P_P         ]
                     + lbModel.fTemp[lbModel.DV_M_P_P         ]
                     + lbModel.fTemp[lbModel.DV_M_M_P         ]
                     + lbModel.fTemp[lbModel.DV_P_M_P         ]
                     + lbModel.fTemp[lbModel.DV_P_P_M         ]
                     + lbModel.fTemp[lbModel.DV_M_P_M         ]
                     + lbModel.fTemp[lbModel.DV_M_M_M         ]
                     + lbModel.fTemp[lbModel.DV_P_M_M         ]
                     + lbModel.fTemp[lbModel.DV_P1_P1_P1      ]
                     + lbModel.fTemp[lbModel.DV_M1_P1_P1      ]
                     + lbModel.fTemp[lbModel.DV_M1_M1_P1      ]
                     + lbModel.fTemp[lbModel.DV_P1_M1_P1      ]
                     + lbModel.fTemp[lbModel.DV_P1_P1_M1      ]
                     + lbModel.fTemp[lbModel.DV_M1_P1_M1      ]
                     + lbModel.fTemp[lbModel.DV_M1_M1_M1      ]
                     + lbModel.fTemp[lbModel.DV_P1_M1_M1      ];
                
                
                
                double rhoRatio = rho[i]/rho1 ;  
                
                lbModel.fTemp0 [i][lbModel.CENTER_DV_ZERO_ZERO_ZERO] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_ZERO_ZERO];
                lbModel.fTemp1 [i][lbModel.G1_DV_ZERO_ZERO_P1      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_ZERO_P1  ];
                lbModel.fTemp1 [i][lbModel.G1_DV_ZERO_ZERO_P2      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_ZERO_P2  ];
                lbModel.fTemp1 [i][lbModel.G1_DV_ZERO_P2_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_P2_ZERO  ];
                lbModel.fTemp1 [i][lbModel.G1_DV_P2_ZERO_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO  ];
                lbModel.fTemp2 [i][lbModel.G2_DV_ZERO_ZERO_M1      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_ZERO_M1  ];
                lbModel.fTemp2 [i][lbModel.G2_DV_ZERO_ZERO_M2      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_ZERO_M2  ];
                lbModel.fTemp2 [i][lbModel.G2_DV_ZERO_M2_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_M2_ZERO  ];
                lbModel.fTemp2 [i][lbModel.G2_DV_M2_ZERO_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_M2_ZERO_ZERO  ];
                lbModel.fTemp3 [i][lbModel.G3_DV_ZERO_P1_P1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_P1_P1    ];
                lbModel.fTemp3 [i][lbModel.G3_DV_ZERO_M1_P1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_M1_P1    ];
                lbModel.fTemp3 [i][lbModel.G3_DV_P1_ZERO_P1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_ZERO_P1    ];
                lbModel.fTemp3 [i][lbModel.G3_DV_M1_ZERO_P1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_ZERO_P1    ];
                lbModel.fTemp4 [i][lbModel.G4_DV_ZERO_P1_M1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_P1_M1    ];
                lbModel.fTemp4 [i][lbModel.G4_DV_ZERO_M1_M1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_M1_M1    ];
                lbModel.fTemp4 [i][lbModel.G4_DV_P1_ZERO_M1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_ZERO_M1    ];
                lbModel.fTemp4 [i][lbModel.G4_DV_M1_ZERO_M1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_ZERO_M1    ];
                lbModel.fTemp5 [i][lbModel.G5_DV_P1_P1_ZERO        ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_P1_ZERO    ];
                lbModel.fTemp5 [i][lbModel.G5_DV_M1_P1_ZERO        ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_P1_ZERO    ];
                lbModel.fTemp5 [i][lbModel.G5_DV_ZERO_P1_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_P1_ZERO  ];
                lbModel.fTemp5 [i][lbModel.G5_DV_P1_ZERO_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_ZERO_ZERO  ];
                lbModel.fTemp6 [i][lbModel.G6_DV_P1_M1_ZERO        ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_M1_ZERO    ];
                lbModel.fTemp6 [i][lbModel.G6_DV_M1_M1_ZERO        ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_M1_ZERO    ];
                lbModel.fTemp6 [i][lbModel.G6_DV_ZERO_M1_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_M1_ZERO  ];
                lbModel.fTemp6 [i][lbModel.G6_DV_M1_ZERO_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_ZERO_ZERO  ];
                lbModel.fTemp7 [i][lbModel.G7_DV_P_P_P             ] = rhoRatio * lbModel.fTemp[lbModel.DV_P_P_P         ];
                lbModel.fTemp7 [i][lbModel.G7_DV_M_P_P             ] = rhoRatio * lbModel.fTemp[lbModel.DV_M_P_P         ];
                lbModel.fTemp7 [i][lbModel.G7_DV_M_M_P             ] = rhoRatio * lbModel.fTemp[lbModel.DV_M_M_P         ];
                lbModel.fTemp7 [i][lbModel.G7_DV_P_M_P             ] = rhoRatio * lbModel.fTemp[lbModel.DV_P_M_P         ];
                lbModel.fTemp8 [i][lbModel.G8_DV_P_P_M             ] = rhoRatio * lbModel.fTemp[lbModel.DV_P_P_M         ];
                lbModel.fTemp8 [i][lbModel.G8_DV_M_P_M             ] = rhoRatio * lbModel.fTemp[lbModel.DV_M_P_M         ];
                lbModel.fTemp8 [i][lbModel.G8_DV_M_M_M             ] = rhoRatio * lbModel.fTemp[lbModel.DV_M_M_M         ];
                lbModel.fTemp8 [i][lbModel.G8_DV_P_M_M             ] = rhoRatio * lbModel.fTemp[lbModel.DV_P_M_M         ];
                lbModel.fTemp9 [i][lbModel.G9_DV_P1_P1_P1          ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_P1_P1      ];
                lbModel.fTemp9 [i][lbModel.G9_DV_M1_P1_P1          ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_P1_P1      ];
                lbModel.fTemp9 [i][lbModel.G9_DV_M1_M1_P1          ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_M1_P1      ];
                lbModel.fTemp9 [i][lbModel.G9_DV_P1_M1_P1          ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_M1_P1      ];
                lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1         ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_P1_M1      ];
                lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1         ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_P1_M1      ];
                lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1         ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_M1_M1      ];
                lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1         ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_M1_M1      ];             
                
        }
    }




/*
    template<typename dataType1>
    void getIterativeFEq(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta)
    {
        dataType1 uSq, E;

        dataType1 rho1 ,mxy ,mxz ,myz ,bx ,by ,bz ;
        dataType1 mxx ,myy ,mzz ,Vx ,Vy ,Vz ,qx ,qy ,qz ,qm ,q2m ;
        dataType1 dBx ,dBy ,dBz ,dga , del_theta ,ga ;

        dataType1 dot;
        
        for(int i=0;i<VECT_LENGTH;i++)
        {
            uSq = uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i];
            E   = uSq + 3.0*theta[i];
            bx  = uX[i]/theta[i];
            by  = uY[i]/theta[i];
            bz  = uZ[i]/theta[i];
            del_theta = (theta[i]-lbModel.theta0)*lbModel.oneByTheta0 ;
            ga  = 0.5*lbModel.oneByTheta0*del_theta*( 1.0-del_theta*del_theta ) / (1.0+del_theta); 
        
            lbModel.fTemp[lbModel.DV_ZERO_ZERO_ZERO]    = lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_ZERO] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_ZERO] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_ZERO] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_ZERO]) ;
            lbModel.fTemp[lbModel.DV_ZERO_ZERO_P1  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_P1  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_P1  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_P1  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_P1  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_P1  ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_ZERO_P2  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_P2  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_P2  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_P2  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_P2  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_P2  ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_P2_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_P2_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P2_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_P2_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_P2_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_P2_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_P2_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_P2_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_P2_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_P2_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_P2_ZERO_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_ZERO_M1  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_M1  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_M1  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_M1  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_M1  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_M1  ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_ZERO_M2  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_M2  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_M2  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_M2  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_M2  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_M2  ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_M2_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_M2_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M2_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_M2_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_M2_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_M2_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_M2_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_M2_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_M2_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_M2_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_M2_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_M2_ZERO_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_P1_P1    ]    = lbModel.wt[lbModel.DV_ZERO_P1_P1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P1_P1    ] + by*lbModel.cy[lbModel.DV_ZERO_P1_P1    ] + bz*lbModel.cz[lbModel.DV_ZERO_P1_P1    ] + ga*lbModel.cc[lbModel.DV_ZERO_P1_P1    ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_M1_P1    ]    = lbModel.wt[lbModel.DV_ZERO_M1_P1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M1_P1    ] + by*lbModel.cy[lbModel.DV_ZERO_M1_P1    ] + bz*lbModel.cz[lbModel.DV_ZERO_M1_P1    ] + ga*lbModel.cc[lbModel.DV_ZERO_M1_P1    ]) ;
            lbModel.fTemp[lbModel.DV_P1_ZERO_P1    ]    = lbModel.wt[lbModel.DV_P1_ZERO_P1    ]* exp( bx*lbModel.cx[lbModel.DV_P1_ZERO_P1    ] + by*lbModel.cy[lbModel.DV_P1_ZERO_P1    ] + bz*lbModel.cz[lbModel.DV_P1_ZERO_P1    ] + ga*lbModel.cc[lbModel.DV_P1_ZERO_P1    ]) ;
            lbModel.fTemp[lbModel.DV_M1_ZERO_P1    ]    = lbModel.wt[lbModel.DV_M1_ZERO_P1    ]* exp( bx*lbModel.cx[lbModel.DV_M1_ZERO_P1    ] + by*lbModel.cy[lbModel.DV_M1_ZERO_P1    ] + bz*lbModel.cz[lbModel.DV_M1_ZERO_P1    ] + ga*lbModel.cc[lbModel.DV_M1_ZERO_P1    ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_P1_M1    ]    = lbModel.wt[lbModel.DV_ZERO_P1_M1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P1_M1    ] + by*lbModel.cy[lbModel.DV_ZERO_P1_M1    ] + bz*lbModel.cz[lbModel.DV_ZERO_P1_M1    ] + ga*lbModel.cc[lbModel.DV_ZERO_P1_M1    ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_M1_M1    ]    = lbModel.wt[lbModel.DV_ZERO_M1_M1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M1_M1    ] + by*lbModel.cy[lbModel.DV_ZERO_M1_M1    ] + bz*lbModel.cz[lbModel.DV_ZERO_M1_M1    ] + ga*lbModel.cc[lbModel.DV_ZERO_M1_M1    ]) ;
            lbModel.fTemp[lbModel.DV_P1_ZERO_M1    ]    = lbModel.wt[lbModel.DV_P1_ZERO_M1    ]* exp( bx*lbModel.cx[lbModel.DV_P1_ZERO_M1    ] + by*lbModel.cy[lbModel.DV_P1_ZERO_M1    ] + bz*lbModel.cz[lbModel.DV_P1_ZERO_M1    ] + ga*lbModel.cc[lbModel.DV_P1_ZERO_M1    ]) ;
            lbModel.fTemp[lbModel.DV_M1_ZERO_M1    ]    = lbModel.wt[lbModel.DV_M1_ZERO_M1    ]* exp( bx*lbModel.cx[lbModel.DV_M1_ZERO_M1    ] + by*lbModel.cy[lbModel.DV_M1_ZERO_M1    ] + bz*lbModel.cz[lbModel.DV_M1_ZERO_M1    ] + ga*lbModel.cc[lbModel.DV_M1_ZERO_M1    ]) ;
            lbModel.fTemp[lbModel.DV_P1_P1_ZERO    ]    = lbModel.wt[lbModel.DV_P1_P1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_P1_P1_ZERO    ] + by*lbModel.cy[lbModel.DV_P1_P1_ZERO    ] + bz*lbModel.cz[lbModel.DV_P1_P1_ZERO    ] + ga*lbModel.cc[lbModel.DV_P1_P1_ZERO    ]) ;
            lbModel.fTemp[lbModel.DV_M1_P1_ZERO    ]    = lbModel.wt[lbModel.DV_M1_P1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_M1_P1_ZERO    ] + by*lbModel.cy[lbModel.DV_M1_P1_ZERO    ] + bz*lbModel.cz[lbModel.DV_M1_P1_ZERO    ] + ga*lbModel.cc[lbModel.DV_M1_P1_ZERO    ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_P1_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_P1_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P1_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_P1_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_P1_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_P1_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_P1_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_P1_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_P1_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_P1_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_P1_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_P1_ZERO_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_P1_M1_ZERO    ]    = lbModel.wt[lbModel.DV_P1_M1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_P1_M1_ZERO    ] + by*lbModel.cy[lbModel.DV_P1_M1_ZERO    ] + bz*lbModel.cz[lbModel.DV_P1_M1_ZERO    ] + ga*lbModel.cc[lbModel.DV_P1_M1_ZERO    ]) ;
            lbModel.fTemp[lbModel.DV_M1_M1_ZERO    ]    = lbModel.wt[lbModel.DV_M1_M1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_M1_M1_ZERO    ] + by*lbModel.cy[lbModel.DV_M1_M1_ZERO    ] + bz*lbModel.cz[lbModel.DV_M1_M1_ZERO    ] + ga*lbModel.cc[lbModel.DV_M1_M1_ZERO    ]) ;
            lbModel.fTemp[lbModel.DV_ZERO_M1_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_M1_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M1_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_M1_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_M1_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_M1_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_M1_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_M1_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_M1_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_M1_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_M1_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_M1_ZERO_ZERO  ]) ;
            lbModel.fTemp[lbModel.DV_P_P_P         ]    = lbModel.wt[lbModel.DV_P_P_P         ]* exp( bx*lbModel.cx[lbModel.DV_P_P_P         ] + by*lbModel.cy[lbModel.DV_P_P_P         ] + bz*lbModel.cz[lbModel.DV_P_P_P         ] + ga*lbModel.cc[lbModel.DV_P_P_P         ]) ;
            lbModel.fTemp[lbModel.DV_M_P_P         ]    = lbModel.wt[lbModel.DV_M_P_P         ]* exp( bx*lbModel.cx[lbModel.DV_M_P_P         ] + by*lbModel.cy[lbModel.DV_M_P_P         ] + bz*lbModel.cz[lbModel.DV_M_P_P         ] + ga*lbModel.cc[lbModel.DV_M_P_P         ]) ;
            lbModel.fTemp[lbModel.DV_M_M_P         ]    = lbModel.wt[lbModel.DV_M_M_P         ]* exp( bx*lbModel.cx[lbModel.DV_M_M_P         ] + by*lbModel.cy[lbModel.DV_M_M_P         ] + bz*lbModel.cz[lbModel.DV_M_M_P         ] + ga*lbModel.cc[lbModel.DV_M_M_P         ]) ;
            lbModel.fTemp[lbModel.DV_P_M_P         ]    = lbModel.wt[lbModel.DV_P_M_P         ]* exp( bx*lbModel.cx[lbModel.DV_P_M_P         ] + by*lbModel.cy[lbModel.DV_P_M_P         ] + bz*lbModel.cz[lbModel.DV_P_M_P         ] + ga*lbModel.cc[lbModel.DV_P_M_P         ]) ;
            lbModel.fTemp[lbModel.DV_P_P_M         ]    = lbModel.wt[lbModel.DV_P_P_M         ]* exp( bx*lbModel.cx[lbModel.DV_P_P_M         ] + by*lbModel.cy[lbModel.DV_P_P_M         ] + bz*lbModel.cz[lbModel.DV_P_P_M         ] + ga*lbModel.cc[lbModel.DV_P_P_M         ]) ;
            lbModel.fTemp[lbModel.DV_M_P_M         ]    = lbModel.wt[lbModel.DV_M_P_M         ]* exp( bx*lbModel.cx[lbModel.DV_M_P_M         ] + by*lbModel.cy[lbModel.DV_M_P_M         ] + bz*lbModel.cz[lbModel.DV_M_P_M         ] + ga*lbModel.cc[lbModel.DV_M_P_M         ]) ;
            lbModel.fTemp[lbModel.DV_M_M_M         ]    = lbModel.wt[lbModel.DV_M_M_M         ]* exp( bx*lbModel.cx[lbModel.DV_M_M_M         ] + by*lbModel.cy[lbModel.DV_M_M_M         ] + bz*lbModel.cz[lbModel.DV_M_M_M         ] + ga*lbModel.cc[lbModel.DV_M_M_M         ]) ;
            lbModel.fTemp[lbModel.DV_P_M_M         ]    = lbModel.wt[lbModel.DV_P_M_M         ]* exp( bx*lbModel.cx[lbModel.DV_P_M_M         ] + by*lbModel.cy[lbModel.DV_P_M_M         ] + bz*lbModel.cz[lbModel.DV_P_M_M         ] + ga*lbModel.cc[lbModel.DV_P_M_M         ]) ;
            lbModel.fTemp[lbModel.DV_P1_P1_P1      ]    = lbModel.wt[lbModel.DV_P1_P1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_P1_P1      ] + by*lbModel.cy[lbModel.DV_P1_P1_P1      ] + bz*lbModel.cz[lbModel.DV_P1_P1_P1      ] + ga*lbModel.cc[lbModel.DV_P1_P1_P1      ]) ;
            lbModel.fTemp[lbModel.DV_M1_P1_P1      ]    = lbModel.wt[lbModel.DV_M1_P1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_P1_P1      ] + by*lbModel.cy[lbModel.DV_M1_P1_P1      ] + bz*lbModel.cz[lbModel.DV_M1_P1_P1      ] + ga*lbModel.cc[lbModel.DV_M1_P1_P1      ]) ;
            lbModel.fTemp[lbModel.DV_M1_M1_P1      ]    = lbModel.wt[lbModel.DV_M1_M1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_M1_P1      ] + by*lbModel.cy[lbModel.DV_M1_M1_P1      ] + bz*lbModel.cz[lbModel.DV_M1_M1_P1      ] + ga*lbModel.cc[lbModel.DV_M1_M1_P1      ]) ;
            lbModel.fTemp[lbModel.DV_P1_M1_P1      ]    = lbModel.wt[lbModel.DV_P1_M1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_M1_P1      ] + by*lbModel.cy[lbModel.DV_P1_M1_P1      ] + bz*lbModel.cz[lbModel.DV_P1_M1_P1      ] + ga*lbModel.cc[lbModel.DV_P1_M1_P1      ]) ;
            lbModel.fTemp[lbModel.DV_P1_P1_M1      ]    = lbModel.wt[lbModel.DV_P1_P1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_P1_M1      ] + by*lbModel.cy[lbModel.DV_P1_P1_M1      ] + bz*lbModel.cz[lbModel.DV_P1_P1_M1      ] + ga*lbModel.cc[lbModel.DV_P1_P1_M1      ]) ;
            lbModel.fTemp[lbModel.DV_M1_P1_M1      ]    = lbModel.wt[lbModel.DV_M1_P1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_P1_M1      ] + by*lbModel.cy[lbModel.DV_M1_P1_M1      ] + bz*lbModel.cz[lbModel.DV_M1_P1_M1      ] + ga*lbModel.cc[lbModel.DV_M1_P1_M1      ]) ;
            lbModel.fTemp[lbModel.DV_M1_M1_M1      ]    = lbModel.wt[lbModel.DV_M1_M1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_M1_M1      ] + by*lbModel.cy[lbModel.DV_M1_M1_M1      ] + bz*lbModel.cz[lbModel.DV_M1_M1_M1      ] + ga*lbModel.cc[lbModel.DV_M1_M1_M1      ]) ;
            lbModel.fTemp[lbModel.DV_P1_M1_M1      ]    = lbModel.wt[lbModel.DV_P1_M1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_M1_M1      ] + by*lbModel.cy[lbModel.DV_P1_M1_M1      ] + bz*lbModel.cz[lbModel.DV_P1_M1_M1      ] + ga*lbModel.cc[lbModel.DV_P1_M1_M1      ]) ;

            for (int iter=0;iter<1;iter++)
            {
                rho1= 0.0;
                Vx  = 0.0;
                Vy  = 0.0;
                Vz  = 0.0;
                mxx = 0.0;
                myy = 0.0;
                mzz = 0.0;
                mxy = 0.0;
                mxz = 0.0;
                myz = 0.0;
                qx  = 0.0;
                qy  = 0.0;
                qz  = 0.0;
                qm  = 0.0;
                q2m = 0.0;

                mat <double,4,4> matA;
                mat <double,4,1> matB;

                for (int dv=0;dv<41;dv++) 
                {
                    rho1+= lbModel.fTemp[dv];
                    Vx  += lbModel.fTemp[dv]*lbModel.cx[dv] ;
                    Vy  += lbModel.fTemp[dv]*lbModel.cy[dv] ;
                    Vz  += lbModel.fTemp[dv]*lbModel.cz[dv] ;
                    mxx += lbModel.fTemp[dv]*lbModel.cx[dv]*lbModel.cx[dv] ;
                    myy += lbModel.fTemp[dv]*lbModel.cy[dv]*lbModel.cy[dv] ;
                    mzz += lbModel.fTemp[dv]*lbModel.cz[dv]*lbModel.cz[dv] ;
                    mxy += lbModel.fTemp[dv]*lbModel.cx[dv]*lbModel.cy[dv] ;
                    mxz += lbModel.fTemp[dv]*lbModel.cx[dv]*lbModel.cz[dv] ;
                    myz += lbModel.fTemp[dv]*lbModel.cz[dv]*lbModel.cy[dv] ;
                    qx  += lbModel.fTemp[dv]*lbModel.cx[dv]*lbModel.cc[dv];
                    qy  += lbModel.fTemp[dv]*lbModel.cy[dv]*lbModel.cc[dv];
                    qz  += lbModel.fTemp[dv]*lbModel.cz[dv]*lbModel.cc[dv];
                    qm  += lbModel.fTemp[dv]*lbModel.cc[dv];
                    q2m += lbModel.fTemp[dv]*lbModel.cc[dv]*lbModel.cc[dv];
                }

                mat <double,4,1> matX;

                matA(0,0) = mxx - uX[i]*Vx;
                matA(0,1) = mxy - uX[i]*Vy;
                matA(0,2) = mxz - uX[i]*Vz;
                matA(0,3) = qx  - uX[i]*qm;

                matA(1,0) = mxy - uY[i]*Vx;
                matA(1,1) = myy - uY[i]*Vy;
                matA(1,2) = myz - uY[i]*Vz;
                matA(1,3) = qy  - uY[i]*qm;

                matA(2,0) = mxz - uZ[i]*Vx;
                matA(2,1) = myz - uZ[i]*Vy;
                matA(2,2) = mzz - uZ[i]*Vz;
                matA(2,3) = qz  - uZ[i]*qm;

                matA(3,0) = qx - E*Vx;
                matA(3,1) = qy - E*Vy;
                matA(3,2) = qz - E*Vz;
                matA(3,3) = q2m- E*qm;

                matB(0,0) = rho1*uX[i] - Vx;
                matB(1,0) = rho1*uY[i] - Vy;
                matB(2,0) = rho1*uZ[i] - Vz;
                matB(3,0) = rho1*E - qm;

                matX = gaussElimination(matA,matB);

                dBx = matX(0,0);
                dBy = matX(1,0);
                dBz = matX(2,0);
                dga = matX(3,0);

                bx += dBx;
                by += dBy;
                bz += dBz;
                ga += dga;
            
                lbModel.fTemp[lbModel.DV_ZERO_ZERO_ZERO]     = lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_ZERO] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_ZERO] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_ZERO] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_ZERO]) ;
                lbModel.fTemp[lbModel.DV_ZERO_ZERO_P1  ]     = lbModel.wt[lbModel.DV_ZERO_ZERO_P1  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_P1  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_P1  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_P1  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_P1  ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_ZERO_P2  ]     = lbModel.wt[lbModel.DV_ZERO_ZERO_P2  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_P2  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_P2  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_P2  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_P2  ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_P2_ZERO  ]     = lbModel.wt[lbModel.DV_ZERO_P2_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P2_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_P2_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_P2_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_P2_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO  ]     = lbModel.wt[lbModel.DV_P2_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_P2_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_P2_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_P2_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_P2_ZERO_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_ZERO_M1  ]     = lbModel.wt[lbModel.DV_ZERO_ZERO_M1  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_M1  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_M1  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_M1  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_M1  ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_ZERO_M2  ]     = lbModel.wt[lbModel.DV_ZERO_ZERO_M2  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_ZERO_M2  ] + by*lbModel.cy[lbModel.DV_ZERO_ZERO_M2  ] + bz*lbModel.cz[lbModel.DV_ZERO_ZERO_M2  ] + ga*lbModel.cc[lbModel.DV_ZERO_ZERO_M2  ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_M2_ZERO  ]     = lbModel.wt[lbModel.DV_ZERO_M2_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M2_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_M2_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_M2_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_M2_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_M2_ZERO_ZERO  ]     = lbModel.wt[lbModel.DV_M2_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_M2_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_M2_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_M2_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_M2_ZERO_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_P1_P1    ]     = lbModel.wt[lbModel.DV_ZERO_P1_P1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P1_P1    ] + by*lbModel.cy[lbModel.DV_ZERO_P1_P1    ] + bz*lbModel.cz[lbModel.DV_ZERO_P1_P1    ] + ga*lbModel.cc[lbModel.DV_ZERO_P1_P1    ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_M1_P1    ]     = lbModel.wt[lbModel.DV_ZERO_M1_P1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M1_P1    ] + by*lbModel.cy[lbModel.DV_ZERO_M1_P1    ] + bz*lbModel.cz[lbModel.DV_ZERO_M1_P1    ] + ga*lbModel.cc[lbModel.DV_ZERO_M1_P1    ]) ;
                lbModel.fTemp[lbModel.DV_P1_ZERO_P1    ]     = lbModel.wt[lbModel.DV_P1_ZERO_P1    ]* exp( bx*lbModel.cx[lbModel.DV_P1_ZERO_P1    ] + by*lbModel.cy[lbModel.DV_P1_ZERO_P1    ] + bz*lbModel.cz[lbModel.DV_P1_ZERO_P1    ] + ga*lbModel.cc[lbModel.DV_P1_ZERO_P1    ]) ;
                lbModel.fTemp[lbModel.DV_M1_ZERO_P1    ]     = lbModel.wt[lbModel.DV_M1_ZERO_P1    ]* exp( bx*lbModel.cx[lbModel.DV_M1_ZERO_P1    ] + by*lbModel.cy[lbModel.DV_M1_ZERO_P1    ] + bz*lbModel.cz[lbModel.DV_M1_ZERO_P1    ] + ga*lbModel.cc[lbModel.DV_M1_ZERO_P1    ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_P1_M1    ]     = lbModel.wt[lbModel.DV_ZERO_P1_M1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P1_M1    ] + by*lbModel.cy[lbModel.DV_ZERO_P1_M1    ] + bz*lbModel.cz[lbModel.DV_ZERO_P1_M1    ] + ga*lbModel.cc[lbModel.DV_ZERO_P1_M1    ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_M1_M1    ]     = lbModel.wt[lbModel.DV_ZERO_M1_M1    ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M1_M1    ] + by*lbModel.cy[lbModel.DV_ZERO_M1_M1    ] + bz*lbModel.cz[lbModel.DV_ZERO_M1_M1    ] + ga*lbModel.cc[lbModel.DV_ZERO_M1_M1    ]) ;
                lbModel.fTemp[lbModel.DV_P1_ZERO_M1    ]     = lbModel.wt[lbModel.DV_P1_ZERO_M1    ]* exp( bx*lbModel.cx[lbModel.DV_P1_ZERO_M1    ] + by*lbModel.cy[lbModel.DV_P1_ZERO_M1    ] + bz*lbModel.cz[lbModel.DV_P1_ZERO_M1    ] + ga*lbModel.cc[lbModel.DV_P1_ZERO_M1    ]) ;
                lbModel.fTemp[lbModel.DV_M1_ZERO_M1    ]     = lbModel.wt[lbModel.DV_M1_ZERO_M1    ]* exp( bx*lbModel.cx[lbModel.DV_M1_ZERO_M1    ] + by*lbModel.cy[lbModel.DV_M1_ZERO_M1    ] + bz*lbModel.cz[lbModel.DV_M1_ZERO_M1    ] + ga*lbModel.cc[lbModel.DV_M1_ZERO_M1    ]) ;
                lbModel.fTemp[lbModel.DV_P1_P1_ZERO    ]     = lbModel.wt[lbModel.DV_P1_P1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_P1_P1_ZERO    ] + by*lbModel.cy[lbModel.DV_P1_P1_ZERO    ] + bz*lbModel.cz[lbModel.DV_P1_P1_ZERO    ] + ga*lbModel.cc[lbModel.DV_P1_P1_ZERO    ]) ;
                lbModel.fTemp[lbModel.DV_M1_P1_ZERO    ]     = lbModel.wt[lbModel.DV_M1_P1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_M1_P1_ZERO    ] + by*lbModel.cy[lbModel.DV_M1_P1_ZERO    ] + bz*lbModel.cz[lbModel.DV_M1_P1_ZERO    ] + ga*lbModel.cc[lbModel.DV_M1_P1_ZERO    ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_P1_ZERO  ]     = lbModel.wt[lbModel.DV_ZERO_P1_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_P1_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_P1_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_P1_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_P1_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_P1_ZERO_ZERO  ]     = lbModel.wt[lbModel.DV_P1_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_P1_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_P1_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_P1_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_P1_ZERO_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_P1_M1_ZERO    ]     = lbModel.wt[lbModel.DV_P1_M1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_P1_M1_ZERO    ] + by*lbModel.cy[lbModel.DV_P1_M1_ZERO    ] + bz*lbModel.cz[lbModel.DV_P1_M1_ZERO    ] + ga*lbModel.cc[lbModel.DV_P1_M1_ZERO    ]) ;
                lbModel.fTemp[lbModel.DV_M1_M1_ZERO    ]     = lbModel.wt[lbModel.DV_M1_M1_ZERO    ]* exp( bx*lbModel.cx[lbModel.DV_M1_M1_ZERO    ] + by*lbModel.cy[lbModel.DV_M1_M1_ZERO    ] + bz*lbModel.cz[lbModel.DV_M1_M1_ZERO    ] + ga*lbModel.cc[lbModel.DV_M1_M1_ZERO    ]) ;
                lbModel.fTemp[lbModel.DV_ZERO_M1_ZERO  ]     = lbModel.wt[lbModel.DV_ZERO_M1_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_ZERO_M1_ZERO  ] + by*lbModel.cy[lbModel.DV_ZERO_M1_ZERO  ] + bz*lbModel.cz[lbModel.DV_ZERO_M1_ZERO  ] + ga*lbModel.cc[lbModel.DV_ZERO_M1_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_M1_ZERO_ZERO  ]     = lbModel.wt[lbModel.DV_M1_ZERO_ZERO  ]* exp( bx*lbModel.cx[lbModel.DV_M1_ZERO_ZERO  ] + by*lbModel.cy[lbModel.DV_M1_ZERO_ZERO  ] + bz*lbModel.cz[lbModel.DV_M1_ZERO_ZERO  ] + ga*lbModel.cc[lbModel.DV_M1_ZERO_ZERO  ]) ;
                lbModel.fTemp[lbModel.DV_P_P_P         ]     = lbModel.wt[lbModel.DV_P_P_P         ]* exp( bx*lbModel.cx[lbModel.DV_P_P_P         ] + by*lbModel.cy[lbModel.DV_P_P_P         ] + bz*lbModel.cz[lbModel.DV_P_P_P         ] + ga*lbModel.cc[lbModel.DV_P_P_P         ]) ;
                lbModel.fTemp[lbModel.DV_M_P_P         ]     = lbModel.wt[lbModel.DV_M_P_P         ]* exp( bx*lbModel.cx[lbModel.DV_M_P_P         ] + by*lbModel.cy[lbModel.DV_M_P_P         ] + bz*lbModel.cz[lbModel.DV_M_P_P         ] + ga*lbModel.cc[lbModel.DV_M_P_P         ]) ;
                lbModel.fTemp[lbModel.DV_M_M_P         ]     = lbModel.wt[lbModel.DV_M_M_P         ]* exp( bx*lbModel.cx[lbModel.DV_M_M_P         ] + by*lbModel.cy[lbModel.DV_M_M_P         ] + bz*lbModel.cz[lbModel.DV_M_M_P         ] + ga*lbModel.cc[lbModel.DV_M_M_P         ]) ;
                lbModel.fTemp[lbModel.DV_P_M_P         ]     = lbModel.wt[lbModel.DV_P_M_P         ]* exp( bx*lbModel.cx[lbModel.DV_P_M_P         ] + by*lbModel.cy[lbModel.DV_P_M_P         ] + bz*lbModel.cz[lbModel.DV_P_M_P         ] + ga*lbModel.cc[lbModel.DV_P_M_P         ]) ;
                lbModel.fTemp[lbModel.DV_P_P_M         ]     = lbModel.wt[lbModel.DV_P_P_M         ]* exp( bx*lbModel.cx[lbModel.DV_P_P_M         ] + by*lbModel.cy[lbModel.DV_P_P_M         ] + bz*lbModel.cz[lbModel.DV_P_P_M         ] + ga*lbModel.cc[lbModel.DV_P_P_M         ]) ;
                lbModel.fTemp[lbModel.DV_M_P_M         ]     = lbModel.wt[lbModel.DV_M_P_M         ]* exp( bx*lbModel.cx[lbModel.DV_M_P_M         ] + by*lbModel.cy[lbModel.DV_M_P_M         ] + bz*lbModel.cz[lbModel.DV_M_P_M         ] + ga*lbModel.cc[lbModel.DV_M_P_M         ]) ;
                lbModel.fTemp[lbModel.DV_M_M_M         ]     = lbModel.wt[lbModel.DV_M_M_M         ]* exp( bx*lbModel.cx[lbModel.DV_M_M_M         ] + by*lbModel.cy[lbModel.DV_M_M_M         ] + bz*lbModel.cz[lbModel.DV_M_M_M         ] + ga*lbModel.cc[lbModel.DV_M_M_M         ]) ;
                lbModel.fTemp[lbModel.DV_P_M_M         ]     = lbModel.wt[lbModel.DV_P_M_M         ]* exp( bx*lbModel.cx[lbModel.DV_P_M_M         ] + by*lbModel.cy[lbModel.DV_P_M_M         ] + bz*lbModel.cz[lbModel.DV_P_M_M         ] + ga*lbModel.cc[lbModel.DV_P_M_M         ]) ;
                lbModel.fTemp[lbModel.DV_P1_P1_P1      ]     = lbModel.wt[lbModel.DV_P1_P1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_P1_P1      ] + by*lbModel.cy[lbModel.DV_P1_P1_P1      ] + bz*lbModel.cz[lbModel.DV_P1_P1_P1      ] + ga*lbModel.cc[lbModel.DV_P1_P1_P1      ]) ;
                lbModel.fTemp[lbModel.DV_M1_P1_P1      ]     = lbModel.wt[lbModel.DV_M1_P1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_P1_P1      ] + by*lbModel.cy[lbModel.DV_M1_P1_P1      ] + bz*lbModel.cz[lbModel.DV_M1_P1_P1      ] + ga*lbModel.cc[lbModel.DV_M1_P1_P1      ]) ;
                lbModel.fTemp[lbModel.DV_M1_M1_P1      ]     = lbModel.wt[lbModel.DV_M1_M1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_M1_P1      ] + by*lbModel.cy[lbModel.DV_M1_M1_P1      ] + bz*lbModel.cz[lbModel.DV_M1_M1_P1      ] + ga*lbModel.cc[lbModel.DV_M1_M1_P1      ]) ;
                lbModel.fTemp[lbModel.DV_P1_M1_P1      ]     = lbModel.wt[lbModel.DV_P1_M1_P1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_M1_P1      ] + by*lbModel.cy[lbModel.DV_P1_M1_P1      ] + bz*lbModel.cz[lbModel.DV_P1_M1_P1      ] + ga*lbModel.cc[lbModel.DV_P1_M1_P1      ]) ;
                lbModel.fTemp[lbModel.DV_P1_P1_M1      ]     = lbModel.wt[lbModel.DV_P1_P1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_P1_M1      ] + by*lbModel.cy[lbModel.DV_P1_P1_M1      ] + bz*lbModel.cz[lbModel.DV_P1_P1_M1      ] + ga*lbModel.cc[lbModel.DV_P1_P1_M1      ]) ;
                lbModel.fTemp[lbModel.DV_M1_P1_M1      ]     = lbModel.wt[lbModel.DV_M1_P1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_P1_M1      ] + by*lbModel.cy[lbModel.DV_M1_P1_M1      ] + bz*lbModel.cz[lbModel.DV_M1_P1_M1      ] + ga*lbModel.cc[lbModel.DV_M1_P1_M1      ]) ;
                lbModel.fTemp[lbModel.DV_M1_M1_M1      ]     = lbModel.wt[lbModel.DV_M1_M1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_M1_M1_M1      ] + by*lbModel.cy[lbModel.DV_M1_M1_M1      ] + bz*lbModel.cz[lbModel.DV_M1_M1_M1      ] + ga*lbModel.cc[lbModel.DV_M1_M1_M1      ]) ;
                lbModel.fTemp[lbModel.DV_P1_M1_M1      ]     = lbModel.wt[lbModel.DV_P1_M1_M1      ]* exp( bx*lbModel.cx[lbModel.DV_P1_M1_M1      ] + by*lbModel.cy[lbModel.DV_P1_M1_M1      ] + bz*lbModel.cz[lbModel.DV_P1_M1_M1      ] + ga*lbModel.cc[lbModel.DV_P1_M1_M1      ]) ;

            }    

                rho1 = lbModel.fTemp[lbModel.DV_ZERO_ZERO_ZERO]
                     + lbModel.fTemp[lbModel.DV_ZERO_ZERO_P1  ]
                     + lbModel.fTemp[lbModel.DV_ZERO_ZERO_P2  ]
                     + lbModel.fTemp[lbModel.DV_ZERO_P2_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_ZERO_ZERO_M1  ]
                     + lbModel.fTemp[lbModel.DV_ZERO_ZERO_M2  ]
                     + lbModel.fTemp[lbModel.DV_ZERO_M2_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_M2_ZERO_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_ZERO_P1_P1    ]
                     + lbModel.fTemp[lbModel.DV_ZERO_M1_P1    ]
                     + lbModel.fTemp[lbModel.DV_P1_ZERO_P1    ]
                     + lbModel.fTemp[lbModel.DV_M1_ZERO_P1    ]
                     + lbModel.fTemp[lbModel.DV_ZERO_P1_M1    ]
                     + lbModel.fTemp[lbModel.DV_ZERO_M1_M1    ]
                     + lbModel.fTemp[lbModel.DV_P1_ZERO_M1    ]
                     + lbModel.fTemp[lbModel.DV_M1_ZERO_M1    ]
                     + lbModel.fTemp[lbModel.DV_P1_P1_ZERO    ]
                     + lbModel.fTemp[lbModel.DV_M1_P1_ZERO    ]
                     + lbModel.fTemp[lbModel.DV_ZERO_P1_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_P1_ZERO_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_P1_M1_ZERO    ]
                     + lbModel.fTemp[lbModel.DV_M1_M1_ZERO    ]
                     + lbModel.fTemp[lbModel.DV_ZERO_M1_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_M1_ZERO_ZERO  ]
                     + lbModel.fTemp[lbModel.DV_P_P_P         ]
                     + lbModel.fTemp[lbModel.DV_M_P_P         ]
                     + lbModel.fTemp[lbModel.DV_M_M_P         ]
                     + lbModel.fTemp[lbModel.DV_P_M_P         ]
                     + lbModel.fTemp[lbModel.DV_P_P_M         ]
                     + lbModel.fTemp[lbModel.DV_M_P_M         ]
                     + lbModel.fTemp[lbModel.DV_M_M_M         ]
                     + lbModel.fTemp[lbModel.DV_P_M_M         ]
                     + lbModel.fTemp[lbModel.DV_P1_P1_P1      ]
                     + lbModel.fTemp[lbModel.DV_M1_P1_P1      ]
                     + lbModel.fTemp[lbModel.DV_M1_M1_P1      ]
                     + lbModel.fTemp[lbModel.DV_P1_M1_P1      ]
                     + lbModel.fTemp[lbModel.DV_P1_P1_M1      ]
                     + lbModel.fTemp[lbModel.DV_M1_P1_M1      ]
                     + lbModel.fTemp[lbModel.DV_M1_M1_M1      ]
                     + lbModel.fTemp[lbModel.DV_P1_M1_M1      ];
                
                
                
                double rhoRatio = rho[i]/rho1 ;  
                
                lbModel.fTemp0 [i][lbModel.CENTER_DV_ZERO_ZERO_ZERO] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_ZERO_ZERO];
                lbModel.fTemp1 [i][lbModel.G1_DV_ZERO_ZERO_P1      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_ZERO_P1  ];
                lbModel.fTemp1 [i][lbModel.G1_DV_ZERO_ZERO_P2      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_ZERO_P2  ];
                lbModel.fTemp1 [i][lbModel.G1_DV_ZERO_P2_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_P2_ZERO  ];
                lbModel.fTemp1 [i][lbModel.G1_DV_P2_ZERO_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO  ];
                lbModel.fTemp2 [i][lbModel.G2_DV_ZERO_ZERO_M1      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_ZERO_M1  ];
                lbModel.fTemp2 [i][lbModel.G2_DV_ZERO_ZERO_M2      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_ZERO_M2  ];
                lbModel.fTemp2 [i][lbModel.G2_DV_ZERO_M2_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_M2_ZERO  ];
                lbModel.fTemp2 [i][lbModel.G2_DV_M2_ZERO_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_M2_ZERO_ZERO  ];
                lbModel.fTemp3 [i][lbModel.G3_DV_ZERO_P1_P1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_P1_P1    ];
                lbModel.fTemp3 [i][lbModel.G3_DV_ZERO_M1_P1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_M1_P1    ];
                lbModel.fTemp3 [i][lbModel.G3_DV_P1_ZERO_P1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_ZERO_P1    ];
                lbModel.fTemp3 [i][lbModel.G3_DV_M1_ZERO_P1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_ZERO_P1    ];
                lbModel.fTemp4 [i][lbModel.G4_DV_ZERO_P1_M1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_P1_M1    ];
                lbModel.fTemp4 [i][lbModel.G4_DV_ZERO_M1_M1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_M1_M1    ];
                lbModel.fTemp4 [i][lbModel.G4_DV_P1_ZERO_M1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_ZERO_M1    ];
                lbModel.fTemp4 [i][lbModel.G4_DV_M1_ZERO_M1        ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_ZERO_M1    ];
                lbModel.fTemp5 [i][lbModel.G5_DV_P1_P1_ZERO        ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_P1_ZERO    ];
                lbModel.fTemp5 [i][lbModel.G5_DV_M1_P1_ZERO        ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_P1_ZERO    ];
                lbModel.fTemp5 [i][lbModel.G5_DV_ZERO_P1_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_P1_ZERO  ];
                lbModel.fTemp5 [i][lbModel.G5_DV_P1_ZERO_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_ZERO_ZERO  ];
                lbModel.fTemp6 [i][lbModel.G6_DV_P1_M1_ZERO        ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_M1_ZERO    ];
                lbModel.fTemp6 [i][lbModel.G6_DV_M1_M1_ZERO        ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_M1_ZERO    ];
                lbModel.fTemp6 [i][lbModel.G6_DV_ZERO_M1_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_ZERO_M1_ZERO  ];
                lbModel.fTemp6 [i][lbModel.G6_DV_M1_ZERO_ZERO      ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_ZERO_ZERO  ];
                lbModel.fTemp7 [i][lbModel.G7_DV_P_P_P             ] = rhoRatio * lbModel.fTemp[lbModel.DV_P_P_P         ];
                lbModel.fTemp7 [i][lbModel.G7_DV_M_P_P             ] = rhoRatio * lbModel.fTemp[lbModel.DV_M_P_P         ];
                lbModel.fTemp7 [i][lbModel.G7_DV_M_M_P             ] = rhoRatio * lbModel.fTemp[lbModel.DV_M_M_P         ];
                lbModel.fTemp7 [i][lbModel.G7_DV_P_M_P             ] = rhoRatio * lbModel.fTemp[lbModel.DV_P_M_P         ];
                lbModel.fTemp8 [i][lbModel.G8_DV_P_P_M             ] = rhoRatio * lbModel.fTemp[lbModel.DV_P_P_M         ];
                lbModel.fTemp8 [i][lbModel.G8_DV_M_P_M             ] = rhoRatio * lbModel.fTemp[lbModel.DV_M_P_M         ];
                lbModel.fTemp8 [i][lbModel.G8_DV_M_M_M             ] = rhoRatio * lbModel.fTemp[lbModel.DV_M_M_M         ];
                lbModel.fTemp8 [i][lbModel.G8_DV_P_M_M             ] = rhoRatio * lbModel.fTemp[lbModel.DV_P_M_M         ];
                lbModel.fTemp9 [i][lbModel.G9_DV_P1_P1_P1          ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_P1_P1      ];
                lbModel.fTemp9 [i][lbModel.G9_DV_M1_P1_P1          ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_P1_P1      ];
                lbModel.fTemp9 [i][lbModel.G9_DV_M1_M1_P1          ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_M1_P1      ];
                lbModel.fTemp9 [i][lbModel.G9_DV_P1_M1_P1          ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_M1_P1      ];
                lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1         ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_P1_M1      ];
                lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1         ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_P1_M1      ];
                lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1         ] = rhoRatio * lbModel.fTemp[lbModel.DV_M1_M1_M1      ];
                lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1         ] = rhoRatio * lbModel.fTemp[lbModel.DV_P1_M1_M1      ];             
                
        }
    }



*/


// Routines from Mr.Nilesh Sawant, Sankhyasutra labs, Bangalore
template<class T>
void solve3By3(T (&A)[3][3],T (&b)[3],T (&x)[3])
{
 T detInv,cA[3][3];
 
 cA[0][0]= A[1][1]*A[2][2] - A[1][2]*A[2][1]; 
 cA[0][1]=-A[1][0]*A[2][2] + A[1][2]*A[2][0]; 
 cA[0][2]= A[1][0]*A[2][1] - A[1][1]*A[2][0];
           
 cA[1][0]=-A[0][1]*A[2][2] + A[0][2]*A[2][1]; 
 cA[1][1]= A[0][0]*A[2][2] - A[0][2]*A[2][0];
 cA[1][2]=-A[0][0]*A[2][1] + A[0][1]*A[2][0];
           
 cA[2][0]= A[0][1]*A[1][2] - A[0][2]*A[1][1];
 cA[2][1]=-A[0][0]*A[1][2] + A[0][2]*A[1][0];
 cA[2][2]= A[0][0]*A[1][1] - A[0][1]*A[1][0];
 
 detInv = 1.0/(A[0][0]*cA[0][0] + A[0][1]*cA[0][1] + A[0][2]*cA[0][2]);
 
 x[0]=(cA[0][0]*b[0]+cA[1][0]*b[1]+cA[2][0]*b[2])*detInv;
 x[1]=(cA[0][1]*b[0]+cA[1][1]*b[1]+cA[2][1]*b[2])*detInv;
 x[2]=(cA[0][2]*b[0]+cA[1][2]*b[1]+cA[2][2]*b[2])*detInv;
}

template<class T>
void solve3By3Vec(T (&A)[3][3],T (&b)[3],T (&x)[3])
{
T detInv,cA[3][3] __attribute__ ((aligned(__BIGGEST_ALIGNMENT__)));

#pragma GCC ivdep
for (int i=0;i<3;i++) for (int j=0;j<3;j++) cA[i][j] = A[(i+1)%3][(j+1)%3]*A[(i+2)%3][(j+2)%3]-A[(i+1)%3][(j+2)%3]*A[(i+2)%3][(j+1)%3];

detInv = 1.0/(A[0][0]*cA[0][0] + A[0][1]*cA[0][1] + A[0][2]*cA[0][2]);

#pragma GCC ivdep
for (int i=0;i<3;i++) x[i]=(cA[0][i]*b[0]+cA[1][i]*b[1]+cA[2][i]*b[2])*detInv;

}



// template<typename dType,int dvN>
// void calculateWeightdot(lbmRD3Q41<dType> &lbModel,dType (&fEq)[dvN],dType &bx,dType &by,dType &bz,dType &ga)
// {
//   dType ex(exp(bx)),ey(exp(by)),ez(exp(bz)),eg(exp(ga));
//   dType sqrt_eg(sqrt(eg));
//   
//   dType exSC1(ex),eySC1(ey),ezSC1(ez),egSC1(eg);
//   dType exSC2(exSC1*exSC1),eySC2(eySC1*eySC1),ezSC2(ezSC1*ezSC1),egSC2(eg*eg*eg*eg);
//   
//   dType exFCC1(exSC1),eyFCC1(eySC1),ezFCC1(ezSC1),egFCC1(eg*eg);
//   
//   dType exBCC(sqrt(ex)),eyBCC(sqrt(ey)),ezBCC(sqrt(ez)),egBCC(sqrt(sqrt_eg)*sqrt_eg);
//   dType exBCC1(exSC1),eyBCC1(eySC1),ezBCC1(ezSC1),egBCC1(eg*egFCC1);
//   
//   lbModel.fTemp[lbModel.DV_ZERO_ZERO_ZERO]    = lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO] ;
//   lbModel.fTemp[lbModel.DV_ZERO_ZERO_P1  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_P1  ] * ezSC1   * egSC1;
//   lbModel.fTemp[lbModel.DV_ZERO_ZERO_P2  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_P2  ] * ezSC2   * egSC2;
//   lbModel.fTemp[lbModel.DV_ZERO_P2_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_P2_ZERO  ] * eySC2   * egSC2;
//   lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_P2_ZERO_ZERO  ] * exSC2   * egSC2;
//   lbModel.fTemp[lbModel.DV_ZERO_ZERO_M1  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_M1  ] / ezSC1   * egSC1;
//   lbModel.fTemp[lbModel.DV_ZERO_ZERO_M2  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_M2  ] / ezSC2   * egSC2;
//   lbModel.fTemp[lbModel.DV_ZERO_M2_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_M2_ZERO  ] / eySC2   * egSC2;
//   lbModel.fTemp[lbModel.DV_M2_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_M2_ZERO_ZERO  ] / exSC2   * egSC2;
//   lbModel.fTemp[lbModel.DV_ZERO_P1_P1    ]    = lbModel.wt[lbModel.DV_ZERO_P1_P1    ] * eyFCC1  * ezFCC1  * egFCC1;
//   lbModel.fTemp[lbModel.DV_ZERO_M1_P1    ]    = lbModel.wt[lbModel.DV_ZERO_M1_P1    ] / eyFCC1  * ezFCC1  * egFCC1;  
//   lbModel.fTemp[lbModel.DV_P1_ZERO_P1    ]    = lbModel.wt[lbModel.DV_P1_ZERO_P1    ] * exFCC1  * ezFCC1  * egFCC1;   
//   lbModel.fTemp[lbModel.DV_M1_ZERO_P1    ]    = lbModel.wt[lbModel.DV_M1_ZERO_P1    ] / exFCC1  * ezFCC1  * egFCC1;    
//   lbModel.fTemp[lbModel.DV_ZERO_P1_M1    ]    = lbModel.wt[lbModel.DV_ZERO_P1_M1    ] * eyFCC1  / ezFCC1  * egFCC1;     
//   lbModel.fTemp[lbModel.DV_ZERO_M1_M1    ]    = lbModel.wt[lbModel.DV_ZERO_M1_M1    ] / eyFCC1  / ezFCC1  * egFCC1;      
//   lbModel.fTemp[lbModel.DV_P1_ZERO_M1    ]    = lbModel.wt[lbModel.DV_P1_ZERO_M1    ] * exFCC1  / ezFCC1  * egFCC1;  
//   lbModel.fTemp[lbModel.DV_M1_ZERO_M1    ]    = lbModel.wt[lbModel.DV_M1_ZERO_M1    ] / exFCC1  / ezFCC1  * egFCC1;   
//   lbModel.fTemp[lbModel.DV_P1_P1_ZERO    ]    = lbModel.wt[lbModel.DV_P1_P1_ZERO    ] * exFCC1  * eyFCC1  * egFCC1;   
//   lbModel.fTemp[lbModel.DV_M1_P1_ZERO    ]    = lbModel.wt[lbModel.DV_M1_P1_ZERO    ] / exFCC1  * eyFCC1  * egFCC1;    
//   lbModel.fTemp[lbModel.DV_ZERO_P1_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_P1_ZERO  ] * eyFCC1  * egFCC1;     
//   lbModel.fTemp[lbModel.DV_P1_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_P1_ZERO_ZERO  ] * exFCC1  * egFCC1;     
//   lbModel.fTemp[lbModel.DV_P1_M1_ZERO    ]    = lbModel.wt[lbModel.DV_P1_M1_ZERO    ] * exFCC1  / eyFCC1  * egFCC1;     
//   lbModel.fTemp[lbModel.DV_M1_M1_ZERO    ]    = lbModel.wt[lbModel.DV_M1_M1_ZERO    ] / exFCC1  / eyFCC1  * egFCC1;      
//   lbModel.fTemp[lbModel.DV_ZERO_M1_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_M1_ZERO  ] / eyFCC1  * egSC1;       
//   lbModel.fTemp[lbModel.DV_M1_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_M1_ZERO_ZERO  ] / exFCC1  * egSC1; 
//   lbModel.fTemp[lbModel.DV_P_P_P         ]    = lbModel.wt[lbModel.DV_P_P_P         ] * exBCC   * eyBCC   * ezBCC   * egBCC; 
//   lbModel.fTemp[lbModel.DV_M_P_P         ]    = lbModel.wt[lbModel.DV_M_P_P         ] / exBCC   * eyBCC   * ezBCC   * egBCC; 
//   lbModel.fTemp[lbModel.DV_M_M_P         ]    = lbModel.wt[lbModel.DV_M_M_P         ] / exBCC   / eyBCC   * ezBCC   * egBCC; 
//   lbModel.fTemp[lbModel.DV_P_M_P         ]    = lbModel.wt[lbModel.DV_P_M_P         ] * exBCC   / eyBCC   * ezBCC   * egBCC; 
//   lbModel.fTemp[lbModel.DV_P_P_M         ]    = lbModel.wt[lbModel.DV_P_P_M         ] * exBCC   * eyBCC   / ezBCC   * egBCC;  
//   lbModel.fTemp[lbModel.DV_M_P_M         ]    = lbModel.wt[lbModel.DV_M_P_M         ] / exBCC   * eyBCC   * ezBCC   * egBCC;  
//   lbModel.fTemp[lbModel.DV_M_M_M         ]    = lbModel.wt[lbModel.DV_M_M_M         ] / exBCC   / eyBCC   / ezBCC   * egBCC;  
//   lbModel.fTemp[lbModel.DV_P_M_M         ]    = lbModel.wt[lbModel.DV_P_M_M         ] * exBCC   / eyBCC   / ezBCC   * egBCC;  
//   lbModel.fTemp[lbModel.DV_P1_P1_P1      ]    = lbModel.wt[lbModel.DV_P1_P1_P1      ] * exBCC1  * eyBCC1  * ezBCC1  * egBCC1; 
//   lbModel.fTemp[lbModel.DV_M1_P1_P1      ]    = lbModel.wt[lbModel.DV_M1_P1_P1      ] / exBCC1  * eyBCC1  * ezBCC1  * egBCC1; 
//   lbModel.fTemp[lbModel.DV_M1_M1_P1      ]    = lbModel.wt[lbModel.DV_M1_M1_P1      ] / exBCC1  / eyBCC1  * ezBCC1  * egBCC1; 
//   lbModel.fTemp[lbModel.DV_P1_M1_P1      ]    = lbModel.wt[lbModel.DV_P1_M1_P1      ] * exBCC1  / eyBCC1  * ezBCC1  * egBCC1; 
//   lbModel.fTemp[lbModel.DV_P1_P1_M1      ]    = lbModel.wt[lbModel.DV_P1_P1_M1      ] * exBCC1  * eyBCC1  / ezBCC1  * egBCC1;  
//   lbModel.fTemp[lbModel.DV_M1_P1_M1      ]    = lbModel.wt[lbModel.DV_M1_P1_M1      ] / exBCC1  * eyBCC1  * ezBCC1  * egBCC1;  
//   lbModel.fTemp[lbModel.DV_M1_M1_M1      ]    = lbModel.wt[lbModel.DV_M1_M1_M1      ] / exBCC1  / eyBCC1  / ezBCC1  * egBCC1;  
//   lbModel.fTemp[lbModel.DV_P1_M1_M1      ]    = lbModel.wt[lbModel.DV_P1_M1_M1      ] * exBCC1  / eyBCC1  / ezBCC1  * egBCC1;  
// };
// 


template<typename dType>
void getFEqNewIterative(lbmRD3Q41<dType> &lbModel,int VECT_LENGTH, dType *rho, dType *uX, dType *uY, dType *uZ, dType *theta)
{
 dType matA[3][3]__attribute__ ((aligned(__BIGGEST_ALIGNMENT__))),matB[3]__attribute__ ((aligned(__BIGGEST_ALIGNMENT__))),matX[3]__attribute__ ((aligned(__BIGGEST_ALIGNMENT__)));
 dType qxx,qyy,qzz,qxxE,qyyE,qzzE,EE,rEInv;
 dType rho1(0.0),Pxy(0.0),Pxz(0.0),Pyz(0.0),Pxx(0.0),Pyy(0.0),Pzz(0.0),jxc(0.0),jyc(0.0),jzc(0.0),qx(0.0),qy(0.0),qz(0.0),ec(0.0),r(0.0);

 for(int i=0;i<VECT_LENGTH;i++)
 {  
  dType qxx,qyy,qzz,qxxE,qyyE,qzzE,EE,rEInv;
  dType rho1(0.0),Pxy(0.0),Pxz(0.0),Pyz(0.0),Pxx(0.0),Pyy(0.0),Pzz(0.0),jxc(0.0),jyc(0.0),jzc(0.0),qx(0.0),qy(0.0),qz(0.0),ec(0.0),r(0.0);

  dType E          = uX[i]*uX[i]+uY[i]*uY[i]+uZ[i]*uZ[i] + 3.0*theta[i];
  dType oneByTheta = 1.0/theta[i];
  dType bx         = uX[i]*oneByTheta;
  dType by         = uY[i]*oneByTheta;
  dType bz         = uZ[i]*oneByTheta;
  dType del_theta  = (theta[i]-lbModel.theta0)*lbModel.oneByTheta0 ;
  dType ga         = (0.5*lbModel.oneByTheta0)*del_theta*(1.0/(1.0+del_theta));

  dType ex         = exp(bx);
  dType ey         = exp(by);
  dType ez         = exp(bz);
  dType exSC2      = ex*ex;
  dType eySC2      = ey*ey;
  dType ezSC2      = ez*ez;
  dType exBCC      = sqrt(ex);
  dType eyBCC      = sqrt(ey);
  dType ezBCC      = sqrt(ez);

  dType eg         = exp(ga);
  dType sqrt_eg    = sqrt(eg);
  dType egSC1      = eg;
  dType egFCC1     = eg*eg;
  dType egBCC1     = eg*egFCC1;
  dType egSC2      = egFCC1*egFCC1;
  dType egBCC      = sqrt(sqrt_eg)*sqrt_eg;

  lbModel.fTemp[lbModel.DV_ZERO_ZERO_ZERO]    = lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO] ;
  lbModel.fTemp[lbModel.DV_ZERO_ZERO_P1  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_P1  ] * ez                        * egSC1;
  lbModel.fTemp[lbModel.DV_ZERO_ZERO_P2  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_P2  ] * ezSC2                     * egSC2;
  lbModel.fTemp[lbModel.DV_ZERO_P2_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_P2_ZERO  ] * eySC2                     * egSC2;
  lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_P2_ZERO_ZERO  ] * exSC2                     * egSC2;
  lbModel.fTemp[lbModel.DV_ZERO_ZERO_M1  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_M1  ] / ez                        * egSC1;
  lbModel.fTemp[lbModel.DV_ZERO_ZERO_M2  ]    = lbModel.wt[lbModel.DV_ZERO_ZERO_M2  ] / ezSC2                     * egSC2;
  lbModel.fTemp[lbModel.DV_ZERO_M2_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_M2_ZERO  ] / eySC2                     * egSC2;
  lbModel.fTemp[lbModel.DV_M2_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_M2_ZERO_ZERO  ] / exSC2                     * egSC2;           
  lbModel.fTemp[lbModel.DV_ZERO_P1_P1    ]    = lbModel.wt[lbModel.DV_ZERO_P1_P1    ] * ey     * ez               * egFCC1;
  lbModel.fTemp[lbModel.DV_ZERO_M1_P1    ]    = lbModel.wt[lbModel.DV_ZERO_M1_P1    ] / ey     * ez               * egFCC1;  
  lbModel.fTemp[lbModel.DV_P1_ZERO_P1    ]    = lbModel.wt[lbModel.DV_P1_ZERO_P1    ] * ex     * ez               * egFCC1;   
  lbModel.fTemp[lbModel.DV_M1_ZERO_P1    ]    = lbModel.wt[lbModel.DV_M1_ZERO_P1    ] / ex     * ez               * egFCC1;    
  lbModel.fTemp[lbModel.DV_ZERO_P1_M1    ]    = lbModel.wt[lbModel.DV_ZERO_P1_M1    ] * ey     / ez               * egFCC1;     
  lbModel.fTemp[lbModel.DV_ZERO_M1_M1    ]    = lbModel.wt[lbModel.DV_ZERO_M1_M1    ] / ey     / ez               * egFCC1;      
  lbModel.fTemp[lbModel.DV_P1_ZERO_M1    ]    = lbModel.wt[lbModel.DV_P1_ZERO_M1    ] * ex     / ez               * egFCC1;  
  lbModel.fTemp[lbModel.DV_M1_ZERO_M1    ]    = lbModel.wt[lbModel.DV_M1_ZERO_M1    ] / ex     / ez               * egFCC1;   
  lbModel.fTemp[lbModel.DV_P1_P1_ZERO    ]    = lbModel.wt[lbModel.DV_P1_P1_ZERO    ] * ex     * ey               * egFCC1;   
  lbModel.fTemp[lbModel.DV_M1_P1_ZERO    ]    = lbModel.wt[lbModel.DV_M1_P1_ZERO    ] / ex     * ey               * egFCC1;    
  lbModel.fTemp[lbModel.DV_ZERO_P1_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_P1_ZERO  ] * ey                        * egSC1;               
  lbModel.fTemp[lbModel.DV_P1_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_P1_ZERO_ZERO  ] * ex                        * egSC1;               
  lbModel.fTemp[lbModel.DV_P1_M1_ZERO    ]    = lbModel.wt[lbModel.DV_P1_M1_ZERO    ] * ex     / ey               * egFCC1;     
  lbModel.fTemp[lbModel.DV_M1_M1_ZERO    ]    = lbModel.wt[lbModel.DV_M1_M1_ZERO    ] / ex     / ey               * egFCC1;      
  lbModel.fTemp[lbModel.DV_ZERO_M1_ZERO  ]    = lbModel.wt[lbModel.DV_ZERO_M1_ZERO  ] / ey                        * egSC1;                 
  lbModel.fTemp[lbModel.DV_M1_ZERO_ZERO  ]    = lbModel.wt[lbModel.DV_M1_ZERO_ZERO  ] / ex                        * egSC1; 
  lbModel.fTemp[lbModel.DV_P_P_P         ]    = lbModel.wt[lbModel.DV_P_P_P         ] * exBCC  * eyBCC  * ezBCC   * egBCC; 
  lbModel.fTemp[lbModel.DV_M_P_P         ]    = lbModel.wt[lbModel.DV_M_P_P         ] / exBCC  * eyBCC  * ezBCC   * egBCC; 
  lbModel.fTemp[lbModel.DV_M_M_P         ]    = lbModel.wt[lbModel.DV_M_M_P         ] / exBCC  / eyBCC  * ezBCC   * egBCC; 
  lbModel.fTemp[lbModel.DV_P_M_P         ]    = lbModel.wt[lbModel.DV_P_M_P         ] * exBCC  / eyBCC  * ezBCC   * egBCC; 
  lbModel.fTemp[lbModel.DV_P_P_M         ]    = lbModel.wt[lbModel.DV_P_P_M         ] * exBCC  * eyBCC  / ezBCC   * egBCC;  
  lbModel.fTemp[lbModel.DV_M_P_M         ]    = lbModel.wt[lbModel.DV_M_P_M         ] / exBCC  * eyBCC  * ezBCC   * egBCC;  
  lbModel.fTemp[lbModel.DV_M_M_M         ]    = lbModel.wt[lbModel.DV_M_M_M         ] / exBCC  / eyBCC  / ezBCC   * egBCC;  
  lbModel.fTemp[lbModel.DV_P_M_M         ]    = lbModel.wt[lbModel.DV_P_M_M         ] * exBCC  / eyBCC  / ezBCC   * egBCC;  
  lbModel.fTemp[lbModel.DV_P1_P1_P1      ]    = lbModel.wt[lbModel.DV_P1_P1_P1      ] * ex     * ey     * ez      * egBCC1; 
  lbModel.fTemp[lbModel.DV_M1_P1_P1      ]    = lbModel.wt[lbModel.DV_M1_P1_P1      ] / ex     * ey     * ez      * egBCC1; 
  lbModel.fTemp[lbModel.DV_M1_M1_P1      ]    = lbModel.wt[lbModel.DV_M1_M1_P1      ] / ex     / ey     * ez      * egBCC1; 
  lbModel.fTemp[lbModel.DV_P1_M1_P1      ]    = lbModel.wt[lbModel.DV_P1_M1_P1      ] * ex     / ey     * ez      * egBCC1; 
  lbModel.fTemp[lbModel.DV_P1_P1_M1      ]    = lbModel.wt[lbModel.DV_P1_P1_M1      ] * ex     * ey     / ez      * egBCC1;  
  lbModel.fTemp[lbModel.DV_M1_P1_M1      ]    = lbModel.wt[lbModel.DV_M1_P1_M1      ] / ex     * ey     * ez      * egBCC1;  
  lbModel.fTemp[lbModel.DV_M1_M1_M1      ]    = lbModel.wt[lbModel.DV_M1_M1_M1      ] / ex     / ey     / ez      * egBCC1;  
  lbModel.fTemp[lbModel.DV_P1_M1_M1      ]    = lbModel.wt[lbModel.DV_P1_M1_M1      ] * ex     / ey     / ez      * egBCC1;  

  
  for (int dv=1;dv<lbModel.dvN;dv++) 
  {
   rho1 += lbModel.fTemp[dv] ;
   jxc  += lbModel.fTemp[dv]*lbModel.cx[dv] ;
   jyc  += lbModel.fTemp[dv]*lbModel.cy[dv] ;
   jzc  += lbModel.fTemp[dv]*lbModel.cz[dv] ;
   Pxx  += lbModel.fTemp[dv]*lbModel.cx2[dv] ;
   Pyy  += lbModel.fTemp[dv]*lbModel.cy2[dv] ;
   Pzz  += lbModel.fTemp[dv]*lbModel.cz2[dv] ;
   Pxy  += lbModel.fTemp[dv]*lbModel.cxcy[dv] ;
   Pxz  += lbModel.fTemp[dv]*lbModel.czcx[dv] ;
   Pyz  += lbModel.fTemp[dv]*lbModel.cycz[dv] ;
   qx   += lbModel.fTemp[dv]*lbModel.cxcsq[dv];
   qy   += lbModel.fTemp[dv]*lbModel.cycsq[dv];
   qz   += lbModel.fTemp[dv]*lbModel.czcsq[dv];
   r    += lbModel.fTemp[dv]*lbModel.csq2[dv] ;
  }
  
  rho1 += lbModel.fTemp[0];
  ec    = Pxx+Pyy+Pzz;
        
  EE    = (rho1*E-ec);
  rEInv = 1.0/(r-ec*E);
  qxx   = (qx-ec*uX[i])*rEInv;
  qyy   = (qy-ec*uY[i])*rEInv;
  qzz   = (qz-ec*uZ[i])*rEInv;
  qxxE  = (qx-jxc*E);
  qyyE  = (qy-jyc*E);
  qzzE  = (qz-jzc*E);
  
  matA[0][0] = Pxx -uX[i]*jxc-qxxE*qxx;
  matA[0][1] = Pxy -uX[i]*jyc-qyyE*qxx;
  matA[0][2] = Pxz -uX[i]*jzc-qzzE*qxx;
               
  matA[1][0] = Pxy -uY[i]*jxc-qxxE*qyy;
  matA[1][1] = Pyy -uY[i]*jyc-qyyE*qyy;
  matA[1][2] = Pyz -uY[i]*jzc-qzzE*qyy;
               
  matA[2][0] = Pxz -uZ[i]*jxc-qxxE*qzz;
  matA[2][1] = Pyz -uZ[i]*jyc-qyyE*qzz;
  matA[2][2] = Pzz -uZ[i]*jzc-qzzE*qzz;
               
  matB[0]    = rho1*uX[i]-jxc-EE*qxx;
  matB[1]    = rho1*uY[i]-jyc-EE*qyy;
  matB[2]    = rho1*uZ[i]-jzc-EE*qzz;
  
  solve3By3Vec(matA,matB,matX);
  
  dType dga=(EE-matX[0]*qxxE-matX[1]*qyyE-matX[2]*qzzE)*rEInv;
  
  rho1=rho[i]/(rho1+matX[0]*jxc+matX[1]*jyc+matX[2]*jzc+dga*ec);

  dType dgaSC1Plus1  = 1.0 + dga;
  dType dgaSC2Plus1  = 1.0 + dga*4.0;
  dType dgaFCC1Plus1 = 1.0 + dga*2.0;
  dType dgaBCCPlus1  = 1.0 + dga*0.75;
  dType dgaBCC1Plus1 = 1.0 + dga*3.0;
    
  lbModel.fTemp0 [i][lbModel.CENTER_DV_ZERO_ZERO_ZERO] = lbModel.fTemp[lbModel.DV_ZERO_ZERO_ZERO]   * rho1 ;
  lbModel.fTemp1 [i][lbModel.G1_DV_ZERO_ZERO_P1      ] = lbModel.fTemp[lbModel.DV_ZERO_ZERO_P1  ]   * rho1 *(                                    matX[2]      + dgaSC1Plus1 )  ;
  lbModel.fTemp1 [i][lbModel.G1_DV_ZERO_ZERO_P2      ] = lbModel.fTemp[lbModel.DV_ZERO_ZERO_P2  ]   * rho1 *(                               +2.0*matX[2]      + dgaSC2Plus1 )  ;
  lbModel.fTemp1 [i][lbModel.G1_DV_ZERO_P2_ZERO      ] = lbModel.fTemp[lbModel.DV_ZERO_P2_ZERO  ]   * rho1 *(               2.0* matX[1]                      + dgaSC2Plus1 )  ;
  lbModel.fTemp1 [i][lbModel.G1_DV_P2_ZERO_ZERO      ] = lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO  ]   * rho1 *( 2.0*matX[0]                                     + dgaSC2Plus1 )  ;
  lbModel.fTemp2 [i][lbModel.G2_DV_ZERO_ZERO_M1      ] = lbModel.fTemp[lbModel.DV_ZERO_ZERO_M1  ]   * rho1 *(                               -    matX[2]      + dgaSC1Plus1 )  ;
  lbModel.fTemp2 [i][lbModel.G2_DV_ZERO_ZERO_M2      ] = lbModel.fTemp[lbModel.DV_ZERO_ZERO_M2  ]   * rho1 *(                               -2.0*matX[2]      + dgaSC2Plus1 )  ;
  lbModel.fTemp2 [i][lbModel.G2_DV_ZERO_M2_ZERO      ] = lbModel.fTemp[lbModel.DV_ZERO_M2_ZERO  ]   * rho1 *(              -2.0* matX[1]                      + dgaSC2Plus1 )  ;
  lbModel.fTemp2 [i][lbModel.G2_DV_M2_ZERO_ZERO      ] = lbModel.fTemp[lbModel.DV_M2_ZERO_ZERO  ]   * rho1 *( -2.0*matX[0]                                    + dgaSC2Plus1 )  ;
  lbModel.fTemp3 [i][lbModel.G3_DV_ZERO_P1_P1        ] = lbModel.fTemp[lbModel.DV_ZERO_P1_P1    ]   * rho1 *(                    matX[1]    +    matX[2]      + dgaFCC1Plus1)  ;
  lbModel.fTemp3 [i][lbModel.G3_DV_ZERO_M1_P1        ] = lbModel.fTemp[lbModel.DV_ZERO_M1_P1    ]   * rho1 *(             -      matX[1]    +    matX[2]      + dgaFCC1Plus1)  ;
  lbModel.fTemp3 [i][lbModel.G3_DV_P1_ZERO_P1        ] = lbModel.fTemp[lbModel.DV_P1_ZERO_P1    ]   * rho1 *(  matX[0]                      +    matX[2]      + dgaFCC1Plus1)  ;
  lbModel.fTemp3 [i][lbModel.G3_DV_M1_ZERO_P1        ] = lbModel.fTemp[lbModel.DV_M1_ZERO_P1    ]   * rho1 *( -matX[0]                      +    matX[2]      + dgaFCC1Plus1)  ;
  lbModel.fTemp4 [i][lbModel.G4_DV_ZERO_P1_M1        ] = lbModel.fTemp[lbModel.DV_ZERO_P1_M1    ]   * rho1 *(                    matX[1]    -    matX[2]      + dgaFCC1Plus1)  ;
  lbModel.fTemp4 [i][lbModel.G4_DV_ZERO_M1_M1        ] = lbModel.fTemp[lbModel.DV_ZERO_M1_M1    ]   * rho1 *(             -      matX[1]    -    matX[2]      + dgaFCC1Plus1)  ;
  lbModel.fTemp4 [i][lbModel.G4_DV_P1_ZERO_M1        ] = lbModel.fTemp[lbModel.DV_P1_ZERO_M1    ]   * rho1 *(  matX[0]                      -    matX[2]      + dgaFCC1Plus1)  ;
  lbModel.fTemp4 [i][lbModel.G4_DV_M1_ZERO_M1        ] = lbModel.fTemp[lbModel.DV_M1_ZERO_M1    ]   * rho1 *( -matX[0]                      -    matX[2]      + dgaFCC1Plus1)  ;
  lbModel.fTemp5 [i][lbModel.G5_DV_P1_P1_ZERO        ] = lbModel.fTemp[lbModel.DV_P1_P1_ZERO    ]   * rho1 *(  matX[0]    +      matX[1]                      + dgaFCC1Plus1)  ;
  lbModel.fTemp5 [i][lbModel.G5_DV_M1_P1_ZERO        ] = lbModel.fTemp[lbModel.DV_M1_P1_ZERO    ]   * rho1 *( -matX[0]    +      matX[1]                      + dgaFCC1Plus1)  ;
  lbModel.fTemp5 [i][lbModel.G5_DV_ZERO_P1_ZERO      ] = lbModel.fTemp[lbModel.DV_ZERO_P1_ZERO  ]   * rho1 *(                    matX[1]                      + dgaSC1Plus1 )  ;
  lbModel.fTemp5 [i][lbModel.G5_DV_P1_ZERO_ZERO      ] = lbModel.fTemp[lbModel.DV_P1_ZERO_ZERO  ]   * rho1 *(  matX[0]                                        + dgaSC1Plus1 )  ;
  lbModel.fTemp6 [i][lbModel.G6_DV_P1_M1_ZERO        ] = lbModel.fTemp[lbModel.DV_P1_M1_ZERO    ]   * rho1 *(  matX[0]    -      matX[1]                      + dgaFCC1Plus1)  ;
  lbModel.fTemp6 [i][lbModel.G6_DV_M1_M1_ZERO        ] = lbModel.fTemp[lbModel.DV_M1_M1_ZERO    ]   * rho1 *( -matX[0]    -      matX[1]                      + dgaFCC1Plus1)  ;
  lbModel.fTemp6 [i][lbModel.G6_DV_ZERO_M1_ZERO      ] = lbModel.fTemp[lbModel.DV_ZERO_M1_ZERO  ]   * rho1 *(             -      matX[1]                      + dgaSC1Plus1 )  ;
  lbModel.fTemp6 [i][lbModel.G6_DV_M1_ZERO_ZERO      ] = lbModel.fTemp[lbModel.DV_M1_ZERO_ZERO  ]   * rho1 *( -matX[0]                                        + dgaSC1Plus1 )  ;
  lbModel.fTemp7 [i][lbModel.G7_DV_P_P_P             ] = lbModel.fTemp[lbModel.DV_P_P_P         ]   * rho1 *(( matX[0]    +      matX[1]    +    matX[2])*0.5 + dgaBCCPlus1 )  ;
  lbModel.fTemp7 [i][lbModel.G7_DV_M_P_P             ] = lbModel.fTemp[lbModel.DV_M_P_P         ]   * rho1 *((-matX[0]    +      matX[1]    +    matX[2])*0.5 + dgaBCCPlus1 )  ;
  lbModel.fTemp7 [i][lbModel.G7_DV_M_M_P             ] = lbModel.fTemp[lbModel.DV_M_M_P         ]   * rho1 *((-matX[0]    -      matX[1]    +    matX[2])*0.5 + dgaBCCPlus1 )  ;
  lbModel.fTemp7 [i][lbModel.G7_DV_P_M_P             ] = lbModel.fTemp[lbModel.DV_P_M_P         ]   * rho1 *(( matX[0]    -      matX[1]    +    matX[2])*0.5 + dgaBCCPlus1 )  ;
  lbModel.fTemp8 [i][lbModel.G8_DV_P_P_M             ] = lbModel.fTemp[lbModel.DV_P_P_M         ]   * rho1 *(( matX[0]    +      matX[1]    -    matX[2])*0.5 + dgaBCCPlus1 )  ;
  lbModel.fTemp8 [i][lbModel.G8_DV_M_P_M             ] = lbModel.fTemp[lbModel.DV_M_P_M         ]   * rho1 *((-matX[0]    +      matX[1]    -    matX[2])*0.5 + dgaBCCPlus1 )  ;
  lbModel.fTemp8 [i][lbModel.G8_DV_M_M_M             ] = lbModel.fTemp[lbModel.DV_M_M_M         ]   * rho1 *((-matX[0]    -      matX[1]    -    matX[2])*0.5 + dgaBCCPlus1 )  ;
  lbModel.fTemp8 [i][lbModel.G8_DV_P_M_M             ] = lbModel.fTemp[lbModel.DV_P_M_M         ]   * rho1 *(( matX[0]    -      matX[1]    -    matX[2])*0.5 + dgaBCCPlus1 )  ;
  lbModel.fTemp9 [i][lbModel.G9_DV_P1_P1_P1          ] = lbModel.fTemp[lbModel.DV_P1_P1_P1      ]   * rho1 *(  matX[0]    +      matX[1]    +    matX[2]      + dgaBCC1Plus1)  ;
  lbModel.fTemp9 [i][lbModel.G9_DV_M1_P1_P1          ] = lbModel.fTemp[lbModel.DV_M1_P1_P1      ]   * rho1 *( -matX[0]    +      matX[1]    +    matX[2]      + dgaBCC1Plus1)  ;
  lbModel.fTemp9 [i][lbModel.G9_DV_M1_M1_P1          ] = lbModel.fTemp[lbModel.DV_M1_M1_P1      ]   * rho1 *( -matX[0]    -      matX[1]    +    matX[2]      + dgaBCC1Plus1)  ;
  lbModel.fTemp9 [i][lbModel.G9_DV_P1_M1_P1          ] = lbModel.fTemp[lbModel.DV_P1_M1_P1      ]   * rho1 *(  matX[0]    -      matX[1]    +    matX[2]      + dgaBCC1Plus1)  ;
  lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1         ] = lbModel.fTemp[lbModel.DV_P1_P1_M1      ]   * rho1 *(  matX[0]    +      matX[1]    -    matX[2]      + dgaBCC1Plus1)  ;
  lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1         ] = lbModel.fTemp[lbModel.DV_M1_P1_M1      ]   * rho1 *( -matX[0]    +      matX[1]    -    matX[2]      + dgaBCC1Plus1)  ;
  lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1         ] = lbModel.fTemp[lbModel.DV_M1_M1_M1      ]   * rho1 *( -matX[0]    -      matX[1]    -    matX[2]      + dgaBCC1Plus1)  ;
  lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1         ] = lbModel.fTemp[lbModel.DV_P1_M1_M1      ]   * rho1 *(  matX[0]    -      matX[1]    -    matX[2]      + dgaBCC1Plus1)  ; 
 }
};













template<typename dataType1>
void getGrad(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx)   
{
    dataType1 oneByTheta0             = 1.0/lbModel.theta0;
    dataType1 oneByTheta0Sqrd         = oneByTheta0*oneByTheta0;     
    dataType1 oneByTwoTheta0          = 0.5*oneByTheta0;   
    dataType1 oneByTwoTheta0Sqrd      = 0.50*oneByTheta0Sqrd;
    dataType1 oneByFourTheta0Sqrd     = 0.25*oneByTheta0Sqrd;   
    dataType1 twoByTheta0Sqrd         = 2.0*oneByTheta0Sqrd;
    dataType1 factorBCC               = (0.125*oneByTheta0Sqrd - oneByTwoTheta0);
    dataType1 factorBCC1              = (oneByTwoTheta0Sqrd - oneByTwoTheta0);
    dataType1 jX(0.0),jY(0.0),jZ(0.0),pTrace(0.0),oneByTwoTheta0TimesPTrace(0.0);
    dataType1 term0,term1,term2;
    
    for(int i=0;i<VECT_LENGTH;i++)
    {
        jX         =  uX[i]*rho[i]*oneByTheta0;
        jY         =  uY[i]*rho[i]*oneByTheta0;
        jZ         =  uZ[i]*rho[i]*oneByTheta0; 
        pTrace     =  Pxx[i] + Pyy[i] + Pzz[i];
        
        oneByTwoTheta0TimesPTrace =  oneByTwoTheta0*pTrace;
        
        lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO]   =   lbModel.w0    * (rho[i] - oneByTwoTheta0TimesPTrace);
        term0 = rho[i]- oneByTwoTheta0TimesPTrace;
        
        term1 = lbModel.wSC1  * (term0 + Pzz[i]*oneByTwoTheta0Sqrd);
        term2 = lbModel.wSC1  * jZ;
        lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P1]         =   term1 + term2 ;
        lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M1]         =   term1 - term2 ; 
        
        
        term1 = lbModel.wSC2  * (term0 + Pzz[i]*twoByTheta0Sqrd);
        term2 = lbModel.wSC2  * 2.0*jZ;
        lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P2]         =   term1 + term2 ;
        lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M2]         =   term1 - term2 ;
        
        
        
        term1 =  lbModel.wSC2  * (term0 + Pyy[i]*twoByTheta0Sqrd);
        term2 =  lbModel.wSC2  * 2.0*jY;
        lbModel.fTemp1[i][lbModel.G1_DV_ZERO_P2_ZERO]         =   term1 + term2;
        lbModel.fTemp2[i][lbModel.G2_DV_ZERO_M2_ZERO]         =   term1 - term2;
        
        
        term1 =  lbModel.wSC2  * (term0 + Pxx[i]*twoByTheta0Sqrd);
        term2 =  lbModel.wSC2  * 2.0*jX;
        lbModel.fTemp1[i][lbModel.G1_DV_P2_ZERO_ZERO]         =  term1 + term2 ;
        lbModel.fTemp2[i][lbModel.G2_DV_M2_ZERO_ZERO]         =  term1 - term2 ;
        
        
        term1 =  lbModel.wFCC  * (term0 + (Pyy[i]+Pzz[i])*oneByTwoTheta0Sqrd + Pyz[i]*oneByTheta0Sqrd);
        term2 =  lbModel.wFCC  * (jY + jZ);
        lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]           =   term1 + term2 ;
        lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1]           =   term1 - term2 ;
        
        
        term1 =  lbModel.wFCC  * (term0 + (Pyy[i]+Pzz[i])*oneByTwoTheta0Sqrd - Pyz[i]*oneByTheta0Sqrd);
        term2 =  lbModel.wFCC  * (jY - jZ);            
        lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1]           =   term1 + term2 ;
        lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1]           =   term1 - term2 ;
        
        
        term1 =  lbModel.wFCC  * (term0 + (Pxx[i]+Pzz[i])*oneByTwoTheta0Sqrd - Pzx[i]*oneByTheta0Sqrd);
        term2 =  lbModel.wFCC  * (jX - jZ);            
        lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1]           =   term1 + term2 ;
        lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1]           =   term1 - term2 ;
        
        term1 =  lbModel.wFCC  * (term0 + (Pxx[i]+Pzz[i])*oneByTwoTheta0Sqrd + Pzx[i]*oneByTheta0Sqrd);
        term2 =  lbModel.wFCC  * (jX + jZ);              
        lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1]           =   term1 + term2 ;
        lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1]           =   term1 - term2 ;
        
        term1 =  lbModel.wFCC  * (term0 + (Pxx[i]+Pyy[i])*oneByTwoTheta0Sqrd + Pxy[i]*oneByTheta0Sqrd);
        term2 =  lbModel.wFCC  * (jX + jY);             
        lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO  ]         =   term1 + term2 ;
        lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO  ]         =   term1 - term2 ;
        
        term1 =  lbModel.wFCC  * (term0 + (Pxx[i]+Pyy[i])*oneByTwoTheta0Sqrd - Pxy[i]*oneByTheta0Sqrd);
        term2 =  lbModel.wFCC  * (jX - jY);             
        lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO  ]         =   term1 + term2  ;
        lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO  ]         =   term1 - term2  ;
        
        term1 =  lbModel.wSC1  * (term0 + Pyy[i]*oneByTwoTheta0Sqrd);
        term2 =  lbModel.wSC1  * jY;              
        lbModel.fTemp5[i][lbModel.G5_DV_ZERO_P1_ZERO]         =   term1 + term2 ;
        lbModel.fTemp6[i][lbModel.G6_DV_ZERO_M1_ZERO]         =   term1 - term2 ;
        
        term1 =  lbModel.wSC1  * (term0 + Pxx[i]*oneByTwoTheta0Sqrd);
        term2 =  lbModel.wSC1  * jX;                
        lbModel.fTemp5[i][lbModel.G5_DV_P1_ZERO_ZERO]         =   term1 + term2 ;
        lbModel.fTemp6[i][lbModel.G6_DV_M1_ZERO_ZERO]         =   term1 - term2 ;
        
        term0 =  lbModel.wBCC  * (rho[i] + factorBCC*pTrace);
        term1 =  lbModel.wBCC  * ((Pxy[i] + Pyz[i] + Pzx[i] )*oneByFourTheta0Sqrd);
        term2 =  lbModel.wBCC  * 0.5 *(  jX + jY + jZ );                            
        lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]                =   term0 + term1 + term2;
        lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]                =   term0 + term1 - term2;
        
        term1 =  lbModel.wBCC  * (( -Pxy[i] + Pyz[i] - Pzx[i] )*oneByFourTheta0Sqrd);
        term2 =  lbModel.wBCC  * 0.5 *( -jX + jY + jZ );             
        lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]                =   term0 + term1 + term2;
        lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]                =   term0 + term1 - term2;
        
        
        term1 =  lbModel.wBCC  * ((Pxy[i] - Pyz[i] - Pzx[i])*oneByFourTheta0Sqrd);
        term2 =  lbModel.wBCC  * 0.5 *(-jX - jY + jZ);            
        lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]                =   term0 + term1 + term2;
        lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]                =   term0 + term1 - term2;
        
        term1 =  lbModel.wBCC  * (( -Pxy[i] - Pyz[i] + Pzx[i])*oneByFourTheta0Sqrd);
        term2 =  lbModel.wBCC  * 0.5 *(jX - jY + jZ);              
        lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]                =   term0 + term1 + term2;
        lbModel.fTemp8[i][lbModel.G8_DV_M_P_M]                =   term0 + term1 - term2;
        
        
        term0 =  lbModel.wBCC1  * (rho[i] + factorBCC1*pTrace);
        term1 =  lbModel.wBCC1  * ((Pxy[i] + Pyz[i] + Pzx[i] )*oneByTheta0Sqrd);
        term2 =  lbModel.wBCC1  * (  jX + jY + jZ );                            
        lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]             =   term0 + term1 + term2;
        lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]           =   term0 + term1 - term2;
        
        term1 =  lbModel.wBCC1  * (( -Pxy[i] + Pyz[i] - Pzx[i] )*oneByTheta0Sqrd);
        term2 =  lbModel.wBCC1  * ( -jX + jY + jZ );             
        lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]             =   term0 + term1 + term2;
        lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]           =   term0 + term1 - term2;
        
        
        term1 =  lbModel.wBCC1  * ((Pxy[i] - Pyz[i] - Pzx[i])*oneByTheta0Sqrd);
        term2 =  lbModel.wBCC1  * (-jX - jY + jZ);            
        lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]             =   term0 + term1 + term2;
        lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]           =   term0 + term1 - term2;
        
        term1 =  lbModel.wBCC1  * (( -Pxy[i] - Pyz[i] + Pzx[i])*oneByTheta0Sqrd);
        term2 =  lbModel.wBCC1  * (jX - jY + jZ);              
        lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]             =   term0 + term1 + term2;
        lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]           =   term0 + term1 - term2;            
    }
    
}    

    template<typename dataType1>
    void getFEqSinglePointIntoArray(lbmRD3Q41<dataType1> &lbModel , dataType1 rho, dataType1 uX, dataType1 uY, dataType1 uZ, dataType1 theta,dataType1*  array)
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
            array[lbModel.DV_ZERO_ZERO_ZERO]=  f0V*twoMinusU2;    
            
            // G1 and G2
            // ZZP1 and ZZM1
            order1       =  delThetaBy2     *lbModel.yMinus3_SC1;        
            order2       =  delSqTheta0p125 *lbModel.yTenMinusYSqrPlus15_SC1;         
            f0V          =  rho*lbModel.wt[lbModel.DV_ZERO_ZERO_P1]*(1.0 + order1+ order2);        
            uDotC        =  uZ;
            array[lbModel.DV_ZERO_ZERO_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC)); 
            array[lbModel.DV_ZERO_ZERO_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
            
            order1       =  delThetaBy2     *lbModel.yMinus3_SC2;        
            order2       =  delSqTheta0p125 *lbModel.yTenMinusYSqrPlus15_SC2;         
            f0V          =  rho*lbModel.wt[lbModel.DV_ZERO_ZERO_P2]*(1.0 + order1+ order2);   
            // ZZP2 and ZZM2
            uDotC        =  2.0*uZ;
            array[lbModel.DV_ZERO_ZERO_P2] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC)); 
            array[lbModel.DV_ZERO_ZERO_M2] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC)); 
            // ZP2Z and ZM2Z
            uDotC        =  2.0*uY;
            array[lbModel.DV_ZERO_P2_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC)); 
            array[lbModel.DV_ZERO_M2_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));         
            // P2ZZ and M2ZZ
            uDotC        =  2.0*uX;
            array[lbModel.DV_P2_ZERO_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC)); 
            array[lbModel.DV_M2_ZERO_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));   


      
            // G3 and G4
            order1    =  delThetaBy2    *lbModel.yMinus3_FCC;        
            order2    =  delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_FCC;         
            f0V       =  rho*lbModel.wt[lbModel.DV_ZERO_P1_P1]*(1.0 + order1+ order2);        
            //ZPP and ZMM         
            uDotC     =  (uY + uZ);        
            array[lbModel.DV_ZERO_P1_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));        
            array[lbModel.DV_ZERO_M1_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
            //ZMP and ZPM                                                                                            
            uDotC     =  (uZ - uY);                                                                       
            array[lbModel.DV_ZERO_M1_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));        
            array[lbModel.DV_ZERO_P1_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
            //PZP and MZM                                                                                            
            uDotC     =  (uX + uZ);                                                                       
            array[lbModel.DV_P1_ZERO_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));        
            array[lbModel.DV_M1_ZERO_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
            //MZP and PZM                                                                                            
            uDotC     =  (uZ - uX);                                                                   
            array[lbModel.DV_M1_ZERO_P1] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));    
            array[lbModel.DV_P1_ZERO_M1] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
            
      
            order1    =  delThetaBy2    *lbModel.yMinus3_FCC;
            order2    =  delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_FCC; 
            f0V       =  rho*lbModel.wt[lbModel.DV_P1_P1_ZERO]*(1.0 + order1+ order2);
            //PPZ and MMZ
            uDotC     =  uX + uY;
           array[lbModel.DV_P1_P1_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
           array[lbModel.DV_M1_M1_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
            //MPZ and PMZ                                                                               
            uDotC     =  (uY - uX);                                                             
            array[lbModel.DV_M1_P1_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            array[lbModel.DV_P1_M1_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
            //ZPZ and ZMZ
            order1    =  delThetaBy2    *lbModel.yMinus3_SC1;
            order2    =  delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_SC1; 
            f0V       =  rho*lbModel.wt[lbModel.DV_ZERO_P1_ZERO]*(1.0 + order1+ order2);    
            uDotC     =  uY;
            array[lbModel.DV_ZERO_P1_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            array[lbModel.DV_ZERO_M1_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));
            //PZZ aqnd MZZ                                                                                  
            uDotC     =   uX;                                                                          
            array[lbModel.DV_P1_ZERO_ZERO] = f0V*(twoMinusU2 + uDotC * (2.0 + uDotC));
            array[lbModel.DV_M1_ZERO_ZERO] = f0V*(twoMinusU2 - uDotC * (2.0 - uDotC));        
  
            
         
            order1    =  1.0 + delThetaBy2*lbModel.yMinus3_BCC +  delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_BCC; 
            f0V       =  rho*lbModel.wt[lbModel.DV_P_P_P]*(order1);
            order2    =  1.0 + delThetaBy2*lbModel.yMinus3_BCC1 + delSqTheta0p125*lbModel.yTenMinusYSqrPlus15_BCC1;         
            f0V1      =  rho*lbModel.wt[lbModel.DV_P1_P1_P1]*(order2);        
            //PPP1 and MMM1        
            uDotC     = (uX + uY + uZ);      
            uDotCSq   = uDotC*uDotC;        
            array[lbModel.DV_P1_P1_P1]   = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);
            array[lbModel.DV_M1_M1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);   
            //PPP and MMM
            //uDotC     =   uDotC*0.5 ;
            array[lbModel.DV_P_P_P]   = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);
            array[lbModel.DV_M_M_M]   = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);
            
            //MPP1 and PMM1                                                                             
            uDotC     = (uY - uX + uZ); 
            uDotCSq   = uDotC*uDotC;
            array[lbModel.DV_M1_P1_P1]  = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);
            array[lbModel.DV_P1_M1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);  
            //MPP and PMM                                                                           
            // uDotC     =  uDotC*0.5 ;                                           
            array[lbModel.DV_M_P_P] = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);     
            array[lbModel.DV_P_M_M] = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);       
            
            //MMP1 and PPM 1                                                                                     
            uDotC     = (uZ - uX - uY);           
            uDotCSq   = uDotC*uDotC;       
            array[lbModel.DV_M1_M1_P1] = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);       
            array[lbModel.DV_P1_P1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);        
            //MMP and PPM                                                                                  
            // uDotC     =   uDotC*0.5 ;                                                  
            array[lbModel.DV_M_M_P] = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);       
            array[lbModel.DV_P_P_M] = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);       
            
            //PMP and MPM                                                                                      
            uDotC     = (uX - uY + uZ);        
            uDotCSq   = uDotC*uDotC;       
            array[lbModel.DV_P1_M1_P1]  = f0V1*(twoMinusU2 + uDotC * 2.0 + uDotCSq);       
            array[lbModel.DV_M1_P1_M1] = f0V1*(twoMinusU2 - uDotC * 2.0 + uDotCSq);         
            //PMP and MPM                                                                                  
            // uDotC     =  uDotC *0.5 ;                                                  
            array[lbModel.DV_P_M_P] = f0V*(twoMinusU2 + uDotC  + 0.25*uDotCSq);       
            array[lbModel.DV_M_P_M] = f0V*(twoMinusU2 - uDotC  + 0.25*uDotCSq);  
            
        
    }
  
        
    template<typename dataType1>
    void getGradSinglePointIntoArray(lbmRD3Q41<dataType1> &lbModel , dataType1 rho, dataType1 uX, dataType1 uY, dataType1 uZ, dataType1 Pxx, dataType1 Pyy, dataType1 Pzz, dataType1 Pxy, dataType1 Pyz, dataType1 Pzx ,dataType1*  array)
    {
    dataType1 oneByTheta0             = 1.0/lbModel.theta0;
    dataType1 oneByTheta0Sqrd         = oneByTheta0*oneByTheta0;     
    dataType1 oneByTwoTheta0          = 0.5*oneByTheta0;   
    dataType1 oneByTwoTheta0Sqrd      = 0.50*oneByTheta0Sqrd;
    dataType1 oneByFourTheta0Sqrd     = 0.25*oneByTheta0Sqrd;   
    dataType1 twoByTheta0Sqrd         = 2.0*oneByTheta0Sqrd;
    dataType1 factorBCC               = (0.125*oneByTheta0Sqrd - oneByTwoTheta0);
    dataType1 factorBCC1              = (oneByTwoTheta0Sqrd - oneByTwoTheta0);
    dataType1 jX(0.0),jY(0.0),jZ(0.0),pTrace(0.0),oneByTwoTheta0TimesPTrace(0.0);
    dataType1 term0,term1,term2;
    
        jX         =  uX*rho*oneByTheta0;
        jY         =  uY*rho*oneByTheta0;
        jZ         =  uZ*rho*oneByTheta0; 
        pTrace     =  Pxx + Pyy + Pzz;
        
        oneByTwoTheta0TimesPTrace =  oneByTwoTheta0*pTrace;
        
        array[lbModel.DV_ZERO_ZERO_ZERO]   =   lbModel.w0    * (rho - oneByTwoTheta0TimesPTrace);
        term0 = rho- oneByTwoTheta0TimesPTrace;
        
        term1 = lbModel.wSC1  * (term0 + Pzz*oneByTwoTheta0Sqrd);
        term2 = lbModel.wSC1  * jZ;
        array[lbModel.DV_ZERO_ZERO_P1]         =   term1 + term2 ;
        array[lbModel.DV_ZERO_ZERO_M1]         =   term1 - term2 ; 
        
        
        term1 = lbModel.wSC2  * (term0 + Pzz*twoByTheta0Sqrd);
        term2 = lbModel.wSC2  * 2.0*jZ;
        array[lbModel.DV_ZERO_ZERO_P2]         =   term1 + term2 ;
        array[lbModel.DV_ZERO_ZERO_M2]         =   term1 - term2 ;
        
        
        
        term1 =  lbModel.wSC2  * (term0 + Pyy*twoByTheta0Sqrd);
        term2 =  lbModel.wSC2  * 2.0*jY;
        array[lbModel.DV_ZERO_P2_ZERO]         =   term1 + term2;
        array[lbModel.DV_ZERO_M2_ZERO]         =   term1 - term2;
        
        
        term1 =  lbModel.wSC2  * (term0 + Pxx*twoByTheta0Sqrd);
        term2 =  lbModel.wSC2  * 2.0*jX;
        array[lbModel.DV_P2_ZERO_ZERO]         =  term1 + term2 ;
        array[lbModel.DV_M2_ZERO_ZERO]         =  term1 - term2 ;
        
        
        term1 =  lbModel.wFCC  * (term0 + (Pyy+Pzz)*oneByTwoTheta0Sqrd + Pyz*oneByTheta0Sqrd);
        term2 =  lbModel.wFCC  * (jY + jZ);
        array[lbModel.DV_ZERO_P1_P1]           =   term1 + term2 ;
        array[lbModel.DV_ZERO_M1_M1]           =   term1 - term2 ;
        
        
        term1 =  lbModel.wFCC  * (term0 + (Pyy+Pzz)*oneByTwoTheta0Sqrd - Pyz*oneByTheta0Sqrd);
        term2 =  lbModel.wFCC  * (jY - jZ);            
        array[lbModel.DV_ZERO_P1_M1]           =   term1 + term2 ;
        array[lbModel.DV_ZERO_M1_P1]           =   term1 - term2 ;
        
        
        term1 =  lbModel.wFCC  * (term0 + (Pxx+Pzz)*oneByTwoTheta0Sqrd - Pzx*oneByTheta0Sqrd);
        term2 =  lbModel.wFCC  * (jX - jZ);            
        array[lbModel.DV_P1_ZERO_M1]           =   term1 + term2 ;
        array[lbModel.DV_M1_ZERO_P1]           =   term1 - term2 ;
        
        term1 =  lbModel.wFCC  * (term0 + (Pxx+Pzz)*oneByTwoTheta0Sqrd + Pzx*oneByTheta0Sqrd);
        term2 =  lbModel.wFCC  * (jX + jZ);              
        array[lbModel.DV_P1_ZERO_P1]           =   term1 + term2 ;
        array[lbModel.DV_M1_ZERO_M1]           =   term1 - term2 ;
        
        term1 =  lbModel.wFCC  * (term0 + (Pxx+Pyy)*oneByTwoTheta0Sqrd + Pxy*oneByTheta0Sqrd);
        term2 =  lbModel.wFCC  * (jX + jY);             
        array[lbModel.DV_P1_P1_ZERO  ]         =   term1 + term2 ;
        array[lbModel.DV_M1_M1_ZERO  ]         =   term1 - term2 ;
        
        term1 =  lbModel.wFCC  * (term0 + (Pxx+Pyy)*oneByTwoTheta0Sqrd - Pxy*oneByTheta0Sqrd);
        term2 =  lbModel.wFCC  * (jX - jY);             
        array[lbModel.DV_P1_M1_ZERO  ]         =   term1 + term2  ;
        array[lbModel.DV_M1_P1_ZERO  ]         =   term1 - term2  ;
        
        term1 =  lbModel.wSC1  * (term0 + Pyy*oneByTwoTheta0Sqrd);
        term2 =  lbModel.wSC1  * jY;              
        array[lbModel.DV_ZERO_P1_ZERO]         =   term1 + term2 ;
        array[lbModel.DV_ZERO_M1_ZERO]         =   term1 - term2 ;
        
        term1 =  lbModel.wSC1  * (term0 + Pxx*oneByTwoTheta0Sqrd);
        term2 =  lbModel.wSC1  * jX;                
        array[lbModel.DV_P1_ZERO_ZERO]         =   term1 + term2 ;
        array[lbModel.DV_M1_ZERO_ZERO]         =   term1 - term2 ;
        
        term0 =  lbModel.wBCC  * (rho + factorBCC*pTrace);
        term1 =  lbModel.wBCC  * ((Pxy + Pyz + Pzx )*oneByFourTheta0Sqrd);
        term2 =  lbModel.wBCC  * 0.5 *(  jX + jY + jZ );                            
        array[lbModel.DV_P_P_P]                =   term0 + term1 + term2;
        array[lbModel.DV_M_M_M]                =   term0 + term1 - term2;
        
        term1 =  lbModel.wBCC  * (( -Pxy + Pyz - Pzx )*oneByFourTheta0Sqrd);
        term2 =  lbModel.wBCC  * 0.5 *( -jX + jY + jZ );             
        array[lbModel.DV_M_P_P]                =   term0 + term1 + term2;
        array[lbModel.DV_P_M_M]                =   term0 + term1 - term2;
        
        
        term1 =  lbModel.wBCC  * ((Pxy - Pyz - Pzx)*oneByFourTheta0Sqrd);
        term2 =  lbModel.wBCC  * 0.5 *(-jX - jY + jZ);            
        array[lbModel.DV_M_M_P]                =   term0 + term1 + term2;
        array[lbModel.DV_P_P_M]                =   term0 + term1 - term2;
        
        term1 =  lbModel.wBCC  * (( -Pxy - Pyz + Pzx)*oneByFourTheta0Sqrd);
        term2 =  lbModel.wBCC  * 0.5 *(jX - jY + jZ);              
        array[lbModel.DV_P_M_P]                =   term0 + term1 + term2;
        array[lbModel.DV_M_P_M]                =   term0 + term1 - term2;
        
        
        term0 =  lbModel.wBCC1  * (rho + factorBCC1*pTrace);
        term1 =  lbModel.wBCC1  * ((Pxy + Pyz + Pzx )*oneByTheta0Sqrd);
        term2 =  lbModel.wBCC1  * (  jX + jY + jZ );                            
        array[lbModel.DV_P1_P1_P1]             =   term0 + term1 + term2;
        array[lbModel.DV_M1_M1_M1]           =   term0 + term1 - term2;
        
        term1 =  lbModel.wBCC1  * (( -Pxy + Pyz - Pzx )*oneByTheta0Sqrd);
        term2 =  lbModel.wBCC1  * ( -jX + jY + jZ );             
        array[lbModel.DV_M1_P1_P1]             =   term0 + term1 + term2;
        array[lbModel.DV_P1_M1_M1]           =   term0 + term1 - term2;
        
        
        term1 =  lbModel.wBCC1  * ((Pxy - Pyz - Pzx)*oneByTheta0Sqrd);
        term2 =  lbModel.wBCC1  * (-jX - jY + jZ);            
        array[lbModel.DV_M1_M1_P1]             =   term0 + term1 + term2;
        array[lbModel.DV_P1_P1_M1]           =   term0 + term1 - term2;
        
        term1 =  lbModel.wBCC1  * (( -Pxy - Pyz + Pzx)*oneByTheta0Sqrd);
        term2 =  lbModel.wBCC1  * (jX - jY + jZ);              
        array[lbModel.DV_P1_M1_P1]             =   term0 + term1 + term2;
        array[lbModel.DV_M1_P1_M1]           =   term0 + term1 - term2;            
    }
        
        
        



        
// Explicit definition
template void getFEq<double>(lbmRD3Q41<double> &,int , double *, double *, double *, double *, double *);
template void getFEqSIMD<double>(lbmRD3Q41<double> &,int , double *, double *, double *, double *, double *, int);
template void getIterativeFEq<double>(lbmRD3Q41<double> &,int , double *, double *, double *, double *, double *);
template void getFEqNewIterative<double>(lbmRD3Q41<double> &,int,double *,double *,double *,double *,double *);
template void getFEqSinglePointIntoArray<double>(lbmRD3Q41<double> &, double, double , double,double, double,double*  );
template void getGrad<double>(lbmRD3Q41<double> &lbModel,int VECT_LENGTH, double *rho, double *uX, double *uY, double *uZ, double *Pxx, double *Pyy, double *Pzz, double *Pxy, double *Pyz, double *Pzx); 
template void getGradSinglePointIntoArray<double>(lbmRD3Q41<double> &lbModel , double rho, double uX, double uY, double uZ, double Pxx, double Pyy, double Pzz, double Pxy, double Pyz, double Pzx ,double*  array);
