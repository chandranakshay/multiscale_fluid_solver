#include"advection41.h"

template <int N,int numblock, typename dataType1>
void advectionG1(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)            // ZZP,ZZP2,ZP2Z,P2ZZ
{     
 unsigned long long int indexMul4i1;   
 for(int nodeType =0; nodeType<2;nodeType++)
 {
  for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
  {                                                                
   for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                  
   {                    
    unsigned long long int index = myGrid.getIndex(0,i2,i3,lbModel.G1,nodeType, 0);
    for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
    {
     indexMul4i1 = index+i1*4;    
     myGrid.value(indexMul4i1)                              = myGrid(i1  ,i2  ,i3-1,lbModel.G1,nodeType, lbModel.G1_DV_ZERO_ZERO_P1);
     myGrid.value(indexMul4i1 + lbModel.G1_DV_ZERO_ZERO_P2) = myGrid(i1  ,i2  ,i3-2,lbModel.G1,nodeType, lbModel.G1_DV_ZERO_ZERO_P2); 
     myGrid.value(indexMul4i1 + lbModel.G1_DV_ZERO_P2_ZERO) = myGrid(i1  ,i2-2,i3  ,lbModel.G1,nodeType, lbModel.G1_DV_ZERO_P2_ZERO);
     myGrid.value(indexMul4i1 + lbModel.G1_DV_P2_ZERO_ZERO) = myGrid.value(index+(i1-2)*4 + lbModel.G1_DV_P2_ZERO_ZERO);    
    }
   }
  }
 }
}
 
template <int N,int numblock, typename dataType1>
void advectionG2(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)            // ZZM,ZZM2,ZM2Z,M2ZZ
{         
 unsigned long long int  indexMul4i1;   
 for(int nodeType =0; nodeType<2;nodeType++)
 {
  for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)
  {
   for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)
   {
    unsigned long long int index = myGrid.getIndex(0,i2,i3,lbModel.G2,nodeType, 0);
    for (int i1 = myGrid.ndB1; i1 <= myGrid.ndE1; i1++)
    {    
     indexMul4i1 = index+i1*4;    
     myGrid.value(indexMul4i1)                              = myGrid(i1  ,i2  ,i3+1,lbModel.G2,nodeType, lbModel.G2_DV_ZERO_ZERO_M1);
     myGrid.value(indexMul4i1 + lbModel.G2_DV_ZERO_ZERO_M2) = myGrid(i1  ,i2  ,i3+2,lbModel.G2,nodeType, lbModel.G2_DV_ZERO_ZERO_M2); 
     myGrid.value(indexMul4i1 + lbModel.G2_DV_ZERO_M2_ZERO) = myGrid(i1  ,i2+2,i3  ,lbModel.G2,nodeType, lbModel.G2_DV_ZERO_M2_ZERO);
     myGrid.value(indexMul4i1 + lbModel.G2_DV_M2_ZERO_ZERO) = myGrid.value(index+(i1+2)*4 + lbModel.G2_DV_M2_ZERO_ZERO);     
    }
   }
  }
 }
}
 
template <int N,int numblock, typename dataType1>                                                
void advectionG3(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)           // ZPP,ZMP,PZP,MZP 
{
 unsigned long long int indexMul4i1;   
 for(int nodeType =0; nodeType<2;nodeType++)
 {
  for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)                   
  {                                                                
   for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)                  
   {                                                               
    unsigned long long int  index = myGrid.getIndex(0,i2,i3,lbModel.G3,nodeType, 0);
    for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)                 
    {  
     indexMul4i1 = index+i1*4;    
     myGrid.value(indexMul4i1)                            = myGrid(i1  ,i2-1,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1);
     myGrid.value(indexMul4i1 + lbModel.G3_DV_ZERO_M1_P1) = myGrid(i1  ,i2+1,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1);
     myGrid.value(indexMul4i1 + lbModel.G3_DV_P1_ZERO_P1) = myGrid(i1-1,i2  ,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1);     
     myGrid.value(indexMul4i1 + lbModel.G3_DV_M1_ZERO_P1) = myGrid(i1+1,i2  ,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1);  
    }                                                            
   }
  }
 }
}

template <int N,int numblock, typename dataType1>
void advectionG4(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)           // ZPM,ZMM,PZM,MZM
{
 unsigned long long int  indexMul4i1;   
 for(int nodeType =0; nodeType<2;nodeType++)
 {
  for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)
  {
   for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)
   {
    unsigned long long int  index = myGrid.getIndex(0,i2,i3,lbModel.G4,nodeType, 0);
    for (int i1 = myGrid.ndB1; i1 <= myGrid.ndE1; i1++)
    {
     indexMul4i1 = index+i1*4;    
     myGrid.value(indexMul4i1 + lbModel.G4_DV_M1_ZERO_M1) = myGrid(i1+1,i2  ,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1);     
     myGrid.value(indexMul4i1 + lbModel.G4_DV_ZERO_P1_M1) = myGrid(i1  ,i2-1,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1);
     myGrid.value(indexMul4i1)                            = myGrid(i1  ,i2+1,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1);
     myGrid.value(indexMul4i1 + lbModel.G4_DV_P1_ZERO_M1) = myGrid(i1-1,i2  ,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1); 
    }                                              
   } 
  }
 }
}  

template <int N,int numblock, typename dataType1>
void advectionG5(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)          // ZPM,ZMM,PZM,MZM
{
 unsigned long long int indexMul4i1;   
 for(int nodeType =0; nodeType<2;nodeType++)
 {
  for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)  
  {
   for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)      
   {       
    unsigned long long int index = myGrid.getIndex(0,i2,i3,lbModel.G5,nodeType, 0);
    for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)     
    {                                                     
     indexMul4i1 = index+i1*4;    
     myGrid.value(indexMul4i1)                              = myGrid(i1-1 ,i2-1,i3,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO  );
     myGrid.value(indexMul4i1 + lbModel.G5_DV_M1_P1_ZERO  ) = myGrid(i1+1 ,i2-1,i3,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO  );
     myGrid.value(indexMul4i1 + lbModel.G5_DV_ZERO_P1_ZERO) = myGrid(i1   ,i2-1,i3,lbModel.G5,nodeType,lbModel.G5_DV_ZERO_P1_ZERO);     
     myGrid.value(indexMul4i1 + lbModel.G5_DV_P1_ZERO_ZERO) = myGrid(i1-1 ,i2  ,i3,lbModel.G5,nodeType,lbModel.G5_DV_P1_ZERO_ZERO); 
    }                                              
   }  
  }
 }
}
 
template <int N,int numblock, typename dataType1>
void advectionG6(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)            // ZZP,ZZT,ZTZ,TZZ
{
 unsigned long long int indexMul4i1;   
 for(int nodeType =0; nodeType<2;nodeType++)
 { 
  for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)  
  {
   for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)
   {
    unsigned long long int index = myGrid.getIndex(0,i2,i3,lbModel.G6,nodeType, 0);
    for (int i1 = myGrid.ndB1; i1 <= myGrid.ndE1; i1++)
    {
     indexMul4i1 = index+i1*4;    
     myGrid.value(indexMul4i1)                              = myGrid(i1+1 ,i2+1 ,i3,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO  );
     myGrid.value(indexMul4i1 + lbModel.G6_DV_P1_M1_ZERO  ) = myGrid(i1-1 ,i2+1 ,i3,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO  );
     myGrid.value(indexMul4i1 + lbModel.G6_DV_ZERO_M1_ZERO) = myGrid(i1   ,i2+1 ,i3,lbModel.G6,nodeType,lbModel.G6_DV_ZERO_M1_ZERO); 
     myGrid.value(indexMul4i1 + lbModel.G6_DV_M1_ZERO_ZERO) = myGrid.value(index+(i1+1)*4+lbModel.G6_DV_M1_ZERO_ZERO);     
    }
   }
  }
 }
}
  
template <int N,int numblock, typename dataType1>
void advectionG7(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)          // ZPM,ZMM,PZM,MZM
{
 unsigned long long int indexMul4i1;   
 for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)
 {
  for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)
  {
   unsigned long long int index = myGrid.getIndex(0,i2,i3,lbModel.G7,myGrid.node[lbModel.G7], 0);
   for (int i1 = myGrid.ndB1; i1 <= myGrid.ndE1; i1++)
   {    
     indexMul4i1 = index+i1*4;    
    myGrid.value(indexMul4i1 + lbModel.G7_DV_M_P_P) = myGrid(i1+1,i2  ,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);
    myGrid.value(indexMul4i1 + lbModel.G7_DV_M_M_P) = myGrid(i1+1,i2+1,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);
    myGrid.value(indexMul4i1 + lbModel.G7_DV_P_M_P) = myGrid(i1  ,i2+1,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);   
   }
  }
 }

 for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)  
 {
  for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)      
  {                                                      
   unsigned long long int  index = myGrid.getIndex(0,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7], 0);
   for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)     
   { 
    indexMul4i1 = index+i1*4;    
    myGrid.value(indexMul4i1)                       = myGrid(i1-1,i2-1,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);
    myGrid.value(indexMul4i1 + lbModel.G7_DV_M_P_P) = myGrid(i1  ,i2-1,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);
    myGrid.value(indexMul4i1 + lbModel.G7_DV_M_M_P) = myGrid(i1  ,i2  ,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);
    myGrid.value(indexMul4i1 + lbModel.G7_DV_P_M_P) = myGrid(i1-1,i2  ,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);       
   }
  }
 }
  int               temp  = myGrid.node[lbModel.G7];
  myGrid.node[lbModel.G7] = myGrid.cell[lbModel.G7];
  myGrid.cell[lbModel.G7] = temp;
}
 
template <int N,int numblock, typename dataType1>
void advectionG8(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)           // ZPM,ZMM,PZM,MZM
{
 unsigned long long int indexMul4i1;   
 for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)  
 {
  for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)      
  {                                                      
   unsigned long long int index = myGrid.getIndex(0,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8], 0);
   for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)     
   {    
    indexMul4i1 = index+i1*4;    
    myGrid.value(indexMul4i1 + lbModel.G8_DV_M_P_M) = myGrid(i1  ,i2-1,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);
    myGrid.value(indexMul4i1 + lbModel.G8_DV_P_M_M) = myGrid(i1-1,i2  ,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);
    myGrid.value(indexMul4i1 + lbModel.G8_DV_P_P_M) = myGrid(i1-1,i2-1,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);
   }
  }
 }
  
 for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)
 {
  for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)
  {
   unsigned long long int index = myGrid.getIndex(0,i2,i3,lbModel.G8,myGrid.node[lbModel.G8], 0);
   for (int i1 = myGrid.ndB1; i1 <= myGrid.ndE1; i1++)
   {   
    indexMul4i1 = index+i1*4;    
    myGrid.value(indexMul4i1)                       = myGrid(i1+1,i2+1,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);
    myGrid.value(indexMul4i1 + lbModel.G8_DV_M_P_M) = myGrid(i1+1,i2  ,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);
    myGrid.value(indexMul4i1 + lbModel.G8_DV_P_M_M) = myGrid(i1  ,i2+1,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);
    myGrid.value(indexMul4i1 + lbModel.G8_DV_P_P_M) = myGrid(i1  ,i2  ,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);
   }
  }
 }
  int               temp  = myGrid.node[lbModel.G8];
  myGrid.node[lbModel.G8] = myGrid.cell[lbModel.G8];
  myGrid.cell[lbModel.G8] = temp;
}
  
template <int N,int numblock, typename dataType1>
void advectionG9(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)           //PPP,MPP,MMP,PMP
{
 unsigned long long int indexMul4i1;   
 for(int nodeType =0; nodeType<2;nodeType++)
 {
  for (int i3 = myGrid.ndE3; i3 >= myGrid.ndB3; i3--)  
  {
   for (int i2 = myGrid.ndE2; i2 >= myGrid.ndB2; i2--)      
   {        
    unsigned long long int index  = myGrid.getIndex(0,i2  ,i3,  lbModel.G9,nodeType, 0);
    unsigned long long int index1 = myGrid.getIndex(0,i2+1,i3-1,lbModel.G9,nodeType, 0);
    unsigned long long int index2 = myGrid.getIndex(0,i2-1,i3-1,lbModel.G9,nodeType, 0);
    for (int i1 = myGrid.ndE1; i1 >= myGrid.ndB1; i1--)     
    {
     indexMul4i1 = index+i1*4;    
     myGrid.value(indexMul4i1)                          = myGrid.value(index2+ (i1-1)*4);
     myGrid.value(indexMul4i1 + lbModel.G9_DV_M1_P1_P1) = myGrid.value(index2+ (i1+1)*4 +lbModel.G9_DV_M1_P1_P1);
     myGrid.value(indexMul4i1 + lbModel.G9_DV_M1_M1_P1) = myGrid.value(index1+ (i1+1)*4 +lbModel.G9_DV_M1_M1_P1);
     myGrid.value(indexMul4i1 + lbModel.G9_DV_P1_M1_P1) = myGrid.value(index1+ (i1-1)*4 +lbModel.G9_DV_P1_M1_P1);  
    }
   }
  } 
 }
}
  
template <int N,int numblock, typename dataType1>
void advectionG10(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)           //PPP,MPP,MMP,PMP
{
 unsigned long long int indexMul4i1;   
 for(int nodeType =0; nodeType<2;nodeType++)
 {
  for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)
  {
   for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)
   {
    unsigned long long int index  = myGrid.getIndex(0,i2  ,i3  ,lbModel.G10,nodeType, 0);
    unsigned long long int index1 = myGrid.getIndex(0,i2+1,i3+1,lbModel.G10,nodeType, 0);
    unsigned long long int index2 = myGrid.getIndex(0,i2-1,i3+1,lbModel.G10,nodeType, 0);
    for (int i1 = myGrid.ndB1; i1 <= myGrid.ndE1; i1++)
    {   
     indexMul4i1 = index+i1*4;    
     myGrid.value(indexMul4i1)                           = myGrid.value(index1 + (i1+1)*4);
     myGrid.value(indexMul4i1 + lbModel.G10_DV_M1_P1_M1) = myGrid.value(index2 + (i1+1)*4 +lbModel.G10_DV_M1_P1_M1);
     myGrid.value(indexMul4i1 + lbModel.G10_DV_P1_M1_M1) = myGrid.value(index1 + (i1-1)*4 +lbModel.G10_DV_P1_M1_M1);
     myGrid.value(indexMul4i1 + lbModel.G10_DV_P1_P1_M1) = myGrid.value(index2 + (i1-1)*4 +lbModel.G10_DV_P1_P1_M1);  
    }
   }
  }
 }
}
 
template <int N,int numblock, typename dataType1>
void advection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)
{
 advectionG1 (lbModel,myGrid); 
 advectionG2 (lbModel,myGrid);
 advectionG3 (lbModel,myGrid); 
 advectionG4 (lbModel,myGrid);
 advectionG5 (lbModel,myGrid); 
 advectionG6 (lbModel,myGrid);   
 advectionG7 (lbModel,myGrid); 
 advectionG8 (lbModel,myGrid);
 advectionG9 (lbModel,myGrid); 
 advectionG10(lbModel,myGrid); 
}
 

 
 //Explicit function declarations
 template void advectionG1<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &) ;
 template void advectionG2<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &) ;
 template void advectionG3<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &) ;
 template void advectionG4<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &) ;
 template void advectionG5<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &) ;
 template void advectionG6<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &) ;
 template void advectionG7<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &) ;
 template void advectionG8<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &) ;
 template void advectionG9<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &) ;
 template void advectionG10<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &) ;
 template void advection<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &) ;

 
