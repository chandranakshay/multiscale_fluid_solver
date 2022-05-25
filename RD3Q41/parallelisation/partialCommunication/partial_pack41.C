/*********************************************************************************************
 * Copyright (c) <2015>, <Santosh Ansumali@JNCASR>                                           *
 *  All rights reserved.                                                                     *
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
#include"partial_pack41.h"


PANINI_REAL* partial_allocateForPacking(int m1, int m2)
{
 PANINI_REAL *temp;  
 temp = new PANINI_REAL[m1*m2*30];
 return temp;
} 

template <int N,int numblock, typename dataType1>
void partial_packTo1P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n2*myGrid.n3;  
 int stride = 1; 

 probeComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;
 probeComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 probeComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 probeComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
 probeComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1,lbModel.G7_DV_P_M_P);  
 index += size;       
 probeComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
////////////////// 
 probeComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 probeComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
////////////////// 
 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;
 probeComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 probeComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 probeComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
 probeComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1,lbModel.G7_DV_P_M_P);  
 index += size;       
 probeComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
////////////////// 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 probeComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
}

template <int N,int numblock, typename dataType1>
void partial_unpackFrom1P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n2*myGrid.n3;  
 int stride = 1; 
 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-2,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;
 modifyComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 modifyComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 modifyComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
 modifyComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G7_DV_P_M_P);  
 index += size;       
 modifyComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
////////////////// 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 modifyComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
//////////////////  
 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-2,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G7_DV_P_M_P);  
 index += size;       
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
////////////////// 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
//////////////////  
}

template <int N,int numblock, typename dataType1>
void partial_packTo1M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n2*myGrid.n3;  
 int stride = 1; 
  
 probeComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1+1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;
 probeComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 probeComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 probeComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
 probeComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1,lbModel.G8_DV_M_P_M);  
 index += size;       
 probeComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
///////////////////////
 probeComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
 probeComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
/////////////////////// 
 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1+1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;
 probeComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 probeComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 probeComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
 probeComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1,lbModel.G8_DV_M_P_M);  
 index += size;       
 probeComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
///////////////////////
 probeComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
 probeComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
/////////////////////// 
}

template <int N,int numblock, typename dataType1>
void partial_unpackFrom1M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n2*myGrid.n3;  
 int stride = 1; 

 modifyComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+2,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;
 modifyComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 modifyComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 modifyComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
 modifyComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G8_DV_M_P_M);  
 index += size;       
 modifyComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
///////////////////////
 modifyComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
/////////////////////// 
 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+2,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G8_DV_M_P_M);  
 index += size;       
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
///////////////////////
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
/////////////////////// 
}

template <int N,int numblock, typename dataType1>
void partial_packMSgTo1Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp1Psend, dataType1 *temp1Msend)
{
 int index = 0;
 partial_packTo1P(lbModel,myGrid,temp1Psend,index);  
 index = 0;
 partial_packTo1M(lbModel,myGrid,temp1Msend,index); 
}

template <int N,int numblock, typename dataType1>
void partial_recvMSgFrom1Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp1Preceive, dataType1 *temp1Mreceive)
{
 int index = 0;
 partial_unpackFrom1P(lbModel,myGrid,temp1Preceive,index);  
 index = 0;
 partial_unpackFrom1M(lbModel,myGrid,temp1Mreceive,index); 
}

template <int N,int numblock, typename dataType1>
void partial_packTo2P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n3;  
 int stride = 1; 
 
 probeComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;
 probeComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2,lbModel.G8_DV_M_P_M);  
 index += size;       
 probeComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2,lbModel.G8_DV_P_P_M);  
 index += size;     
///////////////////////
 probeComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
 probeComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
////////////////////////
 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;
 probeComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2,lbModel.G8_DV_M_P_M);  
 index += size;       
 probeComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2,lbModel.G8_DV_P_P_M);  
 index += size;     
///////////////////////
 probeComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
 probeComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
////////////////////////
}

template <int N,int numblock, typename dataType1>
void partial_unpackFrom2P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n3;  
 int stride = 1; 
 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-2,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;
 modifyComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G8_DV_M_P_M);  
 index += size;       
 modifyComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G8_DV_P_P_M);                                            
 index += size;                                                                                                                                            
///////////////////
 modifyComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
/////////////////// 
 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-2,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G8_DV_M_P_M);  
 index += size;       
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G8_DV_P_P_M);                                            
 index += size;                                                                                                                                            
///////////////////
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
/////////////////// 
}


template <int N,int numblock, typename dataType1>
void partial_packTo2M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n3;  
 int stride = 1; 
 
 probeComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2+1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;
 probeComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2,lbModel.G7_DV_P_M_P);  
 index += size;       
 probeComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
//////////////////////
 probeComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 probeComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
///////////////////// 
 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2+1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;
 probeComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2,lbModel.G7_DV_P_M_P);  
 index += size;       
 probeComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
//////////////////////
 probeComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 probeComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
///////////////////// 
}

template <int N,int numblock, typename dataType1>
void partial_unpackFrom2M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n3;  
 int stride = 1; 
 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+2,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;
 modifyComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G7_DV_P_M_P);  
 index += size;       
 modifyComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
////////////////////
 modifyComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 modifyComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
/////////////////// 

 modifyComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+2,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G7_DV_P_M_P);  
 index += size;       
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
////////////////////
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
/////////////////// 
}


template <int N,int numblock, typename dataType1>
void partial_packMSgTo2Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp2Psend, dataType1 *temp2Msend)
{
 int index = 0;
 partial_packTo2P(lbModel,myGrid,temp2Psend,index);  
 index = 0;
 partial_packTo2M(lbModel,myGrid,temp2Msend,index); 
}

template <int N,int numblock, typename dataType1>
void partial_recvMSgFrom2Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp2Preceive, dataType1 *temp2Mreceive)
{
 int index = 0;
 partial_unpackFrom2P(lbModel,myGrid,temp2Preceive,index);  
 index = 0;
 partial_unpackFrom2M(lbModel,myGrid,temp2Mreceive,index); 
}


template <int N,int numblock, typename dataType1>
void partial_packTo3P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n2;  
 int stride = 1; 
 
 probeComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;
 probeComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 probeComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_M_P_P);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_P_M_P);  
 index += size;       
//////////////////////////
 probeComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
////////////////////////////// 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;
 probeComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 probeComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_M_P_P);  
 index += size;     
 probeComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_P_M_P);  
 index += size;       
//////////////////////////
 probeComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
}

template <int N,int numblock, typename dataType1>
void partial_unpackFrom3P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n2;  
 int stride = 1; 
 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-2,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;
 modifyComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 modifyComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_P_M_P);  
 index += size;       
 ///////////////
 modifyComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
/////////////////// 
 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-2,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_P_M_P);  
 index += size;       
 ///////////////
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
/////////////////// 
}

template <int N,int numblock, typename dataType1>
void partial_packTo3M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n2;  
 int stride = 1; 
  
 probeComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3+1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;
 probeComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;                                                                        
 probeComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size;                                                                        
 probeComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 probeComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_M_P_M);  
 index += size;       
 probeComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_P_P_M);  
 index += size;     
///////////////
 probeComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
 probeComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
////////////////// 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3+1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;
 probeComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;                                                                        
 probeComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size;                                                                        
 probeComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 probeComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_M_P_M);  
 index += size;       
 probeComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_P_P_M);  
 index += size;     
///////////////
 probeComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
 probeComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
////////////////// 
}

template <int N,int numblock, typename dataType1>
void partial_unpackFrom3M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n2;  
 int stride = 1; 
 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+2,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;
 modifyComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;                                                                         
 modifyComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size;                                                                         
 modifyComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 modifyComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_M_P_M);  
 index += size;       
 modifyComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_P_P_M);  
 index += size;     
///////////////////////
 modifyComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
////////////////////// 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+2,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;                                                                         
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size;                                                                         
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_M_P_M);  
 index += size;       
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_P_P_M);  
 index += size;     
///////////////////////
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
}

template <int N,int numblock, typename dataType1>
void partial_packMSgTo3Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp3Psend, dataType1 *temp3Msend)
{
 int index = 0;
 partial_packTo3P(lbModel,myGrid,temp3Psend,index);   
 index = 0;
 partial_packTo3M(lbModel,myGrid,temp3Msend,index); 
}

template <int N,int numblock, typename dataType1>
void partial_recvMSgFrom3Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp3Preceive, dataType1 *temp3Mreceive)
{
 int index = 0;
 partial_unpackFrom3P(lbModel,myGrid,temp3Preceive,index);  
 index = 0;
 partial_unpackFrom3M(lbModel,myGrid,temp3Mreceive,index); 
}


// Explicit declaration
template void partial_packTo1P<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void partial_unpackFrom1P<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void partial_packTo1M<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void partial_unpackFrom1M<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void partial_packMSgTo1Neb<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *, double *);
template void partial_recvMSgFrom1Neb<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *, double *);

template void partial_packTo2P<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void partial_unpackFrom2P<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void partial_packTo2M<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void partial_unpackFrom2M<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void partial_packMSgTo2Neb<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *, double *);
template void partial_recvMSgFrom2Neb<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *, double *);

template void partial_packTo3P<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void partial_unpackFrom3P<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void partial_packTo3M<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void partial_unpackFrom3M<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void partial_packMSgTo3Neb<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *, double *);
template void partial_recvMSgFrom3Neb<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *, double *);
