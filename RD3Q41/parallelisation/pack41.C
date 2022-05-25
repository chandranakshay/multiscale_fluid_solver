#include"pack41.h"

PANINI_REAL* allocateForPacking(int m1, int m2)
{
 PANINI_REAL *temp;  
 temp = new PANINI_REAL[m1*m2*84];
 return temp;
} 

template <int N,int numblock, typename dataType1>
void packTo1P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n2*myGrid.n3;  
 int stride = 1; 

 /////////////
 //   NODE  //
 ///////////// 
 
 probeComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;

 
 
 ////////////////
 //   Group 0  //
 ////////////////  
 probeComputationalGridPlane23Node(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 probeComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 probeComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 probeComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
 
 

 /////////////
 //   CELL  //
 ///////////// 

 probeComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size; 

 ////////////////
 //   Group 0  //
 ////////////////  
 probeComputationalGridPlane23Cell(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 probeComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
}

template <int N,int numblock, typename dataType1> 
void unpackFrom1P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n2*myGrid.n3;  
 int stride = 1; 
 /////////////
 //   NODE  //
 ///////////// 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-2,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 modifyComputationalGridPlane23Node(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB1-1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 modifyComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     


 /////////////
 //   CELL  //
 ///////////// 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-2,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size; 

 ////////////////
 //   Group 0  //
 ////////////////  
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB1-1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1-1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     
}

template <int N,int numblock, typename dataType1>
void packTo1M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n2*myGrid.n3;  
 int stride = 1; 

 /////////////
 //   NODE  //
 ///////////// 
 
 probeComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1+1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;

 
 
 ////////////////
 //   Group 0  //
 ////////////////  
 probeComputationalGridPlane23Node(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 probeComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 probeComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 probeComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
 
 

 /////////////
 //   CELL  //
 ///////////// 

 probeComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1+1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size; 

 ////////////////
 //   Group 0  //
 ////////////////  
 probeComputationalGridPlane23Cell(myGrid,lbModel.G0,tempArray,index,stride,myGrid.nB1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 probeComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB1,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
}

template <int N,int numblock, typename dataType1> 
void unpackFrom1M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n2*myGrid.n3;  
 int stride = 1; 
 /////////////
 //   NODE  //
 ///////////// 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+2,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 modifyComputationalGridPlane23Node(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE1+1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane23Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 modifyComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane23Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane23Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane23Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane23Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     


 /////////////
 //   CELL  //
 ///////////// 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+2,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE1+1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane23Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE1+1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     
}

template <int N,int numblock, typename dataType1> 
void packMSgTo1Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp1Psend, dataType1 *temp1Msend)
{
 int index = 0;
 packTo1P(lbModel,myGrid,temp1Psend,index);  
 index = 0;
 packTo1M(lbModel,myGrid,temp1Msend,index); 
}

template <int N,int numblock, typename dataType1> 
void recvMSgFrom1Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp1Preceive, dataType1 *temp1Mreceive)
{
 int index = 0;
 unpackFrom1P(lbModel,myGrid,temp1Preceive,index);  
 index = 0;
 unpackFrom1M(lbModel,myGrid,temp1Mreceive,index); 
}




template <int N,int numblock, typename dataType1>
void packTo2P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n3;  
 int stride = 1; 

 /////////////
 //   NODE  //
 ///////////// 
 
 probeComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 probeComputationalGridPlane31Node(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE2,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 probeComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 probeComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 probeComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
 
 

 /////////////
 //   CELL  //
 ///////////// 

 probeComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size; 

 ////////////////
 //   Group 0  //
 ////////////////  
 probeComputationalGridPlane31Cell(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE2,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 probeComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
}

template <int N,int numblock, typename dataType1> 
void unpackFrom2P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n3;  
 int stride = 1; 
 /////////////
 //   NODE  //
 ///////////// 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-2,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 modifyComputationalGridPlane31Node(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB2-1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 modifyComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     


 /////////////
 //   CELL  //
 ///////////// 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-2,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size; 

 ////////////////
 //   Group 0  //
 ////////////////  
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB2-1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2-1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     
}

template <int N,int numblock, typename dataType1>
void packTo2M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n3;  
 int stride = 1; 

 /////////////
 //   NODE  //
 ///////////// 
 
 probeComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2+1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 probeComputationalGridPlane31Node(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB2,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 probeComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 probeComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 probeComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
 
 

 /////////////
 //   CELL  //
 ///////////// 

 probeComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2+1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 probeComputationalGridPlane31Cell(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB2,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB2,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 probeComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB2,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB2,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB2,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB2,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB2,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB2,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB2,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB2,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB2,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
}

template <int N,int numblock, typename dataType1> 
void unpackFrom2M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n3;  
 int stride = 1; 
 /////////////
 //   NODE  //
 ///////////// 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+2,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 modifyComputationalGridPlane31Node(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE2+1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 modifyComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane31Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane31Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane31Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     


 /////////////
 //   CELL  //
 ///////////// 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+2,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE2+1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane31Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE2+1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     
}



template <int N,int numblock, typename dataType1> 
void packMSgTo2Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp2Psend, dataType1 *temp2Msend)
{
 int index = 0;
 packTo2P(lbModel,myGrid,temp2Psend,index);  
 index = 0;
 packTo2M(lbModel,myGrid,temp2Msend,index); 
}

template <int N,int numblock, typename dataType1> 
void recvMSgFrom2Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp2Preceive, dataType1 *temp2Mreceive)
{
 int index = 0;
 unpackFrom2P(lbModel,myGrid,temp2Preceive,index);  
 index = 0;
 unpackFrom2M(lbModel,myGrid,temp2Mreceive,index); 
}



template <int N,int numblock, typename dataType1>
void packTo3P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n2;  
 int stride = 1; 

 /////////////
 //   NODE  //
 ///////////// 
 
 probeComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 probeComputationalGridPlane12Node(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE3,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 probeComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 probeComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 probeComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
 
 

 /////////////
 //   CELL  //
 ///////////// 

 probeComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size; 

 ////////////////
 //   Group 0  //
 ////////////////  
 probeComputationalGridPlane12Cell(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE3,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 probeComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 probeComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
}

template <int N,int numblock, typename dataType1> 
void unpackFrom3P(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n2;  
 int stride = 1; 
 /////////////
 //   NODE  //
 ///////////// 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-2,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 modifyComputationalGridPlane12Node(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB3-1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 modifyComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     


 /////////////
 //   CELL  //
 ///////////// 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-2,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB3-1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size; 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;  
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;  
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3-1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     
}

template <int N,int numblock, typename dataType1>
void packTo3M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n2;  
 int stride = 1; 

 /////////////
 //   NODE  //
 ///////////// 
 
 probeComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3+1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 probeComputationalGridPlane12Node(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB3,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size;                                                                        
 probeComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;                                                                        
 probeComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 probeComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size;                                                                        
 probeComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;                                                                        
 probeComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 probeComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3,lbModel.G3_DV_ZERO_M1_P1);
 index += size;                                                                        
 probeComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 probeComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 probeComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 probeComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 probeComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 probeComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 probeComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 probeComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 probeComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 probeComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 probeComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3,lbModel.G7_DV_M_M_P);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 probeComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 probeComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_P_P_M);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 probeComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 probeComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 probeComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 probeComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 probeComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_M1_P1_M1);  
 index += size;    
 
 

 /////////////
 //   CELL  //
 ///////////// 

 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3+1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 probeComputationalGridPlane12Cell(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nB3,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size; 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;  
 probeComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 probeComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nB3,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
  probeComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3,lbModel.G2_DV_ZERO_ZERO_M1);
  index += size; 
  probeComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3,lbModel.G2_DV_ZERO_ZERO_M2);
  index += size;  
  probeComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3,lbModel.G2_DV_ZERO_M2_ZERO);
  index += size;   
  probeComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nB3,lbModel.G2_DV_M2_ZERO_ZERO);
  index += size;   
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
  probeComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3,lbModel.G3_DV_ZERO_P1_P1);                                        
  index += size;                                                                                                                        
  probeComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3,lbModel.G3_DV_ZERO_M1_P1);
  index += size;                                                                        
  probeComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3,lbModel.G3_DV_P1_ZERO_P1);
  index += size;   
  probeComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nB3,lbModel.G3_DV_M1_ZERO_P1);
  index += size;   
 
 ////////////////                                                                       
 //   Group 4  //                                                                       
 ////////////////                                                                                                                                        
  probeComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_ZERO_M1_M1);                                        
  index += size;                                                                                                                        
  probeComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_ZERO_P1_M1);                               
  index += size;                                                                                                              
  probeComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_M1_ZERO_M1);
  index += size;   
  probeComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nB3,lbModel.G4_DV_P1_ZERO_M1);
  index += size;   
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
  probeComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
  index += size;                                                                                        
  probeComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3,lbModel.G5_DV_M1_P1_ZERO);          
  index += size;                                                                                        
  probeComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3,lbModel.G5_DV_ZERO_P1_ZERO);  
  index += size;   
  probeComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nB3,lbModel.G5_DV_P1_ZERO_ZERO);  
  index += size;   
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
  probeComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
  index += size;                                                                                                                         
  probeComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3,lbModel.G6_DV_P1_M1_ZERO);                                         
  index += size;                                                                                                                        
  probeComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3,lbModel.G6_DV_ZERO_M1_ZERO);  
  index += size;   
  probeComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nB3,lbModel.G6_DV_M1_ZERO_ZERO);  
  index += size;   
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
  probeComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
  index += size;                                                                                                                          
  probeComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3,lbModel.G7_DV_M_P_P);                                          
  index += size;                                                                                                                         
  probeComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3,lbModel.G7_DV_M_M_P);  
  index += size;   
  probeComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nB3,lbModel.G7_DV_P_M_P);  
  index += size;   
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
  probeComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
  index += size;                                                                                                                            
  probeComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_P_M_M);                                            
  index += size;                                                                                                                           
  probeComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_P_P_M);  
  index += size;   
  probeComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nB3,lbModel.G8_DV_M_P_M);  
  index += size;   
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
  probeComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
  index += size;                                                                                                                             
  probeComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3,lbModel.G9_DV_M1_P1_P1);                                             
  index += size;                                                                                                                            
  probeComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3,lbModel.G9_DV_M1_M1_P1);  
  index += size;   
  probeComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nB3,lbModel.G9_DV_P1_M1_P1);  
  index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 probeComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                              
 probeComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                             
 probeComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_P1_P1_M1);  
 index += size;   
 probeComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nB3,lbModel.G10_DV_M1_P1_M1);  
 index += size;   
}

template <int N,int numblock, typename dataType1> 
void unpackFrom3M(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *tempArray,int &index)
{ 
 int size = myGrid.n1*myGrid.n2;  
 int stride = 1; 
 /////////////
 //   NODE  //
 ///////////// 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+2,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 modifyComputationalGridPlane12Node(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE3+1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size;                                                                         
 modifyComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;                                                                         
 modifyComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane12Node(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////                                                                       
 //   Group 2  //                                                                       
 ////////////////                                                                       
 modifyComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size;                                                                         
 modifyComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;                                                                         
 modifyComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane12Node(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;                                                                         
 modifyComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane12Node(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////                                                                       
 //   Group 4  //                                                                       
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane12Node(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane12Node(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     


 /////////////
 //   CELL  //
 ///////////// 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+2,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;

 ////////////////
 //   Group 0  //
 ////////////////  
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G0 ,tempArray,index,stride,myGrid.nE3+1,lbModel.CENTER_DV_ZERO_ZERO_ZERO);
 index += size;  
 
 ////////////////
 //   Group 1  //
 //////////////// 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G1_DV_ZERO_ZERO_P1);
 index += size;                                                                         
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G1_DV_ZERO_ZERO_P2);
 index += size;                                                                         
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G1_DV_ZERO_P2_ZERO);
 index += size;   
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G1 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G1_DV_P2_ZERO_ZERO);
 index += size;   
 
 ////////////////
 //   Group 2  //
 //////////////// 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G2_DV_ZERO_ZERO_M1);
 index += size;                                                                         
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G2_DV_ZERO_ZERO_M2);
 index += size;                                                                         
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G2_DV_ZERO_M2_ZERO);
 index += size;   
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G2 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G2_DV_M2_ZERO_ZERO);
 index += size;    
 
 ////////////////
 //   Group 3  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G3_DV_ZERO_P1_P1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G3_DV_ZERO_M1_P1);
 index += size;                                                                         
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G3_DV_P1_ZERO_P1);
 index += size;   
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G3 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G3_DV_M1_ZERO_P1);
 index += size;    
 
 ////////////////
 //   Group 4  //
 ////////////////                                                                                                                                        
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_ZERO_M1_M1);                                        
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_ZERO_P1_M1);                               
 index += size;                                                                                                                               
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_M1_ZERO_M1);
 index += size;   
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G4 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G4_DV_P1_ZERO_M1);
 index += size;    
 
 ////////////////                                                                                                                                                            
 //   Group 5  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G5_DV_P1_P1_ZERO);                                                                                                                                                                   
 index += size;                                                                                                          
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G5_DV_M1_P1_ZERO);          
 index += size;                                                                                                         
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G5_DV_ZERO_P1_ZERO);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G5 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G5_DV_P1_ZERO_ZERO);  
 index += size;      
   
 ////////////////                                                                                                                                                            
 //   Group 6  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G6_DV_M1_M1_ZERO);                                                                                                                                                                       
 index += size;                                                                                                                                           
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G6_DV_P1_M1_ZERO);                                         
 index += size;                                                                                                                                         
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G6_DV_ZERO_M1_ZERO);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G6 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G6_DV_M1_ZERO_ZERO);  
 index += size;       
  
  
 ////////////////                                                                                                                                                            
 //   Group 7  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G7_DV_P_P_P);                                                                                                                                                                        
 index += size;                                                                                                                                            
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G7_DV_M_P_P);                                          
 index += size;                                                                                                                                          
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G7_DV_M_M_P);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G7 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G7_DV_P_M_P);  
 index += size;       
    
 ////////////////                                                                                                                                                            
 //   Group 8  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_M_M_M);                                                                                                                                                                          
 index += size;                                                                                                                                              
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_P_M_M);                                            
 index += size;                                                                                                                                            
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_P_P_M);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G8 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G8_DV_M_P_M);  
 index += size;       
      
 ////////////////                                                                                                                                                            
 //   Group 9  //                                                                                                                                                            
 ////////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G9_DV_P1_P1_P1);                                                                                                                                                                           
 index += size;                                                                                                                                               
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G9_DV_M1_P1_P1);                                             
 index += size;                                                                                                                                             
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G9_DV_M1_M1_P1);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G9 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G9_DV_P1_M1_P1);  
 index += size;   
 
 ///////////////                                                                                                                                                            
 //  Group 10 //                                                                                                                                                            
 ///////////////                                                                                                                                                                                                                                                                 
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_M1_M1_M1);                                                                                                                                                                            
 index += size;                                                                                                                                                
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_P1_M1_M1);                                              
 index += size;                                                                                                                                              
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_P1_P1_M1);  
 index += size;     
 modifyComputationalGridPlane12Cell(myGrid,lbModel.G10 ,tempArray,index,stride,myGrid.nE3+1,lbModel.G10_DV_M1_P1_M1);  
 index += size;     
}

template <int N,int numblock, typename dataType1> 
void packMSgTo3Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp3Psend, dataType1 *temp3Msend)
{
 int index = 0;
 packTo3P(lbModel,myGrid,temp3Psend,index);  
 index = 0;
 packTo3M(lbModel,myGrid,temp3Msend,index); 
}

template <int N,int numblock, typename dataType1> 
void recvMSgFrom3Neb(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 *temp3Preceive, dataType1 *temp3Mreceive)
{
 int index = 0;
 unpackFrom3P(lbModel,myGrid,temp3Preceive,index);  
 index = 0;
 unpackFrom3M(lbModel,myGrid,temp3Mreceive,index); 
}



// Explicit declaration
template void packTo1P<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void unpackFrom1P<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void packTo1M<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void unpackFrom1M<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void packMSgTo1Neb<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *, double *);
template void recvMSgFrom1Neb<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *, double *);

template void packTo2P<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void unpackFrom2P<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void packTo2M<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void unpackFrom2M<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void packMSgTo2Neb<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *, double *);
template void recvMSgFrom2Neb<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *, double *);

template void packTo3P<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void unpackFrom3P<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void packTo3M<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void unpackFrom3M<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *,int &);
template void packMSgTo3Neb<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *, double *);
template void recvMSgFrom3Neb<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double *, double *);
