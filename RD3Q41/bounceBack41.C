#include"bounceBack41.h"
namespace panini {
    
template <typename dataType1> 
inline void SWAP(dataType1 &a,dataType1 &b)
{
 dataType1 temp;
 temp = a;
 a = b;
 b = temp; 
}  

    template <int N,int numblock, typename dataType1>
    inline void bounceBackG12(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID)            
    {                                                                 
        for(int nodeType=0; nodeType<2;nodeType++)
        {
            for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
            {                                                                
                for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
                {                    
                    for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
                    {
                        if(marker(i1,i2,i3,0,nodeType,0) != FLUID)
                        {
                            
                            //////////////////
                            ////  Group 1 ////
                            //////////////////
                            // ZZP
                            if( marker(i1,i2,i3-1,0,nodeType,0)==FLUID )
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P1),
                                     myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M1));
                            }
                            
                            // ZZP2
                            if( (marker(i1,i2,i3-1,0,nodeType,0)==FLUID ) && (marker(i1,i2,i3-2,0,nodeType,0)==FLUID) )
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2),
                                     myGrid(i1  ,i2  ,i3-2,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2));
                            }                            
                            
                            // Single point correction
                            if( ((marker(i1,i2,i3-1,0,nodeType,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeType,0)==FLUID)) || ( (marker(i1,i2,i3-1,0,nodeType,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeType,0)!=FLUID) ) )
                            {
                                SWAP(myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2),
                                     myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2));
                            }   
                            
                            
                            // ZP2Z
                            if( (marker(i1,i2-1,i3,0,nodeType,0)==FLUID ) && (marker(i1,i2-2,i3,0,nodeType,0)==FLUID) )
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO),
                                     myGrid(i1  ,i2-2,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO));
                            }                            
                            
                            // Single point correction
                            if( ((marker(i1,i2-1,i3,0,nodeType,0)==FLUID ) && (marker(i1,i2+1,i3,0,nodeType,0)==FLUID)) || ( (marker(i1,i2-1,i3,0,nodeType,0)==FLUID ) && (marker(i1,i2+1,i3,0,nodeType,0)!=FLUID) ) )
                            {
                                SWAP(myGrid(i1  ,i2+1,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO),
                                     myGrid(i1  ,i2-1,i3 ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO));
                            }   
                            
                            // P2ZZ
                            if( (marker(i1-1,i2,i3,0,nodeType,0)==FLUID ) && (marker(i1-2,i2,i3,0,nodeType,0)==FLUID) )
                            {
                                SWAP(myGrid(i1  ,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO),
                                     myGrid(i1-2,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO));
                            }                            
                            
                            // Single point correction
                            if( ((marker(i1-1,i2,i3,0,nodeType,0)==FLUID ) && (marker(i1+1,i2,i3,0,nodeType,0)==FLUID)) || ( (marker(i1-1,i2,i3,0,nodeType,0)==FLUID ) && (marker(i1+1,i2,i3,0,nodeType,0)!=FLUID) ) )
                            {
                                SWAP(myGrid(i1+1,i2 ,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO),
                                     myGrid(i1-1,i2 ,i3 ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO));
                            }   
                            
                            
                            //////////////////
                            ////  Group 2 ////
                            //////////////////                            
                            // ZZM
                            if( marker(i1,i2,i3+1,0,nodeType,0)==FLUID )
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M1),
                                     myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P1));
                            }                            
                            
                            // ZZM2
                            if( (marker(i1,i2,i3+1,0,nodeType,0)==FLUID ) && (marker(i1,i2,i3+2,0,nodeType,0)==FLUID) )
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2),
                                     myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2));
                            }                            
                            
                            // Single point correction 
                            // ((marker(i1,i2,i3-1,0,nodeType,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeType,0)==FLUID)) ----------------> Not done                             
                            // I think this is not necessary as this was done for ZZP2 already, those populations are already swapped 
                            
                            
                            if( (marker(i1,i2,i3-1,0,nodeType,0)!=FLUID ) && (marker(i1,i2,i3+1,0,nodeType,0)==FLUID) )
                            {
                                SWAP(myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2),
                                     myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2));
                            }                              
                            
                            
                            // ZM2Z
                            if( (marker(i1,i2+1,i3,0,nodeType,0)==FLUID) && (marker(i1,i2+2,i3,0,nodeType,0)==FLUID) )
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO),
                                     myGrid(i1  ,i2+2,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO));
                            }                            
                            
                            // Single point correction  
                            if( (marker(i1,i2-1,i3,0,nodeType,0)!=FLUID) && (marker(i1,i2+1,i3,0,nodeType,0)==FLUID) )
                            {
                                SWAP(myGrid(i1  ,i2-1,i3 ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO),
                                     myGrid(i1  ,i2+1,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO));
                            }  
                            
                            // M2ZZ
                            if( (marker(i1+1,i2,i3,0,nodeType,0)==FLUID) && (marker(i1+2,i2,i3,0,nodeType,0)==FLUID) )
                            {
                                SWAP(myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO),
                                     myGrid(i1+2,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO));
                            }                            
                            
                            // Single point correction  
                            if( (marker(i1-1,i2,i3,0,nodeType,0)!=FLUID) && (marker(i1+1,i2,i3,0,nodeType,0)==FLUID) )
                            {
                                SWAP(myGrid(i1-1,i2 ,i3 ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO),
                                     myGrid(i1+1,i2 ,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO));
                            }          
                        }
                    }
                }
            }
        }
    }
    
    template <int N,int numblock, typename dataType1>
    inline void bounceBackG34(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID)                  
    {                                                                                                                                                       
        
        for(int nodeType =0; nodeType<2;nodeType++)                                                                                                         
        {                                                                                                                                                   
            for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                                                                                               
            {                                                                                                                                               
                for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                                                                                           
                {                                                                                                                                           
                    for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                                                                           
                    {                                                                                                                           
                        if(marker(i1,i2,i3,0,nodeType,0) != FLUID)
                        {
                            //////////////////                                                                            
                            ////  Group 3 ////                                                                            
                            //////////////////                                                                           
                            //  G3_DV_ZERO_P1_P1 //  G3_DV_ZERO_M1_P1    //  G3_DV_P1_ZERO_P1   //  G3_DV_M1_ZERO_P1
                            if(marker(i1  ,i2-1,i3-1,0,nodeType,0)==FLUID)
                            {
                                SWAP( myGrid(i1 ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1) ,
                                      myGrid(i1 ,i2-1,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1));    
                            }
                            
                            if(marker(i1  ,i2+1,i3-1,0,nodeType,0)==FLUID)
                            {
                                SWAP( myGrid(i1 ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1),
                                      myGrid(i1 ,i2+1,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1));    
                            }
                            
                            if(marker(i1-1,i2  ,i3-1,0,nodeType,0)==FLUID)
                            {
                                SWAP( myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1),
                                      myGrid(i1-1,i2  ,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1));
                            }
                            
                            if(marker(i1+1,i2  ,i3-1,0,nodeType,0)==FLUID)
                            {
                                SWAP( myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1),
                                      myGrid(i1+1,i2  ,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1));
                            }
                            
                            //////////////////                                                                                                     
                            ////  Group 4 ////                                                                                                     
                            //////////////////                                                                                                      
                            //  G4_DV_ZERO_M1_M1  //  G4_DV_ZERO_P1_M1  //  G4_DV_M1_ZERO_M1    //  G4_DV_P1_ZERO_M1                            
                            if(marker(i1  ,i2+1,i3+1,0,nodeType,0)==FLUID)                                                                                    
                            {                                                                                                                             
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1),                                                 
                                     myGrid(i1  ,i2+1,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1));                                                
                            }                                                                                                                             
                            
                            if(marker(i1  ,i2-1,i3+1,0,nodeType,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1),
                                     myGrid(i1  ,i2-1,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1));
                            }
                            
                            if(marker(i1+1,i2  ,i3+1,0,nodeType,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1),
                                     myGrid(i1+1,i2  ,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1));
                            }
                            
                            
                            if(marker(i1-1,i2  ,i3+1,0,nodeType,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1),
                                     myGrid(i1-1,i2  ,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1));
                                
                            }
                        }                                                            
                    }
                }
            }
        }
    }  
    
    template <int N,int numblock, typename dataType1>
    inline void bounceBackG56(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID)            
    {                                                                                                                                                
        for(int nodeType =0; nodeType<2;nodeType++)                                                                                                  
        {                                                                                                                                            
            for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                                                                                         
            {                                                                                                                                         
                for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                                                                                     
                {                                                                                                                                     
                    for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
                    {
                        if(marker(i1,i2,i3,0,nodeType,0)!=FLUID)
                        {                                                                                               
                            //////////////////                                                                           
                            ////  Group 5 ////                                                                          
                            //////////////////                                                                           
                            //G5_DV_P1_P1_ZERO   //G5_DV_M1_P1_ZERO   //G5_DV_ZERO_P1_ZERO   //G5_DV_P1_ZERO_ZERO
                            if(marker(i1-1,i2-1,i3  ,0,nodeType,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO),
                                     myGrid(i1-1,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO));
                            }
                            
                            if(marker(i1+1,i2-1,i3  ,0,nodeType,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO),
                                     myGrid(i1+1,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO));
                            }
                            
                            if(marker(i1  ,i2-1,i3  ,0,nodeType,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_ZERO_P1_ZERO),
                                     myGrid(i1  ,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_ZERO_M1_ZERO));
                            }
                            
                            if(marker(i1-1,i2  ,i3  ,0,nodeType,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_ZERO_ZERO),
                                     myGrid(i1-1,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_ZERO_ZERO));
                            } 
                            
                            //////////////////                                                                           
                            ////  Group 6 ////                                                                           
                            //////////////////                                                                           
                            //G6_DV_M1_M1_ZERO   //G6_DV_P1_M1_ZERO    //G6_DV_ZERO_M1_ZERO   //G6_DV_M1_ZERO_ZERO
                            if(marker(i1+1,i2+1,i3  ,0,nodeType,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO),
                                     myGrid(i1+1,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO));
                            }
                            
                            if(marker(i1-1,i2+1,i3  ,0,nodeType,0) ==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO),
                                     myGrid(i1-1,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO));
                            }
                            
                            if(marker(i1  ,i2+1,i3  ,0,nodeType,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_ZERO_M1_ZERO),
                                     myGrid(i1  ,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_ZERO_P1_ZERO));
                            }
                            
                            if(marker(i1+1,i2  ,i3  ,0,nodeType,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_ZERO_ZERO),
                                     myGrid(i1+1,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_ZERO_ZERO));
                            }
                        }                                                            
                    }
                }
            }
        }
    }  
    
    template <int N,int numblock, typename dataType1>                                                                                                          
    inline void bounceBackG78(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID)                      
    {    
        ///////////////////
        // Node --> Cell //
        ///////////////////
        for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                                                                                           
        {                                                                                                                                             
            for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                                                                                         
            {                                                                                                                                            
                for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                                                                                        
                {                                                                                                                                        
                    if(marker(i1,i2,i3,0,nodeTYPE::CELL,0)!=FLUID)
                    {                                                      
                        //////////////////                                                                         
                        ////  Group 7 ////                                                                         
                        //////////////////  
                        // G7_DV_P_P_P    // G7_DV_M_P_P   // G7_DV_M_M_P   // G7_DV_P_M_P 
                        
                        // PPP
                        if(marker(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0)==FLUID)
                        { 
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P),
                                 myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M));
                            
                            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
                            {
                                SWAP(myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9 ], lbModel.G9_DV_P1_P1_P1),
                                     myGrid(i1  ,i2  ,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1));     
                            }                               
                            
                        }                                                         
                        
                        // MPP
                        if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0)==FLUID)
                        {
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P),
                                 myGrid(i1+1,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M));
                            
                            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2+1,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9 ], lbModel.G9_DV_M1_P1_P1),
                                     myGrid(i1+1,i2  ,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1));  
                            }      
                        }
                        
                        // MMP
                        if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
                        {
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P),
                                 myGrid(i1+1,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M));
                            
                            if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_M1_M1_P1),
                                     myGrid(i1+1,i2+1,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1));
                            }
                        }
                        
                        // PMP
                        if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
                        {
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P),
                                 myGrid(i1  ,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M));
                            
                            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
                            {
                                SWAP(myGrid(i1+1,i2  ,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_P1_M1_P1),
                                     myGrid(i1  ,i2+1,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1));
                            }    
                        }
                        
                        //////////////////                                                                                   
                        ////  Group 8 ////                                                                                   
                        //////////////////                                                                                   
                        
                        // Node --> Cell      // G8_DV_M_M_M  // G8_DV_P_M_M    // G8_DV_P_P_M // G8_DV_M_P_M         
                        
                        // MMM
                        if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
                        {
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M),
                                 myGrid(i1+1,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P));
                            
//                             if(marker(i1  ,i2  ,i3  ,0,nodeTYPE::NODE,0)==FLUID)
//                             {
//                                 SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1),
//                                      myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_P1_P1_P1));
//                             }
                        }                  
                        
                        // PMM
                        if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
                        {
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M),
                                 myGrid(i1  ,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P));
                            
//                             if(marker(i1+1,i2  ,i3   ,0,nodeTYPE::NODE,0)==FLUID)
//                             {
//                                 SWAP(myGrid(i1+1,i2  ,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1),
//                                      myGrid(i1  ,i2+1,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_M1_P1_P1));
//                             }
                            
                        }                  
                        
                        // PPM
                        if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
                        {
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M),
                                 myGrid(i1  ,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P));
                            
//                             if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
//                             {
//                                 SWAP(myGrid(i1+1,i2+1,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1),
//                                      myGrid(i1  ,i2  ,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_M1_M1_P1));
//                             }                           
                        }                  
                        
                        // MPM
                        if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
                        {
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M),
                                 myGrid(i1+1,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P));
                            
//                             if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
//                             {
//                                 SWAP(myGrid(i1  ,i2+1,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1),
//                                      myGrid(i1+1,i2  ,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_P1_M1_P1));
//                             }                           
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
                        //////////////////                                                                                                 
                        ////  Group 7 ////                                                                                                  
                        //////////////////  
                        // G7_DV_P_P_P           // G7_DV_P_M_P         // G7_DV_M_P_P         // G7_DV_M_M_P
                        
                        // PPP
                        if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
                        { 
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P),
                                 myGrid(i1-1,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M));
                            
                            if(marker(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_P1_P1_P1),
                                     myGrid(i1-1,i2-1,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1));
                            }                             
                        }                  
                        // MPP
                        if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
                        { 
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P),
                                 myGrid(i1  ,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M));
                            
                            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
                            {
                                SWAP(myGrid(i1-1,i2  ,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_M1_P1_P1),
                                     myGrid(i1  ,i2-1,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1));
                                
                            }             
                        }       
                        // MMP
                        if(marker(i1  ,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)
                        { 
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P),
                                 myGrid(i1  ,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M));
                            
                            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
                            {
                                SWAP(myGrid(i1-1,i2-1,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_M1_M1_P1),
                                     myGrid(i1  ,i2  ,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1));
                            }                             
                        }        
                        // PMP
                        if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)                                                                                    
                        {                                                                                                                                   
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P),                                             
                                 myGrid(i1-1,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M)); 
                            
                            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
                            {
                                SWAP(myGrid(i1  ,i2-1,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_P1_M1_P1),
                                     myGrid(i1-1,i2  ,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1));
                            }      
                            
                            
                        }                              
                        
                        //////////////////                                                                                                 
                        ////  Group 8 ////                                                                                                  
                        //////////////////   
                        // G8_DV_M_M_M           // G8_DV_M_P_M         // G8_DV_P_M_M         // G8_DV_P_P_M 
                        
                        // MMM
                        if(marker(i1  ,i2   ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
                        { 
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M),
                                 myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P));
                            
//                             if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
//                             {
//                                 SWAP(myGrid(i1-1,i2-1,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1),
//                                      myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] , lbModel.G9_DV_P1_P1_P1));
//                             }                                 
                        }                  
                        // PMM
                        if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
                        { 
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M),
                                 myGrid(i1-1,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P));
                            
//                             if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
//                             {
//                                 SWAP(myGrid(i1  ,i2-1,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1),
//                                      myGrid(i1-1,i2  ,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] , lbModel.G9_DV_M1_P1_P1));
//                             }                      
                            
                        }       
                        // PPM
                        if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
                        { 
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M),
                                 myGrid(i1-1,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P));
                            
//                             if(marker(i1  ,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)
//                             {
//                                 SWAP(myGrid(i1  ,i2  ,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1),
//                                      myGrid(i1-1,i2-1,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] , lbModel.G9_DV_M1_M1_P1));
//                             }                                 
                        }        
                        // MPM
                        if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
                        { 
                            SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M),
                                 myGrid(i1  ,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P));
                            
//                             if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)
//                             {
//                                 SWAP(myGrid(i1-1,i2  ,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1),
//                                      myGrid(i1  ,i2-1,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] , lbModel.G9_DV_P1_M1_P1));
//                             }
                        }                           
                    }
                }
            }
        }
  
        
    }           
    
    template <int N,int numblock, typename dataType1>                                                                                                          
    inline void bounceBackG910(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID)                      
    {                                                                                                                                                          
        
        for(int nodeType=0; nodeType<2;nodeType++)
        {        
            for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                                                                                    
            {                                                                                                                                    
                for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                                                                                
                {                                                                                                                                            
                    for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                                                                                        
                    {                                                                                                                                        
                        if(marker(i1  ,i2  ,i3  ,0,nodeType,0)!=FLUID )
                        {            
                            //////////////////                                                                                   
                            ////  Group 9 ////                                                                                   
                            //////////////////                             
                            // G9_DV_P1_P1_P1      // G9_DV_M1_P1_P1   // G9_DV_M1_M1_P1 // G9_DV_P1_M1_P1
                            
                            // PPP
                            if(marker(i1-1,i2-1,i3-1,0,nodeType,0)==FLUID)  
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1),
                                     myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1));        
                            }
                            
                            //MPP
                            if(marker(i1+1,i2-1,i3-1,0,nodeType,0)==FLUID)  
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1),
                                     myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1));         
                            }
                            
                            //MMP
                            if(marker(i1+1,i2+1,i3-1,0,nodeType,0)==FLUID)  
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1),
                                     myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1));       
                            }
                            
                            //PMP
                            if(marker(i1-1,i2+1,i3-1,0,nodeType,0)==FLUID)  
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1),
                                     myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1));       
                            }  
                            
                            ///////////////////                                                                                   
                            ////  Group 10 ////                                                                                   
                            ///////////////////                             
                            // G10_DV_M1_M1_M1     // G10_DV_P1_M1_M1  // G10_DV_P1_P1_M1 // G10_DV_M1_P1_M1
                            
                            //MMM
                            if(marker(i1+1,i2+1,i3+1,0,nodeType,0)==FLUID)  
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1),
                                     myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1));        
                            }                       
                            
                            //PMM
                            if(marker(i1-1,i2+1,i3+1,0,nodeType,0)==FLUID)  
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1),
                                     myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1));         
                            }                        
                            
                            //PPM
                            if(marker(i1-1,i2-1,i3+1,0,nodeType,0)==FLUID)  
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1),
                                     myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1));       
                            }                         
                            
                            //MPM
                            if(marker(i1+1,i2-1,i3+1,0,nodeType,0)==FLUID)  
                            {
                                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1),
                                     myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1));       
                            }                          
                        }    
                    }
                }
            }
        }
    }
    
    template <int N,int numblock, typename dataType1>                                                                                                          
    inline void bounceBack(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID)                      
    {
        bounceBackG12(lbModel,myGrid,marker,FLUID);   
        bounceBackG34(lbModel,myGrid,marker,FLUID);
        bounceBackG56(lbModel,myGrid,marker,FLUID); 
        bounceBackG78(lbModel,myGrid,marker,FLUID);
        bounceBackG910(lbModel,myGrid,marker,FLUID);     
    }
    
}
