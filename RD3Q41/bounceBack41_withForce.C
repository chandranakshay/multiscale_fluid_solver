#include<bounceBack41_withForce.h>
  
  template <typename dataType1> 
  inline void SWAP(dataType1 &a,dataType1 &b)
  {
   dataType1 temp;
   temp = a;
   a = b;
   b = temp; 
  }  
  
  
  template <int N,int numblock, typename dataType1>
  void bounceBackG12(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)        
  {                                                                 
    for(int nodeType=0; nodeType<2;nodeType++)
    {
      //          int nodeType=0; 
      for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
      {                                                                                                                                         
        for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
        {                                                                                                                                     
          for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                         
          {
            if(marker(i1,i2,i3,0,nodeType,0) != FLUID)
            {
              
              //////////////////
              ////  Group 1 ////
              //////////////////
              //ZZP
              if( marker(i1,i2,i3-1,0,nodeType,0)==FLUID )
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P1),
                     myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M1));
                ForceZ += myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M1);
              }
              
              // ZZP2
              if( (marker(i1,i2,i3-1,0,nodeType,0)==FLUID ) && (marker(i1,i2,i3-2,0,nodeType,0)==FLUID) )
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2),
                     myGrid(i1  ,i2  ,i3-2,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2));
                ForceZ += 2.0*myGrid(i1  ,i2  ,i3-2,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2);
              }                            
              
              // Single point correction
              if((marker(i1,i2,i3-1,0,nodeType,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeType,0)==FLUID)) 
              {
                SWAP(myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2),
                     myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2));
                ForceZ += 2.0*myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2);
                ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2);
              }   
              
              if( (marker(i1,i2,i3-1,0,nodeType,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeType,0)!=FLUID) )
              {
                SWAP(myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2),
                     myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2));
                ForceZ += 2.0*myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2);
              }
              
              // ZP2Z
              if( (marker(i1,i2-1,i3,0,nodeType,0)==FLUID ) && (marker(i1,i2-2,i3,0,nodeType,0)==FLUID) )
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO),
                     myGrid(i1  ,i2-2,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO));
                ForceY  += 2.0*myGrid(i1  ,i2-2,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO);
              }                            
              
              // Single point correction
              if((marker(i1,i2-1,i3,0,nodeType,0)==FLUID ) && (marker(i1,i2+1,i3,0,nodeType,0)==FLUID)) 
              {
                SWAP(myGrid(i1  ,i2+1,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO),
                     myGrid(i1  ,i2-1,i3 ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO));
                ForceY += 2.0*myGrid(i1  ,i2-1,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO);
                ForceY -= 2.0*myGrid(i1  ,i2+1,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO);
              }   
              
              if( (marker(i1,i2-1,i3,0,nodeType,0)==FLUID ) && (marker(i1,i2+1,i3,0,nodeType,0)!=FLUID) )
              {
                SWAP(myGrid(i1  ,i2+1,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO),
                     myGrid(i1  ,i2-1,i3 ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO));
                ForceY += 2.0*myGrid(i1  ,i2-1,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO);
              }
              
              // P2ZZ
              if( (marker(i1-1,i2,i3,0,nodeType,0)==FLUID ) && (marker(i1-2,i2,i3,0,nodeType,0)==FLUID) )
              {
                SWAP(myGrid(i1  ,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO),
                     myGrid(i1-2,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO));
                ForceX += 2.0*myGrid(i1-2,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO);
              }                            
              
              // Single point correction
              if( (marker(i1-1,i2,i3,0,nodeType,0)==FLUID ) && (marker(i1+1,i2,i3,0,nodeType,0)==FLUID) )
              {
                SWAP(myGrid(i1+1,i2 ,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO),
                     myGrid(i1-1,i2 ,i3 ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO));
                ForceX += 2.0*myGrid(i1-1,i2 ,i3 ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO);
                ForceX -= 2.0*myGrid(i1+1,i2 ,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO);                                
              }   
              
              
              if( (marker(i1-1,i2,i3,0,nodeType,0)==FLUID ) && (marker(i1+1,i2,i3,0,nodeType,0)!=FLUID) ) 
              {
                SWAP(myGrid(i1+1,i2 ,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO),
                     myGrid(i1-1,i2 ,i3 ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO)); 
                ForceX += 2.0*myGrid(i1-1,i2 ,i3 ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO);                                
              }
              
              
              //////////////////
              ////  Group 2 ////
              //////////////////                            
              //                             ZZM
              if( marker(i1,i2,i3+1,0,nodeType,0)==FLUID )
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M1),
                     myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P1));
                ForceZ -= myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P1);
              }                            
              
              // ZZM2
              if( (marker(i1,i2,i3+1,0,nodeType,0)==FLUID ) && (marker(i1,i2,i3+2,0,nodeType,0)==FLUID) )
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2),
                     myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2));
                ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2);
              }                            
              
              // Single point correction 
              // ((marker(i1,i2,i3-1,0,nodeType,0)==FLUID ) && (marker(i1,i2,i3+1,0,nodeType,0)==FLUID)) ----------------> Not done                             
              // I think this is not necessary as this was done for ZZP2 already, those populations are already swapped 
              
              
              if( (marker(i1,i2,i3-1,0,nodeType,0)!=FLUID ) && (marker(i1,i2,i3+1,0,nodeType,0)==FLUID) )
              {
                SWAP(myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2),
                     myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2));
                ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+1,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2);
              }                              
              
              
              // ZM2Z
              if( (marker(i1,i2+1,i3,0,nodeType,0)==FLUID) && (marker(i1,i2+2,i3,0,nodeType,0)==FLUID) )
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO),
                     myGrid(i1  ,i2+2,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO));
                ForceY -= 2.0*myGrid(i1  ,i2+2,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO);
              }                            
              
              // Single point correction  
              if( (marker(i1,i2-1,i3,0,nodeType,0)!=FLUID) && (marker(i1,i2+1,i3,0,nodeType,0)==FLUID) )
              {
                SWAP(myGrid(i1  ,i2-1,i3 ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO),
                     myGrid(i1  ,i2+1,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO));
                ForceY -= 2.0*myGrid(i1  ,i2+1,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO);
              }  
              
              // M2ZZ
              if( (marker(i1+1,i2,i3,0,nodeType,0)==FLUID) && (marker(i1+2,i2,i3,0,nodeType,0)==FLUID) )
              {
                SWAP(myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO),
                     myGrid(i1+2,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO));
                ForceX -= 2.0*myGrid(i1+2,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO);
              }                            
              
              // Single point correction  
              if( (marker(i1-1,i2,i3,0,nodeType,0)!=FLUID) && (marker(i1+1,i2,i3,0,nodeType,0)==FLUID) )
              {
                SWAP(myGrid(i1-1,i2 ,i3 ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO),
                     myGrid(i1+1,i2 ,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO));
                ForceX -= 2.0*myGrid(i1+1,i2 ,i3 ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO);
              }          
            }
          }
        }
      }
    }
  }

   ///////////////////////////////////////////////////////////////////////
   /// \brief bounceBackG12Solid
   /////////////////////////////////////////////////////////////////////////  
  template <int N,int numblock, typename dataType1>
  void bounceBackG12Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)        
  {                                                                 
    for(int nodeType=0; nodeType<2;nodeType++)
    {
      for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
      {                                                                                                                                         
        for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
        {                                                                                                                                     
          for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                            
          {
            if(marker(i1,i2,i3,0,nodeType,0) == FLUID)
            {
              
              //////////////////
              ////  Group 1 ////
              //////////////////
              //ZZP
              if( marker(i1,i2,i3+1,0,nodeType,0) != FLUID )
              {
                SWAP(myGrid(i1  ,i2  ,i3+1  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P1),
                     myGrid(i1  ,i2  ,i3    ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M1);
              }
              
              // ZZP2
              if( marker(i1,i2,i3+2,0,nodeType,0) != FLUID) 
              {
                SWAP(myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2));
               
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceZ += 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2);
              }                            
              
              if( (marker(i1,i2,i3+1,0,nodeType,0) != FLUID && marker(i1,i2,i3+2,0,nodeType,0) == FLUID) )
              {
                SWAP(myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2));
               
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceZ += 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2);
                  ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2);
                }
              }  
              
              // ZP2Z
              if( (marker(i1,i2+2,i3,0,nodeType,0)!=FLUID) )
              {
                SWAP(myGrid(i1  ,i2+2,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceY  += 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO);
              }                            
              
              if( (marker(i1,i2+1,i3,0,nodeType,0)!=FLUID) && (marker(i1,i2+2,i3,0,nodeType,0)==FLUID))
              {
                SWAP(myGrid(i1  ,i2+2,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceY  += 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO);
                  ForceY  -= 2.0*myGrid(i1  ,i2+2,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO);
                }  
              }    
              
              // P2ZZ
              if( (marker(i1+2,i2,i3,0,nodeType,0)!=FLUID) )
              {
                SWAP(myGrid(i1+2,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO),
                     myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceX += 2.0*myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO);
              }                            
              
              if( (marker(i1+1,i2,i3,0,nodeType,0)!=FLUID) &&  (marker(i1+2,i2,i3,0,nodeType,0)==FLUID) )
              {
                SWAP(myGrid(i1+2,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO),
                     myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX += 2.0*myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO);
                  ForceX -= 2.0*myGrid(i1+2,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO);
                }
              }                                
              
              //////////////////
              ////  Group 2 ////
              //////////////////                            
              //                             ZZM
              if( marker(i1,i2,i3-1,0,nodeType,0)!=FLUID )
              {
                SWAP(myGrid(i1  ,i2  ,i3-1,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M1),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P1);
              }                            
              
              // ZZM2
              if(marker(i1,i2,i3-2,0,nodeType,0)!=FLUID) 
              {
                SWAP(myGrid(i1  ,i2  ,i3-2,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceZ -= 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2);
              }                            
              
              // ZM2Z
              if((marker(i1,i2-2,i3,0,nodeType,0)!=FLUID) )
              {
                SWAP(myGrid(i1  ,i2-2,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceY -= 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO);
              }                            
              
              // M2ZZ
              if((marker(i1-2,i2,i3,0,nodeType,0)!=FLUID) )
              {
                SWAP(myGrid(i1-2,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO),
                     myGrid(i1  ,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceX -= 2.0*myGrid(i1  ,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO);
              } 
            }
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void bounceBackG34(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                     
  {                                                                                                                                                       
    
    for(int nodeType =0; nodeType<2;nodeType++)                                                                                                         
    {                                           
      //       int nodeType(0);
      for (int i3 = myGrid.nB3; i3 <= myGrid.nE3; i3++)                                                                                         
      {                                                                                                                                         
        for (int i2 = myGrid.nB2; i2 <= myGrid.nE2; i2++)                                                                                         
        {                                                                                                                                     
          for (int i1 = myGrid.nB1; i1<= myGrid.nE1; i1++)                                                                                         
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
                ForceY += myGrid(i1 ,i2-1,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1);
                ForceZ += myGrid(i1 ,i2-1,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1);
              }
              
              if(marker(i1  ,i2+1,i3-1,0,nodeType,0)==FLUID)
              {
                SWAP( myGrid(i1 ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1),
                      myGrid(i1 ,i2+1,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1));    
                ForceY -= myGrid(i1 ,i2+1,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1);
                ForceZ += myGrid(i1 ,i2+1,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1);
              }
              
              if(marker(i1-1,i2  ,i3-1,0,nodeType,0)==FLUID)
              {
                SWAP( myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1),
                      myGrid(i1-1,i2  ,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1));
                ForceX += myGrid(i1-1,i2  ,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1);
                ForceZ += myGrid(i1-1,i2  ,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1);
              }
              
              if(marker(i1+1,i2  ,i3-1,0,nodeType,0)==FLUID)
              {
                SWAP( myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1),
                      myGrid(i1+1,i2  ,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1));
                ForceX -= myGrid(i1+1,i2  ,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1);
                ForceZ += myGrid(i1+1,i2  ,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1);
              }
              
              //////////////////                                                                                                     
              ////  Group 4 ////                                                                                                     
              //////////////////                                                                                                      
              //  G4_DV_ZERO_M1_M1  //  G4_DV_ZERO_P1_M1  //  G4_DV_M1_ZERO_M1    //  G4_DV_P1_ZERO_M1                            
              if(marker(i1  ,i2+1,i3+1,0,nodeType,0)==FLUID)                                                                                    
              {                                                                                                                             
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1),                                                 
                     myGrid(i1  ,i2+1,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1));                                                
                ForceY -= myGrid(i1  ,i2+1,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1);
                ForceZ -= myGrid(i1  ,i2+1,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1);
              }                                                                                                                             
              
              if(marker(i1  ,i2-1,i3+1,0,nodeType,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1),
                     myGrid(i1  ,i2-1,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1));
                ForceY += myGrid(i1  ,i2-1,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1);
                ForceZ -= myGrid(i1  ,i2-1,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1);
              }
              
              if(marker(i1+1,i2  ,i3+1,0,nodeType,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1),
                     myGrid(i1+1,i2  ,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1));
                ForceX -= myGrid(i1+1,i2  ,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1);
                ForceZ -= myGrid(i1+1,i2  ,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1);
              }
              
              
              if(marker(i1-1,i2  ,i3+1,0,nodeType,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1),
                     myGrid(i1-1,i2  ,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1));
                ForceX += myGrid(i1-1,i2  ,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1);
                ForceZ -= myGrid(i1-1,i2  ,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1);
              }
            }                                                            
          }
        }
      }
    }
  }  
  
  template <int N,int numblock, typename dataType1>
  void bounceBackG34Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                     
  {                                                                                                                                                       
    
    for(int nodeType =0; nodeType<2;nodeType++)                                                                                                         
    {                                           
      //       int nodeType(0);
      for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
      {                                                                                                                                         
        for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
        {                                                                                                                                     
          for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                             
          {                                                                                                                       
            if(marker(i1,i2,i3,0,nodeType,0) == FLUID)
            {
              //////////////////                                                                            
              ////  Group 3 ////                                                                            
              //////////////////                                                                           
              //  G3_DV_ZERO_P1_P1 //  G3_DV_ZERO_M1_P1    //  G3_DV_P1_ZERO_P1   //  G3_DV_M1_ZERO_P1
              if(marker(i1  ,i2+1,i3+1,0,nodeType,0) != FLUID)
              {
                SWAP( myGrid(i1 ,i2+1,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1) ,
                      myGrid(i1 ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1)  );    
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceY += myGrid(i1 ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1);
                  ForceZ += myGrid(i1 ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1);
                }
              }
              
              if(marker(i1  ,i2-1,i3+1,0,nodeType,0)!=FLUID)
              {
                SWAP( myGrid(i1 ,i2-1,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1),
                      myGrid(i1 ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1));    
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceY -= myGrid(i1 ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1);
                  ForceZ += myGrid(i1 ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1);
                }
              }
              
              if(marker(i1+1,i2  ,i3+1,0,nodeType,0)!=FLUID)
              {
                SWAP( myGrid(i1+1,i2  ,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1),
                      myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1);
                  ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1);
                }
              }
              
              if(marker(i1-1,i2  ,i3+1,0,nodeType,0)!=FLUID)
              {
                SWAP( myGrid(i1-1,i2  ,i3+1,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1),
                      myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1);
                  ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1);
                }
              }
              
              //////////////////                                                                                                     
              ////  Group 4 ////                                                                                                     
              //////////////////                                                                                                      
              //  G4_DV_ZERO_M1_M1  //  G4_DV_ZERO_P1_M1  //  G4_DV_M1_ZERO_M1    //  G4_DV_P1_ZERO_M1                            
              if(marker(i1  ,i2-1,i3-1,0,nodeType,0)!=FLUID)                                                                                    
              {                                                                                                                             
                SWAP(myGrid(i1  ,i2-1,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1),                                                 
                     myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1));                                                
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1);
                  ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1);
                }                                                                                                                             
              }                            
              if(marker(i1  ,i2+1,i3-1,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2+1,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1);
                  ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1);
                }
              }                            
              if(marker(i1-1,i2  ,i3-1,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1-1,i2  ,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1);
                  ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1);
                }
              }                            
              
              if(marker(i1+1,i2  ,i3-1,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1+1,i2  ,i3-1,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1);
                  ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1);
                }
              }
            }                                                            
          }
        }
      }
    }
  }  
  
  template <int N,int numblock, typename dataType1>
  void bounceBackG56(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)               
  {                                                                                                                                                
    for(int nodeType =0; nodeType<2;nodeType++)                                                                                                  
    {       
      //    int nodeType =0;
      for (int i3 = myGrid.nB3; i3 <= myGrid.nE3; i3++)                                                                                         
      {                                                                                                                                         
        for (int i2 = myGrid.nB2; i2 <= myGrid.nE2; i2++)                                                                                         
        {                                                                                                                                     
          for (int i1 = myGrid.nB1; i1<= myGrid.nE1; i1++)                                                                                         
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
                ForceX += myGrid(i1-1,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO);
                ForceY += myGrid(i1-1,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO);
              }
              
              if(marker(i1+1,i2-1,i3  ,0,nodeType,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO),
                     myGrid(i1+1,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO));
                ForceX -= myGrid(i1+1,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO);
                ForceY += myGrid(i1+1,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO);
              }
              
              if(marker(i1  ,i2-1,i3  ,0,nodeType,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_ZERO_P1_ZERO),
                     myGrid(i1  ,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_ZERO_M1_ZERO));
                ForceY += myGrid(i1  ,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_ZERO_M1_ZERO);
              }
              
              if(marker(i1-1,i2  ,i3  ,0,nodeType,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_ZERO_ZERO),
                     myGrid(i1-1,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_ZERO_ZERO));
                ForceX += myGrid(i1-1,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_ZERO_ZERO);
              } 
              
              //////////////////                                                                           
              ////  Group 6 ////                                                                           
              //////////////////                                                                           
              //G6_DV_M1_M1_ZERO   //G6_DV_P1_M1_ZERO    //G6_DV_ZERO_M1_ZERO   //G6_DV_M1_ZERO_ZERO
              if(marker(i1+1,i2+1,i3  ,0,nodeType,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO),
                     myGrid(i1+1,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO));
                ForceX -= myGrid(i1+1,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO);
                ForceY -= myGrid(i1+1,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO); 
              }
              
              if(marker(i1-1,i2+1,i3  ,0,nodeType,0) ==FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO),
                     myGrid(i1-1,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO));
                ForceX += myGrid(i1-1,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO);
                ForceY -= myGrid(i1-1,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO); 
              }
              
              if(marker(i1  ,i2+1,i3  ,0,nodeType,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_ZERO_M1_ZERO),
                     myGrid(i1  ,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_ZERO_P1_ZERO));
                ForceY -= myGrid(i1  ,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_ZERO_P1_ZERO);
              }
              
              if(marker(i1+1,i2  ,i3  ,0,nodeType,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_ZERO_ZERO),
                     myGrid(i1+1,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_ZERO_ZERO));
                ForceX -= myGrid(i1+1,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_ZERO_ZERO);
              }
            }                                                            
          }
        }
      }
    }
  }  
  
  template <int N,int numblock, typename dataType1>
  void bounceBackG56Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)               
  {                                                                                                                                                
    for(int nodeType =0; nodeType<2;nodeType++)                                                                                                  
    {       
      //    int nodeType =0;
      for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
      {                                                                                                                                         
        for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
        {                                                                                                                                     
          for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                           
          {
            if(marker(i1,i2,i3,0,nodeType,0)==FLUID)
            {                                                                                               
              //////////////////                                                                           
              ////  Group 5 ////                                                                          
              //////////////////                                                                           
              //G5_DV_P1_P1_ZERO   //G5_DV_M1_P1_ZERO   //G5_DV_ZERO_P1_ZERO   //G5_DV_P1_ZERO_ZERO
              if(marker(i1+1,i2+1,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1+1,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO);
                  ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO);
                }
              }
              
              if(marker(i1-1,i2+1,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1-1,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO);
                  ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO);
                }
              }
              
              if(marker(i1  ,i2+1,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2+1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_ZERO_P1_ZERO),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_ZERO_M1_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_ZERO_M1_ZERO);
                }
              }
              
              if(marker(i1+1,i2  ,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1+1,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_ZERO_ZERO),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_ZERO_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_ZERO_ZERO);
                } 
              }
              
              //////////////////                                                                           
              ////  Group 6 ////                                                                           
              //////////////////                                                                           
              //G6_DV_M1_M1_ZERO   //G6_DV_P1_M1_ZERO    //G6_DV_ZERO_M1_ZERO   //G6_DV_M1_ZERO_ZERO
              if(marker(i1-1,i2-1,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1-1,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO);
                  ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO); 
                }
              }
              
              if(marker(i1+1,i2-1,i3  ,0,nodeType,0) !=FLUID)
              {
                SWAP(myGrid(i1+1,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO));

                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO);
                  ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO); 
                }
              }
              
              if(marker(i1  ,i2-1,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2-1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_ZERO_M1_ZERO),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_ZERO_P1_ZERO));

                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_ZERO_P1_ZERO);
              }
              
              if(marker(i1-1,i2  ,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1-1,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_ZERO_ZERO),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_ZERO_ZERO));

                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_ZERO_ZERO);
              }
            }                                                            
          }
        }
      }
    }
  }  
  
  template <int N,int numblock, typename dataType1>                                                                                                          
  void bounceBackG78(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                        
  {    
    ///////////////////
    // Node --> Cell //
    ///////////////////
    for (int i3 = myGrid.nB3; i3 <= myGrid.nE3; i3++)                                                                                         
    {                                                                                                                                         
      for (int i2 = myGrid.nB2; i2 <= myGrid.nE2; i2++)                                                                                         
      {                                                                                                                                     
        for (int i1 = myGrid.nB1; i1<= myGrid.nE1; i1++)                                                                                         
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
              ForceX += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);
              ForceY += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);                       
              ForceZ += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M); 
              
              if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
                SWAP(myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9 ], lbModel.G9_DV_P1_P1_P1),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1));  
                
                ForceX -=  myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9 ], lbModel.G9_DV_P1_P1_P1);
                ForceY -=  myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9 ], lbModel.G9_DV_P1_P1_P1);
                ForceZ -=  myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9 ], lbModel.G9_DV_P1_P1_P1);
                
                ForceX +=  myGrid(i1  ,i2  ,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1);
                ForceY +=  myGrid(i1  ,i2  ,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1);
                ForceZ +=  myGrid(i1  ,i2  ,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1);                           
              }                               
              
            }                                                         
            
            // MPP
            if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P),
                   myGrid(i1+1,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M));
              ForceX -= 0.5*myGrid(i1+1,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);
              ForceY += 0.5*myGrid(i1+1,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);                       
              ForceZ += 0.5*myGrid(i1+1,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M); 
              
              if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2+1,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9 ], lbModel.G9_DV_M1_P1_P1),
                     myGrid(i1+1,i2  ,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1));  
                
                ForceX +=  myGrid(i1  ,i2+1,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9 ], lbModel.G9_DV_M1_P1_P1);
                ForceY -=  myGrid(i1  ,i2+1,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9 ], lbModel.G9_DV_M1_P1_P1);
                ForceZ -=  myGrid(i1  ,i2+1,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9 ], lbModel.G9_DV_M1_P1_P1);
                
                ForceX -=  myGrid(i1+1,i2  ,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1);
                ForceY +=  myGrid(i1+1,i2  ,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1);
                ForceZ +=  myGrid(i1+1,i2  ,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1);                           
                
              }      
            }
            
            // MMP
            if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P),
                   myGrid(i1+1,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M));
              ForceX -= 0.5*myGrid(i1+1,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);
              ForceY -= 0.5*myGrid(i1+1,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);                       
              ForceZ += 0.5*myGrid(i1+1,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M); 
              
              if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_M1_M1_P1),
                     myGrid(i1+1,i2+1,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1));
                
                ForceX +=  myGrid(i1  ,i2  ,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_M1_M1_P1);
                ForceY +=  myGrid(i1  ,i2  ,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_M1_M1_P1);
                ForceZ -=  myGrid(i1  ,i2  ,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_M1_M1_P1);
                
                ForceX -=  myGrid(i1+1,i2+1,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1);
                ForceY -=  myGrid(i1+1,i2+1,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1);
                ForceZ +=  myGrid(i1+1,i2+1,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1);                           
                
              }
            }
            
            // PMP
            if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0)==FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P),
                   myGrid(i1  ,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M));
              ForceX += 0.5*myGrid(i1  ,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);
              ForceY -= 0.5*myGrid(i1  ,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);                       
              ForceZ += 0.5*myGrid(i1  ,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M); 
              
              if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
              {
                SWAP(myGrid(i1+1,i2  ,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_P1_M1_P1),
                     myGrid(i1  ,i2+1,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1));
                
                ForceX -=  myGrid(i1+1,i2  ,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_P1_M1_P1);
                ForceY +=  myGrid(i1+1,i2  ,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_P1_M1_P1);
                ForceZ -=  myGrid(i1+1,i2  ,i3+1,lbModel.G9 ,myGrid.node[lbModel.G9] , lbModel.G9_DV_P1_M1_P1);
                
                ForceX +=  myGrid(i1  ,i2+1,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1);
                ForceY -=  myGrid(i1  ,i2+1,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1);
                ForceZ +=  myGrid(i1  ,i2+1,i3  ,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1);                               
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
              
              ForceX -= 0.5*myGrid(i1+1,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);
              ForceY -= 0.5*myGrid(i1+1,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);                       
              ForceZ -= 0.5*myGrid(i1+1,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P); 
            }                  
            
            // PMM
            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M),
                   myGrid(i1  ,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P));
              ForceX += 0.5*myGrid(i1  ,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);
              ForceY -= 0.5*myGrid(i1  ,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);                       
              ForceZ -= 0.5*myGrid(i1  ,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);                             
              
            }                  
            
            // PPM
            if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M),
                   myGrid(i1  ,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P));
              ForceX += 0.5*myGrid(i1  ,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);
              ForceY += 0.5*myGrid(i1  ,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);                       
              ForceZ -= 0.5*myGrid(i1  ,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);                         
            }                  
            
            // MPM
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)==FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M),
                   myGrid(i1+1,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P));
              ForceX -= 0.5*myGrid(i1+1,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);
              ForceY += 0.5*myGrid(i1+1,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);                       
              ForceZ -= 0.5*myGrid(i1+1,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);                       
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
              ForceX += 0.5*myGrid(i1-1,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);
              ForceY += 0.5*myGrid(i1-1,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);                       
              ForceZ += 0.5*myGrid(i1-1,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);                       
              
              if(marker(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_P1_P1_P1),
                     myGrid(i1-1,i2-1,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1));
                
                ForceX -=  myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_P1_P1_P1);
                ForceY -=  myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_P1_P1_P1);
                ForceZ -=  myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_P1_P1_P1); 
                
                ForceX +=  myGrid(i1-1,i2-1,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1);
                ForceY +=  myGrid(i1-1,i2-1,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1);
                ForceZ +=  myGrid(i1-1,i2-1,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1);                            
              }                             
            }                  
            // MPP
            if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P),
                   myGrid(i1  ,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M));
              ForceX -= 0.5*myGrid(i1  ,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);
              ForceY += 0.5*myGrid(i1  ,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);                       
              ForceZ += 0.5*myGrid(i1  ,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);                       
              
              if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
              {
                SWAP(myGrid(i1-1,i2  ,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_M1_P1_P1),
                     myGrid(i1  ,i2-1,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1));
                
                ForceX +=  myGrid(i1-1,i2  ,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_M1_P1_P1);
                ForceY -=  myGrid(i1-1,i2  ,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_M1_P1_P1);
                ForceZ -=  myGrid(i1-1,i2  ,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_M1_P1_P1); 
                
                ForceX -=  myGrid(i1  ,i2-1,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1);
                ForceY +=  myGrid(i1  ,i2-1,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1);
                ForceZ +=  myGrid(i1  ,i2-1,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1); 
                
              }             
            }       
            // MMP
            if(marker(i1  ,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P),
                   myGrid(i1  ,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M));
              ForceX -= 0.5*myGrid(i1  ,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);
              ForceY -= 0.5*myGrid(i1  ,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);                       
              ForceZ += 0.5*myGrid(i1  ,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);                       
              
              if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
              {
                SWAP(myGrid(i1-1,i2-1,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_M1_M1_P1),
                     myGrid(i1  ,i2  ,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1));
                
                ForceX +=  myGrid(i1-1,i2-1,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_M1_M1_P1);
                ForceY +=  myGrid(i1-1,i2-1,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_M1_M1_P1);
                ForceZ -=  myGrid(i1-1,i2-1,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_M1_M1_P1); 
                
                ForceX -=  myGrid(i1  ,i2  ,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1);
                ForceY -=  myGrid(i1  ,i2  ,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1);
                ForceZ +=  myGrid(i1  ,i2  ,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1); 
                
              }                             
            }        
            // PMP
            if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0)==FLUID)                                                                                    
            {                                                                                                                                   
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P),                                             
                   myGrid(i1-1,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M)); 
              ForceX += 0.5*myGrid(i1-1,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);
              ForceY -= 0.5*myGrid(i1-1,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);                       
              ForceZ += 0.5*myGrid(i1-1,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);                       
              
              if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
              {
                SWAP(myGrid(i1  ,i2-1,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_P1_M1_P1),
                     myGrid(i1-1,i2  ,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1));
                
                ForceX -=  myGrid(i1  ,i2-1,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_P1_M1_P1);
                ForceY +=  myGrid(i1  ,i2-1,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_P1_M1_P1);
                ForceZ -=  myGrid(i1  ,i2-1,i3  ,lbModel.G9 ,myGrid.cell[lbModel.G9] ,lbModel.G9_DV_P1_M1_P1); 
                
                ForceX +=  myGrid(i1-1,i2  ,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1);
                ForceY -=  myGrid(i1-1,i2  ,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1);
                ForceZ +=  myGrid(i1-1,i2  ,i3-1,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1); 
                
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
              ForceX -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);
              ForceY -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);                       
              ForceZ -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);                       
            }                  
            // PMM
            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M),
                   myGrid(i1-1,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P));
              ForceX += 0.5*myGrid(i1-1,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);
              ForceY -= 0.5*myGrid(i1-1,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);                       
              ForceZ -= 0.5*myGrid(i1-1,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);                       
              
            }       
            // PPM
            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M),
                   myGrid(i1-1,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P));
              ForceX += 0.5*myGrid(i1-1,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);
              ForceY += 0.5*myGrid(i1-1,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);                       
              ForceZ -= 0.5*myGrid(i1-1,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);                       
            }        
            // MPM
            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)==FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M),
                   myGrid(i1  ,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P));
              ForceX -= 0.5*myGrid(i1  ,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);
              ForceY += 0.5*myGrid(i1  ,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);                       
              ForceZ -= 0.5*myGrid(i1  ,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);                       
            }                           
          }
        }
      }
    }
    
    
  }           
  
  template <int N,int numblock, typename dataType1>                                                                                                          
  void bounceBackG78Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                        
  {    
    ///////////////////
    // Node --> Cell //
    ///////////////////
    for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
    {                                                                                                                                         
      for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
      {                                                                                                                                     
        for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                             
        {                                                                                                                                       
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0)==FLUID)
          {                                                      
            //////////////////                                                                         
            ////  Group 7 ////                                                                         
            //////////////////  
            // G7_DV_P_P_P    // G7_DV_M_P_P   // G7_DV_M_M_P   // G7_DV_P_M_P 
            
            // PPP
            if(marker(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);
                ForceY += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);                       
                ForceZ += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M); 
              }                            
            }                                                         
            
            // MPP
            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1-1,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);
                ForceY += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);                       
                ForceZ += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M); 
              }                            
            }
            
            // MMP
            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1-1,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);
                ForceY -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);                       
                ForceZ += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M); 
              }                            
            }
            
            // PMP
            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1  ,i2-1,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);
                ForceY -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);                       
                ForceZ += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M); 
              }                            
            }
            
            //////////////////                                                                                   
            ////  Group 8 ////                                                                                   
            //////////////////                                                                                   
            
            // Node --> Cell      // G8_DV_M_M_M  // G8_DV_P_M_M    // G8_DV_P_P_M // G8_DV_M_P_M         
            
            // MMM
            if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1-1,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                
                ForceX -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);
                ForceY -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);                       
                ForceZ -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P); 
              }                  
            }                        
            // PMM
            if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1  ,i2-1,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);
                ForceY -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);                       
                ForceZ -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);                             
              }                            
            }                  
            
            // PPM
            if(marker(i1  ,i2  ,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);
                ForceY += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);                       
                ForceZ -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);                         
              }                  
            }                        
            // MPM
            if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1-1,i2  ,i3-1,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);
                ForceY += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);                       
                ForceZ -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);                       
              }                      
            }                  
          }
        }
      }
    }
    
    ////////////////////
    // Cell --> Node  //
    ////////////////////
    for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
    {                                                                                                                                         
      for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
      {                                                                                                                                     
        for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                            
        {                                                                                                                                        
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0)==FLUID)
          {                                                                        
            //////////////////                                                                                                 
            ////  Group 7 ////                                                                                                  
            //////////////////  
            // G7_DV_P_P_P           // G7_DV_P_M_P         // G7_DV_M_P_P         // G7_DV_M_M_P
            
            // PPP
            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1+1,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);
                ForceY += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);                       
                ForceZ += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);                       
              }                  
            }                        // MPP
            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2+1,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);
                ForceY += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);                       
                ForceZ += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);                       
              }                       
              
            }       
            // MMP
            if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);
                ForceY -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);                       
                ForceZ += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);                       
              }                      
              
            }        
            // PMP
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID)                                                                                    
            {                                                                                                                                   
              SWAP(myGrid(i1+1,i2  ,i3+1,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P),                                             
                   myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M)); 

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);
                ForceY -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);                       
                ForceZ += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);                       
              }                      
              
            }                              
            
            //////////////////                                                                                                 
            ////  Group 8 ////                                                                                                  
            //////////////////   
            // G8_DV_M_M_M           // G8_DV_M_P_M         // G8_DV_P_M_M         // G8_DV_P_P_M 
            
            // MMM
            if(marker(i1  ,i2   ,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);
                ForceY -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);                       
                ForceZ -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);                       
              }                  
            }                  
            // PMM
            if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1+1,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);
                ForceY -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);                       
                ForceZ -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);                       
              }                            
            }       
            // PPM
            if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1+1,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);
                ForceY += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);                       
                ForceZ -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);                       
              }                       
            }        
            // MPM
            if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2+1,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);
                ForceY += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);                       
                ForceZ -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);                       
              }                      
            }                           
          }
        }
      }
    }
    
  }           
  
  template <int N,int numblock, typename dataType1>                                                                                                          
  void bounceBackG910(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                       
  {                                                                                                                                                          
    
    for(int nodeType=0; nodeType<2;nodeType++)
    {        
      //          int nodeType=0;
      for (int i3 = myGrid.nB3; i3 <= myGrid.nE3; i3++)                                                                                         
      {                                                                                                                                         
        for (int i2 = myGrid.nB2; i2 <= myGrid.nE2; i2++)                                                                                         
        {                                                                                                                                     
          for (int i1 = myGrid.nB1; i1<= myGrid.nE1; i1++)                                                                                         
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
                ForceX += myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
                ForceY += myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
                ForceZ += myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
              }
              
              //MPP
              if(marker(i1+1,i2-1,i3-1,0,nodeType,0)==FLUID)  
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1),
                     myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1));         
                ForceX -= myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);
                ForceY += myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
                ForceZ += myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
              }
              
              //MMP
              if(marker(i1+1,i2+1,i3-1,0,nodeType,0)==FLUID)  
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1),
                     myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1));       
                ForceX -= myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);
                ForceY -= myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
                ForceZ += myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
              }
              
              //PMP
              if(marker(i1-1,i2+1,i3-1,0,nodeType,0)==FLUID)  
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1),
                     myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1));       
                ForceX += myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);
                ForceY -= myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
                ForceZ += myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
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
                ForceX -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);
                ForceY -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
                ForceZ -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
              }                       
              
              //PMM
              if(marker(i1-1,i2+1,i3+1,0,nodeType,0)==FLUID)  
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1),
                     myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1));         
                ForceX += myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);
                ForceY -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
                ForceZ -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
              }                        
              
              //PPM
              if(marker(i1-1,i2-1,i3+1,0,nodeType,0)==FLUID)  
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1),
                     myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1));       
                ForceX += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);
                ForceY += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
                ForceZ -= myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
              }                         
              
              //MPM
              if(marker(i1+1,i2-1,i3+1,0,nodeType,0)==FLUID)  
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1),
                     myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1));       
                ForceX -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);
                ForceY += myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
                ForceZ -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
              }                          
            }    
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>                                                                                                          
  void bounceBackG910Solid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                       
  {                                                                                                                                                          
    //////////
    // Node //
    //////////     
    int nodeType=nodeTYPE::NODE;
    for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
    {                                                                                                                                         
      for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
      {                                                                                                                                     
        for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                         
        {                                                                                                                                       
          if(marker(i1  ,i2  ,i3  ,0,nodeType,0)==FLUID )
          {            
            //////////////////                                                                                   
            ////  Group 9 ////                                                                                   
            //////////////////                             
            // G9_DV_P1_P1_P1      // G9_DV_M1_P1_P1   // G9_DV_M1_M1_P1 // G9_DV_P1_M1_P1
            
            // PPP
            if(marker(i1+1,i2+1,i3+1,0,nodeType,0)!=FLUID)  
            {
              SWAP(myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1)); 

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
                ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
              }
              
            }
              
              //MPP
              if(marker(i1-1,i2+1,i3+1,0,nodeType,0)!=FLUID) 
              {
                SWAP(myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1),
                     myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1));         

                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);
                  ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
                  ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
                }   
                
              }
                //MMP
                if(marker(i1-1,i2-1,i3+1,0,nodeType,0)!=FLUID)   
                {
                  SWAP(myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1),
                       myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1));       

                  if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  {
                    ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);
                    ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
                    ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
                  }
                  
                }
                  //PMP
                  if(marker(i1+1,i2-1,i3+1,0,nodeType,0)!=FLUID)  
                  {
                    SWAP(myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1),
                         myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1));       

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);
                      ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
                      ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
                    }  
                  }                        
                  if (marker(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID && marker(i1+1,i2+1,i3+1,0,nodeType,0)==FLUID)
                  {
                    SWAP(myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1),
                         myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1)); 

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
                      ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
                      ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
                      ForceX -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);
                      ForceY -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
                      ForceZ -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                
                    }                       
                  }
                  if((marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID && marker(i1-1,i2+1,i3+1,0,nodeType,0)==FLUID)) 
                  {
                    SWAP(myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1),
                         myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1));         

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);
                      ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
                      ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);  
                      ForceX += myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);
                      ForceY -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
                      ForceZ -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);  
                    }                        
                  }             
                  if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID && marker(i1-1,i2-1,i3+1,0,nodeType,0)==FLUID)
                  {
                    SWAP(myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1),
                         myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1));       

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);
                      ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
                      ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                 
                      ForceX += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);
                      ForceY += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
                      ForceZ -= myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1); 
                    }                       
                  }
                  if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID && marker(i1+1,i2-1,i3+1,0,nodeType,0)==FLUID)
                  {
                    SWAP(myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1),
                         myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1));       

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);
                      ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
                      ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1); 
                      ForceX -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);
                      ForceY += myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
                      ForceZ -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1); 
                    }                        
                  }
                  ///////////////////                                                                                   
                  ////  Group 10 ////                                                                                   
                  ///////////////////                             
                  // G10_DV_M1_M1_M1     // G10_DV_P1_M1_M1  // G10_DV_P1_P1_M1 // G10_DV_M1_P1_M1
                  
                  //MMM
                  if(marker(i1-1,i2-1,i3-1,0,nodeType,0)!=FLUID)  
                  {
                    SWAP(myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1),
                         myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1));        

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);
                      ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
                      ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
                    }                       
                  }                       
                  
                  //PMM
                  if(marker(i1+1,i2-1,i3-1,0,nodeType,0)!=FLUID)  
                  {
                    SWAP(myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1),
                         myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1));         

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);
                      ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
                      ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
                    }                       
                  }                        
                  
                  //PPM
                  if(marker(i1+1,i2+1,i3-1,0,nodeType,0)!=FLUID)  
                  {
                    SWAP(myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1),
                         myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1));       

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);
                      ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
                      ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
                    }                       
                  }                         
                  
                  //MPM
                  if(marker(i1-1,i2+1,i3-1,0,nodeType,0)!=FLUID)  
                  {
                    SWAP(myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1),
                         myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1));       

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);
                      ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
                      ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
                    }                      
                  }                          
          }
        }
      }
    }
    //////////
    // Cell //
    //////////     
    nodeType=nodeTYPE::CELL;
    for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
    {                                                                                                                                         
      for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
      {                                                                                                                                     
        for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                             
        {                                                                                                                                       
          if(marker(i1  ,i2  ,i3  ,0,nodeType,0)==FLUID )
          {            
            //////////////////                                                                                   
            ////  Group 9 ////                                                                                   
            //////////////////                             
            // G9_DV_P1_P1_P1      // G9_DV_M1_P1_P1   // G9_DV_M1_M1_P1 // G9_DV_P1_M1_P1
            // PPP
            if(marker(i1+1,i2+1,i3+1,0,nodeType,0)!=FLUID) 
            {
              SWAP(myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1)); 

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
                ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
              }                     
            }
            //MPP
            if(marker(i1-1,i2+1,i3+1,0,nodeType,0)!=FLUID)   
            {
              SWAP(myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1));         

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);
                ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
              }                       
            }
            //MMP
            if(marker(i1-1,i2-1,i3+1,0,nodeType,0)!=FLUID)   
            {
              SWAP(myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1));       

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);
                ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
              }                       
            }
            //PMP
            if(marker(i1+1,i2-1,i3+1,0,nodeType,0)!=FLUID)  
            {
              SWAP(myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1));       

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);
                ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
              }                      
            }  
            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID && marker(i1+1,i2+1,i3+1,0,nodeType,0)==FLUID)
            {
              SWAP(myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1));         

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
                ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
                ForceX -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);
                ForceY -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
                ForceZ -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                
              }                       
            }
            if( (marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID && marker(i1-1,i2+1,i3+1,0,nodeType,0)==FLUID))
            {
              SWAP(myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1));         

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);
                ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);  
                ForceX += myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);
                ForceY -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
                ForceZ -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);  
              }                 
            }
            if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID && marker(i1-1,i2-1,i3+1,0,nodeType,0)==FLUID )
            {
              SWAP(myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1));       

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);
                ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                 
                ForceX += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);
                ForceY += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
                ForceZ -= myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1); 
              }                     
            }
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID && marker(i1+1,i2-1,i3+1,0,nodeType,0)==FLUID)
            {
              SWAP(myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1));       

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);
                ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
                ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1); 
                ForceX -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);
                ForceY += myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
                ForceZ -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1); 
              }                      
            }
            ///////////////////                                                                                   
            ////  Group 10 ////                                                                                   
            ///////////////////                             
            // G10_DV_M1_M1_M1     // G10_DV_P1_M1_M1  // G10_DV_P1_P1_M1 // G10_DV_M1_P1_M1
            //MMM
            if(marker(i1-1,i2-1,i3-1,0,nodeType,0)!=FLUID)  
            {
              SWAP(myGrid(i1-1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1));        

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);
                ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
                ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
              }                     
            }                       
            //PMM
            if(marker(i1+1,i2-1,i3-1,0,nodeType,0)!=FLUID) 
            {
              SWAP(myGrid(i1+1,i2-1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1));         

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);
                ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
                ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
              }                      
            }                        
            //PPM
            if(marker(i1+1,i2+1,i3-1,0,nodeType,0)!=FLUID)
            {
              SWAP(myGrid(i1+1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1));       

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);
                ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
                ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
              }                       
            }                         
            //MPM
            if(marker(i1-1,i2+1,i3-1,0,nodeType,0)!=FLUID)  
            {
              SWAP(myGrid(i1-1,i2+1,i3-1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1));       

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);
                ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
                ForceZ -= myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
              }                        
            }                          
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>                                                                                                          
  void bounceBackWithForces(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ,int step)                      
  {
    ForceX= 0.0;
    ForceY= 0.0;
    ForceZ= 0.0;
    
    bounceBackG12Solid(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ);   
    bounceBackG34Solid(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ); 
    bounceBackG56Solid(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ); 
    bounceBackG78Solid(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ);
    bounceBackG910Solid(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ);   
    /*        bounceBackG12(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ);   
     *        bounceBackG34(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ); 
     *        bounceBackG56(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ); 
     *        bounceBackG78(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ);
     *        bounceBackG910(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ); */         
    
  }

  template <int N,int numblock, typename dataType1>
  void bounceBackG12Solid_beforeAdvection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)        
  {                                                                 
    for(int nodeType=0; nodeType<2;nodeType++)
    {
      for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
      {                                                                                                                                         
        for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
        {                                                                                                                                     
          for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                            
          {
            if(marker(i1,i2,i3,0,nodeType,0) == FLUID)
            {
              
              //////////////////
              ////  Group 1 ////
              //////////////////
              //ZZP
              if( marker(i1,i2,i3+1,0,nodeType,0) != FLUID )
              {
                SWAP(myGrid(i1  ,i2  ,i3    ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P1),
                     myGrid(i1  ,i2  ,i3+1  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceZ += myGrid(i1  ,i2  ,i3+1  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M1);
              }
              
              // ZZP2
              if( marker(i1,i2,i3+2,0,nodeType,0) != FLUID) 
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2),
                     myGrid(i1  ,i2  ,i3+2,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2));
               
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceZ += 2.0*myGrid(i1  ,i2  ,i3+2,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2);
              }                            
              
//               if( (marker(i1,i2,i3+1,0,nodeType,0) != FLUID && marker(i1,i2,i3+2,0,nodeType,0) == FLUID) )
//               {
//                 SWAP(myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2),
//                      myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2));
//                
//                 if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
//                 {
//                   ForceZ += 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2);
//                   ForceZ -= 2.0*myGrid(i1  ,i2  ,i3+2,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2);
//                 }
//               }  
              
              // ZP2Z
              if( (marker(i1,i2+2,i3,0,nodeType,0)!=FLUID) )
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO),
                     myGrid(i1  ,i2+2,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceY  += 2.0*myGrid(i1  ,i2+2,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO);
              }                            
              
//               if( (marker(i1,i2+1,i3,0,nodeType,0)!=FLUID) && (marker(i1,i2+2,i3,0,nodeType,0)==FLUID))
//               {
//                 SWAP(myGrid(i1  ,i2+2,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO),
//                      myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO));
//                 
//                 if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
//                 {
//                   ForceY  += 2.0*myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO);
//                   ForceY  -= 2.0*myGrid(i1  ,i2+2,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO);
//                 }  
//               }    
              
              // P2ZZ
              if( (marker(i1+2,i2,i3,0,nodeType,0)!=FLUID) )
              {
                SWAP(myGrid(i1  ,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO),
                     myGrid(i1+2,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceX += 2.0*myGrid(i1+2,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO);
              }                            
              
//               if( (marker(i1+1,i2,i3,0,nodeType,0)!=FLUID) &&  (marker(i1+2,i2,i3,0,nodeType,0)==FLUID) )
//               {
//                 SWAP(myGrid(i1+2,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO),
//                      myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO));
//                 
//                 if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
//                 {
//                   ForceX += 2.0*myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO);
//                   ForceX -= 2.0*myGrid(i1+2,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO);
//                 }
//               }                                
              
              //////////////////
              ////  Group 2 ////
              //////////////////                            
              //                             ZZM
              if( marker(i1,i2,i3-1,0,nodeType,0)!=FLUID )
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M1),
                     myGrid(i1  ,i2  ,i3-1,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceZ -= myGrid(i1  ,i2  ,i3-1,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P1);
              }                            
              
              // ZZM2
              if(marker(i1,i2,i3-2,0,nodeType,0)!=FLUID) 
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_ZERO_M2),
                     myGrid(i1  ,i2  ,i3-2,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceZ -= 2.0*myGrid(i1  ,i2  ,i3-2,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_ZERO_P2);
              }                            
              
              // ZM2Z
              if((marker(i1,i2-2,i3,0,nodeType,0)!=FLUID) )
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_ZERO_M2_ZERO),
                     myGrid(i1  ,i2-2,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceY -= 2.0*myGrid(i1  ,i2-2,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_ZERO_P2_ZERO);
              }                            
              
              // M2ZZ
              if((marker(i1-2,i2,i3,0,nodeType,0)!=FLUID) )
              {
                SWAP(myGrid(i1  ,i2 ,i3  ,lbModel.G2,nodeType,lbModel.G2_DV_M2_ZERO_ZERO),
                     myGrid(i1-2,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceX -= 2.0*myGrid(i1-2,i2 ,i3  ,lbModel.G1,nodeType,lbModel.G1_DV_P2_ZERO_ZERO);
              } 
            }
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>
  void bounceBackG34Solid_beforeAdvection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                     
  {                                                                                                                                                       
    
    for(int nodeType =0; nodeType<2;nodeType++)                                                                                                         
    {                                           
      //       int nodeType(0);
      for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
      {                                                                                                                                         
        for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
        {                                                                                                                                     
          for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                             
          {                                                                                                                       
            if(marker(i1,i2,i3,0,nodeType,0) == FLUID)
            {
              //////////////////                                                                            
              ////  Group 3 ////                                                                            
              //////////////////                                                                           
              //  G3_DV_ZERO_P1_P1 //  G3_DV_ZERO_M1_P1    //  G3_DV_P1_ZERO_P1   //  G3_DV_M1_ZERO_P1
              if(marker(i1  ,i2+1,i3+1,0,nodeType,0) != FLUID)
              {
                SWAP( myGrid(i1 ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1) ,
                      myGrid(i1 ,i2+1,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1)  );    
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceY += myGrid(i1 ,i2+1,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1);
                  ForceZ += myGrid(i1 ,i2+1,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1);
                }
              }
              
              if(marker(i1  ,i2-1,i3+1,0,nodeType,0)!=FLUID)
              {
                SWAP( myGrid(i1 ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1),
                      myGrid(i1 ,i2-1,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1));    
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceY -= myGrid(i1 ,i2-1,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1);
                  ForceZ += myGrid(i1 ,i2-1,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1);
                }
              }
              
              if(marker(i1+1,i2  ,i3+1,0,nodeType,0)!=FLUID)
              {
                SWAP( myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1),
                      myGrid(i1+1,i2  ,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX += myGrid(i1+1,i2  ,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1);
                  ForceZ += myGrid(i1+1,i2  ,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1);
                }
              }
              
              if(marker(i1-1,i2  ,i3+1,0,nodeType,0)!=FLUID)
              {
                SWAP( myGrid(i1  ,i2  ,i3  ,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1),
                      myGrid(i1-1,i2  ,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX -= myGrid(i1-1,i2  ,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1);
                  ForceZ += myGrid(i1-1,i2  ,i3+1,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1);
                }
              }
              
              //////////////////                                                                                                     
              ////  Group 4 ////                                                                                                     
              //////////////////                                                                                                      
              //  G4_DV_ZERO_M1_M1  //  G4_DV_ZERO_P1_M1  //  G4_DV_M1_ZERO_M1    //  G4_DV_P1_ZERO_M1                            
              if(marker(i1  ,i2-1,i3-1,0,nodeType,0)!=FLUID)                                                                                    
              {                                                                                                                             
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_M1_M1),                                                 
                     myGrid(i1  ,i2-1,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1));                                                
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceY -= myGrid(i1  ,i2-1,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1);
                  ForceZ -= myGrid(i1  ,i2-1,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_P1_P1);
                }                                                                                                                             
              }                            
              if(marker(i1  ,i2+1,i3-1,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_ZERO_P1_M1),
                     myGrid(i1  ,i2+1,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceY += myGrid(i1  ,i2+1,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1);
                  ForceZ -= myGrid(i1  ,i2+1,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_ZERO_M1_P1);
                }
              }                            
              if(marker(i1-1,i2  ,i3-1,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_M1_ZERO_M1),
                     myGrid(i1-1,i2  ,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX -= myGrid(i1-1,i2  ,i3-1 ,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1);
                  ForceZ -= myGrid(i1-1,i2  ,i3-1 ,lbModel.G3,nodeType,lbModel.G3_DV_P1_ZERO_P1);
                }
              }                            
              
              if(marker(i1+1,i2  ,i3-1,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G4,nodeType,lbModel.G4_DV_P1_ZERO_M1),
                     myGrid(i1+1,i2  ,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX += myGrid(i1+1,i2  ,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1);
                  ForceZ -= myGrid(i1+1,i2  ,i3-1,lbModel.G3,nodeType,lbModel.G3_DV_M1_ZERO_P1);
                }
              }
            }                                                            
          }
        }
      }
    }
  }  
  
  template <int N,int numblock, typename dataType1>
  void bounceBackG56Solid_beforeAdvection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)               
  {                                                                                                                                                
    for(int nodeType =0; nodeType<2;nodeType++)                                                                                                  
    {       
      //    int nodeType =0;
      for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
      {                                                                                                                                         
        for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
        {                                                                                                                                     
          for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                           
          {
            if(marker(i1,i2,i3,0,nodeType,0)==FLUID)
            {                                                                                               
              //////////////////                                                                           
              ////  Group 5 ////                                                                          
              //////////////////                                                                           
              //G5_DV_P1_P1_ZERO   //G5_DV_M1_P1_ZERO   //G5_DV_ZERO_P1_ZERO   //G5_DV_P1_ZERO_ZERO
              if(marker(i1+1,i2+1,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO),
                     myGrid(i1+1,i2+1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX += myGrid(i1+1,i2+1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO);
                  ForceY += myGrid(i1+1,i2+1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO);
                }
              }
              
              if(marker(i1-1,i2+1,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO),
                     myGrid(i1-1,i2+1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX -= myGrid(i1-1,i2+1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO);
                  ForceY += myGrid(i1-1,i2+1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO);
                }
              }
              
              if(marker(i1  ,i2+1,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_ZERO_P1_ZERO),
                     myGrid(i1  ,i2+1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_ZERO_M1_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceY += myGrid(i1  ,i2+1,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_ZERO_M1_ZERO);
                }
              }
              
              if(marker(i1+1,i2  ,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_ZERO_ZERO),
                     myGrid(i1+1,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_ZERO_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX += myGrid(i1+1,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_ZERO_ZERO);
                } 
              }
              
              //////////////////                                                                           
              ////  Group 6 ////                                                                           
              //////////////////                                                                           
              //G6_DV_M1_M1_ZERO   //G6_DV_P1_M1_ZERO    //G6_DV_ZERO_M1_ZERO   //G6_DV_M1_ZERO_ZERO
              if(marker(i1-1,i2-1,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_M1_ZERO),
                     myGrid(i1-1,i2-1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO));
                
                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX -= myGrid(i1-1,i2-1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO);
                  ForceY -= myGrid(i1-1,i2-1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_P1_ZERO); 
                }
              }
              
              if(marker(i1+1,i2-1,i3  ,0,nodeType,0) !=FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_P1_M1_ZERO),
                     myGrid(i1+1,i2-1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO));

                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX += myGrid(i1+1,i2-1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO);
                  ForceY -= myGrid(i1+1,i2-1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_M1_P1_ZERO); 
                }
              }
              
              if(marker(i1  ,i2-1,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_ZERO_M1_ZERO),
                     myGrid(i1  ,i2-1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_ZERO_P1_ZERO));

                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceY -= myGrid(i1  ,i2-1,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_ZERO_P1_ZERO);
              }
              
              if(marker(i1-1,i2  ,i3  ,0,nodeType,0)!=FLUID)
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G6,nodeType,lbModel.G6_DV_M1_ZERO_ZERO),
                     myGrid(i1-1,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_ZERO_ZERO));

                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  ForceX -= myGrid(i1-1,i2  ,i3  ,lbModel.G5,nodeType,lbModel.G5_DV_P1_ZERO_ZERO);
              }
            }                                                            
          }
        }
      }
    }
  }  
    
  template <int N,int numblock, typename dataType1>                                                                                                          
  void bounceBackG910Solid_beforeAdvection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                       
  {                                                                                                                                                          
    //////////
    // Node //
    //////////     
    int nodeType=nodeTYPE::NODE;
    for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
    {                                                                                                                                         
      for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
      {                                                                                                                                     
        for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                         
        {                                                                                                                                       
          if(marker(i1  ,i2  ,i3  ,0,nodeType,0)==FLUID )
          {            
            //////////////////                                                                                   
            ////  Group 9 ////                                                                                   
            //////////////////                             
            // G9_DV_P1_P1_P1      // G9_DV_M1_P1_P1   // G9_DV_M1_M1_P1 // G9_DV_P1_M1_P1
            
            // PPP
            if(marker(i1+1,i2+1,i3+1,0,nodeType,0)!=FLUID)  
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1),
                   myGrid(i1+1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1)); 

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += myGrid(i1+1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
                ForceY += myGrid(i1+1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
                ForceZ += myGrid(i1+1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
              }
              
            }
              
              //MPP
              if(marker(i1-1,i2+1,i3+1,0,nodeType,0)!=FLUID) 
              {
                SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1),
                     myGrid(i1-1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1));         

                if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                {
                  ForceX -= myGrid(i1-1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);
                  ForceY += myGrid(i1-1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
                  ForceZ += myGrid(i1-1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
                }   
                
              }
                //MMP
                if(marker(i1-1,i2-1,i3+1,0,nodeType,0)!=FLUID)   
                {
                  SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1),
                       myGrid(i1-1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1));       

                  if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                  {
                    ForceX -= myGrid(i1-1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);
                    ForceY -= myGrid(i1-1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
                    ForceZ += myGrid(i1-1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
                  }
                  
                }
                  //PMP
                  if(marker(i1+1,i2-1,i3+1,0,nodeType,0)!=FLUID)  
                  {
                    SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1),
                         myGrid(i1+1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1));       

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX += myGrid(i1+1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);
                      ForceY -= myGrid(i1+1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
                      ForceZ += myGrid(i1+1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
                    }  
                  }                        
//                   if (marker(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID && marker(i1+1,i2+1,i3+1,0,nodeType,0)==FLUID)
//                   {
//                     SWAP(myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1),
//                          myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1)); 
// 
//                     if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
//                     {
//                       ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
//                       ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
//                       ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
//                       ForceX -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);
//                       ForceY -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
//                       ForceZ -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                
//                     }                       
//                   }
//                   if((marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID && marker(i1-1,i2+1,i3+1,0,nodeType,0)==FLUID)) 
//                   {
//                     SWAP(myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1),
//                          myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1));         
// 
//                     if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
//                     {
//                       ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);
//                       ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
//                       ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);  
//                       ForceX += myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);
//                       ForceY -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
//                       ForceZ -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);  
//                     }                        
//                   }             
//                   if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID && marker(i1-1,i2-1,i3+1,0,nodeType,0)==FLUID)
//                   {
//                     SWAP(myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1),
//                          myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1));       
// 
//                     if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
//                     {
//                       ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);
//                       ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
//                       ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                 
//                       ForceX += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);
//                       ForceY += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
//                       ForceZ -= myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1); 
//                     }                       
//                   }
//                   if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID && marker(i1+1,i2-1,i3+1,0,nodeType,0)==FLUID)
//                   {
//                     SWAP(myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1),
//                          myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1));       
// 
//                     if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
//                     {
//                       ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);
//                       ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
//                       ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1); 
//                       ForceX -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);
//                       ForceY += myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
//                       ForceZ -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1); 
//                     }                        
//                   }
                  ///////////////////                                                                                   
                  ////  Group 10 ////                                                                                   
                  ///////////////////                             
                  // G10_DV_M1_M1_M1     // G10_DV_P1_M1_M1  // G10_DV_P1_P1_M1 // G10_DV_M1_P1_M1
                  
                  //MMM
                  if(marker(i1-1,i2-1,i3-1,0,nodeType,0)!=FLUID)  
                  {
                    SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1),
                         myGrid(i1-1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1));        

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX -= myGrid(i1-1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);
                      ForceY -= myGrid(i1-1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
                      ForceZ -= myGrid(i1-1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
                    }                       
                  }                       
                  
                  //PMM
                  if(marker(i1+1,i2-1,i3-1,0,nodeType,0)!=FLUID)  
                  {
                    SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1),
                         myGrid(i1+1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1));         

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX += myGrid(i1+1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);
                      ForceY -= myGrid(i1+1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
                      ForceZ -= myGrid(i1+1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
                    }                       
                  }                        
                  
                  //PPM
                  if(marker(i1+1,i2+1,i3-1,0,nodeType,0)!=FLUID)  
                  {
                    SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1),
                         myGrid(i1+1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1));       

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX += myGrid(i1+1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);
                      ForceY += myGrid(i1+1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
                      ForceZ -= myGrid(i1+1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
                    }                       
                  }                         
                  
                  //MPM
                  if(marker(i1-1,i2+1,i3-1,0,nodeType,0)!=FLUID)  
                  {
                    SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1),
                         myGrid(i1-1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1));       

                    if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
                    {
                      ForceX -= myGrid(i1-1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);
                      ForceY += myGrid(i1-1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
                      ForceZ -= myGrid(i1-1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
                    }                      
                  }                          
          }
        }
      }
    }
    //////////
    // Cell //
    //////////     
    nodeType=nodeTYPE::CELL;
    for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
    {                                                                                                                                         
      for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
      {                                                                                                                                     
        for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                             
        {                                                                                                                                       
          if(marker(i1  ,i2  ,i3  ,0,nodeType,0)==FLUID )
          {            
            //////////////////                                                                                   
            ////  Group 9 ////                                                                                   
            //////////////////                             
            // G9_DV_P1_P1_P1      // G9_DV_M1_P1_P1   // G9_DV_M1_M1_P1 // G9_DV_P1_M1_P1
            // PPP
            if(marker(i1+1,i2+1,i3+1,0,nodeType,0)!=FLUID) 
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1),
                   myGrid(i1+1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1)); 

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += myGrid(i1+1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
                ForceY += myGrid(i1+1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
                ForceZ += myGrid(i1+1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
              }                     
            }
            //MPP
            if(marker(i1-1,i2+1,i3+1,0,nodeType,0)!=FLUID)   
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1),
                   myGrid(i1-1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1));         

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= myGrid(i1-1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);
                ForceY += myGrid(i1-1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
                ForceZ += myGrid(i1-1,i2+1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
              }                       
            }
            //MMP
            if(marker(i1-1,i2-1,i3+1,0,nodeType,0)!=FLUID)   
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1),
                   myGrid(i1-1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1));       

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= myGrid(i1-1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);
                ForceY -= myGrid(i1-1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
                ForceZ += myGrid(i1-1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
              }                       
            }
            //PMP
            if(marker(i1+1,i2-1,i3+1,0,nodeType,0)!=FLUID)  
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1),
                   myGrid(i1+1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1));       

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += myGrid(i1+1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);
                ForceY -= myGrid(i1+1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
                ForceZ += myGrid(i1+1,i2-1,i3+1,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
              }                      
            }  
//             if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID && marker(i1+1,i2+1,i3+1,0,nodeType,0)==FLUID)
//             {
//               SWAP(myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1),
//                    myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1));         
// 
//               if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
//               {
//                 ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
//                 ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);                       
//                 ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1);
//                 ForceX -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);
//                 ForceY -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
//                 ForceZ -= myGrid(i1+1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                
//               }                       
//             }
//             if( (marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID && marker(i1-1,i2+1,i3+1,0,nodeType,0)==FLUID))
//             {
//               SWAP(myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1),
//                    myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1));         
// 
//               if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
//               {
//                 ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);
//                 ForceY += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);                       
//                 ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1);  
//                 ForceX += myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);
//                 ForceY -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
//                 ForceZ -= myGrid(i1-1,i2+1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);  
//               }                 
//             }
//             if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID && marker(i1-1,i2-1,i3+1,0,nodeType,0)==FLUID )
//             {
//               SWAP(myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1),
//                    myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1));       
// 
//               if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
//               {
//                 ForceX -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);
//                 ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                       
//                 ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1);                 
//                 ForceX += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);
//                 ForceY += myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
//                 ForceZ -= myGrid(i1-1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1); 
//               }                     
//             }
//             if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID && marker(i1+1,i2-1,i3+1,0,nodeType,0)==FLUID)
//             {
//               SWAP(myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1),
//                    myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1));       
// 
//               if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
//               {
//                 ForceX += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);
//                 ForceY -= myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1);                       
//                 ForceZ += myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1); 
//                 ForceX -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);
//                 ForceY += myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
//                 ForceZ -= myGrid(i1+1,i2-1,i3+1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1); 
//               }                      
//             }
            ///////////////////                                                                                   
            ////  Group 10 ////                                                                                   
            ///////////////////                             
            // G10_DV_M1_M1_M1     // G10_DV_P1_M1_M1  // G10_DV_P1_P1_M1 // G10_DV_M1_P1_M1
            //MMM
            if(marker(i1-1,i2-1,i3-1,0,nodeType,0)!=FLUID)  
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_M1_M1),
                   myGrid(i1-1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1));        

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= myGrid(i1-1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);
                ForceY -= myGrid(i1-1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
                ForceZ -= myGrid(i1-1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_P1_P1);                       
              }                     
            }                       
            //PMM
            if(marker(i1+1,i2-1,i3-1,0,nodeType,0)!=FLUID) 
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_M1_M1),
                   myGrid(i1+1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1));         

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += myGrid(i1+1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);
                ForceY -= myGrid(i1+1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
                ForceZ -= myGrid(i1+1,i2-1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_P1_P1);                       
              }                      
            }                        
            //PPM
            if(marker(i1+1,i2+1,i3-1,0,nodeType,0)!=FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_P1_P1_M1),
                   myGrid(i1+1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1));       

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += myGrid(i1+1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);
                ForceY += myGrid(i1+1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
                ForceZ -= myGrid(i1+1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_M1_M1_P1);                       
              }                       
            }                         
            //MPM
            if(marker(i1-1,i2+1,i3-1,0,nodeType,0)!=FLUID)  
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G10,nodeType,lbModel.G10_DV_M1_P1_M1),
                   myGrid(i1-1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1));       

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= myGrid(i1-1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);
                ForceY += myGrid(i1-1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
                ForceZ -= myGrid(i1-1,i2+1,i3-1,lbModel.G9 ,nodeType, lbModel.G9_DV_P1_M1_P1);                       
              }                        
            }                          
          }
        }
      }
    }
  }
  
  template <int N,int numblock, typename dataType1>                                                                                                          
  void bounceBackG78Solid_beforeAdvection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ)                        
  {    
    ///////////////////
    // Node --> Cell //
    ///////////////////
    for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
    {                                                                                                                                         
      for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
      {                                                                                                                                     
        for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                             
        {                                                                                                                                       
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0)==FLUID)
          {                                                      
            //////////////////                                                                         
            ////  Group 7 ////                                                                         
            //////////////////  
            // G7_DV_P_P_P    // G7_DV_M_P_P   // G7_DV_M_M_P   // G7_DV_P_M_P 
            
            // PPP
            if(marker(i1  ,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);
                ForceY += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);                       
                ForceZ += 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M); 
              }                            
            }                                                         
            
            // MPP
            if(marker(i1-1,i2  ,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P),
                   myGrid(i1-1,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1-1,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);
                ForceY += 0.5*myGrid(i1-1,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M);                       
                ForceZ += 0.5*myGrid(i1-1,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M); 
              }                            
            }
            
            // MMP
            if(marker(i1-1,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P),
                   myGrid(i1-1,i2-1,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1-1,i2-1,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);
                ForceY -= 0.5*myGrid(i1-1,i2-1,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M);                       
                ForceZ += 0.5*myGrid(i1-1,i2-1,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M); 
              }                            
            }
            
            // PMP
            if(marker(i1  ,i2-1,i3  ,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P),
                   myGrid(i1  ,i2-1,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1  ,i2-1,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);
                ForceY -= 0.5*myGrid(i1  ,i2-1,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M);                       
                ForceZ += 0.5*myGrid(i1  ,i2-1,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M); 
              }                            
            }
            
            //////////////////                                                                                   
            ////  Group 8 ////                                                                                   
            //////////////////                                                                                   
            
            // Node --> Cell      // G8_DV_M_M_M  // G8_DV_P_M_M    // G8_DV_P_P_M // G8_DV_M_P_M         
            
            // MMM
            if(marker(i1-1,i2-1,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M),
                   myGrid(i1-1,i2-1,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                
                ForceX -= 0.5*myGrid(i1-1,i2-1,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);
                ForceY -= 0.5*myGrid(i1-1,i2-1,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);                       
                ForceZ -= 0.5*myGrid(i1-1,i2-1,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P); 
              }                  
            }                        
            // PMM
            if(marker(i1  ,i2-1,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M),
                   myGrid(i1  ,i2-1,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1  ,i2-1,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);
                ForceY -= 0.5*myGrid(i1  ,i2-1,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);                       
                ForceZ -= 0.5*myGrid(i1  ,i2-1,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P);                             
              }                            
            }                  
            
            // PPM
            if(marker(i1  ,i2  ,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M),
                   myGrid(i1  ,i2  ,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1  ,i2  ,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);
                ForceY += 0.5*myGrid(i1  ,i2  ,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);                       
                ForceZ -= 0.5*myGrid(i1  ,i2  ,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P);                         
              }                  
            }                        
            // MPM
            if(marker(i1-1,i2  ,i3-1,0,nodeTYPE::CELL,0)!=FLUID)
            {
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M),
                   myGrid(i1-1,i2  ,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1-1,i2  ,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);
                ForceY += 0.5*myGrid(i1-1,i2  ,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);                       
                ForceZ -= 0.5*myGrid(i1-1,i2  ,i3-1,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P);                       
              }                      
            }                  
          }
        }
      }
    }
    
    ////////////////////
    // Cell --> Node  //
    ////////////////////
    for (int i3 = myGrid.ndB3; i3 <= myGrid.ndE3; i3++)                                                                                         
    {                                                                                                                                         
      for (int i2 = myGrid.ndB2; i2 <= myGrid.ndE2; i2++)                                                                                         
      {                                                                                                                                     
        for (int i1 = myGrid.ndB1; i1<= myGrid.ndE1; i1++)                                                                                            
        {                                                                                                                                        
          if(marker(i1,i2,i3,0,nodeTYPE::CELL,0)==FLUID)
          {                                                                        
            //////////////////                                                                                                 
            ////  Group 7 ////                                                                                                  
            //////////////////  
            // G7_DV_P_P_P           // G7_DV_P_M_P         // G7_DV_M_P_P         // G7_DV_M_M_P
            
            // PPP
            if(marker(i1+1,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P),
                   myGrid(i1+1,i2+1,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1+1,i2+1,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);
                ForceY += 0.5*myGrid(i1+1,i2+1,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);                       
                ForceZ += 0.5*myGrid(i1+1,i2+1,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);                       
              }                  
            }                        // MPP
            if(marker(i1  ,i2+1,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P),
                   myGrid(i1  ,i2+1,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1  ,i2+1,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);
                ForceY += 0.5*myGrid(i1  ,i2+1,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);                       
                ForceZ += 0.5*myGrid(i1  ,i2+1,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M);                       
              }                       
              
            }       
            // MMP
            if(marker(i1  ,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P),
                   myGrid(i1  ,i2  ,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1  ,i2  ,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);
                ForceY -= 0.5*myGrid(i1  ,i2  ,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);                       
                ForceZ += 0.5*myGrid(i1  ,i2  ,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M);                       
              }                      
              
            }        
            // PMP
            if(marker(i1+1,i2  ,i3+1,0,nodeTYPE::NODE,0)!=FLUID)                                                                                    
            {                                                                                                                                   
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P),                                             
                   myGrid(i1+1,i2  ,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M)); 

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1+1,i2  ,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);
                ForceY -= 0.5*myGrid(i1+1,i2  ,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);                       
                ForceZ += 0.5*myGrid(i1+1,i2  ,i3+1,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M);                       
              }                      
              
            }                              
            
            //////////////////                                                                                                 
            ////  Group 8 ////                                                                                                  
            //////////////////   
            // G8_DV_M_M_M           // G8_DV_M_P_M         // G8_DV_P_M_M         // G8_DV_P_P_M 
            
            // MMM
            if(marker(i1  ,i2   ,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M),
                   myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);
                ForceY -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);                       
                ForceZ -= 0.5*myGrid(i1  ,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);                       
              }                  
            }                  
            // PMM
            if(marker(i1+1,i2  ,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M),
                   myGrid(i1+1,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1+1,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);
                ForceY -= 0.5*myGrid(i1+1,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);                       
                ForceZ -= 0.5*myGrid(i1+1,i2  ,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P);                       
              }                            
            }       
            // PPM
            if(marker(i1+1,i2+1,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M),
                   myGrid(i1+1,i2+1,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX += 0.5*myGrid(i1+1,i2+1,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);
                ForceY += 0.5*myGrid(i1+1,i2+1,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);                       
                ForceZ -= 0.5*myGrid(i1+1,i2+1,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P);                       
              }                       
            }        
            // MPM
            if(marker(i1  ,i2+1,i3  ,0,nodeTYPE::NODE,0)!=FLUID)
            { 
              SWAP(myGrid(i1  ,i2  ,i3  ,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M),
                   myGrid(i1  ,i2+1,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P));

              if(i1 >= myGrid.nB1 && i2 >= myGrid.nB2 && i3 >= myGrid.nB3 && i1 <= myGrid.nE1 && i2 <= myGrid.nE2 && i3 <= myGrid.nE3)
              {
                ForceX -= 0.5*myGrid(i1  ,i2+1,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);
                ForceY += 0.5*myGrid(i1  ,i2+1,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);                       
                ForceZ -= 0.5*myGrid(i1  ,i2+1,i3  ,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P);                       
              }                      
            }                           
          }
        }
      }
    }
  }           
   
  template <int N,int numblock, typename dataType1>                                                                                                          
  void bounceBackWithForces_beforeAdvection(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int FLUID,dataType1 &ForceX, dataType1 &ForceY, dataType1 &ForceZ,int step)                      
  {
    ForceX= 0.0;
    ForceY= 0.0;
    ForceZ= 0.0;
    
    bounceBackG12Solid_beforeAdvection(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ);   
    bounceBackG34Solid_beforeAdvection(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ); 
    bounceBackG56Solid_beforeAdvection(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ); 
    bounceBackG78Solid_beforeAdvection(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ);
    bounceBackG910Solid_beforeAdvection(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ);   
    /*        bounceBackG12(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ);   
     *        bounceBackG34(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ); 
     *        bounceBackG56(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ); 
     *        bounceBackG78(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ);
     *        bounceBackG910(lbModel,myGrid,marker,FLUID,ForceX,ForceY,ForceZ); */         
    
  }

  
// Explicit declarations
  template void SWAP<double>(double &,double &);
  template void bounceBackG12<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG34<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG56<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG78<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG910<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG12Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG34Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG56Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG78Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG910Solid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackWithForces<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &,int );      
  template void bounceBackG12Solid_beforeAdvection<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG34Solid_beforeAdvection<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG56Solid_beforeAdvection<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG78Solid_beforeAdvection<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackG910Solid_beforeAdvection<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &);
  template void bounceBackWithForces_beforeAdvection<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,gridBCC3D<1, 1, int> &,int ,double &, double &, double &,int );      
  
  
  

