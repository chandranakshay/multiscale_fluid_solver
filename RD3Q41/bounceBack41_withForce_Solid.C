  template <int N,int numblock, typename dataType1>
  inline void getRhoForDiffuseBounceBack(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,gridBCC3D<1, 1, dataType1> &rhoGrid ,int FLUID, int SOLID)            
  {
    dataType1 rho(0.0);
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) == FLUID)
          {
            getRhoNodeSinglePoint(lbModel,myGrid,rho,i1,i2,i3);
            rhoGrid(i1,i2,i3,0,nodeTYPE::NODE,0) = rho; 
          }
        }
      }
    }
    for (int i3 = myGrid.nE3; i3 >= myGrid.nB3; i3--)                   
    {                                                                
      for (int i2 = myGrid.nE2; i2 >= myGrid.nB2; i2--)                  
      {                    
        for (int i1 = myGrid.nE1; i1 >= myGrid.nB1; i1--)                 
        {
          if(marker(i1,i2,i3,0,nodeTYPE::NODE,0) == FLUID)
          {
            getRhoCellSinglePoint(lbModel,myGrid,rho,i1,i2,i3);
            rhoGrid(i1,i2,i3,0,nodeTYPE::CELL,0) = rho; 
          }
        }
      }
    }    
    
  }
