     
namespace panini { 
    
     template <int N,int numblock, typename dataType1>        
     inline void printVelocity(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int SOLID,int FLUID,int VECT_LENGTH, int step)        
     {        
      dataType1 rho[VECT_LENGTH],uX[VECT_LENGTH],uY[VECT_LENGTH],uZ[VECT_LENGTH],theta[VECT_LENGTH];
      std::ofstream file;        
      char fileName[150];        
      sprintf(fileName,"./results/velocity_%d.vtk",step);        
      file.open(fileName);        
     //   vtk file header    
      file<<"# vtk DataFile Version 3.0"<<std::endl<<"Velocity"<<std::endl<<"ASCII"<<std::endl<<"DATASET STRUCTURED_GRID"<<std::endl;    
      file<<"DIMENSIONS "<<myGrid.m1<<" "<<myGrid.m2<<" "<<myGrid.m3<<std::endl;
      file<<"POINTS "<<myGrid.actualGridSize<<" double"<<std::endl; 
        for(int k=myGrid.nB3; k<=myGrid.nE3;k++)
            for(int j=myGrid.nB2; j<=myGrid.nE2;j++)
                for(int i=myGrid.nB1; i<=myGrid.nE1;i=i+4)
                     for(int index=0;index<VECT_LENGTH;index++)
                      file<<i+index<<" "<<j<<" "<<k<<std::endl;
       
      file<<"POINT_DATA "<<myGrid.actualGridSize<<std::endl;
      file<<"VECTORS"<<" "<<"velocity"<<" "<<"double"<<std::endl;
        for(int k=myGrid.nB3; k<=myGrid.nE3;k++)
            for(int j=myGrid.nB2; j<=myGrid.nE2;j++)
                for(int i=myGrid.nB1; i<=myGrid.nE1;i=i+4)
                {
//                  copyFromNode(lbModel,myGrid,VECT_LENGTH,i,j,k);    
                 getHydroMoment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);                     
                 for(int index=0;index<VECT_LENGTH;index++)
                 {
                  if(marker(i+index,j,k,0,0,0) == FLUID)   
                   file << uX[index] << " " << uY[index] << " " << uZ[index] << std::endl;
                  else
                   file << "0.0  0.0  0.0" << std::endl;
                 }
               }
     } 
     
template <typename dataType1>  
inline void printMarkerNode(gridBCC3D<1,1,dataType1> &marker,int step)        
     {        
      std::ofstream file;        
      char fileName[150];        
      sprintf(fileName,"./results/markerNode_%d.vtk",step);        
      file.open(fileName);        
     //   vtk file header    
      file<<"# vtk DataFile Version 3.0"<<std::endl<<"Marker"<<std::endl<<"ASCII"<<std::endl<<"DATASET STRUCTURED_GRID"<<std::endl;    
      file<<"DIMENSIONS "<<marker.m1<<" "<<marker.m2<<" "<<marker.m3<<std::endl;
      file<<"POINTS "<<marker.actualGridSize<<" double"<<std::endl; 
        for(int k=marker.nB3; k<=marker.nE3;k++)
            for(int j=marker.nB2; j<=marker.nE2;j++)
                for(int i=marker.nB1; i<=marker.nE1;i++)
                      file<<i<<" "<<j<<" "<<k<<std::endl;
       
      file<<"POINT_DATA "<<marker.actualGridSize<<std::endl;
      file<<"SCALARS"<<" "<<"Marker"<<" "<<"double"<<std::endl;
      file<<"LOOKUP_TABLE default\n"<<std::endl;
        for(int k=marker.nB3; k<=marker.nE3;k++)
            for(int j=marker.nB2; j<=marker.nE2;j++)
                for(int i=marker.nB1; i<=marker.nE1;i++)
                     file <<marker(i,j,k,0,0,0) << std::endl;
     }      

template <typename dataType1>  
inline void printMarkerCell(gridBCC3D<1,1,dataType1> &marker,int step)        
     {        
      std::ofstream file;        
      char fileName[150];        
      sprintf(fileName,"./results/markerCell_%d.vtk",step);        
      file.open(fileName);        
     //   vtk file header    
      file<<"# vtk DataFile Version 3.0"<<std::endl<<"Marker"<<std::endl<<"ASCII"<<std::endl<<"DATASET STRUCTURED_GRID"<<std::endl;    
      file<<"DIMENSIONS "<<marker.m1<<" "<<marker.m2<<" "<<marker.m3<<std::endl;
      file<<"POINTS "<<marker.actualGridSize<<" double"<<std::endl; 
        for(int k=marker.nB3; k<=marker.nE3;k++)
            for(int j=marker.nB2; j<=marker.nE2;j++)
                for(int i=marker.nB1; i<=marker.nE1;i++)
                      file<<i<<" "<<j<<" "<<k<<std::endl;
       
      file<<"POINT_DATA "<<marker.actualGridSize<<std::endl;
      file<<"SCALARS"<<" "<<"Marker"<<" "<<"double"<<std::endl;
      file<<"LOOKUP_TABLE default\n"<<std::endl;
        for(int k=marker.nB3; k<=marker.nE3;k++)
            for(int j=marker.nB2; j<=marker.nE2;j++)
                for(int i=marker.nB1; i<=marker.nE1;i++)
                     file <<marker(i,j,k,0,1,0) << std::endl;
     }      
}
