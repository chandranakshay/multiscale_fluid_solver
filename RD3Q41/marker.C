#include"marker.h"

void initializeMarker(gridBCC3D<1, 1, int> &marker, PANINI_INT FLUID)
{
    // Initisalise everything as Fluid --> Including Ghost Markers
    for(int k=0; k<marker.n3;k++)
        for(int j=0; j<marker.n2;j++)
            for(int i=0; i<marker.n1;i++)
            {
                marker(i,j,k,0,0,0) = FLUID;  
                marker(i,j,k,0,1,0) = FLUID;  
            }
}

void markTopAndBottomWalls      (gridBCC3D<1, 1, int> &marker, PANINI_INT FLUID, PANINI_INT SOLID, int *coord,int nP1,int nP2,int nP3)
{
    int topWall    = nP2*marker.m2 + marker.nB2 - 1  ;
    int bottomWall = marker.nB2;
    
    for(int k=0; k<marker.n3;k++)
        for(int j=0; j<marker.n2;j++)
            for(int i=0; i<marker.n1;i++)
            {
          int jGlobal = j + marker.m2*coord[1] ; 
                
          if( (jGlobal == bottomWall) || (jGlobal == topWall) )
            marker(i,j,k,0,0,0) = SOLID;
          if( (jGlobal == bottomWall) || (jGlobal == topWall) )
            marker(i,j,k,0,1,0) = SOLID;
            }    
}


void markOnePercentBumpOnTheWall(gridBCC3D<1, 1, int> &marker, PANINI_INT FLUID, PANINI_INT SOLID, int *coord,int nP1,int nP2,int nP3,int nX,int nY,int nZ)
{
    int bottomWall = marker.nB2;
    
    int XCENTER = (int)(0.5*(double)nP1*(double)nX) + marker.nB1 + 5;
    int ZCENTER = (int)(0.5*(double)nP3*(double)nZ) + marker.nB3 + 5;

    int onepercentBump = (int)(((marker.m2*nP2)-2.5)*0.01)+1; 
    
    for(int k=0; k<marker.n3;k++)
        for(int j=0; j<marker.n2;j++)
            for(int i=0; i<marker.n1;i++)
            {
             int iGlobal = i + marker.m1*coord[0] ; 
             int jGlobal = j + marker.m2*coord[1] ; 
             int kGlobal = k + marker.m3*coord[2] ; 
             
               if((jGlobal == bottomWall) && (iGlobal == XCENTER) && (kGlobal == ZCENTER ))
                {
                 std::cout<<"Bump size: "<<onepercentBump<<"\n the following points are marked \n"<<std::endl;
                   
                for(int index=1;index<=onepercentBump;index++)
                {    
                 marker(i,j+index,k,0,0,0) = SOLID;
                 marker(i,j+index,k,0,1,0) = SOLID;
                 std::cout<<iGlobal<<"  "<<jGlobal+index<<"  "<<kGlobal<<std::endl;
                }
               }
            } 
}   





