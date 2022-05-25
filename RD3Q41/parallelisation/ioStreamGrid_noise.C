#include"ioStreamGrid_noise.h"
 
template <int N,int numblock, typename dataType1>
void probeComputationalGridPlane12Node_Noise(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i3, int k)
{  
 int i1, i2 ; 
 int index = beginIndex; 
 for(i2= 0; i2< myGrid.n2 ;  i2++)
 {
  for(i1= 0; i1< myGrid.n1;  i1++)
  {
   tmpData[index] = myGrid(i1,i2,i3,groupNumber,myGrid.node[groupNumber],k) ;
   index +=stride;  
  }
 } 
}

template <int N,int numblock, typename dataType1>
void probeComputationalGridPlane12Cell_Noise(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i3, int k)
{  
 int i1, i2 ; 
 int index = beginIndex; 
 for(i2= 0; i2< myGrid.n2;  i2++)
 {
  for(i1= 0; i1< myGrid.n1;  i1++)
  {
   tmpData[index] = myGrid(i1,i2,i3,groupNumber,myGrid.cell[groupNumber],k) ;
   index +=stride;  
  }
 } 
}

template <int N,int numblock, typename dataType1>
void probeComputationalGridPlane23Node_Noise(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i1, int k)
{  
 int i3, i2 ; 
 int index = beginIndex; 
 for(i3= 0  ; i3< myGrid.n3 ;  i3++)
 {
  for(i2= 0; i2< myGrid.n2 ;  i2++)
  {
   tmpData[index] = myGrid(i1,i2,i3,groupNumber,myGrid.node[groupNumber],k) ;
   index +=stride;  
  }
 } 
}


template <int N,int numblock, typename dataType1>
void probeComputationalGridPlane23Cell_Noise(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i1, int k)
{  
 int i3, i2 ; 
 int index = beginIndex; 
 for(i3= 0  ; i3< myGrid.n3  ;  i3++)
 {
  for(i2= 0; i2< myGrid.n2;  i2++)
  {
   tmpData[index] = myGrid(i1,i2,i3,groupNumber,myGrid.cell[groupNumber],k) ;
   index +=stride;  
  }
 } 
}

template <int N,int numblock, typename dataType1>
void probeComputationalGridPlane31Node_Noise(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i2, int k)
{  
 int i1, i3 ; 
 int index = beginIndex; 
 for(i3= 0; i3< myGrid.n3 ;  i3++)
 {
  for(i1= 0; i1< myGrid.n1 ;  i1++)
  {
   tmpData[index] = myGrid(i1,i2,i3,groupNumber,myGrid.node[groupNumber],k) ;
   index +=stride;  
  }
 } 
}

template <int N,int numblock, typename dataType1>
void probeComputationalGridPlane31Cell_Noise(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i2, int k)
{  
 int i1, i3 ; 
 int index = beginIndex; 
 for(i3= 0; i3< myGrid.n3;  i3++)
 {
  for(i1= 0; i1< myGrid.n1;  i1++)
  {
   tmpData[index] = myGrid(i1,i2,i3,groupNumber,myGrid.cell[groupNumber],k) ;
   index +=stride;  
  }
 } 
}

template <int N,int numblock, typename dataType1>
void modifyComputationalGridPlane12Node_Noise(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i3, int k)
{  
 int i1, i2 ; 
 int index = beginIndex; 
 for(i2= myGrid.nB2 - myGrid.nB2 ; i2<= myGrid.nE2 + myGrid.nB2;  i2++)
 {
  for(i1= myGrid.nB1 - myGrid.nB1; i1<= myGrid.nE1 + myGrid.nB1;  i1++)
  {           
    myGrid(i1,i2,i3,groupNumber,myGrid.node[groupNumber],k) = tmpData[index];   
   index +=stride;  
  }
 } 
}

template <int N,int numblock, typename dataType1>
void modifyComputationalGridPlane12Cell_Noise(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i3, int k)
{  
 int i1, i2 ; 
 int index = beginIndex; 
 for(i2= myGrid.nB2 - myGrid.nB2 ; i2<= myGrid.nE2 + myGrid.nB2;  i2++)
 {
  for(i1= myGrid.nB1 - myGrid.nB1; i1<= myGrid.nE1 + myGrid.nB1;  i1++)
  {           
    myGrid(i1,i2,i3,groupNumber,myGrid.cell[groupNumber],k) = tmpData[index];   
   index +=stride;  
  }
 } 
}

template <int N,int numblock, typename dataType1>
void modifyComputationalGridPlane23Node_Noise(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i1, int k)
{  
 int i3, i2 ; 
 int index = beginIndex; 
 for(i3= myGrid.nB3 -  myGrid.nB3  ; i3<= myGrid.nE3 +  myGrid.nB3 ;  i3++)
 {
  for(i2= myGrid.nB2 - myGrid.nB2; i2<= myGrid.nE2 + myGrid.nB2;  i2++)
  {           
    myGrid(i1,i2,i3,groupNumber,myGrid.node[groupNumber],k) = tmpData[index];   
   index +=stride;  
  }
 } 
}

template <int N,int numblock, typename dataType1>
void modifyComputationalGridPlane23Cell_Noise(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i1, int k)
{  
 int i3, i2 ; 
 int index = beginIndex; 
 for(i3= myGrid.nB3 -  myGrid.nB3  ; i3<= myGrid.nE3 +  myGrid.nB3 ;  i3++)
 {
  for(i2= myGrid.nB2 - myGrid.nB2; i2<= myGrid.nE2 + myGrid.nB2;  i2++)
  {           
    myGrid(i1,i2,i3,groupNumber,myGrid.cell[groupNumber],k) = tmpData[index];   
   index +=stride;  
  }
 } 
}

template <int N,int numblock, typename dataType1>
void modifyComputationalGridPlane31Node_Noise(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i2, int k)
{  
 int i1, i3 ; 
 int index = beginIndex; 
 for(i3= myGrid.nB3 - myGrid.nB3; i3<= myGrid.nE3 + myGrid.nB3;  i3++)
 {
  for(i1= myGrid.nB1 - myGrid.nB1; i1<= myGrid.nE1 + myGrid.nB1;  i1++)
  {           
    myGrid(i1,i2,i3,groupNumber,myGrid.node[groupNumber],k) = tmpData[index];   
   index +=stride;  
  }
 } 
}

template <int N,int numblock, typename dataType1>
void modifyComputationalGridPlane31Cell_Noise(gridBCC3D<N, numblock, dataType1> &myGrid, int groupNumber, dataType1 *tmpData, int beginIndex,   int stride, int  i2, int k)
{  
 int i1, i3 ; 
 int index = beginIndex; 
 for(i3= myGrid.nB3 - myGrid.nB3; i3<= myGrid.nE3 + myGrid.nB3;  i3++)
 {
  for(i1= myGrid.nB1 - myGrid.nB1; i1<= myGrid.nE1 + myGrid.nB1;  i1++)
  {           
    myGrid(i1,i2,i3,groupNumber,myGrid.cell[groupNumber],k) = tmpData[index];   
   index +=stride;  
  }
 } 
}



// Explicit declarations
template void probeComputationalGridPlane12Node_Noise<4,1,double>(gridBCC3D<4,1,double> &, int , double *, int ,   int , int  , int );
template void probeComputationalGridPlane12Cell_Noise<4,1,double>(gridBCC3D<4,1,double> &, int , double *, int ,   int , int  , int );
template void probeComputationalGridPlane23Node_Noise<4,1,double>(gridBCC3D<4,1,double> &, int , double *, int ,   int , int  , int );
template void probeComputationalGridPlane23Cell_Noise<4,1,double>(gridBCC3D<4,1,double> &, int , double *, int ,   int , int  , int );
template void probeComputationalGridPlane31Node_Noise<4,1,double>(gridBCC3D<4,1,double> &, int , double *, int ,   int , int  , int );
template void probeComputationalGridPlane31Cell_Noise<4,1,double>(gridBCC3D<4,1,double> &, int , double *, int ,   int , int  , int );

template void modifyComputationalGridPlane12Node_Noise<4,1,double>(gridBCC3D<4,1,double> &, int , double *, int ,   int , int  , int );
template void modifyComputationalGridPlane12Cell_Noise<4,1,double>(gridBCC3D<4,1,double> &, int , double *, int ,   int , int  , int );
template void modifyComputationalGridPlane23Node_Noise<4,1,double>(gridBCC3D<4,1,double> &, int , double *, int ,   int , int  , int );
template void modifyComputationalGridPlane23Cell_Noise<4,1,double>(gridBCC3D<4,1,double> &, int , double *, int ,   int , int  , int );
template void modifyComputationalGridPlane31Node_Noise<4,1,double>(gridBCC3D<4,1,double> &, int , double *, int ,   int , int  , int );
template void modifyComputationalGridPlane31Cell_Noise<4,1,double>(gridBCC3D<4,1,double> &, int , double *, int ,   int , int  , int );

















