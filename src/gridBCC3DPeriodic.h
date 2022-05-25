template <int numfield, int numBlock, typename dataType>
void periodicX1(gridBCC3D<numfield, numBlock,dataType> &myGrid)
{
 // Left 
   for (int block=0; block<numBlock; block++)
     for(int nodeType =0; nodeType<2;nodeType++)
 for (int i3=0; i3<myGrid.n3; i3++) 
  for (int i2=0; i2<myGrid.n2; i2++)
    for (int temp=0; temp<myGrid.nB1; temp++)      
   for(int dv=0; dv<numfield; dv++)
      myGrid(temp,i2,i3,block,nodeType, dv)=myGrid(myGrid.nE1-(myGrid.nB1-1-temp),i2,i3,block,nodeType, dv);         
    
 //Right
    for (int block=0; block<numBlock; block++)
     for(int nodeType =0; nodeType<2;nodeType++)
 for (int i3=0; i3<myGrid.n3; i3++) 
  for (int i2=0; i2<myGrid.n2; i2++)
    for(int temp=myGrid.nB1; temp>0; temp--)     
   for(int dv=0; dv<numfield; dv++)
       myGrid(myGrid.nE1+temp,i2,i3,block,nodeType, dv)=myGrid(myGrid.nB1-1+temp,i2,i3,block,nodeType, dv); 
    
} 


template <int numfield, int numBlock, typename dataType>
void periodicX2(gridBCC3D<numfield, numBlock,dataType> &myGrid)
{
 // Left 
   for (int block=0; block<numBlock; block++)
     for(int nodeType =0; nodeType<2;nodeType++)
 for (int i3=0; i3<myGrid.n3; i3++) 
    for(int temp=0;temp<myGrid.nB2;temp++)
 for (int i1=0; i1< myGrid.n1;  i1++)
   for(int dv=0; dv<numfield; dv++)
      myGrid(i1,temp,i3,block,nodeType, dv)=myGrid(i1,myGrid.nE2-(myGrid.nB2-1-temp),i3,block,nodeType, dv);         
    
 //Right
    for (int block=0; block<numBlock; block++)
     for(int nodeType =0; nodeType<2;nodeType++)
 for (int i3=0; i3<myGrid.n3; i3++) 
    for(int temp=myGrid.nB2;temp>0;temp--)
  for (int i1=0; i1< myGrid.n1;  i1++)
   for(int dv=0; dv<numfield; dv++)
       myGrid(i1,myGrid.nE2+temp,i3,block,nodeType, dv)=myGrid(i1,myGrid.nB2-1+temp,i3,block,nodeType, dv); 
    
} 


template <int numfield, int numBlock, typename dataType>
inline void periodicX3(gridBCC3D<numfield, numBlock,dataType> &myGrid)
{
 // Left 
   for (int block=0; block<numBlock; block++)
     for(int nodeType =0; nodeType<2;nodeType++)
for(int temp=0;temp<myGrid.nB3;temp++)    
  for (int i2=0; i2<myGrid.n2; i2++) 
 for (int i1=0; i1< myGrid.n1;  i1++)
   for(int dv=0; dv<numfield; dv++)
      myGrid(i1,i2,temp,block,nodeType, dv)=myGrid(i1, i2,myGrid.nE3-(myGrid.nB3-1-temp),block,nodeType, dv);         
    
 //Right
    for (int block=0; block<numBlock; block++)
     for(int nodeType =0; nodeType<2;nodeType++)
       for(int temp=myGrid.nB3;temp>0;temp--)
 for (int i2=0; i2<myGrid.n2; i2++) 
  for (int i1=0; i1< myGrid.n1;  i1++)
   for(int dv=0; dv<numfield; dv++)
       myGrid(i1,i2, myGrid.nE2+temp, block,nodeType, dv)=myGrid(i1,i2, myGrid.nB3-1+temp,block,nodeType, dv); 
    
} 
 
  
template <int numfield, int numBlock, typename dataType>
inline void periodicX12(gridBCC3D<numfield, numBlock,dataType> &myGrid)
{
 periodicX1(myGrid);
 periodicX2(myGrid);
}

template <int numfield, int numBlock, typename dataType>
inline void periodicX23(gridBCC3D<numfield, numBlock,dataType> &myGrid)
{
 periodicX2(myGrid);
 periodicX3(myGrid);
}
template <int numfield, int numBlock, typename dataType>
inline void periodicX13(gridBCC3D<numfield, numBlock,dataType> &myGrid)
{
 periodicX1(myGrid);
 periodicX3(myGrid);
}
template <int numfield, int numBlock, typename dataType>
inline void periodicAll(gridBCC3D<numfield, numBlock,dataType> &myGrid)
 
{
 periodicX1(myGrid); 
 periodicX2(myGrid);
 periodicX3(myGrid);
}
