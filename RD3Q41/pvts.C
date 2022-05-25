#include"pvts.h"



template <int N,int numblock, typename dataType1>
void dumpResultsNodeToPVTI(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,procInfo3D  CartersianGrid3D,int myRank,int x1,int x2,int x3,int nP1, int nP2, int nP3,int counter,dataType1 F1, dataType1 F2, dataType1 F3,dataType1 dt)
{
 std::ofstream file;
 char fileName[150];
 sprintf(fileName,"./results/result_node_%d_%d.vti",counter,myRank);
 file.open(fileName);
 int point(0);
 dataType1 rho1(0.0),uX1(0.0),uY1(0.0),uZ1(0.0),theta1(0.0);

 int ix(0),jy(0),kz(0);
//  int ix1(0),jy1(0),kz1(0);
 if(CartersianGrid3D.coordinates[0]==nP1-1)
  ix = 1;
 else
  ix = 0;
 if(CartersianGrid3D.coordinates[1]==nP2-1)
  jy = 1;
 else
  jy = 0;
 if(CartersianGrid3D.coordinates[2]==nP3-1)
  kz = 1;
 else
  kz = 0;

//   int VECT_LENGTH = 4;
 file<<"<?xml version=\"1.0\"?>   "<<std::endl;
 file<<"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<std::endl;


 file<<"  <ImageData WholeExtent=\""<<x1*CartersianGrid3D.coordinates[0]<<" "<<x1-1+x1*CartersianGrid3D.coordinates[0]+1-ix<<" "<<x2*CartersianGrid3D.coordinates[1]<<" "<<x2-1+x2*CartersianGrid3D.coordinates[1]+1-jy<<" "<<x3*CartersianGrid3D.coordinates[2]<<" "<<x3-1+x3*CartersianGrid3D.coordinates[2]+1-kz<<"\" Origin=\"0 0 0\" Spacing=\"1 1 1\">"<<std::endl;



 file<<"    <Piece Extent=\""<<x1*CartersianGrid3D.coordinates[0]<<" "<<x1-1+x1*CartersianGrid3D.coordinates[0]+1-ix<<" "<<x2*CartersianGrid3D.coordinates[1]<<" "<<x2-1+x2*CartersianGrid3D.coordinates[1]+1-jy<<" "<<x3*CartersianGrid3D.coordinates[2]<<" "<<x3-1+x3*CartersianGrid3D.coordinates[2]+1-kz<<"\">"<<std::endl;


 file<<"       <PointData Vectors=\"Velocity\">"<<std::endl;
 file<<"        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\"  format=\"ascii\">"<<std::endl;  //Change to int


 for(int i3= myGrid.nB3;  i3<= myGrid.nE3+(1-kz);  i3++)
  for(int i2= myGrid.nB2;  i2<= myGrid.nE2+(1-jy);  i2++)
   for(int i1= myGrid.nB1;  i1<= myGrid.nE1+(1-ix);  i1++)
   {
     copyFromNodeSinglePoint(lbModel,myGrid,i1,i2,i3,point);
//      F1 = 0.0; F2 = 0.0; F3=0.0;
      getHydroMomentSinglePointWithForce(lbModel,rho1,uX1,uY1,uZ1,theta1,F1,F2,F3,dt); ;//getHydroMomentSinglePoint(lbModel,rho1,uX1,uY1,uZ1,theta1);//
     file<<uX1<<" "<<uY1<<" "<<uZ1<<std::endl;
    }
 file<<"         </DataArray>"<<std::endl;


 file<<"        <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\"  format=\"ascii\">"<<std::endl;  //Change to int

 for(int i3= myGrid.nB3;  i3<= myGrid.nE3+(1-kz);  i3++)
  for(int i2= myGrid.nB2;  i2<= myGrid.nE2+(1-jy);  i2++)
   for(int i1= myGrid.nB1;  i1<= myGrid.nE1+(1-ix);  i1++)
   {
     copyFromNodeSinglePoint(lbModel,myGrid,i1,i2,i3,point);
//      F1 = 0.0; F2 = 0.0; F3=0.0;
     getHydroMomentSinglePointWithForce(lbModel,rho1,uX1,uY1,uZ1,theta1,F1,F2,F3,dt); //getHydroMomentSinglePoint(lbModel,rho1,uX1,uY1,uZ1,theta1);//
     file<<rho1<<std::endl;
    }
 file<<"         </DataArray>"<<std::endl;


 file<<"        <DataArray type=\"Float32\" Name=\"Temperature\" NumberOfComponents=\"1\"  format=\"ascii\">"<<std::endl;  //Change to int

 for(int i3= myGrid.nB3;  i3<= myGrid.nE3+(1-kz);  i3++)
  for(int i2= myGrid.nB2;  i2<= myGrid.nE2+(1-jy);  i2++)
   for(int i1= myGrid.nB1;  i1<= myGrid.nE1+(1-ix);  i1++)
   {
     copyFromNodeSinglePoint(lbModel,myGrid,i1,i2,i3,point);
//      F1 = 0.0; F2 = 0.0; F3=0.0;
     getHydroMomentSinglePointWithForce(lbModel,rho1,uX1,uY1,uZ1,theta1,F1,F2,F3,dt); //getHydroMomentSinglePoint(lbModel,rho1,uX1,uY1,uZ1,theta1);//
     file<<theta1<<std::endl;
    }
 file<<"         </DataArray>"<<std::endl;

 file<<"       </PointData>"<<std::endl;



 file<<"    </Piece>"<<std::endl;
 file<<"  </ImageData>"<<std::endl;
 file<<"</VTKFile>"<<std::endl;
 file.close();
}

void writeResultsNodeMasterPVTI(int x1,int x2,int x3,int nP1, int nP2, int nP3,int counter)
{
 std::ofstream file;
 char fileName[150];
 sprintf(fileName,"./results/resultsNode%d.pvti",counter);
 file.open(fileName);
 file<<"<?xml version=\"1.0\"?>"<<std::endl;
 file<<"<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\"> "<<std::endl;
 file<<"  <PImageData WholeExtent=\""<<"0 "<<nP1*x1-1<<" 0 "<<nP2*x2-1<<" 0 "<<nP3*x3-1<<"\" Origin=\"0 0 0\" Spacing=\"1 1 1\" GhostLevel=\"0\">"<<std::endl;

 file<<"    <PPointData Vectors=\"Velocity\" >"<<std::endl;
 file<<"      <PDataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\"/>"<<std::endl;
 file<<"      <PDataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\"/>"<<std::endl;
 file<<"      <PDataArray type=\"Float32\" Name=\"Temperature\" NumberOfComponents=\"1\"/>"<<std::endl;
 file<<"   </PPointData>"<<std::endl;


 int count =0;
 int ix(0),jy(0),kz(0);//,kz1(0);
 for(int k=0;k<nP3;k++)
  for(int j=0;j<nP2;j++)
   for(int i=0;i<nP1;i++)
   {
    if(i==nP1-1)
     ix = 1;
    else
     ix = 0;
    if(j==nP2-1)
     jy = 1;
    else
     jy = 0;
    if(k==nP3-1)
     kz = 1;
    else
     kz = 0;


    file<<"    <Piece Extent=\""<<i*x1<<" "<<(i+1)*x1-1+(1-ix)<<" "<<j*x2<<" "<<(j+1)*x2-1+(1-jy)<<" "<<k*x3<<" "<<(k+1)*x3-1+(1-kz)<<" "   <<"\" Source=\"result_node_"<<counter<<"_"<<count<<".vti\"/>"<<std::endl;
    count++;


   }
 file<<"  </PImageData>"<<std::endl;
 file<<"</VTKFile>"<<std::endl;
 file.close();
}

template <typename dataType1>
void dumpMarkerNodeToPVTI(lbmRD3Q41<dataType1> &lbModel,gridBCC3D<1, 1, int> &marker,procInfo3D  CartersianGrid3D,int myRank,int x1,int x2,int x3,int nP1, int nP2, int nP3,int counter)
{
 std::ofstream file;
 char fileName[150];
 sprintf(fileName,"./results/marker_node_%d_%d.vti",counter,myRank);
 file.open(fileName);
//  int point(0);
//  dataType1 rho1(0.0),uX1(0.0),uY1(0.0),uZ1(0.0),theta1(0.0);

 int ix(0),jy(0),kz(0);
//  int ix1(0),jy1(0),kz1(0);
 if(CartersianGrid3D.coordinates[0]==nP1-1)
  ix = 1;
 else
  ix = 0;
 if(CartersianGrid3D.coordinates[1]==nP2-1)
  jy = 1;
 else
  jy = 0;
 if(CartersianGrid3D.coordinates[2]==nP3-1)
  kz = 1;
 else
  kz = 0;


 file<<"<?xml version=\"1.0\"?>   "<<std::endl;
 file<<"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<std::endl;


 file<<"  <ImageData WholeExtent=\""<<x1*CartersianGrid3D.coordinates[0]<<" "<<x1-1+x1*CartersianGrid3D.coordinates[0]+1-ix<<" "<<x2*CartersianGrid3D.coordinates[1]<<" "<<x2-1+x2*CartersianGrid3D.coordinates[1]+1-jy<<" "<<x3*CartersianGrid3D.coordinates[2]<<" "<<x3-1+x3*CartersianGrid3D.coordinates[2]+1-kz<<"\" Origin=\"0 0 0\" Spacing=\"1 1 1\">"<<std::endl;



 file<<"    <Piece Extent=\""<<x1*CartersianGrid3D.coordinates[0]<<" "<<x1-1+x1*CartersianGrid3D.coordinates[0]+1-ix<<" "<<x2*CartersianGrid3D.coordinates[1]<<" "<<x2-1+x2*CartersianGrid3D.coordinates[1]+1-jy<<" "<<x3*CartersianGrid3D.coordinates[2]<<" "<<x3-1+x3*CartersianGrid3D.coordinates[2]+1-kz<<"\">"<<std::endl;


 file<<"       <PointData Vectors=\"Velocity\">"<<std::endl;
 file<<"        <DataArray type=\"Float32\" Name=\"Marker\" NumberOfComponents=\"1\"  format=\"ascii\">"<<std::endl;  //Change to int

 for(int i3= marker.nB3;  i3<= marker.nE3+(1-kz);  i3++)
  for(int i2= marker.nB2;  i2<= marker.nE2+(1-jy);  i2++)
   for(int i1= marker.nB1;  i1<= marker.nE1+(1-ix);  i1++)
   {
     file<<marker(i1,i2,i3,0,0,0)<<std::endl;
   }
 file<<"         </DataArray>"<<std::endl;
 file<<"       </PointData>"<<std::endl;



 file<<"    </Piece>"<<std::endl;
 file<<"  </ImageData>"<<std::endl;
 file<<"</VTKFile>"<<std::endl;
 file.close();
}

void writeResultsNodeMarkerPVTI(int x1,int x2,int x3,int nP1, int nP2, int nP3,int counter)
{
 std::ofstream file;
 char fileName[150];
 sprintf(fileName,"./results/markerNode%d.pvti",counter);
 file.open(fileName);
 file<<"<?xml version=\"1.0\"?>"<<std::endl;
 file<<"<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\"> "<<std::endl;
 file<<"  <PImageData WholeExtent=\""<<"0 "<<nP1*x1-1<<" 0 "<<nP2*x2-1<<" 0 "<<nP3*x3-1<<"\" Origin=\"0 0 0\" Spacing=\"1 1 1\" GhostLevel=\"0\">"<<std::endl;

 file<<"    <PPointData Vectors=\"Velocity\" >"<<std::endl;
 file<<"      <PDataArray type=\"Float32\" Name=\"Marker\" NumberOfComponents=\"1\"/>"<<std::endl;
 file<<"   </PPointData>"<<std::endl;


 int count =0;
 int ix(0),jy(0),kz(0);//,kz1(0);
 for(int k=0;k<nP3;k++)
  for(int j=0;j<nP2;j++)
   for(int i=0;i<nP1;i++)
   {
    if(i==nP1-1)
     ix = 1;
    else
     ix = 0;
    if(j==nP2-1)
     jy = 1;
    else
     jy = 0;
    if(k==nP3-1)
     kz = 1;
    else
     kz = 0;


    file<<"    <Piece Extent=\""<<i*x1<<" "<<(i+1)*x1-1+(1-ix)<<" "<<j*x2<<" "<<(j+1)*x2-1+(1-jy)<<" "<<k*x3<<" "<<(k+1)*x3-1+(1-kz)<<" "   <<"\" Source=\"marker_node_"<<counter<<"_"<<count<<".vti\"/>"<<std::endl;
    count++;


   }
 file<<"  </PImageData>"<<std::endl;
 file<<"</VTKFile>"<<std::endl;
 file.close();
}

template <int N,int numblock, typename dataType1>
inline void averageVelocity(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step,int size,int myRank,dataType1 F1,dataType1 F2,dataType1 F3,dataType1 uRef,dataType1 &averageVelocity1,dataType1 dt)
{
    dataType1 nodeMass(0.0),cellMass(0.0);
    int nodeCount(0),cellCount(0);
    dataType1 rho1,uX1,uY1,uZ1,theta1;
    int nodeCountGlobal(0.0);
    int cellCountGlobal(0.0);
    dataType1 nodeMassGlobal(0.0);
    dataType1 cellMassGlobal(0.0);
    dataType1 nodeVelocity(0.0),cellVelocity(0.0);
    int point = 0;
    dataType1 nodeVelocityGlobal(0.0);
    dataType1 cellVelocityGlobal(0.0);
    for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
            {
                if(marker(i1,i2,i3,0,0,0) == FLUID)
                {
                   copyFromNodeSinglePoint(lbModel,myGrid,i1,i2,i3,point);
                   getHydroMomentSinglePointWithForce(lbModel,rho1,uX1,uY1,uZ1,theta1,F1,F2,F3,dt);
                   nodeVelocity += uX1;
                   nodeCount++;
                }
            }

            for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
                for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
                    for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
                    {
                        if(marker(i1,i2,i3,0,1,0) == FLUID)
                        {
                         copyFromCellSinglePoint(lbModel,myGrid,i1,i2,i3,point);
                         getHydroMomentSinglePointWithForce(lbModel,rho1,uX1,uY1,uZ1,theta1,F1,F2,F3,dt);
                         cellVelocity += uX1;
                         cellCount++;
                        }
                    }

                    //             nodeMass = nodeMass;
                    //             cellMass = cellMass;

                    MPI_Reduce(&nodeCount,&nodeCountGlobal, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
                    MPI_Reduce(&cellCount,&cellCountGlobal, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
                    MPI_Reduce(&nodeVelocity ,&nodeVelocityGlobal, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
                    MPI_Reduce(&cellVelocity ,&cellVelocityGlobal, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
                    //             int density = 1;


                    if(myRank == 0)
                    {
                        averageVelocity1 =  (cellVelocityGlobal+nodeVelocityGlobal)/((cellCountGlobal+nodeCountGlobal));
//                         std::cout<<"FLUID: "<< std::fixed<<std::setw(14)<<step<<"    "<<std::setw(16)<<averageVelocity1<<"      "<<std::setw(16)<<averageVelocity1/uRef<<std::setw(16)<< ((averageVelocity1/uRef * (3.0/2.0))-1.0)*100.0 <<std::endl;
                    }
                                            MPI_Bcast(&averageVelocity1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}


// Explicit definitions


template void dumpMarkerNodeToPVTI<double>(lbmRD3Q41<double> &,gridBCC3D<1, 1, int> &,procInfo3D  ,int ,int ,int ,int ,int , int , int ,int );

template void dumpResultsNodeToPVTI<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,procInfo3D ,int ,int ,int ,int ,int , int , int ,int ,double , double , double, double);

template void averageVelocity<4,11,double>(lbmRD3Q41<double> &lbModel, gridBCC3D<4,11,double> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step,int size,int myRank,double F1,double F2,double F3,double uRef,double &averageVelocity1,double dt);
