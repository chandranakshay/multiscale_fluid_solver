#include<restart.h>

template <int N,int numblock, typename dataType1>
void dumpRestart(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,char* restartFile)
{
	std::ofstream file2;
	file2.open(restartFile,std::ios::out|std::ios::binary);

    for(int k=myGrid.nB3; k<=myGrid.nE3;k++)
        for(int j=myGrid.nB2; j<=myGrid.nE2;j++)
            for(int i=myGrid.nB1; i<=myGrid.nE1;i++)
	    {
			   file2.write((char *)&myGrid(i,j,k,lbModel.G0 ,myGrid.node[lbModel.G0] ,0),   sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G1 ,myGrid.node[lbModel.G1] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G2 ,myGrid.node[lbModel.G2] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G3 ,myGrid.node[lbModel.G3] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G4 ,myGrid.node[lbModel.G4] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G5 ,myGrid.node[lbModel.G5] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G6 ,myGrid.node[lbModel.G6] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],0), 4*sizeof(dataType1));

			   file2.write((char *)&myGrid(i,j,k,lbModel.G0 ,myGrid.cell[lbModel.G0] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G1 ,myGrid.cell[lbModel.G1] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G2 ,myGrid.cell[lbModel.G2] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G3 ,myGrid.cell[lbModel.G3] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G4 ,myGrid.cell[lbModel.G4] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G5 ,myGrid.cell[lbModel.G5] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G6 ,myGrid.cell[lbModel.G6] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9] ,0), 4*sizeof(dataType1));
			   file2.write((char *)&myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],0), 4*sizeof(dataType1));

	   }

	file2.close();

}

template <int N,int numblock, typename dataType1>
void restartIC(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,char* restartFile)
{

	std::ifstream file2;
	file2.open(restartFile,std::ios::in|std::ios::binary);

    for(int k=myGrid.nB3; k<=myGrid.nE3;k++)
        for(int j=myGrid.nB2; j<=myGrid.nE2;j++)
            for(int i=myGrid.nB1; i<=myGrid.nE1;i++)
	    {
			   file2.read((char *)&myGrid(i,j,k,lbModel.G0 ,myGrid.node[lbModel.G0] ,0),   sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G1 ,myGrid.node[lbModel.G1] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G2 ,myGrid.node[lbModel.G2] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G3 ,myGrid.node[lbModel.G3] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G4 ,myGrid.node[lbModel.G4] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G5 ,myGrid.node[lbModel.G5] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G6 ,myGrid.node[lbModel.G6] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],0), 4*sizeof(dataType1));

			   file2.read((char *)&myGrid(i,j,k,lbModel.G0 ,myGrid.cell[lbModel.G0] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G1 ,myGrid.cell[lbModel.G1] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G2 ,myGrid.cell[lbModel.G2] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G3 ,myGrid.cell[lbModel.G3] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G4 ,myGrid.cell[lbModel.G4] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G5 ,myGrid.cell[lbModel.G5] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G6 ,myGrid.cell[lbModel.G6] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9] ,0), 4*sizeof(dataType1));
			   file2.read((char *)&myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],0), 4*sizeof(dataType1));
			   
			   /*std::cout<<myGrid(i,j,k,lbModel.G0 ,myGrid.node[lbModel.G0] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G1 ,myGrid.node[lbModel.G1] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G2 ,myGrid.node[lbModel.G2] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G3 ,myGrid.node[lbModel.G3] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G4 ,myGrid.node[lbModel.G4] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G5 ,myGrid.node[lbModel.G5] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G6 ,myGrid.node[lbModel.G6] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],0)<<std::endl;

			   std::cout<<myGrid(i,j,k,lbModel.G0 ,myGrid.cell[lbModel.G0] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G1 ,myGrid.cell[lbModel.G1] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G2 ,myGrid.cell[lbModel.G2] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G3 ,myGrid.cell[lbModel.G3] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G4 ,myGrid.cell[lbModel.G4] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G5 ,myGrid.cell[lbModel.G5] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G6 ,myGrid.cell[lbModel.G6] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9] ,0)<<std::endl;
			   std::cout<<myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],0)<<std::endl;*/

            }

	file2.close();

}

template void dumpRestart<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,char *) ;
template void restartIC<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,char *) ;
