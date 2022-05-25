#include"inletOutlet.h"

template <int N,int numblock, typename dataType1>
void X1BeginInlet(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 u_inlet, dataType1 cs,int step,int VECT_LENGTH)
{
    dataType1 *rho,*uX,*uY,*uZ,*theta,*Pxx,*Pyy,*Pzz,*Pxy,*Pyz,*Pzx; 
    dataType1 *rhoI,*uXI,*uYI,*uZI,*thetaI;//,*PxxI,*PyyI,*PzzI,*PxyI,*PyzI,*PzxI; 
    rho   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              rhoI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    uX    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              uXI    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    uY    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              uYI    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    uZ    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              uZI    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    theta = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              thetaI = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32); 
    Pxx   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              //PxxI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    Pyy   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              //PyyI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    Pzz   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              //PzzI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    Pxy   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              //PxyI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    Pyz   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              //PyzI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32); 
    Pzx   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              //PzxI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);    
    
    dataType1 grad[lbModel.dvN];       
    
    for (int k =0; k <  myGrid.n3; k++) 
    {
        for (int j =0; j <  myGrid.n2; j++)  
        {        
/*                for(int i=0;i<VECT_LENGTH;i++)
                {
                    thetaI[i] = lbModel.theta0;           
                    rhoI  [i] = 1.0;
                    uXI   [i] = u_inlet;
                    uYI   [i] = 0.0;
                    uZI   [i] = 0.0;
                }
//                         getHydroMomentsFromNode(lbModel,myGrid,myGrid.nB1,j,k,rhoI,uXI,uYI,uZI,thetaI);
                        getFEq(lbModel,VECT_LENGTH,rhoI,uXI,uYI,uZI,thetaI);  
                        copyToNodeSinglePoint(lbModel,myGrid,myGrid.nB1-1,j,k,0);
                        copyToNodeSinglePoint(lbModel,myGrid,myGrid.nB1-2,j,k,0);
                        copyToNodeSinglePoint(lbModel,myGrid,myGrid.nB1-3,j,k,0);
                        copyToNodeSinglePoint(lbModel,myGrid,myGrid.nB1-4,j,k,0);
                        copyToCellSinglePoint(lbModel,myGrid,myGrid.nB1-1,j,k,0);
                        copyToCellSinglePoint(lbModel,myGrid,myGrid.nB1-2,j,k,0);
                        copyToCellSinglePoint(lbModel,myGrid,myGrid.nB1-3,j,k,0);   
                        copyToCellSinglePoint(lbModel,myGrid,myGrid.nB1-4,j,k,0);  */          
            copyFromNode(lbModel,myGrid,VECT_LENGTH,myGrid.nB1,j,k);
            getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);  
            
            for(int i=0;i<VECT_LENGTH;i++)
            {
                thetaI[i] = lbModel.theta0;//+0.02;           
                rhoI  [i] = 1.0;//rho[i];
                uXI   [i] = u_inlet;///rhoI[i];
                uYI   [i] = 0.0;
                uZI   [i] = 0.0;
//                 PxxI  [i] = Pxx[i] - rho[i]*(theta[i] + uX[i]*uX[i]) + lbModel.theta0 + u_inlet*u_inlet;// + (u_inlet - uX[i])*cs ;
                //                  PxxI[0] = 3.0*rhoI[0]*thetaI[0]+rhoI[0]*(uXI[0]*uXI[0]+uYI[0]*uYI[0]+uZI[0]*uZI[0])-Pyy[0]-Pzz[0];   
            }
//             
            getGradThermalSinglePoint(lbModel, rhoI[0], uXI[0],uYI[0],uZI[0],thetaI[0],Pxx[0],Pyy[0],Pzz[0],Pxy[0],Pyz[0],Pzx[0],grad);
            
            copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nB1-1,j,k);
            copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nB1-2,j,k);
            copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nB1-3,j,k);
            copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nB1-4,j,k);
            copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-1,j,k);
            copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-2,j,k);
            copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-3,j,k);
            copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-4,j,k);
        }
    }
} 

template <int N,int numblock, typename dataType1>
void X1EndOutlet(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH)
{
    //Outlet
    dataType1 *rho,*uX,*uY,*uZ,*theta,*Pxx,*Pyy,*Pzz,*Pxy,*Pyz,*Pzx; 
//     dataType1 *rhoI,*uXI,*uYI,*uZI,*thetaI,*PxxI,*PyyI,*PzzI,*PxyI,*PyzI,*PzxI; 
    rho   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);             // rhoI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    uX    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);             // uXI    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    uY    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);             // uYI    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    uZ    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);             // uZI    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    theta = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);             // thetaI = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32); 
    Pxx   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);             // PxxI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    Pyy   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);             // PyyI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    Pzz   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);             // PzzI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    Pxy   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);             // PxyI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
    Pyz   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);             // PyzI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32); 
    Pzx   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);             // PzxI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);   
    
    dataType1 grad[lbModel.dvN];       
    
    for (int k =0; k <  myGrid.n3; k++) 
    {
        for (int j =0; j <  myGrid.n2; j++) 
        {           
            copyFromCell(lbModel,myGrid,VECT_LENGTH,myGrid.nE1,j,k);
            getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx); 
            
            getGradThermalSinglePoint(lbModel, rho[0], uX[0],uY[0],uZ[0],theta[0],Pxx[0],Pyy[0],Pzz[0],Pxy[0],Pyz[0],Pzx[0],grad);
            
            copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+1,j,k);
            copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+2,j,k);
            copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+3,j,k);
            copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+4,j,k);
            copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+1,j,k);
            copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+2,j,k);
            copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+3,j,k);
            copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+4,j,k);              
        }
    }
}     



// template <typename dataType1> 
// inline void X1BeginInlet(lbmRD3Q41<dataType1> &lbModel, lbGridD3Q41<dataType1> &compGrid, dataType1 u_inlet, dataType1 cs,int step)
// {
// 	//Inlet
// 	dataType1 pXX,pYY,pZZ,pXY,pYZ,pZX;
// 	dataType1 pXXI,pYYI,pZZI,pXYI,pYZI,pZXI, rhoI, uXI, uYI, uZI, thetaI;
// 	dataType1 u_inlet1=u_inlet;
// 	//  if(step<10)
// 	//    u_inlet1 *= step/10.0;
// 	for (int k = 0; k < compGrid.CENTER.n3; k++) 
// 	{
// 		for (int j = 0; j < compGrid.CENTER.n2; j++) 
// 		{
// 			//Node-1
//                         copyFromNode(lbModel,myGrid,VECT_LENGTH,myGrid.nB1,j,k);
//                         getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);  
// 
// 
//                         
// 			thetaI = lbModel.theta0;//lbModel.theta0;// 
// 			rhoI   = lbModel.rho;
// 			uXI    = u_inlet1/rhoI;
// 			uYI    = 0.0;//lbModel.uY;//
// 			uZI    = 0.0;//lbModel.uZ;//
// 			pXXI = pXX - lbModel.rho*(lbModel.uX*lbModel.uX + lbModel.theta) + rhoI*(uXI*uXI + thetaI);
// 			pYYI = pYY - lbModel.rho*(lbModel.uY*lbModel.uY + lbModel.theta) + rhoI*(uYI*uYI + thetaI);  
// 			pZZI = pZZ - lbModel.rho*(lbModel.uZ*lbModel.uZ + lbModel.theta) + rhoI*(uZI*uZI + thetaI); 
// 			pXYI = pXY - lbModel.rho* lbModel.uX*lbModel.uY + rhoI*uXI*uYI;
// 			pYZI = pYZ - lbModel.rho* lbModel.uZ*lbModel.uY + rhoI*uZI*uYI; 
// 			pZXI = pZX - lbModel.rho* lbModel.uX*lbModel.uZ + rhoI*uXI*uZI; 
// 			lbModel.rho = rhoI;
// 			lbModel.uX  = uXI;
// 			lbModel.uY  = uYI;
// 			lbModel.uZ  = uZI;
// 			lbModel.theta = thetaI;
// 			getGrad(lbModel,lbModel.fTemp,pXXI,pYYI,pZZI,pXYI,pYZI,pZXI);      
// 		//	 getFEq(lbModel,lbModel.fTemp);
// 			copyToNode(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nB1,j,k);
// 
// 			copyFromNode(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nB1+3,j,k);
// 			getPMoments(lbModel,lbModel.fTemp,pXX,pYY,pZZ,pXY,pYZ,pZX);
// 			thetaI = lbModel.theta0;//lbModel.theta0;// 
// 			rhoI   = lbModel.rho;
// 			uXI    = u_inlet1/rhoI;
// 			uYI    = 0.0;//lbModel.uY;//
// 			uZI    = 0.0;//lbModel.uZ;//
// 			pXXI = pXX - lbModel.rho*(lbModel.uX*lbModel.uX + lbModel.theta) + rhoI*(uXI*uXI + thetaI);
// 			pYYI = pYY - lbModel.rho*(lbModel.uY*lbModel.uY + lbModel.theta) + rhoI*(uYI*uYI + thetaI);  
// 			pZZI = pZZ - lbModel.rho*(lbModel.uZ*lbModel.uZ + lbModel.theta) + rhoI*(uZI*uZI + thetaI); 
// 			pXYI = pXY - lbModel.rho* lbModel.uX*lbModel.uY + rhoI*uXI*uYI;
// 			pYZI = pYZ - lbModel.rho* lbModel.uZ*lbModel.uY + rhoI*uZI*uYI; 
// 			pZXI = pZX - lbModel.rho* lbModel.uX*lbModel.uZ + rhoI*uXI*uZI; 
// 			lbModel.rho = rhoI;
// 			lbModel.uX  = uXI;
// 			lbModel.uY  = uYI;
// 			lbModel.uZ  = uZI;
// 			lbModel.theta = thetaI;
// 			getGrad(lbModel,lbModel.fTemp,pXXI,pYYI,pZZI,pXYI,pYZI,pZXI);      
// 	//		getFEq(lbModel,lbModel.fTemp);
// 			copyToNode(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nB1+1,j,k);
// 
// 			// Cell-1
// 			copyFromCell(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nB1+2,j,k);
// 			getPMoments(lbModel,lbModel.fTemp,pXX,pYY,pZZ,pXY,pYZ,pZX);
// 			thetaI = lbModel.theta0;//lbModel.theta0;// 
// 			rhoI   = 1.0;//lbModel.rho;
// 			uXI    = u_inlet1/rhoI;
// 			uYI    = 0.0;//lbModel.uY;//
// 			uZI    = 0.0;//lbModel.uZ;//
// 			pXXI = pXX - lbModel.rho*(lbModel.uX*lbModel.uX + lbModel.theta) + rhoI*(uXI*uXI + thetaI);
// 			pYYI = pYY - lbModel.rho*(lbModel.uY*lbModel.uY + lbModel.theta) + rhoI*(uYI*uYI + thetaI);  
// 			pZZI = pZZ - lbModel.rho*(lbModel.uZ*lbModel.uZ + lbModel.theta) + rhoI*(uZI*uZI + thetaI); 
// 			pXYI = pXY - lbModel.rho* lbModel.uX*lbModel.uY + rhoI*uXI*uYI;
// 			pYZI = pYZ - lbModel.rho* lbModel.uZ*lbModel.uY + rhoI*uZI*uYI; 
// 			pZXI = pZX - lbModel.rho* lbModel.uX*lbModel.uZ + rhoI*uXI*uZI; 
// 			lbModel.rho = rhoI;
// 			lbModel.uX  = uXI;
// 			lbModel.uY  = uYI;
// 			lbModel.uZ  = uZI;
// 			lbModel.theta = thetaI;
// 			//getFEq(lbModel,lbModel.fTemp);
// 			getGrad(lbModel,lbModel.fTemp,pXXI,pYYI,pZZI,pXYI,pYZI,pZXI);      
// 			copyToCell(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nB1,j,k);
// 
// 			copyFromCell(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nB1+3,j,k);
// 			getPMoments(lbModel,lbModel.fTemp,pXX,pYY,pZZ,pXY,pYZ,pZX);
// 			thetaI = lbModel.theta0;//lbModel.theta0;// 
// 			rhoI   = 1.0;//lbModel.rho;
// 			uXI    = u_inlet1/rhoI;
// 			uYI    = 0.0;//lbModel.uY;//
// 			uZI    = 0.0;//lbModel.uZ;//
// 			pXXI = pXX - lbModel.rho*(lbModel.uX*lbModel.uX + lbModel.theta) + rhoI*(uXI*uXI + thetaI);
// 			pYYI = pYY - lbModel.rho*(lbModel.uY*lbModel.uY + lbModel.theta) + rhoI*(uYI*uYI + thetaI);  
// 			pZZI = pZZ - lbModel.rho*(lbModel.uZ*lbModel.uZ + lbModel.theta) + rhoI*(uZI*uZI + thetaI); 
// 			pXYI = pXY - lbModel.rho* lbModel.uX*lbModel.uY + rhoI*uXI*uYI;
// 			pYZI = pYZ - lbModel.rho* lbModel.uZ*lbModel.uY + rhoI*uZI*uYI; 
// 			pZXI = pZX - lbModel.rho* lbModel.uX*lbModel.uZ + rhoI*uXI*uZI; 
// 			lbModel.rho = rhoI;
// 			lbModel.uX  = uXI;
// 			lbModel.uY  = uYI;
// 			lbModel.uZ  = uZI;
// 			lbModel.theta = thetaI;
// 			//getFEq(lbModel,lbModel.fTemp);
// 			getGrad(lbModel,lbModel.fTemp,pXXI,pYYI,pZZI,pXYI,pYZI,pZXI);      
// 			copyToCell(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nB1+1,j,k);
// 
// 		}
// 	}
// }

/*
	template <typename dataType1> 
inline void X1EndOutlet(lbmRD3Q41<dataType1> &lbModel, lbGridD3Q41<dataType1> &compGrid)
{
	//Inlet
	dataType1 pXX,pYY,pZZ,pXY,pYZ,pZX;
	dataType1 pXXI,pYYI,pZZI,pXYI,pYZI,pZXI, rhoI, uXI, uYI, uZI, thetaI;
	//  if(step<10)
	//    u_inlet1 *= step/10.0;
	for (int k = 0; k < compGrid.CENTER.n3; k++) 
	{
		for (int j = 0; j < compGrid.CENTER.n2; j++) 
		{
			//Node-1
			copyFromNode(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nE1-2,j,k);
			getPMoments(lbModel,lbModel.fTemp,pXX,pYY,pZZ,pXY,pYZ,pZX);
			thetaI = lbModel.theta0;//lbModel.theta0;// 
			rhoI   = lbModel.rho;
			uXI    = lbModel.uX;
			uYI    = lbModel.uY;//
			uZI    = lbModel.uZ;//
			pXXI = pXX - lbModel.rho*(lbModel.uX*lbModel.uX + lbModel.theta) + rhoI*(uXI*uXI + thetaI);
			pYYI = pYY ;//- lbModel.rho*(lbModel.uY*lbModel.uY + lbModel.theta) + rhoI*(uYI*uYI + thetaI);  
			pZZI = pZZ ;//- lbModel.rho*(lbModel.uZ*lbModel.uZ + lbModel.theta) + rhoI*(uZI*uZI + thetaI); 
			pXYI = pXY ;//- lbModel.rho* lbModel.uX*lbModel.uY + rhoI*uXI*uYI;
			pYZI = pYZ ;//- lbModel.rho* lbModel.uZ*lbModel.uY + rhoI*uZI*uYI; 
			pZXI = pZX ;//- lbModel.rho* lbModel.uX*lbModel.uZ + rhoI*uXI*uZI; 
			lbModel.rho = rhoI;
			lbModel.uX  = uXI;
			lbModel.uY  = uYI;
			lbModel.uZ  = uZI;
			lbModel.theta = thetaI;
			getGrad(lbModel,lbModel.fTemp,pXXI,pYYI,pZZI,pXYI,pYZI,pZXI);      
			//getFEq(lbModel,lbModel.fTemp);
			copyToNode(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nE1-1,j,k);
			copyToNode(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nE1,j,k);
			//Node-2
			copyFromNode(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nE1-3,j,k);
			getPMoments(lbModel,lbModel.fTemp,pXX,pYY,pZZ,pXY,pYZ,pZX);
			thetaI = lbModel.theta0;//lbModel.theta0;// 
			rhoI   = lbModel.rho;
			uXI    = lbModel.uX;
			uYI    = lbModel.uY;//
			uZI    = lbModel.uZ;//
			pXXI = pXX - lbModel.rho*(lbModel.uX*lbModel.uX + lbModel.theta) + rhoI*(uXI*uXI + thetaI);
			pYYI = pYY ;//- lbModel.rho*(lbModel.uY*lbModel.uY + lbModel.theta) + rhoI*(uYI*uYI + thetaI);  
			pZZI = pZZ ;//- lbModel.rho*(lbModel.uZ*lbModel.uZ + lbModel.theta) + rhoI*(uZI*uZI + thetaI); 
			pXYI = pXY ;//- lbModel.rho* lbModel.uX*lbModel.uY + rhoI*uXI*uYI;
			pYZI = pYZ ;//- lbModel.rho* lbModel.uZ*lbModel.uY + rhoI*uZI*uYI; 
			pZXI = pZX ;//- lbModel.rho* lbModel.uX*lbModel.uZ + rhoI*uXI*uZI; 
			lbModel.rho = rhoI;
			lbModel.uX  = uXI;
			lbModel.uY  = uYI;
			lbModel.uZ  = uZI;
			lbModel.theta = thetaI;
			//getFEq(lbModel,lbModel.fTemp);
			getGrad(lbModel,lbModel.fTemp,pXXI,pYYI,pZZI,pXYI,pYZI,pZXI);      
			compGrid.SC1.Node(compGrid.CENTER.nE1-1,j,k,lbModel.SC1_DV_P2_ZERO_ZERO) = lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO]; 
			//
			// Cell-1
			copyFromCell(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nE1-2,j,k);
			getPMoments(lbModel,lbModel.fTemp,pXX,pYY,pZZ,pXY,pYZ,pZX);
			thetaI = lbModel.theta0;//lbModel.theta0;// 
			rhoI   = lbModel.rho;
			uXI    = lbModel.uX;
			uYI    = lbModel.uY;//
			uZI    = lbModel.uZ;//
			pXXI = pXX - lbModel.rho*(lbModel.uX*lbModel.uX + lbModel.theta) + rhoI*(uXI*uXI + thetaI);
			pYYI = pYY ;//- lbModel.rho*(lbModel.uY*lbModel.uY + lbModel.theta) + rhoI*(uYI*uYI + thetaI);  
			pZZI = pZZ ;//- lbModel.rho*(lbModel.uZ*lbModel.uZ + lbModel.theta) + rhoI*(uZI*uZI + thetaI); 
			pXYI = pXY ;//- lbModel.rho* lbModel.uX*lbModel.uY + rhoI*uXI*uYI;
			pYZI = pYZ ;//- lbModel.rho* lbModel.uZ*lbModel.uY + rhoI*uZI*uYI; 
			pZXI = pZX ;//- lbModel.rho* lbModel.uX*lbModel.uZ + rhoI*uXI*uZI; 
			lbModel.rho = rhoI;
			lbModel.uX  = uXI;
			lbModel.uY  = uYI;
			lbModel.uZ  = uZI;
			lbModel.theta = thetaI;
			//getFEq(lbModel,lbModel.fTemp);
			getGrad(lbModel,lbModel.fTemp,pXXI,pYYI,pZZI,pXYI,pYZI,pZXI);      
			copyToCell(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nE1,j,k);
			copyToCell(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nE1-1,j,k);
			// //Cell-2   
			copyFromCell(lbModel,compGrid,lbModel.fTemp,compGrid.CENTER.nE1-3,j,k);
			getPMoments(lbModel,lbModel.fTemp,pXX,pYY,pZZ,pXY,pYZ,pZX);
			thetaI = lbModel.theta0;//lbModel.theta0;// 
			rhoI   = lbModel.rho;
			uXI    = lbModel.uX;
			uYI    = lbModel.uY;//
			uZI    = lbModel.uZ;//
			pXXI = pXX - lbModel.rho*(lbModel.uX*lbModel.uX + lbModel.theta) + rhoI*(uXI*uXI + thetaI);
			pYYI = pYY ;//- lbModel.rho*(lbModel.uY*lbModel.uY + lbModel.theta) + rhoI*(uYI*uYI + thetaI);  
			pZZI = pZZ ;//- lbModel.rho*(lbModel.uZ*lbModel.uZ + lbModel.theta) + rhoI*(uZI*uZI + thetaI); 
			pXYI = pXY ;//- lbModel.rho* lbModel.uX*lbModel.uY + rhoI*uXI*uYI;
			pYZI = pYZ ;//- lbModel.rho* lbModel.uZ*lbModel.uY + rhoI*uZI*uYI; 
			pZXI = pZX ;//- lbModel.rho* lbModel.uX*lbModel.uZ + rhoI*uXI*uZI; 
			lbModel.rho = rhoI;
			lbModel.uX  = uXI;
			lbModel.uY  = uYI;
			lbModel.uZ  = uZI;
			lbModel.theta = thetaI;
			getGrad(lbModel,lbModel.fTemp,pXXI,pYYI,pZZI,pXYI,pYZI,pZXI);      
			//getFEq(lbModel,lbModel.fTemp);
			compGrid.SC1.Cell(compGrid.CENTER.nE1-1,j,k,lbModel.SC1_DV_P2_ZERO_ZERO) = lbModel.fTemp[lbModel.DV_P2_ZERO_ZERO]; 

		}
	}
}*/
















// template <int N,int numblock, typename dataType1>
// void X1BeginInlet(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, dataType1 u_inlet, dataType1 cs,int step,int VECT_LENGTH)
// {
//     dataType1 *rho,*uX,*uY,*uZ,*theta,*Pxx,*Pyy,*Pzz,*Pxy,*Pyz,*Pzx; 
//     dataType1 *rhoI,*uXI,*uYI,*uZI,*thetaI,*PxxI,*PyyI,*PzzI,*PxyI,*PyzI,*PzxI; 
//     rho   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              rhoI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     uX    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              uXI    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     uY    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              uYI    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     uZ    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              uZI    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     theta = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              thetaI = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32); 
//     Pxx   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              PxxI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     Pyy   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              PyyI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     Pzz   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              PzzI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     Pxy   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              PxyI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     Pyz   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              PyzI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32); 
//     Pzx   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              PzxI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);    
//     
//     dataType1 grad[lbModel.dvN];       
//     
//     //     for (int k = myGrid.ndB3; k <= myGrid.ndE3; k++) 
//     //     {
//     //         for (int j =  myGrid.ndB2; j <= myGrid.ndE2; j++) 
//     //         {        
//     for (int k =0; k <  myGrid.n3; k++) 
//     {
//         for (int j =0; j <  myGrid.n2; j++) 
//         {           
//             //         copyFromNode(lbModel,myGrid,VECT_LENGTH,myGrid.nB1,j,k);
//             //         getFEq(lbModel,VECT_LENGTH,rhoI,uXI,uYI,uZI,thetaI);  
//             copyFromNode(lbModel,myGrid,VECT_LENGTH,myGrid.nB1,j,k);
//             getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);  
//             //             getGrad10Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx); 
//             for(int i=0;i<VECT_LENGTH;i++)
//             {
//                 //                 if(k==4 && j==4)
//                 //                  std::cout<<rho[i]<<"  "<<uX[i]<<"  "<<uY[i]<<"  "<<uZ[i]<<"  "<<theta[i]<<"   "<<Pxx[i]<<"  "<<Pyy[i]<<"  "<<Pzz[i]<<"  "<<Pxy[i]<<"  "<<Pyz[i]<<"   "<<Pzx[i]<<std::endl;
//                 
//                 thetaI[i] = lbModel.theta0;//+0.02;           
//                 rhoI  [i] = 1.0;
//                 uXI   [i] = u_inlet;
//                 uYI   [i] = 0.0;
//                 uZI   [i] = 0.0;
//                 //                 PxxI  [i] = Pxx[i] ;//- rho[i]*theta[i] - rho[i]*uX[i]*uX[i] + lbModel.theta0 + u_inlet*u_inlet;// + (u_inlet - uX[i])*cs ;
// //                 PxxI[0] = 3.0*rhoI[0]*thetaI[0]+rhoI[0]*(uXI[0]*uXI[0]+uYI[0]*uYI[0]+uZI[0]*uZI[0])-Pyy[0]-Pzz[0];   
//                 
//             }
//             
//             getGradThermalSinglePoint(lbModel, rhoI[0], uXI[0],uYI[0],uZI[0],thetaI[0],Pxx[0],Pyy[0],Pzz[0],Pxy[0],Pyz[0],Pzx[0],grad);
//             
//             copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nB1-1,j,k);
//             copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nB1-2,j,k);
//             copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nB1-3,j,k);
//             copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nB1-4,j,k);
//             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-1,j,k);
//             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-2,j,k);
//             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-3,j,k);
//             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-4,j,k);
//             
//             //                         copyFromNode(lbModel,myGrid,VECT_LENGTH,myGrid.nB1-4,j,k);
//             //                         getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);  
//             
//             //             if(k==4 && j==4)
//             //             for(int i=0;i<4;i++)
//             //              std::cout<<rho[i]<<"  "<<uX[i]<<"  "<<uY[i]<<"  "<<uZ[i]<<"  "<<theta[i]<<"   "<<Pxx[i]<<"  "<<Pyy[i]<<"  "<<Pzz[i]<<"  "<<Pxy[i]<<"  "<<Pyz[i]<<"   "<<Pzx[i]<<std::endl;
//             
//             // for(int i=0;i<4;i++)
//             // {
//             // //              std::cout<<rho[i]<<std::endl;
//             // //              std::cout<<uX[i]/cs<<std::endl;
//             // //              std::cout<<uY[i]<<std::endl;
//             // //              std::cout<<uZ[i]<<std::endl;
//             //              std::cout<<theta[i]-lbModel.theta0<<std::endl;
//             // }
//             
//             
//             //             copyToNodeSinglePoint(lbModel,myGrid,myGrid.nB1-1,j,k,0);
//             //             copyToNodeSinglePoint(lbModel,myGrid,myGrid.nB1-2,j,k,0);
//             //             copyToNodeSinglePoint(lbModel,myGrid,myGrid.nB1-3,j,k,0);
//             //             copyToNodeSinglePoint(lbModel,myGrid,myGrid.nB1-4,j,k,0);
//             
//             //         getHydroMomentsFromNode(lbModel,myGrid,myGrid.nB1-4,j,k,rho,uX,uY,uZ,theta);
//             //        std::cout<<rho[0]<<"   "<<rho[1]<<"   "<<rho[2]<<"   "<<rho[3]<<std::endl;
//             //        std::cout<<uX[0]<<"   "<<uX[1]<<"   "<<uX[2]<<"   "<<uX[3]<<std::endl;
//             //        std::cout<<uY[0]<<"   "<<uY[1]<<"   "<<uY[2]<<"   "<<uY[3]<<std::endl;
//             //        std::cout<<uZ[0]<<"   "<<uZ[1]<<"   "<<uZ[2]<<"   "<<uZ[3]<<std::endl;
//             //        std::cout<<theta[0]<<"   "<<theta[1]<<"   "<<theta[2]<<"   "<<theta[3]<<std::endl;
//             
//         }
//     }
//     
//     //     for (int k = myGrid.nB3; k <= myGrid.nE3; k++) 
//     //     {
//     //         for (int j =  myGrid.nB2; j <= myGrid.nE2; j++) 
//     //         {2
//     //             //         copyFromCell(lbModel,myGrid,VECT_LENGTH,myGrid.nB1,j,k);
//     // //             for(int i=0;i<VECT_LENGTH;i++)
//     // //             {
//     // //                 thetaI[i] = lbModel.theta0;            
//     // //                 rhoI  [i] = 1.0;
//     // //                 uXI   [i] = u_inlet;
//     // //                 uYI   [i] = 0.0;
//     // //                 uZI   [i] = 0.0;
//     // //             }
//     // //             getFEq(lbModel,VECT_LENGTH,rhoI,uXI,uYI,uZI,thetaI);      
//     //              copyFromNode(lbModel,myGrid,VECT_LENGTH,myGrid.nB1,j,k);
//     //              getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);  
//     //      
//     //             for(int i=0;i<VECT_LENGTH;i++)
//     //             {
//     //                 thetaI[i] = lbModel.theta0;           
//     //                 rhoI  [i] = 1.0;
//     //                 uXI   [i] = u_inlet;
//     //                 uYI   [i] = 0.0;
//     //                 uZI   [i] = 0.0;
//     //             }
//     //             getGradThermalSinglePoint(lbModel, rhoI[3], uXI[3],uYI[3],uZI[3],thetaI[3],Pxx[3],Pyy[3],Pzz[3],Pxy[3],Pyz[3],Pzx[3],grad);
//     //             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-1,j,k);
//     //             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nB1-2,j,k);
//     // //             copyToCellSinglePoint(lbModel,myGrid,myGrid.nB1-1,j,k,0);
//     // //             copyToCellSinglePoint(lbModel,myGrid,myGrid.nB1-2,j,k,0);
//     // //             copyToCellSinglePoint(lbModel,myGrid,myGrid.nB1-3,j,k,0);
//     // //             copyToCellSinglePoint(lbModel,myGrid,myGrid.nB1-4,j,k,0);
//     //             /*        getHydroMomentsFromCell(lbModel,myGrid,myGrid.nB1-4,j,k,rho,uX,uY,uZ,theta);
//     //              *       std::cout<<rho[0]<<"   "<<rho[1]<<"   "<<rho[2]<<"   "<<rho[3]<<std::endl;
//     //              *       std::cout<<uX[0]<<"   "<<uX[1]<<"   "<<uX[2]<<"   "<<uX[3]<<std::endl;
//     //              *       std::cout<<uY[0]<<"   "<<uY[1]<<"   "<<uY[2]<<"   "<<uY[3]<<std::endl;
//     //              *       std::cout<<uZ[0]<<"   "<<uZ[1]<<"   "<<uZ[2]<<"   "<<uZ[3]<<std::endl;
//     //              *       std::cout<<theta[0]<<"   "<<theta[1]<<"   "<<theta[2]<<"   "<<theta[3]<<std::endl;  */      
//     //         }
//     //     }
//     
// } 
// 
// template <int N,int numblock, typename dataType1>
// void X1EndOutlet(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH)
// {
//     //Outlet
//     dataType1 *rho,*uX,*uY,*uZ,*theta,*Pxx,*Pyy,*Pzz,*Pxy,*Pyz,*Pzx; 
//     dataType1 *rhoI,*uXI,*uYI,*uZI,*thetaI,*PxxI,*PyyI,*PzzI,*PxyI,*PyzI,*PzxI; 
//     rho   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              rhoI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     uX    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              uXI    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     uY    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              uYI    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     uZ    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              uZI    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     theta = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              thetaI = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32); 
//     Pxx   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              PxxI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     Pyy   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              PyyI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     Pzz   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              PzzI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     Pxy   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              PxyI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);       
//     Pyz   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              PyzI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32); 
//     Pzx   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);              PzxI   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);   
//     
//     dataType1 grad[lbModel.dvN];       
//     
//     for (int k =0; k <  myGrid.n3; k++) 
//     {
//         for (int j =0; j <  myGrid.n2; j++) 
//         {           
//  
//             //             getHydroMomentsFromNode(lbModel,myGrid,myGrid.nE1-3,j,k,rhoI,uXI,uYI,uZI,thetaI);
//             //             getFEq(lbModel,VECT_LENGTH,rhoI,uXI,uYI,uZI,thetaI);  
//             //             copyToNodeSinglePoint(lbModel,myGrid,myGrid.nE1+1,j,k,3);
//             //             copyToNodeSinglePoint(lbModel,myGrid,myGrid.nE1+2,j,k,3);
//             //             copyToNodeSinglePoint(lbModel,myGrid,myGrid.nE1+3,j,k,3);
//             //             copyToNodeSinglePoint(lbModel,myGrid,myGrid.nE1+4,j,k,3);
//             
//             
//             copyFromCell(lbModel,myGrid,VECT_LENGTH,myGrid.nE1,j,k);
//             getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx); 
//  /*           std::cout<<"Before: "<<std::endl;
//             for(int i=0;i<4;i++)
//              std::cout<<rho[i]<<"  "<<uX[i]<<"  "<<uY[i]<<"  "<<uZ[i]<<"  "<<theta[i]<<"   "<<Pxx[i]<<"  "<<Pyy[i]<<"  "<<Pzz[i]<<"  "<<Pxy[i]<<"  "<<Pyz[i]<<"   "<<Pzx[i]<<std::endl;
//  */           
//             //              PxxI  [i] = Pxx[i] ;//- rho[i]*theta[i] - rho[i]*uX[i]*uX[i] + lbModel.theta0 + u_inlet*u_inlet;// + (u_inlet - uX[i])*cs ;
//             //  if(k==4 && j==4)
//             // for(int i=0;i<4;i++)
//             //               std::cout<<rho[i]<<"  "<<uX[i]<<"  "<<uY[i]<<"  "<<uZ[i]<<"  "<<theta[i]<<"   "<<Pxx[i]<<"  "<<Pyy[i]<<"  "<<Pzz[i]<<"  "<<Pxy[i]<<"  "<<Pyz[i]<<"   "<<Pzx[i]<<std::endl;
//             
//             getGradThermalSinglePoint(lbModel, rho[0], uX[0],uY[0],uZ[0],theta[0],Pxx[0],Pyy[0],Pzz[0],Pxy[0],Pyz[0],Pzx[0],grad);
// //             
//             copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+1,j,k);
//             copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+2,j,k);
//             copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+3,j,k);
//             copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+4,j,k);
// //    getGradThermalSinglePoint(lbModel, rho[0], uX[0],uY[0],uZ[0],theta[0],Pxx[0],Pyy[0],Pzz[0],Pxy[0],Pyz[0],Pzx[0],grad);
// //    copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+1,j,k);
// // 
// //    getGradThermalSinglePoint(lbModel, rho[1], uX[1],uY[1],uZ[1],theta[1],Pxx[1],Pyy[1],Pzz[1],Pxy[1],Pyz[1],Pzx[1],grad);
// //    copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+2,j,k);
// // 
// //    getGradThermalSinglePoint(lbModel, rho[2], uX[2],uY[2],uZ[2],theta[2],Pxx[2],Pyy[2],Pzz[2],Pxy[2],Pyz[2],Pzx[2],grad);
// //    copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+3,j,k);
// // 
// //    getGradThermalSinglePoint(lbModel, rho[3], uX[3],uY[3],uZ[3],theta[3],Pxx[3],Pyy[3],Pzz[3],Pxy[3],Pyz[3],Pzx[3],grad);
// //    copyFromArrayToNode(lbModel,grad, myGrid,myGrid.nE1+4,j,k);
// //             
//             
// //             copyFromNode(lbModel,myGrid,VECT_LENGTH,myGrid.nE1+1,j,k);
// //             getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx); 
// //             std::cout<<"After: "<<std::endl;
// //             for(int i=0;i<4;i++)
// //              std::cout<<rho[i]<<"  "<<uX[i]<<"  "<<uY[i]<<"  "<<uZ[i]<<"  "<<theta[i]<<"   "<<Pxx[i]<<"  "<<Pyy[i]<<"  "<<Pzz[i]<<"  "<<Pxy[i]<<"  "<<Pyz[i]<<"   "<<Pzx[i]<<std::endl;
// 
// 
//             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+1,j,k);
//             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+2,j,k);
//             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+3,j,k);
//             copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+4,j,k);              
//         }
//     }
//     
//     //     for (int k = myGrid.ndB3; k <= myGrid.ndE3; k++) 
//     //     {
//     //         for (int j =  myGrid.ndB2; j <= myGrid.ndE2; j++) 
//     //         {
//     // //             //Cell
//     // //             getHydroMomentsFromCell(lbModel,myGrid,myGrid.nE1-3,j,k,rhoI,uXI,uYI,uZI,thetaI);
//     // //             getFEq(lbModel,VECT_LENGTH,rhoI,uXI,uYI,uZI,thetaI);      
//     // //             copyToCellSinglePoint(lbModel,myGrid,myGrid.nE1+1,j,k,3);
//     // //             copyToCellSinglePoint(lbModel,myGrid,myGrid.nE1+2,j,k,3);
//     // //             copyToCellSinglePoint(lbModel,myGrid,myGrid.nE1+3,j,k,3);
//     // //             copyToCellSinglePoint(lbModel,myGrid,myGrid.nE1+4,j,k,3);
//     //              copyFromCell(lbModel,myGrid,VECT_LENGTH,myGrid.nE1-3,j,k);
//     //              getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);  
//     //              getGradThermalSinglePoint(lbModel, rho[3], uX[3],uY[3],uZ[3],theta[3],Pxx[3],Pyy[3],Pzz[3],Pxy[3],Pyz[3],Pzx[3],grad);
//     //              copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+1,j,k);
//     //              copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+2,j,k);
//     //              copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+3,j,k);
//     //              copyFromArrayToCell(lbModel,grad, myGrid,myGrid.nE1+4,j,k);            
//     //         }
//     //     }
//     
//     
// }   

// Explicit declarations
template void X1BeginInlet<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &, double , double ,int ,int );
template void X1EndOutlet<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,int );

