#include"collide41.h"

void timestamp ( void )
{
 # define TIME_SIZE 40
 static char time_buffer[TIME_SIZE] ;
 const struct tm *tm;
 size_t len;
 time_t now;
 now = time ( NULL ) ;
 tm = localtime ( &now ) ;
 len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm ) ;
 printf ( "%s\n", time_buffer ) ;
 return;
 # undef TIME_SIZE
}

template <int N,int numblock, typename dataType1>
void collide(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 twoBeta,dataType1 oneMinustwoBeta)
{
    dataType1 *rho,*uX,*uY,*uZ,*theta;
    rho   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uX    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uY    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uZ    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    theta = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);

    dataType1 *pointer1,*pointer2;

    SIMD_REG _temp1,_temp2,_temp3,_temp4;
    _temp3 = SET1_PD(oneMinustwoBeta);    // oneMinustwoBeta
    _temp4 = SET1_PD(twoBeta);            // twoBeta


    for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1=i1+4)
            {
                getHydroMomentsFromNode(lbModel,myGrid,i1,i2,i3,rho,uX,uY,uZ,theta);

                pointer1 = &rho[0];
                _temp1 = LOAD_PD(pointer1);
                _temp1 = MUL_PD(_temp1,_temp4);                // Multiplying rho by two Beta
                STORE_PD(rho,_temp1);                          // write it back to rho
                getFEqSIMD(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);


                for(int i=0;i<VECT_LENGTH;i++)
                    myGrid(i1+i,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) =  oneMinustwoBeta * myGrid(i1+i,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) +  lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO];

                /////////
                // SC1 //
                /////////
                //SC1 - i
                pointer1 = &lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1);

                _temp1 = LOAD_PD(pointer1);           // fEq
                _temp2 = LOAD_PD(pointer2);           // G1
                _temp2 = MUL_PD(_temp2,_temp3);         // oneMinustwoBeta*SC1
                _temp1 = ADD_PD(_temp1,_temp2);         // oneMinustwoBeta*SC1 + fEq
                STORE_PD(pointer2,_temp1);            // write it back to SC1

                //SC1 - i+1
                pointer1 = &lbModel.fTemp1[1][lbModel.G1_DV_ZERO_ZERO_P1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //SC1 - i+2
                pointer1  = &lbModel.fTemp1[2][lbModel.G1_DV_ZERO_ZERO_P1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //SC1 - i+3
                pointer1  = &lbModel.fTemp1[3][lbModel.G1_DV_ZERO_ZERO_P1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);


                /////////
                // SC2 //
                /////////
                //SC2 - i
                pointer1 = &lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1);

                _temp1 = LOAD_PD(pointer1);           // fEq
                _temp2 = LOAD_PD(pointer2);           // G1
                _temp2 = MUL_PD(_temp2,_temp3);         // oneMinustwoBeta*SC1
                _temp1 = ADD_PD(_temp1,_temp2);         // oneMinustwoBeta*SC1 + fEq
                STORE_PD(pointer2,_temp1);            // write it back to SC1

                //SC2 - i+1
                pointer1 = &lbModel.fTemp2[1][lbModel.G2_DV_ZERO_ZERO_M1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //SC2 - i+2
                pointer1  = &lbModel.fTemp2[2][lbModel.G2_DV_ZERO_ZERO_M1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //SC2 - i+3
                pointer1  = &lbModel.fTemp2[3][lbModel.G2_DV_ZERO_ZERO_M1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);



                /////////
                // G3 //
                /////////
                //G3 - i
                pointer1 = &lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1);

                _temp1 = LOAD_PD(pointer1);           // fEq
                _temp2 = LOAD_PD(pointer2);           // G1
                _temp2 = MUL_PD(_temp2,_temp3);         // oneMinustwoBeta*SC1
                _temp1 = ADD_PD(_temp1,_temp2);         // oneMinustwoBeta*SC1 + fEq
                STORE_PD(pointer2,_temp1);            // write it back to SC1

                //G3 - i+1
                pointer1 = &lbModel.fTemp3[1][lbModel.G3_DV_ZERO_P1_P1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G3 - i+2
                pointer1  = &lbModel.fTemp3[2][lbModel.G3_DV_ZERO_P1_P1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G3 - i+3
                pointer1  = &lbModel.fTemp3[3][lbModel.G3_DV_ZERO_P1_P1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);


                /////////
                // G4 //
                /////////
                //G4 - i
                pointer1 = &lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1);

                _temp1 = LOAD_PD(pointer1);           // fEq
                _temp2 = LOAD_PD(pointer2);           // G1
                _temp2 = MUL_PD(_temp2,_temp3);         // oneMinustwoBeta*SC1
                _temp1 = ADD_PD(_temp1,_temp2);         // oneMinustwoBeta*SC1 + fEq
                STORE_PD(pointer2,_temp1);            // write it back to SC1

                //G4 - i+1
                pointer1 = &lbModel.fTemp4[1][lbModel.G4_DV_ZERO_M1_M1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G4 - i+2
                pointer1  = &lbModel.fTemp4[2][lbModel.G4_DV_ZERO_M1_M1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G4 - i+3
                pointer1  = &lbModel.fTemp4[3][lbModel.G4_DV_ZERO_M1_M1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);



                /////////
                // G5 //
                /////////
                //G5 - i
                pointer1 = &lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO);

                _temp1 = LOAD_PD(pointer1);           // fEq
                _temp2 = LOAD_PD(pointer2);           // G1
                _temp2 = MUL_PD(_temp2,_temp3);         // oneMinustwoBeta*SC1
                _temp1 = ADD_PD(_temp1,_temp2);         // oneMinustwoBeta*SC1 + fEq
                STORE_PD(pointer2,_temp1);            // write it back to SC1

                //G5 - i+1
                pointer1 = &lbModel.fTemp5[1][lbModel.G5_DV_P1_P1_ZERO];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G5 - i+2
                pointer1  = &lbModel.fTemp5[2][lbModel.G5_DV_P1_P1_ZERO];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G5 - i+3
                pointer1  = &lbModel.fTemp5[3][lbModel.G5_DV_P1_P1_ZERO];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);





                /////////
                // G6 //
                /////////
                //G6 - i
                pointer1 = &lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO);

                _temp1 = LOAD_PD(pointer1);           // fEq
                _temp2 = LOAD_PD(pointer2);           // G1
                _temp2 = MUL_PD(_temp2,_temp3);       // oneMinustwoBeta*SC1
                _temp1 = ADD_PD(_temp1,_temp2);       // oneMinustwoBeta*SC1 + fEq
                STORE_PD(pointer2,_temp1);            // write it back to SC1

                //G6 - i+1
                pointer1 = &lbModel.fTemp6[1][lbModel.G6_DV_M1_M1_ZERO];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G6 - i+2
                pointer1  = &lbModel.fTemp6[2][lbModel.G6_DV_M1_M1_ZERO];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G6 - i+3
                pointer1  = &lbModel.fTemp6[3][lbModel.G6_DV_M1_M1_ZERO];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);






                /////////
                // G7 //
                /////////
                //G7 - i
                pointer1 = &lbModel.fTemp7[0][lbModel.G7_DV_P_P_P];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);

                _temp1 = LOAD_PD(pointer1);           // fEq
                _temp2 = LOAD_PD(pointer2);           // G1
                _temp2 = MUL_PD(_temp2,_temp3);       // oneMinustwoBeta*SC1
                _temp1 = ADD_PD(_temp1,_temp2);       // oneMinustwoBeta*SC1 + fEq
                STORE_PD(pointer2,_temp1);            // write it back to SC1

                //G7 - i+1
                pointer1 = &lbModel.fTemp7[1][lbModel.G7_DV_P_P_P];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G7 - i+2
                pointer1  = &lbModel.fTemp7[2][lbModel.G7_DV_P_P_P];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G7 - i+3
                pointer1  = &lbModel.fTemp7[3][lbModel.G7_DV_P_P_P];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);


                /////////
                // G8 //
                /////////
                //G8 - i
                pointer1 = &lbModel.fTemp8[0][lbModel.G8_DV_M_M_M];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);

                _temp1 = LOAD_PD(pointer1);           // fEq
                _temp2 = LOAD_PD(pointer2);           // G1
                _temp2 = MUL_PD(_temp2,_temp3);       // oneMinustwoBeta*SC1
                _temp1 = ADD_PD(_temp1,_temp2);       // oneMinustwoBeta*SC1 + fEq
                STORE_PD(pointer2,_temp1);            // write it back to SC1

                //G8 - i+1
                pointer1 = &lbModel.fTemp8[1][lbModel.G8_DV_M_M_M];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G8 - i+2
                pointer1  = &lbModel.fTemp8[2][lbModel.G8_DV_M_M_M];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G8 - i+3
                pointer1  = &lbModel.fTemp8[3][lbModel.G8_DV_M_M_M];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);




                /////////
                // G9 //
                /////////
                //G9 - i
                pointer1 = &lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_P1_P1);

                _temp1 = LOAD_PD(pointer1);           // fEq
                _temp2 = LOAD_PD(pointer2);           // G1
                _temp2 = MUL_PD(_temp2,_temp3);       // oneMinustwoBeta*SC1
                _temp1 = ADD_PD(_temp1,_temp2);       // oneMinustwoBeta*SC1 + fEq
                STORE_PD(pointer2,_temp1);            // write it back to SC1

                //G9 - i+1
                pointer1 = &lbModel.fTemp9[1][lbModel.G9_DV_P1_P1_P1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G9 - i+2
                pointer1  = &lbModel.fTemp9[2][lbModel.G9_DV_P1_P1_P1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G9 - i+3
                pointer1  = &lbModel.fTemp9[3][lbModel.G9_DV_P1_P1_P1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);



                /////////
                // G10 //
                /////////
                //G10 - i
                pointer1 = &lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1);

                _temp1 = LOAD_PD(pointer1);           // fEq
                _temp2 = LOAD_PD(pointer2);           // G1
                _temp2 = MUL_PD(_temp2,_temp3);       // oneMinustwoBeta*SC1
                _temp1 = ADD_PD(_temp1,_temp2);       // oneMinustwoBeta*SC1 + fEq
                STORE_PD(pointer2,_temp1);            // write it back to SC1

                //G10 - i+1
                pointer1 = &lbModel.fTemp10[1][lbModel.G10_DV_M1_M1_M1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G10 - i+2
                pointer1  = &lbModel.fTemp10[2][lbModel.G10_DV_M1_M1_M1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);

                //G10 - i+3
                pointer1  = &lbModel.fTemp10[3][lbModel.G10_DV_M1_M1_M1];
                pointer2 += 4;
                _temp1 = LOAD_PD(pointer1);
                _temp2 = LOAD_PD(pointer2);
                _temp2 = MUL_PD(_temp2,_temp3);
                _temp1 = ADD_PD(_temp1,_temp2);
                STORE_PD(pointer2,_temp1);
            }


            for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
                for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
                    for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1=i1+4)
                    {
                        getHydroMomentsFromCell(lbModel,myGrid,i1,i2,i3,rho,uX,uY,uZ,theta);

                        pointer1 = &rho[0];
                        _temp1 = LOAD_PD(pointer1);
                        _temp1 = MUL_PD(_temp1,_temp4);                // Multiplying rho by two Beta
                        STORE_PD(rho,_temp1);                          // write it back to rho  */

                        getFEqSIMD(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);

                        for(int i=0;i<VECT_LENGTH;i++)
                            myGrid(i1+i,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) =  oneMinustwoBeta * myGrid(i1+i,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) +  lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO];

                        /////////
                        // SC1 //
                        /////////
                        //SC1 - i
                        pointer1 = &lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1];
                        pointer2 = &myGrid(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1);

                        _temp1 = LOAD_PD(pointer1);           // fEq
                        _temp2 = LOAD_PD(pointer2);           // G1
                        _temp2 = MUL_PD(_temp2,_temp3);         // oneMinustwoBeta*SC1
                        _temp1 = ADD_PD(_temp1,_temp2);         // oneMinustwoBeta*SC1 + fEq
                        STORE_PD(pointer2,_temp1);            // write it back to SC1

                        //SC1 - i+1
                        pointer1 = &lbModel.fTemp1[1][lbModel.G1_DV_ZERO_ZERO_P1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //SC1 - i+2
                        pointer1  = &lbModel.fTemp1[2][lbModel.G1_DV_ZERO_ZERO_P1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //SC1 - i+3
                        pointer1  = &lbModel.fTemp1[3][lbModel.G1_DV_ZERO_ZERO_P1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);


                        /////////
                        // SC2 //
                        /////////
                        //SC2 - i
                        pointer1 = &lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];
                        pointer2 = &myGrid(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1);

                        _temp1 = LOAD_PD(pointer1);           // fEq
                        _temp2 = LOAD_PD(pointer2);           // G1
                        _temp2 = MUL_PD(_temp2,_temp3);         // oneMinustwoBeta*SC1
                        _temp1 = ADD_PD(_temp1,_temp2);         // oneMinustwoBeta*SC1 + fEq
                        STORE_PD(pointer2,_temp1);            // write it back to SC1

                        //SC2 - i+1
                        pointer1 = &lbModel.fTemp2[1][lbModel.G2_DV_ZERO_ZERO_M1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //SC2 - i+2
                        pointer1  = &lbModel.fTemp2[2][lbModel.G2_DV_ZERO_ZERO_M1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //SC2 - i+3
                        pointer1  = &lbModel.fTemp2[3][lbModel.G2_DV_ZERO_ZERO_M1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);



                        /////////
                        // G3 //
                        /////////
                        //G3 - i
                        pointer1 = &lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1];
                        pointer2 = &myGrid(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1);

                        _temp1 = LOAD_PD(pointer1);           // fEq
                        _temp2 = LOAD_PD(pointer2);           // G1
                        _temp2 = MUL_PD(_temp2,_temp3);         // oneMinustwoBeta*SC1
                        _temp1 = ADD_PD(_temp1,_temp2);         // oneMinustwoBeta*SC1 + fEq
                        STORE_PD(pointer2,_temp1);            // write it back to SC1

                        //G3 - i+1
                        pointer1 = &lbModel.fTemp3[1][lbModel.G3_DV_ZERO_P1_P1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G3 - i+2
                        pointer1  = &lbModel.fTemp3[2][lbModel.G3_DV_ZERO_P1_P1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G3 - i+3
                        pointer1  = &lbModel.fTemp3[3][lbModel.G3_DV_ZERO_P1_P1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);


                        /////////
                        // G4 //
                        /////////
                        //G4 - i
                        pointer1 = &lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1];
                        pointer2 = &myGrid(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1);

                        _temp1 = LOAD_PD(pointer1);           // fEq
                        _temp2 = LOAD_PD(pointer2);           // G1
                        _temp2 = MUL_PD(_temp2,_temp3);         // oneMinustwoBeta*SC1
                        _temp1 = ADD_PD(_temp1,_temp2);         // oneMinustwoBeta*SC1 + fEq
                        STORE_PD(pointer2,_temp1);            // write it back to SC1

                        //G4 - i+1
                        pointer1 = &lbModel.fTemp4[1][lbModel.G4_DV_ZERO_M1_M1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G4 - i+2
                        pointer1  = &lbModel.fTemp4[2][lbModel.G4_DV_ZERO_M1_M1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G4 - i+3
                        pointer1  = &lbModel.fTemp4[3][lbModel.G4_DV_ZERO_M1_M1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);



                        /////////
                        // G5 //
                        /////////
                        //G5 - i
                        pointer1 = &lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO];
                        pointer2 = &myGrid(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO);

                        _temp1 = LOAD_PD(pointer1);           // fEq
                        _temp2 = LOAD_PD(pointer2);           // G1
                        _temp2 = MUL_PD(_temp2,_temp3);         // oneMinustwoBeta*SC1
                        _temp1 = ADD_PD(_temp1,_temp2);         // oneMinustwoBeta*SC1 + fEq
                        STORE_PD(pointer2,_temp1);            // write it back to SC1

                        //G5 - i+1
                        pointer1 = &lbModel.fTemp5[1][lbModel.G5_DV_P1_P1_ZERO];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G5 - i+2
                        pointer1  = &lbModel.fTemp5[2][lbModel.G5_DV_P1_P1_ZERO];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G5 - i+3
                        pointer1  = &lbModel.fTemp5[3][lbModel.G5_DV_P1_P1_ZERO];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);





                        /////////
                        // G6 //
                        /////////
                        //G6 - i
                        pointer1 = &lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO];
                        pointer2 = &myGrid(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO);

                        _temp1 = LOAD_PD(pointer1);           // fEq
                        _temp2 = LOAD_PD(pointer2);           // G1
                        _temp2 = MUL_PD(_temp2,_temp3);       // oneMinustwoBeta*SC1
                        _temp1 = ADD_PD(_temp1,_temp2);       // oneMinustwoBeta*SC1 + fEq
                        STORE_PD(pointer2,_temp1);            // write it back to SC1

                        //G6 - i+1
                        pointer1 = &lbModel.fTemp6[1][lbModel.G6_DV_M1_M1_ZERO];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G6 - i+2
                        pointer1  = &lbModel.fTemp6[2][lbModel.G6_DV_M1_M1_ZERO];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G6 - i+3
                        pointer1  = &lbModel.fTemp6[3][lbModel.G6_DV_M1_M1_ZERO];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);






                        /////////
                        // G7 //
                        /////////
                        //G7 - i
                        pointer1 = &lbModel.fTemp7[0][lbModel.G7_DV_P_P_P];
                        pointer2 = &myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);

                        _temp1 = LOAD_PD(pointer1);           // fEq
                        _temp2 = LOAD_PD(pointer2);           // G1
                        _temp2 = MUL_PD(_temp2,_temp3);       // oneMinustwoBeta*SC1
                        _temp1 = ADD_PD(_temp1,_temp2);       // oneMinustwoBeta*SC1 + fEq
                        STORE_PD(pointer2,_temp1);            // write it back to SC1

                        //G7 - i+1
                        pointer1 = &lbModel.fTemp7[1][lbModel.G7_DV_P_P_P];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G7 - i+2
                        pointer1  = &lbModel.fTemp7[2][lbModel.G7_DV_P_P_P];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G7 - i+3
                        pointer1  = &lbModel.fTemp7[3][lbModel.G7_DV_P_P_P];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);


                        /////////
                        // G8 //
                        /////////
                        //G8 - i
                        pointer1 = &lbModel.fTemp8[0][lbModel.G8_DV_M_M_M];
                        pointer2 = &myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);

                        _temp1 = LOAD_PD(pointer1);           // fEq
                        _temp2 = LOAD_PD(pointer2);           // G1
                        _temp2 = MUL_PD(_temp2,_temp3);       // oneMinustwoBeta*SC1
                        _temp1 = ADD_PD(_temp1,_temp2);       // oneMinustwoBeta*SC1 + fEq
                        STORE_PD(pointer2,_temp1);            // write it back to SC1

                        //G8 - i+1
                        pointer1 = &lbModel.fTemp8[1][lbModel.G8_DV_M_M_M];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G8 - i+2
                        pointer1  = &lbModel.fTemp8[2][lbModel.G8_DV_M_M_M];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G8 - i+3
                        pointer1  = &lbModel.fTemp8[3][lbModel.G8_DV_M_M_M];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);




                        /////////
                        // G9 //
                        /////////
                        //G9 - i
                        pointer1 = &lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1];
                        pointer2 = &myGrid(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_P1_P1);

                        _temp1 = LOAD_PD(pointer1);           // fEq
                        _temp2 = LOAD_PD(pointer2);           // G1
                        _temp2 = MUL_PD(_temp2,_temp3);       // oneMinustwoBeta*SC1
                        _temp1 = ADD_PD(_temp1,_temp2);       // oneMinustwoBeta*SC1 + fEq
                        STORE_PD(pointer2,_temp1);            // write it back to SC1

                        //G9 - i+1
                        pointer1 = &lbModel.fTemp9[1][lbModel.G9_DV_P1_P1_P1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G9 - i+2
                        pointer1  = &lbModel.fTemp9[2][lbModel.G9_DV_P1_P1_P1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G9 - i+3
                        pointer1  = &lbModel.fTemp9[3][lbModel.G9_DV_P1_P1_P1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);



                        /////////
                        // G10 //
                        /////////
                        //G10 - i
                        pointer1 = &lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
                        pointer2 = &myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1);

                        _temp1 = LOAD_PD(pointer1);           // fEq
                        _temp2 = LOAD_PD(pointer2);           // G1
                        _temp2 = MUL_PD(_temp2,_temp3);       // oneMinustwoBeta*SC1
                        _temp1 = ADD_PD(_temp1,_temp2);       // oneMinustwoBeta*SC1 + fEq
                        STORE_PD(pointer2,_temp1);            // write it back to SC1

                        //G10 - i+1
                        pointer1 = &lbModel.fTemp10[1][lbModel.G10_DV_M1_M1_M1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G10 - i+2
                        pointer1  = &lbModel.fTemp10[2][lbModel.G10_DV_M1_M1_M1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                        //G10 - i+3
                        pointer1  = &lbModel.fTemp10[3][lbModel.G10_DV_M1_M1_M1];
                        pointer2 += 4;
                        _temp1 = LOAD_PD(pointer1);
                        _temp2 = LOAD_PD(pointer2);
                        _temp2 = MUL_PD(_temp2,_temp3);
                        _temp1 = ADD_PD(_temp1,_temp2);
                        STORE_PD(pointer2,_temp1);

                    }



}

template <typename dataType1>
void getMomentsWithForce(lbmRD3Q41<dataType1> &lbModel,dataType1 *f, double F1, double F2, double F3)
{
    lbModel.rho = 0.0;  lbModel.uX = 0.0;  lbModel.uY = 0.0;  lbModel.uZ = 0.0;  lbModel.theta = 0.0;
    double rhoSC1= 0.0, rhoSC2= 0.0, rhoFCC= 0.0, rhoBCC= 0.0, rhoBCC1= 0.0;


    dataType1 csq(0.0),dot(0.0);

    lbModel.rho = f[lbModel.DV_ZERO_ZERO_ZERO];

    rhoSC1 = f[lbModel.DV_ZERO_ZERO_P1] +  f[lbModel.DV_ZERO_P1_ZERO]  + f[lbModel.DV_P1_ZERO_ZERO] + f[lbModel.DV_ZERO_ZERO_M1] + f[lbModel.DV_ZERO_M1_ZERO] + f[lbModel.DV_M1_ZERO_ZERO] ;

    rhoSC2 = f[lbModel.DV_ZERO_ZERO_P2] +  f[lbModel.DV_ZERO_P2_ZERO]  + f[lbModel.DV_P2_ZERO_ZERO] + f[lbModel.DV_ZERO_ZERO_M2] + f[lbModel.DV_ZERO_M2_ZERO] + f[lbModel.DV_M2_ZERO_ZERO] ;

    rhoFCC = f[lbModel.DV_ZERO_P1_P1 ]  + f[lbModel.DV_ZERO_M1_P1 ] + f [lbModel.DV_ZERO_M1_M1 ] + f[lbModel.DV_ZERO_P1_M1 ] +  f[lbModel.DV_M1_ZERO_M1 ]  + f[lbModel.DV_M1_ZERO_P1 ] + f [lbModel.DV_P1_ZERO_M1 ] + f[lbModel.DV_P1_ZERO_P1 ] +  f[lbModel.DV_M1_P1_ZERO ]  + f[lbModel.DV_P1_P1_ZERO ] + f [lbModel.DV_M1_M1_ZERO ] + f[lbModel.DV_P1_M1_ZERO ] ;

    rhoBCC   = f[lbModel.DV_P_P_P ] + f[lbModel.DV_M_P_P ] + f[lbModel.DV_M_M_P ] + f[lbModel.DV_P_M_P ] + f[lbModel.DV_M_M_M ] + f[lbModel.DV_P_M_M ] + f[lbModel.DV_P_P_M ] + f[lbModel.DV_M_P_M ] ;

    rhoBCC1 = f[lbModel.DV_P1_P1_P1] + f[lbModel.DV_M1_P1_P1] + f[lbModel.DV_M1_M1_P1]+ f[lbModel.DV_P1_M1_P1]+ f[lbModel.DV_M1_M1_M1]+ f[lbModel.DV_P1_M1_M1]+ f[lbModel.DV_P1_P1_M1]+ f[lbModel.DV_M1_P1_M1] ;

    lbModel.rho +=   rhoSC1+ rhoSC2+ rhoFCC +rhoBCC+ rhoBCC1;

    lbModel.uX =  f[lbModel.DV_P1_ZERO_ZERO] - f[lbModel.DV_M1_ZERO_ZERO] ;
    lbModel.uY =  f[lbModel.DV_ZERO_P1_ZERO] - f[lbModel.DV_ZERO_M1_ZERO] ;
    lbModel.uZ =  f[lbModel.DV_ZERO_ZERO_P1] - f[lbModel.DV_ZERO_ZERO_M1] ;

    lbModel.uX += 2.0*( f[lbModel.DV_P2_ZERO_ZERO] - f[lbModel.DV_M2_ZERO_ZERO]) ;
    lbModel.uY += 2.0*( f[lbModel.DV_ZERO_P2_ZERO] - f[lbModel.DV_ZERO_M2_ZERO]) ;
    lbModel.uZ += 2.0*( f[lbModel.DV_ZERO_ZERO_P2] - f[lbModel.DV_ZERO_ZERO_M2]) ;

    lbModel.uX += f[lbModel.DV_P1_ZERO_M1] - f[lbModel.DV_M1_ZERO_M1] - f[lbModel.DV_M1_ZERO_P1] + f[lbModel.DV_P1_ZERO_P1] - f[lbModel.DV_M1_P1_ZERO] + f[lbModel.DV_P1_P1_ZERO] - f[lbModel.DV_M1_M1_ZERO] + f[lbModel.DV_P1_M1_ZERO ] ;

    lbModel.uY += f[lbModel.DV_ZERO_P1_P1] - f[lbModel.DV_ZERO_M1_P1] - f[lbModel.DV_ZERO_M1_M1] + f[lbModel.DV_ZERO_P1_M1] + f[lbModel.DV_M1_P1_ZERO] + f[lbModel.DV_P1_P1_ZERO] - f[lbModel.DV_M1_M1_ZERO] - f[lbModel.DV_P1_M1_ZERO ] ;

    lbModel.uZ += f[lbModel.DV_ZERO_P1_P1] + f[lbModel.DV_ZERO_M1_P1] - f[lbModel.DV_ZERO_M1_M1] - f[lbModel.DV_ZERO_P1_M1] - f[lbModel.DV_M1_ZERO_M1] + f[lbModel.DV_M1_ZERO_P1] - f[lbModel.DV_P1_ZERO_M1] + f[lbModel.DV_P1_ZERO_P1 ] ;

    lbModel.uX += 0.5*( f[lbModel.DV_P_P_P ] - f[lbModel.DV_M_P_P ] - f[lbModel.DV_M_M_P ] + f[lbModel.DV_P_M_P ] - f[lbModel.DV_M_M_M ] + f[lbModel.DV_P_M_M ] + f[lbModel.DV_P_P_M ] - f[lbModel.DV_M_P_M ]) ;
    lbModel.uY += 0.5*( f[lbModel.DV_P_P_P ] + f[lbModel.DV_M_P_P ] - f[lbModel.DV_M_M_P ] - f[lbModel.DV_P_M_P ] - f[lbModel.DV_M_M_M ] - f[lbModel.DV_P_M_M ] + f[lbModel.DV_P_P_M ] + f[lbModel.DV_M_P_M ]) ;
    lbModel.uZ += 0.5*( f[lbModel.DV_P_P_P ] + f[lbModel.DV_M_P_P ] + f[lbModel.DV_M_M_P ] + f[lbModel.DV_P_M_P ] - f[lbModel.DV_M_M_M ] - f[lbModel.DV_P_M_M ] - f[lbModel.DV_P_P_M ] - f[lbModel.DV_M_P_M ]) ;

    lbModel.uX +=  f[lbModel.DV_P1_P1_P1] - f[lbModel.DV_M1_P1_P1] - f[lbModel.DV_M1_M1_P1] + f[lbModel.DV_P1_M1_P1] - f[lbModel.DV_M1_M1_M1] + f[lbModel.DV_P1_M1_M1] + f[lbModel.DV_P1_P1_M1] - f[lbModel.DV_M1_P1_M1] ;
    lbModel.uY +=  f[lbModel.DV_P1_P1_P1] + f[lbModel.DV_M1_P1_P1] - f[lbModel.DV_M1_M1_P1] - f[lbModel.DV_P1_M1_P1] - f[lbModel.DV_M1_M1_M1] - f[lbModel.DV_P1_M1_M1] + f[lbModel.DV_P1_P1_M1] + f[lbModel.DV_M1_P1_M1] ;
    lbModel.uZ +=  f[lbModel.DV_P1_P1_P1] + f[lbModel.DV_M1_P1_P1] + f[lbModel.DV_M1_M1_P1] + f[lbModel.DV_P1_M1_P1] - f[lbModel.DV_M1_M1_M1] - f[lbModel.DV_P1_M1_M1] - f[lbModel.DV_P1_P1_M1] - f[lbModel.DV_M1_P1_M1] ;

    lbModel.theta  = rhoSC1 + rhoSC2* 4.0 + rhoFCC* 2.0 + rhoBCC* 0.75 + rhoBCC1 * 3.0 ;


    dataType1  oneByRho = 1.0/lbModel.rho;
    dataType1  oneByThree = 1.0/3.0;
    lbModel.uX    = lbModel.uX * oneByRho;
    lbModel.uY    = lbModel.uY * oneByRho;
    lbModel.uZ    = lbModel.uZ * oneByRho;

    lbModel.uX += 0.5*F1;
    lbModel.uY += 0.5*F2;
    lbModel.uZ += 0.5*F3;

    lbModel.theta = oneByThree * (lbModel.theta - lbModel.rho*(lbModel.uX*lbModel.uX + lbModel.uY*lbModel.uY + lbModel.uZ*lbModel.uZ));
    lbModel.theta = lbModel.theta * oneByRho;

    for (int dv=0;dv<lbModel.dvN;dv++)
    {
        csq = lbModel.cx[dv]*lbModel.cx[dv] + lbModel.cy[dv]*lbModel.cy[dv] + lbModel.cz[dv]*lbModel.cz[dv];
        dot = csq*(F1*lbModel.cx[dv]+F2*lbModel.cy[dv]+F3*lbModel.cz[dv]) ;
        lbModel.theta +=  (0.5*lbModel.wt[dv]*dot) ;
    }
}

template <int N,int numblock, typename dataType1>
void collideWithForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH, dataType1 twoBeta, dataType1 tau , dataType1 F1, dataType1 F2, dataType1 F3, dataType1 dt)
{
    dataType1 fTempF[4][lbModel.dvN] ;

    dataType1 rho[VECT_LENGTH] __attribute__ ((aligned(32)));
    dataType1 uX[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 uY[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 uZ[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 theta[VECT_LENGTH] __attribute__ ((aligned(32)));

    dataType1 twoBetaTau =  twoBeta*tau  ;
    dataType1 oneByTheta0 = 1.0/lbModel.theta0;

   for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1=i1+4)
            {
                getHydroMomentsFromNodeWithForce(lbModel,myGrid,VECT_LENGTH,i1,i2,i3,rho,uX,uY,uZ,theta,F1,F2,F3,dt);

                for(int index=0;index<VECT_LENGTH;index++)
                    theta[index] = lbModel.theta0;

                getFEq(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);

                for(int index=0;index<VECT_LENGTH;index++)
                {
                    fTempF[index][lbModel.DV_ZERO_ZERO_ZERO] = lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO]*(F1*lbModel.cx[lbModel.DV_ZERO_ZERO_ZERO] + F2*lbModel.cy[lbModel.DV_ZERO_ZERO_ZERO] + F3*lbModel.cz[lbModel.DV_ZERO_ZERO_ZERO])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_ZERO_ZERO_P1  ] = lbModel.wt[lbModel.DV_ZERO_ZERO_P1  ]*(F1*lbModel.cx[lbModel.DV_ZERO_ZERO_P1  ] + F2*lbModel.cy[lbModel.DV_ZERO_ZERO_P1  ] + F3*lbModel.cz[lbModel.DV_ZERO_ZERO_P1  ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_ZERO_ZERO_P2  ] = lbModel.wt[lbModel.DV_ZERO_ZERO_P2  ]*(F1*lbModel.cx[lbModel.DV_ZERO_ZERO_P2  ] + F2*lbModel.cy[lbModel.DV_ZERO_ZERO_P2  ] + F3*lbModel.cz[lbModel.DV_ZERO_ZERO_P2  ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_ZERO_P2_ZERO  ] = lbModel.wt[lbModel.DV_ZERO_P2_ZERO  ]*(F1*lbModel.cx[lbModel.DV_ZERO_P2_ZERO  ] + F2*lbModel.cy[lbModel.DV_ZERO_P2_ZERO  ] + F3*lbModel.cz[lbModel.DV_ZERO_P2_ZERO  ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P2_ZERO_ZERO  ] = lbModel.wt[lbModel.DV_P2_ZERO_ZERO  ]*(F1*lbModel.cx[lbModel.DV_P2_ZERO_ZERO  ] + F2*lbModel.cy[lbModel.DV_P2_ZERO_ZERO  ] + F3*lbModel.cz[lbModel.DV_P2_ZERO_ZERO  ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_ZERO_ZERO_M1  ] = lbModel.wt[lbModel.DV_ZERO_ZERO_M1  ]*(F1*lbModel.cx[lbModel.DV_ZERO_ZERO_M1  ] + F2*lbModel.cy[lbModel.DV_ZERO_ZERO_M1  ] + F3*lbModel.cz[lbModel.DV_ZERO_ZERO_M1  ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_ZERO_ZERO_M2  ] = lbModel.wt[lbModel.DV_ZERO_ZERO_M2  ]*(F1*lbModel.cx[lbModel.DV_ZERO_ZERO_M2  ] + F2*lbModel.cy[lbModel.DV_ZERO_ZERO_M2  ] + F3*lbModel.cz[lbModel.DV_ZERO_ZERO_M2  ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_ZERO_M2_ZERO  ] = lbModel.wt[lbModel.DV_ZERO_M2_ZERO  ]*(F1*lbModel.cx[lbModel.DV_ZERO_M2_ZERO  ] + F2*lbModel.cy[lbModel.DV_ZERO_M2_ZERO  ] + F3*lbModel.cz[lbModel.DV_ZERO_M2_ZERO  ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M2_ZERO_ZERO  ] = lbModel.wt[lbModel.DV_M2_ZERO_ZERO  ]*(F1*lbModel.cx[lbModel.DV_M2_ZERO_ZERO  ] + F2*lbModel.cy[lbModel.DV_M2_ZERO_ZERO  ] + F3*lbModel.cz[lbModel.DV_M2_ZERO_ZERO  ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_ZERO_P1_P1    ] = lbModel.wt[lbModel.DV_ZERO_P1_P1    ]*(F1*lbModel.cx[lbModel.DV_ZERO_P1_P1    ] + F2*lbModel.cy[lbModel.DV_ZERO_P1_P1    ] + F3*lbModel.cz[lbModel.DV_ZERO_P1_P1    ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_ZERO_M1_P1    ] = lbModel.wt[lbModel.DV_ZERO_M1_P1    ]*(F1*lbModel.cx[lbModel.DV_ZERO_M1_P1    ] + F2*lbModel.cy[lbModel.DV_ZERO_M1_P1    ] + F3*lbModel.cz[lbModel.DV_ZERO_M1_P1    ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P1_ZERO_P1    ] = lbModel.wt[lbModel.DV_P1_ZERO_P1    ]*(F1*lbModel.cx[lbModel.DV_P1_ZERO_P1    ] + F2*lbModel.cy[lbModel.DV_P1_ZERO_P1    ] + F3*lbModel.cz[lbModel.DV_P1_ZERO_P1    ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M1_ZERO_P1    ] = lbModel.wt[lbModel.DV_M1_ZERO_P1    ]*(F1*lbModel.cx[lbModel.DV_M1_ZERO_P1    ] + F2*lbModel.cy[lbModel.DV_M1_ZERO_P1    ] + F3*lbModel.cz[lbModel.DV_M1_ZERO_P1    ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_ZERO_M1_M1    ] = lbModel.wt[lbModel.DV_ZERO_M1_M1    ]*(F1*lbModel.cx[lbModel.DV_ZERO_M1_M1    ] + F2*lbModel.cy[lbModel.DV_ZERO_M1_M1    ] + F3*lbModel.cz[lbModel.DV_ZERO_M1_M1    ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_ZERO_P1_M1    ] = lbModel.wt[lbModel.DV_ZERO_P1_M1    ]*(F1*lbModel.cx[lbModel.DV_ZERO_P1_M1    ] + F2*lbModel.cy[lbModel.DV_ZERO_P1_M1    ] + F3*lbModel.cz[lbModel.DV_ZERO_P1_M1    ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M1_ZERO_M1    ] = lbModel.wt[lbModel.DV_M1_ZERO_M1    ]*(F1*lbModel.cx[lbModel.DV_M1_ZERO_M1    ] + F2*lbModel.cy[lbModel.DV_M1_ZERO_M1    ] + F3*lbModel.cz[lbModel.DV_M1_ZERO_M1    ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P1_ZERO_M1    ] = lbModel.wt[lbModel.DV_P1_ZERO_M1    ]*(F1*lbModel.cx[lbModel.DV_P1_ZERO_M1    ] + F2*lbModel.cy[lbModel.DV_P1_ZERO_M1    ] + F3*lbModel.cz[lbModel.DV_P1_ZERO_M1    ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P1_P1_ZERO    ] = lbModel.wt[lbModel.DV_P1_P1_ZERO    ]*(F1*lbModel.cx[lbModel.DV_P1_P1_ZERO    ] + F2*lbModel.cy[lbModel.DV_P1_P1_ZERO    ] + F3*lbModel.cz[lbModel.DV_P1_P1_ZERO    ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M1_P1_ZERO    ] = lbModel.wt[lbModel.DV_M1_P1_ZERO    ]*(F1*lbModel.cx[lbModel.DV_M1_P1_ZERO    ] + F2*lbModel.cy[lbModel.DV_M1_P1_ZERO    ] + F3*lbModel.cz[lbModel.DV_M1_P1_ZERO    ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_ZERO_P1_ZERO  ] = lbModel.wt[lbModel.DV_ZERO_P1_ZERO  ]*(F1*lbModel.cx[lbModel.DV_ZERO_P1_ZERO  ] + F2*lbModel.cy[lbModel.DV_ZERO_P1_ZERO  ] + F3*lbModel.cz[lbModel.DV_ZERO_P1_ZERO  ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P1_ZERO_ZERO  ] = lbModel.wt[lbModel.DV_P1_ZERO_ZERO  ]*(F1*lbModel.cx[lbModel.DV_P1_ZERO_ZERO  ] + F2*lbModel.cy[lbModel.DV_P1_ZERO_ZERO  ] + F3*lbModel.cz[lbModel.DV_P1_ZERO_ZERO  ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M1_M1_ZERO    ] = lbModel.wt[lbModel.DV_M1_M1_ZERO    ]*(F1*lbModel.cx[lbModel.DV_M1_M1_ZERO    ] + F2*lbModel.cy[lbModel.DV_M1_M1_ZERO    ] + F3*lbModel.cz[lbModel.DV_M1_M1_ZERO    ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P1_M1_ZERO    ] = lbModel.wt[lbModel.DV_P1_M1_ZERO    ]*(F1*lbModel.cx[lbModel.DV_P1_M1_ZERO    ] + F2*lbModel.cy[lbModel.DV_P1_M1_ZERO    ] + F3*lbModel.cz[lbModel.DV_P1_M1_ZERO    ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_ZERO_M1_ZERO  ] = lbModel.wt[lbModel.DV_ZERO_M1_ZERO  ]*(F1*lbModel.cx[lbModel.DV_ZERO_M1_ZERO  ] + F2*lbModel.cy[lbModel.DV_ZERO_M1_ZERO  ] + F3*lbModel.cz[lbModel.DV_ZERO_M1_ZERO  ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M1_ZERO_ZERO  ] = lbModel.wt[lbModel.DV_M1_ZERO_ZERO  ]*(F1*lbModel.cx[lbModel.DV_M1_ZERO_ZERO  ] + F2*lbModel.cy[lbModel.DV_M1_ZERO_ZERO  ] + F3*lbModel.cz[lbModel.DV_M1_ZERO_ZERO  ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P_P_P         ] = lbModel.wt[lbModel.DV_P_P_P         ]*(F1*lbModel.cx[lbModel.DV_P_P_P         ] + F2*lbModel.cy[lbModel.DV_P_P_P         ] + F3*lbModel.cz[lbModel.DV_P_P_P         ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M_P_P         ] = lbModel.wt[lbModel.DV_M_P_P         ]*(F1*lbModel.cx[lbModel.DV_M_P_P         ] + F2*lbModel.cy[lbModel.DV_M_P_P         ] + F3*lbModel.cz[lbModel.DV_M_P_P         ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M_M_P         ] = lbModel.wt[lbModel.DV_M_M_P         ]*(F1*lbModel.cx[lbModel.DV_M_M_P         ] + F2*lbModel.cy[lbModel.DV_M_M_P         ] + F3*lbModel.cz[lbModel.DV_M_M_P         ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P_M_P         ] = lbModel.wt[lbModel.DV_P_M_P         ]*(F1*lbModel.cx[lbModel.DV_P_M_P         ] + F2*lbModel.cy[lbModel.DV_P_M_P         ] + F3*lbModel.cz[lbModel.DV_P_M_P         ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M_M_M         ] = lbModel.wt[lbModel.DV_M_M_M         ]*(F1*lbModel.cx[lbModel.DV_M_M_M         ] + F2*lbModel.cy[lbModel.DV_M_M_M         ] + F3*lbModel.cz[lbModel.DV_M_M_M         ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P_M_M         ] = lbModel.wt[lbModel.DV_P_M_M         ]*(F1*lbModel.cx[lbModel.DV_P_M_M         ] + F2*lbModel.cy[lbModel.DV_P_M_M         ] + F3*lbModel.cz[lbModel.DV_P_M_M         ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P_P_M         ] = lbModel.wt[lbModel.DV_P_P_M         ]*(F1*lbModel.cx[lbModel.DV_P_P_M         ] + F2*lbModel.cy[lbModel.DV_P_P_M         ] + F3*lbModel.cz[lbModel.DV_P_P_M         ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M_P_M         ] = lbModel.wt[lbModel.DV_M_P_M         ]*(F1*lbModel.cx[lbModel.DV_M_P_M         ] + F2*lbModel.cy[lbModel.DV_M_P_M         ] + F3*lbModel.cz[lbModel.DV_M_P_M         ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P1_P1_P1      ] = lbModel.wt[lbModel.DV_P1_P1_P1      ]*(F1*lbModel.cx[lbModel.DV_P1_P1_P1      ] + F2*lbModel.cy[lbModel.DV_P1_P1_P1      ] + F3*lbModel.cz[lbModel.DV_P1_P1_P1      ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M1_P1_P1      ] = lbModel.wt[lbModel.DV_M1_P1_P1      ]*(F1*lbModel.cx[lbModel.DV_M1_P1_P1      ] + F2*lbModel.cy[lbModel.DV_M1_P1_P1      ] + F3*lbModel.cz[lbModel.DV_M1_P1_P1      ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M1_M1_P1      ] = lbModel.wt[lbModel.DV_M1_M1_P1      ]*(F1*lbModel.cx[lbModel.DV_M1_M1_P1      ] + F2*lbModel.cy[lbModel.DV_M1_M1_P1      ] + F3*lbModel.cz[lbModel.DV_M1_M1_P1      ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P1_M1_P1      ] = lbModel.wt[lbModel.DV_P1_M1_P1      ]*(F1*lbModel.cx[lbModel.DV_P1_M1_P1      ] + F2*lbModel.cy[lbModel.DV_P1_M1_P1      ] + F3*lbModel.cz[lbModel.DV_P1_M1_P1      ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M1_M1_M1      ] = lbModel.wt[lbModel.DV_M1_M1_M1      ]*(F1*lbModel.cx[lbModel.DV_M1_M1_M1      ] + F2*lbModel.cy[lbModel.DV_M1_M1_M1      ] + F3*lbModel.cz[lbModel.DV_M1_M1_M1      ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P1_M1_M1      ] = lbModel.wt[lbModel.DV_P1_M1_M1      ]*(F1*lbModel.cx[lbModel.DV_P1_M1_M1      ] + F2*lbModel.cy[lbModel.DV_P1_M1_M1      ] + F3*lbModel.cz[lbModel.DV_P1_M1_M1      ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_P1_P1_M1      ] = lbModel.wt[lbModel.DV_P1_P1_M1      ]*(F1*lbModel.cx[lbModel.DV_P1_P1_M1      ] + F2*lbModel.cy[lbModel.DV_P1_P1_M1      ] + F3*lbModel.cz[lbModel.DV_P1_P1_M1      ])*oneByTheta0*rho[index];
                    fTempF[index][lbModel.DV_M1_P1_M1      ] = lbModel.wt[lbModel.DV_M1_P1_M1      ]*(F1*lbModel.cx[lbModel.DV_M1_P1_M1      ] + F2*lbModel.cy[lbModel.DV_M1_P1_M1      ] + F3*lbModel.cz[lbModel.DV_M1_P1_M1      ])*oneByTheta0*rho[index];
                }


                for(int index=0;index<VECT_LENGTH;index++)
                {
                    myGrid(i1+index,i2,i3,lbModel.G0 ,myGrid.node[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO) = myGrid(i1+index,i2,i3,lbModel.G0 ,myGrid.node[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO) + twoBeta*(lbModel.fTemp0 [index][lbModel.CENTER_DV_ZERO_ZERO_ZERO] - myGrid(i1+index,i2,i3,lbModel.G0 ,myGrid.node[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO))  +   fTempF[index][lbModel.DV_ZERO_ZERO_ZERO]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P1      ) = myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P1      ) + twoBeta*(lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_ZERO_P1      ] - myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P1      ))  +   fTempF[index][lbModel.DV_ZERO_ZERO_P1  ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P2      ) = myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P2      ) + twoBeta*(lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_ZERO_P2      ] - myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P2      ))  +   fTempF[index][lbModel.DV_ZERO_ZERO_P2  ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_ZERO_P2_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_ZERO_P2_ZERO      ) + twoBeta*(lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_P2_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_ZERO_P2_ZERO      ))  +   fTempF[index][lbModel.DV_ZERO_P2_ZERO  ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_P2_ZERO_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_P2_ZERO_ZERO      ) + twoBeta*(lbModel.fTemp1 [index][lbModel.G1_DV_P2_ZERO_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.node[lbModel.G1 ],lbModel.G1_DV_P2_ZERO_ZERO      ))  +   fTempF[index][lbModel.DV_P2_ZERO_ZERO  ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M1      ) = myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M1      ) + twoBeta*(lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_ZERO_M1      ] - myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M1      ))  +   fTempF[index][lbModel.DV_ZERO_ZERO_M1  ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M2      ) = myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M2      ) + twoBeta*(lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_ZERO_M2      ] - myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M2      ))  +   fTempF[index][lbModel.DV_ZERO_ZERO_M2  ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_ZERO_M2_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_ZERO_M2_ZERO      ) + twoBeta*(lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_M2_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_ZERO_M2_ZERO      ))  +   fTempF[index][lbModel.DV_ZERO_M2_ZERO  ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_M2_ZERO_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_M2_ZERO_ZERO      ) + twoBeta*(lbModel.fTemp2 [index][lbModel.G2_DV_M2_ZERO_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.node[lbModel.G2 ],lbModel.G2_DV_M2_ZERO_ZERO      ))  +   fTempF[index][lbModel.DV_M2_ZERO_ZERO  ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_ZERO_P1_P1        ) = myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_ZERO_P1_P1        ) + twoBeta*(lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_P1_P1        ] - myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_ZERO_P1_P1        ))  +   fTempF[index][lbModel.DV_ZERO_P1_P1    ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_ZERO_M1_P1        ) = myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_ZERO_M1_P1        ) + twoBeta*(lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_M1_P1        ] - myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_ZERO_M1_P1        ))  +   fTempF[index][lbModel.DV_ZERO_M1_P1    ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_P1_ZERO_P1        ) = myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_P1_ZERO_P1        ) + twoBeta*(lbModel.fTemp3 [index][lbModel.G3_DV_P1_ZERO_P1        ] - myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_P1_ZERO_P1        ))  +   fTempF[index][lbModel.DV_P1_ZERO_P1    ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_M1_ZERO_P1        ) = myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_M1_ZERO_P1        ) + twoBeta*(lbModel.fTemp3 [index][lbModel.G3_DV_M1_ZERO_P1        ] - myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.node[lbModel.G3 ],lbModel.G3_DV_M1_ZERO_P1        ))  +   fTempF[index][lbModel.DV_M1_ZERO_P1    ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_ZERO_M1_M1        ) = myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_ZERO_M1_M1        ) + twoBeta*(lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_M1_M1        ] - myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_ZERO_M1_M1        ))  +   fTempF[index][lbModel.DV_ZERO_M1_M1    ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_ZERO_P1_M1        ) = myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_ZERO_P1_M1        ) + twoBeta*(lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_P1_M1        ] - myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_ZERO_P1_M1        ))  +   fTempF[index][lbModel.DV_ZERO_P1_M1    ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_M1_ZERO_M1        ) = myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_M1_ZERO_M1        ) + twoBeta*(lbModel.fTemp4 [index][lbModel.G4_DV_M1_ZERO_M1        ] - myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_M1_ZERO_M1        ))  +   fTempF[index][lbModel.DV_M1_ZERO_M1    ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_P1_ZERO_M1        ) = myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_P1_ZERO_M1        ) + twoBeta*(lbModel.fTemp4 [index][lbModel.G4_DV_P1_ZERO_M1        ] - myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.node[lbModel.G4 ],lbModel.G4_DV_P1_ZERO_M1        ))  +   fTempF[index][lbModel.DV_P1_ZERO_M1    ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_P1_P1_ZERO        ) = myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_P1_P1_ZERO        ) + twoBeta*(lbModel.fTemp5 [index][lbModel.G5_DV_P1_P1_ZERO        ] - myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_P1_P1_ZERO        ))  +   fTempF[index][lbModel.DV_P1_P1_ZERO    ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_M1_P1_ZERO        ) = myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_M1_P1_ZERO        ) + twoBeta*(lbModel.fTemp5 [index][lbModel.G5_DV_M1_P1_ZERO        ] - myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_M1_P1_ZERO        ))  +   fTempF[index][lbModel.DV_M1_P1_ZERO    ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_ZERO_P1_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_ZERO_P1_ZERO      ) + twoBeta*(lbModel.fTemp5 [index][lbModel.G5_DV_ZERO_P1_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_ZERO_P1_ZERO      ))  +   fTempF[index][lbModel.DV_ZERO_P1_ZERO  ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_P1_ZERO_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_P1_ZERO_ZERO      ) + twoBeta*(lbModel.fTemp5 [index][lbModel.G5_DV_P1_ZERO_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.node[lbModel.G5 ],lbModel.G5_DV_P1_ZERO_ZERO      ))  +   fTempF[index][lbModel.DV_P1_ZERO_ZERO  ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_M1_M1_ZERO        ) = myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_M1_M1_ZERO        ) + twoBeta*(lbModel.fTemp6 [index][lbModel.G6_DV_M1_M1_ZERO        ] - myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_M1_M1_ZERO        ))  +   fTempF[index][lbModel.DV_M1_M1_ZERO    ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_P1_M1_ZERO        ) = myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_P1_M1_ZERO        ) + twoBeta*(lbModel.fTemp6 [index][lbModel.G6_DV_P1_M1_ZERO        ] - myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_P1_M1_ZERO        ))  +   fTempF[index][lbModel.DV_P1_M1_ZERO    ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_ZERO_M1_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_ZERO_M1_ZERO      ) + twoBeta*(lbModel.fTemp6 [index][lbModel.G6_DV_ZERO_M1_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_ZERO_M1_ZERO      ))  +   fTempF[index][lbModel.DV_ZERO_M1_ZERO  ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_M1_ZERO_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_M1_ZERO_ZERO      ) + twoBeta*(lbModel.fTemp6 [index][lbModel.G6_DV_M1_ZERO_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.node[lbModel.G6 ],lbModel.G6_DV_M1_ZERO_ZERO      ))  +   fTempF[index][lbModel.DV_M1_ZERO_ZERO  ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_P_P_P             ) = myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_P_P_P             ) + twoBeta*(lbModel.fTemp7 [index][lbModel.G7_DV_P_P_P             ] - myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_P_P_P             ))  +   fTempF[index][lbModel.DV_P_P_P         ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_M_P_P             ) = myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_M_P_P             ) + twoBeta*(lbModel.fTemp7 [index][lbModel.G7_DV_M_P_P             ] - myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_M_P_P             ))  +   fTempF[index][lbModel.DV_M_P_P         ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_M_M_P             ) = myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_M_M_P             ) + twoBeta*(lbModel.fTemp7 [index][lbModel.G7_DV_M_M_P             ] - myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_M_M_P             ))  +   fTempF[index][lbModel.DV_M_M_P         ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_P_M_P             ) = myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_P_M_P             ) + twoBeta*(lbModel.fTemp7 [index][lbModel.G7_DV_P_M_P             ] - myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.node[lbModel.G7 ],lbModel.G7_DV_P_M_P             ))  +   fTempF[index][lbModel.DV_P_M_P         ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_M_M_M             ) = myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_M_M_M             ) + twoBeta*(lbModel.fTemp8 [index][lbModel.G8_DV_M_M_M             ] - myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_M_M_M             ))  +   fTempF[index][lbModel.DV_M_M_M         ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_P_M_M             ) = myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_P_M_M             ) + twoBeta*(lbModel.fTemp8 [index][lbModel.G8_DV_P_M_M             ] - myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_P_M_M             ))  +   fTempF[index][lbModel.DV_P_M_M         ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_P_P_M             ) = myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_P_P_M             ) + twoBeta*(lbModel.fTemp8 [index][lbModel.G8_DV_P_P_M             ] - myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_P_P_M             ))  +   fTempF[index][lbModel.DV_P_P_M         ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_M_P_M             ) = myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_M_P_M             ) + twoBeta*(lbModel.fTemp8 [index][lbModel.G8_DV_M_P_M             ] - myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.node[lbModel.G8 ],lbModel.G8_DV_M_P_M             ))  +   fTempF[index][lbModel.DV_M_P_M         ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_P1_P1_P1          ) = myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_P1_P1_P1          ) + twoBeta*(lbModel.fTemp9 [index][lbModel.G9_DV_P1_P1_P1          ] - myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_P1_P1_P1          ))  +   fTempF[index][lbModel.DV_P1_P1_P1      ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_M1_P1_P1          ) = myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_M1_P1_P1          ) + twoBeta*(lbModel.fTemp9 [index][lbModel.G9_DV_M1_P1_P1          ] - myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_M1_P1_P1          ))  +   fTempF[index][lbModel.DV_M1_P1_P1      ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_M1_M1_P1          ) = myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_M1_M1_P1          ) + twoBeta*(lbModel.fTemp9 [index][lbModel.G9_DV_M1_M1_P1          ] - myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_M1_M1_P1          ))  +   fTempF[index][lbModel.DV_M1_M1_P1      ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_P1_M1_P1          ) = myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_P1_M1_P1          ) + twoBeta*(lbModel.fTemp9 [index][lbModel.G9_DV_P1_M1_P1          ] - myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.node[lbModel.G9 ],lbModel.G9_DV_P1_M1_P1          ))  +   fTempF[index][lbModel.DV_P1_M1_P1      ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1         ) = myGrid(i1+index,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1         ) + twoBeta*(lbModel.fTemp10[index][lbModel.G10_DV_M1_M1_M1         ] - myGrid(i1+index,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1         ))  +   fTempF[index][lbModel.DV_M1_M1_M1      ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1         ) = myGrid(i1+index,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1         ) + twoBeta*(lbModel.fTemp10[index][lbModel.G10_DV_P1_M1_M1         ] - myGrid(i1+index,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1         ))  +   fTempF[index][lbModel.DV_P1_M1_M1      ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1         ) = myGrid(i1+index,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1         ) + twoBeta*(lbModel.fTemp10[index][lbModel.G10_DV_P1_P1_M1         ] - myGrid(i1+index,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1         ))  +   fTempF[index][lbModel.DV_P1_P1_M1      ]*twoBetaTau  ;
                    myGrid(i1+index,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1         ) = myGrid(i1+index,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1         ) + twoBeta*(lbModel.fTemp10[index][lbModel.G10_DV_M1_P1_M1         ] - myGrid(i1+index,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1         ))  +   fTempF[index][lbModel.DV_M1_P1_M1      ]*twoBetaTau  ;
                }
            }


            for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
                for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
                    for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1=i1+4)
                    {
                        getHydroMomentsFromCellWithForce(lbModel,myGrid,VECT_LENGTH,i1,i2,i3,rho,uX,uY,uZ,theta,F1,F2,F3,dt);

                        for(int index=0;index<VECT_LENGTH;index++)
                            theta[index] = lbModel.theta0;

                        getFEq(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);


                        for(int index=0;index<VECT_LENGTH;index++)
                        {
                            fTempF[index][lbModel.DV_ZERO_ZERO_ZERO] = lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO]*(F1*lbModel.cx[lbModel.DV_ZERO_ZERO_ZERO] + F2*lbModel.cy[lbModel.DV_ZERO_ZERO_ZERO] + F3*lbModel.cz[lbModel.DV_ZERO_ZERO_ZERO])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_ZERO_ZERO_P1  ] = lbModel.wt[lbModel.DV_ZERO_ZERO_P1  ]*(F1*lbModel.cx[lbModel.DV_ZERO_ZERO_P1  ] + F2*lbModel.cy[lbModel.DV_ZERO_ZERO_P1  ] + F3*lbModel.cz[lbModel.DV_ZERO_ZERO_P1  ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_ZERO_ZERO_P2  ] = lbModel.wt[lbModel.DV_ZERO_ZERO_P2  ]*(F1*lbModel.cx[lbModel.DV_ZERO_ZERO_P2  ] + F2*lbModel.cy[lbModel.DV_ZERO_ZERO_P2  ] + F3*lbModel.cz[lbModel.DV_ZERO_ZERO_P2  ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_ZERO_P2_ZERO  ] = lbModel.wt[lbModel.DV_ZERO_P2_ZERO  ]*(F1*lbModel.cx[lbModel.DV_ZERO_P2_ZERO  ] + F2*lbModel.cy[lbModel.DV_ZERO_P2_ZERO  ] + F3*lbModel.cz[lbModel.DV_ZERO_P2_ZERO  ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P2_ZERO_ZERO  ] = lbModel.wt[lbModel.DV_P2_ZERO_ZERO  ]*(F1*lbModel.cx[lbModel.DV_P2_ZERO_ZERO  ] + F2*lbModel.cy[lbModel.DV_P2_ZERO_ZERO  ] + F3*lbModel.cz[lbModel.DV_P2_ZERO_ZERO  ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_ZERO_ZERO_M1  ] = lbModel.wt[lbModel.DV_ZERO_ZERO_M1  ]*(F1*lbModel.cx[lbModel.DV_ZERO_ZERO_M1  ] + F2*lbModel.cy[lbModel.DV_ZERO_ZERO_M1  ] + F3*lbModel.cz[lbModel.DV_ZERO_ZERO_M1  ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_ZERO_ZERO_M2  ] = lbModel.wt[lbModel.DV_ZERO_ZERO_M2  ]*(F1*lbModel.cx[lbModel.DV_ZERO_ZERO_M2  ] + F2*lbModel.cy[lbModel.DV_ZERO_ZERO_M2  ] + F3*lbModel.cz[lbModel.DV_ZERO_ZERO_M2  ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_ZERO_M2_ZERO  ] = lbModel.wt[lbModel.DV_ZERO_M2_ZERO  ]*(F1*lbModel.cx[lbModel.DV_ZERO_M2_ZERO  ] + F2*lbModel.cy[lbModel.DV_ZERO_M2_ZERO  ] + F3*lbModel.cz[lbModel.DV_ZERO_M2_ZERO  ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M2_ZERO_ZERO  ] = lbModel.wt[lbModel.DV_M2_ZERO_ZERO  ]*(F1*lbModel.cx[lbModel.DV_M2_ZERO_ZERO  ] + F2*lbModel.cy[lbModel.DV_M2_ZERO_ZERO  ] + F3*lbModel.cz[lbModel.DV_M2_ZERO_ZERO  ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_ZERO_P1_P1    ] = lbModel.wt[lbModel.DV_ZERO_P1_P1    ]*(F1*lbModel.cx[lbModel.DV_ZERO_P1_P1    ] + F2*lbModel.cy[lbModel.DV_ZERO_P1_P1    ] + F3*lbModel.cz[lbModel.DV_ZERO_P1_P1    ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_ZERO_M1_P1    ] = lbModel.wt[lbModel.DV_ZERO_M1_P1    ]*(F1*lbModel.cx[lbModel.DV_ZERO_M1_P1    ] + F2*lbModel.cy[lbModel.DV_ZERO_M1_P1    ] + F3*lbModel.cz[lbModel.DV_ZERO_M1_P1    ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P1_ZERO_P1    ] = lbModel.wt[lbModel.DV_P1_ZERO_P1    ]*(F1*lbModel.cx[lbModel.DV_P1_ZERO_P1    ] + F2*lbModel.cy[lbModel.DV_P1_ZERO_P1    ] + F3*lbModel.cz[lbModel.DV_P1_ZERO_P1    ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M1_ZERO_P1    ] = lbModel.wt[lbModel.DV_M1_ZERO_P1    ]*(F1*lbModel.cx[lbModel.DV_M1_ZERO_P1    ] + F2*lbModel.cy[lbModel.DV_M1_ZERO_P1    ] + F3*lbModel.cz[lbModel.DV_M1_ZERO_P1    ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_ZERO_M1_M1    ] = lbModel.wt[lbModel.DV_ZERO_M1_M1    ]*(F1*lbModel.cx[lbModel.DV_ZERO_M1_M1    ] + F2*lbModel.cy[lbModel.DV_ZERO_M1_M1    ] + F3*lbModel.cz[lbModel.DV_ZERO_M1_M1    ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_ZERO_P1_M1    ] = lbModel.wt[lbModel.DV_ZERO_P1_M1    ]*(F1*lbModel.cx[lbModel.DV_ZERO_P1_M1    ] + F2*lbModel.cy[lbModel.DV_ZERO_P1_M1    ] + F3*lbModel.cz[lbModel.DV_ZERO_P1_M1    ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M1_ZERO_M1    ] = lbModel.wt[lbModel.DV_M1_ZERO_M1    ]*(F1*lbModel.cx[lbModel.DV_M1_ZERO_M1    ] + F2*lbModel.cy[lbModel.DV_M1_ZERO_M1    ] + F3*lbModel.cz[lbModel.DV_M1_ZERO_M1    ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P1_ZERO_M1    ] = lbModel.wt[lbModel.DV_P1_ZERO_M1    ]*(F1*lbModel.cx[lbModel.DV_P1_ZERO_M1    ] + F2*lbModel.cy[lbModel.DV_P1_ZERO_M1    ] + F3*lbModel.cz[lbModel.DV_P1_ZERO_M1    ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P1_P1_ZERO    ] = lbModel.wt[lbModel.DV_P1_P1_ZERO    ]*(F1*lbModel.cx[lbModel.DV_P1_P1_ZERO    ] + F2*lbModel.cy[lbModel.DV_P1_P1_ZERO    ] + F3*lbModel.cz[lbModel.DV_P1_P1_ZERO    ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M1_P1_ZERO    ] = lbModel.wt[lbModel.DV_M1_P1_ZERO    ]*(F1*lbModel.cx[lbModel.DV_M1_P1_ZERO    ] + F2*lbModel.cy[lbModel.DV_M1_P1_ZERO    ] + F3*lbModel.cz[lbModel.DV_M1_P1_ZERO    ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_ZERO_P1_ZERO  ] = lbModel.wt[lbModel.DV_ZERO_P1_ZERO  ]*(F1*lbModel.cx[lbModel.DV_ZERO_P1_ZERO  ] + F2*lbModel.cy[lbModel.DV_ZERO_P1_ZERO  ] + F3*lbModel.cz[lbModel.DV_ZERO_P1_ZERO  ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P1_ZERO_ZERO  ] = lbModel.wt[lbModel.DV_P1_ZERO_ZERO  ]*(F1*lbModel.cx[lbModel.DV_P1_ZERO_ZERO  ] + F2*lbModel.cy[lbModel.DV_P1_ZERO_ZERO  ] + F3*lbModel.cz[lbModel.DV_P1_ZERO_ZERO  ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M1_M1_ZERO    ] = lbModel.wt[lbModel.DV_M1_M1_ZERO    ]*(F1*lbModel.cx[lbModel.DV_M1_M1_ZERO    ] + F2*lbModel.cy[lbModel.DV_M1_M1_ZERO    ] + F3*lbModel.cz[lbModel.DV_M1_M1_ZERO    ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P1_M1_ZERO    ] = lbModel.wt[lbModel.DV_P1_M1_ZERO    ]*(F1*lbModel.cx[lbModel.DV_P1_M1_ZERO    ] + F2*lbModel.cy[lbModel.DV_P1_M1_ZERO    ] + F3*lbModel.cz[lbModel.DV_P1_M1_ZERO    ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_ZERO_M1_ZERO  ] = lbModel.wt[lbModel.DV_ZERO_M1_ZERO  ]*(F1*lbModel.cx[lbModel.DV_ZERO_M1_ZERO  ] + F2*lbModel.cy[lbModel.DV_ZERO_M1_ZERO  ] + F3*lbModel.cz[lbModel.DV_ZERO_M1_ZERO  ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M1_ZERO_ZERO  ] = lbModel.wt[lbModel.DV_M1_ZERO_ZERO  ]*(F1*lbModel.cx[lbModel.DV_M1_ZERO_ZERO  ] + F2*lbModel.cy[lbModel.DV_M1_ZERO_ZERO  ] + F3*lbModel.cz[lbModel.DV_M1_ZERO_ZERO  ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P_P_P         ] = lbModel.wt[lbModel.DV_P_P_P         ]*(F1*lbModel.cx[lbModel.DV_P_P_P         ] + F2*lbModel.cy[lbModel.DV_P_P_P         ] + F3*lbModel.cz[lbModel.DV_P_P_P         ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M_P_P         ] = lbModel.wt[lbModel.DV_M_P_P         ]*(F1*lbModel.cx[lbModel.DV_M_P_P         ] + F2*lbModel.cy[lbModel.DV_M_P_P         ] + F3*lbModel.cz[lbModel.DV_M_P_P         ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M_M_P         ] = lbModel.wt[lbModel.DV_M_M_P         ]*(F1*lbModel.cx[lbModel.DV_M_M_P         ] + F2*lbModel.cy[lbModel.DV_M_M_P         ] + F3*lbModel.cz[lbModel.DV_M_M_P         ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P_M_P         ] = lbModel.wt[lbModel.DV_P_M_P         ]*(F1*lbModel.cx[lbModel.DV_P_M_P         ] + F2*lbModel.cy[lbModel.DV_P_M_P         ] + F3*lbModel.cz[lbModel.DV_P_M_P         ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M_M_M         ] = lbModel.wt[lbModel.DV_M_M_M         ]*(F1*lbModel.cx[lbModel.DV_M_M_M         ] + F2*lbModel.cy[lbModel.DV_M_M_M         ] + F3*lbModel.cz[lbModel.DV_M_M_M         ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P_M_M         ] = lbModel.wt[lbModel.DV_P_M_M         ]*(F1*lbModel.cx[lbModel.DV_P_M_M         ] + F2*lbModel.cy[lbModel.DV_P_M_M         ] + F3*lbModel.cz[lbModel.DV_P_M_M         ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P_P_M         ] = lbModel.wt[lbModel.DV_P_P_M         ]*(F1*lbModel.cx[lbModel.DV_P_P_M         ] + F2*lbModel.cy[lbModel.DV_P_P_M         ] + F3*lbModel.cz[lbModel.DV_P_P_M         ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M_P_M         ] = lbModel.wt[lbModel.DV_M_P_M         ]*(F1*lbModel.cx[lbModel.DV_M_P_M         ] + F2*lbModel.cy[lbModel.DV_M_P_M         ] + F3*lbModel.cz[lbModel.DV_M_P_M         ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P1_P1_P1      ] = lbModel.wt[lbModel.DV_P1_P1_P1      ]*(F1*lbModel.cx[lbModel.DV_P1_P1_P1      ] + F2*lbModel.cy[lbModel.DV_P1_P1_P1      ] + F3*lbModel.cz[lbModel.DV_P1_P1_P1      ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M1_P1_P1      ] = lbModel.wt[lbModel.DV_M1_P1_P1      ]*(F1*lbModel.cx[lbModel.DV_M1_P1_P1      ] + F2*lbModel.cy[lbModel.DV_M1_P1_P1      ] + F3*lbModel.cz[lbModel.DV_M1_P1_P1      ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M1_M1_P1      ] = lbModel.wt[lbModel.DV_M1_M1_P1      ]*(F1*lbModel.cx[lbModel.DV_M1_M1_P1      ] + F2*lbModel.cy[lbModel.DV_M1_M1_P1      ] + F3*lbModel.cz[lbModel.DV_M1_M1_P1      ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P1_M1_P1      ] = lbModel.wt[lbModel.DV_P1_M1_P1      ]*(F1*lbModel.cx[lbModel.DV_P1_M1_P1      ] + F2*lbModel.cy[lbModel.DV_P1_M1_P1      ] + F3*lbModel.cz[lbModel.DV_P1_M1_P1      ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M1_M1_M1      ] = lbModel.wt[lbModel.DV_M1_M1_M1      ]*(F1*lbModel.cx[lbModel.DV_M1_M1_M1      ] + F2*lbModel.cy[lbModel.DV_M1_M1_M1      ] + F3*lbModel.cz[lbModel.DV_M1_M1_M1      ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P1_M1_M1      ] = lbModel.wt[lbModel.DV_P1_M1_M1      ]*(F1*lbModel.cx[lbModel.DV_P1_M1_M1      ] + F2*lbModel.cy[lbModel.DV_P1_M1_M1      ] + F3*lbModel.cz[lbModel.DV_P1_M1_M1      ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_P1_P1_M1      ] = lbModel.wt[lbModel.DV_P1_P1_M1      ]*(F1*lbModel.cx[lbModel.DV_P1_P1_M1      ] + F2*lbModel.cy[lbModel.DV_P1_P1_M1      ] + F3*lbModel.cz[lbModel.DV_P1_P1_M1      ])*oneByTheta0*rho[index];
                            fTempF[index][lbModel.DV_M1_P1_M1      ] = lbModel.wt[lbModel.DV_M1_P1_M1      ]*(F1*lbModel.cx[lbModel.DV_M1_P1_M1      ] + F2*lbModel.cy[lbModel.DV_M1_P1_M1      ] + F3*lbModel.cz[lbModel.DV_M1_P1_M1      ])*oneByTheta0*rho[index];
                        }


                        for(int index=0;index<VECT_LENGTH;index++)
                        {
                            myGrid(i1+index,i2,i3,lbModel.G0 ,myGrid.cell[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO) = myGrid(i1+index,i2,i3,lbModel.G0 ,myGrid.cell[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO) + twoBeta*(lbModel.fTemp0 [index][lbModel.CENTER_DV_ZERO_ZERO_ZERO] - myGrid(i1+index,i2,i3,lbModel.G0 ,myGrid.cell[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO))  +   fTempF[index][lbModel.DV_ZERO_ZERO_ZERO]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.cell[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P1      ) = myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.cell[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P1      ) + twoBeta*(lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_ZERO_P1      ] - myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.cell[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P1      ))  +   fTempF[index][lbModel.DV_ZERO_ZERO_P1  ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.cell[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P2      ) = myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.cell[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P2      ) + twoBeta*(lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_ZERO_P2      ] - myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.cell[lbModel.G1 ],lbModel.G1_DV_ZERO_ZERO_P2      ))  +   fTempF[index][lbModel.DV_ZERO_ZERO_P2  ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.cell[lbModel.G1 ],lbModel.G1_DV_ZERO_P2_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.cell[lbModel.G1 ],lbModel.G1_DV_ZERO_P2_ZERO      ) + twoBeta*(lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_P2_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.cell[lbModel.G1 ],lbModel.G1_DV_ZERO_P2_ZERO      ))  +   fTempF[index][lbModel.DV_ZERO_P2_ZERO  ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.cell[lbModel.G1 ],lbModel.G1_DV_P2_ZERO_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.cell[lbModel.G1 ],lbModel.G1_DV_P2_ZERO_ZERO      ) + twoBeta*(lbModel.fTemp1 [index][lbModel.G1_DV_P2_ZERO_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G1 ,myGrid.cell[lbModel.G1 ],lbModel.G1_DV_P2_ZERO_ZERO      ))  +   fTempF[index][lbModel.DV_P2_ZERO_ZERO  ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.cell[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M1      ) = myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.cell[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M1      ) + twoBeta*(lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_ZERO_M1      ] - myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.cell[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M1      ))  +   fTempF[index][lbModel.DV_ZERO_ZERO_M1  ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.cell[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M2      ) = myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.cell[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M2      ) + twoBeta*(lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_ZERO_M2      ] - myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.cell[lbModel.G2 ],lbModel.G2_DV_ZERO_ZERO_M2      ))  +   fTempF[index][lbModel.DV_ZERO_ZERO_M2  ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.cell[lbModel.G2 ],lbModel.G2_DV_ZERO_M2_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.cell[lbModel.G2 ],lbModel.G2_DV_ZERO_M2_ZERO      ) + twoBeta*(lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_M2_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.cell[lbModel.G2 ],lbModel.G2_DV_ZERO_M2_ZERO      ))  +   fTempF[index][lbModel.DV_ZERO_M2_ZERO  ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.cell[lbModel.G2 ],lbModel.G2_DV_M2_ZERO_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.cell[lbModel.G2 ],lbModel.G2_DV_M2_ZERO_ZERO      ) + twoBeta*(lbModel.fTemp2 [index][lbModel.G2_DV_M2_ZERO_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G2 ,myGrid.cell[lbModel.G2 ],lbModel.G2_DV_M2_ZERO_ZERO      ))  +   fTempF[index][lbModel.DV_M2_ZERO_ZERO  ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_ZERO_P1_P1        ) = myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_ZERO_P1_P1        ) + twoBeta*(lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_P1_P1        ] - myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_ZERO_P1_P1        ))  +   fTempF[index][lbModel.DV_ZERO_P1_P1    ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_ZERO_M1_P1        ) = myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_ZERO_M1_P1        ) + twoBeta*(lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_M1_P1        ] - myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_ZERO_M1_P1        ))  +   fTempF[index][lbModel.DV_ZERO_M1_P1    ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_P1_ZERO_P1        ) = myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_P1_ZERO_P1        ) + twoBeta*(lbModel.fTemp3 [index][lbModel.G3_DV_P1_ZERO_P1        ] - myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_P1_ZERO_P1        ))  +   fTempF[index][lbModel.DV_P1_ZERO_P1    ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_M1_ZERO_P1        ) = myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_M1_ZERO_P1        ) + twoBeta*(lbModel.fTemp3 [index][lbModel.G3_DV_M1_ZERO_P1        ] - myGrid(i1+index,i2,i3,lbModel.G3 ,myGrid.cell[lbModel.G3 ],lbModel.G3_DV_M1_ZERO_P1        ))  +   fTempF[index][lbModel.DV_M1_ZERO_P1    ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_ZERO_M1_M1        ) = myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_ZERO_M1_M1        ) + twoBeta*(lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_M1_M1        ] - myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_ZERO_M1_M1        ))  +   fTempF[index][lbModel.DV_ZERO_M1_M1    ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_ZERO_P1_M1        ) = myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_ZERO_P1_M1        ) + twoBeta*(lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_P1_M1        ] - myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_ZERO_P1_M1        ))  +   fTempF[index][lbModel.DV_ZERO_P1_M1    ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_M1_ZERO_M1        ) = myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_M1_ZERO_M1        ) + twoBeta*(lbModel.fTemp4 [index][lbModel.G4_DV_M1_ZERO_M1        ] - myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_M1_ZERO_M1        ))  +   fTempF[index][lbModel.DV_M1_ZERO_M1    ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_P1_ZERO_M1        ) = myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_P1_ZERO_M1        ) + twoBeta*(lbModel.fTemp4 [index][lbModel.G4_DV_P1_ZERO_M1        ] - myGrid(i1+index,i2,i3,lbModel.G4 ,myGrid.cell[lbModel.G4 ],lbModel.G4_DV_P1_ZERO_M1        ))  +   fTempF[index][lbModel.DV_P1_ZERO_M1    ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_P1_P1_ZERO        ) = myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_P1_P1_ZERO        ) + twoBeta*(lbModel.fTemp5 [index][lbModel.G5_DV_P1_P1_ZERO        ] - myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_P1_P1_ZERO        ))  +   fTempF[index][lbModel.DV_P1_P1_ZERO    ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_M1_P1_ZERO        ) = myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_M1_P1_ZERO        ) + twoBeta*(lbModel.fTemp5 [index][lbModel.G5_DV_M1_P1_ZERO        ] - myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_M1_P1_ZERO        ))  +   fTempF[index][lbModel.DV_M1_P1_ZERO    ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_ZERO_P1_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_ZERO_P1_ZERO      ) + twoBeta*(lbModel.fTemp5 [index][lbModel.G5_DV_ZERO_P1_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_ZERO_P1_ZERO      ))  +   fTempF[index][lbModel.DV_ZERO_P1_ZERO  ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_P1_ZERO_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_P1_ZERO_ZERO      ) + twoBeta*(lbModel.fTemp5 [index][lbModel.G5_DV_P1_ZERO_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G5 ,myGrid.cell[lbModel.G5 ],lbModel.G5_DV_P1_ZERO_ZERO      ))  +   fTempF[index][lbModel.DV_P1_ZERO_ZERO  ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_M1_M1_ZERO        ) = myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_M1_M1_ZERO        ) + twoBeta*(lbModel.fTemp6 [index][lbModel.G6_DV_M1_M1_ZERO        ] - myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_M1_M1_ZERO        ))  +   fTempF[index][lbModel.DV_M1_M1_ZERO    ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_P1_M1_ZERO        ) = myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_P1_M1_ZERO        ) + twoBeta*(lbModel.fTemp6 [index][lbModel.G6_DV_P1_M1_ZERO        ] - myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_P1_M1_ZERO        ))  +   fTempF[index][lbModel.DV_P1_M1_ZERO    ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_ZERO_M1_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_ZERO_M1_ZERO      ) + twoBeta*(lbModel.fTemp6 [index][lbModel.G6_DV_ZERO_M1_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_ZERO_M1_ZERO      ))  +   fTempF[index][lbModel.DV_ZERO_M1_ZERO  ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_M1_ZERO_ZERO      ) = myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_M1_ZERO_ZERO      ) + twoBeta*(lbModel.fTemp6 [index][lbModel.G6_DV_M1_ZERO_ZERO      ] - myGrid(i1+index,i2,i3,lbModel.G6 ,myGrid.cell[lbModel.G6 ],lbModel.G6_DV_M1_ZERO_ZERO      ))  +   fTempF[index][lbModel.DV_M1_ZERO_ZERO  ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_P_P_P             ) = myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_P_P_P             ) + twoBeta*(lbModel.fTemp7 [index][lbModel.G7_DV_P_P_P             ] - myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_P_P_P             ))  +   fTempF[index][lbModel.DV_P_P_P         ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_M_P_P             ) = myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_M_P_P             ) + twoBeta*(lbModel.fTemp7 [index][lbModel.G7_DV_M_P_P             ] - myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_M_P_P             ))  +   fTempF[index][lbModel.DV_M_P_P         ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_M_M_P             ) = myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_M_M_P             ) + twoBeta*(lbModel.fTemp7 [index][lbModel.G7_DV_M_M_P             ] - myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_M_M_P             ))  +   fTempF[index][lbModel.DV_M_M_P         ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_P_M_P             ) = myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_P_M_P             ) + twoBeta*(lbModel.fTemp7 [index][lbModel.G7_DV_P_M_P             ] - myGrid(i1+index,i2,i3,lbModel.G7 ,myGrid.cell[lbModel.G7 ],lbModel.G7_DV_P_M_P             ))  +   fTempF[index][lbModel.DV_P_M_P         ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_M_M_M             ) = myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_M_M_M             ) + twoBeta*(lbModel.fTemp8 [index][lbModel.G8_DV_M_M_M             ] - myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_M_M_M             ))  +   fTempF[index][lbModel.DV_M_M_M         ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_P_M_M             ) = myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_P_M_M             ) + twoBeta*(lbModel.fTemp8 [index][lbModel.G8_DV_P_M_M             ] - myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_P_M_M             ))  +   fTempF[index][lbModel.DV_P_M_M         ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_P_P_M             ) = myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_P_P_M             ) + twoBeta*(lbModel.fTemp8 [index][lbModel.G8_DV_P_P_M             ] - myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_P_P_M             ))  +   fTempF[index][lbModel.DV_P_P_M         ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_M_P_M             ) = myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_M_P_M             ) + twoBeta*(lbModel.fTemp8 [index][lbModel.G8_DV_M_P_M             ] - myGrid(i1+index,i2,i3,lbModel.G8 ,myGrid.cell[lbModel.G8 ],lbModel.G8_DV_M_P_M             ))  +   fTempF[index][lbModel.DV_M_P_M         ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_P1_P1_P1          ) = myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_P1_P1_P1          ) + twoBeta*(lbModel.fTemp9 [index][lbModel.G9_DV_P1_P1_P1          ] - myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_P1_P1_P1          ))  +   fTempF[index][lbModel.DV_P1_P1_P1      ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_M1_P1_P1          ) = myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_M1_P1_P1          ) + twoBeta*(lbModel.fTemp9 [index][lbModel.G9_DV_M1_P1_P1          ] - myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_M1_P1_P1          ))  +   fTempF[index][lbModel.DV_M1_P1_P1      ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_M1_M1_P1          ) = myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_M1_M1_P1          ) + twoBeta*(lbModel.fTemp9 [index][lbModel.G9_DV_M1_M1_P1          ] - myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_M1_M1_P1          ))  +   fTempF[index][lbModel.DV_M1_M1_P1      ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_P1_M1_P1          ) = myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_P1_M1_P1          ) + twoBeta*(lbModel.fTemp9 [index][lbModel.G9_DV_P1_M1_P1          ] - myGrid(i1+index,i2,i3,lbModel.G9 ,myGrid.cell[lbModel.G9 ],lbModel.G9_DV_P1_M1_P1          ))  +   fTempF[index][lbModel.DV_P1_M1_P1      ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1         ) = myGrid(i1+index,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1         ) + twoBeta*(lbModel.fTemp10[index][lbModel.G10_DV_M1_M1_M1         ] - myGrid(i1+index,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1         ))  +   fTempF[index][lbModel.DV_M1_M1_M1      ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1         ) = myGrid(i1+index,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1         ) + twoBeta*(lbModel.fTemp10[index][lbModel.G10_DV_P1_M1_M1         ] - myGrid(i1+index,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1         ))  +   fTempF[index][lbModel.DV_P1_M1_M1      ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1         ) = myGrid(i1+index,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1         ) + twoBeta*(lbModel.fTemp10[index][lbModel.G10_DV_P1_P1_M1         ] - myGrid(i1+index,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1         ))  +   fTempF[index][lbModel.DV_P1_P1_M1      ]*twoBetaTau  ;
                            myGrid(i1+index,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1         ) = myGrid(i1+index,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1         ) + twoBeta*(lbModel.fTemp10[index][lbModel.G10_DV_M1_P1_M1         ] - myGrid(i1+index,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1         ))  +   fTempF[index][lbModel.DV_M1_P1_M1      ]*twoBetaTau  ;
                        }
                    }




}

template<typename dataType1>
void getfTempFSIMD(lbmRD3Q41<dataType1> &lbModel,
                   dataType1 *rho,
                   dataType1 F1,
                   dataType1 F2,
                   dataType1 F3) {
  SIMD_REG _RHOTHETA0, _temp1, _temp2, _temp3;
  SIMD_REG _F1, _F2, _F3, _twoReg, _temp4, _temp5;
  SIMD_REG _temp6;
  dataType1 *pointer1;
  alignas(32) dataType1 tempArray[4];

  pointer1 = &rho[0];
  _RHOTHETA0 = LOAD_PD(pointer1);
  _temp1 = SET1_PD(lbModel.oneByTheta0);
  _RHOTHETA0 = MUL_PD(_RHOTHETA0, _temp1);

  _F1 = SET1_PD(F1);
  _F2 = SET1_PD(F2);
  _F3 = SET1_PD(F3);
  _twoReg = SET1_PD(2.0);

  lbModel.forceTemp0[0][0] = 0.0;
  lbModel.forceTemp0[1][0] = 0.0;
  lbModel.forceTemp0[2][0] = 0.0;
  lbModel.forceTemp0[3][0] = 0.0;

  // 1 Mag P and 1 Mag M
  _temp1 = SET1_PD(lbModel.wt[lbModel.DV_ZERO_ZERO_P1]);
  _temp1 = MUL_PD(_RHOTHETA0, _temp1);    // wt*rho*oneByTheta0

  pointer1 = &tempArray[0];

  _temp2 = MUL_PD(_temp1, _F3);
  STORE_PD(pointer1, _temp2);
  lbModel.forceTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1] = tempArray[0];
  lbModel.forceTemp1[1][lbModel.G1_DV_ZERO_ZERO_P1] = tempArray[1];
  lbModel.forceTemp1[2][lbModel.G1_DV_ZERO_ZERO_P1] = tempArray[2];
  lbModel.forceTemp1[3][lbModel.G1_DV_ZERO_ZERO_P1] = tempArray[3];

  _temp2 = MUL_PD(_temp1, _F2);
  STORE_PD(pointer1, _temp2);
  lbModel.forceTemp5[0][lbModel.G5_DV_ZERO_P1_ZERO] = tempArray[0];
  lbModel.forceTemp5[1][lbModel.G5_DV_ZERO_P1_ZERO] = tempArray[1];
  lbModel.forceTemp5[2][lbModel.G5_DV_ZERO_P1_ZERO] = tempArray[2];
  lbModel.forceTemp5[3][lbModel.G5_DV_ZERO_P1_ZERO] = tempArray[3];

  _temp2 = MUL_PD(_temp1, _F1);
  STORE_PD(pointer1, _temp2);
  lbModel.forceTemp5[0][lbModel.G5_DV_P1_ZERO_ZERO] = tempArray[0];
  lbModel.forceTemp5[1][lbModel.G5_DV_P1_ZERO_ZERO] = tempArray[1];
  lbModel.forceTemp5[2][lbModel.G5_DV_P1_ZERO_ZERO] = tempArray[2];
  lbModel.forceTemp5[3][lbModel.G5_DV_P1_ZERO_ZERO] = tempArray[3];

  _temp3 = SET1_PD(-1.0);
  _temp1 = MUL_PD(_temp3, _temp1);     // -1.0 * wt*rho*oneByTheta0

  _temp2 = MUL_PD(_temp1, _F3);
  STORE_PD(pointer1, _temp2);// G2_DV_ZERO_ZERO_M1
  lbModel.forceTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1] = tempArray[0];
  lbModel.forceTemp2[1][lbModel.G2_DV_ZERO_ZERO_M1] = tempArray[1];
  lbModel.forceTemp2[2][lbModel.G2_DV_ZERO_ZERO_M1] = tempArray[2];
  lbModel.forceTemp2[3][lbModel.G2_DV_ZERO_ZERO_M1] = tempArray[3];

  _temp2 = MUL_PD(_temp1, _F2);
  STORE_PD(pointer1, _temp2);//G6_DV_ZERO_M1_ZERO
  lbModel.forceTemp6[0][lbModel.G6_DV_ZERO_M1_ZERO] = tempArray[0];
  lbModel.forceTemp6[1][lbModel.G6_DV_ZERO_M1_ZERO] = tempArray[1];
  lbModel.forceTemp6[2][lbModel.G6_DV_ZERO_M1_ZERO] = tempArray[2];
  lbModel.forceTemp6[3][lbModel.G6_DV_ZERO_M1_ZERO] = tempArray[3];

  _temp2 = MUL_PD(_temp1, _F1);
  STORE_PD(pointer1, _temp2);//G6_DV_M1_ZERO_ZERO
  lbModel.forceTemp6[0][lbModel.G6_DV_M1_ZERO_ZERO] = tempArray[0];
  lbModel.forceTemp6[1][lbModel.G6_DV_M1_ZERO_ZERO] = tempArray[1];
  lbModel.forceTemp6[2][lbModel.G6_DV_M1_ZERO_ZERO] = tempArray[2];
  lbModel.forceTemp6[3][lbModel.G6_DV_M1_ZERO_ZERO] = tempArray[3];

  // 2 Mag P and 2 Mag M
  _temp1 = SET1_PD(lbModel.wt[lbModel.DV_ZERO_ZERO_P2]);
  _temp1 = MUL_PD(_RHOTHETA0, _temp1);    // wt*rho*oneByTheta0
  _temp1 = MUL_PD(_twoReg, _temp1);       // 2.0*wt*rho*oneByTheta0

  _temp2 = MUL_PD(_temp1, _F3);
  STORE_PD(pointer1, _temp2);
  lbModel.forceTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2] = tempArray[0];
  lbModel.forceTemp1[1][lbModel.G1_DV_ZERO_ZERO_P2] = tempArray[1];
  lbModel.forceTemp1[2][lbModel.G1_DV_ZERO_ZERO_P2] = tempArray[2];
  lbModel.forceTemp1[3][lbModel.G1_DV_ZERO_ZERO_P2] = tempArray[3];

  _temp2 = MUL_PD(_temp1, _F2);
  STORE_PD(pointer1, _temp2);
  lbModel.forceTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO] = tempArray[0];
  lbModel.forceTemp1[1][lbModel.G1_DV_ZERO_P2_ZERO] = tempArray[1];
  lbModel.forceTemp1[2][lbModel.G1_DV_ZERO_P2_ZERO] = tempArray[2];
  lbModel.forceTemp1[3][lbModel.G1_DV_ZERO_P2_ZERO] = tempArray[3];

  _temp2 = MUL_PD(_temp1, _F1);
  STORE_PD(pointer1, _temp2);
  lbModel.forceTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO] = tempArray[0];
  lbModel.forceTemp1[1][lbModel.G1_DV_P2_ZERO_ZERO] = tempArray[1];
  lbModel.forceTemp1[2][lbModel.G1_DV_P2_ZERO_ZERO] = tempArray[2];
  lbModel.forceTemp1[3][lbModel.G1_DV_P2_ZERO_ZERO] = tempArray[3];

  _temp3 = SET1_PD(-1.0);
  _temp1 = MUL_PD(_temp3, _temp1);     // -2.0 * wt*rho*oneByTheta0

  _temp2 = MUL_PD(_temp1, _F3);
  STORE_PD(pointer1, _temp2); //G2_DV_ZERO_ZERO_M2
  lbModel.forceTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2] = tempArray[0];
  lbModel.forceTemp2[1][lbModel.G2_DV_ZERO_ZERO_M2] = tempArray[1];
  lbModel.forceTemp2[2][lbModel.G2_DV_ZERO_ZERO_M2] = tempArray[2];
  lbModel.forceTemp2[3][lbModel.G2_DV_ZERO_ZERO_M2] = tempArray[3];

  _temp2 = MUL_PD(_temp1, _F2);
  STORE_PD(pointer1, _temp2); // G2_DV_ZERO_M2_ZERO
  lbModel.forceTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO] = tempArray[0];
  lbModel.forceTemp2[1][lbModel.G2_DV_ZERO_M2_ZERO] = tempArray[1];
  lbModel.forceTemp2[2][lbModel.G2_DV_ZERO_M2_ZERO] = tempArray[2];
  lbModel.forceTemp2[3][lbModel.G2_DV_ZERO_M2_ZERO] = tempArray[3];

  _temp2 = MUL_PD(_temp1, _F1);
  STORE_PD(pointer1, _temp2); // G2_DV_M2_ZERO_ZERO
  lbModel.forceTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO] = tempArray[0];
  lbModel.forceTemp2[1][lbModel.G2_DV_M2_ZERO_ZERO] = tempArray[1];
  lbModel.forceTemp2[2][lbModel.G2_DV_M2_ZERO_ZERO] = tempArray[2];
  lbModel.forceTemp2[3][lbModel.G2_DV_M2_ZERO_ZERO] = tempArray[3];

  // FCC
  _temp1 = SET1_PD(lbModel.wt[lbModel.DV_ZERO_P1_P1]);
  _temp1 = MUL_PD(_RHOTHETA0, _temp1);    // wt*rho*oneByTheta0

  _temp2 = ADD_PD(_F2, _F3);
  _temp2 = MUL_PD(_temp1, _temp2);
  STORE_PD(pointer1, _temp2); //G3_DV_ZERO_P1_P1
  lbModel.forceTemp3[0][lbModel.G3_DV_ZERO_P1_P1] = tempArray[0];
  lbModel.forceTemp3[1][lbModel.G3_DV_ZERO_P1_P1] = tempArray[1];
  lbModel.forceTemp3[2][lbModel.G3_DV_ZERO_P1_P1] = tempArray[2];
  lbModel.forceTemp3[3][lbModel.G3_DV_ZERO_P1_P1] = tempArray[3];

  _temp2 = SUB_PD(_F3, _F2);
  _temp2 = MUL_PD(_temp1, _temp2);
  STORE_PD(pointer1, _temp2); //G3_DV_ZERO_M1_P1
  lbModel.forceTemp3[0][lbModel.G3_DV_ZERO_M1_P1] = tempArray[0];
  lbModel.forceTemp3[1][lbModel.G3_DV_ZERO_M1_P1] = tempArray[1];
  lbModel.forceTemp3[2][lbModel.G3_DV_ZERO_M1_P1] = tempArray[2];
  lbModel.forceTemp3[3][lbModel.G3_DV_ZERO_M1_P1] = tempArray[3];

  _temp2 = ADD_PD(_F1, _F3);
  _temp2 = MUL_PD(_temp1, _temp2);
  STORE_PD(pointer1, _temp2); // G3_DV_P1_ZERO_P1
  lbModel.forceTemp3[0][lbModel.G3_DV_P1_ZERO_P1] = tempArray[0];
  lbModel.forceTemp3[1][lbModel.G3_DV_P1_ZERO_P1] = tempArray[1];
  lbModel.forceTemp3[2][lbModel.G3_DV_P1_ZERO_P1] = tempArray[2];
  lbModel.forceTemp3[3][lbModel.G3_DV_P1_ZERO_P1] = tempArray[3];

  _temp2 = SUB_PD(_F3, _F1);
  _temp2 = MUL_PD(_temp1, _temp2);
  STORE_PD(pointer1, _temp2); // G3_DV_M1_ZERO_P1
  lbModel.forceTemp3[0][lbModel.G3_DV_M1_ZERO_P1] = tempArray[0];
  lbModel.forceTemp3[1][lbModel.G3_DV_M1_ZERO_P1] = tempArray[1];
  lbModel.forceTemp3[2][lbModel.G3_DV_M1_ZERO_P1] = tempArray[2];
  lbModel.forceTemp3[3][lbModel.G3_DV_M1_ZERO_P1] = tempArray[3];

  _temp2 = ADD_PD(_F1, _F2);
  _temp2 = MUL_PD(_temp1, _temp2);
  STORE_PD(pointer1, _temp2); // G5_DV_P1_P1_ZERO
  lbModel.forceTemp5[0][lbModel.G5_DV_P1_P1_ZERO] = tempArray[0];
  lbModel.forceTemp5[1][lbModel.G5_DV_P1_P1_ZERO] = tempArray[1];
  lbModel.forceTemp5[2][lbModel.G5_DV_P1_P1_ZERO] = tempArray[2];
  lbModel.forceTemp5[3][lbModel.G5_DV_P1_P1_ZERO] = tempArray[3];

  _temp2 = SUB_PD(_F2, _F1);
  _temp2 = MUL_PD(_temp1, _temp2);
  STORE_PD(pointer1, _temp2); // G5_DV_M1_P1_ZERO
  lbModel.forceTemp5[0][lbModel.G5_DV_M1_P1_ZERO] = tempArray[0];
  lbModel.forceTemp5[1][lbModel.G5_DV_M1_P1_ZERO] = tempArray[1];
  lbModel.forceTemp5[2][lbModel.G5_DV_M1_P1_ZERO] = tempArray[2];
  lbModel.forceTemp5[3][lbModel.G5_DV_M1_P1_ZERO] = tempArray[3];

  _temp1 = MUL_PD(_temp3, _temp1);     // -1.0 * wt*rho*oneByTheta0

  _temp2 = ADD_PD(_F2, _F3);
  _temp2 = MUL_PD(_temp1, _temp2);
  STORE_PD(pointer1, _temp2); //G4_DV_ZERO_M1_M1
  lbModel.forceTemp4[0][lbModel.G4_DV_ZERO_M1_M1] = tempArray[0];
  lbModel.forceTemp4[1][lbModel.G4_DV_ZERO_M1_M1] = tempArray[1];
  lbModel.forceTemp4[2][lbModel.G4_DV_ZERO_M1_M1] = tempArray[2];
  lbModel.forceTemp4[3][lbModel.G4_DV_ZERO_M1_M1] = tempArray[3];

  _temp2 = SUB_PD(_F3, _F2);
  _temp2 = MUL_PD(_temp1, _temp2);
  STORE_PD(pointer1, _temp2); //G4_DV_ZERO_P1_M1
  lbModel.forceTemp4[0][lbModel.G4_DV_ZERO_P1_M1] = tempArray[0];
  lbModel.forceTemp4[1][lbModel.G4_DV_ZERO_P1_M1] = tempArray[1];
  lbModel.forceTemp4[2][lbModel.G4_DV_ZERO_P1_M1] = tempArray[2];
  lbModel.forceTemp4[3][lbModel.G4_DV_ZERO_P1_M1] = tempArray[3];

  _temp2 = ADD_PD(_F1, _F3);
  _temp2 = MUL_PD(_temp1, _temp2);
  STORE_PD(pointer1, _temp2); // G4_DV_M1_ZERO_M1
  lbModel.forceTemp4[0][lbModel.G4_DV_M1_ZERO_M1] = tempArray[0];
  lbModel.forceTemp4[1][lbModel.G4_DV_M1_ZERO_M1] = tempArray[1];
  lbModel.forceTemp4[2][lbModel.G4_DV_M1_ZERO_M1] = tempArray[2];
  lbModel.forceTemp4[3][lbModel.G4_DV_M1_ZERO_M1] = tempArray[3];

  _temp2 = SUB_PD(_F3, _F1);
  _temp2 = MUL_PD(_temp1, _temp2);
  STORE_PD(pointer1, _temp2); // G4_DV_P1_ZERO_M1
  lbModel.forceTemp4[0][lbModel.G4_DV_P1_ZERO_M1] = tempArray[0];
  lbModel.forceTemp4[1][lbModel.G4_DV_P1_ZERO_M1] = tempArray[1];
  lbModel.forceTemp4[2][lbModel.G4_DV_P1_ZERO_M1] = tempArray[2];
  lbModel.forceTemp4[3][lbModel.G4_DV_P1_ZERO_M1] = tempArray[3];

  _temp2 = ADD_PD(_F1, _F2);
  _temp2 = MUL_PD(_temp1, _temp2);
  STORE_PD(pointer1, _temp2); // G6_DV_M1_M1_ZERO
  lbModel.forceTemp6[0][lbModel.G6_DV_M1_M1_ZERO] = tempArray[0];
  lbModel.forceTemp6[1][lbModel.G6_DV_M1_M1_ZERO] = tempArray[1];
  lbModel.forceTemp6[2][lbModel.G6_DV_M1_M1_ZERO] = tempArray[2];
  lbModel.forceTemp6[3][lbModel.G6_DV_M1_M1_ZERO] = tempArray[3];

  _temp2 = SUB_PD(_F2, _F1);
  _temp2 = MUL_PD(_temp1, _temp2);
  STORE_PD(pointer1, _temp2); // G6_DV_P1_M1_ZERO
  lbModel.forceTemp6[0][lbModel.G6_DV_P1_M1_ZERO] = tempArray[0];
  lbModel.forceTemp6[1][lbModel.G6_DV_P1_M1_ZERO] = tempArray[1];
  lbModel.forceTemp6[2][lbModel.G6_DV_P1_M1_ZERO] = tempArray[2];
  lbModel.forceTemp6[3][lbModel.G6_DV_P1_M1_ZERO] = tempArray[3];

  // BCC and BCC1
  _temp5 = SET1_PD(0.5);
  _temp1 = SET1_PD(lbModel.wt[lbModel.DV_P_P_P]);
  _temp1 = MUL_PD(_RHOTHETA0, _temp1);//  wt*rho*oneByTheta0 - BCC
  _temp1 = MUL_PD(_temp5, _temp1);//  0.5*wt*rho*oneByTheta0 - BCC

  _temp4 = SET1_PD(lbModel.wt[lbModel.DV_P1_P1_P1]);
  _temp4 = MUL_PD(_RHOTHETA0, _temp4);//  wt*rho*oneByTheta0 - BCC1

  _temp2 = ADD_PD(_F1, _F2);
  _temp2 = ADD_PD(_temp2, _F3);

  _temp6 = MUL_PD(_temp2, _temp4);
  STORE_PD(pointer1, _temp6); // G9_DV_P1_P1_P1
  lbModel.forceTemp9[0][lbModel.G9_DV_P1_P1_P1] = tempArray[0];
  lbModel.forceTemp9[1][lbModel.G9_DV_P1_P1_P1] = tempArray[1];
  lbModel.forceTemp9[2][lbModel.G9_DV_P1_P1_P1] = tempArray[2];
  lbModel.forceTemp9[3][lbModel.G9_DV_P1_P1_P1] = tempArray[3];

  _temp2 = MUL_PD(_temp2, _temp1);
  STORE_PD(pointer1, _temp2); // G7_DV_P_P_P
  lbModel.forceTemp7[0][lbModel.G7_DV_P_P_P] = tempArray[0];
  lbModel.forceTemp7[1][lbModel.G7_DV_P_P_P] = tempArray[1];
  lbModel.forceTemp7[2][lbModel.G7_DV_P_P_P] = tempArray[2];
  lbModel.forceTemp7[3][lbModel.G7_DV_P_P_P] = tempArray[3];

  _temp2 = SUB_PD(_F2, _F1);
  _temp2 = ADD_PD(_temp2, _F3);

  _temp6 = MUL_PD(_temp2, _temp4);
  STORE_PD(pointer1, _temp6); // G9_DV_M1_P1_P1
  lbModel.forceTemp9[0][lbModel.G9_DV_M1_P1_P1] = tempArray[0];
  lbModel.forceTemp9[1][lbModel.G9_DV_M1_P1_P1] = tempArray[1];
  lbModel.forceTemp9[2][lbModel.G9_DV_M1_P1_P1] = tempArray[2];
  lbModel.forceTemp9[3][lbModel.G9_DV_M1_P1_P1] = tempArray[3];

  _temp2 = MUL_PD(_temp2, _temp1);
  STORE_PD(pointer1, _temp2); //  G7_DV_M_P_P
  lbModel.forceTemp7[0][lbModel.G7_DV_M_P_P] = tempArray[0];
  lbModel.forceTemp7[1][lbModel.G7_DV_M_P_P] = tempArray[1];
  lbModel.forceTemp7[2][lbModel.G7_DV_M_P_P] = tempArray[2];
  lbModel.forceTemp7[3][lbModel.G7_DV_M_P_P] = tempArray[3];

  _temp2 = SUB_PD(_F3, _F2);
  _temp2 = SUB_PD(_temp2, _F1);

  _temp6 = MUL_PD(_temp2, _temp4);
  STORE_PD(pointer1, _temp6); // G9_DV_M1_M1_P1
  lbModel.forceTemp9[0][lbModel.G9_DV_M1_M1_P1] = tempArray[0];
  lbModel.forceTemp9[1][lbModel.G9_DV_M1_M1_P1] = tempArray[1];
  lbModel.forceTemp9[2][lbModel.G9_DV_M1_M1_P1] = tempArray[2];
  lbModel.forceTemp9[3][lbModel.G9_DV_M1_M1_P1] = tempArray[3];

  _temp2 = MUL_PD(_temp2, _temp1);
  STORE_PD(pointer1, _temp2); // G7_DV_M_M_P
  lbModel.forceTemp7[0][lbModel.G7_DV_M_M_P] = tempArray[0];
  lbModel.forceTemp7[1][lbModel.G7_DV_M_M_P] = tempArray[1];
  lbModel.forceTemp7[2][lbModel.G7_DV_M_M_P] = tempArray[2];
  lbModel.forceTemp7[3][lbModel.G7_DV_M_M_P] = tempArray[3];

  _temp2 = SUB_PD(_F1, _F2);
  _temp2 = ADD_PD(_temp2, _F3);

  _temp6 = MUL_PD(_temp2, _temp4);
  STORE_PD(pointer1, _temp6);  // G9_DV_P1_M1_P1
  lbModel.forceTemp9[0][lbModel.G9_DV_P1_M1_P1] = tempArray[0];
  lbModel.forceTemp9[1][lbModel.G9_DV_P1_M1_P1] = tempArray[1];
  lbModel.forceTemp9[2][lbModel.G9_DV_P1_M1_P1] = tempArray[2];
  lbModel.forceTemp9[3][lbModel.G9_DV_P1_M1_P1] = tempArray[3];

  _temp2 = MUL_PD(_temp2, _temp1);
  STORE_PD(pointer1, _temp2);   // G7_DV_P_M_P
  lbModel.forceTemp7[0][lbModel.G7_DV_P_M_P] = tempArray[0];
  lbModel.forceTemp7[1][lbModel.G7_DV_P_M_P] = tempArray[1];
  lbModel.forceTemp7[2][lbModel.G7_DV_P_M_P] = tempArray[2];
  lbModel.forceTemp7[3][lbModel.G7_DV_P_M_P] = tempArray[3];

  _temp1 = MUL_PD(_temp3, _temp1);     // -0.5 * wt*rho*oneByTheta0
  _temp4 = MUL_PD(_temp3, _temp4);     // -1.0 * wt*rho*oneByTheta0

  _temp2 = ADD_PD(_F1, _F2);
  _temp2 = ADD_PD(_temp2, _F3);

  _temp6 = MUL_PD(_temp2, _temp4);
  STORE_PD(pointer1, _temp6);    //G10_DV_M1_M1_M1
  lbModel.forceTemp10[0][lbModel.G10_DV_M1_M1_M1] = tempArray[0];
  lbModel.forceTemp10[1][lbModel.G10_DV_M1_M1_M1] = tempArray[1];
  lbModel.forceTemp10[2][lbModel.G10_DV_M1_M1_M1] = tempArray[2];
  lbModel.forceTemp10[3][lbModel.G10_DV_M1_M1_M1] = tempArray[3];

  _temp2 = MUL_PD(_temp2, _temp1);
  STORE_PD(pointer1, _temp2);  // G8_DV_M_M_M
  lbModel.forceTemp8[0][lbModel.G8_DV_M_M_M] = tempArray[0];
  lbModel.forceTemp8[1][lbModel.G8_DV_M_M_M] = tempArray[1];
  lbModel.forceTemp8[2][lbModel.G8_DV_M_M_M] = tempArray[2];
  lbModel.forceTemp8[3][lbModel.G8_DV_M_M_M] = tempArray[3];

  _temp2 = SUB_PD(_F2, _F1);
  _temp2 = ADD_PD(_temp2, _F3);

  _temp6 = MUL_PD(_temp2, _temp4);
  STORE_PD(pointer1, _temp6); // G10_DV_P1_M1_M1
  lbModel.forceTemp10[0][lbModel.G10_DV_P1_M1_M1] = tempArray[0];
  lbModel.forceTemp10[1][lbModel.G10_DV_P1_M1_M1] = tempArray[1];
  lbModel.forceTemp10[2][lbModel.G10_DV_P1_M1_M1] = tempArray[2];
  lbModel.forceTemp10[3][lbModel.G10_DV_P1_M1_M1] = tempArray[3];

  _temp2 = MUL_PD(_temp2, _temp1);
  STORE_PD(pointer1, _temp2); // G8_DV_P_M_M
  lbModel.forceTemp8[0][lbModel.G8_DV_P_M_M] = tempArray[0];
  lbModel.forceTemp8[1][lbModel.G8_DV_P_M_M] = tempArray[1];
  lbModel.forceTemp8[2][lbModel.G8_DV_P_M_M] = tempArray[2];
  lbModel.forceTemp8[3][lbModel.G8_DV_P_M_M] = tempArray[3];

  _temp2 = SUB_PD(_F3, _F2);
  _temp2 = SUB_PD(_temp2, _F1);

  _temp6 = MUL_PD(_temp2, _temp4);
  STORE_PD(pointer1, _temp6); // G10_DV_P1_P1_M1
  lbModel.forceTemp10[0][lbModel.G10_DV_P1_P1_M1] = tempArray[0];
  lbModel.forceTemp10[1][lbModel.G10_DV_P1_P1_M1] = tempArray[1];
  lbModel.forceTemp10[2][lbModel.G10_DV_P1_P1_M1] = tempArray[2];
  lbModel.forceTemp10[3][lbModel.G10_DV_P1_P1_M1] = tempArray[3];

  _temp2 = MUL_PD(_temp2, _temp1);
  STORE_PD(pointer1, _temp2);  // G8_DV_P_P_M
  lbModel.forceTemp8[0][lbModel.G8_DV_P_P_M] = tempArray[0];
  lbModel.forceTemp8[1][lbModel.G8_DV_P_P_M] = tempArray[1];
  lbModel.forceTemp8[2][lbModel.G8_DV_P_P_M] = tempArray[2];
  lbModel.forceTemp8[3][lbModel.G8_DV_P_P_M] = tempArray[3];

  _temp2 = SUB_PD(_F1, _F2);
  _temp2 = ADD_PD(_temp2, _F3);

  _temp6 = MUL_PD(_temp2, _temp4);
  STORE_PD(pointer1, _temp6); // G10_DV_M1_P1_M1
  lbModel.forceTemp10[0][lbModel.G10_DV_M1_P1_M1] = tempArray[0];
  lbModel.forceTemp10[1][lbModel.G10_DV_M1_P1_M1] = tempArray[1];
  lbModel.forceTemp10[2][lbModel.G10_DV_M1_P1_M1] = tempArray[2];
  lbModel.forceTemp10[3][lbModel.G10_DV_M1_P1_M1] = tempArray[3];

  _temp2 = MUL_PD(_temp2, _temp1);
  STORE_PD(pointer1, _temp2); // G8_DV_M_P_M
  lbModel.forceTemp8[0][lbModel.G8_DV_M_P_M] = tempArray[0];
  lbModel.forceTemp8[1][lbModel.G8_DV_M_P_M] = tempArray[1];
  lbModel.forceTemp8[2][lbModel.G8_DV_M_P_M] = tempArray[2];
  lbModel.forceTemp8[3][lbModel.G8_DV_M_P_M] = tempArray[3];
}

template <int N,int numblock, typename dataType1>
void collideWithForceSIMD(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH, dataType1 twoBeta, dataType1 tau , dataType1 F1, dataType1 F2, dataType1 F3,dataType1 dt)
{
    dataType1 fTempF[4][41] ;

    dataType1 rho[VECT_LENGTH] __attribute__ ((aligned(32)));
    dataType1 uX[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 uY[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 uZ[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 theta[VECT_LENGTH] __attribute__ ((aligned(32)));

    dataType1 twoBetaTau =  twoBeta*tau  ;
    dataType1 oneByTheta0 = 1.0/lbModel.theta0;

    SIMD_REG _temp1,_temp2,_temp3;
    SIMD_REG _temp4,_temp5,_temp6;
    SIMD_REG _temp7,_temp8;

    dataType1 oneMinustwoBeta = 1.0-twoBeta;


    _temp1 = SET1_PD(twoBetaTau);         // twoBetaTau
    _temp2 = SET1_PD(twoBeta);            // twoBeta
    _temp3 = SET1_PD(oneMinustwoBeta);    // oneMinustwoBeta

    dataType1 *pointer1,*pointer2,*pointer3;

    for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1=i1+4)
            {
                getHydroMomentsFromNodeWithForce(lbModel,myGrid,VECT_LENGTH,i1,i2,i3,rho,uX,uY,uZ,theta,F1,F2,F3,dt);

  //                 for(int index=0;index<VECT_LENGTH;index++)
  //                     theta[index] = lbModel.theta0;

                getFEqSIMD(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);
                getfTempFSIMD(lbModel,rho,F1,F2,F3);

                ///////////
                // CENTER /
                ///////////

                for(int index=0;index<VECT_LENGTH;index++)
                 myGrid(i1+index,i2,i3,lbModel.G0 ,myGrid.node[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO) = myGrid(i1+index,i2,i3,lbModel.G0 ,myGrid.node[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO) + twoBeta*(lbModel.fTemp0 [index][lbModel.CENTER_DV_ZERO_ZERO_ZERO] - myGrid(i1+index,i2,i3,lbModel.G0 ,myGrid.node[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO)) ;

                /////////
                // G1 //
                /////////
                //G1 - i
                pointer1 = &lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1);
                pointer3 = &fTempF[0][lbModel.DV_ZERO_ZERO_P1];
                _temp5 = LOAD_PD(pointer1);           // fEq
                _temp6 = LOAD_PD(pointer2);           // G1
                _temp7 = LOAD_PD(pointer3);           // fTempF
                _temp5 = MUL_PD(_temp2,_temp5);       // twoBeta*fEq
                _temp6 = MUL_PD(_temp3,_temp6);       // oneMinustwoBeta*SC1
                _temp7 = MUL_PD(_temp1,_temp7);       // twoBetaTau*fTempF
                _temp8 = ADD_PD(_temp5,_temp6);       // oneMinustwoBeta*SC1  + twoBeta*fEq
                _temp8 = ADD_PD(_temp7,_temp8);       // oneMinustwoBeta*SC1  + twoBeta*fEq + twoBetaTau*fTempF
                STORE_PD(pointer2,_temp8);            // write it back to G1

                //G1-i+1
                pointer1 = &lbModel.fTemp1[1][lbModel.G1_DV_ZERO_ZERO_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_ZERO_ZERO_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                //G1-i+2
                pointer1 = &lbModel.fTemp1[2][lbModel.G1_DV_ZERO_ZERO_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_ZERO_ZERO_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                //G1-i+2
                pointer1 = &lbModel.fTemp1[3][lbModel.G1_DV_ZERO_ZERO_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_ZERO_ZERO_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);



                /////////
                // G2 //
                /////////
                //G2 - i
                pointer1 = &lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1);
                pointer3 = &fTempF[0][lbModel.DV_ZERO_ZERO_M1];
                _temp5 = LOAD_PD(pointer1);           // fEq
                _temp6 = LOAD_PD(pointer2);           // G2
                _temp7 = LOAD_PD(pointer3);           // fTempF
                _temp5 = MUL_PD(_temp2,_temp5);       // twoBeta*fEq
                _temp6 = MUL_PD(_temp3,_temp6);       // oneMinustwoBeta*SC2
                _temp7 = MUL_PD(_temp1,_temp7);       // twoBetaTau*fTempF
                _temp8 = ADD_PD(_temp5,_temp6);       // oneMinustwoBeta*SC2  + twoBeta*fEq
                _temp8 = ADD_PD(_temp7,_temp8);       // oneMinustwoBeta*SC2  + twoBeta*fEq + twoBetaTau*fTempF
                STORE_PD(pointer2,_temp8);            // write it back to G2

                //G2-i+1
                pointer1 = &lbModel.fTemp2[1][lbModel.G2_DV_ZERO_ZERO_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_ZERO_ZERO_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                //G2-i+2
                pointer1 = &lbModel.fTemp2[2][lbModel.G2_DV_ZERO_ZERO_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_ZERO_ZERO_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                //G2-i+3
                pointer1 = &lbModel.fTemp2[3][lbModel.G2_DV_ZERO_ZERO_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_ZERO_ZERO_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);


                /////////
                // G3 //
                /////////
                //G3 - i
                pointer1 = &lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1);
                pointer3 = &fTempF[0][lbModel.DV_ZERO_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp3[1][lbModel.G3_DV_ZERO_P1_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_ZERO_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp3[2][lbModel.G3_DV_ZERO_P1_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_ZERO_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp3[3][lbModel.G3_DV_ZERO_P1_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_ZERO_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);


                /////////
                // G4 //
                /////////
                //G4 - i
                pointer1 = &lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1);
                pointer3 = &fTempF[0][lbModel.DV_ZERO_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp4[1][lbModel.G4_DV_ZERO_M1_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_ZERO_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp4[2][lbModel.G4_DV_ZERO_M1_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_ZERO_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp4[3][lbModel.G4_DV_ZERO_M1_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_ZERO_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                /////////
                // G5 //
                /////////
                //G5 - i
                pointer1 = &lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO);
                pointer3 = &fTempF[0][lbModel.DV_P1_P1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp5[1][lbModel.G5_DV_P1_P1_ZERO];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_P1_P1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp5[2][lbModel.G5_DV_P1_P1_ZERO];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_P1_P1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp5[3][lbModel.G5_DV_P1_P1_ZERO];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_P1_P1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                /////////
                // G6 //
                /////////
                //G6 - i
                pointer1 = &lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO);
                pointer3 = &fTempF[0][lbModel.DV_M1_M1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp6[1][lbModel.G6_DV_M1_M1_ZERO];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_M1_M1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp6[2][lbModel.G6_DV_M1_M1_ZERO];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_M1_M1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp6[3][lbModel.G6_DV_M1_M1_ZERO];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_M1_M1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                /////////
                // G7 //
                /////////
                //G7 - i
                pointer1 = &lbModel.fTemp7[0][lbModel.G7_DV_P_P_P];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);
                pointer3 = &fTempF[0][lbModel.DV_P_P_P];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp7[1][lbModel.G7_DV_P_P_P];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_P_P_P];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp7[2][lbModel.G7_DV_P_P_P];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_P_P_P];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp7[3][lbModel.G7_DV_P_P_P];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_P_P_P];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                /////////
                // G8 //
                /////////
                //G8 - i
                pointer1 = &lbModel.fTemp8[0][lbModel.G8_DV_M_M_M];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);
                pointer3 = &fTempF[0][lbModel.DV_M_M_M];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp8[1][lbModel.G8_DV_M_M_M];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_M_M_M];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp8[2][lbModel.G8_DV_M_M_M];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_M_M_M];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp8[3][lbModel.G8_DV_M_M_M];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_M_M_M];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);


                /////////
                // G9 //
                /////////
                //G9 - i
                pointer1 = &lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_P1_P1);
                pointer3 = &fTempF[0][lbModel.DV_P1_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp9[1][lbModel.G9_DV_P1_P1_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_P1_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp9[2][lbModel.G9_DV_P1_P1_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_P1_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp9[3][lbModel.G9_DV_P1_P1_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_P1_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);


                /////////
                // G10 //
                /////////
                //G10 - i
                pointer1 = &lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1);
                pointer3 = &fTempF[0][lbModel.DV_M1_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp10[1][lbModel.G10_DV_M1_M1_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_M1_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp10[2][lbModel.G10_DV_M1_M1_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_M1_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp10[3][lbModel.G10_DV_M1_M1_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_M1_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);


            }


    for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1=i1+4)
            {
                getHydroMomentsFromCellWithForce(lbModel,myGrid,VECT_LENGTH,i1,i2,i3,rho,uX,uY,uZ,theta,F1,F2,F3,dt);

  /*                for(int index=0;index<VECT_LENGTH;index++)
                    theta[index] = lbModel.theta0;*/

                getFEqSIMD(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta);
                getfTempFSIMD(lbModel,rho,F1,F2,F3);

                ///////////
                // CENTER /
                ///////////

                for(int index=0;index<VECT_LENGTH;index++)
                 myGrid(i1+index,i2,i3,lbModel.G0 ,myGrid.cell[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO) = myGrid(i1+index,i2,i3,lbModel.G0 ,myGrid.cell[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO) + twoBeta*(lbModel.fTemp0 [index][lbModel.CENTER_DV_ZERO_ZERO_ZERO] - myGrid(i1+index,i2,i3,lbModel.G0 ,myGrid.cell[lbModel.G0 ],lbModel.CENTER_DV_ZERO_ZERO_ZERO)) ;

                /////////
                // G1 //
                /////////
                //G1 - i
                pointer1 = &lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1);
                pointer3 = &fTempF[0][lbModel.DV_ZERO_ZERO_P1];
                _temp5 = LOAD_PD(pointer1);           // fEq
                _temp6 = LOAD_PD(pointer2);           // G1
                _temp7 = LOAD_PD(pointer3);           // fTempF
                _temp5 = MUL_PD(_temp2,_temp5);       // twoBeta*fEq
                _temp6 = MUL_PD(_temp3,_temp6);       // oneMinustwoBeta*SC1
                _temp7 = MUL_PD(_temp1,_temp7);       // twoBetaTau*fTempF
                _temp8 = ADD_PD(_temp5,_temp6);       // oneMinustwoBeta*SC1  + twoBeta*fEq
                _temp8 = ADD_PD(_temp7,_temp8);       // oneMinustwoBeta*SC1  + twoBeta*fEq + twoBetaTau*fTempF
                STORE_PD(pointer2,_temp8);            // write it back to G1

                //G1-i+1
                pointer1 = &lbModel.fTemp1[1][lbModel.G1_DV_ZERO_ZERO_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_ZERO_ZERO_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                //G1-i+2
                pointer1 = &lbModel.fTemp1[2][lbModel.G1_DV_ZERO_ZERO_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_ZERO_ZERO_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                //G1-i+2
                pointer1 = &lbModel.fTemp1[3][lbModel.G1_DV_ZERO_ZERO_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_ZERO_ZERO_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);



                /////////
                // G2 //
                /////////
                //G2 - i
                pointer1 = &lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1);
                pointer3 = &fTempF[0][lbModel.DV_ZERO_ZERO_M1];
                _temp5 = LOAD_PD(pointer1);           // fEq
                _temp6 = LOAD_PD(pointer2);           // G2
                _temp7 = LOAD_PD(pointer3);           // fTempF
                _temp5 = MUL_PD(_temp2,_temp5);       // twoBeta*fEq
                _temp6 = MUL_PD(_temp3,_temp6);       // oneMinustwoBeta*SC2
                _temp7 = MUL_PD(_temp1,_temp7);       // twoBetaTau*fTempF
                _temp8 = ADD_PD(_temp5,_temp6);       // oneMinustwoBeta*SC2  + twoBeta*fEq
                _temp8 = ADD_PD(_temp7,_temp8);       // oneMinustwoBeta*SC2  + twoBeta*fEq + twoBetaTau*fTempF
                STORE_PD(pointer2,_temp8);            // write it back to G2

                //G2-i+1
                pointer1 = &lbModel.fTemp2[1][lbModel.G2_DV_ZERO_ZERO_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_ZERO_ZERO_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                //G2-i+2
                pointer1 = &lbModel.fTemp2[2][lbModel.G2_DV_ZERO_ZERO_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_ZERO_ZERO_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                //G2-i+3
                pointer1 = &lbModel.fTemp2[3][lbModel.G2_DV_ZERO_ZERO_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_ZERO_ZERO_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);


                /////////
                // G3 //
                /////////
                //G3 - i
                pointer1 = &lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1);
                pointer3 = &fTempF[0][lbModel.DV_ZERO_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp3[1][lbModel.G3_DV_ZERO_P1_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_ZERO_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp3[2][lbModel.G3_DV_ZERO_P1_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_ZERO_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp3[3][lbModel.G3_DV_ZERO_P1_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_ZERO_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);


                /////////
                // G4 //
                /////////
                //G4 - i
                pointer1 = &lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1);
                pointer3 = &fTempF[0][lbModel.DV_ZERO_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp4[1][lbModel.G4_DV_ZERO_M1_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_ZERO_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp4[2][lbModel.G4_DV_ZERO_M1_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_ZERO_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp4[3][lbModel.G4_DV_ZERO_M1_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_ZERO_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                /////////
                // G5 //
                /////////
                //G5 - i
                pointer1 = &lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO);
                pointer3 = &fTempF[0][lbModel.DV_P1_P1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp5[1][lbModel.G5_DV_P1_P1_ZERO];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_P1_P1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp5[2][lbModel.G5_DV_P1_P1_ZERO];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_P1_P1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp5[3][lbModel.G5_DV_P1_P1_ZERO];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_P1_P1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                /////////
                // G6 //
                /////////
                //G6 - i
                pointer1 = &lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO);
                pointer3 = &fTempF[0][lbModel.DV_M1_M1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp6[1][lbModel.G6_DV_M1_M1_ZERO];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_M1_M1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp6[2][lbModel.G6_DV_M1_M1_ZERO];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_M1_M1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp6[3][lbModel.G6_DV_M1_M1_ZERO];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_M1_M1_ZERO];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                /////////
                // G7 //
                /////////
                //G7 - i
                pointer1 = &lbModel.fTemp7[0][lbModel.G7_DV_P_P_P];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);
                pointer3 = &fTempF[0][lbModel.DV_P_P_P];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp7[1][lbModel.G7_DV_P_P_P];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_P_P_P];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp7[2][lbModel.G7_DV_P_P_P];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_P_P_P];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp7[3][lbModel.G7_DV_P_P_P];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_P_P_P];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                /////////
                // G8 //
                /////////
                //G8 - i
                pointer1 = &lbModel.fTemp8[0][lbModel.G8_DV_M_M_M];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);
                pointer3 = &fTempF[0][lbModel.DV_M_M_M];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp8[1][lbModel.G8_DV_M_M_M];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_M_M_M];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp8[2][lbModel.G8_DV_M_M_M];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_M_M_M];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp8[3][lbModel.G8_DV_M_M_M];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_M_M_M];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);


                /////////
                // G9 //
                /////////
                //G9 - i
                pointer1 = &lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_P1_P1);
                pointer3 = &fTempF[0][lbModel.DV_P1_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp9[1][lbModel.G9_DV_P1_P1_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_P1_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp9[2][lbModel.G9_DV_P1_P1_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_P1_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp9[3][lbModel.G9_DV_P1_P1_P1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_P1_P1_P1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);


                /////////
                // G10 //
                /////////
                //G10 - i
                pointer1 = &lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
                pointer2 = &myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1);
                pointer3 = &fTempF[0][lbModel.DV_M1_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp10[1][lbModel.G10_DV_M1_M1_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[1][lbModel.DV_M1_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp10[2][lbModel.G10_DV_M1_M1_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[2][lbModel.DV_M1_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

                pointer1 = &lbModel.fTemp10[3][lbModel.G10_DV_M1_M1_M1];
                pointer2 = pointer2+4;
                pointer3 = &fTempF[3][lbModel.DV_M1_M1_M1];
                _temp5 = LOAD_PD(pointer1);
                _temp6 = LOAD_PD(pointer2);
                _temp7 = LOAD_PD(pointer3);
                _temp5 = MUL_PD(_temp2,_temp5);
                _temp6 = MUL_PD(_temp3,_temp6);
                _temp7 = MUL_PD(_temp1,_temp7);
                _temp8 = ADD_PD(_temp5,_temp6);
                _temp8 = ADD_PD(_temp7,_temp8);
                STORE_PD(pointer2,_temp8);

            }


}

template <typename dataType1>
void copyFrom4To41Array(lbmRD3Q41<dataType1> &lbModel,double feq[][41])
{
  feq[0][lbModel.DV_ZERO_ZERO_ZERO]  = lbModel.fTemp0[0][lbModel.CENTER_DV_ZERO_ZERO_ZERO];
  feq[1][lbModel.DV_ZERO_ZERO_ZERO]  = lbModel.fTemp0[1][lbModel.CENTER_DV_ZERO_ZERO_ZERO];
  feq[2][lbModel.DV_ZERO_ZERO_ZERO]  = lbModel.fTemp0[2][lbModel.CENTER_DV_ZERO_ZERO_ZERO];
  feq[3][lbModel.DV_ZERO_ZERO_ZERO]  = lbModel.fTemp0[3][lbModel.CENTER_DV_ZERO_ZERO_ZERO];

  feq[0][lbModel.DV_ZERO_ZERO_P1]    = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1];
  feq[0][lbModel.DV_ZERO_ZERO_P2]    = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2];
  feq[0][lbModel.DV_ZERO_P2_ZERO]    = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO];
  feq[0][lbModel.DV_P2_ZERO_ZERO]    = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO];

  feq[1][lbModel.DV_ZERO_ZERO_P1]  = lbModel.fTemp1[1][lbModel.G1_DV_ZERO_ZERO_P1];
  feq[1][lbModel.DV_ZERO_ZERO_P2]  = lbModel.fTemp1[1][lbModel.G1_DV_ZERO_ZERO_P2];
  feq[1][lbModel.DV_ZERO_P2_ZERO]  = lbModel.fTemp1[1][lbModel.G1_DV_ZERO_P2_ZERO];
  feq[1][lbModel.DV_P2_ZERO_ZERO]  = lbModel.fTemp1[1][lbModel.G1_DV_P2_ZERO_ZERO];

  feq[2][lbModel.DV_ZERO_ZERO_P1]  = lbModel.fTemp1[2][lbModel.G1_DV_ZERO_ZERO_P1];
  feq[2][lbModel.DV_ZERO_ZERO_P2]  = lbModel.fTemp1[2][lbModel.G1_DV_ZERO_ZERO_P2];
  feq[2][lbModel.DV_ZERO_P2_ZERO]  = lbModel.fTemp1[2][lbModel.G1_DV_ZERO_P2_ZERO];
  feq[2][lbModel.DV_P2_ZERO_ZERO]  = lbModel.fTemp1[2][lbModel.G1_DV_P2_ZERO_ZERO];

  feq[3][lbModel.DV_ZERO_ZERO_P1]  = lbModel.fTemp1[3][lbModel.G1_DV_ZERO_ZERO_P1];
  feq[3][lbModel.DV_ZERO_ZERO_P2]  = lbModel.fTemp1[3][lbModel.G1_DV_ZERO_ZERO_P2];
  feq[3][lbModel.DV_ZERO_P2_ZERO]  = lbModel.fTemp1[3][lbModel.G1_DV_ZERO_P2_ZERO];
  feq[3][lbModel.DV_P2_ZERO_ZERO]  = lbModel.fTemp1[3][lbModel.G1_DV_P2_ZERO_ZERO];


  feq[0][lbModel.DV_ZERO_ZERO_M1]  = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];
  feq[0][lbModel.DV_ZERO_ZERO_M2]  = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
  feq[0][lbModel.DV_ZERO_M2_ZERO]  = lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
  feq[0][lbModel.DV_M2_ZERO_ZERO]  = lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];

  feq[1][lbModel.DV_ZERO_ZERO_M1]  = lbModel.fTemp2[1][lbModel.G2_DV_ZERO_ZERO_M1];
  feq[1][lbModel.DV_ZERO_ZERO_M2]  = lbModel.fTemp2[1][lbModel.G2_DV_ZERO_ZERO_M2];
  feq[1][lbModel.DV_ZERO_M2_ZERO]  = lbModel.fTemp2[1][lbModel.G2_DV_ZERO_M2_ZERO];
  feq[1][lbModel.DV_M2_ZERO_ZERO]  = lbModel.fTemp2[1][lbModel.G2_DV_M2_ZERO_ZERO];

  feq[2][lbModel.DV_ZERO_ZERO_M1]  = lbModel.fTemp2[2][lbModel.G2_DV_ZERO_ZERO_M1];
  feq[2][lbModel.DV_ZERO_ZERO_M2]  = lbModel.fTemp2[2][lbModel.G2_DV_ZERO_ZERO_M2];
  feq[2][lbModel.DV_ZERO_M2_ZERO]  = lbModel.fTemp2[2][lbModel.G2_DV_ZERO_M2_ZERO];
  feq[2][lbModel.DV_M2_ZERO_ZERO]  = lbModel.fTemp2[2][lbModel.G2_DV_M2_ZERO_ZERO];

  feq[3][lbModel.DV_ZERO_ZERO_M1]  = lbModel.fTemp2[3][lbModel.G2_DV_ZERO_ZERO_M1];
  feq[3][lbModel.DV_ZERO_ZERO_M2]  = lbModel.fTemp2[3][lbModel.G2_DV_ZERO_ZERO_M2];
  feq[3][lbModel.DV_ZERO_M2_ZERO]  = lbModel.fTemp2[3][lbModel.G2_DV_ZERO_M2_ZERO];
  feq[3][lbModel.DV_M2_ZERO_ZERO]  = lbModel.fTemp2[3][lbModel.G2_DV_M2_ZERO_ZERO];

  feq[0][lbModel.DV_ZERO_P1_P1]    = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1];
  feq[0][lbModel.DV_ZERO_M1_P1]    = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_M1_P1];
  feq[0][lbModel.DV_P1_ZERO_P1]    = lbModel.fTemp3[0][lbModel.G3_DV_P1_ZERO_P1];
  feq[0][lbModel.DV_M1_ZERO_P1]    = lbModel.fTemp3[0][lbModel.G3_DV_M1_ZERO_P1];

  feq[1][lbModel.DV_ZERO_P1_P1]    = lbModel.fTemp3[1][lbModel.G3_DV_ZERO_P1_P1];
  feq[1][lbModel.DV_ZERO_M1_P1]    = lbModel.fTemp3[1][lbModel.G3_DV_ZERO_M1_P1];
  feq[1][lbModel.DV_P1_ZERO_P1]    = lbModel.fTemp3[1][lbModel.G3_DV_P1_ZERO_P1];
  feq[1][lbModel.DV_M1_ZERO_P1]    = lbModel.fTemp3[1][lbModel.G3_DV_M1_ZERO_P1];

  feq[2][lbModel.DV_ZERO_P1_P1]    = lbModel.fTemp3[2][lbModel.G3_DV_ZERO_P1_P1];
  feq[2][lbModel.DV_ZERO_M1_P1]    = lbModel.fTemp3[2][lbModel.G3_DV_ZERO_M1_P1];
  feq[2][lbModel.DV_P1_ZERO_P1]    = lbModel.fTemp3[2][lbModel.G3_DV_P1_ZERO_P1];
  feq[2][lbModel.DV_M1_ZERO_P1]    = lbModel.fTemp3[2][lbModel.G3_DV_M1_ZERO_P1];

  feq[3][lbModel.DV_ZERO_P1_P1]    = lbModel.fTemp3[3][lbModel.G3_DV_ZERO_P1_P1];
  feq[3][lbModel.DV_ZERO_M1_P1]    = lbModel.fTemp3[3][lbModel.G3_DV_ZERO_M1_P1];
  feq[3][lbModel.DV_P1_ZERO_P1]    = lbModel.fTemp3[3][lbModel.G3_DV_P1_ZERO_P1];
  feq[3][lbModel.DV_M1_ZERO_P1]    = lbModel.fTemp3[3][lbModel.G3_DV_M1_ZERO_P1];

  feq[0][lbModel.DV_M1_ZERO_M1]    = lbModel.fTemp4[0][lbModel.G4_DV_M1_ZERO_M1];
  feq[0][lbModel.DV_ZERO_P1_M1]    = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_P1_M1];
  feq[0][lbModel.DV_ZERO_M1_M1]    = lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1];
  feq[0][lbModel.DV_P1_ZERO_M1]    = lbModel.fTemp4[0][lbModel.G4_DV_P1_ZERO_M1];

  feq[1][lbModel.DV_M1_ZERO_M1]    = lbModel.fTemp4[1][lbModel.G4_DV_M1_ZERO_M1];
  feq[1][lbModel.DV_ZERO_P1_M1]    = lbModel.fTemp4[1][lbModel.G4_DV_ZERO_P1_M1];
  feq[1][lbModel.DV_ZERO_M1_M1]    = lbModel.fTemp4[1][lbModel.G4_DV_ZERO_M1_M1];
  feq[1][lbModel.DV_P1_ZERO_M1]    = lbModel.fTemp4[1][lbModel.G4_DV_P1_ZERO_M1];

  feq[2][lbModel.DV_M1_ZERO_M1]    = lbModel.fTemp4[2][lbModel.G4_DV_M1_ZERO_M1];
  feq[2][lbModel.DV_ZERO_P1_M1]    = lbModel.fTemp4[2][lbModel.G4_DV_ZERO_P1_M1];
  feq[2][lbModel.DV_ZERO_M1_M1]    = lbModel.fTemp4[2][lbModel.G4_DV_ZERO_M1_M1];
  feq[2][lbModel.DV_P1_ZERO_M1]    = lbModel.fTemp4[2][lbModel.G4_DV_P1_ZERO_M1];

  feq[3][lbModel.DV_M1_ZERO_M1]    = lbModel.fTemp4[3][lbModel.G4_DV_M1_ZERO_M1];
  feq[3][lbModel.DV_ZERO_P1_M1]    = lbModel.fTemp4[3][lbModel.G4_DV_ZERO_P1_M1];
  feq[3][lbModel.DV_ZERO_M1_M1]    = lbModel.fTemp4[3][lbModel.G4_DV_ZERO_M1_M1];
  feq[3][lbModel.DV_P1_ZERO_M1]    = lbModel.fTemp4[3][lbModel.G4_DV_P1_ZERO_M1];

  feq[0][lbModel.DV_P1_P1_ZERO]    = lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO];
  feq[0][lbModel.DV_M1_P1_ZERO]    = lbModel.fTemp5[0][lbModel.G5_DV_M1_P1_ZERO];
  feq[0][lbModel.DV_ZERO_P1_ZERO]  = lbModel.fTemp5[0][lbModel.G5_DV_ZERO_P1_ZERO];
  feq[0][lbModel.DV_P1_ZERO_ZERO]  = lbModel.fTemp5[0][lbModel.G5_DV_P1_ZERO_ZERO];

  feq[1][lbModel.DV_P1_P1_ZERO]    = lbModel.fTemp5[1][lbModel.G5_DV_P1_P1_ZERO];
  feq[1][lbModel.DV_M1_P1_ZERO]    = lbModel.fTemp5[1][lbModel.G5_DV_M1_P1_ZERO];
  feq[1][lbModel.DV_ZERO_P1_ZERO]  = lbModel.fTemp5[1][lbModel.G5_DV_ZERO_P1_ZERO];
  feq[1][lbModel.DV_P1_ZERO_ZERO]  = lbModel.fTemp5[1][lbModel.G5_DV_P1_ZERO_ZERO];

  feq[2][lbModel.DV_P1_P1_ZERO]    = lbModel.fTemp5[2][lbModel.G5_DV_P1_P1_ZERO];
  feq[2][lbModel.DV_M1_P1_ZERO]    = lbModel.fTemp5[2][lbModel.G5_DV_M1_P1_ZERO];
  feq[2][lbModel.DV_ZERO_P1_ZERO]  = lbModel.fTemp5[2][lbModel.G5_DV_ZERO_P1_ZERO];
  feq[2][lbModel.DV_P1_ZERO_ZERO]  = lbModel.fTemp5[2][lbModel.G5_DV_P1_ZERO_ZERO];

  feq[3][lbModel.DV_P1_P1_ZERO]    = lbModel.fTemp5[3][lbModel.G5_DV_P1_P1_ZERO];
  feq[3][lbModel.DV_M1_P1_ZERO]    = lbModel.fTemp5[3][lbModel.G5_DV_M1_P1_ZERO];
  feq[3][lbModel.DV_ZERO_P1_ZERO]  = lbModel.fTemp5[3][lbModel.G5_DV_ZERO_P1_ZERO];
  feq[3][lbModel.DV_P1_ZERO_ZERO]  = lbModel.fTemp5[3][lbModel.G5_DV_P1_ZERO_ZERO];

  feq[0][lbModel.DV_M1_M1_ZERO]    = lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO];
  feq[0][lbModel.DV_P1_M1_ZERO]    = lbModel.fTemp6[0][lbModel.G6_DV_P1_M1_ZERO];
  feq[0][lbModel.DV_ZERO_M1_ZERO]  = lbModel.fTemp6[0][lbModel.G6_DV_ZERO_M1_ZERO];
  feq[0][lbModel.DV_M1_ZERO_ZERO]  = lbModel.fTemp6[0][lbModel.G6_DV_M1_ZERO_ZERO];

  feq[1][lbModel.DV_M1_M1_ZERO]    = lbModel.fTemp6[1][lbModel.G6_DV_M1_M1_ZERO];
  feq[1][lbModel.DV_P1_M1_ZERO]    = lbModel.fTemp6[1][lbModel.G6_DV_P1_M1_ZERO];
  feq[1][lbModel.DV_ZERO_M1_ZERO]  = lbModel.fTemp6[1][lbModel.G6_DV_ZERO_M1_ZERO];
  feq[1][lbModel.DV_M1_ZERO_ZERO]  = lbModel.fTemp6[1][lbModel.G6_DV_M1_ZERO_ZERO];

  feq[2][lbModel.DV_M1_M1_ZERO]    = lbModel.fTemp6[2][lbModel.G6_DV_M1_M1_ZERO];
  feq[2][lbModel.DV_P1_M1_ZERO]    = lbModel.fTemp6[2][lbModel.G6_DV_P1_M1_ZERO];
  feq[2][lbModel.DV_ZERO_M1_ZERO]  = lbModel.fTemp6[2][lbModel.G6_DV_ZERO_M1_ZERO];
  feq[2][lbModel.DV_M1_ZERO_ZERO]  = lbModel.fTemp6[2][lbModel.G6_DV_M1_ZERO_ZERO];

  feq[3][lbModel.DV_M1_M1_ZERO]    = lbModel.fTemp6[3][lbModel.G6_DV_M1_M1_ZERO];
  feq[3][lbModel.DV_P1_M1_ZERO]    = lbModel.fTemp6[3][lbModel.G6_DV_P1_M1_ZERO];
  feq[3][lbModel.DV_ZERO_M1_ZERO]  = lbModel.fTemp6[3][lbModel.G6_DV_ZERO_M1_ZERO];
  feq[3][lbModel.DV_M1_ZERO_ZERO]  = lbModel.fTemp6[3][lbModel.G6_DV_M1_ZERO_ZERO];

  feq[0][lbModel.DV_P_P_P]         = lbModel.fTemp7[0][lbModel.G7_DV_P_P_P];
  feq[0][lbModel.DV_M_P_P]         = lbModel.fTemp7[0][lbModel.G7_DV_M_P_P];
  feq[0][lbModel.DV_M_M_P]         = lbModel.fTemp7[0][lbModel.G7_DV_M_M_P];
  feq[0][lbModel.DV_P_M_P]         = lbModel.fTemp7[0][lbModel.G7_DV_P_M_P];

  feq[1][lbModel.DV_P_P_P]         = lbModel.fTemp7[1][lbModel.G7_DV_P_P_P];
  feq[1][lbModel.DV_M_P_P]         = lbModel.fTemp7[1][lbModel.G7_DV_M_P_P];
  feq[1][lbModel.DV_M_M_P]         = lbModel.fTemp7[1][lbModel.G7_DV_M_M_P];
  feq[1][lbModel.DV_P_M_P]         = lbModel.fTemp7[1][lbModel.G7_DV_P_M_P];

  feq[2][lbModel.DV_P_P_P]         = lbModel.fTemp7[2][lbModel.G7_DV_P_P_P];
  feq[2][lbModel.DV_M_P_P]         = lbModel.fTemp7[2][lbModel.G7_DV_M_P_P];
  feq[2][lbModel.DV_M_M_P]         = lbModel.fTemp7[2][lbModel.G7_DV_M_M_P];
  feq[2][lbModel.DV_P_M_P]         = lbModel.fTemp7[2][lbModel.G7_DV_P_M_P];

  feq[3][lbModel.DV_P_P_P]         = lbModel.fTemp7[3][lbModel.G7_DV_P_P_P];
  feq[3][lbModel.DV_M_P_P]         = lbModel.fTemp7[3][lbModel.G7_DV_M_P_P];
  feq[3][lbModel.DV_M_M_P]         = lbModel.fTemp7[3][lbModel.G7_DV_M_M_P];
  feq[3][lbModel.DV_P_M_P]         = lbModel.fTemp7[3][lbModel.G7_DV_P_M_P];

  feq[0][lbModel.DV_M_M_M]         = lbModel.fTemp8[0][lbModel.G8_DV_M_M_M];
  feq[0][lbModel.DV_M_P_M]         = lbModel.fTemp8[0][lbModel.G8_DV_M_P_M];
  feq[0][lbModel.DV_P_M_M]         = lbModel.fTemp8[0][lbModel.G8_DV_P_M_M];
  feq[0][lbModel.DV_P_P_M]         = lbModel.fTemp8[0][lbModel.G8_DV_P_P_M];

  feq[1][lbModel.DV_M_M_M]         = lbModel.fTemp8[1][lbModel.G8_DV_M_M_M];
  feq[1][lbModel.DV_M_P_M]         = lbModel.fTemp8[1][lbModel.G8_DV_M_P_M];
  feq[1][lbModel.DV_P_M_M]         = lbModel.fTemp8[1][lbModel.G8_DV_P_M_M];
  feq[1][lbModel.DV_P_P_M]         = lbModel.fTemp8[1][lbModel.G8_DV_P_P_M];

  feq[2][lbModel.DV_M_M_M]         = lbModel.fTemp8[2][lbModel.G8_DV_M_M_M];
  feq[2][lbModel.DV_M_P_M]         = lbModel.fTemp8[2][lbModel.G8_DV_M_P_M];
  feq[2][lbModel.DV_P_M_M]         = lbModel.fTemp8[2][lbModel.G8_DV_P_M_M];
  feq[2][lbModel.DV_P_P_M]         = lbModel.fTemp8[2][lbModel.G8_DV_P_P_M];

  feq[3][lbModel.DV_M_M_M]         = lbModel.fTemp8[3][lbModel.G8_DV_M_M_M];
  feq[3][lbModel.DV_M_P_M]         = lbModel.fTemp8[3][lbModel.G8_DV_M_P_M];
  feq[3][lbModel.DV_P_M_M]         = lbModel.fTemp8[3][lbModel.G8_DV_P_M_M];
  feq[3][lbModel.DV_P_P_M]         = lbModel.fTemp8[3][lbModel.G8_DV_P_P_M];

  feq[0][lbModel.DV_P1_P1_P1]      = lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1];
  feq[0][lbModel.DV_M1_P1_P1]      = lbModel.fTemp9[0][lbModel.G9_DV_M1_P1_P1];
  feq[0][lbModel.DV_M1_M1_P1]      = lbModel.fTemp9[0][lbModel.G9_DV_M1_M1_P1];
  feq[0][lbModel.DV_P1_M1_P1]      = lbModel.fTemp9[0][lbModel.G9_DV_P1_M1_P1];

  feq[1][lbModel.DV_P1_P1_P1]      = lbModel.fTemp9[1][lbModel.G9_DV_P1_P1_P1];
  feq[1][lbModel.DV_M1_P1_P1]      = lbModel.fTemp9[1][lbModel.G9_DV_M1_P1_P1];
  feq[1][lbModel.DV_M1_M1_P1]      = lbModel.fTemp9[1][lbModel.G9_DV_M1_M1_P1];
  feq[1][lbModel.DV_P1_M1_P1]      = lbModel.fTemp9[1][lbModel.G9_DV_P1_M1_P1];

  feq[2][lbModel.DV_P1_P1_P1]      = lbModel.fTemp9[2][lbModel.G9_DV_P1_P1_P1];
  feq[2][lbModel.DV_M1_P1_P1]      = lbModel.fTemp9[2][lbModel.G9_DV_M1_P1_P1];
  feq[2][lbModel.DV_M1_M1_P1]      = lbModel.fTemp9[2][lbModel.G9_DV_M1_M1_P1];
  feq[2][lbModel.DV_P1_M1_P1]      = lbModel.fTemp9[2][lbModel.G9_DV_P1_M1_P1];

  feq[3][lbModel.DV_P1_P1_P1]      = lbModel.fTemp9[3][lbModel.G9_DV_P1_P1_P1];
  feq[3][lbModel.DV_M1_P1_P1]      = lbModel.fTemp9[3][lbModel.G9_DV_M1_P1_P1];
  feq[3][lbModel.DV_M1_M1_P1]      = lbModel.fTemp9[3][lbModel.G9_DV_M1_M1_P1];
  feq[3][lbModel.DV_P1_M1_P1]      = lbModel.fTemp9[3][lbModel.G9_DV_P1_M1_P1];

  feq[0][lbModel.DV_M1_M1_M1]      = lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
  feq[0][lbModel.DV_M1_P1_M1]      = lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1];
  feq[0][lbModel.DV_P1_M1_M1]      = lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1];
  feq[0][lbModel.DV_P1_P1_M1]      = lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1];

  feq[1][lbModel.DV_M1_M1_M1]      = lbModel.fTemp10[1][lbModel.G10_DV_M1_M1_M1];
  feq[1][lbModel.DV_M1_P1_M1]      = lbModel.fTemp10[1][lbModel.G10_DV_M1_P1_M1];
  feq[1][lbModel.DV_P1_M1_M1]      = lbModel.fTemp10[1][lbModel.G10_DV_P1_M1_M1];
  feq[1][lbModel.DV_P1_P1_M1]      = lbModel.fTemp10[1][lbModel.G10_DV_P1_P1_M1];

  feq[2][lbModel.DV_M1_M1_M1]      = lbModel.fTemp10[2][lbModel.G10_DV_M1_M1_M1];
  feq[2][lbModel.DV_M1_P1_M1]      = lbModel.fTemp10[2][lbModel.G10_DV_M1_P1_M1];
  feq[2][lbModel.DV_P1_M1_M1]      = lbModel.fTemp10[2][lbModel.G10_DV_P1_M1_M1];
  feq[2][lbModel.DV_P1_P1_M1]      = lbModel.fTemp10[2][lbModel.G10_DV_P1_P1_M1];

  feq[3][lbModel.DV_M1_M1_M1]      = lbModel.fTemp10[3][lbModel.G10_DV_M1_M1_M1];
  feq[3][lbModel.DV_M1_P1_M1]      = lbModel.fTemp10[3][lbModel.G10_DV_M1_P1_M1];
  feq[3][lbModel.DV_P1_M1_M1]      = lbModel.fTemp10[3][lbModel.G10_DV_P1_M1_M1];
  feq[3][lbModel.DV_P1_P1_M1]      = lbModel.fTemp10[3][lbModel.G10_DV_P1_P1_M1];
}

template <int N,int numblock, typename dataType1>
void copyFromNodeTo4Pt41Array(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3, dataType1 feq[][4])
{
  unsigned long long int index;

  for(int i=0;i<VECT_LENGTH;i++)
  {
    feq[lbModel.DV_ZERO_ZERO_ZERO][i] = myGrid(i1+i,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],  0);
    for(int dv=0;dv<4;dv++)
      feq[1+dv][i]  =  myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],  0);
    for(int dv=0;dv<4;dv++)
      feq[5+dv][i]  =  myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],  0);
    for(int dv=0;dv<4;dv++)
      feq[9+dv][i]  =  myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],  0);
    for(int dv=0;dv<4;dv++)
      feq[13+dv][i]  = myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],  0);
    for(int dv=0;dv<4;dv++)
      feq[17+dv][i]  = myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],  0);
    for(int dv=0;dv<4;dv++)
      feq[21+dv][i]  = myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],  0);
    for(int dv=0;dv<4;dv++)
      feq[25+dv][i]  = myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],  0);
    for(int dv=0;dv<4;dv++)
      feq[29+dv][i]  = myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],  0);
    for(int dv=0;dv<4;dv++)
      feq[33+dv][i]  = myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],  0);
    for(int dv=0;dv<4;dv++)
      feq[37+dv][i]  = myGrid.value(index+dv) ;
  }
}

template <int N,int numblock, typename dataType1>
void copyFromCellTo4Pt41Array(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3, dataType1 feq[][4])
{
  unsigned long long int index;

  for(int i=0;i<VECT_LENGTH;i++)
  {
    feq[lbModel.DV_ZERO_ZERO_ZERO][i] = myGrid(i1+i,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],  0);
    for(int dv=0;dv<4;dv++)
      feq[1+dv][i]  =  myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],  0);
    for(int dv=0;dv<4;dv++)
      feq[5+dv][i]  =  myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],  0);
    for(int dv=0;dv<4;dv++)
      feq[9+dv][i] =  myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],  0);
    for(int dv=0;dv<4;dv++)
      feq[13+dv][i]  = myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],  0);
    for(int dv=0;dv<4;dv++)
      feq[17+dv][i]  = myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],  0);
    for(int dv=0;dv<4;dv++)
      feq[21+dv][i]  = myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],  0);
    for(int dv=0;dv<4;dv++)
      feq[25+dv][i]  = myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],  0);
    for(int dv=0;dv<4;dv++)
      feq[29+dv][i]  = myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],  0);
    for(int dv=0;dv<4;dv++)
      feq[33+dv][i]  = myGrid.value(index+dv) ;

    index = myGrid.getIndex(i1+i,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],  0);
    for(int dv=0;dv<4;dv++)
      feq[37+dv][i] = myGrid.value(index+dv) ;
  }
}

template<typename dataType1>
void calculateAlpha(lbmRD3Q41<dataType1> &lbModel,dataType1* x_i, double f_i[][41],dataType1 &alpha,int index)
{
 double a(0.), b(0.), c(0.);
 alpha = 2.0;
 for(int dv = 0;dv<lbModel.dvN;dv++)
 {
     if(x_i[dv]<0.0)
     {
         a += f_i[index][dv]*x_i[dv]*x_i[dv]*x_i[dv]*0.5;
         c += f_i[index][dv]*x_i[dv]*x_i[dv];
     }

     else
     {
         c += f_i[index][dv]*x_i[dv]*x_i[dv]/(1.0+x_i[dv]);
     }

     b += f_i[index][dv]*x_i[dv]*x_i[dv]*0.5;
 }
 alpha = (b-sqrt(b*b - 4.0*a*c))/(2.0*a);
}



template<typename dataType1>
void calculateAlphaSIMD(lbmRD3Q41<dataType1> &lbModel,
                        dataType1 x_i[][4],
                        dataType1 f_i[][4],
                        dataType1 beta,
                        dataType1* alpha)
                        {
  alignas(32) dataType1 ximin[4];
  alignas(32) dataType1 ximax[4];
  alignas(32) dataType1 k[4];
  alignas(32) dataType1 xSqr[41][4];
  alignas(32) dataType1 xCube[41][4];
  alignas(32) dataType1 f_xSqr[41][4];
  alignas(32) dataType1 f_xCube[41][4];

  SIMD_REG REG1_ = SET_ZERO();
  SIMD_REG REG2_ = SET_ZERO();
  double *pointer;
  SIMD_REG REG3_, REG4_, REG5_, REG6_, REG7_, REG8_; // 8

  for (int dv = 0; dv < 41; dv++) {
    pointer= &x_i[dv][0];
    REG3_ = LOAD_PD(pointer); // x_i
    REG2_ = MAX_PD(REG2_, REG3_); // max
    REG1_ = MIN_PD(REG1_, REG3_); // min
  }

  pointer = &ximin[0];
  STORE_PD(pointer,REG1_);
  pointer = &ximax[0];
  STORE_PD(pointer,REG2_);

  SIMD_REG A_ = SET_ZERO(); // 10
  SIMD_REG B_ = SET_ZERO(); // 11
  SIMD_REG C_ = SET_ZERO(); // 12

  SIMD_REG MASK_REG_; ///// _13th Register_

  for (int dv = 0; dv < 41; dv++) {
    pointer = &x_i[dv][0];
    REG1_ = LOAD_PD(pointer); // x_i[dv][]
    REG2_ = MUL_PD(REG1_,REG1_);  // xSqr[dv][4]
    pointer = &xSqr[dv][0];
    STORE_PD(pointer,REG2_);
    REG3_ = MUL_PD(REG2_,REG1_);  // xCube[dv][4]
    pointer = &xCube[dv][0];
    STORE_PD(pointer,REG3_);

    pointer = &f_i[dv][0];
    REG4_ = LOAD_PD(pointer);     // f_i[dv][index]
    REG4_ = MUL_PD(REG4_,REG2_);  // f_xSqr[dv][index]
    pointer = &f_xSqr[dv][0];
    STORE_PD(pointer,REG4_);
    REG5_ = MUL_PD(REG4_,REG1_);  // f_xCube[dv][index]
    pointer = &f_xCube[dv][0];
    STORE_PD(pointer,REG5_);

    REG6_ = SET_ZERO();

    MASK_REG_ = _mm256_cmp_pd(REG1_, REG6_, _CMP_LT_OQ);
    REG7_ = SET1_PD(0.5); //  0.5
    REG6_ = MUL_PD(REG5_,REG7_);  //  0.5 * f_xCube[dv][index]
    REG6_ = _mm256_and_pd(MASK_REG_, REG6_);  //  if x_i[dv][] < 0.0

    A_ = ADD_PD(A_,REG6_);  // a + 0.5 * f_xCube[dv][index] (if)

    REG8_ = MUL_PD(REG4_, REG7_); // f_xSqr[dv][index]* 0.5
    B_ = ADD_PD(REG8_, B_); // b = f_xSqr[dv][index]*0.5 + b

    REG7_ = MUL_PD(REG7_,REG1_); // 0.5*x_i[dv][]
    REG8_ = SET1_PD(1.0);
    REG8_ = ADD_PD(REG8_,REG7_); // 1.0 + 0.5*x_i[dv][]
    REG7_ = DIV_PD(REG4_,REG8_);  //  f_xSqr[index] [dv] / (0.5*x_i[dv][] + 1.0)
    C_ = ADD_PD(REG7_,C_); // c = c + f_xSqr[index] [dv] / (0.5*x_i[dv][] + 1.0)
  }

  REG5_ = SET_ZERO();
  REG6_ = _mm256_cmp_pd(A_, REG5_, _CMP_LT_OQ); // a < 0
  REG7_ = _mm256_cmp_pd(B_, REG5_, _CMP_GT_OQ); // b > 0
  REG8_ = _mm256_cmp_pd(C_, REG5_, _CMP_GT_OQ); // c > 0
  MASK_REG_ = AND_PD(REG6_,REG7_);
  MASK_REG_ = AND_PD(MASK_REG_,REG8_); // a<0 && b>0 && c>0

  REG6_ = MUL_PD(B_,B_); //  b*b
  REG7_ = MUL_PD(A_,C_); //  ac
  REG8_ = SET1_PD(4.0);  //  4
  REG8_ = MUL_PD(REG7_,REG8_); //  4ac
  REG7_ = SUB_PD(REG6_,REG8_); //  b*b-4ac
  REG7_ = SQRT_PD(REG7_); // sqrt(b*b-4ac)
  REG8_ = SET1_PD(2.0);  // 2
  REG8_ = MUL_PD(REG8_,A_);  //  2a
  REG6_ = SUB_PD(B_,REG7_);  // b-sqrt(b*b-4ac)
  REG7_ = DIV_PD(REG6_,REG8_); // (b-sqrt(b*b-4ac))/2a

  REG8_ = SET1_PD(1.5);  // 1.5

  REG1_ = AND_PD(MASK_REG_,REG7_); // if (true) (b-sqrt(b*b-4ac))/2a
  REG2_ = ANDNOT_PD(MASK_REG_, REG8_);  //  else 1.5
  SIMD_REG K_ = OR_PD(REG1_, REG2_); ///// _14th Register_
  pointer = &k[0];
  STORE_PD(pointer,K_);
  A_ = SET_ZERO();  //  a = 0
  B_ = SET_ZERO();  //  b = 0
  C_ = SET_ZERO();  //  c = 0

  SIMD_REG REG9_; //  ///// _15th Register_
  SIMD_REG REG10_; //  ///// _16th Register_

  for (int dv = 0; dv < 41; dv++) {
    pointer = &k[0];
    SIMD_REG K_ = LOAD_PD(pointer);
    REG1_ = SET1_PD(4.0);
    REG1_ = DIV_PD(REG1_,K_); // 4/k

    pointer = &x_i[dv][0];
    REG3_ = LOAD_PD(pointer);  // x_i[dv][0]

    pointer = &f_xSqr[dv][0];
    REG7_ =  LOAD_PD(pointer); // f_xSqr

    pointer = &f_xCube[dv][0];
    REG9_ = LOAD_PD(pointer); // f_xCube

    REG10_ = SET1_PD(beta*beta);

    REG2_ = SET_ZERO();
    MASK_REG_ = _mm256_cmp_pd(REG3_, REG2_, _CMP_LT_OQ); // x_i[dv][0] < 0

    REG5_ = SET1_PD(1.0/6.0);
    REG5_ = MUL_PD(REG5_,REG10_); // beta2/6
    REG4_ = MUL_PD(REG9_, REG5_); // f_xCube*beta2/6
    REG5_ = SET1_PD(0.5);
    REG8_ = MUL_PD(REG7_,REG5_); // f_xSqr*0.5
    REG4_ = AND_PD(MASK_REG_, REG4_); // if (true) f_xCube*beta2/6
    REG5_ = AND_PD(MASK_REG_, REG8_); // if (true) f_xSqr*0.5

    A_ = ADD_PD(A_,REG4_); // a = a + if (true) f_xCube*beta2/6
    B_ = ADD_PD(B_,REG5_); // b = b + if (true) f_xSqr*0.5

    REG2_ = ADD_PD(REG1_,REG3_); // 4/k + x_i[index][dv]    // temp1
    REG6_ = ADD_PD(REG2_,REG3_); // 4/k + 2*x_i[index][dv]  // temp2
    REG5_ = ADD_PD(REG6_,REG3_); // 4/k + 3*x_i[index][dv]  // temp3

    REG4_ = SET1_PD(4.0);
    REG2_ = DIV_PD(REG4_,REG2_); // 4/temp1

    REG5_ = DIV_PD(REG4_,REG5_); // 4/temp3
    REG4_ = SET1_PD(2.0);
    REG6_ = DIV_PD(REG4_,REG6_);  //  2/temp2
    REG4_ = ADD_PD(REG2_,REG6_);  //  4/temp1 + 2/temp2
    REG4_ = ADD_PD(REG4_,REG5_);  //  4/temp1 + 2/temp2 + 4/temp3

    REG6_ = MUL_PD(REG10_,REG4_); // beta2*(4/temp1 + 2/temp2 + 4/temp3)
    REG4_ = SET1_PD(15.0);
    REG10_ = DIV_PD(REG9_,REG4_);  // f_xCube/15
    REG10_ = MUL_PD(REG10_,REG6_);//(f_xCube/15)*beta2*(4/temp1+2/temp2+4/temp3)
    REG8_ = SUB_PD(REG8_,REG10_);//f_xSqr*0.5 -["Above term "]

    REG2_ = SET_ZERO();
    MASK_REG_ = _mm256_cmp_pd(REG3_, REG2_, _CMP_GT_OQ); // x_i[dv][0] > 0
    REG8_ = AND_PD(MASK_REG_,REG8_); // if x_i[dv][0] > 0

    B_ = ADD_PD(B_,REG8_);
    REG1_ = SET1_PD(60.0);
    REG2_ = MUL_PD(REG1_,REG7_); // 60*f_xSqr
    REG4_ = MUL_PD(REG1_,REG9_); // 60*f_xCube
    REG5_ = MUL_PD(REG9_,REG3_); // f_xCube * x_i[dv][0]
    REG7_ = ADD_PD(REG2_,REG4_); // 60*f_xSqr + 60*f_xCube
    REG2_ = SET1_PD(11.0);
    REG5_ = MUL_PD(REG5_,REG2_); // 11 * f_xCube * x_i[dv][0]
    REG7_ = ADD_PD(REG7_,REG5_); // 60*f_xSqr+60*f_xCube+11*f_xCube*x_i[dv][0]

    pointer = &xSqr[dv][0];
    REG2_ = LOAD_PD(pointer); // xSqr
    pointer = &xCube[dv][0];
    REG4_ = LOAD_PD(pointer); // xCube
    REG6_ = SET1_PD(90.0);
    REG5_ = MUL_PD(REG6_,REG3_);  // 90*x1
    REG1_ = ADD_PD(REG1_,REG5_);  // 60 + 90*x1
    REG6_ = SET1_PD(36.0);
    REG2_ = MUL_PD(REG2_,REG6_); // 36*xSqr
    REG1_ = ADD_PD(REG1_,REG2_);  // 60 + 90*x1 + 36*xSqr
    REG6_ = SET1_PD(3.0);
    REG3_ = MUL_PD(REG4_,REG6_); // 3*xCube
    REG1_ = ADD_PD(REG1_,REG3_);  // 60 + 90*x1 + 36*xSqr + 3*xCube

    REG7_ = DIV_PD(REG7_,REG1_);
    C_ = ADD_PD(C_,REG7_); // c += (60*f_xSqr+60*f_xCube+11*f_xCube*x_i[dv][0])
    //        /(60 + 90*x1 + 36*xSqr + 3*xCube)
  }

  REG1_ = SET_ZERO();
  REG2_ = _mm256_cmp_pd(A_, REG1_, _CMP_LT_OQ); // a < 0
  REG3_ = _mm256_cmp_pd(B_, REG1_, _CMP_GT_OQ); // b > 0
  REG4_ = _mm256_cmp_pd(C_, REG1_, _CMP_GT_OQ); // c > 0

  MASK_REG_ = AND_PD(REG2_,REG3_);
  MASK_REG_ = AND_PD(MASK_REG_,REG4_);

  REG6_ = MUL_PD(B_,B_); //  b*b
  REG7_ = MUL_PD(A_,C_); //  ac
  REG8_ = SET1_PD(4.0);  //  4
  REG8_ = MUL_PD(REG7_,REG8_); //  4ac
  REG7_ = SUB_PD(REG6_,REG8_); //  b*b-4ac
  REG7_ = SQRT_PD(REG7_); // sqrt(b*b-4ac)
  REG8_ = SET1_PD(2.0);  // 2
  REG8_ = MUL_PD(REG8_,A_);  //  2a
  REG6_ = SUB_PD(B_,REG7_);  // b-sqrt(b*b-4ac)
  REG7_ = DIV_PD(REG6_,REG8_); // (b-sqrt(b*b-4ac))/2a

  REG8_ = SET1_PD(2.1);  // 2.1

  REG1_ = AND_PD(MASK_REG_,REG7_); // if (true) (b-sqrt(b*b-4ac))/2a
  REG2_ = ANDNOT_PD(MASK_REG_, REG8_);  //  else 2.1
  REG10_ = OR_PD(REG1_, REG2_); // h

  A_ = SET_ZERO();

  for (int dv = 0; dv < 41; dv++) {
    REG2_ = SET1_PD(beta);
    REG9_ = MUL_PD(REG10_, REG2_); // h*beta

    pointer = &x_i[dv][0];
    REG3_ = LOAD_PD(pointer);  // x_i[dv][0]

    REG4_ = MUL_PD(REG9_,REG3_);  // hxBeta1 = hBeta * x_i[index][dv];
    REG5_ = MUL_PD(REG4_,REG4_);  // hxBeta1 *  hxBeta1
    REG6_ = MUL_PD(REG5_,REG4_);  // hxBeta1 *  hxBeta1 *  hxBeta1

    pointer = &f_xCube[dv][0];
    REG8_ = LOAD_PD(pointer); // f_xCube
    REG2_ = SET1_PD(beta*beta);
    REG8_ = MUL_PD(REG8_,REG2_); // f_xCube*beta2
    REG2_ = SET1_PD(12.0);
    REG4_ = DIV_PD(REG4_,REG2_); // hxBeta1/12
    REG2_ = SET1_PD(20.0);
    REG5_ = DIV_PD(REG5_,REG2_); // hxBeta2/20
    REG2_ = SET1_PD(5.0);
    REG6_ = DIV_PD(REG6_,REG2_); // hxBeta3/5
    REG2_ = SET1_PD(1.0/6.0);
    REG2_ = SUB_PD(REG2_,REG4_); // 1/6 - hxBeta1/12.0
    REG2_ = ADD_PD(REG2_,REG5_); // 1/6 - hxBeta1/12.0 + hxBeta2/20
    REG2_ = SUB_PD(REG2_,REG6_); // 1/6-hxBeta1/12.0+hxBeta2/20-hxBeta3/5
    REG8_ = MUL_PD(REG8_,REG2_); // f_xCube*beta2*["above term"]

    REG2_ = SET1_PD(0.0);
    MASK_REG_ = _mm256_cmp_pd(REG3_, REG2_, _CMP_LT_OQ); // x_i[dv][0] < 0

    REG8_ = AND_PD(MASK_REG_,REG8_);
    A_ = ADD_PD(A_, REG8_);
  }

  REG1_ = SET_ZERO();
  REG2_ = _mm256_cmp_pd(A_, REG1_, _CMP_LT_OQ); // a < 0
  REG3_ = _mm256_cmp_pd(B_, REG1_, _CMP_GT_OQ); // b > 0
  REG4_ = _mm256_cmp_pd(C_, REG1_, _CMP_GT_OQ); // c > 0

  MASK_REG_ = AND_PD(REG2_,REG3_);
  MASK_REG_ = AND_PD(MASK_REG_,REG4_);

  REG6_ = MUL_PD(B_,B_); //  b*b
  REG7_ = MUL_PD(A_,C_); //  ac
  REG8_ = SET1_PD(4.0);  //  4
  REG8_ = MUL_PD(REG7_,REG8_); //  4ac
  REG7_ = SUB_PD(REG6_,REG8_); //  b*b-4ac
  REG7_ = SQRT_PD(REG7_); // sqrt(b*b-4ac)
  REG8_ = SET1_PD(2.0);  // 2
  REG8_ = MUL_PD(REG8_,C_);  //  2c
  REG6_ = ADD_PD(B_,REG7_);  // b+sqrt(b*b-4ac)
  REG7_ = DIV_PD(REG8_,REG6_); // 2c/(b+sqrt(b*b-4ac))

  REG1_ = AND_PD(MASK_REG_,REG7_); // alpha: - if (true) (b-sqrt(b*b-4ac))/2a

  pointer = &ximin[0];
  REG2_ = LOAD_PD(pointer); // xi_min
  REG3_ = SET1_PD(beta); // beta
  REG3_ = MUL_PD(REG2_,REG3_); // beta * xi_min
  REG2_ = SET1_PD(-1.0);
  REG2_ = DIV_PD(REG2_,REG3_); // alpha_max:  -1/(beta*xi_min)


  //   REG4_ = SET1_PD(0.95);
  //   REG7_ = MUL_PD(REG4_,REG2_); // 0.95*alpha_max // if NOT(alphaMax[index]>1.0)

  REG7_ = SET1_PD(2.0);

  REG4_ = SET1_PD(1.0);
  REG6_ = ADD_PD(REG4_,REG2_); // 1 + alpha_max
  REG3_ = SET1_PD(0.5);
  REG6_ = MUL_PD(REG3_,REG6_); // 0.5 *(1+alpha_max) // if(alphaMax[index]>1.0)

  REG5_ = _mm256_cmp_pd(REG2_, REG4_, _CMP_GT_OQ); // alpha_max > 1
  REG3_ = AND_PD(REG5_,REG6_); // true elements
  REG4_ = ANDNOT_PD(REG5_,REG7_); // false
  REG3_ = OR_PD(REG3_,REG4_);

  MASK_REG_ = _mm256_cmp_pd(REG1_, REG2_, _CMP_GT_OQ); // alpha > alpha_max
  REG7_ = AND_PD(MASK_REG_,REG3_);
  REG8_ = ANDNOT_PD(MASK_REG_,REG1_);
  REG9_ = OR_PD(REG7_,REG8_);

  //   MASK_REG_ = _mm256_cmp_pd(REG1_, REG2_, _CMP_GT_OQ); // alpha > alpha_max
  //
  //   REG4_ = SET1_PD(1.0);
  //   REG5_ = _mm256_cmp_pd(REG2_, REG4_, _CMP_GT_OQ); // alpha_max > 1
  //
  //   REG6_ = ADD_PD(REG4_,REG2_); // 1 + alpha_max
  //   REG4_ = SET1_PD(0.5);
  //   REG6_ = MUL_PD(REG4_,REG6_); // 0.5 *(1+alpha_max)
  //
  //   REG4_ = SET1_PD(0.95);
  //   REG7_ = MUL_PD(REG4_,REG2_); // 0.95*alpha_max
  //
  //   REG8_ = AND_PD(MASK_REG_,REG5_); // (alpha > alpha_max) && (alpha_max > 1)
  //
  //   REG1_ = AND_PD(REG8_,REG6_); // if (above true) alpha = (1 + alpha_max)
  //
  //   REG8_ = ANDNOT_PD(REG5_,MASK_REG_); // if(alpha>alpha_max) and
  //   // NOT(alpha_max > 1)
  //   REG1_ = AND_PD(REG8_,REG7_); // if (above true) alpha = 0.95*alpha_max

  pointer = &alpha[0];
  STORE_PD(pointer, REG9_);
}

template<int N, int numblock, typename dataType1>
void get_xi_SIMD_Node(lbmRD3Q41<dataType1> &lbModel,
                      gridBCC3D<N, numblock, dataType1> &myGrid,
                      int VECT_LENGTH,
                      int i1, int i2, int i3,
                      dataType1 x_SIMD_i[][4]
                      ){
  for (int index = 0; index < VECT_LENGTH; index++) {
    x_SIMD_i[lbModel.DV_ZERO_ZERO_ZERO][index] =
      lbModel.fTemp0 [index][lbModel.CENTER_DV_ZERO_ZERO_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G0, myGrid.node[lbModel.G0],
             lbModel.CENTER_DV_ZERO_ZERO_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_ZERO_P1][index] =
      lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_ZERO_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G1, myGrid.node[lbModel.G1],
             lbModel.G1_DV_ZERO_ZERO_P1) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_ZERO_P2][index] =
      lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_ZERO_P2]/
      myGrid(i1 + index, i2, i3, lbModel.G1, myGrid.node[lbModel.G1],
             lbModel.G1_DV_ZERO_ZERO_P2)  - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_P2_ZERO  ][index] =
      lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_P2_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G1, myGrid.node[lbModel.G1],
             lbModel.G1_DV_ZERO_P2_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_P2_ZERO_ZERO  ][index] =
      lbModel.fTemp1 [index][lbModel.G1_DV_P2_ZERO_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G1, myGrid.node[lbModel.G1],
             lbModel.G1_DV_P2_ZERO_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_ZERO_M1][index] =
      lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_ZERO_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G2, myGrid.node[lbModel.G2],
             lbModel.G2_DV_ZERO_ZERO_M1) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_ZERO_M2][index] =
      lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_ZERO_M2]/
      myGrid(i1 + index, i2, i3, lbModel.G2, myGrid.node[lbModel.G2],
             lbModel.G2_DV_ZERO_ZERO_M2) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_M2_ZERO][index] =
      lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_M2_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G2, myGrid.node[lbModel.G2],
             lbModel.G2_DV_ZERO_M2_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_M2_ZERO_ZERO][index] =
      lbModel.fTemp2 [index][lbModel.G2_DV_M2_ZERO_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G2, myGrid.node[lbModel.G2],
             lbModel.G2_DV_M2_ZERO_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_P1_P1][index] =
      lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_P1_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G3, myGrid.node[lbModel.G3],
             lbModel.G3_DV_ZERO_P1_P1) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_M1_P1][index] =
      lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_M1_P1]/
       myGrid(i1 + index, i2, i3, lbModel.G3, myGrid.node[lbModel.G3],
             lbModel.G3_DV_ZERO_M1_P1) - 1.0;

    x_SIMD_i[lbModel.DV_P1_ZERO_P1][index] =
      lbModel.fTemp3 [index][lbModel.G3_DV_P1_ZERO_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G3, myGrid.node[lbModel.G3],
             lbModel.G3_DV_P1_ZERO_P1) - 1.0;

    x_SIMD_i[lbModel.DV_M1_ZERO_P1][index] =
      lbModel.fTemp3 [index][lbModel.G3_DV_M1_ZERO_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G3, myGrid.node[lbModel.G3],
             lbModel.G3_DV_M1_ZERO_P1) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_M1_M1][index] =
      lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_M1_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G4, myGrid.node[lbModel.G4],
             lbModel.G4_DV_ZERO_M1_M1) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_P1_M1][index] =
      lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_P1_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G4, myGrid.node[lbModel.G4],
             lbModel.G4_DV_ZERO_P1_M1) - 1.0;

    x_SIMD_i[lbModel.DV_M1_ZERO_M1][index] =
      lbModel.fTemp4 [index][lbModel.G4_DV_M1_ZERO_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G4, myGrid.node[lbModel.G4],
             lbModel.G4_DV_M1_ZERO_M1) - 1.0;

    x_SIMD_i[lbModel.DV_P1_ZERO_M1][index] =
      lbModel.fTemp4 [index][lbModel.G4_DV_P1_ZERO_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G4, myGrid.node[lbModel.G4],
             lbModel.G4_DV_P1_ZERO_M1) - 1.0;

    x_SIMD_i[lbModel.DV_P1_P1_ZERO][index] =
      lbModel.fTemp5 [index][lbModel.G5_DV_P1_P1_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G5, myGrid.node[lbModel.G5],
             lbModel.G5_DV_P1_P1_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_M1_P1_ZERO][index] =
      lbModel.fTemp5 [index][lbModel.G5_DV_M1_P1_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G5, myGrid.node[lbModel.G5],
             lbModel.G5_DV_M1_P1_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_P1_ZERO][index] =
      lbModel.fTemp5 [index][lbModel.G5_DV_ZERO_P1_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G5, myGrid.node[lbModel.G5],
             lbModel.G5_DV_ZERO_P1_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_P1_ZERO_ZERO][index] =
      lbModel.fTemp5 [index][lbModel.G5_DV_P1_ZERO_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G5, myGrid.node[lbModel.G5],
             lbModel.G5_DV_P1_ZERO_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_M1_M1_ZERO][index] =
      lbModel.fTemp6 [index][lbModel.G6_DV_M1_M1_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G6, myGrid.node[lbModel.G6],
             lbModel.G6_DV_M1_M1_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_P1_M1_ZERO][index] =
      lbModel.fTemp6 [index][lbModel.G6_DV_P1_M1_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G6, myGrid.node[lbModel.G6],
             lbModel.G6_DV_P1_M1_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_M1_ZERO][index] =
      lbModel.fTemp6 [index][lbModel.G6_DV_ZERO_M1_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G6, myGrid.node[lbModel.G6],
             lbModel.G6_DV_ZERO_M1_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_M1_ZERO_ZERO][index] =
      lbModel.fTemp6 [index][lbModel.G6_DV_M1_ZERO_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G6, myGrid.node[lbModel.G6],
             lbModel.G6_DV_M1_ZERO_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_P_P_P][index] =
      lbModel.fTemp7 [index][lbModel.G7_DV_P_P_P]/
      myGrid(i1 + index, i2, i3, lbModel.G7, myGrid.node[lbModel.G7],
             lbModel.G7_DV_P_P_P) - 1.0;

    x_SIMD_i[lbModel.DV_M_P_P][index] =
      lbModel.fTemp7 [index][lbModel.G7_DV_M_P_P]/
      myGrid(i1 + index, i2, i3, lbModel.G7, myGrid.node[lbModel.G7],
             lbModel.G7_DV_M_P_P) - 1.0;

    x_SIMD_i[lbModel.DV_M_M_P][index] =
      lbModel.fTemp7 [index][lbModel.G7_DV_M_M_P]/
      myGrid(i1 + index, i2, i3, lbModel.G7, myGrid.node[lbModel.G7],
             lbModel.G7_DV_M_M_P) - 1.0;

    x_SIMD_i[lbModel.DV_P_M_P][index] =
      lbModel.fTemp7 [index][lbModel.G7_DV_P_M_P]/
      myGrid(i1 + index, i2, i3, lbModel.G7, myGrid.node[lbModel.G7],
             lbModel.G7_DV_P_M_P) - 1.0;

    x_SIMD_i[lbModel.DV_M_M_M][index] =
      lbModel.fTemp8 [index][lbModel.G8_DV_M_M_M]/
      myGrid(i1 + index, i2, i3, lbModel.G8, myGrid.node[lbModel.G8],
             lbModel.G8_DV_M_M_M) - 1.0;

    x_SIMD_i[lbModel.DV_P_M_M][index] =
      lbModel.fTemp8 [index][lbModel.G8_DV_P_M_M]/
      myGrid(i1 + index, i2, i3, lbModel.G8, myGrid.node[lbModel.G8],
             lbModel.G8_DV_P_M_M) - 1.0;

    x_SIMD_i[lbModel.DV_P_P_M][index] =
      lbModel.fTemp8 [index][lbModel.G8_DV_P_P_M]/
      myGrid(i1 + index, i2, i3, lbModel.G8, myGrid.node[lbModel.G8],
             lbModel.G8_DV_P_P_M) - 1.0;

    x_SIMD_i[lbModel.DV_M_P_M][index] =
      lbModel.fTemp8 [index][lbModel.G8_DV_M_P_M]/
      myGrid(i1 + index, i2, i3, lbModel.G8, myGrid.node[lbModel.G8],
             lbModel.G8_DV_M_P_M) - 1.0;

    x_SIMD_i[lbModel.DV_P1_P1_P1][index] =
      lbModel.fTemp9 [index][lbModel.G9_DV_P1_P1_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G9, myGrid.node[lbModel.G9],
             lbModel.G9_DV_P1_P1_P1) - 1.0;

    x_SIMD_i[lbModel.DV_M1_P1_P1][index] =
      lbModel.fTemp9 [index][lbModel.G9_DV_M1_P1_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G9, myGrid.node[lbModel.G9],
             lbModel.G9_DV_M1_P1_P1) - 1.0;

    x_SIMD_i[lbModel.DV_M1_M1_P1][index] =
      lbModel.fTemp9 [index][lbModel.G9_DV_M1_M1_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G9, myGrid.node[lbModel.G9],
             lbModel.G9_DV_M1_M1_P1) - 1.0;

    x_SIMD_i[lbModel.DV_P1_M1_P1][index] =
      lbModel.fTemp9 [index][lbModel.G9_DV_P1_M1_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G9, myGrid.node[lbModel.G9],
             lbModel.G9_DV_P1_M1_P1) - 1.0;

    x_SIMD_i[lbModel.DV_M1_M1_M1][index] =
      lbModel.fTemp10[index][lbModel.G10_DV_M1_M1_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G10, myGrid.node[lbModel.G10],
             lbModel.G10_DV_M1_M1_M1) - 1.0;

    x_SIMD_i[lbModel.DV_P1_M1_M1][index] =
      lbModel.fTemp10[index][lbModel.G10_DV_P1_M1_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G10, myGrid.node[lbModel.G10],
             lbModel.G10_DV_P1_M1_M1) - 1.0;

    x_SIMD_i[lbModel.DV_P1_P1_M1][index] =
      lbModel.fTemp10[index][lbModel.G10_DV_P1_P1_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G10, myGrid.node[lbModel.G10],
             lbModel.G10_DV_P1_P1_M1) - 1.0;

    x_SIMD_i[lbModel.DV_M1_P1_M1][index] =
      lbModel.fTemp10[index][lbModel.G10_DV_M1_P1_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G10, myGrid.node[lbModel.G10],
             lbModel.G10_DV_M1_P1_M1) - 1.0;

  }

}

template<int N, int numblock, typename dataType1>
void get_xi_SIMD_Cell(lbmRD3Q41<dataType1> &lbModel,
                      gridBCC3D<N, numblock, dataType1> &myGrid,
                      int VECT_LENGTH,
                      int i1, int i2, int i3,
                      dataType1 x_SIMD_i[][4]){
  for (int index = 0; index < VECT_LENGTH; index++) {
    x_SIMD_i[lbModel.DV_ZERO_ZERO_ZERO][index] =
      lbModel.fTemp0 [index][lbModel.CENTER_DV_ZERO_ZERO_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G0, myGrid.cell[lbModel.G0],
             lbModel.CENTER_DV_ZERO_ZERO_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_ZERO_P1][index] =
      lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_ZERO_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G1, myGrid.cell[lbModel.G1],
             lbModel.G1_DV_ZERO_ZERO_P1) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_ZERO_P2][index] =
      lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_ZERO_P2]/
      myGrid(i1 + index, i2, i3, lbModel.G1, myGrid.cell[lbModel.G1],
             lbModel.G1_DV_ZERO_ZERO_P2)  - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_P2_ZERO  ][index] =
      lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_P2_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G1, myGrid.cell[lbModel.G1],
             lbModel.G1_DV_ZERO_P2_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_P2_ZERO_ZERO  ][index] =
      lbModel.fTemp1 [index][lbModel.G1_DV_P2_ZERO_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G1, myGrid.cell[lbModel.G1],
             lbModel.G1_DV_P2_ZERO_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_ZERO_M1][index] =
      lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_ZERO_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G2, myGrid.cell[lbModel.G2],
             lbModel.G2_DV_ZERO_ZERO_M1) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_ZERO_M2][index] =
      lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_ZERO_M2]/
      myGrid(i1 + index, i2, i3, lbModel.G2, myGrid.cell[lbModel.G2],
             lbModel.G2_DV_ZERO_ZERO_M2) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_M2_ZERO][index] =
      lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_M2_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G2, myGrid.cell[lbModel.G2],
             lbModel.G2_DV_ZERO_M2_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_M2_ZERO_ZERO][index] =
      lbModel.fTemp2 [index][lbModel.G2_DV_M2_ZERO_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G2, myGrid.cell[lbModel.G2],
             lbModel.G2_DV_M2_ZERO_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_P1_P1][index] =
      lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_P1_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G3, myGrid.cell[lbModel.G3],
             lbModel.G3_DV_ZERO_P1_P1) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_M1_P1][index] =
      lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_M1_P1]/
       myGrid(i1 + index, i2, i3, lbModel.G3, myGrid.cell[lbModel.G3],
             lbModel.G3_DV_ZERO_M1_P1) - 1.0;

    x_SIMD_i[lbModel.DV_P1_ZERO_P1][index] =
      lbModel.fTemp3 [index][lbModel.G3_DV_P1_ZERO_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G3, myGrid.cell[lbModel.G3],
             lbModel.G3_DV_P1_ZERO_P1) - 1.0;

    x_SIMD_i[lbModel.DV_M1_ZERO_P1][index] =
      lbModel.fTemp3 [index][lbModel.G3_DV_M1_ZERO_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G3, myGrid.cell[lbModel.G3],
             lbModel.G3_DV_M1_ZERO_P1) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_M1_M1][index] =
      lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_M1_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G4, myGrid.cell[lbModel.G4],
             lbModel.G4_DV_ZERO_M1_M1) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_P1_M1][index] =
      lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_P1_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G4, myGrid.cell[lbModel.G4],
             lbModel.G4_DV_ZERO_P1_M1) - 1.0;

    x_SIMD_i[lbModel.DV_M1_ZERO_M1][index] =
      lbModel.fTemp4 [index][lbModel.G4_DV_M1_ZERO_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G4, myGrid.cell[lbModel.G4],
             lbModel.G4_DV_M1_ZERO_M1) - 1.0;

    x_SIMD_i[lbModel.DV_P1_ZERO_M1][index] =
      lbModel.fTemp4 [index][lbModel.G4_DV_P1_ZERO_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G4, myGrid.cell[lbModel.G4],
             lbModel.G4_DV_P1_ZERO_M1) - 1.0;

    x_SIMD_i[lbModel.DV_P1_P1_ZERO][index] =
      lbModel.fTemp5 [index][lbModel.G5_DV_P1_P1_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G5, myGrid.cell[lbModel.G5],
             lbModel.G5_DV_P1_P1_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_M1_P1_ZERO][index] =
      lbModel.fTemp5 [index][lbModel.G5_DV_M1_P1_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G5, myGrid.cell[lbModel.G5],
             lbModel.G5_DV_M1_P1_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_P1_ZERO][index] =
      lbModel.fTemp5 [index][lbModel.G5_DV_ZERO_P1_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G5, myGrid.cell[lbModel.G5],
             lbModel.G5_DV_ZERO_P1_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_P1_ZERO_ZERO][index] =
      lbModel.fTemp5 [index][lbModel.G5_DV_P1_ZERO_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G5, myGrid.cell[lbModel.G5],
             lbModel.G5_DV_P1_ZERO_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_M1_M1_ZERO][index] =
      lbModel.fTemp6 [index][lbModel.G6_DV_M1_M1_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G6, myGrid.cell[lbModel.G6],
             lbModel.G6_DV_M1_M1_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_P1_M1_ZERO][index] =
      lbModel.fTemp6 [index][lbModel.G6_DV_P1_M1_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G6, myGrid.cell[lbModel.G6],
             lbModel.G6_DV_P1_M1_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_ZERO_M1_ZERO][index] =
      lbModel.fTemp6 [index][lbModel.G6_DV_ZERO_M1_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G6, myGrid.cell[lbModel.G6],
             lbModel.G6_DV_ZERO_M1_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_M1_ZERO_ZERO][index] =
      lbModel.fTemp6 [index][lbModel.G6_DV_M1_ZERO_ZERO]/
      myGrid(i1 + index, i2, i3, lbModel.G6, myGrid.cell[lbModel.G6],
             lbModel.G6_DV_M1_ZERO_ZERO) - 1.0;

    x_SIMD_i[lbModel.DV_P_P_P][index] =
      lbModel.fTemp7 [index][lbModel.G7_DV_P_P_P]/
      myGrid(i1 + index, i2, i3, lbModel.G7, myGrid.cell[lbModel.G7],
             lbModel.G7_DV_P_P_P) - 1.0;

    x_SIMD_i[lbModel.DV_M_P_P][index] =
      lbModel.fTemp7 [index][lbModel.G7_DV_M_P_P]/
      myGrid(i1 + index, i2, i3, lbModel.G7, myGrid.cell[lbModel.G7],
             lbModel.G7_DV_M_P_P) - 1.0;

    x_SIMD_i[lbModel.DV_M_M_P][index] =
      lbModel.fTemp7 [index][lbModel.G7_DV_M_M_P]/
      myGrid(i1 + index, i2, i3, lbModel.G7, myGrid.cell[lbModel.G7],
             lbModel.G7_DV_M_M_P) - 1.0;

    x_SIMD_i[lbModel.DV_P_M_P][index] =
      lbModel.fTemp7 [index][lbModel.G7_DV_P_M_P]/
      myGrid(i1 + index, i2, i3, lbModel.G7, myGrid.cell[lbModel.G7],
             lbModel.G7_DV_P_M_P) - 1.0;

    x_SIMD_i[lbModel.DV_M_M_M][index] =
      lbModel.fTemp8 [index][lbModel.G8_DV_M_M_M]/
      myGrid(i1 + index, i2, i3, lbModel.G8, myGrid.cell[lbModel.G8],
             lbModel.G8_DV_M_M_M) - 1.0;

    x_SIMD_i[lbModel.DV_P_M_M][index] =
      lbModel.fTemp8 [index][lbModel.G8_DV_P_M_M]/
      myGrid(i1 + index, i2, i3, lbModel.G8, myGrid.cell[lbModel.G8],
             lbModel.G8_DV_P_M_M) - 1.0;

    x_SIMD_i[lbModel.DV_P_P_M][index] =
      lbModel.fTemp8 [index][lbModel.G8_DV_P_P_M]/
      myGrid(i1 + index, i2, i3, lbModel.G8, myGrid.cell[lbModel.G8],
             lbModel.G8_DV_P_P_M) - 1.0;

    x_SIMD_i[lbModel.DV_M_P_M][index] =
      lbModel.fTemp8 [index][lbModel.G8_DV_M_P_M]/
      myGrid(i1 + index, i2, i3, lbModel.G8, myGrid.cell[lbModel.G8],
             lbModel.G8_DV_M_P_M) - 1.0;

    x_SIMD_i[lbModel.DV_P1_P1_P1][index] =
      lbModel.fTemp9 [index][lbModel.G9_DV_P1_P1_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G9, myGrid.cell[lbModel.G9],
             lbModel.G9_DV_P1_P1_P1) - 1.0;

    x_SIMD_i[lbModel.DV_M1_P1_P1][index] =
      lbModel.fTemp9 [index][lbModel.G9_DV_M1_P1_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G9, myGrid.cell[lbModel.G9],
             lbModel.G9_DV_M1_P1_P1) - 1.0;

    x_SIMD_i[lbModel.DV_M1_M1_P1][index] =
      lbModel.fTemp9 [index][lbModel.G9_DV_M1_M1_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G9, myGrid.cell[lbModel.G9],
             lbModel.G9_DV_M1_M1_P1) - 1.0;

    x_SIMD_i[lbModel.DV_P1_M1_P1][index] =
      lbModel.fTemp9 [index][lbModel.G9_DV_P1_M1_P1]/
      myGrid(i1 + index, i2, i3, lbModel.G9, myGrid.cell[lbModel.G9],
             lbModel.G9_DV_P1_M1_P1) - 1.0;

    x_SIMD_i[lbModel.DV_M1_M1_M1][index] =
      lbModel.fTemp10[index][lbModel.G10_DV_M1_M1_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G10, myGrid.cell[lbModel.G10],
             lbModel.G10_DV_M1_M1_M1) - 1.0;

    x_SIMD_i[lbModel.DV_P1_M1_M1][index] =
      lbModel.fTemp10[index][lbModel.G10_DV_P1_M1_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G10, myGrid.cell[lbModel.G10],
             lbModel.G10_DV_P1_M1_M1) - 1.0;

    x_SIMD_i[lbModel.DV_P1_P1_M1][index] =
      lbModel.fTemp10[index][lbModel.G10_DV_P1_P1_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G10, myGrid.cell[lbModel.G10],
             lbModel.G10_DV_P1_P1_M1) - 1.0;

    x_SIMD_i[lbModel.DV_M1_P1_M1][index] =
      lbModel.fTemp10[index][lbModel.G10_DV_M1_P1_M1]/
      myGrid(i1 + index, i2, i3, lbModel.G10, myGrid.cell[lbModel.G10],
             lbModel.G10_DV_M1_P1_M1) - 1.0;

  }

}

template<int N, int numblock, typename dataType1>
void collideWithForceEntropicSIMD(lbmRD3Q41<dataType1> &lbModel,
                                  gridBCC3D<N, numblock, dataType1> &myGrid,
                                  int VECT_LENGTH,
                                  dataType1 beta,
                                  dataType1 betaTau,
                                  dataType1 F1,
                                  dataType1 F2,
                                  dataType1 F3,
                                  dataType1 dt,
                                  int step,
                                  int procRank) {
  alignas(32) dataType1 f_i[4][41];
  dataType1 x_i[41];
  dataType1 f[41];
  dataType1 f_eq[41];

  alignas(32) dataType1 rho[4];
  alignas(32) dataType1 uX[4];
  alignas(32) dataType1 uY[4];
  alignas(32) dataType1 uZ[4];
  alignas(32) dataType1 theta[4];
  alignas(32) dataType1 alpha[4];
  alignas(32) dataType1 oneMinusalphaBeta[4];
  alignas(32) dataType1 alphaBeta[4];
  alignas(32) dataType1 rhoByTheta0[4];
  alignas(32) dataType1 x_SIMD_i[41][4];
  alignas(32) dataType1 f_SIMD_i[41][4];

  SIMD_REG _temp1, _temp2, _temp3;
  SIMD_REG _temp4, _temp5, _temp6;
  SIMD_REG _temp7, _temp8;

  dataType1 *pointer1, *pointer2, *pointer3;
  dataType1 bodyForceFactor = dt;
  dataType1 oneByTheta0 = 1.0 / lbModel.theta0;

  for (int i3 = myGrid.nB3; i3 <= myGrid.nE3; i3++) {
    for (int i2 = myGrid.nB2; i2 <= myGrid.nE2; i2++) {
      for (int i1 = myGrid.nB1; i1 <= myGrid.nE1; i1 = i1 + 4) {
        getHydroMomentsFromNodeWithForce(lbModel, myGrid, VECT_LENGTH,
                                         i1, i2, i3, rho, uX, uY, uZ, theta,
                                         F1, F2, F3, dt);

        for(int index=0; index<VECT_LENGTH;index++)
          theta[index] = lbModel.theta0;

        getFEqSIMD(lbModel, VECT_LENGTH, rho, uX, uY, uZ, theta);

        copyFromNodeTo4Pt41Array(lbModel, myGrid, VECT_LENGTH,
                                 i1, i2, i3, f_SIMD_i);

        getfTempFSIMD(lbModel, rho, F1, F2, F3);

        get_xi_SIMD_Node(lbModel, myGrid, VECT_LENGTH, i1, i2, i3, x_SIMD_i);

        for (int index = 0; index < VECT_LENGTH; index++) {
         alpha[index] = 2.0;
        }

        bool entropic = false;
        for (int index = 0; index < VECT_LENGTH; index++) {
          for( int dv = 0; dv < 41; dv++){
            if(std::fabs(x_SIMD_i[dv][index]) > 0.0001){
              entropic = true;
              break;
            }
            if(entropic)
              break;
          }
        }

        if(entropic){
          calculateAlphaSIMD(lbModel,x_SIMD_i,f_SIMD_i,beta,alpha);
        }

        for (int index = 0; index < VECT_LENGTH; index++) {
          alphaBeta[index] = alpha[index] * beta;
          oneMinusalphaBeta[index] = 1.0 - alphaBeta[index];
        }

        _temp1 = SET1_PD(bodyForceFactor); // dt


        for (int index = 0; index < VECT_LENGTH; index++) {
          myGrid(i1 + index, i2, i3, lbModel.G0,
                 myGrid.node[lbModel.G0], 0)
          = myGrid(i1 + index, i2, i3, lbModel.G0, myGrid.node[lbModel.G0], 0)
                  * oneMinusalphaBeta[index]
            + alphaBeta[index] * lbModel.fTemp0[index][0];

          _temp2 = SET1_PD(alphaBeta[index]);
          _temp3 = SET1_PD(oneMinusalphaBeta[index]);
          /////////
          // index //
          /////////
          ////////
          //G1
          ////////
          pointer1 = &lbModel.fTemp1[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G1,
                             myGrid.node[lbModel.G1], 0);
          pointer3 = &lbModel.forceTemp1[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G2
          ////////
          pointer1 = &lbModel.fTemp2[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G2,
                             myGrid.node[lbModel.G2], 0);
          pointer3 = &lbModel.forceTemp2[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G3
          ////////
          pointer1 = &lbModel.fTemp3[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G3,
                             myGrid.node[lbModel.G3], 0);
          pointer3 = &lbModel.forceTemp3[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G4
          ////////
          pointer1 = &lbModel.fTemp4[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G4,
                             myGrid.node[lbModel.G4], 0);
          pointer3 = &lbModel.forceTemp4[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G5
          ////////
          pointer1 = &lbModel.fTemp5[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G5,
                             myGrid.node[lbModel.G5], 0);
          pointer3 = &lbModel.forceTemp5[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G6
          ////////
          pointer1 = &lbModel.fTemp6[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G6,
                             myGrid.node[lbModel.G6], 0);
          pointer3 = &lbModel.forceTemp6[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G7
          ////////
          pointer1 = &lbModel.fTemp7[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G7,
                             myGrid.node[lbModel.G7], 0);
          pointer3 = &lbModel.forceTemp7[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G8
          ////////
          pointer1 = &lbModel.fTemp8[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G8,
                             myGrid.node[lbModel.G8], 0);
          pointer3 = &lbModel.forceTemp8[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G9
          ////////
          pointer1 = &lbModel.fTemp9[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G9,
                             myGrid.node[lbModel.G9], 0);
          pointer3 = &lbModel.forceTemp9[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G10
          ////////
          pointer1 = &lbModel.fTemp10[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G10,
                             myGrid.node[lbModel.G10], 0);
          pointer3 = &lbModel.forceTemp10[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
        }
      }
    }
  }

  for (int i3 = myGrid.nB3; i3 <= myGrid.nE3; i3++) {
    for (int i2 = myGrid.nB2; i2 <= myGrid.nE2; i2++) {
      for (int i1 = myGrid.nB1; i1 <= myGrid.nE1; i1 = i1 + 4) {
        getHydroMomentsFromCellWithForce(lbModel, myGrid, VECT_LENGTH,
                                         i1, i2, i3, rho, uX, uY, uZ, theta,
                                         F1, F2, F3, dt);
        for(int index=0; index<VECT_LENGTH;index++)
          theta[index] = lbModel.theta0;

        getFEqSIMD(lbModel, VECT_LENGTH, rho, uX, uY, uZ, theta);

        copyFromCellTo4Pt41Array(lbModel, myGrid, VECT_LENGTH,
                                 i1, i2, i3, f_SIMD_i);

        getfTempFSIMD(lbModel, rho, F1, F2, F3);

        get_xi_SIMD_Cell(lbModel,myGrid, VECT_LENGTH, i1, i2, i3, x_SIMD_i);

        for (int index = 0; index < VECT_LENGTH; index++) {
         alpha[index] = 2.0;
        }

        bool entropic = false;
        for (int index = 0; index < VECT_LENGTH; index++) {
          for( int dv = 0; dv < 41; dv++){
            if(std::fabs(x_SIMD_i[dv][index]) > 0.0001){
              entropic = true;
              break;
            }
            if(entropic)
              break;
          }
        }

        if(entropic){
          calculateAlphaSIMD(lbModel,x_SIMD_i,f_SIMD_i,beta,alpha);
        }

        for (int index = 0; index < VECT_LENGTH; index++) {
          alphaBeta[index] = alpha[index] * beta;
          oneMinusalphaBeta[index] = 1.0 - alphaBeta[index];
        }

        _temp1 = SET1_PD(bodyForceFactor);           // to dt

        for (int index = 0; index < VECT_LENGTH; index++) {
          myGrid(i1 + index, i2, i3, lbModel.G0,
                 myGrid.cell[lbModel.G0], 0)
          = myGrid(i1 + index, i2, i3, lbModel.G0, myGrid.cell[lbModel.G0], 0)
          * oneMinusalphaBeta[index]
          + alphaBeta[index] * lbModel.fTemp0[index][0];

          _temp2 = SET1_PD(alphaBeta[index]);
          _temp3 = SET1_PD(oneMinusalphaBeta[index]);

          ////////
          //G1
          ////////
          pointer1 = &lbModel.fTemp1[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G1,
                             myGrid.cell[lbModel.G1], 0);
          pointer3 = &lbModel.forceTemp1[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G2
          ////////
          pointer1 = &lbModel.fTemp2[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G2,
                             myGrid.cell[lbModel.G2], 0);
          pointer3 = &lbModel.forceTemp2[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G3
          ////////
          pointer1 = &lbModel.fTemp3[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G3,
                             myGrid.cell[lbModel.G3], 0);
          pointer3 = &lbModel.forceTemp3[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G4
          ////////
          pointer1 = &lbModel.fTemp4[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G4,
                             myGrid.cell[lbModel.G4], 0);
          pointer3 = &lbModel.forceTemp4[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G5
          ////////
          pointer1 = &lbModel.fTemp5[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G5,
                             myGrid.cell[lbModel.G5], 0);
          pointer3 = &lbModel.forceTemp5[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G6
          ////////
          pointer1 = &lbModel.fTemp6[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G6,
                             myGrid.cell[lbModel.G6], 0);
          pointer3 = &lbModel.forceTemp6[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G7
          ////////
          pointer1 = &lbModel.fTemp7[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G7,
                             myGrid.cell[lbModel.G7], 0);
          pointer3 = &lbModel.forceTemp7[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G8
          ////////
          pointer1 = &lbModel.fTemp8[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G8,
                             myGrid.cell[lbModel.G8], 0);
          pointer3 = &lbModel.forceTemp8[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G9
          ////////
          pointer1 = &lbModel.fTemp9[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G9,
                             myGrid.cell[lbModel.G9], 0);
          pointer3 = &lbModel.forceTemp9[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
          ////////
          //G10
          ////////
          pointer1 = &lbModel.fTemp10[index][0];
          pointer2 = &myGrid(i1 + index, i2, i3, lbModel.G10,
                             myGrid.cell[lbModel.G10], 0);
          pointer3 = &lbModel.forceTemp10[index][0];
          _temp5 = LOAD_PD(pointer1);      // feq - i dv-0-3
          _temp6 = LOAD_PD(pointer2);      // SC1 - i dv-0-3
          _temp7 = LOAD_PD(pointer3);      // fTempF
          _temp5 = MUL_PD(_temp2, _temp5);  // ab*fEq
          _temp6 = MUL_PD(_temp3, _temp6);  // (1-ab)*SC1
          _temp7 = MUL_PD(_temp1, _temp7);  // dt*fTempF
          _temp8 = ADD_PD(_temp5, _temp6);  // (1-ab)*SC1 + ab*fEq
          _temp8 = ADD_PD(_temp7, _temp8);  // (1-ab)*SC1 + ab*fEq + fTempF
          STORE_PD(pointer2, _temp8);       // write it back to G1
        }
      }
    }
  }
}


// Explicit declaration
template void collide<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,int ,double ,double );
template void getMomentsWithForce<double>(lbmRD3Q41<double> &,double *, double , double , double );
template void collideWithForce<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,int , double , double  , double , double , double  , double);
template void getfTempFSIMD<double>(lbmRD3Q41<double> &, double *, double , double , double);
template void copyFromNodeTo4Pt41Array<4,11,double>(lbmRD3Q41<double> &lbModel, gridBCC3D<4,11,double> &myGrid,int VECT_LENGTH,int i1, int i2, int i3, double feq[][4]);
template void copyFromCellTo4Pt41Array<4,11,double>(lbmRD3Q41<double> &lbModel, gridBCC3D<4,11,double> &myGrid,int VECT_LENGTH,int i1, int i2, int i3, double feq[][4]);
template void collideWithForceSIMD<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,int , double , double  , double , double , double, double );
template void collideWithForceEntropicSIMD<4, 11, double>(lbmRD3Q41<double>&, gridBCC3D<4, 11, double>&, int, double, double, double, double, double, double, int, int);
template void calculateAlpha<double>(lbmRD3Q41<double> &,double* , double [][41],double &,int );
template void calculateAlphaSIMD<double>(lbmRD3Q41<double>&, double [][4], double [][4], double, double*);

template void get_xi_SIMD_Node<4, 11, double>(lbmRD3Q41<double> &,
 gridBCC3D<4, 11, double> &,
 int,
 int, int , int,
 double [][4]);

template void get_xi_SIMD_Cell<4, 11, double>(lbmRD3Q41<double> &,
 gridBCC3D<4, 11, double> &,
 int,
 int, int , int,
 double [][4]);
