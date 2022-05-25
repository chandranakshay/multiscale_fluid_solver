#include<diffuse41.h>


template <int N,int numblock, typename dataType1>
void prepareDiffuseTopWall(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)
{
 for (int k = myGrid.nB3; k<= myGrid.nE3; k++)
  for (int i = myGrid.nB1; i <= myGrid.nE1; i++)
  {
   myGrid(i,myGrid.nE2+2,k,lbModel.G1 ,1,lbModel.G1_DV_ZERO_P2_ZERO) = myGrid(i,myGrid.nE2-1,k,lbModel.G1 ,myGrid.cell[lbModel.G1]  ,lbModel.G1_DV_ZERO_P2_ZERO);

   myGrid(i,myGrid.nE2+1,k,lbModel.G1 ,1,lbModel.G1_DV_ZERO_P2_ZERO) = myGrid(i,myGrid.nE2  ,k,lbModel.G1 ,myGrid.cell[lbModel.G1]  ,lbModel.G1_DV_ZERO_P2_ZERO);
   myGrid(i,myGrid.nE2+1,k,lbModel.G3 ,1,lbModel.G3_DV_ZERO_P1_P1  ) = myGrid(i,myGrid.nE2  ,k,lbModel.G3 ,myGrid.cell[lbModel.G3]  ,lbModel.G3_DV_ZERO_P1_P1  );
   myGrid(i,myGrid.nE2+1,k,lbModel.G4 ,1,lbModel.G4_DV_ZERO_P1_M1  ) = myGrid(i,myGrid.nE2  ,k,lbModel.G4 ,myGrid.cell[lbModel.G4]  ,lbModel.G4_DV_ZERO_P1_M1  );
   myGrid(i,myGrid.nE2+1,k,lbModel.G5 ,1,lbModel.G5_DV_P1_P1_ZERO  ) = myGrid(i,myGrid.nE2  ,k,lbModel.G5 ,myGrid.cell[lbModel.G5]  ,lbModel.G5_DV_P1_P1_ZERO  );
   myGrid(i,myGrid.nE2+1,k,lbModel.G5 ,1,lbModel.G5_DV_M1_P1_ZERO  ) = myGrid(i,myGrid.nE2  ,k,lbModel.G5 ,myGrid.cell[lbModel.G5]  ,lbModel.G5_DV_M1_P1_ZERO  );
   myGrid(i,myGrid.nE2+1,k,lbModel.G5 ,1,lbModel.G5_DV_ZERO_P1_ZERO) = myGrid(i,myGrid.nE2  ,k,lbModel.G5 ,myGrid.cell[lbModel.G5]  ,lbModel.G5_DV_ZERO_P1_ZERO);
   myGrid(i,myGrid.nE2+1,k,lbModel.G9 ,1,lbModel.G9_DV_P1_P1_P1    ) = myGrid(i,myGrid.nE2  ,k,lbModel.G9 ,myGrid.cell[lbModel.G9]  ,lbModel.G9_DV_P1_P1_P1    );
   myGrid(i,myGrid.nE2+1,k,lbModel.G9 ,1,lbModel.G9_DV_M1_P1_P1    ) = myGrid(i,myGrid.nE2  ,k,lbModel.G9 ,myGrid.cell[lbModel.G9]  ,lbModel.G9_DV_M1_P1_P1    );
   myGrid(i,myGrid.nE2+1,k,lbModel.G10,1,lbModel.G10_DV_P1_P1_M1   ) = myGrid(i,myGrid.nE2  ,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_P1_P1_M1   );
   myGrid(i,myGrid.nE2+1,k,lbModel.G10,1,lbModel.G10_DV_M1_P1_M1   ) = myGrid(i,myGrid.nE2  ,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_M1_P1_M1   );

   myGrid(i,myGrid.nE2+1,k,lbModel.G7 ,1,lbModel.G7_DV_P_P_P       ) = myGrid(i,myGrid.nE2  ,k,lbModel.G7 ,myGrid.cell[lbModel.G7]  ,lbModel.G7_DV_P_P_P       );
   myGrid(i,myGrid.nE2+1,k,lbModel.G7 ,1,lbModel.G7_DV_M_P_P       ) = myGrid(i,myGrid.nE2  ,k,lbModel.G7 ,myGrid.cell[lbModel.G7]  ,lbModel.G7_DV_M_P_P       );
   myGrid(i,myGrid.nE2+1,k,lbModel.G8 ,1,lbModel.G8_DV_P_P_M       ) = myGrid(i,myGrid.nE2  ,k,lbModel.G8 ,myGrid.cell[lbModel.G8]  ,lbModel.G8_DV_P_P_M       );
   myGrid(i,myGrid.nE2+1,k,lbModel.G8 ,1,lbModel.G8_DV_M_P_M       ) = myGrid(i,myGrid.nE2  ,k,lbModel.G8 ,myGrid.cell[lbModel.G8]  ,lbModel.G8_DV_M_P_M       );
  }

 for (int k = myGrid.nB3; k<= myGrid.nE3; k++)
  for (int i = myGrid.nB1; i <= myGrid.nE1; i++)
  {
   myGrid(i,myGrid.nE2+2,k,lbModel.G1 ,0,lbModel.G1_DV_ZERO_P2_ZERO) = myGrid(i,myGrid.nE2-1,k,lbModel.G1 ,myGrid.node[lbModel.G1]  ,lbModel.G1_DV_ZERO_P2_ZERO);

   myGrid(i,myGrid.nE2+1,k,lbModel.G1 ,0,lbModel.G1_DV_ZERO_P2_ZERO) = myGrid(i,myGrid.nE2  ,k,lbModel.G1 ,myGrid.node[lbModel.G1]  ,lbModel.G1_DV_ZERO_P2_ZERO) ;
   myGrid(i,myGrid.nE2+1,k,lbModel.G3 ,0,lbModel.G3_DV_ZERO_P1_P1  ) = myGrid(i,myGrid.nE2  ,k,lbModel.G3 ,myGrid.node[lbModel.G3]  ,lbModel.G3_DV_ZERO_P1_P1  ) ;
   myGrid(i,myGrid.nE2+1,k,lbModel.G4 ,0,lbModel.G4_DV_ZERO_P1_M1  ) = myGrid(i,myGrid.nE2  ,k,lbModel.G4 ,myGrid.node[lbModel.G4]  ,lbModel.G4_DV_ZERO_P1_M1  ) ;
   myGrid(i,myGrid.nE2+1,k,lbModel.G5 ,0,lbModel.G5_DV_P1_P1_ZERO  ) = myGrid(i,myGrid.nE2  ,k,lbModel.G5 ,myGrid.node[lbModel.G5]  ,lbModel.G5_DV_P1_P1_ZERO  ) ;
   myGrid(i,myGrid.nE2+1,k,lbModel.G5 ,0,lbModel.G5_DV_M1_P1_ZERO  ) = myGrid(i,myGrid.nE2  ,k,lbModel.G5 ,myGrid.node[lbModel.G5]  ,lbModel.G5_DV_M1_P1_ZERO  ) ;
   myGrid(i,myGrid.nE2+1,k,lbModel.G5 ,0,lbModel.G5_DV_ZERO_P1_ZERO) = myGrid(i,myGrid.nE2  ,k,lbModel.G5 ,myGrid.node[lbModel.G5]  ,lbModel.G5_DV_ZERO_P1_ZERO) ;
   myGrid(i,myGrid.nE2+1,k,lbModel.G9 ,0,lbModel.G9_DV_P1_P1_P1    ) = myGrid(i,myGrid.nE2  ,k,lbModel.G9 ,myGrid.node[lbModel.G9]  ,lbModel.G9_DV_P1_P1_P1    ) ;
   myGrid(i,myGrid.nE2+1,k,lbModel.G9 ,0,lbModel.G9_DV_M1_P1_P1    ) = myGrid(i,myGrid.nE2  ,k,lbModel.G9 ,myGrid.node[lbModel.G9]  ,lbModel.G9_DV_M1_P1_P1    ) ;
   myGrid(i,myGrid.nE2+1,k,lbModel.G10,0,lbModel.G10_DV_P1_P1_M1   ) = myGrid(i,myGrid.nE2  ,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_P1_P1_M1   ) ;
   myGrid(i,myGrid.nE2+1,k,lbModel.G10,0,lbModel.G10_DV_M1_P1_M1   ) = myGrid(i,myGrid.nE2  ,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_M1_P1_M1   ) ;
  }
}

template <int N,int numblock, typename dataType1>
void correctDiffuseF0massManipTopWall(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int VECT_LENGTH, dataType1 tempTop, double uRef)
{
    dataType1 *rhoWall,*uX,*uY,*uZ,*thetaTop;
    rhoWall   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uX        = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uY        = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uZ        = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    thetaTop  = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);

 dataType1 wallFlux, wallFluxEq, probability, massLeak ;


 for(int i=0; i<VECT_LENGTH; i++)
 {
  thetaTop[i]    = tempTop;
  uX[i]          = uRef;
  uY[i]          = 0.0;
  uZ[i]          = 0.0;
  rhoWall[i]     = 1.0;
 }

                getFEqSIMD(lbModel,VECT_LENGTH,rhoWall,uX,uY,uZ,thetaTop,10);
//    getIterativeFEq(lbModel,VECT_LENGTH,rhoWall,uX,uY,uZ,thetaTop);
//               getFEqNewIterative(lbModel,VECT_LENGTH,rhoWall,uX,uY,uZ,thetaTop);


 for (int k = myGrid.nB3; k<= myGrid.nE3; k++)
  for (int i = myGrid.nB1; i <= myGrid.nE1; i=i+VECT_LENGTH)
  {
   for(int index=0; index<VECT_LENGTH; index++)
   {
    wallFlux = 0.0;
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_P2_ZERO] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G1 ,1,lbModel.G1_DV_ZERO_P2_ZERO));
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_P1_P1  ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G3 ,1,lbModel.G3_DV_ZERO_P1_P1  ));
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_P1_M1  ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G4 ,1,lbModel.G4_DV_ZERO_P1_M1  ));
    wallFlux += (lbModel.cy[lbModel.DV_P1_P1_ZERO  ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G5 ,1,lbModel.G5_DV_P1_P1_ZERO  ));
    wallFlux += (lbModel.cy[lbModel.DV_M1_P1_ZERO  ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G5 ,1,lbModel.G5_DV_M1_P1_ZERO  ));
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_P1_ZERO] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G5 ,1,lbModel.G5_DV_ZERO_P1_ZERO));
    wallFlux += (lbModel.cy[lbModel.DV_P1_P1_P1    ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G9 ,1,lbModel.G9_DV_P1_P1_P1    ));
    wallFlux += (lbModel.cy[lbModel.DV_M1_P1_P1    ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G9 ,1,lbModel.G9_DV_M1_P1_P1    ));
    wallFlux += (lbModel.cy[lbModel.DV_P1_P1_M1    ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G10,1,lbModel.G10_DV_P1_P1_M1   ));
    wallFlux += (lbModel.cy[lbModel.DV_M1_P1_M1    ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G10,1,lbModel.G10_DV_M1_P1_M1   ));
    wallFlux += (lbModel.cy[lbModel.DV_P_P_P       ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G7 ,1,lbModel.G7_DV_P_P_P       ));
    wallFlux += (lbModel.cy[lbModel.DV_M_P_P       ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G7 ,1,lbModel.G7_DV_M_P_P       ));
    wallFlux += (lbModel.cy[lbModel.DV_P_P_M       ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G8 ,1,lbModel.G8_DV_P_P_M       ));
    wallFlux += (lbModel.cy[lbModel.DV_M_P_M       ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G8 ,1,lbModel.G8_DV_M_P_M       ));

    massLeak = 0.0;
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G1 ,1,lbModel.G1_DV_ZERO_P2_ZERO);
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G3 ,1,lbModel.G3_DV_ZERO_P1_P1  );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G4 ,1,lbModel.G4_DV_ZERO_P1_M1  );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G5 ,1,lbModel.G5_DV_P1_P1_ZERO  );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G5 ,1,lbModel.G5_DV_M1_P1_ZERO  );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G5 ,1,lbModel.G5_DV_ZERO_P1_ZERO);
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G9 ,1,lbModel.G9_DV_P1_P1_P1    );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G9 ,1,lbModel.G9_DV_M1_P1_P1    );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G10,1,lbModel.G10_DV_P1_P1_M1   );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G10,1,lbModel.G10_DV_M1_P1_M1   );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G7 ,1,lbModel.G7_DV_P_P_P       );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G7 ,1,lbModel.G7_DV_M_P_P       );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G8 ,1,lbModel.G8_DV_P_P_M       );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G8 ,1,lbModel.G8_DV_M_P_M       );


    wallFluxEq  = 0.0;
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_P2_ZERO] * lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_P2_ZERO]);
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_P1_P1  ] * lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_P1_P1  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_P1_M1  ] * lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_P1_M1  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P1_P1_ZERO  ] * lbModel.fTemp5 [index][lbModel.G5_DV_P1_P1_ZERO  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M1_P1_ZERO  ] * lbModel.fTemp5 [index][lbModel.G5_DV_M1_P1_ZERO  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_P1_ZERO] * lbModel.fTemp5 [index][lbModel.G5_DV_ZERO_P1_ZERO]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P1_P1_P1    ] * lbModel.fTemp9 [index][lbModel.G9_DV_P1_P1_P1    ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M1_P1_P1    ] * lbModel.fTemp9 [index][lbModel.G9_DV_M1_P1_P1    ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P1_P1_M1    ] * lbModel.fTemp10[index][lbModel.G10_DV_P1_P1_M1  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M1_P1_M1    ] * lbModel.fTemp10[index][lbModel.G10_DV_M1_P1_M1  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P_P_P       ] * lbModel.fTemp7 [index][lbModel.G7_DV_P_P_P       ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M_P_P       ] * lbModel.fTemp7 [index][lbModel.G7_DV_M_P_P       ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P_P_M       ] * lbModel.fTemp8 [index][lbModel.G8_DV_P_P_M       ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M_P_M       ] * lbModel.fTemp8 [index][lbModel.G8_DV_M_P_M       ]);

    probability = wallFlux/wallFluxEq ;

    myGrid(i+index,myGrid.nE2,k,lbModel.G2 ,myGrid.cell[lbModel.G2 ] ,lbModel.G2_DV_ZERO_M2_ZERO)   =  probability * lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_P2_ZERO]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G4 ,myGrid.cell[lbModel.G4 ] ,lbModel.G4_DV_ZERO_M1_M1  )   =  probability * lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_P1_P1  ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G3 ,myGrid.cell[lbModel.G3 ] ,lbModel.G3_DV_ZERO_M1_P1  )   =  probability * lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_P1_M1  ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.cell[lbModel.G6 ] ,lbModel.G6_DV_M1_M1_ZERO  )   =  probability * lbModel.fTemp5 [index][lbModel.G5_DV_P1_P1_ZERO  ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.cell[lbModel.G6 ] ,lbModel.G6_DV_P1_M1_ZERO  )   =  probability * lbModel.fTemp5 [index][lbModel.G5_DV_M1_P1_ZERO  ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.cell[lbModel.G6 ] ,lbModel.G6_DV_ZERO_M1_ZERO)   =  probability * lbModel.fTemp5 [index][lbModel.G5_DV_ZERO_P1_ZERO]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_M1_M1_M1   )   =  probability * lbModel.fTemp9 [index][lbModel.G9_DV_P1_P1_P1    ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_P1_M1_M1   )   =  probability * lbModel.fTemp9 [index][lbModel.G9_DV_M1_P1_P1    ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G9 ,myGrid.cell[lbModel.G9 ] ,lbModel.G9_DV_M1_M1_P1    )   =  probability * lbModel.fTemp10[index][lbModel.G10_DV_P1_P1_M1  ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G9 ,myGrid.cell[lbModel.G9 ] ,lbModel.G9_DV_P1_M1_P1    )   =  probability * lbModel.fTemp10[index][lbModel.G10_DV_M1_P1_M1  ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G8 ,myGrid.cell[lbModel.G8 ] ,lbModel.G8_DV_M_M_M       )   =  probability * lbModel.fTemp7 [index][lbModel.G7_DV_P_P_P       ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G8 ,myGrid.cell[lbModel.G8 ] ,lbModel.G8_DV_P_M_M       )   =  probability * lbModel.fTemp7 [index][lbModel.G7_DV_M_P_P       ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G7 ,myGrid.cell[lbModel.G7 ] ,lbModel.G7_DV_M_M_P       )   =  probability * lbModel.fTemp8 [index][lbModel.G8_DV_P_P_M       ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G7 ,myGrid.cell[lbModel.G7 ] ,lbModel.G7_DV_P_M_P       )   =  probability * lbModel.fTemp8 [index][lbModel.G8_DV_M_P_M       ]  ;

    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G2 ,myGrid.cell[lbModel.G2 ] ,lbModel.G2_DV_ZERO_M2_ZERO) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G4 ,myGrid.cell[lbModel.G4 ] ,lbModel.G4_DV_ZERO_M1_M1  ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G3 ,myGrid.cell[lbModel.G3 ] ,lbModel.G3_DV_ZERO_M1_P1  ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.cell[lbModel.G6 ] ,lbModel.G6_DV_M1_M1_ZERO  ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.cell[lbModel.G6 ] ,lbModel.G6_DV_P1_M1_ZERO  ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.cell[lbModel.G6 ] ,lbModel.G6_DV_ZERO_M1_ZERO) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_M1_M1_M1   ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_P1_M1_M1   ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G9 ,myGrid.cell[lbModel.G9 ] ,lbModel.G9_DV_M1_M1_P1    ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G9 ,myGrid.cell[lbModel.G9 ] ,lbModel.G9_DV_P1_M1_P1    ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G8 ,myGrid.cell[lbModel.G8 ] ,lbModel.G8_DV_M_M_M       ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G8 ,myGrid.cell[lbModel.G8 ] ,lbModel.G8_DV_P_M_M       ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G7 ,myGrid.cell[lbModel.G7 ] ,lbModel.G7_DV_M_M_P       ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G7 ,myGrid.cell[lbModel.G7 ] ,lbModel.G7_DV_P_M_P       ) ;


    myGrid(i+index,myGrid.nE2,k, lbModel.G0, myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) += massLeak;

//     myGrid(i+index,myGrid.nE2,k,lbModel.G2 ,myGrid.cell[lbModel.G2 ] ,lbModel.G2_DV_ZERO_M2_ZERO)  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G4 ,myGrid.cell[lbModel.G4 ] ,lbModel.G4_DV_ZERO_M1_M1  )  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G3 ,myGrid.cell[lbModel.G3 ] ,lbModel.G3_DV_ZERO_M1_P1  )  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.cell[lbModel.G6 ] ,lbModel.G6_DV_M1_M1_ZERO  )  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.cell[lbModel.G6 ] ,lbModel.G6_DV_P1_M1_ZERO  )  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.cell[lbModel.G6 ] ,lbModel.G6_DV_ZERO_M1_ZERO)  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_M1_M1_M1   )  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_P1_M1_M1   )  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G9 ,myGrid.cell[lbModel.G9 ] ,lbModel.G9_DV_M1_M1_P1    )  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G9 ,myGrid.cell[lbModel.G9 ] ,lbModel.G9_DV_P1_M1_P1    )  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G8 ,myGrid.cell[lbModel.G8 ] ,lbModel.G8_DV_M_M_M       )  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G8 ,myGrid.cell[lbModel.G8 ] ,lbModel.G8_DV_P_M_M       )  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G7 ,myGrid.cell[lbModel.G7 ] ,lbModel.G7_DV_M_M_P       )  += (massLeak/14.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G7 ,myGrid.cell[lbModel.G7 ] ,lbModel.G7_DV_P_M_P       )  += (massLeak/14.0) ;
   }



   // This can be repalced with normal BB later
   for(int index=0; index<VECT_LENGTH; index++)
   {
    double wallFlux = 0.0;
    wallFlux   += (lbModel.cy[lbModel.DV_ZERO_P2_ZERO] * myGrid(i+index,myGrid.nE2+2,k,lbModel.G1 ,1  ,lbModel.G1_DV_ZERO_P2_ZERO));

    massLeak = 0.0;
    massLeak   += myGrid(i+index,myGrid.nE2+2,k,lbModel.G1 ,1 ,lbModel.G1_DV_ZERO_P2_ZERO);

    wallFluxEq = 0.0;
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_P2_ZERO] * lbModel.fTemp1[index][lbModel.G1_DV_ZERO_P2_ZERO]);

    probability = wallFlux/wallFluxEq ;

    myGrid(i+index,myGrid.nE2-1,k,lbModel.G2 ,myGrid.cell[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO)   =  probability * lbModel.fTemp1[index][lbModel.G1_DV_ZERO_P2_ZERO]  ;
    massLeak   -=  myGrid(i+index,myGrid.nE2-1,k,lbModel.G2 ,myGrid.cell[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO) ;

    myGrid(i+index,myGrid.nE2-1,k, lbModel.G0, myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) += massLeak;
//     myGrid(i+index,myGrid.nE2-1,k,lbModel.G2 ,myGrid.cell[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO) += massLeak;
   }
 }



 // Node
 for (int k = myGrid.nB3; k<= myGrid.nE3; k++)
  for (int i = myGrid.nB1; i <= myGrid.nE1; i=i+VECT_LENGTH)
  {
   for(int index=0; index<VECT_LENGTH; index++)
   {
    wallFlux = 0.0;

    wallFlux += (lbModel.cy[lbModel.DV_ZERO_P2_ZERO] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G1 ,0,lbModel.G1_DV_ZERO_P2_ZERO));
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_P1_P1  ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G3 ,0,lbModel.G3_DV_ZERO_P1_P1  ));
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_P1_M1  ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G4 ,0,lbModel.G4_DV_ZERO_P1_M1  ));
    wallFlux += (lbModel.cy[lbModel.DV_P1_P1_ZERO  ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G5 ,0,lbModel.G5_DV_P1_P1_ZERO  ));
    wallFlux += (lbModel.cy[lbModel.DV_M1_P1_ZERO  ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G5 ,0,lbModel.G5_DV_M1_P1_ZERO  ));
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_P1_ZERO] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G5 ,0,lbModel.G5_DV_ZERO_P1_ZERO));
    wallFlux += (lbModel.cy[lbModel.DV_P1_P1_P1    ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G9 ,0,lbModel.G9_DV_P1_P1_P1    ));
    wallFlux += (lbModel.cy[lbModel.DV_M1_P1_P1    ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G9 ,0,lbModel.G9_DV_M1_P1_P1    ));
    wallFlux += (lbModel.cy[lbModel.DV_P1_P1_M1    ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G10,0,lbModel.G10_DV_P1_P1_M1   ));
    wallFlux += (lbModel.cy[lbModel.DV_M1_P1_M1    ] * myGrid(i+index,myGrid.nE2+1,k,lbModel.G10,0,lbModel.G10_DV_M1_P1_M1   ));

    massLeak = 0.0;

    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G1 ,0 ,lbModel.G1_DV_ZERO_P2_ZERO);
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G3 ,0 ,lbModel.G3_DV_ZERO_P1_P1  );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G4 ,0 ,lbModel.G4_DV_ZERO_P1_M1  );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G5 ,0 ,lbModel.G5_DV_P1_P1_ZERO  );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G5 ,0 ,lbModel.G5_DV_M1_P1_ZERO  );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G5 ,0 ,lbModel.G5_DV_ZERO_P1_ZERO);
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G9 ,0 ,lbModel.G9_DV_P1_P1_P1    );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G9 ,0 ,lbModel.G9_DV_M1_P1_P1    );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G10,0 ,lbModel.G10_DV_P1_P1_M1   );
    massLeak += myGrid(i+index,myGrid.nE2+1,k,lbModel.G10,0 ,lbModel.G10_DV_M1_P1_M1   );


    wallFluxEq  = 0.0;
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_P2_ZERO] * lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_P2_ZERO]);
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_P1_P1  ] * lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_P1_P1  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_P1_M1  ] * lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_P1_M1  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P1_P1_ZERO  ] * lbModel.fTemp5 [index][lbModel.G5_DV_P1_P1_ZERO  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M1_P1_ZERO  ] * lbModel.fTemp5 [index][lbModel.G5_DV_M1_P1_ZERO  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_P1_ZERO] * lbModel.fTemp5 [index][lbModel.G5_DV_ZERO_P1_ZERO]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P1_P1_P1    ] * lbModel.fTemp9 [index][lbModel.G9_DV_P1_P1_P1    ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M1_P1_P1    ] * lbModel.fTemp9 [index][lbModel.G9_DV_M1_P1_P1    ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P1_P1_M1    ] * lbModel.fTemp10[index][lbModel.G10_DV_P1_P1_M1  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M1_P1_M1    ] * lbModel.fTemp10[index][lbModel.G10_DV_M1_P1_M1  ]);

    probability = wallFlux/wallFluxEq ;

    myGrid(i+index,myGrid.nE2,k,lbModel.G2 ,myGrid.node[lbModel.G2 ] ,lbModel.G2_DV_ZERO_M2_ZERO)   =  probability * lbModel.fTemp1 [index][lbModel.G1_DV_ZERO_P2_ZERO]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G4 ,myGrid.node[lbModel.G4 ] ,lbModel.G4_DV_ZERO_M1_M1  )   =  probability * lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_P1_P1  ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G3 ,myGrid.node[lbModel.G3 ] ,lbModel.G3_DV_ZERO_M1_P1  )   =  probability * lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_P1_M1  ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.node[lbModel.G6 ] ,lbModel.G6_DV_M1_M1_ZERO  )   =  probability * lbModel.fTemp5 [index][lbModel.G5_DV_P1_P1_ZERO  ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.node[lbModel.G6 ] ,lbModel.G6_DV_P1_M1_ZERO  )   =  probability * lbModel.fTemp5 [index][lbModel.G5_DV_M1_P1_ZERO  ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.node[lbModel.G6 ] ,lbModel.G6_DV_ZERO_M1_ZERO)   =  probability * lbModel.fTemp5 [index][lbModel.G5_DV_ZERO_P1_ZERO]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_M1_M1_M1   )   =  probability * lbModel.fTemp9 [index][lbModel.G9_DV_P1_P1_P1    ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_P1_M1_M1   )   =  probability * lbModel.fTemp9 [index][lbModel.G9_DV_M1_P1_P1    ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G9 ,myGrid.node[lbModel.G9 ] ,lbModel.G9_DV_M1_M1_P1    )   =  probability * lbModel.fTemp10[index][lbModel.G10_DV_P1_P1_M1  ]  ;
    myGrid(i+index,myGrid.nE2,k,lbModel.G9 ,myGrid.node[lbModel.G9 ] ,lbModel.G9_DV_P1_M1_P1    )   =  probability * lbModel.fTemp10[index][lbModel.G10_DV_M1_P1_M1  ]  ;

    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G2 ,myGrid.node[lbModel.G2 ] ,lbModel.G2_DV_ZERO_M2_ZERO) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G4 ,myGrid.node[lbModel.G4 ] ,lbModel.G4_DV_ZERO_M1_M1  ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G3 ,myGrid.node[lbModel.G3 ] ,lbModel.G3_DV_ZERO_M1_P1  ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.node[lbModel.G6 ] ,lbModel.G6_DV_M1_M1_ZERO  ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.node[lbModel.G6 ] ,lbModel.G6_DV_P1_M1_ZERO  ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.node[lbModel.G6 ] ,lbModel.G6_DV_ZERO_M1_ZERO) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_M1_M1_M1   ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_P1_M1_M1   ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G9 ,myGrid.node[lbModel.G9 ] ,lbModel.G9_DV_M1_M1_P1    ) ;
    massLeak -=  myGrid(i+index,myGrid.nE2,k,lbModel.G9 ,myGrid.node[lbModel.G9 ] ,lbModel.G9_DV_P1_M1_P1    ) ;


    myGrid(i+index,myGrid.nE2,k, lbModel.G0, myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) += massLeak;

//     myGrid(i+index,myGrid.nE2,k,lbModel.G2 ,myGrid.node[lbModel.G2 ] ,lbModel.G2_DV_ZERO_M2_ZERO)   +=  (massLeak/10.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G4 ,myGrid.node[lbModel.G4 ] ,lbModel.G4_DV_ZERO_M1_M1  )   +=  (massLeak/10.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G3 ,myGrid.node[lbModel.G3 ] ,lbModel.G3_DV_ZERO_M1_P1  )   +=  (massLeak/10.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.node[lbModel.G6 ] ,lbModel.G6_DV_M1_M1_ZERO  )   +=  (massLeak/10.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.node[lbModel.G6 ] ,lbModel.G6_DV_P1_M1_ZERO  )   +=  (massLeak/10.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G6 ,myGrid.node[lbModel.G6 ] ,lbModel.G6_DV_ZERO_M1_ZERO)   +=  (massLeak/10.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_M1_M1_M1   )   +=  (massLeak/10.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_P1_M1_M1   )   +=  (massLeak/10.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G9 ,myGrid.node[lbModel.G9 ] ,lbModel.G9_DV_M1_M1_P1    )   +=  (massLeak/10.0) ;
//     myGrid(i+index,myGrid.nE2,k,lbModel.G9 ,myGrid.node[lbModel.G9 ] ,lbModel.G9_DV_P1_M1_P1    )   +=  (massLeak/10.0) ;
   }

   // This can be repalced with normal BB later
   for(int index=0; index<VECT_LENGTH; index++)
   {
    wallFlux = 0.0;
    wallFlux   += (lbModel.cy[lbModel.DV_ZERO_P2_ZERO] * myGrid(i+index,myGrid.nE2+2,k,lbModel.G1 ,0 ,lbModel.G1_DV_ZERO_P2_ZERO));

    massLeak = 0.0;
    massLeak   += myGrid(i+index,myGrid.nE2+2,k,lbModel.G1 ,0 ,lbModel.G1_DV_ZERO_P2_ZERO);

    wallFluxEq = 0.0;
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_P2_ZERO] * lbModel.fTemp1[index][lbModel.G1_DV_ZERO_P2_ZERO]);

    probability = wallFlux/wallFluxEq ;

    myGrid(i+index,myGrid.nE2-1,k,lbModel.G2 ,myGrid.node[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO)   =  probability * lbModel.fTemp1[index][lbModel.G1_DV_ZERO_P2_ZERO]  ;
    massLeak   -=  myGrid(i+index,myGrid.nE2-1,k,lbModel.G2 ,myGrid.node[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO) ;

    myGrid(i+index,myGrid.nE2-1,k, lbModel.G0, myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) += massLeak;
//     myGrid(i+index,myGrid.nE2-1,k,lbModel.G2 ,myGrid.node[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO) += massLeak;
   }
  }


}

template <int N,int numblock, typename dataType1>
void prepareDiffuseBottomWall(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid)
{
 for (int k = myGrid.nB3; k<= myGrid.nE3; k++)
  for (int i = myGrid.nB1; i <= myGrid.nE1; i++)
  {
   myGrid(i,myGrid.nB2-1,k,lbModel.G2 ,0,lbModel.G2_DV_ZERO_M2_ZERO) = myGrid(i,myGrid.nB2  ,k,lbModel.G2 ,myGrid.node[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO) ;
   myGrid(i,myGrid.nB2-2,k,lbModel.G2 ,0,lbModel.G2_DV_ZERO_M2_ZERO) = myGrid(i,myGrid.nB2+1,k,lbModel.G2 ,myGrid.node[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G4 ,0,lbModel.G4_DV_ZERO_M1_M1  ) = myGrid(i,myGrid.nB2  ,k,lbModel.G4 ,myGrid.node[lbModel.G4]  ,lbModel.G4_DV_ZERO_M1_M1  ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G3 ,0,lbModel.G3_DV_ZERO_M1_P1  ) = myGrid(i,myGrid.nB2  ,k,lbModel.G3 ,myGrid.node[lbModel.G3]  ,lbModel.G3_DV_ZERO_M1_P1  ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G6 ,0,lbModel.G6_DV_M1_M1_ZERO  ) = myGrid(i,myGrid.nB2  ,k,lbModel.G6 ,myGrid.node[lbModel.G6]  ,lbModel.G6_DV_M1_M1_ZERO  ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G6 ,0,lbModel.G6_DV_P1_M1_ZERO  ) = myGrid(i,myGrid.nB2  ,k,lbModel.G6 ,myGrid.node[lbModel.G6]  ,lbModel.G6_DV_P1_M1_ZERO  ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G6 ,0,lbModel.G6_DV_ZERO_M1_ZERO) = myGrid(i,myGrid.nB2  ,k,lbModel.G6 ,myGrid.node[lbModel.G6]  ,lbModel.G6_DV_ZERO_M1_ZERO) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G10,0,lbModel.G10_DV_M1_M1_M1   ) = myGrid(i,myGrid.nB2  ,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_M1_M1_M1   ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G10,0,lbModel.G10_DV_P1_M1_M1   ) = myGrid(i,myGrid.nB2  ,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_P1_M1_M1   ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G9 ,0,lbModel.G9_DV_M1_M1_P1    ) = myGrid(i,myGrid.nB2  ,k,lbModel.G9 ,myGrid.node[lbModel.G9]  ,lbModel.G9_DV_M1_M1_P1    ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G9 ,0,lbModel.G9_DV_P1_M1_P1    ) = myGrid(i,myGrid.nB2  ,k,lbModel.G9 ,myGrid.node[lbModel.G9]  ,lbModel.G9_DV_P1_M1_P1    ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G8 ,0,lbModel.G8_DV_M_M_M       ) = myGrid(i,myGrid.nB2  ,k,lbModel.G8,myGrid.node[lbModel.G8]   ,lbModel.G8_DV_M_M_M       ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G8 ,0,lbModel.G8_DV_P_M_M       ) = myGrid(i,myGrid.nB2  ,k,lbModel.G8,myGrid.node[lbModel.G8]   ,lbModel.G8_DV_P_M_M       ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G7 ,0,lbModel.G7_DV_M_M_P       ) = myGrid(i,myGrid.nB2  ,k,lbModel.G7,myGrid.node[lbModel.G7]   ,lbModel.G7_DV_M_M_P       ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G7 ,0,lbModel.G7_DV_P_M_P       ) = myGrid(i,myGrid.nB2  ,k,lbModel.G7,myGrid.node[lbModel.G7]   ,lbModel.G7_DV_P_M_P       ) ;
  }

 for (int k = myGrid.nB3; k<= myGrid.nE3; k++)
  for (int i = myGrid.nB1; i <= myGrid.nE1; i++)
  {
   myGrid(i,myGrid.nB2-1,k,lbModel.G2 ,1,lbModel.G2_DV_ZERO_M2_ZERO) =  myGrid(i,myGrid.nB2  ,k,lbModel.G2 ,myGrid.cell[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO) ;
   myGrid(i,myGrid.nB2-2,k,lbModel.G2 ,1,lbModel.G2_DV_ZERO_M2_ZERO) =  myGrid(i,myGrid.nB2+1,k,lbModel.G2 ,myGrid.cell[lbModel.G2]  ,lbModel.G2_DV_ZERO_M2_ZERO) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G4 ,1,lbModel.G4_DV_ZERO_M1_M1  ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G4 ,myGrid.cell[lbModel.G4]  ,lbModel.G4_DV_ZERO_M1_M1  ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G3 ,1,lbModel.G3_DV_ZERO_M1_P1  ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G3 ,myGrid.cell[lbModel.G3]  ,lbModel.G3_DV_ZERO_M1_P1  ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G6 ,1,lbModel.G6_DV_M1_M1_ZERO  ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G6 ,myGrid.cell[lbModel.G6]  ,lbModel.G6_DV_M1_M1_ZERO  ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G6 ,1,lbModel.G6_DV_P1_M1_ZERO  ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G6 ,myGrid.cell[lbModel.G6]  ,lbModel.G6_DV_P1_M1_ZERO  ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G6 ,1,lbModel.G6_DV_ZERO_M1_ZERO) =  myGrid(i,myGrid.nB2  ,k,lbModel.G6 ,myGrid.cell[lbModel.G6]  ,lbModel.G6_DV_ZERO_M1_ZERO) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G10,1,lbModel.G10_DV_M1_M1_M1   ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_M1_M1_M1   ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G10,1,lbModel.G10_DV_P1_M1_M1   ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_P1_M1_M1   ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G9 ,1,lbModel.G9_DV_M1_M1_P1    ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G9 ,myGrid.cell[lbModel.G9]  ,lbModel.G9_DV_M1_M1_P1    ) ;
   myGrid(i,myGrid.nB2-1,k,lbModel.G9 ,1,lbModel.G9_DV_P1_M1_P1    ) =  myGrid(i,myGrid.nB2  ,k,lbModel.G9 ,myGrid.cell[lbModel.G9]  ,lbModel.G9_DV_P1_M1_P1    ) ;
 }
}


template <int N,int numblock, typename dataType1>
void correctDiffuseF0massManipBottomWall(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int VECT_LENGTH, dataType1 tempBottom, double uRef)
{
    dataType1 *rhoWall,*uX,*uY,*uZ,*thetaBottom;
    rhoWall   = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uX    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uY    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    uZ    = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);
    thetaBottom = (dataType1*) _mm_malloc(VECT_LENGTH*sizeof(dataType1), 32);

 dataType1 wallFlux, wallFluxEq, probability, massLeak ;


 for(int i=0; i<VECT_LENGTH; i++)
 {
  thetaBottom[i] = tempBottom;
  uX[i]          = -uRef;
  uY[i]          = 0.0;
  uZ[i]          = 0.0;
  rhoWall[i]     = 1.0;
 }

                getFEqSIMD(lbModel,VECT_LENGTH,rhoWall,uX,uY,uZ,thetaBottom,10);
//    getIterativeFEq(lbModel,VECT_LENGTH,rhoWall,uX,uY,uZ,thetaBottom);
//               getFEqNewIterative(lbModel,VECT_LENGTH,rhoWall,uX,uY,uZ,thetaBottom);


 for (int k = myGrid.nB3; k<= myGrid.nE3; k++)
  for (int i = myGrid.nB1; i <= myGrid.nE1; i=i+VECT_LENGTH)
  {
   for(int index=0; index<VECT_LENGTH; index++)
   {
    wallFlux = 0.0;
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_M2_ZERO] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G2 ,0,lbModel.G2_DV_ZERO_M2_ZERO));
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_M1_M1  ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G4 ,0,lbModel.G4_DV_ZERO_M1_M1  ));
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_M1_P1  ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G3 ,0,lbModel.G3_DV_ZERO_M1_P1  ));
    wallFlux += (lbModel.cy[lbModel.DV_M1_M1_ZERO  ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G6 ,0,lbModel.G6_DV_M1_M1_ZERO  ));
    wallFlux += (lbModel.cy[lbModel.DV_P1_M1_ZERO  ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G6 ,0,lbModel.G6_DV_P1_M1_ZERO  ));
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_M1_ZERO] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G6 ,0,lbModel.G6_DV_ZERO_M1_ZERO));
    wallFlux += (lbModel.cy[lbModel.DV_M1_M1_M1    ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G10,0,lbModel.G10_DV_M1_M1_M1   ));
    wallFlux += (lbModel.cy[lbModel.DV_P1_M1_M1    ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G10,0,lbModel.G10_DV_P1_M1_M1   ));
    wallFlux += (lbModel.cy[lbModel.DV_M1_M1_P1    ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G9 ,0,lbModel.G9_DV_M1_M1_P1    ));
    wallFlux += (lbModel.cy[lbModel.DV_P1_M1_P1    ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G9 ,0,lbModel.G9_DV_P1_M1_P1    ));
    wallFlux += (lbModel.cy[lbModel.DV_M_M_M       ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G8 ,0,lbModel.G8_DV_M_M_M       ));
    wallFlux += (lbModel.cy[lbModel.DV_P_M_M       ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G8 ,0,lbModel.G8_DV_P_M_M       ));
    wallFlux += (lbModel.cy[lbModel.DV_M_M_P       ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G7 ,0,lbModel.G7_DV_M_M_P       ));
    wallFlux += (lbModel.cy[lbModel.DV_P_M_P       ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G7 ,0,lbModel.G7_DV_P_M_P       ));

    massLeak = 0.0;
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G2 ,0,lbModel.G2_DV_ZERO_M2_ZERO);
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G4 ,0,lbModel.G4_DV_ZERO_M1_M1  );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G3 ,0,lbModel.G3_DV_ZERO_M1_P1  );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G6 ,0,lbModel.G6_DV_M1_M1_ZERO  );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G6 ,0,lbModel.G6_DV_P1_M1_ZERO  );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G6 ,0,lbModel.G6_DV_ZERO_M1_ZERO);
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G10,0,lbModel.G10_DV_M1_M1_M1   );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G10,0,lbModel.G10_DV_P1_M1_M1   );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G9 ,0,lbModel.G9_DV_M1_M1_P1    );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G9 ,0,lbModel.G9_DV_P1_M1_P1    );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G8 ,0,lbModel.G8_DV_M_M_M       );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G8 ,0,lbModel.G8_DV_P_M_M       );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G7 ,0,lbModel.G7_DV_M_M_P       );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G7 ,0,lbModel.G7_DV_P_M_P       );




    wallFluxEq  = 0.0;
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_M2_ZERO]  *  lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_M2_ZERO]);
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_M1_M1  ]  *  lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_M1_M1  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_M1_P1  ]  *  lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_M1_P1  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M1_M1_ZERO  ]  *  lbModel.fTemp6 [index][lbModel.G6_DV_M1_M1_ZERO  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P1_M1_ZERO  ]  *  lbModel.fTemp6 [index][lbModel.G6_DV_P1_M1_ZERO  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_M1_ZERO]  *  lbModel.fTemp6 [index][lbModel.G6_DV_ZERO_M1_ZERO]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M1_M1_M1    ]  *  lbModel.fTemp10[index][lbModel.G10_DV_M1_M1_M1   ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P1_M1_M1    ]  *  lbModel.fTemp10[index][lbModel.G10_DV_P1_M1_M1   ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M1_M1_P1    ]  *  lbModel.fTemp9 [index][lbModel.G9_DV_M1_M1_P1    ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P1_M1_P1    ]  *  lbModel.fTemp9 [index][lbModel.G9_DV_P1_M1_P1    ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M_M_M       ]  *  lbModel.fTemp8 [index][lbModel.G8_DV_M_M_M       ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P_M_M       ]  *  lbModel.fTemp8 [index][lbModel.G8_DV_P_M_M       ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M_M_P       ]  *  lbModel.fTemp7 [index][lbModel.G7_DV_M_M_P       ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P_M_P       ]  *  lbModel.fTemp7 [index][lbModel.G7_DV_P_M_P       ]);

    probability = wallFlux/wallFluxEq ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G1 ,myGrid.node[lbModel.G1 ] ,lbModel.G1_DV_ZERO_P2_ZERO)   =  probability * lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_M2_ZERO]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G3 ,myGrid.node[lbModel.G3 ] ,lbModel.G3_DV_ZERO_P1_P1  )   =  probability * lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_M1_M1  ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G4 ,myGrid.node[lbModel.G4 ] ,lbModel.G4_DV_ZERO_P1_M1  )   =  probability * lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_M1_P1  ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.node[lbModel.G5 ] ,lbModel.G5_DV_P1_P1_ZERO  )   =  probability * lbModel.fTemp6 [index][lbModel.G6_DV_M1_M1_ZERO  ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.node[lbModel.G5 ] ,lbModel.G5_DV_M1_P1_ZERO  )   =  probability * lbModel.fTemp6 [index][lbModel.G6_DV_P1_M1_ZERO  ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.node[lbModel.G5 ] ,lbModel.G5_DV_ZERO_P1_ZERO)   =  probability * lbModel.fTemp6 [index][lbModel.G6_DV_ZERO_M1_ZERO]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G9 ,myGrid.node[lbModel.G9 ] ,lbModel.G9_DV_P1_P1_P1    )   =  probability * lbModel.fTemp10[index][lbModel.G10_DV_M1_M1_M1   ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G9 ,myGrid.node[lbModel.G9 ] ,lbModel.G9_DV_M1_P1_P1    )   =  probability * lbModel.fTemp10[index][lbModel.G10_DV_P1_M1_M1   ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_P1_P1_M1   )   =  probability * lbModel.fTemp9 [index][lbModel.G9_DV_M1_M1_P1    ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_M1_P1_M1   )   =  probability * lbModel.fTemp9 [index][lbModel.G9_DV_P1_M1_P1    ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G7 ,myGrid.node[lbModel.G7 ] ,lbModel.G7_DV_P_P_P       )   =  probability * lbModel.fTemp8 [index][lbModel.G8_DV_M_M_M       ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G7 ,myGrid.node[lbModel.G7 ] ,lbModel.G7_DV_M_P_P       )   =  probability * lbModel.fTemp8 [index][lbModel.G8_DV_P_M_M       ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G8 ,myGrid.node[lbModel.G8 ] ,lbModel.G8_DV_P_P_M       )   =  probability * lbModel.fTemp7 [index][lbModel.G7_DV_M_M_P       ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G8 ,myGrid.node[lbModel.G8 ] ,lbModel.G8_DV_M_P_M       )   =  probability * lbModel.fTemp7 [index][lbModel.G7_DV_P_M_P       ]  ;

    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G1 ,myGrid.node[lbModel.G1 ] ,lbModel.G1_DV_ZERO_P2_ZERO);
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G3 ,myGrid.node[lbModel.G3 ] ,lbModel.G3_DV_ZERO_P1_P1  );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G4 ,myGrid.node[lbModel.G4 ] ,lbModel.G4_DV_ZERO_P1_M1  );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.node[lbModel.G5 ] ,lbModel.G5_DV_P1_P1_ZERO  );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.node[lbModel.G5 ] ,lbModel.G5_DV_M1_P1_ZERO  );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.node[lbModel.G5 ] ,lbModel.G5_DV_ZERO_P1_ZERO);
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G9 ,myGrid.node[lbModel.G9 ] ,lbModel.G9_DV_P1_P1_P1    );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G9 ,myGrid.node[lbModel.G9 ] ,lbModel.G9_DV_M1_P1_P1    );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_P1_P1_M1   );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_M1_P1_M1   );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G7 ,myGrid.node[lbModel.G7 ] ,lbModel.G7_DV_P_P_P       );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G7 ,myGrid.node[lbModel.G7 ] ,lbModel.G7_DV_M_P_P       );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G8 ,myGrid.node[lbModel.G8 ] ,lbModel.G8_DV_P_P_M       );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G8 ,myGrid.node[lbModel.G8 ] ,lbModel.G8_DV_M_P_M       );

    myGrid(i+index,myGrid.nB2,k, lbModel.G0, myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) += massLeak;

//     myGrid(i+index,myGrid.nB2,k,lbModel.G1 ,myGrid.node[lbModel.G1 ] ,lbModel.G1_DV_ZERO_P2_ZERO)   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G3 ,myGrid.node[lbModel.G3 ] ,lbModel.G3_DV_ZERO_P1_P1  )   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G4 ,myGrid.node[lbModel.G4 ] ,lbModel.G4_DV_ZERO_P1_M1  )   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.node[lbModel.G5 ] ,lbModel.G5_DV_P1_P1_ZERO  )   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.node[lbModel.G5 ] ,lbModel.G5_DV_M1_P1_ZERO  )   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.node[lbModel.G5 ] ,lbModel.G5_DV_ZERO_P1_ZERO)   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G9 ,myGrid.node[lbModel.G9 ] ,lbModel.G9_DV_P1_P1_P1    )   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G9 ,myGrid.node[lbModel.G9 ] ,lbModel.G9_DV_M1_P1_P1    )   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_P1_P1_M1   )   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G10,myGrid.node[lbModel.G10] ,lbModel.G10_DV_M1_P1_M1   )   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G7 ,myGrid.node[lbModel.G7 ] ,lbModel.G7_DV_P_P_P       )   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G7 ,myGrid.node[lbModel.G7 ] ,lbModel.G7_DV_M_P_P       )   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G8 ,myGrid.node[lbModel.G8 ] ,lbModel.G8_DV_P_P_M       )   += (massLeak/14.0);
//     myGrid(i+index,myGrid.nB2,k,lbModel.G8 ,myGrid.node[lbModel.G8 ] ,lbModel.G8_DV_M_P_M       )   += (massLeak/14.0);
//
   }



   // This can be repalced with normal BB later
   for(int index=0; index<VECT_LENGTH; index++)
   {
    wallFlux = 0.0;
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_M2_ZERO] * myGrid(i+index,myGrid.nB2-2,k,lbModel.G2 ,0 ,lbModel.G2_DV_ZERO_M2_ZERO));

    massLeak = 0.0;
    massLeak += myGrid(i+index,myGrid.nB2-2,k,lbModel.G2 ,0 ,lbModel.G2_DV_ZERO_M2_ZERO);

    wallFluxEq = 0.0;
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_M2_ZERO]  *  lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_M2_ZERO]);

    probability = wallFlux/wallFluxEq ;

    myGrid(i+index,myGrid.nB2+1,k,lbModel.G1 ,myGrid.node[lbModel.G1 ] ,lbModel.G1_DV_ZERO_P2_ZERO)   =  probability * lbModel.fTemp1[index][lbModel.G1_DV_ZERO_P2_ZERO]  ;
    massLeak   -=  myGrid(i+index,myGrid.nB2+1,k,lbModel.G1 ,myGrid.node[lbModel.G1 ] ,lbModel.G1_DV_ZERO_P2_ZERO) ;

    myGrid(i+index,myGrid.nB2+1,k, lbModel.G0, myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) += massLeak;
//     myGrid(i+index,myGrid.nB2+1,k,lbModel.G1 ,myGrid.node[lbModel.G1 ] ,lbModel.G1_DV_ZERO_P2_ZERO)  += massLeak;
   }
 }


 for (int k = myGrid.nB3; k<= myGrid.nE3; k++)
  for (int i = myGrid.nB1; i <= myGrid.nE1; i=i+VECT_LENGTH)
  {
   for(int index=0; index<VECT_LENGTH; index++)
   {
    wallFlux = 0.0;
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_M2_ZERO] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G2 ,1,lbModel.G2_DV_ZERO_M2_ZERO));
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_M1_M1  ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G4 ,1,lbModel.G4_DV_ZERO_M1_M1  ));
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_M1_P1  ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G3 ,1,lbModel.G3_DV_ZERO_M1_P1  ));
    wallFlux += (lbModel.cy[lbModel.DV_M1_M1_ZERO  ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G6 ,1,lbModel.G6_DV_M1_M1_ZERO  ));
    wallFlux += (lbModel.cy[lbModel.DV_P1_M1_ZERO  ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G6 ,1,lbModel.G6_DV_P1_M1_ZERO  ));
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_M1_ZERO] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G6 ,1,lbModel.G6_DV_ZERO_M1_ZERO));
    wallFlux += (lbModel.cy[lbModel.DV_M1_M1_M1    ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G10,1,lbModel.G10_DV_M1_M1_M1   ));
    wallFlux += (lbModel.cy[lbModel.DV_P1_M1_M1    ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G10,1,lbModel.G10_DV_P1_M1_M1   ));
    wallFlux += (lbModel.cy[lbModel.DV_M1_M1_P1    ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G9 ,1,lbModel.G9_DV_M1_M1_P1    ));
    wallFlux += (lbModel.cy[lbModel.DV_P1_M1_P1    ] * myGrid(i+index,myGrid.nB2-1,k,lbModel.G9 ,1,lbModel.G9_DV_P1_M1_P1    ));

    massLeak = 0.0;
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G2 ,1,lbModel.G2_DV_ZERO_M2_ZERO);
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G4 ,1,lbModel.G4_DV_ZERO_M1_M1  );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G3 ,1,lbModel.G3_DV_ZERO_M1_P1  );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G6 ,1,lbModel.G6_DV_M1_M1_ZERO  );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G6 ,1,lbModel.G6_DV_P1_M1_ZERO  );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G6 ,1,lbModel.G6_DV_ZERO_M1_ZERO);
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G10,1,lbModel.G10_DV_M1_M1_M1   );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G10,1,lbModel.G10_DV_P1_M1_M1   );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G9 ,1,lbModel.G9_DV_M1_M1_P1    );
    massLeak += myGrid(i+index,myGrid.nB2-1,k,lbModel.G9 ,1,lbModel.G9_DV_P1_M1_P1    );



    wallFluxEq  = 0.0;
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_M2_ZERO]  *  lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_M2_ZERO]);
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_M1_M1  ]  *  lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_M1_M1  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_M1_P1  ]  *  lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_M1_P1  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M1_M1_ZERO  ]  *  lbModel.fTemp6 [index][lbModel.G6_DV_M1_M1_ZERO  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P1_M1_ZERO  ]  *  lbModel.fTemp6 [index][lbModel.G6_DV_P1_M1_ZERO  ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_M1_ZERO]  *  lbModel.fTemp6 [index][lbModel.G6_DV_ZERO_M1_ZERO]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M1_M1_M1    ]  *  lbModel.fTemp10[index][lbModel.G10_DV_M1_M1_M1   ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P1_M1_M1    ]  *  lbModel.fTemp10[index][lbModel.G10_DV_P1_M1_M1   ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_M1_M1_P1    ]  *  lbModel.fTemp9 [index][lbModel.G9_DV_M1_M1_P1    ]);
    wallFluxEq += (lbModel.cy[lbModel.DV_P1_M1_P1    ]  *  lbModel.fTemp9 [index][lbModel.G9_DV_P1_M1_P1    ]);

    probability = wallFlux/wallFluxEq ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G1 ,myGrid.cell[lbModel.G1 ] ,lbModel.G1_DV_ZERO_P2_ZERO)   =  probability * lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_M2_ZERO]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G3 ,myGrid.cell[lbModel.G3 ] ,lbModel.G3_DV_ZERO_P1_P1  )   =  probability * lbModel.fTemp4 [index][lbModel.G4_DV_ZERO_M1_M1  ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G4 ,myGrid.cell[lbModel.G4 ] ,lbModel.G4_DV_ZERO_P1_M1  )   =  probability * lbModel.fTemp3 [index][lbModel.G3_DV_ZERO_M1_P1  ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.cell[lbModel.G5 ] ,lbModel.G5_DV_P1_P1_ZERO  )   =  probability * lbModel.fTemp6 [index][lbModel.G6_DV_M1_M1_ZERO  ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.cell[lbModel.G5 ] ,lbModel.G5_DV_M1_P1_ZERO  )   =  probability * lbModel.fTemp6 [index][lbModel.G6_DV_P1_M1_ZERO  ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.cell[lbModel.G5 ] ,lbModel.G5_DV_ZERO_P1_ZERO)   =  probability * lbModel.fTemp6 [index][lbModel.G6_DV_ZERO_M1_ZERO]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G9 ,myGrid.cell[lbModel.G9 ] ,lbModel.G9_DV_P1_P1_P1    )   =  probability * lbModel.fTemp10[index][lbModel.G10_DV_M1_M1_M1   ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G9 ,myGrid.cell[lbModel.G9 ] ,lbModel.G9_DV_M1_P1_P1    )   =  probability * lbModel.fTemp10[index][lbModel.G10_DV_P1_M1_M1   ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_P1_P1_M1   )   =  probability * lbModel.fTemp9 [index][lbModel.G9_DV_M1_M1_P1    ]  ;
    myGrid(i+index,myGrid.nB2,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_M1_P1_M1   )   =  probability * lbModel.fTemp9 [index][lbModel.G9_DV_P1_M1_P1    ]  ;

    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G1 ,myGrid.cell[lbModel.G1 ] ,lbModel.G1_DV_ZERO_P2_ZERO);
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G3 ,myGrid.cell[lbModel.G3 ] ,lbModel.G3_DV_ZERO_P1_P1  );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G4 ,myGrid.cell[lbModel.G4 ] ,lbModel.G4_DV_ZERO_P1_M1  );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.cell[lbModel.G5 ] ,lbModel.G5_DV_P1_P1_ZERO  );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.cell[lbModel.G5 ] ,lbModel.G5_DV_M1_P1_ZERO  );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.cell[lbModel.G5 ] ,lbModel.G5_DV_ZERO_P1_ZERO);
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G9 ,myGrid.cell[lbModel.G9 ] ,lbModel.G9_DV_P1_P1_P1    );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G9 ,myGrid.cell[lbModel.G9 ] ,lbModel.G9_DV_M1_P1_P1    );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_P1_P1_M1   );
    massLeak -= myGrid(i+index,myGrid.nB2,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_M1_P1_M1   );

    myGrid(i+index,myGrid.nB2,k, lbModel.G0, myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) += massLeak;

/*    myGrid(i+index,myGrid.nB2,k,lbModel.G1 ,myGrid.cell[lbModel.G1 ] ,lbModel.G1_DV_ZERO_P2_ZERO) += (massLeak/10.0);
    myGrid(i+index,myGrid.nB2,k,lbModel.G3 ,myGrid.cell[lbModel.G3 ] ,lbModel.G3_DV_ZERO_P1_P1  ) += (massLeak/10.0);
    myGrid(i+index,myGrid.nB2,k,lbModel.G4 ,myGrid.cell[lbModel.G4 ] ,lbModel.G4_DV_ZERO_P1_M1  ) += (massLeak/10.0);
    myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.cell[lbModel.G5 ] ,lbModel.G5_DV_P1_P1_ZERO  ) += (massLeak/10.0);
    myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.cell[lbModel.G5 ] ,lbModel.G5_DV_M1_P1_ZERO  ) += (massLeak/10.0);
    myGrid(i+index,myGrid.nB2,k,lbModel.G5 ,myGrid.cell[lbModel.G5 ] ,lbModel.G5_DV_ZERO_P1_ZERO) += (massLeak/10.0);
    myGrid(i+index,myGrid.nB2,k,lbModel.G9 ,myGrid.cell[lbModel.G9 ] ,lbModel.G9_DV_P1_P1_P1    ) += (massLeak/10.0);
    myGrid(i+index,myGrid.nB2,k,lbModel.G9 ,myGrid.cell[lbModel.G9 ] ,lbModel.G9_DV_M1_P1_P1    ) += (massLeak/10.0);
    myGrid(i+index,myGrid.nB2,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_P1_P1_M1   ) += (massLeak/10.0);
    myGrid(i+index,myGrid.nB2,k,lbModel.G10,myGrid.cell[lbModel.G10] ,lbModel.G10_DV_M1_P1_M1   ) += (massLeak/10.0);  */
   }

   // This can be repalced with normal BB later
   for(int index=0; index<VECT_LENGTH; index++)
   {
    wallFlux = 0.0;
    wallFlux += (lbModel.cy[lbModel.DV_ZERO_M2_ZERO] * myGrid(i+index,myGrid.nB2-2,k,lbModel.G2 ,1 ,lbModel.G2_DV_ZERO_M2_ZERO));

    massLeak = 0.0;
    massLeak += myGrid(i+index,myGrid.nB2-2,k,lbModel.G2 ,1 ,lbModel.G2_DV_ZERO_M2_ZERO);


    wallFluxEq = 0.0;
    wallFluxEq += (lbModel.cy[lbModel.DV_ZERO_M2_ZERO]  *  lbModel.fTemp2 [index][lbModel.G2_DV_ZERO_M2_ZERO]);

    probability = wallFlux/wallFluxEq ;

    myGrid(i+index,myGrid.nB2+1,k,lbModel.G1 ,myGrid.cell[lbModel.G1 ] ,lbModel.G1_DV_ZERO_P2_ZERO)   =  probability * lbModel.fTemp1[index][lbModel.G1_DV_ZERO_P2_ZERO]  ;
    massLeak -= myGrid(i+index,myGrid.nB2+1,k,lbModel.G1 ,myGrid.cell[lbModel.G1 ] ,lbModel.G1_DV_ZERO_P2_ZERO);
    myGrid(i+index,myGrid.nB2+1,k, lbModel.G0, myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) += massLeak;
//     myGrid(i+index,myGrid.nB2+1,k,lbModel.G1 ,myGrid.cell[lbModel.G1 ] ,lbModel.G1_DV_ZERO_P2_ZERO) += massLeak;


   }
 }
}


 template<typename dataType1>
 void getHydroMomentSinglePoint(lbmRD3Q41<dataType1> &lbModel, dataType1 &rho, dataType1 &uX, dataType1 &uY, dataType1 &uZ, dataType1 &theta)
 {
    dataType1 sum(0.0);
    rho=0.0;uX=0.0;uY=0.0;uZ=0.0;theta=0.0;

      rho = lbModel.fTemp0[0][lbModel.CENTER_DV_ZERO_ZERO_ZERO];


    // G1 and G2

      sum    = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO]+lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
      rho    += sum;
      theta  = 4.0*sum;
      uX     = 2.0*(lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO]-lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO]);

      sum       = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO]+lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO];
      rho   += sum;
      theta += 4.0*sum;
      uY     = 2.0*(lbModel.fTemp1[0][lbModel.G1_DV_ZERO_P2_ZERO]-lbModel.fTemp2[0][lbModel.G2_DV_ZERO_M2_ZERO]);

      sum       = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2]+lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2];
      rho   += sum;
      theta += 4.0*sum;
      uZ     = 2.0*(lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P2]-lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M2]);

      sum       = lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1]+lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];
      rho   += sum;
      theta += sum;
      uZ    += lbModel.fTemp1[0][lbModel.G1_DV_ZERO_ZERO_P1]-lbModel.fTemp2[0][lbModel.G2_DV_ZERO_ZERO_M1];

    // G3 and G4
      sum       = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1]+lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1];
      rho   += sum;
      theta += 2.0*sum;
      uY    += (lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1]-lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1]);
      uZ    += (lbModel.fTemp3[0][lbModel.G3_DV_ZERO_P1_P1]-lbModel.fTemp4[0][lbModel.G4_DV_ZERO_M1_M1]);

      sum       = lbModel.fTemp3[0][lbModel.G3_DV_ZERO_M1_P1 ]+lbModel.fTemp4[0][lbModel.G4_DV_ZERO_P1_M1];
      rho   += sum;
      theta += 2.0*sum;
      uY    -= (lbModel.fTemp3[0][lbModel.G3_DV_ZERO_M1_P1]-lbModel.fTemp4[0][lbModel.G4_DV_ZERO_P1_M1]);
      uZ    += (lbModel.fTemp3[0][lbModel.G3_DV_ZERO_M1_P1]-lbModel.fTemp4[0][lbModel.G4_DV_ZERO_P1_M1]);

      sum       = lbModel.fTemp3[0][lbModel.G3_DV_P1_ZERO_P1]+lbModel.fTemp4[0][lbModel.G4_DV_M1_ZERO_M1];
      rho   += sum;
      theta += 2.0*sum;
      uX    += (lbModel.fTemp3[0][lbModel.G3_DV_P1_ZERO_P1]-lbModel.fTemp4[0][lbModel.G4_DV_M1_ZERO_M1]);
      uZ    += (lbModel.fTemp3[0][lbModel.G3_DV_P1_ZERO_P1]-lbModel.fTemp4[0][lbModel.G4_DV_M1_ZERO_M1]);

      sum       = lbModel.fTemp3[0][lbModel.G3_DV_M1_ZERO_P1]+lbModel.fTemp4[0][lbModel.G4_DV_P1_ZERO_M1];
      rho   += sum;
      theta += 2.0*sum;
      uX    -= (lbModel.fTemp3[0][lbModel.G3_DV_M1_ZERO_P1]-lbModel.fTemp4[0][lbModel.G4_DV_P1_ZERO_M1]);
      uZ    += (lbModel.fTemp3[0][lbModel.G3_DV_M1_ZERO_P1]-lbModel.fTemp4[0][lbModel.G4_DV_P1_ZERO_M1]);


    // G5 and G6
      sum       = lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO]+lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO];
      rho   += sum;
      theta += 2.0*sum;
      uX    += (lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO]-lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO]);
      uY    += (lbModel.fTemp5[0][lbModel.G5_DV_P1_P1_ZERO]-lbModel.fTemp6[0][lbModel.G6_DV_M1_M1_ZERO]);

      sum       = lbModel.fTemp5[0][lbModel.G5_DV_M1_P1_ZERO]+lbModel.fTemp6[0][lbModel.G6_DV_P1_M1_ZERO];
      rho   += sum;
      theta += 2.0*sum;
      uX    += (lbModel.fTemp6[0][lbModel.G6_DV_P1_M1_ZERO]-lbModel.fTemp5[0][lbModel.G5_DV_M1_P1_ZERO]);
      uY    -= (lbModel.fTemp6[0][lbModel.G6_DV_P1_M1_ZERO]-lbModel.fTemp5[0][lbModel.G5_DV_M1_P1_ZERO]);

      sum      = lbModel.fTemp5[0][lbModel.G5_DV_ZERO_P1_ZERO]+lbModel.fTemp6[0][lbModel.G6_DV_ZERO_M1_ZERO];
      rho   += sum;
      theta += sum;
      uY    += (lbModel.fTemp5[0][lbModel.G5_DV_ZERO_P1_ZERO]-lbModel.fTemp6[0][lbModel.G6_DV_ZERO_M1_ZERO]);

      sum       = lbModel.fTemp5[0][lbModel.G5_DV_P1_ZERO_ZERO]+lbModel.fTemp6[0][lbModel.G6_DV_M1_ZERO_ZERO];
      rho   += sum;
      theta += sum;
      uX    += (lbModel.fTemp5[0][lbModel.G5_DV_P1_ZERO_ZERO]-lbModel.fTemp6[0][lbModel.G6_DV_M1_ZERO_ZERO]);


    // G7 and G8
    const dataType1 threeByFour=3.0/4.0;

      sum       = lbModel.fTemp7[0][lbModel.G7_DV_P_P_P]+lbModel.fTemp8[0][lbModel.G8_DV_M_M_M];
      rho   += sum;
      theta += threeByFour*sum ;
      uX    += 0.5*(lbModel.fTemp7[0][lbModel.G7_DV_P_P_P]-lbModel.fTemp8[0][lbModel.G8_DV_M_M_M]);
      uY    += 0.5*(lbModel.fTemp7[0][lbModel.G7_DV_P_P_P]-lbModel.fTemp8[0][lbModel.G8_DV_M_M_M]);
      uZ    += 0.5*(lbModel.fTemp7[0][lbModel.G7_DV_P_P_P]-lbModel.fTemp8[0][lbModel.G8_DV_M_M_M]);

      sum       = lbModel.fTemp7[0][lbModel.G7_DV_M_P_P]+lbModel.fTemp8[0][lbModel.G8_DV_P_M_M];
      rho   += sum;
      theta += threeByFour*sum ;
      uX    -= 0.5*(lbModel.fTemp7[0][lbModel.G7_DV_M_P_P]-lbModel.fTemp8[0][lbModel.G8_DV_P_M_M]);
      uY    += 0.5*(lbModel.fTemp7[0][lbModel.G7_DV_M_P_P]-lbModel.fTemp8[0][lbModel.G8_DV_P_M_M]);
      uZ    += 0.5*(lbModel.fTemp7[0][lbModel.G7_DV_M_P_P]-lbModel.fTemp8[0][lbModel.G8_DV_P_M_M]);


      sum       = lbModel.fTemp7[0][lbModel.G7_DV_M_M_P]+lbModel.fTemp8[0][lbModel.G8_DV_P_P_M];
      rho   += sum;
      theta += threeByFour*sum ;
      uX    -= 0.5*(lbModel.fTemp7[0][lbModel.G7_DV_M_M_P]-lbModel.fTemp8[0][lbModel.G8_DV_P_P_M]);
      uY    -= 0.5*(lbModel.fTemp7[0][lbModel.G7_DV_M_M_P]-lbModel.fTemp8[0][lbModel.G8_DV_P_P_M]);
      uZ    += 0.5*(lbModel.fTemp7[0][lbModel.G7_DV_M_M_P]-lbModel.fTemp8[0][lbModel.G8_DV_P_P_M]);

      sum       = lbModel.fTemp7[0][lbModel.G7_DV_P_M_P]+lbModel.fTemp8[0][lbModel.G8_DV_M_P_M];
      rho   += sum;
      theta += threeByFour*sum ;
      uX    += 0.5*(lbModel.fTemp7[0][lbModel.G7_DV_P_M_P]-lbModel.fTemp8[0][lbModel.G8_DV_M_P_M ]);
      uY    -= 0.5*(lbModel.fTemp7[0][lbModel.G7_DV_P_M_P]-lbModel.fTemp8[0][lbModel.G8_DV_M_P_M ]);
      uZ    += 0.5*(lbModel.fTemp7[0][lbModel.G7_DV_P_M_P]-lbModel.fTemp8[0][lbModel.G8_DV_M_P_M ]);

    // G9 and G10

      sum       = lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1]+lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1];
      rho   += sum;
      theta += 3.0*sum ;
      uX    += (lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1]-lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1]);
      uY    += (lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1]-lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1]);
      uZ    += (lbModel.fTemp9[0][lbModel.G9_DV_P1_P1_P1]-lbModel.fTemp10[0][lbModel.G10_DV_M1_M1_M1]);

      sum       = lbModel.fTemp9[0][lbModel.G9_DV_M1_P1_P1]+lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1];
      rho   += sum;
      theta += 3.0*sum ;
      uX    -= (lbModel.fTemp9[0][lbModel.G9_DV_M1_P1_P1]-lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1]);
      uY    += (lbModel.fTemp9[0][lbModel.G9_DV_M1_P1_P1]-lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1]);
      uZ    += (lbModel.fTemp9[0][lbModel.G9_DV_M1_P1_P1]-lbModel.fTemp10[0][lbModel.G10_DV_P1_M1_M1]);


      sum       = lbModel.fTemp9[0][lbModel.G9_DV_M1_M1_P1]+lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1];
      rho   += sum;
      theta += 3.0*sum ;
      uX    -= ( lbModel.fTemp9[0][lbModel.G9_DV_M1_M1_P1]-lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1]);
      uY    -= ( lbModel.fTemp9[0][lbModel.G9_DV_M1_M1_P1]-lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1]);
      uZ    += ( lbModel.fTemp9[0][lbModel.G9_DV_M1_M1_P1]-lbModel.fTemp10[0][lbModel.G10_DV_P1_P1_M1]);

      sum       = lbModel.fTemp9[0][lbModel.G9_DV_P1_M1_P1 ]+lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1];
      rho   += sum;
      theta += 3.0*sum ;
      uX    += (lbModel.fTemp9[0][lbModel.G9_DV_P1_M1_P1]-lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1]);
      uY    -= (lbModel.fTemp9[0][lbModel.G9_DV_P1_M1_P1]-lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1]);
      uZ    += (lbModel.fTemp9[0][lbModel.G9_DV_P1_M1_P1]-lbModel.fTemp10[0][lbModel.G10_DV_M1_P1_M1]);

      dataType1  oneByRho = 1.0/rho;
      dataType1  oneByThree = 1.0/3.0;
      uX  *= oneByRho;
      uY  *= oneByRho;
      uZ  *= oneByRho;
      theta = oneByThree * (theta - rho*(uX*uX + uY*uY + uZ*uZ));
      theta *= oneByRho;
  }


template <int N,int numblock, typename dataType1>
void getTempProfile(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int step, double visc, double delT,  int myRank,  int time)
{
  dataType1 rho, uX, uZ, uY, theta, y;

  int XCENTER = (myGrid.nB1 + myGrid.nE1)*0.5 ;
  int ZCENTER = (myGrid.nB3 + myGrid.nE3)*0.5 ;
  if(myRank==0)
   std::cout << "\n y and theta at diffusion time=" << time;

 std::ofstream file;
 char fileName[150];
 sprintf(fileName,"./results/temperature_%d.txt",step);
 file.open(fileName);

  for(int i2=myGrid.nB2; i2<=myGrid.nE2-1;i2++)
  {
    copyFromNodeSinglePoint(lbModel,myGrid,XCENTER,i2,ZCENTER,0);
   getHydroMomentSinglePoint(lbModel,rho,uX,uY,uZ,theta);
   file<<((double)(i2-myGrid.nB2)/(double)(myGrid.m2))<<"   "<<rho<<"   "<<uY/std::sqrt((5.0/3.0)*lbModel.theta0)<<"  "<<theta/(lbModel.theta0)<<"   "<<rho*theta/(lbModel.theta0)<<"   "<<(theta-lbModel.theta0)/(lbModel.theta0)<<"  "<<i2-myGrid.nB2 <<std::endl;
   //                             1                              2                3                                4                              5                                  6
  }

  file.close();
}

template <int N,int numblock, typename dataType1>
void globalMass(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int step, int myRank)
 {
        dataType1 nodeMass(0.0),cellMass(0.0);

        int nodeCount(0),cellCount(0);
        int nodeCountGlobal(0);
        int cellCountGlobal(0);
        dataType1 nodeMassGlobal(0.0);
        dataType1 cellMassGlobal(0.0);

        for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
            for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
                for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
                    {
                        nodeMass +=  ( myGrid(i1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO)
                                     + myGrid(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1)
                                     + myGrid(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P2)
                                     + myGrid(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_P2_ZERO)
                                     + myGrid(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO)
                                     + myGrid(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1)
                                     + myGrid(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M2)
                                     + myGrid(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_M2_ZERO)
                                     + myGrid(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO)
                                     + myGrid(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1)
                                     + myGrid(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],lbModel.G3_DV_ZERO_M1_P1)
                                     + myGrid(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],lbModel.G3_DV_P1_ZERO_P1)
                                     + myGrid(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],lbModel.G3_DV_M1_ZERO_P1)
                                     + myGrid(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],lbModel.G4_DV_M1_ZERO_M1)
                                     + myGrid(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],lbModel.G4_DV_ZERO_P1_M1)
                                     + myGrid(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1)
                                     + myGrid(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],lbModel.G4_DV_P1_ZERO_M1)
                                     + myGrid(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO  )
                                     + myGrid(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],lbModel.G5_DV_M1_P1_ZERO  )
                                     + myGrid(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],lbModel.G5_DV_ZERO_P1_ZERO)
                                     + myGrid(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],lbModel.G5_DV_P1_ZERO_ZERO)
                                     + myGrid(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO  )
                                     + myGrid(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],lbModel.G6_DV_P1_M1_ZERO  )
                                     + myGrid(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],lbModel.G6_DV_ZERO_M1_ZERO)
                                     + myGrid(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],lbModel.G6_DV_M1_ZERO_ZERO)
                                     + myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P)
                                     + myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P)
                                     + myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P)
                                     + myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P)
                                     + myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M)
                                     + myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M)
                                     + myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M)
                                     + myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M)
                                     + myGrid(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_P1_P1)
                                     + myGrid(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],lbModel.G9_DV_M1_P1_P1)
                                     + myGrid(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],lbModel.G9_DV_M1_M1_P1)
                                     + myGrid(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_M1_P1)
                                     + myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1)
                                     + myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1)
                                     + myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1)
                                     + myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1));
                                     nodeCount++;
                    }

                    for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
                        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
                            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
                                {
                                    cellMass += (    myGrid(i1,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO)
                                                   + myGrid(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1)
                                                   + myGrid(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P2)
                                                   + myGrid(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_P2_ZERO)
                                                   + myGrid(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO)
                                                   + myGrid(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1)
                                                   + myGrid(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M2)
                                                   + myGrid(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_M2_ZERO)
                                                   + myGrid(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO)
                                                   + myGrid(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1)
                                                   + myGrid(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],lbModel.G3_DV_ZERO_M1_P1)
                                                   + myGrid(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],lbModel.G3_DV_P1_ZERO_P1)
                                                   + myGrid(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],lbModel.G3_DV_M1_ZERO_P1)
                                                   + myGrid(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],lbModel.G4_DV_M1_ZERO_M1)
                                                   + myGrid(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],lbModel.G4_DV_ZERO_P1_M1)
                                                   + myGrid(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1)
                                                   + myGrid(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],lbModel.G4_DV_P1_ZERO_M1)
                                                   + myGrid(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO  )
                                                   + myGrid(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],lbModel.G5_DV_M1_P1_ZERO  )
                                                   + myGrid(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],lbModel.G5_DV_ZERO_P1_ZERO)
                                                   + myGrid(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],lbModel.G5_DV_P1_ZERO_ZERO)
                                                   + myGrid(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO  )
                                                   + myGrid(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],lbModel.G6_DV_P1_M1_ZERO  )
                                                   + myGrid(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],lbModel.G6_DV_ZERO_M1_ZERO)
                                                   + myGrid(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],lbModel.G6_DV_M1_ZERO_ZERO)
                                                   + myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P)
                                                   + myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P)
                                                   + myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P)
                                                   + myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P)
                                                   + myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M)
                                                   + myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M)
                                                   + myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M)
                                                   + myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M)
                                                   + myGrid(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_P1_P1)
                                                   + myGrid(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],lbModel.G9_DV_M1_P1_P1)
                                                   + myGrid(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],lbModel.G9_DV_M1_M1_P1)
                                                   + myGrid(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_M1_P1)
                                                   + myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1)
                                                   + myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1)
                                                   + myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1)
                                                   + myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1));
                                                   cellCount++;
                                }

        MPI_Reduce(&nodeCount,&nodeCountGlobal, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&cellCount,&cellCountGlobal, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&nodeMass ,&nodeMassGlobal, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&cellMass ,&cellMassGlobal, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);

        if(myRank==0)
        printf("\n %d \t %.16lf \t %.16lf \t %.16lf",step,nodeMassGlobal,cellMassGlobal,nodeMassGlobal+cellMassGlobal);

    }


template void prepareDiffuseTopWall<4, 11, double>              (lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &);
template void correctDiffuseF0massManipTopWall<4, 11, double>   (lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,int,double,double);
template void prepareDiffuseBottomWall<4, 11, double>           (lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &);
template void correctDiffuseF0massManipBottomWall<4, 11, double>(lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,int,double,double);
template void getHydroMomentSinglePoint<double>                 (lbmRD3Q41<double> &, double &, double &, double &, double &, double &);
template void getTempProfile<4, 11, double>                     (lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,int ,int , double , double ,  int ,  int );
template void globalMass<4, 11, double>                         (lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,int ,int , int );
