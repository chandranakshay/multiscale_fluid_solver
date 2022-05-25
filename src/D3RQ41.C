#include"D3RQ41.h"

template <typename dataType1>
void setModelParameters(lbmRD3Q41<dataType1> &lbModel)
{
    lbModel.theta0      = 0.2948964908710636;
    lbModel.oneByTheta0 =  1.0/lbModel.theta0;
    lbModel.wSC1        = 0.04743040745116578;
    lbModel.wSC2        = 0.0016568766450157603 ;
    lbModel.wFCC        = 0.006511753278324646;
    lbModel.wBCC        = 0.04917980624482672;
    lbModel.wBCC1       = 0.004540878011544409 ;
    lbModel.w0          = 1.0 - (6.0*(lbModel.wSC1+lbModel.wSC2) + 12.0*lbModel.wFCC + 8.0*lbModel.wBCC + 8.0*lbModel.wBCC1);
    lbModel.yMinus3_CENTER               =  -3.0;
    lbModel.yTenMinusYSqrPlus15_CENTER   =  15.0;
    lbModel.yMinus3_SC2                  =  (4.0*lbModel.oneByTheta0)-3.0;
    lbModel.yTenMinusYSqrPlus15_SC2      =  (16.0*lbModel.oneByTheta0*lbModel.oneByTheta0)-10.0*(4.0*lbModel.oneByTheta0) + 15.0;
    lbModel.yMinus3_SC1                  =  (1.0*lbModel.oneByTheta0)-3.0;
    lbModel.yTenMinusYSqrPlus15_SC1      =  (lbModel.oneByTheta0*lbModel.oneByTheta0)-10.0*(1.0*lbModel.oneByTheta0) + 15.0;
    lbModel.yMinus3_FCC                  =  (2.0*lbModel.oneByTheta0)-3.0;
    lbModel.yTenMinusYSqrPlus15_FCC      =  (4.0*lbModel.oneByTheta0*lbModel.oneByTheta0)-10.0*(2.0*lbModel.oneByTheta0) + 15.0;
    lbModel.yMinus3_BCC                  =  (0.75*lbModel.oneByTheta0)-3.0;
    lbModel.yTenMinusYSqrPlus15_BCC      =  (0.5625*lbModel.oneByTheta0*lbModel.oneByTheta0)-10.0*(0.75*lbModel.oneByTheta0) + 15.0;
    lbModel.yMinus3_BCC1                 =  (3.0*lbModel.oneByTheta0)-3.0;
    lbModel.yTenMinusYSqrPlus15_BCC1     =  (9.0*lbModel.oneByTheta0*lbModel.oneByTheta0)-10.0*(3.0*lbModel.oneByTheta0) + 15.0;
}


template <typename dataType1>
void getLatticeParameter(lbmRD3Q41<dataType1> &lbModel)
{
    lbModel.wt[lbModel.DV_ZERO_ZERO_ZERO]  = lbModel.w0;     lbModel.cx[lbModel.DV_ZERO_ZERO_ZERO]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_ZERO_ZERO]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_ZERO_ZERO]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_ZERO_ZERO_P1  ]  = lbModel.wSC1;   lbModel.cx[lbModel.DV_ZERO_ZERO_P1  ]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_ZERO_P1  ]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_ZERO_P1  ]  =  1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_ZERO_ZERO_P2  ]  = lbModel.wSC2;   lbModel.cx[lbModel.DV_ZERO_ZERO_P2  ]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_ZERO_P2  ]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_ZERO_P2  ]  =  2.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_ZERO_P2_ZERO  ]  = lbModel.wSC2;   lbModel.cx[lbModel.DV_ZERO_P2_ZERO  ]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_P2_ZERO  ]  =  2.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_P2_ZERO  ]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P2_ZERO_ZERO  ]  = lbModel.wSC2;   lbModel.cx[lbModel.DV_P2_ZERO_ZERO  ]  =  2.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P2_ZERO_ZERO  ]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_P2_ZERO_ZERO  ]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_ZERO_ZERO_M1  ]  = lbModel.wSC1;   lbModel.cx[lbModel.DV_ZERO_ZERO_M1  ]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_ZERO_M1  ]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_ZERO_M1  ]  = -1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_ZERO_ZERO_M2  ]  = lbModel.wSC2;   lbModel.cx[lbModel.DV_ZERO_ZERO_M2  ]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_ZERO_M2  ]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_ZERO_M2  ]  = -2.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_ZERO_M2_ZERO  ]  = lbModel.wSC2;   lbModel.cx[lbModel.DV_ZERO_M2_ZERO  ]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_M2_ZERO  ]  = -2.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_M2_ZERO  ]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M2_ZERO_ZERO  ]  = lbModel.wSC2;   lbModel.cx[lbModel.DV_M2_ZERO_ZERO  ]  = -2.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M2_ZERO_ZERO  ]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_M2_ZERO_ZERO  ]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_ZERO_P1_P1    ]  = lbModel.wFCC;   lbModel.cx[lbModel.DV_ZERO_P1_P1    ]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_P1_P1    ]  =  1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_P1_P1    ]  =  1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_ZERO_M1_P1    ]  = lbModel.wFCC;   lbModel.cx[lbModel.DV_ZERO_M1_P1    ]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_M1_P1    ]  = -1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_M1_P1    ]  =  1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P1_ZERO_P1    ]  = lbModel.wFCC;   lbModel.cx[lbModel.DV_P1_ZERO_P1    ]  =  1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P1_ZERO_P1    ]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_P1_ZERO_P1    ]  =  1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M1_ZERO_P1    ]  = lbModel.wFCC;   lbModel.cx[lbModel.DV_M1_ZERO_P1    ]  = -1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M1_ZERO_P1    ]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_M1_ZERO_P1    ]  =  1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_ZERO_P1_M1    ]  = lbModel.wFCC;   lbModel.cx[lbModel.DV_ZERO_P1_M1    ]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_P1_M1    ]  =  1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_P1_M1    ]  = -1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_ZERO_M1_M1    ]  = lbModel.wFCC;   lbModel.cx[lbModel.DV_ZERO_M1_M1    ]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_M1_M1    ]  = -1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_M1_M1    ]  = -1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P1_ZERO_M1    ]  = lbModel.wFCC;   lbModel.cx[lbModel.DV_P1_ZERO_M1    ]  =  1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P1_ZERO_M1    ]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_P1_ZERO_M1    ]  = -1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M1_ZERO_M1    ]  = lbModel.wFCC;   lbModel.cx[lbModel.DV_M1_ZERO_M1    ]  = -1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M1_ZERO_M1    ]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_M1_ZERO_M1    ]  = -1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P1_P1_ZERO    ]  = lbModel.wFCC;   lbModel.cx[lbModel.DV_P1_P1_ZERO    ]  =  1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P1_P1_ZERO    ]  =  1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_P1_P1_ZERO    ]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M1_P1_ZERO    ]  = lbModel.wFCC;   lbModel.cx[lbModel.DV_M1_P1_ZERO    ]  = -1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M1_P1_ZERO    ]  =  1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_M1_P1_ZERO    ]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_ZERO_P1_ZERO  ]  = lbModel.wSC1;   lbModel.cx[lbModel.DV_ZERO_P1_ZERO  ]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_P1_ZERO  ]  =  1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_P1_ZERO  ]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P1_ZERO_ZERO  ]  = lbModel.wSC1;   lbModel.cx[lbModel.DV_P1_ZERO_ZERO  ]  =  1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P1_ZERO_ZERO  ]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_P1_ZERO_ZERO  ]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P1_M1_ZERO    ]  = lbModel.wFCC;   lbModel.cx[lbModel.DV_P1_M1_ZERO    ]  =  1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P1_M1_ZERO    ]  = -1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_P1_M1_ZERO    ]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M1_M1_ZERO    ]  = lbModel.wFCC;   lbModel.cx[lbModel.DV_M1_M1_ZERO    ]  = -1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M1_M1_ZERO    ]  = -1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_M1_M1_ZERO    ]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_ZERO_M1_ZERO  ]  = lbModel.wSC1;   lbModel.cx[lbModel.DV_ZERO_M1_ZERO  ]  =  0.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_ZERO_M1_ZERO  ]  = -1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_ZERO_M1_ZERO  ]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M1_ZERO_ZERO  ]  = lbModel.wSC1;   lbModel.cx[lbModel.DV_M1_ZERO_ZERO  ]  = -1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M1_ZERO_ZERO  ]  =  0.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_M1_ZERO_ZERO  ]  =  0.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P_P_P         ]  = lbModel.wBCC;   lbModel.cx[lbModel.DV_P_P_P         ]  =  0.5*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P_P_P         ]  =  0.5*lbModel.cellSize;   lbModel.cz[lbModel.DV_P_P_P         ]  =  0.5*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M_P_P         ]  = lbModel.wBCC;   lbModel.cx[lbModel.DV_M_P_P         ]  = -0.5*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M_P_P         ]  =  0.5*lbModel.cellSize;   lbModel.cz[lbModel.DV_M_P_P         ]  =  0.5*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M_M_P         ]  = lbModel.wBCC;   lbModel.cx[lbModel.DV_M_M_P         ]  = -0.5*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M_M_P         ]  = -0.5*lbModel.cellSize;   lbModel.cz[lbModel.DV_M_M_P         ]  =  0.5*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P_M_P         ]  = lbModel.wBCC;   lbModel.cx[lbModel.DV_P_M_P         ]  =  0.5*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P_M_P         ]  = -0.5*lbModel.cellSize;   lbModel.cz[lbModel.DV_P_M_P         ]  =  0.5*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P_P_M         ]  = lbModel.wBCC;   lbModel.cx[lbModel.DV_P_P_M         ]  =  0.5*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P_P_M         ]  =  0.5*lbModel.cellSize;   lbModel.cz[lbModel.DV_P_P_M         ]  = -0.5*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M_P_M         ]  = lbModel.wBCC;   lbModel.cx[lbModel.DV_M_P_M         ]  = -0.5*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M_P_M         ]  =  0.5*lbModel.cellSize;   lbModel.cz[lbModel.DV_M_P_M         ]  = -0.5*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M_M_M         ]  = lbModel.wBCC;   lbModel.cx[lbModel.DV_M_M_M         ]  = -0.5*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M_M_M         ]  = -0.5*lbModel.cellSize;   lbModel.cz[lbModel.DV_M_M_M         ]  = -0.5*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P_M_M         ]  = lbModel.wBCC;   lbModel.cx[lbModel.DV_P_M_M         ]  =  0.5*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P_M_M         ]  = -0.5*lbModel.cellSize;   lbModel.cz[lbModel.DV_P_M_M         ]  = -0.5*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P1_P1_P1      ]  = lbModel.wBCC1;  lbModel.cx[lbModel.DV_P1_P1_P1      ]  =  1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P1_P1_P1      ]  =  1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_P1_P1_P1      ]  =  1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M1_P1_P1      ]  = lbModel.wBCC1;  lbModel.cx[lbModel.DV_M1_P1_P1      ]  = -1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M1_P1_P1      ]  =  1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_M1_P1_P1      ]  =  1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M1_M1_P1      ]  = lbModel.wBCC1;  lbModel.cx[lbModel.DV_M1_M1_P1      ]  = -1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M1_M1_P1      ]  = -1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_M1_M1_P1      ]  =  1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P1_M1_P1      ]  = lbModel.wBCC1;  lbModel.cx[lbModel.DV_P1_M1_P1      ]  =  1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P1_M1_P1      ]  = -1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_P1_M1_P1      ]  =  1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P1_P1_M1      ]  = lbModel.wBCC1;  lbModel.cx[lbModel.DV_P1_P1_M1      ]  =  1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P1_P1_M1      ]  =  1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_P1_P1_M1      ]  = -1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M1_P1_M1      ]  = lbModel.wBCC1;  lbModel.cx[lbModel.DV_M1_P1_M1      ]  = -1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M1_P1_M1      ]  =  1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_M1_P1_M1      ]  = -1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_M1_M1_M1      ]  = lbModel.wBCC1;  lbModel.cx[lbModel.DV_M1_M1_M1      ]  = -1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_M1_M1_M1      ]  = -1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_M1_M1_M1      ]  = -1.0*lbModel.cellSize;
    lbModel.wt[lbModel.DV_P1_M1_M1      ]  = lbModel.wBCC1;  lbModel.cx[lbModel.DV_P1_M1_M1      ]  =  1.0*lbModel.cellSize ;   lbModel.cy[lbModel.DV_P1_M1_M1      ]  = -1.0*lbModel.cellSize;   lbModel.cz[lbModel.DV_P1_M1_M1      ]  = -1.0*lbModel.cellSize;

    for(int dv=0;dv<lbModel.dvN;dv++)
    {
      lbModel.cxcy [dv]  = lbModel.cx[dv]*lbModel.cy[dv]  ;
      lbModel.cycz [dv]  = lbModel.cy[dv]*lbModel.cz[dv]  ;
      lbModel.czcx [dv]  = lbModel.cz[dv]*lbModel.cx[dv]  ;
      lbModel.cx2  [dv]  = lbModel.cx[dv]*lbModel.cx[dv]  ;
      lbModel.cy2  [dv]  = lbModel.cy[dv]*lbModel.cy[dv]  ;
      lbModel.cz2  [dv]  = lbModel.cz[dv]*lbModel.cz[dv]  ;
      lbModel.cc   [dv]  = lbModel.cx2[dv]+lbModel.cy2[dv]+lbModel.cz2[dv];
      lbModel.cxcsq[dv]  = lbModel.cc[dv]*lbModel.cx[dv];
      lbModel.cycsq[dv]  = lbModel.cc[dv]*lbModel.cy[dv];
      lbModel.czcsq[dv]  = lbModel.cc[dv]*lbModel.cz[dv];
      lbModel.csq2 [dv]  = lbModel.cc   [dv];
    }



}



template <int N,int numblock, typename dataType1>
void copyFromNode(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3)
{
    for(int i=0;i<VECT_LENGTH;i++)
        lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO] = myGrid(i1+i,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp1[i][dv]   =    myGrid.value(index+dv) ;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp2[i][dv]  =    myGrid.value(index+dv) ;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp3[i][dv]   =    myGrid.value(index+dv) ;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp4[i][dv]   =    myGrid.value(index+dv) ;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp5[i][dv]   =    myGrid.value(index+dv) ;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp6[i][dv]   =    myGrid.value(index+dv) ;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp7[i][dv]   =    myGrid.value(index+dv) ;
    }
    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp8[i][dv]   =    myGrid.value(index+dv) ;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp9[i][dv]   =    myGrid.value(index+dv) ;
    }
    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp10[i][dv]   =    myGrid.value(index+dv) ;
    }
}


template <int N,int numblock, typename dataType1>
void copyToNode(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3)
{
    for(int i=0;i<VECT_LENGTH;i++)
        myGrid(i1+i,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO)=   lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO] ;

    for(int i=0;i<VECT_LENGTH;i++)
        for(int dv=0;dv<4;dv++)
            myGrid(i1+i,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],dv)     =    lbModel.fTemp1[i][dv] ;

        for(int i=0;i<VECT_LENGTH;i++)
            for(int dv=0;dv<4;dv++)
                myGrid(i1+i,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],dv)     =    lbModel.fTemp2[i][dv] ;

            for(int i=0;i<VECT_LENGTH;i++)
                for(int dv=0;dv<4;dv++)
                    myGrid(i1+i,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],dv)     =    lbModel.fTemp3[i][dv] ;

                for(int i=0;i<VECT_LENGTH;i++)
                    for(int dv=0;dv<4;dv++)
                        myGrid(i1+i,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],dv)     =    lbModel.fTemp4[i][dv] ;

                    for(int i=0;i<VECT_LENGTH;i++)
                        for(int dv=0;dv<4;dv++)
                            myGrid(i1+i,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],dv)     =    lbModel.fTemp5[i][dv] ;

                        for(int i=0;i<VECT_LENGTH;i++)
                            for(int dv=0;dv<4;dv++)
                                myGrid(i1+i,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],dv)     =    lbModel.fTemp6[i][dv] ;

                            for(int i=0;i<VECT_LENGTH;i++)
                                for(int dv=0;dv<4;dv++)
                                    myGrid(i1+i,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],dv)     =    lbModel.fTemp7[i][dv] ;


                                for(int i=0;i<VECT_LENGTH;i++)
                                    for(int dv=0;dv<4;dv++)
                                        myGrid(i1+i,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],dv)     =    lbModel.fTemp8[i][dv] ;

                                    for(int i=0;i<VECT_LENGTH;i++)
                                        for(int dv=0;dv<4;dv++)
                                            myGrid(i1+i,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],dv)     =    lbModel.fTemp9[i][dv] ;

                                        for(int i=0;i<VECT_LENGTH;i++)
                                            for(int dv=0;dv<4;dv++)
                                                myGrid(i1+i,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],dv)   =    lbModel.fTemp10[i][dv];
}

template <int N,int numblock, typename dataType1>
inline void copyToCell(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3)
{
    for(int i=0;i<VECT_LENGTH;i++)
        myGrid(i1+i,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO)=   lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO] ;

    for(int i=0;i<VECT_LENGTH;i++)
        for(int dv=0;dv<4;dv++)
            myGrid(i1+i,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],dv)     =    lbModel.fTemp1[i][dv] ;

        for(int i=0;i<VECT_LENGTH;i++)
            for(int dv=0;dv<4;dv++)
                myGrid(i1+i,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],dv)     =    lbModel.fTemp2[i][dv] ;

            for(int i=0;i<VECT_LENGTH;i++)
                for(int dv=0;dv<4;dv++)
                    myGrid(i1+i,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],dv)     =    lbModel.fTemp3[i][dv] ;

                for(int i=0;i<VECT_LENGTH;i++)
                    for(int dv=0;dv<4;dv++)
                        myGrid(i1+i,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],dv)     =    lbModel.fTemp4[i][dv] ;

                    for(int i=0;i<VECT_LENGTH;i++)
                        for(int dv=0;dv<4;dv++)
                            myGrid(i1+i,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],dv)     =    lbModel.fTemp5[i][dv] ;

                        for(int i=0;i<VECT_LENGTH;i++)
                            for(int dv=0;dv<4;dv++)
                                myGrid(i1+i,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],dv)     =    lbModel.fTemp6[i][dv] ;

                            for(int i=0;i<VECT_LENGTH;i++)
                                for(int dv=0;dv<4;dv++)
                                    myGrid(i1+i,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],dv)     =    lbModel.fTemp7[i][dv] ;

                                for(int i=0;i<VECT_LENGTH;i++)
                                    for(int dv=0;dv<4;dv++)
                                        myGrid(i1+i,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],dv)     =    lbModel.fTemp8[i][dv] ;

                                    for(int i=0;i<VECT_LENGTH;i++)
                                        for(int dv=0;dv<4;dv++)
                                            myGrid(i1+i,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],dv)     =    lbModel.fTemp9[i][dv] ;

                                        for(int i=0;i<VECT_LENGTH;i++)
                                            for(int dv=0;dv<4;dv++)
                                                myGrid(i1+i,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],dv)   =    lbModel.fTemp10[i][dv];
}

template <int N,int numblock, typename dataType1>
inline void copyFromCell(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3)
{
    for(int i=0;i<VECT_LENGTH;i++)
    {
        lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO] = myGrid(i1+i,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    }
    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp1[i][dv]   =    myGrid.value(index+dv) ;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp2[i][dv]   =    myGrid.value(index+dv) ;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp3[i][dv]   =    myGrid.value(index+dv) ;
    }



    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp4[i][dv]   =    myGrid.value(index+dv) ;
    }


    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp5[i][dv]   =    myGrid.value(index+dv) ;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp6[i][dv]   =    myGrid.value(index+dv) ;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp7[i][dv]   =    myGrid.value(index+dv) ;
    }
    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp8[i][dv]   =    myGrid.value(index+dv) ;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp9[i][dv]   =    myGrid.value(index+dv) ;
    }
    for(int i=0;i<VECT_LENGTH;i++)
    {
        unsigned long long int index = myGrid.getIndex(i1+i,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp10[i][dv]   =    myGrid.value(index+dv) ;
    }
}


    template <int N,int numblock, typename dataType1>
    void copyToNodeSinglePoint(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3,int point)
    {
       myGrid(i1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) = lbModel.fTemp0[point][lbModel.CENTER_DV_ZERO_ZERO_ZERO];

        int index = myGrid.getIndex(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],  0);
        for(int dv=0;dv<4;dv++)
           myGrid.value(index+dv)  =  lbModel.fTemp1[point][dv];

        index = myGrid.getIndex(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],  0);
        for(int dv=0;dv<4;dv++)
           myGrid.value(index+dv) =  lbModel.fTemp2[point][dv] ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],  0);
        for(int dv=0;dv<4;dv++)
           myGrid.value(index+dv) = lbModel.fTemp3[point][dv];

        index = myGrid.getIndex(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],  0);
        for(int dv=0;dv<4;dv++)
            myGrid.value(index+dv) = lbModel.fTemp4[point][dv] ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],  0);
        for(int dv=0;dv<4;dv++)
           myGrid.value(index+dv) = lbModel.fTemp5[point][dv] ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],  0);
        for(int dv=0;dv<4;dv++)
            myGrid.value(index+dv) = lbModel.fTemp6[point][dv]   ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],  0);
        for(int dv=0;dv<4;dv++)
            myGrid.value(index+dv) = lbModel.fTemp7[point][dv]    ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],  0);
        for(int dv=0;dv<4;dv++)
            myGrid.value(index+dv) = lbModel.fTemp8[point][dv]     ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],  0);
        for(int dv=0;dv<4;dv++)
          myGrid.value(index+dv) = lbModel.fTemp9[point][dv]    ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],  0);
        for(int dv=0;dv<4;dv++)
           myGrid.value(index+dv) = lbModel.fTemp10[point][dv]  ;

    }

    template <int N,int numblock, typename dataType1>
    void copyToCellSinglePoint(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3,int point)
    {
       myGrid(i1,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) = lbModel.fTemp0[point][lbModel.CENTER_DV_ZERO_ZERO_ZERO];

        int index = myGrid.getIndex(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],  0);
        for(int dv=0;dv<4;dv++)
           myGrid.value(index+dv)  =  lbModel.fTemp1[point][dv];

        index = myGrid.getIndex(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],  0);
        for(int dv=0;dv<4;dv++)
           myGrid.value(index+dv) =  lbModel.fTemp2[point][dv] ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],  0);
        for(int dv=0;dv<4;dv++)
           myGrid.value(index+dv) = lbModel.fTemp3[point][dv];

        index = myGrid.getIndex(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],  0);
        for(int dv=0;dv<4;dv++)
            myGrid.value(index+dv) = lbModel.fTemp4[point][dv] ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],  0);
        for(int dv=0;dv<4;dv++)
           myGrid.value(index+dv) = lbModel.fTemp5[point][dv] ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],  0);
        for(int dv=0;dv<4;dv++)
            myGrid.value(index+dv) = lbModel.fTemp6[point][dv]   ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],  0);
        for(int dv=0;dv<4;dv++)
            myGrid.value(index+dv) = lbModel.fTemp7[point][dv]    ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],  0);
        for(int dv=0;dv<4;dv++)
            myGrid.value(index+dv) = lbModel.fTemp8[point][dv]     ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],  0);
        for(int dv=0;dv<4;dv++)
          myGrid.value(index+dv) = lbModel.fTemp9[point][dv]    ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],  0);
        for(int dv=0;dv<4;dv++)
           myGrid.value(index+dv) = lbModel.fTemp10[point][dv]  ;
    }


    template <int N,int numblock, typename dataType1>
    void copyFromNodeSinglePoint(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3,int point)
    {
       lbModel.fTemp0[point][lbModel.CENTER_DV_ZERO_ZERO_ZERO] = myGrid(i1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) ;

        int index = myGrid.getIndex(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp1[point][dv] = myGrid.value(index+dv) ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],  0);
        for(int dv=0;dv<4;dv++)
          lbModel.fTemp2[point][dv] = myGrid.value(index+dv) ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],  0);
        for(int dv=0;dv<4;dv++)
           lbModel.fTemp3[point][dv] = myGrid.value(index+dv) ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],  0);
        for(int dv=0;dv<4;dv++)
           lbModel.fTemp4[point][dv] = myGrid.value(index+dv)  ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],  0);
        for(int dv=0;dv<4;dv++)
           lbModel.fTemp5[point][dv] = myGrid.value(index+dv)  ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp6[point][dv] = myGrid.value(index+dv)  ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp7[point][dv] = myGrid.value(index+dv)   ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],  0);
        for(int dv=0;dv<4;dv++)
           lbModel.fTemp8[point][dv] =  myGrid.value(index+dv)   ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],  0);
        for(int dv=0;dv<4;dv++)
          lbModel.fTemp9[point][dv]   = myGrid.value(index+dv)   ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],  0);
        for(int dv=0;dv<4;dv++)
           lbModel.fTemp10[point][dv]  = myGrid.value(index+dv)  ;

    }

    template <int N,int numblock, typename dataType1>
    void copyFromCellSinglePoint(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int i1,int i2,int i3,int point)
    {
       lbModel.fTemp0[point][lbModel.CENTER_DV_ZERO_ZERO_ZERO] = myGrid(i1,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO) ;

        int index = myGrid.getIndex(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp1[point][dv] = myGrid.value(index+dv) ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],  0);
        for(int dv=0;dv<4;dv++)
          lbModel.fTemp2[point][dv] = myGrid.value(index+dv) ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],  0);
        for(int dv=0;dv<4;dv++)
           lbModel.fTemp3[point][dv] = myGrid.value(index+dv) ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],  0);
        for(int dv=0;dv<4;dv++)
           lbModel.fTemp4[point][dv] = myGrid.value(index+dv)  ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],  0);
        for(int dv=0;dv<4;dv++)
           lbModel.fTemp5[point][dv] = myGrid.value(index+dv)  ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp6[point][dv] = myGrid.value(index+dv)  ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],  0);
        for(int dv=0;dv<4;dv++)
            lbModel.fTemp7[point][dv] = myGrid.value(index+dv)   ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],  0);
        for(int dv=0;dv<4;dv++)
           lbModel.fTemp8[point][dv] =  myGrid.value(index+dv)   ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],  0);
        for(int dv=0;dv<4;dv++)
          lbModel.fTemp9[point][dv]   = myGrid.value(index+dv)   ;

        index = myGrid.getIndex(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],  0);
        for(int dv=0;dv<4;dv++)
           lbModel.fTemp10[point][dv]  = myGrid.value(index+dv)  ;

    }



template<typename dataType1>
void getHydroMoment(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta)
{
    dataType1 sum(0.0);

    for(int i=0;i<VECT_LENGTH;i++)
    {
        rho[i] = lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO];
    }

    // G1 and G2
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp1[i][lbModel.G1_DV_P2_ZERO_ZERO]+lbModel.fTemp2[i][lbModel.G2_DV_M2_ZERO_ZERO];
        rho[i]    += sum;
        theta[i]  = 4.0*sum;
        uX[i]     = 2.0*(lbModel.fTemp1[i][lbModel.G1_DV_P2_ZERO_ZERO]-lbModel.fTemp2[i][lbModel.G2_DV_M2_ZERO_ZERO]);

        sum       = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_P2_ZERO]+lbModel.fTemp2[i][lbModel.G2_DV_ZERO_M2_ZERO];
        rho[i]   += sum;
        theta[i] += 4.0*sum;
        uY[i]     = 2.0*(lbModel.fTemp1[i][lbModel.G1_DV_ZERO_P2_ZERO]-lbModel.fTemp2[i][lbModel.G2_DV_ZERO_M2_ZERO]);

        sum       = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P2]+lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M2];
        rho[i]   += sum;
        theta[i] += 4.0*sum;
        uZ[i]     = 2.0*(lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P2]-lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M2]);

        sum       = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P1]+lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M1];
        rho[i]   += sum;
        theta[i] += sum;
        uZ[i]    += lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P1]-lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M1];
    }

    // G3 and G4
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]+lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        uY[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1]);

        sum       = lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1 ]+lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        uY[i]    -= (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1]-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1]-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1]);

        sum       = lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1]+lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        uX[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1]-lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1]-lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1]);

        sum       = lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1]+lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        uX[i]    -= (lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1]-lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1]-lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1]);
    }

    // G5 and G6
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO]+lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        uX[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO]-lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO]);
        uY[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO]-lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO]);

        sum       = lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO]+lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        uX[i]    += (lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO]-lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO]);
        uY[i]    -= (lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO]-lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO]);

        sum      = lbModel.fTemp5[i][lbModel.G5_DV_ZERO_P1_ZERO]+lbModel.fTemp6[i][lbModel.G6_DV_ZERO_M1_ZERO];
        rho[i]   += sum;
        theta[i] += sum;
        uY[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_ZERO_P1_ZERO]-lbModel.fTemp6[i][lbModel.G6_DV_ZERO_M1_ZERO]);

        sum       = lbModel.fTemp5[i][lbModel.G5_DV_P1_ZERO_ZERO]+lbModel.fTemp6[i][lbModel.G6_DV_M1_ZERO_ZERO];
        rho[i]   += sum;
        theta[i] += sum;
        uX[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_P1_ZERO_ZERO]-lbModel.fTemp6[i][lbModel.G6_DV_M1_ZERO_ZERO]);
    }

    // G7 and G8
    const dataType1 threeByFour=3.0/4.0;
    for(int i=0;i< VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]+lbModel.fTemp8[i][lbModel.G8_DV_M_M_M];
        rho[i]   += sum;
        theta[i] += threeByFour*sum ;
        uX[i]    += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]);
        uY[i]    += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]);
        uZ[i]    += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]);

        sum       = lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]+lbModel.fTemp8[i][lbModel.G8_DV_P_M_M];
        rho[i]   += sum;
        theta[i] += threeByFour*sum ;
        uX[i]    -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]);
        uY[i]    += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]);
        uZ[i]    += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]);


        sum       = lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]+lbModel.fTemp8[i][lbModel.G8_DV_P_P_M];
        rho[i]   += sum;
        theta[i] += threeByFour*sum ;
        uX[i]    -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]);
        uY[i]    -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]);
        uZ[i]    += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]);

        sum       = lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]+lbModel.fTemp8[i][lbModel.G8_DV_M_P_M];
        rho[i]   += sum;
        theta[i] += threeByFour*sum ;
        uX[i]    += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_P_M ]);
        uY[i]    -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_P_M ]);
        uZ[i]    += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_P_M ]);
    }

    // G9 and G10
    for(int i=0;i< VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]+lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        uX[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]);
        uY[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]);

        sum       = lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]+lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        uX[i]    -= (lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]);
        uY[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]);


        sum       = lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]+lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        uX[i]    -= ( lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]);
        uY[i]    -= ( lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]);
        uZ[i]    += ( lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]);

        sum       = lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1 ]+lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        uX[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]);
        uY[i]    -= (lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]);
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        dataType1  oneByRho = 1.0/rho[i];
        dataType1  oneByThree = 1.0/3.0;
        uX[i]  *= oneByRho;
        uY[i]  *= oneByRho;
        uZ[i]  *= oneByRho;
        theta[i] = oneByThree * (theta[i] - rho[i]*(uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i]));
        theta[i] *= oneByRho;
    }
}


template <int N,int numblock, typename dataType1>
void getHydroMomentsFromNode(lbmRD3Q41<dataType1> &lbModel,gridBCC3D<N, numblock, dataType1> &myGrid,int i1, int i2, int i3, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta)
{
    //         dataType1 sum(0.0);
    SIMD_REG _RHO, _UX, _UY, _UZ, _THETA;              //4
    SIMD_REG _temp1,_temp2,_temp3,_temp4;              //4
    SIMD_REG _temp5,_temp6,_temp7,_temp8;              //4
    SIMD_REG _temp9,_twoReg,_fourReg;                  //3   = 16 Registers in Total
    dataType1 tempArray[4] __attribute__ ((aligned(32)));
    dataType1  *pointer1;

    //////////////
    //  CENTER  //
    //////////////
    tempArray[0]=myGrid(i1  ,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    tempArray[1]=myGrid(i1+1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    tempArray[2]=myGrid(i1+2,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    tempArray[3]=myGrid(i1+3,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    pointer1 = &tempArray[0];
    _RHO = LOAD_PD(pointer1);


    ///////////////////
    //  SC1 and SC2  //
    ///////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1);
    _temp1 = LOAD_PD(pointer1);      // All SC1 at  i
    _temp2 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains ZZP1 from i,i+1,i+2,i+3
    // temp2 now contains ZZP2 from i,i+1,i+2,i+3
    // temp3 now contains ZP2Z from i,i+1,i+2,i+3
    // temp4 now contains P2ZZ from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1);
    _UZ    = LOAD_PD(pointer1);       // All SC2 at  i
    _temp6 = LOAD_PD(pointer1+4);     //     "      i+1
    _UY    = LOAD_PD(pointer1+8);     //     "      i+2
    _UX    = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_UZ,_temp6,_UY,_UX);
    // _UZ   now contains ZZM1 from i,i+1,i+2,i+3
    // temp6 now contains ZZM2 from i,i+1,i+2,i+3
    // _UY   now contains ZM2Z from i,i+1,i+2,i+3
    // _UX   now contains M2ZZ from i,i+1,i+2,i+3

    _THETA = ADD_PD(_temp1,_UZ)   ;  _UZ    = SUB_PD(_temp1,_UZ);
    _temp1 = ADD_PD(_temp2,_temp6);  _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_UY)   ;  _UY    = SUB_PD(_temp3,_UY);
    _temp3 = ADD_PD(_temp4,_UX)   ;  _UX    = SUB_PD(_temp4,_UX);

    _RHO = ADD_PD(_RHO,_THETA); // ZZP1 + ZZM1
    _RHO = ADD_PD(_RHO,_temp1); // ZZP1 + ZZM1 + ZZP2 + ZZM2
    _RHO = ADD_PD(_RHO,_temp2); // ZZP1 + ZZM1 + ZZP2 + ZZM2 + ZP2Z + ZM2Z
    _RHO = ADD_PD(_RHO,_temp3); // ZZP1 + ZZM1 + ZZP2 + ZZM2 + ZP2Z + ZM2Z + P2ZZ + M2ZZ  at i,i+1,i+2,i+3


    _twoReg = SET1_PD(2.0);
    _temp6  = MUL_PD(_temp6,_twoReg);    // (ZZP2-ZZM2)*2.0
    _UY     = MUL_PD(_UY   ,_twoReg);    // (ZP2Z-ZM2Z)*2.0
    _UX     = MUL_PD(_UX   ,_twoReg);    // (P2ZZ-M2ZZ)*2.0

    _UZ     = ADD_PD(_UZ,_temp6);        // (ZZP1-ZZM1) + (ZZP2-ZZM2)*2.0

    _fourReg = SET1_PD(4.0);
    _temp1   = MUL_PD(_temp1,_fourReg);    // (ZZP2+ZZM2)*4.0
    _temp2   = MUL_PD(_temp2,_fourReg);    // (ZP2Z+ZM2Z)*4.0
    _temp3   = MUL_PD(_temp3,_fourReg);    // (P2ZZ+M2ZZ)*4.0

    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);



    /////////////////////
    //  FCC1 and FCC2  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1);

    _temp1 = LOAD_PD(pointer1);    // All FCC1 at  i
    _temp2 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains ZPP  from i,i+1,i+2,i+3
    // temp2 now contains ZMP  from i,i+1,i+2,i+3
    // temp3 now contains PZP  from i,i+1,i+2,i+3
    // temp4 now contains MZP  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1);

    _temp5 = LOAD_PD(pointer1);      // All FCC2 at  i
    _temp6 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains ZMM from i,i+1,i+2,i+3
    // temp6 now contains ZPM from i,i+1,i+2,i+3
    // temp7 now contains MZM from i,i+1,i+2,i+3
    // temp8 now contains PZM from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9); // ZPP + ZMM
    _RHO = ADD_PD(_RHO,_temp1); // ZPP + ZMM + ZMP + ZPM
    _RHO = ADD_PD(_RHO,_temp2); // ZPP + ZMM + ZMP + ZPM + PZP + MZM
    _RHO = ADD_PD(_RHO,_temp3); // ZPP + ZMM + ZMP + ZPM + PZP + MZM + MZP + PZM at i,i+1,i+2,i+3

    //  temp9 =  ZPP + ZMM    temp5 =  ZPP - ZMM
    //  temp1 =  ZMP + ZPM    temp6 =  ZMP - ZPM
    //  temp2 =  PZP + MZM    temp7 =  PZP - MZM
    //  temp3 =  MZP + PZM    temp8 =  MZP - PZM

    _UZ  = ADD_PD(_UZ,_temp5);
    _UZ  = ADD_PD(_UZ,_temp6);
    _UZ  = ADD_PD(_UZ,_temp7);
    _UZ  = ADD_PD(_UZ,_temp8);
    _UY  = ADD_PD(_UY,_temp5);
    _UY  = SUB_PD(_UY,_temp6);
    _UX  = ADD_PD(_UX,_temp7);
    _UX  = SUB_PD(_UX,_temp8); //at i,i+1,i+2,i+3



    _temp9 = MUL_PD(_temp9,_twoReg);
    _temp1 = MUL_PD(_temp1,_twoReg);
    _temp2 = MUL_PD(_temp2,_twoReg);
    _temp3 = MUL_PD(_temp3,_twoReg);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);



    /////////////////////
    //  FCC3 and FCC4  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO);

    _temp1 = LOAD_PD(pointer1);      // All FCC3 at  i
    _temp2 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains PPZ  from i,i+1,i+2,i+3
    // temp2 now contains MPZ  from i,i+1,i+2,i+3
    // temp3 now contains ZPZ  from i,i+1,i+2,i+3
    // temp4 now contains PZZ  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO);

    _temp5 = LOAD_PD(pointer1);      // All FCC4 at  i
    _temp6 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains MMZ from i,i+1,i+2,i+3
    // temp6 now contains PMZ from i,i+1,i+2,i+3
    // temp7 now contains ZMZ from i,i+1,i+2,i+3
    // temp8 now contains MZZ from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9);
    _RHO = ADD_PD(_RHO,_temp1);
    _RHO = ADD_PD(_RHO,_temp2);
    _RHO = ADD_PD(_RHO,_temp3); //at i,i+1,i+2,i+3

    //  temp1 =  PPZ + MMZ    temp5 =  PPZ - MMZ
    //  temp2 =  MPZ + PMZ    temp6 =  MPZ - PMZ
    //  temp3 =  ZPZ + ZMZ    temp7 =  ZPZ - ZMZ
    //  temp4 =  PZZ + MZZ    temp8 =  PZZ - MZZ

    _UY  = ADD_PD(_UY,_temp5);
    _UY  = ADD_PD(_UY,_temp6);
    _UY  = ADD_PD(_UY,_temp7);
    _UX  = ADD_PD(_UX,_temp5);
    _UX  = SUB_PD(_UX,_temp6);
    _UX  = ADD_PD(_UX,_temp8); //at i,i+1,i+2,i+3


    _temp9 = MUL_PD(_temp9,_twoReg);
    _temp1 = MUL_PD(_temp1,_twoReg);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);


    /////////////////////
    // BCCP and BCCM  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);

    _temp1 = LOAD_PD(pointer1);       // All BCCP at  i
    _temp2 = LOAD_PD(pointer1+4);     //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);     //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains PPP  from i,i+1,i+2,i+3
    // temp2 now contains MPP  from i,i+1,i+2,i+3
    // temp3 now contains MMP  from i,i+1,i+2,i+3
    // temp4 now contains PMP  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);

    _temp5 = LOAD_PD(pointer1);       // All BCCM at  i
    _temp6 = LOAD_PD(pointer1+4);     //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);     //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains MMM from i,i+1,i+2,i+3
    // temp6 now contains PMM from i,i+1,i+2,i+3
    // temp7 now contains PPM from i,i+1,i+2,i+3
    // temp8 now contains MPM from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9);
    _RHO = ADD_PD(_RHO,_temp1);
    _RHO = ADD_PD(_RHO,_temp2);
    _RHO = ADD_PD(_RHO,_temp3); //at i,i+1,i+2,i+3

    //  temp1 =  PPP + MMM    temp5 =  PPP - MMM
    //  temp2 =  MPP + PMM    temp6 =  MPP - PMM
    //  temp3 =  MMP + PPM    temp7 =  MMP - PPM
    //  temp4 =  PMP + MPM    temp8 =  PMP - MPM

    _temp4 = SET1_PD(0.5);
    _temp5 = MUL_PD(_temp5,_temp4);    // (PPP - MMM)*0.5
    _temp6 = MUL_PD(_temp6,_temp4);    // (MPP - PMM)*0.5
    _temp7 = MUL_PD(_temp7,_temp4);    // (MMP - PPM)*0.5
    _temp8 = MUL_PD(_temp8,_temp4);    // (PMP - MPM)*0.5

    _UZ  = ADD_PD(_UZ,_temp5);
    _UZ  = ADD_PD(_UZ,_temp6);
    _UZ  = ADD_PD(_UZ,_temp7);
    _UZ  = ADD_PD(_UZ,_temp8);
    _UY  = ADD_PD(_UY,_temp5);
    _UY  = ADD_PD(_UY,_temp6);
    _UY  = SUB_PD(_UY,_temp7);
    _UY  = SUB_PD(_UY,_temp8);
    _UX  = ADD_PD(_UX,_temp5);
    _UX  = SUB_PD(_UX,_temp6);
    _UX  = SUB_PD(_UX,_temp7);
    _UX  = ADD_PD(_UX,_temp8);  //at i,i+1,i+2,i+3

    _temp4 = SET1_PD(0.75);
    _temp9 = MUL_PD(_temp9,_temp4);
    _temp1 = MUL_PD(_temp1,_temp4);
    _temp2 = MUL_PD(_temp2,_temp4);
    _temp3 = MUL_PD(_temp3,_temp4);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);


    /////////////////////
    // BCCP1 and BCCM1  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_P1_P1);

    _temp1 = LOAD_PD(pointer1);    // All BCCP1 at  i
    _temp2 = LOAD_PD(pointer1+4);      //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);      //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);     //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains PPP  from i,i+1,i+2,i+3
    // temp2 now contains MPP  from i,i+1,i+2,i+3
    // temp3 now contains MMP  from i,i+1,i+2,i+3
    // temp4 now contains PMP  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1);

    _temp5 = LOAD_PD(pointer1);    // All BCCM1 at  i
    _temp6 = LOAD_PD(pointer1+4);      //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);      //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);     //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains MMM from i,i+1,i+2,i+3
    // temp6 now contains PMM from i,i+1,i+2,i+3
    // temp7 now contains PPM from i,i+1,i+2,i+3
    // temp8 now contains MPM from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9);
    _RHO = ADD_PD(_RHO,_temp1);
    _RHO = ADD_PD(_RHO,_temp2);
    _RHO = ADD_PD(_RHO,_temp3); //at i,i+1,i+2,i+3

    //  temp9 =  PPP1 + MMM1    temp5 =  PPP1 - MMM1
    //  temp1 =  MPP1 + PMM1    temp6 =  MPP1 - PMM1
    //  temp2 =  MMP1 + PPM1    temp7 =  MMP1 - PPM1
    //  temp3 =  PMP1 + MPM1    temp8 =  PMP1 - MPM1

    _UZ  = ADD_PD(_UZ,_temp5);
    _UZ  = ADD_PD(_UZ,_temp6);
    _UZ  = ADD_PD(_UZ,_temp7);
    _UZ  = ADD_PD(_UZ,_temp8);
    _UY  = ADD_PD(_UY,_temp5);
    _UY  = ADD_PD(_UY,_temp6);
    _UY  = SUB_PD(_UY,_temp7);
    _UY  = SUB_PD(_UY,_temp8);
    _UX  = ADD_PD(_UX,_temp5);
    _UX  = SUB_PD(_UX,_temp6);
    _UX  = SUB_PD(_UX,_temp7);
    _UX  = ADD_PD(_UX,_temp8);  //at i,i+1,i+2,i+3

    _UX = DIV_PD(_UX,_RHO);
    _UY = DIV_PD(_UY,_RHO);
    _UZ = DIV_PD(_UZ,_RHO);

    _temp4 = SET1_PD(3.0);
    _temp9 = MUL_PD(_temp9,_temp4);
    _temp1 = MUL_PD(_temp1,_temp4);
    _temp2 = MUL_PD(_temp2,_temp4);
    _temp3 = MUL_PD(_temp3,_temp4);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);



    _temp1  = MUL_PD(_UX,_UX);             //temp9  = uX*uX
    _temp2  = MUL_PD(_UY,_UY);             //temp10 = uY*uY
    _temp3  = ADD_PD(_temp1,_temp2);       //temp9  = uX*uX + uY*uY
    _temp4  = MUL_PD(_UZ,_UZ);             //temp10 = uZ*uZ
    _temp1  = ADD_PD(_temp3,_temp4);       //temp9  = uX*uX + uY*uY + uZ*uZ
    _temp2  = MUL_PD(_RHO,_temp1);         //temp9  = rho*(uX*uX + uY*uY + uZ*uZ)
    _THETA  = SUB_PD(_THETA,_temp2);       //theta  = theta - rho*(uX*uX + uY*uY + uZ*uZ)
    _temp1  = SET1_PD(1.0/3.0);
    _THETA  = MUL_PD(_THETA,_temp1);  //theta  = 1.0/3.0 *(theta - rho*(uX*uX + uY*uY + uZ*uZ))
    _THETA = DIV_PD(_THETA,_RHO);         // oneByRho


    STORE_PD(rho  ,_RHO  );
    STORE_PD(uX   ,_UX   );
    STORE_PD(uY   ,_UY   );
    STORE_PD(uZ   ,_UZ   );
    STORE_PD(theta,_THETA);
}


template <int N,int numblock, typename dataType1>
void getHydroMomentsFromCell(lbmRD3Q41<dataType1> &lbModel,gridBCC3D<N, numblock, dataType1> &myGrid,int i1, int i2, int i3, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta)
{
    //         dataType1 sum(0.0);
    SIMD_REG _RHO, _UX, _UY, _UZ, _THETA;              //4
    SIMD_REG _temp1,_temp2,_temp3,_temp4;              //4
    SIMD_REG _temp5,_temp6,_temp7,_temp8;              //4
    SIMD_REG _temp9,_twoReg,_fourReg;                  //3   = 16 Registers in Total
    dataType1 tempArray[4] __attribute__ ((aligned(32)));
    dataType1  *pointer1;

    //////////////
    //  CENTER  //
    //////////////
    tempArray[0]=myGrid(i1  ,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    tempArray[1]=myGrid(i1+1,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    tempArray[2]=myGrid(i1+2,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    tempArray[3]=myGrid(i1+3,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    pointer1 = &tempArray[0];
    _RHO = LOAD_PD(pointer1);


    ///////////////////
    //  SC1 and SC2  //
    ///////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1);
    _temp1 = LOAD_PD(pointer1);      // All SC1 at  i
    _temp2 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains ZZP1 from i,i+1,i+2,i+3
    // temp2 now contains ZZP2 from i,i+1,i+2,i+3
    // temp3 now contains ZP2Z from i,i+1,i+2,i+3
    // temp4 now contains P2ZZ from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1);

    _UZ    = LOAD_PD(pointer1);       // All SC2 at  i
    _temp6 = LOAD_PD(pointer1+4);     //     "      i+1
    _UY    = LOAD_PD(pointer1+8);     //     "      i+2
    _UX    = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_UZ,_temp6,_UY,_UX);
    // temp5 now contains ZZM1 from i,i+1,i+2,i+3
    // temp6 now contains ZZM2 from i,i+1,i+2,i+3
    // temp7 now contains ZM2Z from i,i+1,i+2,i+3
    // temp8 now contains M2ZZ from i,i+1,i+2,i+3

    _THETA = ADD_PD(_temp1,_UZ)   ;  _UZ    = SUB_PD(_temp1,_UZ);
    _temp1 = ADD_PD(_temp2,_temp6);  _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_UY)   ;  _UY    = SUB_PD(_temp3,_UY);
    _temp3 = ADD_PD(_temp4,_UX)   ;  _UX    = SUB_PD(_temp4,_UX);

    _RHO = ADD_PD(_RHO,_THETA); // ZZP1 + ZZM1
    _RHO = ADD_PD(_RHO,_temp1); // ZZP1 + ZZM1 + ZZP2 + ZZM2
    _RHO = ADD_PD(_RHO,_temp2); // ZZP1 + ZZM1 + ZZP2 + ZZM2 + ZP2Z + ZM2Z
    _RHO = ADD_PD(_RHO,_temp3); // ZZP1 + ZZM1 + ZZP2 + ZZM2 + ZP2Z + ZM2Z + P2ZZ + M2ZZ  at i,i+1,i+2,i+3


    _twoReg = SET1_PD(2.0);
    _temp6  = MUL_PD(_temp6,_twoReg);    // (ZZP2-ZZM2)*2.0
    _UY     = MUL_PD(_UY   ,_twoReg);    // (ZP2Z-ZM2Z)*2.0
    _UX     = MUL_PD(_UX   ,_twoReg);    // (P2ZZ-M2ZZ)*2.0

    _UZ     = ADD_PD(_UZ,_temp6);        // (ZZP1-ZZM1) + (ZZP2-ZZM2)*2.0

    _fourReg = SET1_PD(4.0);
    _temp1   = MUL_PD(_temp1,_fourReg);    // (ZZP2+ZZM2)*4.0
    _temp2   = MUL_PD(_temp2,_fourReg);    // (ZP2Z+ZM2Z)*4.0
    _temp3   = MUL_PD(_temp3,_fourReg);    // (P2ZZ+M2ZZ)*4.0

    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);



    /////////////////////
    //  FCC1 and FCC2  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1);

    _temp1 = LOAD_PD(pointer1);    // All FCC1 at  i
    _temp2 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains ZPP  from i,i+1,i+2,i+3
    // temp2 now contains ZMP  from i,i+1,i+2,i+3
    // temp3 now contains PZP  from i,i+1,i+2,i+3
    // temp4 now contains MZP  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1);

    _temp5 = LOAD_PD(pointer1);      // All FCC2 at  i
    _temp6 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains ZMM from i,i+1,i+2,i+3
    // temp6 now contains ZPM from i,i+1,i+2,i+3
    // temp7 now contains MZM from i,i+1,i+2,i+3
    // temp8 now contains PZM from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9); // ZPP + ZMM
    _RHO = ADD_PD(_RHO,_temp1); // ZPP + ZMM + ZMP + ZPM
    _RHO = ADD_PD(_RHO,_temp2); // ZPP + ZMM + ZMP + ZPM + PZP + MZM
    _RHO = ADD_PD(_RHO,_temp3); // ZPP + ZMM + ZMP + ZPM + PZP + MZM + MZP + PZM at i,i+1,i+2,i+3

    //  temp9 =  ZPP + ZMM    temp5 =  ZPP - ZMM
    //  temp1 =  ZMP + ZPM    temp6 =  ZMP - ZPM
    //  temp2 =  PZP + MZM    temp7 =  PZP - MZM
    //  temp3 =  MZP + PZM    temp8 =  MZP - PZM

    _UZ  = ADD_PD(_UZ,_temp5);
    _UZ  = ADD_PD(_UZ,_temp6);
    _UZ  = ADD_PD(_UZ,_temp7);
    _UZ  = ADD_PD(_UZ,_temp8);
    _UY  = ADD_PD(_UY,_temp5);
    _UY  = SUB_PD(_UY,_temp6);
    _UX  = ADD_PD(_UX,_temp7);
    _UX  = SUB_PD(_UX,_temp8); //at i,i+1,i+2,i+3



    _temp9 = MUL_PD(_temp9,_twoReg);
    _temp1 = MUL_PD(_temp1,_twoReg);
    _temp2 = MUL_PD(_temp2,_twoReg);
    _temp3 = MUL_PD(_temp3,_twoReg);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);



    /////////////////////
    //  FCC3 and FCC4  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO);

    _temp1 = LOAD_PD(pointer1);      // All FCC3 at  i
    _temp2 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains PPZ  from i,i+1,i+2,i+3
    // temp2 now contains MPZ  from i,i+1,i+2,i+3
    // temp3 now contains ZPZ  from i,i+1,i+2,i+3
    // temp4 now contains PZZ  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO);

    _temp5 = LOAD_PD(pointer1);      // All FCC4 at  i
    _temp6 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains MMZ from i,i+1,i+2,i+3
    // temp6 now contains PMZ from i,i+1,i+2,i+3
    // temp7 now contains ZMZ from i,i+1,i+2,i+3
    // temp8 now contains MZZ from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9);
    _RHO = ADD_PD(_RHO,_temp1);
    _RHO = ADD_PD(_RHO,_temp2);
    _RHO = ADD_PD(_RHO,_temp3); //at i,i+1,i+2,i+3

    //  temp1 =  PPZ + MMZ    temp5 =  PPZ - MMZ
    //  temp2 =  MPZ + PMZ    temp6 =  MPZ - PMZ
    //  temp3 =  ZPZ + ZMZ    temp7 =  ZPZ - ZMZ
    //  temp4 =  PZZ + MZZ    temp8 =  PZZ - MZZ

    _UY  = ADD_PD(_UY,_temp5);
    _UY  = ADD_PD(_UY,_temp6);
    _UY  = ADD_PD(_UY,_temp7);
    _UX  = ADD_PD(_UX,_temp5);
    _UX  = SUB_PD(_UX,_temp6);
    _UX  = ADD_PD(_UX,_temp8); //at i,i+1,i+2,i+3


    _temp9 = MUL_PD(_temp9,_twoReg);
    _temp1 = MUL_PD(_temp1,_twoReg);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);


    /////////////////////
    // BCCP and BCCM  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);

    _temp1 = LOAD_PD(pointer1);       // All BCCP at  i
    _temp2 = LOAD_PD(pointer1+4);     //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);     //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains PPP  from i,i+1,i+2,i+3
    // temp2 now contains MPP  from i,i+1,i+2,i+3
    // temp3 now contains MMP  from i,i+1,i+2,i+3
    // temp4 now contains PMP  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);

    _temp5 = LOAD_PD(pointer1);       // All BCCM at  i
    _temp6 = LOAD_PD(pointer1+4);     //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);     //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains MMM from i,i+1,i+2,i+3
    // temp6 now contains PMM from i,i+1,i+2,i+3
    // temp7 now contains PPM from i,i+1,i+2,i+3
    // temp8 now contains MPM from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9);
    _RHO = ADD_PD(_RHO,_temp1);
    _RHO = ADD_PD(_RHO,_temp2);
    _RHO = ADD_PD(_RHO,_temp3); //at i,i+1,i+2,i+3

    //  temp1 =  PPP + MMM    temp5 =  PPP - MMM
    //  temp2 =  MPP + PMM    temp6 =  MPP - PMM
    //  temp3 =  MMP + PPM    temp7 =  MMP - PPM
    //  temp4 =  PMP + MPM    temp8 =  PMP - MPM

    _temp4 = SET1_PD(0.5);
    _temp5 = MUL_PD(_temp5,_temp4);    // (PPP - MMM)*0.5
    _temp6 = MUL_PD(_temp6,_temp4);    // (MPP - PMM)*0.5
    _temp7 = MUL_PD(_temp7,_temp4);    // (MMP - PPM)*0.5
    _temp8 = MUL_PD(_temp8,_temp4);    // (PMP - MPM)*0.5

    _UZ  = ADD_PD(_UZ,_temp5);
    _UZ  = ADD_PD(_UZ,_temp6);
    _UZ  = ADD_PD(_UZ,_temp7);
    _UZ  = ADD_PD(_UZ,_temp8);
    _UY  = ADD_PD(_UY,_temp5);
    _UY  = ADD_PD(_UY,_temp6);
    _UY  = SUB_PD(_UY,_temp7);
    _UY  = SUB_PD(_UY,_temp8);
    _UX  = ADD_PD(_UX,_temp5);
    _UX  = SUB_PD(_UX,_temp6);
    _UX  = SUB_PD(_UX,_temp7);
    _UX  = ADD_PD(_UX,_temp8);  //at i,i+1,i+2,i+3

    _temp4 = SET1_PD(0.75);
    _temp9 = MUL_PD(_temp9,_temp4);
    _temp1 = MUL_PD(_temp1,_temp4);
    _temp2 = MUL_PD(_temp2,_temp4);
    _temp3 = MUL_PD(_temp3,_temp4);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);


    /////////////////////
    // BCCP1 and BCCM1  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_P1_P1);

    _temp1 = LOAD_PD(pointer1);    // All BCCP1 at  i
    _temp2 = LOAD_PD(pointer1+4);      //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);      //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);     //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains PPP  from i,i+1,i+2,i+3
    // temp2 now contains MPP  from i,i+1,i+2,i+3
    // temp3 now contains MMP  from i,i+1,i+2,i+3
    // temp4 now contains PMP  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1);

    _temp5 = LOAD_PD(pointer1);    // All BCCM1 at  i
    _temp6 = LOAD_PD(pointer1+4);      //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);      //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);     //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains MMM from i,i+1,i+2,i+3
    // temp6 now contains PMM from i,i+1,i+2,i+3
    // temp7 now contains PPM from i,i+1,i+2,i+3
    // temp8 now contains MPM from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9);
    _RHO = ADD_PD(_RHO,_temp1);
    _RHO = ADD_PD(_RHO,_temp2);
    _RHO = ADD_PD(_RHO,_temp3); //at i,i+1,i+2,i+3

    //  temp9 =  PPP1 + MMM1    temp5 =  PPP1 - MMM1
    //  temp1 =  MPP1 + PMM1    temp6 =  MPP1 - PMM1
    //  temp2 =  MMP1 + PPM1    temp7 =  MMP1 - PPM1
    //  temp3 =  PMP1 + MPM1    temp8 =  PMP1 - MPM1

    _UZ  = ADD_PD(_UZ,_temp5);
    _UZ  = ADD_PD(_UZ,_temp6);
    _UZ  = ADD_PD(_UZ,_temp7);
    _UZ  = ADD_PD(_UZ,_temp8);
    _UY  = ADD_PD(_UY,_temp5);
    _UY  = ADD_PD(_UY,_temp6);
    _UY  = SUB_PD(_UY,_temp7);
    _UY  = SUB_PD(_UY,_temp8);
    _UX  = ADD_PD(_UX,_temp5);
    _UX  = SUB_PD(_UX,_temp6);
    _UX  = SUB_PD(_UX,_temp7);
    _UX  = ADD_PD(_UX,_temp8);  //at i,i+1,i+2,i+3

    _UX = DIV_PD(_UX,_RHO);
    _UY = DIV_PD(_UY,_RHO);
    _UZ = DIV_PD(_UZ,_RHO);

    _temp4 = SET1_PD(3.0);
    _temp9 = MUL_PD(_temp9,_temp4);
    _temp1 = MUL_PD(_temp1,_temp4);
    _temp2 = MUL_PD(_temp2,_temp4);
    _temp3 = MUL_PD(_temp3,_temp4);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);

    _temp1  = MUL_PD(_UX,_UX);             //temp9  = uX*uX
    _temp2  = MUL_PD(_UY,_UY);             //temp10 = uY*uY
    _temp3  = ADD_PD(_temp1,_temp2);       //temp9  = uX*uX + uY*uY
    _temp4  = MUL_PD(_UZ,_UZ);             //temp10 = uZ*uZ
    _temp1  = ADD_PD(_temp3,_temp4);       //temp9  = uX*uX + uY*uY + uZ*uZ
    _temp2  = MUL_PD(_RHO,_temp1);         //temp9  = rho*(uX*uX + uY*uY + uZ*uZ)
    _THETA  = SUB_PD(_THETA,_temp2);       //theta  = theta - rho*(uX*uX + uY*uY + uZ*uZ)
    _temp1  = SET1_PD(1.0/3.0);
    _THETA  = MUL_PD(_THETA,_temp1);  //theta  = 1.0/3.0 *(theta - rho*(uX*uX + uY*uY + uZ*uZ))
    _THETA = DIV_PD(_THETA,_RHO);         // oneByRho

    STORE_PD(rho  ,_RHO  );
    STORE_PD(uX   ,_UX   );
    STORE_PD(uY   ,_UY   );
    STORE_PD(uZ   ,_UZ   );
    STORE_PD(theta,_THETA);
}


template <int N,int numblock, typename dataType1>
void getHydroMomentsFromNodeWithForce(lbmRD3Q41<dataType1> &lbModel,gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta,dataType1 F1, dataType1 F2, dataType1 F3,dataType1 dt)
{
    //         dataType1 sum(0.0);
    SIMD_REG _RHO, _UX, _UY, _UZ, _THETA;              //4
    SIMD_REG _temp1,_temp2,_temp3,_temp4;              //4
    SIMD_REG _temp5,_temp6,_temp7,_temp8;              //4
    SIMD_REG _temp9,_twoReg,_fourReg;                  //3   = 16 Registers in Total
    dataType1 tempArray[4] __attribute__ ((aligned(32)));
    dataType1  *pointer1;

    dataType1 pt5F1 = 0.5* dt * F1;
    dataType1 pt5F2 = 0.5* dt * F2;
    dataType1 pt5F3 = 0.5* dt * F3;
    dataType1 dot(0.0),csq(0.0);
    //////////////
    //  CENTER  //
    //////////////
    tempArray[0]=myGrid(i1  ,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    tempArray[1]=myGrid(i1+1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    tempArray[2]=myGrid(i1+2,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    tempArray[3]=myGrid(i1+3,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    pointer1 = &tempArray[0];
    _RHO = LOAD_PD(pointer1);


    ///////////////////
    //  SC1 and SC2  //
    ///////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1);
    _temp1 = LOAD_PD(pointer1);      // All SC1 at  i
    _temp2 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains ZZP1 from i,i+1,i+2,i+3
    // temp2 now contains ZZP2 from i,i+1,i+2,i+3
    // temp3 now contains ZP2Z from i,i+1,i+2,i+3
    // temp4 now contains P2ZZ from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1);
    _UZ    = LOAD_PD(pointer1);       // All SC2 at  i
    _temp6 = LOAD_PD(pointer1+4);     //     "      i+1
    _UY    = LOAD_PD(pointer1+8);     //     "      i+2
    _UX    = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_UZ,_temp6,_UY,_UX);
    // _UZ   now contains ZZM1 from i,i+1,i+2,i+3
    // temp6 now contains ZZM2 from i,i+1,i+2,i+3
    // _UY   now contains ZM2Z from i,i+1,i+2,i+3
    // _UX   now contains M2ZZ from i,i+1,i+2,i+3

    _THETA = ADD_PD(_temp1,_UZ)   ;  _UZ    = SUB_PD(_temp1,_UZ);
    _temp1 = ADD_PD(_temp2,_temp6);  _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_UY)   ;  _UY    = SUB_PD(_temp3,_UY);
    _temp3 = ADD_PD(_temp4,_UX)   ;  _UX    = SUB_PD(_temp4,_UX);

    _RHO = ADD_PD(_RHO,_THETA); // ZZP1 + ZZM1
    _RHO = ADD_PD(_RHO,_temp1); // ZZP1 + ZZM1 + ZZP2 + ZZM2
    _RHO = ADD_PD(_RHO,_temp2); // ZZP1 + ZZM1 + ZZP2 + ZZM2 + ZP2Z + ZM2Z
    _RHO = ADD_PD(_RHO,_temp3); // ZZP1 + ZZM1 + ZZP2 + ZZM2 + ZP2Z + ZM2Z + P2ZZ + M2ZZ  at i,i+1,i+2,i+3


    _twoReg = SET1_PD(2.0);
    _temp6  = MUL_PD(_temp6,_twoReg);    // (ZZP2-ZZM2)*2.0
    _UY     = MUL_PD(_UY   ,_twoReg);    // (ZP2Z-ZM2Z)*2.0
    _UX     = MUL_PD(_UX   ,_twoReg);    // (P2ZZ-M2ZZ)*2.0

    _UZ     = ADD_PD(_UZ,_temp6);        // (ZZP1-ZZM1) + (ZZP2-ZZM2)*2.0

    _fourReg = SET1_PD(4.0);
    _temp1   = MUL_PD(_temp1,_fourReg);    // (ZZP2+ZZM2)*4.0
    _temp2   = MUL_PD(_temp2,_fourReg);    // (ZP2Z+ZM2Z)*4.0
    _temp3   = MUL_PD(_temp3,_fourReg);    // (P2ZZ+M2ZZ)*4.0

    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);



    /////////////////////
    //  FCC1 and FCC2  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1);

    _temp1 = LOAD_PD(pointer1);    // All FCC1 at  i
    _temp2 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains ZPP  from i,i+1,i+2,i+3
    // temp2 now contains ZMP  from i,i+1,i+2,i+3
    // temp3 now contains PZP  from i,i+1,i+2,i+3
    // temp4 now contains MZP  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1);

    _temp5 = LOAD_PD(pointer1);      // All FCC2 at  i
    _temp6 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains ZMM from i,i+1,i+2,i+3
    // temp6 now contains ZPM from i,i+1,i+2,i+3
    // temp7 now contains MZM from i,i+1,i+2,i+3
    // temp8 now contains PZM from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9); // ZPP + ZMM
    _RHO = ADD_PD(_RHO,_temp1); // ZPP + ZMM + ZMP + ZPM
    _RHO = ADD_PD(_RHO,_temp2); // ZPP + ZMM + ZMP + ZPM + PZP + MZM
    _RHO = ADD_PD(_RHO,_temp3); // ZPP + ZMM + ZMP + ZPM + PZP + MZM + MZP + PZM at i,i+1,i+2,i+3

    //  temp9 =  ZPP + ZMM    temp5 =  ZPP - ZMM
    //  temp1 =  ZMP + ZPM    temp6 =  ZMP - ZPM
    //  temp2 =  PZP + MZM    temp7 =  PZP - MZM
    //  temp3 =  MZP + PZM    temp8 =  MZP - PZM

    _UZ  = ADD_PD(_UZ,_temp5);
    _UZ  = ADD_PD(_UZ,_temp6);
    _UZ  = ADD_PD(_UZ,_temp7);
    _UZ  = ADD_PD(_UZ,_temp8);
    _UY  = ADD_PD(_UY,_temp5);
    _UY  = SUB_PD(_UY,_temp6);
    _UX  = ADD_PD(_UX,_temp7);
    _UX  = SUB_PD(_UX,_temp8); //at i,i+1,i+2,i+3



    _temp9 = MUL_PD(_temp9,_twoReg);
    _temp1 = MUL_PD(_temp1,_twoReg);
    _temp2 = MUL_PD(_temp2,_twoReg);
    _temp3 = MUL_PD(_temp3,_twoReg);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);



    /////////////////////
    //  FCC3 and FCC4  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO);

    _temp1 = LOAD_PD(pointer1);      // All FCC3 at  i
    _temp2 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains PPZ  from i,i+1,i+2,i+3
    // temp2 now contains MPZ  from i,i+1,i+2,i+3
    // temp3 now contains ZPZ  from i,i+1,i+2,i+3
    // temp4 now contains PZZ  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO);

    _temp5 = LOAD_PD(pointer1);      // All FCC4 at  i
    _temp6 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains MMZ from i,i+1,i+2,i+3
    // temp6 now contains PMZ from i,i+1,i+2,i+3
    // temp7 now contains ZMZ from i,i+1,i+2,i+3
    // temp8 now contains MZZ from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9);
    _RHO = ADD_PD(_RHO,_temp1);
    _RHO = ADD_PD(_RHO,_temp2);
    _RHO = ADD_PD(_RHO,_temp3); //at i,i+1,i+2,i+3

    //  temp1 =  PPZ + MMZ    temp5 =  PPZ - MMZ
    //  temp2 =  MPZ + PMZ    temp6 =  MPZ - PMZ
    //  temp3 =  ZPZ + ZMZ    temp7 =  ZPZ - ZMZ
    //  temp4 =  PZZ + MZZ    temp8 =  PZZ - MZZ

    _UY  = ADD_PD(_UY,_temp5);
    _UY  = ADD_PD(_UY,_temp6);
    _UY  = ADD_PD(_UY,_temp7);
    _UX  = ADD_PD(_UX,_temp5);
    _UX  = SUB_PD(_UX,_temp6);
    _UX  = ADD_PD(_UX,_temp8); //at i,i+1,i+2,i+3


    _temp9 = MUL_PD(_temp9,_twoReg);
    _temp1 = MUL_PD(_temp1,_twoReg);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);


    /////////////////////
    // BCCP and BCCM  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P);

    _temp1 = LOAD_PD(pointer1);       // All BCCP at  i
    _temp2 = LOAD_PD(pointer1+4);     //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);     //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains PPP  from i,i+1,i+2,i+3
    // temp2 now contains MPP  from i,i+1,i+2,i+3
    // temp3 now contains MMP  from i,i+1,i+2,i+3
    // temp4 now contains PMP  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M);

    _temp5 = LOAD_PD(pointer1);       // All BCCM at  i
    _temp6 = LOAD_PD(pointer1+4);     //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);     //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains MMM from i,i+1,i+2,i+3
    // temp6 now contains PMM from i,i+1,i+2,i+3
    // temp7 now contains PPM from i,i+1,i+2,i+3
    // temp8 now contains MPM from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9);
    _RHO = ADD_PD(_RHO,_temp1);
    _RHO = ADD_PD(_RHO,_temp2);
    _RHO = ADD_PD(_RHO,_temp3); //at i,i+1,i+2,i+3

    //  temp1 =  PPP + MMM    temp5 =  PPP - MMM
    //  temp2 =  MPP + PMM    temp6 =  MPP - PMM
    //  temp3 =  MMP + PPM    temp7 =  MMP - PPM
    //  temp4 =  PMP + MPM    temp8 =  PMP - MPM

    _temp4 = SET1_PD(0.5);
    _temp5 = MUL_PD(_temp5,_temp4);    // (PPP - MMM)*0.5
    _temp6 = MUL_PD(_temp6,_temp4);    // (MPP - PMM)*0.5
    _temp7 = MUL_PD(_temp7,_temp4);    // (MMP - PPM)*0.5
    _temp8 = MUL_PD(_temp8,_temp4);    // (PMP - MPM)*0.5

    _UZ  = ADD_PD(_UZ,_temp5);
    _UZ  = ADD_PD(_UZ,_temp6);
    _UZ  = ADD_PD(_UZ,_temp7);
    _UZ  = ADD_PD(_UZ,_temp8);
    _UY  = ADD_PD(_UY,_temp5);
    _UY  = ADD_PD(_UY,_temp6);
    _UY  = SUB_PD(_UY,_temp7);
    _UY  = SUB_PD(_UY,_temp8);
    _UX  = ADD_PD(_UX,_temp5);
    _UX  = SUB_PD(_UX,_temp6);
    _UX  = SUB_PD(_UX,_temp7);
    _UX  = ADD_PD(_UX,_temp8);  //at i,i+1,i+2,i+3

    _temp4 = SET1_PD(0.75);
    _temp9 = MUL_PD(_temp9,_temp4);
    _temp1 = MUL_PD(_temp1,_temp4);
    _temp2 = MUL_PD(_temp2,_temp4);
    _temp3 = MUL_PD(_temp3,_temp4);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);


    /////////////////////
    // BCCP1 and BCCM1  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_P1_P1);

    _temp1 = LOAD_PD(pointer1);    // All BCCP1 at  i
    _temp2 = LOAD_PD(pointer1+4);      //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);      //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);     //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains PPP  from i,i+1,i+2,i+3
    // temp2 now contains MPP  from i,i+1,i+2,i+3
    // temp3 now contains MMP  from i,i+1,i+2,i+3
    // temp4 now contains PMP  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1);

    _temp5 = LOAD_PD(pointer1);    // All BCCM1 at  i
    _temp6 = LOAD_PD(pointer1+4);      //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);      //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);     //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains MMM from i,i+1,i+2,i+3
    // temp6 now contains PMM from i,i+1,i+2,i+3
    // temp7 now contains PPM from i,i+1,i+2,i+3
    // temp8 now contains MPM from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9);
    _RHO = ADD_PD(_RHO,_temp1);
    _RHO = ADD_PD(_RHO,_temp2);
    _RHO = ADD_PD(_RHO,_temp3); //at i,i+1,i+2,i+3

    //  temp9 =  PPP1 + MMM1    temp5 =  PPP1 - MMM1
    //  temp1 =  MPP1 + PMM1    temp6 =  MPP1 - PMM1
    //  temp2 =  MMP1 + PPM1    temp7 =  MMP1 - PPM1
    //  temp3 =  PMP1 + MPM1    temp8 =  PMP1 - MPM1

    _UZ  = ADD_PD(_UZ,_temp5);
    _UZ  = ADD_PD(_UZ,_temp6);
    _UZ  = ADD_PD(_UZ,_temp7);
    _UZ  = ADD_PD(_UZ,_temp8);
    _UY  = ADD_PD(_UY,_temp5);
    _UY  = ADD_PD(_UY,_temp6);
    _UY  = SUB_PD(_UY,_temp7);
    _UY  = SUB_PD(_UY,_temp8);
    _UX  = ADD_PD(_UX,_temp5);
    _UX  = SUB_PD(_UX,_temp6);
    _UX  = SUB_PD(_UX,_temp7);
    _UX  = ADD_PD(_UX,_temp8);  //at i,i+1,i+2,i+3

    _UX = DIV_PD(_UX,_RHO);
    _UY = DIV_PD(_UY,_RHO);
    _UZ = DIV_PD(_UZ,_RHO);

    _temp5 = SET1_PD(pt5F1);        // LINES FOR FORCE METHOD ARE ADDED HERE
    _temp6 = SET1_PD(pt5F2);        // LINES FOR FORCE METHOD ARE ADDED HERE
    _temp7 = SET1_PD(pt5F3);        // LINES FOR FORCE METHOD ARE ADDED HERE
    // LINES FOR FORCE METHOD ARE ADDED HERE
    _UX = ADD_PD(_UX,_temp5);       // LINES FOR FORCE METHOD ARE ADDED HERE
    _UY = ADD_PD(_UY,_temp6);       // LINES FOR FORCE METHOD ARE ADDED HERE
    _UZ = ADD_PD(_UZ,_temp7);       // LINES FOR FORCE METHOD ARE ADDED HERE


    _temp4 = SET1_PD(3.0);
    _temp9 = MUL_PD(_temp9,_temp4);
    _temp1 = MUL_PD(_temp1,_temp4);
    _temp2 = MUL_PD(_temp2,_temp4);
    _temp3 = MUL_PD(_temp3,_temp4);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);



    _temp1  = MUL_PD(_UX,_UX);             //temp9  = uX*uX
    _temp2  = MUL_PD(_UY,_UY);             //temp10 = uY*uY
    _temp3  = ADD_PD(_temp1,_temp2);       //temp9  = uX*uX + uY*uY
    _temp4  = MUL_PD(_UZ,_UZ);             //temp10 = uZ*uZ
    _temp1  = ADD_PD(_temp3,_temp4);       //temp9  = uX*uX + uY*uY + uZ*uZ
    _temp2  = MUL_PD(_RHO,_temp1);         //temp9  = rho*(uX*uX + uY*uY + uZ*uZ)
    _THETA  = SUB_PD(_THETA,_temp2);       //theta  = theta - rho*(uX*uX + uY*uY + uZ*uZ)
    _temp1  = SET1_PD(1.0/3.0);
    _THETA  = MUL_PD(_THETA,_temp1);  //theta  = 1.0/3.0 *(theta - rho*(uX*uX + uY*uY + uZ*uZ))
    _THETA = DIV_PD(_THETA,_RHO);         // oneByRho

    STORE_PD(rho  ,_RHO  );
    STORE_PD(uX   ,_UX   );
    STORE_PD(uY   ,_UY   );
    STORE_PD(uZ   ,_UZ   );
    STORE_PD(theta,_THETA);

    for (int index=0;index<VECT_LENGTH;index++)                                                        // LINES FOR FORCE METHOD ARE ADDED HERE
    {
        for (int dv=0;dv<lbModel.dvN;dv++)                                                                        // LINES FOR FORCE METHOD ARE ADDED HERE
        {                                                                                                         // LINES FOR FORCE METHOD ARE ADDED HERE
            csq = lbModel.cx[dv]*lbModel.cx[dv] + lbModel.cy[dv]*lbModel.cy[dv] + lbModel.cz[dv]*lbModel.cz[dv];     // LINES FOR FORCE METHOD ARE ADDED HERE
            dot = csq*(F1*lbModel.cx[dv]+F2*lbModel.cy[dv]+F3*lbModel.cz[dv])*1.68103*1.68103 ;                                      // LINES FOR FORCE METHOD ARE ADDED HERE
            // LINES FOR FORCE METHOD ARE ADDED HERE
            theta[index] +=  (0.5*lbModel.wt[dv]*dot) ;                                                             // LINES FOR FORCE METHOD ARE ADDED HERE
        }                                                                                                        // LINES FOR FORCE METHOD ARE ADDED HERE
    }                                                                                                         // LINES FOR FORCE METHOD ARE ADDED HERE



}



template <int N,int numblock, typename dataType1>
void getHydroMomentsFromCellWithForce(lbmRD3Q41<dataType1> &lbModel,gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,int i1, int i2, int i3, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta,dataType1 F1, dataType1 F2, dataType1 F3,dataType1 dt)
{
    //         dataType1 sum(0.0);
    SIMD_REG _RHO, _UX, _UY, _UZ, _THETA;              //4
    SIMD_REG _temp1,_temp2,_temp3,_temp4;              //4
    SIMD_REG _temp5,_temp6,_temp7,_temp8;              //4
    SIMD_REG _temp9,_twoReg,_fourReg;                  //3   = 16 Registers in Total
    dataType1 tempArray[4] __attribute__ ((aligned(32)));
    dataType1  *pointer1;

    dataType1 pt5F1 = 0.5* dt * F1;
    dataType1 pt5F2 = 0.5* dt * F2;
    dataType1 pt5F3 = 0.5* dt * F3;
    dataType1 dot(0.0),csq(0.0);

    //////////////
    //  CENTER  //
    //////////////
    tempArray[0]=myGrid(i1  ,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    tempArray[1]=myGrid(i1+1,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    tempArray[2]=myGrid(i1+2,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    tempArray[3]=myGrid(i1+3,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO);
    pointer1 = &tempArray[0];
    _RHO = LOAD_PD(pointer1);


    ///////////////////
    //  SC1 and SC2  //
    ///////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1);
    _temp1 = LOAD_PD(pointer1);      // All SC1 at  i
    _temp2 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains ZZP1 from i,i+1,i+2,i+3
    // temp2 now contains ZZP2 from i,i+1,i+2,i+3
    // temp3 now contains ZP2Z from i,i+1,i+2,i+3
    // temp4 now contains P2ZZ from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1);

    _UZ    = LOAD_PD(pointer1);       // All SC2 at  i
    _temp6 = LOAD_PD(pointer1+4);     //     "      i+1
    _UY    = LOAD_PD(pointer1+8);     //     "      i+2
    _UX    = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_UZ,_temp6,_UY,_UX);
    // temp5 now contains ZZM1 from i,i+1,i+2,i+3
    // temp6 now contains ZZM2 from i,i+1,i+2,i+3
    // temp7 now contains ZM2Z from i,i+1,i+2,i+3
    // temp8 now contains M2ZZ from i,i+1,i+2,i+3

    _THETA = ADD_PD(_temp1,_UZ)   ;  _UZ    = SUB_PD(_temp1,_UZ);
    _temp1 = ADD_PD(_temp2,_temp6);  _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_UY)   ;  _UY    = SUB_PD(_temp3,_UY);
    _temp3 = ADD_PD(_temp4,_UX)   ;  _UX    = SUB_PD(_temp4,_UX);

    _RHO = ADD_PD(_RHO,_THETA); // ZZP1 + ZZM1
    _RHO = ADD_PD(_RHO,_temp1); // ZZP1 + ZZM1 + ZZP2 + ZZM2
    _RHO = ADD_PD(_RHO,_temp2); // ZZP1 + ZZM1 + ZZP2 + ZZM2 + ZP2Z + ZM2Z
    _RHO = ADD_PD(_RHO,_temp3); // ZZP1 + ZZM1 + ZZP2 + ZZM2 + ZP2Z + ZM2Z + P2ZZ + M2ZZ  at i,i+1,i+2,i+3


    _twoReg = SET1_PD(2.0);
    _temp6  = MUL_PD(_temp6,_twoReg);    // (ZZP2-ZZM2)*2.0
    _UY     = MUL_PD(_UY   ,_twoReg);    // (ZP2Z-ZM2Z)*2.0
    _UX     = MUL_PD(_UX   ,_twoReg);    // (P2ZZ-M2ZZ)*2.0

    _UZ     = ADD_PD(_UZ,_temp6);        // (ZZP1-ZZM1) + (ZZP2-ZZM2)*2.0

    _fourReg = SET1_PD(4.0);
    _temp1   = MUL_PD(_temp1,_fourReg);    // (ZZP2+ZZM2)*4.0
    _temp2   = MUL_PD(_temp2,_fourReg);    // (ZP2Z+ZM2Z)*4.0
    _temp3   = MUL_PD(_temp3,_fourReg);    // (P2ZZ+M2ZZ)*4.0

    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);



    /////////////////////
    //  FCC1 and FCC2  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1);

    _temp1 = LOAD_PD(pointer1);    // All FCC1 at  i
    _temp2 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains ZPP  from i,i+1,i+2,i+3
    // temp2 now contains ZMP  from i,i+1,i+2,i+3
    // temp3 now contains PZP  from i,i+1,i+2,i+3
    // temp4 now contains MZP  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1);

    _temp5 = LOAD_PD(pointer1);      // All FCC2 at  i
    _temp6 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains ZMM from i,i+1,i+2,i+3
    // temp6 now contains ZPM from i,i+1,i+2,i+3
    // temp7 now contains MZM from i,i+1,i+2,i+3
    // temp8 now contains PZM from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9); // ZPP + ZMM
    _RHO = ADD_PD(_RHO,_temp1); // ZPP + ZMM + ZMP + ZPM
    _RHO = ADD_PD(_RHO,_temp2); // ZPP + ZMM + ZMP + ZPM + PZP + MZM
    _RHO = ADD_PD(_RHO,_temp3); // ZPP + ZMM + ZMP + ZPM + PZP + MZM + MZP + PZM at i,i+1,i+2,i+3

    //  temp9 =  ZPP + ZMM    temp5 =  ZPP - ZMM
    //  temp1 =  ZMP + ZPM    temp6 =  ZMP - ZPM
    //  temp2 =  PZP + MZM    temp7 =  PZP - MZM
    //  temp3 =  MZP + PZM    temp8 =  MZP - PZM

    _UZ  = ADD_PD(_UZ,_temp5);
    _UZ  = ADD_PD(_UZ,_temp6);
    _UZ  = ADD_PD(_UZ,_temp7);
    _UZ  = ADD_PD(_UZ,_temp8);
    _UY  = ADD_PD(_UY,_temp5);
    _UY  = SUB_PD(_UY,_temp6);
    _UX  = ADD_PD(_UX,_temp7);
    _UX  = SUB_PD(_UX,_temp8); //at i,i+1,i+2,i+3



    _temp9 = MUL_PD(_temp9,_twoReg);
    _temp1 = MUL_PD(_temp1,_twoReg);
    _temp2 = MUL_PD(_temp2,_twoReg);
    _temp3 = MUL_PD(_temp3,_twoReg);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);



    /////////////////////
    //  FCC3 and FCC4  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO);

    _temp1 = LOAD_PD(pointer1);      // All FCC3 at  i
    _temp2 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains PPZ  from i,i+1,i+2,i+3
    // temp2 now contains MPZ  from i,i+1,i+2,i+3
    // temp3 now contains ZPZ  from i,i+1,i+2,i+3
    // temp4 now contains PZZ  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO);

    _temp5 = LOAD_PD(pointer1);      // All FCC4 at  i
    _temp6 = LOAD_PD(pointer1+4);    //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);    //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);   //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains MMZ from i,i+1,i+2,i+3
    // temp6 now contains PMZ from i,i+1,i+2,i+3
    // temp7 now contains ZMZ from i,i+1,i+2,i+3
    // temp8 now contains MZZ from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9);
    _RHO = ADD_PD(_RHO,_temp1);
    _RHO = ADD_PD(_RHO,_temp2);
    _RHO = ADD_PD(_RHO,_temp3); //at i,i+1,i+2,i+3

    //  temp1 =  PPZ + MMZ    temp5 =  PPZ - MMZ
    //  temp2 =  MPZ + PMZ    temp6 =  MPZ - PMZ
    //  temp3 =  ZPZ + ZMZ    temp7 =  ZPZ - ZMZ
    //  temp4 =  PZZ + MZZ    temp8 =  PZZ - MZZ

    _UY  = ADD_PD(_UY,_temp5);
    _UY  = ADD_PD(_UY,_temp6);
    _UY  = ADD_PD(_UY,_temp7);
    _UX  = ADD_PD(_UX,_temp5);
    _UX  = SUB_PD(_UX,_temp6);
    _UX  = ADD_PD(_UX,_temp8); //at i,i+1,i+2,i+3


    _temp9 = MUL_PD(_temp9,_twoReg);
    _temp1 = MUL_PD(_temp1,_twoReg);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);


    /////////////////////
    // BCCP and BCCM  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P);

    _temp1 = LOAD_PD(pointer1);       // All BCCP at  i
    _temp2 = LOAD_PD(pointer1+4);     //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);     //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains PPP  from i,i+1,i+2,i+3
    // temp2 now contains MPP  from i,i+1,i+2,i+3
    // temp3 now contains MMP  from i,i+1,i+2,i+3
    // temp4 now contains PMP  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M);

    _temp5 = LOAD_PD(pointer1);       // All BCCM at  i
    _temp6 = LOAD_PD(pointer1+4);     //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);     //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);    //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains MMM from i,i+1,i+2,i+3
    // temp6 now contains PMM from i,i+1,i+2,i+3
    // temp7 now contains PPM from i,i+1,i+2,i+3
    // temp8 now contains MPM from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9);
    _RHO = ADD_PD(_RHO,_temp1);
    _RHO = ADD_PD(_RHO,_temp2);
    _RHO = ADD_PD(_RHO,_temp3); //at i,i+1,i+2,i+3

    //  temp1 =  PPP + MMM    temp5 =  PPP - MMM
    //  temp2 =  MPP + PMM    temp6 =  MPP - PMM
    //  temp3 =  MMP + PPM    temp7 =  MMP - PPM
    //  temp4 =  PMP + MPM    temp8 =  PMP - MPM

    _temp4 = SET1_PD(0.5);
    _temp5 = MUL_PD(_temp5,_temp4);    // (PPP - MMM)*0.5
    _temp6 = MUL_PD(_temp6,_temp4);    // (MPP - PMM)*0.5
    _temp7 = MUL_PD(_temp7,_temp4);    // (MMP - PPM)*0.5
    _temp8 = MUL_PD(_temp8,_temp4);    // (PMP - MPM)*0.5

    _UZ  = ADD_PD(_UZ,_temp5);
    _UZ  = ADD_PD(_UZ,_temp6);
    _UZ  = ADD_PD(_UZ,_temp7);
    _UZ  = ADD_PD(_UZ,_temp8);
    _UY  = ADD_PD(_UY,_temp5);
    _UY  = ADD_PD(_UY,_temp6);
    _UY  = SUB_PD(_UY,_temp7);
    _UY  = SUB_PD(_UY,_temp8);
    _UX  = ADD_PD(_UX,_temp5);
    _UX  = SUB_PD(_UX,_temp6);
    _UX  = SUB_PD(_UX,_temp7);
    _UX  = ADD_PD(_UX,_temp8);  //at i,i+1,i+2,i+3

    _temp4 = SET1_PD(0.75);
    _temp9 = MUL_PD(_temp9,_temp4);
    _temp1 = MUL_PD(_temp1,_temp4);
    _temp2 = MUL_PD(_temp2,_temp4);
    _temp3 = MUL_PD(_temp3,_temp4);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);


    /////////////////////
    // BCCP1 and BCCM1  //
    /////////////////////
    pointer1 = &myGrid(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_P1_P1);

    _temp1 = LOAD_PD(pointer1);    // All BCCP1 at  i
    _temp2 = LOAD_PD(pointer1+4);      //     "      i+1
    _temp3 = LOAD_PD(pointer1+8);      //     "      i+2
    _temp4 = LOAD_PD(pointer1+12);     //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp1,_temp2,_temp3,_temp4);
    // temp1 now contains PPP  from i,i+1,i+2,i+3
    // temp2 now contains MPP  from i,i+1,i+2,i+3
    // temp3 now contains MMP  from i,i+1,i+2,i+3
    // temp4 now contains PMP  from i,i+1,i+2,i+3

    pointer1 = &myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1);

    _temp5 = LOAD_PD(pointer1);    // All BCCM1 at  i
    _temp6 = LOAD_PD(pointer1+4);      //     "      i+1
    _temp7 = LOAD_PD(pointer1+8);      //     "      i+2
    _temp8 = LOAD_PD(pointer1+12);     //     "      i+3

    //Transpose them
    AVX_TRANSPOSE4_PD(_temp5,_temp6,_temp7,_temp8);
    // temp5 now contains MMM from i,i+1,i+2,i+3
    // temp6 now contains PMM from i,i+1,i+2,i+3
    // temp7 now contains PPM from i,i+1,i+2,i+3
    // temp8 now contains MPM from i,i+1,i+2,i+3

    _temp9 = ADD_PD(_temp1,_temp5); _temp5 = SUB_PD(_temp1,_temp5);
    _temp1 = ADD_PD(_temp2,_temp6); _temp6 = SUB_PD(_temp2,_temp6);
    _temp2 = ADD_PD(_temp3,_temp7); _temp7 = SUB_PD(_temp3,_temp7);
    _temp3 = ADD_PD(_temp4,_temp8); _temp8 = SUB_PD(_temp4,_temp8);

    _RHO = ADD_PD(_RHO,_temp9);
    _RHO = ADD_PD(_RHO,_temp1);
    _RHO = ADD_PD(_RHO,_temp2);
    _RHO = ADD_PD(_RHO,_temp3); //at i,i+1,i+2,i+3

    //  temp9 =  PPP1 + MMM1    temp5 =  PPP1 - MMM1
    //  temp1 =  MPP1 + PMM1    temp6 =  MPP1 - PMM1
    //  temp2 =  MMP1 + PPM1    temp7 =  MMP1 - PPM1
    //  temp3 =  PMP1 + MPM1    temp8 =  PMP1 - MPM1

    _UZ  = ADD_PD(_UZ,_temp5);
    _UZ  = ADD_PD(_UZ,_temp6);
    _UZ  = ADD_PD(_UZ,_temp7);
    _UZ  = ADD_PD(_UZ,_temp8);
    _UY  = ADD_PD(_UY,_temp5);
    _UY  = ADD_PD(_UY,_temp6);
    _UY  = SUB_PD(_UY,_temp7);
    _UY  = SUB_PD(_UY,_temp8);
    _UX  = ADD_PD(_UX,_temp5);
    _UX  = SUB_PD(_UX,_temp6);
    _UX  = SUB_PD(_UX,_temp7);
    _UX  = ADD_PD(_UX,_temp8);  //at i,i+1,i+2,i+3

    _UX = DIV_PD(_UX,_RHO);
    _UY = DIV_PD(_UY,_RHO);
    _UZ = DIV_PD(_UZ,_RHO);

    _temp5 = SET1_PD(pt5F1);        // LINES FOR FORCE METHOD ARE ADDED HERE
    _temp6 = SET1_PD(pt5F2);        // LINES FOR FORCE METHOD ARE ADDED HERE
    _temp7 = SET1_PD(pt5F3);        // LINES FOR FORCE METHOD ARE ADDED HERE
    // LINES FOR FORCE METHOD ARE ADDED HERE
    _UX = ADD_PD(_UX,_temp5);       // LINES FOR FORCE METHOD ARE ADDED HERE
    _UY = ADD_PD(_UY,_temp6);       // LINES FOR FORCE METHOD ARE ADDED HERE
    _UZ = ADD_PD(_UZ,_temp7);       // LINES FOR FORCE METHOD ARE ADDED HERE

    _temp4 = SET1_PD(3.0);
    _temp9 = MUL_PD(_temp9,_temp4);
    _temp1 = MUL_PD(_temp1,_temp4);
    _temp2 = MUL_PD(_temp2,_temp4);
    _temp3 = MUL_PD(_temp3,_temp4);

    _THETA = ADD_PD(_THETA,_temp9);
    _THETA = ADD_PD(_THETA,_temp1);
    _THETA = ADD_PD(_THETA,_temp2);
    _THETA = ADD_PD(_THETA,_temp3);

    _temp1  = MUL_PD(_UX,_UX);             //temp9  = uX*uX
    _temp2  = MUL_PD(_UY,_UY);             //temp10 = uY*uY
    _temp3  = ADD_PD(_temp1,_temp2);       //temp9  = uX*uX + uY*uY
    _temp4  = MUL_PD(_UZ,_UZ);             //temp10 = uZ*uZ
    _temp1  = ADD_PD(_temp3,_temp4);       //temp9  = uX*uX + uY*uY + uZ*uZ
    _temp2  = MUL_PD(_RHO,_temp1);         //temp9  = rho*(uX*uX + uY*uY + uZ*uZ)
    _THETA  = SUB_PD(_THETA,_temp2);       //theta  = theta - rho*(uX*uX + uY*uY + uZ*uZ)
    _temp1  = SET1_PD(1.0/3.0);
    _THETA  = MUL_PD(_THETA,_temp1);  //theta  = 1.0/3.0 *(theta - rho*(uX*uX + uY*uY + uZ*uZ))
    _THETA = DIV_PD(_THETA,_RHO);         // oneByRho

    STORE_PD(rho  ,_RHO  );
    STORE_PD(uX   ,_UX   );
    STORE_PD(uY   ,_UY   );
    STORE_PD(uZ   ,_UZ   );
    STORE_PD(theta,_THETA);

    for (int dv=0;dv<lbModel.dvN;dv++)                                                                        // LINES FOR FORCE METHOD ARE ADDED HERE
    {                                                                                                         // LINES FOR FORCE METHOD ARE ADDED HERE
        csq = lbModel.cx[dv]*lbModel.cx[dv] + lbModel.cy[dv]*lbModel.cy[dv] + lbModel.cz[dv]*lbModel.cz[dv];     // LINES FOR FORCE METHOD ARE ADDED HERE
        dot = csq*(F1*lbModel.cx[dv]+F2*lbModel.cy[dv]+F3*lbModel.cz[dv])*1.68103*1.68103 ;                                      // LINES FOR FORCE METHOD ARE ADDED HERE
        for (int index=0;index<VECT_LENGTH;index++)                                                              // LINES FOR FORCE METHOD ARE ADDED HERE
        {                                                                                                        // LINES FOR FORCE METHOD ARE ADDED HERE
            theta[index] +=  (0.5*lbModel.wt[dv]*dot) ;                                                             // LINES FOR FORCE METHOD ARE ADDED HERE
        }                                                                                                        // LINES FOR FORCE METHOD ARE ADDED HERE
    }                                                                                                         // LINES FOR FORCE METHOD ARE ADDED HERE
}

  template<typename dataType1>
  void getHydroMomentSinglePoint(lbmRD3Q41<dataType1> &lbModel, dataType1 &rho, dataType1 &uX, dataType1 &uY, dataType1 &uZ, dataType1 &theta)
  {
    dataType1 sum(0.0);
    rho=0.0;uX=0.0;uY=0.0;uZ=0.0;theta=0.0;

      rho = lbModel.fTemp0[0][lbModel.CENTER_DV_ZERO_ZERO_ZERO];


    // G1 and G2

      sum       = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO]+lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
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

  template<typename dataType1>
  void getHydroMomentSinglePointWithForce(lbmRD3Q41<dataType1> &lbModel, dataType1 &rho, dataType1 &uX, dataType1 &uY, dataType1 &uZ, dataType1 &theta,dataType1 F1, dataType1 F2, dataType1 F3,dataType1 dt )
  {
    dataType1 sum(0.0),csq(0.0),dot(0.0);
    rho=0.0;uX=0.0;uY=0.0;uZ=0.0;theta=0.0;

      rho = lbModel.fTemp0[0][lbModel.CENTER_DV_ZERO_ZERO_ZERO];


    // G1 and G2

      sum       = lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO]+lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO];
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


      uX += 0.5*dt*F1 ;
      uY += 0.5*dt*F2 ;
      uZ += 0.5*dt*F3 ;

//       uX *= lbModel.cellSize;
//       uY *= lbModel.cellSize;
//       uZ *= lbModel.cellSize;

      theta = oneByThree * (theta - rho*(uX*uX + uY*uY + uZ*uZ));
      theta *= oneByRho;

//       theta *= lbModel.cellSize*lbModel.cellSize;

      for (int dv=0;dv<lbModel.dvN;dv++)
      {
       csq = (lbModel.cx[dv]*lbModel.cx[dv] + lbModel.cy[dv]*lbModel.cy[dv] + lbModel.cz[dv]*lbModel.cz[dv]);
       dot = csq*(F1*lbModel.cx[dv]+F2*lbModel.cy[dv]+F3*lbModel.cz[dv]) ;
       theta +=  (0.5*lbModel.wt[dv]*dot) ;
      }

  }



template<typename dataType1>
void getGrad10Moment(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx)
{


    dataType1 sum(0.0),point25sum(0.0);

    for(int i=0;i<VECT_LENGTH;i++)
    {
        rho[i] = lbModel.fTemp0[lbModel.CENTER_DV_ZERO_ZERO_ZERO][i];
    }

    // G1 and G2
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp1[lbModel.G1_DV_P2_ZERO_ZERO][i]+lbModel.fTemp2[lbModel.G2_DV_M2_ZERO_ZERO][i];
        rho[i]   += sum;
        Pxx[i]    = 4.0*sum;
        uX[i]     = 2.0*(lbModel.fTemp1[lbModel.G1_DV_P2_ZERO_ZERO][i]-lbModel.fTemp2[lbModel.G2_DV_M2_ZERO_ZERO][i]);

        sum       = lbModel.fTemp1[lbModel.G1_DV_ZERO_P2_ZERO][i]+lbModel.fTemp2[lbModel.G2_DV_ZERO_M2_ZERO][i];
        rho[i]   += sum;
        Pyy[i]    = 4.0*sum;
        uY[i]     = 2.0*(lbModel.fTemp1[lbModel.G1_DV_ZERO_P2_ZERO][i]-lbModel.fTemp2[lbModel.G2_DV_ZERO_M2_ZERO][i]);

        sum       = lbModel.fTemp1[lbModel.G1_DV_ZERO_ZERO_P2][i]+lbModel.fTemp2[lbModel.G2_DV_ZERO_ZERO_M2][i];
        rho[i]   += sum;
        Pzz[i]    = 4.0*sum;
        uZ[i]     = 2.0*(lbModel.fTemp1[lbModel.G1_DV_ZERO_ZERO_P2][i]-lbModel.fTemp2[lbModel.G2_DV_ZERO_ZERO_M2][i]);

        sum       = lbModel.fTemp1[lbModel.G1_DV_ZERO_ZERO_P1][i]+lbModel.fTemp2[lbModel.G2_DV_ZERO_ZERO_M1][i];
        rho[i]   += sum;
        Pzz[i]   += sum;
        uZ[i]    += lbModel.fTemp1[lbModel.G1_DV_ZERO_ZERO_P1][i]-lbModel.fTemp2[lbModel.G2_DV_ZERO_ZERO_M1][i];
    }

    // G3 and G4
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp3[lbModel.G3_DV_ZERO_P1_P1][i]+lbModel.fTemp4[lbModel.G4_DV_ZERO_M1_M1][i];
        rho[i]   += sum;
        Pyy[i]   += sum;
        Pyz[i]    = sum;
        Pzz[i]   += sum;
        uY[i]    += (lbModel.fTemp3[lbModel.G3_DV_ZERO_P1_P1][i]-lbModel.fTemp4[lbModel.G4_DV_ZERO_M1_M1][i]);
        uZ[i]    += (lbModel.fTemp3[lbModel.G3_DV_ZERO_P1_P1][i]-lbModel.fTemp4[lbModel.G4_DV_ZERO_M1_M1][i]);

        sum       = lbModel.fTemp3[lbModel.G3_DV_ZERO_M1_P1 ][i]+lbModel.fTemp4[lbModel.G4_DV_ZERO_P1_M1][i];
        rho[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pyz[i]   -= sum;
        uY[i]    -= (lbModel.fTemp3[lbModel.G3_DV_ZERO_M1_P1][i]-lbModel.fTemp4[lbModel.G4_DV_ZERO_P1_M1][i]);
        uZ[i]    += (lbModel.fTemp3[lbModel.G3_DV_ZERO_M1_P1][i]-lbModel.fTemp4[lbModel.G4_DV_ZERO_P1_M1][i]);

        sum       = lbModel.fTemp3[lbModel.G3_DV_P1_ZERO_P1][i]+lbModel.fTemp4[lbModel.G4_DV_M1_ZERO_M1][i];
        rho[i]   += sum;
        Pxx[i]   += sum;
        Pzz[i]   += sum;
        Pzx[i]    = sum;
        uX[i]    += (lbModel.fTemp3[lbModel.G3_DV_P1_ZERO_P1][i]-lbModel.fTemp4[lbModel.G4_DV_M1_ZERO_M1][i]);
        uZ[i]    += (lbModel.fTemp3[lbModel.G3_DV_P1_ZERO_P1][i]-lbModel.fTemp4[lbModel.G4_DV_M1_ZERO_M1][i]);

        sum       = lbModel.fTemp3[lbModel.G3_DV_M1_ZERO_P1][i]+lbModel.fTemp4[lbModel.G4_DV_P1_ZERO_M1][i];
        rho[i]   += sum;
        Pxx[i]   += sum;
        Pzz[i]   += sum;
        Pzx[i]   -= sum;
        uX[i]    -= (lbModel.fTemp3[lbModel.G3_DV_M1_ZERO_P1][i]-lbModel.fTemp4[lbModel.G4_DV_P1_ZERO_M1][i]);
        uZ[i]    += (lbModel.fTemp3[lbModel.G3_DV_M1_ZERO_P1][i]-lbModel.fTemp4[lbModel.G4_DV_P1_ZERO_M1][i]);
    }

    // G5 and G6
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp5[lbModel.G5_DV_P1_P1_ZERO][i]+lbModel.fTemp6[lbModel.G6_DV_M1_M1_ZERO][i];
        rho[i]   += sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pxy[i]   += sum;
        uX[i]    += (lbModel.fTemp5[lbModel.G5_DV_P1_P1_ZERO][i]-lbModel.fTemp6[lbModel.G6_DV_M1_M1_ZERO][i]);
        uY[i]    += (lbModel.fTemp5[lbModel.G5_DV_P1_P1_ZERO][i]-lbModel.fTemp6[lbModel.G6_DV_M1_M1_ZERO][i]);

        sum       = lbModel.fTemp5[lbModel.G5_DV_M1_P1_ZERO][i]+lbModel.fTemp6[lbModel.G6_DV_P1_M1_ZERO][i];
        rho[i]   += sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pxy[i]   -= sum;
        uX[i]    += (lbModel.fTemp6[lbModel.G6_DV_P1_M1_ZERO][i]-lbModel.fTemp5[lbModel.G5_DV_M1_P1_ZERO][i]);
        uY[i]    -= (lbModel.fTemp6[lbModel.G6_DV_P1_M1_ZERO][i]-lbModel.fTemp5[lbModel.G5_DV_M1_P1_ZERO][i]);

        sum      = lbModel.fTemp5[lbModel.G5_DV_ZERO_P1_ZERO][i]+lbModel.fTemp6[lbModel.G6_DV_ZERO_M1_ZERO][i];
        rho[i]   += sum;
        Pyy[i]   += sum;
        uY[i]    += (lbModel.fTemp5[lbModel.G5_DV_ZERO_P1_ZERO][i]-lbModel.fTemp6[lbModel.G6_DV_ZERO_M1_ZERO][i]);

        sum       = lbModel.fTemp5[lbModel.G5_DV_P1_ZERO_ZERO][i]+lbModel.fTemp6[lbModel.G6_DV_M1_ZERO_ZERO][i];
        rho[i]   += sum;
        Pxx[i]   += sum;
        uX[i]    += (lbModel.fTemp5[lbModel.G5_DV_P1_ZERO_ZERO][i]-lbModel.fTemp6[lbModel.G6_DV_M1_ZERO_ZERO][i]);
    }

    // G7 and G8
    const dataType1 threeByFour=3.0/4.0;
    for(int i=0;i< VECT_LENGTH;i++)
    {
        sum        = lbModel.fTemp7[lbModel.G7_DV_P_P_P][i]+lbModel.fTemp8[lbModel.G8_DV_M_M_M][i];
        rho[i]    += sum;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    += point25sum;
        Pyz[i]    += point25sum;
        Pzx[i]    += point25sum;
        uX[i]     += 0.5*(lbModel.fTemp7[lbModel.G7_DV_P_P_P][i]-lbModel.fTemp8[lbModel.G8_DV_M_M_M][i]);
        uY[i]     += 0.5*(lbModel.fTemp7[lbModel.G7_DV_P_P_P][i]-lbModel.fTemp8[lbModel.G8_DV_M_M_M][i]);
        uZ[i]     += 0.5*(lbModel.fTemp7[lbModel.G7_DV_P_P_P][i]-lbModel.fTemp8[lbModel.G8_DV_M_M_M][i]);

        sum        = lbModel.fTemp7[lbModel.G7_DV_M_P_P][i]+lbModel.fTemp8[lbModel.G8_DV_P_M_M][i];
        rho[i]    += sum;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    -= point25sum;
        Pyz[i]    += point25sum;
        Pzx[i]    -= point25sum;
        uX[i]     -= 0.5*(lbModel.fTemp7[lbModel.G7_DV_M_P_P][i]-lbModel.fTemp8[lbModel.G8_DV_P_M_M][i]);
        uY[i]     += 0.5*(lbModel.fTemp7[lbModel.G7_DV_M_P_P][i]-lbModel.fTemp8[lbModel.G8_DV_P_M_M][i]);
        uZ[i]     += 0.5*(lbModel.fTemp7[lbModel.G7_DV_M_P_P][i]-lbModel.fTemp8[lbModel.G8_DV_P_M_M][i]);


        sum        = lbModel.fTemp7[lbModel.G7_DV_M_M_P][i]+lbModel.fTemp8[lbModel.G8_DV_P_P_M][i];
        rho[i]    += sum;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    += point25sum;
        Pyz[i]    -= point25sum;
        Pzx[i]    -= point25sum;
        uX[i]     -= 0.5*(lbModel.fTemp7[lbModel.G7_DV_M_M_P][i]-lbModel.fTemp8[lbModel.G8_DV_P_P_M][i]);
        uY[i]     -= 0.5*(lbModel.fTemp7[lbModel.G7_DV_M_M_P][i]-lbModel.fTemp8[lbModel.G8_DV_P_P_M][i]);
        uZ[i]     += 0.5*(lbModel.fTemp7[lbModel.G7_DV_M_M_P][i]-lbModel.fTemp8[lbModel.G8_DV_P_P_M][i]);

        sum        = lbModel.fTemp7[lbModel.G7_DV_P_M_P][i]+lbModel.fTemp8[lbModel.G8_DV_M_P_M][i];
        rho[i]    += sum;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    -= point25sum;
        Pyz[i]    -= point25sum;
        Pzx[i]    += point25sum;
        uX[i]     += 0.5*(lbModel.fTemp7[lbModel.G7_DV_P_M_P][i]-lbModel.fTemp8[lbModel.G8_DV_M_P_M ][i]);
        uY[i]     -= 0.5*(lbModel.fTemp7[lbModel.G7_DV_P_M_P][i]-lbModel.fTemp8[lbModel.G8_DV_M_P_M ][i]);
        uZ[i]     += 0.5*(lbModel.fTemp7[lbModel.G7_DV_P_M_P][i]-lbModel.fTemp8[lbModel.G8_DV_M_P_M ][i]);
    }

    // G9 and G10
    for(int i=0;i< VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp9[lbModel.G9_DV_P1_P1_P1][i]+lbModel.fTemp10[lbModel.G10_DV_M1_M1_M1][i];
        rho[i]   += sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   += sum;
        Pyz[i]   += sum;
        Pzx[i]   += sum;
        uX[i]    += (lbModel.fTemp9[lbModel.G9_DV_P1_P1_P1][i]-lbModel.fTemp10[lbModel.G10_DV_M1_M1_M1][i]);
        uY[i]    += (lbModel.fTemp9[lbModel.G9_DV_P1_P1_P1][i]-lbModel.fTemp10[lbModel.G10_DV_M1_M1_M1][i]);
        uZ[i]    += (lbModel.fTemp9[lbModel.G9_DV_P1_P1_P1][i]-lbModel.fTemp10[lbModel.G10_DV_M1_M1_M1][i]);

        sum       = lbModel.fTemp9[lbModel.G9_DV_M1_P1_P1][i]+lbModel.fTemp10[lbModel.G10_DV_P1_M1_M1][i];
        rho[i]   += sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   -= sum;
        Pyz[i]   += sum;
        Pzx[i]   -= sum;
        uX[i]    -= (lbModel.fTemp9[lbModel.G9_DV_M1_P1_P1][i]-lbModel.fTemp10[lbModel.G10_DV_P1_M1_M1][i]);
        uY[i]    += (lbModel.fTemp9[lbModel.G9_DV_M1_P1_P1][i]-lbModel.fTemp10[lbModel.G10_DV_P1_M1_M1][i]);
        uZ[i]    += (lbModel.fTemp9[lbModel.G9_DV_M1_P1_P1][i]-lbModel.fTemp10[lbModel.G10_DV_P1_M1_M1][i]);


        sum       = lbModel.fTemp9[lbModel.G9_DV_M1_M1_P1][i]+lbModel.fTemp10[lbModel.G10_DV_P1_P1_M1][i];
        rho[i]   += sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   += sum;
        Pyz[i]   -= sum;
        Pzx[i]   -= sum;
        uX[i]    -= ( lbModel.fTemp9[lbModel.G9_DV_M1_M1_P1][i]-lbModel.fTemp10[lbModel.G10_DV_P1_P1_M1][i]);
        uY[i]    -= ( lbModel.fTemp9[lbModel.G9_DV_M1_M1_P1][i]-lbModel.fTemp10[lbModel.G10_DV_P1_P1_M1][i]);
        uZ[i]    += ( lbModel.fTemp9[lbModel.G9_DV_M1_M1_P1][i]-lbModel.fTemp10[lbModel.G10_DV_P1_P1_M1][i]);

        sum       = lbModel.fTemp9[lbModel.G9_DV_P1_M1_P1 ][i]+lbModel.fTemp10[lbModel.G10_DV_M1_P1_M1][i];
        rho[i]   += sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   -= sum;
        Pyz[i]   -= sum;
        Pzx[i]   += sum;
        uX[i]    += (lbModel.fTemp9[lbModel.G9_DV_P1_M1_P1][i]-lbModel.fTemp10[lbModel.G10_DV_M1_P1_M1][i]);
        uY[i]    -= (lbModel.fTemp9[lbModel.G9_DV_P1_M1_P1][i]-lbModel.fTemp10[lbModel.G10_DV_M1_P1_M1][i]);
        uZ[i]    += (lbModel.fTemp9[lbModel.G9_DV_P1_M1_P1][i]-lbModel.fTemp10[lbModel.G10_DV_M1_P1_M1][i]);
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
        dataType1  oneByRho = 1.0/rho[i];
        dataType1  oneByThree = 1.0/3.0;
        uX[i]    *= oneByRho;
        uY[i]    *= oneByRho;
        uZ[i]    *= oneByRho;
    }
}



template<typename dataType1>
void getGrad11Moment(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx)
{
    dataType1 sum(0.0),point25sum(0.0);

    for(int i=0;i<VECT_LENGTH;i++)
    {
        rho[i] = lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO];
    }

    // G1 and G2
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp1[i][lbModel.G1_DV_P2_ZERO_ZERO]+lbModel.fTemp2[i][lbModel.G2_DV_M2_ZERO_ZERO];
        rho[i]   += sum;
        theta[i]  = 4.0*sum;
        Pxx[i]    = 4.0*sum;
        uX[i]     = 2.0*(lbModel.fTemp1[i][lbModel.G1_DV_P2_ZERO_ZERO]-lbModel.fTemp2[i][lbModel.G2_DV_M2_ZERO_ZERO]);

        sum       = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_P2_ZERO]+lbModel.fTemp2[i][lbModel.G2_DV_ZERO_M2_ZERO];
        rho[i]   += sum;
        theta[i] += 4.0*sum;
        Pyy[i]    = 4.0*sum;
        uY[i]     = 2.0*(lbModel.fTemp1[i][lbModel.G1_DV_ZERO_P2_ZERO]-lbModel.fTemp2[i][lbModel.G2_DV_ZERO_M2_ZERO]);

        sum       = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P2]+lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M2];
        rho[i]   += sum;
        theta[i] += 4.0*sum;
        Pzz[i]    = 4.0*sum;
        uZ[i]     = 2.0*(lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P2]-lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M2]);

        sum       = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P1]+lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M1];
        rho[i]   += sum;
        theta[i] += sum;
        Pzz[i]   += sum;
        uZ[i]    += lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P1]-lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M1];
    }

    // G3 and G4
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]+lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pyz[i]    = sum;
        uY[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1]);

        sum       = lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1 ]+lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pyz[i]   -= sum;
        uY[i]    -= (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1]-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1]-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1]);

        sum       = lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1]+lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pxx[i]   += sum;
        Pzz[i]   += sum;
        Pzx[i]    = sum;
        uX[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1]-lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1]-lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1]);

        sum       = lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1]+lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1];
        theta[i] += 2.0*sum;
        rho[i]   += sum;
        Pxx[i]   += sum;
        Pzz[i]   += sum;
        Pzx[i]   -= sum;
        uX[i]    -= (lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1]-lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1]-lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1]);
    }

    // G5 and G6
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO]+lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pxy[i]    = sum;
        uX[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO]-lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO]);
        uY[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO]-lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO]);

        sum       = lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO]+lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pxy[i]   -= sum;
        uX[i]    += (lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO]-lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO]);
        uY[i]    -= (lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO]-lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO]);

        sum      = lbModel.fTemp5[i][lbModel.G5_DV_ZERO_P1_ZERO]+lbModel.fTemp6[i][lbModel.G6_DV_ZERO_M1_ZERO];
        rho[i]   += sum;
        theta[i] += sum;
        Pyy[i]   += sum;
        uY[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_ZERO_P1_ZERO]-lbModel.fTemp6[i][lbModel.G6_DV_ZERO_M1_ZERO]);

        sum       = lbModel.fTemp5[i][lbModel.G5_DV_P1_ZERO_ZERO]+lbModel.fTemp6[i][lbModel.G6_DV_M1_ZERO_ZERO];
        rho[i]   += sum;
        theta[i] += sum;
        Pxx[i]   += sum;
        uX[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_P1_ZERO_ZERO]-lbModel.fTemp6[i][lbModel.G6_DV_M1_ZERO_ZERO]);
    }

    // G7 and G8
    const dataType1 threeByFour=3.0/4.0;
    for(int i=0;i< VECT_LENGTH;i++)
    {
        sum        = lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]+lbModel.fTemp8[i][lbModel.G8_DV_M_M_M];
        rho[i]    += sum;
        point25sum = 0.25*sum;
        theta[i]  += threeByFour*sum ;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    += point25sum;
        Pyz[i]    += point25sum;
        Pzx[i]    += point25sum;
        uX[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]);
        uY[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]);
        uZ[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]);

        sum        = lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]+lbModel.fTemp8[i][lbModel.G8_DV_P_M_M];
        rho[i]    += sum;
        theta[i]  += threeByFour*sum ;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    -= point25sum;
        Pyz[i]    += point25sum;
        Pzx[i]    -= point25sum;
        uX[i]     -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]);
        uY[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]);
        uZ[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]);


        sum        = lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]+lbModel.fTemp8[i][lbModel.G8_DV_P_P_M];
        rho[i]    += sum;
        theta[i]  += threeByFour*sum ;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    += point25sum;
        Pyz[i]    -= point25sum;
        Pzx[i]    -= point25sum;
        uX[i]     -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]);
        uY[i]     -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]);
        uZ[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]);

        sum        = lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]+lbModel.fTemp8[i][lbModel.G8_DV_M_P_M];
        rho[i]    += sum;
        theta[i] += threeByFour*sum ;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    -= point25sum;
        Pyz[i]    -= point25sum;
        Pzx[i]    += point25sum;
        uX[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_P_M ]);
        uY[i]     -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_P_M ]);
        uZ[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_P_M ]);
    }

    // G9 and G10
    for(int i=0;i< VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]+lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   += sum;
        Pyz[i]   += sum;
        Pzx[i]   += sum;
        uX[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]);
        uY[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]);

        sum       = lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]+lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   -= sum;
        Pyz[i]   += sum;
        Pzx[i]   -= sum;
        uX[i]    -= (lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]);
        uY[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]);


        sum       = lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]+lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   += sum;
        Pyz[i]   -= sum;
        Pzx[i]   -= sum;
        uX[i]    -= ( lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]);
        uY[i]    -= ( lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]);
        uZ[i]    += ( lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]);

        sum       = lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1 ]+lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   -= sum;
        Pyz[i]   -= sum;
        Pzx[i]   += sum;
        uX[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]);
        uY[i]    -= (lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]);
    }
    dataType1  oneByThree = 1.0/3.0;
    for(int i=0;i<VECT_LENGTH;i++)
    {
        dataType1  oneByRho = 1.0/rho[i];
        uX[i]    *= oneByRho;
        uY[i]    *= oneByRho;
        uZ[i]    *= oneByRho;
        theta[i] = oneByThree * (theta[i] - rho[i]*(uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i]));
        theta[i] *= oneByRho;
    }
}

template<typename dataType1>
void getGrad11Moment_heatAlpha(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx, dataType1 *q1, dataType1 *q2, dataType1 *q3)
{
    dataType1 sum(0.0),point25sum(0.0);

    for(int i=0;i<VECT_LENGTH;i++)
    {
     q1[i] = 0.0;
     q2[i] = 0.0;
     q3[i] = 0.0;
    }


    for(int i=0;i<VECT_LENGTH;i++)
    {
        rho[i] = lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO];
    }

    // G1 and G2
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp1[i][lbModel.G1_DV_P2_ZERO_ZERO]+lbModel.fTemp2[i][lbModel.G2_DV_M2_ZERO_ZERO];
        rho[i]   += sum;
        theta[i]  = 4.0*sum;
        Pxx[i]    = 4.0*sum;
        uX[i]     = 2.0*(lbModel.fTemp1[i][lbModel.G1_DV_P2_ZERO_ZERO]-lbModel.fTemp2[i][lbModel.G2_DV_M2_ZERO_ZERO]);

        sum       = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_P2_ZERO]+lbModel.fTemp2[i][lbModel.G2_DV_ZERO_M2_ZERO];
        rho[i]   += sum;
        theta[i] += 4.0*sum;
        Pyy[i]    = 4.0*sum;
        uY[i]     = 2.0*(lbModel.fTemp1[i][lbModel.G1_DV_ZERO_P2_ZERO]-lbModel.fTemp2[i][lbModel.G2_DV_ZERO_M2_ZERO]);

        sum       = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P2]+lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M2];
        rho[i]   += sum;
        theta[i] += 4.0*sum;
        Pzz[i]    = 4.0*sum;
        uZ[i]     = 2.0*(lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P2]-lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M2]);

        sum       = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P1]+lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M1];
        rho[i]   += sum;
        theta[i] += sum;
        Pzz[i]   += sum;
        uZ[i]    += lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P1]-lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M1];





    }

    // G3 and G4
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]+lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pyz[i]    = sum;
        uY[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1]);

        sum       = lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1 ]+lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pyz[i]   -= sum;
        uY[i]    -= (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1]-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1]-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1]);

        sum       = lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1]+lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pxx[i]   += sum;
        Pzz[i]   += sum;
        Pzx[i]    = sum;
        uX[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1]-lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1]-lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1]);

        sum       = lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1]+lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1];
        theta[i] += 2.0*sum;
        rho[i]   += sum;
        Pxx[i]   += sum;
        Pzz[i]   += sum;
        Pzx[i]   -= sum;
        uX[i]    -= (lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1]-lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1]-lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1]);
    }

    // G5 and G6
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO]+lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pxy[i]    = sum;
        uX[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO]-lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO]);
        uY[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO]-lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO]);

        sum       = lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO]+lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pxy[i]   -= sum;
        uX[i]    += (lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO]-lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO]);
        uY[i]    -= (lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO]-lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO]);

        sum      = lbModel.fTemp5[i][lbModel.G5_DV_ZERO_P1_ZERO]+lbModel.fTemp6[i][lbModel.G6_DV_ZERO_M1_ZERO];
        rho[i]   += sum;
        theta[i] += sum;
        Pyy[i]   += sum;
        uY[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_ZERO_P1_ZERO]-lbModel.fTemp6[i][lbModel.G6_DV_ZERO_M1_ZERO]);

        sum       = lbModel.fTemp5[i][lbModel.G5_DV_P1_ZERO_ZERO]+lbModel.fTemp6[i][lbModel.G6_DV_M1_ZERO_ZERO];
        rho[i]   += sum;
        theta[i] += sum;
        Pxx[i]   += sum;
        uX[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_P1_ZERO_ZERO]-lbModel.fTemp6[i][lbModel.G6_DV_M1_ZERO_ZERO]);
    }

    // G7 and G8
    const dataType1 threeByFour=3.0/4.0;
    for(int i=0;i< VECT_LENGTH;i++)
    {
        sum        = lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]+lbModel.fTemp8[i][lbModel.G8_DV_M_M_M];
        rho[i]    += sum;
        point25sum = 0.25*sum;
        theta[i]  += threeByFour*sum ;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    += point25sum;
        Pyz[i]    += point25sum;
        Pzx[i]    += point25sum;
        uX[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]);
        uY[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]);
        uZ[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]);

        sum        = lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]+lbModel.fTemp8[i][lbModel.G8_DV_P_M_M];
        rho[i]    += sum;
        theta[i]  += threeByFour*sum ;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    -= point25sum;
        Pyz[i]    += point25sum;
        Pzx[i]    -= point25sum;
        uX[i]     -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]);
        uY[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]);
        uZ[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]);


        sum        = lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]+lbModel.fTemp8[i][lbModel.G8_DV_P_P_M];
        rho[i]    += sum;
        theta[i]  += threeByFour*sum ;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    += point25sum;
        Pyz[i]    -= point25sum;
        Pzx[i]    -= point25sum;
        uX[i]     -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]);
        uY[i]     -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]);
        uZ[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]);

        sum        = lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]+lbModel.fTemp8[i][lbModel.G8_DV_M_P_M];
        rho[i]    += sum;
        theta[i] += threeByFour*sum ;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    -= point25sum;
        Pyz[i]    -= point25sum;
        Pzx[i]    += point25sum;
        uX[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_P_M ]);
        uY[i]     -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_P_M ]);
        uZ[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]-lbModel.fTemp8[i][lbModel.G8_DV_M_P_M ]);
    }

    // G9 and G10
    for(int i=0;i< VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]+lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   += sum;
        Pyz[i]   += sum;
        Pzx[i]   += sum;
        uX[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]);
        uY[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]);

        sum       = lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]+lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   -= sum;
        Pyz[i]   += sum;
        Pzx[i]   -= sum;
        uX[i]    -= (lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]);
        uY[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]);


        sum       = lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]+lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   += sum;
        Pyz[i]   -= sum;
        Pzx[i]   -= sum;
        uX[i]    -= ( lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]);
        uY[i]    -= ( lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]);
        uZ[i]    += ( lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]);

        sum       = lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1 ]+lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   -= sum;
        Pyz[i]   -= sum;
        Pzx[i]   += sum;
        uX[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]);
        uY[i]    -= (lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]-lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]);
    }

    dataType1  oneByThree = 1.0/3.0;
    for(int i=0;i<VECT_LENGTH;i++)
    {
        dataType1  oneByRho = 1.0/rho[i];
        uX[i]    *= oneByRho;
        uY[i]    *= oneByRho;
        uZ[i]    *= oneByRho;
        theta[i] = oneByThree * (theta[i] - rho[i]*(uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i]));
        theta[i] *= oneByRho;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
      for(int dv=0;dv<4;dv++)
      {
        q1[i]    += lbModel.fTemp1[i][dv]*lbModel.cxcsq[1+dv];
        q2[i]    += lbModel.fTemp1[i][dv]*lbModel.cycsq[1+dv];
        q3[i]    += lbModel.fTemp1[i][dv]*lbModel.czcsq[1+dv];
      }

      for(int dv=0;dv<4;dv++)
      {
        q1[i]    += lbModel.fTemp2[i][dv]*lbModel.cxcsq[5+dv];
        q2[i]    += lbModel.fTemp2[i][dv]*lbModel.cycsq[5+dv];
        q3[i]    += lbModel.fTemp2[i][dv]*lbModel.czcsq[5+dv];
      }

      for(int dv=0;dv<4;dv++)
      {
        q1[i]    += lbModel.fTemp3[i][dv]*lbModel.cxcsq[9+dv];
        q2[i]    += lbModel.fTemp3[i][dv]*lbModel.cycsq[9+dv];
        q3[i]    += lbModel.fTemp3[i][dv]*lbModel.czcsq[9+dv];
      }
      for(int dv=0;dv<4;dv++)
      {
        q1[i]    += lbModel.fTemp4[i][dv]*lbModel.cxcsq[13+dv];
        q2[i]    += lbModel.fTemp4[i][dv]*lbModel.cycsq[13+dv];
        q3[i]    += lbModel.fTemp4[i][dv]*lbModel.czcsq[13+dv];
      }
      for(int dv=0;dv<4;dv++)
      {
        q1[i]    += lbModel.fTemp5[i][dv]*lbModel.cxcsq[17+dv];
        q2[i]    += lbModel.fTemp5[i][dv]*lbModel.cycsq[17+dv];
        q3[i]    += lbModel.fTemp5[i][dv]*lbModel.czcsq[17+dv];
      }
      for(int dv=0;dv<4;dv++)
      {
        q1[i]    += lbModel.fTemp6[i][dv]*lbModel.cxcsq[21+dv];
        q2[i]    += lbModel.fTemp6[i][dv]*lbModel.cycsq[21+dv];
        q3[i]    += lbModel.fTemp6[i][dv]*lbModel.czcsq[21+dv];
      }
      for(int dv=0;dv<4;dv++)
      {
        q1[i]    += lbModel.fTemp7[i][dv]*lbModel.cxcsq[25+dv];
        q2[i]    += lbModel.fTemp7[i][dv]*lbModel.cycsq[25+dv];
        q3[i]    += lbModel.fTemp7[i][dv]*lbModel.czcsq[25+dv];
      }
      for(int dv=0;dv<4;dv++)
      {
        q1[i]    += lbModel.fTemp8[i][dv]*lbModel.cxcsq[29+dv];
        q2[i]    += lbModel.fTemp8[i][dv]*lbModel.cycsq[29+dv];
        q3[i]    += lbModel.fTemp8[i][dv]*lbModel.czcsq[29+dv];
      }
      for(int dv=0;dv<4;dv++)
      {
        q1[i]    += lbModel.fTemp9[i][dv]*lbModel.cxcsq[33+dv];
        q2[i]    += lbModel.fTemp9[i][dv]*lbModel.cycsq[33+dv];
        q3[i]    += lbModel.fTemp9[i][dv]*lbModel.czcsq[33+dv];
      }
      for(int dv=0;dv<4;dv++)
      {
        q1[i]    += lbModel.fTemp10[i][dv]*lbModel.cxcsq[37+dv];
        q2[i]    += lbModel.fTemp10[i][dv]*lbModel.cycsq[37+dv];
        q3[i]    += lbModel.fTemp10[i][dv]*lbModel.czcsq[37+dv];
      }
    }

}

template<typename dataType1>
void getGrad11Moment_heatAlpha_YplusFlux(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx, dataType1 *q1, dataType1 *q2, dataType1 *q3)
{
    dataType1 sum(0.0),point25sum(0.0);

    for(int i=0;i<VECT_LENGTH;i++)
    {
     q1[i] = 0.0;
     q2[i] = 0.0;
     q3[i] = 0.0;
    }


    for(int i=0;i<VECT_LENGTH;i++)
    {
        rho[i] = 0.;
    }

    // G1 and G2
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = 0.;
        rho[i]   += sum;
        theta[i]  = 4.0*sum;
        Pxx[i]    = 4.0*sum;
        uX[i]     = 2.0*0.;

        sum       = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_P2_ZERO];
        rho[i]   += sum;
        theta[i] += 4.0*sum;
        Pyy[i]    = 4.0*sum;
        uY[i]     = 2.0*(lbModel.fTemp1[i][lbModel.G1_DV_ZERO_P2_ZERO]);

        sum       = 0.;
        rho[i]   += sum;
        theta[i] += 4.0*sum;
        Pzz[i]    = 4.0*sum;
        uZ[i]     = 2.0*0.;

        sum       = 0.;
        rho[i]   += sum;
        theta[i] += sum;
        Pzz[i]   += sum;
        uZ[i]    += 0.;





    }

    // G3 and G4
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pyz[i]    = sum;
        uY[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1]);

        sum       = lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pyz[i]   -= sum;
        uY[i]    -= (-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1]);
        uZ[i]    += (-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1]);

        sum       = 0.;
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pxx[i]   += sum;
        Pzz[i]   += sum;
        Pzx[i]    = sum;
        uX[i]    += 0.;
        uZ[i]    += 0.;

        sum       = 0.;
        theta[i] += 2.0*sum;
        rho[i]   += sum;
        Pxx[i]   += sum;
        Pzz[i]   += sum;
        Pzx[i]   -= sum;
        uX[i]    -= 0.;
        uZ[i]    += 0.;
    }

    // G5 and G6
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pxy[i]    = sum;
        uX[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO]);
        uY[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO]);

        sum       = lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pxy[i]   -= sum;
        uX[i]    += (-lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO]);
        uY[i]    -= (-lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO]);

        sum      = lbModel.fTemp5[i][lbModel.G5_DV_ZERO_P1_ZERO];
        rho[i]   += sum;
        theta[i] += sum;
        Pyy[i]   += sum;
        uY[i]    += (lbModel.fTemp5[i][lbModel.G5_DV_ZERO_P1_ZERO]);

        sum       = 0.;
        rho[i]   += sum;
        theta[i] += sum;
        Pxx[i]   += sum;
        uX[i]    += 0.;
    }

    // G7 and G8
    const dataType1 threeByFour=3.0/4.0;
    for(int i=0;i< VECT_LENGTH;i++)
    {
        sum        = lbModel.fTemp7[i][lbModel.G7_DV_P_P_P];
        rho[i]    += sum;
        point25sum = 0.25*sum;
        theta[i]  += threeByFour*sum ;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    += point25sum;
        Pyz[i]    += point25sum;
        Pzx[i]    += point25sum;
        uX[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]);
        uY[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]);
        uZ[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_P_P]);

        sum        = lbModel.fTemp7[i][lbModel.G7_DV_M_P_P];
        rho[i]    += sum;
        theta[i]  += threeByFour*sum ;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    -= point25sum;
        Pyz[i]    += point25sum;
        Pzx[i]    -= point25sum;
        uX[i]     -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]);
        uY[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]);
        uZ[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_P_P]);


        sum        = lbModel.fTemp8[i][lbModel.G8_DV_P_P_M];
        rho[i]    += sum;
        theta[i]  += threeByFour*sum ;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    += point25sum;
        Pyz[i]    -= point25sum;
        Pzx[i]    -= point25sum;
        uX[i]     -= 0.5*(-lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]);
        uY[i]     -= 0.5*(-lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]);
        uZ[i]     += 0.5*(-lbModel.fTemp8[i][lbModel.G8_DV_P_P_M]);

        sum        = lbModel.fTemp8[i][lbModel.G8_DV_M_P_M];
        rho[i]    += sum;
        theta[i] += threeByFour*sum ;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    -= point25sum;
        Pyz[i]    -= point25sum;
        Pzx[i]    += point25sum;
        uX[i]     += 0.5*(-lbModel.fTemp8[i][lbModel.G8_DV_M_P_M]);
        uY[i]     -= 0.5*(-lbModel.fTemp8[i][lbModel.G8_DV_M_P_M]);
        uZ[i]     += 0.5*(-lbModel.fTemp8[i][lbModel.G8_DV_M_P_M]);
    }

    // G9 and G10
    for(int i=0;i< VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   += sum;
        Pyz[i]   += sum;
        Pzx[i]   += sum;
        uX[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]);
        uY[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1]);

        sum       = lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   -= sum;
        Pyz[i]   += sum;
        Pzx[i]   -= sum;
        uX[i]    -= (lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]);
        uY[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1]);


        sum       = lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   += sum;
        Pyz[i]   -= sum;
        Pzx[i]   -= sum;
        uX[i]    -= (-lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]);
        uY[i]    -= (-lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]);
        uZ[i]    += (-lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1]);

        sum       = lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   -= sum;
        Pyz[i]   -= sum;
        Pzx[i]   += sum;
        uX[i]    += (-lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]);
        uY[i]    -= (-lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]);
        uZ[i]    += (-lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1]);
    }

    dataType1  oneByThree = 1.0/3.0;
    for(int i=0;i<VECT_LENGTH;i++)
    {
        dataType1  oneByRho = 1.0/rho[i];
        uX[i]    *= oneByRho;
        uY[i]    *= oneByRho;
        uZ[i]    *= oneByRho;
        theta[i] = oneByThree * (theta[i] - rho[i]*(uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i]));
        theta[i] *= oneByRho;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
      int dvYminus = lbModel.G1_DV_ZERO_P2_ZERO;
      q1[i]    += lbModel.fTemp1[i][dvYminus]*lbModel.cxcsq[1+dvYminus];
      q2[i]    += lbModel.fTemp1[i][dvYminus]*lbModel.cycsq[1+dvYminus];
      q3[i]    += lbModel.fTemp1[i][dvYminus]*lbModel.czcsq[1+dvYminus];

      dvYminus = lbModel.G3_DV_ZERO_P1_P1;
      q1[i]    += lbModel.fTemp3[i][dvYminus]*lbModel.cxcsq[9+dvYminus];
      q2[i]    += lbModel.fTemp3[i][dvYminus]*lbModel.cycsq[9+dvYminus];
      q3[i]    += lbModel.fTemp3[i][dvYminus]*lbModel.czcsq[9+dvYminus];

      dvYminus = lbModel.G4_DV_ZERO_P1_M1;
      q1[i]    += lbModel.fTemp4[i][dvYminus]*lbModel.cxcsq[13+dvYminus];
      q2[i]    += lbModel.fTemp4[i][dvYminus]*lbModel.cycsq[13+dvYminus];
      q3[i]    += lbModel.fTemp4[i][dvYminus]*lbModel.czcsq[13+dvYminus];

      int dvYminus1 = lbModel.G5_DV_P1_P1_ZERO;
      q1[i]    += lbModel.fTemp5[i][dvYminus1]*lbModel.cxcsq[17+dvYminus1];
      q2[i]    += lbModel.fTemp5[i][dvYminus1]*lbModel.cycsq[17+dvYminus1];
      q3[i]    += lbModel.fTemp5[i][dvYminus1]*lbModel.czcsq[17+dvYminus1];

      int dvYminus2 = lbModel.G5_DV_M1_P1_ZERO;
      q1[i]    += lbModel.fTemp5[i][dvYminus2]*lbModel.cxcsq[17+dvYminus2];
      q2[i]    += lbModel.fTemp5[i][dvYminus2]*lbModel.cycsq[17+dvYminus2];
      q3[i]    += lbModel.fTemp5[i][dvYminus2]*lbModel.czcsq[17+dvYminus2];

      int dvYminus3 = lbModel.G5_DV_ZERO_P1_ZERO;
      q1[i]    += lbModel.fTemp5[i][dvYminus3]*lbModel.cxcsq[17+dvYminus3];
      q2[i]    += lbModel.fTemp5[i][dvYminus3]*lbModel.cycsq[17+dvYminus3];
      q3[i]    += lbModel.fTemp5[i][dvYminus3]*lbModel.czcsq[17+dvYminus3];

      dvYminus1 = lbModel.G7_DV_P_P_P;
      q1[i]    += lbModel.fTemp7[i][dvYminus1]*lbModel.cxcsq[25+dvYminus1];
      q2[i]    += lbModel.fTemp7[i][dvYminus1]*lbModel.cycsq[25+dvYminus1];
      q3[i]    += lbModel.fTemp7[i][dvYminus1]*lbModel.czcsq[25+dvYminus1];

      dvYminus2 = lbModel.G7_DV_M_P_P;
      q1[i]    += lbModel.fTemp7[i][dvYminus2]*lbModel.cxcsq[25+dvYminus2];
      q2[i]    += lbModel.fTemp7[i][dvYminus2]*lbModel.cycsq[25+dvYminus2];
      q3[i]    += lbModel.fTemp7[i][dvYminus2]*lbModel.czcsq[25+dvYminus2];

      dvYminus1 = lbModel.G8_DV_P_P_M;
      q1[i]    += lbModel.fTemp8[i][dvYminus1]*lbModel.cxcsq[29+dvYminus1];
      q2[i]    += lbModel.fTemp8[i][dvYminus1]*lbModel.cycsq[29+dvYminus1];
      q3[i]    += lbModel.fTemp8[i][dvYminus1]*lbModel.czcsq[29+dvYminus1];

      dvYminus2 = lbModel.G8_DV_M_P_M;
      q1[i]    += lbModel.fTemp8[i][dvYminus2]*lbModel.cxcsq[29+dvYminus2];
      q2[i]    += lbModel.fTemp8[i][dvYminus2]*lbModel.cycsq[29+dvYminus2];
      q3[i]    += lbModel.fTemp8[i][dvYminus2]*lbModel.czcsq[29+dvYminus2];

      dvYminus1 = lbModel.G9_DV_P1_P1_P1;
      q1[i]    += lbModel.fTemp9[i][dvYminus1]*lbModel.cxcsq[33+dvYminus1];
      q2[i]    += lbModel.fTemp9[i][dvYminus1]*lbModel.cycsq[33+dvYminus1];
      q3[i]    += lbModel.fTemp9[i][dvYminus1]*lbModel.czcsq[33+dvYminus1];

      dvYminus2 = lbModel.G9_DV_M1_P1_P1;
      q1[i]    += lbModel.fTemp9[i][dvYminus2]*lbModel.cxcsq[33+dvYminus2];
      q2[i]    += lbModel.fTemp9[i][dvYminus2]*lbModel.cycsq[33+dvYminus2];
      q3[i]    += lbModel.fTemp9[i][dvYminus2]*lbModel.czcsq[33+dvYminus2];

      dvYminus1 = lbModel.G10_DV_P1_P1_M1;
      q1[i]    += lbModel.fTemp10[i][dvYminus1]*lbModel.cxcsq[37+dvYminus1];
      q2[i]    += lbModel.fTemp10[i][dvYminus1]*lbModel.cycsq[37+dvYminus1];
      q3[i]    += lbModel.fTemp10[i][dvYminus1]*lbModel.czcsq[37+dvYminus1];

      dvYminus2 = lbModel.G10_DV_M1_P1_M1;
      q1[i]    += lbModel.fTemp10[i][dvYminus2]*lbModel.cxcsq[37+dvYminus2];
      q2[i]    += lbModel.fTemp10[i][dvYminus2]*lbModel.cycsq[37+dvYminus2];
      q3[i]    += lbModel.fTemp10[i][dvYminus2]*lbModel.czcsq[37+dvYminus2];
    }

}


template<typename dataType1>
void getGrad11Moment_heatAlpha_YminusFlux(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx, dataType1 *q1, dataType1 *q2, dataType1 *q3)
{
    dataType1 sum(0.0),point25sum(0.0);

    for(int i=0;i<VECT_LENGTH;i++)
    {
     q1[i] = 0.0;
     q2[i] = 0.0;
     q3[i] = 0.0;
    }


    for(int i=0;i<VECT_LENGTH;i++)
    {
        rho[i] = 0.;
    }

    // G1 and G2
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = 0.;
        rho[i]   += sum;
        theta[i]  = 4.0*sum;
        Pxx[i]    = 4.0*sum;
        uX[i]     = 2.0*0.;

        sum       = lbModel.fTemp2[i][lbModel.G2_DV_ZERO_M2_ZERO];
        rho[i]   += sum;
        theta[i] += 4.0*sum;
        Pyy[i]    = 4.0*sum;
        uY[i]     = 2.0*(-lbModel.fTemp2[i][lbModel.G2_DV_ZERO_M2_ZERO]);

        sum       = 0.;
        rho[i]   += sum;
        theta[i] += 4.0*sum;
        Pzz[i]    = 4.0*sum;
        uZ[i]     = 2.0*0.;

        sum       = 0.;
        rho[i]   += sum;
        theta[i] += sum;
        Pzz[i]   += sum;
        uZ[i]    += 0.;





    }

    // G3 and G4
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pyz[i]    = sum;
        uY[i]    += (-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1]);
        uZ[i]    += (-lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1]);

        sum       = lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1 ];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pyz[i]   -= sum;
        uY[i]    -= (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1]);
        uZ[i]    += (lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1]);

        sum       = 0.;
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pxx[i]   += sum;
        Pzz[i]   += sum;
        Pzx[i]    = sum;
        uX[i]    += 0.;
        uZ[i]    += 0.;

        sum       = 0.;
        theta[i] += 2.0*sum;
        rho[i]   += sum;
        Pxx[i]   += sum;
        Pzz[i]   += sum;
        Pzx[i]   -= sum;
        uX[i]    -= 0.;
        uZ[i]    += 0.;
    }

    // G5 and G6
    for(int i=0;i<VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pxy[i]    = sum;
        uX[i]    += (-lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO]);
        uY[i]    += (-lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO]);

        sum       = lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO];
        rho[i]   += sum;
        theta[i] += 2.0*sum;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pxy[i]   -= sum;
        uX[i]    += (lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO]);
        uY[i]    -= (lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO]);

        sum       = lbModel.fTemp6[i][lbModel.G6_DV_ZERO_M1_ZERO];
        rho[i]   += sum;
        theta[i] += sum;
        Pyy[i]   += sum;
        uY[i]    += (-lbModel.fTemp6[i][lbModel.G6_DV_ZERO_M1_ZERO]);

        sum       = 0.;
        rho[i]   += sum;
        theta[i] += sum;
        Pxx[i]   += sum;
        uX[i]    += 0.;
    }

    // G7 and G8
    const dataType1 threeByFour=3.0/4.0;
    for(int i=0;i< VECT_LENGTH;i++)
    {
        sum        = lbModel.fTemp8[i][lbModel.G8_DV_M_M_M];
        rho[i]    += sum;
        point25sum = 0.25*sum;
        theta[i]  += threeByFour*sum ;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    += point25sum;
        Pyz[i]    += point25sum;
        Pzx[i]    += point25sum;
        uX[i]     += 0.5*(-lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]);
        uY[i]     += 0.5*(-lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]);
        uZ[i]     += 0.5*(-lbModel.fTemp8[i][lbModel.G8_DV_M_M_M]);

        sum        = lbModel.fTemp8[i][lbModel.G8_DV_P_M_M];
        rho[i]    += sum;
        theta[i]  += threeByFour*sum ;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    -= point25sum;
        Pyz[i]    += point25sum;
        Pzx[i]    -= point25sum;
        uX[i]     -= 0.5*(-lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]);
        uY[i]     += 0.5*(-lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]);
        uZ[i]     += 0.5*(-lbModel.fTemp8[i][lbModel.G8_DV_P_M_M]);


        sum        = lbModel.fTemp7[i][lbModel.G7_DV_M_M_P];
        rho[i]    += sum;
        theta[i]  += threeByFour*sum ;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    += point25sum;
        Pyz[i]    -= point25sum;
        Pzx[i]    -= point25sum;
        uX[i]     -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]);
        uY[i]     -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]);
        uZ[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_M_M_P]);

        sum        = lbModel.fTemp7[i][lbModel.G7_DV_P_M_P];
        rho[i]    += sum;
        theta[i] += threeByFour*sum ;
        point25sum = 0.25*sum;
        Pxx[i]    += point25sum;
        Pyy[i]    += point25sum;
        Pzz[i]    += point25sum;
        Pxy[i]    -= point25sum;
        Pyz[i]    -= point25sum;
        Pzx[i]    += point25sum;
        uX[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]);
        uY[i]     -= 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]);
        uZ[i]     += 0.5*(lbModel.fTemp7[i][lbModel.G7_DV_P_M_P]);
    }

    // G9 and G10
    for(int i=0;i< VECT_LENGTH;i++)
    {
        sum       = lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   += sum;
        Pyz[i]   += sum;
        Pzx[i]   += sum;
        uX[i]    += (-lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]);
        uY[i]    += (-lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]);
        uZ[i]    += (-lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1]);

        sum       = lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   -= sum;
        Pyz[i]   += sum;
        Pzx[i]   -= sum;
        uX[i]    -= (-lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]);
        uY[i]    += (-lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]);
        uZ[i]    += (-lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1]);


        sum       = lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   += sum;
        Pyz[i]   -= sum;
        Pzx[i]   -= sum;
        uX[i]    -= (lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]);
        uY[i]    -= (lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1]);

        sum       = lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1];
        rho[i]   += sum;
        theta[i] += 3.0*sum ;
        Pxx[i]   += sum;
        Pyy[i]   += sum;
        Pzz[i]   += sum;
        Pxy[i]   -= sum;
        Pyz[i]   -= sum;
        Pzx[i]   += sum;
        uX[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]);
        uY[i]    -= (lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]);
        uZ[i]    += (lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1]);
    }

    dataType1  oneByThree = 1.0/3.0;
    for(int i=0;i<VECT_LENGTH;i++)
    {
        dataType1  oneByRho = 1.0/rho[i];
        uX[i]    *= oneByRho;
        uY[i]    *= oneByRho;
        uZ[i]    *= oneByRho;
        theta[i] = oneByThree * (theta[i] - rho[i]*(uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i]));
        theta[i] *= oneByRho;
    }

    for(int i=0;i<VECT_LENGTH;i++)
    {
      int dvYminus = lbModel.G2_DV_ZERO_M2_ZERO;
      q1[i]    += lbModel.fTemp2[i][dvYminus]*lbModel.cxcsq[5+dvYminus];
      q2[i]    += lbModel.fTemp2[i][dvYminus]*lbModel.cycsq[5+dvYminus];
      q3[i]    += lbModel.fTemp2[i][dvYminus]*lbModel.czcsq[5+dvYminus];

      dvYminus = lbModel.G3_DV_ZERO_M1_P1;
      q1[i]    += lbModel.fTemp3[i][dvYminus]*lbModel.cxcsq[9+dvYminus];
      q2[i]    += lbModel.fTemp3[i][dvYminus]*lbModel.cycsq[9+dvYminus];
      q3[i]    += lbModel.fTemp3[i][dvYminus]*lbModel.czcsq[9+dvYminus];

      dvYminus = lbModel.G4_DV_ZERO_M1_M1;
      q1[i]    += lbModel.fTemp4[i][dvYminus]*lbModel.cxcsq[13+dvYminus];
      q2[i]    += lbModel.fTemp4[i][dvYminus]*lbModel.cycsq[13+dvYminus];
      q3[i]    += lbModel.fTemp4[i][dvYminus]*lbModel.czcsq[13+dvYminus];

      int dvYminus1 = lbModel.G6_DV_M1_M1_ZERO;
      q1[i]    += lbModel.fTemp6[i][dvYminus1]*lbModel.cxcsq[21+dvYminus1];
      q2[i]    += lbModel.fTemp6[i][dvYminus1]*lbModel.cycsq[21+dvYminus1];
      q3[i]    += lbModel.fTemp6[i][dvYminus1]*lbModel.czcsq[21+dvYminus1];

      int dvYminus2 = lbModel.G6_DV_P1_M1_ZERO;
      q1[i]    += lbModel.fTemp6[i][dvYminus2]*lbModel.cxcsq[21+dvYminus2];
      q2[i]    += lbModel.fTemp6[i][dvYminus2]*lbModel.cycsq[21+dvYminus2];
      q3[i]    += lbModel.fTemp6[i][dvYminus2]*lbModel.czcsq[21+dvYminus2];

      int dvYminus3 = lbModel.G6_DV_ZERO_M1_ZERO;
      q1[i]    += lbModel.fTemp6[i][dvYminus3]*lbModel.cxcsq[21+dvYminus3];
      q2[i]    += lbModel.fTemp6[i][dvYminus3]*lbModel.cycsq[21+dvYminus3];
      q3[i]    += lbModel.fTemp6[i][dvYminus3]*lbModel.czcsq[21+dvYminus3];

      dvYminus1 = lbModel.G7_DV_M_M_P;
      q1[i]    += lbModel.fTemp7[i][dvYminus1]*lbModel.cxcsq[25+dvYminus1];
      q2[i]    += lbModel.fTemp7[i][dvYminus1]*lbModel.cycsq[25+dvYminus1];
      q3[i]    += lbModel.fTemp7[i][dvYminus1]*lbModel.czcsq[25+dvYminus1];

      dvYminus2 = lbModel.G7_DV_P_M_P;
      q1[i]    += lbModel.fTemp7[i][dvYminus2]*lbModel.cxcsq[25+dvYminus2];
      q2[i]    += lbModel.fTemp7[i][dvYminus2]*lbModel.cycsq[25+dvYminus2];
      q3[i]    += lbModel.fTemp7[i][dvYminus2]*lbModel.czcsq[25+dvYminus2];

      dvYminus1 = lbModel.G8_DV_M_M_M;
      q1[i]    += lbModel.fTemp8[i][dvYminus1]*lbModel.cxcsq[29+dvYminus1];
      q2[i]    += lbModel.fTemp8[i][dvYminus1]*lbModel.cycsq[29+dvYminus1];
      q3[i]    += lbModel.fTemp8[i][dvYminus1]*lbModel.czcsq[29+dvYminus1];

      dvYminus2 = lbModel.G8_DV_P_M_M;
      q1[i]    += lbModel.fTemp8[i][dvYminus2]*lbModel.cxcsq[29+dvYminus2];
      q2[i]    += lbModel.fTemp8[i][dvYminus2]*lbModel.cycsq[29+dvYminus2];
      q3[i]    += lbModel.fTemp8[i][dvYminus2]*lbModel.czcsq[29+dvYminus2];

      dvYminus1 = lbModel.G9_DV_M1_M1_P1;
      q1[i]    += lbModel.fTemp9[i][dvYminus1]*lbModel.cxcsq[33+dvYminus1];
      q2[i]    += lbModel.fTemp9[i][dvYminus1]*lbModel.cycsq[33+dvYminus1];
      q3[i]    += lbModel.fTemp9[i][dvYminus1]*lbModel.czcsq[33+dvYminus1];

      dvYminus2 = lbModel.G9_DV_P1_M1_P1;
      q1[i]    += lbModel.fTemp9[i][dvYminus2]*lbModel.cxcsq[33+dvYminus2];
      q2[i]    += lbModel.fTemp9[i][dvYminus2]*lbModel.cycsq[33+dvYminus2];
      q3[i]    += lbModel.fTemp9[i][dvYminus2]*lbModel.czcsq[33+dvYminus2];

      dvYminus1 = lbModel.G10_DV_M1_M1_M1;
      q1[i]    += lbModel.fTemp10[i][dvYminus1]*lbModel.cxcsq[37+dvYminus1];
      q2[i]    += lbModel.fTemp10[i][dvYminus1]*lbModel.cycsq[37+dvYminus1];
      q3[i]    += lbModel.fTemp10[i][dvYminus1]*lbModel.czcsq[37+dvYminus1];

      dvYminus2 = lbModel.G10_DV_P1_M1_M1;
      q1[i]    += lbModel.fTemp10[i][dvYminus2]*lbModel.cxcsq[37+dvYminus2];
      q2[i]    += lbModel.fTemp10[i][dvYminus2]*lbModel.cycsq[37+dvYminus2];
      q3[i]    += lbModel.fTemp10[i][dvYminus2]*lbModel.czcsq[37+dvYminus2];
    }

}

template <int N,int numblock, typename dataType1>
void calcLBTopFlux(lbmRD3Q41<dataType1> &lbModel41,int VECT_LENGTH, gridBCC3D<N, numblock, dataType1> &gridLBM41, int LBpointsX, int LBpointsY, int numMoments, int numSpaceAvgCells_LB_X, int numSpaceAvgCells_LB_Y, int yIndex, dataType1* momentsLBFlux, int myRank)
{
    for(int cellX = 0; cellX < LBpointsX; cellX++)
    {
        for(int cellY = 0; cellY < LBpointsY; cellY++)
        {
            int xIndex = cellY + gridLBM41.nB1;
            int zIndex = cellX + gridLBM41.nB3;

            int cellNumberLocal           = ((int)LBpointsX)*((int)(cellY/numSpaceAvgCells_LB_Y)) + ((int)(cellX/numSpaceAvgCells_LB_X));
            copyFromNode(lbModel41,gridLBM41,VECT_LENGTH,xIndex,yIndex,zIndex);
            getGrad11Moment_heatAlpha_YplusFlux(lbModel41,VECT_LENGTH,&momentsLBFlux[numMoments*cellNumberLocal + 0],&momentsLBFlux[numMoments*cellNumberLocal + 1],&momentsLBFlux[numMoments*cellNumberLocal + 2],&momentsLBFlux[numMoments*cellNumberLocal + 3],&momentsLBFlux[numMoments*cellNumberLocal + 4],&momentsLBFlux[numMoments*cellNumberLocal + 5],&momentsLBFlux[numMoments*cellNumberLocal + 6],&momentsLBFlux[numMoments*cellNumberLocal + 7],&momentsLBFlux[numMoments*cellNumberLocal + 8],&momentsLBFlux[numMoments*cellNumberLocal + 9],&momentsLBFlux[numMoments*cellNumberLocal + 10],&momentsLBFlux[numMoments*cellNumberLocal + 11],&momentsLBFlux[numMoments*cellNumberLocal + 12],&momentsLBFlux[numMoments*cellNumberLocal + 13]);

            if(myRank == 0)
                std::cout<<"cellNumberLocal: "<<cellNumberLocal<<std::endl;
                std::cout<<"fluxLB: "<<cellNumberLocal<<"\t"<<momentsLBFlux[numMoments*cellNumberLocal + 0]<<std::endl;
        }
    }
}

template <int N,int numblock, typename dataType1>
void calcLBBotFlux(lbmRD3Q41<dataType1> &lbModel41,int VECT_LENGTH, gridBCC3D<N, numblock, dataType1> &gridLBM41, int LBpointsX, int LBpointsY, int numMoments, int numSpaceAvgCells_LB_X, int numSpaceAvgCells_LB_Y, int yIndex, dataType1* momentsLBFlux)
{
    for(int cellX = 0; cellX < LBpointsX; cellX++)
    {
        for(int cellY = 0; cellY < LBpointsY; cellY++)
        {
            int xIndex = cellY + gridLBM41.nB1;
            int zIndex = cellX + gridLBM41.nB3;

            int cellNumberLocal           = ((int)LBpointsX)*((int)(cellY/numSpaceAvgCells_LB_Y)) + ((int)(cellX/numSpaceAvgCells_LB_X));
            copyFromNode(lbModel41,gridLBM41,VECT_LENGTH,xIndex,yIndex,zIndex);
            getGrad11Moment_heatAlpha_YminusFlux(lbModel41,VECT_LENGTH,&momentsLBFlux[numMoments*cellNumberLocal + 0],&momentsLBFlux[numMoments*cellNumberLocal + 1],&momentsLBFlux[numMoments*cellNumberLocal + 2],&momentsLBFlux[numMoments*cellNumberLocal + 3],&momentsLBFlux[numMoments*cellNumberLocal + 4],&momentsLBFlux[numMoments*cellNumberLocal + 5],&momentsLBFlux[numMoments*cellNumberLocal + 6],&momentsLBFlux[numMoments*cellNumberLocal + 7],&momentsLBFlux[numMoments*cellNumberLocal + 8],&momentsLBFlux[numMoments*cellNumberLocal + 9],&momentsLBFlux[numMoments*cellNumberLocal + 10],&momentsLBFlux[numMoments*cellNumberLocal + 11],&momentsLBFlux[numMoments*cellNumberLocal + 12],&momentsLBFlux[numMoments*cellNumberLocal + 13]);
        }
    }
}



template<typename dataType1>
void getGrad11MomentSIMD(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx)
{
    SIMD_REG _RHO, _UX, _UY, _UZ;              //4
    SIMD_REG _PXX,_PYY,_PZZ,_PXY;              //4
    SIMD_REG _PYZ,_PZX,_THETA,_twoReg;         //4
    SIMD_REG _fourReg,_temp1,_temp2,_temp4;    //4

    dataType1 tempArray[4] __attribute__ ((aligned(32)));
    dataType1  *pointer1;

    _twoReg  = SET1_PD(2.0);
    _fourReg = SET1_PD(4.0);

    _UX      =  SET1_PD(0.0);
    _UY      =  SET1_PD(0.0);
    _UZ      =  SET1_PD(0.0);
    _PXX     =  SET1_PD(0.0);
    _PYY     =  SET1_PD(0.0);
    _PZZ     =  SET1_PD(0.0);
    _PXY     =  SET1_PD(0.0);
    _PYZ     =  SET1_PD(0.0);
    _PZX     =  SET1_PD(0.0);
    _THETA   =  SET1_PD(0.0);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp0[i][lbModel.CENTER_DV_ZERO_ZERO_ZERO];
    pointer1 = &tempArray[0];

    _RHO = LOAD_PD(pointer1);


    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp1[i][lbModel.G1_DV_P2_ZERO_ZERO];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp2[i][lbModel.G2_DV_M2_ZERO_ZERO];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4 = ADD_PD(_temp1,_temp2);
    _RHO   = ADD_PD(_RHO,_temp4);
    _temp4 = MUL_PD(_temp4,_fourReg);
    _THETA = ADD_PD(_THETA,_temp4);
    _PXX   = ADD_PD(_PXX,_temp4);
    _temp4 = SUB_PD(_temp1,_temp2);
    _temp4 = MUL_PD(_temp4,_twoReg);
    _UX    = ADD_PD(_UX,_temp4);


    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_P2_ZERO];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp2[i][lbModel.G2_DV_ZERO_M2_ZERO];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4 = ADD_PD(_temp1,_temp2);
    _RHO   = ADD_PD(_RHO,_temp4);
    _temp4 = MUL_PD(_temp4,_fourReg);
    _THETA = ADD_PD(_THETA,_temp4);
    _PYY   = ADD_PD(_PYY,_temp4);
    _temp4 = SUB_PD(_temp1,_temp2);
    _temp4 = MUL_PD(_temp4,_twoReg);
    _UY    = ADD_PD(_UY,_temp4);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P2];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M2];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4 = ADD_PD(_temp1,_temp2);
    _RHO   = ADD_PD(_RHO,_temp4);
    _temp4 = MUL_PD(_temp4,_fourReg);
    _THETA = ADD_PD(_THETA,_temp4);
    _PZZ   = ADD_PD(_PZZ,_temp4);
    _temp4 = SUB_PD(_temp1,_temp2);
    _temp4 = MUL_PD(_temp4,_twoReg);
    _UZ    = ADD_PD(_UZ,_temp4);


    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp1[i][lbModel.G1_DV_ZERO_ZERO_P1];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp2[i][lbModel.G2_DV_ZERO_ZERO_M1];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4 = ADD_PD(_temp1,_temp2);
    _RHO   = ADD_PD(_RHO,_temp4);
    _THETA = ADD_PD(_THETA,_temp4);
    _PZZ   = ADD_PD(_PZZ,_temp4);
    _temp4 = SUB_PD(_temp1,_temp2);
    _UZ    = ADD_PD(_UZ,_temp4);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp3[i][lbModel.G3_DV_ZERO_P1_P1];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp4[i][lbModel.G4_DV_ZERO_M1_M1];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4 = ADD_PD(_temp1,_temp2);
    _RHO   = ADD_PD(_RHO,_temp4);
    _PYY   = ADD_PD(_PYY,_temp4);
    _PZZ   = ADD_PD(_PZZ,_temp4);
    _PYZ   = ADD_PD(_PYZ,_temp4);
    _temp4 = MUL_PD(_twoReg,_temp4);
    _THETA = ADD_PD(_THETA,_temp4);
    _temp4 = SUB_PD(_temp1,_temp2);
    _UY    = ADD_PD(_UY,_temp4);
    _UZ    = ADD_PD(_UZ,_temp4);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp3[i][lbModel.G3_DV_ZERO_M1_P1];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp4[i][lbModel.G4_DV_ZERO_P1_M1];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4 = ADD_PD(_temp1,_temp2);
    _RHO   = ADD_PD(_RHO,_temp4);
    _PYY   = ADD_PD(_PYY,_temp4);
    _PZZ   = ADD_PD(_PZZ,_temp4);
    _PYZ   = SUB_PD(_PYZ,_temp4);
    _temp4 = MUL_PD(_twoReg,_temp4);
    _THETA = ADD_PD(_THETA,_temp4);
    _temp4 = SUB_PD(_temp1,_temp2);
    _UY    = SUB_PD(_UY,_temp4);
    _UZ    = ADD_PD(_UZ,_temp4);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp3[i][lbModel.G3_DV_P1_ZERO_P1];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp4[i][lbModel.G4_DV_M1_ZERO_M1];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4 = ADD_PD(_temp1,_temp2);
    _RHO   = ADD_PD(_RHO,_temp4);
    _PXX   = ADD_PD(_PXX,_temp4);
    _PZZ   = ADD_PD(_PZZ,_temp4);
    _PZX   = ADD_PD(_PZX,_temp4);
    _temp4 = MUL_PD(_twoReg,_temp4);
    _THETA = ADD_PD(_THETA,_temp4);
    _temp4 = SUB_PD(_temp1,_temp2);
    _UX    = ADD_PD(_UX,_temp4);
    _UZ    = ADD_PD(_UZ,_temp4);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp3[i][lbModel.G3_DV_M1_ZERO_P1];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp4[i][lbModel.G4_DV_P1_ZERO_M1];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4 = ADD_PD(_temp1,_temp2);
    _RHO   = ADD_PD(_RHO,_temp4);
    _PXX   = ADD_PD(_PXX,_temp4);
    _PZZ   = ADD_PD(_PZZ,_temp4);
    _PZX   = SUB_PD(_PZX,_temp4);
    _temp4 = MUL_PD(_twoReg,_temp4);
    _THETA = ADD_PD(_THETA,_temp4);
    _temp4 = SUB_PD(_temp1,_temp2);
    _UX    = SUB_PD(_UX,_temp4);
    _UZ    = ADD_PD(_UZ,_temp4);


    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp5[i][lbModel.G5_DV_P1_P1_ZERO];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp6[i][lbModel.G6_DV_M1_M1_ZERO];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4 = ADD_PD(_temp1,_temp2);
    _RHO   = ADD_PD(_RHO,_temp4);
    _PXX   = ADD_PD(_PXX,_temp4);
    _PYY   = ADD_PD(_PYY,_temp4);
    _PXY   = ADD_PD(_PXY,_temp4);
    _temp4 = MUL_PD(_twoReg,_temp4);
    _THETA = ADD_PD(_THETA,_temp4);
    _temp4 = SUB_PD(_temp1,_temp2);
    _UX    = ADD_PD(_UX,_temp4);
    _UY    = ADD_PD(_UY,_temp4);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp5[i][lbModel.G5_DV_M1_P1_ZERO];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp6[i][lbModel.G6_DV_P1_M1_ZERO];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4 = ADD_PD(_temp1,_temp2);
    _RHO   = ADD_PD(_RHO,_temp4);
    _PXX   = ADD_PD(_PXX,_temp4);
    _PYY   = ADD_PD(_PYY,_temp4);
    _PXY   = SUB_PD(_PXY,_temp4);
    _temp4 = MUL_PD(_twoReg,_temp4);
    _THETA = ADD_PD(_THETA,_temp4);
    _temp4 = SUB_PD(_temp2,_temp1);
    _UX    = ADD_PD(_UX,_temp4);
    _UY    = SUB_PD(_UY,_temp4);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp5[i][lbModel.G5_DV_ZERO_P1_ZERO];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp6[i][lbModel.G6_DV_ZERO_M1_ZERO];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4 = ADD_PD(_temp1,_temp2);
    _RHO   = ADD_PD(_RHO,_temp4);
    _PYY   = ADD_PD(_PYY,_temp4);
    _THETA = ADD_PD(_THETA,_temp4);
    _temp4 = SUB_PD(_temp1,_temp2);
    _UY    = ADD_PD(_UY,_temp4);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp5[i][lbModel.G5_DV_P1_ZERO_ZERO];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp6[i][lbModel.G6_DV_M1_ZERO_ZERO];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4 = ADD_PD(_temp1,_temp2);
    _RHO   = ADD_PD(_RHO,_temp4);
    _PXX   = ADD_PD(_PXX,_temp4);
    _THETA = ADD_PD(_THETA,_temp4);
    _temp4 = SUB_PD(_temp1,_temp2);
    _UX    = ADD_PD(_UX,_temp4);



    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp7[i][lbModel.G7_DV_P_P_P];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp8[i][lbModel.G8_DV_M_M_M];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4   = ADD_PD(_temp1,_temp2);
    _RHO     = ADD_PD(_RHO,_temp4);
    _twoReg  = SET1_PD(0.25);
    _temp4   = MUL_PD(_twoReg,_temp4); // 0.25*SUM
    _PXX     = ADD_PD(_PXX,_temp4);
    _PYY     = ADD_PD(_PYY,_temp4);
    _PZZ     = ADD_PD(_PZZ,_temp4);
    _PXY     = ADD_PD(_PXY,_temp4);
    _PYZ     = ADD_PD(_PYZ,_temp4);
    _PZX     = ADD_PD(_PZX,_temp4);
    _twoReg  = SET1_PD(3.0);
    _temp4   = MUL_PD(_twoReg,_temp4); // 0.75*SUM
    _THETA   = ADD_PD(_THETA,_temp4);
    _temp4   = SUB_PD(_temp1,_temp2);
    _twoReg  = SET1_PD(0.5);
    _temp4   = MUL_PD(_twoReg,_temp4);
    _UX      = ADD_PD(_UX,_temp4);
    _UY      = ADD_PD(_UY,_temp4);
    _UZ      = ADD_PD(_UZ,_temp4);


    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp7[i][lbModel.G7_DV_M_P_P];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp8[i][lbModel.G8_DV_P_M_M];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4   = ADD_PD(_temp1,_temp2);
    _RHO     = ADD_PD(_RHO,_temp4);
    _twoReg  = SET1_PD(0.25);
    _temp4   = MUL_PD(_twoReg,_temp4); // 0.25*SUM
    _PXX     = ADD_PD(_PXX,_temp4);
    _PYY     = ADD_PD(_PYY,_temp4);
    _PZZ     = ADD_PD(_PZZ,_temp4);
    _PXY     = SUB_PD(_PXY,_temp4);
    _PYZ     = ADD_PD(_PYZ,_temp4);
    _PZX     = SUB_PD(_PZX,_temp4);
    _twoReg  = SET1_PD(3.0);
    _temp4   = MUL_PD(_twoReg,_temp4); // 0.75*SUM
    _THETA   = ADD_PD(_THETA,_temp4);
    _temp4   = SUB_PD(_temp1,_temp2);
    _twoReg  = SET1_PD(0.5);
    _temp4   = MUL_PD(_twoReg,_temp4);
    _UX      = SUB_PD(_UX,_temp4);
    _UY      = ADD_PD(_UY,_temp4);
    _UZ      = ADD_PD(_UZ,_temp4);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp7[i][lbModel.G7_DV_M_M_P];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp8[i][lbModel.G8_DV_P_P_M];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4   = ADD_PD(_temp1,_temp2);
    _RHO     = ADD_PD(_RHO,_temp4);
    _twoReg  = SET1_PD(0.25);
    _temp4   = MUL_PD(_twoReg,_temp4); // 0.25*SUM
    _PXX     = ADD_PD(_PXX,_temp4);
    _PYY     = ADD_PD(_PYY,_temp4);
    _PZZ     = ADD_PD(_PZZ,_temp4);
    _PXY     = ADD_PD(_PXY,_temp4);
    _PYZ     = SUB_PD(_PYZ,_temp4);
    _PZX     = SUB_PD(_PZX,_temp4);
    _twoReg  = SET1_PD(3.0);
    _temp4   = MUL_PD(_twoReg,_temp4); // 0.75*SUM
    _THETA   = ADD_PD(_THETA,_temp4);
    _temp4   = SUB_PD(_temp1,_temp2);
    _twoReg  = SET1_PD(0.5);
    _temp4   = MUL_PD(_twoReg,_temp4);
    _UX      = SUB_PD(_UX,_temp4);
    _UY      = SUB_PD(_UY,_temp4);
    _UZ      = ADD_PD(_UZ,_temp4);


    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp7[i][lbModel.G7_DV_P_M_P];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp8[i][lbModel.G8_DV_M_P_M];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4   = ADD_PD(_temp1,_temp2);
    _RHO     = ADD_PD(_RHO,_temp4);
    _twoReg  = SET1_PD(0.25);
    _temp4   = MUL_PD(_twoReg,_temp4); // 0.25*SUM
    _PXX     = ADD_PD(_PXX,_temp4);
    _PYY     = ADD_PD(_PYY,_temp4);
    _PZZ     = ADD_PD(_PZZ,_temp4);
    _PXY     = SUB_PD(_PXY,_temp4);
    _PYZ     = SUB_PD(_PYZ,_temp4);
    _PZX     = ADD_PD(_PZX,_temp4);
    _twoReg  = SET1_PD(3.0);
    _temp4   = MUL_PD(_twoReg,_temp4); // 0.75*SUM
    _THETA   = ADD_PD(_THETA,_temp4);
    _temp4   = SUB_PD(_temp1,_temp2);
    _twoReg  = SET1_PD(0.5);
    _temp4   = MUL_PD(_twoReg,_temp4);
    _UX      = ADD_PD(_UX,_temp4);
    _UY      = SUB_PD(_UY,_temp4);
    _UZ      = ADD_PD(_UZ,_temp4);


    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp9[i][lbModel.G9_DV_P1_P1_P1];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp10[i][lbModel.G10_DV_M1_M1_M1];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4   = ADD_PD(_temp1,_temp2);
    _RHO     = ADD_PD(_RHO,_temp4);
    _PXX     = ADD_PD(_PXX,_temp4);
    _PYY     = ADD_PD(_PYY,_temp4);
    _PZZ     = ADD_PD(_PZZ,_temp4);
    _PXY     = ADD_PD(_PXY,_temp4);
    _PYZ     = ADD_PD(_PYZ,_temp4);
    _PZX     = ADD_PD(_PZX,_temp4);
    _twoReg  = SET1_PD(3.0);
    _temp4   = MUL_PD(_twoReg,_temp4);
    _THETA   = ADD_PD(_THETA,_temp4);
    _temp4   = SUB_PD(_temp1,_temp2);
    _UX      = ADD_PD(_UX,_temp4);
    _UY      = ADD_PD(_UY,_temp4);
    _UZ      = ADD_PD(_UZ,_temp4);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp9[i][lbModel.G9_DV_M1_P1_P1];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp10[i][lbModel.G10_DV_P1_M1_M1];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4   = ADD_PD(_temp1,_temp2);
    _RHO     = ADD_PD(_RHO,_temp4);
    _PXX     = ADD_PD(_PXX,_temp4);
    _PYY     = ADD_PD(_PYY,_temp4);
    _PZZ     = ADD_PD(_PZZ,_temp4);
    _PXY     = SUB_PD(_PXY,_temp4);
    _PYZ     = ADD_PD(_PYZ,_temp4);
    _PZX     = SUB_PD(_PZX,_temp4);
    _temp4   = MUL_PD(_twoReg,_temp4);
    _THETA   = ADD_PD(_THETA,_temp4);
    _temp4   = SUB_PD(_temp1,_temp2);
    _UX      = SUB_PD(_UX,_temp4);
    _UY      = ADD_PD(_UY,_temp4);
    _UZ      = ADD_PD(_UZ,_temp4);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp9[i][lbModel.G9_DV_M1_M1_P1];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp10[i][lbModel.G10_DV_P1_P1_M1];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4   = ADD_PD(_temp1,_temp2);
    _RHO     = ADD_PD(_RHO,_temp4);
    _PXX     = ADD_PD(_PXX,_temp4);
    _PYY     = ADD_PD(_PYY,_temp4);
    _PZZ     = ADD_PD(_PZZ,_temp4);
    _PXY     = ADD_PD(_PXY,_temp4);
    _PYZ     = SUB_PD(_PYZ,_temp4);
    _PZX     = SUB_PD(_PZX,_temp4);
    _temp4   = MUL_PD(_twoReg,_temp4);
    _THETA   = ADD_PD(_THETA,_temp4);
    _temp4   = SUB_PD(_temp1,_temp2);
    _UX      = SUB_PD(_UX,_temp4);
    _UY      = SUB_PD(_UY,_temp4);
    _UZ      = ADD_PD(_UZ,_temp4);

    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp9[i][lbModel.G9_DV_P1_M1_P1];
    pointer1 = &tempArray[0];
    _temp1 = LOAD_PD(pointer1);
    for(int i=0;i<VECT_LENGTH;i++)
     tempArray[i] = lbModel.fTemp10[i][lbModel.G10_DV_M1_P1_M1];
    pointer1 = &tempArray[0];
    _temp2 = LOAD_PD(pointer1);

    _temp4   = ADD_PD(_temp1,_temp2);
    _RHO     = ADD_PD(_RHO,_temp4);
    _PXX     = ADD_PD(_PXX,_temp4);
    _PYY     = ADD_PD(_PYY,_temp4);
    _PZZ     = ADD_PD(_PZZ,_temp4);
    _PXY     = SUB_PD(_PXY,_temp4);
    _PYZ     = SUB_PD(_PYZ,_temp4);
    _PZX     = ADD_PD(_PZX,_temp4);
    _temp4   = MUL_PD(_twoReg,_temp4);
    _THETA   = ADD_PD(_THETA,_temp4);
    _temp4   = SUB_PD(_temp1,_temp2);
    _UX      = ADD_PD(_UX,_temp4);
    _UY      = SUB_PD(_UY,_temp4);
    _UZ      = ADD_PD(_UZ,_temp4);

    _twoReg  = SET1_PD(1.0/3.0);
    _fourReg = SET1_PD(1.0);
    _temp4   = DIV_PD(_fourReg,_RHO);
    _UX      = MUL_PD(_UX,_temp4);
    _UY      = MUL_PD(_UY,_temp4);
    _UZ      = MUL_PD(_UZ,_temp4);
    _temp1   = MUL_PD(_UX,_UX);
    _temp2   = MUL_PD(_UY,_UY);
    _temp1   = ADD_PD(_temp1,_temp2);
    _temp2   = MUL_PD(_UZ,_UZ);
    _temp1   = ADD_PD(_temp1,_temp2);
    _temp1   = MUL_PD(_RHO,_temp1);
    _temp1   = SUB_PD(_THETA,_temp1);
    _THETA   = MUL_PD(_twoReg,_temp1);
    _THETA   = MUL_PD(_temp4,_THETA);


    STORE_PD(rho  ,_RHO  );
    STORE_PD(uX   ,_UX   );
    STORE_PD(uY   ,_UY   );
    STORE_PD(uZ   ,_UZ   );
    STORE_PD(theta,_THETA);
    STORE_PD(Pxx  ,_PXX  );
    STORE_PD(Pyy  ,_PYY  );
    STORE_PD(Pzz  ,_PZZ  );
    STORE_PD(Pxy  ,_PXY  );
    STORE_PD(Pyz  ,_PYZ  );
    STORE_PD(Pzx  ,_PZX  );
}


template<typename dataType1>
void getEnstrophyMoments(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx, dataType1 *Sxx, dataType1 *Syy, dataType1 *Szz, dataType1 *Sxy, dataType1 *Syz, dataType1 *Szx,dataType1 dt, dataType1 beta)
{
    getGrad11Moment(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);//Sxx,Syy,Szz,Sxy,Syz,Szx
    dataType1 factor[VECT_LENGTH];
    for(int i=0;i<VECT_LENGTH;i++)
    {
        factor[i] = beta/(dt*rho[i]*theta[i]);
        Sxx[i] = (rho[i]*uX[i]*uX[i] + rho[i]*theta[i] - Pxx[i]) * factor[i];
        Syy[i] = (rho[i]*uY[i]*uY[i] + rho[i]*theta[i] - Pyy[i]) * factor[i];
        Szz[i] = (rho[i]*uZ[i]*uZ[i] + rho[i]*theta[i] - Pzz[i]) * factor[i];
        Sxy[i] = (rho[i]*uX[i]*uY[i] - Pxy[i]) * factor[i];
        Syz[i] = (rho[i]*uY[i]*uZ[i] - Pyz[i]) * factor[i];
        Szx[i] = (rho[i]*uZ[i]*uX[i] - Pzx[i]) * factor[i];
    }
}


template<typename dataType1>
void getEnstrophyMomentsSIMD(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH, dataType1 *rho, dataType1 *uX, dataType1 *uY, dataType1 *uZ, dataType1 *theta, dataType1 *Pxx, dataType1 *Pyy, dataType1 *Pzz, dataType1 *Pxy, dataType1 *Pyz, dataType1 *Pzx, dataType1 *Sxx, dataType1 *Syy, dataType1 *Szz, dataType1 *Sxy, dataType1 *Syz, dataType1 *Szx,dataType1 dt, dataType1 beta)
{
    getGrad11MomentSIMD(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx);
    dataType1 factor[VECT_LENGTH];
    SIMD_REG _RHO, _UX, _UY, _UZ;              //4
    SIMD_REG _PXX,_PYY,_PZZ,_PXY;              //4
    SIMD_REG _PYZ,_PZX,_THETA,_twoReg;         //4
    SIMD_REG _fourReg,_temp1,_temp2,_temp4;    //4
    dataType1 *pointer1;

    pointer1   = &rho[0];
    _RHO       = LOAD_PD(pointer1);
    pointer1   = &uX[0];
    _UX        = LOAD_PD(pointer1);
    pointer1   = &uY[0];
    _UY        = LOAD_PD(pointer1);
    pointer1   = &uZ[0];
    _UZ        = LOAD_PD(pointer1);
    pointer1   = &theta[0];
    _THETA     = LOAD_PD(pointer1);
    pointer1   = &Pxx[0];
    _PXX       = LOAD_PD(pointer1);
    pointer1   = &Pyy[0];
    _PYY       = LOAD_PD(pointer1);
    pointer1   = &Pzz[0];
    _PZZ       = LOAD_PD(pointer1);
    pointer1   = &Pxy[0];
    _PXY       = LOAD_PD(pointer1);
    pointer1   = &Pyz[0];
    _PYZ       = LOAD_PD(pointer1);
    pointer1   = &Pzx[0];
    _PZX       = LOAD_PD(pointer1);

    _temp1 = SET1_PD(beta);
    _temp2 = SET1_PD(dt);
    _temp1 = DIV_PD(_temp1,_temp2);   // beta/dt
    _temp2 = MUL_PD(_RHO,_THETA);     // rho*theta
    _temp1 = DIV_PD(_temp1,_temp2);   // beta/(dt*rho*theta)

    _PXX   = SUB_PD(_temp2,_PXX);     // rho*theta - PXX
    _PYY   = SUB_PD(_temp2,_PYY);
    _PZZ   = SUB_PD(_temp2,_PZZ);

    _temp2 = MUL_PD(_UX,_UX);
    _temp2 = MUL_PD(_temp2,_RHO);    // rho*Ux^2
    _temp2 = ADD_PD(_temp2,_PXX);    // rho*Ux^2 + rho*theta - PXX
    _PXX   = MUL_PD(_temp1,_temp2);  // SXX = (rho*Ux^2 + rho*theta - PXX)*factor

    _temp2 = MUL_PD(_UY,_UY);
    _temp2 = MUL_PD(_temp2,_RHO);    // rho*Uy^2
    _temp2 = ADD_PD(_temp2,_PYY);    // rho*Uy^2 + rho*theta - PYY
    _PYY   = MUL_PD(_temp1,_temp2);

    _temp2 = MUL_PD(_UZ,_UZ);
    _temp2 = MUL_PD(_temp2,_RHO);    // rho*Uy^2
    _temp2 = ADD_PD(_temp2,_PZZ);    // rho*Uy^2 + rho*theta - PYY
    _PZZ   = MUL_PD(_temp1,_temp2);

    _temp2 = MUL_PD(_UX,_UY);
    _temp2 = MUL_PD(_temp2,_RHO);    // rho*Ux*Uy
    _temp2 = SUB_PD(_temp2,_PXY);    // rho*Ux*Uy
    _PXY   = MUL_PD(_temp1,_temp2);

    _temp2 = MUL_PD(_UY,_UZ);
    _temp2 = MUL_PD(_temp2,_RHO);
    _temp2 = SUB_PD(_temp2,_PYZ);
    _PYZ   = MUL_PD(_temp1,_temp2);

    _temp2 = MUL_PD(_UZ,_UX);
    _temp2 = MUL_PD(_temp2,_RHO);
    _temp2 = SUB_PD(_temp2,_PZX);
    _PZX   = MUL_PD(_temp1,_temp2);


    STORE_PD(Sxx  ,_PXX  );
    STORE_PD(Syy  ,_PYY  );
    STORE_PD(Szz  ,_PZZ  );
    STORE_PD(Sxy  ,_PXY  );
    STORE_PD(Syz  ,_PYZ  );
    STORE_PD(Szx  ,_PZX  );
}



template <int N,int numblock, typename dataType1>
void printEnstrophy(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 u_inlet,dataType1 dt, dataType1 beta,int step,dataType1 convectionTime,int size,int myRank, std::ofstream &file)
{
    dataType1 rho[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 uX[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 uY[VECT_LENGTH] __attribute__ ((aligned(32))) ;
    dataType1 uZ[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 theta[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 Pxx[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 Pyy[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 Pzz[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 Pxy[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 Pyz[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 Pzx[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 Sxx[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 Syy[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 Szz[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 Sxy[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 Syz[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 Szx[VECT_LENGTH]  __attribute__ ((aligned(32)));
    dataType1 enstrophy[VECT_LENGTH]  __attribute__ ((aligned(32)));

    dataType1 rhoSum(0.0),thetasum(0.0),internalEnergy(0.0),kineticEnergy(0.0),enstrophySum(0.0),enstrophyMax(0.0),global_EnstrophyMax(0.0);
    for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1=i1+4)
            {
                copyFromNode(lbModel,myGrid,VECT_LENGTH,i1,i2,i3);
                getEnstrophyMomentsSIMD(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx,Sxx,Syy,Szz,Sxy,Syz,Szx,dt,beta);
                for(int i=0;i<VECT_LENGTH;i++)
                {
                    rhoSum           += rho[i];
                    thetasum         += theta[i];
                    internalEnergy   += 1.5*rho[i]*theta[i];
                    kineticEnergy    += 0.5*rho[i]*(uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i]);
                    enstrophy[i]      = Sxx[i]*Sxx[i] + Syy[i]*Syy[i] + Szz[i]*Szz[i] + 2.0*(Sxy[i]*Sxy[i] + Syz[i]*Syz[i] + Szx[i]*Szx[i]);
                    if(enstrophy[i] > enstrophyMax)
                        enstrophyMax = enstrophy[i];
                    enstrophySum     += enstrophy[i];
                }
            }


     for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
         for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
             for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1=i1+4)
             {
                 copyFromCell(lbModel,myGrid,VECT_LENGTH,i1,i2,i3);
                 getEnstrophyMomentsSIMD(lbModel,VECT_LENGTH,rho,uX,uY,uZ,theta,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx,Sxx,Syy,Szz,Sxy,Syz,Szx,dt,beta);
                 for(int i=0;i<VECT_LENGTH;i++)
                 {
                     rhoSum           += rho[i];
                     thetasum         += theta[i];
                     internalEnergy   += 1.5*rho[i]*theta[i];
                     kineticEnergy    += 0.5*rho[i]*(uX[i]*uX[i] + uY[i]*uY[i] + uZ[i]*uZ[i]);
                     enstrophy[i]      = Sxx[i]*Sxx[i] + Syy[i]*Syy[i] + Szz[i]*Szz[i] + 2.0*(Sxy[i]*Sxy[i] + Syz[i]*Syz[i] + Szx[i]*Szx[i]);
                     if(enstrophy[i] > enstrophyMax)
                         enstrophyMax = enstrophy[i];
                     enstrophySum     += enstrophy[i];
                 }
             }


                    dataType1  global_rho,global_theta,global_KE,global_IE,global_TE,global_Enstrophy;

                    dataType1 oneByDen      = 1.0/(2.0*myGrid.actualGridSize);
                    dataType1 oneByU_inlet2 = 1.0/(u_inlet*u_inlet);
                    kineticEnergy           = kineticEnergy  * oneByDen * oneByU_inlet2;
                    internalEnergy          = internalEnergy * oneByDen * oneByU_inlet2;
                    dataType1 totalEnergy   = internalEnergy + kineticEnergy;
                    rhoSum                  = rhoSum*oneByDen;
                    thetasum                = thetasum*oneByDen;
                    enstrophySum            = enstrophySum*oneByDen*oneByU_inlet2;

                    MPI_Reduce(&rhoSum         ,&global_rho         , 1, MPI_DOUBLE, MPI_SUM   ,0,MPI_COMM_WORLD);
                    MPI_Reduce(&thetasum       ,&global_theta       , 1, MPI_DOUBLE, MPI_SUM   ,0,MPI_COMM_WORLD);
                    MPI_Reduce(&kineticEnergy  ,&global_KE          , 1, MPI_DOUBLE, MPI_SUM   ,0,MPI_COMM_WORLD);
                    MPI_Reduce(&internalEnergy ,&global_IE          , 1, MPI_DOUBLE, MPI_SUM   ,0,MPI_COMM_WORLD);
                    MPI_Reduce(&totalEnergy    ,&global_TE          , 1, MPI_DOUBLE, MPI_SUM   ,0,MPI_COMM_WORLD);
                    MPI_Reduce(&enstrophySum   ,&global_Enstrophy   , 1, MPI_DOUBLE, MPI_SUM   ,0,MPI_COMM_WORLD);
                    MPI_Reduce(&enstrophyMax   ,&global_EnstrophyMax, 1, MPI_DOUBLE, MPI_MAX   ,0,MPI_COMM_WORLD);
                    int space(16);
                    if(myRank==0)
                        file<<static_cast<dataType1>(step/convectionTime) << std::setw(space) << global_rho/size <<std::setw(space)<< global_theta/size<<std::setw(space)<< global_KE/size<<std::setw(space)<< global_IE/size<<std::setw(space)<< global_TE/size<< std::setw(space)<< global_Enstrophy/size<< /*std::setw(space)<< global_EnstrophyMax <<*/std::endl;
}


template <int N,int numblock, typename dataType1>
void printTransverseKineticEnergy(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,int VECT_LENGTH,dataType1 u_inlet,dataType1 dt, dataType1 beta,int step,dataType1 convectionTime,int size,int myRank, std::ofstream &file,int yLocation, int *coord,dataType1 F1, dataType1 F2, dataType1 F3, int coresYtot)
{
    dataType1 rho,uX,uY,uZ,theta;
    int sizeGrid = 1;

    dataType1 rhoSum(0.0),thetasum(0.0),transverseKineticEnergy(0.0);
    for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
            {
             int jGlobal = i2 + myGrid.m2*coord[1];

             if(jGlobal == yLocation)
             {
              copyFromNodeSinglePoint(lbModel,myGrid,i1,i2,i3,0);
              getHydroMomentSinglePointWithForce(lbModel,rho,uX,uY,uZ,theta,F1,F2,F3,dt);

              rhoSum                     += rho;
              thetasum                   += theta;
              transverseKineticEnergy    += 0.5*rho*(uY*uY + uZ*uZ);

	      //sizeGrid++;
             }
            }

            dataType1  global_rho,global_theta,global_TKE;

            dataType1 oneByDen      = 1.0/(myGrid.m1*myGrid.m3);
            dataType1 oneByU_inlet2 = 1.0/(u_inlet*u_inlet);
            transverseKineticEnergy = transverseKineticEnergy  * oneByDen * oneByU_inlet2;
            rhoSum                  = rhoSum   * oneByDen;
            thetasum                = thetasum * oneByDen;

            MPI_Reduce(&rhoSum                   ,&global_rho         , 1, MPI_DOUBLE, MPI_SUM   ,0,MPI_COMM_WORLD);
            MPI_Reduce(&thetasum                 ,&global_theta       , 1, MPI_DOUBLE, MPI_SUM   ,0,MPI_COMM_WORLD);
            MPI_Reduce(&transverseKineticEnergy  ,&global_TKE         , 1, MPI_DOUBLE, MPI_SUM   ,0,MPI_COMM_WORLD);
            int space(16);
            if(myRank==0)
                file<<static_cast<dataType1>(step/convectionTime) << std::setw(space) << global_rho/(sizeGrid*coresYtot) <<std::setw(space)<< global_theta/(sizeGrid*coresYtot)<<std::setw(space)<< global_TKE/(sizeGrid*coresYtot)<<std::endl;
}

template <int N,int numblock, typename dataType1>
void printVelocityAtTheCenterOfDomain(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int *coord, dataType1 F1, dataType1 F2, dataType1 F3,int centerX, int centerY, int centerZ, std::ofstream &file,int step,dataType1 convectionTime,dataType1 dt)
{
    dataType1 rho,uX,uY,uZ,theta;

    for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
            {
             int iGlobal = i1 + myGrid.m1*coord[0];
             int jGlobal = i2 + myGrid.m2*coord[1];
             int kGlobal = i3 + myGrid.m3*coord[2];

             if((iGlobal==centerX) && (jGlobal == centerY) && (kGlobal == centerZ))
             {
              copyFromNodeSinglePoint(lbModel,myGrid,i1,i2,i3,0);
              getHydroMomentSinglePointWithForce(lbModel,rho,uX,uY,uZ,theta,F1,F2,F3,dt);
              file<<dataType1(step)/convectionTime<<"  "<<iGlobal<<"  "<<jGlobal<<"  "<<kGlobal<<"  "<<rho<<"  "<<uX<<"  "<<uY<<"  "<<uZ<<"  "<<theta<<std::endl;
             }
            }
}





template <typename dataType1>
inline void printFvalues(lbmRD3Q41<dataType1> &lbModel,int VECT_LENGTH,int index)
{
    for(int i=0;i<VECT_LENGTH;i++)
        std::cout<<"Group 1 : "<< lbModel.fTemp1[index][i] <<"   ";
    std::cout<<std::endl;

    for(int i=0;i<VECT_LENGTH;i++)
        std::cout<<"Group 2 : "<< lbModel.fTemp2[index][i] <<"   ";
    std::cout<<std::endl;

    for(int i=0;i<VECT_LENGTH;i++)
        std::cout<<"Group 3 : "<< lbModel.fTemp3[index][i] <<"   ";
    std::cout<<std::endl;

    for(int i=0;i<VECT_LENGTH;i++)
        std::cout<<"Group 4 : "<< lbModel.fTemp4[index][i] <<"   ";
    std::cout<<std::endl;

    for(int i=0;i<VECT_LENGTH;i++)
        std::cout<<"Group 5 : "<< lbModel.fTemp5[index][i] <<"   ";
    std::cout<<std::endl;

    for(int i=0;i<VECT_LENGTH;i++)
        std::cout<<"Group 6 : "<< lbModel.fTemp6[index][i] <<"   ";
    std::cout<<std::endl;

    for(int i=0;i<VECT_LENGTH;i++)
        std::cout<<"Group 7 : "<< lbModel.fTemp7[index][i] <<"   ";
    std::cout<<std::endl;

    for(int i=0;i<VECT_LENGTH;i++)
        std::cout<<"Group 8 : "<< lbModel.fTemp8[index][i] <<"   ";
    std::cout<<std::endl;

    for(int i=0;i<VECT_LENGTH;i++)
        std::cout<<"Group 9 : "<< lbModel.fTemp9[index][i] <<"   ";
    std::cout<<std::endl;

    for(int i=0;i<VECT_LENGTH;i++)
        std::cout<<"Group 10 : "<< lbModel.fTemp10[index][i] <<"   ";
    std::cout<<std::endl;

}


template <int N,int numblock, typename dataType1>
inline void getForce(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int SOLID,dataType1 *forceNode,dataType1 *forceCell)
{
    forceNode[0] = forceCell[0] = 0.0;
    forceNode[1] = forceCell[1] = 0.0;
    forceNode[2] = forceCell[2] = 0.0;


    for(int k=myGrid.nB3; k<=myGrid.nE3;k++)
        for(int j=myGrid.nB2; j<=myGrid.nE2;j++)
            for(int i=myGrid.nB1; i<=myGrid.nE1;i++)
            {
                if(marker(i,j,k,0,0,0) == SOLID)
                {

                    forceNode[0]  +=  2.0*(  myGrid(i,j,k,lbModel.G1 ,myGrid.node[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO)
                        - myGrid(i,j,k,lbModel.G2 ,myGrid.node[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO))
                        + myGrid(i,j,k,lbModel.G3 ,myGrid.node[lbModel.G3],lbModel.G3_DV_P1_ZERO_P1)
                        - myGrid(i,j,k,lbModel.G3 ,myGrid.node[lbModel.G3],lbModel.G3_DV_M1_ZERO_P1)
                        - myGrid(i,j,k,lbModel.G4 ,myGrid.node[lbModel.G4],lbModel.G4_DV_M1_ZERO_M1)
                        + myGrid(i,j,k,lbModel.G4 ,myGrid.node[lbModel.G4],lbModel.G4_DV_P1_ZERO_M1)
                        + myGrid(i,j,k,lbModel.G5 ,myGrid.node[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO)
                        - myGrid(i,j,k,lbModel.G5 ,myGrid.node[lbModel.G5],lbModel.G5_DV_M1_P1_ZERO)
                        + myGrid(i,j,k,lbModel.G5 ,myGrid.node[lbModel.G5],lbModel.G5_DV_P1_ZERO_ZERO)
                        - myGrid(i,j,k,lbModel.G6 ,myGrid.node[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO)
                        + myGrid(i,j,k,lbModel.G6 ,myGrid.node[lbModel.G6],lbModel.G6_DV_P1_M1_ZERO)
                        - myGrid(i,j,k,lbModel.G6 ,myGrid.node[lbModel.G6],lbModel.G6_DV_M1_ZERO_ZERO)
                        +0.5*(   myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P)
                        - myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P)
                        - myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P)
                        + myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P)
                        - myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M)
                        + myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M)
                        + myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M)
                        - myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M) )
                        + myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_P1_P1)
                        - myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9],lbModel.G9_DV_M1_P1_P1)
                        - myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9],lbModel.G9_DV_M1_M1_P1)
                        + myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_M1_P1)
                        - myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1)
                        + myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1)
                        + myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1)
                        - myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1) ;


                        forceNode[1]  +=  2.0*(  myGrid(i,j,k,lbModel.G1 ,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_P2_ZERO)
                            - myGrid(i,j,k,lbModel.G2 ,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_M2_ZERO))
                            + myGrid(i,j,k,lbModel.G3 ,myGrid.node[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1)
                            - myGrid(i,j,k,lbModel.G3 ,myGrid.node[lbModel.G3],lbModel.G3_DV_ZERO_M1_P1)
                            - myGrid(i,j,k,lbModel.G4 ,myGrid.node[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1)
                            + myGrid(i,j,k,lbModel.G4 ,myGrid.node[lbModel.G4],lbModel.G4_DV_ZERO_P1_M1)
                            + myGrid(i,j,k,lbModel.G5 ,myGrid.node[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO)
                            + myGrid(i,j,k,lbModel.G5 ,myGrid.node[lbModel.G5],lbModel.G5_DV_M1_P1_ZERO)
                            + myGrid(i,j,k,lbModel.G5 ,myGrid.node[lbModel.G5],lbModel.G5_DV_ZERO_P1_ZERO)
                            - myGrid(i,j,k,lbModel.G6 ,myGrid.node[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO)
                            - myGrid(i,j,k,lbModel.G6 ,myGrid.node[lbModel.G6],lbModel.G6_DV_P1_M1_ZERO)
                            - myGrid(i,j,k,lbModel.G6 ,myGrid.node[lbModel.G6],lbModel.G6_DV_ZERO_M1_ZERO)
                            +0.5*(   myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P)
                            + myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P)
                            - myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P)
                            - myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P)
                            - myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M)
                            - myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M)
                            + myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M)
                            + myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M) )
                            + myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_P1_P1)
                            + myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9],lbModel.G9_DV_M1_P1_P1)
                            - myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9],lbModel.G9_DV_M1_M1_P1)
                            - myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_M1_P1)
                            - myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1)
                            - myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1)
                            + myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1)
                            + myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1) ;


                            forceNode[2]  +=  2.0*(  myGrid(i,j,k,lbModel.G1 ,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P2)
                                - myGrid(i,j,k,lbModel.G2 ,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M2))
                                + myGrid(i,j,k,lbModel.G1 ,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1)
                                - myGrid(i,j,k,lbModel.G2 ,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1)
                                + myGrid(i,j,k,lbModel.G3 ,myGrid.node[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1)
                                + myGrid(i,j,k,lbModel.G3 ,myGrid.node[lbModel.G3],lbModel.G3_DV_ZERO_M1_P1)
                                + myGrid(i,j,k,lbModel.G3 ,myGrid.node[lbModel.G3],lbModel.G3_DV_P1_ZERO_P1)
                                + myGrid(i,j,k,lbModel.G3 ,myGrid.node[lbModel.G3],lbModel.G3_DV_M1_ZERO_P1)
                                - myGrid(i,j,k,lbModel.G4 ,myGrid.node[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1)
                                - myGrid(i,j,k,lbModel.G4 ,myGrid.node[lbModel.G4],lbModel.G4_DV_ZERO_P1_M1)
                                - myGrid(i,j,k,lbModel.G4 ,myGrid.node[lbModel.G4],lbModel.G4_DV_M1_ZERO_M1)
                                - myGrid(i,j,k,lbModel.G4 ,myGrid.node[lbModel.G4],lbModel.G4_DV_P1_ZERO_M1)
                                +0.5*(   myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P)
                                + myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P)
                                + myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P)
                                + myGrid(i,j,k,lbModel.G7 ,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P)
                                - myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M)
                                - myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M)
                                - myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M)
                                - myGrid(i,j,k,lbModel.G8 ,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M) )
                                + myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_P1_P1)
                                + myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9],lbModel.G9_DV_M1_P1_P1)
                                + myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9],lbModel.G9_DV_M1_M1_P1)
                                + myGrid(i,j,k,lbModel.G9 ,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_M1_P1)
                                - myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1)
                                - myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1)
                                - myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1)
                                - myGrid(i,j,k,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1) ;
                }
                if(marker(i,j,k,0,1,0) == SOLID)
                {

                    forceCell[0]  +=  2.0*(  myGrid(i,j,k,lbModel.G1 ,myGrid.cell[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO)
                        - myGrid(i,j,k,lbModel.G2 ,myGrid.cell[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO))
                        + myGrid(i,j,k,lbModel.G3 ,myGrid.cell[lbModel.G3],lbModel.G3_DV_P1_ZERO_P1)
                        - myGrid(i,j,k,lbModel.G3 ,myGrid.cell[lbModel.G3],lbModel.G3_DV_M1_ZERO_P1)
                        - myGrid(i,j,k,lbModel.G4 ,myGrid.cell[lbModel.G4],lbModel.G4_DV_M1_ZERO_M1)
                        + myGrid(i,j,k,lbModel.G4 ,myGrid.cell[lbModel.G4],lbModel.G4_DV_P1_ZERO_M1)
                        + myGrid(i,j,k,lbModel.G5 ,myGrid.cell[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO)
                        - myGrid(i,j,k,lbModel.G5 ,myGrid.cell[lbModel.G5],lbModel.G5_DV_M1_P1_ZERO)
                        + myGrid(i,j,k,lbModel.G5 ,myGrid.cell[lbModel.G5],lbModel.G5_DV_P1_ZERO_ZERO)
                        - myGrid(i,j,k,lbModel.G6 ,myGrid.cell[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO)
                        + myGrid(i,j,k,lbModel.G6 ,myGrid.cell[lbModel.G6],lbModel.G6_DV_P1_M1_ZERO)
                        - myGrid(i,j,k,lbModel.G6 ,myGrid.cell[lbModel.G6],lbModel.G6_DV_M1_ZERO_ZERO)
                        +0.5*(   myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P)
                        - myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P)
                        - myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P)
                        + myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P)
                        - myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M)
                        + myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M)
                        + myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M)
                        - myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M) )
                        + myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_P1_P1)
                        - myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9],lbModel.G9_DV_M1_P1_P1)
                        - myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9],lbModel.G9_DV_M1_M1_P1)
                        + myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_M1_P1)
                        - myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1)
                        + myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1)
                        + myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1)
                        - myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1) ;


                        forceCell[1]  +=  2.0*(  myGrid(i,j,k,lbModel.G1 ,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_P2_ZERO)
                            - myGrid(i,j,k,lbModel.G2 ,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_M2_ZERO))
                            + myGrid(i,j,k,lbModel.G3 ,myGrid.cell[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1)
                            - myGrid(i,j,k,lbModel.G3 ,myGrid.cell[lbModel.G3],lbModel.G3_DV_ZERO_M1_P1)
                            - myGrid(i,j,k,lbModel.G4 ,myGrid.cell[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1)
                            + myGrid(i,j,k,lbModel.G4 ,myGrid.cell[lbModel.G4],lbModel.G4_DV_ZERO_P1_M1)
                            + myGrid(i,j,k,lbModel.G5 ,myGrid.cell[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO)
                            + myGrid(i,j,k,lbModel.G5 ,myGrid.cell[lbModel.G5],lbModel.G5_DV_M1_P1_ZERO)
                            + myGrid(i,j,k,lbModel.G5 ,myGrid.cell[lbModel.G5],lbModel.G5_DV_ZERO_P1_ZERO)
                            - myGrid(i,j,k,lbModel.G6 ,myGrid.cell[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO)
                            - myGrid(i,j,k,lbModel.G6 ,myGrid.cell[lbModel.G6],lbModel.G6_DV_P1_M1_ZERO)
                            - myGrid(i,j,k,lbModel.G6 ,myGrid.cell[lbModel.G6],lbModel.G6_DV_ZERO_M1_ZERO)
                            +0.5*(   myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P)
                            + myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P)
                            - myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P)
                            - myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P)
                            - myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M)
                            - myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M)
                            + myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M)
                            + myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M) )
                            + myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_P1_P1)
                            + myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9],lbModel.G9_DV_M1_P1_P1)
                            - myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9],lbModel.G9_DV_M1_M1_P1)
                            - myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_M1_P1)
                            - myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1)
                            - myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1)
                            + myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1)
                            + myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1) ;


                            forceCell[2]  +=  2.0*(  myGrid(i,j,k,lbModel.G1 ,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P2)
                                - myGrid(i,j,k,lbModel.G2 ,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M2))
                                + myGrid(i,j,k,lbModel.G1 ,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1)
                                - myGrid(i,j,k,lbModel.G2 ,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1)
                                + myGrid(i,j,k,lbModel.G3 ,myGrid.cell[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1)
                                + myGrid(i,j,k,lbModel.G3 ,myGrid.cell[lbModel.G3],lbModel.G3_DV_ZERO_M1_P1)
                                + myGrid(i,j,k,lbModel.G3 ,myGrid.cell[lbModel.G3],lbModel.G3_DV_P1_ZERO_P1)
                                + myGrid(i,j,k,lbModel.G3 ,myGrid.cell[lbModel.G3],lbModel.G3_DV_M1_ZERO_P1)
                                - myGrid(i,j,k,lbModel.G4 ,myGrid.cell[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1)
                                - myGrid(i,j,k,lbModel.G4 ,myGrid.cell[lbModel.G4],lbModel.G4_DV_ZERO_P1_M1)
                                - myGrid(i,j,k,lbModel.G4 ,myGrid.cell[lbModel.G4],lbModel.G4_DV_M1_ZERO_M1)
                                - myGrid(i,j,k,lbModel.G4 ,myGrid.cell[lbModel.G4],lbModel.G4_DV_P1_ZERO_M1)
                                +0.5*(   myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P)
                                + myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P)
                                + myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P)
                                + myGrid(i,j,k,lbModel.G7 ,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P)
                                - myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M)
                                - myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M)
                                - myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M)
                                - myGrid(i,j,k,lbModel.G8 ,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M) )
                                + myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_P1_P1)
                                + myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9],lbModel.G9_DV_M1_P1_P1)
                                + myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9],lbModel.G9_DV_M1_M1_P1)
                                + myGrid(i,j,k,lbModel.G9 ,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_M1_P1)
                                - myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1)
                                - myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1)
                                - myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1)
                                - myGrid(i,j,k,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1) ;
                }
            }

}


template <int N,int numblock, typename dataType1>
inline void globalMassTotalFluid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step,int size,int myRank)
{
    dataType1 nodeMass(0.0),cellMass(0.0);
    int nodeCount(0),cellCount(0);

    int nodeCountGlobal(0.0);
    int cellCountGlobal(0.0);
    dataType1 nodeMassGlobal(0.0);
    dataType1 cellMassGlobal(0.0);

    for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
            {
                if(marker(i1,i2,i3,0,0,0) == FLUID)
                {
                    nodeMass += (  myGrid(i1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO)
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

                    // std::cout<<"Grid Node Section: "<<myRank<<std::endl<<std::endl;

                    // if(myRank == 0)
                    // {
                    //   std::cout<<"Grid Values Node: "<<myGrid(i1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P2)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_ZERO_P2_ZERO)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G1,myGrid.node[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M2)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_ZERO_M2_ZERO)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G2,myGrid.node[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],lbModel.G3_DV_ZERO_M1_P1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],lbModel.G3_DV_P1_ZERO_P1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G3,myGrid.node[lbModel.G3],lbModel.G3_DV_M1_ZERO_P1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],lbModel.G4_DV_M1_ZERO_M1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],lbModel.G4_DV_ZERO_P1_M1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G4,myGrid.node[lbModel.G4],lbModel.G4_DV_P1_ZERO_M1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO  )
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],lbModel.G5_DV_M1_P1_ZERO  )
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],lbModel.G5_DV_ZERO_P1_ZERO)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G5,myGrid.node[lbModel.G5],lbModel.G5_DV_P1_ZERO_ZERO)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO  )
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],lbModel.G6_DV_P1_M1_ZERO  )
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],lbModel.G6_DV_ZERO_M1_ZERO)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G6,myGrid.node[lbModel.G6],lbModel.G6_DV_M1_ZERO_ZERO)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_P_P)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_P_P)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_M_M_P)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G7,myGrid.node[lbModel.G7],lbModel.G7_DV_P_M_P)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_M_M)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_M_P_M)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_M_M)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G8,myGrid.node[lbModel.G8],lbModel.G8_DV_P_P_M)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_P1_P1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],lbModel.G9_DV_M1_P1_P1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],lbModel.G9_DV_M1_M1_P1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G9,myGrid.node[lbModel.G9],lbModel.G9_DV_P1_M1_P1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_M1_M1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_M1_P1_M1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_M1_M1)
                    //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G10,myGrid.node[lbModel.G10],lbModel.G10_DV_P1_P1_M1)<<std::endl;
                    // }

                    nodeCount++;
                }
            }
            for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
                for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
                    for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
                    {
                        if(marker(i1,i2,i3,0,1,0) == FLUID)
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
                            + myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1) );

                            // std::cout<<"Grid Cell Section: "<<myRank<<std::endl<<std::endl;

                            // if(myRank == 0)
                            // {
                            //   std::cout<<"Grid Values Cell: "<<myGrid(i1,i2,i3,lbModel.G0,myGrid.cell[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_ZERO_P2)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_ZERO_P2_ZERO)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G1,myGrid.cell[lbModel.G1],lbModel.G1_DV_P2_ZERO_ZERO)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_ZERO_M2)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_ZERO_M2_ZERO)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G2,myGrid.cell[lbModel.G2],lbModel.G2_DV_M2_ZERO_ZERO)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],lbModel.G3_DV_ZERO_P1_P1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],lbModel.G3_DV_ZERO_M1_P1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],lbModel.G3_DV_P1_ZERO_P1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G3,myGrid.cell[lbModel.G3],lbModel.G3_DV_M1_ZERO_P1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],lbModel.G4_DV_M1_ZERO_M1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],lbModel.G4_DV_ZERO_P1_M1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],lbModel.G4_DV_ZERO_M1_M1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G4,myGrid.cell[lbModel.G4],lbModel.G4_DV_P1_ZERO_M1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],lbModel.G5_DV_P1_P1_ZERO  )
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],lbModel.G5_DV_M1_P1_ZERO  )
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],lbModel.G5_DV_ZERO_P1_ZERO)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G5,myGrid.cell[lbModel.G5],lbModel.G5_DV_P1_ZERO_ZERO)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],lbModel.G6_DV_M1_M1_ZERO  )
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],lbModel.G6_DV_P1_M1_ZERO  )
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],lbModel.G6_DV_ZERO_M1_ZERO)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G6,myGrid.cell[lbModel.G6],lbModel.G6_DV_M1_ZERO_ZERO)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_P_P)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_P_P)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_M_M_P)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G7,myGrid.cell[lbModel.G7],lbModel.G7_DV_P_M_P)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_M_M)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_M_P_M)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_M_M)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G8,myGrid.cell[lbModel.G8],lbModel.G8_DV_P_P_M)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_P1_P1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],lbModel.G9_DV_M1_P1_P1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],lbModel.G9_DV_M1_M1_P1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G9,myGrid.cell[lbModel.G9],lbModel.G9_DV_P1_M1_P1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_M1_M1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_M1_P1_M1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_M1_M1)
                            //   <<"\n"<<myGrid(i1,i2,i3,lbModel.G10,myGrid.cell[lbModel.G10],lbModel.G10_DV_P1_P1_M1)<<std::endl;
                            // }

                            cellCount++;

                        }
                    }

                    //             nodeMass = nodeMass;
                    //             cellMass = cellMass;
                    //
                    // MPI_Reduce(&nodeCount,&nodeCountGlobal, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
                    // MPI_Reduce(&cellCount,&cellCountGlobal, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
                    // MPI_Reduce(&nodeMass ,&nodeMassGlobal, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
                    // MPI_Reduce(&cellMass ,&cellMassGlobal, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
                    //             int density = 1;

                    // if(myRank == 0)
                    {
                        std::cout<<"FLUID: "<< std::fixed<<std::setw(14)<<myRank<<std::fixed<<std::setw(14)<<step<<"   "<<std::setw(26)<<(nodeMass/*/(nodeCountGlobal)*/)<<"  "<<std::setw(26)<<(cellMass/*/(cellCountGlobal)*/)<<"    "<<std::setw(16)<<(nodeMass+cellMass)<<"      "<<std::setw(16)<<(nodeMass+cellMass)/((cellCount+nodeCount))<<"      "<<std::setw(16)<<nodeCount<<"        "<<cellCount<<std::endl;//<<density*(cellCount+nodeCount)<<std::endl;
                    }

                    // if((nodeMassGlobal > 1.0e6) || (cellMassGlobal > 1.0e6))
                    // {
                    //     std::cout<<"-nan here: "<<step<<" "<<myRank<<" "<<nodeMass<<" "<<cellMass<<std::endl;
                    // }
                    // if(myRank == 0)
                    // {
                    //     std::cout<<"FLUID: "<< std::fixed<<std::setw(14)<<step<<"   "<<std::setw(26)<<(nodeMassGlobal/*/(nodeCountGlobal)*/)<<"  "<<std::setw(26)<<(cellMassGlobal/*/(cellCountGlobal)*/)<<"    "<<std::setw(16)<<(nodeMassGlobal+cellMassGlobal)<<"      "<<std::setw(16)<<(nodeMassGlobal+cellMassGlobal)/((cellCountGlobal+nodeCountGlobal))<<"      "<<std::setw(16)<<nodeCountGlobal<<"        "<<cellCountGlobal<<std::endl;//<<density*(cellCountGlobal+nodeCountGlobal)<<std::endl;
                    // }
}

template <int N,int numblock, typename dataType1>
inline void globalMassTotalSolid(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step, int size,int myRank)
{
    dataType1 nodeMass(0.0),cellMass(0.0);
    int nodeCount(0),cellCount(0);

    int nodeCountGlobal(0.0);
    int cellCountGlobal(0.0);
    dataType1 nodeMassGlobal(0.0);
    dataType1 cellMassGlobal(0.0);

    for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
            {
                if(marker(i1,i2,i3,0,0,0) == SOLID)
                {
                    nodeMass += (  myGrid(i1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO)
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
            }
            for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
                for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
                    for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
                    {
                        if(marker(i1,i2,i3,0,1,0) == SOLID)
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
                    }

                    //             nodeMass = nodeMass;
                    //             cellMass = cellMass;

                    MPI_Reduce(&nodeCount,&nodeCountGlobal, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
                    MPI_Reduce(&cellCount,&cellCountGlobal, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
                    MPI_Reduce(&nodeMass ,&nodeMassGlobal, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
                    MPI_Reduce(&cellMass ,&cellMassGlobal, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);

                    if(myRank == 0)
                    {
                        std::cout<<"SOLID: "<< std::fixed<<std::setw(14)<<step<<"   "<<std::setw(26)<<(nodeMassGlobal/*/(nodeCountGlobal)*/)<<"  "<<std::setw(26)<<(cellMassGlobal/*/(cellCountGlobal)*/)<<"    "<<std::setw(16)<<(nodeMassGlobal+cellMassGlobal)<<"      "<<std::setw(16)<<(nodeMassGlobal+cellMassGlobal)/((cellCountGlobal+nodeCountGlobal))<<"      "<<std::setw(16)<<nodeCountGlobal<<"        "<<cellCountGlobal<<std::endl;//<<density*(cellCountGlobal+nodeCountGlobal)<<std::endl;
                    }
}

template <int N,int numblock, typename dataType1>
inline void globalMassInternal(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step, std::ofstream &file,int size,int myRank, int *coord)
{
    dataType1 nodeMass(0.0),cellMass(0.0);
    int nodeCount(0),cellCount(0);

    int nodeCountGlobal(0.0);
    int cellCountGlobal(0.0);
    dataType1 nodeMassGlobal(0.0);
    dataType1 cellMassGlobal(0.0);

    int XWALL1 = 30;
    int XWALL2 = 80;
    int YWALL1 = 30;
    int YWALL2 = 80;
    int ZWALL1 = 30;
    int ZWALL2 = 80;



    for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
            {

                int iGlobal = (int) i1 + (int)(marker.m1) * coord[0] ;
                int jGlobal = (int) i2 + (int)(marker.m2) * coord[1] ;
                int kGlobal = (int) i3 + (int)(marker.m3) * coord[2] ;

                if( ( (iGlobal>=XWALL1 && iGlobal<=XWALL2) && (jGlobal>=YWALL1 && jGlobal<=YWALL2) && (kGlobal>=ZWALL1 && kGlobal<=ZWALL2) ) && (marker(i1,i2,i3,0,0,0) == FLUID)  )
                {
                    nodeMass += (  myGrid(i1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO)
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
            }
            for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
                for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
                    for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
                    {


                        int iGlobal = (int) i1 + (int)(marker.m1) * coord[0] ;
                        int jGlobal = (int) i2 + (int)(marker.m2) * coord[1] ;
                        int kGlobal = (int) i3 + (int)(marker.m3) * coord[2] ;

                        if( ( (jGlobal>=YWALL1 && jGlobal<YWALL2) && (kGlobal>=ZWALL1 && kGlobal<ZWALL2) && (iGlobal >= XWALL1 && iGlobal < XWALL2) ) && (marker(i1,i2,i3,0,1,0) == FLUID)  )
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
                    }

                    //             nodeMass = nodeMass/nodeCount;
                    //             cellMass = cellMass/cellCount;

                    MPI_Reduce(&nodeCount,&nodeCountGlobal, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
                    MPI_Reduce(&cellCount,&cellCountGlobal, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
                    MPI_Reduce(&nodeMass ,&nodeMassGlobal, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
                    MPI_Reduce(&cellMass ,&cellMassGlobal, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);

                    if(myRank == 0)
                    {
                        file<< std::fixed<<std::setw(4)<<step<<"   "<<std::setw(16)<<(nodeMassGlobal/(nodeCountGlobal))<<std::setw(16)<<(cellMassGlobal/(cellCountGlobal))<<std::setw(16)<<(nodeMassGlobal+cellMassGlobal)/((cellCountGlobal+nodeCountGlobal))<<std::setw(16)<<nodeCountGlobal<<"        "<<cellCountGlobal<<std::endl;//<<density*(cellCountGlobal+nodeCountGlobal)<<std::endl;
                    }
}


template <int N,int numblock, typename dataType1>
inline void globalMassExternal(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid,gridBCC3D<1, 1, int> &marker,int VECT_LENGTH,int SOLID,int FLUID,int step, std::ofstream &file,int size,int myRank, int *coord)
{
    dataType1 nodeMass(0.0),cellMass(0.0);
    int nodeCount(0),cellCount(0);

    int nodeCountGlobal(0.0);
    int cellCountGlobal(0.0);
    dataType1 nodeMassGlobal(0.0);
    dataType1 cellMassGlobal(0.0);

    int XWALL1 = 30;
    int XWALL2 = 80;
    int YWALL1 = 30;
    int YWALL2 = 80;
    int ZWALL1 = 30;
    int ZWALL2 = 80;



    for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
        for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
            for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
            {

                int iGlobal = (int) i1 + (int)(marker.m1) * coord[0] ;
                int jGlobal = (int) i2 + (int)(marker.m2) * coord[1] ;
                int kGlobal = (int) i3 + (int)(marker.m3) * coord[2] ;

                if( ( (iGlobal<=XWALL1 || iGlobal>=XWALL2) && (jGlobal<=YWALL1 || jGlobal>=YWALL2) && (kGlobal<=ZWALL1 || kGlobal>=ZWALL2) ) && (marker(i1,i2,i3,0,0,0) == FLUID)  )
                {
                    nodeMass += (  myGrid(i1,i2,i3,lbModel.G0,myGrid.node[lbModel.G0],lbModel.CENTER_DV_ZERO_ZERO_ZERO)
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
            }
            for(int i3=myGrid.nB3; i3<=myGrid.nE3;i3++)
                for(int i2=myGrid.nB2; i2<=myGrid.nE2;i2++)
                    for(int i1=myGrid.nB1; i1<=myGrid.nE1;i1++)
                    {

                        int iGlobal = (int) i1 + (int)(marker.m1) * coord[0] ;
                        int jGlobal = (int) i2 + (int)(marker.m2) * coord[1] ;
                        int kGlobal = (int) i3 + (int)(marker.m3) * coord[2] ;

                        if( ( (jGlobal<=YWALL1-1 || jGlobal>=YWALL2) && (kGlobal<=ZWALL1-1 || kGlobal>=ZWALL2) && (iGlobal<=XWALL1-1 || iGlobal>=XWALL2) ) && (marker(i1,i2,i3,0,1,0) == FLUID)  )
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
                    }

                    //             nodeMass = nodeMass/nodeCount;
                    //             cellMass = cellMass/cellCount;

                    MPI_Reduce(&nodeCount,&nodeCountGlobal, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
                    MPI_Reduce(&cellCount,&cellCountGlobal, 1, MPI_INT, MPI_SUM,0,MPI_COMM_WORLD);
                    MPI_Reduce(&nodeMass ,&nodeMassGlobal, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
                    MPI_Reduce(&cellMass ,&cellMassGlobal, 1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);

                    if(myRank == 0)
                    {
                        file<< std::fixed<<std::setw(4)<<step<<"   "<<std::setw(16)<<(nodeMassGlobal/(nodeCountGlobal))<<std::setw(16)<<(cellMassGlobal/(cellCountGlobal))<<std::setw(16)<<(nodeMassGlobal+cellMassGlobal)/((cellCountGlobal+nodeCountGlobal))<<std::setw(16)<<nodeCountGlobal<<"        "<<cellCountGlobal<<std::endl;//<<density*(cellCountGlobal+nodeCountGlobal)<<std::endl;
                    }
}


template <int N,int numblock, typename dataType1>
inline void printGridZplane(lbmRD3Q41<dataType1> &lbModel, gridBCC3D<N, numblock, dataType1> &myGrid, int *coord,int n1,int n2,int n3,int ZPlaneValue)
{
    //     for(int k=myGrid.nB3; k<=myGrid.nE3;k++)
    //     {
    for(int j=myGrid.ndB2 ;j<=myGrid.ndE2;j++)
    {
        for(int i=myGrid.ndB1; i<=myGrid.ndE1;i++)
        {
            //           copyFromNode(lbModel,compGrid,lbModel.fTemp,i,j,k);
            copyFromNodeSinglePoint(lbModel,myGrid,i,j,ZPlaneValue,0);
            std::cout << "(" << i+coord[0]*(myGrid.m1) << "," << j+coord[1]*(myGrid.m2) << "," << ZPlaneValue << ") " <<"("<< lbModel.fTemp1[0][lbModel.G1_DV_P2_ZERO_ZERO] << ", " << lbModel.fTemp2[0][lbModel.G2_DV_M2_ZERO_ZERO] <<")"<< '\t';
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}















// Explicit template declarations
template void getLatticeParameter<double>(lbmRD3Q41<double> &);
template void setModelParameters<double>(lbmRD3Q41<double> &);
template void getHydroMoment<double>(lbmRD3Q41<double> &,int , double *, double *, double *, double *, double *);
template void copyFromNode<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,int ,int , int , int );
template void copyToNode<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,int ,int , int , int );
template void copyFromCell<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,int ,int , int , int );
template void copyToCell<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,int ,int , int , int );
template void copyToNodeSinglePoint<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,int ,int ,int ,int );
template void copyToCellSinglePoint<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,int ,int ,int ,int );
template void copyFromNodeSinglePoint<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,int ,int ,int ,int );
template void copyFromCellSinglePoint<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4,11,double> &,int ,int ,int ,int );

template void getHydroMomentsFromNode<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,int , int , int , double *, double *, double *, double *, double *);
template void getHydroMomentsFromCell<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,int , int , int , double *, double *, double *, double *, double *);
template void getHydroMomentsFromNodeWithForce<4,11,double>(lbmRD3Q41<double> &,gridBCC3D<4,11,double> &,int ,int , int , int , double *, double *, double *, double *, double *,double , double , double, double );
template void getHydroMomentsFromCellWithForce<4,11,double>(lbmRD3Q41<double> &,gridBCC3D<4,11,double> &,int ,int , int , int , double *, double *, double *, double *, double *,double , double , double, double );
template void globalMassTotalFluid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,gridBCC3D<1, 1, int> &,int ,int ,int ,int ,int ,int );
template void globalMassTotalSolid<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,gridBCC3D<1, 1, int> &,int ,int ,int ,int ,int ,int );
template void globalMassInternal<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,gridBCC3D<1, 1, int> &,int ,int ,int ,int , std::ofstream &,int ,int ,int *);
template void globalMassExternal<4,11,double>(lbmRD3Q41<double> &, gridBCC3D<4, 11, double> &,gridBCC3D<1, 1, int> &,int ,int ,int ,int , std::ofstream &,int ,int ,int *);
template void getGrad11Moment<double>(lbmRD3Q41<double> &,int , double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
template void getEnstrophyMoments(lbmRD3Q41<double> &,int , double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,double , double );
template void getEnstrophyMomentsSIMD(lbmRD3Q41<double> &,int , double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,double , double );
template void printEnstrophy<4,11,double>(lbmRD3Q41<double> &lbModel, gridBCC3D<4,11,double> &myGrid,int VECT_LENGTH,double u_inlet,double dt, double beta,int step,double convectionTime,int size,int myRank, std::ofstream &file);
template void getGrad11MomentSIMD<double>(lbmRD3Q41<double> &,int , double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

template void getHydroMomentSinglePoint<double>(lbmRD3Q41<double> &, double &, double &, double &, double &, double &);
template void getHydroMomentSinglePointWithForce<double>(lbmRD3Q41<double> &, double &, double &, double &, double &, double &,double , double , double ,double );

template void printTransverseKineticEnergy<4,11,double>(lbmRD3Q41<double> &lbModel, gridBCC3D<4,11,double> &myGrid,int VECT_LENGTH,double u_inlet,double dt, double beta,int step,double convectionTime,int size,int myRank, std::ofstream &file,int yLocation, int *coord,double F1, double F2, double F3, int coresYtot);



template void printVelocityAtTheCenterOfDomain<4,11,double>(lbmRD3Q41<double> &lbModel, gridBCC3D<4,11,double> &myGrid, int *coord, double F1, double F2, double F3,int centerX, int centerY, int centerZ, std::ofstream &file,int step,double convectionTime,double dt);


template void getGrad11Moment_heatAlpha<double>(lbmRD3Q41<double> &,int , double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

template void getGrad11Moment_heatAlpha_YminusFlux<double>(lbmRD3Q41<double> &,int , double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

template void getGrad11Moment_heatAlpha_YplusFlux<double>(lbmRD3Q41<double> &,int , double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

template void calcLBTopFlux(lbmRD3Q41<double> &,int , gridBCC3D<4, 11, double> &, int , int , int , int , int , int , double* , int);
template void calcLBBotFlux(lbmRD3Q41<double> &,int , gridBCC3D<4, 11, double> &, int , int , int , int , int , int , double* );
