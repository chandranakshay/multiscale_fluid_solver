#include<iostream>
#include<math.h>

void projectiveEuler(double* momentsDSMCZero, double* momentsDSMCEnd, double* momentsDSMCProj, int DSMC_Count, int DSMC_Count_Actual, int numMomentsAllSpaceCells) {

  // initTimeInv ----- Inverse of the difference between time at moment calculation at momentsDSMCEnd and momentsDSMCZero

  double deltaTProj  = (double)DSMC_Count_Actual - (double)(0.25*DSMC_Count);
  double initTimeInv = 1./((double)(0.5*DSMC_Count));

  for(int currMoment = 0; currMoment < numMomentsAllSpaceCells; currMoment++) {

    momentsDSMCProj[currMoment] = momentsDSMCZero[currMoment] + deltaTProj*(momentsDSMCEnd[currMoment] - momentsDSMCZero[currMoment])*initTimeInv;
    // std::cout<<"momentsProjected: "<<currMoment<<"\t"<<momentsDSMCZero[currMoment]<<"\t"<<momentsDSMCEnd[currMoment]<<"\t"<<momentsDSMCProj[currMoment]<<std::endl;
  }
}
