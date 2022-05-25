#ifndef _MARKER_H_
#define _MARKER_H_

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>
#include"procInfo3D.h"
#include"gridBCC3DBasic.h"
#include"paniniAliases.h"
#include"D3RQ41.h"
#include"simdHelper.h"

void initializeMarker(gridBCC3D<1, 1, int> &marker, PANINI_INT FLUID);
void markSphereIntParallel(gridBCC3D<1, 1, int> &marker, PANINI_INT SOLID, int XCENTER,int YCENTER, int ZCENTER, int RADIUS, int *coord);
void markTopAndBottomWalls      (gridBCC3D<1, 1, int> &marker, PANINI_INT FLUID, PANINI_INT SOLID, int *coord,int nP1,int nP2,int nP3);
void markOnePercentBumpOnTheWall(gridBCC3D<1, 1, int> &marker, PANINI_INT FLUID, PANINI_INT SOLID, int *coord,int nP1,int nP2,int nP3,int nX,int nY,int nZ);

#endif
