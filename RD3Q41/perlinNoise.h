#ifndef _PERLINNOISE_H_
#define _PERLINNOISE_H_

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<mpi.h>
#include"gridBCC3DBasic.h"
#include"paniniAliases.h"
#include"D3RQ41.h"
#include"procInfo3D.h"
#include <time.h>


double noise3(double x1, double x2, double x3);
static void normalize3(double v[3]);
static void init(void);


#endif
