/*********************************************************************************************
 *   Copyright (c) <2022>, <Santosh Ansumali@JNCASR>                                         *
 *   All rights reserved.                                                                    *
 *   Redistribution and use in source and binary forms, with or without modification, are    *
 *   permitted provided that the following conditions are met:                               *
 *                                                                                           *
 *    1. Redistributions of source code must retain the above copyright notice, this list of *
 *       conditions and the following disclaimer.                                            *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list *
 *       of conditions and the following disclaimer in the documentation and/or other        *
 *       materials provided with the distribution.                                           *
 *    3. Neither the name of the <JNCASR> nor the names of its contributors may be used to   *
 *       endorse or promote products derived from this software without specific prior       *
 *       written permission.                                                                 *
 *                                                                                           *
 *       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND     *
 *       ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED       *
 *       WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  *
 *       IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,    *
 *       INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,      *
 *       BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,       *
 *       DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF     *
 *       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE     *
 *       OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED   *
 *       OF THE POSSIBILITY OF SUCH DAMAGE.                                                  *
 *                                                                                           *
 *       Suggestions:          ansumali@jncasr.ac.in                                         *
 *       Bugs:                 ansumali@jncasr.ac.in                                         *
 *                                                                                           *
 *********************************************************************************************/


/*********************************************************************************************

 *  @Author: Akshay Chandran and Santosh Ansumali                                            *

 *********************************************************************************************/


#ifndef PARTICLE_H_INCLUDED
#define PARTICLE_H_INCLUDED

#include "dof.h"
#include <mm_malloc.h>

template<int DIM,typename real,typename dof>
class particle
{
    public:
        particle()
        {}

        ~particle()
        {
            _mm_free(data);
            _mm_free(cellCoord);
        }

        void initializeParticles(int Xsize, int YSize, int ZSize, int nCx, int nCy, int nCz, int nPCell)
        {
            numParticles = ZSize*(YSize*Xsize/(nCx*nCy))*nPCell*2.;
            numElements  = DIM*numParticles;
            data         = (dof*) _mm_malloc(numElements*sizeof(dof), 4096);
            cellCoord    = (int*) _mm_malloc(numParticles*sizeof(int), 4096);
        }

        dof& getData(int particlenum, int dim) {return data[(int)numParticles*(dim - 1) + particlenum];}

        int& getCellCoord(int particlenum) {return cellCoord[particlenum];}

        int getnumParticles() const {return numParticles;}

        int& setnumParticles() {return numParticles;}

        int getDim() {return DIM;}

        real getnumElements() {return numElements;}

        void swapParticles(int, int);

	      dof& operator()(int particlenum, int dim) {return data[(int)numParticles*(dim - 1) + particlenum];}

    private:
        int numParticles;
        int numElements;
        int cellSizeZ;
        dof* data;
        int* cellCoord;
};

#endif // PARTICLE_H_INCLUDED
