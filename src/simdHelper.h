/*********************************************************************************************
 * Copyright (c) <2010>, <Santosh Ansumali@JNCASR>                                           *
 *  All rights reserved.                                                                     *
 *   Redistribution and use in source and binary forms, with or without modification, are    *
 *   permitted provided that the following conditions are met:                               *
 *                                                                                           *
 *    1. Redistributions of source code must retain the above copyright notice, this list of *
 *       conditions and the following disclaimer.                                            *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list *
 *       of conditions and the following disclaimer in the documentation and/or other        *
 *       materials provided with the distribution.                                           *
 *    3. Neither the name of the <JNCASR> nor the names of its contributors may be used to   *
 *       endorse or promote products derived from this softwTSCre without specific prior       *
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
// namespace panini{
#ifndef _SIMD_HELPER_PANINI_
#define  _SIMD_HELPER_PANINI_

#ifdef SSE
#include <xmmintrin.h>
#define SIMD_REG __m128d
#define SIMD_WIDTH 2
#define ALIGN_LEN  16
#define SET_ZERO() _mm_setzero_pd()
#define SET1_PD(d) _mm_set1_pd(d)
#define LOAD_PD(p) _mm_load_pd(p)
#define ADD_PD(r0,r1) _mm_add_pd(r0,r1)
#define MUL_PD(r0,r1) _mm_mul_pd(r0,r1)
#define FMADD_PD(r0,r1,r2) r1 = _mm_mul_pd(r0,r1); r2 = _mm_add_pd(r1,r2)
#define SUB_PD(r0,r1) _mm_sub_pd(r0,r1)
#define DIV_PD(r0,r1) _mm_div_pd(r0,r1)
#define STORE_PD(p,r)  _mm_store_pd(p,r)
#elif AVX
#include "immintrin.h"
#define SIMD_REG __m256d
#define SIMD_WIDTH 4
#define ALIGN_LEN  32
#define SET_ZERO() _mm256_setzero_pd()
#define SET1_PD(d) _mm256_set1_pd(d)
#define LOAD_PD(p) _mm256_load_pd(p)
#define ADD_PD(r0,r1) _mm256_add_pd(r0,r1)
#define MUL_PD(r0,r1) _mm256_mul_pd(r0,r1)
#define SUB_PD(r0,r1) _mm256_sub_pd(r0,r1)
#define DIV_PD(r0,r1) _mm256_div_pd(r0,r1)
#define MAX_PD(r0,r1) _mm256_max_pd(r0,r1)
#define MIN_PD(r0,r1) _mm256_min_pd(r0,r1)
#define STORE_PD(p,r) _mm256_store_pd(p,r)
#define AND_PD(r0,r1) _mm256_and_pd(r0,r1)
#define ANDNOT_PD(r0,r1) _mm256_andnot_pd(r0,r1)
#define OR_PD(r0,r1) _mm256_or_pd(r0,r1)
#define SQRT_PD(r0) _mm256_sqrt_pd(r0)
#define AVX_TRANSPOSE4_PD(row0, row1, row2, row3) {              \
    __m256d tmp0, tmp1, tmp2, tmp3;                                 \
                                                                    \
    tmp0   = _mm256_unpacklo_pd((row0), (row1));                    \
    tmp1   = _mm256_unpackhi_pd((row0), (row1));                    \
    tmp2   = _mm256_unpacklo_pd((row2), (row3));                    \
    tmp3   = _mm256_unpackhi_pd((row2), (row3));                    \
                                                                    \
    (row0) = _mm256_permute2f128_pd(tmp0, tmp2, 32);                \
    (row1) = _mm256_permute2f128_pd(tmp1, tmp3, 32);                \
    (row2) = _mm256_permute2f128_pd(tmp0, tmp2, 49);                \
    (row3) = _mm256_permute2f128_pd(tmp1, tmp3, 49);                \
}
#endif

#endif
// }
