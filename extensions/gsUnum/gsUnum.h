
#pragma once

#ifndef NDEBUG
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#endif

#include <universal/posit/posit>

typedef sw::unum::posit<256,5> posit_256_5;
typedef sw::unum::posit<128,4> posit_128_4;
typedef sw::unum::posit< 64,3> posit_64_3;
typedef sw::unum::posit< 32,2> posit_32_2;
typedef sw::unum::posit< 16,1> posit_16_1;
typedef sw::unum::posit<  8,1> posit_8_1;
typedef sw::unum::posit<  8,0> posit_8_0;
typedef sw::unum::posit<  4,0> posit_4_0;
typedef sw::unum::posit<  3,1> posit_3_1;
typedef sw::unum::posit<  3,0> posit_3_0;
typedef sw::unum::posit<  2,0> posit_2_0;
