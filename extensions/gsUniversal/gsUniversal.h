
#pragma once

#ifndef NDEBUG
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#endif

#include <universal/utility/directives.hpp>
#include <universal/number/posit/posit.hpp>

typedef sw::universal::posit<256,5> posit_256_5;
typedef sw::universal::posit<128,4> posit_128_4;
typedef sw::universal::posit< 64,3> posit_64_3;
typedef sw::universal::posit< 32,2> posit_32_2;
typedef sw::universal::posit< 16,1> posit_16_1;
typedef sw::universal::posit<  8,1> posit_8_1;
typedef sw::universal::posit<  8,0> posit_8_0;
typedef sw::universal::posit<  4,0> posit_4_0;
typedef sw::universal::posit<  3,1> posit_3_1;
typedef sw::universal::posit<  3,0> posit_3_0;
typedef sw::universal::posit<  2,0> posit_2_0;
