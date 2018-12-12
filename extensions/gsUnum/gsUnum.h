
#pragma once

#include <gsCore/gsConfig.h>

#ifndef NDEBUG
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#endif

#include <posit>

typedef sw::unum::posit<32,2> posit_32_2;
