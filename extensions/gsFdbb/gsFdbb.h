
#pragma once

#include <gsCore/gsConfig.h>

#ifndef SUPPORT_EIGEN
#define SUPPORT_EIGEN
#endif

#if !defined(SUPPORT_CXX11) && !defined(SUPPORT_CXX14)
#if __cplusplus == 201103L
#define SUPPORT_CXX11
#elif __cplusplus == 201402L
#define SUPPORT_CXX14
#else
#error "FDBB requires C++11 or C++14 support enabled"
#endif
#endif

#include <fdbb.h>
