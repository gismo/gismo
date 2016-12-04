# Precompiled Header creation for C++
#
# Author: Adam Strzelecki <ono@java.pl>
# Copyright (c) 2014-2015 Adam Strzelecki. All rights reserved.
# This code is licensed under the MIT License, see README.md.
#
# Defines compiler command templates based on existing C++ compiler.

include(CMakePCHCompiler)

__define_pch_compiler(CXX)
