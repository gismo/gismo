# Precompiled Header creation for C++
#
# Author: Adam Strzelecki <ono@java.pl>
# Copyright (c) 2014-2015 Adam Strzelecki. All rights reserved.
# This code is licensed under the MIT License, see README.md.
#
# Main entry point for new compiler. Here it just proxies to C++ compiler.

include(CMakePCHCompiler)

__configure_pch_compiler(CXX)
