/** @file gsBinaryFunctions.h

    @brief Various helper classes derived from from std::binary_function for
    stl-style functional programming


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C Hofer -- taken from DUNE
    Created on: 2016-03-23
*/


#pragma once

// #include <functional>
// #include <algorithm>

namespace gismo
{
template<typename Type>
struct Min
    : std::binary_function<Type,Type,Type>
{
    Type operator()(const Type& t1, const Type& t2) const
    {
        return std::min(t1,t2);
    }
};

template<typename Type>
struct Max
    : std::binary_function<Type,Type,Type>
{
    Type operator()(const Type& t1, const Type& t2) const
    {
        return std::max(t1,t2);
    }
};



template<typename Type, typename BinaryFunction>
class Generic_MPI_Op
{

public:
    static MPI_Op get ()
    {
        if (!op) // if op is null then create new op
        {
            op = memory::shared_ptr<MPI_Op>(new MPI_Op);
            MPI_Op_create((void (*)(void*, void*, int*, MPI_Datatype*))&operation,true,op.get());
        }
        return *op;
    }
private:
    static void operation (Type *in, Type *inout, int *len, MPI_Datatype*)
    {
        BinaryFunction func;

        for (int i=0; i< *len; ++i, ++in, ++inout) {
            Type temp;
            temp = func(*in, *inout);
            *inout = temp;
        }
    }
    Generic_MPI_Op () {}
    Generic_MPI_Op (const Generic_MPI_Op& ) {}
    static memory::shared_ptr<MPI_Op> op;     // to check: can do with static MPI_Op ?
};


template<typename Type, typename BinaryFunction>
memory::shared_ptr<MPI_Op> Generic_MPI_Op<Type,BinaryFunction>::op = memory::shared_ptr<MPI_Op>(static_cast<MPI_Op*>(0));

#define ComposeMPIOp(type,func,op)              \
    template<>                                  \
    class Generic_MPI_Op<type, func<type> >{    \
    public:                                     \
    static MPI_Op get(){                        \
        return op;                              \
    }                                           \
    private:                                    \
    Generic_MPI_Op () {}                        \
    Generic_MPI_Op (const Generic_MPI_Op & ) {} \
    }


ComposeMPIOp(char, std::plus, MPI_SUM);
ComposeMPIOp(unsigned char, std::plus, MPI_SUM);
ComposeMPIOp(short, std::plus, MPI_SUM);
ComposeMPIOp(unsigned short, std::plus, MPI_SUM);
ComposeMPIOp(int, std::plus, MPI_SUM);
ComposeMPIOp(unsigned int, std::plus, MPI_SUM);
ComposeMPIOp(long, std::plus, MPI_SUM);
ComposeMPIOp(unsigned long, std::plus, MPI_SUM);
ComposeMPIOp(float, std::plus, MPI_SUM);
ComposeMPIOp(double, std::plus, MPI_SUM);
ComposeMPIOp(long double, std::plus, MPI_SUM);

ComposeMPIOp(char, std::multiplies, MPI_PROD);
ComposeMPIOp(unsigned char, std::multiplies, MPI_PROD);
ComposeMPIOp(short, std::multiplies, MPI_PROD);
ComposeMPIOp(unsigned short, std::multiplies, MPI_PROD);
ComposeMPIOp(int, std::multiplies, MPI_PROD);
ComposeMPIOp(unsigned int, std::multiplies, MPI_PROD);
ComposeMPIOp(long, std::multiplies, MPI_PROD);
ComposeMPIOp(unsigned long, std::multiplies, MPI_PROD);
ComposeMPIOp(float, std::multiplies, MPI_PROD);
ComposeMPIOp(double, std::multiplies, MPI_PROD);
ComposeMPIOp(long double, std::multiplies, MPI_PROD);

ComposeMPIOp(char, Min, MPI_MIN);
ComposeMPIOp(unsigned char, Min, MPI_MIN);
ComposeMPIOp(short, Min, MPI_MIN);
ComposeMPIOp(unsigned short, Min, MPI_MIN);
ComposeMPIOp(int, Min, MPI_MIN);
ComposeMPIOp(unsigned int, Min, MPI_MIN);
ComposeMPIOp(long, Min, MPI_MIN);
ComposeMPIOp(unsigned long, Min, MPI_MIN);
ComposeMPIOp(float, Min, MPI_MIN);
ComposeMPIOp(double, Min, MPI_MIN);
ComposeMPIOp(long double, Min, MPI_MIN);

ComposeMPIOp(char, Max, MPI_MAX);
ComposeMPIOp(unsigned char, Max, MPI_MAX);
ComposeMPIOp(short, Max, MPI_MAX);
ComposeMPIOp(unsigned short, Max, MPI_MAX);
ComposeMPIOp(int, Max, MPI_MAX);
ComposeMPIOp(unsigned int, Max, MPI_MAX);
ComposeMPIOp(long, Max, MPI_MAX);
ComposeMPIOp(unsigned long, Max, MPI_MAX);
ComposeMPIOp(float, Max, MPI_MAX);
ComposeMPIOp(double, Max, MPI_MAX);
ComposeMPIOp(long double, Max, MPI_MAX);

#undef ComposeMPIOp


}
