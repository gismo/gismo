/** @file gsSortedVector_test.cpp

    @brief test gsUtils/gsSortedVector

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris,  H. Weiner
 **/

#include "gismo_unittest.h"

// general declarations/definitions
void push_unique(std::vector<std::string>& vector, std::string& string)
{
    bool contains = std::find(vector.begin(), vector.end(), string) != vector.end();
    if (!contains)
    {
        vector.push_back(string);
    }
}

int getIndex(std::vector<std::string> & vector, std::string& string)
{
    int index = std::find(vector.begin(), vector.end(), string) - vector.begin();
    return index;
}

// declarations/definitions for testing custom sorting of gsSortedVector
class my_clazz
{
public:
    int first;
    int second;

    // We need the operator== to compare two my_clazz objects.
    bool operator==(const my_clazz &other) const
    {
        return (first == other.first) && (second == other.second);
    }
    // We need the operator< to comparare two my_clazz objects.
    /*
    bool operator<(const my_clazz &other) const
    {
        return first < other.first;
    }
    */
};

int getIndex(gismo::gsSortedVector<my_clazz> & vector,
        my_clazz ptr)
{
    int index = std::find(vector.begin(), vector.end(), ptr) - vector.begin();
    return index;
}

/*
typedef struct my_clazz_comparable_first : gismo::gsSortedVectorComparable<my_clazz>
{
    bool operator() (my_clazz a, my_clazz b)
    {
            return (a.first < b.first);
    }
} my_comparable_first;

typedef struct my_clazz_comparable_second : gismo::gsSortedVectorComparable<my_clazz>
{
    bool operator() (my_clazz a, my_clazz b)
    {
            return (a.second < b.second);
    }
} my_comparable_second;
*/

// declarations/definitions for testing gsBoundaryConditions
/*
#define FUNC_PTR gismo::boundary_condition<real_t>::function_ptr

void push_unique(std::vector<gismo::boundary_condition<real_t>::function_ptr> & vector,
        gismo::boundary_condition<real_t>::function_ptr & ptr)
{
    bool contains = std::find(vector.begin(), vector.end(), ptr) != vector.end();
    if (!contains)
    {
        vector.push_back(ptr);
    }
}

int getIndex(std::vector<gismo::boundary_condition<real_t>::function_ptr> & vector,
        gismo::boundary_condition<real_t>::function_ptr & ptr)
{
    int index = std::find(vector.begin(), vector.end(), ptr) - vector.begin();
    return index;
}

int getIndex(std::vector<gismo::boundary_condition<real_t>::function_ptr> & vector,
        gismo::boundary_condition<real_t>::function_ptr & ptr)
{
    int index = std::find(vector.begin(), vector.end(), ptr) - vector.begin();
    return index;
}

FUNC_PTR getFunctionPtrFor(gsFunctionExpr<real_t> & func)
{
    FUNC_PTR result = memory::make_shared(new gsFunctionExpr<real_t>);
    *result = func;
    return result;
}
*/

SUITE(gsSortedVector_test)
{

    TEST(gsSortedVector_with_string)
    {
        std::string gismo = "Hello G+smo";
        std::string world = "Hello World!";
        std::string string1a = gismo;
        std::string string1b = gismo;
        std::string string2a = world;
        std::string string2b = world;
        gsSortedVector<std::string> vector;
        vector.push_sorted_unique(gismo);
        vector.push_sorted_unique(world);
        vector.push_sorted_unique(string1a);
        vector.push_sorted_unique(string1b);
        vector.push_sorted_unique(string2a);
        vector.push_sorted_unique(string2b);
        CHECK_EQUAL(2, (int )vector.size());
        CHECK_EQUAL(0, (int )vector.getIndex(gismo));
        CHECK_EQUAL(1, (int )vector.getIndex(world));
        CHECK_EQUAL(0, (int )vector.getIndex(string1a));
        CHECK_EQUAL(0, (int )vector.getIndex(string1b));
        CHECK_EQUAL(1, (int )vector.getIndex(string2a));
        CHECK_EQUAL(1, (int )vector.getIndex(string2b));
    }

    TEST(gsSortedVector_with_string_2)
    {
        std::string world = "Hello World!";
        std::string gismo = "Hello G+smo";
        std::string worldA = world;
        std::string worldB = world;
        std::string gismoA = gismo;
        std::string gismoB = gismo;
        gsSortedVector<std::string> vector;
        vector.push_sorted_unique(world);
        vector.push_sorted_unique(gismo);
        vector.push_sorted_unique(worldA);
        vector.push_sorted_unique(worldB);
        vector.push_sorted_unique(gismoA);
        vector.push_sorted_unique(gismoB);
        CHECK_EQUAL(2, (int )vector.size());
        CHECK_EQUAL(0, (int )vector.getIndex(gismo));
        CHECK_EQUAL(1, (int )vector.getIndex(world));
        CHECK_EQUAL(0, (int )vector.getIndex(gismoA));
        CHECK_EQUAL(0, (int )vector.getIndex(gismoB));
        CHECK_EQUAL(1, (int )vector.getIndex(worldA));
        CHECK_EQUAL(1, (int )vector.getIndex(worldB));
    }

    TEST(gsSortedVector_normal_vector_with_string)
    {
        std::string gismo = "Hello G+smo!";
        std::string world = "Hello World";
        std::string string1a = gismo;
        std::string string1b = gismo;
        std::string string2a = world;
        std::string string2b = world;
        std::vector<std::string> vector;
        push_unique(vector, gismo);
        push_unique(vector, world);
        push_unique(vector, string1a);
        push_unique(vector, string1b);
        push_unique(vector, string2a);
        push_unique(vector, string2b);
        CHECK_EQUAL(2, (int )vector.size());
        CHECK_EQUAL(0, getIndex(vector, gismo));
        CHECK_EQUAL(1, getIndex(vector, world));
        CHECK_EQUAL(0, getIndex(vector, string1a));
        CHECK_EQUAL(0, getIndex(vector, string1b));
        CHECK_EQUAL(1, getIndex(vector, string2a));
        CHECK_EQUAL(1, getIndex(vector, string2b));
    }

    /*TEST(gsSortedVector_with_my_clazz_first)
    {
        my_comparable_first* comparable_first = new my_comparable_first();
        my_clazz clazzA;
        clazzA.first = 1;
        clazzA.second = 100;
        my_clazz clazzB;
        clazzB.first = 4000;
        clazzB.second = 3;
        my_clazz clazzC;
        clazzC.first = 2;
        clazzC.second = 5;
        gismo::gsSortedVector<my_clazz> vector;
        vector.setComparable(comparable_first);
        vector.push_sorted(clazzA);
        vector.push_sorted(clazzB);
        vector.push_sorted(clazzC);
        CHECK_EQUAL(3, (int )vector.size());
        CHECK_EQUAL(0, getIndex(vector, clazzA));
        CHECK_EQUAL(2, getIndex(vector, clazzB));
        CHECK_EQUAL(1, getIndex(vector, clazzC));
        delete(comparable_first);
    }

    TEST(gsSortedVector_with_my_clazz_second)
    {
        my_comparable_second* comparable_second = new my_comparable_second();

        my_clazz clazzA;
        clazzA.first = 1;
        clazzA.second = 100;
        my_clazz clazzB;
        clazzB.first = 4000;
        clazzB.second = 3;
        my_clazz clazzC;
        clazzC.first = 2;
        clazzC.second = 5;
        gismo::gsSortedVector<my_clazz> vector;
        vector.setComparable(comparable_second);
        vector.push_sorted(clazzA);
        vector.push_sorted(clazzB);
        vector.push_sorted(clazzC);
        CHECK_EQUAL(3, (int )vector.size());
        CHECK_EQUAL(2, getIndex(vector, clazzA));
        CHECK_EQUAL(0, getIndex(vector, clazzB));
        CHECK_EQUAL(1, getIndex(vector, clazzC));
        delete(comparable_second);
    }*/

    /*
    TEST(gsSortedVector_gsFunctionPtr_push_unique_and_get)
    {
        int dim0 = 1;
        int dim1 = 2;
        std::string funcName0 = "sin(x)*sin(y)";
        std::string funcName1 = "tan(x)";
        gsFunctionExpr<real_t> func0 = gsFunctionExpr<real_t>(funcName0, dim0);
        gsFunctionExpr<real_t> func1 = gsFunctionExpr<real_t>(funcName1, dim1);
        FUNC_PTR funcPtr0 = getFunctionPtrFor(func0);
        FUNC_PTR funcPtr1 = getFunctionPtrFor(func1);
        FUNC_PTR funcPtr0a = funcPtr0;
        FUNC_PTR funcPtr0b = funcPtr0;
        FUNC_PTR funcPtr1a = funcPtr1;
        FUNC_PTR funcPtr1b = funcPtr1;
        gsSortedVector<gismo::boundary_condition<real_t>::function_ptr> vector;
        vector.push_sorted_unique(funcPtr0);
        vector.push_sorted_unique(funcPtr1);
        vector.push_sorted_unique(funcPtr0a);
        vector.push_sorted_unique(funcPtr0b);
        vector.push_sorted_unique(funcPtr1a);
        vector.push_sorted_unique(funcPtr1b);

        // type-casting is safe - prevent compiler warnings
        CHECK_EQUAL(2, (int )vector.size());
        CHECK_EQUAL(0, (int )vector.getIndex(funcPtr0));
        CHECK_EQUAL(1, (int )vector.getIndex(funcPtr1));
        CHECK_EQUAL(0, (int )vector.getIndex(funcPtr0a));
        CHECK_EQUAL(0, (int )vector.getIndex(funcPtr0b));
        CHECK_EQUAL(1, (int )vector.getIndex(funcPtr1a));
        CHECK_EQUAL(1, (int )vector.getIndex(funcPtr1b));
    }

    TEST(gsSortedVector_gsFunctionPtr_with_normal_vector)
    {
        int dim0 = 1;
        int dim1 = 2;
        std::string funcName0 = "sin(x)*sin(y)";
        std::string funcName1 = "tan(x)";
        gsFunctionExpr<real_t> func0 = gsFunctionExpr<real_t>(funcName0, dim0);
        gsFunctionExpr<real_t> func1 = gsFunctionExpr<real_t>(funcName1, dim1);
        FUNC_PTR funcPtr0 = getFunctionPtrFor(func0);
        FUNC_PTR funcPtr1 = getFunctionPtrFor(func1);
        FUNC_PTR funcPtr0a = funcPtr0;
        FUNC_PTR funcPtr0b = funcPtr0;
        FUNC_PTR funcPtr1a = funcPtr1;
        FUNC_PTR funcPtr1b = funcPtr1;

        std::vector<gismo::boundary_condition<real_t>::function_ptr> vector;
        push_unique(vector, funcPtr0);
        push_unique(vector, funcPtr1);
        push_unique(vector, funcPtr0a);
        push_unique(vector, funcPtr0b);
        push_unique(vector, funcPtr1a);
        push_unique(vector, funcPtr1b);

        CHECK_EQUAL(2, (int )vector.size());
        CHECK_EQUAL(0, getIndex(vector, funcPtr0));
        CHECK_EQUAL(1, getIndex(vector, funcPtr1));
        CHECK_EQUAL(0, getIndex(vector, funcPtr0a));
        CHECK_EQUAL(0, getIndex(vector, funcPtr0b));
        CHECK_EQUAL(1, getIndex(vector, funcPtr1a));
        CHECK_EQUAL(1, getIndex(vector, funcPtr1b));
    }

    TEST(gsSortedVector_test_sorting_of_gsFunctionExprs)
    {
        int dim = 2;
        std::string tanName = "tan(x)*sin(y)";
        std::string sinName = "sin(x)*sin(y)";
        gsFunctionExpr<real_t> tanFunc = gsFunctionExpr<real_t>(tanName, dim);
        gsFunctionExpr<real_t> sinFunc = gsFunctionExpr<real_t>(sinName, dim);
        gsFunctionExpr<real_t> tanFuncA = tanFunc;
        gsFunctionExpr<real_t> tanFuncB = tanFunc;
        gsFunctionExpr<real_t> sinFuncA = sinFunc;
        gsFunctionExpr<real_t> sinFuncB = sinFunc;
        gsSortedVector<gsFunctionExpr<real_t> > vector;
        vector.push_sorted_unique(tanFunc);
        vector.push_sorted_unique(sinFunc);
        vector.push_sorted_unique(tanFuncA);
        vector.push_sorted_unique(tanFuncB);
        vector.push_sorted_unique(sinFuncA);
        vector.push_sorted_unique(sinFuncB);

        // type-casting is safe - prevent compiler warnings
        CHECK_EQUAL(2, (int )vector.size());
        CHECK_EQUAL(0, (int )vector.getIndex(sinFunc));
        CHECK_EQUAL(1, (int )vector.getIndex(tanFunc));
        CHECK_EQUAL(0, (int )vector.getIndex(sinFuncA));
        CHECK_EQUAL(0, (int )vector.getIndex(sinFuncB));
        CHECK_EQUAL(1, (int )vector.getIndex(tanFuncA));
        CHECK_EQUAL(1, (int )vector.getIndex(tanFuncB));
    }

    TEST(gsSortedVector_test_sorting_of_gsFunctionPtrs)
    {
        int dim = 2;
        std::string tanName = "tan(x)*sin(y)";
        std::string sinName = "sin(x)*sin(y)";
        gsFunctionExpr<real_t> tanFunc = gsFunctionExpr<real_t>(tanName, dim);
        gsFunctionExpr<real_t> sinFunc = gsFunctionExpr<real_t>(sinName, dim);
        FUNC_PTR funcPtrTan = getFunctionPtrFor(tanFunc);
        FUNC_PTR funcPtrSin = getFunctionPtrFor(sinFunc);
        FUNC_PTR funcPtrTanA = funcPtrTan;
        FUNC_PTR funcPtrTanB = funcPtrTan;
        FUNC_PTR funcPtrSinA = funcPtrSin;
        FUNC_PTR funcPtrSinB = funcPtrSin;
        gsSortedVector<gismo::boundary_condition<real_t>::function_ptr> vector;
        vector.push_sorted_unique(funcPtrTan);
        vector.push_sorted_unique(funcPtrSin);
        vector.push_sorted_unique(funcPtrTanA);
        vector.push_sorted_unique(funcPtrTanB);
        vector.push_sorted_unique(funcPtrSinA);
        vector.push_sorted_unique(funcPtrSinB);

        // type-casting is safe - prevent compiler warnings
        CHECK_EQUAL(2, (int )vector.size());
        CHECK_EQUAL(0, (int )vector.getIndex(funcPtrSin));
        CHECK_EQUAL(1, (int )vector.getIndex(funcPtrTan));
        CHECK_EQUAL(0, (int )vector.getIndex(funcPtrSinA));
        CHECK_EQUAL(0, (int )vector.getIndex(funcPtrSinB));
        CHECK_EQUAL(1, (int )vector.getIndex(funcPtrTanA));
        CHECK_EQUAL(1, (int )vector.getIndex(funcPtrTanB));
    }
    */
}
