/** @file gsBoundaryConditions_test.cpp

    @brief test gsPde/gsBoundaryConditions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris,  H. Weiner

**/

#include "gismo_unittest.h"
#include <typeinfo>
#include <memory>

std::string getBoxSideExpectedString(std::string string, int index)
{
    std::ostringstream stream;
    std::string result;
    stream << string << " (" << index << ")";
    result = stream.str();
    return result;
}

std::string getBoxSideActualString(gismo::boxSide boxSide)
{
    std::ostringstream stream;
    std::string result;
    stream << boxSide;
    result = stream.str();
    return result;
}

std::string getFunctionExprExpectedString(std::string string)
{
    std::ostringstream stream;
    std::string result;
    stream << "[ " << string << " ]";
    result = stream.str();
    return result;
}

std::string getFunctionExprActualString(gsFunctionExpr<real_t> func)
{
    std::ostringstream stream;
    std::string result;
    func.print(stream);
    result = stream.str();
    return result;
}

gsBoundaryConditions<real_t> gsBoundaryConditions_loadFromFile(std::string path)
{
    gsBoundaryConditions<real_t> result;
    gsReadFile<>(path, result);
    return result;
}

void gsBoundaryConditions_saveToFile(std::string path,
        gsBoundaryConditions<real_t> &sut)
{
    // clean-up before using file
    remove(path.c_str());
    // now save...
    gsFileData<> write;
    write << sut;
    write.dump(path);
}

// create the same gsBoundary condition as in test-data file "bc.xml"
gsBoundaryConditions<real_t> createBcGsBoundaryConditiions()
{
    // set-up
    gsBoundaryConditions<real_t> result;
    return result;
}

gsBoundaryConditions<real_t> createSimpleGsBoundaryConditiions()
{
    gsBoundaryConditions<real_t> result;
    return result;
}

gsFunctionExpr<real_t> getFunctionExpr(boundary_condition<real_t> bc)
{
    CHECK( 0 != dynamic_cast<gsFunctionExpr<real_t>*>(bc.m_function.get()) );
    gsFunction<real_t>::Ptr ptr = bc.m_function;
    gsFunctionExpr<real_t> * ptr2 =
        dynamic_cast<gsFunctionExpr<real_t> *>(ptr.get());
    gsFunctionExpr<real_t> result = *ptr2;
    return result;
}

void checkBoundaryCondition(boundary_condition<real_t> bc, bool parametric,
        std::string label, gismo::condition_type::type conditionType, index_t patch,
        index_t index, int unknown, int unkcomp, short_t domainDim,
        std::string funcName)
{
    // check boundary condition itself
    CHECK_EQUAL(parametric, bc.m_parametric);
    CHECK_EQUAL(label, bc.m_label);
    CHECK_EQUAL(conditionType, bc.m_type);
    CHECK_EQUAL(patch, bc.ps.patch);
    CHECK_EQUAL(index, bc.ps.m_index);
    CHECK_EQUAL(unknown, bc.m_unknown);
    CHECK_EQUAL(unkcomp, bc.m_unkcomp);
    // check functionExpr
    gismo::gsFunctionExpr<real_t> func = getFunctionExpr(bc);
    CHECK_EQUAL(domainDim, func.domainDim());
    std::string expectedName = getFunctionExprExpectedString(funcName);
    std::string actualName = getFunctionExprActualString(func);
    CHECK_EQUAL(expectedName, actualName);
}

void checkGsBoundaryCondition(const gsBoundaryConditions<real_t> & sut)
{
    // check corner value(s)
    gismo::gsBoundaryConditions<real_t>::cornerContainer c1 =
            sut.cornerValues();
    unsigned elems1 = 1;
    CHECK_EQUAL(elems1, c1.size());
    gismo::corner_value<real_t> cv1 = c1[0];
    CHECK_EQUAL(0, cv1.corner);
    CHECK_EQUAL(0, cv1.corner.m_index);
    CHECK_EQUAL(0, cv1.patch);
    CHECK_EQUAL(0.0, cv1.value);
    CHECK_EQUAL(0, cv1.unknown);

    // check dirichlet conditions
    std::string label1 = "Dirichlet";
    gismo::gsBoundaryConditions<real_t>::bcContainer bcc0 = sut.container(
            label1);
    gismo::gsBoundaryConditions<real_t>::bcContainer bcc1 =
            sut.dirichletSides();
    unsigned elems2 = 6;
    CHECK_EQUAL(elems2, bcc0.size());
    CHECK_EQUAL(elems2, bcc1.size());
    // check boundary conditions themselves
    std::string funcName0 = "0";
    int patches[4] =
    { 0, 0, 1, 1 };
    int indizes[4] =
    { 1, 3, 1, 2 };
    int unknown = 0;
    int unkcomp = -1; //option not present in the XML yields default value -1
    short_t domainDim = 2;
    for (int i = 0; i < 4; i++)
    {
        checkBoundaryCondition(bcc0[i], false, label1,
                gismo::condition_type::dirichlet, patches[i], indizes[i],
                unknown, unkcomp, domainDim, funcName0);
    }

    std::string funcName1 = "sin(pi*x)*sin(pi*y)";
    int patch4 = 0;
    int index4 = 2;
    int patch5 = 0;
    int index5 = 4;
    checkBoundaryCondition(bcc0[4], false, label1,
            gismo::condition_type::dirichlet, patch4, index4, unknown, unkcomp,
            domainDim, funcName1);
    checkBoundaryCondition(bcc0[5], false, label1,
            gismo::condition_type::dirichlet, patch5, index5, unknown, unkcomp,
            domainDim, funcName1);

    // check neumann conditions
    std::string label2 = "Neumann";
    gismo::gsBoundaryConditions<real_t>::bcContainer bcc2 = sut.container(
            label2);
    gismo::gsBoundaryConditions<real_t>::bcContainer bcc3 = sut.neumannSides();
    unsigned elems3 = 2;
    CHECK_EQUAL(elems3, bcc2.size());
    CHECK_EQUAL(elems3, bcc3.size());
    int patch6 = 1;
    int index6 = 3;
    int patch7 = 1;
    int index7 = 4;
    checkBoundaryCondition(bcc2[0], false, label2,
            gismo::condition_type::neumann, patch6, index6, unknown, unkcomp,
            domainDim, funcName0);
    checkBoundaryCondition(bcc3[0], false, label2,
            gismo::condition_type::neumann, patch6, index6, unknown, unkcomp,
            domainDim, funcName0);
    checkBoundaryCondition(bcc2[1], false, label2,
            gismo::condition_type::neumann, patch7, index7, unknown, unkcomp,
            domainDim, funcName0);
    checkBoundaryCondition(bcc3[1], false, label2,
            gismo::condition_type::neumann, patch7, index7, unknown, unkcomp,
            domainDim, funcName0);

    // check robin conditions
    std::string label3 = "Robin";
    gismo::gsBoundaryConditions<real_t>::bcContainer bcc4 = sut.container(
            label3);
    gismo::gsBoundaryConditions<real_t>::bcContainer bcc5 = sut.robinSides();
    unsigned elems4 = 0;
    CHECK_EQUAL(elems4, bcc4.size());
    CHECK_EQUAL(elems4, bcc5.size());
    // check bctype_iterator
    typedef gismo::gsBoundaryConditions<real_t>::const_bciterator bctype_it;
    int c = 0;
    for (bctype_it it = sut.beginAll(); it != sut.endAll(); ++it)
    {
        c++;
    }
    CHECK_EQUAL(3, c);
}

SUITE(gsBoundaryConditions_test)
{
    TEST(test_condition_type)
    {
        gismo::condition_type::type myType = gismo::condition_type::neumann;
        CHECK_EQUAL(gismo::condition_type::neumann, myType);
        CHECK_EQUAL(1, myType);
    }

    TEST(test_function_expr)
    {
        int dim1 = 1;
        std::string funcName1 = "tan(x)";
        gsFunctionExpr<real_t> func1 = gsFunctionExpr<real_t>(funcName1, dim1);
        short_t domainDim1 = func1.domainDim();
        CHECK_EQUAL(dim1, domainDim1);
        std::string expectedName1 = getFunctionExprExpectedString(funcName1);
        std::string actualName1 = getFunctionExprActualString(func1);
        CHECK_EQUAL(expectedName1, actualName1);
        int dim2 = 2;
        std::string funcName2 = "sin(x)*sin(y)";
        gsFunctionExpr<real_t> func2 = gsFunctionExpr<real_t>(funcName2, dim2);
        short_t domainDim2 = func2.domainDim();
        CHECK_EQUAL(dim2, domainDim2);
        std::string expectedName2 = getFunctionExprExpectedString(funcName2);
        std::string actualName2 = getFunctionExprActualString(func2);
        CHECK_EQUAL(expectedName2, actualName2);
    }

    TEST(test_box_side)
    {
        int index1 = 1;
        gismo::boxSide boxSide1 = gismo::boxSide(index1);
        CHECK_EQUAL(index1, boxSide1.m_index);
        std::string expect1 = getBoxSideExpectedString("west", index1);
        std::string actual1 = getBoxSideActualString(boxSide1);
        CHECK_EQUAL(expect1, actual1);
        int index2 = 7;
        gismo::boxSide boxSide2 = gismo::boxSide(index2);
        CHECK_EQUAL(index2, boxSide2.m_index);
        std::string expect2 = getBoxSideExpectedString("side", index2);
        std::string actual2 = getBoxSideActualString(boxSide2);
        CHECK_EQUAL(expect2, actual2);
    }

    TEST(test_boundary_condition)
    {
        int dim1 = 1;
        std::string funcName1 = "tan(x)";
        index_t index1 = 1; // Eigen Index
        index_t index2 = 2; // index_t
        gismo::boxSide boxSide1 = gismo::boxSide(index1);
        gismo::boundary_condition<real_t>::function_ptr funcPtr1 =
                gismo::memory::make_shared(
                        new gsFunctionExpr<real_t>(funcName1, dim1));
        std::string label1 = "Dirichlet";
        bool parametric1 = true;
        int unknown1 = 1;
        int unkcomp1 = 2;
        gismo::boundary_condition<real_t> bound = gismo::boundary_condition<real_t>(
                index2, boxSide1, funcPtr1, label1, unknown1, unkcomp1,
                parametric1);
        CHECK_EQUAL(parametric1, bound.m_parametric);
        CHECK_EQUAL(label1, bound.m_label);
        CHECK_EQUAL(gismo::condition_type::dirichlet, bound.m_type);
        CHECK_EQUAL(index2, bound.ps.patch);
        CHECK_EQUAL(index2, bound.patch());
        CHECK_EQUAL(index1, bound.ps.m_index);
        CHECK_EQUAL(funcPtr1, bound.m_function);
        CHECK_EQUAL(unknown1, bound.m_unknown);
        CHECK_EQUAL(unkcomp1, bound.m_unkcomp);
    }

    TEST(test_box_corner)
    {
        int index1 = 3;
        gismo::boxCorner c1 = gismo::boxCorner(index1);
        CHECK_EQUAL(index1, c1.m_index);
    }

    TEST(test_corner_value)
    {
        int index1 = 3;
        gismo::boxCorner c1 = gismo::boxCorner(index1);
        index_t p1 = 2;
        real_t v1 = 3.0;
        int u1 = 4;
        gismo::corner_value<real_t> cornerVal1 = gismo::corner_value<real_t>(p1, c1,
                v1, u1);
        CHECK_EQUAL(c1, cornerVal1.corner);
        CHECK_EQUAL(index1, cornerVal1.corner.m_index);
        CHECK_EQUAL(p1, cornerVal1.patch);
        CHECK_EQUAL(v1, cornerVal1.value);
        CHECK_EQUAL(u1, cornerVal1.unknown);
    }

    /***
     * test loading from bc.xml
     */
    TEST(load_from_bc_xml)
    {
        std::string path = gsFileManager::findInDataDir( "gsBoundaryConditions/bc.xml" );
        gsBoundaryConditions<real_t> sut = gsBoundaryConditions_loadFromFile(path);
        checkGsBoundaryCondition(sut);
    }

    /***
     * test loading from bc.xml and
     * saving to bc2.xml
     * and ensure that bc.xml and bc2.xml have
     * the same content
     */
    TEST(save_load_bc_xml)
    {
        std::string path1 = gsFileManager::findInDataDir( "gsBoundaryConditions/bc.xml" );
        std::string path2 = gsFileManager::getTempPath() + "/bc2.xml";

        gsBoundaryConditions<real_t> sut = gsBoundaryConditions_loadFromFile(path1);
        gsBoundaryConditions_saveToFile(path2, sut);
        gsBoundaryConditions<real_t> sut2 = gsBoundaryConditions_loadFromFile(path2);
        checkGsBoundaryCondition(sut2);
    }

}
