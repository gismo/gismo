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
#include <fstream>
#include <sstream>
#include <gsAssembler/gsDofMapperCreator.h>

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
    gsFunctionSet<real_t>::Ptr ptr = bc.m_function;
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

namespace
{

// Two unit-square patches with a degree-2 basis, 2 knots inserted per span
// (3x3 elements, 5x5 functions per patch): small, with boundary dofs on every side.
gsMultiBasis<real_t> traceTestPatches(gsMultiPatch<real_t> & mp)
{
    mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
    gsMultiBasis<real_t> mb(mp);
    mb.degreeElevate(1);
    mb.uniformRefine(2);
    return mb;
}

// Condition set with Dirichlet data g on (patch 0, west) and (patch 1, east),
// both with the trace tr
gsBoundaryConditions<real_t> traceTestConditions(gsFunctionExpr<real_t> & g,
                                                 const gsMultiPatch<real_t> & mp,
                                                 trace_type::type tr)
{
    gsBoundaryConditions<real_t> bcs;
    bcs.addCondition(0, boundary::west, condition_type::dirichlet, tr, &g);
    bcs.addCondition(1, boundary::east, condition_type::dirichlet, tr, &g);
    bcs.setGeoMap(mp);
    return bcs;
}

// Number of conditions in the container on (patch, side) with the trace tr
int countTrace(const gsBoundaryConditions<real_t>::bcContainer & list,
               index_t patch, boxSide side, trace_type::type tr)
{
    int n = 0;
    for (auto it = list.begin(); it != list.end(); ++it)
        if (it->patch() == patch && it->side() == side && it->trace() == tr)
            ++n;
    return n;
}

// The two mappers have the same sizes and the same flat-local to global numbering
void checkSameMapper(const gsDofMapper & a, const gsDofMapper & b)
{
    CHECK_EQUAL(a.size(), b.size());
    CHECK_EQUAL(a.freeSize(), b.freeSize());
    CHECK_EQUAL(a.boundarySize(), b.boundarySize());
    const gsVector<index_t> va = a.asVector(0), vb = b.asVector(0);
    CHECK_EQUAL(va.size(), vb.size());
    if (va.size() == vb.size())
        CHECK(va == vb);
}

} // namespace

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

    // Property: a condition built without a trace selector selects the
    // essential trace. Expected values come from the declaration of
    // trace_type (essential = 0), not from running the code.
    TEST(trace_default_is_essential)
    {
        CHECK_EQUAL(0, trace_type::essential);

        // Direct construction with the (patch, side, function, type, unknown, parametric) constructor
        gsFunctionExpr<real_t> g("0", 2);
        boundary_condition<real_t>::function_ptr fp =
            memory::make_shared_not_owned(&g);
        boundary_condition<real_t> bc(0, boundary::west, fp,
                                      condition_type::dirichlet, 0, false);
        CHECK_EQUAL(trace_type::essential, bc.trace());

        // Conditions added through the overload that has no trace argument
        gsBoundaryConditions<real_t> bcs;
        bcs.addCondition(0, boundary::west, condition_type::dirichlet, &g);
        bcs.addCondition(0, boundary::east, condition_type::neumann, &g);
        CHECK_EQUAL(1u, bcs.dirichletSides().size());
        CHECK_EQUAL(1u, bcs.neumannSides().size());
        CHECK_EQUAL(trace_type::essential, bcs.dirichletSides()[0].trace());
        CHECK_EQUAL(trace_type::essential, bcs.neumannSides()[0].trace());
    }

    // Property: a data file that predates the trace selector (no "trace"
    // attribute on any <bc>) reads as essential everywhere. bc.xml holds 6
    // Dirichlet and 2 Neumann conditions (see checkGsBoundaryCondition); the
    // count check makes sure the loops below are not vacuous.
    TEST(trace_xml_missing_attribute_reads_essential)
    {
        std::string path = gsFileManager::findInDataDir("gsBoundaryConditions/bc.xml");
        gsBoundaryConditions<real_t> bcs = gsBoundaryConditions_loadFromFile(path);
        CHECK_EQUAL(6u, bcs.dirichletSides().size());
        CHECK_EQUAL(2u, bcs.neumannSides().size());

        // Every condition read carries the default trace
        for (size_t i = 0; i < bcs.dirichletSides().size(); ++i)
            CHECK_EQUAL(trace_type::essential, bcs.dirichletSides()[i].trace());
        for (size_t i = 0; i < bcs.neumannSides().size(); ++i)
            CHECK_EQUAL(trace_type::essential, bcs.neumannSides()[i].trace());
    }

    // Property: traces survive a write/read cycle and the "trace" attribute
    // appears in the file. The traces written are the ones assigned below, so
    // the reference is the input itself.
    TEST(trace_xml_round_trip)
    {
        // The normal and the tangential condition on (0, west) share one
        // function, so the writer must keep them apart by their trace; the
        // other two have functions of their own.
        gsFunctionExpr<real_t> g0("0", 2), g1("1", 2), g3("3", 2);
        gsBoundaryConditions<real_t> bcs;
        bcs.addCondition(0, boundary::west, condition_type::dirichlet,
                         trace_type::normal, &g1);
        bcs.addCondition(0, boundary::west, condition_type::dirichlet,
                         trace_type::tangential, &g1);
        bcs.addCondition(0, boundary::south, condition_type::dirichlet, &g0);
        bcs.addCondition(0, boundary::north, condition_type::neumann,
                         trace_type::all, &g3);

        const std::string path = gsFileManager::getTempPath() + "/bc_trace.xml";
        gsBoundaryConditions_saveToFile(path, bcs);

        // Text of the dumped file: a non-essential trace must be written out
        std::ifstream in(path.c_str());
        std::stringstream ss;
        ss << in.rdbuf();
        CHECK(std::string::npos != ss.str().find("trace="));

        // Each condition comes back once, with its trace
        gsBoundaryConditions<real_t> back = gsBoundaryConditions_loadFromFile(path);
        CHECK_EQUAL(3u, back.dirichletSides().size());
        CHECK_EQUAL(1u, back.neumannSides().size());
        CHECK_EQUAL(1, countTrace(back.dirichletSides(), 0, boundary::west, trace_type::normal));
        CHECK_EQUAL(1, countTrace(back.dirichletSides(), 0, boundary::west, trace_type::tangential));
        CHECK_EQUAL(1, countTrace(back.dirichletSides(), 0, boundary::south, trace_type::essential));
        CHECK_EQUAL(1, countTrace(back.neumannSides(), 0, boundary::north, trace_type::all));
    }

    // Property: a set with only default traces writes no "trace" attribute
    // (files stay readable by older versions) and reads back as essential.
    // The set is built both through the old overload and through the new one
    // with trace_type::essential, which must be indistinguishable.
    TEST(trace_xml_default_not_written)
    {
        gsFunctionExpr<real_t> g0("0", 2), g1("1", 2);
        gsBoundaryConditions<real_t> bcs;
        bcs.addCondition(0, boundary::west, condition_type::dirichlet, &g0);
        bcs.addCondition(0, boundary::east, condition_type::dirichlet,
                         trace_type::essential, &g1);
        CHECK_EQUAL(2u, bcs.dirichletSides().size());

        const std::string path = gsFileManager::getTempPath() + "/bc_trace_default.xml";
        gsBoundaryConditions_saveToFile(path, bcs);
        std::ifstream in(path.c_str());
        std::stringstream ss;
        ss << in.rdbuf();
        CHECK(!ss.str().empty());
        CHECK(std::string::npos == ss.str().find("trace="));

        gsBoundaryConditions<real_t> back = gsBoundaryConditions_loadFromFile(path);
        CHECK_EQUAL(2u, back.dirichletSides().size());
        for (size_t i = 0; i < back.dirichletSides().size(); ++i)
            CHECK_EQUAL(trace_type::essential, back.dirichletSides()[i].trace());
    }

    // Property: isSame distinguishes conditions that differ only in the
    // trace, and still identifies conditions with equal traces. Expected
    // truth values follow from the definition (same patch side, unknown,
    // component and trace).
    TEST(trace_isSame)
    {
        gsFunctionExpr<real_t> g("0", 2), h("1", 2);
        gsBoundaryConditions<real_t> d, n1c, n2c;
        d  .addCondition(0, boundary::west, condition_type::dirichlet, &g);
        n1c.addCondition(0, boundary::west, condition_type::dirichlet,
                         trace_type::normal, &g);
        n2c.addCondition(0, boundary::west, condition_type::dirichlet,
                         trace_type::normal, &h);

        // isSame is non-const, so work on copies; n1 and n2 differ in the
        // function only, a and n1 in the trace only
        boundary_condition<real_t> a  = d.dirichletSides()[0];
        boundary_condition<real_t> n1 = n1c.dirichletSides()[0];
        boundary_condition<real_t> n2 = n2c.dirichletSides()[0];
        CHECK_EQUAL(trace_type::normal, n1.trace());
        CHECK(a.isSame(a));
        CHECK(!a.isSame(n1));
        CHECK(!n1.isSame(a));
        CHECK(n1.isSame(n2));
    }

    // Property: conditions of one type on one patch side for one unknown are
    // kept only if they constrain disjoint parts: different components, or
    // the normal and the tangential part. Every other combination, an exact
    // duplicate included, throws and leaves the set as it was. The expected
    // outcomes are the rule itself.
    TEST(overlapping_conditions)
    {
        gsFunctionExpr<real_t> g("0", 2);

        // Accepted: different components; normal with tangential; the same
        // side for another unknown or another condition type
        gsBoundaryConditions<real_t> ok;
        ok.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, false, 0);
        ok.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, false, 1);
        ok.addCondition(0, boundary::east, condition_type::dirichlet, trace_type::normal, &g);
        ok.addCondition(0, boundary::east, condition_type::dirichlet, trace_type::tangential, &g);
        ok.addCondition(0, boundary::east, condition_type::dirichlet, &g, 1);
        ok.addCondition(0, boundary::east, condition_type::neumann, trace_type::tangential, &g);
        CHECK_EQUAL(5u, ok.dirichletSides().size());
        CHECK_EQUAL(1u, ok.neumannSides().size());

        // Refused: each addCondition below overlaps a condition of ok
        CHECK_THROW(ok.addCondition(0, boundary::west, condition_type::dirichlet, &g, 0, false, 1),
                    std::runtime_error); // exact duplicate
        CHECK_THROW(ok.addCondition(0, boundary::west, condition_type::dirichlet, &g),
                    std::runtime_error); // all components over component 0
        CHECK_THROW(ok.addCondition(0, boundary::east, condition_type::dirichlet, trace_type::all, &g),
                    std::runtime_error); // all over normal
        CHECK_THROW(ok.addCondition(0, boundary::east, condition_type::dirichlet, trace_type::normal, &g),
                    std::runtime_error); // normal twice
        CHECK_THROW(ok.addCondition(0, boundary::east, condition_type::dirichlet, &g, 0, false, 0),
                    std::runtime_error); // component 0 over normal
        CHECK_EQUAL(5u, ok.dirichletSides().size());

        // Refused: essential over all, which constrain the same components
        gsBoundaryConditions<real_t> full;
        full.addCondition(0, boundary::west, condition_type::dirichlet, trace_type::all, &g);
        CHECK_THROW(full.addCondition(0, boundary::west, condition_type::dirichlet, &g),
                    std::runtime_error);

        // Refused: a component together with a trace
        gsBoundaryConditions<real_t> comp;
        CHECK_THROW(comp.addCondition(0, boundary::west, condition_type::dirichlet,
                                      trace_type::normal, &g, 0, false, 1),
                    std::runtime_error);
        CHECK_EQUAL(0u, comp.dirichletSides().size());
    }

    // Property: getConditionsForPatch keeps the component, the trace and the
    // parametric flag of each condition, so two conditions on different
    // components of one side are copied as two disjoint conditions.
    TEST(conditions_for_patch_keep_component)
    {
        gsFunctionExpr<real_t> g("0", 2);
        gsBoundaryConditions<real_t> bcs, onPatch;
        bcs.addCondition(1, boundary::west, condition_type::dirichlet, &g, 0, true, 0);
        bcs.addCondition(1, boundary::west, condition_type::dirichlet, &g, 0, false, 1);
        bcs.addCondition(1, boundary::east, condition_type::dirichlet, trace_type::normal, &g);
        bcs.getConditionsForPatch(1, onPatch);

        // Patch 1 becomes patch 0 of the copy
        CHECK_EQUAL(3u, onPatch.dirichletSides().size());
        CHECK_EQUAL(0, onPatch.dirichletSides()[0].patch());
        CHECK_EQUAL(0, onPatch.dirichletSides()[0].unkComponent());
        CHECK(onPatch.dirichletSides()[0].parametric());
        CHECK_EQUAL(1, onPatch.dirichletSides()[1].unkComponent());
        CHECK(!onPatch.dirichletSides()[1].parametric());
        CHECK_EQUAL(trace_type::normal, onPatch.dirichletSides()[2].trace());
    }

    // Property: createMapper accepts the traces essential and all, which
    // both constrain every component, and refuses normal and tangential,
    // which it cannot impose. The reference for the accepted traces is the
    // mapper of the set with the default trace.
    TEST(trace_createMapper)
    {
        gsMultiPatch<real_t> mp;
        gsMultiBasis<real_t> mb = traceTestPatches(mp);
        gsFunctionExpr<real_t> g("x*y", 2);
        gsBoundaryConditions<real_t> E = traceTestConditions(g, mp, trace_type::essential);
        gsBoundaryConditions<real_t> A = traceTestConditions(g, mp, trace_type::all);
        gsBoundaryConditions<real_t> N = traceTestConditions(g, mp, trace_type::normal);
        gsBoundaryConditions<real_t> Tn = traceTestConditions(g, mp, trace_type::tangential);

        gsDofMapper mE = createMapper(mb, E, 1, 0, true, true);
        gsDofMapper mA = createMapper(mb, A, 1, 0, true, true);
        CHECK(mE.boundarySize() > 0);
        checkSameMapper(mE, mA);

        CHECK_THROW(createMapper(mb, N, 1, 0, true, true), std::runtime_error);
        CHECK_THROW(createMapper(mb, Tn, 1, 0, true, true), std::runtime_error);
    }

    // Property: an expression-assembler space gives the same mapper and the
    // same interpolated Dirichlet values for the traces essential and all,
    // and refuses normal in its setup and in gsDirichletValues.
    TEST(trace_exprAssembler)
    {
        gsMultiPatch<real_t> mp;
        gsMultiBasis<real_t> mb = traceTestPatches(mp);
        // Non-zero on both Dirichlet sides (x = 0 and x = 2), so the
        // interpolated values are non-trivial
        gsFunctionExpr<real_t> g("1 + x*y", 2);
        gsBoundaryConditions<real_t> E = traceTestConditions(g, mp, trace_type::essential);
        gsBoundaryConditions<real_t> A = traceTestConditions(g, mp, trace_type::all);
        gsBoundaryConditions<real_t> N = traceTestConditions(g, mp, trace_type::normal);

        gsExprAssembler<real_t> eaE(1,1), eaA(1,1), eaN(1,1);
        eaE.setIntegrationElements(mb);
        eaA.setIntegrationElements(mb);
        eaN.setIntegrationElements(mb);
        auto uE = eaE.getSpace(mb);
        auto uA = eaA.getSpace(mb);
        auto uN = eaN.getSpace(mb);
        uE.setup(E, dirichlet::interpolation, 0);
        uA.setup(A, dirichlet::interpolation, 0);
        CHECK(uE.mapper().boundarySize() > 0);
        checkSameMapper(uE.mapper(), uA.mapper());

        // The interpolation does not read the trace, so the same operations
        // on the same data give the same values, bit for bit
        const gsMatrix<real_t> & fE = uE.fixedPart();
        const gsMatrix<real_t> & fA = uA.fixedPart();
        CHECK_EQUAL(fE.rows(), fA.rows());
        CHECK_EQUAL(fE.cols(), fA.cols());
        CHECK(fE.rows() > 0);
        CHECK(fE.cwiseAbs().minCoeff() > 0);
        if (fE.rows() == fA.rows() && fE.cols() == fA.cols())
            CHECK(fE == fA);

        // Normal is refused by the setup, and by gsDirichletValues on a
        // space set up without it
        CHECK_THROW(uN.setup(N, dirichlet::interpolation, 0), std::runtime_error);
        CHECK_THROW(gsDirichletValues(N, dirichlet::interpolation, uE), std::runtime_error);
    }

    // Property: gsPoissonAssembler builds the same dof mapper with the
    // traces essential and all, and refuses tangential.
    TEST(trace_poissonAssembler)
    {
        gsMultiPatch<real_t> mp;
        gsMultiBasis<real_t> mb = traceTestPatches(mp);
        gsFunctionExpr<real_t> g("x*y", 2), rhs("1", 2);
        gsBoundaryConditions<real_t> E = traceTestConditions(g, mp, trace_type::essential);
        gsBoundaryConditions<real_t> A = traceTestConditions(g, mp, trace_type::all);
        gsBoundaryConditions<real_t> Tn = traceTestConditions(g, mp, trace_type::tangential);

        // Same number of dofs and same numbering for essential and all
        gsPoissonAssembler<real_t> pE(mp, mb, E, rhs);
        gsPoissonAssembler<real_t> pA(mp, mb, A, rhs);
        CHECK_EQUAL(pE.numDofs(), pA.numDofs());
        const gsDofMapper & mE = pE.system().colMapper(0);
        const gsDofMapper & mA = pA.system().colMapper(0);
        CHECK(mE.boundarySize() > 0);
        checkSameMapper(mE, mA);

        CHECK_THROW(gsPoissonAssembler<real_t> pT(mp, mb, Tn, rhs), std::runtime_error);
    }
}
