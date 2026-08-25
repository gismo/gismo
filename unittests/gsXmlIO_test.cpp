/** @file gsXmlIO_test.cpp

    @brief Tests writer/reader round-trip of labelled boundary groups
    (gsXml<gsBoundaryConditions>::get_into's name="..." selection) and of
    gsFileData::reserveIds.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G+Smo agent
**/

#include "gismo_unittest.h"
#include <fstream>
#include <set>
#include <tuple>

using namespace gismo;

namespace
{

/// Two BSplineSquares sharing the edge x==1 (patch0 east / patch1 west),
/// with topology computed. Every side with side()==boundary::west gets the
/// label "west", every other boundary side gets "other" -- so exactly one
/// side (patch 0, west) carries "west", the remaining boundary sides carry
/// "other" (at least two distinct label groups).
gsMultiPatch<real_t> labelledFixture()
{
    gsMultiPatch<real_t> mp;
    mp.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1.0, 0.0, 0.0));
    mp.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1.0, 1.0, 0.0));
    mp.computeTopology();

    for (gsBoxTopology::biterator it = mp.bBegin(); it != mp.bEnd(); ++it)
    {
        if (it->side() == boundary::west)
            it->setLabel("west");
        else
            it->setLabel("other");
    }
    return mp;
}

typedef std::tuple<index_t,int,std::string> LabelledSide;

std::multiset<LabelledSide> labelledSides(const gsBoxTopology & topo)
{
    std::multiset<LabelledSide> result;
    for (gsBoxTopology::const_biterator it = topo.bBegin(); it != topo.bEnd(); ++it)
        result.insert( LabelledSide(it->patch, it->side().index(), it->label()) );
    return result;
}

std::multiset<std::pair<index_t,int> > unlabelledSides(const gsBoxTopology & topo)
{
    std::multiset<std::pair<index_t,int> > result;
    for (gsBoxTopology::const_biterator it = topo.bBegin(); it != topo.bEnd(); ++it)
        result.insert( std::make_pair(it->patch, it->side().index()) );
    return result;
}

std::string readWholeFile(const std::string & path)
{
    std::ifstream in(path.c_str());
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::string tmpPath(const std::string & name)
{
    return gsFileManager::getTempPath()
        + gsFileManager::getNativePathSeparator() + name;
}

} // anonymous namespace

SUITE(gsXmlIO_test)
{

// Label round-trip. FAILS on unmodified master: appendBoxTopology dropped
// labels, so every read-back side comes back with an empty label.
TEST(label_roundtrip)
{
    gsMultiPatch<real_t> mp = labelledFixture();
    std::multiset<LabelledSide> expected = labelledSides(mp);

    std::string path = tmpPath("gsXmlIO_test_label_roundtrip.xml");
    gsWrite(mp, path);

    gsMultiPatch<real_t> mp2;
    gsReadFile<real_t>(path, mp2);

    CHECK( expected == labelledSides(mp2) );
}

// Boundary patch indices under a non-zero id base (reserveIds(500)).
// FAILS on unmodified master: the boundary branch of appendBoxTopology
// writes the raw (un-mapped) patch index, so every boundary reads back
// on patch 0.
TEST(boundary_patch_indices_survive_id_offset)
{
    gsMultiPatch<real_t> mp = labelledFixture();
    std::set<index_t> expectedPatches;
    for (gsBoxTopology::const_biterator it = mp.bBegin(); it != mp.bEnd(); ++it)
        expectedPatches.insert(it->patch);
    CHECK( expectedPatches.size() > 1 ); // fixture sanity: >1 patch has a boundary side

    std::string path = tmpPath("gsXmlIO_test_id_offset.xml");
    gsFileData<real_t> fd;
    fd.reserveIds(500);
    fd.add(mp, 0);
    fd.save(path);

    gsFileData<real_t> rd(path);
    gsMultiPatch<real_t> mp2;
    rd.getId(0, mp2);

    std::set<index_t> actualPatches;
    for (gsBoxTopology::const_biterator it = mp2.bBegin(); it != mp2.bEnd(); ++it)
        actualPatches.insert(it->patch);

    CHECK( expectedPatches == actualPatches );
}

// Id layout: after reserveIds(500), id 0 is free for the MultiPatch, the
// saved <patches type="id_range"> starts at 500, and no stray <string> node
// (the reserveIds bookkeeping dummy) survives in the saved file.
TEST(reserveIds_id_layout)
{
    gsMultiPatch<real_t> mp = labelledFixture();

    std::string path = tmpPath("gsXmlIO_test_id_layout.xml");
    gsFileData<real_t> fd;
    fd.reserveIds(500);
    fd.add(mp, 0);
    fd.save(path);

    gsFileData<real_t> rd(path);
    CHECK( rd.hasId(0) );
    gsMultiPatch<real_t> mp2;
    rd.getId(0, mp2);
    CHECK_EQUAL( mp.nPatches(), mp2.nPatches() );

    std::string text = readWholeFile(path);
    CHECK( text.find("<patches type=\"id_range\">500 ") != std::string::npos );
    CHECK( text.find("<string") == std::string::npos );
}

// reserveIds must refuse to go backwards but must be idempotent at the
// current boundary.
TEST(reserveIds_backwards_throws_forwards_idempotent)
{
    gsFileData<real_t> fd;
    gsMultiPatch<real_t> mp = labelledFixture();
    fd.reserveIds(500);
    fd.add(mp, 0); // consumes id 0, geometries/topology take further auto ids

    // Going below the current max must throw (GISMO_ENSURE -> std::runtime_error).
    CHECK_THROW( fd.reserveIds(1), std::runtime_error );

    // On a fresh tree, calling reserveIds(500) twice must succeed: the
    // second call is a no-op (max_Id already == 499), not a throw.
    gsFileData<real_t> fd2;
    fd2.reserveIds(500);
    bool threw = false;
    try { fd2.reserveIds(500); }
    catch (const std::exception &) { threw = true; }
    CHECK( !threw );

    // The next auto-assigned id after the idempotent pair is still 500.
    gsMultiPatch<real_t> mp2 = labelledFixture();
    fd2.add(mp2, 0);
    std::string path = tmpPath("gsXmlIO_test_reserveIds_idempotent.xml");
    fd2.save(path);
    std::string text = readWholeFile(path);
    CHECK( text.find("<patches type=\"id_range\">500 ") != std::string::npos );
}

// End-to-end: a <bc name="west"> picks exactly the "west"-labelled sides,
// and none of the "other"-labelled sides.
TEST(bc_selects_by_label)
{
    gsMultiPatch<real_t> mp = labelledFixture();

    std::set<std::pair<index_t,int> > westSides, otherSides;
    for (gsBoxTopology::const_biterator it = mp.bBegin(); it != mp.bEnd(); ++it)
    {
        if (it->label() == "west")
            westSides.insert( std::make_pair(it->patch, it->side().index()) );
        else
            otherSides.insert( std::make_pair(it->patch, it->side().index()) );
    }
    CHECK( westSides.size() >= 1 );
    CHECK( otherSides.size() >= 1 );

    std::string path = tmpPath("gsXmlIO_test_bc_by_label.xml");
    gsFileData<real_t> fd;
    fd.reserveIds(500);
    fd.add(mp, 0);
    fd.save(path);

    // Inject a <boundaryConditions> block selecting the "west" label,
    // just before the closing </xml> tag.
    std::string text = readWholeFile(path);
    const std::string bcBlock =
        "<boundaryConditions id=\"2\" multipatch=\"0\">\n"
        "<Function type=\"FunctionExpr\" dim=\"2\" index=\"0\">0</Function>\n"
        "<bc type=\"Dirichlet\" function=\"0\" unknown=\"0\" name=\"west\"></bc>\n"
        "</boundaryConditions>\n";
    const std::string closing = "</xml>";
    size_t pos = text.rfind(closing);
    CHECK( pos != std::string::npos );
    text.insert(pos, bcBlock);

    std::ofstream out(path.c_str());
    out << text;
    out.close();

    gsFileData<real_t> rd(path);
    gsBoundaryConditions<real_t> bc;
    rd.getId(2, bc);

    std::set<std::pair<index_t,int> > selected;
    typedef gsBoundaryConditions<real_t>::bcContainer bcContainer;
    bcContainer all = bc.allConditions();
    for (bcContainer::const_iterator it = all.begin(); it != all.end(); ++it)
        selected.insert( std::make_pair(it->patch(), it->side().index()) );

    CHECK( westSides == selected );
    // discriminating half: none of the "other" sides were picked
    for (std::set<std::pair<index_t,int> >::const_iterator it = otherSides.begin();
         it != otherSides.end(); ++it)
        CHECK( selected.find(*it) == selected.end() );
}

// No churn on unlabelled input: same interfaces, same boundaries, same
// patch indices as before this task's changes.
TEST(no_churn_on_unlabelled_input)
{
    gsMultiPatch<real_t> mp;
    mp.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1.0, 0.0, 0.0));
    mp.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1.0, 1.0, 0.0));
    mp.computeTopology();
    // labels intentionally left empty

    std::multiset<std::pair<index_t,int> > expectedBoundaries = unlabelledSides(mp);
    const size_t expectedInterfaces = mp.nInterfaces();

    std::string path = tmpPath("gsXmlIO_test_no_churn.xml");
    gsWrite(mp, path);

    gsMultiPatch<real_t> mp2;
    gsReadFile<real_t>(path, mp2);

    CHECK( expectedBoundaries == unlabelledSides(mp2) );
    CHECK_EQUAL( expectedInterfaces, mp2.nInterfaces() );
    CHECK_EQUAL( mp.nPatches(), mp2.nPatches() );
}

// gsMultiBasis labelled round-trip: the Edit-4 fail-before. With Edit 1
// applied but Edit 4 omitted, this fails by COUNT -- every <boundary>
// group after the first is silently dropped by the reader's single-node
// read. (On the fully unmodified library it passes vacuously: one group,
// no labels ever written.)
TEST(multibasis_labelled_roundtrip)
{
    gsMultiPatch<real_t> mp = labelledFixture();
    gsMultiBasis<real_t> mb(mp);

    std::multiset<std::pair<index_t,int> > expected = unlabelledSides(mb.topology());
    // fixture sanity: at least two distinct labels present
    std::set<std::string> labels;
    for (gsBoxTopology::const_biterator it = mb.topology().bBegin();
         it != mb.topology().bEnd(); ++it)
        labels.insert(it->label());
    CHECK( labels.size() >= 2 );

    std::string path = tmpPath("gsXmlIO_test_multibasis.xml");
    gsWrite(mb, path);

    gsMultiBasis<real_t> mb2;
    gsReadFile<real_t>(path, mb2);

    CHECK( expected == unlabelledSides(mb2.topology()) );
}

} // SUITE(gsXmlIO_test)
