/** @file gsModuleManager_test.cpp

    @brief Tests gsRegistry, gsModuleManager, gsFileIo and the 3dm file
    format served by the gsOpennurbs module (built in or loaded dynamically)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Angelos Mantzaflaris
**/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework

// gsRegistry<Base>::get() is defined in the .hpp, which a GISMO_BUILD_LIB
// build does not include from the .h; the test-local registry needs it.
#include <gsSystem/gsRegistry.hpp>
#include <gsSystem/gsModuleConfig.h>

#include <algorithm>
#include <cstdio>
#include <fstream>

using namespace gismo;

SUITE(gsModuleManager_test)
{
namespace
{

// Test-local hierarchy for the registry tests
struct MMBase
{
    virtual ~MMBase() {}
    virtual int id() const = 0;
};

struct MMOne : public MMBase { int id() const { return 1; } };
struct MMTwo : public MMBase { int id() const { return 2; } };

typedef gsRegistry<MMBase> MMRegistry;

// true if @a v holds @a s
bool contains(const std::vector<std::string> & v, const std::string & s)
{
    return std::find(v.begin(), v.end(), s) != v.end();
}

// Position of @a s in @a v, v.size() if absent
size_t position(const std::vector<std::string> & v, const std::string & s)
{
    return std::find(v.begin(), v.end(), s) - v.begin();
}

// true when an io for ".3dm" can be made; otherwise prints why the caller returns early
bool have3dm(const char * test)
{
    gsFileIo::uPtr io = gsFileIo::find("3dm");   // dynamic mode: gsOpennurbs was loaded with the library
    if (io) return true;
    gsInfo << "gsModuleManager_test: 3dm io not available, skipping " << test << "\n";
    return false;
}

// Cubic B-spline curve with 6 basis functions (two interior knots) in
// dimension @a dim; coefficient row i is (i, 0.25 i^2, 1 - 0.5 i) cut to dim columns.
gsBSpline<> makeCurve(index_t dim)
{
    gsKnotVector<> kv(0.0, 1.0, 2, 4, 1, 3);
    gsMatrix<> c(6, 3);
    for (index_t i = 0; i < 6; ++i)
        c.row(i) << i, 0.25 * i * i, 1.0 - 0.5 * i;
    gsMatrix<> cd = c.leftCols(dim);
    return gsBSpline<>(kv, cd);
}

// Tensor B-spline surface with 4 x 3 basis functions of degrees (2,1);
// coefficient row i + 4 j is (i, j, 0.1 i + 0.3 j^2 + 0.05 i j).
gsTensorBSpline<2> makeSurface()
{
    gsKnotVector<> ku(0.0, 1.0, 1, 3, 1, 2);
    gsKnotVector<> kv(0.0, 1.0, 1, 2, 1, 1);
    gsMatrix<> c(12, 3);
    for (index_t j = 0; j < 3; ++j)
        for (index_t i = 0; i < 4; ++i)
            c.row(i + 4 * j) << i, j, 0.1 * i + 0.3 * j * j + 0.05 * i * j;
    return gsTensorBSpline<2>(ku, kv, c);
}

// Writes a curve and a surface to a .3dm file, reads it back and compares
// degrees, knots and coefficients. The tolerance EPSILON (about 5.6e-12 for
// double) exceeds the round-off of the double values the writer stores and
// the reader prints with 16 significant digits (relative error about 5e-16,
// values of magnitude at most 6.25).
void roundTripCurveAndSurface()
{
    // Reference: the objects themselves; nothing is recomputed from the file.
    gsBSpline<> crv = makeCurve(3);
    gsTensorBSpline<2> srf = makeSurface();

    // Write both objects into one 3dm file.
    const std::string fn = gsFileManager::getTempPath() + "gsModuleManager_test_roundtrip.3dm";
    {
        gsFileData<> out;
        out << crv;
        out << srf;
        out.save(fn);
    }
    CHECK(gsFileManager::fileExists(fn));

    // Read it back; both kinds must be present.
    gsFileData<> in(fn);
    CHECK(in.has< gsBSpline<> >());
    CHECK(in.has< gsTensorBSpline<2> >());

    // Curve: degree, dimension, size, knots and coefficients are preserved.
    gsBSpline<>::uPtr c2 = in.getFirst< gsBSpline<> >();
    CHECK(c2.get() != NULL);
    if (c2)
    {
        CHECK_EQUAL(crv.degree(), c2->degree());
        CHECK_EQUAL(3, c2->geoDim());
        CHECK_EQUAL(crv.coefsSize(), c2->coefsSize());
        gsMatrix<real_t> k1 = crv.knots().asMatrix();
        gsMatrix<real_t> k2 = c2->knots().asMatrix();
        CHECK_MATRIX_CLOSE(k1, k2, EPSILON);
        gsMatrix<real_t> a1 = crv.coefs();
        gsMatrix<real_t> a2 = c2->coefs();
        CHECK_MATRIX_CLOSE(a1, a2, EPSILON);
    }

    // Surface: sizes and degrees differ per direction, so a transposed
    // read-back cannot pass.
    gsTensorBSpline<2>::uPtr s2 = in.getFirst< gsTensorBSpline<2> >();
    CHECK(s2.get() != NULL);
    if (s2)
    {
        CHECK_EQUAL(srf.degree(0), s2->degree(0));
        CHECK_EQUAL(srf.degree(1), s2->degree(1));
        CHECK_EQUAL(3, s2->geoDim());
        CHECK_EQUAL(srf.coefsSize(), s2->coefsSize());
        for (short_t d = 0; d < 2; ++d)
        {
            gsMatrix<real_t> k1 = srf.knots(d).asMatrix();
            gsMatrix<real_t> k2 = s2->knots(d).asMatrix();
            CHECK_MATRIX_CLOSE(k1, k2, EPSILON);
        }
        gsMatrix<real_t> a1 = srf.coefs();
        gsMatrix<real_t> a2 = s2->coefs();
        CHECK_MATRIX_CLOSE(a1, a2, EPSILON);
    }

    std::remove(fn.c_str());
}

} // namespace

// Property: a new key is stored, found and produces the registered type.
// Expected values: the ids 1 and 2 returned by MMOne and MMTwo by construction.
TEST(Registry_addHasMake)
{
    MMRegistry & r = MMRegistry::get();

    // Register two keys of this test.
    CHECK(r.add("mm_add_a", &gsRegistryFactory<MMBase, MMOne>));
    CHECK(r.add("mm_add_b", &gsRegistryFactory<MMBase, MMTwo>));
    CHECK(r.has("mm_add_a"));
    CHECK(r.has("mm_add_b"));

    // Each key makes the type it was registered with.
    MMRegistry::BasePtr a = r.make("mm_add_a");
    MMRegistry::BasePtr b = r.make("mm_add_b");
    CHECK(a.get() != NULL);
    CHECK(b.get() != NULL);
    if (a) CHECK_EQUAL(1, a->id());
    if (b) CHECK_EQUAL(2, b->id());

    // A key that was never added is absent and makes a null pointer.
    CHECK(!r.has("mm_add_never"));
    CHECK(r.make("mm_add_never").get() == NULL);

    // keys() lists the added keys.
    const std::vector<std::string> k = r.keys();
    CHECK(contains(k, "mm_add_a"));
    CHECK(contains(k, "mm_add_b"));
}

// Property: the first registration of a key wins and a null factory is
// refused. Expected values: the documented contract of gsRegistry::add.
TEST(Registry_firstWins)
{
    MMRegistry & r = MMRegistry::get();

    // A second registration under the same key is refused and changes nothing.
    CHECK(r.add("mm_first", &gsRegistryFactory<MMBase, MMOne>));
    CHECK(!r.add("mm_first", &gsRegistryFactory<MMBase, MMTwo>));
    MMRegistry::BasePtr p = r.make("mm_first");
    CHECK(p.get() != NULL);
    if (p) CHECK_EQUAL(1, p->id());

    // A null factory is refused and leaves no entry.
    CHECK(!r.add("mm_null", NULL));
    CHECK(!r.has("mm_null"));
}

// Property: removeModule erases exactly the entries of the named owner.
// Expected values: entries added from a test have owner "" (no register
// function is running), so removeModule("") removes all entries of this
// suite's registry, removeModule of another name none.
TEST(Registry_removeModule)
{
    MMRegistry & r = MMRegistry::get();

    // Outside a register function the owner is empty.
    CHECK(gsModuleManager::get().registering().empty());
    CHECK(r.add("mm_rm_a", &gsRegistryFactory<MMBase, MMOne>));
    CHECK(r.add("mm_rm_b", &gsRegistryFactory<MMBase, MMTwo>));

    // An unknown module owns nothing: nothing is removed.
    CHECK_EQUAL(0u, r.removeModule("gsNoSuchModule"));
    CHECK(r.has("mm_rm_a"));
    CHECK(r.has("mm_rm_b"));

    // Owner "": every entry of this registry goes.
    const size_t n = r.keys().size();
    CHECK_EQUAL(n, r.removeModule(""));
    CHECK(r.keys().empty());
    CHECK(r.make("mm_rm_a").get() == NULL);
    CHECK(r.make("mm_rm_b").get() == NULL);

    // A removed key is free again.
    CHECK(r.add("mm_rm_a", &gsRegistryFactory<MMBase, MMOne>));
}

// Property: a path to no file is refused, recorded in rejected() and never
// throws. Expected values: the documented behaviour of gsModuleManager::load.
TEST(Manager_loadMissingPath)
{
    gsModuleManager & mm = gsModuleManager::get();
    const std::string p = gsFileManager::getTempPath() + "gsModuleManager_test_missing.so";

    // Precondition: the file really does not exist.
    CHECK(!gsFileManager::fileExists(p));

    // Loading fails, and the path is listed as rejected.
    CHECK(!mm.load(p));
    CHECK(contains(mm.rejected(), p));
    CHECK(!mm.has("gsModuleManager_test_missing"));
    CHECK(!mm.isLoaded("gsModuleManager_test_missing"));
}

// Property: an existing file that is not a shared library is refused (dlopen
// fails) and recorded in rejected(). Expected values: documented behaviour.
TEST(Manager_loadNotALibrary)
{
    gsModuleManager & mm = gsModuleManager::get();

    // A text file with a library-like name; the temp directory is not a search
    // path, so no later loadAll() touches it.
    const std::string p = gsFileManager::getTempPath() + "gsModuleManager_test_notalib.so";
    {
        std::ofstream f(p.c_str());
        f << "this is not an ELF file\n";
    }
    CHECK(gsFileManager::fileExists(p));

    // Loading fails, and the path is listed as rejected.
    CHECK(!mm.load(p));
    CHECK(contains(mm.rejected(), p));

    std::remove(p.c_str());
}

// Property: a bare name found in no search directory is refused and recorded
// under the name; unload of it is false. Expected values: documented behaviour.
TEST(Manager_loadUnknownName)
{
    gsModuleManager & mm = gsModuleManager::get();

    // The unknown name is refused, listed under the name, and cannot be unloaded.
    CHECK(!mm.load("gsNoSuchModule"));
    CHECK(contains(mm.rejected(), "gsNoSuchModule"));
    CHECK(!mm.unload("gsNoSuchModule"));
    CHECK(!mm.has("gsNoSuchModule"));
}

// Property: addSearchPath keeps order, ignores empty strings and duplicates,
// and the configured module directory follows the added ones.
// Expected values: the order documented for searchPaths().
TEST(Manager_searchPaths)
{
    gsModuleManager & mm = gsModuleManager::get();

    // A directory that does not exist; never the temp directory itself, as
    // loadAll() would dlopen the library files lying there.
    const std::string d = gsFileManager::getTempPath() + "gsModuleManager_test_nodir";
    mm.addSearchPath(d);
    mm.addSearchPath("");

    // d is listed, before the built-in module directory, and no entry is empty.
    const std::vector<std::string> sp = mm.searchPaths();
    CHECK(contains(sp, d));
    CHECK(contains(sp, GISMO_MODULE_DIR));
    CHECK(position(sp, d) < position(sp, GISMO_MODULE_DIR));
    CHECK(!contains(sp, ""));

    // Adding d again changes nothing.
    mm.addSearchPath(d);
    CHECK_EQUAL(sp.size(), mm.searchPaths().size());
}

// Property: gsOpennurbs is available (3dm format) and the manager reports it
// as built in or as loaded dynamically consistently with the build.
// Expected values: with gsOpennurbs_ENABLED the module is compiled into
// libgismo (has, not isLoaded); without it an available module can only have
// been loaded from a shared library (isLoaded).
TEST(Manager_hasOpennurbs)
{
    if (!have3dm("Manager_hasOpennurbs")) return;
    gsModuleManager & mm = gsModuleManager::get();

    // The module is known and listed among the modules.
    CHECK(mm.has("gsOpennurbs"));
    CHECK(contains(mm.modules(), "gsOpennurbs"));

    // Mode as the expected answer of the build configuration.
#ifdef gsOpennurbs_ENABLED
    CHECK(!mm.isLoaded("gsOpennurbs"));
#else
    CHECK(mm.isLoaded("gsOpennurbs"));
#endif

    // The log shows which mode this run exercised.
    if (mm.isLoaded("gsOpennurbs"))
        gsInfo << "gsModuleManager_test: gsOpennurbs is loaded dynamically\n";
    else
        gsInfo << "gsModuleManager_test: gsOpennurbs is built in\n";
}

// Property: find normalises the extension and returns null for none or an
// unknown one; the 3dm io can read and write and is registered under "3dm".
// Expected values: documented behaviour of gsFileIo::find.
TEST(FileIo_find)
{
    // Holds in every build: no extension, and an extension nobody registers.
    CHECK(gsFileIo::find("").get() == NULL);
    CHECK(gsFileIo::find("gsmmtest_noext").get() == NULL);

    if (!have3dm("FileIo_find")) return;

    // A leading dot and upper case are dropped.
    gsFileIo::uPtr io = gsFileIo::find(".3DM");
    CHECK(io.get() != NULL);
    if (io)
    {
        CHECK(io->canRead());
        CHECK(io->canWrite());
    }
    CHECK(gsRegistry<gsFileIo>::get().has("3dm"));
}

// Property: a 3-D cubic curve and a 3-D bi-(2,1) surface survive a write and
// read of a .3dm file. Reference: the in-memory objects (analytic
// coefficients), see roundTripCurveAndSurface for the tolerance.
TEST(FileIo_roundTrip3dm)
{
    if (!have3dm("FileIo_roundTrip3dm")) return;
    roundTripCurveAndSurface();
}

// Property: a planar curve is stored as a 3-D curve with zero third
// coordinate. Expected values: the input's two columns and zeros, since the
// 3dm writer always makes 3-D objects and pads with 0.
TEST(FileIo_roundTrip3dmPlanarCurve)
{
    if (!have3dm("FileIo_roundTrip3dmPlanarCurve")) return;

    // A planar curve alone in the file.
    gsBSpline<> crv = makeCurve(2);
    const std::string fn = gsFileManager::getTempPath() + "gsModuleManager_test_planar.3dm";
    {
        gsFileData<> out;
        out << crv;
        out.save(fn);
    }
    CHECK(gsFileManager::fileExists(fn));

    gsFileData<> in(fn);
    gsBSpline<>::uPtr c2 = in.getFirst< gsBSpline<> >();
    CHECK(c2.get() != NULL);
    if (c2)
    {
        // Third dimension added; the first two columns are the input.
        CHECK_EQUAL(3, c2->geoDim());
        gsMatrix<real_t> a1 = crv.coefs();
        gsMatrix<real_t> a2 = c2->coefs().leftCols(2);
        CHECK_MATRIX_CLOSE(a1, a2, EPSILON);
        // Zero column, exact by construction; EPSILON for round-off safety.
        gsMatrix<real_t> z = c2->coefs().col(2);
        gsMatrix<real_t> zero = gsMatrix<real_t>::Zero(6, 1);
        CHECK_MATRIX_CLOSE(zero, z, EPSILON);
    }

    std::remove(fn.c_str());
}

// Property: a dynamically loaded module can be unloaded (its registrations
// disappear) and can be loaded again; a built-in module refuses unload.
// Expected values: documented behaviour of gsModuleManager::unload and of
// gsFileIo::find (which loads nothing). Last test of the suite, as it changes
// the module state.
TEST(Manager_unloadReload)
{
    if (!have3dm("Manager_unloadReload")) return;
    gsModuleManager & mm = gsModuleManager::get();

    if (!mm.isLoaded("gsOpennurbs"))
    {
        // Built in: unload refuses and nothing changes.
        CHECK(!mm.unload("gsOpennurbs"));
        CHECK(mm.has("gsOpennurbs"));
        CHECK(gsFileIo::find("3dm").get() != NULL);
        return;
    }

    // Dynamic: no gsFileIo object is alive here (have3dm's was destroyed).
    CHECK(mm.unload("gsOpennurbs"));
    CHECK(!mm.isLoaded("gsOpennurbs"));
    CHECK(!mm.has("gsOpennurbs"));
    CHECK(!gsRegistry<gsFileIo>::get().has("3dm"));
    CHECK(!mm.unload("gsOpennurbs"));

    // A lookup of the extension does not load the module again.
    CHECK(gsFileIo::find("3dm").get() == NULL);
    CHECK(!mm.isLoaded("gsOpennurbs"));

    // An explicit load does; the io object is destroyed before the scope ends.
    CHECK(mm.load("gsOpennurbs"));
    {
        gsFileIo::uPtr io = gsFileIo::find("3dm");
        CHECK(io.get() != NULL);
    }
    CHECK(mm.isLoaded("gsOpennurbs"));
    CHECK(gsRegistry<gsFileIo>::get().has("3dm"));

    // The reloaded module writes and reads.
    roundTripCurveAndSurface();
}

} // SUITE(gsModuleManager_test)
