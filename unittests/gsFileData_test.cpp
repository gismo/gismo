/** @file gsFileData_test.cpp

    @brief Typed XML lookups in gsFileData in the presence of nodes without a type attribute

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include "gismo_unittest.h"

SUITE(gsFileData_test)
{
    typedef gsTensorBSpline<2, real_t> Geo;

    // A typed geometry whose second coefficient is (c0+3, 0)
    std::string geoNode(int id, int c0)
    {
        std::ostringstream os;
        os << "<Geometry type=\"TensorBSpline2\" id=\"" << id << "\">\n"
           << " <Basis type=\"TensorBSplineBasis2\">\n"
           << "  <Basis type=\"BSplineBasis\" index=\"0\"><KnotVector degree=\"1\">0 0 1 1</KnotVector></Basis>\n"
           << "  <Basis type=\"BSplineBasis\" index=\"1\"><KnotVector degree=\"1\">0 0 1 1</KnotVector></Basis>\n"
           << " </Basis>\n"
           << " <coefs geoDim=\"2\">" << c0 << " 0  " << c0 + 3 << " 0  "
           << c0 << " 5  " << c0 + 3 << " 5</coefs>\n"
           << "</Geometry>\n";
        return os.str();
    }

    const std::string nodeA = geoNode(7, 0);   // coefs()(1,0) == 3
    const std::string nodeB = geoNode(8, 10);  // coefs()(1,0) == 13

    // Writes <xml>body</xml> to a temporary file, reads it into fd and removes it
    void load(gsFileData<> & fd, const std::string & name, const std::string & body)
    {
        const std::string tmp = gsFileManager::getTempPath();
        CHECK(!tmp.empty());
        const std::string fn = tmp + name + ".xml";
        {
            std::ofstream os(fn.c_str());
            os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<xml>\n" << body << "</xml>\n";
        }
        CHECK(fd.read(fn));
        std::remove(fn.c_str());
    }

    void checkIsA(const Geo & g)
    {
        CHECK_EQUAL(4, g.coefs().rows());
        CHECK_EQUAL(3, g.coefs()(1, 0));
    }

    void checkAnyFindsA(const std::string & name, const std::string & body)
    {
        gsFileData<> fd;
        load(fd, name, body);
        CHECK(fd.hasAny<Geo>());
        memory::unique_ptr<Geo> g = fd.getAnyFirst<Geo>();
        CHECK(g.get() != NULL);
        if (g) checkIsA(*g);
        Geo obj;
        CHECK(fd.getAnyFirst(obj));
        checkIsA(obj);
    }

    TEST(fdxAnyLevel1)
    {
        checkAnyFindsA("gsFileData_test_any1", "<Geometry id=\"99\"/>\n" + nodeA);
    }

    TEST(fdxAnyLevel2)
    {
        checkAnyFindsA("gsFileData_test_any2",
                       "<Wrap1><Geometry id=\"99\"/></Wrap1>\n" + nodeA);
    }

    TEST(fdxAnyLevel3)
    {
        checkAnyFindsA("gsFileData_test_any3",
                       "<Wrap1><Wrap2><Geometry id=\"99\"/></Wrap2></Wrap1>\n" + nodeA);
    }

    TEST(fdxTopLevel)
    {
        gsFileData<> fd;
        load(fd, "gsFileData_test_top", "<Geometry id=\"99\"/>\n" + nodeA);
        CHECK(fd.has<Geo>());
        memory::unique_ptr<Geo> g = fd.getFirst<Geo>();
        CHECK(g.get() != NULL);
        if (g) checkIsA(*g);
        Geo obj;
        CHECK(fd.getFirst(obj));
        checkIsA(obj);
        CHECK_EQUAL(1, fd.count<Geo>());
    }

    TEST(fdxSiblings)
    {
        gsFileData<> fd;
        load(fd, "gsFileData_test_sib", nodeA + "<Geometry id=\"99\"/>\n" + nodeB);
        CHECK_EQUAL(2, fd.count<Geo>());
        std::vector<memory::unique_ptr<Geo> > all = fd.getAll<Geo>();
        CHECK_EQUAL((size_t)2, all.size());
        if (all.size() == 2)
        {
            CHECK_EQUAL(3, all[0]->coefs()(1, 0));
            CHECK_EQUAL(13, all[1]->coefs()(1, 0));
        }
    }

    TEST(fdxOnlyUntyped)
    {
        gsFileData<> fd;
        load(fd, "gsFileData_test_untyped",
             "<Geometry id=\"98\"/>\n<Wrap1><Wrap2><Geometry id=\"99\"/></Wrap2></Wrap1>\n");
        CHECK(!fd.hasAny<Geo>());
        CHECK(fd.getAnyFirst<Geo>().get() == NULL);
        CHECK(!fd.has<Geo>());
        CHECK(fd.getFirst<Geo>().get() == NULL);
        CHECK_EQUAL(0, fd.count<Geo>());
        CHECK(fd.getAll<Geo>().empty());
    }
}
