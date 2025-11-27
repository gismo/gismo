/** @file gsIO_test.cpp

    @brief Comprehensive tests for gsIO module (input/output operations)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Copilot
*/

#include "gismo_unittest.h"

SUITE(gsIO_test)
{
    TEST(gsIO_BSpline_XML)
    {
        // Create a simple B-spline curve
        gsKnotVector<> kv(0, 1, 1, 3);
        gsBSplineBasis<> basis(kv);
        gsMatrix<> coefs(3, 2);
        coefs << 0, 0,
                 0.5, 1,
                 1, 0;
        gsBSpline<> curve(kv, coefs);

        // Write to XML
        std::string filename = "/tmp/test_curve.xml";
        gsWrite(curve, filename);

        // Verify file exists
        std::ifstream file(filename);
        CHECK(file.good());
        file.close();

        // Read it back
        gsFileData<> fd(filename);
        auto geom = fd.getFirst<gsGeometry<>>();
        
        CHECK(geom.get() != nullptr);
        CHECK_EQUAL(geom->domainDim(), 1);
    }

    TEST(gsWrite_MultiPatch_XML)
    {
        gsKnotVector<> kv(0, 1, 1, 2);
        gsBSplineBasis<> basis(kv);
        gsMatrix<> coefs(2, 2);
        coefs << 0, 0,
                 1, 0;
        gsBSpline<> curve(kv, coefs);
        
        gsMultiPatch<> mp;
        mp.addPatch(curve);
        mp.addPatch(curve);
        
        std::string filename = "/tmp/test_multipatch.xml";
        gsWrite(mp, filename);
        
        std::ifstream file(filename);
        CHECK(file.good());
        file.close();
    }

    TEST(gsIO_Matrix_XML)
    {
        // Write a matrix to XML
        gsMatrix<> mat(3, 4);
        mat << 1, 2, 3, 4,
               5, 6, 7, 8,
               9, 10, 11, 12;
        std::string filename = "/tmp/test_matrix.xml";
        gsWrite(mat, filename);

        // Verify file exists
        std::ifstream file(filename);
        CHECK(file.good());
        file.close();

        // Read the matrix back
        gsFileData<> fd(filename);
        auto readMat = fd.getFirst<gsMatrix<>>();
        CHECK_EQUAL(readMat->rows(), 3);
        CHECK_EQUAL(readMat->cols(), 4);
        CHECK_CLOSE((*readMat)(0, 0), 1, 1e-10);
    }

    TEST(gsFileData_HasId)
    {
        gsMatrix<> mat(2, 2);
        mat << 1, 2, 3, 4;
        
        std::string filename = "/tmp/test_hasid.xml";
        gsWrite(mat, filename);
        
        gsFileData<> fd(filename);
        CHECK(fd.has<gsMatrix<>>());
    }

    TEST(gsFileData_Count)
    {
        gsKnotVector<> kv(0, 1, 1, 2);
        gsBSplineBasis<> basis(kv);
        gsMatrix<> coefs(2, 1);
        coefs << 0, 1;
        gsBSpline<> curve(kv, coefs);
        
        gsMultiPatch<> mp;
        mp.addPatch(curve);
        mp.addPatch(curve);
        mp.addPatch(curve);
        
        std::string filename = "/tmp/test_count.xml";
        gsWrite(mp, filename);
        
        gsFileData<> fd(filename);
        int count = fd.count<gsGeometry<>>();
        CHECK(count > 0);
    }

    /* 
     * Step-by-step instructions for additional complex gsIO tests:
     * 
     * 1. gsWriteParaview tests:
     *    - Create geometry and solution field
     *    - Write VTK/VTU files for visualization
     *    - Test different cell types (quads, hexas, triangles)
     *    - Test multi-component fields (vector, tensor)
     *    - Verify file format correctness
     * 
     * 2. gsWrite/gsRead for different formats:
     *    - Test IGES format (CAD exchange)
     *    - Test OBJ format (mesh)
     *    - Test OFF format (surface mesh)
     *    - Test STL format (triangulated surface)
     *    - Test CSV format for point clouds
     * 
     * 3. gsMesh writing tests:
     *    - Create gsMesh object
     *    - Write to various mesh formats
     *    - Test boundary marking in output
     *    - Test element data attachment
     * 
     * 4. gsXml advanced tests:
     *    - Test nested XML structures
     *    - Test attribute reading/writing
     *    - Test XML validation
     *    - Test error handling for malformed XML
     * 
     * 5. gsFileManager tests:
     *    - Test file path resolution
     *    - Test data directory management
     *    - Test file search in multiple directories
     *    - Test file extension handling
     * 
     * 6. Binary I/O tests:
     *    - Test binary matrix writing
     *    - Test binary geometry writing
     *    - Compare binary vs XML file sizes
     *    - Test endianness handling
     * 
     * 7. gsOptionList I/O tests:
     *    - Write option list to XML
     *    - Read option list from XML
     *    - Test default value handling
     *    - Test type conversion
     */
}