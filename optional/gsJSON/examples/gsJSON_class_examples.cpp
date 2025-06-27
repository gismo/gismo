/** @file gsJSON_class_examples.cpp

    @brief Examples using gsJSON class methods

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, J. Li
*/

#include <gismo.h>
#include <gsJSON/gsJSON.h>

using namespace gismo;

void exampleBasicUsage()
{
    gsInfo << "\n--- Basic gsJSON usage ---\n";
    
    // Create empty gsJSON object
    gsJSON j;
    
    // Add data using add() method
    j.add("name", "test");
    j.add("value", 42);
    j.add("pi", 3.14159);
    
    // Access using operator[]
    j["array"] = {1, 2, 3, 4};
    
    // Get data using get() method
    std::string name = j.get<std::string>("name");
    int val = j.get<int>("value");
    double pi = j.get<double>("pi");
    
    gsInfo << "name: " << name << "\n";
    gsInfo << "value: " << val << "\n";
    gsInfo << "pi: " << pi << "\n";
    
    // Check if key exists
    if (j.contains("name"))
        gsInfo << "Found 'name' key\n";
    
    // Size and empty check
    gsInfo << "Size: " << j.size() << "\n";
    gsInfo << "Empty: " << (j.empty() ? "yes" : "no") << "\n";
}

void exampleWithGismoTypes()
{
    gsInfo << "\n--- gsJSON with G+Smo types ---\n";
    
    gsJSON j;
    
    // Vector
    gsVector<> vec(3);
    vec << 1.0, 2.0, 3.0;
    j.add("vector", vec);
    
    // Matrix  
    gsMatrix<> mat(2, 2);
    mat << 1, 2, 3, 4;
    j.add("matrix", mat);
    
    // Option list
    gsOptionList opt;
    opt.addInt("iter", "iterations", 100);
    opt.addReal("tol", "tolerance", 1e-6);
    j.add("options", opt);
    
    // Knot vector
    gsKnotVector<> kv(0, 1, 1, 2);
    j.add("knots", kv);
    
    // Retrieve them back
    gsVector<> vecBack = j.get<gsVector<>>("vector");
    gsMatrix<> matBack = j.get<gsMatrix<>>("matrix");
    gsOptionList optBack = j.get<gsOptionList>("options");
    gsKnotVector<> kvBack = j.get<gsKnotVector<>>("knots");
    
    gsInfo << "Vector: " << vecBack.transpose() << "\n";
    gsInfo << "Matrix:\n" << matBack << "\n";
    gsInfo << "Options: " << optBack << "\n";
    gsInfo << "Knots: " << kvBack << "\n";
}

void exampleIterators()
{
    gsInfo << "\n--- gsJSON iterators ---\n";
    
    gsJSON j;
    j.add("a", 1);
    j.add("b", 2.5);
    j.add("c", "three");
    
    // Iterate over all items
    for (auto it = j.begin(); it != j.end(); ++it)
    {
        gsInfo << it.key() << ": " << it.value() << "\n";
    }
    
    // Find specific key
    auto found = j.find("b");
    if (found != j.end())
        gsInfo << "Found key 'b' with value: " << found.value() << "\n";
    
    // Count occurrences
    gsInfo << "Count of 'a': " << j.count("a") << "\n";
    
    // Erase a key
    j.erase("c");
    gsInfo << "After erasing 'c', size: " << j.size() << "\n";
}

void exampleFileOperations()
{
    gsInfo << "\n--- File operations ---\n";
    
    // Create gsJSON and populate it
    gsJSON j;
    
    // Add some geometry data
    gsKnotVector<> kv(0, 1, 2, 3);
    gsBSplineBasis<> basis(kv);
    j.add("basis", basis);
    
    gsMatrix<> coefs(basis.size(), 2);
    coefs.setRandom();
    j.add("coefficients", coefs);
    
    // Save to file
    std::string filename = "test_gsjson.json";
    j.save(filename);
    gsInfo << "Saved to " << filename << "\n";
    
    // Read from file
    gsJSON j2(filename);
    gsInfo << "Loaded from file\n";
    
    // Print contents
    gsInfo << "Contents:\n" << j2 << "\n";
    
    // Get data back
    gsBSplineBasis<> basisRead = j2.get<gsBSplineBasis<>>("basis");
    gsInfo << "Basis degree: " << basisRead.degree() << "\n";
}

void exampleFromOptionList()
{
    gsInfo << "\n--- Create from OptionList ---\n";
    
    // Create option list
    gsOptionList opt;
    opt.addInt("maxIter", "", 1000);
    opt.addReal("tol", "", 1e-8);
    opt.addString("solver", "", "CG");
    opt.addSwitch("verbose", "", true);
    
    // Create gsJSON from option list
    gsJSON j(opt);
    
    gsInfo << "JSON from options:\n" << j << "\n";
    
    // Modify it
    j.add("extraParam", 123);
    
    // Get back as option list
    gsOptionList optBack = j.get<gsOptionList>();
    gsInfo << "Options back:\n" << optBack << "\n";
}

void exampleClearAndModify()
{
    gsInfo << "\n--- Clear and modify ---\n";
    
    gsJSON j;
    j.add("a", 1);
    j.add("b", 2);
    j.add("c", 3);
    
    gsInfo << "Initial size: " << j.size() << "\n";
    
    // Clear all
    j.clear();
    gsInfo << "After clear, empty: " << (j.empty() ? "yes" : "no") << "\n";
    
    // Add new data
    j["new_data"]["nested"] = "value";
    j["new_data"]["number"] = 42;
    
    gsInfo << "New contents:\n" << j << "\n";
}

void exampleGetTo()
{
    gsInfo << "\n--- Using get_to ---\n";
    
    gsJSON j;
    
    // Add a matrix
    gsMatrix<> mat(3, 2);
    mat << 1, 2, 3, 4, 5, 6;
    j.add("data", mat);
    
    // For specific field, use operator[] with get_to
    gsMatrix<> matOut;
    j["data"].get_to(matOut);
    
    gsInfo << "Matrix from get_to:\n" << matOut << "\n";
    
    // Also works with gsJSON's get method
    gsMatrix<> matOut2 = j.get<gsMatrix<>>("data");
    gsInfo << "Same matrix from get:\n" << matOut2 << "\n";
}

void exampleMultiPatchCreation()
{
    gsInfo << "\n--- MultiPatch from JSON ---\n";
    
    // Create JSON with patch creation info
    gsJSON j;
    j["patch"]["create"]["type"] = "BSplineSquare";
    j["patch"]["create"]["length"] = 2.0;
    j["patch"]["create"]["x"] = 0.0;
    j["patch"]["create"]["y"] = 0.0;
    
    // This will create a patch when read
    gsMultiPatch<> mp = j["patch"].get<gsMultiPatch<>>();
    gsInfo << "Created multipatch with " << mp.nPatches() << " patches\n";
    
    // Read from XML file using path functionality
    gsJSON j2;
    // Use relative path to the filedata directory
    std::string xmlPath = "../optional/gsJSON/filedata/thbs_face_3levels.xml";
    j2["geometry"]["path"] = xmlPath;
    
    gsMultiPatch<> mp2 = j2["geometry"].get<gsMultiPatch<>>();
    gsInfo << "Loaded from XML: " << mp2.nPatches() << " patches\n";
    gsInfo << "First patch domain dim: " << mp2.patch(0).domainDim() << "\n";
    gsInfo << "First patch target dim: " << mp2.patch(0).targetDim() << "\n";
    
    // Check if it's a THB-spline
    if (mp2.patch(0).targetDim() == 3) 
    {
        gsInfo << "This is a 3D surface geometry\n";
        // Get bounding box
        gsMatrix<> bbox;
        mp2.boundingBox(bbox);
        gsInfo << "Bounding box:\n";
        gsInfo << "  X: [" << bbox(0,0) << ", " << bbox(0,1) << "]\n";
        gsInfo << "  Y: [" << bbox(1,0) << ", " << bbox(1,1) << "]\n";
        gsInfo << "  Z: [" << bbox(2,0) << ", " << bbox(2,1) << "]\n";
    }
    
    // Visualize the geometry (save as VTK for visualization)
    std::string vtkFile = "thbs_face_visualization.vts";
    gsWriteParaview(mp2, vtkFile, 1000);
    gsInfo << "Visualization saved to: " << vtkFile << "\n";
}

int main(int argc, char *argv[])
{
    gsCmdLine cmd("Examples for gsJSON class usage.");
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
    
    gsInfo << "=== gsJSON Class Examples ===\n";
    
    exampleBasicUsage();
    exampleWithGismoTypes();
    exampleIterators();
    exampleFileOperations();
    exampleFromOptionList();
    exampleClearAndModify();
    exampleGetTo();
    exampleMultiPatchCreation();
    
    gsInfo << "\nAll examples completed.\n";
    
    return EXIT_SUCCESS;
}