#include <iostream>

#include <gismo.h>

#include <fstream>


using namespace gismo;

using std::cout;
using std::endl;
using std::size_t;

int main(int argc, char *argv[])
{ 
    bool plot = 1; // If user gives --plot as argiment, paraview file is generated and launched on exit
    bool toxml =1;
    bool noSmooth = false;
    bool writePatchNumbers = false;
    std::string filename;
    filename= GISMO_SOURCE_DIR;
    filename += "/filedata/off/mushroom_triangulated.off";
    double cutoffAngle = 40.0;
    double innerAngle = 15.0;
    double patchAreaWeight = 0.2;
    double mergeSmallPatches = 2;
    int degree = 4;
    int interiorPts = 5;
    int wEdge = 10;
    int wInterior = 1;

    try 
    {
        
        gsCmdLine cmd("Recover the features of a triangulated surface.");
        
        gsArgValPlain<std::string> a1("filename","File containing the input mesh.",
                                      false,filename,"string",cmd );
        gsArgSwitch ap("", "plot", "Output mesh in ParaView format", cmd);
        gsArgSwitch ax("", "xml", "Output solid to xml file", cmd);
        gsArgSwitch ans("", "nosmooth", "Do not smooth corners of faces", cmd);
        gsArgSwitch wpn("", "writenums", "Write patch numbers to patchnumbers.txt", cmd);
        gsArgVal<int> aWInterior("f", "interiorweight","Surface fitting: weight for interior points.",
                                      false,wInterior,"double precision float",cmd );
        gsArgVal<int> aWEdge("e", "edgeweight","Surface fitting: weight for edges.",
                                      false,wEdge,"double precision float",cmd );
        gsArgVal<int> aInteriorPts("i", "interior","Interior knot points of each dimension of trimmed surface.",
                                      false,interiorPts,"double precision float",cmd );
        gsArgVal<int> aDegree("d", "degree","Degree of each dimension of trimmed surface.",
                                      false,degree,"double precision float",cmd );
        gsArgVal<> aPAW("p", "paw","Patch area weight.",
                                      false,patchAreaWeight,"double precision float",cmd );
        gsArgVal<> aMSP("m", "msp","Merge small patches.",
                                      false,mergeSmallPatches,"double precision float",cmd );
        gsArgVal<> aInnerAngle("y", "innerangle","Cutoff angle (degrees) for second pass.",
                                      false,innerAngle,"double precision float",cmd );
        gsArgVal<> aCutoff("c", "cutoff","Cutoff angle (degrees).",
                                      false,cutoffAngle,"double precision float",cmd );
        
        
        
        cmd.parse(argc,argv);
        filename = a1.getValue();
        if (filename.empty() )
        {
            std::cout<< "Waiting for file input.\n";
            return 0;
        }
        plot = ap.getValue();
        toxml = ax.getValue();
        noSmooth = ans.getValue();
        writePatchNumbers = wpn.getValue();
        cutoffAngle = aCutoff.getValue();
        innerAngle = aInnerAngle.getValue();
        patchAreaWeight = aPAW.getValue();
        mergeSmallPatches = aMSP.getValue();
        degree = aDegree.getValue();
        interiorPts = aInteriorPts.getValue();
        wEdge = aWEdge.getValue();
        wInterior = aWInterior.getValue();

    }
    catch ( gsArgException& e )
    {
        cout << "Error: " << e.error() << " " << e.argId() << endl;
        return 1;
    }

    // decide on the base filename
    size_t nameStartIdx = filename.rfind('/');
    if ( nameStartIdx == std::string::npos) 
        nameStartIdx = 0;
    else
        nameStartIdx+= 1;
    size_t nameEndIdx = filename.rfind('.');
    std::string baseName = filename.substr(nameStartIdx, nameEndIdx - nameStartIdx);
    cout << "--- Settings ---\n";
    cout << "Processing " << baseName << "\n";
    cout << "Cutoff angle is " << cutoffAngle << "\n";
    cout << "Master surfaces have degree " << degree << " and " << interiorPts << " interior knot points.\n";
    cout << "Surface fit: edge weighting is " << wEdge << " and interior point weighting is " << wInterior << ".\n";
    if(noSmooth) cout << "Will NOT smooth corners of faces\n";
    else cout << "WILL smooth corners of faces\n";
    if(plot) cout << "WILL create paraview plot\n";
    else cout << "Will NOT create paraview plot\n";
    cout << "----------------\n";

    gsMesh<> * m = gsReadFile<>(filename);
    if (m)
      cout<< "Got "<< *m <<endl;
    else
    {
      cout<< "Problem encountered in file "<<filename<<", quitting." <<endl;
      return 0;
    }
    std::vector< gsMesh<> * > paraMeshes;// filled inside toSolid
    std::vector< gsMesh<> * > fitMeshes; // filled inside toSolid
    std::vector< gsMesh<> * > patchMeshes; // filled inside toSolid

    std::vector< std::vector< gsVertex<>* > > iPoints;
    std::vector< std::vector< gsVertex<>* > > oPoints;
    std::vector< std::vector< std::vector<gsVertex<>* > > >  innerBdrys;
    std::vector< std::vector< gsVertex<> > > innerBdrysMassP;
    std::vector< std::vector< bool > > oPointsConvexFlag;
    
    // set up the gsTriMeshToSolid class to perform the feature detection
    gsTriMeshToSolid<> tmts(m);
    bool non_manifold, warning_borders;
    
    // compute the features
    tmts.getFeatures(cutoffAngle, non_manifold, warning_borders);
    m->cleanStlMesh();
    
    // give every face a patch number
    tmts.calcPatchNumbers();
    
    // improve quality by further dividing and merging patches according to more complex rules
    tmts.storeNeighboringFaces();
    tmts.divideAndMergePatches(innerAngle, patchAreaWeight, mergeSmallPatches);
    
    // recompute patch numbers
    tmts.calcPatchNumbers();
    
    // write patch numbers
    cout << "Writing patch numbers...\n";
    if(writePatchNumbers)
    {
        std::ofstream pn("patchnumbers.txt");
        
        for(std::vector<gsMeshElement<>::gsFaceHandle >::iterator it(m->face.begin());it!=m->face.end();++it)
        {
            pn << (**it).faceIdentity << "\n";
        }
    }
    
    // get the patches
    tmts.getFaces(iPoints, oPoints, innerBdrys, innerBdrysMassP, oPointsConvexFlag);
    
    gsSolid<> * sl = new gsSolid<>();
    //if you need a higher number of interior points when fitting the surface increase the 8th input variiable(now 5)
    tmts.toSolid(*sl,iPoints,oPoints,innerBdrys,innerBdrysMassP,oPointsConvexFlag,paraMeshes,fitMeshes,patchMeshes,degree,interiorPts,1,300,1,wEdge,wInterior,1, noSmooth);
    gsInfo<<*sl<<'\n';
    
    if (toxml)
    {
        cout << "Writing xml file..." << endl;
        
        gsFileData<> newdata;
        newdata << *sl;
        newdata.dump(baseName);
    }

    //write sharpness information to a file:
    // first line: #patches
    // next #patches lines: #outer points of patch i, sharpness of first point ... sharpness of last point of patch i
    // next #patches lines: #holes of patch i
    // next #holes of patch i lines: #points of hole j of patch i, sharpness of first point ... sharpness of last point of hole j of patch i
    std::string outputFilename;
    outputFilename = "corner_data.txt";
    std::ofstream myfile (outputFilename.c_str());
    if (myfile.is_open())
    {
      myfile << oPoints.size()<<"\n";

      for(size_t i=0;i<oPoints.size();i++)
      {
          myfile << oPoints[i].size()<<" ";
          for(size_t j=0;j<oPoints[i].size();j++)
             myfile<<(oPoints[i][j]->numEdges>2)<<" ";
          myfile << "\n";

      }
      myfile << "\n";

      for(size_t i=0;i<innerBdrys.size();i++)
      {
          myfile << innerBdrys[i].size()<<" ";
          for(size_t j=0;j<innerBdrys[i].size();j++)
          {
              myfile << innerBdrys[i][j].size()<<" ";
              for(size_t k=0;k<innerBdrys[i][j].size();k++)
                  myfile<<(innerBdrys[i][j][k]->numEdges>2)<<" ";
          }
          myfile << "\n";
      }


      myfile.close();
    }
    if (plot)
    {
        // Write a paraview file
        cout<<"Writing paraview file..." << endl;

        gsWriteParaview( *m, "output");
    }
    delete m;
    delete sl;

    // free meshes
    freeAll(fitMeshes);
    freeAll(paraMeshes);
    freeAll(patchMeshes);

    return 0;
}
