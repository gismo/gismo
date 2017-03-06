#include "stdafx.h"
#include "TestCommandImplementation.h"
#include "ThbSurfaceObject.h"

CTestCommandImplementation::CTestCommandImplementation()
{
}


CTestCommandImplementation::~CTestCommandImplementation()
{
}

CRhinoCommand::result CTestCommandImplementation::RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand)
{
    gsKnotVector<> kU(0, 1, 3, 3);
    gsKnotVector<> kV(0.0, 1.0, 3, 3);

    gsTensorBSplineBasis<2> basis(kU, kV);

    // the same refinement boxes as in the example
    std::vector<unsigned int> box;
    box.push_back(1);
    box.push_back(2);
    box.push_back(2);
    box.push_back(6);
    box.push_back(6);
    //box.push_back(2);
    //box.push_back(4);
    //box.push_back(4);
    //box.push_back(10);
    //box.push_back(10);


    gsMatrix<> coefs(basis.size(), 4);

    // fill the coefficients in the unit square.
    for (int i = 0; i < basis.size(0); ++i)
    {
        double x = (double)i / (basis.size(0) - 1);
        for (int j = 0; j < basis.size(1); ++j)
        {
            double y = (double)j / (basis.size(1) - 1);
            // note the orientation of the basis functions is required to be 
            // like this: i-runs-faster-than-j. If not, the refinement does not work.
            coefs(j*basis.size(0) + i, 0) = x;
            coefs(j*basis.size(0) + i, 1) = y;
            coefs(j*basis.size(0) + i, 2) = .2* math::sin(x * 2 * ON_PI) * math::sin(y * 2 * ON_PI);
            coefs(j*basis.size(0) + i, 3) = 1.0; // weight
        }
    }

    gsTHBSpline2* th = new gsTHBSpline2(basis, coefs);
    th->refineElements(box);

    CThbSurfaceObject* newObj = new CThbSurfaceObject();
    newObj->SetHierarchicalSurface(th);

    if (!context.m_doc.AddObject(newObj))
    {
        RhinoApp().Print(L"Failed to add object \n");
    }

    context.m_doc.Redraw();
    return CRhinoCommand::success;
}