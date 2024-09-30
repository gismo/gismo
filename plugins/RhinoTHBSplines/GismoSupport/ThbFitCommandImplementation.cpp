#include "stdafx.h"
#include "ThbFitCommandImplementation.h"
#include "ThbSurfaceObject.h"
#include "gsStreamData.h"

CThbFitCommandImplementation::CThbFitCommandImplementation()
{
}


CThbFitCommandImplementation::~CThbFitCommandImplementation()
{
}

CRhinoCommand::result CThbFitCommandImplementation::RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand)
{
    CRhinoGetObject go;
    go.SetGeometryFilter(CRhinoGetObject::mesh_object);
    go.SetCommandPrompt(L"Select mesh with UV-mapping (texture coordinates)");

    CRhinoGet::result gr = go.GetObjects();
    if (gr != CRhinoGet::object)
        return CRhinoCommand::success;

    CRhinoObjRef oRef = go.Object(0);
    const ON_Mesh* m = oRef.Mesh();
    if(!m) 
        return CRhinoCommand::success;

    if (m->m_V.Count() != m->m_T.Count())
    {
        RhinoApp().Print("No texture coordinates defined\n");
        return CRhinoCommand::success;
    }

    const ON_MeshTopology& top = m->Topology();
    gsMatrix<> xyz(3, top.m_topv.Count());
    gsMatrix<>  uv(2, top.m_topv.Count());

    for (int topvi = 0; topvi < top.m_topv.Count(); ++topvi)
    {
        if (top.m_topv[topvi].m_v_count <= 0) continue;
        
        ON_3dPoint pt = top.TopVertexPoint(topvi);
        xyz(0, topvi) = pt.x;
        xyz(1, topvi) = pt.y;
        xyz(2, topvi) = pt.z;

        int vi = top.m_topv[topvi].m_vi[0];
        ON_2dPoint tc = m->m_T[vi];

        uv(0, topvi) = tc.x;
        uv(1, topvi) = tc.y;
    }

    real_t u_min = uv.row(0).minCoeff(),
        u_max = uv.row(0).maxCoeff(),
        v_min = uv.row(1).minCoeff(),
        v_max = uv.row(1).maxCoeff();

    const int iter = 2;
    const int deg_x = 3;
    const int deg_y = 3;
    const int numURef = 3;
    const int extension = 2;
    real_t refPercent = 0.1;
    real_t lambda = 1e-07;
    real_t threshold = 1e-02;
    real_t tolerance = 1e-02;


    // Create knot-vectors without interior knots
    gsKnotVector<> u_knots(u_min, u_max, 0, deg_x + 1);
    gsKnotVector<> v_knots(v_min, v_max, 0, deg_y + 1);

    // Create a tensor-basis and apply initial uniform refinement    
    gsTensorBSplineBasis<2> T_tbasis(u_knots, v_knots);
    T_tbasis.uniformRefine((1 << numURef) - 1);

    // Create Initial hierarchical basis
    gsTHBSplineBasis<2>  THB(T_tbasis);

    // Specify extension size in u and v cells
    std::vector<unsigned> ext;
    ext.push_back(extension);
    ext.push_back(extension);

    // Create hierarchical refinement object
    gsHFitting<2, real_t> ref(uv, xyz, THB, refPercent, ext, lambda);

    for (int i = 0; i <= iter; i++)
    {
        ref.nextIteration(tolerance, threshold);
        if (ref.maxPointError() < tolerance)
        {
            RhinoApp().Print(L"Error tolerance achieved after %i iterations.\n", i);
            break;
        }
    }
    
    gsStreamData sd;
    gsGeometry<>* res = ref.result();
    gsTHBSpline2* resThb = static_cast<gsTHBSpline2*>(res);
    
    // vital: make a copy of the result, because when the fitting
    // object 'ref' goes out of scope, the result object will be deleted, leading to an
    // already deleted thb spline definition. crash guaranteed.
    gsTHBSpline2* fitted = new gsTHBSpline2(*resThb);
    CThbSurfaceObject* obj = new CThbSurfaceObject();
    obj->SetHierarchicalSurface(fitted);
    context.m_doc.AddObject(obj);
    context.m_doc.Redraw();

    return CRhinoCommand::success;
}
