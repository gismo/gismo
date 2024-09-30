#include "stdafx.h"
#include "ThbSurfaceUtils.h"
/// @cond
int ON_GismoUtils::NurbForm(const gsBSpline<>& b, ON_NurbsCurve& crv, double tolerance)
{
    crv = ON_NurbsCurve(3, true, b.degree() + 1, b.coefs().rows());

    const gsKnotVector<>& kts = b.knots();
    for (int i = 1; i < kts.size() - 1; ++i)
    {
        crv.SetKnot(i - 1, kts[i]);
    }

    bool inputRational = b.coefs().cols() == 4;
    for (int i = 0; i < b.coefs().rows(); ++i)
    {
        if (!inputRational)
        {
            ON_4dPoint pt(b.coef(i, 0), b.coef(i, 1), b.coef(i, 2), 1);
            crv.SetCV(i, pt);
        }
        else
        {
            ON_4dPoint pt(b.coef(i, 0), b.coef(i, 1), b.coef(i, 2), b.coef(i, 3));
            crv.SetCV(i, pt);
        }
    }
    return 1;
}

int ON_GismoUtils::NurbForm(const gsTensorBSpline2& nurbs, ON_NurbsSurface& ns, double tolerance)
{
    int rc = ON_GismoUtils::NurbFormImpl(nurbs, ns, tolerance); 
    return rc;
}
/// @endcond

int ON_GismoUtils::NurbForm(const gsTensorNurbs2& nurbs, ON_NurbsSurface& ns, double tolerance)
{
    int rc = ON_GismoUtils::NurbFormImpl(nurbs, ns, tolerance);
    // set weights
    int i = 0;

    for (int col = 0; col < ns.CVCount(1); ++col)
    {
        for (int row = 0; row < ns.CVCount(0); ++row)
        {
            ON_4dPoint cv;
            ns.GetCV(row, col, cv);
            RhinoApp().ActiveDoc()->AddPointObject(cv);

            double w = nurbs.basis().weight(i);
            ++i;
            cv.x = cv.x * w;
            cv.y = cv.y * w;
            cv.z = cv.z * w;
            cv.w = w;
            ns.SetCV(row, col, cv);
        }
    }

    return rc;
}

int ON_GismoUtils::NurbForm(const gsTHBSpline2& thb, ON_NurbsSurface& ns, double tolerance)
{
    gsTensorBSpline2 nurbs;

    gsTHBSpline2 copy(thb);
    copy.convertToBSpline(nurbs);
    
    // note that the orientation is transposed in the 
    // thbspline definition, so we compensate for this here.
    return NurbForm(nurbs, ns);
}

bool ON_GismoUtils::FromSurface(const ON_Surface& srf, gsTHBSpline2& thb)
{
    ON_NurbsSurface ns;
    if (0 == srf.GetNurbForm(ns))
        return false;

    std::vector<double> uKnots;
    uKnots.reserve(2 + ns.KnotCount(0));
    uKnots.push_back(ns.Knot(0, 0));
    for (int i = 0; i < ns.KnotCount(0); ++i)
    {
        uKnots.push_back(ns.Knot(0, i));
    }
    uKnots.push_back(ns.Knot(0, ns.KnotCount(0) - 1));

    std::vector<double> vKnots;
    vKnots.reserve(2 + ns.KnotCount(1));
    vKnots.push_back(ns.Knot(1, 0));
    for (int i = 0; i < ns.KnotCount(1); ++i)
    {
        vKnots.push_back(ns.Knot(1, i));
    }
    vKnots.push_back(ns.Knot(1, ns.KnotCount(1) - 1));

    gsKnotVector<> kU(uKnots, ns.Degree(0));
    gsKnotVector<> kV(vKnots, ns.Degree(1));

    //always use a rational approach with homogeneous coordinates
    gsMatrix<> coefs(ns.CVCount(), 4);
    for (int i = 0; i < ns.CVCount(0); ++i)
    {
        for (int j = 0; j < ns.CVCount(1); ++j)
        {
            ON_4dPoint rcv;
            ns.GetCV(i, j, rcv);

            coefs(j*ns.CVCount(0) + i, 0) = rcv.x;
            coefs(j*ns.CVCount(0) + i, 1) = rcv.y;
            coefs(j*ns.CVCount(0) + i, 2) = rcv.z;
            coefs(j*ns.CVCount(0) + i, 3) = rcv.w;
            
        }
    }
    gsTensorBSplineBasis<2> basis(kU, kV);        

  thb = gsTHBSpline2(basis, coefs);
    
    return true;
}

bool ON_GismoUtils::FromSurface(const ON_Surface& srf, const std::vector<unsigned int>& boxes, gsTHBSpline2& thb)
{
    if (!FromSurface(srf, thb))
        return false;

    thb.refineElements(boxes);
    return true;
}


int ON_GismoUtils::BrepForm(const gsTHBSpline2& thb, ON_Brep& b, double tolerance)
{
    b = ON_Brep();// re-initialize to empty brep. 

    gsMatrix<unsigned> b1, b2;
    gsVector<unsigned> level;
    thb.basis().tree().getBoxes(b1, b2, level); // splitting based on the quadtree

    const int nboxes = level.size();
    gsVector<unsigned> p1, p2;
    gsMatrix<> temp1;
    gsKnotVector<> cku, ckv;
    
    for (int i = 0; i < nboxes; i++) // for all boxes
    {
        p1 = b1.row(i).transpose();
        p2 = b2.row(i).transpose();

        thb.basis().getBsplinePatchGlobal(p1, p2, level[i], thb.coefs(), temp1, cku, ckv);
        gsTensorBSpline2 tbspline(cku, ckv, give(temp1));

        ON_NurbsSurface ns;
        if (0 != ON_GismoUtils::NurbForm(tbspline, ns))
        {
            b.NewFace(ns);
        }
    }

  return 1;
}