/** @file ThbSurfaceUtils.h 
    @brief some general utility functions to convert between G+SMO and OpenNurbs
*/
#pragma once

class ON_GismoUtils
{
public:
    static int NurbForm(const gsTHBSpline2& thb, ON_NurbsSurface& ns, double tolerance = 0);
    static int BrepForm(const gsTHBSpline2& thb, ON_Brep& b, double tolerance = 0);
    static int NurbForm(const gsTensorNurbs2& nurbs, ON_NurbsSurface& ns, double tolerance = 0);
    static int NurbForm(const gsTensorBSpline2& nurbs, ON_NurbsSurface& ns, double tolerance = 0);
    static int NurbForm(const gsBSpline<>& b, ON_NurbsCurve& crv, double tolerance = 0);

    static bool FromSurface(const ON_Surface& srf, gsTHBSpline2& thb);
    static bool FromSurface(const ON_Surface& srf, const std::vector<unsigned int>& boxes, gsTHBSpline2& thb);

private:
    template<class T>
    static int NurbFormImpl(const T& nurbs, ON_NurbsSurface& ns, double tolerance = 0)
    {
        int cvc0 = nurbs.basis().component(0).size();
        int cvc1 = nurbs.basis().component(1).size();
        
        // always use a rational surface, even if all weights are 1.0
        ns = ON_NurbsSurface(3, true, nurbs.degree(0) + 1, nurbs.degree(1) + 1, cvc0, cvc1);
        
        for (int dir = 0; dir < 2; ++dir)
        {
            gsKnotVector<> kts = nurbs.knots(dir);
            // ignore superfluous knots
            for (int i = 1; i < kts.size() - 1; ++i)
            {
                ns.SetKnot(dir, i - 1, kts[i]);
            }
        }
        
        int c = 0;
        bool inputRational = nurbs.coefs().cols() == 4;
        for (int j = 0; j < cvc1; ++j)
        {
            for (int i = 0; i < cvc0; ++i)
            {
                if (!inputRational)
                {
                    // again, i-runs-faster-than-j
                    ON_4dPoint pt(nurbs.coef(c, 0), nurbs.coef(c, 1), nurbs.coef(c, 2), 1.0);
                    ns.SetCV(i, j, pt);
                    ++c;
                }
                else
                {
                    ON_4dPoint pt(nurbs.coef(c, 0), nurbs.coef(c, 1), nurbs.coef(c, 2), nurbs.coef(c, 3));
                    ns.SetCV(i, j, pt);
                    ++c;
                }
        
            }
        }
        ON_wString err;
    ON_TextLog log(err);
    if (!ns.IsValid(&log))
    {
      RhinoApp().Print(err);
      return 0;
    }
        return 1;
        
    }
};
