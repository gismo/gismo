#include "stdafx.h"
//#include "THBSpline.h"
//
//ON_OBJECT_IMPLEMENT(ON_ThbSurface, ON_NurbsSurface, "218397bf-1bde-4d03-9b9e-f2c3130edd29");
//
//ON_ThbSurface::ON_ThbSurface()
//	: ON_NurbsSurface(3, false, 1, 1, 0, 0)
//	, m_thb(nullptr)	
//	, m_deleteThb(false)
//{
//	
//}
//
//ON_ThbSurface::ON_ThbSurface(const ON_ThbSurface& src)
//	: ON_NurbsSurface(src)
//{
//	if (m_thb) delete m_thb;
//	m_thb = new gsTHBSpline<2>(*src.m_thb);
//	m_deleteThb = true; // I made it
//}
//
//ON_ThbSurface& ON_ThbSurface::operator=(const ON_ThbSurface& src)
//{
//	if (&src != this)
//	{
//		ON_NurbsSurface::operator=(src);
//		if (m_thb)
//		{
//			delete m_thb;
//			m_thb = nullptr;
//		}
//		if (src.m_thb)
//		{
//			m_thb = new gsTHBSpline<2>(*src.m_thb);
//			m_deleteThb = true; // I made it
//
//		}
//	}
//	return *this;
//}
//
//bool ON_ThbSurface::GiveTHBSpline(gsTHBSpline<2>* p)
//{
//	m_thb = p;
//	m_deleteThb = true; // I am now responsible for it
//	return 0 != GetNurbForm(*this);
//}
//
//bool ON_ThbSurface::LoanTHBSpline(gsTHBSpline<2>* p)
//{
//	m_thb = p;
//	m_deleteThb = false; // I am not responsible for it
//	return 0 != GetNurbForm(*this);
//}
//
//bool ON_ThbSurface::FromSurface(const ON_Surface& srf)
//{
//	if (m_thb && m_deleteThb)
//		delete m_thb;
//	m_thb = nullptr;
//
//	ON_NurbsSurface ns;
//	if (0 == srf.GetNurbForm(ns))
//		return false;
//
//	std::vector<double> uKnots;
//	uKnots.reserve(2 + ns.KnotCount(0));
//	uKnots.push_back(ns.Knot(0, 0));
//	for (int i = 0; i < ns.KnotCount(0); ++i)
//	{
//		uKnots.push_back(ns.Knot(0, i));
//	}
//	uKnots.push_back(ns.Knot(0, ns.KnotCount(0) - 1));
//	
//	std::vector<double> vKnots;
//	vKnots.reserve(2 + ns.KnotCount(1));
//	vKnots.push_back(ns.Knot(1, 0));
//	for (int i = 0; i < ns.KnotCount(1); ++i)
//	{
//		vKnots.push_back(ns.Knot(1, i));
//	}
//	vKnots.push_back(ns.Knot(1, ns.KnotCount(1) - 1));
//
//
//	gsKnotVector<> kU(uKnots, ns.Degree(0));
//	gsKnotVector<> kV(vKnots, ns.Degree(1));
//
//	gsMatrix<> coefs(ns.CVCount(), 3);
//
//	for (int i = 0; i < ns.CVCount(0); ++i)
//	{
//		for (int j = 0; j < ns.CVCount(1); ++j)
//		{
//			ON_3dPoint cv;
//			ns.GetCV(i, j, cv);
//			coefs(j*ns.CVCount(0) + i, 0) = cv.x;
//			coefs(j*ns.CVCount(0) + i, 1) = cv.y;
//			coefs(j*ns.CVCount(0) + i, 2) = cv.z;
//		}
//	}
//
//	gsTensorBSplineBasis<2> basis(kU, kV);
//	gsTHBSpline<2>* thb = new gsTHBSpline<2>(basis, coefs);
//	return GiveTHBSpline(thb);	
//}
//
//bool ON_ThbSurface::FromSurface(const ON_Surface& srf, const std::vector<unsigned int>& boxes)
//{
//	if (!FromSurface(srf))
//		return false;
//
//	m_thb->refineElements(boxes);
//	return 0 != GetNurbForm(*this);
//}
//
//ON_ThbSurface::~ON_ThbSurface()
//{
//	if (m_deleteThb && m_thb) delete m_thb;
//	m_thb = nullptr;
//}
//
//ON_Surface* ON_ThbSurface::DuplicateSurface() const
//{
//	gsTHBSpline<2>* copy(nullptr);
//	if (m_thb) 
//		copy = new gsTHBSpline<2>(*m_thb);	
//	
//	ON_ThbSurface* toReturn = new ON_ThbSurface();
//	toReturn->GiveTHBSpline(copy);
//	return toReturn;
//}
//
//ON_BOOL32 ON_ThbSurface::SetDomain(int dir, double t0, double t1)
//{
//	if (!m_thb) return false;
//	// not sure how to do this, may need to re-create the basis
//	return false;
//}
//
//ON_Interval ON_ThbSurface::Domain(int dir) const
//{
//	ON_Interval iv = ON_Interval::EmptyInterval;
//	if (!m_thb) return iv;
//
//	const gsTHBSplineBasis<2>& basis = m_thb->basis();
//	gsMatrix<> supp = basis.supportInterval(dir);
//	iv.m_t[0] = supp(0, 0);
//	iv.m_t[1] = supp(0, 1);
//
//	return iv;
//}
//
//int ON_ThbSurface::SpanCount(int dir) const
//{
//	if (!m_thb) return 0;
//	int nb = m_thb->basis().numBreaks(m_thb->basis().maxLevel(), dir);
//	return nb;
//}
//
//int ON_ThbSurface::Degree(int dir) const
//{
//	if (!m_thb) return 0;
//	return m_thb->basis().degree(dir);
//}
//
//ON_BOOL32 ON_ThbSurface::GetSpanVector(int dir, double* sv) const
//{
//	if (!m_thb) return false;
//	int nKnots = m_thb->basis().numKnots(m_thb->basis().maxLevel(), dir);
//	int svLength = 1 + m_thb->basis().numBreaks(m_thb->basis().maxLevel(), dir);
//	for (int i = 0, n = 0; i < nKnots && n < svLength; ++i)
//	{
//		double k = m_thb->basis().knot(m_thb->basis().maxLevel(), dir, i);
//		if (n > 0 && sv[n - 1] == k) continue;
//		sv[n++] = k;
//	}
//	return true;
//}
//
//ON_BOOL32 ON_ThbSurface::Reverse(int dir)
//{
//	RhinoApp().Print("%s\n", __FUNCTION__);
//
//	if (!m_thb) return false;
//	return false; // for now
//	//gsTHBSplineBasis<2>& basis = m_thb->basis();
//	//basis.basis(dir).reverse();
//
//	//return true;
//}
//
//ON_BOOL32 ON_ThbSurface::Transpose()
//{
//	RhinoApp().Print("%s\n", __FUNCTION__);
//
//	if (!m_thb) return false;
//	return false; // for now
//	//gsTHBSplineBasis<2>& basis = m_thb->basis();
//	//basis.transpose();
//
//	//return true;
//}
//
//ON_BOOL32 ON_ThbSurface::Evaluate( // returns false if unable to evaluate
//	double u, double v,   // evaluation parameters
//	int num_der,          // number of derivatives (>=0)
//	int array_stride,     // array stride (>=Dimension())
//	double* der_array,    // array of length stride*(ndir+1)*(ndir+2)/2
//	int quadrant,         // optional - determines which quadrant to evaluate from
//						  //         0 = default
//						  //         1 from NE quadrant
//						  //         2 from NW quadrant
//						  //         3 from SW quadrant
//						  //         4 from SE quadrant
//	int* hint             // optional - evaluation hint (int[2]) used to speed
//						  //            repeated evaluations
//) const
//{
//	if (!m_thb) return false;
//
//	gsMatrix<> uv(2, 1); 
//	uv(0, 0) = u;
//	uv(1, 0) = v;
//
//	std::vector<gsMatrix<> > results;
//	m_thb->evalAllDers_into(uv, num_der, results);
//
//	int n = 0;
//	for (int i = 0; i <= num_der; ++i)
//	{
//		const gsMatrix<>& result = results[i];
//		auto nCols = result.cols();
//		for (int j = 0; j < nCols; ++j)
//		{
//#ifdef _DEBUG
//			if (result.rows() != array_stride)
//			{
//				return false;
//			}
//#endif
//			for (int k = 0; k < array_stride; ++k)
//			{
//				der_array[n++] = result(k, j);
//			}
//		}
//	}
//
//	return true;
//}
//
//int ON_ThbSurface::GetNurbForm(ON_NurbsSurface& ns, double tolerance) const
//{
//	if (!m_thb) return 0;
//
//	gsTensorBSpline<2> nurbs;
//
//	gsTHBSpline<2> copy(*m_thb);
//	copy.convertToBSpline(nurbs);
//
//	// note that the orientation is transposed in the 
//	// thbspline definition, so we compensate for this here.
//
//	// I'm not sure if this is really necessary, it could be achieved
//	// by rotating the U- and V- of the ON_NUrbsSurface too?!
//	int cvc0 = nurbs.basis().component(0).size();
//	int cvc1 = nurbs.basis().component(1).size();
//
//	// nice: I did not know this was possibble, but the surface can be initialized like this.
//	ns = ON_NurbsSurface(3, false, nurbs.degree(0) + 1, nurbs.degree(1) + 1, cvc0, cvc1);
//
//	for (int dir = 0; dir < 2; ++dir)
//	{
//		gsKnotVector<> kts = nurbs.basis().component(dir).knots();//nurbs.knots(dir);
//																  // ignore superfluous knots
//		for (int i = 1; i < kts.size() - 1; ++i)
//		{
//			ns.SetKnot(dir, i - 1, kts[i]);
//		}
//	}
//
//	int c = 0;
//	for (int j = 0; j < cvc1; ++j)
//	{
//		for (int i = 0; i < cvc0; ++i)
//		{
//			// again, i-runs-faster-than-j
//			ON_3dPoint pt(nurbs.coef(c, 0), nurbs.coef(c, 1), nurbs.coef(c, 2));
//			ns.SetCV(i, j, pt);
//			++c;
//		}
//	}
//
//	return 1;
//
//}
//
//int ON_ThbSurface::HasNurbForm() const
//{
//	if (!m_thb) return 0;
//	return 1;
//}
//
//ON_BOOL32 ON_ThbSurface::Transform(const ON_Xform& x)
//{
//	if (!m_thb) return false;
//	const gsMatrix<>& coefs = m_thb->coefs();
//	gsMatrix<> transformed(coefs);
//
//	ON_3dPoint cv;
//	for (int c = 0; c < coefs.rows(); ++c)
//	{
//		cv.x = coefs(c, 0);
//		cv.y = coefs(c, 1);
//		cv.z = coefs(c, 2);
//		cv.Transform(x);
//		transformed(c, 0) = cv.x;
//		transformed(c, 1) = cv.y;
//		transformed(c, 2) = cv.z;
//
//	}
//
//	m_thb->setCoefs(transformed);
//	return true;
//}
//
//bool ON_ThbSurface::Morph(const ON_SpaceMorph& morph)
//{
//	if (!m_thb) return false;
//	const gsMatrix<>& coefs = m_thb->coefs();
//	gsMatrix<> morphed(coefs);
//
//	ON_3dPoint cv;
//	for (int c = 0; c < coefs.rows(); ++c)
//	{
//		cv.x = coefs(c, 0);
//		cv.y = coefs(c, 1);
//		cv.z = coefs(c, 2);
//		cv = morph.MorphPoint(cv);
//		morphed(c, 0) = cv.x;
//		morphed(c, 1) = cv.y;
//		morphed(c, 2) = cv.z;
//
//	}
//
//	m_thb->setCoefs(morphed);
//	return true;
//}
//
//
//bool ON_ThbSurface::IsMorphable() const
//{
//	return true;
//}
//bool ON_ThbSurface::IsDeformable() const
//{
//	return true;
//}