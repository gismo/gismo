// stdafx.cpp : source file that includes just the standard includes
// GismoSupport.pch will be the pre-compiled header
// stdafx.obj will contain the pre-compiled type information

#include "stdafx.h"


//CLASS_TEMPLATE_INST gsTHBSpline<2>;
//CLASS_TEMPLATE_INST gsTensorBSpline<2>;
//CLASS_TEMPLATE_INST gsTensorNurbs<2>;

//gsTensorBSpline2 a;
//gsTensorNurbs2 b;
//gsTHBSpline2 c;

gsTHBSpline2::gsTHBSpline2()
  : gsTHBSpline<2>()
{}
gsTHBSpline2::gsTHBSpline2(const gsTensorBSplineBasis<2>& b, const gsMatrix<>& c)
  
    : gsTHBSpline<2>(b,c)
{}
gsTHBSpline2::gsTHBSpline2(gsTHBSpline<2>& temp)
  : gsTHBSpline<2>(temp)
{}

gsTensorNurbs2::gsTensorNurbs2(gsTensorNurbs<2>& temp)
  : gsTensorNurbs<2>(temp){}