/** @file HTBasisDegreeElevate.cpp

    @brief Test that

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/


#include "gismo_unittest.h"


template <int dim>
gsTensorBSplineBasis<dim> *makeTBSBasis(gsKnotVector<real_t> knots)
{
    return new gsTensorBSplineBasis<dim>(knots);
}
template <>
gsTensorBSplineBasis<2> *makeTBSBasis<2>(gsKnotVector<real_t> knots)
{
    return new gsTensorBSplineBasis<2>(knots,knots);
}
template <>
gsTensorBSplineBasis<3> *makeTBSBasis<3>(gsKnotVector<real_t> knots)
{
    return new gsTensorBSplineBasis<3>(knots,knots,knots);
}

template <typename Btype>
void insertBoxes(Btype &basis)
{
    const int dim=basis.dim();
    std::vector<unsigned> boxes;

    boxes.push_back(2);
    for (index_t i=0;i<dim;++i)
        boxes.push_back(2*i);
    for (index_t i=0;i<dim;++i)
        boxes.push_back(2*i+4);

    boxes.push_back(3);
    for (index_t i=0;i<dim;++i)
        boxes.push_back(4*i);
    for (index_t i=0;i<dim;++i)
        boxes.push_back(4*i+2);

    basis.refineElements(boxes);
}

template <typename Btype>
void setMultiplicity(Btype &basis,index_t mult)
{
    bool allChecked=true;
    typedef typename Btype::tensorBasis tensorT;
    std::vector<tensorT*> tensorB=basis.getBases();
    do {
        allChecked=true;
        for (size_t l=0; l<tensorB.size();++l)
        {
            for (int d=0;d<Btype::Dim;++d)
            {
                index_t m_min=INT_MAX;
                gsKnotVector<real_t> &knots = tensorB[l]->knots(d);
                std::vector<real_t> toRise;
                typedef typename gsKnotVector<real_t>::iterator iter;
                for (iter i=knots.begin();i<knots.end();++i)
                {
                    if (knots.multiplicity(*i)<mult)
                    {
                        toRise.push_back(*i);
                        m_min=math::min<int>(m_min,mult-knots.multiplicity(*i));
                        allChecked=false;
                    }
                }
                if (m_min==0) break;
                basis.increaseMultiplicity(l,d,toRise,m_min);
            }
        }
    } while (!allChecked);
}


template <typename Btype>
void testImplementation(index_t deg,index_t mult)
{
    const index_t dim=Btype::Dim;
    const real_t eps=2*deg*deg*std::numeric_limits<real_t>::epsilon();

    gsKnotVector<> knots(0, 1, 3, deg+1, mult);
    gsVector<index_t> numPoints;
    numPoints.setConstant(dim,(deg+mult+1)*3);

    // prepare basis
    gsTensorBSplineBasis<Btype::Dim>*temp=makeTBSBasis<Btype::Dim>(knots);
    Btype                             basis  (*temp);
    temp->degreeElevate();
    Btype                             basisE (*temp);
    delete temp;

    insertBoxes(basis);
    insertBoxes(basisE);
    setMultiplicity(basis, mult);
    basis.degreeElevate();
    setMultiplicity(basisE, mult+1);

    // computes values
    gsGridIterator<real_t,CUBE> ptIter(basis.support(), numPoints);
    gsMatrix<real_t>  pt(dim, numPoints.prod());

    for (index_t c = 0; ptIter; ++ptIter, ++c )
        pt.col(c) = *ptIter;

    gsMatrix<real_t> baseEval=basis.eval(pt);
    gsMatrix<real_t> baseEvalE=basisE.eval(pt);
    gsMatrix<unsigned> baseAct  =basis.active(pt);
    gsMatrix<unsigned> baseActE =basisE.active(pt);
    CHECK( gsAllCloseAbsolute(baseEval,baseEvalE, eps ));
    CHECK( (baseAct.array()==baseActE.array()).all() );
}



TEST(HBasisDegreeElevate)
{
// Disabled: expensive test (due to the function setMultiplicity which
// is highly inefficient) which does not help: it computes the same
// thing with 2 identical ways.

// TO DO: a good test would be to test the partition of unity for THB,
// or test (for HB and THB) that interpolation of a polynomial
// function of the same degree is accurate upto machine precision

/* 
    for (index_t deg=0; deg<4;++deg)
    {
        for (index_t mult=1; mult<=deg; ++mult)
        {
            try
            {
                testImplementation<gsHBSplineBasis<2> >(deg,mult);
                testImplementation<gsHBSplineBasis<3> >(deg,mult);

                testImplementation<gsTHBSplineBasis<2> >(deg,mult);
                testImplementation<gsTHBSplineBasis<3> >(deg,mult);
            }
            catch(...)
            {
                std::cout<<"error at deg"<<deg<<" mult "<<mult<<std::endl;
                throw;
            }
        }
    }
//*/
}
