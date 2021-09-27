/** @file gsPatchRule.hpp

    @brief Provides implementation of the Gauss-Legendre quadrature rule

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

template<class T>
gsPatchRule<T>::gsPatchRule(const gsBasis<T> & basis,
                            const index_t degree,
                            const index_t regularity,
                            const bool overintegrate,
                            const short_t fixDir
                            )
                            :
                            m_basis(&basis),
                            m_deg(degree),
                            m_reg(regularity),
                            m_over(overintegrate),
                            m_fixDir(fixDir)
{
    GISMO_ENSURE(m_reg<m_deg,"regularity cannot be greater or equal to the order!");
    // Initialize some stuff
    m_dim = m_basis->dim();

    GISMO_ASSERT( m_fixDir < short_t(m_dim) && m_fixDir>-2, "Invalid input fixDir = "<<m_fixDir);

    m_nodes.resize(m_dim);
    m_weights.resize(m_dim);
    m_maps.resize(m_dim);

    gsKnotVector<T> knots;
    gsMatrix<T> greville;
    gsVector<T> integral;
    gsBSplineBasis<T> * Bbasis;

    // Loop over dimensions of the basis and store the nodes and weights for each dimension
    for (size_t d = 0; d != m_dim; d++)
    {
        m_end = m_basis->support().col(1);
        if (short_t(d)==m_fixDir && m_fixDir!=-1)
        {
            m_nodes[d].resize(2);
            m_nodes[d]<<0,1;
            m_weights[d].resize(2);
            m_weights[d]<<1,1;
        }
        else
        {
            // Construct temporary basis (must be B-spline because we use knots!)
            Bbasis = const_cast<gsBSplineBasis<T> *>(static_cast<const gsBSplineBasis<T> * >(&m_basis->component(d)));

            // Find the knots
            knots = this->_init(Bbasis);
            // Compute exact integrals
            #if __cplusplus >= 201103L || _MSC_VER >= 1600
                std::tie(greville,integral) = this->_integrate(knots);
            #else
                std::pair< gsMatrix<T>,gsVector<T> > tmp1 = this->_integrate(knots);
                tmp1.first.swap(greville);
                tmp1.second.swap(integral);
            #endif

            // Compute quadrule
            #if __cplusplus >= 201103L || _MSC_VER >= 1600
                std::tie(m_nodes[d],m_weights[d]) = this->_compute(knots,greville,integral);
            #else
                std::pair< gsVector<T>,gsVector<T> > tmp2 = this->_compute(knots,greville,integral);
                tmp2.first.swap(m_nodes[d]);
                tmp2.second.swap(m_weights[d]);
            #endif
        }

        // Construct a map with the nodes and the weights
        for (index_t k=0; k!=m_nodes[d].size(); k++)
            m_maps[d][m_nodes[d].at(k)] = m_weights[d].at(k);

    }
}



/*
    Maps as follows:

    COORDINATES
    d=1:[x1 x2 x3 x4]

    d=2:[x1 x2 x3 x4][x1 x2 x3 x4][x1 x2 x3 x4][x1 x2 x3 x4]
        [y1 y1 y1 y1][y2 y2 y2 y2][y3 y3 y3 y3][y4 y4 y4 y4]

    d=3:[x1 x2 x3 x4][x1 x2 x3 x4][x1 x2 x3 x4][x1 x2 x3 x4]...
        [y1 y1 y1 y1][y2 y2 y2 y2][y3 y3 y3 y3][y4 y4 y4 y4]...
        [z1 z1 z1 z1][z1 z1 z1 z1][z1 z1 z1 z1][z1 z1 z1 z1]...

        [x1 x2 x3 x4][x1 x2 x3 x4][x1 x2 x3 x4][x1 x2 x3 x4]...
        [y1 y1 y1 y1][y2 y2 y2 y2][y3 y3 y3 y3][y4 y4 y4 y4]...
        [z2 z2 z2 z2][z2 z2 z2 z2][z2 z2 z2 z2][z2 z2 z2 z2]...

        [x1 x2 x3 x4][x1 x2 x3 x4][x1 x2 x3 x4][x1 x2 x3 x4]...
        [y1 y1 y1 y1][y2 y2 y2 y2][y3 y3 y3 y3][y4 y4 y4 y4]...
        [z3 z3 z3 z3][z3 z3 z3 z3][z3 z3 z3 z3][z3 z3 z3 z3]...

        [x1 x2 x3 x4][x1 x2 x3 x4][x1 x2 x3 x4][x1 x2 x3 x4]
        [y1 y1 y1 y1][y2 y2 y2 y2][y3 y3 y3 y3][y4 y4 y4 y4]
        [z4 z4 z4 z4][z4 z4 z4 z4][z4 z4 z4 z4][z4 z4 z4 z4]

    WEIGHTS
    d=1:[wx1 wx2 wx3 wx4]

    d=2: [wy1*[wx1 wx2 wx3 wx4]] [wy2*[wx1 wx2 wx3 wx4]] [wy3*[wx1 wx2 wx3 wx4]] [wy4*[wx1 wx2 wx3 wx4]]

    d=3: wz1*[[wy1*[wx1 wx2 wx3 wx4]] [wy2*[wx1 wx2 wx3 wx4]] [wy3*[wx1 wx2 wx3 wx4]] [wy4*[wx1 wx2 wx3 wx4]]]...
         wz2*[[wy1*[wx1 wx2 wx3 wx4]] [wy2*[wx1 wx2 wx3 wx4]] [wy3*[wx1 wx2 wx3 wx4]] [wy4*[wx1 wx2 wx3 wx4]]]...
         wz3*[[wy1*[wx1 wx2 wx3 wx4]] [wy2*[wx1 wx2 wx3 wx4]] [wy3*[wx1 wx2 wx3 wx4]] [wy4*[wx1 wx2 wx3 wx4]]]...
         wz4*[[wy1*[wx1 wx2 wx3 wx4]] [wy2*[wx1 wx2 wx3 wx4]] [wy3*[wx1 wx2 wx3 wx4]] [wy4*[wx1 wx2 wx3 wx4]]]
*/
template<class T>
void gsPatchRule<T>::mapTo( const gsVector<T>& lower,
                            const gsVector<T>& upper,
                            gsMatrix<T> & nodes,
                            gsVector<T> & weights ) const
{
    // First, we compute the nodes and weights that are between the lower and upper corners.
    // This number can be different per element
    // we also compute the total number of points in the tensor product (size)
    index_t size = 1;
    std::vector<gsVector<T> > elNodes(m_dim), elWeights(m_dim);
    index_t k = 0; // counter per direction d
    for (size_t d = 0; d!=m_dim; d++)
    {
        elNodes[d].resize(m_nodes[d].size());
        elWeights[d].resize(m_weights[d].size());
        for (typename std::map<T,T>::const_iterator it = m_maps[d].lower_bound(lower[d]); it!= ( upper[d]==m_end[d] ? m_maps[d].end() : m_maps[d].lower_bound(upper[d]) ); it++, k++) // lower_bound = geq, upper_bound= greather than
        {
            elNodes[d].at(k) = it->first;
            elWeights[d].at(k) = it->second;
        }
        elNodes[d].conservativeResize(k);
        elWeights[d].conservativeResize(k);
        size *= elNodes[d].size(); // compute tensor size
        k = 0; // reset counter
    }

    // This could work. However, in cases when there are elements which have no quadPoint, it goes wrong
    // this->computeTensorProductRule_into(elNodes,elWeights,nodes,weights);

    // initialize the number of nodes and weights
    nodes.resize(m_dim,size);
    weights.resize(size);
    if (size==0)
        return;

    // Now we fill the matrix with the points and we construct the tensor product (according to the scheme on top)
    gsMatrix<T> tmpNodes, tmpWeights;
    gsVector<T> ones;
    tmpNodes = elNodes[0].transpose();
    tmpWeights = elWeights[0].transpose();
    size = 1;
    for (size_t d = 1; d!=m_dim; d++)
    {
        nodes.block( 0, 0, d, tmpNodes.cols()*elNodes[d].size() ) = tmpNodes.replicate(1,elNodes[d].size());

        ones.setOnes(tmpNodes.cols());
        for (k = 0; k != elNodes[d].size(); k++)
        {
            nodes.block(d,k*ones.size(),1,ones.size()) = ones.transpose()*elNodes[d].at(k);
            weights.segment(k*ones.size(),ones.size()) = (ones.transpose()*elWeights[d].at(k)).cwiseProduct(tmpWeights);
        }

        tmpWeights.transpose() = weights.segment(0,tmpWeights.cols()*elWeights[d].size());
        tmpNodes = nodes.block( 0, 0, d+1, tmpNodes.cols()*elNodes[d].size() );
    }
};

template<class T> void
gsPatchRule<T>::mapTo( T startVal, T endVal,
                      gsMatrix<T> & nodes, gsVector<T> & weights ) const
{
    GISMO_ASSERT( 1 == m_nodes.size(), "Inconsistent quadrature mapping (dimension != 1)");

    nodes.resize(1,m_nodes[0].size());
    weights.resize(m_weights[0].size());
    index_t k=0;
    for (typename std::map<T,T>::const_iterator it = m_maps[0].lower_bound(startVal); it!= ( endVal==m_end[0] ? m_maps[0].end() : m_maps[0].lower_bound(endVal) ); it++, k++) // lower_bound = geq, upper_bound= greather than
    {
        nodes.at(k) = it->first;
        weights.at(k) = it->second;
    }
    nodes.conservativeResize(1,k);
    weights.conservativeResize(k);
}

template<class T>
gsKnotVector<T> gsPatchRule<T>::_init(const gsBSplineBasis<T> * Bbasis) const
{
    // get the knot vector and the size of the basis
    gsKnotVector<T> knots = Bbasis->knots();
    index_t size = Bbasis->size();

    // check the difference in the current order and the desired order
    index_t pdiff = m_deg - knots.degree();
    // check the difference in the current regularity and the desired regularity
    std::vector<index_t> multiplicities = knots.multiplicities();
    index_t rmin = *std::min_element(multiplicities.begin(), multiplicities.end());
    index_t rdiff = (m_deg-m_reg)-rmin ;

    // Increase order and regularity
    if (pdiff>0)
        knots.degreeIncrease(pdiff);
    else if (pdiff<0)
        knots.degreeDecrease(-pdiff);
    if (rdiff>0)
        knots.increaseMultiplicity(rdiff);
    else if (rdiff>0)
        knots.reduceMultiplicity(-rdiff);

    gsBSplineBasis<T> basis = gsBSplineBasis<T>(knots);
    size = basis.size();

    // If basis should be over-integrated, then add extra knots in the boundary elements
    if (m_over)
    {
        index_t numOver = knots.degree()-1;

        T lowerLength = (knots(1)-knots.first())/(numOver+1);
        T upperLength = (knots.last()-knots(knots.uSize()-2))/(numOver+1);
        for (index_t k=0; k!=numOver; k++)
        {
            knots.insert(knots.first()+(k+1)*lowerLength);
            knots.insert(knots.last()-(k+1)*upperLength);
        }

        size += 2*numOver;
    }

    // Add a middle knot if the size of the knot vector is odd
    if (size % 2 == 1)
    {
        typedef typename gsKnotVector<T>::const_iterator knotIterator;
        knotIterator prevKnot = knots.begin();
        std::vector<T> diff;
        for (knotIterator it = knots.begin()+1; it!=knots.end(); it++ )
        {
            diff.push_back(*it-*prevKnot);
            prevKnot = it;
        }

        T max = *std::max_element(diff.begin(),diff.end());
        std::vector<index_t> maxIdx;
        index_t k=0;
        for (typename std::vector<T>::iterator it = diff.begin(); it!=diff.end(); it++,k++)
        {
            if (math::abs(*it-max)/(max)<1e-15)//warning: fixed acuracy
                maxIdx.push_back(k);
        }

        const index_t i = maxIdx.at(
            static_cast<index_t>(std::ceil(maxIdx.size()/2.))-1);
        knots.insert((knots.at(i) + knots.at(i+1))/2.) ;

        /*  ALTERNATIVE (does not work well)
            T first = knots.first();
            T last  = knots.last();
            T knot = (last - first) / 2.0;

            knots.insert( knot );

            size++;
        */

    }
    return knots;
};

template <class T>
std::pair<gsMatrix<T>,gsVector<T> > gsPatchRule<T>::_integrate(const gsKnotVector<T> & knots ) const
{
    // Obtain a temporary bspline basis and quadrule
    gsBSplineBasis<T> basis = gsBSplineBasis<T>(knots);
    gsQuadRule<T> quRule = gsGaussRule<T>(basis,1,1);


    // obtain the greville points
    gsMatrix<T> greville;
    knots.greville_into(greville);

    gsMatrix<T> integrals(basis.size(),1);
    integrals.setZero();
    gsMatrix<T> nodes, tmp;
    gsMatrix<index_t> actives;
    gsVector<T> weights;

    // obtain the exact integrals of the functions
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    for (; domIt->good(); domIt->next() )
    {
        quRule.mapTo(domIt->lowerCorner(),domIt->upperCorner(),nodes,weights);
        basis.active_into(nodes,actives);
        basis.eval_into(nodes,tmp);
        tmp *= weights;
        for (index_t k = 0; k!=weights.size(); k++)
            integrals(actives.at(k),0) += tmp(k,0);
    }

    return std::make_pair(greville,integrals);
};

template <class T>
std::pair<gsVector<T>,gsVector<T> > gsPatchRule<T>::_compute(const gsKnotVector<T> & knots,
                                                            const gsMatrix<T> & greville,
                                                            const gsVector<T> & integrals,
                                                            const T tol) const
{
    // Initialization
    gsVector<T> nodes, weights;
    gsBSplineBasis<T> basis = gsBSplineBasis<T>(knots);
    gsQuadRule<T> quRule = gsGaussRule<T>(basis,1,1);
    index_t size = basis.size();
    GISMO_ENSURE((size % 2 == 0),"Number of points should be even!");

    // Construct initial guess
    nodes.resize(size/2);
    weights.resize(size/2);
    for (index_t k = 0; k!=size/2; k++)
    {
        nodes.at(k)  = 0.5 * (greville.at(2*k) + greville.at(2*k+1));
        weights.at(k) = integrals(2*k,0) + integrals(2*k+1,0);
    }

    // Initialize Newton iterations
    index_t itMax = 100;
    gsMatrix<index_t> actives;
    gsMatrix<T> vals,dvals,vals_tmp,dvals_tmp,res,dres;
    gsVector<T> update(size);
    vals.resize(size,size/2);   vals.setZero();                 // Jacobian: dF/dw
    dvals.resize(size,size/2);  dvals.setZero();                // Jacobian: dF/dxi
    index_t it;

    // Newton Iterations
    for (it = 0; it != itMax; it++)
    {
        // Compute matrix with values and derivatives of the basis on the quadrature point estimates
        basis.active_into(nodes.transpose(),actives);
        basis.eval_into(nodes.transpose(),vals_tmp);
        basis.deriv_into(nodes.transpose(),dvals_tmp);
        for (index_t act=0; act!=actives.rows(); act++)
            for (index_t pt=0; pt!=actives.cols(); pt++)
            {
                vals(actives(act,pt),pt)    = vals_tmp(act,pt); // Jacobian: dF/dw
                dvals(actives(act,pt),pt)   = dvals_tmp(act,pt);// Jacobian: dF/dxi
            }

        // Compute the residual (res) and the Jacobian (dres)
        res = vals * weights - integrals;                       // Residual
        dres = gsMatrix<T>::Zero(size,size);
        dres.block(0,0,size,size/2) = vals;
        dres.block(0,size/2,size,size/2) = dvals * weights.asDiagonal();

        // Compute the Newton update
        m_solver.compute(dres.sparseView());
        update = m_solver.solve(-res);

        // Update the weights and the nodes
        weights += update.head(size/2);
        nodes  += update.tail(size/2);

        GISMO_ENSURE(nodes.minCoeff() > knots.first(),"Construction failed: min(nodes) < min(knots) minCoef = "<<nodes.minCoeff()<<"; min(knots) = "<<knots.first());
        GISMO_ENSURE(nodes.maxCoeff() < knots.last(),"Construction failed: max(nodes) > max(knots) maxCoef = "<<nodes.maxCoeff()<<"; min(knots) = "<<knots.last());

        // Check convergence
        if (res.norm() < tol)
            break;
    }
    GISMO_ENSURE(it+1!=itMax,"Maximum iterations reached");
    return std::make_pair(nodes,weights);
}

} // namespace gismo
