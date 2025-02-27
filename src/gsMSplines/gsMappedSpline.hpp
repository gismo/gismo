/** @file gsMappedSpline.hpp

    @brief implementation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsMSplines/gsMappedSpline.h>
#include <gsCore/gsDofMapper.h>

namespace gismo
{

template<short_t d,class T>
gsMappedSpline<d,T>::gsMappedSpline( const gsMultiPatch<T> & mp, const gsSparseMatrix<T> & m )
{
    GISMO_ASSERT(mp.nPatches()>0,"MultiPatch is empty?");
    m_mbases = new gsMappedBasis<d,T>(gsMultiBasis<T>(mp),m);

    // collect and transform the coefficients
    gsMatrix<T> local = mp.coefs();
    m_mbases->local_coef_to_global_coef(local,m_global);

    init(*m_mbases);
}

template<short_t d,class T>
gsMappedSpline<d,T>::gsMappedSpline( const gsMappedBasis<d,T> & mbases, const gsMatrix<T> & coefs )
:
m_global(coefs)
{
    m_mbases=mbases.clone().release();
    init(mbases);
}

template<short_t d,class T>
gsMappedSpline<d,T>::gsMappedSpline( const gsMappedSpline& other )
: gsFunctionSet<T>(), m_global(other.m_global)
{
    m_mbases=other.m_mbases->clone().release();
}

template<short_t d,class T>
gsMappedSpline<d,T> & gsMappedSpline<d,T>::operator=( const gsMappedSpline& other )
{
    delete m_mbases;
    m_mbases=other.m_mbases->clone().release();
    m_global = other.m_global;
    m_ss = other.m_ss;
    for (auto & s : m_ss) s.setSource(*this);
    return *this;
}

template<short_t d,class T>
void gsMappedSpline<d,T>::eval_into(const unsigned patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    gsMatrix<index_t> actives;
    gsMatrix<T> evals;

    index_t N = targetDim();
    // gsC1Basis<d,t> * basis = dynamic_cast<gsC1Basis<d,t> *>(this->getBase(patch));
    // if (basis==NULL)
    // {
        // m_mbases->active_into(patch,u,actives);
        // m_mbases->eval_into(patch,u,evals);
        // m_mbases->getBase(patch).linearCombination_into(m_coefs,actives,evals,result);
    // }
    // else
    // {
        gsMatrix<T> tmp;
        result.resize( N,u.cols());
        // This loop enables that the number of actives can be different for each column in u
        for (index_t k = 0; k!=u.cols(); k++)
        {
            m_mbases->active_into(patch,u.col(k),actives);
            m_mbases->eval_into(patch,u.col(k),evals);
            m_mbases->getBase(patch).linearCombination_into(m_global,actives,evals,tmp);
            result.col(k) = tmp;
        }
    // }
}

template<short_t d,class T>
void gsMappedSpline<d,T>::deriv_into(const unsigned patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    gsMatrix<index_t> actives;
    gsMatrix<T> evals;

    index_t N = targetDim();
    index_t M = domainDim();
    result.resize( N * M,u.cols());

    gsMatrix<T> tmp;
    // This loop enables that the number of actives can be different for each column in u
    for (index_t k = 0; k!=u.cols(); k++)
    {
        m_mbases->active_into(patch,u.col(k),actives);
        m_mbases->deriv_into(patch,u.col(k),evals);
        m_mbases->getBase(patch).linearCombination_into(m_global,actives,evals,tmp);
        result.col(k) = tmp;
    }
}

template<short_t d,class T>
void gsMappedSpline<d,T>::deriv2_into(const unsigned patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    gsMatrix<index_t> actives;
    gsMatrix<T> evals;

    index_t N = targetDim();
    index_t M = domainDim();
    index_t S = M*(M+1)/2;
    result.resize( S*N,u.cols());
    gsMatrix<T> tmp;
    // This loop enables that the number of actives can be different for each column in u
    for (index_t k = 0; k!=u.cols(); k++)
    {
        m_mbases->active_into(patch,u.col(k),actives);
        m_mbases->deriv2_into(patch,u.col(k),evals);
        m_mbases->getBase(patch).linearCombination_into(m_global,actives,evals,tmp);
        result.col(k) = tmp;
    }
}

template<short_t d,class T>
void gsMappedSpline<d,T>::evalAllDers_into(const unsigned patch, const gsMatrix<T> & u,
                                           const int n, std::vector<gsMatrix<T> >& result,
                                           bool sameElement) const
{
    result.resize(n+1);

    gsMatrix<index_t> actives;
    std::vector< gsMatrix<T> > evals;

    index_t N = targetDim();
    index_t m = domainDim();
    index_t S = m*(m+1)/2;

    std::vector<index_t> blocksizes(3);
    blocksizes[0] = N;
    blocksizes[1] = N*m;
    blocksizes[2] = N*S;

    gsMatrix<T> tmp;
    // todo: change the loop over i and the loop over k
    for( int i = 0; i <= n; i++)
    {
        result[i].resize(blocksizes[i],u.cols());
        // This loop enables that the number of actives can be different for each column in u
        for (index_t k = 0; k!=u.cols(); k++)
        {
            m_mbases->active_into(patch,u.col(k),actives);//..
            m_mbases->evalAllDers_into(patch,u.col(k),n,evals,sameElement);
            m_mbases->getBase(patch).linearCombination_into(m_global,actives,evals[i],tmp);
            result[i].col(k) = tmp;

        }
    }

    // for( int i = 0; i <= n; i++)
    //     result[i].resize(blocksizes[i],u.cols());

    // // todo: change the loop over i and the loop over k
    // for (index_t k = 0; k!=u.cols(); k++)
    // {
    //     // This loop enables that the number of actives can be different for each column in u
    //     m_mbases->active_into(patch,u.col(k),actives);
    //     m_mbases->evalAllDers_into(patch,u.col(k),n,evals);
    //     for( int i = 0; i <= n; i++)
    //     {
    //         m_mbases->getBase(patch).linearCombination_into(m_global,actives,evals[i],tmp);
    //         result[i].col(k) = tmp;
    //     }
    // }
}

template<short_t d,class T>
gsMultiPatch<T> gsMappedSpline<d,T>::exportToPatches() const
{
    gsMatrix<T> local;
    m_mbases->global_coef_to_local_coef(m_global,local);
    return m_mbases->exportToPatches(local);
}

template<short_t d,class T>
gsGeometry<T> * gsMappedSpline<d,T>::exportPatch(int i,gsMatrix<T> const & localCoef) const
{
    return m_mbases->exportPatch(i,localCoef);
}

template<short_t d,class T>
std::map<std::array<size_t, 4>, internal::ElementBlock> gsMappedSpline<d,T>::BezierOperator() const
{
    GISMO_ENSURE( 2==domainDim(), "Anything other than bivariate splines is not yet supported!");

    // Loop over all the elements of the given Mapped Spline and collect all relevant
    // information in ElementBlocks. These will be grouped in a std::map
    // with respect to the number of active basis functions ( = NNj/NCV )
    // of each Bezier element
    std::map<std::array<size_t, 4>, internal::ElementBlock> ElementBlocks;

    // index_t NNj; // Number of control points of the Bezier elements of block j
    gsMatrix<index_t> mappedActives, globalActives, localActives; // Active basis functions

    gsMatrix<T> coefVectors, center_point(2,1);
    // Center point of the Bezier element, all basis functions are active here
    center_point.setConstant( (T)(0.5) );

    // The control points of the MappedSpline
    gsMatrix<T> globalCoefs = this->getMappedCoefs();

    // The transpose of the global Bezier Operator
    gsWeightMapper<T> mapper = this->getMapper();

    // The Bezier operator, transpose is used because we are constructing
    // the basis functions, the mapper.asMatrix() is used for the control points
    gsMatrix<> bezOperator = mapper.asMatrix().transpose();

// TODO: Delete if implementation is OK as is
/*     gsMultiPatch<> mp = this->exportToPatches();
    mp.computeTopology();
    // gsMultiBasis<> mb(mp);
    gsDofMapper dofMapper = mp.getMapper(1e-4);
    // dofMapper.finalize(); */

    // The mapped spline's basis functions
    gsMappedBasis<d,T>  mappedBasis = this->getMappedBasis();

    std::array<size_t, 4> key;
    for (size_t p=0; p<this->nPatches(); ++p)
    {
        // Mapped active basis functions
        mappedBasis.active_into(p, center_point,mappedActives);

        // Patch-local active basis functions
        this->getBase(p).active_into(center_point,localActives);
        // gsInfo << "Local ("<< localActives.rows()<<"):" << localActives.transpose() << "\n";

        globalActives.resizeLike(localActives);
        globalActives.setZero();

        // OPTION 1: Global, all considered uncoupled
        // Get the global indices of the active basis functions
        globalActives = mappedBasis.getGlobalIndex(p,localActives);

        // OPTION 2: Global, considering coupling (i.e. as in the multi patch)
        // globalActives.resizeLike(localActives);
        // for (index_t i=0; i<localActives.rows(); ++i)
        // {
        //     globalActives(i) = dofMapper.index( localActives(i), p);
        // }

        //OPTION 3: Global, considering coupling and mapping
        std::vector<index_t> sourceID, mappedID;
        // mappedActives.resizeLike(globalActives);
        // for (index_t i=0; i<localActives.rows(); ++i)
        // {
        //     sourceID.clear();
        //     mappedID.clear();
        //     sourceID.push_back( globalActives(i) );
        //     mapper.sourceToTarget(sourceID, mappedID);
        //     mappedActives(i) = mappedID.front();
        // }


        // gsInfo << "Global("<< globalActives.rows()<<"):" << globalActives.transpose() << "\n";
        // gsInfo << "Mapped("<< mappedActives.rows()<<"):" << mappedActives.transpose() <<"\n";


        // coefVectors is the patch-local Bezier operator, size (NNj x NCVCj), NCVCj = (PR+1)*(PS+1)*(PT+1)
        coefVectors.resize(mappedActives.rows(), globalActives.rows());
        coefVectors.setZero();

        // std::vector<std::pair<index_t,index_t> > preImage;
        // Extract the local bezier operator from the global one
        for (index_t i=0; i<mappedActives.rows(); ++i)
        {
            for (index_t j=0; j<globalActives.rows(); ++j)
            {
                // if (dofMapper.is_coupled_index(globalActives(j)))
                // {
                //     dofMapper.preImage(globalActives(j), preImage);
                //     // gsInfo << globalActives(j) << " is coupled with ";
                //     for (size_t k = 0; k != preImage.size(); ++k)
                //         // gsInfo << preImage.at(k).first << "," << preImage.at(k).second << " ";
                //         coefVectors(i,j) += bezOperator(mappedActives(i),  mappedBasis.getGlobalIndex(preImage.at(k).first, preImage.at(k).second) );
                //     // gsInfo << "\n";
                // }
                // else
                // {
                //     // gsInfo << globalActives(j) << " is free.\n";
                //     coefVectors(i,j) = bezOperator(mappedActives(i), mappedBasis.getGlobalIndex(p, localActives(j))) ;
                // }

                // sourceID.clear();
                // mappedID.clear();
                // sourceID.push_back( globalActives(j) );
                // mapper.sourceToTarget(sourceID, mappedID);
                coefVectors(i,j) = bezOperator(mappedActives(i), globalActives(j)) ;
            }

        }
        // gsInfo << "Coefs size:" << coefVectors.rows() << "x" << coefVectors.cols() << "\n\n";

        // Put everything in the ElementBlock
        key[0] = mappedActives.rows();
        key[1] = this->getBase(p).degree(0);
        key[2] = this->getBase(p).degree(1);
        key[3] = 0; // TODO: if implemented for trivariates fix this
        // NNj = mappedActives.size();             // Number of control points (nodes) of the Bezier element
        ElementBlocks[key].numElements += 1;    // Increment the Number of Elements contained in the ElementBlock
        ElementBlocks[key].actives.push_back(mappedActives);  // Append the active basis functions ( = the Node IDs ) for this element.
        ElementBlocks[key].PR = this->getBase(p).degree(0);   // Degree of the Bezier element in the r-direction
        ElementBlocks[key].PS = this->getBase(p).degree(1);   // Degree of the Bezier element in the s-direction
        ElementBlocks[key].PT = 0; // TODO: if implemented for trivariates fix this
        ElementBlocks[key].coefVectors.push_back( coefVectors ); //(!) If it is not  bezier
    }

    return ElementBlocks;
}


} //namespace gismo
