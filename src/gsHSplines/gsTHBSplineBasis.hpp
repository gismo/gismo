/** @file gsTHBSplineBasis.hpp

    @brief Provides implementation of THBSplineBasis class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once 

#include <gsCore/gsMultiPatch.h>

#include <gsNurbs/gsBoehm.h>
#include <gsNurbs/gsDeboor.hpp>

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>

#include <gsTensor/gsTensorTools.h>

namespace gismo
{

template<short_t d, class T>
typename gsTHBSplineBasis<d,T>::BoundaryBasisType * gsTHBSplineBasis<d,T>::basisSlice(index_t dir_fixed,T par ) const
{
    GISMO_ASSERT(d-1>=0,"d must be greater or equal than 1");
    GISMO_ASSERT(dir_fixed>=0 && static_cast<index_t>(dir_fixed)<d,"cannot fix a dir greater than dim or smaller than 0");
    const boxSide side(dir_fixed,0);
    const typename gsTensorBSplineBasis<d,T>::BoundaryBasisType::uPtr bBSplineBasis =
        this->m_bases[0]->boundaryBasis(side);
    typename gsTHBSplineBasis<d,T>::BoundaryBasisType* bBasis =
        new typename gsTHBSplineBasis<d,T>::BoundaryBasisType(*bBSplineBasis);//,this->m_tree.getMaxInsLevel()+1);

    if(d!=1)
    {
        std::vector<index_t> boxes;
        this->getBoxesAlongSlice(dir_fixed,par,boxes);
        bBasis->refineElements(boxes);
    }
    return bBasis;
}

template<short_t d, class T>
void gsTHBSplineBasis<d,T>::representBasis()
{
    // Cleanup previous basis
    this->m_is_truncated.resize(this->size());
    m_presentation.clear();

    for (index_t j = 0; j < this->size(); ++j)
    {
        unsigned level = this->levelOf(j);
        index_t tensor_index = this->flatTensorIndexOf(j, level);

        // element indices
        gsMatrix<index_t, d, 2> element_ind(d, 2);
        this->m_bases[level]->elementSupport_into(tensor_index, element_ind);

        // I tried with block, I can not trick the compiler to use references
        gsVector<index_t, d> low = element_ind.col(0); //block<d, 1>(0, 0);
        gsVector<index_t, d> high = element_ind.col(1); //block<d, 1>(0, 1);gsMatrix<index_t> element_ind =

        // Finds coarsest level that function, with supports given with
        // support indices of the coarsest level (low & high), has presentation
        // based only on B-Splines (and not THB-Splines).
        // this is not the same as query 3
        unsigned clevel = this->m_tree.query4(low, high, level);

        if (level != clevel) // we must compute its presentation
        {
            this->m_tree.computeFinestIndex(low, level, low);
            this->m_tree.computeFinestIndex(high, level, high);

            this->m_is_truncated[j] = clevel;
            _representBasisFunction(j, clevel, low, high);
        }
        else
        {
            this->m_is_truncated[j] = -1;
        }
    }
}

template<short_t d, class T>
void gsTHBSplineBasis<d,T>::_representBasisFunction(
    const unsigned j,
    const unsigned pres_level,
    const gsVector<index_t, d>& finest_low,
    const gsVector<index_t, d>& finest_high)
{
    const unsigned cur_level = this->levelOf(j);

    // actual size of the coefficients
    gsVector<index_t, d> act_size_of_coefs(d);
    act_size_of_coefs.fill(1);

    // number of new coefficients
    unsigned nmb_of_coefs = _updateSizeOfCoefs(cur_level, pres_level,
                                               finest_low, finest_high,
                                               act_size_of_coefs);
    gsVector<T> one(1);
    one(0) = 1.0;
    gsMatrix<T> coefs(nmb_of_coefs, 1);
    coefs.fill(0);
    coefs.row(0) = one;

    // vector of the numbers of the coefficients (in each dimension)
    // stored in coefs
    gsVector<index_t, d> vec_nmb_of_coefs(d);
    vec_nmb_of_coefs.fill(1);

    unsigned tensor_index = this->flatTensorIndexOf(j, cur_level);

    // B-Spline vector tensor index
    gsVector<index_t, d> bspl_vec_ti =
        this->m_bases[cur_level]->tensorIndex(tensor_index);


    // we need to separately save knot vectors because we will modify
    // them, when we proceed from one level on another
    std::vector<gsKnotVector<T> > vector_of_kv(d);

    // size of the coefficients that are affected in individual iteration
    gsVector<index_t, d> cur_size_of_coefs(d);
    cur_size_of_coefs.fill(1);

    for (unsigned level = cur_level; level < pres_level; ++level)
    {
        _updateSizeOfCoefs(level, level + 1, finest_low,
                           finest_high, cur_size_of_coefs);


        // index of a support of the j-th basis function (l_low, l_high
        // on level, and l1_high, l1_low on level + 1)
        gsVector<index_t, d> clow, chigh, fhigh, flow;

        this->m_tree.computeLevelIndex(finest_low, level, clow);
        this->m_tree.computeLevelIndex(finest_high, level, chigh);

        this->m_tree.computeLevelIndex(finest_low, level + 1, flow);
        this->m_tree.computeLevelIndex(finest_high, level + 1, fhigh);

        std::vector<T> knots;

        for (unsigned dim = 0; dim < d; ++dim)
        {
            const gsKnotVector<T>& ckv = m_bases[level  ]->knots(dim);
            const gsKnotVector<T>& fkv = m_bases[level+1]->knots(dim);

            if (level == cur_level)
                vector_of_kv[dim] = ckv;

            knots.clear();
            std::set_symmetric_difference(ckv.beginAt(clow[dim]), ckv.endAt(chigh[dim]),
                                          fkv.beginAt(flow[dim]), fkv.endAt(fhigh[dim]),
                                          std::back_inserter(knots));

            gsTensorBoehmRefineLocal<d,
                                     gsKnotVector<T>,
                                     gsMatrix<T>,
                                     typename std::vector<T>::const_iterator>
                (vector_of_kv[dim], bspl_vec_ti[dim], coefs, vec_nmb_of_coefs,
                 act_size_of_coefs,
                 cur_size_of_coefs, dim, knots.begin(), knots.end(),
                 true);
        }

        _truncate(coefs, act_size_of_coefs, cur_size_of_coefs,
                  level + 1, bspl_vec_ti, cur_level, finest_low);
    }

    _saveNewBasisFunPresentation(coefs, act_size_of_coefs,
                                 j, pres_level, finest_low);
}


template<short_t d, class T>
void gsTHBSplineBasis<d,T>::_saveNewBasisFunPresentation(
    const gsMatrix<T>& coefs,
    const gsVector<index_t, d>& act_size_of_coefs,
    const unsigned j,
    const unsigned pres_level,
    const gsVector<index_t, d>& finest_low)
{
    const unsigned level = this->levelOf(j);
    const unsigned tensor_index = this->flatTensorIndexOf(j, level);

    gsVector<index_t, d> bspl_vec_ti =
        this->m_bases[level]->tensorIndex(tensor_index);

    // finer tensor index
    const unsigned f_ten_index = _basisFunIndexOnLevel(bspl_vec_ti, level,
                                                       finest_low, pres_level);

    gsVector<index_t, d> act_coefs_strides(d);
    bspline::buildCoeffsStrides<d>(act_size_of_coefs, act_coefs_strides);


    gsVector<index_t, d> position(d);
    position.fill(0);


    gsVector<index_t, d> first_point(position);
    gsVector<index_t, d> last_point(d);
    bspline::getLastIndexLocal<d>(act_size_of_coefs, last_point);

    this->m_presentation[j] =
        gsSparseVector<T>(this->m_bases[pres_level]->size());

    gsSparseVector<T>& presentation = this->m_presentation[j];
    //presentation.reserve(coefs.rows()); // reserve memory ?

    do
    {
        // ten_index - (tensor) index of a bspline function with respect to
        //             the coef at position "position"
        // coef_index - (local) index of a ceofficient at position

        unsigned ten_index = f_ten_index;
        for (unsigned dim = 0; dim < d; dim++)
        {
            ten_index += position(dim) *
                this->m_bases[pres_level]->stride(static_cast<int>(dim));
        }

        unsigned coef_index = bspline::getIndex<d>(act_coefs_strides, position);

        if (coefs(coef_index) != 0)
            presentation(ten_index) = coefs(coef_index);

    } while(nextCubePoint<gsVector<index_t, d> > (position, first_point,
                                                   last_point));
}


template<short_t d, class T>
unsigned gsTHBSplineBasis<d,T>::_basisFunIndexOnLevel(
    const gsVector<index_t, d>& index,
    const unsigned level,
    const gsVector<index_t, d>& fin_low,
    const unsigned new_level)
{
    gsVector<index_t, d> low(d);
    this->m_tree.computeLevelIndex(fin_low, level, low);

    gsVector<index_t, d> flow(d);
    this->m_tree.computeLevelIndex(fin_low, new_level, flow);

    gsVector<index_t, d> new_index(d);

    for (unsigned dim = 0; dim < d; dim++)
    {
        const gsKnotVector<T>& ckv =
            this->m_bases[level]->knots(dim);

        const gsKnotVector<T>& fkv =
            this->m_bases[new_level]->knots(dim);


        unsigned mult = index[dim] - ckv.firstKnotIndex(low[dim]);

        new_index(dim) = fkv.firstKnotIndex(flow[dim]) + mult;
    }

    return this->m_bases[new_level]->index(new_index);
}


template<short_t d, class T>
void gsTHBSplineBasis<d,T>::_truncate(
    gsMatrix<T>& coefs,
    const gsVector<index_t, d>& act_size_of_coefs,
    const gsVector<index_t, d>& size_of_coefs,
    const unsigned level,
    const gsVector<index_t, d>& bspl_vec_ti,
    const unsigned bspl_vec_ti_level,
    const gsVector<index_t, d>& finest_low)
{
    // if we dont have any active function in this level, we do not truncate
    if (this->m_xmatrix[level].size() == 0)
        return;


    // global tensor index
    const unsigned const_ten_index = _basisFunIndexOnLevel(bspl_vec_ti,
                                                           bspl_vec_ti_level, finest_low, level);
    gsVector<index_t, d> act_coefs_strides(d);
    bspline::buildCoeffsStrides<d>(act_size_of_coefs, act_coefs_strides);


    gsVector<index_t, d> last_point(d);
    bspline::getLastIndexLocal<d>(size_of_coefs, last_point);
    last_point(0) = 0;


    gsVector<index_t, d> position(d);
    position.fill(0);

    gsVector<index_t, d> first_point(position);

    unsigned xmatrix_index = 0;
    unsigned tensor_active_index = this->m_xmatrix[level][0];

    unsigned numb_of_point = size_of_coefs[0];

    do
    {
        // indices
        // ten_index - (tensor) index of a bspline function with respect to
        //             the coef at position "position"
        // coef_index - (local) index of a ceofficient at position
        // tensor_active_index - (tensor) index of a active functions


        unsigned ten_index = const_ten_index;
        for (unsigned dim = 1; dim < d; dim++)
        {
            ten_index += position(dim) *
                this->m_bases[level]->stride(static_cast<int>(dim));
        }

        unsigned coef_index = bspline::getIndex<d>(act_coefs_strides, position);

        for (unsigned index = 0; index < numb_of_point; ++index)
        {
            if (tensor_active_index < ten_index)
            {
                while(tensor_active_index < ten_index)
                {
                    xmatrix_index++;
                    if (xmatrix_index == this->m_xmatrix[level].size())
                    {
                        // we don't have any active basis function,
                        // so all the rest basis function in our
                        // representation should not be truncated
                        return;
                    }

                    tensor_active_index =
                        this->m_xmatrix[level][xmatrix_index];
                }
                // ten_index <= tensor_active_index holds
            }

            if (ten_index == tensor_active_index) // truncate
                coefs(coef_index + index, 0) = 0;

            ten_index++;
        }

    } while(nextCubePoint<gsVector<index_t, d> >(position, first_point,
                                                  last_point));
}


template<short_t d, class T>
unsigned gsTHBSplineBasis<d,T>::_updateSizeOfCoefs(
    const unsigned clevel,
    const unsigned flevel,
    const gsVector<index_t, d>& finest_low,
    const gsVector<index_t, d>& finest_high,
    gsVector<index_t, d>& size_of_coefs)
{
    gsVector<index_t, d> clow, chigh;

    this->m_tree.computeLevelIndex(finest_low, clevel, clow);
    this->m_tree.computeLevelIndex(finest_high, clevel, chigh);

    gsVector<index_t, d> flow, fhigh;
    this->m_tree.computeLevelIndex(finest_low, flevel, flow);
    this->m_tree.computeLevelIndex(finest_high, flevel, fhigh);

    // number of new coefficients
    unsigned nmb_of_coefs = 1;

    for (unsigned dim = 0; dim < d; ++dim)
    {
        const gsKnotVector<T>& ckv =
            this->m_bases[clevel]->knots(dim);
        const gsKnotVector<T>& fkv =
            this->m_bases[flevel]->knots(dim);

        unsigned cnmb_knts = ckv.knotsUntilSpan(chigh[dim]) -
            ckv.knotsUntilSpan(clow[dim]);

        unsigned fnmb_knts = fkv.knotsUntilSpan(fhigh[dim]) -
            fkv.knotsUntilSpan(flow[dim]);

        size_of_coefs(dim) += fnmb_knts - cnmb_knts;
        nmb_of_coefs *= size_of_coefs(dim);
    }

    return nmb_of_coefs;
}

// return the B-spline representation of a THB-spline subpatch
template<short_t d, class T>
void gsTHBSplineBasis<d,T>::getBsplinePatchGlobal(gsVector<index_t> b1,
                                                  gsVector<index_t> b2,
                                                  unsigned level, 
                                                  const gsMatrix<T>& geom_coef,
                                                  gsMatrix<T>& cp,
                                                  gsKnotVector<T>& k1,
                                                  gsKnotVector<T>& k2) const
{    
    // check if the indices in b1, and b2 are correct with respect to the given level    
    const unsigned loc2glob = ( 1<< (this->maxLevel() - level) );
    if( b1[0]%loc2glob != 0 ) b1[0] -= b1[0]%loc2glob;
    if( b1[1]%loc2glob != 0 ) b1[1] -= b1[1]%loc2glob;
    if( b2[0]%loc2glob != 0 ) b2[0] += loc2glob -(b2[0]%loc2glob);
    if( b2[1]%loc2glob != 0 ) b2[1] += loc2glob -(b2[1]%loc2glob);

    // select the indices of all B-splines of the given level acting on the given box
    gsVector<index_t,d> b1_outputs, b2_outputs;
    this->m_tree.computeLevelIndex( b1, level, b1_outputs );
    this->m_tree.computeLevelIndex( b2, level, b2_outputs );
    int i0 = b1_outputs(0);
    int i1 = b2_outputs(0);
    int j0 = b1_outputs(1);
    int j1 = b2_outputs(1);
    i0 = m_bases[level]->knots(0).lastKnotIndex(i0) - m_deg[0];
    i1 = m_bases[level]->knots(0).firstKnotIndex(i1) - 1;
    j0 = m_bases[level]->knots(1).lastKnotIndex(j0) - m_deg[1];
    j1 = m_bases[level]->knots(1).firstKnotIndex(j1) - 1;

    const index_t sz0   = m_bases[level]->size(0);
    const index_t newSz = (i1 - i0 + 1)*(j1 - j0 + 1);
    cp.resize(newSz, geom_coef.cols());

    gsMatrix<T> temp;
    globalRefinement(geom_coef, level, temp);

    index_t cc = 0;
    for(int j = j0; j <= j1; j++)
        for(int k = i0; k <= i1; k++)
            cp.row(cc++) = temp.row(j*sz0+k);
    
    // compute the new vectors for the B-spline patch
    k1 = gsKnotVector<T>(m_deg[0], m_bases[level]->knots(0).begin() + i0 , 
                         m_bases[level]->knots(0).begin() + i1 + m_deg[0] + 2);
    k2 = gsKnotVector<T>(m_deg[1], m_bases[level]->knots(1).begin() + j0 , 
                         m_bases[level]->knots(1).begin() + j1 + m_deg[1] + 2);
}

// returns the list of B-spline patches to represent a THB-spline geometry
template<short_t d, class T>
void gsTHBSplineBasis<d,T>::getBsplinePatches(const gsMatrix<T>& geom_coef, gsMatrix<T>& cp,
                                              gsMatrix<index_t>& b1, gsMatrix<index_t>& b2,
                                              gsVector<index_t>& level, gsMatrix<index_t>& nvertices) const
{ 
    this->m_tree.getBoxes(b1,b2,level); // splitting based on the quadtree
    int nboxes = level.size();
    //------------------------------------------------------------------------------------------------------------------------------
    // iteration on the boxes to call getBsplinePatchGlobal()
    //------------------------------------------------------------------------------------------------------------------------------
    gsVector<index_t> p1, p2;
    p1.resize(this->dim());
    p2.resize(this->dim());
    gsMatrix<T> temp1, temp2;
    gsKnotVector<T> cku, ckv;
    nvertices.resize(nboxes,this->dim());

    for (int i = 0; i < nboxes; i++)
    {
        p1 = b1.row(i).transpose();
        p2 = b2.row(i).transpose();

        this->getBsplinePatchGlobal(p1, p2, level[i], geom_coef, temp1, cku, ckv);        

        if (i == 0)
        {
            cp = temp1;
        }
        else
        {
            int cprows = cp.rows();
            temp2.resize(cp.rows()+temp1.rows(),cp.cols());

            for (int j = 0; j < cp.rows(); j++)
            {
                for (int k = 0; k < cp.cols(); k++)
                {
                    temp2(j,k) = cp(j,k);
                }
            }

            for (int j = 0; j < temp1.rows(); j++)
            {
                for (int k = 0; k < temp1.cols(); k++)
                {
                    temp2(cprows+j,k) = temp1(j,k);
                }
            }
            cp = temp2;
        }

        nvertices(i,0) = cku.size()-cku.degree()-1;
        nvertices(i,1) = ckv.size()-ckv.degree()-1;
    }
}

// returns the list of B-spline patches to represent a THB-spline geometry
template<short_t d, class T>
gsMultiPatch<T> gsTHBSplineBasis<d,T>::getBsplinePatchesToMultiPatch(const gsMatrix<T>& geom_coef) const
{
    GISMO_ASSERT(d==2,"Dim must be 2 for now");

    gsMultiPatch<T> result;

    gsMatrix<index_t> b1, b2;
    gsVector<index_t> level;
    this->m_tree.getBoxes(b1,b2,level); // splitting based on the quadtree

    const int nboxes = level.size();
    gsVector<index_t> p1, p2;
    gsMatrix<T> temp1;
    gsKnotVector<T> cku, ckv;

    for (int i = 0; i < nboxes; i++) // for all boxes
    {
        p1 = b1.row(i).transpose();
        p2 = b2.row(i).transpose();

        this->getBsplinePatchGlobal(p1, p2, level[i], geom_coef, temp1, cku, ckv);
        result.addPatch(typename gsTensorBSpline<2, T>::uPtr(new gsTensorBSpline<2, T>(cku, ckv, give(temp1))));
    }

    return result;
}

// /*
template<short_t d, class T>
void gsTHBSplineBasis<d,T>::getConnectedComponents(
    std::vector<std::vector<std::vector< std::vector<index_t> > > >& connectedComponents, gsVector<index_t>& level) const
{
    //identify the outer polylines- conected components
    int first_level = 0;
    for(unsigned int i = 0; i < this->m_xmatrix.size(); i++)
    {
        if(this->m_xmatrix[i].size()>0)
        {
            first_level = i;
            break;
        }
    }
    gsDebug<<"new min level"<<"\n";
    std::vector< std::vector< std::vector< std::vector<index_t> > > > res; //things to assign to trim_curves
    std::vector< std::vector< std::vector<index_t > > > aabb;//axis aligned bounding box
    std::vector< std::vector<index_t > > boxes;
    aabb = this->domainBoundariesIndices(res);
    for(unsigned int i = 0; i < aabb.size(); i++)//level
    {
        //compare every aabb with the others
        for(unsigned int j = 0; j < aabb[i].size(); j++)//aabb
        {
            bool is_boundary = true;
            for(unsigned int k = 0; k < aabb[i].size(); k++)//aabb
            {
                if(j!=k)
                {
                    if(     (aabb[i][j][0] > aabb[i][k][0]) &&
                            (aabb[i][j][2] < aabb[i][k][2]) &&
                            (aabb[i][j][1] > aabb[i][k][1]) &&
                            (aabb[i][j][3] < aabb[i][k][3]) )
                    {
                        is_boundary= !is_boundary;
                    }
                }
            }
            if(is_boundary)
            {
                boxes.push_back(aabb[i][j]);
                boxes[boxes.size()-1].push_back(i+first_level);//adding the level information
                connectedComponents.resize(connectedComponents.size()+1);
                connectedComponents[connectedComponents.size()-1].push_back(res[i][j]);//trim_curves
            }
        }
    }
    //create the matrices b1,b2 and vector level

    level.resize(boxes.size());
    for(size_t i = 0; i < boxes.size(); i++){
        level[i] = boxes[i][2*d];
    }

    // identify holes
    for(size_t l = 0; l < aabb.size();l++) //level
    {
        for(size_t i = 0; i < aabb[l].size(); i++) //bb
        {
            int closest_box = -1;
            for(size_t j = 0; j < boxes.size(); j++)
            {
                if(l == static_cast<size_t>(boxes[j][4]) )
                {
                    if(     (aabb[l][i][0] > boxes[j][0]) &&
                            (aabb[l][i][1] > boxes[j][1]) &&
                            (aabb[l][i][2] < boxes[j][2]) &&
                            (aabb[l][i][3] < boxes[j][3]))
                    {
                        if(closest_box!=-1){
                            //test with previous closest_box
                            if(!(
                                   (boxes[closest_box][0] > boxes[j][0]) &&
                                   (boxes[closest_box][1] > boxes[j][1]) &&
                                   (boxes[closest_box][2] < boxes[j][2]) &&
                                   (boxes[closest_box][3] < boxes[j][3])))
                            {
                                closest_box = j;
                            }
                        }
                        else
                        {
                            closest_box = j;
                        }

                    }
                    else if(    (aabb[l][i][0] == boxes[j][0]) &&
                                (aabb[l][i][1] == boxes[j][1]) &&
                                (aabb[l][i][2] == boxes[j][2]) &&
                                (aabb[l][i][3] == boxes[j][3]))
                    {
                        //boxes[j] == aabb[l][i]
                        closest_box = -1;
                        break;
                    }
                }

            }
            //prirad s res do trim_curves
            if(closest_box>-1)
            {
                connectedComponents[closest_box].push_back(res[l][i]);
            }
        }
    }


}
//*/


//return data for trimming in parasolid
template<short_t d, class T>
void gsTHBSplineBasis<d,T>::getBsplinePatches_trimming(
    const gsMatrix<T>& geom_coef,
    gsMatrix<T>& cp,
    gsMatrix<index_t>& b1,
    gsMatrix<index_t>& b2,
    gsVector<index_t>& level,
    gsMatrix<index_t>& nvertices,
    std::vector<std::vector<std::vector< std::vector<T> > > >& trim_curves) const
{
    //identify the outer polylines- conected components
//     int first_level = 0;
//     for(unsigned int i = 0; i < this->m_xmatrix.size(); i++){
//         if(this->m_xmatrix[i].size()>0){
//             first_level = i - 1; 
//             break;
//         }
//     }
    
    // gsDebug<<"new min level"<< first_level << "\n";
    std::vector< std::vector< std::vector< std::vector< T > > > > res; //things to assign to trim_curves
    std::vector< std::vector< std::vector<index_t> > > aabb;//axis aligned bounding box
    std::vector< std::vector<index_t> > boxes;
    aabb = this->domainBoundariesParams(res);
    for(unsigned int i = 0; i < aabb.size(); i++)//level
    {
        //compare every aabb with the others
        for(unsigned int j = 0; j < aabb[i].size(); j++)//aabb
        {
            bool is_boundary = true;
            for(unsigned int k = 0; k < aabb[i].size(); k++)//aabb
            {
                if(j!=k)
                {
                    if(     (aabb[i][j][0] > aabb[i][k][0]) && // aabb[i][j] is completely inside aabb[i][k]
                            (aabb[i][j][2] < aabb[i][k][2]) &&
                            (aabb[i][j][1] > aabb[i][k][1]) &&
                            (aabb[i][j][3] < aabb[i][k][3]) )
                    {
                        is_boundary= !is_boundary;
                    }
                }
            }
            if(is_boundary)
            {
                boxes.push_back(aabb[i][j]);
                boxes[boxes.size()-1].push_back(i);//adding the level information
                trim_curves.resize(trim_curves.size()+1);
                trim_curves[trim_curves.size()-1].push_back(res[i][j]);//trim_curves
            }
        }
    }
    //create the matrices b1,b2 and vector level
    b1.resize(boxes.size(),d);
    b2.resize(boxes.size(),d);
    level.resize(boxes.size());
    for(size_t i = 0; i < boxes.size(); i++){
        for(unsigned j = 0; j < d; j++){
            //convert from param spece to index space in highest level
            // (!) next lines: Conversion form double to int, possible loss of data
            b1(i,j) = boxes[i][j];
            b2(i,j) = boxes[i][j+d];
        }
        level[i] = boxes[i][2*d];
    }
    //use getBsplinePatches
    int nboxes = level.size();
    //------------------------------------------------------------------------------------------------------------------------------
    // iteration on the boxes to call getBsplinePatchGlobal()
    //------------------------------------------------------------------------------------------------------------------------------

    nvertices.resize(nboxes,this->dim());

    for (int i = 0; i < nboxes; i++)
    {
        gsMatrix<T> new_cp;
        gsVector<index_t> p1, p2;
        gsKnotVector<T> cku, ckv;

        p1 = b1.row(i).transpose();
        p2 = b2.row(i).transpose();

        this->getBsplinePatchGlobal(p1, p2, level[i], geom_coef, new_cp, cku, ckv);

        if (i == 0)
        {
            cp = new_cp;
        }
        else
        {
            gsMatrix<T> bigger;
            int cprows = cp.rows();
            /*cp.conservativeResize( cp.rows() + temp_cp.rows(), Eigen::NoChange );
              for( int j=cprows; j < cp.rows(); j++ )
              cp.row(j) = temp2.row(j-cprows);*/
            bigger.resize(cp.rows()+new_cp.rows(), cp.cols());
            for(int j=0; j< bigger.rows(); j++)
            {
                if(j<cprows)
                {
                    bigger.row(j) = cp.row(j);
                }
                else
                {
                    bigger.row(j) = new_cp.row(j-cprows);
                }
            }
            cp=bigger;
        }

        nvertices(i,0) = cku.size()-cku.degree()-1;
        nvertices(i,1) = ckv.size()-ckv.degree()-1;

    }
    // identify holes
    for(size_t l = 0; l < aabb.size();l++) //level
    {
        for(size_t i = 0; i < aabb[l].size(); i++) //bb
        {
            int closest_box = -1;
            for(size_t j = 0; j < boxes.size(); j++)
            {
                if(l == static_cast<size_t>(boxes[j][4]) )
                {
                    if(     (aabb[l][i][0] > boxes[j][0]) &&
                            (aabb[l][i][1] > boxes[j][1]) &&
                            (aabb[l][i][2] < boxes[j][2]) &&
                            (aabb[l][i][3] < boxes[j][3]))
                    {
                        if(closest_box!=-1)
                        {
                            //test with previous closest_box
                            if(!(
                                   (boxes[closest_box][0] > boxes[j][0]) &&
                                   (boxes[closest_box][1] > boxes[j][1]) &&
                                   (boxes[closest_box][2] < boxes[j][2]) &&
                                   (boxes[closest_box][3] < boxes[j][3])))
                            {
                                closest_box = j;
                            }
                        }
                        else
                        {
                            closest_box = j;
                        }

                    }
                    else if(    (aabb[l][i][0] == boxes[j][0]) &&
                                (aabb[l][i][1] == boxes[j][1]) &&
                                (aabb[l][i][2] == boxes[j][2]) &&
                                (aabb[l][i][3] == boxes[j][3]))
                    {
                        //boxes[j] == aabb[l][i]
                        closest_box = -1;
                        break;
                    }
                }

            }
            //prirad s res do trim_curves
            if(closest_box>-1)
            {
                trim_curves[closest_box].push_back(res[l][i]);
            }
        }
    }
}


//return data for trimming in parasolid
template<short_t d, class T>
gsMultiPatch<T> gsTHBSplineBasis<d,T>::getBsplinePatchesToMultiPatch_trimming(
    const gsMatrix<T>& geom_coef,
    std::vector<std::vector<std::vector< std::vector<T> > > >& trim_curves) const
{
    gsMatrix<index_t> b1;
    gsMatrix<index_t> b2;
    gsVector<index_t> level;
    gsMultiPatch<T> result;
    //identify the outer polylines- conected components
    int first_level = 0;
    for(unsigned int i = 0; i < this->m_xmatrix.size(); i++)
    {
        if(this->m_xmatrix[i].size()>0)
        {
            first_level = i;
            break;
        }
    }
    gsDebug<<"new min level"<<"\n";
    std::vector< std::vector< std::vector< std::vector< T > > > > res; //things to assign to trim_curves
    std::vector< std::vector< std::vector<index_t > > > aabb;//axis aligned bounding box
    std::vector< std::vector<index_t > > boxes;
    aabb = this->domainBoundariesParams(res);
    for(size_t i = 0; i < aabb.size(); i++)//level
    {
        //compare every aabb with the others
        for(size_t j = 0; j < aabb[i].size(); j++)//aabb
        {
            bool is_boundary = true;
            for(size_t k = 0; k < aabb[i].size(); k++)//aabb
            {
                if(j!=k)
                {
                    if(     (aabb[i][j][0] > aabb[i][k][0]) &&
                            (aabb[i][j][2] < aabb[i][k][2]) &&
                            (aabb[i][j][1] > aabb[i][k][1]) &&
                            (aabb[i][j][3] < aabb[i][k][3]) )
                    {
                        is_boundary= !is_boundary;
                    }
                }
            }
            if(is_boundary)
            {
                boxes.push_back(aabb[i][j]);
                boxes[boxes.size()-1].push_back(i+first_level);//adding the level information
                trim_curves.resize(trim_curves.size()+1);
                trim_curves[trim_curves.size()-1].push_back(res[i][j]);//trim_curves
            }
        }
    }
    //create the matrices b1,b2 and vector level
    b1.resize(boxes.size(),d);
    b2.resize(boxes.size(),d);
    level.resize(boxes.size());
    for(size_t i = 0; i < boxes.size(); i++)
    {
        for(unsigned j = 0; j < d; j++)
        {
            //convert from param spece to index space in highest level
            // (!) next lines: Conversion form double to int, possible loss of data
            b1(i,j) = boxes[i][j];
            b2(i,j) = boxes[i][j+d];
        }
        level[i] = boxes[i][2*d];
    }
    //use getBsplinePatches

    int nboxes = level.size();
    //------------------------------------------------------------------------------------------------------------------------------
    // iteration on the boxes to call getBsplinePatchGlobal()
    //------------------------------------------------------------------------------------------------------------------------------
    gsVector<index_t> p1, p2;
    p1.resize(this->dim());
    p2.resize(this->dim());
    gsMatrix<T> temp1;
    gsKnotVector<T> cku, ckv;


    for (int i = 0; i < nboxes; i++)
    {
        p1 = b1.row(i).transpose();
        p2 = b2.row(i).transpose();

        this->getBsplinePatchGlobal(p1, p2, level[i], geom_coef, temp1, cku, ckv);
        gsTensorBSplineBasis<2, T> tbasis(cku, ckv);
        result.addPatch(typename gsTensorBSpline<2, T>::uPtr(new gsTensorBSpline<2, T>(tbasis, give(temp1))));

    }
    // identify holes
    for(size_t l = 0; l < aabb.size();l++) //level
    {
        for(size_t i = 0; i < aabb[l].size(); i++) //bb
        {
            int closest_box = -1;
            for(size_t j = 0; j < boxes.size(); j++)
            {
                if(l == static_cast<size_t>(boxes[j][4]) )
                {
                    if(     (aabb[l][i][0] > boxes[j][0]) &&
                            (aabb[l][i][1] > boxes[j][1]) &&
                            (aabb[l][i][2] < boxes[j][2]) &&
                            (aabb[l][i][3] < boxes[j][3]))
                    {
                        if(closest_box!=-1){
                            //test with previous closest_box
                            if(!(
                                   (boxes[closest_box][0] > boxes[j][0]) &&
                                   (boxes[closest_box][1] > boxes[j][1]) &&
                                   (boxes[closest_box][2] < boxes[j][2]) &&
                                   (boxes[closest_box][3] < boxes[j][3])))
                            {
                                closest_box = j;
                            }
                        }
                        else
                        {
                            closest_box = j;
                        }

                    }
                    else if(    (aabb[l][i][0] == boxes[j][0]) &&
                                (aabb[l][i][1] == boxes[j][1]) &&
                                (aabb[l][i][2] == boxes[j][2]) &&
                                (aabb[l][i][3] == boxes[j][3]))
                    {
                        //boxes[j] == aabb[l][i]
                        closest_box = -1;
                        break;
                    }
                }

            }
            //prirad s res do trim_curves
            if(closest_box>-1)
            {
                trim_curves[closest_box].push_back(res[l][i]);
            }
        }
    }
    return result;
}


/*
  Private functions
*/


template<short_t d, class T>
void gsTHBSplineBasis<d,T>::globalRefinement(const gsMatrix<T> & thbCoefs,
                                             int level, gsMatrix<T> & lvlCoefs) const
{
    const index_t n = thbCoefs.cols();

    // Initialize level 0 coefficients
    lvlCoefs.setZero(m_bases[0]->size(), n);
    for(cmatIterator it = m_xmatrix[0].begin(); it != m_xmatrix[0].end(); ++it)
    {
        const int hIndex = m_xmatrix_offset[0] + (it - m_xmatrix[0].begin());
        lvlCoefs.row(*it) = thbCoefs.row(hIndex);
    }

    gsKnotVector<T> k1, k2;// fixme: boehm refine needs non-const kv
    std::vector<T> knots_x, knots_y;

    for(int l = 1; l <=level; l++)
    {
        k1 = m_bases[l-1]->knots(0);
        k2 = m_bases[l-1]->knots(1);

        // global dyadic refinement with respect to previous level
        k1.getUniformRefinementKnots(1,knots_x);
        k2.getUniformRefinementKnots(1,knots_y);

        // refine direction 0
        lvlCoefs.resize(m_bases[l-1]->size(0), n * m_bases[l-1]->size(1));
        gsBoehmRefine(k1, lvlCoefs, m_deg[0], knots_x.begin(), knots_x.end(), false);
        
        // refine direction 1
        lvlCoefs.blockTransposeInPlace(m_bases[l-1]->size(1));
        gsBoehmRefine(k2, lvlCoefs, m_deg[1], knots_y.begin(), knots_y.end(), false);
        lvlCoefs.blockTransposeInPlace(m_bases[l]->size(0));
        lvlCoefs.resize(m_bases[l]->size(), n); //lvlCoefs: control points at level \a l

        // overwrite with the THB coefficients of level \a l
        for(cmatIterator it = m_xmatrix[l].begin(); it != m_xmatrix[l].end(); ++it)
        {
            const int hIndex = m_xmatrix_offset[l] + (it - m_xmatrix[l].begin());
            lvlCoefs.row(*it) = thbCoefs.row(hIndex);
        }
    }
}


template<short_t d, class T>
void gsTHBSplineBasis<d,T>::evalSingle_into(index_t i,
                                            const gsMatrix<T>& u,
                                            gsMatrix<T>& result) const
{    
    if (this->m_is_truncated[i] == -1)  // basis function not truncated
    {
        unsigned level = this->levelOf(i);
        unsigned tensor_index = flatTensorIndexOf(i, level);           
        this->m_bases[level]->evalSingle_into(tensor_index, u, result);
    }
    else
    {
        
        unsigned level = this->m_is_truncated[i];
        
        const gsSparseVector<T>& coefs = getCoefs(i);
        
        const gsTensorBSplineBasis<d, T>& base =
            *this->m_bases[level];
        
        gsTensorDeboor<d, T, gsKnotVector<T>, gsSparseVector<T> >
            (u, base, coefs, result);
    }
}

template<short_t d, class T>
void gsTHBSplineBasis<d,T>::deriv2Single_into(index_t i,
                                              const gsMatrix<T>& u,
                                              gsMatrix<T>& result) const
{
    
    if (this->m_is_truncated[i] == -1) // basis function not truncated
    {
        const unsigned level = this->levelOf(i);
        const unsigned fl_tensor_index = flatTensorIndexOf(i, level);
        this->m_bases[level]->deriv2Single_into(fl_tensor_index, u, result);
    }
    else
    {
        const unsigned level = this->m_is_truncated[i];
        const gsSparseVector<T>& coefs = this->getCoefs(i);
        const gsTensorBSplineBasis<d, T> & base =
            *this->m_bases[level];
        
        gsTensorDeriv2_into<d, T, gsKnotVector<T>,
                            gsSparseVector<T> >(u, base, coefs, result);
    }
}

template<short_t d, class T>
void gsTHBSplineBasis<d,T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    gsMatrix<index_t> indices;
    gsMatrix<T> res(1, 1);
    this->active_into(u, indices);

    result.setZero(indices.rows(), u.cols());

    for (index_t i = 0; i < indices.cols(); i++)
    {
        for (index_t j = 0; j < indices.rows(); j++)
        {
            const index_t index = indices(j, i);
            if (j != 0 && index == 0)
                break;
            
            evalSingle_into(index, u.col(i), res);
            result(j, i) = res.value();
        }
    }
}


template<short_t d, class T>
void gsTHBSplineBasis<d,T>::deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result)const
{
    gsMatrix<index_t> indices;
    this->active_into(u, indices);

    static const short_t numDers = (d * (d + 1)) / 2;
    gsMatrix<T> res(numDers, 1); // result of deriv2Single_into

    result.setZero(indices.rows() * numDers, u.cols());

    for (int i = 0; i < indices.cols(); i++)
    {
        for (int j = 0; j < indices.rows(); j++)
        {
            const index_t index = indices(j, i);
            if (j != 0 && index == 0)
                break;

            this->deriv2Single_into(index, u.col(i), res);

            result.template block<numDers, 1>(j * numDers, i) = res;
        }
    }
}


template<short_t d, class T>
void gsTHBSplineBasis<d,T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{

    gsMatrix<index_t> indices;
    this->active_into(u, indices);
    gsMatrix<T> res(d, 1);

    result.setZero(indices.rows() * d, u.cols());

    for (int i = 0; i < indices.cols(); i++)
    {
        for (int j = 0; j < indices.rows(); j++)
        {

            const unsigned index = indices(j, i);
            if (j != 0 && index == 0)
                break;

            this->derivSingle_into(index, u.col(i), res);
            result.template block<d,1>(j * d, i) = res;
        }
    }
}


template<short_t d, class T>
void gsTHBSplineBasis<d,T>::derivSingle_into(index_t i,
                                             const gsMatrix<T> & u,
                                             gsMatrix<T>& result) const
{

    if (this->m_is_truncated[i] == -1) // basis function not truncated
    {
        unsigned level = this->levelOf(i);
        unsigned fl_tensor_index = flatTensorIndexOf(i, level);
        this->m_bases[level]->derivSingle_into(fl_tensor_index, u, result);
    }
    else
    {
        unsigned level = this->m_is_truncated[i];
        const gsSparseVector<T>& coefs = this->getCoefs(i);
        const gsTensorBSplineBasis<d,T>& base =
            *this->m_bases[level];
        gsTensorDeriv_into<d, T, gsKnotVector<T>,
                           gsSparseVector<T> >(u, base, coefs, result);
    }
}


template<short_t d, class T>
void gsTHBSplineBasis<d, T>::decomposeDomain(
    typename gsTHBSplineBasis<d, T>::AxisAlignedBoundingBox& boundaryAABB,
    typename gsTHBSplineBasis<d, T>::TrimmingCurves& trimCurves) const
{
    Polylines polylines;
    AxisAlignedBoundingBox aabb;
    
    aabb = this->domainBoundariesParams(polylines);
    breakCycles(aabb, polylines);

    int numBoundaryBoxes = 0;
    for (unsigned level = 0; level != aabb.size(); level++)
    {
        boundaryAABB.push_back(std::vector< std::vector<index_t> >());
        trimCurves.push_back(std::vector< std::vector< std::vector< std::vector<T> > > >());

        //compare every aabb with the others
        for (unsigned boxI = 0; boxI != aabb[level].size(); boxI++)
        {
            bool isBoundaryBox = true; 
            for (unsigned boxJ = 0; boxJ != aabb[level].size(); boxJ++)
            {
                if (boxI != boxJ)
                {
                    if (isFirstBoxCompletelyInsideSecond(aabb[level][boxI], aabb[level][boxJ]))
                    {
                        isBoundaryBox = !isBoundaryBox;
                    }
                }
            }
	    
            if (isBoundaryBox)
            {
                numBoundaryBoxes++;
                boundaryAABB[level].push_back(aabb[level][boxI]);
		
                // make new componenet
                trimCurves[level].push_back(std::vector< std::vector< std::vector<T> > >());
                trimCurves[level][trimCurves[level].size() - 1].push_back(polylines[level][boxI]);
            }
        }	
    }
    
    for (unsigned level = 0; level != aabb.size(); level++)
    {
        for (unsigned box = 0; box != aabb[level].size(); box++)
        {
            int closestBox = -1;
            for (unsigned boundBox = 0; 
                 boundBox != boundaryAABB[level].size(); 
                 boundBox++)
            {
                if (isFirstBoxCompletelyInsideSecond(aabb[level][box], boundaryAABB[level][boundBox]))
                {
                    if (closestBox == -1 ||  
                        !isFirstBoxCompletelyInsideSecond(boundaryAABB[level][closestBox],
                                                          boundaryAABB[level][boundBox]))
                    {
                        closestBox = boundBox;
                    }
                }
                else if (areBoxesTheSame(aabb[level][box], boundaryAABB[level][boundBox]))
                {
                    closestBox = -1;
                    break;
                }
            }
	    
            if (-1 < closestBox)
            {
                trimCurves[level][closestBox].push_back(polylines[level][box]);
            }
        }
    }
}


template<short_t d, class T>
gsTensorBSpline<d, T> 
gsTHBSplineBasis<d,T>::getBSplinePatch(const std::vector<index_t>& boundingBox,
                                       const unsigned level,
                                       const gsMatrix<T>& geomCoefs) const
{
    gsVector<index_t, d> low, upp;
    for (unsigned dim = 0; dim != d; dim++)
    {
        low(dim) = boundingBox[dim];
        upp(dim) = boundingBox[d + dim];
    }
    this->m_tree.computeLevelIndex(low, level, low);
    this->m_tree.computeLevelIndex(upp, level, upp);
    
    const gsKnotVector<T> & knots0 = m_bases[level]->knots(0);
    const gsKnotVector<T> & knots1 = m_bases[level]->knots(1);

    const int lowIndex0 = knots0.lastKnotIndex (low(0)) - m_deg[0];
    const int uppIndex0 = knots0.firstKnotIndex(upp(0)) - 1;
    const int lowIndex1 = knots1.lastKnotIndex (low(1)) - m_deg[1];
    const int uppIndex1 = knots1.firstKnotIndex(upp(1)) - 1;

    const int numDirection0 = uppIndex0 - lowIndex0 + 1;
    const int numDirection1 = uppIndex1 - lowIndex1 + 1;
    const int numNewCoefs = numDirection0 * numDirection1;
    gsMatrix<T> newCoefs(numNewCoefs, geomCoefs.cols());
    const index_t sz0 = m_bases[level]->size(0);

    gsMatrix<T> coefs;
    globalRefinement(geomCoefs, level, coefs);
    index_t cc = 0;
    for (int j = lowIndex1; j <= uppIndex1; j++)
        for (int i = lowIndex0; i <= uppIndex0; i++)
            newCoefs.row(cc++) =  coefs.row(j*sz0+i);

    std::vector<gsKnotVector<T> > kv(2);

    kv[0] = gsKnotVector<T>(m_deg[0], knots0.begin() + lowIndex0, 
                            knots0.begin() + uppIndex0 + m_deg[0] + 2);

    kv[1] = gsKnotVector<T>(m_deg[1], knots1.begin() + lowIndex1, 
                            knots1.begin() + uppIndex1 + m_deg[1] + 2);

    tensorBasis basis(kv);

    return gsTensorBSpline<d, T> (basis, newCoefs);
}

template<short_t d, class T>
void gsTHBSplineBasis<d, T>::breakCycles(
    typename gsTHBSplineBasis<d, T>::AxisAlignedBoundingBox& aabb,
    typename gsTHBSplineBasis<d, T>::Polylines& polylines) const
{
    for (size_t level = 0; level != polylines.size(); level++)
    {
        for (size_t line = 0; line != polylines[level].size(); line++)
        {
            std::pair<T, T> pt; // point
            index_t segment = identifyCycle(polylines[level][line], pt);
	    
            if (-1 < segment)
            {
                std::vector< std::vector<T> > part1, part2;
                breakPolylineIntoTwoParts(polylines[level][line], segment, pt,
                                          part1, part2);
		
                polylines[level][line] = part1;
                polylines[level].push_back(part2);
		
                std::vector<index_t> aabb1, aabb2;
                findNewAABB(part1, aabb1);
                findNewAABB(part2, aabb2);
		
                aabb[level][line] = aabb1;
                aabb[level].push_back(aabb2);
		
                // very important, this will check current line again if it has more cycles
                line--; 
            }
        }
    }
}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// utility funcitions for breakCycles
// ......................................................................

// utility funcition for breakCycles
template<short_t d, class T>
index_t gsTHBSplineBasis<d, T>::identifyCycle(const std::vector< std::vector< T> >& line,
                                              std::pair<T, T>& pt) const
{
    std::map< std::pair<T, T>, index_t > times;
    std::map< std::pair<T, T>, index_t > index;
    
    for (size_t seg = 0; seg != line.size(); seg++)
    {
        const size_t seg1 = (seg + 1) % line.size();
	
        std::pair<T, T> currentPt( line[seg][0], line[seg][1] );
        if (!((currentPt.first == line[seg1][0] && currentPt.second == line[seg1][1]) ||
              (currentPt.first == line[seg1][2] && currentPt.second == line[seg1][3])))
        {
            currentPt.first = line[seg][2];
            currentPt.second = line[seg][3];
        }
	
        size_t count = times.count(currentPt);
        if (0 < count)
        {
            times[currentPt] += 1;
        }
        else
        {
            times[currentPt] = 1;
            index[currentPt] = seg1;
        }
    }

    typedef typename std::map< std::pair<T, T>, index_t>::iterator iterator;
    for (iterator it = times.begin(); it != times.end(); it++)
    {
        if (it->second == 2)
        {
            pt = it->first;
            return index[it->first];
        }
        else
        {
            GISMO_ENSURE(2 > it->second, "Internal error. Check the polylines from the domainBoundariesParam." );
        }
    }
    return -1;
}

// utility funcition for breakCycles
template<short_t d, class T>
void gsTHBSplineBasis<d, T>::breakPolylineIntoTwoParts(
    const std::vector< std::vector< T> >& line,
    const index_t segment, 
    const std::pair<T, T>& meetingPt,
    std::vector< std::vector<T> >& part1,
    std::vector< std::vector<T> >& part2) const
{
    bool p1 = false; // inside part 1
    bool p2 = false; // inside part 2
    
    index_t length = static_cast<index_t> (line.size());
    for (index_t i = 0; i != length; i++)
    {
        const index_t seg = (i + segment) % length;
	
        if (!p1 && !p2) // start
        {
            p1 = true;
            part1.push_back(line[seg]);
        }
        else // not start
        {
            // we hit the meeting point again
            if ((meetingPt.first == line[seg][0] && meetingPt.second == line[seg][1]) ||
                (meetingPt.first == line[seg][2] && meetingPt.second == line[seg][3]))
            {
                if (p1) // end of part 1
                {
                    part1.push_back(line[seg]);
                    p1 = false;
                    p2 = true;
                }
                else if (p2) // start or finish
                {
                    part2.push_back(line[seg]);
                }
            }
            else
            {
                if (p1)
                {
                    part1.push_back(line[seg]);
                }
                else if (p2)
                {
                    part2.push_back(line[seg]);
                }
            }
        }
    }
}

// utility funcition for breakCycles
template<short_t d, class T>
void gsTHBSplineBasis<d, T>::findNewAABB(const std::vector< std::vector<T> >& polyline,
                                         std::vector<index_t>& aabb) const
{
    T minX = polyline[0][0];
    T minY = polyline[0][1];
    T maxX = polyline[0][2];
    T maxY = polyline[0][3];
    

    for (size_t seg = 0; seg != polyline.size(); seg++)
    {
        if (polyline[seg][0] < minX)
        {
            minX = polyline[seg][0];
        }
        if (polyline[seg][1] < minY)
        {
            minY = polyline[seg][1];
        }
        if (maxX < polyline[seg][2])
        {
            maxX = polyline[seg][2];
        }
        if (maxY < polyline[seg][3])
        {
            maxY = polyline[seg][3];
        }
    }
    
    unsigned maxLevel = this->maxLevel();
    const gsKnotVector<T>& kv0 = this->m_bases[maxLevel]->knots(0);
    const gsKnotVector<T>& kv1 = this->m_bases[maxLevel]->knots(1);
    
    aabb.resize(4);
    for (unsigned i = 0; i != kv0.uSize(); i++)
    {
        if (kv0.uValue(i) <= minX)
        {
            aabb[0] = i;
        }
        if (maxX <= kv0.uValue(i))
        {
            aabb[2] = i;
            break;
        }
    }
    
    for (unsigned i = 0; i != kv1.uSize(); i++)
    {
        if (kv1.uValue(i) <= minY)
        {
            aabb[1] = i;
        }
        if (maxY <= kv1.uValue(i))
        {
            aabb[3] = i;
            break;
        }
    }
}



// --------------------------------------------------------------------------------
// Code for hierarchical coarsening
// --------------------------------------------------------------------------------

template<short_t d, class T>
void gsTHBSplineBasis<d,T>::transferbyLvl (std::vector<gsSparseMatrix<T> >& result)
{
    result.clear();
    gsVector<index_t> level;
    gsMatrix<index_t> b1, b2;//boxes in highes level numbering
    this->m_tree.getBoxesInLevelIndex(b1,b2,level);//return boxes in level indices
    tensorBasis T_0_copy = this->tensorLevel(0);
    std::vector< gsSparseMatrix<T,RowMajor> > transfer;
    transfer.resize(this->maxLevel() );
    std::vector<std::vector<T> > knots(d);

    for(unsigned i = 0; i < this->maxLevel(); ++i)
    {
        //T_0_copy.uniformRefine_withTransfer(transfer[i], 1);
        for(short_t dim = 0; dim < d; dim++)
        {
            const gsKnotVector<T> & ckv = m_bases[i]->knots(dim);
            const gsKnotVector<T> & fkv = m_bases[i + 1]->knots(dim);
            ckv.symDifference(fkv, knots[dim]);

            //gsDebug << "level: " << i << "\n"
            //        << "direction: " << dim << "\n"
            //        << "dirknots:\n"<<gsAsMatrix<T>(dirKnots)<<std::endl;
        }
        T_0_copy.refine_withTransfer(transfer[i], knots);

        // Must use refine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, const std::vector<T>& knots)
        // must correctly find the knots to insert
        // T_0_copy.uniformRefine_withTransfer(transfer[i], 1);
    }
    std::vector< gsSparseMatrix<T,RowMajor> > temp_transf;
    for(unsigned j = 0; j < this->maxLevel(); ++j)
    {
        std::vector<CMatrix> x_mat_old_0, x_matrix_lvl;
        this->setActiveToLvl(j,x_mat_old_0);
        this->setActiveToLvl(j+1,x_matrix_lvl);
        temp_transf.push_back(transfer[j]);

        gsSparseMatrix<T> crs = this->coarsening_direct(x_mat_old_0, x_matrix_lvl, temp_transf);
        result.push_back(crs);
    }
}

//todo remove
template<short_t d, class T>
gsSparseMatrix<T> gsTHBSplineBasis<d,T>::coarsening( const std::vector<gsSortedVector<index_t> >& old, const std::vector<gsSortedVector<index_t> >& n, const gsSparseMatrix<T,RowMajor> & transfer) const
{
    int size1= 0, size2 = 0;
    int glob_numb = 0;//continous numbering of hierarchical basis
    for(unsigned int i =0; i< old.size();i++)
    {//count the number of basis functions in old basis
        size1 += old[i].size();
    }
    for(unsigned int i =0; i< n.size();i++)
    {//count the number of basis functions in new basis
        size2 += n[i].size();
    }
    gsSparseMatrix<T> result(size2,size1);

    gsMatrix<T> transferDense = transfer;

    for (unsigned int i = 0; i < old.size(); i++)//iteration through the levels of the old basis
    {
        // find starting index of level i in new basis
        int start_lv_i = 0;
        for(unsigned int l =0; l < i; l++)
        {
            start_lv_i += n[l].size();
        }

        for (unsigned int j = 0; j < old[i].size();j++)//iteration through the basis functions in the given level
        {
            start_lv_i = 0;
            for(unsigned int l =0; l < i; l++)
            {
                start_lv_i += n[l].size();
            }
            const unsigned old_ij = old[i][j];  // tensor product index

            if( n[i].bContains(old_ij) )//it he basis function was not refined
            {
                result(start_lv_i + std::distance(n[i].begin(), n[i].find_it_or_fail(old_ij) ),glob_numb ) = 1;//settign the coefficient of the not refined basis function to 1
                for(int k = 0; k < transferDense.rows(); k++)//basis function was refined->looking for the coeficinets from the global transfer matrix
                {
                    //TODO test if the basis function is in the basis
                    if(transferDense(k, old_ij) != 0)//if the coefficient is non zero we find the coresponding function in n
                    {
                        if(n[i+1].bContains(k))
                        {
                            const int pos = start_lv_i + n[i].size() + std::distance(n[i+1].begin(), n[i+1].find_it_or_fail(k));
                            result(pos,glob_numb) = transferDense(k, old_ij);
                        }
                    }
                }
            }
            else
            {
                for(int k = 0; k < transferDense.rows(); k++)//basis function was refined->looking for the coeficinets from the global transfer matrix
                {
                    if(transferDense(k, old_ij) != 0)//if the coefficient is non zero we find the coresponding function in n
                    {
                        const int pos = start_lv_i + n[i].size() + std::distance(n[i+1].begin(), n[i+1].find_it_or_fail(k));
                        result(pos,glob_numb) = transferDense(k, old_ij);
                    }
                }
            }
            glob_numb++;
        }
    }
    return result;
}
template<short_t d, class T>
gsSparseMatrix<T> gsTHBSplineBasis<d,T>::coarsening_direct2( const std::vector<gsSortedVector<index_t> >& old,
                                                       const std::vector<gsSortedVector<index_t> >& n,
                                                       const std::vector<gsSparseMatrix<T,RowMajor> >& transfer) const
{
    int size1= 0;int size2 = 0;
    int glob_numb = 0;//continous numbering of hierarchical basis
    for(unsigned int i =0; i< old.size();i++){//count the number of basis functions in old basis
        size1 += old[i].size();
    }
    for(unsigned int i =0; i< n.size();i++){//count the number of basis functions in new basis
        size2 += n[i].size();
    }
    gsSparseMatrix<T> result(size2,size1);


    //std::vector<gsMatrix<T> > transferDense;// = transfer;
    //transferDense.resize(transfer.size());
    //for (unsigned int i = 0; i < transfer.size();i++){
    //    transferDense[i] = transfer[i];
    //}
    std::vector<gsSparseMatrix<T,ColMajor> > temptransfer;// = transfer;
    temptransfer.resize(transfer.size());
    for (unsigned int i = 0; i < transfer.size();i++){
        temptransfer[i] = transfer[i];
        //gsDebug<<"transfer"<<i<<"\n"<<transfer[i]<<std::endl;
    }

    for (unsigned int i = 0; i < old.size(); i++)//iteration through the levels of the old basis
    {
        // find starting index of level i in new basis
        int start_lv_i = 0;
        for(unsigned int l =0; l < i; l++)
        {
            start_lv_i += n[l].size();
        }

        //gsDebug<<old[i].size()<<std::endl;
        for (unsigned int j = 0; j < old[i].size();j++)//iteration through the basis functions in the given level
        {
            //gsDebug<<"j = "<<j<<std::endl;
            start_lv_i = 0;
            for(unsigned int l =0; l < i; l++)
            {
                start_lv_i += n[l].size();
            }
            const unsigned old_ij = old[i][j];  // tensor product index
            gsMatrix<index_t, d, 2> supp(d, 2);
            this->m_bases[i]->elementSupport_into(old_ij, supp);//this->support(start_lv_i+old_ij);
            //gsDebug<<"supp "<< supp<<std::endl;
            //unsigned max_lvl =
            //    math::min<index_t>( this->m_tree.query4(supp.col(0),supp.col(1), i), transfer.size() ) ;
            gsSparseVector<T,RowMajor> t(this->m_bases[i]->size());
            t.setZero();
            t[old_ij] = 1;
            for(unsigned int k = i; k < n.size();k++){
                //gsDebug<<"i: "<<i<<" j: "<<j<<" k: "<<k<<std::endl;
                if(k > i)
                {
                    //compare with old matrix
                    for(int l = 0 ; l < t.size();l++){
                        if(t[l]!=0)
                            if( (k) < (old.size()))
                                if(old[k].bContains(l)){
                                    t[l]=0;
                                }

                    }
                }
                //gsDebug<<"ksize:"<<n[k].size()<<std::endl;
                if(k!=0)
                {
                    start_lv_i = 0;
                    for(unsigned int l =0; l < k-1; l++)
                    {
                        start_lv_i += n[l].size();
                    }
                }
                //for all non zero in t comapre with new
                for(int l = 0 ; l < t.size();l++)
                {
                    //gsDebug<<"i: "<<i<<" j: "<<j<<" l: "<<l<<"nsize"<<n.size()<<std::endl;
                    if(t[l]!=0)
                        if(n[k].bContains(l))
                        {
                            //gsDebug<<"j: "<<j<<" "<<"oldij"<<old_ij<<" "<<"l:"<<l<<"    ";
                            //gsDebug<<"k "<<k<<" j: "<<j<<" "<<"globnumb "<<glob_numb<<" oldij "<<old_ij<<" "<<"l:"<<l<<" "<<t[l]<<"  ";
                            int p = 0;
                            if(k!=0){
                                p = n[k-1].size();
                            }
                            const int pos = start_lv_i + p + std::distance(n[k].begin(), n[k].find_it_or_fail(l));
                            //const int pos = start_lv_i + n[k].size() + std::distance(n[k+1].begin(), n[k+1].find_it_or_fail(l));
                            //gsDebug<<"pos:"<<pos<<std::endl;
                            //double ppp =  transferDense[coeffs[0].lvl](k, coeffs[0].pos);
                            result(pos,glob_numb) = t[l];//transferDense[coeffs[0].lvl](k, coeffs[0].pos);
                            //t[l] = 0;
                            //const int pos = start_lv_i + n[coeffs[0].lvl].size() + std::distance(n[coeffs[0].lvl+1].begin(), n[coeffs[0].lvl+1].find_it_or_fail(k.row()));
                            //double ppp =  transferDense[coeffs[0].lvl](k, coeffs[0].pos);
                            //result(pos,glob_numb) += coeffs[0].coef * temptransfer[coeffs[0].lvl](k.row(), coeffs[0].pos);//transferDense[coeffs[0].lvl](k, coeffs[0].pos);
                        }
                }
                //gsDebug<<"b"<<std::endl;
                gsSparseVector<T,RowMajor> temp(this->m_bases[i]->size());
                //gsDebug<<"c"<<std::endl;
                temp.setZero();
                if(k<temptransfer.size())
                    temp = temptransfer[k] * t.transpose();
                t = temp;
                //refine
            }


//            if( n[i].bContains(old_ij) )//it he basis function was not refined
//            {
//                result(start_lv_i + std::distance(n[i].begin(), n[i].find_it_or_fail(old_ij) ),glob_numb ) = 1;//settign the coefficient of the not refined basis function to 1
//            }
//            else
//            {

//                gsMatrix<index_t, d, 2> supp(d, 2);
//                this->m_bases[i]->elementSupport_into(old_ij, supp);//this->support(start_lv_i+old_ij);
//                //gsDebug<<<"supp "<< supp<<std::endl;
//                unsigned max_lvl =
//                    math::min<index_t>( this->m_tree.query4(supp.col(0),supp.col(1), i), transfer.size() ) ;
//                gsSparseVector<T,RowMajor> t(this->m_bases[i]->size());
//                t.setZero();
//                t[old_ij] = 1;
//                //gsDebug<<"j = "<<j<<std::endl;
//                //gsDebug<<"nsize"<<n.size()<<std::endl;
//                for(unsigned int k = i+1; k < n.size();k++){
//                    start_lv_i = 0;
//                    for(unsigned int l =0; l < k-1; l++)
//                    {
//                        start_lv_i += n[l].size();
//                    }
// //                    gsDebug<<"k "<<k<<std::endl;
// //                    gsDebug<<"nk"<<std::endl;
// //                    for(int a = 0; a < n[k].size();a++){
// //                        gsDebug<<n[k][a]<<" ";
// //                    }
//                    //gsDebug<<std::endl;
//                    gsSparseVector<T,RowMajor> M;
//                    M.setZero();
//                    M = temptransfer[k-1] * t.transpose();

//                    //gsDebug<<"M\n"<<M<<std::endl;
//                    for(int l = 0 ; l < M.size();l++)
//                    //for(typename gsSparseVector<T,RowMajor>::InnerIterator l(M); l; ++l)//basis function was refined->looking for the coeficinets from the global transfer matrix
//                    {
//                        if(M[l]!=0)
//                            if(n[k].bContains(l))
//                            {
//                                //gsDebug<<"j: "<<j<<" "<<"oldij"<<old_ij<<" "<<"l:"<<l<<"    ";
//                                //gsDebug<<"k "<<k<<" j: "<<j<<" "<<"globnumb "<<glob_numb<<" "<<"l:"<<l<<" "<<M[l]<<"  ";
//                                const int pos = start_lv_i + n[k-1].size() + std::distance(n[k].begin(), n[k].find_it_or_fail(l));
//                                //gsDebug<<"pos:"<<pos<<std::endl;
//                                //double ppp =  transferDense[coeffs[0].lvl](k, coeffs[0].pos);
//                                result(pos,glob_numb) = M[l];//transferDense[coeffs[0].lvl](k, coeffs[0].pos);
//                                M[l] = 0;
//                                //const int pos = start_lv_i + n[coeffs[0].lvl].size() + std::distance(n[coeffs[0].lvl+1].begin(), n[coeffs[0].lvl+1].find_it_or_fail(k.row()));
//                                //double ppp =  transferDense[coeffs[0].lvl](k, coeffs[0].pos);
//                                //result(pos,glob_numb) += coeffs[0].coef * temptransfer[coeffs[0].lvl](k.row(), coeffs[0].pos);//transferDense[coeffs[0].lvl](k, coeffs[0].pos);
//                            }
//                    }
//                    t = M;
//                    //gsDebug<<"M after \n"<<M<<std::endl;
//                }
//            }
            glob_numb++;
        }

    }

    return result;
}

template<short_t d, class T>
gsSparseMatrix<T> gsTHBSplineBasis<d,T>::coarsening_direct( const std::vector<gsSortedVector<index_t> >& old,
                                                      const std::vector<gsSortedVector<index_t> >& n,
                                                      const std::vector<gsSparseMatrix<T,RowMajor> >& transfer) const
{
    GISMO_ASSERT(old.size() < n.size(), "old,n problem in coarsening.");

    int size1= 0;int size2 = 0;
    int glob_numb = 0;//continous numbering of hierarchical basis
    for(unsigned int i =0; i< old.size();i++)
    {//count the number of basis functions in old basis
        size1 += old[i].size();
    }
    for(unsigned int i =0; i< n.size();i++)
    {//count the number of basis functions in new basis
        size2 += n[i].size();
    }
    gsSparseMatrix<T> result(size2,size1);


//    std::vector<gsMatrix<T> > transferDense;// = transfer;
//    transferDense.resize(transfer.size());
//    for (unsigned int i = 0; i < transfer.size();i++)
//    {
//        transferDense[i] = transfer[i];
//    }
    std::vector<gsSparseMatrix<T,ColMajor> > temptransfer;// = transfer;
    temptransfer.resize(transfer.size());
    for (unsigned int i = 0; i < transfer.size();i++)
    {
        temptransfer[i] = transfer[i];
    }
    //gsDebug<<"temp:\n"<< temptransfer[0]<<std::endl;

    for (unsigned int i = 0; i < old.size(); i++)//iteration through the levels of the old basis
    {
        // find starting index of level i in new basis
        int start_lv_i = 0;
        for(unsigned int l =0; l < i; l++)
        {
            start_lv_i += n[l].size();
        }

        for (unsigned int j = 0; j < old[i].size();j++)//iteration through the basis functions in the given level
        {
            //gsDebug<<"j...."<<j<<endl;
            //gsDebug<<"i = "<< i<< " j= "<<j<<std::endl;
            start_lv_i = 0;
            for(unsigned int l =0; l < i; l++)
            {
                start_lv_i += n[l].size();
            }
            const unsigned old_ij = old[i][j];  // tensor product index
            //gsDebug<<"oldij = "<< old_ij<<std::endl;
            if( n[i].bContains(old_ij) )//it he basis function was not refined
            {
                result(start_lv_i + std::distance(n[i].begin(), n[i].find_it_or_fail(old_ij) ),glob_numb ) = 1;//settign the coefficient of the not refined basis function to 1
                std::vector<lvl_coef> coeffs;
                gsMatrix<index_t, d, 2> supp(d, 2);
                this->m_bases[i]->elementSupport_into(old_ij, supp);//this->support(start_lv_i+old_ij);
                unsigned max_lvl = math::min<index_t>( this->m_tree.query4(supp.col(0),supp.col(1), i), transfer.size()) ;//transfer.size();//
                //gsDebug<<"transfer size "<< transfer.size()<<" max lvl"<< this->m_tree.query4(supp.col(0),supp.col(1), this->levelOf(start_lv_i+j))<<"   lvl of"<< this->levelOf(start_lv_i+j)<<" support\n"<<supp<<std::endl;
                lvl_coef temp;
                temp.pos = old_ij;
                temp.coef = 1;
                temp.lvl = i;
                if(temp.lvl<n.size()-1)//in case of no refinement
                {
                    coeffs.push_back(temp);
                }
                for (size_t coeff_index = 0; coeff_index < coeffs.size(); ++coeff_index)
                {
                    const lvl_coef coeff = coeffs[coeff_index];

                    start_lv_i = 0;
                    for(unsigned int l =0; l < coeff.lvl; l++)
                    {
                        start_lv_i += n[l].size();
                    }

                    for(typename gsSparseMatrix<T,ColMajor>::InnerIterator k(temptransfer[coeff.lvl],coeff.pos); k; ++k)//basis function was refined->looking for the coeficinets from the global transfer matrix
                    {
//                        if(transferDense[coeff.lvl](k, coeff.pos) != 0)
//                        {
                        //gsDebug<<"first while "<< coeffs.size()<<std::endl;
                        bool p = true;
                        if( (coeff.lvl+1) < old.size())
                        {
                            if(old[coeff.lvl+1].bContains(k.row()))
                            {
                                p = false;
                            }
                        }
                        if(n[coeff.lvl+1].bContains(k.row()))
                        {
                            if(p)
                            {
                                const int pos = start_lv_i + n[coeff.lvl].size() + std::distance(n[coeff.lvl+1].begin(), n[coeff.lvl+1].find_it_or_fail(k.row()));

                                result(pos,glob_numb) += coeff.coef * temptransfer[coeff.lvl](k.row(), coeff.pos);//transferDense[coeff.lvl](k, coeff.pos);
                                if(coeff.lvl + 1 < max_lvl)//transfer.size()
                                {
                                    temp.pos = k.row();
                                    temp.coef = temptransfer[coeff.lvl](k.row(), coeff.pos) * coeff.coef;
                                    temp.lvl = coeff.lvl+1;
                                    coeffs.push_back(temp);
                                }
                            }
                        }else
                        {
                            //TODO test if level would not be higher than the q4 for this function- or max inserted level
                            if(p)
                            {
                                if( coeff.lvl + 1< max_lvl)
                                {
                                    temp.pos = k.row();
                                    temp.coef = temptransfer[coeff.lvl](k.row(), coeff.pos) * coeff.coef;
                                    temp.lvl = coeff.lvl+1;
                                    coeffs.push_back(temp);
                                }
                            }
                        }
                        //}
                    }

//                        for(unsigned int ii = 0; ii < coeffs.size();ii++){
//                            gsDebug<<"( "<< coeffs[ii].pos<<" , "<<coeffs[ii].coef<<" , "<<coeffs[ii].lvl<<" )"<<std::endl;
//                        }
//                        gsDebug<<"j ="<<j<<std::endl;
                }
            }
            else
            {
                gsMatrix<index_t, d, 2> supp(d, 2);
                this->m_bases[i]->elementSupport_into(old_ij, supp);//this->support(start_lv_i+old_ij);
                //gsDebug<<"supp "<< supp<<std::endl;
                unsigned max_lvl =
                    math::min<index_t>( this->m_tree.query4(supp.col(0),supp.col(1), i), transfer.size() ) ;
                std::vector<lvl_coef> coeffs;
                lvl_coef temp;
                temp.pos = old_ij;
                temp.coef = 1;
                temp.lvl = i;
                //gsDebug<<"temp.pos: "<<temp.pos<<" temp.coef: "<<temp.coef<<" temp.lvl "<< temp.lvl<<std::endl;
                coeffs.push_back(temp);
                for (size_t coeff_index = 0; coeff_index < coeffs.size(); ++coeff_index)
                {
                    const lvl_coef coeff = coeffs[coeff_index];

                    start_lv_i = 0;
                    for(unsigned int l =0; l < coeff.lvl; l++)
                    {
                        start_lv_i += n[l].size();
                    }

                    for(typename gsSparseMatrix<T,ColMajor>::InnerIterator k(temptransfer[coeff.lvl],coeff.pos); k; ++k)//basis function was refined->looking for the coeficinets from the global transfer matrix
                    {
                        //if(transferDense[coeff.lvl](k, coeff.pos) != 0)
                        //{
                        //gsDebug<<"second while "<< coeffs.size()<<std::endl;
                        bool p = true;
                        if( (coeff.lvl+1) < (old.size()))
                        {
                            if(old[coeff.lvl+1].bContains(k.row()))
                            {
                                p = false;
                            }
                        }
                        //gsDebug<<"k.row "<<k.row()<<std::endl;
                        if(n[coeff.lvl+1].bContains(k.row()))
                        {
                            if(p)
                            {
                                const int pos = start_lv_i + n[coeff.lvl].size() + std::distance(n[coeff.lvl+1].begin(), n[coeff.lvl+1].find_it_or_fail(k.row()));
                                //T ppp =  transferDense[coeff.lvl](k, coeff.pos);
                                //gsDebug<<"pos "<<pos<<" oldij "<<old_ij<<" "<< "coeflvl "<<coeff.lvl<<" ";
                                //gsDebug<<"inserted coef "<< coeff.coef<<"*"<< temptransfer[coeff.lvl](k.row(), coeff.pos)<<std::endl;
                                result(pos,glob_numb) += coeff.coef * temptransfer[coeff.lvl](k.row(), coeff.pos);//transferDense[coeff.lvl](k, coeff.pos);
                                if(coeff.lvl < max_lvl-1)
                                {
                                    temp.pos = k.row();
                                    temp.coef = temptransfer[coeff.lvl](k.row(), coeff.pos) * coeff.coef;
                                    temp.lvl = coeff.lvl+1;
                                    //gsDebug<<"temp.pos: "<<temp.pos<<" temp.coef: "<<temp.coef<<" temp.lvl "<< temp.lvl<<std::endl;
                                    coeffs.push_back(temp);
                                }
                            }
                        }
                        else
                        {
                            if(p)
                            {
                                if(coeff.lvl < max_lvl-1)
                                {
                                    temp.pos = k.row();
                                    temp.coef = temptransfer[coeff.lvl](k.row(), coeff.pos) * coeff.coef;
                                    temp.lvl = coeff.lvl+1;
                                    //gsDebug<<"temp.pos: "<<temp.pos<<" temp.coef: "<<temp.coef<<" temp.lvl "<< temp.lvl<<std::endl;
                                    coeffs.push_back(temp);
                                }
                            }
                        }
                        //}
                    }

//                        for(unsigned int ii = 0; ii < coeffs.size();ii++){
//                            gsDebug<<"( "<< coeffs[ii].pos<<" , "<<coeffs[ii].coef<<" , "<<coeffs[ii].lvl<<" )"<<std::endl;
//                        }
                }
            }
            glob_numb++;

        }
    }
    return result;
}


namespace internal
{

/// Get a Truncated Hierarchical B-spline basis from XML data
template<short_t d, class T>
class gsXml< gsTHBSplineBasis<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTHBSplineBasis<TMPLA2(d,T)>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "THBSplineBasis"+ (d>1 ? to_string(d):""); }

    static gsTHBSplineBasis<d,T> * get (gsXmlNode * node)
    {
        return getHTensorBasisFromXml< gsTHBSplineBasis<d,T> > (node);
    }

    static gsXmlNode * put (const gsTHBSplineBasis<d,T> & obj,
                            gsXmlTree & data )
    {
        return putHTensorBasisToXml< gsTHBSplineBasis<d,T> > (obj, data);
    }
};


}

} // namespace gismo
