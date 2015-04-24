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

namespace gismo
{

template<unsigned d, class T>
typename gsTHBSplineBasis<d,T>::BoundaryBasisType * gsTHBSplineBasis<d,T>::basisSlice(index_t dir_fixed,T par ) const
{
    const boxSide side(dir_fixed,0);
    const typename gsTensorBSplineBasis<d,T, gsCompactKnotVector<T> >::BoundaryBasisType * bBSplineBasis =
        this->m_bases[0]->boundaryBasis(side);
    typename gsTHBSplineBasis<d,T>::BoundaryBasisType* bBasis =
        new typename gsTHBSplineBasis<d,T>::BoundaryBasisType(*bBSplineBasis);//,this->m_tree.getMaxInsLevel()+1);

    std::vector<unsigned> boxes;
    this->getBoxesAlongSlice(dir_fixed,par,boxes);
    bBasis->refineElements(boxes);

    delete bBSplineBasis;
    return bBasis;
}

template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::representBasis()
{
    // Cleanup previous basis
    this->m_is_truncated.resize(this->size());
    m_presentation.clear();

    for (unsigned j = 0; j < static_cast<unsigned>(this->size()); ++j)
    {
        unsigned level = static_cast<unsigned>(this->levelOf(j));
        unsigned tensor_index = this->flatTensorIndexOf(j, level);

        // element indices
        gsMatrix<unsigned, d, 2> element_ind(d, 2);
        this->m_bases[level]->elementSupport_into(tensor_index, element_ind);

        // I tried with block, I can not trick the compiler to use references
        gsVector<unsigned, d> low = element_ind.col(0); //block<d, 1>(0, 0);
        gsVector<unsigned, d> high = element_ind.col(1); //block<d, 1>(0, 1);gsMatrix<unsigned> element_ind =

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

template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::_representBasisFunction(
    const unsigned j,
    const unsigned pres_level,
    const gsVector<unsigned, d>& finest_low,
    const gsVector<unsigned, d>& finest_high)
{
    const unsigned cur_level = this->levelOf(j);

    // actual size of the coefficients
    gsVector<unsigned, d> act_size_of_coefs(d);
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
    gsVector<unsigned, d> vec_nmb_of_coefs(d);
    vec_nmb_of_coefs.fill(1);

    unsigned tensor_index = this->flatTensorIndexOf(j, cur_level);

    // B-Spline vector tensor index
    gsVector<unsigned, d> bspl_vec_ti =
        this->m_bases[cur_level]->tensorIndex(tensor_index);


    // we need to separately save knot vectors because we will modify
    // them, when we proceed from one level on another
    std::vector<gsCompactKnotVector<T> > vector_of_kv(d);

    // size of the coefficients that are affected in individual iteration
    gsVector<unsigned, d> cur_size_of_coefs(d);
    cur_size_of_coefs.fill(1);

    for (unsigned level = cur_level; level < pres_level; ++level)
    {
        _updateSizeOfCoefs(level, level + 1, finest_low,
                           finest_high, cur_size_of_coefs);


        // index of a support of the j-th basis function (l_low, l_high
        // on level, and l1_high, l1_low on level + 1)
        gsVector<unsigned, d> clow, chigh, fhigh, flow;

        this->m_tree.computeLevelIndex(finest_low, level, clow);
        this->m_tree.computeLevelIndex(finest_high, level, chigh);

        this->m_tree.computeLevelIndex(finest_low, level + 1, flow);
        this->m_tree.computeLevelIndex(finest_high, level + 1, fhigh);


        for (unsigned dim = 0; dim < d; ++dim)
        {
            std::vector<T> knots;
            gsCompactKnotVector<T>& ckv =
                this->m_bases[level]->knots(dim);
            gsCompactKnotVector<T>& fkv =
                this->m_bases[level + 1]->knots(dim);


            if (level == cur_level)
                vector_of_kv[dim] = ckv;

            this->_differenceBetweenKnotVectors(ckv, clow[dim], chigh[dim],
                                                fkv, flow[dim], fhigh[dim],
                                                knots);

            gsTensorBoehmRefineLocal<d,
                                     gsCompactKnotVector<T>,
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


template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::_saveNewBasisFunPresentation(
    const gsMatrix<T>& coefs,
    const gsVector<unsigned, d>& act_size_of_coefs,
    const unsigned j,
    const unsigned pres_level,
    const gsVector<unsigned, d>& finest_low)
{
    const unsigned level = this->levelOf(j);
    const unsigned tensor_index = this->flatTensorIndexOf(j, level);

    gsVector<unsigned, d> bspl_vec_ti =
        this->m_bases[level]->tensorIndex(tensor_index);

    // finer tensor index
    const unsigned f_ten_index = _basisFunIndexOnLevel(bspl_vec_ti, level,
                                                       finest_low, pres_level);

    gsVector<unsigned, d> act_coefs_strides(d);
    bspline::buildCoeffsStrides<d>(act_size_of_coefs, act_coefs_strides);


    gsVector<unsigned, d> position(d);
    position.fill(0);


    gsVector<unsigned, d> first_point(position);
    gsVector<unsigned, d> last_point(d);
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

    } while(nextCubePoint<gsVector<unsigned, d> > (position, first_point,
                                                   last_point));
}


template<unsigned d, class T>
unsigned gsTHBSplineBasis<d,T>::_basisFunIndexOnLevel(
    const gsVector<unsigned, d>& index,
    const unsigned level,
    const gsVector<unsigned, d>& fin_low,
    const unsigned new_level)
{
    gsVector<unsigned, d> low(d);
    this->m_tree.computeLevelIndex(fin_low, level, low);

    gsVector<unsigned, d> flow(d);
    this->m_tree.computeLevelIndex(fin_low, new_level, flow);

    gsVector<unsigned, d> new_index(d);

    for (unsigned dim = 0; dim < d; dim++)
    {
        gsCompactKnotVector<T>& ckv =
            this->m_bases[level]->knots(dim);

        gsCompactKnotVector<T>& fkv =
            this->m_bases[new_level]->knots(dim);


        unsigned mult = index[dim] - ckv.firstKnotIndex(low[dim]);

        new_index(dim) = fkv.firstKnotIndex(flow[dim]) + mult;
    }

    return this->m_bases[new_level]->index(new_index);
}


template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::_truncate(
    gsMatrix<T>& coefs,
    const gsVector<unsigned, d>& act_size_of_coefs,
    const gsVector<unsigned, d>& size_of_coefs,
    const unsigned level,
    const gsVector<unsigned, d>& bspl_vec_ti,
    const unsigned bspl_vec_ti_level,
    const gsVector<unsigned, d>& finest_low)
{
    // if we dont have any active function in this level, we do not truncate
    if (this->m_xmatrix[level].size() == 0)
        return;


    // global tensor index
    const unsigned const_ten_index = _basisFunIndexOnLevel(bspl_vec_ti,
                                                           bspl_vec_ti_level, finest_low, level);
    gsVector<unsigned, d> act_coefs_strides(d);
    bspline::buildCoeffsStrides<d>(act_size_of_coefs, act_coefs_strides);


    gsVector<unsigned, d> last_point(d);
    bspline::getLastIndexLocal<d>(size_of_coefs, last_point);
    last_point(0) = 0;


    gsVector<unsigned, d> position(d);
    position.fill(0);

    gsVector<unsigned, d> first_point(position);

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

    } while(nextCubePoint<gsVector<unsigned, d> >(position, first_point,
                                                  last_point));
}


template<unsigned d, class T>
unsigned gsTHBSplineBasis<d,T>::_updateSizeOfCoefs(
    const unsigned clevel,
    const unsigned flevel,
    const gsVector<unsigned, d>& finest_low,
    const gsVector<unsigned, d>& finest_high,
    gsVector<unsigned, d>& size_of_coefs)
{
    gsVector<unsigned, d> clow, chigh;

    this->m_tree.computeLevelIndex(finest_low, clevel, clow);
    this->m_tree.computeLevelIndex(finest_high, clevel, chigh);

    gsVector<unsigned, d> flow, fhigh;
    this->m_tree.computeLevelIndex(finest_low, flevel, flow);
    this->m_tree.computeLevelIndex(finest_high, flevel, fhigh);

    // number of new coefficients
    unsigned nmb_of_coefs = 1;

    for (unsigned dim = 0; dim < d; ++dim)
    {
        gsCompactKnotVector<T>& ckv =
            this->m_bases[clevel]->knots(dim);
        gsCompactKnotVector<T>& fkv =
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
template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::getBsplinePatchGlobal(gsVector<unsigned> b1, gsVector<unsigned> b2, unsigned level, const gsMatrix<T>& geom_coef, gsMatrix<T>& cp, gsCompactKnotVector<T>& k1, gsCompactKnotVector<T>& k2) const
{
    std::vector< std::map<unsigned,T> > cmatrix;
    update_cmatrix(cmatrix);
    
    // check if the indices in b1, and b2 are correct with respect to the given level
    
    // The following should be equivalent to the scary 
    //Qlocal2global(1,0, this->m_tree.getIndexLevel()-level) that was here before.
    const unsigned loc2glob = ( 1<< (this->m_tree.getIndexLevel() - level) );
    if( b1[0]%loc2glob != 0 )
    {
        b1[0] -= b1[0]%loc2glob;
    }
    if( b1[1]%loc2glob != 0 )
    {
        b1[1] -= (b1[1]%loc2glob);
    }
    if( b2[0]%loc2glob != 0 )
    {
        b2[0] += loc2glob -(b2[0]%loc2glob);
    }
    if( b2[1]%loc2glob != 0 )
    {
        b2[1] += loc2glob -(b2[1]%loc2glob);
    }

    // select the indices of all B-splines of the given level acting on the given box

    // The following should be equivalent to the commented Qglobal2locals.
    gsVector<unsigned,d> b1_outputs, b2_outputs;

    this->m_tree.global2localIndex( b1, level, b1_outputs );
    this->m_tree.global2localIndex( b2, level, b2_outputs );

    int i0 = b1_outputs(0); //Qglobal2local(b1[0],level,this->m_tree.getIndexLevel());
    int i1 = b2_outputs(0); //Qglobal2local(b2[0],level,this->m_tree.getIndexLevel());
    int j0 = b1_outputs(1); //Qglobal2local(b1[1],level,this->m_tree.getIndexLevel());
    int j1 = b2_outputs(1); //Qglobal2local(b2[1],level,this->m_tree.getIndexLevel());

    i0 = this->m_bases[level]->knots(0).lastKnotIndex(i0) - this->m_deg[0];
    i1 = this->m_bases[level]->knots(0).firstKnotIndex(i1) - 1;
    j0 = this->m_bases[level]->knots(1).lastKnotIndex(j0) - this->m_deg[1];
    j1 = this->m_bases[level]->knots(1).firstKnotIndex(j1) - 1;

    // It would also be nice to rename b0, b1, i0, i1, j0 and j1 to something more intuitive.
    for(int i = 0; i < geom_coef.cols(); i++)
    {
        update_cmatrix(geom_coef, i, level, cmatrix);
        gsMatrix<T> temp(1,1); // ? (!) ?
        globalRefinement(level, temp, cmatrix);

        for(int k = i0; k <= i1; k++)
        {
            for(int j = j0; j <= j1; j++)
            {
                temp(j-j0,k-i0) = temp(j,k);
            }
        }
        temp.conservativeResize(j1-j0+1,i1-i0+1);

        if(i == 0)
        {
            cp.resize(temp.cols()*temp.rows(), geom_coef.cols());
        }
        return_cp_1D(temp, i, cp);
    }
    // compute the new vectors for the B-spline patch
    k1 = gsCompactKnotVector<T>(this->m_deg[0], this->m_bases[level]->knots(0).begin() + i0 , this->m_bases[level]->knots(0).begin() + i1 + this->m_deg[0] + 2);
    k2 = gsCompactKnotVector<T>(this->m_deg[1], this->m_bases[level]->knots(1).begin() + j0 , this->m_bases[level]->knots(1).begin() + j1 + this->m_deg[1] + 2);
}

// returns the list of B-spline patches to represent a THB-spline geometry
template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::getBsplinePatches(const gsMatrix<T>& geom_coef, gsMatrix<T>& cp,
                                              gsMatrix<unsigned>& b1, gsMatrix<unsigned>& b2,
                                              gsVector<unsigned>& level, gsMatrix<unsigned>& nvertices) const
{ 
    this->m_tree.getBoxes(b1,b2,level); // splitting based on the quadtree
    int nboxes = level.size();
    //------------------------------------------------------------------------------------------------------------------------------
    // iteration on the boxes to call getBsplinePatchGlobal()
    //------------------------------------------------------------------------------------------------------------------------------
    gsVector<unsigned> p1, p2;
    p1.resize(this->dim());
    p2.resize(this->dim());
    gsMatrix<T> temp1, temp2;
    gsCompactKnotVector<T> cku, ckv;
    nvertices.resize(nboxes,this->dim());

    for (int i = 0; i < nboxes; i++){
        p1(0) = b1(i,0); p1(1) = b1(i,1); p2(0) = b2(i,0); p2(1) = b2(i,1);

        this->getBsplinePatchGlobal(p1, p2, level[i], geom_coef, temp1, cku, ckv);        

        if (i == 0){
            cp = temp1;
        }else{
            int cprows = cp.rows();
            temp2.resize(cp.rows()+temp1.rows(),cp.cols());

            for (int j = 0; j < cp.rows(); j++){
                for (int k = 0; k < cp.cols(); k++){
                    temp2(j,k) = cp(j,k);
                }
            }

            for (int j = 0; j < temp1.rows(); j++){
                for (int k = 0; k < temp1.cols(); k++){                    
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
template<unsigned d, class T>
gsMultiPatch<T> gsTHBSplineBasis<d,T>::getBsplinePatchesToMultiPatch(const gsMatrix<T>& geom_coef) const
{
    gsMultiPatch<T> result;
    gsMatrix<unsigned> b1;
    gsMatrix<unsigned> b2;
    gsVector<unsigned> level;
    this->m_tree.getBoxes(b1,b2,level); // splitting based on the quadtree
    int nboxes = level.size();
    //------------------------------------------------------------------------------------------------------------------------------
    // iteration on the boxes to call getBsplinePatchGlobal()
    //------------------------------------------------------------------------------------------------------------------------------
    gsVector<unsigned> p1, p2;
    p1.resize(this->dim());
    p2.resize(this->dim());
    gsMatrix<T> temp1;
    gsCompactKnotVector<T> cku, ckv;

    for (int i = 0; i < nboxes; i++){
        p1(0) = b1(i,0); p1(1) = b1(i,1); p2(0) = b2(i,0); p2(1) = b2(i,1);

        this->getBsplinePatchGlobal(p1, p2, level[i], geom_coef, temp1, cku, ckv);
        gsTensorBSplineBasis<d, T, gsCompactKnotVector<T> > tbasis(cku, ckv);
        gsTensorBSpline<d, T, gsCompactKnotVector<T> > *tbspline = new gsTensorBSpline<d, T, gsCompactKnotVector<T> >(tbasis, give(temp1));
        result.addPatch(tbspline);

    }

    return result;
}

template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::getConnectedComponents(
    std::vector<std::vector<std::vector< std::vector<unsigned int> > > >& connectedComponents, gsVector<unsigned>& level) const
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
    gsInfo<<"new min level"<<"\n";
    std::vector< std::vector< std::vector< std::vector< unsigned int > > > > res; //things to assign to trim_curves
    std::vector< std::vector< std::vector< unsigned int > > > aabb;//axis aligned bounding box
    std::vector< std::vector< unsigned int > > boxes;
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
    for(std::size_t i = 0; i < boxes.size(); i++){
        level[i] = boxes[i][2*d];
    }

    // identify holes
    for(unsigned int l = 0; l < aabb.size();l++) //level
    {
        for(unsigned int i = 0; i < aabb[l].size(); i++) //bb
        {
            int closest_box = -1;
            for(unsigned int j = 0; j < boxes.size(); j++)
            {
                if(l == boxes[j][4])
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




//return data for trimming in parasolid
template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::getBsplinePatches_trimming(
    const gsMatrix<T>& geom_coef,
    gsMatrix<T>& cp,
    gsMatrix<unsigned>& b1,
    gsMatrix<unsigned>& b2,
    gsVector<unsigned>& level,
    gsMatrix<unsigned>& nvertices,
    std::vector<std::vector<std::vector< std::vector<T> > > >& trim_curves) const
{
    //identify the outer polylines- conected components
    int first_level = 0;
    for(unsigned int i = 0; i < this->m_xmatrix.size(); i++){
        if(this->m_xmatrix[i].size()>0){
            first_level = i-1;
            break;
        }
    }
    gsInfo<<"new min level"<<"\n";
    std::vector< std::vector< std::vector< std::vector< T > > > > res; //things to assign to trim_curves
    std::vector< std::vector< std::vector< unsigned int > > > aabb;//axis aligned bounding box
    std::vector< std::vector< unsigned int > > boxes;
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
    for(std::size_t i = 0; i < boxes.size(); i++){
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
        gsVector<unsigned> p1, p2;
        p1.resize(this->dim());
        p2.resize(this->dim());
        gsCompactKnotVector<T> cku, ckv;

        p1(0) = b1(i,0); p1(1) = b1(i,1); p2(0) = b2(i,0); p2(1) = b2(i,1);

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
    for(unsigned int l = 0; l < aabb.size();l++) //level
    {
        for(unsigned int i = 0; i < aabb[l].size(); i++) //bb
        {
            int closest_box = -1;
            for(unsigned int j = 0; j < boxes.size(); j++)
            {
                if(l == boxes[j][4])
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
}


//return data for trimming in parasolid
template<unsigned d, class T>
gsMultiPatch<T> gsTHBSplineBasis<d,T>::getBsplinePatchesToMultiPatch_trimming(
    const gsMatrix<T>& geom_coef,
    std::vector<std::vector<std::vector< std::vector<T> > > >& trim_curves) const
{
    gsMatrix<unsigned> b1;
    gsMatrix<unsigned> b2;
    gsVector<unsigned> level;
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
    gsInfo<<"new min level"<<"\n";
    std::vector< std::vector< std::vector< std::vector< T > > > > res; //things to assign to trim_curves
    std::vector< std::vector< std::vector< unsigned int > > > aabb;//axis aligned bounding box
    std::vector< std::vector< unsigned int > > boxes;
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
    for(std::size_t i = 0; i < boxes.size(); i++)
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
    gsVector<unsigned> p1, p2;
    p1.resize(this->dim());
    p2.resize(this->dim());
    gsMatrix<T> temp1;
    gsCompactKnotVector<T> cku, ckv;


    for (int i = 0; i < nboxes; i++)
    {
        p1(0) = b1(i,0); p1(1) = b1(i,1); p2(0) = b2(i,0); p2(1) = b2(i,1);

        this->getBsplinePatchGlobal(p1, p2, level[i], geom_coef, temp1, cku, ckv);
        gsTensorBSplineBasis<d, T, gsCompactKnotVector<T> > tbasis(cku, ckv);
        gsTensorBSpline<d, T, gsCompactKnotVector<T> > *tbspline = new gsTensorBSpline<d, T, gsCompactKnotVector<T> >(tbasis, give(temp1));
        result.addPatch(tbspline);

    }
    // identify holes
    for(unsigned int l = 0; l < aabb.size();l++) //level
    {
        for(unsigned int i = 0; i < aabb[l].size(); i++) //bb
        {
            int closest_box = -1;
            for(unsigned int j = 0; j < boxes.size(); j++)
            {
                if(l == boxes[j][4])
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




/**
 * Private functions
 */



// convert the coefficient matrix mat in the given direction to a column of the control points matrix
template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::return_cp_1D(const gsMatrix<T> & mat, int direction, gsMatrix<T>& cp)const
{
    GISMO_ASSERT((mat.cols()*mat.rows() == cp.rows()), "Wrong matrix dimension.");
    int counter = 0;
    for(int j = 0; j < mat.rows(); j++){
        for(int i = 0; i < mat.cols(); i++){
            cp(counter, direction) = mat(j,i);
            counter ++;
        }
    }
}

// called by getBsplinePatchGlobal
template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::globalRefinement(int level, gsMatrix<T>& coeffs, 
                                             std::vector< std::map<unsigned,T> > & cmatrix) const
{    
    //coeffs.resize(this->m_bases[0]->component(1).knots().size()-this->m_deg[1]-1,this->m_bases[0]->component(0).knots().size()-this->m_deg[0]-1);
    coeffs.resize(this->m_bases[0]->size(1),this->m_bases[0]->size(0));
    for(int j = 0; j < coeffs.rows(); j++){
        for(int k = 0; k < coeffs.cols(); k++){
            unsigned s  = this->m_bases[0]->index(k,j);//this->fromPair(k,j, 0);
            if(this->m_xmatrix[0].bContains(s) ){
                coeffs(j,k) = cmatrix[0].find( s )->second;
            }else{
                coeffs(j,k) = 0;
            }
        }
    }

    for(int l = 1; l <=level; l++){

        // global dyadic refinement with respect to previous level
        gsCompactKnotVector<T>k1;
        gsCompactKnotVector<T>k2;
        k1 = this->m_bases[l-1]->knots(0);
        k2 = this->m_bases[l-1]->knots(1);
        std::vector<T> knots_x;
        std::vector<T> knots_y;
        for(unsigned int i = 1; i < this->m_bases[l]->knots(0).unique().size(); i = i+2)
	{
            knots_x.push_back(this->m_bases[l]->knots(0).unique()[i]);
        }
        for(unsigned int i = 1; i < this->m_bases[l]->knots(1).unique().size(); i = i+2)
	{
            knots_y.push_back(this->m_bases[l]->knots(1).unique()[i]);
        }

        gsBoehmRefine(k2,coeffs,this->m_deg[1],knots_y.begin(),knots_y.end());
        coeffs.transposeInPlace();

        gsBoehmRefine(k1,coeffs,this->m_deg[0],knots_x.begin(),knots_x.end());
        coeffs.transposeInPlace();

        //overwrite the whole matrix
        for(int j = 0; j < coeffs.rows(); j++)
	{
            for(int k = 0; k < coeffs.cols(); k++)
	    {
                unsigned s  = this->m_bases[l]->index(k,j);
                if(this->m_xmatrix[l].bContains(s) )
		{
                    coeffs(j,k) = cmatrix[l].find( s )->second;
                }
            }
        }
    }
}


template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::evalSingle_into(unsigned i,
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
        
        const gsTensorBSplineBasis<d, T, gsCompactKnotVector<T> >& base =
            *this->m_bases[level];
        
        gsTensorDeboor<d, T, gsCompactKnotVector<T>, gsSparseVector<T> >
            (u, base, coefs, result);
    }
}

template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::deriv2Single_into(unsigned i,
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
        const gsTensorBSplineBasis<d, T, gsCompactKnotVector<T> >& base =
            *this->m_bases[level];
        
        gsTensorDeriv2_into<d, T, gsCompactKnotVector<T>,
                            gsSparseVector<T> >(u, base, coefs, result);
    }
}

template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    gsMatrix<unsigned> indices;
    gsMatrix<T> res(1, 1);
    this->active_into(u, indices);

    result.setZero(indices.rows(), u.cols());

    for (index_t i = 0; i < indices.cols(); i++)
    {
        for (index_t j = 0; j < indices.rows(); j++)
        {
            const unsigned index = indices(j, i);
            if (j != 0 && index == 0)
                break;
            
            evalSingle_into(index, u.col(i), res);
            result(j, i) = res.value();
        }
    }
}


template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result)const
{
    gsMatrix<unsigned> indices;
    this->active_into(u, indices);

    static const unsigned numDers = (d * (d + 1)) / 2;
    gsMatrix<T> res(numDers, 1); // result of deriv2Single_into

    result.setZero(indices.rows() * numDers, u.cols());

    for (int i = 0; i < indices.cols(); i++)
    {
        for (int j = 0; j < indices.rows(); j++)
        {
            const unsigned index = indices(j, i);
            if (j != 0 && index == 0)
                break;

            this->deriv2Single_into(index, u.col(i), res);

            result.template block<numDers, 1>(j * numDers, i) = res;
        }
    }
}


template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{

    gsMatrix<unsigned> indices;
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


template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::derivSingle_into(unsigned i,
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
        const gsTensorBSplineBasis<d, T, gsCompactKnotVector<T> >& base =
            *this->m_bases[level];
        gsTensorDeriv_into<d, T, gsCompactKnotVector<T>,
                           gsSparseVector<T> >(u, base, coefs, result);
    }

}


template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::update_cmatrix(std::vector< std::map<unsigned,T> > & cmatrix) const
{
    //srand((unsigned)time(NULL));//seed the random alg.
    cmatrix.clear();
    //initializing coefficient matrices to 0 for all active basis
    //functions
    for(size_t i = 0; i < this->m_xmatrix.size();i++)
    {
        cmatrix.push_back( std::map< unsigned, T>());
        for(typename CMatrix::const_iterator it = this->m_xmatrix[i].begin(); 
            it != this->m_xmatrix[i].end(); it++)
        {
            cmatrix[i][*it] = T(0.0);
            //cmatrix[i][*it] = 1.0/((rand()%100)+0.01);
        }
    }
}

template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::update_cmatrix(const gsMatrix<T>&geom_coeff, int col, int c_level,
                                           std::vector< std::map<unsigned,T> > & cmatrix) const
{
    int counter = 0;
    for(int i = 0; i <= c_level; i++)
    {
        cmatrix.push_back( std::map <unsigned, T>() );
        for(typename CMatrix::const_iterator it = this->m_xmatrix[i].begin();
            it != this->m_xmatrix[i].end(); it++)
        {
            cmatrix[i][ *it ] = geom_coeff(counter,col);
            counter++;
        }
    }
}


template<unsigned d, class T>
void gsTHBSplineBasis<d,T>::transferbyLvl (std::vector<gsMatrix<T> >& result){
    //std::vector< gsMatrix<T> > result;
    result.clear();
    gsVector<unsigned> level;
    gsMatrix<unsigned> b1, b2;//boxes in highes level numbering
    this->m_tree.getBoxesInLevelIndex(b1,b2,level);//return boxes in level indices
    gsTensorBSplineBasis<d,T, gsCompactKnotVector<T> > T_0_copy = this->tensorLevel(0);
    std::vector< gsSparseMatrix<T,RowMajor> > transfer;
    transfer.resize(this->maxLevel() );
    for(unsigned i = 0; i < this->maxLevel();i++)
    {
        T_0_copy.uniformRefine_withTransfer(transfer[i], 1);
    }
    std::vector< gsSparseMatrix<T,RowMajor> > temp_transf;
    for(unsigned j = 0; j < this->maxLevel();j++)
    {
        std::vector<gsSortedVector<unsigned> >x_mat_old_0;
        this->setActiveToLvl(j,x_mat_old_0);
        std::vector<gsSortedVector<unsigned> > x_matrix_lvl;
        this->setActiveToLvl(j+1,x_matrix_lvl);
        temp_transf.push_back(transfer[j]);

        gsMatrix<T> crs = this->coarsening_direct(x_mat_old_0, x_matrix_lvl, temp_transf);
        //std::cout<<"matrix ready"<<std::endl;
        result.push_back(crs);
    }
    //return result;
}

//todo remove
template<unsigned d, class T>
gsMatrix<T> gsTHBSplineBasis<d,T>::coarsening( const std::vector<gsSortedVector<unsigned> >& old, const std::vector<gsSortedVector<unsigned> >& n, const gsSparseMatrix<T,RowMajor> & transfer){
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
    gsMatrix<T> result(size2,size1);
    result.setZero();

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
                        if(n[i+1].bContains(k)){
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


template<unsigned d, class T>
gsMatrix<T> gsTHBSplineBasis<d,T>::coarsening_direct( const std::vector<gsSortedVector<unsigned> >& old,
                                                      const std::vector<gsSortedVector<unsigned> >& n,
                                                      const std::vector<gsSparseMatrix<T,RowMajor> >& transfer)
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
    gsMatrix<T> result(size2,size1);
    result.setZero();


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
            //cout<<"j...."<<j<<endl;
            start_lv_i = 0;
            for(unsigned int l =0; l < i; l++)
            {
                start_lv_i += n[l].size();
            }
            const unsigned old_ij = old[i][j];  // tensor product index

            if( n[i].bContains(old_ij) )//it he basis function was not refined
            {
                result(start_lv_i + std::distance(n[i].begin(), n[i].find_it_or_fail(old_ij) ),glob_numb ) = 1;//settign the coefficient of the not refined basis function to 1
                std::vector<lvl_coef> coeffs;
                gsMatrix<unsigned, d, 2> supp(d, 2);
                this->m_bases[i]->elementSupport_into(old_ij, supp);//this->support(start_lv_i+old_ij);
                unsigned max_lvl = math::min<unsigned>( this->m_tree.query4(supp.col(0),supp.col(1), i), transfer.size()) ;//transfer.size();//
                //std::cout<<"transfer size "<< transfer.size()<<" max lvl"<< this->m_tree.query4(supp.col(0),supp.col(1), this->levelOf(start_lv_i+j))<<"   lvl of"<< this->levelOf(start_lv_i+j)<<" support\n"<<supp<<std::endl;
                lvl_coef temp;
                temp.pos = old_ij;
                temp.coef = 1;
                temp.lvl = i;
                if(temp.lvl<n.size()-1)//in case of no refinement
                {
                    coeffs.push_back(temp);
                }
                while(!coeffs.empty())
                {
                    start_lv_i = 0;
                    for(unsigned int l =0; l < coeffs[0].lvl; l++)
                    {
                        start_lv_i += n[l].size();
                    }
//                    for(int k = 0; k < transferDense[coeffs[0].lvl].rows(); k++)//basis function was refined->looking for the coeficinets from the global transfer matrix
//                    {
//                        if(transferDense[coeffs[0].lvl](k, coeffs[0].pos) != 0)
//                        {
//                            bool p = true;
//                            if( (coeffs[0].lvl+1) < (old.size()))
//                            {
//                                if(old[coeffs[0].lvl+1].bContains(k))
//                                {
//                                    p = false;
//                                }
//                            }
//                            if(n[coeffs[0].lvl+1].bContains(k))
//                            {
//                                if(p){
//                                    const int pos = start_lv_i + n[coeffs[0].lvl].size() + std::distance(n[coeffs[0].lvl+1].begin(), n[coeffs[0].lvl+1].find_it_or_fail(k));

//                                    result(pos,glob_numb) += coeffs[0].coef * transferDense[coeffs[0].lvl](k, coeffs[0].pos);//transferDense[coeffs[0].lvl](k, coeffs[0].pos);
//                                    if(coeffs[0].lvl < transfer.size()-1)
//                                    {
//                                        temp.pos = k;
//                                        temp.coef = transferDense[coeffs[0].lvl](k, coeffs[0].pos) * coeffs[0].coef;
//                                        temp.lvl = coeffs[0].lvl+1;
//                                        coeffs.push_back(temp);
//                                    }
//                                }
//                            }else
//                            {
//                                //TODO test if level would not be higher than the q4 for this function- or max inserted level
//                                if(p){
//                                    if(coeffs[0].lvl<transfer.size()-1)
//                                    {
//                                        temp.pos = k;
//                                        temp.coef = transferDense[coeffs[0].lvl](k, coeffs[0].pos) * coeffs[0].coef;
//                                        temp.lvl = coeffs[0].lvl+1;
//                                        coeffs.push_back(temp);
//                                    }
//                                }
//                            }
//                        }
//                    }
                    for(typename gsSparseMatrix<T,ColMajor>::InnerIterator k(temptransfer[coeffs[0].lvl],coeffs[0].pos); k; ++k)//basis function was refined->looking for the coeficinets from the global transfer matrix
                    {
//                        if(transferDense[coeffs[0].lvl](k, coeffs[0].pos) != 0)
//                        {
                        //std::cout<<"first while "<< coeffs.size()<<std::endl;
                        bool p = true;
                        if( (coeffs[0].lvl+1) < (old.size()))
                        {
                            if(old[coeffs[0].lvl+1].bContains(k.row()))
                            {
                                p = false;
                            }
                        }
                        if(n[coeffs[0].lvl+1].bContains(k.row()))
                        {
                            if(p){
                                const int pos = start_lv_i + n[coeffs[0].lvl].size() + std::distance(n[coeffs[0].lvl+1].begin(), n[coeffs[0].lvl+1].find_it_or_fail(k.row()));

                                result(pos,glob_numb) += coeffs[0].coef * temptransfer[coeffs[0].lvl](k.row(), coeffs[0].pos);//transferDense[coeffs[0].lvl](k, coeffs[0].pos);
                                if(coeffs[0].lvl + 1 < max_lvl)//transfer.size()
                                {
                                    temp.pos = k.row();
                                    temp.coef = temptransfer[coeffs[0].lvl](k.row(), coeffs[0].pos) * coeffs[0].coef;
                                    temp.lvl = coeffs[0].lvl+1;
                                    coeffs.push_back(temp);
                                }
                            }
                        }else
                        {
                            //TODO test if level would not be higher than the q4 for this function- or max inserted level
                            if(p)
			    {
                                if( coeffs[0].lvl + 1< max_lvl)
                                {
                                    temp.pos = k.row();
                                    temp.coef = temptransfer[coeffs[0].lvl](k.row(), coeffs[0].pos) * coeffs[0].coef;
                                    temp.lvl = coeffs[0].lvl+1;
                                    coeffs.push_back(temp);
                                }
                            }
                        }
                        //}
                    }

//                        for(unsigned int ii = 0; ii < coeffs.size();ii++){
//                            std::cout<<"( "<< coeffs[ii].pos<<" , "<<coeffs[ii].coef<<" , "<<coeffs[ii].lvl<<" )"<<std::endl;
//                        }
//                        std::cout<<"j ="<<j<<std::endl;
                    coeffs.erase(coeffs.begin());
                }
            }
            else
            {
                gsMatrix<unsigned, d, 2> supp(d, 2);
                this->m_bases[i]->elementSupport_into(old_ij, supp);//this->support(start_lv_i+old_ij);
                unsigned max_lvl = 
                    math::min<unsigned>( this->m_tree.query4(supp.col(0),supp.col(1), i), transfer.size() ) ;
                std::vector<lvl_coef> coeffs;
                lvl_coef temp;
                temp.pos = old_ij;
                temp.coef = 1;
                temp.lvl = i;
                coeffs.push_back(temp);
                while(!coeffs.empty())
		{
                    start_lv_i = 0;
                    for(unsigned int l =0; l < coeffs[0].lvl; l++)
                    {
                        start_lv_i += n[l].size();
                    }
//                    for(int k = 0; k < transferDense[coeffs[0].lvl].rows(); k++)//basis function was refined->looking for the coeficinets from the global transfer matrix
//                    {
//                        if(transferDense[coeffs[0].lvl](k, coeffs[0].pos) != 0)
//                        {
//                            bool p = true;
//                            if( (coeffs[0].lvl+1) < (old.size()))
//                            {
//                                if(old[coeffs[0].lvl+1].bContains(k))
//                                {
//                                    p = false;
//                                }
//                            }
//                            if(n[coeffs[0].lvl+1].bContains(k))
//                            {
//                                if(p){
//                                    const int pos = start_lv_i + n[coeffs[0].lvl].size() + std::distance(n[coeffs[0].lvl+1].begin(), n[coeffs[0].lvl+1].find_it_or_fail(k));
//                                    //T ppp =  transferDense[coeffs[0].lvl](k, coeffs[0].pos);
//                                    result(pos,glob_numb) += coeffs[0].coef * transferDense[coeffs[0].lvl](k, coeffs[0].pos);//transferDense[coeffs[0].lvl](k, coeffs[0].pos);
//                                    if(coeffs[0].lvl < transfer.size()-1)
//                                    {
//                                        temp.pos = k;
//                                        temp.coef = transferDense[coeffs[0].lvl](k, coeffs[0].pos) * coeffs[0].coef;
//                                        temp.lvl = coeffs[0].lvl+1;
//                                        coeffs.push_back(temp);
//                                    }
//                                }
//                            }else{
//                                if(p)
//                                {
//                                    if(coeffs[0].lvl < transfer.size()-1)
//                                    {
//                                        temp.pos = k;
//                                        temp.coef = transferDense[coeffs[0].lvl](k, coeffs[0].pos) * coeffs[0].coef;
//                                        temp.lvl = coeffs[0].lvl+1;
//                                        coeffs.push_back(temp);
//                                    }
//                                }
//                            }
//                        }
//                    }
                    for(typename gsSparseMatrix<T,ColMajor>::InnerIterator k(temptransfer[coeffs[0].lvl],coeffs[0].pos); k; ++k)//basis function was refined->looking for the coeficinets from the global transfer matrix
                    {
                        //if(transferDense[coeffs[0].lvl](k, coeffs[0].pos) != 0)
                        //{
                        //std::cout<<"second while "<< coeffs.size()<<std::endl;
                        bool p = true;
                        if( (coeffs[0].lvl+1) < (old.size()))
                        {
                            if(old[coeffs[0].lvl+1].bContains(k.row()))
                            {
                                p = false;
                            }
                        }
                        if(n[coeffs[0].lvl+1].bContains(k.row()))
                        {
                            if(p){
                                const int pos = start_lv_i + n[coeffs[0].lvl].size() + std::distance(n[coeffs[0].lvl+1].begin(), n[coeffs[0].lvl+1].find_it_or_fail(k.row()));
                                //T ppp =  transferDense[coeffs[0].lvl](k, coeffs[0].pos);
                                result(pos,glob_numb) += coeffs[0].coef * temptransfer[coeffs[0].lvl](k.row(), coeffs[0].pos);//transferDense[coeffs[0].lvl](k, coeffs[0].pos);
                                if(coeffs[0].lvl < max_lvl-1)
                                {
                                    temp.pos = k.row();
                                    temp.coef = temptransfer[coeffs[0].lvl](k.row(), coeffs[0].pos) * coeffs[0].coef;
                                    temp.lvl = coeffs[0].lvl+1;
                                    coeffs.push_back(temp);
                                }
                            }
                        }else
			{
                            if(p)
                            {
                                if(coeffs[0].lvl < max_lvl-1)
                                {
                                    temp.pos = k.row();
                                    temp.coef = temptransfer[coeffs[0].lvl](k.row(), coeffs[0].pos) * coeffs[0].coef;
                                    temp.lvl = coeffs[0].lvl+1;
                                    coeffs.push_back(temp);
                                }
                            }
                        }
                        //}
                    }

//                        for(unsigned int ii = 0; ii < coeffs.size();ii++){
//                            std::cout<<"( "<< coeffs[ii].pos<<" , "<<coeffs[ii].coef<<" , "<<coeffs[ii].lvl<<" )"<<std::endl;
//                        }
                    coeffs.erase(coeffs.begin());
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
template<unsigned d, class T>
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
