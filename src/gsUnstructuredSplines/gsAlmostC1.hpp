/** @file gsAlmostC1.hpp

    @brief Creates the D-Patch smoothing matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gsIO/gsWriteParaview.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSpline.h>
#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsMatrixOp.h>
#include <typeinfo>

#include <gsAssembler/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>

namespace gismo
{
    // Constructors
    template<short_t d,class T>
    gsAlmostC1<d,T>::gsAlmostC1(const gsMultiPatch<T> & patches)
    :
    m_patches(patches)
    {
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            for (short_t dim=0; dim!=d; dim++)
            GISMO_ASSERT(m_patches.basis(p).degree(dim)==2,"Degree of the basis ( dimension "<<dim<<" ) of patch "<<p<<" is "<<m_patches.basis(p).degree(dim)<<", but should be 2!");

        this->_initialize();
        this->_computeMapper();
        this->_computeSmoothMatrix();
        // m_RefPatches = m_patches;
        this->_makeTHB();
        this->_computeEVs();
    }

    template<short_t d,class T>
    gsAlmostC1<d,T>::gsAlmostC1(const gsAlmostC1& other)
    :
    m_patches(other.m_patches)
    {
        GISMO_NO_IMPLEMENTATION;
    }

    template<short_t d,class T>
    gsAlmostC1<d,T>::~gsAlmostC1()
    {
        freeAll(m_bases);
    }

    /*=====================================================================================
                                    Output functions
    =====================================================================================*/
    // GIVES SEGFAULT??

    // template<short_t d,class T>
    // gsMultiBasis<T> gsAlmostC1<d,T>::makeGlobalBasis()
    // {
    //     m_MBasis = gsMappedBasis<d,T>(m_bases,m_matrix.transpose());

    //     m_basisContainer = std::vector<gsBasis<T> *>(m_patches.nPatches());
    //     for (size_t p=0; p!=m_patches.nPatches(); p++)
    //     {
    //         gsMappedSingleBasis<2,real_t>::uPtr mbasis = gsMappedSingleBasis<2,real_t>::make(m_MBasis.getMappedSingleBasis(p));
    //         m_basisContainer[p] = static_cast<gsBasis<> *>(mbasis.release());
    //     }

    //     m_globalBasis = gsMultiBasis<T>(m_basisContainer,m_patches);
    //     return m_globalBasis;
    // }

    /*=====================================================================================
                                    Special functions
    =====================================================================================*/
    template<short_t d,class T>
    gsMatrix<T> gsAlmostC1<d,T>::_getNormals(const std::vector<patchCorner> & corners) const
    {
        gsMatrix<T> normals(3,corners.size());

        gsVector<bool> pars;
        gsMatrix<T> mat;

        gsExprEvaluator<T> ev;
        typename gsExprEvaluator<T>::geometryMap Gm = ev.getMap(m_patches);
        index_t k = 0;
        for (typename std::vector<patchCorner>::const_iterator it = corners.begin(); it!=corners.end(); it++, k++)
        {
            it->corner().parameters_into(m_patches.parDim(),pars); // get the parametric coordinates of the corner
            mat = pars.template cast<T>(); // cast to real coordinates
            normals.col(k) = ev.eval(sn(Gm).normalized(),mat,it->patch);
        }
        return normals;
    }


    template<short_t d,class T>
    std::tuple<gsMatrix<T>,gsMatrix<T>,gsMatrix<index_t>> gsAlmostC1<d,T>::_makeTriangle(const patchCorner & corner) const
    {
        GISMO_ASSERT(m_RefPatches.nPatches()!=0,"Are the patches refined?");

        index_t tdim = m_RefPatches.targetDim();

        std::vector<patchCorner> corners;
        m_RefPatches.getCornerList(corner,corners);

        gsVector<bool> pars;
        gsMatrix<T> mat;
        // 1. Get the coordinates of the vertex and set its z coordinate to 0
        gsMatrix<T> um(3,1), midpoint;
        um.setZero();
        corner.corner().parameters_into(m_RefPatches.parDim(),pars); // get the parametric coordinates of the corner
        mat = pars.template cast<T>(); // cast to real coordinates
        um.block(0,0,tdim,1) = m_RefPatches.patch(corner.patch).eval(mat);
        midpoint = um; // store the original midpoint

        // 2. Get the 0,0;0,1; 1,0; 1,1 coordinates
        gsMatrix<T> u(3,corners.size()*4);
        u.setZero();
        gsMatrix<index_t> uind(1,corners.size()*4);
        uind.setZero();

        std::vector<patchSide> csides;
        index_t idx;
        for (size_t c = 0; c!=corners.size(); c++)
        {
            corners[c].getContainingSides(d,csides);
            index_t k=0;
            for (index_t i=0; i!=2; i++)
                for (index_t j=0; j!=2; j++,k++)
                {
                    idx = _indexFromVert(i,corners[c],csides[0],j);
                    uind(0,4*c+k) = m_mapOriginal.index(idx,corners[c].patch);
                    u.block(0,4*c+k,m_RefPatches.targetDim(),1) = m_RefPatches.patch(corners[c].patch).coefs().row(idx).transpose();
                }
        }

        // 3. Translate all points to a coordinate system with origin um
        gsMatrix<T> up = u;
        for (index_t k=0; k!=up.cols(); k++)
            up.col(k) -= um;

        // 4. Rotate the points parallel the xy-plane and set their z-coordinates to 0
        gsMatrix<T,3,3> Rn, Rx;
        Rn.setIdentity();
        if (m_RefPatches.targetDim()==2)
        {
            // do nothing
        }
        else if(m_RefPatches.targetDim()==3)
        {
            // Get the average normal at the corner
            gsVector<T> avgnormal = _getNormals(corners).rowwise().mean();

            // Find the rotation matrix that maps the average normal to the z axis
            gsVector<T,3> ez;
            ez<<0,0,1;
            Rn = _getRotationMatrix(avgnormal.normalized(),ez);

            for (index_t k=0; k!=up.cols(); k++)
                up.col(k).applyOnTheLeft(Rn);

            up.row(2).setZero(); // all points
            um.row(2).setZero();// midpoint
        }
        else
            GISMO_ERROR("Target dimension of the multipatch should be 2 or 3, but is "<<m_RefPatches.targetDim());

        // 5. Find the maximum distance from the midpoint to all points
        T distance, maxDistance = 0;
        gsMatrix<T> umax;
        for (index_t k = 0; k!=up.cols(); k++)
        {
            distance = (up.col(k)).norm();
            if (distance > maxDistance)
            {
                maxDistance = distance;
                umax = up.col(k);
            }
        }

        gsVector<T,3> ex;
        ex<<1,0,0;

        // 6. Rotate all points such that the maximum point is aligned with the x-axis
        Rx = _getRotationMatrix(umax.normalized(),ex);
        for (index_t k=0; k!=up.cols(); k++)
            up.col(k).applyOnTheLeft(Rx);

        // 7. Obtain the coordinates of the triangle that encloses the circle with radius maxDistance in the xy plane
        T r = maxDistance;
        T a = 1. / ( 1./6. * std::sqrt(3) ) * r;
        T rr = 1. / 3. * std::sqrt(3) * a;

        gsMatrix<T> Cp(2,3);
        Cp.col(0)<<rr,0;
        Cp.col(1)<<-r, 0.5*a;
        Cp.col(2)<<-r,-0.5*a;

        // 8. Get the barycentric coordinates of the points
        gsMatrix<T> ub = up;
        up.row(2).setOnes(); // project the parametric points to z=1
        gsMatrix<T> A(3,3);
        A.block(0,0,2,3) = Cp;
        A.row(2).setOnes();

        for (index_t k = 0; k!=ub.cols(); k++)
        {
            ub.col(k) = A.colPivHouseholderQr().solve(up.col(k));
            GISMO_ASSERT((Cp * ub.col(k)-up.col(k).head(2)).norm()<1e-14,"Something went wrong with the computation of the barycentric coordinates");
        }

        // 9. Move the corners of the triangle back to physical coordinates
        gsMatrix<T> Cg(tdim,3);
        Cg.setZero();
        Cg.block(0,0,2,3) = Cp;

        for (index_t k = 0; k!=Cg.cols(); k++)
        {
            Cg.col(k).applyOnTheLeft((Rx).transpose());
            Cg.col(k).applyOnTheLeft((Rn).transpose());
            Cg.col(k) += midpoint;
        }

        if (m_RefPatches.targetDim()==2)
            Cg.conservativeResize(2,Eigen::NoChange);

        return std::make_tuple(Cg,ub,uind);
    }



    template<short_t d,class T>
    void gsAlmostC1<d,T>::_toBarycentricCoordinates(const gsMatrix<T> & Cs, gsMatrix<T> & u) const
    {

    }

    template<short_t d,class T>
    gsMatrix<T,3,3> gsAlmostC1<d,T>::_getRotationMatrix(const gsVector<T,3> & a, const gsVector<T,3> & b) const
    {
        GISMO_ASSERT(std::abs(a.norm()-1)<1e-14,"A must be a unit vector, a.norm() = "<<std::abs(a.norm()-1));
        GISMO_ASSERT(std::abs(b.norm()-1)<1e-14,"A must be a unit vector, b.norm() = "<<std::abs(b.norm()-1));

        gsVector<T,3> v = a.cross(b);
        v.normalize();
        T theta = std::acos( a.dot(b) / ( a.norm() * b.norm() ) );

        T s = std::sin(theta);
        T c = std::cos(theta);
        gsMatrix<T,3,3> R,vx,tmp, I;
        R.setZero();
        vx.setZero();

        vx.row(0)<<0,-v.at(2),v.at(1);
        vx.row(1)<<v.at(2),0,-v.at(0);
        vx.row(2)<<-v.at(1),v.at(0),0;

        I.setIdentity();
        R += I*c;
        R += vx * s;
        tmp = (v*v.transpose()) * (1-c);
        R += tmp;

        GISMO_ASSERT((R * a - b).norm() < 1e-12,"Rotation matrix is wrong, R*a = "<<R*a<<"; b = "<<b);
        return R;
    }

    /*=====================================================================================
                                    Information functions
    =====================================================================================*/

    template<short_t d,class T>
    void gsAlmostC1<d,T>::mapperInfo() const
    {
        for (size_t p=0; p!=m_patches.nPatches(); p++)
        {
            gsInfo<<"----------------------------------\n";
            gsInfo<<"patch "<<p<<"\n";
            index_t size = m_mapModified.patchSize(p);
            for (index_t k=0; k!=size; k++)
            {
                bool free = m_mapModified.is_free(k,p);
                std::string str = free ? "free" : "eliminated";
                gsInfo<<"DoF "<<k<<" is "<<str<<"\n";
            }
        }
    }

    template<short_t d,class T>
    const void gsAlmostC1<d,T>::vertexInfo(patchCorner corner) const
    {
        std::pair<index_t,bool> data = this->_vertexData(corner);
        gsInfo<<"Patch "<<corner.patch<<", corner "<<corner<<" has valence "<<data.first<<" and is "<<(data.second ? "an interior vertex" : "a boundary vertex")<<"\n";

    }

    template<short_t d,class T>
    const void gsAlmostC1<d,T>::sideInfo(patchSide side) const
    {
        gsInfo<<"Patch "<<side.patch<<", side "<<side<<" is "<<(m_patches.isBoundary(side) ? "a boundary side" : "an interface")<<"\n";
    }

    template<short_t d,class T>
    const void gsAlmostC1<d,T>::sideInfo() const
    {
        gsInfo<<"**D-Patch Side info**\n";
        for(size_t i = 0;i<m_patches.nPatches();++i)
            for(int j=1;j<=4;++j)
                sideInfo(patchSide(i,j));
    }

    template<short_t d,class T>
    const void gsAlmostC1<d,T>::cornerInfo() const
    {
        gsInfo<<"**D-Patch Corner info**\n";
        for(size_t i = 0;i<m_patches.nPatches();++i)
            for(int j=1;j<=4;++j)
                vertexInfo(patchCorner(i,j));
    }

    /*=====================================================================================
                                    Coefficients
    =====================================================================================*/
    // ADD THE COEFFICIENTS OF THE TRIANGLES AS EXTRA COEFFICIENTS

    template<short_t d,class T>
    gsMatrix<T> gsAlmostC1<d,T>::preCoefficients()
    {
        GISMO_ASSERT(m_mapModified.isFinalized(),"Mapper is not finalized, run XXXX first");

        gsMatrix<T> coefs(m_mapModified.freeSize(),m_patches.geoDim());

        index_t size;
        for (size_t p=0; p!=m_patches.nPatches(); p++) // patches
        {
            size = m_mapModified.patchSize(p);
            for (index_t k=0; k!=size; k++)
            {
                if (m_mapModified.is_free(k,p))
                    coefs.row(m_mapModified.index(k,p,0)) = m_patches.patch(p).coefs().row(k);
            }
        }

        // Correct the EVs
        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);
        index_t cidx;
        std::vector<patchCorner> pcorners;
        patchCorner pcorner;
        for (size_t p=0; p!=m_patches.nPatches(); p++)
        {
            for (index_t c=1; c<=4; c++)
            {
                cidx = _vertIndex(p,c);
                if (m_vertCheck.at(cidx))
                    continue;

                pcorner = patchCorner(p,c);
                std::pair<index_t,bool> vdata = _vertexData(pcorner); // corner c
                if (vdata.first > 2 && !(vdata.first==4 && vdata.second)) // valence must be 3 or larger, but it must not be an interior vertex with v=4
                {
                    m_RefPatches.getCornerList(pcorner,pcorners);

                    // get the triangle
                    gsMatrix<T> Cg;
                    std::tie(Cg,std::ignore,std::ignore) = _makeTriangle(pcorner);

                    // The barycentric coordinates will be attached to the matrix rows corresponding to the 0,0 corners (and the three lowest patch index corners whenever valence > 3)
                    // We use _getLowestCorners such that the corners are assigned to increasing patch corners
                    _getLowestCorners(pcorners,3);

                    std::vector<index_t> rowIndices;
                    std::vector<patchSide> sides(2);
                    for (size_t k=0; k!=pcorners.size(); k++ )
                    {
                        pcorners[k].getContainingSides(d,sides);
                        rowIndices.push_back(m_mapModified.index(_indexFromVert(m_Bbases,0,pcorners[k],sides[0],0),pcorners[k].patch));
                    }

                    index_t rowIdx;
                    for (index_t j=0; j!=Cg.cols(); j++)
                    {
                        rowIdx = rowIndices[j];
                        coefs.row(rowIdx) = Cg.col(j).transpose();
                    }

                    for (size_t k = 0; k!=pcorners.size(); k++)
                        m_vertCheck[ _vertIndex(pcorners[k].patch, pcorners[k].corner()) ] = true;
                }
                else
                    continue;
            }
        }

        return coefs;
    }

    template<short_t d,class T>
    gsMatrix<T> gsAlmostC1<d,T>::allCoefficients() const
    {
        std::vector<index_t> sizes(m_patches.nPatches());
        index_t totalsize = 0;
        for (size_t p=0; p!=m_patches.nPatches(); p++) // patches
        {
            sizes.at(p) = m_patches.patch(p).coefs().rows();
            totalsize += sizes.at(p);
        }

        gsMultiBasis<T> basis(m_patches);
        gsDofMapper tmpMap(basis);
        tmpMap.finalize();

        gsMatrix<T> coefs(totalsize,m_patches.geoDim());
        index_t offset = 0;
        for (size_t p=0; p!=m_patches.nPatches(); p++) // patches
        {
            for (index_t k=0; k!=sizes.at(p); k++)
            {
                    coefs.row(tmpMap.index(k,p)) = m_patches.patch(p).coefs().row(k);
            }
            offset += sizes.at(p);
        }

        return coefs;
    }

    /*=====================================================================================
                                    Output functions
    =====================================================================================*/
    template<short_t d,class T>
    gsGeometry<T>* gsAlmostC1<d,T>::exportPatch(index_t patch, bool computeCoefs)
    {
        ////////////////////////////////////////////////
        // This can be done more efficient!!
        // Do it once instead of for every patch
        ////////////////////////////////////////////////
        if (computeCoefs)
        {
            m_coefs = this->preCoefficients(); // gets coefficients of the modified size
            m_coefs = m_matrix.transpose() * m_coefs; // maps to local size
        }

        ////////////////////////////////////////////////
        index_t size,offset = 0;
        for (index_t p=0; p!=patch; p++)
            offset += m_mapOriginal.patchSize(p);

        size = m_mapOriginal.patchSize(patch);
        gsMatrix<T> local = m_coefs.block(offset,0,size,m_patches.geoDim());
        return m_bases[patch].makeGeometry( give(local) ).release();
    }

    template<short_t d,class T>
    gsMultiPatch<T> gsAlmostC1<d,T>::exportToPatches()
    {
        m_coefs = this->preCoefficients();
        m_coefs = m_matrix.transpose() * m_coefs;

        std::vector<gsGeometry<T> *> patches(m_patches.nPatches());
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            patches[p]= this->exportPatch(p,false);

        return gsMultiPatch<T>(patches,m_patches.boundaries(),m_patches.interfaces());
    }

    /*=====================================================================================
                                    Output functions
    =====================================================================================*/

    template<short_t d,class T>
    const index_t gsAlmostC1<d,T>::_indexFromSides(index_t index1, const patchSide side1, index_t index2, const patchSide side2)
    {
        /*
            Finds the index index1 away from side1 and index2 from side2
            index 1 is the index parallel to side 1
            index 2 is the index parallel to side 2
        */
        GISMO_ASSERT(side1.patch==side2.patch,"Sides must be from the same patch");
        GISMO_ASSERT(side1.side().direction()!=side2.side().direction(),"Sides must have different direction");
        index_t index;

        gsBasis<T> * basis = &m_bases.basis(side1.patch);

        gsVector<index_t> indices1 = static_cast<gsVector<index_t>>(basis->boundaryOffset(side1.side(),index2));

        index_t n = indices1.rows();
        if (side1.side()==1) //west
            if (side2.side().parameter()==0) //south
                index = indices1(index1);
            else
                index = indices1(n-1-index1); //north
        else if (side1.side()==2) //east
            if (side2.side().parameter()==0) //south
                index = indices1(index1);
            else
                index = indices1(n-1-index1); //north
        else if (side1.side()==3) //south
            if (side2.side().parameter()==0) //west
                index = indices1(index1);
            else
                index = indices1(n-1-index1); //east
        else if (side1.side()==4) //north
            if (side2.side().parameter()==0) //west
                index = indices1(index1);
            else
                index = indices1(n-1-index1); //east
        else
            GISMO_ERROR("Side unknown. index = "<<side1.side());
        return index;
    }

    template<short_t d,class T>
    const gsVector<index_t> gsAlmostC1<d,T>::_indicesFromVert(index_t index, const patchCorner corner, const patchSide side, index_t offset)
    {
        /*
            Finds indices i1,...,in in the direction of side away from the vertex
        */

        gsVector<index_t> result(index);

        gsBasis<T> * basis = &m_bases.basis(side.patch);

        gsVector<index_t> indices = static_cast<gsVector<index_t>>(basis->boundaryOffset(side.side(),offset));
        if (side.side()==1) //west
        {
            if (corner.corner()==1)//southwest
                result = indices.head(index);
            else if (corner.corner()==3) //northwest
                result = indices.tail(index);
            else
                GISMO_ERROR(corner.corner() << "is not on side "<<side.side()<<"!");
        }
        else if (side.side()==2) //east
        {
            if (corner.corner()==2)//southeast
                result = indices.head(index);
            else if (corner.corner()==4) //northeast
                result = indices.tail(index);
            else
                GISMO_ERROR(corner.corner() << "is not on side "<<side.side()<<"!");
        }
        else if (side.side()==3) //south
        {
            if (corner.corner()==1)//southwest
                result = indices.head(index);
            else if (corner.corner()==2) //southeast
                result = indices.tail(index);
            else
                GISMO_ERROR(corner.corner() << "is not adjacent to side "<<side.side()<<"!");
        }
        else if (side.side()==4) //north
        {
            if (corner.corner()==3)//northwest
                result = indices.head(index);
            else if (corner.corner()==4) //northeast
                result = indices.tail(index);
            else
                GISMO_ERROR(corner.corner() << "is not adjacent to side "<<side.side()<<"!");
        }
        return result;
    }

    template<short_t d,class T>
    const index_t gsAlmostC1<d,T>::_indexFromVert(index_t index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        return _indexFromVert(m_bases,index, corner, side, offset, levelOffset);
    }

    template<short_t d,class T>
    const index_t gsAlmostC1<d,T>::_indexFromVert(gsMultiBasis<T> bases, index_t index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        if ((index==0) && (offset==0))
        {
            gsHTensorBasis<d,T> *basis = dynamic_cast<gsHTensorBasis<d,T>*>(&bases.basis(corner.patch));
            index_t idx = basis->functionAtCorner(corner.corner());
            return idx;
        }
        else
        {
            std::vector<index_t> indices(1);
            indices.at(0) = index;
            std::vector<index_t> result = _indexFromVert(bases,indices,corner,side,offset,levelOffset);
            return result.at(0);
        }
    }

    template<short_t d,class T>
    const std::vector<index_t> gsAlmostC1<d,T>::_indexFromVert(std::vector<index_t> index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        return _indexFromVert(m_bases,index,corner,side,offset,levelOffset);
    }

    template<short_t d,class T>
    const std::vector<index_t> gsAlmostC1<d,T>::_indexFromVert(gsMultiBasis<T> bases, std::vector<index_t> index, const patchCorner corner, const patchSide side, index_t offset, index_t levelOffset) const
    {
        /*
            Finds indices i in the direction of side away from the vertex
            if index = 0, the corner index is requested
        */

        std::vector<index_t> result(index.size());

        // gsBasis<T> * basis = &bases.basis(corner.patch);
        gsHTensorBasis<2,T> *basis = dynamic_cast<gsHTensorBasis<2,T>*>(&bases.basis(corner.patch));
        index_t level = basis->levelAtCorner(corner.corner()) + levelOffset;

        gsVector<index_t> indices = static_cast<gsVector<index_t>>(basis->boundaryOffset(side.side(),offset,level));
        index_t end = indices.rows()-1;

        for (size_t k=0; k!=index.size(); k++)
        {
            if (side.side()==1) //west
            {
                if (corner.corner()==1)//southwest
                    result[k] = indices.at(index[k]);
                else if (corner.corner()==3) //northwest
                    result[k] = indices.at(end-index[k]);
                else
                    GISMO_ERROR(corner.corner() << " is not adjacent to side "<<side.side()<<"!");
            }
            else if (side.side()==2) //east
            {
                if (corner.corner()==2)//southeast
                    result[k] = indices.at(index[k]);
                else if (corner.corner()==4) //northeast
                    result[k] = indices.at(end-index[k]);
                else
                    GISMO_ERROR(corner.corner() << " is not adjacent to side "<<side.side()<<"!");
            }
            else if (side.side()==3) //south
            {
                if (corner.corner()==1)//southwest
                    result[k] = indices.at(index[k]);
                else if (corner.corner()==2) //southeast
                    result[k] = indices.at(end-index[k]);
                else
                    GISMO_ERROR(corner.corner() << " is not adjacent to side "<<side.side()<<"!");
            }
            else if (side.side()==4) //north
            {
                if (corner.corner()==3)//northwest
                    result[k] = indices.at(index[k]);
                else if (corner.corner()==4) //northeast
                    result[k] = indices.at(end-index[k]);
                else
                    GISMO_ERROR(corner.corner() << " is not adjacent to side "<<side.side()<<"!");

            }
            if (levelOffset==0) // if levelOffset!=0, the index in the original basis is requested
                result[k] = basis->flatTensorIndexToHierachicalIndex(result[k],level);
        }

        return result;
    }

    template<short_t d,class T>
    const std::pair<index_t,bool> gsAlmostC1<d,T>::_vertexData(const patchCorner corner) const
    {
        std::vector<patchCorner> corners;
        std::pair<index_t,bool> output;
        output.second = m_patches.getCornerList(corner,corners); // bool is true if interior vertex
        output.first = corners.size();
        return output;
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_getLowestCorners(std::vector<patchCorner> & pcorners, index_t n) const
    {
        struct {
            bool operator()(patchCorner a, patchCorner b) const { return a.patch < b.patch; }
        } customLess;

        // Sort
        std::sort(pcorners.begin(), pcorners.end(),customLess);
        // Resize; this vector are the indices we want to keep the DoFs of
        pcorners.resize(n);
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_removeLowestCorners(std::vector<patchCorner> & pcorners, index_t n) const
    {
        struct {
            bool operator()(patchCorner a, patchCorner b) const { return a.patch > b.patch; }
        } customGreater;

        // Sort
        std::sort(pcorners.begin(), pcorners.end(),customGreater);
        // Resize; this vector are the indices we want to keep the DoFs of
        pcorners.resize(n);
    }

    /*=====================================================================================
                                    Construction functions
    =====================================================================================*/


    template<short_t d,class T>
    void gsAlmostC1<d,T>::_makeTHB()
    {
        m_RefPatches = gsMultiPatch<T>(m_patches);
        // gsMultiPatch<T> refPatches(m_patches);
        typename gsBlockOp<T>::Ptr thbMat=gsBlockOp<T>::make(m_RefPatches.nPatches(),m_RefPatches.nPatches());
        // prepare the geometry
        std::vector<std::vector<patchCorner> > cornerLists;
        m_RefPatches.getEVs(cornerLists,true);

        if (cornerLists.size()!=0)
        {
            std::vector< std::vector<index_t> > elVec(m_RefPatches.nPatches());
            for (size_t v =0; v!=cornerLists.size(); v++)
                for (size_t c = 0; c!=cornerLists[v].size(); c++)
                {
                    patchCorner corner = cornerLists[v].at(c);
                    gsVector<bool> pars;
                    corner.corner().parameters_into(m_RefPatches.parDim(),pars); // get the parametric coordinates of the corner
                    gsMatrix<T> mat = pars.template cast<T>(); // cast to real coordinates

                    gsMatrix<T> boxes(m_RefPatches.parDim(),2);
                    boxes.col(0) << mat;
                    boxes.col(1) << mat;

                    gsHTensorBasis<2,T> *basis = dynamic_cast<gsHTensorBasis<2,T>*>(&m_RefPatches.basis(corner.patch));
                    std::vector<index_t> elements = basis->asElements(boxes,0); // 0-ring

                    elVec.at(corner.patch).insert(elVec.at(corner.patch).end(), elements.begin(), elements.end());

                    // gsHTensorBasis<2,T> *basis = dynamic_cast<gsHTensorBasis<2,T>*>(&m_RefPatches.basis(corner.patch));

                    // basis->refineElements(elements, m_tMatrix);
                }

            std::vector<gsSparseMatrix<T>> xmatrices(m_RefPatches.nPatches());
            gsSparseMatrix<T> tmp;
            index_t rows = 0, cols = 0;
            for (size_t p=0; p!=m_RefPatches.nPatches(); p++)
            {
                gsHTensorBasis<2,T> *basis = dynamic_cast<gsHTensorBasis<2,T>*>(&m_RefPatches.basis(p));
                std::vector< gsSortedVector< index_t > > xmat = basis->getXmatrix();

                m_RefPatches.patch(p).refineElements(elVec[p]);

                basis->transfer(xmat,tmp);
                xmatrices[p] = tmp;

                rows += tmp.rows();
                cols += tmp.cols();

                thbMat->addOperator(p,p,makeMatrixOp(tmp) );
            }

            index_t rowOffset = 0, colOffset = 0;
            m_tMatrix.resize(rows,cols);
            for (size_t p=0; p!=m_RefPatches.nPatches(); p++)
            {
                for (index_t i = 0; i!=xmatrices[p].rows(); i++)
                    for (index_t j = 0; j!=xmatrices[p].cols(); j++)
                    {
                        if ( 0 == xmatrices[p](i,j) ) continue;
                        m_tMatrix.coeffRef(rowOffset + i,colOffset + j) = xmatrices[p](i,j);
                    }
                rowOffset += xmatrices[p].rows();
                colOffset += xmatrices[p].cols();
            }

            m_tMatrix.makeCompressed();
            m_bases = gsMultiBasis<T>(m_RefPatches);
        }

        // redefine the mappers
        // m_mapModified = gsDofMapper(m_bases);
        m_mapOriginal = gsDofMapper(m_bases);
        m_mapOriginal.finalize();

        // gsWriteParaview<>(m_RefPatches,"mp_ref",1000,true);
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_computeEVs()
    {
        /*
            Our goal is to create three vectors c11, c12, c21 which all contain the
            c11, c12 and c21 coefficients of the patches around the EV in the right order
            (counter)-clockwise.
        */

        std::vector<std::vector<patchCorner> > cornerLists;
        m_patches.getEVs(cornerLists);

        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);
        if (cornerLists.size()!=0)
        {
            std::vector<patchCorner> pcorners;
            patchCorner pcorner;
            gsMatrix<T> Cg;         // coefficients
            gsMatrix<T> ub;         // baricentric coordinates
            gsMatrix<index_t> uind; // column indices of baricentric coordinates
            index_t cidx;

            m_matrix = m_matrix * m_tMatrix.transpose();

            for (size_t p=0; p!=m_patches.nPatches(); p++)
            {
                for (index_t c=1; c<=4; c++)
                {
                    cidx = _vertIndex(p,c);
                    if (m_vertCheck.at(cidx))
                        continue;

                    pcorner = patchCorner(p,c);
                    std::pair<index_t,bool> vdata = _vertexData(pcorner); // corner c
                    if (vdata.first > 2 && !(vdata.first==4 && vdata.second)) // valence must be 3 or larger, but it must not be an interior vertex with v=4
                    {
                        m_RefPatches.getCornerList(pcorner,pcorners);

                        // get the triangle
                        std::tie(std::ignore,ub,uind) = _makeTriangle(pcorner);

                        // The barycentric coordinates will be attached to the matrix rows corresponding to the 0,0 corners (and the three lowest patch index corners whenever valence > 3)
                        // We use _getLowestCorners such that the corners are assigned to increasing patch corners
                        _getLowestCorners(pcorners,3);

                        std::vector<index_t> rowIndices;
                        rowIndices.reserve(3);
                        index_t idx;
                        std::vector<patchSide> sides(2);
                        for (size_t k=0; k!=pcorners.size(); k++ )
                        {
                            pcorners[k].getContainingSides(d,sides);
                            // We need the index on the old basis (the unrefined basis), because we plug it into the mapModified (which maps the local DoFs to the global ones)
                            idx = _indexFromVert(m_Bbases,0,pcorners[k],sides[0],0); // On the 'old' (unrefined) multibasis
                            GISMO_ASSERT(m_mapModified.is_free(idx,pcorners[k].patch),"This DoF must be free!");
                            rowIndices.push_back(m_mapModified.index(idx,pcorners[k].patch));
                        }

                        index_t rowIdx,colIdx;
                        // set the colums related to the barycentric columsn equal to zero
                        for (index_t j=0; j!=ub.cols(); j++)
                        {
                            colIdx = uind(0,j);
                            m_matrix.prune(
                                            [&colIdx](index_t i, index_t j, T)
                                            { return j!=colIdx; }
                                            );
                        }

                        for (index_t i=0; i!=ub.rows(); i++)
                            for (index_t j=0; j!=ub.cols(); j++)
                            {
                                rowIdx = rowIndices[i];
                                colIdx = uind(0,j);
                                m_matrix(rowIdx,colIdx) = ub(i,j);
                            }

                        for (size_t k = 0; k!=pcorners.size(); k++)
                            m_vertCheck[ _vertIndex(pcorners[k].patch, pcorners[k].corner()) ] = true;
                    }
                    else
                        continue;
                }
            }
        }
        m_matrix.makeCompressed();
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_initialize() // also initialize the mappers!
    {
        m_verbose = false;
        m_patches.checkConsistency();
        size_t nSides = 2*m_patches.nInterfaces() + m_patches.nBoundary();
        size_t nVerts = 4*m_patches.nPatches();

        m_sideCheck.resize(nSides);
        std::fill(m_sideCheck.begin(), m_sideCheck.end(), false);
        m_vertCheck.resize(nVerts);
        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);

        m_bases = gsMultiBasis<T>(m_patches);

        m_mapModified = gsDofMapper(m_bases);
        m_mapOriginal = gsDofMapper(m_bases);

        m_size = -1;

        // Cast all patches of the mp object to THB splines
        gsTHBSpline<d,T> thb;
        for (size_t k=0; k!=m_patches.nPatches(); ++k)
        {
            if (typeid(gsTHBSplineBasis<2,T>)!=typeid(m_bases.basis(k)))
            {
                gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&m_patches.patch(k));

                thb = gsTHBSpline<2,real_t>(*geo);
                m_patches.patch(k) = thb;
            }
        }


        /*

            v   b   b   b   b   b   b   v   |   v   b   b   b   b   b   b   v           o   o   o   o   o   o   o   X   |   X   o   o   o   o   o   o   o
                                            |                                                                           |
            b   x   x   x   x   x   x   i   |   i   x   x   x   x   x   x   b           o   o   o   o   o   o   o   X   |   X   o   o   o   o   o   o   o
                                            |                                                                           |
            b   x   x   x   x   x   x   i   |   i   x   x   x   x   x   x   b           o   o   o   o   o   o   o   X   |   X   o   o   o   o   o   o   o
                                            |                                                                           |
            b   x   x   x   x   x   x   i   |   i   x   x   x   x   x   x   b  ------>  o   o   o   o   o   o   o   X   |   X   o   o   o   o   o   o   o
                                            |                                                                           |
            b   x   x   x   x   x   x   i   |   i   x   x   x   x   x   x   b           o   o   o   o   o   o   o   X   |   X   o   o   o   o   o   o   o
                                            |                                                                           |
            b   x   x   x   x   x   x   i   |   i   x   x   x   x   x   x   b           o   o   o   o   o   o   o   X   |   X   o   o   o   o   o   o   o
                                            |                                                                           |
            v   i   i   i   i   i   i   v   |   v   b   b   b   b   b   b   v           X   X   X   X   X   X   X   o   |   o   o   o   o   o   o   o   o
                                            |                                                                           |
            ------------------------------------------------------------------                                           ------------------------------------------------------------------
            v   i   i   i   i   i   i   v   |                                           X   X   X   X   X   X   X   o   |
                                            |                                                                           |
            b   x   x   x   x   x   x   b   |                                           o   o   o   o   o   o   o   o   |
                                            |                                                                           |
            b   x   x   x   x   x   x   b   |                                           o   o   o   o   o   o   o   o   |
                                            |                                                                           |
            b   x   x   x   x   x   x   b   |                                           o   o   o   o   o   o   o   o   |
                                            |                                                                           |
            b   x   x   x   x   x   x   b   |                                           o   o   o   o   o   o   o   o   |
                                            |                                                                           |
            b   x   x   x   x   x   x   b   |                                           o   o   o   o   o   o   o   o   |
                                            |                                                                           |
            v   b   b   b   b   b   b   v   |                                           o   o   o   o   o   o   o   o   |

         */

        size_t tmp;
        m_size = tmp = 0;

        // number of interior basis functions (denoted by x)
        // (ONLY FOR TENSOR-PRODUCT STRUCTURE!!)
        for (size_t k=0; k!=m_patches.nPatches(); k++)
            tmp += (m_bases.basis(k).component(0).size()-2)*(m_bases.basis(k).component(1).size()-2);
        // gsDebug<<"Number of interior DoFs: "<<tmp<<"\n";
        m_size += tmp;

        // we don't add the DoFs on the interfaces (denoted by i)

        // // interfaces (denoted by i)
        gsBasis<T> * basis1;
        // gsBasis<T> * basis2;
        // gsVector<index_t> indices1,indices2;
        // tmp = 0;
        // for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
        // {
        //     basis1 = &m_bases.basis(iit->first().patch);
        //     basis2 = &m_bases.basis(iit->second().patch);
        //     tmp += basis1->boundary(iit->first().side()).size() - 2;
        //     tmp += basis2->boundary(iit->second().side()).size() - 2;
        // }
        // gsDebug<<"Number of interface DoFs: "<<tmp<<"\n";
        // m_size += tmp;

        // add the DoFs on the boundaries boundaries (denoted by b)
        tmp = 0;
        for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
        {
            basis1 = &m_bases.basis(bit->patch);
            tmp += (basis1->boundary(bit->side()).size() - 2);
        }
        // gsDebug<<"Number of boundary DoFs: "<<tmp<<"\n";
        m_size += tmp;

        // add DoFs for the vertices (denoted by v) if
        // - part of a boundary vertex with valence 1
        // - valence >2 (interior or boundary vertex) [add 3 in total]

        // vertices (denoted by v)
        tmp = 0;
        std::vector<bool> passed(m_patches.nPatches()*4);
        std::fill(passed.begin(), passed.end(), false);

        std::vector<patchCorner> corners;
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            for (index_t c=1; c<=4; c++)
            {
                index_t idx = _vertIndex(p,c);
                if (!passed.at(idx))
                {
                    // If we pass a corner, we get all connected corners and mark them so that we don't pass them again
                    m_patches.getCornerList(patchCorner(p,c),corners);
                    for (size_t k=0; k!=corners.size(); k++)
                        passed.at(_vertIndex(corners[k].patch,corners[k])) = true;

                    if (corners.size()==1) // valence = 1, must be a boundary vertex
                        tmp += 1;
                    else if (corners.size()>2 && corners.size()!=4)
                        tmp += 3;
                }
            }
        // gsDebug<<"Number of vertex DoFs: "<<tmp<<"\n";
        m_size += tmp;

        m_Bbases = m_bases = gsMultiBasis<T>(m_patches);

        m_mapModified = gsDofMapper(m_bases);
        m_mapOriginal = gsDofMapper(m_bases);

        // gsWriteParaview(m_patches,"mp",1000,true);
        // gsWriteParaview<>( m_bases.basis(0), "basis", 1000, true);

        m_matrix.resize(m_size,m_bases.totalSize());
        // m_matrix.resize(m_bases.totalSize(),m_bases.totalSize());

        m_coefs = this->allCoefficients();
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_computeMapper() // also initialize the mappers!
    {
        // interfaces
        std::fill(m_sideCheck.begin(), m_sideCheck.end(), false);
        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);

        gsBasis<T> * basis;

        std::vector<index_t> patches(2);
        std::vector<patchSide> psides(2);
        gsVector<index_t> indices;
        std::vector<patchCorner> pcorners;
        patchCorner pcorner;
        index_t cidx, sidx;
        std::pair<index_t,bool> vdata1, vdata2, vdata;

        // For the interfaces, we eliminate all DoFs located on the interface, except the ones coinciding with the end vertices
        for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
        {
            sidx = _sideIndex( iit->second().patch,iit->second().side());
            if (m_sideCheck.at(sidx))
                continue;

            patches[0] = iit->first().patch;
            patches[1] = iit->second().patch;
            psides[0] = patchSide(iit->first().patch,iit->first().side()); // the interface on the first patch
            psides[1] = patchSide(iit->second().patch,iit->second().side()); // the interface on the second patch

            for (index_t p = 0; p != 2; p++)
            {
                sidx = _sideIndex( patches[p] ,psides[p] );
                if (m_sideCheck.at(sidx))
                    continue;

                basis = &m_bases.basis(patches[p]);
                indices = static_cast<gsVector<index_t>>( basis->boundary(psides[p]) );

                patchSide(patches[p],psides[p]).getContainedCorners(d,pcorners);
                vdata1 = this->_vertexData(pcorners[0]);
                vdata2 = this->_vertexData(pcorners[1]);

                // cast indices to an std::vector
                std::vector<index_t> allIndices(indices.data(), indices.data() + indices.rows() * indices.cols());
                // for(size_t i=0; i < allIndices.size(); i++)
                //     std::cout << allIndices.at(i) << ' ';

                std::vector<index_t> selectedIndices;
                // for both vertices of the side, add the indices at the vertex and one inside
                for (index_t c =0; c!=2; c++)
                    selectedIndices.push_back(_indexFromVert(0,pcorners[c],psides[p],0,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.

                std::sort(selectedIndices.begin(),selectedIndices.end());
                // for(size_t i=0; i < selectedIndices.size(); i++)
                //     std::cout << selectedIndices.at(i) << ' ';

                std::vector<index_t> result(allIndices.size());
                std::vector<index_t>::iterator it=std::set_difference (allIndices.begin(), allIndices.end(), selectedIndices.begin(), selectedIndices.end(), result.begin());
                result.resize(it-result.begin());

                gsAsMatrix<index_t> indices(result,result.size(),1);
                m_mapModified.markBoundary(patches[p], indices);

                // gsDebug<<"Eliminated "<<indices.transpose()<<" of basis "<<patches[p]<<"\n";
                m_sideCheck.at(sidx) = true;
            }
        }

        // On the boundaries, we don't do anything
        for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
        {
            sidx = _sideIndex(bit->patch,bit->side());
            m_sideCheck.at(sidx) = true;
        }

        // For the vertices, we eliminate as follows (where v is the valence):
        // - No elimination when v==1
        // - One on each side when v==2
        // - All but three when v>2
        for (size_t p=0; p!=m_patches.nPatches(); p++)
        {
            for (index_t c=1; c<=4; c++)
            {
                cidx = _vertIndex(p,c);
                if (m_vertCheck.at(cidx))
                    continue;

                pcorner = patchCorner(p,c);

                basis = &m_bases.basis(p);

                // get valence and vertex info of corner
                std::pair<index_t,bool> vdata = _vertexData(pcorner); // corner c
                if (vdata.first==1) //valence = 1
                {  } // do nothing
                else if (vdata.first==2 || vdata.first==4) //valence = 2
                {
                    /*
                    v = 2
                        o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
                        o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
                        o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
                        o o o @ X |e| X @ o o o                 @: modified DoFs by interface rule
                        o o o * x |r| x * o o o                 *: modified DoFs by vertex rule (unique DoFs)
                        -----------------------                 %: modified DoFs by vertex rule (matched DoFs)
                        -boundary-| | -boundary-
                        -----------------------

                    v = 4
                        o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
                        o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
                        o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
                        @ @ @ * X |e| X * @ @ @                 @: modified DoFs by interface rule
                        X X X X x |r| x X X X X                 *: modified DoFs by vertex rule (unique DoFs)
                        -----------------------
                        -boundary-| |-interface
                        -----------------------
                        X X X X x |i| x X X X X
                        @ @ @ * X |n| X * @ @ @
                        o o o @ X |t| X @ o o o
                        o o o @ X |e| X @ o o o
                        o o o @ X |r| X @ o o o

                    */

                    // we mark the nodes belonging to the interface
                    pcorner.getContainingSides(d,psides);
                    for (size_t p=0; p!=psides.size(); p++)
                    {
                        if (m_patches.isInterface(psides[p]))
                        {
                            // the 0,0 vertex should be eliminated
                            m_mapModified.eliminateDof(basis->functionAtCorner(pcorner),pcorner.patch);
                        }
                    }
                }
                else if (vdata.first==3) //valence >2 but not 4
                {
                    // For the valence 3 case, there is not elimination needed (irrespective if it is a boundary vertex or an interior one)
                }
                else if (vdata.first>4)
                {
                    for (std::vector<patchCorner>::iterator it=pcorners.begin(); it!=pcorners.end(); it++)
                    {
                        // mark the vertex as passed
                        m_vertCheck[ _vertIndex(it->patch, it->corner()) ] = true;
                    }

                    _removeLowestCorners(pcorners,3);
                    for (std::vector<patchCorner>::iterator it=pcorners.begin(); it!=pcorners.end(); it++)
                    {
                        basis = &m_bases.basis(it->patch);
                        // Eliminate the 0,0 index
                        m_mapModified.eliminateDof(basis->functionAtCorner(*it),it->patch);
                    }
                }
                else
                    GISMO_ERROR("Something went terribly wrong");

                // label vertex as processed
                m_vertCheck[ cidx ] = true;
            }
        }
        m_mapModified.finalize();
        m_mapOriginal.finalize();

        m_matrix.resize( m_mapModified.freeSize(), m_mapOriginal.freeSize() );

        // gsDebugVar(m_mapModified.coupledSize());
        // gsDebugVar(m_mapModified.boundarySize());
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_handleVertex(patchCorner pcorner)
    {

        if (m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ])
        {
            // gsDebug<<"corner "<<pcorner.corner()<<" ("<<pcorner.patch<<") skipped!\n";
            return;
        }
        // if (std::count(m_C0s.begin(), m_C0s.end(), pcorner))
        // {
        //     gsInfo<<"patch = "<<pcorner.patch<<", corner = "<<pcorner.corner()<<"\n";
        //     m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
        //     return;
        // }

        std::vector<patchSide> psides;
        std::vector<patchCorner> corners;
        std::vector<index_t> indices;
        std::vector<index_t> weights;

        index_t colIdx, rowIdx, cornerIdx;

        std::pair<index_t,bool> vdata = _vertexData(pcorner); // corner c
        pcorner.getContainingSides(d,psides);

        gsBasis<T> * basis;


        if (!vdata.second) // boundary vertices
        {
            // Correct?
            if (vdata.first==1)
            {
                // Only do the coefficient for the 0,0 DoF
                GISMO_ASSERT(psides.size()==2,"Must have 2 adjacent sides");
                basis = &m_bases.basis(pcorner.patch);
                cornerIdx = basis->functionAtCorner(pcorner);
                rowIdx = m_mapModified.index(cornerIdx,pcorner.patch);
                colIdx = m_mapOriginal.index(cornerIdx,pcorner.patch);
                m_matrix(rowIdx,colIdx) = 1.0;

                m_basisCheck[rowIdx] = true;
                m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
                // gsInfo<<"patch = "<<pcorner.patch<<", corner = "<<pcorner.corner()<<"\n";
                return;
            }
            // Correct?
            else if (vdata.first==2)
            {
                // 1. find the interface
                index_t iindex = m_patches.isInterface(psides[0]) ? 0 : 1;

                boundaryInterface iface;
                GISMO_ASSERT(m_patches.getInterface(psides[iindex],iface),"Must be an interface");
                // 2. collect indices
                // If we want C0 at this vertex, we only handle the row k=1.
                indices.resize(3);

                patchSide otherSide = iface.other(psides[iindex]);
                patchCorner otherCorner = iface.mapCorner(pcorner);
                indices[0] = _indexFromVert(0,pcorner,psides[iindex],1); // 1,0 on patch of iface (offset 1)
                basis = &m_bases.basis(pcorner.patch);
                indices[1] = basis->functionAtCorner(pcorner);
                basis = &m_bases.basis(otherCorner.patch);
                indices[2] = basis->functionAtCorner(otherCorner);

                rowIdx = m_mapModified.index(indices[0],pcorner.patch);
                colIdx = m_mapOriginal.index(indices[0],pcorner.patch);
                m_matrix(rowIdx,colIdx) = 1.0;
                colIdx = m_mapOriginal.index(indices[1],psides[iindex].patch);
                m_matrix(rowIdx,colIdx) = 0.5;
                colIdx = m_mapOriginal.index(indices[2],otherSide.patch);
                m_matrix(rowIdx,colIdx) = 0.5;

                m_basisCheck[rowIdx] = true;

                m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
            }
            // Correct?
            else if (vdata.first>2)
            {
                // 2. make container for the interfaces
                std::vector<boundaryInterface> ifaces;
                boundaryInterface iface;
                std::vector<index_t> rowIndices, colIndices, patchIndices;

                index_t tmpIdx;

                // pcorner is the current corner
                m_patches.getCornerList(pcorner,corners);

                // The 1,1 corners of each patch will be given 0.5 weight in the interface handling, but in addition they will receive a 1/(v+2) weight from the 0,0 DoFs on each patch
                pcorner.getContainingSides(d,psides);

                for (typename std::vector<patchCorner>::iterator it = corners.begin(); it!=corners.end(); it++)
                {
                    basis = &m_bases.basis(it->patch);
                    colIndices.push_back(basis->functionAtCorner(*it));
                    patchIndices.push_back(it->patch);
                }

                basis = &m_bases.basis(pcorner.patch);
                tmpIdx = _indexFromVert(1,pcorner,psides[0],1); // 1,1 corner
                rowIndices.push_back(tmpIdx);
                // Check if one of the adjacent interfaces is a boundary; if so, add weight 1.0 to itself and add it to the rowIndices
                for (index_t k = 0; k!=2; k++)
                    if (!m_patches.getInterface(psides[k],iface)) // check if the side is NOT an interface
                    {
                        tmpIdx = _indexFromVert(1,pcorner,psides[k],0); // 1,0 corner (on the boundary)
                        rowIdx = m_mapModified.index(tmpIdx,pcorner.patch);
                        colIdx = m_mapOriginal.index(tmpIdx,pcorner.patch);
                        m_matrix(rowIdx,colIdx) = 1.0;
                        rowIndices.push_back(tmpIdx);
                    }


                for (size_t k=0; k!=colIndices.size(); k++)
                    for (std::vector<index_t>::iterator rit=rowIndices.begin(); rit!=rowIndices.end(); rit++ )
                    {
                        rowIdx = m_mapModified.index(*rit,pcorner.patch);
                        colIdx = m_mapOriginal.index(colIndices.at(k),patchIndices.at(k));
                        m_matrix(rowIdx,colIdx) = 1. / (vdata.first + 2);
                        m_basisCheck[rowIdx] = true;
                    }

                // Furthermore, if the corner is one of the three DoFs that is preserved, we mark the 0,0 DoF as handled (should be a zero-row)
                basis = &m_bases.basis(pcorner.patch);
                tmpIdx = basis->functionAtCorner(pcorner.corner());
                rowIdx = m_mapModified.index(tmpIdx,pcorner.patch);
                m_basisCheck[rowIdx] = true;

                m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
            }
        } // end boundary vertices
        else // interior vertices
        {
            pcorner.getContainingSides(d,psides);
            index_t b11_p1 = _indexFromVert(1,pcorner,psides[0],1); // point 1,1 (does not matter which reference side is taken)
            index_t rowIdx = m_mapModified.index(b11_p1,pcorner.patch);

            m_patches.getCornerList(pcorner,corners);

            // The 1,1 corners of each patch will be given 0.5 weight in the interface handling, but in addition they will receive a 1/v weight from the 0,0 DoFs on each patch
            gsBasis<T> * tmpBasis;
            index_t index;
            for (std::vector<patchCorner>::iterator corn = corners.begin(); corn != corners.end(); ++corn)
            {
                tmpBasis = &m_bases.basis(corn->patch);
                index = tmpBasis->functionAtCorner(corn->corner());
                colIdx = m_mapOriginal.index(index,corn->patch);
                m_matrix(rowIdx,colIdx) = 1./vdata.first;
            }

            m_basisCheck[rowIdx] = true;

            // Furthermore, if the corner is one of the three DoFs that is preserved, we mark the 0,0 DoF as handled (should be a zero-row)
            tmpBasis = &m_bases.basis(pcorner.patch);
            index = tmpBasis->functionAtCorner(pcorner.corner());
            gsDebugVar(index);
            gsDebugVar(pcorner.patch);
            gsDebugVar(m_mapModified.index(index,pcorner.patch));
            gsDebugVar(m_mapModified.is_free(index,pcorner.patch));
            if (m_mapModified.is_free(index,pcorner.patch))
            {
                rowIdx = m_mapModified.index(index,pcorner.patch);
                m_basisCheck[rowIdx] = true;
            }

            // Mark the vertex as processed
            m_vertCheck[ _vertIndex(pcorner.patch, pcorner.corner()) ] = true;
        }
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_handleInterface(boundaryInterface iface)
    {
        if (m_sideCheck[ _sideIndex(iface.first().patch, iface.first().side()) ] || m_sideCheck[ _sideIndex(iface.second().patch, iface.second().side()) ])
        {
            // gsDebug<<"sides "<<iface.first().side()<<" ("<<iface.first().patch<<") and "<<iface.second().side()<<" ("<<iface.second().patch<<") skipped!\n";
            return;
        }

        std::vector<patchCorner> pcorners;

        std::vector<std::vector<index_t>> selectedIndices(2);
        std::vector<std::vector<index_t>> selectedOIndices(2);

        std::vector<gsBasis<T> *> basis(2);
        std::vector<gsMatrix<index_t>> indices(2); // interface indices
        std::vector<gsMatrix<index_t>> oindices(2); // interface indices
        gsVector<bool> dirOr;

        index_t np;
        // for both vertices of the side, add the indices at the vertex and one inside
        for (index_t p =0; p!=2; p++)
        {
            np = abs(p-1); // not index p;
            // get the bases belonging to both patches
            basis[p] = &m_bases.basis(iface[p].patch);
            // now we treat the offset of the interface
            indices[p] = basis[p]->boundaryOffset(iface[p].side(),0);
            oindices[p] = basis[p]->boundaryOffset(iface[p].side(),1);

            iface[p].getContainedCorners(d,pcorners);
            for (index_t c =0; c!=2; c++)
            {
                selectedIndices[p].push_back(_indexFromVert(0,pcorners[c],iface[p],0,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
                selectedOIndices[p].push_back(_indexFromVert(0,pcorners[c],iface[p],1,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
            }

            std::sort(selectedIndices[p].begin(),selectedIndices[p].end());
            std::sort(selectedOIndices[p].begin(),selectedOIndices[p].end());
            std::vector<index_t> allIndices(indices[p].data(), indices[p].data() + indices[p].rows() * indices[p].cols());
            std::vector<index_t> result(allIndices.size());
            std::vector<index_t>::iterator it=std::set_difference (allIndices.begin(), allIndices.end(), selectedIndices[p].begin(), selectedIndices[p].end(), result.begin());
            result.resize(it-result.begin());
            indices[p] = gsAsMatrix<index_t>(result);

            std::vector<index_t> allOIndices(oindices[p].data(), oindices[p].data() + oindices[p].rows() * oindices[p].cols());
            result.resize(allOIndices.size());
            std::vector<index_t>::iterator ito=std::set_difference (allOIndices.begin(), allOIndices.end(), selectedOIndices[p].begin(), selectedOIndices[p].end(), result.begin());
            result.resize(ito-result.begin());
            oindices[p] = gsAsMatrix<index_t>(result);
        }
        // Flip the index vector if the directions of the indices do not match
        dirOr = iface.dirOrientation();
        for (short_t k = 0; k<m_patches.dim(); ++k )
        {
            if ( k == iface[0].side().direction() ) // skip ?
                continue;

            if ( ! dirOr[k] ) // flip ?
            {
                // gsDebug<<"\t\tReversed direction\n";
                indices[0].reverseInPlace();
                oindices[0].reverseInPlace();
            }
        }

        index_t rowIdx,colIdx;
        // loop over adjacent patches and couple the DoFs.
        for (index_t p =0; p!= 2; p++)
        {
            np = abs(p-1); // not index p;
            for (index_t k=0; k!= indices[p].size(); k++ )
            {
                rowIdx = m_mapModified.index(oindices[p].at(k),iface[p].patch);
                // rowIdx1 = m_mapOriginal.index(oindices[p].at(k),patches[p]);
                colIdx = m_mapOriginal.index(oindices[p].at(k),iface[p].patch);
                m_matrix(rowIdx,colIdx) = 1.0;

                colIdx = m_mapOriginal.index(indices[p].at(k),iface[p].patch);
                m_matrix(rowIdx,colIdx) = 0.5;

                colIdx = m_mapOriginal.index(indices[np].at(k),iface[np].patch);
                m_matrix(rowIdx,colIdx) = 0.5;

                m_basisCheck[rowIdx] = true;
            }
            m_sideCheck[ _sideIndex(iface[p].patch, iface[p].side()) ] = true; // side finished
        }

    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_handleBoundary(patchSide side)
    {
            std::vector<patchCorner> pcorners;
            std::vector<index_t> selectedIndices;
            gsBasis<T> * basis = &m_bases.basis(side.patch);
            gsMatrix<index_t> indices = basis->boundaryOffset(side.side(),0);
            side.getContainedCorners(d,pcorners);
            for (index_t c =0; c!=2; c++)
            {
                selectedIndices.push_back(_indexFromVert(0,pcorners[c],side,0,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
                selectedIndices.push_back(_indexFromVert(1,pcorners[c],side,0,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
            }

            std::sort(selectedIndices.begin(),selectedIndices.end());
            std::vector<index_t> allIndices(indices.data(), indices.data() + indices.rows() * indices.cols());
            std::vector<index_t> result(allIndices.size());
            std::vector<index_t>::iterator it=std::set_difference (allIndices.begin(), allIndices.end(), selectedIndices.begin(), selectedIndices.end(), result.begin());
            result.resize(it-result.begin());

            index_t rowIdx,colIdx;
            for (std::vector<index_t>::iterator it = result.begin(); it!=result.end(); ++it)
            {
                rowIdx = m_mapModified.index(*it,side.patch);
                colIdx = m_mapOriginal.index(*it,side.patch);
                m_matrix(rowIdx,colIdx) = 1.0;
                m_basisCheck[rowIdx] = true;
            }

            m_sideCheck.at( _sideIndex(side.patch,side.side()) ) = true;
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_handleInterior()
    {
        index_t rowIdx,colIdx;
        for(size_t p=0; p!=m_patches.nPatches(); p++)
            for (index_t b=0; b!=m_bases.basis(p).size(); b++)
            {
                rowIdx = m_mapModified.index(b,p);
                // rowIdx = m_mapOriginal.index(b,p);
                if ( (!m_mapModified.is_free(b,p)) || (m_basisCheck[rowIdx]) )
                // if ( (m_basisCheck[rowIdx]) )
                    continue;
                colIdx = m_mapOriginal.index(b,p);
                m_matrix(rowIdx,colIdx) = 1;
                m_basisCheck[rowIdx] = true;
                // gsInfo<<"Basis function "<<rowIdx<<"(patch: "<<p<<"; fun: "<<b<<") is "<< (m_basisCheck[rowIdx] ? "" : "not ")<<"processed\n";
            }
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_whichHandled()
    {
        for (size_t p=0; p!=m_patches.nPatches(); p++)
        {
            gsInfo<<"----------------------------------\n";
            gsInfo<<"patch "<<p<<"\n";
            for (size_t b=0; b!=m_mapOriginal.patchSize(p); b++)
            {
                index_t idx = m_mapModified.index(b,p);
                if (m_mapModified.is_free_index(idx))
                    gsInfo<<"basis function "<<b<<" on patch "<<p<<" is "<<(m_basisCheck[idx] ? "":"not ")<<"handled\n";
                else
                    gsInfo<<"basis function "<<b<<" on patch "<<p<<" is "<<"eliminated\n";

            }
        }
    }

    template<short_t d,class T>
    void gsAlmostC1<d,T>::_computeSmoothMatrix()
    {

        m_basisCheck.resize(m_size);
        // m_basisCheck.resize(m_bases.totalSize());
        std::fill(m_basisCheck.begin(), m_basisCheck.end(), false);

        std::fill(m_sideCheck.begin(), m_sideCheck.end(), false);
        std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);

        // iterate over the vertices
        for (size_t p=0; p!=m_patches.nPatches(); p++)
            for (index_t c=1; c<=4; c++)
                _handleVertex(patchCorner(p,c));

        for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
            _handleInterface(*iit);

        // boundaries
        for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
            _handleBoundary(*bit);

        _handleInterior();

        if (m_verbose) { _whichHandled(); }

        bool checkSides = std::all_of(m_sideCheck.begin(), m_sideCheck.end(), [](bool m_sideCheck) { return m_sideCheck; });
        GISMO_ASSERT(checkSides,"Not all sides are checked");
        bool checkVerts = std::all_of(m_vertCheck.begin(), m_vertCheck.end(), [](bool m_vertCheck) { return m_vertCheck; });
        GISMO_ASSERT(checkVerts,"Not all vertices are checked");
        bool checkBasis = std::all_of(m_basisCheck.begin(), m_basisCheck.end(), [](bool m_basisCheck) { return m_basisCheck; });
        GISMO_ASSERT(checkBasis,"Not all vertices are checked");
    }

    // THIS FUNCTION IS LEGACY!!
    // template<short_t d,class T>
    // void gsAlmostC1<d,T>::_computeSmoothMatrix2()
    // {
    //     m_basisCheck.resize(m_size);
    //     // m_basisCheck.resize(m_bases.totalSize());
    //     std::fill(m_basisCheck.begin(), m_basisCheck.end(), false);

    //     std::fill(m_sideCheck.begin(), m_sideCheck.end(), false);
    //     std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);

    //     //

    //     // std::vector<std::vector<patchCorner> > cornerLists;
    //     // m_patches.getEVs(cornerLists);
    //     // if (cornerLists.size()!=0)
    //     //     GISMO_ERROR("Refinement procedure needed since mesh contains an extraordinary vertex! This is currently not implemented....");

    //     m_matrix.resize( m_mapModified.freeSize(), m_bases.totalSize() );
    //     // m_matrix.resize( m_bases.totalSize(), m_bases.totalSize() );

    //     std::vector<gsBasis<T> *> basis(2);
    //     std::vector<gsMatrix<index_t>> indices(2); // interface indices
    //     std::vector<gsMatrix<index_t>> oindices(2); // interface indices (with offset 1)
    //     gsMatrix<index_t> tmpIndices;
    //     gsVector<bool> dirOr;

    //     boundaryInterface iface;
    //     std::vector<patchCorner> pcorners;
    //     patchCorner pcorner;
    //     std::vector<boxSide> bsides;
    //     std::vector<patchSide> psides(2);
    //     std::vector<index_t> patches(2);
    //     std::vector<index_t> minJ(2);
    //     std::vector<index_t> maxJ(2);
    //     index_t colIdx,rowIdx,rowIdx1,rowIdx2;

    //     for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
    //     {
    //         sideInfo(iit->first());
    //         sideInfo(iit->second());

    //         // get the bases belonging to both patches
    //         basis[0] = &m_bases.basis(iit->first().patch);
    //         basis[1] = &m_bases.basis(iit->second().patch);

    //         // get a boundaryInterface object for the interface and obtain the matched indices along this interface
    //         // m_patches.getInterface(patchSide(iit->first().patch,iit->first().side()),iface);
    //         // basis[0]->matchWith(iface,*basis[1],indices[0],indices[1]);

    //         // now we treat the offset of the interface
    //         indices[0] = basis[0]->boundaryOffset(iit->first().side(),0);
    //         indices[1] = basis[1]->boundaryOffset(iit->second().side(),0);
    //         oindices[0] = basis[0]->boundaryOffset(iit->first().side(),1);
    //         oindices[1] = basis[1]->boundaryOffset(iit->second().side(),1);

    //         // store relevant data in vectors with size 2 for the loop over the associated corners
    //         patches[0] = iit->first().patch;
    //         patches[1] = iit->second().patch;
    //         psides[0] = patchSide(iit->first().patch,iit->first().side()); // the interface on the first patch
    //         psides[1] = patchSide(iit->second().patch,iit->second().side()); // the interface on the second patch

    //         index_t idx, np;
    //         for (index_t p =0; p!= 2; p++)
    //         {
    //             psides[p].getContainedCorners(d,pcorners); //corners on first patch
    //             np = abs(p-1); // not index p;
    //             for (size_t c=0; c!=pcorners.size(); c++)
    //             {
    //                 idx = _vertIndex(patches[p],pcorners[c]) ;
    //                 // if(m_vertCheck[ idx] )
    //                 //     continue;

    //                 // get valence and vertex info of corner
    //                 std::pair<index_t,bool> vdata = _vertexData(pcorners[c]); // corner c
    //                 if (!vdata.second) // boundary vertex
    //                 {
    //                     if (vdata.first==1) //valence = 1
    //                         gsInfo<<"this case should not exist... (interface handling)\n";
    //                     else if (vdata.first==2) //valence = 2
    //                     {
    //                         /*
    //                             o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
    //                             o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
    //                             o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
    //                             o o o * x |e| x * o o o                 @: modified DoFs by interface rule
    //                             o o o * x |r| x * o o o                 *: modified DoFs by vertex rule (unique DoFs)
    //                             -----------------------                 %: modified DoFs by vertex rule (matched DoFs)
    //                             -boundary-| | -boundary-
    //                             -----------------------
    //                         */

    //                         /*
    //                             NOTE: the alternative would be to adjust minJ and maxJ
    //                         */
    //                         std::vector<patchCorner> otherCorners;
    //                         m_patches.getCornerList(pcorners[c],otherCorners);
    //                         patchCorner otherCorner = otherCorners[1]; // because the first entry is pcorner[c]
    //                         for (index_t i=0; i!=2; i++)
    //                         {
    //                             // get indices
    //                             index_t bj0_p1 = _indexFromVert(i,pcorners[c],psides[p],0);
    //                             index_t bj1_p1 = _indexFromVert(i,pcorners[c],psides[p],1);
    //                             index_t bj0_p2 = _indexFromVert(i,otherCorner,psides[np],0);
    //                             // index_t bj1_p2 = _indexFromVert(i,otherCorner,psides[np],1);



    //                             rowIdx1 = m_mapModified.index(bj1_p1,patches[p]);
    //                             // rowIdx1 = m_mapOriginal.index(bj1_p1,patches[p]);

    //                             // Factor 1 to itself
    //                             colIdx = m_mapOriginal.index(bj1_p1,patches[p]);

    //                             m_matrix(rowIdx1,colIdx) = 1;

    //                             // Factor 0.5 to eachother
    //                             colIdx = m_mapOriginal.index(bj0_p1,patches[p]);
    //                             m_matrix(rowIdx1,colIdx) = 0.5;
    //                             colIdx = m_mapOriginal.index(bj0_p2,patches[np]);
    //                             m_matrix(rowIdx1,colIdx) = 0.5;

    //                             m_basisCheck[rowIdx1] = true;
    //                         }
    //                     }
    //                     else if (vdata.first==3) //valence = 3
    //                     {
    //                         /*
    //                             o o o @ X | X @ o o o               x: eliminated DoFs by vertex rule
    //                             o o o @ X | X @ o o o               X: eliminated DoFs by interface/boundary rule
    //                             o o o @ X | X @ o o o               o: preserved DoFs (interior)
    //                             @ @ @ * x | x * @ @ @               @: modified DoFs by interface rule
    //                             X X X x % | % x X X X               *: modified DoFs by vertex rule (unique DoFs)
    //                             ----------|----------               %: modified DoFs by vertex rule (matched DoFs)
    //                                         | % x X X X
    //                                         | x * @ @ @
    //                                         | X @ o o o
    //                                         | X @ o o o
    //                                         | X @ o o o

    //                         */
    //                         // the (0,0) basis functions (relative to the corner) of all patches have weight 1/4
    //                         gsBasis<T> * tmpBasis;
    //                         std::vector<patchCorner> otherCorners;
    //                         m_patches.getCornerList(pcorners[c],otherCorners);

    //                         // find the patch corner which shares the interface
    //                         patchCorner otherCorner;
    //                         if (otherCorners[1].patch == patches[np])
    //                             otherCorner = otherCorners[1];
    //                         else if (otherCorners[2].patch == patches[np])
    //                             otherCorner = otherCorners[2];
    //                         else
    //                             GISMO_ERROR("HUH?"<<otherCorners[0].patch<<"; "<<otherCorners[1].patch<<"; "<<otherCorners[2].patch<<"; "<<patches[p]);

    //                         // get the sides corresponding to the corner on patch p (pcorners[c])
    //                         std::vector<patchSide> csides;
    //                         pcorners[c].getContainingSides(d,csides);

    //                         index_t b01_p1 = _indexFromVert(1,pcorners[c],csides[0],0); // index from vertex pcorners[c] along side psides[0] with offset 0.
    //                         index_t b10_p1 = _indexFromVert(1,pcorners[c],csides[1],0); // index from vertex pcorners[c] along side psides[0] with offset 0.
    //                         index_t b11_p1 = _indexFromVert(1,pcorners[c],psides[p],1); // point 1,1
    //                         index_t b10_p2 = _indexFromVert(1,otherCorner,psides[np],0); // point 0,1

    //                         rowIdx1 = m_mapModified.index(b11_p1,patches[p]);
    //                         // rowIdx1 = m_mapOriginal.index(b11_p1,patches[p]);

    //                         // 1. give the basis function a weight 1 from itself
    //                         colIdx = m_mapOriginal.index(b11_p1,patches[p]);
    //                         m_matrix(rowIdx1,colIdx) = 1;

    //                         // 2. give the basis function a weight 1/4 from the other (0,0) basis functions
    //                         index_t index;
    //                         for (size_t k =0; k!=otherCorners.size(); k++)
    //                         {
    //                             tmpBasis = &m_bases.basis(otherCorners[k].patch);
    //                             index = tmpBasis->functionAtCorner(otherCorners[k].corner());
    //                             colIdx = m_mapOriginal.index(index,otherCorners[k].patch);
    //                             m_matrix(rowIdx1,colIdx) = 0.25;
    //                         }

    //                         // 3. add weight 1/2 from the interface and boundary bcs.
    //                         colIdx = m_mapOriginal.index(b10_p1,patches[p]);
    //                         m_matrix(rowIdx1,colIdx) = 0.5;

    //                         colIdx = m_mapOriginal.index(b01_p1,patches[p]);
    //                         m_matrix(rowIdx1,colIdx) = 0.5;

    //                         colIdx = m_mapOriginal.index(b10_p2,patches[np]);
    //                         m_matrix(rowIdx1,colIdx) = 0.5;


    //                         // 4. Weigh the matched - bonus - DoF
    //                         // This DoF is assigned to the (0,0) basis function of the patch with the lowest number

    //                         // find lowest patch index
    //                         patchCorner pcorn = pcorners[c];
    //                         for (size_t k=0; k!=otherCorners.size(); k++)
    //                             pcorn = (otherCorners[k].patch < pcorn.patch ) ? otherCorners[k] : pcorn;

    //                         // the (0,0) corners weigh 1/4
    //                         tmpBasis = &m_bases.basis(pcorn.patch);
    //                         index_t bindex = tmpBasis->functionAtCorner(pcorn.corner());
    //                         rowIdx2 = m_mapModified.index(bindex,pcorn.patch); // all corners point to this index.
    //                         // rowIdx2 = m_mapOriginal.index(bindex,pcorn.patch); // all corners point to this index.

    //                         // tmpBasis = &m_bases.basis(pcorners[c].patch);
    //                         // index_t bindex = tmpBasis->functionAtCorner(pcorners[c].corner());
    //                         // rowIdx2 = m_mapOriginal.index(bindex,pcorners[c].patch); // all corners point to this index.
    //                         // // rowIdx2 = m_mapModified.index(bindex,pcorners[c].patch); // all corners point to this index.
    //                         for (std::vector<patchCorner>::iterator corner = otherCorners.begin(); corner != otherCorners.end(); corner++)
    //                         {
    //                             tmpBasis = &m_bases.basis(corner->patch);
    //                             bindex = tmpBasis->functionAtCorner(corner->corner());
    //                             colIdx = m_mapOriginal.index(bindex,corner->patch);
    //                             m_matrix(rowIdx2,colIdx) = 0.25;
    //                         }

    //                         // determine if both sides are interfaces
    //                         std::vector<patchSide> sides;
    //                         pcorners[c].getContainingSides(d,sides);

    //                         GISMO_ASSERT(sides.size()==2,"There must be 2 adjacent sides, instead of "<<sides.size());
    //                         for (index_t k = 0; k!=2; k++)
    //                             if (m_patches.isBoundary(sides[k]))
    //                             {
    //                                 index_t b01 = _indexFromVert(1,pcorners[c],sides[k],0);
    //                                 colIdx = m_mapOriginal.index(b01,patches[p]);
    //                                 m_matrix(rowIdx2,colIdx) = 0.5;
    //                             }

    //                         // Mark the functions as processed.
    //                         m_basisCheck[rowIdx1] = true;
    //                         // 0,0 points
    //                         for (std::vector<patchCorner>::iterator corner = otherCorners.begin(); corner != otherCorners.end(); corner++)
    //                         {
    //                             tmpBasis = &m_bases.basis(corner->patch);
    //                             bindex = tmpBasis->functionAtCorner(corner->corner());
    //                             rowIdx1 = m_mapModified.index(bindex,corner->patch);
    //                             // rowIdx1 = m_mapOriginal.index(bindex,corner->patch);
    //                             m_basisCheck[rowIdx1] = true;
    //                         }

    //                     }
    //                     else
    //                         GISMO_ERROR("Boundary vertex with valence = "<<vdata.first<<" has no implementation");
    //                 }
    //                 else // interior vertex
    //                 {
    //                     /*
    //                         o o o @ X | X @ o o o               x: eliminated DoFs by vertex rule
    //                         o o o @ X | X @ o o o               X: eliminated DoFs by interface/boundary rule
    //                         o o o @ X | X @ o o o               o: preserved DoFs (interior)
    //                         @ @ @ * x | x * @ @ @               @: modified DoFs by interface rule
    //                         X X X x % | % x X X X               *: modified DoFs by vertex rule (unique DoFs)
    //                         ----------|----------               %: modified DoFs by vertex rule (matched DoFs)
    //                         X X X x % | % x X X X
    //                         @ @ @ * x | x * @ @ @
    //                         o o o @ x | X @ o o o
    //                         o o o @ x | X @ o o o
    //                         o o o @ x | X @ o o o

    //                     */
    //                     gsBasis<T> * tmpBasis;
    //                     std::vector<patchCorner> otherCorners;
    //                     m_patches.getCornerList(pcorners[c],otherCorners);
    //                     // find the patch corner which shares the interface
    //                     patchCorner otherCorner;
    //                     for (std::vector<patchCorner>::iterator corner = otherCorners.begin(); corner != otherCorners.end(); corner++)
    //                         if (corner->patch == patches[np])
    //                             otherCorner = *corner;

    //                     index_t b10_p1 = _indexFromVert(1,pcorners[c],psides[p],0); // index from vertex pcorners[c] along side psides[0] with offset 0.
    //                     index_t b11_p1 = _indexFromVert(1,pcorners[c],psides[p],1); // point 1,1
    //                     index_t b10_p2 = _indexFromVert(1,otherCorner,psides[np],0); // point 0,1

    //                     rowIdx1 = m_mapModified.index(b11_p1,patches[p]);
    //                     // rowIdx1 = m_mapOriginal.index(b11_p1,patches[p]);

    //                     // 1. give the basis function a weight 1 from itself
    //                     colIdx = m_mapOriginal.index(b11_p1,patches[p]);
    //                     m_matrix(rowIdx1,colIdx) = 1;

    //                     // 2. give the basis function a weight 1/4 from the other (0,0) basis functions
    //                     index_t index;
    //                     for (size_t k =0; k!=otherCorners.size(); k++)
    //                     {
    //                         tmpBasis = &m_bases.basis(otherCorners[k].patch);
    //                         index = tmpBasis->functionAtCorner(otherCorners[k].corner());
    //                         colIdx = m_mapOriginal.index(index,otherCorners[k].patch);
    //                         m_matrix(rowIdx1,colIdx) = 1./vdata.first;
    //                     }

    //                     // 3. add weight 1/2 from the interface bcs.
    //                     colIdx = m_mapOriginal.index(b10_p1,patches[p]);
    //                     m_matrix(rowIdx1,colIdx) = 0.5;

    //                     colIdx = m_mapOriginal.index(b10_p2,patches[np]);
    //                     m_matrix(rowIdx1,colIdx) = 0.5;

    //                     m_basisCheck[rowIdx1] = true;
    //                 } // end vertex if

    //                 // label vertex as processed
    //                 m_vertCheck[ idx ] = true;
    //                 m_basisCheck[rowIdx1] = true;
    //             } // end corner loop

    //             // label side as processed
    //             // nothing
    //         } // end patch loop

    //         /*
    //             Now we treat the interior degrees of freedom. For all interface DoFs, a linear combination of basis functions is constructed and the DoFs located on the interface are removed from the mapper.
    //         */

    //         std::vector<std::vector<index_t>> selectedIndices(2);
    //         std::vector<std::vector<index_t>> selectedOIndices(2);
    //         // for both vertices of the side, add the indices at the vertex and one inside
    //         for (index_t p =0; p!=2; p++)
    //         {
    //             psides[p].getContainedCorners(d,pcorners);
    //             for (index_t c =0; c!=2; c++)
    //             {
    //                 selectedIndices[p].push_back(_indexFromVert(0,pcorners[c],psides[p],0,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
    //                 selectedIndices[p].push_back(_indexFromVert(1,pcorners[c],psides[p],0,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.

    //                 selectedOIndices[p].push_back(_indexFromVert(0,pcorners[c],psides[p],1,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
    //                 selectedOIndices[p].push_back(_indexFromVert(1,pcorners[c],psides[p],1,0)); // index from vertex pcorners[c] along side psides[0] with offset 0.
    //             }

    //             std::sort(selectedIndices[p].begin(),selectedIndices[p].end());
    //             std::sort(selectedOIndices[p].begin(),selectedOIndices[p].end());
    //             std::vector<index_t> allIndices(indices[p].data(), indices[p].data() + indices[p].rows() * indices[p].cols());
    //             std::vector<index_t> result(allIndices.size());
    //             std::vector<index_t>::iterator it=std::set_difference (allIndices.begin(), allIndices.end(), selectedIndices[p].begin(), selectedIndices[p].end(), result.begin());
    //             result.resize(it-result.begin());
    //             indices[p] = gsAsMatrix<index_t>(result);

    //             std::vector<index_t> allOIndices(oindices[p].data(), oindices[p].data() + oindices[p].rows() * oindices[p].cols());
    //             result.resize(allOIndices.size());
    //             std::vector<index_t>::iterator ito=std::set_difference (allOIndices.begin(), allOIndices.end(), selectedOIndices[p].begin(), selectedOIndices[p].end(), result.begin());
    //             result.resize(ito-result.begin());
    //             oindices[p] = gsAsMatrix<index_t>(result);
    //         }
    //         // Flip the index vector if the directions of the indices do not match
    //         dirOr = iit->dirOrientation();
    //         for (short_t k = 0; k<m_patches.dim(); ++k )
    //         {
    //             if ( k == iit->first().side().direction() ) // skip ?
    //                 continue;

    //             if ( ! dirOr[k] ) // flip ?
    //             {
    //                 // gsDebug<<"\t\tReversed direction\n";
    //                 indices[0].reverseInPlace();
    //                 oindices[0].reverseInPlace();
    //             }
    //         }

    //         // loop over adjacent patches and couple the DoFs.
    //         for (index_t p =0; p!= 2; p++)
    //         {
    //             np = abs(p-1); // not index p;
    //             for (index_t k=0; k!= indices[p].size(); k++ )
    //             {
    //                 rowIdx1 = m_mapModified.index(oindices[p].at(k),patches[p]);
    //                 // rowIdx1 = m_mapOriginal.index(oindices[p].at(k),patches[p]);
    //                 colIdx = m_mapOriginal.index(oindices[p].at(k),patches[p]);
    //                 m_matrix(rowIdx1,colIdx) = 1.0;

    //                 colIdx = m_mapOriginal.index(indices[p].at(k),patches[p]);
    //                 m_matrix(rowIdx1,colIdx) = 0.5;

    //                 colIdx = m_mapOriginal.index(indices[np].at(k),patches[np]);
    //                 m_matrix(rowIdx1,colIdx) = 0.5;

    //                 m_basisCheck[rowIdx1] = true;
    //             }
    //             m_sideCheck[ _sideIndex(patches[p], psides[p]) ] = true; // side finished
    //         }

    //     } // end interface iterator

    //     // boundaries
    //     for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
    //     {
    //         // Mark sides as done
    //         m_sideCheck[ _sideIndex(bit->patch, bit->side()) ] = true;
    //     }

    //     patchCorner corner;
    //     // iterate over the other vertices
    //     for (size_t p=0; p!=m_patches.nPatches(); p++)
    //         for (index_t c=1; c<=4; c++)
    //         {
    //             if (m_vertCheck[ _vertIndex(p,c) ]) // vertex already covered?
    //                 continue;

    //             std::pair<index_t,bool> vdata = _vertexData(patchCorner(p,c)); // corner c
    //             if (!vdata.second) // boundary vertex
    //             {
    //                 if (vdata.first!=1) //valence = 1
    //                     gsInfo<<"this case should not exist...\n";
    //             }
    //             else
    //                 gsInfo<<"this case should not exist...\n";

    //             // label vertex as processed
    //             m_vertCheck[ _vertIndex(p,c) ] = true;
    //         }

    //     // Mark all basis functions with a row filled in the m_matrix as TRUE (i.e. passed)
    //     for (int i = 0; i < m_matrix.outerSize(); ++i)
    //         for (typename gsSparseMatrix<T>::InnerIterator it(m_matrix,i); it; ++it)
    //             m_basisCheck[it.row()] = true;

    //     for(size_t p=0; p!=m_patches.nPatches(); p++)
    //         for (index_t b=0; b!=m_bases.basis(p).size(); b++)
    //         {
    //             rowIdx = m_mapModified.index(b,p);
    //             // rowIdx = m_mapOriginal.index(b,p);
    //             if ( (!m_mapModified.is_free(b,p)) || (m_basisCheck[rowIdx]) )
    //             // if ( (m_basisCheck[rowIdx]) )
    //                 continue;
    //             colIdx = m_mapOriginal.index(b,p);
    //             m_matrix(rowIdx,colIdx) = 1;
    //             // gsInfo<<"Basis function "<<rowIdx<<"(patch: "<<p<<"; fun: "<<b<<") is "<< (m_basisCheck[rowIdx] ? "" : "not ")<<"processed\n";
    //         }

    //     // for(index_t k=0; k!=m_basisCheck.size(); k++)
    //     // Handle interior basis

    //     bool checkSides = std::all_of(m_sideCheck.begin(), m_sideCheck.end(), [](bool m_sideCheck) { return m_sideCheck; });
    //     GISMO_ASSERT(checkSides,"Not all sides are checked");
    //     bool checkVerts = std::all_of(m_vertCheck.begin(), m_vertCheck.end(), [](bool m_vertCheck) { return m_vertCheck; });
    //     GISMO_ASSERT(checkVerts,"Not all vertices are checked");

    //     // // normalize matrix
    //     // gsVector<T> sums(m_matrix.rows());
    //     // for (index_t k=0; k!=m_matrix.rows(); ++k)
    //     //     sums.at(k) = m_matrix.row(k).sum();

    //     // for (index_t k=0; k<m_matrix.outerSize(); ++k)
    //     //     for (typename gsSparseMatrix<T>::InnerIterator it(m_matrix,k); it; ++it)
    //     //         it.valueRef() = it.valueRef() / sums(it.row());
    // }


} // namespace gismo