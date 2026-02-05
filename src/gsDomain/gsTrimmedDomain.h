/** @file gsTrimmedDomain.h

    @brief TODO

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomain.h>
#include <gsDomain/gsCutTreeData.h>
#include <gsAssembler/gsQuadrature.h>

namespace gismo
{

// Cells with sign 0
struct BoundarySign { bool operator()(short_t s) const { return s==0; } };
// Cells with sign -1
struct InteriorSign { bool operator()(short_t s) const { return s< 0; } };
// Cells with sign  1
struct ExteriorSign { bool operator()(short_t s) const { return s> 0; } };
// Cells with signs 0 or -1
struct AllSign      { bool operator()(short_t s) const { return s<=0; } };
// Cells with any sign
struct AnySign      { bool operator()(short_t s) const { return true; } };

/**
 * @brief TODO
 *
 * \ingroup Domain
 */

template<short_t d, class T, class Z=size_t>
class gsTrimmedDomain : public gsDomain<T>
{
public:
    typedef gsCutTreeData<d,Z> TData_t;
    typedef gsKdTree<d,Z,TData_t > Tree_t;
    typedef typename gsDomain<T>::Ptr BgDomain_t;
    typedef gsDomainIteratorWrapper<T> domainIter;
    typedef typename Tree_t::const_literator leafIterator;
    typedef typename Tree_t::point_t point_t;

    // template <class _T, short_t _d, class _Z>
    // friend class gsTrimmedDomainIterator;
    // template <class _T, short_t _d, class _Z>
    // friend class gsTrimmedDomainBoundaryIterator;

    template <typename SignOp, short_t _d, class _T, class _Z>
    friend class gsTrimmedDomainIterator;

protected:
    Tree_t m_tree; ///< The tree partition characterizing in/out/cut elements
    BgDomain_t m_bgDomain; ///< The background domain // NOTE: this should be a gsMatrix<

public: // virtual interface

    virtual gsMatrix<T> boundingBox() const override
    {
        // note: a better bounding box might be possible
        // return m_bgDomain->boundingBox();
        GISMO_NO_IMPLEMENTATION
    }

    /// Computes the sign at given points \a u in the domain
    virtual gsVector<short_t> sign(const gsMatrix<T> & u)
    {
        GISMO_NO_IMPLEMENTATION
    }

    ///
    // virtual gsVector<short_t> sign( NODE )


public: // non-virtual interface

    const gsDomain<T> & backgroundDomain() const { return *m_bgDomain; }

    inline std::vector<bool> inDomain(const gsMatrix<T> & u)
    {
        std::vector<bool> res(u.cols());
        gsVector<short_t> val = this->sign(u);
        bool * r = res.data();
        for (short_t * a = val.data(); a != val.data()+val.size(); ++a)
            *(r++) = (*a<0 ? false : true);
        return res;
    }

    inline std::vector<bool> onBoundary(const gsMatrix<T> & u)
    {
        std::vector<bool> res(u.cols());
        gsVector<short_t> val = this->sign(u);
        bool * r = res.data();
        for (short_t * a = val.data(); a != val.data()+val.size(); ++a)
            *(r++) = (*a!=0 ? false : true);
        return res;
    }

    domainIter beginAll() const override
    {
        return begin<AllSign>();
    }

    domainIter beginBdr(const boxSide bs) const override
    {
        return begin<BoundarySign>();
    }

    domainIter beginActive() const
    {
        return beginAll();
    }

    domainIter beginExterior() const
    {
        return begin<ExteriorSign>();
    }

    domainIter beginInterior() const
    {
        return begin<InteriorSign>();
    }

    domainIter beginAny() const
    {
        return begin<AnySign>();
    }

    template<typename SignOp>
    domainIter begin() const
    {
        return domainIter(new gsTrimmedDomainIterator<SignOp,d,T,Z>(*this));
    }

    template<typename SignOp>
    domainIter end() const
    {
        return domainIter(new gsDomainIteratorEnd<T>(this->template numElements<SignOp>()));
    }

    template<typename SignOp>
    size_t numElements() const
    {
        leafIterator it = m_tree.beginLeafIterator(); //generic tree leaf iterator
        size_t nel(0);
        point_t lc, uc;
        while (it.good())
        {
            if (SignOp()(it.data().sign()))
            {
                uc = it.data().upperCorner();
                lc = it.data().lowerCorner();
                nel += (uc - lc).prod();
            }
            it.next();
        }
        return nel;
    }

    size_t numElements() const override
    {
        return numElements<AllSign>();
    }

    size_t numElementsBdr(boxSide const & s = boundary::none) const override
    {
        return numElements<BoundarySign>();
    }

    short_t dim() const override { return d; }

    const Tree_t & tree() const { return m_tree; }

protected:



    /*
        // bool SignOp.operator()(gsCutTreeData<XXX> & data) // sets data.sign(), returns false if splitting needed

        struct SignOp // ASSUMES LEAF
        {
            SignOp( const gsTensorBSplineBasis<d,T> & tbasis, index_t samples )
            : tbasis(tbasis)
            , QuRule( gsVector<index_t>::Constant(d,samples) )
            {}

            bool operator()(gsCutTreeData<XXX> & data)
            {
                auto & k1 = curNode->nodeData().lowerCorner();
                auto & k2 = curNode->nodeData().upperCorner();
                for(short_t j = 0; j < d;j++)
                {
                    const gsKnotVector<T> & kv = tbasis.knots(j);
                    u1[j] = kv.uValue(k1[j]);
                    u2[j] = kv.uValue(k2[j]);
                }

                QuRule.mapTo(u1, u2, pts, wts);
                sgn = this->sign(pts);

                if ( (sgn.array() == 1).all() ) // domain active ?
                {
                    data.sign() = 1;
                    return true;
                }
                if ( (sgn.array() == -1).all() ) // domain inactive ?
                {
                    data.sign() = - 1;
                    return true;
                }
                // There are both positive and negative signs
                if ( 1 == (k2-k1).prod() )
                {
                    data.sign() = 0;
                    return true;
                }
                else
                {
                    return false; // splitting needed
                }
            }
        };

        // OTHER SIGN OPS
        // * With sub-grid iterator (now gsGridIterator) over elements in box, using:
        //   TODO: For all elements inside [k1,k2].....
                    // Could be replaced by a gsSubTensorDomainIterator
                    gsGridIterator<Z,CUBE,d> git(k1, k2);
                    for( ; git; ++git)
                    {
                        // get element corners in param. space
                        for(short_t j = 0; j < d;j++)
                        {
                            const gsKnotVector<T> & kv = tbasis.knots(j);
                            u1[j] = kv.uValue( git.current()[j]     );
                            u2[j] = kv.uValue( git.current()[j] + 1 );
                        }

                        QuRule.mapTo(u1, u2, pts, wts);
                        sgn = this->sign(pts);
                        // process sgn ...
                    }
            // One not depending on tensor product basis (for THB)


        template<typename SignOp, short_t d, class T, class Z>
        void init(const gsTensorBSplineBasis<d,T> & tbasis, index_t samples = 3)
        {

            SignOp signOp;
            Tree_t * curNode;
            while ( ! stack.empty() )
            {
                curNode = stack.back(); //top();
                stack.pop_back();       //pop();

                if ( curNode->isLeaf() ) // reached a leaf
                {
                    if (!signOp( curNode->nodeData() ))
                    {
                        // splitting needed
                        curNode->anyMidSplit(1);
                        stack.push_back(curNode->left );
                        stack.push_back(curNode->right);
                    }
                }
                else // roll down the tree
                {
                    stack.push_back(curNode->left );
                    stack.push_back(curNode->right);
                }
            }
        }
    */

    void init(T maxElementSize = T(1), T minElementSize = T(0.1), index_t samples = 10)
    {
        /*
            Generates a binary partitioning tree of the parametric domain bounded
            by a boundingBox. Each leaf node is marked as inside, outside or
            cut, based on the implicit function sign evaluated at quadrature
            points with \a samples samples within the node's bounding box.
            The maximimum and minimum element sizes can be controlled via
            \a maxElementSize and \a minElementSize respectively. Finally,
            all elements will have a size inside [minElementSize, maxElementSize].

            The algorithm will split partitions until one of the following
            conditions is met:
            - The partition is marked as inside or outside and the partition size is smaller than \a maxElementSize.
            - The partition is marked as cut and the partition size is larger than \a minElementSize and smaller than \a maxElementSize.
        */

    }

    void init(const gsHTensorBasis<d,T> & tbasis, index_t samples = 3)
    {
        /* 
            Similar to init with gsTensorBSplineBasis, but for hierarchical
            tensor b-spline bases.

            Starts from the tree defined by the hierarchical basis, and marks
            each leaf node as inside, outside or cut, based on the implicit
            function sign evaluated at quadrature points with \a samples samples
            within the node's bounding box.

            Like for the gsTensorBSplineBasis version, the algorithm will stop
            when a node has 1 element ((k2-k1).prod() == 1), or when the node is
            marked as inside or outside.
        */
    }


    /** 
     *  Initializes `m_tree` and `m_bgDomain` based on the provided tensor B-spline basis \a tbasis.
     * 
    */
    void init(const gsTensorBSplineBasis<d,T> & tbasis, index_t samples = 3)
    {
        gsVector<index_t>  numNodes(d);
        numNodes.setConstant(samples);
        gsLobattoRule<real_t> QuRule(numNodes);

        point_t upperIndex;
        for ( index_t i = 0; i!=d; ++i )
            upperIndex[i] = tbasis.knots(i).numElements();


        // Make a box
        TData_t iData(point_t::Zero(), upperIndex);
        m_tree = Tree_t(iData);

        // Initialize stack
        std::vector<Tree_t*> stack;
        stack.reserve(  std::pow(4,d) );
        stack.push_back(&m_tree);

        gsVector<T,d> u1, u2;
        gsVector<T> wts;
        gsMatrix<T> pts;
        gsVector<short_t> sgn;

        Tree_t * curNode;
        while ( ! stack.empty() )
        {
            curNode = stack.back(); //top();
            stack.pop_back();       //pop();

            if ( curNode->isLeaf() ) // reached a leaf
            {
                auto & k1 = curNode->nodeData().lowerCorner();
                auto & k2 = curNode->nodeData().upperCorner();
                for(short_t j = 0; j < d;j++)
                {
                    const gsKnotVector<T> & kv = tbasis.knots(j);
                    u1[j] = kv.uValue(k1[j]);
                    u2[j] = kv.uValue(k2[j]);
                }
                //question: treat Trivial boxes ?

                // TODO: For all elements inside [k1,k2].....
                /*
                    // Could be replaced by a gsSubTensorDomainIterator
                    gsGridIterator<Z,CUBE,d> git(k1, k2);
                    for( ; git; ++git)
                    {
                        // get element corners in param. space
                        for(short_t j = 0; j < d;j++)
                        {
                            const gsKnotVector<T> & kv = tbasis.knots(j);
                            u1[j] = kv.uValue( git.current()[j]     );
                            u2[j] = kv.uValue( git.current()[j] + 1 );
                        }

                        QuRule.mapTo(u1, u2, pts, wts);
                        sgn = this->sign(pts);
                        // process sgn ...
                    }
                */
                QuRule.mapTo(u1, u2, pts, wts);
                sgn = this->sign(pts);

                if ( (sgn.array() == 1).all() ) // domain active ?
                {
                    curNode->nodeData().sign() = 1;
                    continue;
                }

                if ( (sgn.array() == -1).all() ) // domain inactive ?
                {
                    curNode->nodeData().sign() = - 1;
                    continue;
                }

                // There are both positive and negative signs
                if ( 1 == (k2-k1).prod() )
                {
                    curNode->nodeData().sign() = 0;
                }
                else //more than one element
                {
                    curNode->anyMidSplit(1);
                    stack.push_back(curNode->left );
                    stack.push_back(curNode->right);
                }
            }
            else // roll down the tree
            {
                stack.push_back(curNode->left );
                stack.push_back(curNode->right);
            }
        }

        m_bgDomain = tbasis.domain();

        //----------------------------

        //First pass:
        // mark element position (-1,0,1) -- (inactive, cut-cell, interior)
        // create dofmapper (mark active basis functions)

        //Second pass:
        // assign active elements
        // subdivide into cut-cells upto cutlevel

    }
};

} // namespace gismo
