/** @file immersed_example.cpp

    @brief Tutorial on gsBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
*/

#include <iostream>
#include <gismo.h>

#include <gsDomain/gsKdTree.h> // A binary partition tree for representing hierarchical domains.

// why needed ?
//#include <gsHSplines/gsHDomain.hpp>

namespace gismo
{

struct BoundarySign { bool operator()(short_t s) const { return s==0; } };
struct InteriorSign { bool operator()(short_t s) const { return s>=0; } };

template<short_t d, class Z>
class gsCutTreeData
{
    typedef gsCutTreeData Self_t;
    typedef gsKdTree<d,Z,Self_t> Tree_t;
public:
    typedef          gsAABB<d, Z> kdBox;
    typedef typename kdBox::point point;
private:

    //int m_level;
    kdBox box;
    short_t m_sign;

public:
    gsCutTreeData(point const & k1, point const & k2 /*, int lvl = 0*/)
    : /*m_level(lvl),*/ box(k1,k2)
    { }

    explicit gsCutTreeData(point const & k2)
    : /*m_level(0),*/ box(k2)
    { }

    friend std::ostream & operator<<(std::ostream & os, const Self_t & td)
    {
        os << "Box("/*<<td.m_level<<*/"): ["<< td.box.first.transpose()
           <<"]~["<<td.box.second.transpose() <<"]\n";
        return os;
    }

    short_t sign() const { return m_sign; }

    short_t & sign() { return m_sign; }

    bool check() const
    { return (box.second.array() >= box.first.array()).any(); }

    static void split(const Tree_t & node)
    {
        node.left->nodeData().upperCorner()[node.axis] =
            node.right->nodeData().lowerCorner()[node.axis] = node.pos;
    }

    static void nextMidSplit(Tree_t &  node)
    {
        node.axis = ( node.parent == 0 ? 0 : (node.parent->axis+1)%d );
        node.pos  = node.box->first[node.axis] +
            (node.box->second[node.axis] - node.box->first[node.axis])/2 ;
    }

    static void anyMidSplit(unsigned h, Tree_t &  node,
                            int & doSplit)
    {
        const unsigned mask = ~(h - 1);
        for ( unsigned i = 0; i < d; ++i )
        {
            const unsigned c =
                (node.nodeData().box.first [i] + (node.nodeData().box.second[i] - node.nodeData().box.first[i])/2) & mask ;
            if ( c != (unsigned)node.nodeData().box.first [i] ) // avoid degenerate split
            {
                node.axis = i;
                node.pos  = c;
                doSplit   = 1;
                return;
            }
        }
        doSplit = 0;
    }

    // Merges the data of right child into left child of the node
    //template<short_t d, class Z = index_t>
    static void mergeToLeft(const gsKdTree<d,Z,Self_t> & node)
    {
        node.left ->nodeData().upperCorner()[node.axis] =
            node.right->nodeData().upperCorner()[node.axis];
    }

    void multiplyByTwo()
    {
        box.first .array() *= 2;
        box.second.array() *= 2;
    }

    void divideByTwo()
    {
        box.first .array() /= 2;
        box.second.array() /= 2;
    }

    /*
    int level() const {return m_level;}
    int & level() {return m_level;}
    */

    const point & lowerCorner() const { return box.first; }
          point & lowerCorner()       { return box.first; }

    const point & upperCorner() const { return box.second; }
          point & upperCorner()       { return box.second; }

    const kdBox & aabb() const { return box; }
    kdBox & aabb() { return box; }

};

// Forward declaration
template<typename SignOp, short_t d, class T, class Z=size_t>
class gsTrimmedDomainIterator;

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
    BgDomain_t m_bgDomain; ///< The background domain

public: // virtual interface

    virtual gsMatrix<T> boundingBox() const override
    {
        // note: a better bounding box might be possible
        return m_bgDomain->boundingBox();
    }

    /// Computes the sign at given points \a u in the domain
    virtual gsVector<short_t> sign(const gsMatrix<T> & u)
    {
        GISMO_NO_IMPLEMENTATION
    }

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
        return domainIter(new gsTrimmedDomainIterator<InteriorSign,d,T,Z>(*this));
    }

    domainIter beginBdr(const boxSide bs) const override
    {
        return domainIter(new gsTrimmedDomainIterator<BoundarySign,d,T,Z>(*this));
    }

    size_t numElements() const override
    {
        leafIterator it = m_tree.beginLeafIterator(); //generic tree leaf iterator
        size_t nel(0);
        point_t lc, uc;
        while (it.good())
        {
            if (InteriorSign()(it.data().sign()))
            {
                uc = it.data().upperCorner();
                lc = it.data().lowerCorner();
                nel += (uc - lc).prod();
            }
            it.next();
        }
        return nel;
    }

    size_t numElementsBdr(boxSide const & s = boundary::none) const override
    {
        leafIterator it = m_tree.beginLeafIterator(); //generic tree leaf iterator
        size_t nel(0);
        point_t lc, uc;
        while (it.good())
        {
            if (BoundarySign()(it.data().sign()))
            {
                uc = it.data().upperCorner();
                lc = it.data().lowerCorner();
                nel += (uc - lc).prod();
            }
            it.next();
        }
        return nel;
    }

    short_t dim() const override { return d; }

    const Tree_t & tree() const { return m_tree; }

protected:

    void init(T tol, index_t samples = 10)
    {

    }

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

template<typename SignOp, short_t d, class T, class Z /*= size_t*/>
class gsTrimmedDomainIterator : public gsDomainIterator<T>
{
public:

    typedef gsTrimmedDomain<d,T,Z> Domain_t;
    typedef typename Domain_t::TData_t TData_t;
    typedef typename Domain_t::Tree_t Tree_t;
    typedef typename Tree_t::const_literator leafIterator;
    typedef typename Tree_t::point_t point_t;

    typedef gsGridIterator<Z,CUBE,d> gridIterator;

    typedef typename std::vector<T>::const_iterator  uiter;

    typedef typename gsDomainIterator<T>::uPtr domainIter;
    typedef gsDomainIteratorWrapper<T> domainIterWrapper;

    // elements:
    // bool active (true/false)
    // int position (-1,0,1) (out, cut, in)

private:
    // The trimmed domain being iterated on
    const Domain_t & m_tdomain;

    // The current leaf node of the tree
    leafIterator m_leaf;

    // The current grid point of the leaf
    gridIterator m_gridIter;

    /// The current element
    domainIterWrapper m_current;

public:

    explicit gsTrimmedDomainIterator(const Domain_t & _tdomain)
    :
    gsDomainIterator<T>(),
    m_tdomain(_tdomain), m_current(give(_tdomain.backgroundDomain().beginAll()))
    {
        reset();
    }

    gsTrimmedDomainIterator(const gsTrimmedDomainIterator & other) = default;
    domainIter clone() const override { return domainIter(new gsTrimmedDomainIterator(*this)); }

    // ---> Documentation in gsDomainIterator.h
    void next() override
    {
        ++m_gridIter;
        if (m_gridIter)
            m_current.skipTo( *m_gridIter );//copy gsVector<D> -->gsVector<>
        else  // went through all elements in m_leaf
        {
            nextLeaf();
        }

    }

    // ---> Documentation in gsDomainIterator.h
    void next(index_t increment) override
    {
        // increment through elements
        //todo: better implementation
        // compute the number of elements between curElement and meshEnd
        // use m_leaf.numElements() to skip leaves
        // arrive at the element or end
        bool isGood(m_leaf.good());
        //gsDebug<<"Next with increment "<<increment<<"\n";
        for (index_t i = 0; i < increment && isGood; ++i)
        {
            if (SignOp()(m_leaf.data().sign()))
            {
                // gsDebug<<"Leaf with corners "<<m_leaf.data().lowerCorner().transpose()
                //   <<"~"<<m_leaf.data().upperCorner().transpose()
                //   <<" and sign "<<m_leaf.data().sign()<<"\n";
                ++m_gridIter;
                if (!m_gridIter) // went through all elements in m_leaf
                    isGood = nextLeaf();
            }
            else
                isGood = nextLeaf();
        }

        if (m_gridIter)
            m_current.skipTo( *m_gridIter );//copy gsVector<D> -->gsVector<>
    }

    /// Resets the iterator so that it can be used for another
    /// iteration through all boundary elements.
    void reset() override
    {
        m_leaf = m_tdomain.tree().beginLeafIterator();
        if (m_leaf.good() && !SignOp()(m_leaf.data().sign()))
            nextLeaf();
        updateLeaf();
    }

    gsVector<T> lowerCorner() const override { return m_current.lowerCorner(); }

    gsVector<T> upperCorner() const override { return m_current.upperCorner(); }

    short_t sign() const { return m_leaf.data().sign(); }

private:

    /// returns true if there is a another leaf with a boundary element
    bool nextLeaf()
    {
        bool isGood = m_leaf.next();
        while (isGood)
        {
            if (SignOp()(m_leaf.data().sign()))
            {
                updateLeaf();
                return true;
            }
            isGood = m_leaf.next();
        }
        return false;
    }

    void updateLeaf()
    {
        if(m_leaf.good())
        {
            m_gridIter.reset( m_leaf.data().lowerCorner(), m_leaf.data().upperCorner() );
            m_current.skipTo( *m_gridIter );//copy gsVector<D> -->gsVector<>
        }
    }
};

/// Class representing an implicit trimmed domain
template<short_t d, class T>
class gsImplTrimmedDomain : public gsTrimmedDomain<d,T>
{
    typedef gsTrimmedDomain<d,T> Base;
private:
    typename gsFunction<T>::Ptr m_implFunction; // implicit function

public:
    gsImplTrimmedDomain(const gsFunction<T> & fnc,
                        const gsTensorBSplineBasis<d,T> & tbasis) :
    m_implFunction(memory::make_shared_not_owned(&fnc))
    {
        this->init(tbasis, 3);
    }

    inline gsVector<short_t> sign(const gsMatrix<T> & u)
    {
        gsVector<T> val = m_implFunction->eval(u).row(0);
        return gsVector<short_t>(val.array().sign().template cast<short_t>());
    }

};

template<short_t d, class T>
class gsImmersedGeometry// : public gsFunction<T>
{
private:
    typename gsFunction<T>::Ptr m_bgGeo; // background geometry
    gsImplTrimmedDomain<d,  T> m_trDomain;  // trimmed domain

public:
    gsImmersedGeometry(const gsFunction<T> & bgGeo,
                       const gsImplTrimmedDomain<d, T> & trDomain) :
    m_bgGeo(bgGeo.clone()), m_trDomain(trDomain)
    { }

public://function interface?

public:
    const gsFunction<T> & background() const { return *m_bgGeo; }

    const gsImplTrimmedDomain<d, T> & domain() const { return m_trDomain; }

    /// Return a triangulation of the boundaty of the immersed geometry
    memory::unique_ptr<gsMesh<T> > toBoundaryMesh(int npoints = 50) const
    {
        gsVector<index_t,d>  numNodes; numNodes.setConstant(3);
        gsLobattoRule<real_t> qurule(numNodes);

        auto supp = m_bgGeo.support();
        // gsKdTree
    }

    /// Return a triangulation of the immersed geometry
    memory::unique_ptr<gsMesh<T> > toVolumeMesh(int npoints = 50) const
    {

    }

};

template<class T>
class gsCutCellRule : public gsQuadRule<T>
{
public:

    typedef memory::unique_ptr<gsCutCellRule> uPtr;

    gsCutCellRule(const gsQuadRule<T> & rule, const gsFunction<T> & levelSet)
    :
    m_rule(rule), m_levelSet(levelSet)
    {

    }

    static uPtr make(const gsQuadRule<T> & rule, const gsFunction<T> & levelSet)
    {
        return uPtr(new gsCutCellRule(rule, levelSet));
    }

    inline void mapTo( const gsVector<T>& lower, const gsVector<T>& upper,
                             gsMatrix<T>& nodes,       gsVector<T>& weights ) const override
    {
        gsMatrix<T> vals;
        m_rule.mapTo(lower, upper, nodes, weights);
        m_levelSet.eval_into(nodes,vals);
        for (index_t i = 0; i < vals.cols(); ++i)
            if (vals(0,i) < 0)
            {
                weights[i] = 0;
            }
    }

protected:
    const gsQuadRule<T> & m_rule;
    const gsFunction<T> & m_levelSet;
};

}//namespace gismo

using namespace gismo;


void test_2D(int k)
{
    // defines the background geometry
    gsTensorBSpline<2,real_t> bg = *gsNurbsCreator<>::BSplineSquare(5,0,0);
    bg.knots(0).transform(-1,1);
    bg.knots(1).transform(-1,1);
    // defines the inside or the outside of the (parametric) domain

    // These ones have some problems
    // gsFunctionExpr<> impl_fun("1", 2);
    // gsFunctionExpr<> impl_fun("-1", 2);
    // gsFunctionExpr<> impl_fun("0", 2);

    gsFunctionExpr<> impl_fun("1 - (x-1)^2 - (y-1)^2", 2);
    //gsFunctionExpr<> impl_fun("1 - x^2 - y^2", 2);
    //gsFunctionExpr<> impl_fun("x-y", 2);
    //gsFunctionExpr<> impl_fun("x*y", 2);
    gsMatrix<> bb(2,2);
    bb << -1, 1,
          -1, 1;

    // uses the two objects above to get an immersed geometry
    //gsImmersedGeometry<real_t> igo(bg, inOut);

    // Background basis
    gsKnotVector<> kv(-1,1,k,2);
    gsTensorBSplineBasis<2,real_t> tbs(kv,kv);

    gsImplTrimmedDomain<2,real_t> tr_domain(impl_fun, tbs);
    gsDebugVar(tr_domain.boundingBox());

    tr_domain.tree().printLeaves();

    // test iterator
    gsVector<> InteriorLabel(tr_domain.numElements());
    gsMatrix<> InteriorBoxes(2,tr_domain.numElements()*2);
    for (auto & elem : tr_domain.allElements())
    {
        gsInfo<<"Element with sign "<<static_cast<gsTrimmedDomainIterator<InteriorSign,2,real_t> &>(elem).sign()
              <<" and corners "<<elem.lowerCorner().transpose()<<" , "<<elem.upperCorner().transpose()<<"\n";
        InteriorLabel[elem.id()] = static_cast<gsTrimmedDomainIterator<InteriorSign,2,real_t> &>(elem).sign();
        InteriorBoxes.col(elem.id()*2)   = elem.lowerCorner();
        InteriorBoxes.col(elem.id()*2+1) = elem.upperCorner();
    }
    gsWriteParaview(InteriorBoxes,InteriorLabel,"interior_elements2d");

    gsVector<> BoundaryLabel(tr_domain.numElementsBdr());
    gsMatrix<> BoundaryBoxes(2,tr_domain.numElementsBdr()*2);
    for (auto & elem : gsDomain<real_t>::ElementRange( tr_domain.beginBdr(boundary::none), tr_domain.endBdr(boundary::none) ) )
    {
        gsInfo<<"Element with sign "<<static_cast<gsTrimmedDomainIterator<BoundarySign,2,real_t> &>(elem).sign()
              <<" and corners "<<elem.lowerCorner().transpose()<<" , "<<elem.upperCorner().transpose()<<"\n";

        BoundaryLabel[elem.id()] = static_cast<gsTrimmedDomainIterator<BoundarySign,2,real_t> &>(elem).sign();
        BoundaryBoxes.col(elem.id()*2)   = elem.lowerCorner();
        BoundaryBoxes.col(elem.id()*2+1) = elem.upperCorner();
    }
    gsWriteParaview(BoundaryBoxes,BoundaryLabel,"boundary_elements2d");

    gsWriteParaview(impl_fun, bb, "implicit_function2d");

    gsMatrix<> pts;
    gsVector<> wts;
    gsGaussRule<real_t> rule(gsVector<index_t,2>::Constant(5));
    gsCutCellRule<real_t> ccrule(rule, impl_fun);
    real_t area(0.0);
    for (auto & elem : tr_domain.allElements())
    {
        ccrule.mapTo(elem.lowerCorner(), elem.upperCorner(), pts, wts);
        area += wts.sum();
    }
    gsInfo<<"Area = "<<area<<"\n";
}


void test_3D(int k)
{
    // defines the background geometry
    gsTensorBSpline<3,real_t> bg = *gsNurbsCreator<>::BSplineCube(5,0,0);
    bg.knots(0).transform(-1,1);
    bg.knots(1).transform(-1,1);
    bg.knots(2).transform(-1,1);
    // defines the inside or the outside of the (parametric) domain
    gsFunctionExpr<> impl_fun("1 - x^2 - y^2 - z^2", 3);
    gsMatrix<> bb(3,2);
    bb << -1, 1,
          -1, 1,
          -1, 1;

    // uses the two objects above to get an immersed geometry
    //gsImmersedGeometry<real_t> igo(bg, inOut);

    // Background basis
    gsKnotVector<> kv(-1,1,k,2);
    gsTensorBSplineBasis<3,real_t> tbs(kv,kv,kv);

    gsImplTrimmedDomain<3,real_t> tr_domain(impl_fun, tbs);

    //tr_domain.tree().printLeaves();

    // test iterator
    gsVector<> InteriorLabel(tr_domain.numElements());
    gsMatrix<> InteriorBoxes(3,tr_domain.numElements()*2);
    for (auto & elem : tr_domain.allElements())
    {
        gsInfo<<"Element with sign "<<static_cast<gsTrimmedDomainIterator<InteriorSign,3,real_t> &>(elem).sign()
              <<" and corners "<<elem.lowerCorner().transpose()<<" , "<<elem.upperCorner().transpose()<<"\n";
        InteriorLabel[elem.id()] = static_cast<gsTrimmedDomainIterator<InteriorSign,3,real_t> &>(elem).sign();
        InteriorBoxes.col(elem.id()*2)   = elem.lowerCorner();
        InteriorBoxes.col(elem.id()*2+1) = elem.upperCorner();
    }
    gsWriteParaview(InteriorBoxes,InteriorLabel,"interior_elements3d");

    gsVector<> BoundaryLabel(tr_domain.numElementsBdr());
    gsMatrix<> BoundaryBoxes(3,tr_domain.numElementsBdr()*2);
    for (auto & elem : gsDomain<real_t>::ElementRange( tr_domain.beginBdr(boundary::none), tr_domain.endBdr(boundary::none) ) )
    {
        gsInfo<<"Element with sign "<<static_cast<gsTrimmedDomainIterator<BoundarySign,3,real_t> &>(elem).sign()
              <<" and corners "<<elem.lowerCorner().transpose()<<" , "<<elem.upperCorner().transpose()<<"\n";

        //BoundaryLabel[elem.id()] = static_cast<gsTrimmedDomainIterator<BoundarySign,3,real_t,index_t> &>(elem).sign();
        BoundaryLabel[elem.id()] = elem.id();
        BoundaryBoxes.col(elem.id()*2)   = elem.lowerCorner();
        BoundaryBoxes.col(elem.id()*2+1) = elem.upperCorner();
    }
    gsWriteParaview(BoundaryBoxes,BoundaryLabel,"boundary_elements3d");

    gsWriteParaview(impl_fun, bb, "implicit_function3d");

    gsMatrix<> pts;
    gsVector<> wts;
    gsGaussRule<real_t> rule(gsVector<index_t,3>::Constant(5));
    gsCutCellRule<real_t> ccrule(rule, impl_fun);
    real_t area(0.0);
    for (auto & elem : tr_domain.allElements())
    {
        ccrule.mapTo(elem.lowerCorner(), elem.upperCorner(), pts, wts);
        gsDebugVar(wts.transpose());
        area += wts.sum();
    }
    gsInfo<<"Area = "<<area<<"\n";
}

int main(int argc, char* argv[])
{
    index_t numKnots  = 9;
    gsCmdLine cmd("Immersing the geometry and the basis.");
    cmd.addInt( "k", "knots", "Number of knots per direction",  numKnots );
    // cmd.addPlainString("input", "G+Smo input basis file.", input);
    // cmd.addString("o", "output", "Name of the output file.", output);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    test_2D(numKnots);
    test_3D(numKnots);

    return EXIT_SUCCESS;
}
