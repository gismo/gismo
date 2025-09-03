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
struct InteriorSign { bool operator()(short_t s) const { return s>0; } };

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
            if ( c != node.nodeData().box.first [i] ) // avoid degenerate split
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
template<typename SignOp, short_t d, class T, class Z>
class gsTrimmedDomainIterator;

template<short_t d, class T, class Z=index_t>
class gsTrimmedDomain : public gsDomain<T>
{
public:
    typedef gsCutTreeData<d,Z> TData_t;
    typedef gsKdTree<d,Z,TData_t > Tree_t;

    typedef gsDomainIteratorWrapper<T> domainIter;
    typedef typename Tree_t::const_literator leafIterator;
    typedef typename Tree_t::point_t point_t;

    // template <class _T, short_t _d, class _Z>
    // friend class gsTrimmedDomainIterator;
    // template <class _T, short_t _d, class _Z>
    // friend class gsTrimmedDomainBoundaryIterator;

    template <typename SignOp, short_t _d, class _T, class _Z>
    friend class gsTrimmedDomainIterator;

public: // virtual interface

    virtual gsMatrix<T> boundingBox() const override
    {
        GISMO_NO_IMPLEMENTATION
    }

    virtual gsVector<short_t> sign(const gsMatrix<T> & u)
    {
        GISMO_NO_IMPLEMENTATION
    }

public: // non-virtual interface

    inline std::vector<bool> inDomain(const gsMatrix<T> & u)
    {
        std::vector<bool> res(u.cols());
        gsVector<T> val = this->sign(u);
        bool * r = res.data();
        for (T * a = val.data(); a != val.data()+val.size(); ++a)
            *(r++) = (*a<0 ? false : true);
        return res;
    }

    inline std::vector<bool> onBoundary(const gsMatrix<T> & u)
    {
        std::vector<bool> res(u.cols());
        gsVector<T> val = this->sign(u);
        bool * r = res.data();
        for (T * a = val.data(); a != val.data()+val.size(); ++a)
            *(r++) = (*a!=0 ? false : true);
        return res;
    }

    domainIter beginAll() const override
    {
        // GISMO_NO_IMPLEMENTATION
        return domainIter(new gsTrimmedDomainIterator<InteriorSign,d,T,Z>(m_tree,this->boundingBox(),m_upperIndex));
    }

    domainIter beginBdr(const boxSide bs) const override
    {
        // GISMO_NO_IMPLEMENTATION
       return domainIter(new gsTrimmedDomainIterator<BoundarySign,d,T,Z>(m_tree,this->boundingBox(),m_upperIndex));
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

        for ( index_t i = 0; i!=d; ++i )
            m_upperIndex[i] = tbasis.knots(i).numElements();


        // Make a box
        TData_t iData(point_t::Zero(), m_upperIndex);
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

        //----------------------------

        //First pass:
        // mark element position (-1,0,1) -- (inactive, cut-cell, interior)
        // create dofmapper (mark active basis functions)

        //Second pass:
        // assign active elements
        // subdivide into cut-cells upto cutlevel

    }

protected:
    Tree_t m_tree;
    point_t m_upperIndex;
};


/// Class representing an implicit trimmed domain
template<short_t d, class T, class Z=index_t>
class gsImplTrimmedDomain : public gsTrimmedDomain<d,T,Z>
{
private:
    typename gsFunction<T>::Ptr m_implFunction; // implicit function
    gsMatrix<T> m_boundingBox;

public:
    gsImplTrimmedDomain(const gsFunction<T> & fnc, const gsTensorBSplineBasis<d,T> & tbasis)
    :
    m_implFunction(memory::make_shared_not_owned(&fnc)),
    m_boundingBox(fnc.support())
    {
        this->init(tbasis, 3);
    }

    void setBoundingBox(const gsMatrix<T> & bb)
    {
        m_boundingBox = bb;
    }

    virtual gsMatrix<T> boundingBox() const override
    {
        return m_boundingBox;
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

template<typename SignOp, short_t d, class T, class Z>
class gsTrimmedDomainIterator : public gsDomainIterator<T>
{

public:

    typedef gsTrimmedDomain<d,T,Z> Domain_t;
    typedef typename Domain_t::TData_t TData_t;
    typedef typename Domain_t::Tree_t Tree_t;
    typedef typename Tree_t::const_literator leafIterator;
    typedef typename Tree_t::point_t point_t;


    typedef typename std::vector<T>::const_iterator  uiter;

    typedef typename gsDomainIterator<T>::uPtr domainIter;

    // elements:
    // bool active (true/false)
    // int position (-1,0,1) (out, cut, in)

public:

    gsTrimmedDomainIterator(const Tree_t & tree, const gsMatrix<T> & boundingBox, const point_t & upperIndex)
    :
    gsDomainIterator<T>(),
    m_tree(tree),
    m_boundingBox(boundingBox),
    m_upperIndex(upperIndex)
    {
        m_leaf = init(m_tree);
        if (!SignOp()(m_leaf.data().sign()))
            nextLeaf();
        updateLeaf();
    }

    gsTrimmedDomainIterator(const gsTrimmedDomainIterator & other) = default;
    domainIter clone() const override { return domainIter(new gsTrimmedDomainIterator(*this)); }

    leafIterator init(const Tree_t & tree)
    {
        // Initialize mesh data
        m_meshStart.resize(d);
        m_meshEnd  .resize(d);

        // Initialize cell data
        m_curElement.resize(d);

        // Allocate breaks
        m_breaks = std::vector<std::vector<T> >(d, std::vector<T>());

        return tree.beginLeafIterator();
    }

    // ---> Documentation in gsDomainIterator.h
    void next() override
    {
        gsDebug<<"Next\n";
        bool isGood(m_leaf.good());
        if (m_leaf.good())
        {
            if (SignOp()(m_leaf.data().sign()))
            {
                gsDebug<<"Leaf with corners "<<m_leaf.data().lowerCorner().transpose()
                      <<"~"<<m_leaf.data().upperCorner().transpose()
                      <<" and sign "<<m_leaf.data().sign()<<"\n";

                isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);
                if (!isGood) // went through all elements in m_leaf
                {
                    isGood = nextLeaf();
                    gsDebug<<"Going to next leaf\n";
                }
            }
            else
            {
                isGood = nextLeaf();
                gsDebug<<"Going to next leaf\n";
            }
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
        gsDebug<<"Next with increment "<<increment<<"\n";
        for (index_t i = 0; i < increment && isGood; ++i)
        {
            if (m_leaf.good())
            {
                if (SignOp()(m_leaf.data().sign()))
                {
                    gsDebug<<"Leaf with corners "<<m_leaf.data().lowerCorner().transpose()
                      <<"~"<<m_leaf.data().upperCorner().transpose()
                      <<" and sign "<<m_leaf.data().sign()<<"\n";
                    isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);
                    if (!isGood) // went through all elements in m_leaf
                        isGood = nextLeaf();
                }
                else
                    isGood = nextLeaf();
            }
        }
    }

    /// Resets the iterator so that it can be used for another
    /// iteration through all boundary elements.
    void reset() override
    {
        m_leaf = m_tree.beginLeafIterator();
        updateLeaf();
    }

    gsVector<T> lowerCorner() const override
    {
        gsVector<T> lower;
        lower.resize(d);
        for (short_t i = 0; i < d ; ++i)
            lower[i] = *m_curElement[i];
        return lower;
    }

    gsVector<T> upperCorner() const override
    {
        gsVector<T> upper;
        upper.resize(d);
        for (short_t i = 0; i < d ; ++i)
            upper[i]  = *(m_curElement[i]+1);
        return upper;
    }

    short_t sign() const
    {
        return m_leaf.data().sign();
    }

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

    /// Computes lower, upper and center point of the current element, maps the reference
    /// quadrature nodes and weights to the current element, and computes the
    /// active functions.
    void updateLeaf()
    {
        point_t lower = m_leaf.data().lowerCorner();
        point_t upper = m_leaf.data().upperCorner();

        // Update leaf box
        for (size_t dim = 0; dim < d; ++dim)
        {
            index_t start = lower(dim);
            index_t end  = upper(dim) ;

            m_breaks[dim].clear();
            for (index_t index = start; index <= end; ++index)
            {
                T coord = (T)index / (T)m_upperIndex[dim];
                m_breaks[dim].push_back( m_boundingBox(dim,0)+coord*(m_boundingBox(dim,1)-m_boundingBox(dim,0)) );
            }

            m_curElement(dim) =
            m_meshStart(dim)  = m_breaks[dim].begin();

            // for n breaks, we have n - 1 elements (spans)
            m_meshEnd(dim) =  m_breaks[dim].end() - 1;
        }
    }

private:

    const Tree_t & m_tree;
    const gsMatrix<T> m_boundingBox;
    const point_t m_upperIndex;

        // The current leaf node of the tree
    leafIterator m_leaf;

    // Coordinates of the grid cell boundaries
    // \todo remove this member
    std::vector< std::vector<T> > m_breaks;

    // Extent of the tensor grid
    gsVector<uiter, d> m_meshStart, m_meshEnd;

    // Current element as pointers to it's supporting mesh-lines
    gsVector<uiter, d> m_curElement;

};

template<class T>
class gsCutCellRule
{ }; // todo

}//namespace gismo

using namespace gismo;

int main(int argc, char* argv[])
{
    gsCmdLine cmd("Immersing the geometry and the basis.");
    // cmd.addPlainString("input", "G+Smo input basis file.", input);
    // cmd.addString("o", "output", "Name of the output file.", output);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // defines the background geometry
    gsTensorBSpline<2,real_t> bg = *gsNurbsCreator<>::BSplineSquare(5,0,0);
    bg.knots(0).transform(-1,1);
    bg.knots(1).transform(-1,1);
    // defines the inside or the outside of the (parametric) domain
    gsFunctionExpr<> impl_fun("1 - x^2 - y^2", 2);
    gsMatrix<> bb(2,2);
    bb << -1, 1,
          -1, 1;

    // uses the two objects above to get an immersed geometry
    //gsImmersedGeometry<real_t> igo(bg, inOut);

    // Background basis
    gsKnotVector<> kv(-1,1,9,2);
    gsTensorBSplineBasis<2,real_t> tbs(kv,kv);

    gsImplTrimmedDomain<2,real_t> tr_domain(impl_fun, tbs);
    tr_domain.setBoundingBox(bb);
    gsDebugVar(tr_domain.boundingBox());


    tr_domain.tree().printLeaves();


    // test iterator
    gsVector<> InteriorLabel(tr_domain.numElements());
    gsMatrix<> InteriorBoxes(2,tr_domain.numElements()*2);
    for (auto & elem : tr_domain.allElements())
    {
        gsInfo<<"Element with sign "<<static_cast<gsTrimmedDomainIterator<InteriorSign,2,real_t,index_t> &>(elem).sign()
              <<" and corners "<<elem.lowerCorner().transpose()<<" , "<<elem.upperCorner().transpose()<<"\n";
        InteriorLabel[elem.id()] = static_cast<gsTrimmedDomainIterator<InteriorSign,2,real_t,index_t> &>(elem).sign();
        InteriorBoxes.col(elem.id()*2)   = elem.lowerCorner();
        InteriorBoxes.col(elem.id()*2+1) = elem.upperCorner();
    }
    gsWriteParaview(InteriorBoxes,InteriorLabel,"interior_elements");

    gsVector<> BoundaryLabel(tr_domain.numElementsBdr());
    gsMatrix<> BoundaryBoxes(2,tr_domain.numElementsBdr()*2);
    for (auto & elem : gsDomain<real_t>::ElementRange( tr_domain.beginBdr(boundary::none), tr_domain.endBdr(boundary::none) ) )
    {
        gsInfo<<"Element with sign "<<static_cast<gsTrimmedDomainIterator<BoundarySign,2,real_t,index_t> &>(elem).sign()
              <<" and corners "<<elem.lowerCorner().transpose()<<" , "<<elem.upperCorner().transpose()<<"\n";

        BoundaryLabel[elem.id()] = static_cast<gsTrimmedDomainIterator<BoundarySign,2,real_t,index_t> &>(elem).sign();
        BoundaryBoxes.col(elem.id()*2)   = elem.lowerCorner();
        BoundaryBoxes.col(elem.id()*2+1) = elem.upperCorner();
    }
    gsWriteParaview(BoundaryBoxes,BoundaryLabel,"boundary_elements");

    gsWriteParaview(impl_fun, bb, "implicit_function");

    gsVector<unsigned,2> upp;
    upp.setConstant(kv.numElements());
    //gsTrimmedDomain<2,unsigned> tree;
    //tree.init(upp,10);
    //tree.size();
    //tree.construct(inOut);
    //tree.printLeaves();
    //tree.coordinates(nullptr);

    return EXIT_SUCCESS;
}



