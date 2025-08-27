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
                (node.box.first [i] + (node.box.second[i] - node.box.first[i])/2) & mask ;
            if ( c != node.box.first [i] ) // avoid degenerate split
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

template<short_t d, class T>
class gsTrimmedDomain : public gsDomain<T>
{
public:
    typedef gsCutTreeData<d,index_t> TData_t;
    typedef gsKdTree<d,index_t,TData_t > Tree_t;

    typedef gsDomainIteratorWrapper<T> domainIter;
    typedef typename Tree_t::const_literator leafIterator;
    typedef typename Tree_t::point point_t;

    template <class _T, short_t _d, class _Z>
    friend class gsTrimmedDomainIterator;
    template <class _T, short_t _d, class _Z>
    friend class gsTrimmedDomainBoundaryIterator;

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
        GISMO_NO_IMPLEMENTATION
        //return domainIter(new gsTrimmedDomainIterator<T,d,Z>(m_tree));
    }

    domainIter beginBdr(const boxSide bs) const override
    {
        GISMO_NO_IMPLEMENTATION
       //return domainIter(new gsTrimmedDomainBoundaryIterator<T,d,Z>(m_tree,bs));
    }

    size_t numElements() const override
    {
        leafIterator it = m_tree.beginLeafIterator(); //generic tree leaf iterator
        size_t nel(0);
        point_t lc, uc;
        while (it.good())
        {
            uc = it.data().upperCorner();
            lc = it.data().lowerCorner();
            nel += (uc - lc).prod();
            it.next();
        }
        return nel;
    }

    size_t numElementsBdr(boxSide const & s = boundary::none) const override
    {
        return 0;
    }

    short_t dim() const override { return d; }

    const Tree_t & tree() const { return m_tree; }

private:

    void init(T tol, index_t samples = 10)
    {

    }

    void init(gsTensorBasis<d,T> & tbasis, index_t samples = 3)
    {
        gsVector<index_t>  numNodes(d);
        numNodes.setConstant(samples);
        gsLobattoRule<real_t> QuRule(numNodes);

        for ( index_t i = 0; i!=d; ++i )
            m_upperIndex[i] = tbasis.knots(i).numElements();


        // Make a box
        TData_t iData(point_t::Zeros(), m_upperIndex);
        m_tree = Tree_t(iData);
        
        // Initialize stack
        std::vector<TData_t*> stack;
        stack.reserve(  std::pow(4,d) );
        stack.push_back(&m_tree);

        gsVector<T> wts, u1, u2;
        gsMatrix<T> pts;
        gsVector<short_t> sgn;

        node_t * curNode;
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
                QuRule->mapTo(u1, u2, pts, wts);
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
                    node * newLeaf = curNode->anyMidSplit(1);
                    stack.push_back(newLeaf);
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
    const Tree_t m_tree;
    point_t m_upperIndex;
};


/// Class representing an implicit trimmed domain
template<short_t d, class T>
class gsImplTrimmedDomain : gsTrimmedDomain<T>
{
private:
    typename gsFunction<T>::Ptr m_implFunction; // implicit function
    
public:
    gsImplTrimmedDomain(const gsFunction<T> & fnc, gsTensorBasis<d,T> & tbasis)
    {
        this->init(tbasis, 3);
    }

    virtual gsMatrix<T> boundingBox() const override
    {
        return fnc.support();
    }
    
    inline gsVector<int> sign(const gsMatrix<T> & u)
    {
        gsMatrix<T> val = m_implFunction->eval(u);
        return val.array().sign();
    }

};

template<class T>
class gsImmersedGeometry// : public gsFunction<T>
{
private:
    typename gsFunction<T>::Ptr m_bgGeo; // background geometry
    gsImplTrimmedDomain<T> m_trDomain;  // trimmed domain

public:
    gsImmersedGeometry(const gsFunction<T> & bgGeo,
                       const gsImplTrimmedDomain<T> & trDomain) :
    m_bgGeo(bgGeo.clone()), m_trDomain(trDomain)
    { }

public://function interface?

public:
    const gsFunction<T> & background() const { return *m_bgGeo; }

    const gsImplTrimmedDomain<T> & domain() const { return m_trDomain; }

    /// Return a triangulation of the boundaty of the immersed geometry
    memory::unique_ptr<gsMesh<T> > toBoundaryMesh(int npoints = 50) const
    {
        gsVector<index_t>  numNodes(2); numNodes.setConstant(3);
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
class gsTrimmedDomainIterator
{
    // elements:
    // bool active (true/false)
    // int position (-1,0,1) (out, cut, in)
private:
    gsDofMapper m_mapper;

    // gsKdTree
    
public:
    gsTrimmedDomainIterator(const gsImplTrimmedDomain<T> & trdom,
                            const gsBasis<T> & basis, index_t cutlevel = 0)
    {
        init(trdom,basis,cutlevel);
    }

    const gsDofMapper & mapper() const { return m_mapper; }
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
    gsFunctionExpr<> impl_fun("1 - x^2 - y^2", 2)

    // uses the two objects above to get an immersed geometry
    //gsImmersedGeometry<real_t> igo(bg, inOut);

    // Background basis
    gsKnotVector<> kv(-1,1,9,2);
    gsTensorBSplineBasis<2,real_t> tbs(kv,kv);

    gsImplTrimmedDomain<real_t> tr_domain(impl_fun, tbs);


    tr_domain.tree().printLeaves();
    

    
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



