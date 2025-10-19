/** @file make_g1.cpp

    @brief 



    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

// control point index starting from halfedge (corner0,corner1) and offset (xstep,ystep)
inline index_t cpIndex(index_t HEcorner, index_t xstep, index_t ystep)
{
    static const index_t str[]  = { 1, 6 }; // degree 5
    //HEcorner a is corner id (0,1,2,3) --> (0,5,30,35)
    // Note boundary::corner equivalent id would be: 1,2,4,3
    switch (HEcorner)
    {
    case 0: // HE 0->1
        return 0 + str[0]*xstep + str[1]*ystep;
        break;
    case 1: // HE 1->3
        return 5 + str[1]*xstep - str[0]*ystep;
        break;
    case 3:  // HE 3->2
        return 35 - str[0]*xstep - str[1]*ystep;
        break;
    case 2:  // HE 2->0
        return 30 - str[1]*xstep + str[0]*ystep;
        break;
    default:
        GISMO_ERROR("indexOnPatch error: "<< HEcorner);
        break;
    }
}

//swap [os os2], [ci ci2], b[..]

class BBMapper
{
    class BBMapperHandle //Rename: BBMapperHandle
    {
    public:
        const BBMapper * m_bbm;
        gsSurfMesh::Halfedge m_he;

        gsSparseMatrix<real_t> & M;
        gsDofMapper & cmapper;
        const std::vector<size_t> & offset;
        std::pair<index_t,index_t> ci;
        std::pair<index_t,index_t> pi;
    public:

        // BBMapperHandle(const BBMapper & _bbm,
        //                gsSurfMesh::Halfedge _he)
        // : m_bbm(_bbm). m_he(_he) { }
        
        BBMapperHandle( gsSparseMatrix<real_t> & _M,
                        gsDofMapper & _cmapper, //information on the multipatch
                        const std::vector<size_t> & _offset,
                        const gsSurfMesh & msh,
                        gsSurfMesh::Halfedge h)
        : M(_M), cmapper(_cmapper), offset(_offset)
        {
            this->reset(msh, h);
        }

        void reset(const gsSurfMesh & msh,
                   gsSurfMesh::Halfedge h)
        {
            m_he = h;
            auto v = msh.from_vertex(h);
            pi.second = msh.face(h).idx();
            if (-1 != pi.second)
                cmapper.indexOnPatch(v.idx(), pi.second, ci.second );
            h = msh.opposite_halfedge(h);
            pi.first = msh.face(h).idx();
            if (-1 != pi.first)
                cmapper.indexOnPatch(v.idx(), pi.first, ci.first );
        }

        int numPatches() const {return (-1!=pi.first)+(-1!=pi.second); }
        
        /// Returns the index of the BB point (i,j) on left patch
        inline index_t left(index_t i, index_t j = 0) const
        { return offset[pi.second] + cpIndex(ci.second , i, j); }

        /// Returns the index of the BB point (i,j) on right patch
        inline index_t right(index_t i, index_t j = 0) const
        { return offset[pi.first] + cpIndex(ci.first, j, i); }

        /// Returns the index of the BB point (i,j) on right (pid=0) or left (pid=1) patch
        inline index_t idx(bool pid, index_t i, index_t j = 0) const
        { return pid ? left(i,j) : right(i,j); }

        /// Sets the coefficient (i,j) of function \a fid on the left patch to \a val
        inline void insertLeft(index_t fid, index_t i, index_t j, real_t val)
        { M.insert( left(i,j), fid ) = val; }

        /// Sets the coefficient (i,j) of function \a fid on the left patch to \a val
        inline void insertRight(index_t fid, index_t i, index_t j, real_t val)
        { M.insert( right(i,j), fid ) = val; }

        /// Sets the coefficient (i,j) of function \a fid on the left patch to \a val
        /// Indexing centered around a halfedge:
        /// j==0: BB points along the halfedge 
        /// j>0: BB points on the over (left)  of the halfedge
        /// j<0: BB points on the under (right) of the halfedge
        /// We always need i>=0
        /// It sets at once all instances of the control point
        ///  |(j+)
        ///  |
        ///  |-------------> (i)
        ///  |
        ///  |(j-)
        inline void insert(index_t fid, index_t i, index_t j, real_t val)
        {
            GISMO_ASSERT(!(i<0), "Index cannot be nagative "<< i);
            if  (j<0)
            {
                M.insert( right(i,-j), fid ) = val;
                //if (0==i) // (0,-j) -- on CW halfedge
                //    gsWarn<<"Did not set B[0,-j] twin CP.\n";
            }
            else
            {
                M.insert( left(i,j)  , fid ) = val;
                if (0==j) //on this halfedge
                {
                    if (-1 != pi.first) M.insert( right(i, 0), fid ) = val;
                    //if (0==i) // (0,0) -- on vertex
                    //    gsWarn<<"Did not set all B[0,0] twin CPs.\n";
                }

                if (0==i) // (0,j) -- on CCW halfedge
                {
                    //    gsWarn<<"Did not set B[0,j] twin CP.\n";
                }
            }
        }

        /// Inserts symmetrically control points on both patches
        /// Assumes that both left and right patches are valid
        inline void insertBoth(index_t fid, index_t i, index_t j, real_t val)
        {
            GISMO_ASSERT(!(j<0), "Index cannot be nagative "<< i);
            insert(fid, i, j, val);
            if (0!=j)
                insert(fid, i,-j, val);
        }
        
        std::vector<index_t> both(index_t i, index_t j = 0) const
        {
            if (-1 != pi.first  && -1 != pi.second) 
                return { right(i,j), left(i,j) };
            else if (-1 != pi.first)
                return { right(i,j) };
            else if (-1 != pi.second)
                return { left(i,j) };
            return {};
        }

        };//end class BBMapperHandle

    /*
        class BBMapperHandle2
        {
        gsDofMapper & cmapper;
        const std::vector<size_t> & offset;
        std::vector<index_t> ci;
        std::vector<index_t> pi;
    public:
        BBMapperHandle2(gsDofMapper & _cmapper, //information on the multipatch
              const std::vector<size_t> & _offset,
              const gsSurfMesh & msh,
              gsSurfMesh::Halfedge h)
        : cmapper(_cmapper), offset(_offset)
        {
            // Register (max) 1+8 patches around this one.
            // 0 (this)
            // 1,2,3,4 --> first neighborhoud
            // 5 6 7 8 --> second neighborhood

            pi.reserve(9);
            ci.reserve(9);
            // pi.push_back( msh.face(h).idx() );
            // if (-1 != pi.back())
            //     cmapper.indexOnPatch(msh.from_vertex(h).idx(), pi.back(), ci.back() );

            auto he = h;
            do //CCW
            {
                pi.push_back( msh.face(h).idx() );
                if (-1 != pi.back())
                    cmapper.indexOnPatch(msh.from_vertex(h).idx(), pi.back(), ci.back() );

                h = msh.ccw_rotated_halfedge(h);
                bb = bbmap.atHalfedge(h);
                }
            while (h!=he);

                            
            auto f = msh.face(h);
            for ( auto he : msh.halfedges(v) ) //for all halfedges he in msh
                              
            
            auto v = msh.from_vertex(h);
            pi.second = msh.face(h).idx();
            if (-1 != pi.second)
                cmapper.indexOnPatch(v.idx(), pi.second, ci.second );
            h = msh.opposite_halfedge(h);
            pi.first = msh.face(h).idx();
            if (-1 != pi.first)
                cmapper.indexOnPatch(v.idx(), pi.first, ci.first );
        }

        /// Returns the index of the BB point (i,j) on right (pid=0) or left (pid=1) patch
        inline index_t idx(bool pid, index_t i, index_t j) const
        { return pid ? left(i,j) : right(i,j); }

        };//end class BBMapperHandle2
    */

    std::vector<size_t>  offset;

    gsSurfMesh & msh;

    gsSurfMesh::Halfedge_property<std::pair<index_t,index_t> > ci;
                
    gsSparseMatrix<real_t> & M;//todo: gsFiberMatrix
    gsDofMapper & cmapper;
public:
    BBMapper(gsSurfMesh & _msh,
             gsSparseMatrix<real_t> & _M,
              gsDofMapper & _cmapper)
    : msh(_msh), M(_M), cmapper(_cmapper)
    {
        const size_t nf = msh.n_faces();
        offset.reserve( nf );
        offset.push_back(0);
        for (size_t k = 1; k < nf; ++k)
            offset.push_back(offset.back() + 36);

        ci = _msh.add_halfedge_property<std::pair<index_t,index_t> >("h:cindex",
                                                                     std::make_pair(-1,-1));
        index_t fid;
        gsSurfMesh::Vertex v;
        for (auto h : _msh.halfedges() )
        {
            v = msh.from_vertex(h);
            fid = msh.face(h).idx();//second
            if (-1 != fid)
                _cmapper.indexOnPatch(v.idx(), fid, ci[h].second);
            fid = msh.face(msh.opposite_halfedge(h)).idx();//first
            if (-1 != fid)
                _cmapper.indexOnPatch(v.idx(), fid, ci[h].first);
        }
    }

    ~BBMapper() { msh.remove_halfedge_property(ci); }

    inline index_t lp(gsSurfMesh::Halfedge h) const
    { return msh.face(h).idx(); }
    inline index_t l0(gsSurfMesh::Halfedge h) const
    { return ci[h].second; }
    inline index_t rp(gsSurfMesh::Halfedge h) const
    { return msh.face(msh.opposite_halfedge(h)).idx(); }
    inline index_t r0(gsSurfMesh::Halfedge h) const
    { return ci[h].first; }

    /// Returns the index of the BB point (i,j) on left patch
    inline index_t lcp(gsSurfMesh::Halfedge h, index_t i, index_t j) const
    { return offset[lp(h)] + cpIndex(l0(h), i, j); }
    /// Returns the index of the BB point (i,j) on right patch
    inline index_t rcp(gsSurfMesh::Halfedge h, index_t i, index_t j) const
    { return offset[rp(h)] + cpIndex(r0(h), j, i);/*note: (j,i) on purpose!*/ }

    inline void insert(gsSurfMesh::Halfedge h, index_t fid, index_t i, index_t j, real_t val)
    {
        GISMO_ASSERT(i>=0, "Index cannot be negative "<< i);
        if  (j<0)
        {
            M.insert( rcp(h,i,-j), fid ) = val; // OK
            if (0==i) // (0,-j) -- on CW halfedge
            {
                auto cwh = msh.cw_rotated_halfedge(h);
                if (-1 != rp(cwh))
                    M.insert( rcp(cwh,-j,0), fid ) = val;
            }
        }
        else//j>=0
        {
            if (0==j) //on this halfedge
            {
                if (0==i) // (0,0) -- on vertex
                {
                    for ( auto vh : msh.halfedges(msh.from_vertex(h)) )
                        if (-1 != lp(vh))
                            M.insert( lcp(vh,0,0), fid) = val;
                    return;
                }
                else if (-1 != rp(h)) //edge twin
                {
                    M.insert( rcp(h, i,0), fid ) = val;
                }
            }

            M.insert( lcp(h,i,j), fid ) = val; // OK

            if (0==i) // (0,j) -- on CCW halfedge
            {
                auto ccwh = msh.ccw_rotated_halfedge(h);
                if (-1 != lp(ccwh))
                    M.insert( lcp(ccwh,j,0), fid ) = val;
            }
        }
    }

    inline void insert(gsSurfMesh::Vertex v, index_t fid, real_t val)
    {
        for ( auto h : msh.halfedges(v) )
            M.insert( lcp(h,0,0), fid ) = val;
    }

    /// Inserts symmetrically control points on both patches
    /// Assumes that both left and right patches are valid
    inline void insertBoth(gsSurfMesh::Halfedge h, index_t fid,
                           index_t i, index_t j, real_t val)
    {
        GISMO_ASSERT(!(j<0), "Index cannot be nagative "<< i);
        insert(h, fid, i, j, val);
        if (0!=j)
            insert(h, fid, i,-j, val);
    }

    //--------------
    
    BBMapperHandle atHalfedge(gsSurfMesh::Halfedge h) const
        { return BBMapperHandle(M, cmapper, offset, msh, h); }
    //{ return BBMapperHandle(*this, h); }

    // Returns all corner indices associated to vertex \a v
    std::vector<index_t> allCorners(gsSurfMesh::Vertex v) const
    {
        std::vector<index_t> res;
        index_t ci;
        res.reserve( msh.valence(v));
        for ( auto ff : msh.faces(v) )
        {
            const index_t pi = ff.idx();
            cmapper.indexOnPatch(v.idx(), pi, ci );
            res.push_back( offset[pi] + cpIndex(ci, 0, 0));
        }
        return res;
    }

};// end class BBMapper

gsMultiPatch<> mesh_to_multipatch(const gsSurfMesh & msh)
{
    gsKnotVector<> kv(0,1,0,2);
    gsTensorBSplineBasis<2> bs(kv,kv);
    gsMatrix<> coefs;
    int cc;
    gsMultiPatch<> result;
    auto points = msh.get_vertex_property<gsSurfMesh::Point>("v:point");
    for ( auto ff : msh.faces() ) //for all faces
    {
        coefs.resize(4, 3);
        cc = 0;
        for (auto v : msh.vertices(ff)) //for all vertices
        {
            coefs.row(cc<2 ? cc : 5-cc) = points[v].transpose();
            ++cc;
        }
        result.addPatch( bs.makeGeometry(give(coefs)) );
    }
    result.computeTopology();
    return result;
}


bool on_xy_plane(gsMultiPatch<> & mp)
{
    index_t n = mp.geoDim();
    if ( n<3 )
        return true;
    n = n - 2;
    for (size_t q = 0; q < mp.nPatches(); ++q)
        if ( (mp.patch(q).coefs().rightCols(n).array()!=0).any() )
            return false;
    return true;
}

// Toepliz of (0,1,0,..0)
gsMatrix<> circulant(index_t n)
{
    gsMatrix<> C(n,n);
    C.setZero();
    C(n-1,0) = (real_t)1;
    for(index_t i = 1; i<n;++i)
        C(i-1,i) = (real_t)1;
    return C;
}

bool checkG1(gsMultiPatch<> & mp, bool verbose = false) // BUG on boundary ??
{
    gsMultiBasis<> basis(mp);
    gsExprEvaluator<> ev;
    ev.setIntegrationDomain(basis.domain());
    auto G = ev.getMap(mp);

    if (verbose)
        gsInfo  <<"\nG : ["<<mp.parDim()<<"] --> ["<<mp.geoDim()<< "]\n";
        
    ev.integralInterface( ( G.left() - G.right() ).sqNorm() );
    ev.calcSqrt();
    bool ok = (ev.value()<1e-3);
    if (verbose && ok) gsInfo  <<"\nSurface is C0.";
    if (!ok) gsInfo  <<"\nResult (C0): "<< ev.value() <<"\n";
    if (verbose)
    {
        gsInfo<<"Per edge: "<<ev.allValues().transpose() <<"\n";//<<"\n";
    }

    ev.integralInterface( (sn(G.left()).normalized()-sn(G.right()).normalized()).sqNorm() );
    ev.calcSqrt();
    ok &= (ev.value()<1e-3);
    if (verbose && ok) gsInfo  <<"\nSurface is G1.";
    if (!ok) gsInfo  <<"\nResult (G1): "<< ev.value() <<"\n";
    if (verbose)
    {
        gsInfo<<"Per edge: "<<ev.allValues().transpose() <<"\n";//<<"\n";
    }

    /*
    ev.integralInterface( (sn(G.left())-sn(G.right())).sqNorm() );
    ev.calcSqrt();
    ok &= (ev.value()<1e-3);
    if (verbose &&ok) gsInfo  <<"\nSurface has continuous normal.";
    if (verbose)
        gsInfo  <<"\nResult (non-unit normal): "<< ev.value() <<"\n";
    */

    /*
    ev.maxInterface( abs( (fform(G.left() ).inv()*fform2nd(G.left() )).det() -
                          (fform(G.right()).inv()*fform2nd(G.right())).det() ) );
    ok &= (ev.value()<1e-3);
    if (ok) gsInfo  <<"\nSurface has continuous Gauss curvature.";
    gsInfo  <<"\nResult (Gauss curv. diff): "<< ev.value() <<"\n";
    if (verbose)
        gsInfo<<"Per edge: "<<ev.allValues().transpose() <<"\n";//<<"\n";




    ev.maxInterface( abs( (fform(G.left() ).inv()*fform2nd(G.left() )).trace().val() -
                          (fform(G.right()).inv()*fform2nd(G.right())).trace().val() ) );
    ok &= (ev.value()<1e-3);
    if (ok) gsInfo  <<"\nSurface has continuous mean curvature.";
    gsInfo  <<"\nResult (mean curv. diff): "<< ev.value() <<"\n";
    if (verbose)
        gsInfo<<"Per edge: "<<ev.allValues().transpose() <<"\n";//<<"\n";


    ev.minInterface( (sn(G.left()).normalized()+sn(G.right()).normalized()).norm()/2.0 );
    if (ev.value()<1e-3) gsInfo <<"\nSurface is NOT oriented properly.";
    else gsInfo <<"\nSurface is oriented.";
    gsInfo  <<"\nResult (orient): "<< ev.value() <<"\n";
    if (verbose)
        gsInfo<<"Per edge: "<<ev.allValues().transpose() <<"\n";//<<"\n";

    */
/*
    //test vertices
    ev.options().setReal("quA",0);
    ev.options().setInt ("quB",2);
    ev.options().addInt ("quRule","",gsQuadrature::GaussLobatto); //Gauss-Lobatto

    ev.maxInterface( (G.left()-G.right()).norm() );
    gsInfo  <<"\nResult (vertex values diff): "<< ev.value() <<"\n";
    if (verbose)
        gsInfo<<"Per edge: "<<ev.allValues().transpose() <<"\n";//<<"\n";

    ev.maxInterface( (sn(G.left()).normalized()-sn(G.right()).normalized()).norm());
    gsInfo  <<"Result (vertex unit normals diff): "<< ev.value() <<"\n";
    if (verbose) gsInfo<<"Per edge: "<<ev.allValues().transpose() <<"\n";

    ev.maxInterface( (sn(G.left())-sn(G.right()) ).norm() );
    gsInfo  <<"Result (vertex normals diff): "<< ev.value() <<"\n";
    if (verbose) gsInfo<<"Per edge: "<<ev.allValues().transpose() <<"\n";
*/
    return ok;
}

void PlotBasis(const gsSparseMatrix<> & cf, const gsMultiPatch<> & mp,
               const int & ns,
               const bool pm = false, bool pn = false)
{
    //gsFileData<> xml;
    gsParaviewCollection col("ccbasis");
    gsMultiPatch<> mmp = mp;
    index_t d = mp.geoDim();
    if ( on_xy_plane(mmp) && 3==d) --d;
    gsInfo << "dim " <<d<<"\n";
    const size_t np = mp.nPatches();
    gsVector<index_t> offset(np+1);
    offset[0] = 0;
    mmp.patch(0).embed(d+1);
    for ( size_t i = 1; i< np; ++i)
    {
        offset[i] = offset[i-1] + mp.patch(i-1).basis().size();
        mmp.patch(i).embed(d+1);
    }
    offset[np] = offset[np-1] + mp.patch(np-1).basis().size();

    for ( index_t k = 0; k != cf.cols(); ++k) // for all basis function
    {
        //if (cf.col(k).nonZeros()==1) continue; //skip free functions
        for ( size_t i = 0; i< np; ++i)
            mmp.patch(i).coefs().col(d).setZero();

        for(auto it = cf.begin(k); it; ++it)
        {
            index_t pid = std::upper_bound(offset.begin(),offset.end(),it.row())
                - offset.begin() - 1;
            index_t ii  = it.row() - offset[pid];
            mmp.patch(pid).coefs()(ii, d) = it.value();
        }

        // Export XML
        //xml << mmp;
        
        //--- Check that the function is G1
        if (2==d && !checkG1(mmp, false))
            gsInfo <<"----------- Basis function number "<< k <<" FAILED the test ------------\n";

        if (d>2) pn = false; // no control net in 4D
        
        // Export Paraview
        std::string gcBasisFct = "ccbasis" + util::to_string(k);
        gsWriteParaview(mmp, gcBasisFct, ns, pm, pn);
        for ( size_t i = 0; i< np; ++i)
            col.addTimestep(gcBasisFct, i, k, ".vts");
        if(pm)
        {
            for ( size_t i = 0; i< np; ++i)
                col.addTimestep(gcBasisFct, i, k, "_mesh.vtp");
        }
        if(pn)
        {
            for ( size_t i = 0; i< np; ++i)
            col.addTimestep(gcBasisFct, i, k, "_cnet.vtp");
        }
    }
    col.save();
    //xml.save("representative_basis");
}

inline unsigned int vertex_degree(const gsSurfMesh & msh, gsSurfMesh::Vertex v)
{
    unsigned int k = 0;
    for (auto f: msh.faces(v)) { ++k; GISMO_UNUSED(f); }
    return msh.is_boundary(v) ? 2*k : k;
}

index_t g1Dimension(const gsSurfMesh & msh)
{
    index_t n, c1dim(0);
    for (auto v : msh.vertices())
    {
        n = msh.valence(v);
        if (n<=3 && msh.is_boundary(v) ) // boundary regular (T-junction) or corner
            c1dim += 4;
        if (n==4 && !msh.is_boundary(v)) // inner regular
            c1dim += 4;
        if (n!=4 && !msh.is_boundary(v)) // inner EV
            c1dim += n + 3;
        if (n>3 && msh.is_boundary(v))   // boundary EV
            c1dim += 2;
    }
    gsInfo << "Dimension V  "<< c1dim <<"\n";

    for ( auto ed : msh.edges() ) //for all edges he in msh
    {
        bool bdr = msh.is_boundary(ed);
        // Note for Extraordinary edge on boundary we still count 4)

        bool ee = !( (4==vertex_degree(msh, msh.vertex(ed,0)) ||
                      2==vertex_degree(msh, msh.vertex(ed,0)) ) &&
                     (4==vertex_degree(msh, msh.vertex(ed,1)) ||
                      2==vertex_degree(msh, msh.vertex(ed,1)) )  );
        // 4*boundary + 4*inner_regular + 2*inner_EE
        c1dim += 4 * bdr  + 4*( !bdr && !ee) + 2*(!bdr && ee);
    }
    gsInfo << "Dimension V+E: "<< c1dim <<"\n";
    
    c1dim += 4 * msh.n_faces();
    gsInfo << "Dimension is (final): "<< c1dim <<"\n";
    return c1dim;
}

// Matrices:
// C   = Circ()^T
// C1s = -a*I + C + C^{n-1}
// K1  = Ker(C1s) and this kernel has dimension 2
// C1  = [C1s; K1^T ]
// sol1 = C1 \ b  where b = [ (2-a); .. ;(2-a) ; 0; 0 ]

// For Vertex value BF:
// C2 = [C + I; Ker(C+I)] where kernel has dimension either 0 (for odd) or 1 (for even)
//sol2 = C2 \ b where b = (1/5) * ( vec(a) - 5*(a-2)*sol + 4*a*( (1/2)*sol - 1/10 ) )

// For Derivatives BFs:
//sol2 = C2 \ b where b = (1/5) * ( -5*(a-2)* K1 + 2*a*K1 )
void g1Matrices(index_t n, real_t & a0,
                gsMatrix<> & C1,
                gsMatrix<> & C2)
{
    a0 = 2 * math::cos(2*EIGEN_PI/n);
    gsMatrix<> C = circulant(n);
    C.transposeInPlace(); // Making column-circulant !!
    //gsInfo << "C:\n" << C <<"\n";

    C1.resize(n+2,n);
    auto C1s = C1.topRows(n);
    C1s = C;
    for(index_t i = 2; i<n;++i) // C^{n-1}
        C1s *= C;
    C1s += C;
    C1s.diagonal().array() -= a0;
    typename gsMatrix<>::JacobiSVD svd = C1s.jacobiSvd(gsEigen::ComputeFullV);
    C1.bottomRows(2) = svd.matrixV().rightCols(2).transpose();
    //gsInfo << "C1:\n" << C1 <<"\n";
    //Note: auto K1 = C1.bottomRows(2).transpose();

    index_t sz = n + (n%2==0);
    C2.resize(sz,n);
    C2.topRows(n) = C;
    C2 /*.topRows(n)*/ .diagonal().array() += 1;
    if (sz!=n) // Case: n even
    {
        //add to C2 one last row: (1,N): [-1 1 -1 1 -1 1 ...]
        // C2 has now size: (N+1)xN
        auto RR = C2.row(n);
        index_t tmp(1);
        for(index_t q = 0; q!=n; ++q, tmp *= -1)
            RR[q]= tmp;
    }
    //gsInfo << "C2:\n" << C2 <<"\n";
}

gsSparseMatrix<> basisFromMesh(gsSurfMesh & msh, gsDofMapper & cmapper)
{
    //todo: version without cmapper, to create basis from a pure mesh
    const size_t nf = msh.n_faces();
    static const size_t cpSize = 36;
    // Initialize offsets
    std::vector<size_t> offset;
    offset.reserve( nf );
    offset.push_back(0);
    for (size_t k = 1; k < nf; ++k)
        offset.push_back(offset.back() + cpSize);

    // Get dimension
    int n, c1dim = g1Dimension(msh);

    //Create transfer matrix
    gsSparseMatrix<> M(offset.back()+cpSize, c1dim);
    index_t pi, pi2, os, os2, ci, ci2;
    c1dim = 0;//reset

    //for neighbor
    gsSurfMesh::Vertex v_N;
    real_t a_N;
    int n_N;

    real_t a;
    gsMatrix<> C1, C2, b;

    BBMapper bbmap(msh, M, cmapper);

    for (auto v : msh.vertices()) //for all vertices
    {
        n = msh.valence(v);
        if (n==2) //corner vertex - we assume msh.is_boundary(v) [manifold mesh] $$4
        {
            //note: msh.halfedge(v) returns a boundary half-edge
            auto h = msh.next_halfedge(msh.opposite_halfedge(msh.halfedge(v)));
            auto bb = bbmap.atHalfedge(h);
            // four "corner functions" to add (square)
            // bb.insert(c1dim++, 0, 0, 1.0);
            // bb.insert(c1dim++, 0, 1, 1.0);
            // bb.insert(c1dim++, 1, 0, 1.0);
            // bb.insert(c1dim++, 1, 1, 1.0);

            bbmap.insert(h, c1dim++, 0, 0, 1.0);
            bbmap.insert(h, c1dim++, 0, 1, 1.0);
            bbmap.insert(h, c1dim++, 1, 0, 1.0);
            bbmap.insert(h, c1dim++, 1, 1, 1.0);

        }

        if (n==3 && msh.is_boundary(v)) // boundary regular (T-junction) $$4
        {
            // interior halfedge of the T-junction
            auto h = msh.cw_rotated_halfedge(msh.halfedge(v));
            auto bb = bbmap.atHalfedge(h);

            for (index_t j : {-1,1})
            {
                //Anchor (01)
                bbmap.insert(h, c1dim, 0,j , 1.0);
                bbmap.insert(h, c1dim, 0,0 , 0.5);
                ++c1dim;
                // Anchor (11)
                bbmap.insert(h, c1dim, 1,j , 1.0);
                bbmap.insert(h, c1dim, 1,0 , 0.5);
                ++c1dim;
            }
            /*
            for (bool i : {0,1})
            {
                //Anchor (01)
                M.insert(bb.idx( i, 0, 1), c1dim  ) = (real_t)1.0;
                M.insert(bb.idx( i, 0, 0), c1dim  ) = (real_t)0.5;// 1/2
                M.insert(bb.idx(!i, 0, 0), c1dim  ) = (real_t)0.5;
                //for ( index_t q : bb.both(0) )
                //    M.insert(q, c1dim) = (real_t)0.5;
                ++c1dim;
                // Anchor (11)
                M.insert(bb.idx( i, 1, 1), c1dim  ) = (real_t)1.0;
                M.insert(bb.idx( i, 1, 0), c1dim  ) = (real_t)0.5;// 1/2
                M.insert(bb.idx(!i, 1, 0), c1dim  ) = (real_t)0.5;
                //for ( index_t q : bb.both(1) )
                //    M.insert(q, c1dim) = (real_t)0.5;
                ++c1dim;
            }
            */
            continue;// no modif
            
            v_N = msh.to_vertex(h);
            n_N = msh.valence(v_N);
            if (false && 4 != n_N && !msh.is_boundary(v_N) )
            {
                a_N = 2 * math::cos(2*EIGEN_PI/n_N);
                real_t b20, tmp;
                for (index_t i = c1dim-4; i<c1dim; ++i)
                {
                    b20 = 0.5 * M(bb.left(1,0), i) - 0.1*M(bb.left(0,0), i); // 1/2 * b10 - 1/10 * b00
                    M.insert(bb.left (2,0), i) = b20;
                    M.insert(bb.right(2,0), i) = b20;
                    tmp = 0.5*(2-a_N)*b20 + 0.2*M(bb.left(1, 0),i);//1/2 *(2-a) * b20 + 1/5*a_N *b10
                    M.insert(bb.right(2,1), i) = tmp;
                    M.insert(bb.left (2,1), i) = tmp;
                    tmp = 0.3 * a_N * b20;// 3/10 * a_N * b20
                    M.insert(bb.left (3,1), i) = tmp;
                    M.insert(bb.right(3,1), i) = tmp;
                }
            }
        }

        if (n==4 && !msh.is_boundary(v) ) // inner regular $$4
        {
            for ( auto he : msh.halfedges(v) ) // iterate over all HE (n=4 functions)
            {
                //corner
                bbmap.insert(he,c1dim,0,0, 0.25);
                //for ( index_t q : bbmap.allCorners( msh.from_vertex(he)) )
                //    M.insert(q, c1dim) = (real_t)0.25;
                auto bb = bbmap.atHalfedge(he);
                //first axis
                bbmap.insert(he,c1dim, 1,1 , 1.0);
                bbmap.insert(he,c1dim, 1,0 , 0.5);
                bbmap.insert(he,c1dim, 0,1 , 0.5);

                // M.insert(bb.left (1, 1), c1dim) = (real_t)1.0;
                // M.insert(bb.left (1, 0), c1dim) = (real_t)0.5;
                // M.insert(bb.right(1, 0), c1dim) = (real_t)0.5;     
                auto bb1 = bbmap.atHalfedge(msh.ccw_rotated_halfedge(he));
                //second axis
//                bb1.insert(c1dim, 1,0 , 0.5);
//                bbmap.insert(msh.ccw_rotated_halfedge(he), c1dim, 1,0 , 0.5);
                ++c1dim;
            }
// /*
            index_t hc = 0;
            for ( auto h : msh.halfedges(v) ) // iterate again
            {
                v_N = msh.to_vertex(h);
                n_N = msh.valence(v_N);
                if ( 4 != n_N && !msh.is_boundary(v_N) )
                {
                    a_N = 2 * math::cos(2*EIGEN_PI/n_N);
                    auto bb = bbmap.atHalfedge(h);
                    for (index_t i = 0; i!=4; ++i)
                    {
                        const index_t cc = c1dim - 4 + i;
                       if ( (i==hc) || (i+1)%4==hc ) //HE adjacent function cc
                       {
                           bbmap.insert    (h,cc, 2,0 , 0.225);//=9/40
                           bbmap.insertBoth(h,cc, 2,1 , 0.225 - a_N/80);
                           bbmap.insertBoth(h,cc, 3,1 , 0.0675*a_N);//=27/400*a_N

                            // M.insert(bb.left (2, 0), cc) = (real_t)0.225;//=9/40;
                            // M.insert(bb.right(2, 0), cc) = (real_t)0.225;//=9/40;
                            // M.insert(bb.left (2, 1), cc) = (real_t)(0.225 - a_N/80);
                            // M.insert(bb.right(2, 1), cc) = (real_t)(0.225 - a_N/80);
                            // M.insert(bb.left (3, 1), cc) = (real_t)(0.0675*a_N);//=27/400*a_N
                            // M.insert(bb.right(3, 1), cc) = (real_t)(0.0675*a_N);//=27/400*a_N
                        }
                        else
                        {
                            bbmap.insert    (h,cc, 2,0 , -0.025);//=-1/40
                            bbmap.insertBoth(h,cc, 2,1 , -(2-a_N)/80);
                            bbmap.insertBoth(h,cc, 3,1 , -0.0075*a_N);//=-3/400*a_N
                           
                            // M.insert(bb.left (2,0), cc) = (real_t)(-0.025);//=-1/40
                            // M.insert(bb.right(2,0), cc) = (real_t)(-0.025);//=-1/40
                            // M.insert(bb.left (2,1), cc) = -(2-a_N)/80;
                            // M.insert(bb.right(2,1), cc) = -(2-a_N)/80;
                            // M.insert(bb.left (3,1), cc) = -0.0075*a_N;//=-3/400*a_N
                            // M.insert(bb.right(3,1), cc) = -0.0075*a_N;//=-3/400*a_N
                        }                            
                    }
                }
                ++hc;
            }
        }

        if (n!=4 && !msh.is_boundary(v)) //inner EV $$(n+3 functions)
        {
            //gsInfo << "EV: " << c1dim <<"\n";
            g1Matrices(n, a, C1, C2);
            auto K1 = C1.bottomRows(2).transpose();
            //gsInfo << "K1:\n" << K1 <<"\n";

            b.resize(n+2,1);
            b.topRows(n).setConstant(2 - a);
            b.bottomRows(2).setZero();
            gsMatrix<> sol1 = C1.fullPivLu().solve(b);
            //gsInfo << "sol1: " << sol1.transpose() <<"\n";

            b.setZero(C2.rows(), 1);
            // b.topRows(n).setConstant(a);
            // b +=-5*(a-2)*sol1 + 4 * a * (0.5*sol1 - 0.1*gsMatrix<>::Ones(n,1)) ;
            // b *= 0.2;
            for(index_t q=0;q!=n;++q)
                b.at(q) = 0.2 * ( a-5*(a-2)*sol1.at(q) + 4*a*(0.5*sol1.at(q)-0.1) );

              gsMatrix<> sol2 = C2.fullPivLu().solve(b);
              //gsInfo << "sol2: " << sol2.transpose() <<"\n";
              //gsInfo << "b:\n" << b.transpose() <<"\n";

              // Basis function attached to the vertex value

              // for( auto vv : bbmap.allCorners(v) )
              //    M.insert(vv, c1dim) = (real_t)1.0;
              // for ( auto ff : msh.faces(v) ) //for all faces, set vertex value
              // {
              //     pi = ff.idx();
              //     cmapper.indexOnPatch( v.idx(), pi, ci);
              //     os = offset[pi];
              //     M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)1.0;
              // }
              // NEW: ?
              bbmap.insert(v, c1dim, 1.0);

              index_t k = 0;// (continues) Basis function attached to the vertex value
              for ( auto he : msh.halfedges(v) ) //for all halfedges he in msh
              {
                  real_t R1 = 0.1 * ( -a + 5*a*sol1.at(k) - 10*(a-2)*(0.5*sol1.at(k)-0.1));

// /*
                  auto bb = bbmap.atHalfedge(he);
                  bbmap.insert    (he,c1dim,  1, 1, sol2.at(k)           );
                  bbmap.insert(he,c1dim,  0, 1, sol1.at(k)           );
                  bbmap.insert(he,c1dim,  0, 2, 0.5*sol1.at(k) - 0.1 );
                  bbmap.insertBoth(he,c1dim,  1, 2, 0.5*R1               );

                  // M.insert(bb.right(1, 1), c1dim) = sol2.at(k);
                  // //
                  // M.insert(bb.right(0, 1), c1dim) = sol1.at(k);
                  // M.insert(bb.right(0, 2), c1dim) = 0.5*sol1.at(k) - 0.1;
                  // M.insert(bb.right(1, 2), c1dim) = 0.5*R1;
                  // //
                  // M.insert(bb.left (0, 1), c1dim) = sol1.at(k);
                  // M.insert(bb.left (0, 2), c1dim) = 0.5*sol1.at(k) - 0.1;
                  // M.insert(bb.left (1, 2), c1dim) = 0.5*R1;
                  ++k;
                  continue;
//*/
                  //--------------------------------------- //known to work

                  pi = msh.face(he).idx();
                  os = offset[pi];
                  cmapper.indexOnPatch(v.idx(), pi, ci);
                  M.insert(os+cpIndex(ci, 1, 1), c1dim) = sol2.at(k);
                  //
                  M.insert(os+cpIndex(ci, 1, 0), c1dim) = sol1.at(k);
                  M.insert(os+cpIndex(ci, 2, 0), c1dim) = 0.5*sol1.at(k) - 0.1;
                  M.insert(os+cpIndex(ci, 2, 1), c1dim) = 0.5*R1;
                  pi = msh.face(msh.opposite_halfedge(he)).idx();
                  cmapper.indexOnPatch(v.idx(), pi, ci);
                  os = offset[pi];
                  M.insert(os+cpIndex(ci, 0, 1), c1dim) = sol1.at(k);
                  M.insert(os+cpIndex(ci, 0, 2), c1dim) = 0.5*sol1.at(k) - 0.1;
                  M.insert(os+cpIndex(ci, 1, 2), c1dim) = 0.5*R1;
                  ++k;
              }
              //gsInfo <<"Function "<< c1dim <<" :\n" << M.col(c1dim).toDense().transpose() <<"\n";
              ++c1dim;

              //Start: basis functions attached to the partial derivatives at the vertex
              b.setZero(C2.rows(),2); // reusing memory for b
              //b += -5*(a-2)*K1 + 4 * a * (K1 * 0.5) ;
              //b *= 0.2;
              for(index_t q=0;q!=n;++q)
              {
                  b(q,0) = 0.2 * (-5*(a-2)*K1(q,0) + 2*a*K1(q,0));
                  b(q,1) = 0.2 * (-5*(a-2)*K1(q,1) + 2*a*K1(q,1));
              }

              sol2 = C2.fullPivLu().solve(b);
              //gsInfo << "sol2: "<< sol2.transpose() <<"\n";

              //Basis functions attached to the two partial derivatives at the vertex.
              k = 0;
              for ( auto he : msh.halfedges(v) ) //for all halfedges he in msh
              {
                  real_t R1 = 0.1 * ( 5*a*K1(k,0) - 10*(a-2)*(0.5*K1(k,0)));
                  real_t R3 = 0.1 * (-5*a*K1(k,0) + 10*a*(0.5*K1(k,0)));
                  //
                  real_t R2 = 0.1 * ( 5*a*K1(k,1) - 10*(a-2)*(0.5*K1(k,1)));
                  real_t R4 = 0.1 * (-5*a*K1(k,1) + 10*a*(0.5*K1(k,1)));
// /*
                  auto bb = bbmap.atHalfedge(he);
                  //first
                  bbmap.insert(he,c1dim, 1,1, sol2(k,0) );
                  bbmap.insert(he,c1dim, 1,0, K1(k,0));
                  bbmap.insert(he,c1dim, 2,0, K1(k,0) * 0.5);
                  bbmap.insertBoth(he,c1dim, 2,1, 0.5 * R1);
                  bbmap.insertBoth(he,c1dim, 3,1, 0.5 * R3);

                  //second
                  bbmap.insert(he,1+c1dim,     1,1, sol2(k,1));
                  bbmap.insert(he,1+c1dim,     1,0, K1(k,1));
                  bbmap.insert(he,1+c1dim,     2,0, K1(k,1) * 0.5);
                  bbmap.insertBoth(he,1+c1dim, 2,1, 0.5 * R2);
                  bbmap.insertBoth(he,1+c1dim, 3,1, 0.5 * R4);

//                  bbmap.insert(he,1+c1dim, 0,-1, K1(k,1));

                  ++k;
                  continue;
                  //--------------------------------------- //known to work
//*/
                  pi = msh.face(he).idx();
                  os = offset[pi];
                  cmapper.indexOnPatch(v.idx(), pi, ci);
                  //first
                   M.insert(os+cpIndex(ci, 1, 0), c1dim) = K1(k,0);
                   M.insert(os+cpIndex(ci, 2, 0), c1dim) = K1(k,0) * 0.5;
                   M.insert(os+cpIndex(ci, 2, 1), c1dim) = 0.5 * R1;
                   M.insert(os+cpIndex(ci, 3, 1), c1dim) = 0.5 * R3;//
                  
                   M.insert(os+cpIndex(ci, 1, 1), c1dim) = sol2(k,0);
                  //second
                  //M.insert(os+cpIndex(ci, 1, 0), 1+c1dim) = K1(k,1);
                  //M.insert(os+cpIndex(ci, 2, 0), 1+c1dim) = K1(k,1) * 0.5;
                  //M.insert(os+cpIndex(ci, 2, 1), 1+c1dim) = 0.5 * R2;//..
                  //M.insert(os+cpIndex(ci, 3, 1), 1+c1dim) = 0.5 * R4;//..
                 //M.insert(os+cpIndex(ci, 1, 1), 1+c1dim) = sol2(k,1);
                  pi = msh.face(msh.opposite_halfedge(he)).idx();
                  cmapper.indexOnPatch(v.idx(), pi, ci);
                  os = offset[pi];
                  //first
                   M.insert(os+cpIndex(ci, 0, 1), c1dim) = K1(k,0);
                   M.insert(os+cpIndex(ci, 0, 2), c1dim) = K1(k,0) * 0.5;
                   M.insert(os+cpIndex(ci, 1, 2), c1dim) = 0.5 * R1;
                   M.insert(os+cpIndex(ci, 1, 3), c1dim) = 0.5 * R3;//..
                  //second
                  //M.insert(os+cpIndex(ci, 0, 1), 1+c1dim) = K1(k,1);
                  //M.insert(os+cpIndex(ci, 0, 2), 1+c1dim) = K1(k,1) * 0.5;
                  //M.insert(os+cpIndex(ci, 1, 2), 1+c1dim) = 0.5 * R2;//..
                  //M.insert(os+cpIndex(ci, 1, 3), 1+c1dim) = 0.5 * R4;//..
                  ++k;
              }
              c1dim+=2;

              // Basis attached to the cross derivatives at the vertex
              k = 0;
              for ( auto he : msh.halfedges(v) ) //for all halfedges he in msh
              {
                  real_t R1 = -(a-2)*(5.0/(4*a)) + (3.0/5.0)*a*(5.0/(4*a));
                  real_t R2 =      a*(5.0/(4*a)) - (a-2)*(5.0/(4*a));
// /*
                  auto bb = bbmap.atHalfedge(he);

                  // current f
                  bbmap.insert(he,c1dim, 1, 1, 1.0       );
                  bbmap.insert(he,c1dim, 2, 0, 5.0/(4*a) );
                  bbmap.insert(he,c1dim, 0, 2, 5.0/(4*a) );

                  bbmap.insert(he,c1dim, 3, 0, 5.0/(4*a) );
                  bbmap.insert(he,c1dim, 0, 3, 5.0/(4*a) );

                  bbmap.insert(he,c1dim, 1, 2, 0.5*R1    );
                  bbmap.insert(he,c1dim, 2, 1, 0.5*R1    );
                  bbmap.insert(he,c1dim, 2,-1, 0.5*R1    );//right face

                  bbmap.insert(he,c1dim, 1, 3, 0.5*R2    );
                  bbmap.insert(he,c1dim, 3, 1, 0.5*R2    );
                  bbmap.insert(he,c1dim, 3,-1, 0.5*R2    );//right face
                  // next face
                  auto hh = msh.ccw_rotated_halfedge(he);
                  bbmap.insert(hh,c1dim, 2,1, 0.5*R1    );
                  bbmap.insert(hh,c1dim, 3,1, 0.5*R2    );

                  if (3==n)
                      M.col(c1dim) *= -1;
                      
                  ++k;
                  ++c1dim;
                  continue;
//*/
//--------------------------------------- //known to work

                  pi = msh.face(he).idx();
                  os = offset[pi];
                  cmapper.indexOnPatch(v.idx(), pi, ci); //current face F
                  M.insert(os+cpIndex(ci, 1, 1), c1dim) = 1;
                  M.insert(os+cpIndex(ci, 2, 0), c1dim) =  5.0/(4*a);
                  M.insert(os+cpIndex(ci, 0, 2), c1dim) =  5.0/(4*a);
                  M.insert(os+cpIndex(ci, 3, 0), c1dim) =  5.0/(4*a);
                  M.insert(os+cpIndex(ci, 0, 3), c1dim) =  5.0/(4*a);

                  M.insert(os+cpIndex(ci, 1, 2), c1dim) =  0.5*R1;
                  M.insert(os+cpIndex(ci, 2, 1), c1dim) =  0.5*R1;
                  M.insert(os+cpIndex(ci, 1, 3), c1dim) =  0.5*R2;
                  M.insert(os+cpIndex(ci, 3, 1), c1dim) =  0.5*R2;
                  //
                  pi = msh.face(msh.opposite_halfedge(he)).idx();
                  cmapper.indexOnPatch(v.idx(), pi, ci); // previous face F-1
                  os = offset[pi];
                  M.insert(os+cpIndex(ci, 0, 2), c1dim) =  5.0/(4*a);
                  M.insert(os+cpIndex(ci, 0, 3), c1dim) =  5.0/(4*a);

                  M.insert(os+cpIndex(ci, 1, 2), c1dim) =  0.5*R1;
                  M.insert(os+cpIndex(ci, 1, 3), c1dim) =  0.5*R2;
                  //
                  pi = msh.face(msh.ccw_rotated_halfedge(he)).idx();
                  cmapper.indexOnPatch(v.idx(), pi, ci); // next face F+1
                  os = offset[pi];
                  M.insert(os+cpIndex(ci, 2, 0), c1dim) =  5.0/(4*a);
                  M.insert(os+cpIndex(ci, 3, 0), c1dim) =  5.0/(4*a);

                  M.insert(os+cpIndex(ci, 2, 1), c1dim) = 0.5*R1;
                  M.insert(os+cpIndex(ci, 3, 1), c1dim) = 0.5*R2;

                  if (3==n)
                      M.col(c1dim) *= -1;
                      
                  ++k;
                  ++c1dim;
              }
          }

          if (n>3 && msh.is_boundary(v)) // boundary EV
          {
              GISMO_ERROR("Boundary EV problem.");
          }

    }//end for vertices

    gsInfo << "--> Dimension V  "<< c1dim <<"\n";
        
    for ( auto he : msh.halfedges() ) //for all halfedges he in msh $$ 1(EV) or 2(RV) or 4(bdr))
    {
        if ( msh.is_boundary(he) ) continue;
        // Edge on boundary ?
        bool bdr = msh.touches_boundary(he);
        // Extraordinary edge? (note: no new functions for boundart EHE)
        bool ee = !( (4==vertex_degree(msh, msh.from_vertex(he)) ||
                      2==vertex_degree(msh, msh.from_vertex(he)) ) &&
                     (4==vertex_degree(msh, msh.to_vertex(he)) ||
                      2==vertex_degree(msh, msh.to_vertex(he)) )  );

        pi = msh.face(he).idx(); // patch
        os = offset[pi];
        cmapper.indexOnPatch( msh.from_vertex(he).idx(), pi, ci); // corner id on patch

        if (!bdr && ee) // inner EHE
        {
            pi2 = msh.face(msh.opposite_halfedge(he)).idx();
            cmapper.indexOnPatch( msh.to_vertex(he).idx(), pi2, ci2 );
            os2 = offset[pi2];
            //two functions per edge (or one per halfedge)
            bbmap.insert(he,c1dim, 2, 1,  1.0);
            bbmap.insert(he,c1dim, 2,-1, -1.0);

            // M.insert(os +cpIndex(ci , 2, 1),c1dim) = (real_t)1.0;
            // M.insert(os2+cpIndex(ci2, 3, 1),c1dim) = (real_t)(-1.0);
            ++c1dim;
        }
        if (!bdr && !ee) // inner regular HE
        {
            bbmap.insert(he,c1dim  , 2,0, 0.5);
            bbmap.insert(he,c1dim++, 2,1, 1.0);

            bbmap.insert(he,c1dim  , 3, 0, 0.5);
            bbmap.insert(he,c1dim++, 3, 1, 1.0);

            /*
            pi2 = msh.face(msh.opposite_halfedge(he)).idx();
            cmapper.indexOnPatch( msh.to_vertex(he).idx(), pi2, ci2 );
            os2 = offset[pi2];
            //+++ from_vertex
            M.insert(os +cpIndex(ci , 2, 0),c1dim  ) = (real_t)0.5;
            M.insert(os +cpIndex(ci , 2, 1),c1dim  ) = (real_t)1.0;
            M.insert(os2+cpIndex(ci2, 3, 0),c1dim++) = (real_t)0.5;
            //+++ symmetric to_vertex
            M.insert(os +cpIndex(ci , 3, 0),c1dim  ) = (real_t)0.5;
            M.insert(os +cpIndex(ci , 3, 1),c1dim  ) = (real_t)1.0;
            M.insert(os2+cpIndex(ci2, 2, 0),c1dim++) = (real_t)0.5;
            */
        }
        if (bdr) // boundary (EHE or regular)
        {
            bbmap.insert(he,c1dim++, 2,0, 1.0);
            bbmap.insert(he,c1dim++, 3,0, 1.0);
            bbmap.insert(he,c1dim++, 2,1, 1.0);
            bbmap.insert(he,c1dim++, 3,1, 1.0);
            
            // M.insert(os +cpIndex(ci , 2, 0),c1dim++) = (real_t)1.0;
            // M.insert(os +cpIndex(ci , 3, 0),c1dim++) = (real_t)1.0;
            // M.insert(os +cpIndex(ci , 2, 1),c1dim++) = (real_t)1.0;
            // M.insert(os +cpIndex(ci , 3, 1),c1dim++) = (real_t)1.0;
        }
    } //end for halfedges

    for ( auto ff : msh.faces() ) //for all faces
    {
        auto he = msh.halfedge(ff);
        bbmap.insert(he,c1dim++, 2,2, 1.0);
        bbmap.insert(he,c1dim++, 2,3, 1.0);
        bbmap.insert(he,c1dim++, 3,2, 1.0);
        bbmap.insert(he,c1dim++, 3,3, 1.0);
        
        // os = offset[ff.idx()];
        // M.insert(os+cpIndex(0, 2, 2),c1dim++) = (real_t)1.0;
        // M.insert(os+cpIndex(0, 2, 3),c1dim++) = (real_t)1.0;
        // M.insert(os+cpIndex(0, 3, 2),c1dim++) = (real_t)1.0;
        // M.insert(os+cpIndex(0, 3, 3),c1dim++) = (real_t)1.0;
    }

    //gsInfo << "Matrix:\n"<< M.toDense() <<"\n";
    gsInfo << "Number of functions is "<< c1dim <<"\n";
    GISMO_ASSERT(c1dim==M.cols(), "Not equal "<< c1dim <<"!="<< M.cols() );
    return M;
}

int main(int argc, char *argv[])
{
    std::string pname("gsview"), fn("/home/amantzaf/gitlab/catmull-clark/Gismo/Basis/triangle_planar.xml");
    index_t numSamples(1000), quantize(0);
    bool plot_mesh = false;
    bool plot_net = false;
    bool plot_patchid = false;
    bool get_basis = false;
    bool get_mesh = false;
    bool get_geo = false;
    bool plot = false;

    //! [Parse Command line]
    gsCmdLine cmd("Hi, give me a file (eg: .xml) and I will try to draw it!");

    cmd.addSwitch("geometry", "Try to find and plot a geometry contained in the file", get_geo);
    cmd.addSwitch("mesh"    , "Try to find and plot a mesh contained in the file", get_mesh);
    cmd.addSwitch("basis"   , "Try to find and plot a basis contained in the file", get_basis);
    cmd.addInt   ("s", "samples", "Number of samples to use for viewing", numSamples);
    cmd.addInt   ("q", "quantize", "Quantize the patches", quantize);
    cmd.addSwitch("element"   , "Plot the element mesh (when applicable)", plot_mesh);
    cmd.addSwitch("controlNet", "Plot the control net (when applicable)", plot_net);
    cmd.addSwitch("pid"  , "Plot the ID of each patch and boudanries as color", plot_patchid);
    cmd.addPlainString("filename", "File containing data to draw (.xml or third-party)", fn);
    cmd.addString("o", "oname", "Filename to use for the ParaView output", pname);
    cmd.addSwitch("plot", "Plot basis", plot);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse Command line]

    gsFileData<>  filedata(fn);
    gsMultiPatch<> mp;
    filedata.getFirst(mp);
    gsInfo<< "Got "<< mp <<"\n";

    for (index_t q= 0; q<quantize; ++q)
        mp = mp.uniformSplit();

    mp.computeTopology(1e-3);

    gsDofMapper mapper = mp.getMapper(1e-3); // np*36

    gsMultiBasis<> mb(mp);
    gsMultiPatch<> mp5(mp);
    if ( mb.degree(0)<5)
        mp5.degreeIncrease(5-mb.degree(0));

    mp.degreeReduce( mb.degree(0)-1 );
    gsDofMapper cmapper = mp.getMapper(1e-3); // np*4
    gsWrite(mp,"mp");
    
    //triangle_planar.xml: 7 vertices, 18 halfedges, 3 faces
    gsSurfMesh msh = mp.toMesh();
    msh.write("mp_mesh.off");
    gsWrite( mesh_to_multipatch(msh), "mp_mesh");

    /// TO DO: add indexOnPatch as face property in msh !!!!!!!!!!!!
    
//    msh.read("/user/amantzaf/home/Workspace/gitlab/catmull-clark-Marsala/Julia/Data/OFF_Meshes/triangle_planar.off");

    msh.write("out.off");
    gsInfo << "\nVertices: "<< msh.n_vertices();
    gsInfo << "\nEdges: "<< msh.n_edges();
    gsInfo << "\nFaces: "<< msh.n_faces();
    gsInfo << "\nHalfedges: "<< msh.n_halfedges() <<"\n";

    //Create transfer matrix
    gsSparseMatrix<> M = basisFromMesh(msh, cmapper);
    //gsMappedBasis<2, real_t> mbasis(mb, M);
    //gsFileData<> out;
    //out << mbasis;
    //out.save("basis");
    
    gsMappedSpline<2, real_t> srf(mp5, M); //Create spline using projection
    gsMultiPatch<> pr = srf.exportToPatches();
    gsFileData<> out;
    out << pr;
    out.save("surface_g1");

    if (!checkG1(pr, false))
            gsInfo <<"----------- Surface is not G1---\n";

//    gsWriteParaview(mp, ""e, numSamples, plot_mesh, plot_net);
//    gsFileManager::open(pname+".pvd");

    if (plot)
        PlotBasis(M, mp5, numSamples, false, true); //last false is for the net
    //gsFileManager::open("ccbasis.pvd");

    
    
    return EXIT_SUCCESS;
}
