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


// control point index starting from halfedge (corner0,corner1) and offset (xstep,ystep)
inline index_t cpIndex(index_t HEcorner, index_t xstep, index_t ystep)
{
    static const index_t str[]  = { 1, 6 }; // degree 5
    //HEcorner a is corner id (0,1,2,3) --> (0,5,30,35)
    // Note boundaty::corner equivalent id would be: 1,2,4,3
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
        GISMO_ERROR("indexOnPatch error.");
        break;
    }
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

gsSparseMatrix<> basisFromMesh(const gsSurfMesh & msh, const gsDofMapper & cmapper)
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
    for (auto v : msh.vertices()) //for all vertices
    {
        n = msh.valence(v);
        if (n==2) //corner vertex - we assume msh.is_boundary(v) [manifold mesh]
        {
            //note: msh.halfedge(v) returns a boundary half-edge
            auto h = msh.opposite_halfedge(msh.halfedge(v));
            pi = msh.face(h).idx();                  // patch
            os = offset[pi];
            cmapper.indexOnPatch(v.idx(), pi, ci); // corner id on patch
            // four "corner functions" to add (square)
            M.insert(os+cpIndex(ci, 0, 0),c1dim++) = (real_t)1.0;
            M.insert(os+cpIndex(ci, 0, 1),c1dim++) = (real_t)1.0;
            M.insert(os+cpIndex(ci, 1, 0),c1dim++) = (real_t)1.0;
            M.insert(os+cpIndex(ci, 1, 1),c1dim++) = (real_t)1.0;
        }

        if (n==3 && msh.is_boundary(v)) // boundary regular (T-junction)
        {
            auto h = msh.opposite_halfedge(msh.prev_halfedge(msh.halfedge(v)));
            pi = msh.face(h).idx();
            h = msh.ccw_rotated_halfedge(h);
            pi2 = msh.face(h).idx();
            //gsInfo <<"patches: "<< pi <<", "<< pi2 <<"\n";
            // first
            os = offset[pi];
            cmapper.indexOnPatch(v.idx(), pi, ci);
            M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)0.5;
            M.insert(os+cpIndex(ci, 1, 0),c1dim) = (real_t)1.0;
            os = offset[pi2];// next patch
            cmapper.indexOnPatch(v.idx(), pi2, ci);
            M.insert(os+cpIndex(ci, 0, 0), c1dim++) = (real_t)0.5;
            // second
            os = offset[pi];
            cmapper.indexOnPatch(v.idx(), pi, ci);
            M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)0.5;
            os = offset[pi2];// next patch
            cmapper.indexOnPatch(v.idx(), pi2, ci);
            M.insert(os+cpIndex(ci, 0, 1),c1dim) = (real_t)1.0;
            M.insert(os+cpIndex(ci, 0, 0), c1dim++) = (real_t)0.5;
            // third
            os = offset[pi];
            cmapper.indexOnPatch(v.idx(), pi, ci);
            M.insert(os+cpIndex(ci, 0, 1), c1dim) = (real_t)0.5;
            M.insert(os+cpIndex(ci, 1, 1), c1dim) = (real_t)1.0;
            os = offset[pi2];// next patch
            cmapper.indexOnPatch(v.idx(), pi2, ci);
            M.insert(os+cpIndex(ci, 1, 0), c1dim++) = (real_t)0.5;
            // fourth
            os = offset[pi];
            cmapper.indexOnPatch(v.idx(), pi, ci);
            M.insert(os+cpIndex(ci, 0, 1), c1dim) = (real_t)0.5;
            os = offset[pi2];// next patch
            cmapper.indexOnPatch(v.idx(), pi2, ci);
            M.insert(os+cpIndex(ci, 1, 0), c1dim) = (real_t)0.5;
            M.insert(os+cpIndex(ci, 1, 1), c1dim++) = (real_t)1.0;

            //MODIF BOUNDARY REGULAR neighbors
        }

        if (n==4 && !msh.is_boundary(v) ) // inner regular
        {
            for ( auto he : msh.halfedges(v) ) // iterate over all HE (n=4 functions)
            {
                auto h = he; // First patch
                pi = msh.face(h).idx();
                os = offset[pi];
                cmapper.indexOnPatch(v.idx(), pi, ci);
                M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)0.25;
                M.insert(os+cpIndex(ci, 1, 0), c1dim) = (real_t)0.5;
                M.insert(os+cpIndex(ci, 0, 1), c1dim) = (real_t)0.5;
                M.insert(os+cpIndex(ci, 1, 1), c1dim) = (real_t)1.0;//ANCHOR!

                v_N = msh.to_vertex(h);
                n_N = msh.valence(v_N); //modif first axis
                if ( 4 != n_N && !msh.is_boundary(v_N) )
                {
                    a_N = 2 * math::cos(2*EIGEN_PI/n_N);
                    M.insert(os+cpIndex(ci, 2, 0), c1dim) = (real_t)0.225;//=9/40
                    M.insert(os+cpIndex(ci, 2, 1), c1dim) = (real_t)(0.225 - a_N/80);
                    M.insert(os+cpIndex(ci, 3, 1), c1dim) = (real_t)(0.0675*a_N);//=27/400*a_N
                }
                v_N = msh.to_vertex( msh.ccw_rotated_halfedge(h) );
                n_N = msh.valence(v_N); //modif second axis
                if ( 4 != n_N && !msh.is_boundary(v_N) )
                {
                    a_N = 2 * math::cos(2*EIGEN_PI/n_N);
                    M.insert(os+cpIndex(ci, 0, 2), c1dim) = (real_t)0.225;//=9/40
                    M.insert(os+cpIndex(ci, 1, 2), c1dim) = (real_t)(0.225 - a_N/80);
                    M.insert(os+cpIndex(ci, 1, 3), c1dim) = (real_t)(0.0675*a_N);//=27/400*a_N
                }
                
                h = msh.ccw_rotated_halfedge(h);// Second patch
                pi = msh.face(h).idx();
                os = offset[pi];
                cmapper.indexOnPatch(v.idx(), pi, ci);
                M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)0.25;
                M.insert(os+cpIndex(ci, 1, 0), c1dim) = (real_t)0.5;

                v_N = msh.to_vertex(h);
                n_N = msh.valence(v_N); //modif first axis
                if ( 4 != n_N && !msh.is_boundary(v_N) )
                {
                    a_N = 2 * math::cos(2*EIGEN_PI/n_N);
                    M.insert(os+cpIndex(ci, 2, 0), c1dim) = (real_t)0.225;//=9/40
                    M.insert(os+cpIndex(ci, 2, 1), c1dim) = (real_t)(0.225 - a_N/80);
                    M.insert(os+cpIndex(ci, 3, 1), c1dim) = (real_t)(0.0675*a_N);//=27/400*a_N
                }
                v_N = msh.to_vertex( msh.ccw_rotated_halfedge(h) );
                n_N = msh.valence(v_N); //modif second axis
                if ( 4 != n_N && !msh.is_boundary(v_N) )
                {
                    a_N = 2 * math::cos(2*EIGEN_PI/n_N);
                    M.insert(os+cpIndex(ci, 0, 2), c1dim) = (real_t)(-0.025);//=-1/40
                    M.insert(os+cpIndex(ci, 1, 2), c1dim) = -(2-a_N)/80;
                    M.insert(os+cpIndex(ci, 1, 3), c1dim) = -0.0075*a_N;//=-3/400*a_N
                }

                h = msh.ccw_rotated_halfedge(h);// Third patch
                pi = msh.face(h).idx();
                os = offset[pi];
                cmapper.indexOnPatch(v.idx(), pi, ci);
                M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)0.25;

                v_N = msh.to_vertex(h);
                n_N = msh.valence(v_N); //modif first axis
                if ( 4 != n_N && !msh.is_boundary(v_N) )
                {
                    a_N = 2 * math::cos(2*EIGEN_PI/n_N);
                    M.insert(os+cpIndex(ci, 2, 0), c1dim) = (real_t)(-0.025);//=-1/40
                    M.insert(os+cpIndex(ci, 2, 1), c1dim) = -(2-a_N)/80;
                    M.insert(os+cpIndex(ci, 3, 1), c1dim) = -0.0075*a_N;//=-3/400*a_N
                }
                v_N = msh.to_vertex( msh.ccw_rotated_halfedge(h) );
                n_N = msh.valence(v_N); //modif second axis
                if ( 4 != n_N && !msh.is_boundary(v_N) )
                {
                    a_N = 2 * math::cos(2*EIGEN_PI/n_N);
                    M.insert(os+cpIndex(ci, 0, 2), c1dim) = (real_t)(-0.025);//=-1/40
                    M.insert(os+cpIndex(ci, 1, 2), c1dim) = -(2-a_N)/80; //------------------
                    M.insert(os+cpIndex(ci, 1, 3), c1dim) = -0.0075*a_N;//=-3/400*a_N
                }
                
                h = msh.ccw_rotated_halfedge(h); // Forth patch
                pi = msh.face(h).idx();
                os = offset[pi];
                cmapper.indexOnPatch(v.idx(), pi, ci);
                M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)0.25;
                M.insert(os+cpIndex(ci, 0, 1), c1dim) = (real_t)0.5;

                n_N = msh.valence(v_N); //modif first axis
                if ( 4 != n_N && !msh.is_boundary(v_N) )
                {
                    a_N = 2 * math::cos(2*EIGEN_PI/n_N);
                    M.insert(os+cpIndex(ci, 2, 0), c1dim) = (real_t)(-0.025);//=-1/40
                    M.insert(os+cpIndex(ci, 2, 1), c1dim) = -(2-a_N)/80;
                    M.insert(os+cpIndex(ci, 3, 1), c1dim) = -0.0075*a_N;//=-3/400*a_N
                }
                v_N = msh.to_vertex( msh.ccw_rotated_halfedge(h) );
                n_N = msh.valence(v_N); //modif second axis
                if ( 4 != n_N && !msh.is_boundary(v_N) )
                {
                    a_N = 2 * math::cos(2*EIGEN_PI/n_N);
                    M.insert(os+cpIndex(ci, 0, 2), c1dim) = (real_t)0.225;//=9/40
                    M.insert(os+cpIndex(ci, 1, 2), c1dim) = (real_t)(0.225 - a_N/80);
                    M.insert(os+cpIndex(ci, 1, 3), c1dim) = (real_t)(0.0675*a_N);//=27/400*a_N
                }
                ++c1dim;
            }
        }

          if (n!=4 && !msh.is_boundary(v)) //inner EV (n+3 functions)
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
              for ( auto ff : msh.faces(v) ) //for all faces, set vertex value
              {
                  pi = ff.idx();
                  cmapper.indexOnPatch( v.idx(), pi, ci);
                  os = offset[pi];
                  M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)1.0;
              }
              index_t k = 0;// (continues) Basis function attached to the vertex value
              for ( auto he : msh.halfedges(v) ) //for all halfedges he in msh
              {
                  real_t R1 = 0.1 * ( -a + 5*a*sol1.at(k) - 10*(a-2)*(0.5*sol1.at(k)-0.1));
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
                  M.insert(os+cpIndex(ci, 1, 0), 1+c1dim) = K1(k,1);
                  M.insert(os+cpIndex(ci, 2, 0), 1+c1dim) = K1(k,1) * 0.5;
                  M.insert(os+cpIndex(ci, 2, 1), 1+c1dim) = 0.5 * R2;//..
                  M.insert(os+cpIndex(ci, 3, 1), 1+c1dim) = 0.5 * R4;//..
                  M.insert(os+cpIndex(ci, 1, 1), 1+c1dim) = sol2(k,1);
                  pi = msh.face(msh.opposite_halfedge(he)).idx();
                  cmapper.indexOnPatch(v.idx(), pi, ci);
                  os = offset[pi];
                  //first
                  M.insert(os+cpIndex(ci, 0, 1), c1dim) = K1(k,0);
                  M.insert(os+cpIndex(ci, 0, 2), c1dim) = K1(k,0) * 0.5;
                  M.insert(os+cpIndex(ci, 1, 2), c1dim) = 0.5 * R1;
                  M.insert(os+cpIndex(ci, 1, 3), c1dim) = 0.5 * R3;//..
                  //second
                  M.insert(os+cpIndex(ci, 0, 1), 1+c1dim) = K1(k,1);
                  M.insert(os+cpIndex(ci, 0, 2), 1+c1dim) = K1(k,1) * 0.5;
                  M.insert(os+cpIndex(ci, 1, 2), 1+c1dim) = 0.5 * R2;//..
                  M.insert(os+cpIndex(ci, 1, 3), 1+c1dim) = 0.5 * R4;//..
                  ++k;
              }
              c1dim+=2;

              // Basis attached to the cross derivatives at the vertex
              k = 0;
              for ( auto he : msh.halfedges(v) ) //for all halfedges he in msh
              {
                  real_t R1 = -(a-2)*(5.0/(4*a)) + (3.0/5.0)*a*(5.0/(4*a));
                  real_t R2 =      a*(5.0/(4*a)) - (a-2)*(5.0/(4*a));

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
        
    for ( auto he : msh.halfedges() ) //for all halfedges he in msh
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
            M.insert(os +cpIndex(ci , 2, 1),c1dim) = (real_t)1.0;
            M.insert(os2+cpIndex(ci2, 3, 1),c1dim) = (real_t)(-1.0);
            ++c1dim;
        }
        if (!bdr && !ee) // inner regular HE
        {
            pi2 = msh.face(msh.opposite_halfedge(he)).idx();
            cmapper.indexOnPatch( msh.to_vertex(he).idx(), pi2, ci2 );
            os2 = offset[pi2];

            M.insert(os +cpIndex(ci , 2, 0),c1dim  ) = (real_t)0.5;
            M.insert(os +cpIndex(ci , 2, 1),c1dim  ) = (real_t)1.0;
            M.insert(os2+cpIndex(ci2, 3, 0),c1dim++) = (real_t)0.5;

            M.insert(os +cpIndex(ci , 3, 0),c1dim  ) = (real_t)0.5;
            M.insert(os +cpIndex(ci , 3, 1),c1dim  ) = (real_t)1.0;
            M.insert(os2+cpIndex(ci2, 2, 0),c1dim++) = (real_t)0.5;
        }
        if (bdr) // boundary (EHE or regular)
        {
            M.insert(os +cpIndex(ci , 2, 0),c1dim++) = (real_t)1.0;
            M.insert(os +cpIndex(ci , 3, 0),c1dim++) = (real_t)1.0;
            M.insert(os +cpIndex(ci , 2, 1),c1dim++) = (real_t)1.0;
            M.insert(os +cpIndex(ci , 3, 1),c1dim++) = (real_t)1.0;
        }
    } //end for halfedges

    for ( auto ff : msh.faces() ) //for all faces
    {
        os = offset[ff.idx()];
        M.insert(os+cpIndex(0, 2, 2),c1dim++) = (real_t)1.0;
        M.insert(os+cpIndex(0, 2, 3),c1dim++) = (real_t)1.0;
        M.insert(os+cpIndex(0, 3, 2),c1dim++) = (real_t)1.0;
        M.insert(os+cpIndex(0, 3, 3),c1dim++) = (real_t)1.0;
    }

    //gsInfo << "Matrix:\n"<< M.toDense() <<"\n";
    gsInfo << "Number of functions is "<< c1dim <<"\n";
    GISMO_ASSERT(c1dim==M.cols(), "Not equal "<< c1dim <<"!="<< M.cols() );
    return M;
}

int main(int argc, char *argv[])
{
    std::string pname("gsview"), fn("/home/amantzaf/gitlab/catmull-clark/Gismo/Basis/triangle_planar.xml");
    index_t numSamples(1000);
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
