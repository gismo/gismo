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

void PlotBasis(const gsSparseMatrix<> & cf, const gsMultiPatch<> & mp,
               const int & ns,
               const bool pm = false, const bool pn = false)
{
    //gsFileData<> xml;
    gsParaviewCollection col("ccbasis");
    gsMultiPatch<> mmp = mp;
    const index_t d = mp.parDim();
    gsInfo << "dim " <<d<<"\n";
    const size_t np = mp.nPatches();
    gsVector<index_t> offset(np+1);
    offset[0] = 0;
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

index_t g1Dimension(gsSurfMesh & msh)
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

    for ( auto ee : msh.edges() ) //for all edges he in msh
    {
        bool bdr = msh.is_boundary(ee);
        // Note for Extraordinary edge on boundary we still count 4)
        bool eed = (4!=msh.valence(msh.vertex(ee,0)) || 4!=msh.valence(msh.vertex(ee,1)) );
        c1dim += (!bdr && eed ? 2 : 4);
    }
    gsInfo << "Dimension V+E: "<< c1dim <<"\n";
    
    c1dim += 4 * msh.n_faces();
    gsInfo << "Dimension is (final): "<< c1dim <<"\n";
    return c1dim;
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
    mp.degreeReduce(4);
    gsDofMapper cmapper = mp.getMapper(1e-3); // np*4
    gsWrite(mp,"mp");

    //triangle_planar.xml: 7 vertices, 18 halfedges, 3 faces
    gsSurfMesh msh = mp.toMesh();
//    msh.read("/user/amantzaf/home/Workspace/gitlab/catmull-clark-Marsala/Julia/Data/OFF_Meshes/triangle_planar.off");
    auto pt = msh.get_vertex_property<gismo::Point>("v:point");
    msh.write("out.off");

    gsInfo << "\nVertices: "<< msh.n_vertices();
    gsInfo << "\nEdges: "<< msh.n_edges();
    gsInfo << "\nFaces: "<< msh.n_faces();
    gsInfo << "\nHalfedges: "<< msh.n_halfedges() <<"\n";

    int n, c1dim = g1Dimension(msh);

    //Create transfer matrix
    gsSparseMatrix<> M(mapper.mapSize(), c1dim);
    index_t pi, pi2, os, os2, ci, ci2;

    /*
    cmapper.print();
    for (auto v : msh.vertices()) //for all vertices
    {
        gsInfo <<"On vertex "<<v.idx() <<"\n";
        for ( auto he : msh.halfedges(v) ) // iterate over all HE (n=4 functions)
        {
            if ( msh.is_boundary(he) ) continue;
            pi = msh.face(he).idx(); // patch
            os = mapper.offset(pi);
            //cmapper.indexOnPatch( msh.from_vertex(he).idx(), pi, ci); // corner id on patch
            cmapper.indexOnPatch( v.idx(), pi, ci); // corner id on patch

            gsInfo <<"  [ci="<<ci<<"], CP("<<pi<<"): "<< cpIndex(ci, 0, 0)<< ", u-:"<< cpIndex(ci, 1, 0)<<", v-:"<<cpIndex(ci, 0, 1)<<", uv-: "<<cpIndex(ci, 1, 1) <<" --> "<<cpIndex(ci, 2, 1)<<"\n";
            gsInfo<<"point = "<< mp.patch(pi).coef(ci)<< " ~~ "<< pt[v].transpose() <<"\n";
        }
    }
    //*/     

    c1dim = 0;//reset

    //auto pid = msh.get_vertex_property<index_t>("v:patch");
    //auto anchor = msh.get_vertex_property<index_t>("v:anchor");
    //auto dof = msh.get_vertex_property<index_t>("v:dof");

    for (auto v : msh.vertices()) //for all vertices
    {                    
        n = msh.valence(v);
        if (n==2) //corner vertex - we assume msh.is_boundary(v) [manifold mesh]
        {
            //note: msh.halfedge(v) returns a boundary half-edge
            auto h = msh.opposite_halfedge(msh.halfedge(v));
            pi = msh.face(h).idx();                  // patch
            os = mapper.offset(pi);
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
            os = mapper.offset(pi);
            cmapper.indexOnPatch(v.idx(), pi, ci);
            M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)0.5;
            M.insert(os+cpIndex(ci, 1, 0),c1dim) = (real_t)1.0;
            os = mapper.offset(pi2);// next patch
            cmapper.indexOnPatch(v.idx(), pi2, ci);
            M.insert(os+cpIndex(ci, 0, 0), c1dim++) = (real_t)0.5;
            // second
            os = mapper.offset(pi);
            cmapper.indexOnPatch(v.idx(), pi, ci);
            M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)0.5;
            os = mapper.offset(pi2);// next patch
            cmapper.indexOnPatch(v.idx(), pi2, ci);
            M.insert(os+cpIndex(ci, 0, 1),c1dim) = (real_t)1.0;
            M.insert(os+cpIndex(ci, 0, 0), c1dim++) = (real_t)0.5;
            // third
            os = mapper.offset(pi);
            cmapper.indexOnPatch(v.idx(), pi, ci);
            M.insert(os+cpIndex(ci, 0, 1), c1dim) = (real_t)0.5;
            M.insert(os+cpIndex(ci, 1, 1),c1dim) = (real_t)1.0;
            os = mapper.offset(pi2);// next patch
            cmapper.indexOnPatch(v.idx(), pi2, ci);
            M.insert(os+cpIndex(ci, 1, 0), c1dim++) = (real_t)0.5;
            // fourth
            os = mapper.offset(pi);
            cmapper.indexOnPatch(v.idx(), pi, ci);
            M.insert(os+cpIndex(ci, 0, 1), c1dim) = (real_t)0.5;
            os = mapper.offset(pi2);// next patch
            cmapper.indexOnPatch(v.idx(), pi2, ci);
            M.insert(os+cpIndex(ci, 1, 0),c1dim) = (real_t)0.5;
            M.insert(os+cpIndex(ci, 1, 1), c1dim++) = (real_t)1.0;
        }

        if (n==4 && !msh.is_boundary(v) ) // inner regular
        {
            for ( auto he : msh.halfedges(v) ) // iterate over all HE (n=4 functions)
            {
                auto h = he;
                pi = msh.face(h).idx();
                os = mapper.offset(pi);
                cmapper.indexOnPatch(v.idx(), pi, ci);
                M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)0.25;
                M.insert(os+cpIndex(ci, 1, 0), c1dim) = (real_t)0.5;
                M.insert(os+cpIndex(ci, 0, 1), c1dim) = (real_t)0.5;
                M.insert(os+cpIndex(ci, 1, 1), c1dim) = (real_t)1.0;    
                h = msh.ccw_rotated_halfedge(h);// next patch
                pi = msh.face(h).idx();
                os = mapper.offset(pi);
                cmapper.indexOnPatch(v.idx(), pi, ci);
                M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)0.25;
                M.insert(os+cpIndex(ci, 1, 0), c1dim) = (real_t)0.5;
                h = msh.ccw_rotated_halfedge(h);// next patch
                pi = msh.face(h).idx();
                os = mapper.offset(pi);
                cmapper.indexOnPatch(v.idx(), pi, ci);
                M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)0.25;
                h = msh.ccw_rotated_halfedge(h); // next patch
                pi = msh.face(h).idx();
                os = mapper.offset(pi);
                cmapper.indexOnPatch(v.idx(), pi, ci);
                M.insert(os+cpIndex(ci, 0, 1), c1dim++) = (real_t)0.5;
            }
        }

          if (n!=4 && !msh.is_boundary(v)) //inner EV (n+3 functions)
          {
              gsInfo << "EV: " << c1dim <<"\n";
//              /*
              // Computing n + 3 functions attached to the EV (1+2+n)
              real_t a = 2 * math::cos(2*EIGEN_PI/n);// rename as a0

              // Basis function attached to the vertex value
              for ( auto ff : msh.faces(v) ) //for all faces, set vertex value
              {
                  pi = ff.idx();
                  cmapper.indexOnPatch( v.idx(), pi, ci);
                  os = mapper.offset(pi);
                  M.insert(os+cpIndex(ci, 0, 0), c1dim) = (real_t)1.0;
              }

              gsMatrix<> C = circulant(n);
              C.transposeInPlace(); // Making column-circulant !!
              gsInfo << "C:\n " << C <<"\n";
              gsMatrix<> A(n+2,n), b(n+2,1);
              b.topRows(n).setConstant(2 - a);
              b.bottomRows(2).setZero();

              index_t sz = n + (n%2==0);
              gsMatrix<> C2(sz,n);
              C2.topRows(n) = C;
              C2 /*.topRows(n)*/ .diagonal().array() += 1;
              gsMatrix<> Cpow = C;
              for(index_t i = 2; i<n;++i)
                  Cpow *= C;
              auto C1 = A.topRows(n); //ref
              C1.noalias() = C + Cpow;
              C1.diagonal().array() -= a;
              gsInfo << "C1:\n " << C1 <<"\n";

              typename gsMatrix<>::JacobiSVD svd = C1.jacobiSvd(gsEigen::ComputeFullV /*| gsEigen::ComputeFullU*/ );
              //gsInfo << "U:\n " << svd.matrixU() <<"\n";
              gsInfo << "V:\n " << svd.matrixV() <<"\n";
              gsInfo << "S:\n " << svd.singularValues().transpose() <<"\n";
              gsMatrix<> K = svd.matrixV().rightCols(2);
              // Get nullspace from Julia for n=3:
              //gsMatrix<> K(3,2);
              //K<<-0.499814, -0.645641,
              //   -0.309234,  0.755672,
              //    0.809049, -0.110032;

              //
              //gsInfo << "svalues: " << svd.singularValues().transpose() <<"\n";
              gsInfo << "K:\n " << K <<"\n";
              gsInfo << "K0-->:\n " << K.col(0).norm() <<"\n";
              gsInfo << "K1-->:\n " << K.col(1).norm() <<"\n";
              
              A.bottomRows(2) = K.transpose();
              gsInfo << "A:\n " << A <<"\n";

              gsMatrix<> sol_odd, sol = A.fullPivLu().solve(b);
              gsInfo << "sol: " << sol.transpose() <<"\n";

              b.resize(sz, 1);// reusing memory for b
              // b.array() = -a;
              // b +=-5*(a-2)*sol + 4 * a * (0.5*sol - 0.1*gsMatrix<>::Ones(n,1)) ;
              // b *= 0.2;
              for(index_t q=0;q!=n;++q)
                  b.at(q) = 0.2 * ( a-5*(a-2)*sol.at(q) + 4*a*(0.5*sol.at(q)-0.1));

              if (sz==n) //odd vertex
              {
                  sol_odd = C2.fullPivLu().solve(b);
                  gsInfo << "sol_odd: " << sol_odd.transpose() <<"\n";
              }
              else //even vertex
              {
                  //add to C2 one last row: (1,N): [-1 1 -1 1 -1 1 ...]
                  // C2 has now size: (N+1)xN
                  auto RR = C2.row(n);
                  index_t tmp(1);
                  for(index_t q = 0; q!=n; ++q)
                  {
                      tmp *= -1;
                      RR[q]= tmp;
                  }
                  b.at(n) = 0;
                  // solve using C2 plus one additional row
                  // NOTE: solve OVERDETERMINED system (LS?)
                  sol_odd = C2.fullPivLu().solve(b);
                  gsInfo << "sol_even: " << sol_odd.transpose() <<"\n";
              }
              gsInfo << "C2:\n " << C2 <<"\n";
              //gsInfo << "b:\n " << b.transpose() <<"\n";

              // (continues) Basis function attached to the vertex value
              index_t k = 0;
              for ( auto he : msh.halfedges(v) ) //for all halfedges he in msh
              {
                  real_t R1 = 0.1 * ( -a + 5*a*sol.at(k) - 10*(a-2)*(0.5*sol.at(k)-0.1));
                  pi = msh.face(he).idx();
                  os = mapper.offset(pi);
                  cmapper.indexOnPatch(v.idx(), pi, ci);
                  M.insert(os+cpIndex(ci, 1, 1), c1dim) = sol_odd.at(k);
                  //
                  M.insert(os+cpIndex(ci, 1, 0), c1dim) = sol.at(k);
                  M.insert(os+cpIndex(ci, 2, 0), c1dim) = 0.5*sol.at(k) - 0.1;
                  M.insert(os+cpIndex(ci, 2, 1), c1dim) = 0.5*R1;
                  pi = msh.face(msh.opposite_halfedge(he)).idx();
                  cmapper.indexOnPatch(v.idx(), pi, ci);
                  os = mapper.offset(pi);
                  M.insert(os+cpIndex(ci, 0, 1), c1dim) = sol.at(k);
                  M.insert(os+cpIndex(ci, 0, 2), c1dim) = 0.5*sol.at(k) - 0.1;
                  M.insert(os+cpIndex(ci, 1, 2), c1dim) = 0.5*R1;
                  ++k;
              }
              //gsInfo <<"Function "<< c1dim <<" :\n" << M.col(c1dim).toDense().transpose() <<"\n";
              ++c1dim;

              //Start: basis functions attached to the partial derivatives at the vertex
              b.setZero(sz,2); // reusing memory for b
              //b += -5*(a-2)*K + 4 * a * (K * 0.5) ;
              //b *= 0.2;
              for(index_t q=0;q!=n;++q)
              {
                  b(q,0) = 0.2 * (-5*(a-2)*K(q,0) + 4*a*(0.5*K(q,0)));
                  b(q,1) = 0.2 * (-5*(a-2)*K(q,1) + 4*a*(0.5*K(q,1)));
              }

              if (sz==n) // odd vertex
              {
                  sol_odd = C2.fullPivLu().solve(b);
                  gsInfo << "b_odd:\n" << b.transpose() <<"\n";
                  gsInfo << "sol_odd:\n" << sol_odd.transpose() <<"\n";
              }
              else //even vertex
              {
                  //add to C2 one last row: (1,N): [-1 1 -1 1 -1 1 ...]
                  // C2 has now size: (N+1)xN
                  auto RR = C2.row(n);
                  index_t tmp(1);
                  for(index_t q = 0; q!=n; ++q)
                  {
                      tmp *= -1;
                      RR[q]= tmp;
                  }
                  b.at(n) = 0;
                  // solve using C2 plus one additional row
                  // NOTE: solve OVERDETERMINED system (LS?)
                  sol_odd = C2.fullPivLu().solve(b);
                  gsInfo << "sol_even: " << sol_odd.transpose() <<"\n";
              }

              //Basis functions attached to the two partial derivatives at the vertex.
              k = 0;
              for ( auto he : msh.halfedges(v) ) //for all halfedges he in msh
              {
                  real_t R1 = 0.1 * ( 5*a*K(k,0) - 10*(a-2)*(0.5*K(k,0)));
                  real_t R3 = 0.1 * (-5*a*K(k,0) + 10*a*(0.5*K(k,0)));
                  //
                  real_t R2 = 0.1 * ( 5*a*K(k,1) - 10*(a-2)*(0.5*K(k,1)));
                  real_t R4 = 0.1 * (-5*a*K(k,1) + 10*a*(0.5*K(k,1)));

                  pi = msh.face(he).idx();
                  os = mapper.offset(pi);
                  cmapper.indexOnPatch(v.idx(), pi, ci);
                  //first
                  M.insert(os+cpIndex(ci, 1, 0), c1dim) = K(k,0);
                  M.insert(os+cpIndex(ci, 2, 0), c1dim) = K(k,0) * 0.5;
                  M.insert(os+cpIndex(ci, 2, 1), c1dim) = 0.5 * R1;
                  M.insert(os+cpIndex(ci, 3, 1), c1dim) = 0.5 * R3;//
                  M.insert(os+cpIndex(ci, 1, 1), c1dim) = sol_odd(k,0);
                  //second
                  M.insert(os+cpIndex(ci, 1, 0), 1+c1dim) = K(k,1);
                  M.insert(os+cpIndex(ci, 2, 0), 1+c1dim) = K(k,1) * 0.5;
                  M.insert(os+cpIndex(ci, 2, 1), 1+c1dim) = 0.5 * R2;//..
                  M.insert(os+cpIndex(ci, 3, 1), 1+c1dim) = 0.5 * R4;//..
                  M.insert(os+cpIndex(ci, 1, 1), 1+c1dim) = sol_odd(k,1);
                  pi = msh.face(msh.opposite_halfedge(he)).idx();
                  cmapper.indexOnPatch(v.idx(), pi, ci);
                  os = mapper.offset(pi);
                  //first
                  M.insert(os+cpIndex(ci, 0, 1), c1dim) = K(k,0);
                  M.insert(os+cpIndex(ci, 0, 2), c1dim) = K(k,0) * 0.5;
                  M.insert(os+cpIndex(ci, 1, 2), c1dim) = 0.5 * R1;
                  M.insert(os+cpIndex(ci, 1, 3), c1dim) = 0.5 * R3;//..
                  //second
                  M.insert(os+cpIndex(ci, 0, 1), 1+c1dim) = K(k,1);
                  M.insert(os+cpIndex(ci, 0, 2), 1+c1dim) = K(k,1) * 0.5;
                  M.insert(os+cpIndex(ci, 1, 2), 1+c1dim) = 0.5 * R2;//..
                  M.insert(os+cpIndex(ci, 1, 3), 1+c1dim) = 0.5 * R4;//..
                  ++k;
              }
              gsInfo <<"Function "<< c1dim <<" :\n";
              gsInfo << M.col(c1dim).toDense().transpose() <<"\n";

              c1dim+=2;

              // Basis attached to the cross derivatives at the vertex
              k = 0;
              for ( auto he : msh.halfedges(v) ) //for all halfedges he in msh
              {
                  real_t R1 = -(a-2)*(5.0/(4*a)) + (3.0/5.0)*a*(5.0/(4*a));
                  real_t R2 =      a*(5.0/(4*a)) - (a-2)*(5.0/(4*a));

                  pi = msh.face(he).idx();
                  os = mapper.offset(pi);
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
                  os = mapper.offset(pi);
                  M.insert(os+cpIndex(ci, 0, 2), c1dim) =  5.0/(4*a);
                  M.insert(os+cpIndex(ci, 0, 3), c1dim) =  5.0/(4*a);

                  M.insert(os+cpIndex(ci, 1, 2), c1dim) =  0.5*R1;
                  M.insert(os+cpIndex(ci, 1, 3), c1dim) =  0.5*R2;
                  //
                  pi = msh.face(msh.ccw_rotated_halfedge(he)).idx();
                  cmapper.indexOnPatch(v.idx(), pi, ci); // next face F+1
                  os = mapper.offset(pi);
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

    for ( auto he : msh.halfedges() ) //for all halfedges he in msh
    {
        if ( msh.is_boundary(he) ) continue;
        // Edge on boundary ?
        bool bdr = msh.touches_boundary(he);
        // Extraordinary edge? (note: no new functions for boundart EHE)
        bool ee = (4!=msh.valence(msh.from_vertex(he)) || 4!=msh.valence(msh.from_vertex(he)) );

        pi = msh.face(he).idx(); // patch
        os = mapper.offset(pi);
        cmapper.indexOnPatch( msh.from_vertex(he).idx(), pi, ci); // corner id on patch

        if (!bdr && ee) // inner EHE
        {
            pi2 = msh.face(msh.opposite_halfedge(he)).idx();
            cmapper.indexOnPatch( msh.to_vertex(he).idx(), pi2, ci2 );
            os2 = mapper.offset(pi2);
            //two functions per edge (or one per halfedge)
            M.insert(os +cpIndex(ci , 2, 1),c1dim) = (real_t)1.0;
            M.insert(os2+cpIndex(ci2, 3, 1),c1dim) = (real_t)(-1.0);
            ++c1dim;
        }
        if (!bdr && !ee) // inner regular HE
        {
            //4 
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
        os = mapper.offset(ff.idx());
        M.insert(os+cpIndex(0, 2, 2),c1dim++) = (real_t)1.0;
        M.insert(os+cpIndex(0, 2, 3),c1dim++) = (real_t)1.0;
        M.insert(os+cpIndex(0, 3, 2),c1dim++) = (real_t)1.0;
        M.insert(os+cpIndex(0, 3, 3),c1dim++) = (real_t)1.0;
    }

    //gsInfo << "Matrix:\n"<< M.toDense() <<"\n";
    gsInfo << "Number of functions is "<< c1dim <<"\n";

    //gsMappedBasis<2, real_t> mbasis(mb, M);
    //gsFileData<> out;
    //out << mbasis;
    //out.save("basis");

//    gsWriteParaview(mp, ""e, numSamples, plot_mesh, plot_net);
//    gsFileManager::open(pname+".pvd");

    PlotBasis(M, mp5, 1000, false, true); //last false is for the net
    //gsFileManager::open("ccbasis.pvd");
                
    return EXIT_SUCCESS;
}
