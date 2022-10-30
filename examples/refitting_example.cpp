/** @file refit_patches.cpp

    @brief Computes patches from structured (tensor-product) data samples by fitting.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#include <gismo.h>
#include <gsModeling/gsPatchGenerator.h>

using namespace gismo;


/**
 * @brief      Parameterizes the points for a curve using centripetal (alpha=0)
 *             or chord-length (alpha=0.5) parameterization
 *
 * @param[in]  xyz    The coordinates of the points, each column is a point
 * @param[in]  alpha  Alpha coefficient
 *
 * @tparam     T      Double
 *
 * @return     Parametric coordinates for each point, per column
 */
template<class T>
gsMatrix<T> parameterize_points1D(const gsMatrix<T> & xyz, T alpha = 0.5)
{
    index_t N = xyz.cols();
    gsMatrix<T> uv(1,N);

    uv(0,0)   = 0;
    uv(0,N-1) = 1;
    std::vector<T> distances(N-1);
    for (index_t k=1; k!=N; k++)
        distances.at(k-1) = std::pow((xyz.col(k)-xyz.col(k-1)).norm(),alpha);

    T d = std::accumulate(distances.begin(),distances.end(),0.0);
    for (index_t k=1; k!=N-1; k++)
        uv(0,k) = uv(0,k-1) + distances.at(k) / d;

    return uv;
}


/**
 * @brief      Parameterizes the points for a surface using centripetal
 *             (alpha=0) or chord-length (alpha=0.5) parameterization
 *
 * @param[in]  bbasis  A Tensor B-spline basis to fit on
 * @param[in]  xyz     The coordinates of the points, each column is a point
 * @param[in]  alpha   Alpha coefficient
 *
 * @tparam     T       Double
 *
 * @return     Parametric coordinates for each point, per column
 */
template<class T>
gsMatrix<T> parameterize_points2D(const gsTensorBSplineBasis<2,real_t> & bbasis, const gsMatrix<T> & xyz, T alpha = 0.5)
{
    gsMatrix<T> uv(2,xyz.cols());
    gsMatrix<index_t> slice;
    gsMatrix<T> tmp_coefs;
    gsMatrix<T> u_mat(bbasis.size(0),bbasis.size(1));
    gsMatrix<T> v_mat(bbasis.size(0),bbasis.size(1));
    for(index_t k=0; k!=bbasis.size(0); k++)
    {
        slice = bbasis.coefSlice(0,k);
        tmp_coefs.resize(xyz.rows(),slice.rows());
        for (index_t c=0; c!=slice.rows(); c++)
            tmp_coefs.col(c) = xyz.col(slice(c,0));

        u_mat.row(k) = parameterize_points1D(tmp_coefs,alpha);
    }
    u_mat.array() /= u_mat.rows();
    u_mat = u_mat.colwise().sum(); // takes the sum of the rows for each column

    for(index_t k=0; k!=bbasis.size(1); k++)
    {
        slice = bbasis.coefSlice(1,k);
        tmp_coefs.resize(xyz.rows(),slice.rows());
        for (index_t c=0; c!=slice.rows(); c++)
            tmp_coefs.col(c) = xyz.col(slice(c,0));
        v_mat.col(k) = parameterize_points1D(tmp_coefs,alpha).transpose();
    }
    v_mat.array() /= v_mat.cols();
    v_mat = v_mat.rowwise().sum();

    for (index_t i=0; i!=u_mat.cols(); i++)
        for (index_t j=0; j!=v_mat.rows(); j++)
        {
            uv.col(bbasis.index(j,i))<<v_mat(j,0),u_mat(0,i);
        }
    return uv;
}

int main(int argc, char *argv[])
{
    std::string filename("domain2d/yeti_mp2.xml");
    bool plot = false;
    real_t tol = 1e-5;
    index_t nknots = 5, degree = 3;
    index_t npts = 100;
    real_t lambda_crv = 0;
    real_t lambda_srf = 0;
    // real_t gtol = 1e-6;
    // bool reparam = false, gaps = true;

    gsCmdLine cmd("Computes patches from structured (tensor-product) data samples by fitting. Give a file path to an XML or 3dm file to refit the patches!");
    cmd.addPlainString("filename", "File containing multipatch input (.xml).", filename);
    cmd.addReal  ("t","tolerance","Tolerance for identifing patch interfaces", tol);
    cmd.addReal  ("l","lambda_crv","lambda for the curve fitting", lambda_crv);
    cmd.addReal  ("L","lambda_srf","lambda for the surface fitting", lambda_srf);
    cmd.addInt   ("d", "degree", "Degree of B-splines for reparameterization", degree);
    cmd.addInt   ("k", "knots", "Number of interior knots for reparameterization", nknots);
    cmd.addInt   ("N", "npts", "Number of points for sampling", npts);
    cmd.addSwitch("plot", "plot results", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<>::uPtr mp0 = gsReadFile<>(filename);
    mp0->computeTopology(tol,true);

    gsInfo <<" Got "<< *mp0 <<" \n" ;


    // STEP 0: sample all patches to get linear patches.
    // (already the input is ASSUMED linear for now)
    // sample each patch on a grid
    // create the linear patches with the samples as coefficients + the topology of the initial one   
    gsMultiPatch<> mp = *mp0;
    gsMultiPatch<> mp_par = *mp0;
    gsMatrix<> ab, pts, eval, par_pts;
    gsVector<> a, b;
    gsVector<unsigned> np(mp.parDim());

    for (size_t p=0; p!=mp0->nPatches(); p++)
    {
        ab = mp0->patch(p).support();
        a = ab.col(0);
        b = ab.col(1);
        np.setConstant((std::ceil(std::pow(npts,1./mp.parDim()))));
        // Uniform parameters for evaluation
        pts = gsPointGrid(a, b, np);

        mp0->patch(p).eval_into(pts,eval);

        gsKnotVector<> kv0(0, 1, np[0]-2, 2, 1, 1);
        gsKnotVector<> kv1(0, 1, np[1]-2, 2, 1, 1);
        gsTensorBSplineBasis<2,real_t> bbasis(kv0,kv1);

        // Reparameterize points
        pts = parameterize_points2D(bbasis,eval);

        mp.patch(p) = give(*bbasis.makeGeometry(eval.transpose()));
        // Contains the parametric values of the points
        mp_par.patch(p) = give(*bbasis.makeGeometry(pts.transpose()));
    }
    if (plot) gsWriteParaview(mp,"mp",200,false);
    
    // STEP 1: Get curve network with merged linear interfaces
    gsInfo<<"Loading curve network...";
    mp.constructInterfaceRep();
    mp.constructBoundaryRep();
    auto & irep = mp.interfaceRep();
    auto & brep = mp.boundaryRep();
    gsDebug <<" irep "<< irep.size() <<" \n" ;
    gsDebug <<" brep "<< brep.size() <<" \n" ;

    // outputing...
    gsMultiPatch<> crv_net, iface_net, bnd_net;
    for (auto it = irep.begin(); it!=irep.end(); ++it)
    {
        iface_net.addPatch((*it->second));
        crv_net.addPatch((*it->second));
    }
    for (auto it = brep.begin(); it!=brep.end(); ++it)
    {
        bnd_net.addPatch((*it->second));
        crv_net.addPatch((*it->second));
    }

    if (plot) gsWriteParaview(iface_net,"iface_net",200);
    if (plot) gsWriteParaview(bnd_net,"bnd_net",200);
    if (plot) gsWriteParaview(crv_net,"crv_net",200);

    //end outputing
    gsInfo<<"Finished\n";
    
    //STEP 3: Fit curve network with B-splines of degree \a d and \a k interior knots
    //parametrizePts
    gsInfo<<"Making boundary representation with fitted interfaces...";
    gsKnotVector<> kv(0, 1, nknots, degree+1, 1, degree);
    gsBSplineBasis<> fbasis(kv);

    // pbdr contains the indices of the boundary curves in he first entry of the pair, and the coefficients in the second
    std::vector<std::pair<std::vector<gsMatrix<index_t>>,gsMatrix<real_t> > > pbdr(mp0->nPatches());
    crv_net.clear();
    gsMatrix<> uv, xyz;

    gsTensorBSplineBasis<2> sbasis(kv,kv);
    for (size_t p=0; p!=mp.nPatches(); p++)
    {
        pbdr[p].first.resize(4);
        pbdr[p].second.resize(sbasis.size(),mp.targetDim());
    }

    index_t k=0;
    for (auto it = irep.begin(); it!=irep.end(); ++it, k++)
    {
        gsMatrix<index_t> bndThis,bndOther;

        const gsBSpline<> & crv = static_cast<const gsBSpline<> &>(*it->second);
        xyz = crv.coefs().transpose();
        uv  = parameterize_points1D(xyz);

        gsFitting<> cfit(uv,xyz,fbasis);

        std::vector<index_t> prescribedDoFs;
        std::vector<gsMatrix<>> prescibedCoefs;
        prescribedDoFs.push_back(0);
        prescribedDoFs.push_back(fbasis.size()-1);
        prescibedCoefs.push_back(crv.coefs().row(0));
        prescibedCoefs.push_back(crv.coefs().row(crv.basis().size()-1));

        cfit.setConstraints(prescribedDoFs,prescibedCoefs);
        cfit.compute();

        sbasis.matchWith(it->first,sbasis,bndThis,bndOther);

        GISMO_ASSERT(bndThis .rows()==cfit.result()->coefs().rows(),"Coefficients and indices do not match!");
        GISMO_ASSERT(bndOther.rows()==cfit.result()->coefs().rows(),"Coefficients and indices do not match!");

        pbdr[it->first.first() .patch].first.at(it->first.first() .side()-1) = bndThis;
        for (index_t k=0; k!=bndThis .rows(); k++)
            pbdr[it->first.first() .patch].second.row(bndThis (k,0)) = cfit.result()->coefs().row(k);

        pbdr[it->first.second() .patch].first.at(it->first.second() .side()-1) = bndOther;
        for (index_t k=0; k!=bndOther.rows(); k++)
            pbdr[it->first.second() .patch].second.row(bndOther(k,0)) = cfit.result()->coefs().row(k);

        crv_net.addPatch( *cfit.result() );
    }
    for (auto it = brep.begin(); it!=brep.end(); ++it)
    {
        const gsBSpline<> & crv = static_cast<const gsBSpline<> &>(*it->second);
        xyz = crv.coefs().transpose();
        uv  = parameterize_points1D(xyz);
        gsFitting<> cfit(uv,xyz,fbasis);

        std::vector<index_t> prescribedDoFs;
        std::vector<gsMatrix<>> prescibedCoefs;
        prescribedDoFs.push_back(0);
        prescribedDoFs.push_back(fbasis.size()-1);
        prescibedCoefs.push_back(crv.coefs().row(0));
        prescibedCoefs.push_back(crv.coefs().row(crv.basis().size()-1));
        cfit.setConstraints(prescribedDoFs,prescibedCoefs);

        cfit.compute(lambda_crv);
        // cfit.computeErrors();
        // real_t tol = 1e-1;
        // if (cfit.maxPointError()> tol)
        //     gsWarn<<"Error of curve fit is large: "<<cfit.maxPointError()<<">"<<tol<<"\n";

        gsMatrix<index_t> bndThis = sbasis.boundary(it->first.side());
        pbdr[it->first .patch].first.at(it->first .side()-1) = bndThis;
        for (index_t k=0; k!=bndThis.rows(); k++)
            pbdr[it->first .patch].second.row(bndThis(k,0)) = cfit.result()->coefs().row(k);
        crv_net.addPatch( *cfit.result() );
    }

    if (plot) gsWriteParaview(crv_net,"crv_fit",200);
    gsInfo<<"Finished\n";

    std::vector<gsGeometry<>*> container(mp0->nPatches());
    //STEP 4: fit interior points of each patch with boundary constraints being the curves..
    gsInfo<<"Fitting surface...";
    for (size_t p=0; p!=pbdr.size(); p++)
    {
        std::vector<index_t> prescribedDoFs;
        std::vector<gsMatrix<>> prescibedCoefs;

        for (size_t s=0; s!=pbdr.at(p).first.size(); s++)
        {
            for (index_t k=0; k!=pbdr.at(p).first.at(s).size(); k++)
            {
                index_t index = (pbdr.at(p).first.at(s) )(k,0);
                prescribedDoFs.push_back(index);
                prescibedCoefs.push_back((pbdr.at(p).second).row(index));
            }
        }

        uv  = mp_par.patch(p).coefs().transpose();
        xyz = mp.patch(p).coefs().transpose();
        gsFitting<> sfit(uv,xyz,sbasis);
        sfit.setConstraints(prescribedDoFs,prescibedCoefs);

        sfit.compute(lambda_srf);
        // sfit.computeErrors();
        // real_t tol = 1e-1;
        // if (sfit.maxPointError()> tol)
        //     gsWarn<<"Error of surface fit is large: "<<sfit.maxPointError()<<">"<<tol<<"\n";

        container.at(p) = sfit.result()->clone().release();
    }

    gsMultiPatch<> mp_res(container,mp0->boundaries(),mp0->interfaces());
    if (plot) gsWriteParaview(mp_res,"final",200,true);
    gsWrite<>(mp_res,"final");
    gsInfo<<"Finished\n";

    return EXIT_SUCCESS;
}
