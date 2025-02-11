/** @file embeddedCurve_refinement.cpp

    @brief \a h -refinement of an embedded curve
           based on a 2D superdomain mesh grid. 

    Authors: C. Chianese (UniNa), H.M. Verhelst (TU Delft)
    Year: 2025
*/

#include <iostream>
#include <numeric>
#include <gismo.h>
#include <gsNurbs/gsKnotVector.h>

using namespace gismo;

template <class T>
gsVector<T> findKnot(const gsGeometry<T> &geometry, const gsKnotVector<T> &kv, const gsVector<T> &value,
                     const index_t i=0, const short_t dir=0, const index_t maxiter=50, const T epsilon=1e-8);

int main(int argc, char *argv[])
{
    //! [Embedding surface]
    std::vector<real_t> vect_s = {0,0,0,0.25,0.5,0.75,1,1,1};
    gsKnotVector<real_t> kv_u(vect_s,2); // knot vector along u direction
    gsKnotVector<real_t> kv_v(vect_s,2); // knot vector along v direction
    
    gsTensorBSplineBasis<2, real_t> basis_s(kv_u, kv_v); // B-spline bivariate basis

    gsVector<real_t,6> vec;
    vec << 0, 0.5, 1.5, 2.5, 3.5, 4;
    gsMatrix<real_t> coef_s(basis_s.size(), 3); // specify (x,y,z)- control net
    coef_s.setZero();
    index_t count = -1;

    for (index_t j = 0; j < vec.size(); j++)
    {
        for (index_t i = 0; i < vec.size(); i++)
        {
            count += 1;
            coef_s(count, 0) = vec.at(i);
            coef_s(count, 1) = vec.at(j);
        }
    }
    
    gsTensorBSpline<2, real_t>  surf(basis_s, coef_s); // build B-spline surface

    gsMultiPatch<> mp_surf;
    mp_surf.addPatch(surf); // add surface to multipatch

    gsWriteParaview(mp_surf,  "surf",  1000, true,  false);
    //! [Embedding surface]

    //! [Embedded curve]
    std::vector<real_t> vect_c = {0,0,0,0.25,0.5,0.75,1,1,1};
    gsKnotVector<real_t> kv_c(vect_c,2); // knot vector along curve

    gsBSplineBasis<> basis_c(kv_c); // B-spline univariate basis

    gsMatrix<real_t> coef_c(basis_c.size(), 2); // specify (u,v)- control net
    coef_c << 0, 0.116,
              0.214, 0.107,
              0.646, 0.248,
              0.646, 0.752,
              0.214, 0.892,
              0, 0.883;

    coef_c.col(0) *= kv_u.last();
    coef_c.col(1) *= kv_v.last();

    gsBSpline<> curve(basis_c, give(coef_c)); // build B-spline curve

    gsMultiPatch<> mp_curve;
    mp_curve.addPatch(curve); // add curve to multipatch

    gsWriteParaview(mp_curve,  "curve",  1000, true,  false);
     //! [Embedded curve]

    //! [Search for crossings between embedded curve and surface knot lines based on multiple γ-guesses]
    gsVector<> guess(3);
    guess << 0.25, 0.50, 0.75;
    index_t lu = guess.rows() * (kv_u.numElements() - 1);
    gsVector<real_t> crossing_u (lu);

    gsInfo << "\nIntersections of embedded curve with u-knot lines.\n";
    for (size_t i = 1; i != kv_u.numElements(); ++i) //loop over internal u-knots regardless of multiplicity
    {
        crossing_u.segment((i-1)*guess.rows(), guess.rows()) = findKnot(curve, kv_u, guess, i, 0);
    }

    index_t lv = guess.rows() * (kv_v.numElements() - 1);
    gsVector<real_t> crossing_v (lv);

     gsInfo << "\nIntersections of embedded curve with v-knot lines.\n";
    for (size_t j = 1; j != kv_v.numElements(); ++j) //loop over internal v-knots regardless of multiplicity
    {
        crossing_v.segment((j-1)*guess.rows(), guess.rows()) = findKnot(curve, kv_v, guess, j, 1);
    }
    //! [Search for crossings between embedded curve and surface knot lines starting from multiple γ-guesses]

    //! [Sort and filter out crossings]
    gsVector<> crossing (crossing_u.size() + crossing_v.size());
    crossing << crossing_u, crossing_v;

    std::sort(crossing.data(), crossing.data() + crossing.size());
    
    int frequencyNaN = std::count(crossing.begin(), crossing.end(), std::numeric_limits<real_t>::max());
    crossing.conservativeResize(crossing.size() - frequencyNaN); //remove NaN-like values

    real_t threshold = 1e-2;
    gsVector<real_t> diff (crossing.size());
    std::adjacent_difference(crossing.begin(), crossing.end(), diff.begin()); 
    gsVector<bool> mask(crossing.size());
    mask = diff.array() > threshold;
    //gsDebugVar(crossing);
    //gsDebugVar(diff);
    //gsDebugVar(mask);
    crossing = crossing.array() * mask.cast<real_t>().array();
    int frequency0 = std::count(crossing.begin(), crossing.end(), 0);
    std::remove(crossing.begin(), crossing.end(), 0);
    crossing.conservativeResize(crossing.size() - frequency0); //remove duplicated values within threshold
    
    gsInfo << "\nTotal intersections: " << crossing.transpose() << "\n";
    //! [Sort and filter out crossings]

    //![Add knots along the curve at found intersections]
    gsAsConstMatrix<real_t> kv_c_ori = kv_c.asMatrix();
    for (index_t i = 0; i != crossing.size(); ++i)
    {
        gsMatrix <> difference = kv_c_ori - gsEigen::MatrixXd::Constant(1, kv_c_ori.size(), 1.0) * crossing(i);
        if (std::all_of(difference.row(0).begin(), difference.row(0).end(), [threshold](real_t x) { return std::abs(x) > threshold; }))
        {
            curve.insertKnot(crossing(i));
        }
    }
    gsAsConstMatrix<real_t> kv_c_refined = curve.knots().get();
    gsInfo << "Knot vector of h-refined embedded curve: " << kv_c_refined << "\n";
    //![Add knots along the curve at found intersections]

    //![Get integration points along embedded curve]
    short_t maxdeg = std::max({curve.degree(),surf.degree(0),surf.degree(1)});
    short_t numNodes = static_cast<short_t>(std::ceil(0.5*maxdeg + 1)); //# of Gauss points per curve knot span 
    gsGaussRule<real_t> rule(numNodes);
    gsMatrix<real_t> gamma_ref = rule.referenceNodes(); //coordinates of Gauss points in reference domain
    gsVector<real_t> w_ref = rule.referenceWeights(); //associated weights

    gsMatrix<real_t> gamma (1,(curve.knots().uSize()-1)*numNodes);
    gsVector<real_t> w ((curve.knots().uSize()-1)*numNodes); 
    gsMatrix<> gamma_ks (curve.parDim(),numNodes);
    gsVector<> w_ks (numNodes);
    for (size_t i = 0; i != curve.knots().uSize()-1; ++i) //loop over curve knot spans
    {
        auto lower = kv_c_refined.block(0,curve.knots().multFirst()-1 + i,1,1);
        auto upper = kv_c_refined.block(0,curve.knots().multFirst()   + i,1,1);
        rule.mapTo(lower,upper,gamma_ks,w_ks);
        gamma.block(0,i*numNodes,1,numNodes) = gamma_ks; //map Gauss point coordinates to current knot span
        w.segment(i*numNodes,numNodes) = w_ks; //map associated weights to current knot span
    }
    //![Get integration points along embedded curve]

    gsDebugVar(gamma);
    gsDebugVar(w);

    return EXIT_SUCCESS;
}

template <class T>
gsVector<T> findKnot(const gsGeometry<T> &geometry, const gsKnotVector<T> &kv, const gsVector<T> &value,
                     const index_t i, const short_t dir, const index_t maxiter, const T epsilon)
{
    /* Function to find the crossing of an embedded curve with a knot line of an embedding surface.
    
       geometry:            embedded entity description in superdomain parameter space
       kv:                  superdomain knot vector in relevant parametric direction
       value:               vector of initial γ-guesses (γ is the embedded curve parameter)
       i:                   index of superdomain knot at which crossings are looked for 
       dir:                 superdomain parametric direction (0 ≡ u, 1 ≡ v)
       maxiter, epsilon:    root finding stopping parameters
    */

    GISMO_ASSERT(geometry.parDim() == 1, "The embedded geometry must be a curve (parameter dimension 1)");
    GISMO_ASSERT(dir == 0 || dir == 1, "The embedding geometry must be a surface");

    gsVector<T> result(value.rows());
    result.setConstant(std::numeric_limits<T>::max());
    gsMatrix<T> xn(1,1);
    T Dfxn , fxn;

    for(index_t k=0; k < value.rows(); k++) // loop over initial guesses
    {
        xn(0,0) = value(k);
        for (index_t n=0; n < maxiter; n++) // Newton iterations
        {
             Dfxn = geometry.deriv(xn)(dir,0);
             if (gsClose(Dfxn,0.0,std::numeric_limits<T>::epsilon()))
             {
                 gsInfo << "Guess " << k+1 << " :zero derivative, no crossing found.\n";
                 break;
             }
             fxn = geometry.eval(xn)(dir,0) - kv(i);
             if (math::abs(fxn) <= epsilon)
             {
                gsInfo << "Guess " << k+1 << " :crossing found after " << n << " iterations.\n";
                result(k) = xn(0,0);
                break;
             }
             xn(0,0) -= fxn/Dfxn;
             if (n == maxiter - 1)
             {
                gsInfo << "Guess" << k+1 << " :maxiter exceeded, no crossing found.\n";
             }
        }
    }
    return result;
}
