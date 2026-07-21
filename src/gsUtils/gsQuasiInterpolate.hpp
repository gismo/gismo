/** @file gsQuasiInterpolate.hpp

    @brief Different Quasi-Interpolation Schemes based on the article
    "Spline methods (Lyche Morken)"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Haberleitner, A. Mantzaflaris, H. Verhelst
*/

#pragma once

#include <gsUtils/gsCombinatorics.h>

namespace gismo {

template<typename T>
gsMatrix<T> gsQuasiInterpolate<T>::localIntpl(const gsBasis<T> &bb,
                                              const gsFunction<T> &fun,
                                              index_t i,
                                              const gsMatrix<T> &ab)
{
    gsMatrix<T> bev, fev, pts, tmp;
    gsVector<index_t> nNodes = gsQuadrature::numNodes(bb,(T)1.0,1);
    gsQuadRule<T>  qRule     = gsQuadrature::get<T>(gsQuadrature::GaussLegendre,nNodes);
    qRule.mapTo(ab, pts);//map points on element
    bb .eval_into(pts, bev);//evaluate basis
    fun.eval_into(pts, fev);//evaluate function
    bev.transposeInPlace();
    fev.transposeInPlace();
    tmp = bev.fullPivLu().solve(fev);//solve on element

    gsMatrix<T> interpolatedFev = bev.transpose() * tmp; // interpolated values at quad points
    gsMatrix<T> error = fev - interpolatedFev; // local error
   
    // find the i-th BS:
    gsMatrix<index_t> act = bb.active(pts.col(0));
    index_t c = std::lower_bound(act.data(), act.data()+act.size(), i) - act.data();
    GISMO_ASSERT(c<act.size(), "Problem with basis function index");
    return tmp.row(c);
}

template<typename T>
template<short_t d>
gsMatrix<T> gsQuasiInterpolate<T>::localIntpl(const gsHTensorBasis<d,T> &bb,
                                              const gsFunction<T> &fun,
                                              index_t i)
{
    index_t lvl = bb.levelOf(i);
    index_t j = bb.flatTensorIndexOf(i);
    return localIntpl(bb.tensorLevel(lvl),fun,j, bb.elementInSupportOf(i)); // uses the H-grid element implementation
    //return localIntpl(bb.tensorLevel(lvl),fun,j); // uses the central element implementation
}

template<typename T>
gsMatrix<T> gsQuasiInterpolate<T>::localIntpl(const gsBasis<T> &bb,
                                              const gsFunction<T> &fun,
                                              index_t i)
{
    if (const gsHTensorBasis<1,T>* b = dynamic_cast<const gsHTensorBasis<1,T>* >(&bb))
        return localIntpl(*b,fun,i);
    if (const gsHTensorBasis<2,T>* b = dynamic_cast<const gsHTensorBasis<2,T>* >(&bb))
         return localIntpl(*b,fun,i);
    if (const gsHTensorBasis<3,T>* b = dynamic_cast<const gsHTensorBasis<3,T>* >(&bb))
        return localIntpl(*b,fun,i);
    if (const gsHTensorBasis<4,T>* b = dynamic_cast<const gsHTensorBasis<4,T>* >(&bb))
        return localIntpl(*b,fun,i);        
    // If it is a gsRationalTHBSplineBasis, we check the source
    if (const gsHTensorBasis<1, T>* b = dynamic_cast<const gsHTensorBasis<1,T>* >(&bb.source())) 
        return localIntpl(*b,fun,i);
    if (const gsHTensorBasis<2, T>* b = dynamic_cast<const gsHTensorBasis<2,T>* >(&bb.source())) 
        return localIntpl(*b,fun,i);
    if (const gsHTensorBasis<3, T>* b = dynamic_cast<const gsHTensorBasis<3,T>* >(&bb.source()))
        return localIntpl(*b,fun,i);
    else
        return localIntpl(bb,fun,i,bb.elementInSupportOf(i));
}

template<typename T>
gsMatrix<T> gsQuasiInterpolate<T>::localL2(const gsBasis<T> &bb,
                                            const gsFunction<T>  &fun,
                                            index_t i,                                    
                                            const gsMatrix<T> &ab)
{
    gsMatrix<T> bev, fev, pts, tmp;
    gsVector<T> weights;
    gsVector<index_t> nNodes = gsQuadrature::numNodes(bb,(T)1.0,1);
    gsQuadRule<T>  qRule     = gsQuadrature::get<T>(gsQuadrature::GaussLobatto,nNodes);
    qRule.mapTo(ab.col(0),ab.col(1), pts, weights);//map points and weights on element (for quadrature)

    bb .eval_into(pts, bev);//evaluate basis on quadrature points (numbasis*pts)
    fun.eval_into(pts, fev);//evaluate function on quadrature points (1*pts)
    gsMatrix<T> M = bev * weights.asDiagonal() * bev.transpose(); // Mass matrix
    gsMatrix<T> RHS = bev * weights.asDiagonal() * fev.transpose();
    tmp = M.fullPivLu().solve(RHS); //solve on element


    // find the i-th BS:
    gsMatrix<index_t> act = bb.active(pts.col(0)); // the same basis functions will be active for the element in consideration!
    index_t c = std::lower_bound(act.data(), act.data()+act.size(), i) - act.data();
    GISMO_ASSERT(c<act.size(), "Problem with basis function index");
    return tmp.row(c);
}

template<typename T>
template<short_t d>
gsMatrix<T> gsQuasiInterpolate<T>::localL2(const gsHTensorBasis<d,T> &bb,   
                                            const gsFunction<T>  &fun,
                                            index_t i)
{
    index_t lvl = bb.levelOf(i);
    index_t j = bb.flatTensorIndexOf(i);
    return localL2(bb.tensorLevel(lvl),fun,j,bb.elementInSupportOf(i)); // uses the H-grid element implementation
}

template<typename T>
gsMatrix<T> gsQuasiInterpolate<T>::localL2(const gsBasis<T> &bb,           
                                            const gsFunction<T>  &fun,
                                            index_t i)
{
    if (const gsHTensorBasis<1,T>* b = dynamic_cast<const gsHTensorBasis<1,T>* >(&bb))
        return localL2(*b,fun,i);
    if (const gsHTensorBasis<2,T>* b = dynamic_cast<const gsHTensorBasis<2,T>* >(&bb))
        return localL2(*b,fun,i);
    if (const gsHTensorBasis<3,T>* b = dynamic_cast<const gsHTensorBasis<3,T>* >(&bb))
        return localL2(*b,fun,i);
    if (const gsHTensorBasis<4,T>* b = dynamic_cast<const gsHTensorBasis<4,T>* >(&bb))
        return localL2(*b,fun,i);
    else
        return localL2(bb,fun,i,bb.elementInSupportOf(i));
}

template <typename T>
void gsQuasiInterpolate<T>::localL2(const gsBasis<T> &b,
                                    const gsFunction<T>  &fun,
                                    gsMatrix<T> &result)
{
    GISMO_ASSERT(b.domainDim()==fun.domainDim(),"Domain dimensions should be equal");
    //assert b.domainDim()==fun.domainDim()
    gsMatrix<>  cf;
    index_t n = b.size();
    index_t dim = fun.targetDim();
    result.resize(n,dim);

#   pragma omp parallel for private(cf)
    for (index_t i = 0; i<n; ++i)
    {
        cf = localL2(b,fun,i);
        result.row(i) = cf;
    }
}

template<typename T>
void gsQuasiInterpolate<T>::Taylor(const gsBasis<T> &bb, const gsFunction<T> &fun, const int &r, gsMatrix<T> & coefs)
{
    const gsBSplineBasis<T>*b = dynamic_cast<const gsBSplineBasis<T> *>(&bb); // cast bb to a gsBSplineBasis
    GISMO_ASSERT(b != nullptr, "Basis should be a gsBSplineBasis"); // assertion to ensure that b is a gsBSplineBasis
    // ONLY 1D

    const gsKnotVector<T> & kv = b->knots();
    int deg = b->degree();
    gsMatrix<T> xj = b->anchors();

    int n = xj.size();
    int dim = fun.targetDim();
    coefs.resize(n,dim);

    std::vector<gsMatrix<T> > derivs;
    fun.evalAllDers_into(xj, r, derivs);


    gsMatrix<T> val;
    std::vector<T> knots;
    for(int j=0; j<n; j++)
    {
        val.setZero(1,dim);
        knots.clear();
        for(int i=j+1; i<=j+deg; i++)
            knots.push_back(kv[i]);

        for(int k=0; k<=r; k++) // (r+1) nodes
        {
            const T factor1 = derivProd(knots, deg-k, xj(j)); //coeff
            for(int i=0; i<dim; i++)
            {
                const T factor2 = derivs[k](i,j); //node
                val(i) += std::pow(-1.0,k) * factor1 * factor2;
            }
        }
        val /= factorial(deg);
        coefs.row(j) = val;
    }
}

template<typename T>
index_t gsQuasiInterpolate<T>::derivRow(const gsVector<index_t> &alpha, short_t d, index_t comp)
{
    const index_t m = alpha.sum(); // total derivative order
    if (0==m) // value
        return comp;

    if (1==m) // first derivatives: [d_0, d_1, ..., d_{d-1}]
    {
        for (short_t dir = 0; dir!=d; ++dir)
            if (1==alpha[dir]) return comp*d + dir;
        GISMO_ERROR("Invalid multi-index");
    }

    if (2==m) // second derivatives: pure [d_00,..,d_{d-1,d-1}] then mixed d_ab (a<b) lex.
    {
        const index_t str = d*(d+1)/2; // block size per component
        for (short_t dir = 0; dir!=d; ++dir)
            if (2==alpha[dir]) return comp*str + dir; // pure
        // mixed: find the two directions a<b with alpha==1 and their lex. pair index
        index_t pair = 0;
        for (short_t a = 0; a!=d; ++a)
            for (short_t b = a+1; b!=d; ++b, ++pair)
                if (1==alpha[a] && 1==alpha[b])
                    return comp*str + d + pair;
        GISMO_ERROR("Invalid multi-index");
    }

    // m>=3 : composition (lexicographic) order, see nextComposition
    const index_t str = numCompositions(m, d); // block size per component
    gsVector<index_t> comp_it;
    firstComposition(m, d, comp_it);
    index_t pos = 0;
    do
    {
        if (comp_it==alpha) return comp*str + pos;
        ++pos;
    } while (nextComposition(comp_it));
    GISMO_ERROR("Invalid multi-index");
}

template<typename T>
template<short_t d>
gsMatrix<T> gsQuasiInterpolate<T>::localTaylor(const gsTensorBSplineBasis<d,T> &b,
                                              const gsFunction<T> &fun,
                                              const int &r,
                                              index_t j)
{
    const index_t dim = fun.targetDim();

    // Per-direction data: degree, the deg roots of rho_j, the anchor x_j
    // and the maximal derivative order used in that direction.
    const gsVector<index_t,d> jj = b.tensorIndex(j);
    gsVector<index_t,d> deg, rdir;
    std::vector<std::vector<T> > knots(d);
    gsMatrix<T> point(d,1);
    for (short_t dir = 0; dir!=d; ++dir)
    {
        deg[dir]  = b.degree(dir);
        rdir[dir] = math::min(r, static_cast<int>(deg[dir])); // r <= p per direction
        const gsKnotVector<T> & kv = b.knots(dir);
        for (index_t q = jj[dir]+1; q<=jj[dir]+deg[dir]; ++q)
            knots[dir].push_back(kv[q]);
        point(dir,0) = b.component(dir).anchors()(0, jj[dir]);
    }

    // Mixed partials of the target function up to the max total order needed
    // (sum of the per-direction orders), evaluated once at the single anchor
    // point x_j. Note: this requires \a fun to provide derivatives up to that
    // order. Tensor B-spline (and hierarchical) geometries do (see
    // gsTensorBasis::evalAllDers_into, which handles arbitrary order); an
    // analytic gsFunctionExpr is limited to order 2, hence usable here only when
    // the total order does not exceed 2.
    const index_t maxOrder = rdir.sum();
    std::vector<gsMatrix<T> > derivs;
    fun.evalAllDers_into(point, maxOrder, derivs);

    // Tensor product of the univariate Taylor QIs: sum over the derivative
    // multi-index k, with 0 <= k[dir] <= rdir[dir]. nextLexicographic uses an
    // exclusive upper bound, hence rdir+1.
    const gsVector<index_t,d> bound = rdir + gsVector<index_t,d>::Ones();
    gsMatrix<T> val;
    val.setZero(1,dim);
    gsVector<index_t,d> k = gsVector<index_t,d>::Zero();
    do
    {
        const index_t order = k.sum();

        T factor1 = (T)1;
        for (short_t dir = 0; dir!=d; ++dir)
            factor1 *= derivProd(knots[dir], deg[dir]-k[dir], point(dir,0));

        const T sign = (0==(order%2)) ? (T)1 : (T)(-1);
        for (index_t i = 0; i!=dim; ++i)
        {
            const T factor2 = derivs[order]( derivRow(k, d, i), 0 );
            val(i) += sign * factor1 * factor2;
        }
    } while ( nextLexicographic(k, bound) );

    T denom = (T)1;
    for (short_t dir = 0; dir!=d; ++dir)
        denom *= factorial(deg[dir]);
    val /= denom;
    return val;
}

template<typename T>
template<short_t d>
gsMatrix<T> gsQuasiInterpolate<T>::localTaylor(const gsHTensorBasis<d,T> &bb,   
                                                const gsFunction<T>  &fun,
                                                const int &r,
                                                index_t i)
{
    index_t lvl = bb.levelOf(i);
    index_t j = bb.flatTensorIndexOf(i);
    return localTaylor<d>(bb.tensorLevel(lvl),fun,r,j); // uses the H-grid element implementation
}

template<typename T>
gsMatrix<T> gsQuasiInterpolate<T>::localTaylor(const gsBasis<T> &bb,
                                              const gsFunction<T> &fun,
                                              const int &r,
                                              index_t i)
{
    // Hierarchical tensor bases: dispatch on the parameter dimension.
    if (const gsHTensorBasis<1,T>* b = dynamic_cast<const gsHTensorBasis<1,T>* >(&bb))
        return localTaylor(*b,fun,r,i);
    if (const gsHTensorBasis<2,T>* b = dynamic_cast<const gsHTensorBasis<2,T>* >(&bb))
        return localTaylor(*b,fun,r,i);
    if (const gsHTensorBasis<3,T>* b = dynamic_cast<const gsHTensorBasis<3,T>* >(&bb))
        return localTaylor(*b,fun,r,i);
    if (const gsHTensorBasis<4,T>* b = dynamic_cast<const gsHTensorBasis<4,T>* >(&bb))
        return localTaylor(*b,fun,r,i);
    // Plain tensor B-spline bases: dispatch on the parameter dimension.
    if (const gsTensorBSplineBasis<1,T>* b = dynamic_cast<const gsTensorBSplineBasis<1,T>* >(&bb))
        return localTaylor<1>(*b,fun,r,i);
    if (const gsTensorBSplineBasis<2,T>* b = dynamic_cast<const gsTensorBSplineBasis<2,T>* >(&bb))
        return localTaylor<2>(*b,fun,r,i);
    if (const gsTensorBSplineBasis<3,T>* b = dynamic_cast<const gsTensorBSplineBasis<3,T>* >(&bb))
        return localTaylor<3>(*b,fun,r,i);
    if (const gsTensorBSplineBasis<4,T>* b = dynamic_cast<const gsTensorBSplineBasis<4,T>* >(&bb))
        return localTaylor<4>(*b,fun,r,i);
    GISMO_ERROR("localTaylor: unsupported basis type/dimension");
}


template<typename T>
void gsQuasiInterpolate<T>::localTaylor(const gsBasis<T> &b,
                                       const gsFunction<T> &fun,
                                       const int &r,
                                       gsMatrix<T> & result)
{
    // GISMO_ASSERT(b.domainDim()==fun.domainDim(),"Domain dimensions should be equal");
    // //assert b.domainDim()==fun.domainDim()
    gsMatrix<T> cf;
    index_t n = b.size();
    index_t dim = fun.targetDim();
    result.resize(n,dim);

#   pragma omp parallel for private(cf)
    for (index_t i = 0; i<n; ++i)
    {
        cf = localTaylor(b,fun,r,i);
        result.row(i) = cf;
    }
}

template<typename T>
void gsQuasiInterpolate<T>::Taylor2D(const gsBasis<T> &bb, const gsFunction<T> &fun, const int &r, gsMatrix<T> & coefs)
{
    // Superseded by the dimension-independent localTaylor. Kept for API
    // compatibility; delegates to the general per-coefficient Taylor QI (which
    // also fixes the old 2D aliasing and targetDim-stride bugs).
    localTaylor(bb, fun, r, coefs);
}

template<typename T> gsMatrix<T>
gsQuasiInterpolate<T>::Schoenberg(const gsBasis<T> &b,
                                  const gsFunction<T> &fun,
                                  index_t i)
{
    gsMatrix<T> res;
    fun.eval_into(b.anchor(i), res);
    res.transposeInPlace();
    return res;
}


template<typename T>
void gsQuasiInterpolate<T>::EvalBased(const gsBasis<T> &bb, const gsFunction<T> &fun, const bool specialCase, gsMatrix<T> &coefs)
{
    const gsBSplineBasis<T>* b = dynamic_cast<const gsBSplineBasis<T> *>(&bb); // cast bb to a gsBSplineBasis
    GISMO_ASSERT(b != nullptr, "Basis should be a gsBSplineBasis"); // assertion to ensure that b is a gsBSplineBasis
    // ONLY 1D

    const gsKnotVector<T> & kv = b->knots();
    const int n = b->size();
    //gsDebugVar(kv);

    coefs.resize(n, fun.targetDim());

    gsMatrix<T> knots(1,kv.size());
    for(unsigned int i=0; i<kv.size(); i++)
        knots(i) = kv[i];

    gsMatrix<T> TmpCoefs;

    int type = 0;
    if(specialCase)
    {
        type = b->degree();
        GISMO_ASSERT( (type == 1 || type == 2 || type == 3),
                      "quasiInterpolateEvalBased is implemented for special cases of deg 1, 2 or 3!");
    }

    switch(type)
    {
    case 1: //piecewise linear (section 8.2.1 of "Spline methods (Lyche Morken)")
    {
        fun.eval_into(knots, TmpCoefs);
        for(int i=0; i<n; i++)
            coefs.row(i) = TmpCoefs.col(i+1).transpose();
        break;
    }
    case 2: //3-point quadratic (section 8.2.2 of "Spline methods (Lyche Morken)")
    {
        fun.eval_into(knots, TmpCoefs);
        gsMatrix<T> knotsAvg(1, kv.size()-1);
        for(unsigned int i=0; i<kv.size()-1; i++)
            knotsAvg(i) = (kv[i]+kv[i+1]) / (T)(2);
        gsMatrix<T> TmpCoefsAvg;
        fun.eval_into(knotsAvg, TmpCoefsAvg);

        coefs.row(0) = TmpCoefs.col(0);
        for(int i=1; i<n-1; i++)
        {
            // formula: (-a + 4b - c)/2;
            coefs.row(i).noalias() =
                ( - TmpCoefs.col(i+1)
                  + 4 * TmpCoefsAvg.col(i+1)
                  - TmpCoefs.col(i+2) ) / (T)(2);
        }
        coefs.row(n-1) = TmpCoefs.col(n);
        break;
    }
    case 3: //5-point cubic (section 8.2.3 of "Spline methods (Lyche Morken)")
    {
        fun.eval_into(knots, TmpCoefs);
        gsMatrix<T> knotsAvg(1, kv.size()-1);
        for(unsigned int i=0; i<kv.size()-1; i++)
            knotsAvg(i) = (kv[i]+kv[i+1]) / (T)(2);
        gsMatrix<T> TmpCoefsAvg;
        fun.eval_into(knotsAvg, TmpCoefsAvg);

        coefs.row(0) = TmpCoefs.col(3).transpose();

        // formula: (- 5a + 40b - 24c + 8d - e)/18
        coefs.row(1).noalias() =
                ( -  5 * TmpCoefs.col(3)
                  + 40 * TmpCoefsAvg.col(3)
                  - 24 * TmpCoefs.col(4)
                  +  8 * TmpCoefsAvg.col(4)
                  - TmpCoefs.col(5) ) / (T)(18);

        for(int i=2; i<n-2; i++)
        {
            // formula: (a - 8b + 20c - 8d + e)/6
            coefs.row(i).noalias() =
                    ( TmpCoefs.col(i+1)
                      -  8 * TmpCoefsAvg.col(i+1)
                      + 20 * TmpCoefs.col(i+2)
                      -  8 * TmpCoefsAvg.col(i+2)
                      + TmpCoefs.col(i+3) ) / (T)(6);
        }

        // formula: (- a + 8b - 24c + 40d - 5e)/18
        coefs.row(n-2).noalias() =
                ( - TmpCoefs.col(n-2)
                  +  8 * TmpCoefsAvg.col(n-2)
                  - 24 * TmpCoefs.col(n-1)
                  + 40 * TmpCoefsAvg.col(n-1)
                  -  5 * TmpCoefs.col(n) ) / (T)(18);

        coefs.row(n-1) = TmpCoefs.col(n);

        break;
    }
    default: //if none of the special cases, (Theorem 8.7 and Lemma 9.7 of "Spline methods (Lyche Morken)")
        gsMatrix<T> xik, weights;
        for(int i=0; i<n; i++)
        {
            //look for the greatest subinterval to chose the interpolation points from
            int gsi = greatestSubInterval(kv, i, i+kv.degree());

            //compute equally distributed points in greatest subinterval
            distributePoints(kv[gsi], kv[gsi+1], kv.degree()+1, xik);
            //gsDebugVar(xik);

            //compute the factors omega_{ik}
            computeWeights(xik, kv, i+1, weights);

            //compute the coefficients lambda_i of the quasi-interpolant
            coefs.row(i) = computeControlPoints(weights, fun, xik);
        }
        break;
    }
}


template<typename T>
T gsQuasiInterpolate<T>::derivProd(const std::vector<T> &zeros, const int &order, const T &x)
{
    if(order == 0) // value
        return (x - gsAsConstMatrix<T,1>(zeros).array()).prod();

    std::vector<T> tmpZeros;

    if(order == 1) // first derivative
    {
        const index_t n1 = zeros.size() - 1;
        tmpZeros = zeros;
        T val = 0;
        for(typename std::vector<T>::iterator it = tmpZeros.begin(); it!=tmpZeros.end(); ++it)
        {
            std::iter_swap(it, tmpZeros.end()-1);
            val += (x - gsAsConstMatrix<T,1>(tmpZeros,1,n1).array()).prod(); // eval product
            std::iter_swap(it, tmpZeros.end()-1);
        }
        return val;
    }

    // Reccursion for higher order derivatives
    const int n = zeros.size();
    T val = 0;
    for(int i=0; i!=n; i++)
    {
        tmpZeros = zeros;
        tmpZeros.erase(tmpZeros.begin() + i);
        val += derivProd(tmpZeros, order-1, x);
    }
    return val;
}


template<typename T>
void gsQuasiInterpolate<T>::distributePoints(T a, T b, int n, gsMatrix<T> &points)
{
    points.resize(1,n);
    for(int k=0; k<n; k++)
        points.at(k) = a + (T)k/(T)(n-1) * (b-a);
}


template<typename T>
void gsQuasiInterpolate<T>::computeWeights(const gsMatrix<T> &points, const gsKnotVector<T> &knots, const int &pos, gsMatrix<T> &weights)
{
    const int deg = knots.degree();
    weights.resize(1,deg+1);

    gsMatrix<T> pointsReduced(1,deg);
    gsVector<int> indices(deg);
    for(int k=0; k<deg+1; k++)
    {
        T constant = (T)factorial(deg);

        //get the list of points without 'points(k)' and
        //multiply the values of the numerator together, to get the total constant
        for(int i=0; i<k; i++)
        {
            pointsReduced(i) = points(i);
            constant *= points(k) - points(i);
        }
        for(int i=k+1; i<deg+1; i++)
        {
            pointsReduced(i-1) = points(i);
            constant *= points(k) - points(i);
        }

        //get all permutations of the list
        for(int i=0; i<deg; i++)
            indices[i]=i;

        T sum = 0;
        do
        {
            //compute the product of values for this permutation
            T factor = 1;
            for(int i=0; i<deg; i++)
            {
                factor *= (knots[pos+indices[i]] - pointsReduced(i));
            }
            sum += factor; //sum up the products for each permutation
        }
        while(std::next_permutation(indices.data(), indices.data()+deg));

        weights(k) = sum / constant;
    }
}


template<typename T>
gsMatrix<T> gsQuasiInterpolate<T>::computeControlPoints(const gsMatrix<T> &weights, const gsFunction<T> &fun, const gsMatrix<T> &xik)
{
    gsMatrix<T> funValues;
    fun.eval_into(xik, funValues);
    return (weights.asDiagonal() * funValues.transpose()).colwise().sum();
}


template<typename T>
int gsQuasiInterpolate<T>::greatestSubInterval(const gsKnotVector<T> &knots, const int &posStart, const int &posEnd)       //ToDo: move to gsKnotVector
{
    const int diff = posEnd-posStart;
    T maxDist=0.0;
    int maxInd = posStart;
    for(int i=1; i<diff+1; i++)
    {
        const T dist = knots[posStart+i+1] - knots[posStart+i];
        if(dist > maxDist)
        {
            maxDist = dist;
            maxInd = posStart+i;
        }
    }
    return maxInd;
}


template<typename T>
void gsQuasiInterpolate<T>::localIntpl(const gsBasis<T> &b,
                                       const gsFunction<T> &fun,
                                       gsMatrix<T> & result)
{
    GISMO_ASSERT(b.domainDim()==fun.domainDim(),"Domain dimensions should be equal");
    //assert b.domainDim()==fun.domainDim()
    gsMatrix<T> cf;
    index_t n = b.size();
    index_t dim = fun.targetDim();
    result.resize(n,dim);

#   pragma omp parallel for private(cf)
    for (index_t i = 0; i<n; ++i)
    {
        cf = localIntpl(b,fun,i);
        result.row(i) = cf;
    }
}

template<typename T>
void gsQuasiInterpolate<T>::Schoenberg(const gsBasis<T> &b,
                                       const gsFunction<T> &fun,
                                       gsMatrix<T> & result)
{
    //assert b.domainDim()==fun.domainDim()
    gsMatrix<T> cf;
    index_t n = b.size();
    index_t dim = fun.targetDim();
    result.resize(n,dim);

    for (index_t i = 0; i!=n; ++i)
    {
        cf = Schoenberg(b,fun,i);
        result.row(i) = cf;
    }
}


} // gismo
