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
    tmp = bev.partialPivLu().solve(fev);//solve on element

    // find the i-th BS:
    gsMatrix<index_t> act = bb.active(pts.col(0));
    index_t c = std::lower_bound(act.data(), act.data()+act.size(), i) - act.data();
    GISMO_ASSERT(c<act.size(), "Problem with basis function index");
    return tmp.row(c);
}

/*
gsMatrix<T> gsQuasiInterpolate<T>::localIntpl(const gsTensorBasis<d,T> &bb,
                                              const gsFunction<T> &fun,
                                              index_t i,
                                              const gsMatrix<T> &ab)
{
    gsMatrix<T> bev, fev, pts, tmp;
    gsVector<index_t> nNodes = gsQuadrature::numNodes(bb,(T)1.0,1);
    gsQuadRule<T>  qRule     = gsQuadrature::get<T>(gsQuadrature::GaussLegendre,nNodes); //gsTPQuadRule ..

    // for(pt..)
    //{
    qRule.mapTo(ab, pt);//map point on element
    fun.eval_into(pts, fev);//evaluate function
    //}
    fev.transposeInPlace();

    //solve
    bev.transposeInPlace();// must be cwise
    tmp = bev.partialPivLu().solve(fev);//solve on element

    // find the i-th BS:
    gsMatrix<index_t> act = bb.active(pts.col(0)); //cwise..?
    index_t c = std::lower_bound(act.data(), act.data()+act.size(), i) - act.data();
    GISMO_ASSERT(c<act.size(), "Problem with basis function index");
    return tmp.row(c);
}
*/

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
    else
        return localIntpl(bb,fun,i,bb.elementInSupportOf(i));
}


template<typename T>
void gsQuasiInterpolate<T>::Taylor(const gsBasis<T> &bb, const gsFunction<T> &fun, const int &r, gsMatrix<T> & coefs)
{
    const gsBSplineBasis<T> & b = dynamic_cast<const gsBSplineBasis<T> &>(bb);
    // ONLY 1D

    const gsKnotVector<T> & kv = b.knots();
    int deg = b.degree();
    gsMatrix<T> xj = b.anchors();

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
    const gsBSplineBasis<T> & b = dynamic_cast<const gsBSplineBasis<T> &>(bb);
    // ONLY 1D

    const gsKnotVector<T> & kv = b.knots();
    const int n = b.size();
    //gsDebugVar(kv);

    coefs.resize(n, fun.targetDim());

    gsMatrix<T> knots(1,kv.size());
    for(unsigned int i=0; i<kv.size(); i++)
        knots(i) = kv[i];

    gsMatrix<T> TmpCoefs;

    int type = 0;
    if(specialCase)
    {
        type = b.degree();
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
            knotsAvg(i) = (kv[i]+kv[i+1])/2;
        gsMatrix<T> TmpCoefsAvg;
        fun.eval_into(knotsAvg, TmpCoefsAvg);

        coefs.row(0) = TmpCoefs.col(0);
        for(int i=1; i<n-1; i++)
        {
            // formula: (-a + 4b - c)/2;
            coefs.row(i).noalias() =
                ( - TmpCoefs.col(i+1)
                  + 4 * TmpCoefsAvg.col(i+1)
                  - TmpCoefs.col(i+2) ) / 2;
        }
        coefs.row(n-1) = TmpCoefs.col(n);
        break;
    }
    case 3: //5-point cubic (section 8.2.3 of "Spline methods (Lyche Morken)")
    {
        fun.eval_into(knots, TmpCoefs);
        gsMatrix<T> knotsAvg(1, kv.size()-1);
        for(unsigned int i=0; i<kv.size()-1; i++)
            knotsAvg(i) = (kv[i]+kv[i+1])/2;
        gsMatrix<T> TmpCoefsAvg;
        fun.eval_into(knotsAvg, TmpCoefsAvg);

        coefs.row(0) = TmpCoefs.col(3).transpose();

        // formula: (- 5a + 40b - 24c + 8d - e)/18
        coefs.row(1).noalias() =
                ( -  5 * TmpCoefs.col(3)
                  + 40 * TmpCoefsAvg.col(3)
                  - 24 * TmpCoefs.col(4)
                  +  8 * TmpCoefsAvg.col(4)
                  - TmpCoefs.col(5) ) / 18;

        for(int i=2; i<n-2; i++)
        {
            // formula: (a - 8b + 20c - 8d + e)/6
            coefs.row(i).noalias() =
                    ( TmpCoefs.col(i+1)
                      -  8 * TmpCoefsAvg.col(i+1)
                      + 20 * TmpCoefs.col(i+2)
                      -  8 * TmpCoefsAvg.col(i+2)
                      + TmpCoefs.col(i+3) ) / 6;
        }

        // formula: (- a + 8b - 24c + 40d - 5e)/18
        coefs.row(n-2).noalias() =
                ( - TmpCoefs.col(n-2)
                  +  8 * TmpCoefsAvg.col(n-2)
                  - 24 * TmpCoefs.col(n-1)
                  + 40 * TmpCoefsAvg.col(n-1)
                  -  5 * TmpCoefs.col(n) ) / 18;

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
        points.at(k) = a + (T)k/(n-1) * (b-a);
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
    //assert b.domainDim()==fun.domainDim()
    gsMatrix<T> cf;
    index_t n = b.size();
    index_t dim = fun.targetDim();
    result.resize(n,dim);

    for (index_t i = 0; i!=n; ++i)
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
