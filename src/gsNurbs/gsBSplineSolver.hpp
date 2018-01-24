/** @file gsBSplineSolver.hpp

    @brief Provides implementation of classes and functions to solve
    equations involving B-splines.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, A. Mantzaflaris
*/

#pragma once

# include <gsNurbs/gsBSpline.h>

#include <numeric>

namespace gismo {

template<class T>
int gsBSplineSolver<T>::insertKnot(int mu) 
{
    // find span (in case mu is not tight)
    mu = math::max(mu,m_d);

    while ( x>=m_t[mu+1]) mu++;

    // make room for one coefficient at position mu
    m_c.insert(m_c.begin()+mu,0);
    //for ( int i=m_n; i>mu; i--) m_c[i] = m_c[i-1];

    // Compute new coefficients
    for ( int i=mu; i>=mu-m_d+1; i--)
    {
        const T alpha = (x-m_t[i])/(m_t[i+m_d]-m_t[i]);
        m_c[i] = (1.0-alpha) * m_c[i-1] + alpha * m_c[i];
    }

    // set last knot
    m_t[m_n+m_d] = m_t[m_n+m_d-1];

    //insert the knot in the knot sequence
    m_t.insert(m_t.begin()+mu+1,x);
    // for ( int i=m_n; i>mu+1; --i)
    //     m_t[i] = m_t[i-1];
    // m_t[mu+1] = x;

    m_n++;
    return mu;
}

template<class T>
bool gsBSplineSolver<T>::nextRoot()
{
    // Do until maxn knots are inserted
    while  ( m_n < maxn )
    {
        // Find sign change
        while( m_k<m_n && (m_c[m_k-1]*m_c[m_k] > 0) )
        {
            ++m_k;
            //gsLog<<"m_c[m_k]= "<< m_c[m_k] <<" ( m_k="<<m_k<<")\n";
        }

        // No sign change ?
        if ( m_k >= m_n )
            return false;
        //gsLog<<"Sign change at "<< m_c[m_k]<<", "<< m_c[m_k+1] <<" (position "<<m_k<<")\n";

        // Interval converged ?
        const T diff = m_t[m_k+m_d]-m_t[m_k];
        if (diff<eps)
        {
            x = m_t[m_k];// root found
            //gsDebug <<"diff converged: Root found at "<< m_k<<", val="<< x <<"m_n="<<m_n <<"\n";
            break;
        }

        // Horizontal alignment ?
        const T cdiff = m_c[m_k-1]-m_c[m_k];
        if (math::abs(cdiff)<eps)
        {
            ++m_k;
            continue;
        }

        // compute line slope
        const T lambda = m_c[m_k-1] / cdiff ;

        // compute intersection
        x = (std::accumulate(m_t.begin()+m_k, m_t.begin()+m_k+m_d, T(0.0) )
             + lambda*diff ) / T(m_d);
        //x = ( m_t.segment(m_k, m_d).sum() + lambda*diff ) / T(m_d);

        // Stopping criterion
        const T e = math::max(x, m_t[m_k+m_d-1]) - math::min(x, m_t[m_k+1] );
        if (e < eps)
        {
            //gsDebug <<"eps: Root found after knot "<< m_k <<" is "<<m_t[m_k]<<", val="<< x <<"\n";
            break;
        }

        insertKnot(m_k);
        //gsLog<<"Inserted knot, m_k="<< m_k <<", m_t[m_k]="<<m_t[m_k]<<", m_n="<< m_n<<"\n";

    } // end for


    //gsLog<<"Loop finished with "<< x <<" m_k="<<m_k<<", m_c[m_k]="<<m_c[m_k]<<", m_t[m_k]="<<m_t[m_k]<< ", m_n="<< m_n<<"\n";

    if ( m_n>=maxn )
    {
        gsWarn<<"gsBSplineRoot: Maximum number of knots reached: "<< m_n << "\n";
        //gsDebug<< "m_c:\n"<< m_c.transpose()<<"\n";
        //gsDebug<< "m_t:\n"<< m_t.transpose()<<"\n";
        return false;
    }

    m_k++; // next root

    // buggy:   while ( m_k<m_n && math::fabs( m_t[m_k]-x)<eps ) m_k++;

    //gsDebug<<"next m_k="<< m_k <<" is "<< m_t[m_k]<<"\n";
    return true;
}

enum Position
{
    // new scheme
    greater= 1,
    smaller=-greater,
    equal  = 1<<5,
    undef  = 0
};

template <typename T>
Position relativePosition (const T pos, const T ref, const T tol)
{
    if ( pos < ref-tol )
        return smaller;
    if ( pos > ref+tol )
        return greater;
    if ( pos <= ref+tol && pos >= ref-tol)
        return equal;
    return undef;
}

template <typename T>
unsigned findHyperPlaneIntersections (
        const gsBSpline<T>    &curve,
        const gsVector<T>     &normal,
        T                      reference,
        T                      tolerance,
        std::vector<Root<T> > &roots
        )
{
    // argument check
    GISMO_ASSERT( curve.coefDim() == normal.rows(),
                  "cannot intersect a curve in R^"<< curve.coefDim()<<" with a hyperplane in R^"<< normal.rows() );
    if ( !math::isfinite(reference + normal.sum()) )
    {
        gsWarn<<"No intersection reportet between curve and invalid hyperplane."<<std::endl;
    }

    // environment description
    const int deg = curve.degree();
    const T   tol = tolerance;
    const T   ref = reference;

    // internal copy of the curve used for knot insertion
    gsBSpline<T> crv(curve);

    // position status
    Position curP = relativePosition(crv.coef(0).dot(normal), ref, tol) ; // current Position
    if ( curP == undef )
    {
        gsWarn<<__FUNCTION__<<": first coefficient of curve is NAN, exit without results."<<std::endl;
        return 0;
    }
    Position newP = undef;
    Position oldP = undef;
    bool     is_odd = true;

    // knots in the parameter domain
    T oldK = NAN; // last knot inserted
    T newK = NAN; // knot to be insert

    // counters
    int curC = 1; // current coefficient (Control Point)
    int lstC = 0; // coefficient in which the position became the current one
    int intB = 0; // we insert a new knot between the Greville abscissas of intB and curC
    int newK_rep = 0; // number of time we inserted the same knot
    unsigned rootCount = 0; // total number of roots found

    for ( ; (unsigned)curC<crv.coefsSize(); ++curC )
    {
        newP = relativePosition(crv.coef(curC).dot(normal), ref, tol); // current Position
        if ( newP == undef )
        {
            gsWarn<<__FUNCTION__<<": stopping the search for intersection because the curve contains nan coefficients"<<std::endl;
            return rootCount;
        }
        if(curP==newP)
            continue;

        switch (curP)
        {
        case greater:
        case smaller:
            if (curP+newP!=0)
            {
                // from one side of the hyperplane to inside the hyperplane
                // we need to scan till the end of the intersection to determine
                // the type so we only change the state
                oldP=curP;
                curP=newP;
                lstC=curC;
                continue;
            }
            // otherwise follow the code to knot insertion
            is_odd=true;
            intB=curC-1;
            break;
        case equal:
            if (curC-lstC > deg) // interval in hyperplane
            {
                is_odd=(newP!=oldP);
                gsMatrix<T> params;
                params.resize(1,2);
                params(0,0) = crv.knots()[lstC+deg];
                params(0,1) = crv.knots()[curC];
                gsMatrix<T> points=crv.eval(params);
                roots.push_back(Root<T>(is_odd, params(0,0),params(0,1),points.col(0),points.col(1)));
                ++rootCount;
                curP=newP;
                newK_rep=0;
                continue;
            }
            else if ( (curC-lstC) >=  deg+1-(int)crv.knots().multiplicityIndex(curC)) // interpolatory point in hyperplane
            {
                is_odd=(newP!=oldP);
                gsMatrix<T> params;
                params.resize(1,1);
                params(0,0) = crv.knots()[curC];
                gsMatrix<T> points=crv.eval(params);
                roots.push_back(Root<T>(is_odd, params(0,0),points.col(0)));
                ++rootCount;
                curP=newP;
                newK_rep=0;
                continue;
            }
            else if (newP==oldP)
            {
                curP=newP;
                continue;
            }
            else
            {
                is_odd = true;
                intB=math::max(0,lstC-1);
            }
            break;
        default:
            GISMO_ASSERT(false, __FUNCTION__<<"The impossible happened! Check for B..S!!" );
        }
        const T grevB = crv.knots().greville(intB);
        const T grevE = crv.knots().greville(curC);
        const T coeffB = crv.coef(curC).dot(normal)-ref;
        const T coeffE = ref-crv.coef(intB).dot(normal);
        const T denom  = (crv.coef(curC)-crv.coef(intB)).dot(normal);
        newK = coeffB*(grevB/denom)+coeffE*(grevE/denom);
        // in rare but real cases we go out of the interval due to numerical approximation
        newK = math::max(newK, curve.knots().first());
        newK = math::min(newK, curve.knots().last());

        if ( math::abs(newK-oldK) < tol )
        {
            gsMatrix<T> point=crv.eval(gsMatrix<T>::Constant(1,1,newK));
            if (math::abs(point.col(0).dot(normal)-ref)<tol)
            {
                roots.push_back(Root<T>(is_odd, newK,point.col(0)));
                ++rootCount;
                curP=newP;
                oldK=NAN;
                lstC=curC;
                newK_rep=0;
                continue;
            }
            if (newK==oldK) { newK_rep+=1; }
            // if we inserted the same knot degre times and still do not
            // satisfy the requirements then fallback to bisection
            // as this could be caused by approximations in the computation of newK
            if (newK_rep>=deg) newK=(grevB+grevE)/2;
            // if the newK equals the previously added knot anyway then we
            // resolved the position of the zero up to numerical precision
            // in the domain, add the zero and warn that the image is outside the
            // the hyperplane.
            if (newK==oldK)
            {
                gsWarn<< __FUNCTION__<<": intersection does not satisfy required tolerance with knot insertion up to numerical precision\n";
                roots.push_back(Root<T>(is_odd, newK,point.col(0)));
                ++rootCount;
                curP=newP;
                oldK=NAN;
                lstC=curC;
                newK_rep=0;
                continue;
            }
        }
        crv.insertKnot(newK);
        curC = math::max(math::max(curC-deg, lstC-1),0); // incremented by one on end of loop
        curP = relativePosition(crv.coef(curC).dot(normal), ref, tol);
        oldK = newK;
    }
    // intersections at the end of the curve
    if(curP==equal)
    {
        if (curC-lstC > deg) // interval in hyperplane
        {
            gsMatrix<T> params;
            params.resize(1,2);
            params(0,0) = crv.knots()[lstC+deg];
            params(0,1) = crv.knots()[curC];
            gsMatrix<T> points=crv.eval(params);
            roots.push_back(Root<T>(true, params(0,0),params(0,1),points.col(0),points.col(1)));
            ++rootCount;
        }
        else
        {
            gsMatrix<T> params;
            params.resize(1,1);
            params(0,0) = crv.knots()[curC];
            gsMatrix<T> points=crv.eval(params);
            roots.push_back(Root<T>(true, params(0,0),points.col(0)));
            ++rootCount;
        }
    }
    return rootCount;
}

} //namespace gismo
