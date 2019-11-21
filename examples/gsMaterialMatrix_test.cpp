/** @file gsExprIntegral_test.cpp

    @brief Testing integral computation using the expression evaluator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

#include <gsAssembler/gsExprEvaluator.h>

using namespace gismo;

/*
    gsIntegrantZ has as input
        -
        -
        -

*/
template<class T>
class gsIntegrantZ : public gismo::gsFunction<T>
{
    protected:
        const typename gsFunction<T>::Ptr _fun;
        mutable gsMatrix<T> tmp;
        gsMatrix<T> _surfPts;
    public:
        /// Shared pointer for gsIntegrantZ
        typedef memory::shared_ptr< gsIntegrantZ > Ptr;

        /// Unique pointer for gsIntegrantZ
        typedef memory::unique_ptr< gsIntegrantZ > uPtr;

    // copy constructor
    explicit gsIntegrantZ(const gsIntegrantZ &other) : _fun(other._fun), _surfPts(other._surfPts) {}

    gsIntegrantZ(const gsFunction<T> & fun) : _fun(fun.clone()) { }

    GISMO_CLONE_FUNCTION(gsIntegrantZ)

    void setPoint(const gsMatrix<T>& surfPts) { _surfPts = surfPts; }

    gsMatrix<T> point() { return _surfPts; }

    short_t domainDim() const {return 1;}

    short_t targetDim() const {return _fun->targetDim();}

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows()==1,
            "The number of rows for the 1D coordinate is not 1 but " <<u.rows()<<"!"  );
        GISMO_ASSERT(_fun->domainDim()==_surfPts.rows() + 1,
            "The domain dimensions do not match! fun.domainDim() != surfPts.rows() + 1! (" <<_fun->domainDim()<<"!="<<_surfPts.rows() + 1<<" )"  );
        GISMO_ASSERT(_surfPts.cols()==1,
            "Multiple ("<<_surfPts.cols()<<") parametric points given, accepts only 1... " <<"!"     );

        index_t m = _surfPts.rows();
        index_t N = u.cols();

        tmp.resize(m + 1,N);
        tmp.topRows(m)=_surfPts.replicate(1,N);
        tmp.bottomRows(1) = u;

        _fun->eval_into(tmp,result);
    }

    std::ostream &print(std::ostream &os) const
      { os << "gsIntegrantZ ( " << _fun << " )"; return os; };
};

template <class T>
class gsIntegrateZ : public gismo::gsFunction<T>
{
    protected:
        const typename gsFunction<T>::Ptr _fun;//, _t;
        mutable gsMatrix<T> tmp;
        gsMatrix<T> _surfPts;
        T _t;
        int _mom;
    public:
        /// Shared pointer for gsIntegrateZ
        typedef memory::shared_ptr< gsIntegrateZ > Ptr;

        /// Unique pointer for gsIntegrateZ
        typedef memory::unique_ptr< gsIntegrateZ > uPtr;

    // copy constructor
    explicit gsIntegrateZ(const gsIntegrateZ &other)
    : _fun(other._fun), _t(other._t), _surfPts(other._surfPts) {}

    gsIntegrateZ(const gsFunction<T> & fun, T thickness, int moment = 0)
    : _fun(fun.clone()), _t(thickness), _mom(moment)
    {
        // gsConstantFunction _t(thickness, 3);
    }

    // gsIntegrateZ(const gsFunction<T> & fun, gsFunction<T> & thickFun)
    // : _fun(fun.clone()), _t(thickFun.clone()) { }

    GISMO_CLONE_FUNCTION(gsIntegrateZ)

    short_t domainDim() const {return 1;}

    short_t targetDim() const {return _fun->targetDim();}

    void setPoint(const gsMatrix<T>& surfPts) { _surfPts = surfPts; }

    // u are z-coordinates only!
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // // Compute the thickness
        // _tmp.points = u;
        // _t.eval_into(_tmp.values[0], thickMat);

        result.resize(_fun->targetDim(),u.cols());

        // Define integrator for the z direction
        gsExprEvaluator<real_t> ev;

        // Define integration interval
        int k = 1; // interior knots
        int p = 1; // B-spline order

        // Make 1D domain with a basis
        gsKnotVector<> KV(-_t/2.0, _t/2.0, k, p+1);
        gsMultiBasis<> basis;
        basis.addBasis(gsBSplineBasis<>::make(KV));

        // Set integration elements along basis
        ev.setIntegrationElements(basis);

        // Define integrant variables
        typedef gsExprEvaluator<real_t>::variable    variable;
        variable    integrant = ev.getVariable(*_fun, 1);

        gsFunctionExpr<> height("x",1);
        variable z = ev.getVariable(height);
        for (index_t i=0; i!=_fun->targetDim(); ++i)
            for (index_t j = 0; j != u.cols(); ++j)
            {
                //thickness integral for all components i with moment _mom
                ev.integral(pow(z,_mom)*integrant.tr()[i]);
                result(i,j) = ev.value();
            }
    }

    std::ostream &print(std::ostream &os) const
      { os << "gsIntegrateZ ( " << _fun << " )"; return os; };
};

template <class T>
class gsIntegrate : public gismo::gsFunction<T>
{
    protected:
        const typename gsFunction<T>::Ptr _fun, _t;
        mutable gsMatrix<T> tmp, thickMat;
        gsMatrix<T> _surfPts;
        int _mom;
        // T _t;
    public:
        /// Shared pointer for gsIntegrate
        typedef memory::shared_ptr< gsIntegrate > Ptr;

        /// Unique pointer for gsIntegrate
        typedef memory::unique_ptr< gsIntegrate > uPtr;

    // copy constructor
    explicit gsIntegrate(const gsIntegrate &other)
    : _fun(other._fun), _t(other._t), _surfPts(other._surfPts) {}

    // gsIntegrate(const gsFunction<T> & fun, T thickness)
    // : _fun(fun.clone()), _t(thickness)
    // {
    //     gsConstantFunction _t(thickness, 3);
    // }

    gsIntegrate(const gsFunction<T> & fun, gsFunction<T> & thickFun, int moment = 0)
    : _fun(fun.clone()), _t(thickFun.clone()), _mom(moment) { }

    GISMO_CLONE_FUNCTION(gsIntegrate)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return _fun->targetDim();}

    void setMoment(const int moment ) { _mom = moment;}

    // u are xy-coordinates only; domainDim=2
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // Compute the thickness
        _t->eval_into(u, thickMat);

        result.resize(_fun->targetDim(),u.cols());

        // Define integrator for the z direction
        gsExprEvaluator<real_t> ev;

        // Define integrant variables
        typedef gsExprEvaluator<real_t>::variable    variable;
        gsIntegrantZ integrant(*_fun);

        T tHalf;
        gsFunctionExpr<> height("x",1);
        variable z = ev.getVariable(height);
        for (index_t i=0; i!=_fun->targetDim(); ++i)
            for (index_t j = 0; j != u.cols(); ++j)
            {
                // this part is quite sloppy since a multi-basis is created every iteration..
                tHalf = thickMat(0,j)/2.0;

                gsKnotVector<> KV(-tHalf, tHalf, 2, 2);
                gsMultiBasis<> basis;
                basis.addBasis(gsBSplineBasis<>::make(KV));

                // Set integration elements along basis
                ev.setIntegrationElements(basis);

                variable intfun = ev.getVariable(integrant, 1);

                // set new integration point
                integrant.setPoint(u.col(j));

                //thickness integral for all components i with moment _mom
                ev.integral(pow(z,_mom)*intfun.tr()[i]);
                result(i,j) = ev.value();
            }
    }

    std::ostream &print(std::ostream &os) const
      { os << "gsIntegrateZ ( " << _fun << " )"; return os; };
};

template <class T>
class gsMaterialMatrixLinear : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
  //
protected:
    const gsFunctionSet<T> * _mp;
    const gsFunction<T> * _YoungsModulus;
    const gsFunction<T> * _PoissonRatio;
    mutable gsMapData<T> _tmp;
    mutable gsMatrix<real_t,3,3> F0;
    mutable gsMatrix<T> Emat,Nmat;
    mutable real_t lambda, mu, E, nu, C_constant;

public:
    /// Shared pointer for gsMaterialMatrixLinear
    typedef memory::shared_ptr< gsMaterialMatrixLinear > Ptr;

    /// Unique pointer for gsMaterialMatrixLinear
    typedef memory::unique_ptr< gsMaterialMatrixLinear > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsMaterialMatrixLinear(const gsFunctionSet<T> & mp, const gsFunction<T> & YoungsModulus,
                   const gsFunction<T> & PoissonRatio) :
    _mp(&mp), _YoungsModulus(&YoungsModulus), _PoissonRatio(&PoissonRatio), _mm_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    ~gsMaterialMatrixLinear() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrixLinear)

    short_t domainDim() const {return 3;}

    short_t targetDim() const {return 9;}

    mutable gsMaterialMatrixLinear<T> * _mm_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _mm_piece;
        _mm_piece = new gsMaterialMatrixLinear(_mp->piece(k), *_YoungsModulus, *_PoissonRatio);
        return *_mm_piece;
    }

    //class .. matMatrix_z
    // should contain eval_into(thickness variable)

    // Input is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // NOTE 1: if the input \a u is considered to be in physical coordinates
        // then we first need to invert the points to parameter space
        // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
        // otherwise we just use the input paramteric points

        _tmp.points = u.topRows(2);

        static_cast<const gsFunction<T>&>( _mp->piece(0) ).computeMap(_tmp);
        // _mp->piece(0).computeMap(_tmp);

        // NOTE 2: in the case that parametric value is needed it suffices
        // to evaluate Youngs modulus and Poisson's ratio at
        // \a u instead of _tmp.values[0].
        _YoungsModulus->eval_into(_tmp.values[0], Emat);
        _PoissonRatio->eval_into(_tmp.values[0], Nmat);

        result.resize( targetDim() , u.cols() );
        for( index_t k=0; k< u.cols(); ++k )
        {
            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(k,3,3);

            F0.leftCols(2) = _tmp.jacobian(k);
            F0.col(2)      = _tmp.normal(k).normalized();
            F0 = F0.inverse();
            F0 = F0 * F0.transpose(); //3x3

            // Evaluate material properties on the quadrature point
            E = Emat(0,k);
            nu = Nmat(0,k);
            lambda = E * nu / ( (1. + nu)*(1.-2.*nu)) ;
            mu     = E / (2.*(1. + nu)) ;

            C_constant = 4*lambda*mu/(lambda+2*mu);

            C(0,0) = C_constant*F0(0,0)*F0(0,0) + 2*mu*(2*F0(0,0)*F0(0,0));
            C(1,1) = C_constant*F0(1,1)*F0(1,1) + 2*mu*(2*F0(1,1)*F0(1,1));
            C(2,2) = C_constant*F0(0,1)*F0(0,1) + 2*mu*(F0(0,0)*F0(1,1) + F0(0,1)*F0(0,1));
            C(1,0) =
            C(0,1) = C_constant*F0(0,0)*F0(1,1) + 2*mu*(2*F0(0,1)*F0(0,1));
            C(2,0) =
            C(0,2) = C_constant*F0(0,0)*F0(0,1) + 2*mu*(2*F0(0,0)*F0(0,1));
            C(2,1) = C(1,2) = C_constant*F0(0,1)*F0(1,1) + 2*mu*(2*F0(0,1)*F0(1,1));

            //gsDebugVar(C);
        }
    }

    // std::ostream &print(std::ostream &os) const
    //   { os << "gsMaterialMatrix "; return os; };
};

// To Do:
// * struct for material model
// Input is parametric coordinates of the surface \a mp
template <class T, int mat>
class gsMaterialMatrixIncompressible : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
  //
protected:
    const gsFunctionSet<T> * _mp;
    const gsFunction<T> * _par1;
    // const gsFunction<T> * _par2;
    const gsFunctionSet<T> * _mp_def;
    mutable gsMapData<T> _tmp;
    mutable gsMapData<T> _tmp_def;
    mutable gsMatrix<T> jacGdef, jacGori, a_ori, a_def, b_ori, b_def, g_def, g_ori;
    mutable gsVector<T> n_def, n_ori;
    mutable gsMatrix<T> par1mat,par2mat;
    mutable real_t mu, J0;

public:
    /// Shared pointer for gsMaterialMatrixIncompressible
    typedef memory::shared_ptr< gsMaterialMatrixIncompressible > Ptr;

    /// Unique pointer for gsMaterialMatrixIncompressible
    typedef memory::unique_ptr< gsMaterialMatrixIncompressible > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsMaterialMatrixIncompressible( const gsFunctionSet<T> & mp,
                                    const gsFunction<T> & par1,
                                    // const gsFunction<T> & par2,
                                    const gsFunctionSet<T> & mp_def) : // deformed multipatch
    _mp(&mp),
    _par1(&par1),
    // _par2(&par2),
    _mp_def(&mp_def),
    _mm_piece(nullptr)
    {
        _tmp.flags     = NEED_VALUE | NEED_JACOBIAN | NEED_NORMAL |  NEED_2ND_DER;
        _tmp_def.flags = NEED_VALUE | NEED_JACOBIAN | NEED_NORMAL |  NEED_2ND_DER;
    }

    ~gsMaterialMatrixIncompressible() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrixIncompressible)

    short_t domainDim() const {return 3;}

    short_t targetDim() const
    {
        switch (mat)
        {
            case 0:
                return 3;
                break;
            case 1:
                return 9;
                break;
        }
    }

    mutable gsMaterialMatrixIncompressible<T,mat> * _mm_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _mm_piece;
        // _mm_piece = new gsMaterialMatrixIncompressible(_mp->piece(k), *_par1, *_par2, _mp_def->piece(k) );
        _mm_piece = new gsMaterialMatrixIncompressible(_mp->piece(k), *_par1, _mp_def->piece(k) );
        return *_mm_piece;
    }

    // Input is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // NOTE 1: if the input \a u is considered to be in physical coordinates
        // then we first need to invert the points to parameter space
        // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
        // otherwise we just use the input paramteric points
        _tmp.points     = u.topRows(2);
        _tmp_def.points = u.topRows(2);

        static_cast<const gsFunction<T>&>( _mp->piece(0)     ).computeMap(_tmp);
        static_cast<const gsFunction<T>&>( _mp_def->piece(0) ).computeMap(_tmp_def);

        // NOTE 2: in the case that parametric value is needed it suffices
        // to evaluate Youngs modulus and Poisson's ratio at
        // \a u instead of _tmp.values[0].
        _par1->eval_into(_tmp.values[0], par1mat);
        // _par2->eval_into(_tmp.values[0], par2mat);

        result.resize( targetDim() , u.cols() );
        for( index_t k=0; k< u.cols(); ++k )
        {
            // Construct metric tensor a = [dcd1*dcd1, dcd1*dcd2; dcd2*dcd1, dcd2*dcd2]
            jacGdef = _tmp_def.jacobian(k);
            a_def   = jacGdef.transpose() * jacGdef;
            jacGori = _tmp.jacobian(k);
            a_ori   = jacGori.transpose() * jacGori;

            // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
            gsAsConstMatrix<T,3,3> deriv2def( _tmp_def.deriv2(0).data(),3,3 ); // gives [d11 c1, d11c2, d11c3; d22c1, d22c1, d22c3; d12c1, d12c2, d12c3]
            gsAsConstMatrix<T,3,3> deriv2ori(     _tmp.deriv2(0).data(),3,3 );
            n_def = _tmp_def.normal(k).normalized();
            n_ori = _tmp.normal(k).normalized();
            b_def.resize(2,2);
            b_ori.resize(2,2);

            b_def(0,0) = deriv2def.row(0).dot(n_def);
            b_def(1,1) = deriv2def.row(1).dot(n_def);
            b_def(0,1) = b_def(1,0) = deriv2def.row(2).dot(n_def);

            b_ori(0,0) = deriv2ori.row(0).dot(n_ori);
            b_ori(1,1) = deriv2ori.row(1).dot(n_ori);
            b_ori(0,1) = b_ori(1,0) = deriv2ori.row(2).dot(n_ori);

            // Construct basis of coordinate system g = [a_ij - 2*theta3*b_ij]
            g_def = a_def - 2 * u(2,k) * b_def;
            g_ori = a_ori - 2 * u(2,k) * b_ori;

            // Evaluate material properties on the quadrature point
            mu = par1mat(0,k);
            J0 = math::sqrt( g_def.determinant() / g_ori.determinant() );
            J0 = math::pow( J0, -2 );

            switch (mat)
            {
                case 0:
                    // Sab(0,0) = mu * (jacGori(0,0) - math::pow(J0,-2.) * jacGdef(0,0) );
                    // Sab(0,1) =
                    // Sab(1,0) = mu * (jacGori(1,0) - math::pow(J0,-2.) * jacGdef(1,0) ); // CHECK SYMMETRIES
                    // Sab(1,1) = mu * (jacGori(1,1) - math::pow(J0,-2.) * jacGdef(1,1) );

                    result(0,k) = mu * (jacGori(0,0) - math::pow(J0,-2.) * jacGdef(0,0) );
                    result(1,k) = mu * (jacGori(1,1) - math::pow(J0,-2.) * jacGdef(1,1) );
                    result(2,k) = mu * (jacGori(1,0) - math::pow(J0,-2.) * jacGdef(1,0) ); // CHECK SYMMETRIES
                    break;

                case 1:
                    /*
                        C =     C1111,  C1122,  C1112
                                symm,   C2222,  C2212
                                symm,   symm,   C1212
                    */
                    gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(k,3,3);
                    C(0,0) = 2*g_ori(0,0)*g_ori(0,0) + g_ori(0,0)*g_ori(0,0) + g_ori(0,0)*g_ori(0,0); // C1111
                    C(1,0) =
                    C(0,1) = 2*g_ori(0,0)*g_ori(1,1) + g_ori(0,1)*g_ori(0,1) + g_ori(0,1)*g_ori(0,1); // C1122
                    C(2,0) =
                    C(0,2) = 2*g_ori(0,0)*g_ori(0,1) + g_ori(0,0)*g_ori(0,1) + g_ori(0,1)*g_ori(0,0); // C1112
                    C(1,1) = 2*g_ori(1,1)*g_ori(1,1) + g_ori(1,1)*g_ori(1,1) + g_ori(1,1)*g_ori(1,1); // C2222
                    C(2,1) =
                    C(1,2) = 2*g_ori(1,1)*g_ori(0,1) + g_ori(1,0)*g_ori(1,1) + g_ori(1,1)*g_ori(1,0); // C2212
                    C(2,2) = 2*g_ori(0,1)*g_ori(0,1) + g_ori(0,0)*g_ori(1,1) + g_ori(0,1)*g_ori(1,0); // C1212

                    C *= mu * J0;
                    break;
            }
        }
    }

};


// To Do:
// * struct for material model
// Input is parametric coordinates of the surface \a mp
template <class T, int mat>
class gsMaterialMatrixCompressible : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
  //
protected:
    const gsFunctionSet<T> * _mp;
    const gsFunction<T> * _par1;
    const gsFunction<T> * _par2;
    const gsFunctionSet<T> * _mp_def;
    mutable gsMapData<T> _tmp;
    mutable gsMapData<T> _tmp_def;
    mutable gsMatrix<T> jacGdef, jacGori, a_ori, a_def, b_ori, b_def, g_def, g_ori;
    mutable gsVector<T> n_def, n_ori;
    mutable gsMatrix<T> par1mat,par2mat;
    mutable real_t mu, K, J0, J;

public:
    /// Shared pointer for gsMaterialMatrixCompressible
    typedef memory::shared_ptr< gsMaterialMatrixCompressible > Ptr;

    /// Unique pointer for gsMaterialMatrixCompressible
    typedef memory::unique_ptr< gsMaterialMatrixCompressible > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsMaterialMatrixCompressible( const gsFunctionSet<T> & mp,
                                    const gsFunction<T> & par1,
                                    const gsFunction<T> & par2,
                                    const gsFunctionSet<T> & mp_def) : // deformed multipatch
    _mp(&mp),
    _par1(&par1),
    _par2(&par2),
    _mp_def(&mp_def),
    _mm_piece(nullptr)
    {
        _tmp.flags     = NEED_VALUE | NEED_JACOBIAN | NEED_NORMAL |  NEED_2ND_DER;
        _tmp_def.flags = NEED_VALUE | NEED_JACOBIAN | NEED_NORMAL |  NEED_2ND_DER;
    }

    ~gsMaterialMatrixCompressible() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrixCompressible)

    short_t domainDim() const {return 3;}

    short_t targetDim() const
    {
        switch (mat)
        {
            case 0:
                return 3;
                break;
            case 1:
                return 9;
                break;
        }
    }

    mutable gsMaterialMatrixCompressible<T,mat> * _mm_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _mm_piece;
        _mm_piece = new gsMaterialMatrixCompressible(_mp->piece(k), *_par1, *_par2, _mp_def->piece(k) );
        // _mm_piece = new gsMaterialMatrixCompressible(_mp->piece(k), *_par1, _mp_def.piece(k) );
        return *_mm_piece;
    }



//class .. matMatrix_z
// should contain eval_into(thickness variable)


    // void computeThickness()
    // {

    // }

    // Input is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // NOTE 1: if the input \a u is considered to be in physical coordinates
        // then we first need to invert the points to parameter space
        // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
        // otherwise we just use the input paramteric points
        _tmp.points     = u.topRows(2);
        _tmp_def.points = u.topRows(2);

        static_cast<const gsFunction<T>&>( _mp->piece(0)     ).computeMap(_tmp);
        static_cast<const gsFunction<T>&>( _mp_def->piece(0) ).computeMap(_tmp_def);

        // NOTE 2: in the case that parametric value is needed it suffices
        // to evaluate Youngs modulus and Poisson's ratio at
        // \a u instead of _tmp.values[0].
        _par1->eval_into(_tmp.values[0], par1mat);
        _par2->eval_into(_tmp.values[0], par2mat);

        result.resize( targetDim() , u.cols() );
        for( index_t k=0; k< u.cols(); ++k )
        {
            // Material parameters
            mu = par1mat(0,k);
            K = par1mat(0,k);

            // Define objects
            gsMatrix<T,3,3> c, cinv;
            T S33, C3333, dc33, traceC;

            // Construct metric tensor a = [dcd1*dcd1, dcd1*dcd2; dcd2*dcd1, dcd2*dcd2]
            jacGdef = _tmp_def.jacobian(k);
            a_def   = jacGdef.transpose() * jacGdef;
            jacGori = _tmp.jacobian(k);
            a_ori   = jacGori.transpose() * jacGori;

            // Construct metric tensor b = [d11c*n, d12c*n ; d21c*n, d22c*n]
            gsAsConstMatrix<T,3,3> deriv2def( _tmp_def.deriv2(0).data(),3,3 ); // gives [d11 c1, d11c2, d11c3; d22c1, d22c1, d22c3; d12c1, d12c2, d12c3]
            gsAsConstMatrix<T,3,3> deriv2ori(     _tmp.deriv2(0).data(),3,3 );
            n_def = _tmp_def.normal(k).normalized();
            n_ori = _tmp.normal(k).normalized();
            b_def.resize(2,2);
            b_ori.resize(2,2);

            b_def(0,0) = deriv2def.row(0).dot(n_def);
            b_def(1,1) = deriv2def.row(1).dot(n_def);
            b_def(0,1) = b_def(1,0) = deriv2def.row(2).dot(n_def);

            b_ori(0,0) = deriv2ori.row(0).dot(n_ori);
            b_ori(1,1) = deriv2ori.row(1).dot(n_ori);
            b_ori(0,1) = b_ori(1,0) = deriv2ori.row(2).dot(n_ori);

            // Construct basis of coordinate system g = [a_ij - 2*theta3*b_ij]
            g_def = g_ori = gsMatrix<T>::Zero(3,3);
            g_def.block(0,0,2,2) = a_def - 2 * u(2,k) * b_def;
            g_ori.block(0,0,2,2) = a_ori - 2 * u(2,k) * b_ori;
            g_def(2,2) = g_ori(2,2) = 1.0;

            // Initialize c
            c.setZero();
            c.block(0,0,2,2) = g_def.block(0,0,2,2);
            c(2,2) = 1.0; // c33
            cinv = c.inverse();
            // note: can also just do c = jacGdef because the normal has length one and hence c(2,2) is 1. CHECK!

            J0 = math::sqrt( g_def.determinant() / g_ori.determinant() );
            J = J0 * math::sqrt( c(2,2) );

            index_t itmax = 20;
            T tol = 1e-6;
            S33 = 0.0;
            C3333 = 1.0;

            // Define lambda function for C
            std::function<T (index_t i, index_t j, index_t k, index_t l)> Cijkl;
            Cijkl = [=](index_t i, index_t j, index_t k, index_t l)
            {
                T res = 1.0 / 9.0 * mu * math::pow( J , -2.0/3.0 ) * ( traceC * ( 2*cinv(i,j)*cinv(k,l) + 3*cinv(i,k)*cinv(j,l) + 3*cinv(i,l)*cinv(j,k) )
                                - 6*g_ori(i,j)*cinv(k,l) + cinv(i,j)*g_ori(k,l) ) + K * ( J*J*cinv(i,j)*cinv(k,l) - 0.5*(J*J-1)*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) ) );
                return res;
            };

            for (index_t it = 0; it < itmax; it++)
            {
                dc33 = -2. * S33 / C3333;
                c(2,2) += dc33;
                cinv(2,2) = 1.0/c(2,2);

                traceC = c.trace();
                J = J0 * math::sqrt( c(2,2) );

                S33     = mu * math::pow( J , -2.0/3.0 ) * ( g_ori(2,2) - 1.0/3.0 * traceC * cinv(2,2) ) + 0.5 * K * ( J*J - 1 ) * cinv(2,2);
                C3333   = Cijkl(2,2,2,2);

                if (S33 < tol)
                {
                    // gsInfo<<"Converged in "<<it<<" iterations, S33 = "<<S33<<" and tolerance = "<<tol<<"\n";
                    switch (mat)
                    {
                        case 0:
                            result(0,k) = mu * math::pow( J , -2.0/3.0 ) * ( g_ori(0,0) - 1.0/3.0 * traceC * cinv(0,0) ) + 0.5 * K * ( J*J - 1 ) * cinv(0,0); // S11
                            result(1,k) = mu * math::pow( J , -2.0/3.0 ) * ( g_ori(1,1) - 1.0/3.0 * traceC * cinv(1,1) ) + 0.5 * K * ( J*J - 1 ) * cinv(1,1); // S22
                            result(2,k) = mu * math::pow( J , -2.0/3.0 ) * ( g_ori(0,1) - 1.0/3.0 * traceC * cinv(0,1) ) + 0.5 * K * ( J*J - 1 ) * cinv(0,1); // S12
                            break;
                        case 1:
                            /*
                                C =     C1111,  C1122,  C1112
                                        symm,   C2222,  C2212
                                        symm,   symm,   C1212
                                Here, Cabcd = Cijkl - Cab33*C33cd / C3333;
                                a,b,c,d = 1,2; i,j,k,l = 1...3;
                            */
                            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(k,3,3);
                            C(0,0) = Cijkl(0,0,0,0) - ( Cijkl(0,0,2,2) * Cijkl(2,2,0,0) ) / (Cijkl(2,2,2,2)); // C1111
                            C(0,1) =
                            C(1,0) = Cijkl(0,0,1,1) - ( Cijkl(0,0,2,2) * Cijkl(2,2,1,1) ) / (Cijkl(2,2,2,2)); // C1122
                            C(0,2) =
                            C(2,0) = Cijkl(0,0,0,1) - ( Cijkl(0,0,2,2) * Cijkl(2,2,0,1) ) / (Cijkl(2,2,2,2)); // C1112
                            C(1,1) = Cijkl(1,1,1,1) - ( Cijkl(1,1,2,2) * Cijkl(2,2,1,1) ) / (Cijkl(2,2,2,2)); // C2222
                            C(1,2) =
                            C(2,1) = Cijkl(1,1,0,1) - ( Cijkl(1,1,2,2) * Cijkl(2,2,0,1) ) / (Cijkl(2,2,2,2)); // C2212
                            C(2,2) = Cijkl(0,1,0,1) - ( Cijkl(0,1,2,2) * Cijkl(2,2,0,1) ) / (Cijkl(2,2,2,2)); // C1212
                            break;
                    }
                    break;
                }
                else if (it == itmax - 1)
                {
                    gsInfo<<"Error: Method did not converge, S33 = "<<S33<<" and tolerance = "<<tol<<"\n";
                    // std::terminate();
                }
            }
        }
    }
};

/*
    Todo:
        * Improve for mu, E, phi as gsFunction instead of reals
*/
template <class T, int mat>
class gsMaterialMatrixD : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
  // NOTE: This material matrix is in local Cartesian coordinates and should be transformed!
  // NOTE: To make this efficient, we can output all matrices stacked and use expressions to pick the 1st, second and third
protected:
    // const gsFunctionSet<T> * _mp;
    const std::vector<std::pair<T,T>> _YoungsModuli;
    const std::vector<T> _ShearModuli;
    const std::vector<std::pair<T,T>> _PoissonRatios;
    const std::vector<T> _thickness;
    const std::vector<T> _phi;
    mutable gsMapData<T> _tmp;
    mutable real_t E1, E2, G12, nu12, nu21, t, t_tot, t_temp, z, z_mid, phi;
    mutable gsMatrix<T> Tmat, Dmat;

public:
    /// Shared pointer for gsMaterialMatrixD
    typedef memory::shared_ptr< gsMaterialMatrixD > Ptr;

    /// Unique pointer for gsMaterialMatrixD
    typedef memory::unique_ptr< gsMaterialMatrixD > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsMaterialMatrixD(  //const gsFunctionSet<T> & mp,
                        const std::vector<std::pair<T,T>> & YoungsModuli,
                        const std::vector<T> & ShearModuli,
                        const std::vector<std::pair<T,T>> & PoissonRatios,
                        const std::vector<T> thickness,
                        const std::vector<T> phi) :
    // _mp(&mp),
    _YoungsModuli(YoungsModuli),
    _ShearModuli(ShearModuli),
    _PoissonRatios(PoissonRatios),
    _thickness(thickness),
    _phi(phi)//,
    // _mm_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    // ~gsMaterialMatrixD() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrixD)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 9;}

    // mutable gsMaterialMatrixD<T> * _mm_piece; // todo: improve the way pieces are accessed

    // const gsFunction<T> & piece(const index_t k) const
    // {
    //     delete _mm_piece;
    //     _mm_piece = new gsMaterialMatrixD(_mp->piece(k), _YoungsModuli, _ShearModuli, _PoissonRatios, _thickness, _phi);
    //     return *_mm_piece;
    // }

    // Input is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // static_cast<const gsFunction<T>*>(_mp)->computeMap(_tmp);
        GISMO_ASSERT(_YoungsModuli.size()==_PoissonRatios.size(),"Size of vectors of Youngs Moduli and Poisson Ratios is not equal: " << _YoungsModuli.size()<<" & "<<_PoissonRatios.size());
        GISMO_ASSERT(_YoungsModuli.size()==_ShearModuli.size(),"Size of vectors of Youngs Moduli and Shear Moduli is not equal: " << _YoungsModuli.size()<<" & "<<_ShearModuli.size());
        GISMO_ASSERT(_thickness.size()==_phi.size(),"Size of vectors of thickness and angles is not equal: " << _thickness.size()<<" & "<<_phi.size());
        GISMO_ASSERT(_YoungsModuli.size()==_thickness.size(),"Size of vectors of material properties and laminate properties is not equal: " << _YoungsModuli.size()<<" & "<<_thickness.size());
        GISMO_ASSERT(_YoungsModuli.size()!=0,"No laminates defined");

        // Compute total thickness (sum of entries)
        t_tot = std::accumulate(_thickness.begin(), _thickness.end(), 0.0);

        // compute mid-plane height of total plate
        z_mid = t_tot / 2.0;

        // now we use t_temp to add the thickness of all plies iteratively
        t_temp = 0.0;

        // Initialize material matrix and result
        Dmat.resize(3,3);
        result.resize( targetDim(), 1 );

        // Initialize transformation matrix
        Tmat.resize(3,3);


        for (size_t i = 0; i != _phi.size(); ++i) // loop over laminates
        {
            // Compute all quantities
            E1 = _YoungsModuli[i].first;
            E2 = _YoungsModuli[i].second;
            G12 = _ShearModuli[i];
            nu12 = _PoissonRatios[i].first;
            nu21 = _PoissonRatios[i].second;
            t = _thickness[i];
            phi = _phi[i];

            GISMO_ASSERT(nu21*E1 == nu12*E2, "No symmetry in material properties for ply "<<i<<". nu12*E2!=nu21*E1:\n"<<
                    "\tnu12 = "<<nu12<<"\t E2 = "<<E2<<"\t nu12*E2 = "<<nu12*E2<<"\n"
                  <<"\tnu21 = "<<nu21<<"\t E1 = "<<E1<<"\t nu21*E1 = "<<nu21*E1);

            // Fill material matrix
            Dmat(0,0) = E1 / (1-nu12*nu21);
            Dmat(1,1) = E2 / (1-nu12*nu21);;
            Dmat(2,2) = G12;
            Dmat(0,1) = nu21*E1 / (1-nu12*nu21);
            Dmat(1,0) = nu12*E2 / (1-nu12*nu21);
            Dmat(2,0) = Dmat(0,2) = Dmat(2,1) = Dmat(1,2) = 0.0;

            // Make transformation matrix
            Tmat(0,0) = Tmat(1,1) = math::pow(math::cos(phi),2);
            Tmat(0,1) = Tmat(1,0) = math::pow(math::sin(phi),2);
            Tmat(2,0) = Tmat(0,2) = Tmat(2,1) = Tmat(1,2) = math::sin(phi) * math::cos(phi);
            Tmat(2,0) *= -2.0;
            Tmat(2,1) *= 2.0;
            Tmat(1,2) *= -1.0;
            Tmat(2,2) = math::pow(math::cos(phi),2) - math::pow(math::sin(phi),2);

            // Compute laminate stiffness matrix
            Dmat = Tmat.transpose() * Dmat * Tmat;

            z = math::abs(z_mid - (t/2.0 + t_temp) ); // distance from mid-plane of plate

            // Make matrices A, B and C
            // [NOTE: HOW TO DO THIS NICELY??]
            // result.reshape(3,3) += Dmat * t; // A
            // result.reshape(3,3) += Dmat * t*z; // B
            // result.reshape(3,3) += Dmat * ( t*z*z + t*t*t/12.0 ); // D
            switch (mat)
            {
                case 0:
                    result.reshape(3,3) += Dmat * t; // A
                    break;
                case 1:
                    result.reshape(3,3) += Dmat * t; // A
                    break;
                case 2:
                    result.reshape(3,3) += Dmat * t; // A
                    break;
            }

            t_temp += t;
        }

        GISMO_ASSERT(t_tot==t_temp,"Total thickness after loop is wrong. t_temp = "<<t_temp<<" and sum(thickness) = "<<t_tot);

        // Replicate for all points since the quantities are equal over the whole domain
        result.replicate(1, u.cols());
    }

    // piece(k) --> for patch k

}; //! [Include namespace]


int main(int argc, char *argv[])
{
    gsCmdLine cmd("Testing expression evaluator.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;

    mp.addPatch( gsNurbsCreator<>::BSplineSquare(4) ); // degree
    mp.addAutoBoundaries();
    mp.embed(3);

    mp.uniformRefine(1);
    // mp.degreeElevate(2);
    gsMultiBasis<> b(mp);

    // Initiate the expression evaluator
    gsExprEvaluator<real_t> ev;

    // Set the parameter mesh as the integration mesh
    ev.setIntegrationElements(b);


    /*
        test gsIntegrantZ function
    */

    gsVector<> pt1D(1); pt1D.setConstant(0.25);
    gsVector<> pt2D(2); pt2D.setConstant(0.25);
    gsVector<> pt3D(3); pt3D.setConstant(0.25);

    // gsMatrix<> points(2,11);
    // points<<0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
    //         0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0;

    gsMatrix<> points(2,3);
    points<<0.0,0.5,1.0,
            0.0,0.5,1.0;

    gsFunctionExpr<> fun("1*x","2*y","x*y*z^2",3);


    gsMatrix<> result;
    fun.eval_into(pt3D,result);
    gsInfo<<"Evaluation of the following function on point (x,y) = ("<<pt3D.at(0)<<","<<pt3D.at(1)<<") and z coordinate "<<pt3D.at(2)<<"\n";
    gsInfo<<fun<<"\n";
    gsInfo<<"result = "<<result.transpose()<<"\n";

    gsInfo<<"Evaluation of the following function on point (x,y) = ("<<pt2D.at(0)<<","<<pt2D.at(1)<<") and z coordinate "<<pt1D.at(0)<<"\n";
    gsInfo<<fun<<"\n";
    gsIntegrantZ fun2(fun);
    fun2.setPoint(pt2D); // if changes to be applied
    fun2.eval_into(pt1D,result);
    gsInfo<<"result = "<<result.transpose()<<"\n";

    pt2D.setConstant(0.1);
    gsInfo<<"Evaluation of the following function on point (x,y) = ("<<pt2D.at(0)<<","<<pt2D.at(1)<<") and z coordinate "<<pt1D.at(0)<<"\n";
    gsInfo<<fun<<"\n";
    fun2.setPoint(pt2D); // if changes to be applied
    fun2.eval_into(pt1D,result);
    gsInfo<<"result = "<<result.transpose()<<"\n";

    /*
        test gsIntegrateZ function
    */
    gsFunctionExpr<> fun3("1","x","x^2",1);

    real_t bound = 1.0;
    gsConstantFunction<> thickFun(bound,2);


    gsIntegrateZ<real_t> integrator(fun3,bound);
    integrator.eval_into(pt1D,result);

    gsInfo<<"Integration of the following function from "<<-bound/2.0<<" to "<<bound/2.0<<": \n";
    gsInfo<<fun3<<"\n";
    gsInfo<<"Result: "<<result.transpose()<<"\n";

    gsFunctionExpr<> fun4("1",3);
    gsInfo<<"Or by directly computing the moments; integral from "<<-bound/2.0<<" to "<<bound/2.0<<": \n";
    // We use the point pt2D since gsIntegrate requires a function fun: R^3->R^n
    gsIntegrate integrator4(fun4,thickFun,0);
    integrator4.setMoment(0);
    integrator4.eval_into(pt2D,result);

    gsInfo<<"Result: "<<result<<"\t";
    integrator4.setMoment(1);
    integrator4.eval_into(pt2D,result);
    gsInfo<<result<<"\t";

    integrator4.setMoment(2);
    integrator4.eval_into(pt2D,result);
    gsInfo<<result<<"\n";


    /*
        test gsIntegrateZ function
    */
    bound = 1.0;
    gsIntegrateZ integrator2(fun2,bound);
    integrator2.eval_into(pt1D,result);

    gsInfo<<"Integration of the third component of the following function from "<<-bound/2.0<<" to "<<bound/2.0<<" on point (x,y) = ("<<pt2D.at(0)<<","<<pt2D.at(1)<<") \n";
    gsInfo<<fun<<"\n";
    gsInfo<<"Result: "<<result.transpose()<<"\n";


    /*
        test gsIntegrate function
    */
    gsIntegrate integrate(fun,thickFun,0);
    integrate.eval_into(points,result);

    gsInfo<<"Integration of the third component of the following function from "<<-bound/2.0<<" to "<<bound/2.0<<"\n";
    gsInfo<<fun<<"\n";
    gsInfo<<"on points (x,y) = \n";
    gsInfo<<points.transpose()<<"\n";
    gsInfo<<"Result: \n"<<result.transpose()<<"\n";


    /*
        Integrate now a material matrix point by point
        NOTE: does not work
    */
    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsMaterialMatrixLinear materialMat(mp, E, nu);
    gsIntegrate integrateMM(materialMat,thickFun,0);

    integrateMM.eval_into(points,result);
    gsInfo<<"Result: \n"<<result<<"\n"; //.reshape(3,3)


    real_t muPar = 1.5*1e6;
    gsFunctionExpr<> mu(std::to_string(muPar),3);
    gsMaterialMatrixIncompressible<real_t,0> materialMatI(mp, mu, mp);
    gsIntegrate integrateMMI(materialMatI,thickFun,0);

    integrateMMI.eval_into(points,result);
    gsInfo<<"Result: \n"<<result<<"\n"; //.reshape(3,3)

    real_t nuPar = 0.45;
    real_t KPar = 2*muPar*(1+nuPar)/(3-6*nuPar);
    gsFunctionExpr<> K(std::to_string(KPar),3);
    gsMaterialMatrixCompressible<real_t,0> materialMatC(mp, mu, K, mp);
    gsIntegrate integrateMMC(materialMatC,thickFun,0);

    integrateMMC.eval_into(points,result);
    gsInfo<<"Result: \n"<<result<<"\n"; //.reshape(3,3)



    // materialMat.eval_into(pt2D, result);
    // materialMat.eval_into(pt3D, result);
    // gsInfo<<"Result: \n"<<result<<"\n"; //.reshape(3,3)




    /*
        make composite material matrix
        NOTE: does not work
    */
    std::vector<std::pair<real_t,real_t>> Emod;
    std::vector<std::pair<real_t,real_t>> Nu;
    std::vector<real_t> G;
    std::vector<real_t> t;
    std::vector<real_t> phi;

    real_t pi = math::atan(1)*4;

    Emod.push_back( std::make_pair(300.0,200.0) );
    Nu.push_back( std::make_pair(0.3,0.2) );
    G.push_back( 100.0 );
    t.push_back( 0.100 );
    phi.push_back( pi/2.0);



    gsInfo<<math::cos(pi/2.0)<<"\n";
    gsInfo<<math::sin(pi/2.0)<<"\n";

    gsMaterialMatrixD<real_t,0> Dmat(Emod,G,Nu,t,phi);
    Dmat.eval_into(pt2D, result);

    gsInfo<<"Result: \n"<<result.reshape(3,3)<<"\n";



    return EXIT_SUCCESS;
}
