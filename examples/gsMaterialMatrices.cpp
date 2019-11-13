// To Do:
// * struct for material model
// Input is parametric coordinates of the surface \a mp
template <class T>
class gsMaterialMatrix : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
  //
protected:
    const gsFunction<T> _mm;
    mutable gsMapData<T> _tmp;
    mutable gsMatrix<T> res, thickMat;
    mutable real_t lambda, mu, E, nu, C_constant;

public:
    /// Shared pointer for gsMaterialMatrix
    typedef memory::shared_ptr< gsMaterialMatrix > Ptr;

    /// Unique pointer for gsMaterialMatrix
    typedef memory::unique_ptr< gsMaterialMatrix > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsMaterialMatrix(   const gsMaterialMatrix & mm,
                        const gsFunction<T> & thickness) :
    _mm(&mm), _thickness(&thickness), _mm_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    ~gsMaterialMatrix() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrix)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 9;}

    mutable gsMaterialMatrix<T> * _mm_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _mm_piece;
        _mm_piece = new gsMaterialMatrix(_mp->piece(k), *_YoungsModulus, *_PoissonRatio);
        return *_mm_piece;
    }

    // Input u is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // Compute the thickness
        _tmp.points = u;
        _thickness.eval_into(_tmp.values[0], thickMat);

        // Define integrator for the z direction
        gsExprEvaluator<real_t> ev;

        // Define integration interval
        int k = 1; // interior knots
        int p = 1; // B-spline order

        // Make 1D domain with a basis
        gsKnotVector<> KV(-thickness/2.0, thickness/2.0, k, p+1);
        gsBasis<>::uPtr tBasis = gsBSplineBasis<>::make(KV);

        // Set integration elements along basis
        ev.setIntegrationElements(tBasis);

        // Define integrant variables
        typedef gsExprEvaluator<real_t>::geometryMap geometryMap;
        typedef gsExprEvaluator<real_t>::variable    variable;
        geometryMap G = ev.getMap(mp);
        gsMaterialMatrix_z materialMat(_mm, u); // make object for material matrix that transforms 3D input to 1D input given a parametric point u
        variable    materialMat = ev.getVariable(materialMat, 1); // material matrix

        //thickness integral
        ev.integral(materialMat);
    }
};

/*
    gsIntegrateThickness has as input
        - gsMaterialMatrixPt: the material matrix point-wise (input: u \in Omega^3, where Omega is parametric domain)
        - Parametric surface coordinates (u \in Omega^2)
        - A vertical coordinate z
*/
template<class T>
class gsIntegrateThickness : public gismo::gsFunction<T>
{
    protected:

    public:
        /// Shared pointer for gsIntegrateThickness
        typedef memory::shared_ptr< gsIntegrateThickness > Ptr;

        /// Unique pointer for gsIntegrateThickness
        typedef memory::unique_ptr< gsIntegrateThickness > uPtr;

    gsIntegrateThickness(const gsMaterialMatrixPt & mm, const gsMatrix<T>& u) : _mm(mm), _u(u) { }
    GISMO_CLONE_FUNCTION(gsIntegrateThickness)

    int domainDim() const {return 1;}

    int targetDim() const {return 9;}

    void eval_into(const gsMatrix<T>& zmat, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(zmat.cols()==1,"The number of columns for the 1D coordinate is not 1 but " <<zmat.cols()<<"!"  );
        GISMO_ASSERT(u.cols()==2,"The number of columns for the 2D surface coordinate is not 2 but " <<u.cols()<<"!");
        GISMO_ASSERT(u.rows()==1, "Multiple ("<<u.cols()<<") parametric points given, accepts only 1... " <<"!"     );

        index_t N = zmat.rows();
        result.resize(N,3);
        result.leftCols(2)=u.replicate(N,1);
        result.col(2) = zmat;
    }
};

// To Do:
// * struct for material model
// Input is parametric coordinates of the surface \a mp
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
    mutable gsMatrix<T> Emat,Nmat, pts;
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
        _tmp.points = u.leftCols(2);

        static_cast<const gsFunction<T>*>(_mp)->computeMap(_tmp);

        // NOTE 2: in the case that parametric value is needed it suffices
        // to evaluate Youngs modulus and Poisson's ratio at
        // \a u instead of _tmp.values[0].
        _YoungsModulus->eval_into(_tmp.values[0], Emat);
        _PoissonRatio->eval_into(_tmp.values[0], Nmat);

        result.resize( targetDim() , u.cols() );
        for( index_t i=0; i< u.cols(); ++i )
        {
            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(i,3,3);

            F0.leftCols(2) = _tmp.jacobian(i);
            F0.col(2)      = _tmp.normal(i).normalized();
            F0 = F0.inverse();
            F0 = F0 * F0.transpose(); //3x3

            // Evaluate material properties on the quadrature point
            E = Emat(0,i);
            nu = Nmat(0,i);
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

    // piece(k) --> for patch k

};


// To Do:
// * struct for material model
// Input is parametric coordinates of the surface \a mp
template <class T>
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
    mutable gsMatrix<real_t,3,3> F0, jacGdef, jacGori;
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
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
        _tmp_def.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    ~gsMaterialMatrixIncompressible() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrixIncompressible)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 9;}

    mutable gsMaterialMatrixIncompressible<T> * _mm_piece; // todo: improve the way pieces are accessed

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
        _tmp.points = u;
        _tmp_def.points = u;

        real_t z = 0; // parametric coordinate in thickness direction

        static_cast<const gsFunction<T>*>(_mp)->computeMap(_tmp);
        static_cast<const gsFunction<T>*>(_mp_def)->computeMap(_tmp_def);

        // NOTE 2: in the case that parametric value is needed it suffices
        // to evaluate Youngs modulus and Poisson's ratio at
        // \a u instead of _tmp.values[0].
        _par1->eval_into(_tmp.values[0], par1mat);
        // _par2->eval_into(_tmp.values[0], par2mat);

        result.resize( targetDim() , u.cols() );
        for( index_t i=0; i< u.cols(); ++i )
        {
            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(i,3,3);

            // Compute metric tensor gab deformed: jacGdef
            jacGdef.leftCols(2) = _tmp_def.jacobian(i);                 // The in-plane surface tangents and...
            jacGdef.col(2)      = _tmp_def.normal(i).normalized();      // The normal vector form the covariant basis, which is stored in the cols of jacGdef
            jacGdef = jacGdef.inverse();                                // Now we have the inverse of this matrix
            jacGdef = jacGdef * jacGdef.transpose(); //3x3              // And this is the contravariant metric tensor

            // Compute metric tensor gab undeformed (original): jacGori
            jacGori.leftCols(2) = _tmp.jacobian(i);                     // The in-plane surface tangents and...
            jacGori.col(2)      = _tmp.normal(i).normalized();          // The normal vector form the covariant basis, which is stored in the cols of jacGori
            jacGori = jacGori.inverse();                                // Now we have the inverse of this matrix
            jacGori = jacGori * jacGori.transpose(); //3x3              // And this is the contravariant metric tensor

            // Evaluate material properties on the quadrature point
            mu = par1mat(0,i);
            J0 = math::sqrt( jacGdef.determinant() / jacGori.determinant() );
            J0 = math::pow( J0, -2 );

            /*
                C =     C1111,  C1122,  C1112
                        symm,   C2222,  C2212
                        symm,   symm,   C1212
            */

            C(0,0) = 2*F0(0,0)*F0(0,0) + F0(0,0)*F0(0,0) + F0(0,0)*F0(0,0);                 // C1111
            C(1,0) =
            C(0,1) = 2*F0(0,0)*F0(0,0) + F0(0,0)*F0(0,0) + F0(0,0)*F0(0,0);                 // C1122
            C(2,0) =
            C(0,2) = 2*F0(0,0)*F0(0,1) + F0(0,0)*F0(0,1) + F0(0,1)*F0(0,0);                 // C1112
            C(1,1) = 2*F0(1,1)*F0(1,1) + F0(1,1)*F0(1,1) + F0(1,1)*F0(1,1);                 // C2222
            C(2,1) =
            C(1,2) = 2*F0(1,1)*F0(0,1) + F0(1,0)*F0(1,1) + F0(1,1)*F0(1,0);                 // C2212
            C(2,2) = 2*F0(0,1)*F0(0,1) + F0(0,0)*F0(1,1) + F0(0,1)*F0(1,0);                 // C1212

            C *= mu * J0;

            //gsDebugVar(C);



            // for completeness, the computation of variable S^\alpha\beta is given here, commented. It should be implemented later in a separate class or function
             // TRANSFORM JACOBIANS!
            // Sab(0,0) = mu * (jacGori(0,0) - math::pow(J0,-2.) * jacGdef(0,0) );
            // Sab(0,1) =
            // Sab(1,0) = mu * (jacGori(1,0) - math::pow(J0,-2.) * jacGdef(1,0) ); // CHECK SYMMETRIES
            // Sab(1,1) = mu * (jacGori(1,1) - math::pow(J0,-2.) * jacGdef(1,1) );

        }
    }

};

// To Do:
// * struct for material model
// Input is parametric coordinates of the surface \a mp
template <class T>
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
    mutable gsMatrix<real_t,3,3> F0, jacGdef, jacGori;
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
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
        _tmp_def.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    ~gsMaterialMatrixCompressible() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrixCompressible)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 9;}

    mutable gsMaterialMatrixCompressible<T> * _mm_piece; // todo: improve the way pieces are accessed

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
        _tmp.points = u;
        _tmp_def.points = u;

        static_cast<const gsFunction<T>*>(_mp)->computeMap(_tmp);
        static_cast<const gsFunction<T>*>(_mp_def)->computeMap(_tmp_def);

        // NOTE 2: in the case that parametric value is needed it suffices
        // to evaluate Youngs modulus and Poisson's ratio at
        // \a u instead of _tmp.values[0].
        _par1->eval_into(_tmp.values[0], par1mat);
        _par2->eval_into(_tmp.values[0], par2mat);

        result.resize( targetDim() , u.cols() );
        for( index_t i=0; i< u.cols(); ++i )
        {
            // Material parameters
            mu = par1mat(0,i);
            K = par1mat(0,i);

            // Define objects
            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(i,3,3);
            gsMatrix<T,3,3> c, cinv;
            T S33, C3333, dc33, traceC;

            // Compute metric tensor gab deformed: jacGdef
            jacGdef.leftCols(2) = _tmp_def.jacobian(i);                 // The in-plane surface tangents and...
            jacGdef.col(2)      = _tmp_def.normal(i).normalized();      // The normal vector form the covariant basis, which is stored in the cols of jacGdef
            jacGdef = jacGdef * jacGdef.transpose(); //3x3              // And this is the covariant metric tensor

            // Compute metric tensor gab undeformed (original): jacGori
            jacGori.leftCols(2) = _tmp.jacobian(i);                     // The in-plane surface tangents and...
            jacGori.col(2)      = _tmp.normal(i).normalized();          // The normal vector form the covariant basis, which is stored in the cols of jacGori
            jacGori = jacGori * jacGori.transpose(); //3x3              // And this is the covariant metric tensor

            // Initialize c
            c.setZero();
            c.block(0,0,2,2) = jacGdef.block(0,0,2,2);
            c(2,2) = 1.0; // c33
            cinv = c.inverse();
            // note: can also just do c = jacGdef because the normal has length one and hence c(2,2) is 1. CHECK!

            J0 = math::sqrt( jacGdef.determinant() / jacGori.determinant() );
            J = J0 * math::sqrt( c(2,2) );

            index_t imax = 20;
            T tol = 1e-6;
            S33 = 0.0;
            C3333 = 1.0;

            // Define lambda function for C
            std::function<T (index_t i, index_t j, index_t k, index_t l)> Cijkl;
            Cijkl = [=](index_t i, index_t j, index_t k, index_t l)
            {
                T res = 1.0 / 9.0 * mu * math::pow( J , -2.0/3.0 ) * ( traceC * ( 2*cinv(i,j)*cinv(k,l) + 3*cinv(i,k)*cinv(j,l) + 3*cinv(i,l)*cinv(j,k) )
                                - 6*jacGori(i,j)*cinv(k,l) + cinv(i,j)*jacGori(k,l) ) + K * ( J*J*cinv(i,j)*cinv(k,l) - 0.5*(J*J-1)*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) ) );
                return res;
            };

            for (index_t i = 0; i < imax; i++)
            {
                dc33 = -2. * S33 / C3333;
                c(2,2) += dc33;
                cinv(2,2) = 1.0/c(2,2);

                traceC = c.trace();
                J = J0 * math::sqrt( c(2,2) );

                S33     = mu * math::pow( J , -2.0/3.0 ) * ( jacGori(2,2) - 1.0/3.0 * traceC * cinv(2,2) ) + 0.5 * K * ( J*J - 1 ) * cinv(2,2);
                C3333   = Cijkl(2,2,2,2);

                if (S33 < tol)
                {
                    gsInfo<<"Converged, S33 = "<<S33<<" and tolerance = "<<tol<<"\n";

/*
                    S.at(0) = mu * math::pow( J , -2.0/3.0 ) * ( jacGori(0,0) - 1.0/3.0 * traceC * cinv(0,0) ) + 0.5 * K * ( J*J - 1 ) * cinv(0,0); // S11
                    S.at(1) = mu * math::pow( J , -2.0/3.0 ) * ( jacGori(1,1) - 1.0/3.0 * traceC * cinv(1,1) ) + 0.5 * K * ( J*J - 1 ) * cinv(1,1); // S22
                    S.at(2) = mu * math::pow( J , -2.0/3.0 ) * ( jacGori(0,1) - 1.0/3.0 * traceC * cinv(0,1) ) + 0.5 * K * ( J*J - 1 ) * cinv(0,1); // S12
                    // S(i,j) = mu * math::pow( J , -2.0/3.0 ) * ( jacGori(i,j) - 1.0/3.0 * traceC * cinv(i,j) ) + 0.5 * K * ( J*J - 1 ) * cinv(i,j);
*/
                    /*
                        C =     C1111,  C1122,  C1112
                                symm,   C2222,  C2212
                                symm,   symm,   C1212
                        Here, Cabcd = Cijkl - Cab33*C33cd / C3333;
                        a,b,c,d = 1,2; i,j,k,l = 1...3;
                    */
                    C(0,0) = Cijkl(0,0,0,0) - ( Cijkl(0,0,2,2) * Cijkl(2,2,0,0) ) / (Cijkl(2,2,2,2)); // C1111
                    C(0,1) =
                    C(1,0) = Cijkl(0,0,1,1) - ( Cijkl(0,0,2,2) * Cijkl(2,2,1,1) ) / (Cijkl(2,2,2,2)); // C1122
                    C(0,2) =
                    C(2,0) = Cijkl(0,0,0,1) - ( Cijkl(0,0,2,2) * Cijkl(2,2,0,1) ) / (Cijkl(2,2,2,2)); // C1112
                    C(1,1) = Cijkl(1,1,1,1) - ( Cijkl(1,1,2,2) * Cijkl(2,2,1,1) ) / (Cijkl(2,2,2,2)); // C2222
                    C(1,2) =
                    C(2,1) = Cijkl(1,1,0,1) - ( Cijkl(1,1,2,2) * Cijkl(2,2,0,1) ) / (Cijkl(2,2,2,2)); // C2212
                    C(2,2) = Cijkl(0,1,0,1) - ( Cijkl(0,1,2,2) * Cijkl(2,2,0,1) ) / (Cijkl(2,2,2,2)); // C1212
                }
                else if (i == imax - 1)
                {
                    gsInfo<<"Error: Method did not converge";
                    std::terminate();
                }
            }
        }
    }
};
