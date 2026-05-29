// /** @file voigt_expr.h

//     @brief Provides an expressions to put another expression in Voigt notation

//     This file is part of the G+Smo library.

//     This Source Code Form is subject to the terms of the Mozilla Public
//     License, v. 2.0. If a copy of the MPL was not distributed with this
//     file, You can obtain one at http://mozilla.org/MPL/2.0/.

//     Author(s):
//         H.M.Verhelst
// */

// #pragma once

// #include <gsPhaseFieldFracture/gsMaterialBase.h>
// #include <gsPhaseFieldFracture/gsMaterialUtils.h>
// #include <gsExpressions/gsExpressions.h>

// namespace gismo
// {
//     namespace expr
//     {
//         /**
//         Transforms a matrix expression into a vector expression by computing the vector
//         [ a b c+d]^T
//         for each matrix block
//         [ a d ]
//         [ c b ]
//         */
//         template<short_t DIM, class E>
//         class voigt_expr  : public _expr<voigt_expr<DIM,E> >
//         {
//         public:
//             typedef typename E::Scalar Scalar;
//             enum {ScalarValued = 0, Space = E::Space, ColBlocks= 0}; // to do: ColBlocks
//         private:
//             typename E::Nested_t _u;
//             mutable gsMatrix<Scalar> tmp, result;

//             template <short_t dim>
//             inline typename std::enable_if<dim == 2, void>::type
//             _voigtStrain_impl(typename gsMatrix<Scalar>::ColXpr res,
//                             const gsMatrix<Scalar> & K) const
//             {
//                 // For 2D: [xx, yy, 2*xy]
//                 res(0) = K(0,0);
//                 res(1) = K(1,1);
//                 res(2) = K(0,1) + K(1,0);
//             }

//             template <short_t dim>
//             inline typename std::enable_if<dim == 3, void>::type
//             _voigtStrain_impl(typename gsMatrix<Scalar>::ColXpr res,
//                             const gsMatrix<Scalar> & K) const
//             {
//                 // For 3D: [xx, yy, zz, 2*xy, 2*xz, 2*yz]
//                 res(0) = K(0,0);
//                 res(1) = K(1,1);
//                 res(2) = K(2,2);
//                 res(3) = K(0,1) + K(1,0);
//                 res(4) = K(0,2) + K(2,0);
//                 res(5) = K(1,2) + K(2,1);
//             }

//         public:

//             voigt_expr(_expr<E> const& u) : _u(u)
//             {
//                 GISMO_ASSERT(_u.rows()==DIM, "Wrong dimension of voigt expression: expected "<<DIM<<" but got "<<_u.rows());
//             }

//             const gsMatrix<Scalar> & eval(const index_t k) const
//             {
//                 tmp = _u.eval(k);
//                 const index_t numActives = _u.cardinality(); // basis functions
//                 result.resize(DIM*(DIM+1)/2,numActives);

//                 for (index_t i = 0; i<numActives; ++i)
//                     _voigtStrain_impl<DIM>(result.col(i),tmp.block(0,i*_u.cols(),_u.rows(),_u.cols()));

//                 // Why is this done????
//                 if ( 1==Space )
//                     result.transposeInPlace();
//                 else if (2!=Space) // if not colSpan and not rowSpan
//                     result.transposeInPlace();


//                 return result;
//             }

//             // Per basis function
//             index_t rows() const { return 1; }
//             index_t cols() const { return DIM*(DIM+1)/2; } // dim*(dim+1)/2

//             void parse(gsExprHelper<Scalar> & evList) const
//             { _u.parse(evList); }

//             const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
//             const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
//             index_t cardinality_impl() const { return _u.cardinality_impl(); }

//             void print(std::ostream &os) const { os << "voigt("; _u.print(os); os<<")"; }
//         };

//         template <short_t DIM, typename E> EIGEN_STRONG_INLINE
//         voigt_expr<DIM,E> const voigt(E const & u)
//         { return voigt_expr<DIM,E>(u); }
//     }
// }
