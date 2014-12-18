
#pragma once

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsHBSpline.h>

namespace gismo
{
    /** 
     * \brief
     * A hierarchical B-spline basis of parametric dimension \em d.
     *
     * See \cite Kraft1997 for the theory behind this kind of basis.
     * 
     * \tparam d the dimension of the parameter domain
     * \tparam T coefficient type
     *
     * \ingroup basis
    */ 
    
template<unsigned d, class T>
class gsHBSplineBasis : public gsHTensorBasis<d,T>
{
public:
    /// Associated geometry type
    typedef gsHBSpline<d,T> GeometryType;
    
    typedef typename gsHTensorBasis<d,T>::CMatrix CMatrix;

    typedef memory::shared_ptr< gsHBSplineBasis > Ptr;

public:
    /// Constructor out of a gsBSplineBasis
    gsHBSplineBasis(gsBSplineBasis<T> &  bsbasis, int nlevels = 10)
        : gsHTensorBasis<d,T>( gsTensorBSplineBasis<d,T>(&bsbasis), nlevels)
    {
        GISMO_ASSERT(d==1, "Wrong dimension");
    }
    
    gsHBSplineBasis( gsBSplineBasis<T> &  bsbasis,
                     std::vector<unsigned> boxes, int nlevels = 10)
        : gsHTensorBasis<d,T>( gsTensorBSplineBasis<d,T>(&bsbasis), nlevels, boxes) 
    {
        GISMO_ASSERT(d==1, "Wrong dimension");
    }
    
    gsHBSplineBasis( gsBSplineBasis<T> &  bsbasis,
                     gsMatrix<T> const & boxes, int nlevels = 10)
        : gsHTensorBasis<d,T>(gsTensorBSplineBasis<d,T>(&bsbasis), nlevels, boxes) 
    {
        GISMO_ASSERT(d==1, "Wrong dimension");
    }
    
    gsHBSplineBasis( gsBSplineBasis<T> &  bsbasis,
                     gsMatrix<T> const & boxes, std::vector<unsigned int> levels, int nlevels = 10)
        : gsHTensorBasis<d,T>(gsTensorBSplineBasis<d,T>(&bsbasis), nlevels, boxes, levels)
    {
        GISMO_ASSERT(d==1, "Wrong dimension");
    }

    /// Constructor out of a tensor BSpline Basis
    gsHBSplineBasis(gsBasis<T> const&  tbasis, int nlevels = 10)
        : gsHTensorBasis<d,T>(tbasis, nlevels) 
    {
        // initialize(); // is done in the base constructor
    }
    
    gsHBSplineBasis( gsTensorBSplineBasis<d,T> const&  tbasis,
                     std::vector<unsigned> boxes, int nlevels = 10)
        : gsHTensorBasis<d,T>(tbasis, nlevels, boxes) 
    {
        // initialize(); // is done in the base constructor
    }
    
    gsHBSplineBasis( gsTensorBSplineBasis<d,T> const&  tbasis,
                     gsMatrix<T> const & boxes, int nlevels = 10)
        : gsHTensorBasis<d,T>(tbasis, nlevels, boxes) 
    {
        // initialize(); // is done in the base constructor
    }
    
    gsHBSplineBasis( gsTensorBSplineBasis<d,T> const&  tbasis,
                     gsMatrix<T> const & boxes, std::vector<unsigned int> levels, int nlevels = 10)
        : gsHTensorBasis<d,T>(tbasis, nlevels, boxes, levels)
    {
        // initialize(); // is done in the base constructor
    }
    
public:
    
    int dim() const { return d; }
    
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;
    
    void evalSingle_into  (unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const;
    
    void derivSingle_into (unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const;
    
    void deriv2Single_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const;
    
    virtual gsHBSplineBasis* clone() const
    { return new gsHBSplineBasis(*this); }
    
    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const;
    ///returns transfer matrices betweend the levels of the given hierarchical spline
    void transferbyLvl(std::vector<gsMatrix<T> >& result);

    GISMO_MAKE_GEOMETRY_NEW
    
private:
    
    /// Initialize the characteristic and coefficient matrices and the
    /// internal bspline representations
    void initialize();

    gsMatrix<T> coarsening(const std::vector<gsSortedVector<unsigned> >& old, const std::vector<gsSortedVector<unsigned> >& n, const gsSparseMatrix<T,RowMajor> & transfer);
    gsMatrix<T> coarsening_direct( const std::vector<gsSortedVector<unsigned> >& old, const std::vector<gsSortedVector<unsigned> >& n, const std::vector<gsSparseMatrix<T,RowMajor> >& transfer);
    
}; // class gsHBSplineBasis


//////////////////////////////////////////////////
//////////////////////////////////////////////////

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHBSplineBasis.hpp)
#endif
