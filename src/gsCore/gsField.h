 /** @file gsField.h

    @brief Provides declaration of the Field class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsGeometry.h>
#include <gsCore/gsMultiPatch.h>
#include <gsUtils/gsNorms.h>
#include <gsCore/gsMultiBasis.h>

namespace gismo
{

  /**
   * \brief A scalar of vector field defined on a parametrized geometry.
   *
   * A gsField is, generally speaking, some mathematical function that is defined on a domain of interest
   * (the name "field" is motivated by, e.g., "scalar field" or "vector field").
   *
   * The gsField combines the following:\n
   * - <b>Geometric information</b> on the domain:\n
   *
   * The domain can be represented as one single patch or as a
   * collection of multiple patches (a.k.a. subdomains).\n This
   * information is stored in a member of the gsMultiPatch class.\n
   *
   * - The <b>function</b> defined on the domain:\n
   *
   * For each patch (a.k.a. subdomain), the gsField contains a member
   * of class gsFunction (which represents the "local field", so to
   * say).  On this, the operations of gsFunction can be carried out
   * (e.g., function evaluation or computation of derivatives).\n
   * Remark: The collection of patch-wise gsFunction is stored in the
   * private member gsField::m_fields.
   *
   * Note that the geometry representation of a single patch can be
   * extracted by calling the member function gsField::patch.
   *
   * The "local field" on a single patch can be extracted by calling gsField::function.
   */


//A field is a function defined on a geometry, or
//      over a multipatch collection
  
template<class T>
class gsField
{
public:
  /// Shared pointer for gsField
  typedef memory::shared_ptr< gsField >  Ptr;

  /// Unique pointer for gsField
  typedef memory::auto_ptr< gsField >   uPtr;

public:
    
    gsField( const gsMultiPatch<T> & mp, const std::vector<gsFunction<T> *>& fs, const bool isparam= true)
        : m_patches(mp), m_fields(fs), parametrized(isparam), m_owning(true)
        { }

    gsField( const gsMultiPatch<T> & mp, gsFunction<T> * f, const bool isparam= true) 
        : m_patches(mp), parametrized(isparam), m_owning(true)
        { 
            for (std::size_t i = 0; i< mp.nPatches(); ++i)
                m_fields.push_back(f); 
        }

    gsField( const gsMultiPatch<T> & mp, gsFunction<T> & f, const bool isparam= true) 
        : m_patches(mp), parametrized(isparam), m_owning(false)
        { 
            for (size_t i = 0; i< mp.nPatches(); ++i)
                m_fields.push_back(&f); 
        }

    ~gsField()
        {
            if (m_owning)
            {
                // avoid deleting twice the same pointer
                typename std::vector<gsFunction<T>*>::iterator itr= 
                    std::unique( m_fields.begin(), m_fields.end() );
                m_fields.resize( itr - m_fields.begin() ) ;

                freeAll( m_fields );
            }
        }
    
public:
    
// TO DO:

// EVAL_physical_position
// Need to solve x = m_geometry(u) for u

    // Return a point in  the physical domain at parameter value u
    /**
     * @brief Maps points \a u from the parameter domain to the physical domain.
     *
     * @param[in] u Evaluation points as gsMatrix of size <em>d</em> x <em>n</em>.\n
     * \a d denotes the dimension of the parameter domain (i.e., d = parDim()).\n
     * \a n denotes the number of evaluation points.\n
     * Each column of \a u corresponds to one evaluation point.
     * @param[in] i Index of the considered patch/subdomain.
     * @param[out] uPtr The <em>j</em>-th column of \a uPtr corresponds
     * to the image of the point \a u_j (which is defined by the \a j-th column of the input parameter \a u).
     */
    typename gsMatrix<T>::uPtr point(const gsMatrix<T>& u, int i = 0) const
    {
        return m_patches[i].eval(u);
    }

    // Return the value of the Field at parameter value u
    // TO DO: rename to evalParam()
    /**
     * @brief Evaluation of the field at points \a u.
     *
     * @param[in] u Evaluation points as gsMatrix of size <em>d</em> x <em>n</em>.\n
     * \a d denotes the dimension of the parameter domain (i.e., d = parDim()).\n
     * \a n denotes the number of evaluation points.\n
     * Each column of \a u corresponds to one evaluation point.
     * @param[in] i Index of the considered patch/subdomain.
     * @param[out] uPtr The <em>j</em>-th column of \a uPtr corresponds
     * to the value of the field at the point \a u_j (which is defined by the \a j-th column of the input parameter \a u).
     */
    typename gsMatrix<T>::uPtr value(const gsMatrix<T>& u, int i = 0)  const
    {
        return parametrized
            ? m_fields[i]->eval(u)
            : m_fields[i]->eval( *point(u, i) );
    }

    // Return the value of the Field at physical value u 
    // TO DO: rename to evalPhys()
    typename gsMatrix<T>::uPtr pvalue(const gsMatrix<T>& u, int i)  const
    { 
        assert( !parametrized );
        return ( m_fields[i]->eval(u) ); 
    }

    /// Computes the L2-distance between the two fields, on the physical domain
    T distanceL2(gsField<T> const & field, int numEvals= 1000) const 
    {
        return computeL2Distance(*this, field, numEvals);
    }

    /// Computes the L2-distance between the field and a function \a func on the physical domain
    T distanceL2(gsFunction<T> const & func, 
                 bool isFunc_param = false,
                 int numEvals=1000) const 
    {
        if ( parametrized ) // isogeometric field
            return igaFieldL2Distance(*this, func, isFunc_param);
        else
            return computeL2Distance(*this, func, isFunc_param,  numEvals);
    }

    /// Computes the L2-distance between the field and a function \a func on the physical domain
    T distanceL2(gsFunction<T> const & func,
                 gsMultiBasis<T> const & B,
                 bool isFunc_param = false,
                 int numEvals=1000) const
    {
        if ( parametrized ) // isogeometric field
            return igaFieldL2Distance(*this, func, B,isFunc_param);
        else
            return computeL2Distance(*this, func, isFunc_param,  numEvals);
    }

    /// Computes the H1-distance between the field and a function \a func on the physical domain
    T distanceH1(gsFunction<T> const & func, 
                 bool isFunc_param = false,
                 int numEvals=1000) const 
    {
        if ( parametrized ) // isogeometric field
            return igaFieldH1Distance(*this, func, isFunc_param);
        else
        {
            GISMO_UNUSED(numEvals);
            gsWarn <<"H1 seminorm not implemented.\n";
            return -1;
        }
    }

    /// Computes the DG-distance between the field and a function \a
    /// func on the physical domain
    T distanceDG(gsFunction<T> const & func, 
                 bool isFunc_param = false,
                 int numEvals=1000) const 
    {
        if ( parametrized ) // isogeometric field
            return igaFieldDGDistance(*this, func, isFunc_param);
        else
        {
            GISMO_UNUSED(numEvals);
            gsWarn <<"DG norm not implemented.\n";
            return -1;
        }
    }
    
    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
        { os << "gsField.\n"; return os; }
    
    /// \brief Returns the dimension of the parameter domain
    /// (e.g., if the domain is a surface in three-dimensional space, it returns 2).
    int parDim() const { return m_patches.parDim(); }

    /// \brief Returns the dimension of the physical domain
    /// (e.g., if the domain is a surface in three-dimensional space, it returns 3).
    int geoDim() const { return m_patches.geoDim(); }

    /// \brief Returns the dimension of the physical domain
    /// (e.g., if the domain is a surface in three-dimensional space, it returns 3).
    int dim() const { return m_fields[0]->targetDim(); }

    /// Returns the number of patches.
    int nPatches()  const { return m_patches.nPatches(); }

    gsGeometry<T>   & geometry() const { assert(m_patches.size()==1); return m_patches[0]; }

    /// Returns gsMultiPatch containing the geometric information on the domain.
    const gsMultiPatch<T> & patches() const    { return m_patches; }

    /// Returns the gsGeometry of patch \a i.
    const gsGeometry<T> & patch(int i=0) const       { return m_patches[i]; }

    /// Returns the gsFunction of patch \a i.
    // const todo
    gsFunction<T>   & function(int i=0) const  { return *m_fields[i]; }

    /// Attempts to return an Isogeometric function for patch i
    const gsGeometry<T> & igaFunction(int i=0) const
    { 
        GISMO_ASSERT( parametrized, "Cannot get an IGA function from non-parametric field.");
        return static_cast<gsGeometry<T> &>(*m_fields[i]);
    }

    bool isParametrized() const         { return parametrized; }

    // S.Kleiss: incomplete documentation.
    /// Returns coefficient vector corresponding to the function on patch \a i.
    gsMatrix<T> & coefficientVector(int i=0) const
    {
        GISMO_ASSERT( parametrized, "Coefficients do not exist.");
        GISMO_ASSERT( i < int( m_patches.nPatches() ) , "Index of patch exceeds number of patches.");
        return (dynamic_cast< gsGeometry<T> *>( m_fields[i] ) )->coefs();
    }


private:
    // disable copying
    gsField(const gsField& other);
    gsField& operator=(const gsField& other);

// Data members
private:

    /// dox to m_patches.
    const gsMultiPatch<T>& m_patches;

    // If there are many patches, one field per patch

    /// \brief Vector containing "local fields" for each patch/subdomain.
    ///
    /// For each patch/subdomain, the "local field" is represented by a gsFunction. This local field can be accessed with gsField::function.
    std::vector< gsFunction<T> *> m_fields;

    /**
     * @brief \a True iff this is an isogeometric field.
     *
     * If \a parametrized is \a true, the evaluation points for calling gsField::value have to be placed in the
     * \a parameter domain.
     *
     * If \a parametrized is \a false, then the evaluation points are in the \a physical domain.
     * This applies to, e.g., given exact solutions which are defined on the physical domain.
     */
    bool parametrized;// True iff this is an Isogeometric field, living on parameter domain

    bool m_owning;      // whether this field owns its function and should destroy it

}; // class gsField


//////////////////////////////////////////////////
//////////////////////////////////////////////////


/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsField<T>& b)
{return b.print(os); }


} // namespace gismo
