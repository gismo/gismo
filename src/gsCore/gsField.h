/** @file gsField.h

    @brief Provides declaration of the Field class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsGeometry.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>
#include <gsUtils/gsPointGrid.h>

namespace gismo
{
/**
 * \brief A scalar of vector field defined on a m_parametric geometry.
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
 *
 * \ingroup Core
 */
template<class T>
class gsField
{

public:
    typedef memory::shared_ptr< gsField >  Ptr;// todo: remove
    typedef memory::unique_ptr< gsField >  uPtr;// todo: remove

    gsField(): m_patches(NULL) { }

    gsField( const gsFunctionSet<T> & mp,
             typename gsFunctionSet<T>::Ptr fs,
             const bool isparam)
    : m_patches(&mp), m_fields(fs), m_parametric(isparam)
    { }

    gsField( const gsGeometry<T> & sp, const gsFunctionSet<T> & pf, const bool isparam = false)
    : m_patches(&sp), m_fields(memory::make_shared_not_owned(&pf)), m_parametric(isparam)
    { }

    gsField( const gsGeometry<T> & sp, const gsGeometry<T> & pf)
    : m_patches(&sp), m_fields(memory::make_shared_not_owned(&pf)), m_parametric(true)
    { }

    gsField( const gsMultiPatch<T> & mp, const gsFunctionSet<T> & f, const bool isparam = false)
    : m_patches(&mp), m_fields(memory::make_shared_not_owned(&f)), m_parametric(isparam)
    { }

    gsField( const gsMultiPatch<T> & mp, const gsMultiPatch<T> & f)
    : m_patches(&mp), m_fields(memory::make_shared_not_owned(&f)), m_parametric(true)
    { }

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
     * @returns uPtr The <em>j</em>-th column of \a uPtr corresponds
     * to the image of the point \a u_j (which is defined by the \a j-th column of the input parameter \a u).
     */
    gsMatrix<T> point(const gsMatrix<T>& u, int i = 0) const
    {
        return m_patches->piece(i).eval(u);
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
     * @returns uPtr The <em>j</em>-th column of \a uPtr corresponds
     * to the value of the field at the point \a u_j (which is defined by the \a j-th column of the input parameter \a u).
     */
    gsMatrix<T> value(const gsMatrix<T>& u, int i = 0)  const
    {
        return m_parametric
            ? m_fields->piece(i).eval(u)
            : m_fields->piece(i).eval( point(u, i) );
    }

    // Return the value of the Field at physical value u
    // TO DO: rename to evalPhys()
    gsMatrix<T> pvalue(const gsMatrix<T>& u, int i)  const
    {
        GISMO_ASSERT(!m_parametric, "Cannot compute physical value");
        return ( m_fields->piece(i).eval(u) );
    }

    /// Computes the L2-distance between the two fields, on the physical domain
    T distanceL2(gsField<T> const & field, int numEvals= 1000) const;

    /// Computes the L2-distance between the field and a function \a func on the physical domain
    T distanceL2(gsFunctionSet<T> const & func,
                 bool isFunc_param = false,
                 int numEvals=1000) const;

    /// Computes the L2-distance between the field and a function \a
    /// func on the physical domain, using mesh from B
    T distanceL2(gsFunctionSet<T> const & func,
                 gsMultiBasis<T> const & B,
                 bool isFunc_param = false,
                 int numEvals=1000) const;

    /// Computes the H1-seminorm of the diff. between the field and a function \a
    /// func on the physical domain
    T distanceH1(gsFunctionSet<T> const & func,
                 bool isFunc_param = false,
                 int = 1000) const;

    /// Computes the H1-seminorm of the diff. between the field and a function \a
    /// func on the physical domain, using mesh from B
    T distanceH1(gsFunctionSet<T> const & func,
                 gsMultiBasis<T> const & B,
                 bool isFunc_param = false,
                 int = 1000) const;

    /// Computes the H2-seminorm of the diff. between the field and a function \a
    /// func on the physical domain, using mesh from B
    T distanceH2(gsFunctionSet<T> const & func,
                 bool isFunc_param = false) const;

    /// Computes the DG-distance between the field and a function \a
    /// func on the physical domain
    T distanceDG(gsFunctionSet<T> const & func,
                 bool isFunc_param = false,
                 int = 1000) const;

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << ( m_parametric ? "Parameterized f" : "F")
           << "unction field.\n Defined on " << m_patches;
        return os;
    }

    /// \brief Returns the dimension of the parameter domain
    /// (e.g., if the domain is a surface in three-dimensional space, it returns 2).
    short_t parDim() const { return m_patches->domainDim(); }

    /// \brief Returns the dimension of the physical domain
    /// (e.g., if the domain is a surface in three-dimensional space, it returns 3).
    short_t geoDim() const { return m_patches->targetDim(); }

    /// \brief Returns the dimension of the physical domain
    /// (e.g., if the domain is a surface in three-dimensional space, it returns 3).
    short_t dim() const { return m_fields->targetDim(); }

    /// Returns the number of patches.
    GISMO_DEPRECATED index_t nPatches()  const { return m_patches->nPieces(); }

    /// Returns the number of pieces.
    int nPieces()  const { return m_patches->nPieces(); }

    const gsGeometry<T> & geometry() const
    {
        GISMO_ASSERT(dynamic_cast<const gsGeometry<T>*>(m_patches),
                     "No geometry in field. The domain is"<< *m_patches);
        return *static_cast<const gsGeometry<T>*>(m_patches);
    }

    /// Returns gsMultiPatch containing the geometric information on the domain.
    const gsMultiPatch<T> & patches() const
    {
        GISMO_ASSERT(dynamic_cast<const gsMultiPatch<T>*>(m_patches),
                     "No patches in field. The field domain is "<< *m_patches);
        return *static_cast<const gsMultiPatch<T>*>(m_patches);
    }

    /// Returns the fields (defined per patch)
    const gsFunctionSet<T> & fields() const { return *m_fields; }

    /// Returns the gsGeometry of patch \a i.
    const gsGeometry<T> & patch(int i=0) const
    {
        GISMO_ASSERT( i<nPieces(),
                      "gsField: Invalid patch index.");
        GISMO_ASSERT(dynamic_cast<const gsGeometry<T>*>(&m_patches->piece(i)),
                     "No geometry in field. The domain is"<< m_patches->piece(i));
        return static_cast<const gsGeometry<T>&>(m_patches->piece(i));
    }

    /// Returns the gsFunction of patch \a i.
    const gsFunction<T> & function(int i=0) const
    {
        GISMO_ASSERT(dynamic_cast<const gsFunction<T>*>(&m_fields->piece(i)),
                     "No function in field. The domain is"<< m_patches->piece(i));
        return static_cast<const gsFunction<T>&>(m_fields->piece(i));
    }

    /// Attempts to return an Isogeometric function for patch i
    const gsGeometry<T> & igaFunction(int i=0) const
    {
        GISMO_ASSERT(i<m_fields->nPieces(),
                     "gsField: Invalid patch index.");
        GISMO_ASSERT(m_parametric,
                     "Cannot get an IGA function from non-parametric field.");
        GISMO_ASSERT(dynamic_cast<const gsGeometry<T>*>(&m_fields->piece(i)),
                     "Cannot return an igaFunction from a function of type: "<< m_fields->piece(i) );
        return static_cast<const gsGeometry<T> &>(m_fields->piece(i));
    }

    /// True if the field function is defined on parametric
    /// coordinates (i.e. the same coordinate system as the patches)
    bool isParametric() const { return m_parametric; }

    /// True if the field function is parametrized by a set of basis
    /// functions in parametric coordinates (i.e. the same coordinate
    /// system as the patches)
    bool isParametrized() const
    { return m_parametric && dynamic_cast<const gsGeometry<T>*>(&m_fields->piece(0));}

    /** \brief Returns the coefficient vector (if it exists)
        corresponding to the function field for patch \a i.

    Returns the coefficients of the field corresponding to the \a i-th
    patch. This is only possible in the case when the field is defined
    in terms of basis functions (ie. it derives from gsGeometry).

    */
    const gsMatrix<T> & coefficientVector(int i=0) const
    {
        return igaFunction(i).coefs();
    }

// Data members
private:

    /// The isogeometric field is defined on this multipatch domain
    const gsFunctionSet<T> * m_patches;

    // If there are many patches, one field per patch

    /// \brief Vector containing "local fields" for each patch/subdomain.
    ///
    /// For each patch/subdomain, the "local field" is represented by
    /// a gsFunction. This local field can be accessed with
    /// gsField::function.
    typename gsFunctionSet<T>::Ptr m_fields;

    /**
     * @brief \a True iff this is an isogeometric field.
     *
     * If \a m_parametric is \a true, the evaluation points for calling gsField::value have to be placed in the
     * \a parameter domain.
     *
     * If \a m_parametric is \a false, then the evaluation points are in the \a physical domain.
     * This applies to, e.g., given exact solutions which are defined on the physical domain.
     */
    bool m_parametric;// True iff this is an Isogeometric field, living on parameter domain

}; // class gsField


/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsField<T>& b)
{return b.print(os); }


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsField.hpp)
#endif
