/** @file gsCPPInterface.h

    @brief Provides a mapping between the corresponding sides of two patches sharing an interface,
    by means of a Closest Point Projection.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):  A. Mantzaflaris, C. Karampatzakis
    Created on: 2022-06-28
*/


#pragma once

#include <gsCore/gsAffineFunction.h>
#include <gsCore/gsBoundary.h>
#include <gsIO/gsOptionList.h>

namespace gismo {

/** @brief Provides a mapping between the corresponding sides of two patches sharing an interface,
  * by means of a closest point projection.  
  *
  *
  * @ingroup Assembler
 **/
template <class T>
class gsCPPInterface : public gsFunction<T>
{
public:

    /// Shared pointer for gsRemapInterface
    typedef memory::shared_ptr< gsCPPInterface > Ptr;

    /// Unique pointer for gsRemapInterface
    typedef memory::unique_ptr< gsCPPInterface > uPtr;

    /// @brief Constructor
    ///
    /// @param  mp                 The multi-patch object
    /// @param  mb                 The multi-basis object
    /// @param  bi                 The boundary interface (specifying the patches and their \a patchSide s)
    /// @param  opt                Options; see below:
    ///
    ///  | Option                  | Meaning
    ///  |------------------------ | --------------------------------------------------------------------------------------
    ///  | CheckAffine             | The number of interior points (per direction) in point grid used to
    ///  |                         | default: 1e-10
    gsCPPInterface(const gsMultiPatch<T>   & mp,
                     const gsMultiBasis<T>   & mb,
                     const boundaryInterface & bi,
                     const gsOptionList      & opt = defaultOptions() );

    static uPtr make (const gsMultiPatch<T>   & mp,
                     const gsMultiBasis<T>   & mb,
                     const boundaryInterface & bi,
                     const gsOptionList      & opt = defaultOptions() )
    { return uPtr(new gsCPPInterface(mp,mb,bi,opt)); }
    
private:
    const gsGeometry<T>* m_slaveGeom, *m_masterGeom; ///< Geometry of first (Slave) patch and second patch (master) 
    typename gsGeometry<T>::Ptr m_masterBdr;  ///< The boundary geometry of second patch -- Master

    const gsBasis<T> * m_slaveBasis ;        ///< Basis on first patch
    const gsBasis<T> * m_masterBasis;        ///< Basis on second patch

    boundaryInterface m_boundaryInterface;   ///< Corresponding boundary interface

    T m_Tolerance;                            ///< Tolerance for closest point algorithm
    
    std::vector<index_t> m_freeDirs;
    index_t m_fixedParam;
    index_t m_fixedDir;

    std::vector< std::vector<T> > m_breakpoints;  ///< Union of breakpoints of both bases

public:
    /// @Returns default options
    ///
    /// See constructor for their meaning
    static gsOptionList defaultOptions()
    {
        gsOptionList defaultOpts;
        defaultOpts.addInt ( "Tied",      "Is it a tied contact?", 0);
        defaultOpts.addReal( "Tolerance", "Tolerance to be used for closest point search.", 1e-5);
        return defaultOpts;
    }


    /// @brief Evaluates the interface map
    ///
    /// Interfaces map \f$ \widehat \Gamma_1 \rightarrow \widehat \Gamma_2 \f$ that represents
    /// \f$ G_2^{-1} \circ G_1 \f$
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const;

    /// Returns parameter dimension of the domains
    virtual short_t domainDim() const { return m_slaveGeom->domainDim(); }

    /// @brief Returns a domain iterator
    ///
    /// The domain iterator lives on \f$ \widehat \Gamma_1 \f$. Its break points are the union of
    /// the brakpoints of the basis on \f$ \widehat \Omega_1 \f$ and the breakpoints of the basis
    /// on \f$ \widehat \Omega_2 \f$, mapped to \f$ \widehat \Omega_1 \f$.
    typename gsDomainIterator<T>::uPtr makeDomainIterator() const;


    /// Returns the break points used in \ref makeDomainIterator
    const std::vector< std::vector<T> > & breakPoints() const { return m_breakpoints; }

    /// Prints the state of the object
    virtual std::ostream & print(std::ostream& os) const;

private:
    /// Computes the box which represents the intersection of sides of incoming patches
    void constructInterfaceBox();

    /// Constructs the breakpoints \a m_breakpoints
    void constructBreaks();
}; // End gsCPPInterface


/// @brief   Prints the state of the object
/// @relates gsCPPInterface
template <class T>
inline std::ostream & operator<<(std::ostream & os, const gsCPPInterface<T> & cppIf)
{ return cppIf.print(os); }

} // End namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCPPInterface.hpp)
#endif
