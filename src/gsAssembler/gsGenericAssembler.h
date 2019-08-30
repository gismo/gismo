/** @file gsGenericAssembler.h

    @brief Provides an assembler for common IGA matrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris
*/

#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsVisitorMass.h>
//#include <gsAssembler/gsVisitorTPmass.h>
#include <gsAssembler/gsVisitorGradGrad.h>
#include <gsAssembler/gsVisitorMoments.h>

#include <gsPde/gsLaplacePde.h>


#include <gsAssembler/gsExpressions.h>
#include <gsAssembler/gsExprHelper.h>
#include <gsAssembler/gsExprAssembler.h>

namespace gismo
{

/**
   @brief Assembles the mass, stiffness matrix on a given domain

   
   \ingroup Assembler
 */
template <class T>
class gsGenericAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:

    /// Constructor with gsMultiBasis
    gsGenericAssembler( const gsMultiPatch<T>    & patches,
                        const gsMultiBasis<T>    & bases,
                        const gsOptionList & opt = Base::defaultOptions(),
                        const gsBoundaryConditions<T> * bc = NULL)
    : m_pde(patches)
    {
        if ( bc != NULL)
        {
            m_pde.boundaryConditions() = *bc;
            this->m_ddof.resize(1);
        }
        
        Base::initialize(m_pde, bases, opt);
        gsGenericAssembler::refresh();
    }

    /* TODO (id geometry)
    // Constructor with gsMultiBasis
    gsGenericAssembler( gsMultiBasis<T> const         & bases,
                        bool conforming = false,
                        const gsBoundaryConditions<T> * bc = NULL)
    : Base(patches)
    {
        m_bases.push_back(bases);

        Base::initialize(m_pde, bases, opt);
        gsGenericAssembler::refresh();
    }
    */

    void refresh()
    {
        // Setup sparse system
        gsDofMapper mapper;
        m_bases[0].getMapper(
            (dirichlet::strategy)(m_options.getInt("DirichletStrategy")),
            (iFace::strategy)(m_options.getInt("InterfaceStrategy")),
            this->pde().bc(), mapper, 0);
        m_system = gsSparseSystem<T>(mapper);
        //note: no allocation here
        //        const index_t nz = m_options.numColNz(m_bases[0][0]);
        //        m_system.reserve(nz, 1);
    }

    /// Mass assembly routine
    const gsSparseMatrix<T> & assembleMass()
    {
        /*
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        
        // Elements used for numerical integration
        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(m_bases.front());
        geometryMap G = A.getMap(m_pde_ptr->patches());
        space u = A.getSpace(m_bases.front());

        A.initSystem();
        //m_system.matrix() = A.assemble( u * u.tr() * meas(G) );
        A.assemble( u * u.tr() * meas(G) );
        m_system.matrix() = A.matrix();
        
        return m_system.matrix();
        */

        // Clean the sparse system
        gsGenericAssembler::refresh();
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
        m_system.matrix().reservePerColumn(nz);
        
        // Assemble mass integrals
        //this->template push<gsVisitorMass<T> >();
        this->template push<gsVisitorMass<T> >();

        // Assembly is done, compress the matrix
        this->finalize();

        return m_system.matrix();
    }

        /// Mass assembly routine
    const gsSparseMatrix<T> & assembleMass2()
    {
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        
        // Elements used for numerical integration
        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(m_bases.front());
        geometryMap G = A.getMap(m_pde_ptr->patches());
        space u = A.getSpace(m_bases.front());

        A.initSystem();
        //m_system.matrix() = A.assemble( u * u.tr() * meas(G) );
        A.assemble( u * u.tr() * meas(G) );
        m_system.matrix() = A.matrix();
        
        return m_system.matrix();
    }

    /*// Mass assembly routine
    const gsSparseMatrix<T> & assembleMass3()
    {
        // Clean the sparse system
        gsGenericAssembler::refresh();
        
        //this->template push<gsVisitorTPmass<T> >();
        this->finalize();
        return m_system.matrix();
    }
    */

    /// Stiffness assembly routine
    const gsSparseMatrix<T> & assembleStiffness()
    {
        // Clean the sparse system
        gsGenericAssembler::refresh();
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
        m_system.matrix().reservePerColumn(nz);

        // Assemble stiffness integrals
        this->template push<gsVisitorGradGrad<T> >();

        // Assembly is done, compress the matrix
        this->finalize();

        return m_system.matrix();
    }

    /// Moments assembly routine
    const gsMatrix<T> & assembleMoments(const gsFunction<T> & func)
    {
        // Reset the right-hand side vector
        m_system.rhs().setZero(m_system.cols(), 1);

        // Assemble moment integrals
        gsVisitorMoments<T> mom(func);
        this->push(mom);
        
        // Assembly is done, compress the matrix
        this->finalize();

        return m_system.rhs();
    }

    /// Stiffness assembly routine on patch \a patchIndex
    const gsSparseMatrix<T> & assembleMass(index_t patchIndex)
    {
        gsGenericAssembler<T> tmp(m_pde.patches().patch(patchIndex), 
                                  m_bases[patchIndex], m_options);
        tmp.assembleMass();
        m_system.matrix().swap(tmp.m_system.matrix());
        return m_system.matrix();
    }

    /// Stiffness assembly routine on patch \a patchIndex
    const gsSparseMatrix<T> & assembleStiffness(index_t patchIndex)
    {
        gsGenericAssembler<T> tmp(m_pde.patches().patch(patchIndex), 
                                  m_bases[patchIndex],  m_options);
        tmp.assembleStiffness();
        m_system.matrix().swap(tmp.m_system.matrix());
        return m_system.matrix();
    }
    
    /// Returns an expression of the "full" assembled sparse
    /// matrix. Note that matrix() might return a lower diagonal
    /// matrix, if we exploit possible symmetry during assembly
    /// (check: m_matrix.symmetry() == true )
    typename gsSparseMatrix<T>::fullView fullMatrix()
    {
        return m_system.matrix().template selfadjointView<Lower>();
    }

    /// Returns an expression of the "full" assembled sparse
    /// matrix. Note that matrix() might return a lower diagonal
    /// matrix, if we exploit possible symmetry during assembly
    /// (check: m_matrix.symmetry() == true )
    const typename gsSparseMatrix<T>::constFullView fullMatrix() const
    {
        return m_system.matrix().template selfadjointView<Lower>();
    }
    

private:

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

private:

    gsLaplacePde<T> m_pde;
};



} // namespace gismo

