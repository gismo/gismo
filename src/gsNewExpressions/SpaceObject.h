/** @file SpaceObject.h

    @brief Defines the SpaceObject expression for finite element spaces

    SpaceObject represents a finite element space (test or trial space) in the
    expression system. It wraps a gsMultiBasis and handles:
    - Basis function evaluation at quadrature points
    - Derivative computation (gradients, hessians, etc.)
    - Boundary condition handling via setup()
    - DOF mapping for assembly
    - Component extraction for vector-valued spaces

    Key features:
    - Space type (Test/Trial) determines role in bilinear/linear forms
    - Order determines tensor rank (0=scalar field, 1=vector field, etc.)
    - setup() method applies boundary conditions and creates DOF mapper
    - isInitialized() checks if space has been set up for assembly

    Usage pattern:
    1. Create space: auto u = space<Test>(basis);
    2. Set up BCs: u.setup(bcInfo, dirichlet);
    3. Assemble: assembler.assemble(u * v);

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsMultiBasis.h>
#include <gsMSplines/gsMappedBasis.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsAssembler/gsAssemblerOptions.h>
#include <gsUtils/gsPointGrid.h>

#pragma once

namespace gismo
{
namespace Expr
{
template<class T>
struct FeSpaceData;

// Forward declaration for component space
template <class T, enum SpaceType _Space>
class ComponentSpaceObject;

template <class T, enum SpaceType _Space, size_t _Order>
struct ExpressionTraits<SpaceObject<T, _Space, _Order>>
{
    typedef T Scalar;
    static constexpr size_t Order = _Order;
    static constexpr SpaceType Space = _Space;
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = false;
};

template <class T, enum SpaceType _Space, size_t _Order>
class SpaceObject : public BaseObject<SpaceObject<T, _Space, _Order>>
{

    friend class NullObject<T, _Space, _Order>; 

    static_assert(_Space==SpaceType::None || _Space==SpaceType::Test || _Space==SpaceType::Trial, "SpaceObject can only be None, Test or Trial space.");

    using Base = BaseObject<SpaceObject<T, _Space, _Order>>;
    using Base::Deriv_;

public:
    // Expose the static traits publicly
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    typedef typename Base::Scalar Scalar;

private:
    const gsFunctionSet<Scalar> * m_fs;
    const gsFuncData<Scalar>    * m_fd;
    FeSpaceData<Scalar>         * m_sd; // Pointer to Space Data for assembly
    size_t m_id;

public:
    SpaceObject(size_t domainDim, const std::array<size_t, Order> & input_sizes, size_t id,
                std::string label=(Space==SpaceType::Test)?"φ":(Space==SpaceType::Trial)?"ψ":"UNKNOWN_SPACE")
    :
    Base(domainDim, input_sizes, label),
    m_fs(NULL), m_fd(NULL), m_sd(NULL), m_id(id)
    {
    }

    // ExpressionResult<Scalar> eval(const index_t k) const
    // {
    //     // Get the number of active basis functions
    //     index_t numActive = m_fd->values[0].rows();

    //     // Test space: cardinality is (numActive, 1) - no runtime branching
    //     ExpressionResult<Scalar> result((_Space==SpaceType::Test) ? numActive : 1, (_Space==SpaceType::Test) ? 1 : numActive);

    //     // Fill in the values for each basis function
    //     for (index_t i = 0; i < numActive; ++i)
    //     {
    //         result((_Space==SpaceType::Test) ? i : 0, (_Space==SpaceType::Test) ? 0 : i) = m_fd->values[0].col(k).segment(i, 1).blockDiag(Order ? this->sizes()[0] : 1);
    //     }
    //     return result;
    // }

    // Specialization for Test space using enable_if
    template <size_t S = _Space>
    typename std::enable_if<S == SpaceType::Test, ExpressionResult<Scalar>>::type
    eval(const index_t k) const
    {
        // Get the number of active basis functions
        index_t numActive = m_fd->values[0].rows();

        // Test space: cardinality is (numActive, 1) - no runtime branching
        ExpressionResult<Scalar> result(numActive, 1);

        // Fill in the values for each basis function
        for (index_t i = 0; i < numActive; ++i)
        {
            result(i, 0) = m_fd->values[0].col(k).segment(i, 1).blockDiag(Order ? this->sizes()[0] : 1);
        }
        return result;
    }

    // Specialization for Trial space using enable_if
    template <size_t S = _Space>
    typename std::enable_if<S == SpaceType::Trial, ExpressionResult<Scalar>>::type
    eval(const index_t k) const
    {
        // Get the number of active basis functions
        index_t numActive = m_fd->values[0].rows();

        // Trial space: cardinality is (1, numActive) - no runtime branching
        ExpressionResult<Scalar> result(1, numActive);

        // Fill in the values for each basis function
        for (index_t i = 0; i < numActive; ++i)
        {
            result(0, i) = m_fd->values[0].col(k).segment(i, 1).blockDiag(Order ? this->sizes()[0] : 1);
        }
        return result;
    }

    // Specialization for None space using enable_if
    template <size_t S = _Space>
    typename std::enable_if<S == SpaceType::None, ExpressionResult<Scalar>>::type
    eval(const index_t k) const
    {
        GISMO_ERROR("THIS CODE SHOULD NEVER BE CALLED. SpaceObject with Space=None should not be evaluated.");
        // Get the number of active basis functions
        index_t numActive = m_fd->values[0].rows();

        // None space: cardinality is (1, 1) - no runtime branching
        ExpressionResult<Scalar> result(1, 1);

        // Fill in the values (single entry)
        for (index_t i = 0; i < numActive; ++i)
        {
            result(0, 0) = m_fd->values[0].col(k).segment(i, 1).blockDiag(Order ? this->sizes()[0] : 1);
        }
        return result;
    }

    const gsFunctionSet<T> & source() const {return *m_fs;}
    const gsFuncData<T>    & data()   const
    {
        GISMO_ASSERT(NULL!=m_fd, "SpaceObject: invalid data "<< m_fs <<","<<m_fd);
        return *m_fd;
    }

    FeSpaceData<Scalar> * spaceData() const { return m_sd; }
    
    gsDofMapper & mapper()
    {
        GISMO_ASSERT(NULL!=m_sd, "Space/mapper not properly initialized.");
        return m_sd->mapper;
    }
    const gsDofMapper & mapper() const
    { return const_cast<SpaceObject*>(this)->mapper(); }
    
    inline const gsMatrix<T> & fixedPart() const { return m_sd->fixedDofs; }
    gsMatrix<T> & fixedPart() { return m_sd->fixedDofs; }

    index_t dim() const { return m_sd ? m_sd->dim : 1; }
    
    index_t interfaceCont() const { return m_sd->cont; }
    index_t & setInterfaceCont(const index_t _r) const
    {
        GISMO_ASSERT(_r>-2 && _r<1, "Invalid or not implemented (r="<<_r<<").");
        return m_sd->cont = _r;
    }

    void setSource(const gsFunctionSet<Scalar> & fs) { m_fs = &fs;}
    void setData(const gsFuncData<Scalar> & val) { m_fd = &val;}
    void setSpaceData(FeSpaceData<Scalar> & sd) { m_sd = &sd; }
    
    /// @brief Check if the space has been properly initialized (setup called)
    bool check() const
    {
        return m_sd != nullptr && m_sd->isInitialized();
    }
    
    /// @brief Get the space id
    size_t id() const { return m_id; }

    /**
     * @brief Component access for vector spaces (Order=1)
     * 
     * Returns a ComponentSpaceObject representing a scalar component of this vector space.
     * Only available for Order=1 (vector) spaces.
     * 
     * @param c Component index (0-based)
     * @return ComponentSpaceObject<T, _Space> representing the c-th component
     */
    template<size_t O = _Order>
    typename std::enable_if<O == 1, ComponentSpaceObject<T, _Space>>::type
    operator[](index_t c) const
    {
        return ComponentSpaceObject<T, _Space>(*this, c);
    }

    /// @brief Simple setup without boundary conditions
    void setup(const index_t _icont = -1) const
    {
        GISMO_ASSERT(m_sd != nullptr, "SpaceData not set. Call setSpaceData first.");
        this->setInterfaceCont(_icont);
        m_sd->mapper = gsDofMapper();

        if (const gsMultiBasis<T> * mb =
            dynamic_cast<const gsMultiBasis<T>*>(&this->source()) )
        {
            m_sd->mapper = gsDofMapper(*mb, this->spaceData()->dim );
            if ( 0 == _icont ) // Conforming boundaries ?
            {
                for ( gsBoxTopology::const_iiterator it = mb->topology().iBegin();
                      it != mb->topology().iEnd(); ++it )
                {
                    mb->matchInterface(*it, m_sd->mapper);
                }
            }
        }
        else if (const gsMappedBasis<2,T> * mb =
            dynamic_cast<const gsMappedBasis<2,T>*>(&this->source()) )
        {
            m_sd->mapper.setIdentity(mb->nPatches(), mb->size() , this->spaceData()->dim);
        }
        
        m_sd->mapper.finalize();
        m_sd->fixedDofs.setZero(m_sd->mapper.boundarySize(), 1);
        m_sd->initialized = true;
    }
    
    /// @brief Setup with boundary conditions (eliminating Dirichlet DOFs)
    void setup(const gsBoundaryConditions<T> & bc, const index_t dir_values,
               const index_t _icont = -1) const
    {
        GISMO_ASSERT(m_sd != nullptr, "SpaceData not set. Call setSpaceData first.");
        this->setInterfaceCont(_icont);
        m_sd->mapper = gsDofMapper();
        const index_t dim = this->dim();
        
        const gsMultiBasis<T> *mb = dynamic_cast<const gsMultiBasis<T> *>(&this->source());
        if (mb != nullptr)
        {
            m_sd->mapper = gsDofMapper(*mb, dim);
            if (0 == _icont) // Conforming boundaries ?
            {
                for (gsBoxTopology::const_iiterator it = mb->topology().iBegin();
                     it != mb->topology().iEnd(); ++it) {
                    if ( it->type() != interaction::contact )
                        mb->matchInterface(*it, m_sd->mapper);
                }
            }

            // Strong Dirichlet conditions - mark boundary DOFs for elimination
            gsMatrix<index_t> bnd;
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
            {
                const index_t cc = it->unkComponent();
                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mb->basis(it->ps.patch).boundary(it->ps.side());
                m_sd->mapper.markBoundary(it->ps.patch, bnd, cc);
            }
            
            // Clamped boundary condition (per DoF)
            gsMatrix<index_t> bnd1;
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Clamped"); it != bc.end("Clamped"); ++it)
            {
                const index_t cc = it->unkComponent();

                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mb->basis(it->ps.patch).boundaryOffset(it->ps.side(), 0);
                bnd1 = mb->basis(it->ps.patch).boundaryOffset(it->ps.side(), 1);
                if (!it->ps.parameter())
                        bnd.swap(bnd1);
                for (index_t c = 0; c!=dim; c++)
                {
                    if (c==cc || cc==-1 )
                        for (index_t k = 0; k < bnd.size(); ++k)
                            m_sd->mapper.matchDof(it->ps.patch, (bnd)(k, 0),
                                                  it->ps.patch, (bnd1)(k, 0), c);
                }
            }

            // Collapsed
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Collapsed"); it != bc.end("Collapsed"); ++it)
            {
                const index_t cc = it->unkComponent();

                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mb->basis(it->ps.patch).boundary(it->ps.side());

                for (index_t c = 0; c!=dim; c++)
                {
                    if (c==cc || cc==-1)
                        for (index_t k = 0; k < bnd.size() - 1; ++k)
                            m_sd->mapper.matchDof(it->ps.patch, (bnd)(0, 0),
                                                  it->ps.patch, (bnd)(k + 1, 0), c);
                }
            }

            // corners
            for (typename gsBoundaryConditions<T>::const_citerator
                     it = bc.cornerBegin(); it != bc.cornerEnd(); ++it)
            {
                for (index_t r = 0; r!=dim; ++r)
                {
                    if (it->component!=-1 && r!=it->component) continue;

                    GISMO_ASSERT(static_cast<size_t>(it->patch) < mb->nBases(),
                                 "Problem: a corner boundary condition is set on a patch id which does not exist.");
                    m_sd->mapper.eliminateDof(mb->basis(it->patch).functionAtCorner(it->corner),
                                              it->patch, it->component);
                }
            }
        }
        else if (const gsBasis<T> *b = dynamic_cast<const gsBasis<T> *>(&this->source()))
        {
            m_sd->mapper = gsDofMapper(*b, dim);
            gsMatrix<index_t> bnd;
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it) {
                GISMO_ASSERT(it->ps.patch == 0,
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = b->boundary(it->ps.side());
                m_sd->mapper.markBoundary(0, bnd, it->unkComponent());
            }
        }
        else if (const gsMappedBasis<2, T> *mapb =
                   dynamic_cast<const gsMappedBasis<2, T> *>(&this->source()))
        {
            m_sd->mapper.setIdentity(mapb->nPatches(), mapb->size(), dim);

            if (0 == _icont)
            {
                gsMatrix<index_t> int1, int2;
                for (gsBoxTopology::const_iiterator it = mapb->getTopol().iBegin();
                     it != mapb->getTopol().iEnd(); ++it) {
                    int1 = mapb->basis(it->first().patch).boundaryOffset(it->first().side(), 0);
                    int2 = mapb->basis(it->second().patch).boundaryOffset(it->second().side(), 0);
                    m_sd->mapper.matchDofs(it->first().patch, int1, it->second().patch, int2);
                }
            }

            gsMatrix<index_t> bnd;
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it) {
                const index_t cc = it->unkComponent();
                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mapb->basis(it->ps.patch).boundary(it->ps.side());
                m_sd->mapper.markBoundary(it->ps.patch, bnd, cc);
            }

            for (typename gsBoundaryConditions<T>::const_citerator
                     it = bc.cornerBegin(); it != bc.cornerEnd(); ++it)
            {
                GISMO_ASSERT(it->patch < mapb->nPieces(),
                             "Problem: a corner boundary condition is set on a patch id which does not exist.");
                m_sd->mapper.eliminateDof(mapb->basis(it->patch).functionAtCorner(it->corner), it->patch, it->component);
            }
        }
        else
        {
            GISMO_ASSERT(0 == bc.size(), "Problem: BCs are ignored.");
            m_sd->mapper.setIdentity(this->source().nPieces(), this->source().size());
        }

        m_sd->mapper.finalize();
        m_sd->fixedDofs.setZero(m_sd->mapper.boundarySize(), 1);
        m_sd->initialized = true;

        // Compute Dirichlet node values
        computeDirichletValues(bc, dir_values);
    }
    
    /// @brief Compute Dirichlet DOF values using the specified method
    void computeDirichletValues(const gsBoundaryConditions<T> & bc, 
                                const index_t dir_values) const
    {
        if (bc.container("Dirichlet").empty() && bc.cornerValues().empty()) return;
        
        const gsDofMapper & mapper = this->mapper();
        gsMatrix<T> & fixedDofs = m_sd->fixedDofs;
        fixedDofs.setZero(mapper.boundarySize(), 1);
        
        switch (dir_values)
        {
        case dirichlet::homogeneous:
        case dirichlet::user:
            // homogeneous: fill with zeros (already done)
            break;
        case dirichlet::interpolation:
            computeDirichletByInterpolation(bc);
            break;
        case dirichlet::l2Projection:
            computeDirichletByL2Projection(bc);
            break;
        default:
            GISMO_ERROR("Unknown Dirichlet value computation method: " << dir_values);
        }
        
        // Corner values
        for (typename gsBoundaryConditions<T>::const_citerator it = bc.cornerBegin(); 
             it != bc.cornerEnd(); ++it)
        {
            if (it->unknown != static_cast<index_t>(this->id())) continue;
            
            const int k = it->patch;
            const gsBasis<T> & basis = this->source().basis(k);
            const int i = basis.functionAtCorner(it->corner);
            const index_t com = it->component;
            
            for (index_t r = 0; r != this->dim(); ++r)
            {
                if (com != -1 && r != com) continue;
                const int ii = mapper.bindex(i, k, r);
                fixedDofs.at(ii) = it->value;
            }
        }
    }

private:
    /// @brief Compute Dirichlet values by interpolation
    void computeDirichletByInterpolation(const gsBoundaryConditions<T> & bc) const
    {
        const index_t parDim = this->source().domainDim();
        std::vector<gsVector<T>> rr;
        gsMatrix<index_t> boundary;
        gsVector<T> b(1);
        gsMatrix<T> fpts, tmp;
        
        gsMatrix<T> & fixedDofs = m_sd->fixedDofs;
        
        for (typename gsBoundaryConditions<T>::const_iterator it = bc.begin("Dirichlet");
             it != bc.end("Dirichlet"); ++it)
        {
            if (it->unknown() != static_cast<index_t>(this->id())) continue;
            
            const index_t com = it->unkComponent();
            const int k = it->patch();
            const gsBasis<T> & basis = this->source().basis(k);
            
            // Get boundary DOF indices
            boundary = basis.boundary(it->ps.side());
            
            // Get the side information
            short_t dir = it->ps.direction();
            index_t param = (it->ps.parameter() ? 1 : 0);
            
            // Build full parametric grid on the face
            rr.clear();
            rr.reserve(parDim);
            for (index_t i = 0; i < parDim; ++i)
            {
                if (i == dir)
                {
                    // Fixed direction: use the boundary parameter value
                    gsVector<T> b_vec(1);
                    b_vec[0] = basis.component(i).support()(0, param);
                    rr.push_back(b_vec);
                }
                else
                {
                    // Free direction: use anchors
                    rr.push_back(basis.component(i).anchors().transpose());
                }
            }
            
            // Create parametric point grid
            gsMatrix<T> paramPts = gsPointGrid<T>(rr);
            GISMO_ASSERT(paramPts.cols() == boundary.size(),
                         "Boundary size mismatch: " << paramPts.cols() << " vs " << boundary.size());
            
            // Evaluate Dirichlet function - map to physical coordinates if needed
            if (it->parametric())
            {
                // Function is defined in parametric coordinates
                it->function()->eval_into(paramPts, tmp);
            }
            else
            {
                // Function is defined in physical coordinates - need geometry map
                gsMatrix<T> physPts;
                bc.geoMap().piece(k).eval_into(paramPts, physPts);
                it->function()->eval_into(physPts, tmp);
            }
            
            for (index_t r = 0; r != this->dim(); ++r)
            {
                if (com != -1 && r != com) continue;
                
                for (index_t l = 0; l != boundary.size(); ++l)
                {
                    const index_t ii = this->mapper().bindex(boundary.at(l), k, r);
                    fixedDofs.at(ii) = tmp(r, l);
                }
            }
        }
    }
    
    /// @brief Compute Dirichlet values by L2 projection
    void computeDirichletByL2Projection(const gsBoundaryConditions<T> & bc) const
    {
        // For simplicity, fallback to interpolation
        // A full implementation would solve a boundary L2 projection problem
        gsWarn << "L2 projection for Dirichlet values not fully implemented, using interpolation.\n";
        computeDirichletByInterpolation(bc);
    }

public:

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        helper.add(*this);
        m_fd->flags |= NEED_VALUE | NEED_ACTIVE;
        if (Deriv_ > 0)
            m_fd->flags |= NEED_DERIV;
        if (Deriv_ > 1)
            m_fd->flags |= NEED_DERIV2;
    }

    void print(std::ostream & os) const
    {
        os<<Base::label_;
        _print_arguments(os);
    }

    const SpaceObject<T, SpaceType::Trial, _Order> & trial() const
    {
        return trial_impl<_Space>();
    }

    const SpaceObject<T, SpaceType::Test, _Order> & test () const
    {
        return test_impl<_Space>();
    }

    template <size_t S = _Space>
    typename std::enable_if<S==SpaceType::Trial, const SpaceObject<T, SpaceType::Trial, _Order> &>::type
    trial_impl() const
    {
        return *this;
    }
    template <size_t S = _Space>
    typename std::enable_if<S!=SpaceType::Trial, const SpaceObject<T, SpaceType::Trial, _Order> &>::type
    trial_impl() const
    {
        return NullObject<T,SpaceType::Trial,_Order>::get();
    }
    
    template <size_t S = _Space>
    typename std::enable_if<S==SpaceType::Test, const SpaceObject<T, SpaceType::Test, _Order> &>::type
    test_impl () const
    {
        return *this;
    }
    template <size_t S = _Space>
    typename std::enable_if<S!=SpaceType::Test, const SpaceObject<T, SpaceType::Test, _Order> &>::type
    test_impl () const
    {
        return NullObject<T,SpaceType::Test,_Order>::get();
    }

protected:
    void _print_arguments(std::ostream & os) const
    {
        os<<"(";
        for (size_t d=0; d!=this->domainDim(); d++)
        {
            os<<"x"<<d;
            if (d!=this->domainDim()-1)
                os<<",";
        }
        os<<")";
    }
};

template <class _LhsScalar, enum SpaceType _LhsSpace, size_t _LhsOrder,
          class _RhsScalar, enum SpaceType _RhsSpace, size_t _RhsOrder> 
inline bool operator==(const SpaceObject<_LhsScalar, _LhsSpace, _LhsOrder>& lhs,
                          const SpaceObject<_RhsScalar, _RhsSpace, _RhsOrder>& rhs)
{
    //TODO add isAcross
    return (&lhs.source() == &rhs.source()) && (lhs.id() == rhs.id()); 
}

}//namespace Expr
}//namespace gismo

