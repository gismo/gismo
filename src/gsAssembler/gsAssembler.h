/** @file gsAssembler.h

    @brief Provides generic assembler routines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris, J. Sogn
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsBasisRefs.h>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsStdVectorRef.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsAffineFunction.h>

#include <gsIO/gsOptionList.h>

#include <gsPde/gsPde.h>
#include <gsPde/gsBoundaryConditions.h>

#include <gsAssembler/gsQuadRule.h>
#include <gsAssembler/gsSparseSystem.h>
#include <gsAssembler/gsRemapInterface.h>



namespace gismo
{

template <class T>
void transformGradients(const gsMapData<T> & md, index_t k, const gsMatrix<T>& allGrads, gsMatrix<T>& trfGradsK)
{
    GISMO_ASSERT(allGrads.rows() % md.dim.first == 0, "Invalid size of gradient matrix");

    const index_t numGrads = allGrads.rows() / md.dim.first;
    const gsAsConstMatrix<T> grads_k(allGrads.col(k).data(), md.dim.first, numGrads);
    trfGradsK.noalias() = md.jacobian(k).cramerInverse().transpose() * grads_k;
}

template <class T>
void transformLaplaceHgrad( const gsMapData<T> & md, index_t k,
                        const gsMatrix<T> & allGrads,
                        const gsMatrix<T> & allHessians,
                        gsMatrix<T> & result)
{
    //todo: check me
    gsMatrix<T> hessians;
    transformDeriv2Hgrad(md, k, allGrads, allHessians, hessians);
    result = hessians.leftCols(md.dim.first).rowwise().sum(); // trace of Hessian
    result.transposeInPlace();
}

template <class T>
void normal(const gsMapData<T> & md, index_t k, gsVector<T> & result)
{
    GISMO_ASSERT( md.dim.first+1 == md.dim.second, "Codimension should be equal to one");
    result.resize(md.dim.first+1);
    const gsMatrix<T> Jk = md.jacobian(k);

    T alt_sgn(1.0);
    typename gsMatrix<T>::RowMinorMatrixType minor;
    for (short_t i = 0; i <= md.dim.first; ++i) // for all components of the normal vector
    {
        Jk.rowMinor(i, minor);
        result[i] = alt_sgn * minor.determinant();
        alt_sgn = -alt_sgn;
    }
}

template <class T>
void outerNormal(const gsMapData<T> & md, index_t k, boxSide s, gsVector<T> & result)
{
    //todo: fix and check me
    short_t m_orientation = md.jacobian(k).determinant() >= 0 ? 1 : -1;

    const T sgn = sideOrientation(s) * m_orientation; // TODO: fix me
    const short_t dir = s.direction();

    // assumes points u on boundary "s"
    result.resize(md.dim.second);
    if (md.dim.first + 1 == md.dim.second) // surface case GeoDim == 3
    {
        const gsMatrix<T,3,1> Jk = md.jacobian(k).col(!dir);
        // fixme: generalize to nD
        normal(md, k, result);
        result = result.template head<3>().normalized().cross(sgn * Jk);

        /*
          gsDebugVar(result.transpose()); // result 1
          normal(k,result);
          Jk.col(dir) = result.normalized();
          gsMatrix<T, ParDim, md.dim.first> minor;
          T alt_sgn = sgn;
          for (short_t i = 0; i != GeoDim; ++i) // for all components of the normal
          {
          Jk.rowMinor(i, minor);
          result[i] = alt_sgn * minor.determinant();
          alt_sgn = -alt_sgn;
          }
          gsDebugVar(result.transpose()); // result 2
        //*/
    }
    else // planar case
    {
        GISMO_ASSERT(md.dim.first == md.dim.second, "Codim different than zero/one");

        if (1 == md.dim.second)
        {
            result[0] = sgn;
            return;
        } // 1D case

        const gsMatrix<T> Jk = md.jacobian(k);

        T alt_sgn = sgn;
        typename gsMatrix<T>::FirstMinorMatrixType minor;
        for (short_t i = 0; i != md.dim.first; ++i) // for all components of the normal
        {
            Jk.firstMinor(i, dir, minor);
            result[i] = alt_sgn * minor.determinant();
            alt_sgn = -alt_sgn;
        }
    }
}

template<typename T>
void secDerToHessian(typename gsMatrix<T>::constRef & secDers,
                     gsMatrix<T> & hessian,
                     short_t parDim)
{
    switch (parDim)
    {
        case 1:
            hessian.resize(1, 1);
            hessian(0, 0) = secDers(0, 0);
            break;
        case 2:
            hessian.resize(2, 2);
            hessian(0, 0) = secDers(0, 0);
            hessian(1, 1) = secDers(1, 0);
            hessian(0, 1) =
            hessian(1, 0) = secDers(2, 0);
            break;
        case 3:
            hessian.resize(3, 3);
            hessian(0, 0) = secDers(0, 0);
            hessian(1, 1) = secDers(1, 0);
            hessian(2, 2) = secDers(2, 0);
            hessian(0, 1) =
            hessian(1, 0) = secDers(3, 0);
            hessian(0, 2) =
            hessian(2, 0) = secDers(4, 0);
            hessian(1, 2) =
            hessian(2, 1) = secDers(5, 0);
            break;
        default:
            GISMO_ERROR("Parametric dimension 1, 2 or 3 was expected.");
    }
}

template<typename T>
void hessianToSecDer (const gsMatrix<T> & hessian,
                      typename gsMatrix<T>::Row secDers,
                      short_t parDim)
{
    switch (parDim)
    {
        case 1:
            secDers(0, 0) = hessian(0, 0);
            break;
        case 2:
            secDers(0, 0) = hessian(0, 0);
            secDers(0, 1) = hessian(1, 1);
            secDers(0, 2) = (hessian(1, 0) + hessian(0, 1)) / 2.0;
            break;
        case 3:
            secDers(0, 0) = hessian(0, 0);
            secDers(0, 1) = hessian(1, 1);
            secDers(0, 2) = hessian(2, 2);
            secDers(0, 3) = (hessian(0, 1) + hessian(1, 0)) / 2.0;
            secDers(0, 4) = (hessian(0, 2) + hessian(2, 0)) / 2.0;
            secDers(0, 5) = (hessian(1, 2) + hessian(2, 1)) / 2.0;
            break;
        default:
            GISMO_ERROR("Parametric dimension 1, 2 or 3 was expected.");
    }
}

template<typename T>
void secDerToTensor(typename Eigen::DenseBase<Eigen::Map<const Eigen::Matrix<T, -1, -1>, 0, Eigen::Stride<0, 0> > >::ConstColXpr & secDers,
                    gsMatrix<T> * a,
                    short_t parDim, short_t geoDim)
{
    const index_t dim = parDim * (parDim + 1) / 2;
    for (short_t i = 0; i < geoDim; ++i)
        secDerToHessian<T>(secDers.segment(i * dim, dim), a[i], parDim);
}

template <class T>
void transformDeriv2Hgrad(const gsMapData<T> & md,
                          index_t              k,
                          const gsMatrix<T> &  funcGrad,
                          const gsMatrix<T> &  funcSecDir,
                          gsMatrix<T> &        result)
{
    //todo: check me
    const short_t ParDim = md.dim.first;
    const short_t GeoDim = md.dim.second;
    GISMO_ASSERT(
        (ParDim == 1 && (GeoDim == 1 || GeoDim == 2 || GeoDim == 3))
        || (ParDim == 2 && (GeoDim == 2 || GeoDim == 3))
        || (ParDim == 3 && GeoDim == 3), "No implementation for this case");

    // important sizes
    const index_t parSecDirSize = ParDim * (ParDim + 1) / 2;
    const index_t fisSecDirSize = GeoDim * (GeoDim + 1) / 2;

    // allgrads
    const index_t numGrads = funcGrad.rows() / ParDim;

    result.resize(numGrads, fisSecDirSize);

    gsMatrix<T> JMT = md.jacobian(k).cramerInverse();  // m_jacInvs.template block<GeoDim,ParDim>(0, k*ParDim).transpose();
    typename gsMatrix<T>::Tr JM1 = JMT.transpose();           // m_jacInvs.template block<GeoDim,ParDim>(0, k*ParDim);

    // First part: J^-T H(u) J^-1
    gsMatrix<T> parFuncHessian;
    for (index_t i = 0; i < numGrads; ++i)
    {
        secDerToHessian<T>(funcSecDir.block(i * parSecDirSize, k, parSecDirSize, 1), parFuncHessian, ParDim);
        hessianToSecDer<T>(JM1 * parFuncHessian * JMT, result.row(i), GeoDim);
    }

    // Second part: sum_i[  J^-T H(G_i) J^-1 ( e_i^T J^-T grad(u) ) ]
    const gsAsConstMatrix<T, -1, -1>  & secDer = md.deriv2(k);
    std::vector<gsMatrix<T> > DDG(GeoDim);
    secDerToTensor<T>(secDer.col(0), DDG.data(), ParDim, GeoDim);
    gsMatrix<T> HGT(GeoDim, fisSecDirSize);
    for (short_t i = 0; i < GeoDim; ++i)
        hessianToSecDer<T>(JM1 * DDG[i] * JMT, HGT.row(i), GeoDim);

    // Lastpart: substract part2 from part1
    const gsAsConstMatrix<T> grads_k(funcGrad.col(k).data(), ParDim, numGrads);
    result.noalias() -= grads_k.transpose() * JMT * HGT;
    // 1 x d * d x d * d x d * d * s -> 1 x s
}

/** @brief The assembler class provides generic routines for volume
  and boundary integrals that are used for for matrix and right-hand
  side generation

  \ingroup Assembler
*/
template <class T>
class gsAssembler
{
private:
    typedef gsStdVectorRef<gsDofMapper> gsDofMappers;

    typedef typename gsBoundaryConditions<T>::bcContainer bcContainer;

protected: // *** Input data members ***

    /// The PDE: contains multi-patch domain, boundary conditions and
    /// coeffcient functions
    memory::shared_ptr<gsPde<T> > m_pde_ptr;

    /// The discretization bases corresponding to the patches and to
    /// the number of solution fields that are to be computed. One
    /// entry of the vector corresponds to one or more unknown (in case
    /// of a system of PDEs)
    std::vector< gsMultiBasis<T> > m_bases;

    /// Options
    gsOptionList m_options;

protected: // *** Output data members ***

    /// Global sparse linear system
    gsSparseSystem<T> m_system;

    /// Fixed DoF values (if applicable, for instance eliminated Dirichlet DoFs). One
    /// entry of the vector corresponds to exactly one unknown, i.e. the size of m_ddof
    /// must fit m_system.colBlocks().
    std::vector<gsMatrix<T> > m_ddof;

public:

    gsAssembler() : m_options(defaultOptions())
    { }

    virtual ~gsAssembler()
    { }

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    /// @brief Create an empty Assembler of the derived type and return a
    /// pointer to it. Call the initialize functions to set the members
    virtual gsAssembler * create() const;

    /// @brief Clone this Assembler, making a deep copy.
    virtual gsAssembler * clone() const;

    /// @brief Intitialize function for, sets data fields
    /// using the pde, a vector of multi-basis and assembler options.
    void initialize(const gsPde<T>                         & pde,
                    const gsStdVectorRef<gsMultiBasis<T> > & bases,
                    //note: merge with initialize(.., const gsMultiBasis<T> ,..) ?
                    const gsOptionList & opt = defaultOptions() )
    {
        typename gsPde<T>::Ptr _pde = memory::make_shared_not_owned(&pde);
        initialize(_pde,bases,opt);
    }

    /// @brief Intitialize function for, sets data fields
    /// using the pde, a vector of multi-basis and assembler options.
    void initialize(typename gsPde<T>::Ptr pde,
                    const gsStdVectorRef<gsMultiBasis<T> > & bases,
                    //note: merge with initialize(.., const gsMultiBasis<T> ,..) ?
                    const gsOptionList & opt = defaultOptions() )
    {
        m_pde_ptr = pde;
        m_bases = bases;
        m_options = opt;
        refresh(); // virtual call to derived
        GISMO_ASSERT( check(), "Something went wrong in assembler initialization");
    }

    /// @brief Intitialize function for, sets data fields
    /// using the pde, a multi-basis and assembler options.
    void initialize(const gsPde<T>           & pde,
                    const gsMultiBasis<T>    & bases,
                    const gsOptionList & opt = defaultOptions() )
    {
        typename gsPde<T>::Ptr _pde = memory::make_shared_not_owned(&pde);
        initialize(_pde,bases,opt);
    }

    void initialize(typename gsPde<T>::Ptr pde,
                    const gsMultiBasis<T>    & bases,
                    const gsOptionList & opt = defaultOptions() )
    {
        m_pde_ptr = pde;
        m_bases.clear();
        m_bases.push_back(bases);
        m_options = opt;
        refresh(); // virtual call to derived
        GISMO_ASSERT( check(), "Something went wrong in assembler initialization");
    }


    /// @brief Intitialize function for, sets data fields
    /// using the pde, a vector of bases and assembler options.
    void initialize(const gsPde<T>           & pde,
                    const gsBasisRefs<T>     & basis,
                    const gsOptionList & opt = defaultOptions() )
    {
        GISMO_ASSERT(pde.domain().nPatches() ==1,
                     "You cannot initialize a multipatch domain with bases on a single patch");
        m_pde_ptr = memory::make_shared_not_owned(&pde);

        m_bases.clear();
        m_bases.reserve(basis.size());
        for(size_t c=0;c<basis.size();c++)
            m_bases.push_back(gsMultiBasis<T>(basis[c]));

        m_options = opt;
        refresh(); // virtual call to derived
        GISMO_ASSERT( check(), "Something went wrong in assembler initialization");
    }

    /// @brief checks for consistency and legal values of the stored members.
    bool check();

    /// @brief finishes the assembling of the system matrix,
    /// i.e. calls its .makeCompressed() method.
    void finalize() { m_system.matrix().makeCompressed(); }

    /** @brief Penalty constant for patch \a k, used for Nitsche and
    / Discontinuous Galerkin methods
    */
    T penalty(index_t k) const
    {
        const short_t deg = m_bases[0][k].maxDegree();
        return (deg + m_bases[0][k].dim()) * (deg + 1) * T(2.0);
    }

    /// @brief Provides an estimation of the number of non-zero matrix
    /// entries per column. This value can be used for sparse matrix
    /// memory allocation
    index_t numColNz() const
    {
        return m_system.numColNz(m_bases.front()[0], m_options);
    }

public:  /* Virtual assembly routines*/

    /// @brief Creates the mappers and setups the sparse system.
    /// to be implemented in derived classes, see scalarProblemGalerkinRefresh()
    /// for a possible implementation
    virtual void refresh();

    /// @brief Main assemble routine, to be implemented in derived classes
    virtual void assemble();

    /// @brief Main non-linear assemble routine with input from
    /// current solution
    virtual void assemble(const gsMultiPatch<T> & curSolution);

    gsOptionList & options() {return m_options;}

public: /* Element visitors */

    /// @brief Iterates over all elements of the domain and applies
    /// the \a ElementVisitor
    template<class ElementVisitor>
    void push()
    {
        for (size_t np = 0; np < m_pde_ptr->domain().nPatches(); ++np)
        {
            ElementVisitor visitor(*m_pde_ptr);
            //Assemble (fill m_matrix and m_rhs) on patch np
            apply(visitor, np);
        }
    }

    /// @brief Iterates over all elements of the boundaries \a BCs and
    /// applies the \a BElementVisitor
    template<class BElementVisitor>
    void push(const bcContainer & BCs)
    {
        for (typename bcContainer::const_iterator it
             = BCs.begin(); it!= BCs.end(); ++it)
        {
            BElementVisitor visitor(*m_pde_ptr, *it);
            //Assemble (fill m_matrix and m_rhs) contribution from this BC
            apply(visitor, it->patch(), it->side());
        }
    }

    /// @brief Iterates over all elements of the domain and applies
    /// the \a ElementVisitor
    template<class ElementVisitor>
    void push(const ElementVisitor & visitor)
    {
        for (size_t np = 0; np < m_pde_ptr->domain().nPatches(); ++np)
        {
            ElementVisitor curVisitor = visitor;
            //Assemble (fill m_matrix and m_rhs) on patch np
            apply(curVisitor, np);
        }
    }

    /// @brief Applies the \a BElementVisitor to the boundary condition \a BC
    template<class BElementVisitor>
    void push(const BElementVisitor & visitor, const boundary_condition<T> & BC)
    {
        BElementVisitor curVisitor = visitor;
        //Assemble (fill m_matrix and m_rhs) contribution from this BC
        apply(curVisitor, BC.patch(), BC.side());
    }

    /// @brief Iterates over all elements of interfaces and
    /// applies the \a InterfaceVisitor
    template<class InterfaceVisitor>
    void pushInterface()
    {
        InterfaceVisitor visitor(*m_pde_ptr);

        const gsMultiPatch<T> & mp = m_pde_ptr->domain();
        for ( typename gsMultiPatch<T>::const_iiterator
                  it = mp.iBegin(); it != mp.iEnd(); ++it )
        {
            const boundaryInterface & iFace = //recover master elemen
                ( m_bases[0][it->first() .patch].numElements(it->first() .side() ) <
                  m_bases[0][it->second().patch].numElements(it->second().side() ) ?
                  it->getInverse() : *it );

            this->apply(visitor, iFace);
        }
    }


public:  /* Dirichlet degrees of freedom computation */

    /// @brief Triggers computation of the Dirichlet dofs
    /// \param[in] unk the considered unknown
    void computeDirichletDofs(short_t unk = 0);

    /// @brief the user can manually set the dirichlet Dofs for a given patch and
    /// unknown, based on the Basis coefficients
    /// \param[in] coefMatrix the coefficients of the function
    /// \param[in] unk the consideren unknown
    /// \param[in] patch the patch index
    void setFixedDofs(const gsMatrix<T> & coefMatrix, short_t unk = 0, size_t patch = 0);

    /// @brief the user can manually set the dirichlet Dofs for a given patch and
    /// unknown.
    /// \param[in] vals the values of the eliminated dofs.
    /// \param[in] unk the considered unknown
    void setFixedDofVector(gsMatrix<T> vals, short_t unk = 0);

    /// Enforce Dirichlet boundary conditions by diagonal penalization
    /// \param[in] unk the considered unknown
    void penalizeDirichletDofs(short_t unk = 0);

    /// @brief Sets any Dirichlet values to homogeneous (if applicable)
    /// \param[in] unk the considered unknown
    void homogenizeFixedDofs(short_t unk = 0)
    {
        if(unk==-1)
        {
            for(size_t i=0;i<m_ddof.size();++i)
                m_ddof[i].setZero();
        }
        else
            m_ddof[unk].setZero();
    }

    // index_t numFixedDofs(short_t unk = 0) {return m_dofMappers[unk].boundarySize();}

    /// @brief Returns all the Dirichlet values (if applicable)
    const std::vector<gsMatrix<T> > & allFixedDofs() const { return m_ddof; }

    /// @brief Returns the Dirichlet values for a unknown (if applicable)
    /// \param[in] unk the considered unknown
    const gsMatrix<T> & fixedDofs(short_t unk=0) const { return m_ddof[unk]; }

    GISMO_DEPRECATED
    const gsMatrix<T> & dirValues(short_t unk=0) const { return m_ddof[unk]; }//remove

protected:  /* Helpers for Dirichlet degrees of freedom computation */

    /// @brief calculates the values of the eliminated dofs based on Interpolation.
    /// \param[in] mapper the dofMapper for the considered unknown
    /// \param[in] mbasis the multipabasis for the considered unknown
    /// \param[in] unk_ the considered unknown
    void computeDirichletDofsIntpl(const gsDofMapper     & mapper,
                                   const gsMultiBasis<T> & mbasis,
                                   const short_t unk_ = 0);

    /// @brief calculates the values of the eliminated dofs based on L2 Projection.
    /// \param[in] mapper the dofMapper for the considered unknown
    /// \param[in] mbasis the multipabasis for the considered unknown
    /// \param[in] unk_ the considered unknown
    void computeDirichletDofsL2Proj(const gsDofMapper     & mapper,
                                    const gsMultiBasis<T> & mbasis,
                                    const short_t unk_ = 0);

public:  /* Solution reconstruction */

    /// @brief Construct solution from computed solution vector for a single unknows
    /// \param[in] solVector the solution vector obtained from the linear system
    /// \param[out] result the solution in form of a gsMultiBasis
    /// \param[in] unk the considered unknown
    virtual void constructSolution(const gsMatrix<T>& solVector,
                                   gsMultiPatch<T>& result, short_t unk = 0) const;



    /// @brief Construct solution from computed solution vector for a set of unknowns.
    /// The result is a vectorfield, where each component is given the corresponding
    /// entry of \par unknowns. This method assumes that the specified unknowns have the
    /// same basis.
    /// \param[in] solVector the solution vector obtained from the linear system
    /// \param[out] result the solution seen as vectorfield in form of a gsMultiBasis
    /// \param[in] unknowns the considered vector of unknowns
    virtual void constructSolution(const gsMatrix<T>& solVector,
                                   gsMultiPatch<T>& result,
                                   const gsVector<index_t>  & unknowns) const;

    gsField<T> constructSolution(const gsMatrix<T>& solVector, short_t unk = 0) const;

    /// @brief Update solution by adding the computed solution vector
    /// to the current solution specified by \par result. This method assumes that all
    /// unknowns have the same basis.
    /// \param[in] solVector the solution vector obtained from the linear system
    /// \param[out] result the solution in form of a gsMultiBasis, \par solVector is added to the
    ///                   coefficients of result.
    /// \param[in] theta damping factor for update, theta = 1 corresponds to a full step.
    virtual void updateSolution(const gsMatrix<T>& solVector,
                                gsMultiPatch<T>& result, T theta = T(1)) const;

public: // *** Accessors ***

    /// @brief Return the Pde
    const gsPde<T> & pde() const { return *m_pde_ptr; }

    /// @brief Return the multipatch
    const gsMultiPatch<T> & patches() const { return m_pde_ptr->patches(); }

    /// @brief Return the multi-basis
    const gsMultiBasis<T> & multiBasis(index_t k = 0) const { return m_bases[k]; }

    /// @brief Return the multi-basis. Note: if the basis is altered,
    /// then refresh() should be called
    gsMultiBasis<T> & multiBasis(index_t k = 0) { return m_bases[k]; }

    /// @brief Returns the number of multi-bases
    size_t numMultiBasis() const {return m_bases.size(); }

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix() const { return m_system.matrix(); }

    /// @brief Returns the left-hand side vector(s)
    /// ( multiple right hand sides possible )
    const gsMatrix<T> & rhs() const { return m_system.rhs(); }
    gsMatrix<T> & rhs() { return m_system.rhs(); }

    /// @brief Returns the left-hand global matrix
    const gsSparseSystem<T> & system() const { return m_system; }
    gsSparseSystem<T> & system() { return m_system; }

    /// @brief Swaps the actual sparse system with the given one
    void setSparseSystem(gsSparseSystem<T> & sys)
    {
        GISMO_ASSERT(sys.initialized(), "Sparse system must be initialized first");
        m_system.swap(sys);
    }

    /// @brief Returns the number of (free) degrees of freedom
    index_t numDofs() const
    {
        index_t sum = 0;
        for (index_t c = 0; c!= m_system.numColBlocks(); ++c)
            sum += m_system.colMapper(c).freeSize();
        return sum;
    }

protected:

    /// @brief A prototype of the refresh function for a "standard"
    /// scalar problem.  Creats one mapper based on the set options
    /// and initializes the sparse system (without allocating memory.
    void scalarProblemGalerkinRefresh();

protected:

    /// @brief Generic assembly routine for volume or boundary integrals
    /// \param[in] visitor The visitor for the boundary or volume integral
    /// \param[in] patchIndex The considered patch
    /// \param[in] side The considered boundary side, only necessary for boundary
    /// integrals.
    template<class ElementVisitor>
    void apply(ElementVisitor & visitor,
               size_t patchIndex = 0,
               boxSide side = boundary::none);

    /// @brief Generic assembly routine for patch-interface integrals
    template<class InterfaceVisitor>
    void apply(InterfaceVisitor & visitor,
               const boundaryInterface & bi);
};

template <class T>
template<class ElementVisitor>
void gsAssembler<T>::apply(ElementVisitor & visitor,
                           size_t patchIndex,
                           boxSide side)
{
    //gsDebug<< "Apply to patch "<< patchIndex <<"("<< side <<")\n";

    const gsBasisRefs<T> bases(m_bases, patchIndex);

#pragma omp parallel
{
    gsQuadRule<T> quRule ; // Quadrature rule
    gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
    gsVector<T> quWeights; // Temp variable for mapped weights

    ElementVisitor
#ifdef _OPENMP
    // Create thread-private visitor
    visitor_(visitor);
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#else
    &visitor_ = visitor;
#endif

    // Initialize reference quadrature rule and visitor data
    visitor_.initialize(bases, patchIndex, m_options, quRule);

    const gsGeometry<T> & patch = m_pde_ptr->patches()[patchIndex];

    // Initialize domain element iterator -- using unknown 0
    typename gsBasis<T>::domainIter domIt = bases[0].makeDomainIterator(side);

    // Start iteration over elements
#ifdef _OPENMP
    for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
    for (; domIt->good(); domIt->next() )
#endif
    {
        // Map the Quadrature rule to the element
        quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

        // Perform required evaluations on the quadrature nodes
        visitor_.evaluate(bases, patch, quNodes);

        // Assemble on element
        visitor_.assemble(*domIt, quWeights);

        // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
        visitor_.localToGlobal(patchIndex, m_ddof, m_system); // omp_locks inside
    }
}//omp parallel

}


template <class T>
template<class InterfaceVisitor>
void gsAssembler<T>::apply(InterfaceVisitor & visitor,
                           const boundaryInterface & bi)
{
    gsRemapInterface<T> interfaceMap(m_pde_ptr->patches(), m_bases[0], bi);

    const index_t patchIndex1      = bi.first().patch;
    const index_t patchIndex2      = bi.second().patch;
    const gsBasis<T> & B1 = m_bases[0][patchIndex1];// (!) unknown 0
    const gsBasis<T> & B2 = m_bases[0][patchIndex2];

    gsQuadRule<T> quRule ; // Quadrature rule
    gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
    gsVector<T> quWeights;         // Mapped weights
    
    // Initialize
    visitor.initialize(B1, B2, bi, m_options, quRule);

    const gsGeometry<T> & patch1 = m_pde_ptr->patches()[patchIndex1];
    const gsGeometry<T> & patch2 = m_pde_ptr->patches()[patchIndex2];

    // Initialize domain element iterators
    typename gsBasis<T>::domainIter domIt = interfaceMap.makeDomainIterator();
    //int count = 0;

    // iterate over all boundary grid cells on the "left"
    for (; domIt->good(); domIt->next() )
    {
        //count++;

        // Compute the quadrature rule on both sides
        quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes1, quWeights);
        interfaceMap.eval_into(quNodes1,quNodes2);

        // Perform required evaluations on the quadrature nodes
        visitor.evaluate(B1, patch1, B2, patch2, quNodes1, quNodes2);

        // Assemble on element
        visitor.assemble(*domIt,*domIt, quWeights);

        // Push to global patch matrix (m_rhs is filled in place)
        visitor.localToGlobal(patchIndex1, patchIndex2, m_ddof, m_system);
    }

}


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAssembler.hpp)
#endif
