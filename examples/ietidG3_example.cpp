/** @file ietidG2_example.cpp

    @brief Provides an example for the ieti solver for a dG setting

    This class uses the gsPoissonAssembler. We use CG to solve for the
    Schur complement formulation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, R. Schneckenleitner
*/
#include <ctime>

#define DEBUGVAR(v) gsInfo << #v << ": " << (v) << " (line nr. " << __LINE__ << ")\n"
#define DEBUGMAT(m) gsInfo << #m << ": " << (m).rows() << "x" << (m).cols() << " (line nr. " << __LINE__ << ")\n"

#include <gismo.h>
#include <gsAssembler/gsVisitorDg.h>
#include <gsIeti/gsIetidGMapper.h>
#include <gsSolver/gsPreconditioner.h>
#include <unsupported/src/gsSolver/gsTimedOp.h>

using namespace gismo;

gsMultiPatch<real_t> approximateQuarterAnnulus(index_t deg)
{
    gsGeometry<>::uPtr quann = gsNurbsCreator<>::NurbsQuarterAnnulus();

    gsKnotVector<> KV1(0,1, 0, deg+1);        // no interior knots in x direction
    gsKnotVector<> KV2(0,1, 0, deg+1);        // no interior knot in y direction

    gsTensorBSplineBasis<2> tbsp (give(KV1), give(KV2));
    gsMatrix<real_t> eval = quann->eval(tbsp.anchors());
    gsGeometry<>::uPtr approxGeom = tbsp.interpolateAtAnchors( eval );
    gsMultiPatch<real_t> mp(*approxGeom);

    return mp;

}

/// struct to measure timings related to the setup of and solving with the IETI method
struct timings
{
    timings() {}

    timings(size_t nPatches): setupInverseInMsD(nPatches), setupPk(nPatches),
                              assemblePrimalBasis(nPatches), totalAssemble(1), totalSolving(1)
    {
        setupInverseInMsD.setZero(); setupPk.setZero(); assemblePrimalBasis.setZero(); totalAssemble.setZero(); totalSolving.setZero();
    }

    gsVector<real_t> setupInverseInMsD;
    //gsVector<real_t> applyMsD;
    gsVector<real_t> setupPk;
    //gsVector<real_t> PkApplication;
    //gsVector<real_t> PApplication;
    gsVector<real_t> assemblePrimalBasis;
    gsVector<real_t> totalAssemble;
    gsVector<real_t> totalSolving;


    void print(std::ostream& out)
    {
        std::ios  state(NULL);
        state.copyfmt(std::cout);

        //out<<std::setprecision(4)<<"\n\n";
        out<<std::setw(12)<<"Assembling time: "<<totalAssemble<<"\n";
        out<<std::setw(12)<<"Setup inverse in scaled Dirichlet Preconditioner: "<<setupInverseInMsD.sum()<<"\n";
        out<<std::setw(12)<<"Setup of P^(k): "<<setupPk.sum()<<"\n";
        out<<std::setw(12)<<"Setup of Psi: " << assemblePrimalBasis.sum() <<"\n";
        out<<std::setw(12)<<"Solving time: "<<totalSolving<<"\n";

        out<<"\n";
        std::cout.copyfmt(state);
    }
};

gsSparseMatrix<> mixedMass( const gsMultiBasis<>& mb, const gsMultiBasis<>&  mb2, const gsOptionList& options )
{
    gsExprAssembler<> mass(1, 1);

    if (mb.basis(0).numElements() < mb2.basis(0).numElements())
        mass.setIntegrationElements(mb2);
    else
        mass.setIntegrationElements(mb);

    gsExprAssembler<>::space u = mass.getSpace(mb2);
    gsExprAssembler<>::space v = mass.getTestSpace(u, mb);
    mass.setOptions(options);
    mass.initMatrix();
    mass.assemble(v * u.tr());

    gsSparseMatrix<> result;
    mass.matrix_into(result);
    return result;
}




/// A preconditioner class for IETI-DP with local inexact preconditioners
template<class T>
class gsInexactIETIPrec : public gsPreconditionerOp<T> {
private:
#if defined(GISMO_WITH_PARDISO)
    typedef typename gsSparseSolver<T>::PardisoLLT Choleskyfac;
#else
    typedef typename gsSparseSolver<T>::SimplicialLDLT Choleskyfac;
#endif

public:

    /// Shared pointer for gsLinearOperator
    typedef memory::shared_ptr <gsInexactIETIPrec> Ptr;

    /// Unique pointer for gsLinearOperator
    typedef memory::unique_ptr <gsInexactIETIPrec> uPtr;

    /// Base class
    typedef gsLinearOperator <T> Base;

    /// Base class
    typedef typename gsLinearOperator<T>::Ptr BasePtr;
    typedef typename gsLinearOperator<T>::uPtr BaseUPtr;

    //gsInexactIETIPrec() : m_originalSize() {}

    gsInexactIETIPrec(const gsSparseMatrix <T> &K,
                      const BasePtr firstprecond,
                      const gsMultiBasis<T> &basis,
                      const std::vector<gsIetidGMapper<>::ArtificialIface> &iFacePairs,
                      const gsDofMapper& augmentedMapper,
                      const T delta,
                      const gsOptionList &opt = gsAssembler<T>::defaultOptions()) :
    m_underlyingOperator(K), m_fastDiag(firstprecond), m_originalSize(m_fastDiag->rows()),
    m_basis(basis), //m_bc1(bc1),
    m_patchIFace(iFacePairs), m_mapper(augmentedMapper), //bc(Bc),
    m_options(opt)
    {
        m_R2T.resize(m_underlyingOperator.cols(), m_underlyingOperator.cols() - m_originalSize);

        gsSparseMatrix<> M;
        for (size_t i = 0; i < m_patchIFace.size(); ++i) {
            M = gsPatchPreconditionersCreator<>::massMatrix(* m_basis[i+1].boundaryBasis(m_patchIFace[i].takenFrom));
            m_mass.push_back(M);
            DEBUGVAR(i);
            DEBUGMAT(m_mass[i]);
            //m_mass.push_back(M.block(1,1,M.rows()-2, M.cols()-2)); //TODO: 3d; w/o corners
        }

        constructRestrictionMatricesforL2();
        index_t colentry = 0;
        for (index_t r = m_originalSize; r < m_underlyingOperator.rows(); r++, colentry++)
            m_R2T.insert(r, colentry) = T(1);
        m_R2T.makeCompressed();

        m_R1T *= 0.7;
        m_R2T *= 3.;

        // do the interface part
        index_t r = 0, c = 0;
        gsSparseEntries<> tripletList;
        for (std::vector<gsIetidGMapper<>::ArtificialIface >::const_iterator it = m_patchIFace.begin(); it != m_patchIFace.end(); ++it)
        {
            //index_t deg = math::max(m_basis[it-m_patchIFace.begin()+1].maxDegree(), m_basis[0].maxDegree());
            gsSparseMatrix <T> iMass = 0.5 * ( delta * (1./m_basis[it - m_patchIFace.begin()+1].getMinCellLength() + 1./m_basis[0].getMinCellLength()) ) * m_mass[it - m_patchIFace.begin()];
            for (index_t k = 0; k < iMass.cols(); ++k) {
                for (index_t i = 0; i < iMass.rows(); ++i)
                    tripletList.add(r + i, c + k, iMass(i, k));
            }
            r += iMass.rows();
            c += iMass.cols();
        }
        gsSparseMatrix <T> edgeMass(m_underlyingOperator.cols() - m_originalSize,
                                    m_underlyingOperator.cols() - m_originalSize);
        edgeMass.setFrom(tripletList);
        gsInfo << "m_underlyingOperator: " << m_underlyingOperator.rows() << "x" << m_underlyingOperator.cols() << "\n";
        gsInfo << "edgeMass: " << edgeMass.rows() << "x" << edgeMass.cols() << "\n";
        m_edgeSolver = memory::unique_ptr<Choleskyfac>(new Choleskyfac(edgeMass));

        m_eig = powerIteration();
        DEBUGVAR(m_eig);
    }

    static BasePtr make(const gsSparseMatrix <T> &K,
                        const BasePtr firstprecond,
                        const gsMultiBasis<T> &basis,
                        //const gsBoundaryConditions<T>& bc1,
                        const std::vector<gsIetidGMapper<>::ArtificialIface> &patchIFace,
                        const gsDofMapper& augmentedMapper,
                        //const std::vector<gsBasis < T> *>& iFaceBasis,
                        //const std::vector<gsBoundaryConditions<> > &Bc,
                        const T delta,
                        const gsOptionList &options = gsAssembler<T>::defaultOptions())
    {
        return BasePtr(new gsInexactIETIPrec(K, firstprecond, basis, patchIFace, augmentedMapper, delta, options));
    }

    /**
     * @brief Apply the method for given right hand side and current iterate
     * @param rhs Right hand side vector
     * @param x   Current iterate vector
     */
    void step(const gsMatrix <T> &rhs, gsMatrix <T> &x) const {
        gsMatrix <T> solExP1;
        ///just extracting the dofs
        m_fastDiag->apply(m_R1T.transpose() * rhs, solExP1);
        x.noalias() = (m_R1T * solExP1 + m_R2T * m_edgeSolver->solve(m_R2T.transpose() * rhs));
    }

    /**
     * @brief Apply the transposed variant of the method for given right hand
     *        side and current iterate
     * @param rhs Right hand side vector
     * @param x   Current iterate vector
     *
     * @warning Derived classes *must* overwrite this virtual function if the
     * preconditioner is not symmetric.
     */
    virtual void stepT(const gsMatrix <T> &rhs, gsMatrix <T> &x) const { step(rhs, x); } // Assume symmetry.

    void apply(const gsMatrix <T> &input, gsMatrix <T> &x) const {
        x.setZero(this->rows(), input.cols()); // we assume quadratic matrices
        step(input, x);
        for (index_t i = 0; i < 0; ++i)
        {
            //x *= (real_t)1./(m_eig);
            gsMatrix<T> update, residual;
            residual = input - m_underlyingOperator*x;
            step(residual,update);
            x += (real_t)0.5/m_eig*update; // TODO: Why only here?
        }
    }

    virtual index_t rows() const { return m_underlyingOperator.rows(); }

    virtual index_t cols() const { return m_underlyingOperator.cols(); }

    /// Return the underlying operator \f$ A \f$.
    virtual BasePtr underlyingOp() const {
        return gsMatrixOp < gsSparseMatrix < T > > ::make(m_underlyingOperator);
    }

    /// Set the number of sweeps to be applied in the member function apply
    void setNumOfSweeps(index_t n) {
        GISMO_ASSERT (n > 0, "Number of sweeps needs to be positive.");
        m_num_of_sweeps = n;
    }

    /// Get the number of sweeps to be applied in the member function \a apply
    index_t numOfSweeps() {
        return m_num_of_sweeps;
    }

    /// Get the default options as gsOptionList object
    static gsOptionList defaultOptions() {
        gsOptionList opt = Base::defaultOptions();
        opt.addInt("NumOfSweeps", "Number of sweeps to be applied in the member function apply", 1);
        return opt;
    }

    /// Set options based on a gsOptionList object
    virtual void setOptions(const gsOptionList &opt) {
        Base::setOptions(opt);
        m_num_of_sweeps = opt.askInt("NumOfSweeps", m_num_of_sweeps);
    }

private:

    // TODO: First only assemble matrices
    // Then, use dof mapper to handle boundary conditions
    // Then, construct the transfer matrix (should this be dense or sparse!?)

    void constructRestrictionMatricesforL2() {
        gsSparseEntries<T> tripletList;
        m_R1T.resize(m_underlyingOperator.rows(), m_fastDiag->rows());
        DEBUGMAT(m_R1T);
        for (int l = 0; l < m_fastDiag->rows(); ++l)
        {
            GISMO_ASSERT( l < m_R1T.rows(), l <<"<"<< m_R1T.rows() );
            GISMO_ASSERT( l < m_R1T.cols(), l <<"<"<< m_R1T.cols() );

            tripletList.add(l, l, T(1));
        }

        std::vector<index_t> dofsPatchIFace;
        std::vector<gsMatrix<T>> values;
        index_t r = 0;
        gsInfo << "m_patchIFace.size(): " << m_patchIFace.size() << "\n";
        for (size_t j = 0; j < m_patchIFace.size(); ++j) {
            gsInfo << "j: " << j << "\n";
            values.clear();
            dofsPatchIFace.clear();

            //for (index_t i = 1; i < m_basis[0].boundary(m_patchIFace[j].assignedTo).rows()-1; ++i) //without the corners //TODO: 3D
            for (index_t i = 0; i < m_basis[0].boundary(m_patchIFace[j].assignedTo).rows(); ++i)
                if (m_mapper.is_free(m_basis[0].boundary(m_patchIFace[j].assignedTo)(i, 0), 0))
                {
                    // subtract for the corner values
                    index_t idx = m_mapper.index(m_basis[0].boundary(m_patchIFace[j].assignedTo)(i, 0));
                    //for(size_t corner = 1; corner <= math::pow(2,m_basis[0].dim()); corner++) {
                    //    if (m_mapper.is_free(m_basis[0].functionAtCorner(corner), 0) && m_mapper.index(m_basis[0].functionAtCorner(corner)) < idx)
                    //    {
                    //        idx--;
                    //    }
                    //}
                    DEBUGVAR(idx);
                    dofsPatchIFace.push_back(idx);
                }
            gsMultiBasis<T> mb2(*m_basis[0].boundaryBasis(m_patchIFace[j].assignedTo));
            gsMultiBasis<T> mb(*m_basis[j+1].boundaryBasis(m_patchIFace[j].takenFrom));

            gsSparseMatrix<T> mixed_mass = mixedMass( mb, mb2, m_options );
            DEBUGMAT(mixed_mass);
            DEBUGVAR(mixed_mass.toDense());
            //TODO: eliminateDirichlet1D(boundaryConditionsForDirection(m_bc1, comp1), result);
            //result = result.block(1,1,result.rows()-2, result.cols()-2); //without the corners //TODO: 3D
            DEBUGMAT(m_mass[j]);

            gsMatrix<index_t> ai_idxes = m_basis[j+1].boundary(m_patchIFace[j].takenFrom);
            gsSparseEntries<> ai_embedding_se;
            index_t free_idx = 0;
            for (index_t i=0; i<ai_idxes.rows(); ++i)
            {
                if (m_mapper.is_free(ai_idxes(i,0),j+1))
                {
                    ai_embedding_se.add(i,free_idx,1);
                    ++free_idx;
                }
            }
            gsSparseMatrix<T> ai_embedding(ai_idxes.rows(), free_idx);
            ai_embedding.setFrom(ai_embedding_se);

            std::vector<gsSparseMatrix<T>> li_embeddings;
            gsMatrix<index_t> li_idxes = m_basis[0].boundary(m_patchIFace[j].assignedTo);
            gsSparseEntries<> li_embedding_se;
            free_idx = 0;
            for (index_t i=0; i<li_idxes.rows(); ++i)
            {
                if (m_mapper.is_free(li_idxes(i,0),0))
                {
                    li_embedding_se.add(i,free_idx,1);
                    ++free_idx;
                }
            }
            gsSparseMatrix<T> li_embedding(li_idxes.rows(), free_idx);
            li_embedding.setFrom(li_embedding_se);

            m_mass[j] = ai_embedding.transpose() * m_mass[j] * ai_embedding;
            DEBUGMAT(m_mass[j]);

            mixed_mass = ai_embedding.transpose() * mixed_mass * li_embedding;
            DEBUGMAT(mixed_mass);

            Choleskyfac solver;
            solver.compute(m_mass[j]);

            gsMatrix<T> projMatrix = solver.solve(mixed_mass.toDense());

            for (index_t k = 0; k < projMatrix.cols(); ++k)
            {
                for (index_t i = 0; i < projMatrix.rows(); ++i)
                {
                    index_t row = m_originalSize+r+i;
                    index_t col = dofsPatchIFace[k];
                    GISMO_ASSERT( row < m_R1T.rows(), row <<"<"<< m_R1T.rows() );
                    GISMO_ASSERT( col < m_R1T.cols(), col <<"<"<< m_R1T.cols() );
                    tripletList.add(row, col, projMatrix(i, k));
                }
            }
            r += projMatrix.rows();
        }

        m_R1T.setFrom(tripletList);
        m_R1T.makeCompressed();
        DEBUGVAR(m_R1T.toDense());
    }

    inline T powerIteration(index_t N = 10)
    {
        GISMO_ASSERT( m_underlyingOperator.rows() == m_underlyingOperator.cols(), "Dimensions do not argee");

        gsMatrix<T> x, y, z;
        z.setRandom(m_underlyingOperator.rows(),1);

        for (index_t i = 0; i < N; ++i)
        {
            z.swap(x);
            GISMO_ASSERT(x.rows()>0, "");
            const T v = x.col(0).dot(x.col(0));
            GISMO_ASSERT( v>=0, v<<"<0");
            x /= math::sqrt( v );
            y = m_underlyingOperator*x;
            step(y, z);
        }

        return y.col(0).dot(z.col(0)) / x.col(0).dot(y.col(0));

    }

    /*
    void eliminateDirichlet1D(const gsBoundaryConditions<T>& bc, gsSparseMatrix<T> & result) const
    {
        dirichlet::strategy ds = (dirichlet::strategy)m_options.askInt("DirichletStrategy",dirichlet::elimination);
        if (ds == dirichlet::elimination)
        {
            patchSide west(0,boundary::west), east(0,boundary::east);
            index_t i = 0;
            if (bc.getConditionFromSide( west ) && bc.getConditionFromSide( west )->type() == condition_type::dirichlet ) i += 1;
            if (bc.getConditionFromSide( east ) && bc.getConditionFromSide( east )->type() == condition_type::dirichlet ) i += 2;
            if (i%2 + i/2 >= result.rows() || i%2 + i/2 >= result.cols())
                result.resize(0,0);
            else switch ( i )
                {
                    case 0: break;
                    case 1: result = result.block( 1, 1, result.rows()-1, result.cols()-1 ); break;
                    case 2: result = result.block( 0, 0, result.rows()-1, result.cols()-1 ); break;
                    case 3: result = result.block( 1, 1, result.rows()-2, result.cols()-2 ); break;
                }

            for (typename gsBoundaryConditions<T>::const_citerator it = bc.cornerValues().begin(); it != bc.cornerValues().end(); ++it) {
                switch ( it->corner )
                {
                    case 1: result = result.block( 1, 1, result.rows()-1, result.cols()-1 ); break;
                    case 2: result = result.block( 0, 0, result.rows()-1, result.cols()-1 ); break;
                    case 3: result = result.block( 1, 1, result.rows()-1, result.cols()-1 ); break;
                    case 4: result = result.block( 0, 0, result.rows()-1, result.cols()-1 ); break;
                }
            }
        }
        else
            GISMO_ERROR("Unknown Dirichlet strategy.");
    }*/

protected:
    using gsPreconditionerOp<T>::m_num_of_sweeps;
    gsSparseMatrix <T> m_underlyingOperator;
    BasePtr m_fastDiag;
    index_t m_originalSize;

    const gsMultiBasis<T> &m_basis;

    //gsBoundaryConditions <T> m_bc1;
    const std::vector<gsIetidGMapper<>::ArtificialIface> m_patchIFace;
    const gsDofMapper& m_mapper;
    const std::vector<gsBasis < T>* > m_ArtiIBasis;
    //const std::vector<gsBoundaryConditions<> > bc;

    const gsOptionList &m_options;

    std::vector<gsSparseMatrix <T> > m_mass;

    gsSparseMatrix <T> m_R1T, m_R2T;
    memory::unique_ptr<Choleskyfac> m_edgeSolver;
    real_t m_eig;

}; // gsInexactIETIPrec

// Assemble routines for inexact dG
template<class T>
class ParameterdomainAssembler
{
public:
    static void apply(gsGenericAssembler<T>& ass,
               const gsMultiPatch<>& patches,
               const boundaryInterface& bi)
    {
        gsRemapInterface<T> interfaceMap(patches, ass.multiBasis(), bi);

        const index_t patchIndex1      = bi.first().patch;
        const index_t patchIndex2      = bi.second().patch;
        const gsBasis<T> & B1 = ass.multiBasis()[patchIndex1];// (!) unknown 0
        const gsBasis<T> & B2 = ass.multiBasis()[patchIndex2];

        gsQuadRule<T> quRule ; // Quadrature rule
        gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
        gsVector<T> quWeights;         // Mapped weights
        T pen;

        // Initialize
        initialize(B1, B2, bi.first(), ass.options(), quRule, pen);

        // Initialize domain element iterators
        typename gsBasis<T>::domainIter domIt = interfaceMap.makeDomainIterator();

        // iterate over all boundary grid cells on the "left"
        for (; domIt->good(); domIt->next() )
        {
            // Compute the quadrature rule on both sides
            quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes1, quWeights);
            interfaceMap.eval_into(quNodes1,quNodes2);

            // Perform required evaluations on the quadrature nodes
            evaluateParameterDomain(B1, B2, pen, quNodes1, quNodes2, quWeights, bi,ass);

        }

    }

    static void initialize(const gsBasis<T> & basis1,
                    const gsBasis<T> & basis2,
                    const patchSide & side1,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule,
                    T & penalty)
    {
        // Setup Quadrature
        rule = gsQuadrature::get(basis1, options, side1.direction());

        penalty     = options.askReal("DG.Penalty",-1);
        // If not given, use default
        if (penalty<0)
        {
            const index_t deg = math::max( basis1.maxDegree(), basis2.maxDegree() );
            penalty = T(4) * (deg + basis1.dim()) * (deg + 1);
        }

    }

    static void evaluateParameterDomain(
                         const gsBasis<T>       & B1,
                         const gsBasis<T>       & B2,
                         const T                & penalty,
                         const gsMatrix<T>      & quNodes1,
                         const gsMatrix<T>      & quNodes2,
                         const gsVector<T>      & quWeights,
                         const boundaryInterface        & bi,
                         gsGenericAssembler<T> & ass)
    {
        gsMapData<T> md1, md2;
        std::vector<gsMatrix<T> > basisData1, basisData2;
        gsMatrix<index_t> actives1, actives2;// active basis functions
        gsMatrix<T> E11, E12, E21, E22;// Matrices for the penalization term
        gsVector<T> unormal = gsVector<T>::Zero(B1.dim()); //unormal.setZero();

        md1.flags = md2.flags = NEED_VALUE|NEED_DERIV;

        md1.points = quNodes1;
        md2.points = quNodes2;
        // Compute the active basis functions
        B1.active_into(md1.points.col(0), actives1);
        B2.active_into(md2.points.col(0), actives2);
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();

        // Evaluate basis functions and their first derivatives
        B1.evalAllDers_into(md1.points, 1, basisData1);
        B2.evalAllDers_into(md2.points, 1, basisData2);

        // Initialize local matrices
        E11.setZero(numActive1, numActive1); E12.setZero(numActive1, numActive2);
        E22.setZero(numActive2, numActive2); E21.setZero(numActive2, numActive1);


        /// Assemble
        T h1    = B1.getMinCellLength();
        T h2    = B2.getMinCellLength();
        patchSide side1 = bi.first();
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector from patch1
            // Compute the outer normal vector of the parameter domain

            unormal( (side1-1)/2 ) = (side1%2) ? -1 : 1;

            /*switch (side1)
            {
                case 1: // west
                    unormal(0) = -1;
                    break;
                case 2: //east
                    unormal(0) = 1;
                    break;
                case 3: //south
                    unormal(1) = -1;
                    break;
                case 4: //north
                    unormal(1) = 1;
                    break;
                default:
                    GISMO_ASSERT(B1.dim()==2, "The dimension is expected to be 2!");
            }*/
            // Multiply quadrature weight by the geometry measure
            // Integral transformation and quadrature weight (patch1)
            // assumed the same on both sides
            const T weight = quWeights[k] * unormal.norm();

            // Take blocks of values and derivatives of basis functions
            const typename gsMatrix<T>::Block val1 = basisData1[0].block(0, k, numActive1, 1);
            const typename gsMatrix<T>::Block val2 = basisData2[0].block(0, k, numActive2, 1);

            // Compute element matrices
            const T c2 = weight * penalty * (1. / h1 + 1. / h2) * (T(1) / 2);

            E11.noalias() += c2 * (val1 * val1.transpose());
            E12.noalias() += c2 * (val1 * val2.transpose());
            E22.noalias() += c2 * (val2 * val2.transpose());
            E21.noalias() += c2 * (val2 * val1.transpose());
        }

        /// Adds the contributions to the sparse system
        // Map patch-local DoFs to global DoFs
        ass.system().mapColIndices(actives1, bi.first().patch, actives1);
        ass.system().mapColIndices(actives2, bi.second().patch, actives2);

        ass.system().pushToMatrix(E11, actives1,actives1,ass.fixedDofs(),0,0);
        ass.system().pushToMatrix(-E21, actives2,actives1,ass.fixedDofs(),0,0);
        ass.system().pushToMatrix(-E12, actives1,actives2,ass.fixedDofs(),0,0);
        ass.system().pushToMatrix(E22, actives2,actives2,ass.fixedDofs(),0,0);
    }
};




int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    std::string geometry("domain2d/square.xml");
    index_t splitPatches = 1;
    real_t stretchGeometry = 1;
    index_t refinements = 1;
    index_t degree = 2;
    bool nonMatching = false;
    real_t alpha = 1;
    real_t beta = 1;
    real_t penalty = 5;
    std::string boundaryConditions("d");
    std::string primals("c");
    bool eliminatePointwiseDofs = !true; // true for MFD
    real_t tolerance = 1.e-6;
    index_t maxIterations = 20000;
    std::string out;
    bool plot = false;

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using an isogeometric tearing and interconnecting (IETI) solver.");
    cmd.addString("g", "Geometry",              "Geometry file", geometry);
    cmd.addInt   ("",  "SplitPatches",          "Split every patch that many times in 2^d patches", splitPatches);
    cmd.addReal  ("",  "StretchGeometry",       "Stretch geometry in x-direction by the given factor", stretchGeometry);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addSwitch("",  "NonMatching",           "Set up a non-matching multi-patch discretization", nonMatching);
    cmd.addReal  ("",  "DG.Alpha",              "Parameter alpha for dG scheme; use 1 for SIPG and NIPG.", alpha );
    cmd.addReal  ("",  "DG.Beta",               "Parameter beta for dG scheme; use 1 for SIPG and -1 for NIPG", beta );
    cmd.addReal  ("",  "DG.Penalty",            "Penalty parameter delta for dG scheme; if negative, default 4(p+d)(p+1) is used.", penalty );
    cmd.addString("b", "BoundaryConditions",    "Boundary conditions", boundaryConditions);
    cmd.addString("c", "Primals",               "Primal constraints (c=corners, e=edges, f=faces)", primals);
    //cmd.addSwitch("e", "EliminateCorners",      "Eliminate corners (if they are primals)", eliminatePointwiseDofs);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Maximum number of iterations for linear solver", maxIterations);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //if ( ! gsFileManager::fileExists(geometry) )
    //{
    //    gsInfo << "Geometry file could not be found.\n";
    //    gsInfo << "I was searching in the current directory and in: " << gsFileManager::getSearchPaths() << "\n";
    //    return EXIT_FAILURE;
    //}

    gsInfo << "Run ietidG2_example with options:\n" << cmd << std::endl;

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;

    //! [Define Geometry]
    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    //! [Define Geometry]
    if (!mpPtr)
    {
        gsInfo << "No geometry found in file " << geometry << ".\n";
        return EXIT_FAILURE;
    }
    //! [Define Geometry2]
    gsMultiPatch<>& mp = *mpPtr;
    {

        /*std::vector<gsGeometry<>*> ptch;
        for(size_t np = 0; np< mp.nPatches(); ++np)
        {
            std::vector<gsGeometry<>* >  ptch_ =  mp.patch(np).uniformSplit(1);
            ptch.insert(ptch.end(),ptch_.begin(),ptch_.end());
        }
        mp = gsMultiPatch<real_t>(ptch);*/

        for (index_t i=0; i<splitPatches; ++i)
        {
            gsInfo << "split patches uniformly... " << std::flush;
            mp = mp.uniformSplit();
        }
        mp.computeTopology();
    }
    timings times(mp.nPatches());
    //! [Define Geometry2]

    if (stretchGeometry!=1)
    {
       gsInfo << "and stretch it... " << std::flush;
       for (size_t i=0; i!=mp.nPatches(); ++i)
           const_cast<gsGeometry<>&>(mp[i]).scale(stretchGeometry,0);
       // Const cast is allowed since the object itself is not const. Stretching the
       // overall domain keeps its topology.
    }

    gsInfo << "done.\n";

    /************** Define boundary conditions **************/

    gsInfo << "Define right-hand-side and boundary conditions... " << std::flush;

    //! [Define Source]
    // Right-hand-side
    gsFunctionExpr<> f( "2*pi^2*sin(pi*x)*sin(pi*y)", mp.geoDim() );

    // Dirichlet function
    gsFunctionExpr<> gD( "0.0", mp.geoDim() );

    // Neumann
    gsConstantFunction<> gN( 1.0, mp.geoDim() );

    gsBoundaryConditions<> bc;
    //! [Define Source]
    {
        const index_t len = boundaryConditions.length();
        index_t i = 0;
        for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
        {
            char b_local;
            if ( len == 1 )
                b_local = boundaryConditions[0];
            else if ( i < len )
                b_local = boundaryConditions[i];
            else
            {
                gsInfo << "\nNot enough boundary conditions given.\n";
                return EXIT_FAILURE;
            }

            if ( b_local == 'd' )
                bc.addCondition( *it, condition_type::dirichlet, &gD );
            else if ( b_local == 'n' )
                bc.addCondition( *it, condition_type::neumann, &gN );
            else
            {
                gsInfo << "\nInvalid boundary condition given; only 'd' (Dirichlet) and 'n' (Neumann) are supported.\n";
                return EXIT_FAILURE;
            }

            ++i;
        }
        if ( len > i )
            gsInfo << "\nToo many boundary conditions have been specified. Ignoring the remaining ones.\n";
        gsInfo << "done. "<<i<<" boundary conditions set.\n";
    }


    /************ Setup bases and adjust degree *************/

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    //! [Define Basis]
    gsMultiBasis<> mb(mp);
    //! [Define Basis]

    //! [Set degree and refine]
    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinements; ++i )
        mb.uniformRefine();

    // We might want to create a non-matching setup such that we are not in the
    // special case where the function spaces on the interfaces actually agree
    // (as this would be the case for the Yeti footprint, the default domain.)
    //gsPiecewiseFunction<> piece(mb.size());
    //gsFunctionExpr<> one("1.0", 2);
    //gsFunctionExpr<> two("2.0", 2);
    //gsFunctionExpr<> zero("0.0", 2);
    if (nonMatching)
    {
        gsInfo << "\n  Option NonMatching: Make uniform refinement for every third patch... " << std::flush;
        for (size_t i = 0; i < mb.nBases(); ++i)
            if ( i%3 == 0 ) {
                mb[i].uniformRefine();
                //piece.addPiece(zero);
            }
        gsInfo << "done.\n";
        --refinements;

        gsInfo << "  Option NonMatching: Increase spline degree for every other third patch... " << std::flush;
        for (size_t i = 0; i < mb.nBases(); ++i)
            if ( i%3 == 1 ) {
                mb[i].setDegreePreservingMultiplicity(degree+1);
                //piece.addPiece(one);
            }
        gsInfo << "done.\n";

        for (size_t i = 0; i < mb.nBases(); ++i)
            if ( i%3 == 2 ) {
                //piece.addPiece(two);
            }
    }
    //! [Set degree and refine]

    gsInfo << "done.\n";

    /************** Compute penalty parameter **************/
    if(penalty < 0)
        penalty *= mb.maxCwiseDegree();

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    const index_t nPatches = mp.nPatches();
    //! [Define global mapper]
    gsIetidGMapper<> ietiMapper;
    {
        // We start by setting up a global assembler that allows us to
        // obtain a dof mapper and the Dirichlet data
        gsOptionList assemblerOptions = gsGenericAssembler<>::defaultOptions();
        assemblerOptions.setInt("DirichletStrategy", dirichlet::elimination);
        assemblerOptions.setInt("InterfaceStrategy", iFace::dg);
        gsGenericAssembler<> assembler(
            mp,
            mb,
            assemblerOptions,
            &bc
        );
        assembler.computeDirichletDofs();
        //! [Define global mapper]

        //! [Define Ieti Mapper]
        ietiMapper.init(
            mb,
            assembler.system().rowMapper(0),
            assembler.fixedDofs(),
            gsIetidGMapper<>::allArtificialIfaces(mb)
        );
    }
    //! [Define Ieti Mapper]

    // Which primal dofs should we choose?
    bool cornersAsPrimals = false, edgesAsPrimals = false, facesAsPrimals = false;
    for (size_t i=0; i<primals.length(); ++i)
        switch (primals[i])
        {
            case 'c': cornersAsPrimals = true;   break;
            case 'e': edgesAsPrimals = true;     break;
            case 'f': facesAsPrimals = true;     break;
            default:
                gsInfo << "\nUnkown type of primal constraint: \"" << primals[i] << "\"\n";
                return EXIT_FAILURE;
        }

    // Compute the jump matrices
    bool fullyRedundant = true,
         noLagrangeMultipliersForCorners = cornersAsPrimals;
    //! [Define jumps]
    ietiMapper.computeJumpMatrices(fullyRedundant, noLagrangeMultipliersForCorners);
    //! [Define jumps]

    // We tell the ieti mapper which primal constraints we want; calling
    // more than one such function is possible.
    //! [Define primals]
    if (cornersAsPrimals)
        ietiMapper.cornersAsPrimals();

    if (edgesAsPrimals)
        ietiMapper.interfaceAveragesAsPrimals(mp,1);

    if (facesAsPrimals)
        ietiMapper.interfaceAveragesAsPrimals(mp,2);
    //! [Define primals]

    //! [Setup]
    // The ieti system does not have a special treatment for the
    // primal dofs. They are just one more subdomain
    gsIetiSystem<> ieti;
    ieti.reserve(nPatches+1);

    // The scaled Dirichlet preconditioner is independent of the
    // primal dofs.
    gsScaledDirichletPrec<> prec;
    prec.reserve(nPatches);

    // Setup the primal system, which needs to know the number of primal dofs.
    gsPrimalSystem<> primal(ietiMapper.nPrimalDofs());

    // Setup of the block-diagonal preconditioner for the saddle point problem
    // First, we need to know its size
    const index_t bdPrecSz = nPatches + 1 + (ietiMapper.nPrimalDofs()>0?1:0);
    gsBlockOp<>::Ptr bdPrec = gsBlockOp<>::make(bdPrecSz,bdPrecSz);


    if (eliminatePointwiseDofs)
        primal.setEliminatePointwiseConstraints(true);
    //! [Setup]


    // We set up the assembler
    gsOptionList assemblerOptions = gsGenericAssembler<>::defaultOptions();
    assemblerOptions.setInt("DirichletStrategy", dirichlet::elimination);
    assemblerOptions.setInt("InterfaceStrategy", iFace::dg);
    assemblerOptions.setSwitch("DG.OneSided", true);
    assemblerOptions.setReal("DG.Alpha", alpha);
    assemblerOptions.setReal("DG.Beta", beta);
    assemblerOptions.setReal("DG.Penalty", penalty);

    gsStopwatch timer;
    gsStopwatch assemblingTimer;
    assemblingTimer.restart();

    //! [Assemble]
    for (index_t k=0; k<nPatches; ++k)
    {
        gsInfo << "k=" << k << " of " << nPatches << "\n";
        gsInfo << "mb[k].numElements()" << mb[k].numElements() << "\n";
        gsInfo << "mb[k].size()" << mb[k].size() << "\n";

        // We use the local variants of everything
        gsBoundaryConditions<> bc_local;
        bc.getConditionsForPatch(k,bc_local);
        gsMultiPatch<> mp_local;
        gsMultiBasis<> mb_local;
        ietiMapper.localSpaces(mp,k,mp_local,mb_local);
/*
        // only required for the yeti footprint
        if(k == 4 || k == 12 || k == 36 || k == 60)
        {
            bc_local.addCornerValue(1,0,0);
            bc_local.addCornerValue(1,0,0);
        } else if(k == 5 || k == 13)
        {
            bc_local.addCornerValue(3, 0, 0);
        } else if(k == 38 || k == 46 || k == 54 || k == 62)
        {
            bc_local.addCornerValue(2, 0, 0);
        } else if(k == 47 || k == 55)
        {
            bc_local.addCornerValue(4, 0, 0);
        }
*/

        gsGenericAssembler<> assembler(
            mp_local,
            mb_local,
            assemblerOptions,
            &bc_local
        );

        // This function provides a new dof mapper and the Dirichlet data
        // This is necessary since it might happen that a 2d-patch touches the
        // Dirichlet boundary just with a corner or that a 3d-patch touches the
        // Dirichlet boundary with a corner or an edge. These cases are not
        // covered by bc.getConditionsForPatch
        assembler.refresh(ietiMapper.augmentedDofMapperLocal(k));
        assembler.setFixedDofVector(ietiMapper.augmentedFixedPart(k));
        assembler.system().reserve(mb_local, assemblerOptions, 1);

        // Assemble
        assembler.assembleStiffness(0,false);

        assembler.assembleMoments(f,0,false);
        gsBoundaryConditions<>::bcContainer neumannSides = bc_local.neumannSides();
        for (gsBoundaryConditions<>::bcContainer::const_iterator it = neumannSides.begin();
                it!= neumannSides.end(); ++it)
            assembler.assembleNeumann(*it,false);

        for (size_t i=0; i<ietiMapper.artificialIfaces(k).size(); ++i)
        {
            patchSide side1(0,ietiMapper.artificialIfaces(k)[i].assignedTo.side());
            patchSide side2(i+1,ietiMapper.artificialIfaces(k)[i].takenFrom.side());
            boundaryInterface bi(side1, side2, mp.geoDim());
            assembler.assembleDG(bi,false);
        }

        // Fetch data
        gsSparseMatrix<real_t, RowMajor> jumpMatrix  = ietiMapper.jumpMatrix(k);
        gsMatrix<>                       localRhs    = assembler.rhs();
        gsSparseMatrix<>                 localMatrix = assembler.matrix();
        //! [Assemble]
        DEBUGMAT(jumpMatrix);
        DEBUGMAT(localRhs);
        DEBUGMAT(localMatrix);

        //! [Patch to preconditioner]
        // Add the patch to the scaled Dirichlet preconditioner
        std::vector<index_t> skeletonDofs = ietiMapper.skeletonDofs(k);

//! inexact stuff

        // TODO: Use parameter domain assembler, which is also much faster!
        gsGenericAssembler<> hatassembler(
                mp.dim()==2 ? gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare()) : gsMultiPatch<>(*gsNurbsCreator<>::BSplineCube()),
                mb_local,
                assemblerOptions,
                &bc_local
        );

        hatassembler.refresh(ietiMapper.augmentedDofMapperLocal(k));
        hatassembler.setFixedDofVector(ietiMapper.augmentedFixedPart(k));
        hatassembler.system().reserve(mb_local, assemblerOptions, 1);

        hatassembler.assembleStiffness(0, false);
        for (size_t i=0; i<ietiMapper.artificialIfaces(k).size(); ++i)
        {
            patchSide side1(0,ietiMapper.artificialIfaces(k)[i].assignedTo.side());
            patchSide side2(i+1,ietiMapper.artificialIfaces(k)[i].takenFrom.side());
            boundaryInterface bi(side1, side2, mp.geoDim());
            ParameterdomainAssembler<real_t>::apply(hatassembler, mp_local, bi);
        }

        gsBoundaryConditions<> bc_A11 = bc_local;
        for (index_t i = 1; i <= 2*mb.dim(); ++i)
        {
            if (!bc_A11.getConditionFromSide(patchSide(0,i)))
                bc_A11.addCondition(0, i, condition_type::dirichlet, &gD);
        }

        gsScaledDirichletPrec<>::Blocks blocks
                = gsScaledDirichletPrec<>::matrixBlocks(hatassembler.matrix(), skeletonDofs);

        DEBUGMAT(hatassembler.matrix());
        DEBUGMAT(blocks.A00);
        DEBUGMAT(blocks.A11);

        timer.restart();
        gsLinearOperator<>::Ptr A11_fd = gsPatchPreconditionersCreator<>::fastDiagonalizationOp(mb_local.basis(0), bc_A11, assemblerOptions);
        times.setupInverseInMsD(k) += timer.stop();
        DEBUGMAT(*A11_fd);

        prec.addSubdomain(
                prec.restrictJumpMatrix(jumpMatrix, skeletonDofs).moveToPtr(),
                gsScaledDirichletPrec<>::schurComplement(blocks, A11_fd)
        );
        //! [Patch to preconditioner]

        // This function writes back to jumpMatrix, localMatrix, and localRhs,
        // so it must be called after prec.addSubdomain().
        //! [Patch to primals]
        gsSparseMatrix<> modifiedLocalMatrix, localEmbedding, embeddingForBasis;
        gsMatrix<> rhsForBasis;
        gsPrimalSystem<>::incorporateConstraints(
                ietiMapper.primalConstraints(k),
                eliminatePointwiseDofs,
                localMatrix,
                modifiedLocalMatrix,
                localEmbedding,
                embeddingForBasis,
                rhsForBasis
        );
        DEBUGMAT(localMatrix);
        DEBUGMAT(modifiedLocalMatrix);
        DEBUGMAT(localEmbedding);
        DEBUGMAT(embeddingForBasis);
        DEBUGMAT(rhsForBasis);

        gsMatrix<>                       modifiedLocalRhs     = localEmbedding * localRhs;
        gsSparseMatrix<real_t, RowMajor> modifiedJumpMatrix   = jumpMatrix * localEmbedding.transpose();

        DEBUGMAT(modifiedLocalRhs);
        DEBUGMAT(modifiedJumpMatrix);

        //! [Patch to primals]
        real_t reg;
        if(bc_local.dirichletSides().size() == 0)
            reg = (real_t)1;
        else
            reg = (real_t)0;

        // Set up the local preconditioners
        if (eliminatePointwiseDofs && cornersAsPrimals)
        {
            for(index_t i = 1; i <= 1<<mb.dim(); i++) // Since corners are always primals
                bc_local.addCornerValue(i, 0, 0);
        }

        timer.restart();
        gsLinearOperator<>::Ptr fastdiagOp = gsPatchPreconditionersCreator<>::fastDiagonalizationOp(
            mb_local.basis(0), bc_local, assemblerOptions, reg, (real_t)1, (real_t)0, true);
        DEBUGMAT(*fastdiagOp);
        {
            // TODO: Scaling of the remainder term
            gsMatrix<> dis;
            fastdiagOp->toMatrix(dis);
            gsInfo << dis << "\n";
        }

        //localPrec does not work propperly now, since it needs to consider that the
        //corners are "there but eliminated" (diag = 1)
        gsLinearOperator<>::Ptr localPrec = gsInexactIETIPrec<real_t>::make(
                                                                       localMatrix,
                                                                       fastdiagOp,
                                                                       mb_local,
                                                                       ietiMapper.artificialIfaces(k),
                                                                       ietiMapper.augmentedDofMapperLocal(k),
                                                                       penalty,
                                                                       assemblerOptions
                                                                       );
        DEBUGMAT(*localPrec);
        {
            gsMatrix<> dis;
            localPrec->toMatrix(dis);
            gsInfo << dis << "\n";
        }

        DEBUGVAR(localEmbedding.rows()-localEmbedding.cols());
        gsLinearOperator<>::Ptr modifiedLocalPrec;
        if (localEmbedding.rows()-localEmbedding.cols())
            modifiedLocalPrec = gsBlockSolverOp<>::make( localPrec,
                modifiedLocalMatrix.bottomRows(localEmbedding.rows()-localEmbedding.cols()).leftCols(localPrec->rows()) );
        else
            modifiedLocalPrec = localPrec;

        DEBUGMAT(*modifiedLocalPrec);

        times.setupPk(k) += timer.stop();

        gsLinearOperator<>::Ptr localSolver = gsIterativeSolverOp<gsConjugateGradient<> >::make(modifiedLocalMatrix, modifiedLocalPrec);
        DEBUGMAT(*localSolver);


        /// "real" elimination
        timer.restart();
        gsSparseMatrix<> basisDel = gsPrimalSystem<>::primalBasis(
                localSolver,
                embeddingForBasis,
                rhsForBasis,
                ietiMapper.primalDofIndices(k),
                primal.nPrimalDofs()
        );
        DEBUGMAT(basisDel);
        DEBUGVAR(basisDel.toDense());

        gsInfo << "ietiMapper.primalDofIndices(k): "; for (index_t idx : ietiMapper.primalDofIndices(k)) { gsInfo << idx << " "; } gsInfo << "\n";

        times.assemblePrimalBasis(k) += timer.stop();

        // TODO: this has to go into gsPrimalSystem<>::primalBasis ?

        /*if(eliminatePointwiseDofs)
        {
            for (size_t i = 0; i < ietiMapper.primalConstraints(k).size(); ++i) {
                gsSparseVector<> constraint = ietiMapper.primalConstraints(k)[i];
                for (int j = 0; j < constraint.size(); ++j) {
                    if(constraint[j])
                    {
                        //TODO: This seems to be wrong
                        GISMO_ASSERT(math::abs(basisDel(j, ietiMapper.primalDofIndices(k)[i]) - 1) < 1e-4,
                              "("<<j<<", ietiMapper.primalDofIndices("<<k<<")["<<i<<"] [aka."<<ietiMapper.primalDofIndices(k)[i]<<"]) = "
                              <<basisDel(j, ietiMapper.primalDofIndices(k)[i])<<"!=1");
                        basisDel(j, ietiMapper.primalDofIndices(k)[i]) = (real_t)1;
                    }
                }
            }
        }*/

        primal.addContribution(
                jumpMatrix,
                localMatrix,
                localRhs,
                basisDel
                //,makeMatrixOp(localEmbedding.moveToPtr())
        );
/*
        timer.restart();
        gsLinearOperator<>::Ptr localSolver = makeSparseLUSolver(localMatrix);
        times.setupPk(k) += timer.stop();
*/
        // Register the local solver to the block preconditioner. We use
        // a sparse LU solver since the local saddle point problem is not
        // positive definite.
        bdPrec->addOperator(k,k, gsTimedOp<>::make("P "+std::to_string(k), modifiedLocalPrec));


        // Add the patch to the Ieti system
        //! [Patch to system]
        ieti.addSubdomain(
            modifiedJumpMatrix.moveToPtr(),
            makeMatrixOp(modifiedLocalMatrix.moveToPtr()),
            give(modifiedLocalRhs)
        );
        gsInfo << "Assembling fin.\n";
        //! [Patch to system]
    //! [End of assembling loop]
    } // end for
    //! [End of assembling loop]
    times.totalAssemble(0) += assemblingTimer.stop();

    // Add the primal problem if there are primal constraints
    //! [Primal to system]
    if (ietiMapper.nPrimalDofs()>0)//TODO
    {
        // It is not required to provide a local solver to .addSubdomain,
        // since a sparse LU solver would be set up on the fly if required.
        // Here, we make use of the fact that we can use a Cholesky solver
        // because the primal problem is symmetric and positive definite:
        bdPrec->addOperator(nPatches, nPatches, makeSparseLUSolver(primal.localMatrix())); //TODO: Cholesky?

        // Add to IETI system
        ieti.addSubdomain(
            primal.jumpMatrix().moveToPtr(),
            makeMatrixOp(primal.localMatrix().moveToPtr()),
            give(primal.localRhs())
        );
    }
    //! [Primal to system]

    gsInfo << "done. " << ietiMapper.nPrimalDofs() << " primal dofs.\n";

    /**************** Setup solver and solve ****************/

    gsInfo << "Setup solver and solve... \n"
        "    Setup multiplicity scaling... " << std::flush;

    // Tell the preconditioner to set up the scaling
    //! [Setup scaling]
    prec.setupMultiplicityScaling();

    // The scaled Dirichlet preconditioner is in the last block
    gsLinearOperator<>::Ptr sdPrec = prec.preconditioner();
    bdPrec->addOperator(bdPrecSz-1,bdPrecSz-1, sdPrec);
    //! [Setup scaling]

    gsInfo << "done.\n    Setup minres solver and solve... \n" << std::flush;
    // Initial guess
    //! [Define initial guess]
    gsMatrix<> x;
    //TODO: Not possible: x.setRandom( bdPrec->rows(), 1 );
    x.setZero( bdPrec->rows(), 1 );

    //! [Define initial guess]

    //gsInfo << ieti.rhsForSaddlePoint().transpose() << "\n";
    //gsInfo << x.transpose() << "\n";

    gsMatrix<> errorHistory;

    // This is the main cg iteration
    //! [Solve]
    timer.restart();
    //gsMinimalResidual
    gsGMRes<>( ieti.saddlePointProblem(), gsTimedOp<>::make("Preconditioner", bdPrec) )
            .setOptions( cmd.getGroup("Solver") )
            .solveDetailed( ieti.rhsForSaddlePoint(), x, errorHistory );
    times.totalSolving(0) = timer.stop();
    //! [Solve]

    gsInfo << "done.\n    Reconstruct solution from Lagrange multipliers... " << std::flush;
    // Now, we want to have the global solution for u
    //! [Recover]
    gsMatrix<> uVec = ietiMapper.constructGlobalSolutionFromLocalSolutions(
        primal.distributePrimalSolution(
            ieti.constructSolutionFromSaddlePoint(x)
        )
    );
    //! [Recover]
    gsInfo << "done.\n\n";

    /******************** Print end Exit ********************/

    gsInfo << "Print timings ... \n";
    times.print(gsInfo);

    const index_t iter = errorHistory.rows()-1;
    const bool success = errorHistory(iter,0) < tolerance;
    if (success)
        gsInfo << "Reached desired tolerance after " << iter << " iterations:\n";
    else
        gsInfo << "Did not reach desired tolerance after " << iter << " iterations:\n";

    if (errorHistory.rows() < 20)
        gsInfo << errorHistory.transpose() << "\n\n";
    else
        gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

    if (!out.empty())
    {
        gsFileData<> fd;
        std::time_t time = std::time(NULL);
        //fd.add(cmd);
        //fd.add(uVec);
        //gsMatrix<> mat; ieti.saddlePointProblem()->toMatrix(mat); fd.add(mat);
        fd.addComment(std::string("ietidG3_example   Timestamp:")+std::ctime(&time));
        fd.save(out);
        gsInfo << "Write solution to file " << out << "\n";
    }

    if (plot)
    {
        gsInfo << "Write Paraview data to file ieti_result.pvd\n";
        gsOptionList assemblerOptions = gsGenericAssembler<>::defaultOptions();
        assemblerOptions.setInt("DirichletStrategy", dirichlet::elimination);
        assemblerOptions.setInt("InterfaceStrategy", iFace::dg);
        gsGenericAssembler<> assembler(
            mp,
            mb,
            assemblerOptions,
            &bc
        );
        assembler.computeDirichletDofs();
        gsMultiPatch<> mpsol;
        assembler.constructSolution(uVec, mpsol);
        //gsField<> sol( assembler.patches(), piece );
        //gsWriteParaview<>(sol, "IETI", 1000);
        system("paraview IETI.pvd  &");

        //gsFileManager::open("ieti_result.pvd");
    }

    if (!plot&&out.empty())
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution or --out to write solution to xml file.\n";
    }
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
