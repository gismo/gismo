/** @file ExprAssembler.h

    @brief Assembler for the new expression module.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsNewExpressions/FeSpaceData.h>
#include <gsNewExpressions/ExpressionHelper.h>
#include <gsMatrix/gsFiberMatrix.h>
#include <gsMatrix/gsSparseMatrix.h>
#include <gsDomain/gsDomain.h>
#include <gsNewExpressions/SpaceObject.h> // Required for SpaceObject definition
#include <gsAssembler/gsQuadrature.h> // For gsQuadrature
#include <gsAssembler/gsGaussRule.h>  // For gsGaussRule (sum factorization)
#include <gsTensor/gsTensorBasis.h>   // For gsTensorBasis (sum factorization)
#include <gsIO/gsOptionList.h> // For gsOptionList

#include <list>
#include <vector>

namespace gismo
{

// Forward declaration
template<class T>
class ExprEvaluator;

template <class T>
class ExprAssembler
{
public:
    // Constructors
    ExprAssembler(const gsDomain<T>& domain, index_t num_test_spaces = 1, index_t num_trial_spaces = 1)
        : m_domain(domain),
          m_test_spaces(num_test_spaces, nullptr),
          m_trial_spaces(num_trial_spaces, nullptr),
          m_options(defaultOptions())
    {
        // Initialize mappers
        resetDimensions();
    }

    // Methods to create expression objects, making the API transparent
    
    // Forwarded methods for Constant and Variable objects
    Expr::ConstantObject<T,0> getConstant(const T s, std::string label="c") { return m_helper.getConstant(s, label); }
    Expr::ConstantObject<T,1> getConstant(const gsVector<T> & v, std::string label="C") { return m_helper.getConstant(v, label); }
    Expr::ConstantObject<T,2> getConstant(const gsMatrix<T> & m, std::string label="C") { return m_helper.getConstant(m, label); }

    Expr::VariableObject<T,0,true> getScalarFunction(const gsConstantFunction<T> & cfunc, std::string label="f") { return m_helper.getScalarFunction(cfunc, label); }
    Expr::VariableObject<T,0,false> getScalarFunction(const gsFunction<T> & func, std::string label="f") { return m_helper.getScalarFunction(func, label); }
    Expr::VariableObject<T,1,true> getVectorFunction(const gsConstantFunction<T> & cfunc, std::string label="F") { return m_helper.getVectorFunction(cfunc, label); }
    Expr::VariableObject<T,1,false> getVectorFunction(const gsFunction<T> & func, std::string label="F") { return m_helper.getVectorFunction(func, label); }
    
    /**
     * @brief Create a geometry map expression
     * @param geom The geometry (typically gsMultiPatch)
     * @return GeometryMap expression
     */
    Expr::GeometryMap<T> getMap(const gsFunctionSet<T> & geom)
    {
        return m_helper.getMap(geom);
    }

    // Space creation methods
    Expr::SpaceObject<T, Expr::SpaceType::Test, 0> getScalarTestSpace(const gsFunctionSet<T>& basis, index_t id = 0, std::string label="φ")
    {
        GISMO_ASSERT(basis.targetDim()==1,"Function is not scalar");
        GISMO_ASSERT(id < m_test_spaces.size(), "Test space ID exceeds pre-allocated size.");

        m_space_data_list.emplace_back(basis, 1, id); // Assuming scalar space, dim=1
        Expr::FeSpaceData<T>* space_data = &m_space_data_list.back();
        space_data->init();
        
        Expr::SpaceObject<T, Expr::SpaceType::Test, 0> expr(basis.domainDim(),std::array<size_t,0>{}, id, label);
        expr.setSource(basis);
        m_helper.add(expr); // Register with helper for eval()
        expr.setSpaceData(*space_data);
        
        m_test_spaces[id] = space_data;
        return expr;
    }

    Expr::SpaceObject<T, Expr::SpaceType::Trial, 0> getScalarTrialSpace(const gsFunctionSet<T>& basis, index_t id = 0, std::string label="ψ")
    {
        GISMO_ASSERT(basis.targetDim()==1,"Function is not scalar");
        GISMO_ASSERT(id < m_trial_spaces.size(), "Trial space ID exceeds pre-allocated size.");
        
        m_space_data_list.emplace_back(basis, 1, id); // Assuming scalar space, dim=1
        Expr::FeSpaceData<T>* space_data = &m_space_data_list.back();
        space_data->init();

        Expr::SpaceObject<T, Expr::SpaceType::Trial, 0> expr(basis.domainDim(),std::array<size_t,0>{}, id, label);
        expr.setSource(basis);
        m_helper.add(expr); // Register with helper for eval()
        expr.setSpaceData(*space_data);
        
        m_trial_spaces[id] = space_data;
        return expr;
    }

    Expr::SpaceObject<T, Expr::SpaceType::Test, 1> getVectorTestSpace(const gsFunctionSet<T>& basis, index_t dim, index_t id = 0, std::string label="v")
    {
        GISMO_ASSERT(id < m_test_spaces.size(), "Test space ID exceeds pre-allocated size.");

        m_space_data_list.emplace_back(basis, dim, id); // Vector space with explicit dimension
        Expr::FeSpaceData<T>* space_data = &m_space_data_list.back();
        space_data->init();
        
        std::array<size_t,1> sizes = {(size_t)dim};
        Expr::SpaceObject<T, Expr::SpaceType::Test, 1> expr(basis.domainDim(), sizes, id, label);
        expr.setSource(basis);
        m_helper.add(expr); // Register with helper for eval()
        expr.setSpaceData(*space_data);
        
        m_test_spaces[id] = space_data;
        return expr;
    }

    Expr::SpaceObject<T, Expr::SpaceType::Trial, 1> getVectorTrialSpace(const gsFunctionSet<T>& basis, index_t dim, index_t id = 0, std::string label="u")
    {
        GISMO_ASSERT(id < m_trial_spaces.size(), "Trial space ID exceeds pre-allocated size.");
        
        m_space_data_list.emplace_back(basis, dim, id); // Vector space with explicit dimension
        Expr::FeSpaceData<T>* space_data = &m_space_data_list.back();
        space_data->init();

        std::array<size_t,1> sizes = {(size_t)dim};
        Expr::SpaceObject<T, Expr::SpaceType::Trial, 1> expr(basis.domainDim(), sizes, id, label);
        expr.setSource(basis);
        m_helper.add(expr); // Register with helper for eval()
        expr.setSpaceData(*space_data);
        
        m_trial_spaces[id] = space_data;
        return expr;
    }

    // Solution object creation - for nonlinear problems
    template<enum Expr::SpaceType _Space, size_t _Order>
    Expr::SolutionObject<T, _Space, _Order> getSolution(
        const Expr::SpaceObject<T, _Space, _Order>& space,
        gsMatrix<T>& solVector,
        std::string label="s")
    {
        return Expr::SolutionObject<T, _Space, _Order>(space, solVector, label);
    }

    // Initialize system
    void initSystem(const index_t numRhs = 1)
    {
        // Check all spaces are initialized
        checkSpacesInitialized();
        
        initMatrix();
        m_rhs.setZero(numTestDofs(), numRhs);
    }
    
    void initMatrix()
    {
        resetDimensions();
        clearMatrix(false);
    }
    
    /// @brief Check that all registered spaces have been set up
    void checkSpacesInitialized() const
    {
        for (size_t i = 0; i < m_test_spaces.size(); ++i)
        {
            if (m_test_spaces[i] != nullptr)
            {
                GISMO_ASSERT(m_test_spaces[i]->isInitialized(),
                    "Test space " << i << " not initialized. Call space.setup() first.");
            }
        }
        for (size_t i = 0; i < m_trial_spaces.size(); ++i)
        {
            if (m_trial_spaces[i] != nullptr)
            {
                GISMO_ASSERT(m_trial_spaces[i]->isInitialized(),
                    "Trial space " << i << " not initialized. Call space.setup() first.");
            }
        }
    }

    // Main assembly method
    template <class E>
    void assemble(const Expr::BaseExpression<E>& expr)
    {
        // Determine if matrix or vector assembly
        bool is_matrix_assembly = false;
        if (expr.Space == Expr::SpaceType::Both) { // Assuming SpaceType::Both for matrix
            is_matrix_assembly = true;
        } else if (expr.Space == Expr::SpaceType::Test) { // Assuming SpaceType::Test for vector
            is_matrix_assembly = false;
        } else {
            GISMO_ERROR("Unsupported expression for assembly. Must involve space objects.");
        }

        expr.parse(m_helper); // Parse expression once

        // Prepare local accumulation
        #pragma omp parallel
        {
            typename gsQuadRule<T>::uPtr quad_rule;
            index_t current_patch = -1;
            
            _assembler_op ee(m_helper, m_fmatrix, m_rhs, is_matrix_assembly); // Pass is_matrix_assembly
            
            #pragma omp for
            for (auto & elem : m_domain.allElements())
            {
                if (current_patch != elem.patch())
                {
                    current_patch = elem.patch();
                    // Using trial space 0 for quadrature basis
                    quad_rule = gsQuadrature::getPtr(m_trial_spaces[0]->fs->basis(current_patch), m_options);
                }

                quad_rule->mapTo(elem.lowerCorner(), elem.upperCorner(), m_helper.points(), m_helper.weights());
                
                m_helper.precompute(current_patch); // Precompute values at quadrature points

                Expr::ExpressionResult<T> eval_val = expr.eval(0);
                
                // Get the space data to compute proper dimensions
                // Cast to derived type to access test() and trial() methods
                const E& derived_expr = static_cast<const E&>(expr);
                auto test_space_data = derived_expr.test().spaceData();
                auto trial_space_data = derived_expr.trial().spaceData();
                
                // Get number of active basis functions
                const gsMatrix<index_t>& test_actives = m_helper.funcData(test_space_data->fs).actives;
                const gsMatrix<index_t>& trial_actives = m_helper.funcData(trial_space_data->fs).actives;
                
                // Initialize localMat with proper size based on number of basis functions
                index_t n_test = test_actives.rows();
                index_t n_trial = trial_actives.rows();
                ee.localMat.setZero(n_test, n_trial);
                
                // Accumulate over all quadrature points
                for (index_t k = 0; k < m_helper.points().cols(); ++k)
                {
                    eval_val = expr.eval(k);
                    T weight = m_helper.weights()[k];
                    
                    // Accumulate weighted evaluation for each basis function pair
                    // Handle both scalar and tensor-valued expressions
                    for (index_t i = 0; i < n_test; ++i)
                    {
                        for (index_t j = 0; j < n_trial; ++j)
                        {
                            // eval_val(i,j) is a gsMatrix<T>
                            // For scalar expressions: (1x1) matrix
                            // For vector expressions: (n x 1) vector
                            // For tensor expressions: (n x m) matrix
                            // We accumulate all entries, weighted by quadrature weight
                            const auto& mat = eval_val(i, j);
                            ee.localMat(i, j) += mat.sum() * weight;
                        }
                    }
                }
                
                // Push local contribution to global system
                // Call test()/trial() on derived type E, not BaseExpression
                ee.push(test_space_data, trial_space_data, elem.patch()); // Pass relevant FeSpaceData and patch ID
            }
        }
        m_fmatrix.toSparseMatrix_into(m_matrix);
    }


    /**
     * @brief Sum factorization assembly for tensor-product bases
     * 
     * This method exploits the tensor-product structure of the basis functions
     * to reduce the computational complexity of matrix assembly.
     * 
     * For a mass matrix in d dimensions with n basis functions per direction
     * and q quadrature points per direction:
     * - Standard assembly: O(n^2d * q^d)
     * - Sum factorization: O(n^(d+1) * q^d)
     * 
     * The key insight is that for tensor-product bases:
     *   B_{ijk}(x,y,z) = B_i(x) * B_j(y) * B_k(z)
     * 
     * So the mass matrix entry M_{ijk,lmn} can be written as:
     *   M_{ijk,lmn} = (∫B_i B_l dx) * (∫B_j B_m dy) * (∫B_k B_n dz)
     *              = M^x_{il} * M^y_{jm} * M^z_{kn}
     * 
     * This is a Kronecker product: M = M^z ⊗ M^y ⊗ M^x
     * 
     * @tparam E Expression type
     * @param expr The bilinear form expression to assemble
     * 
     * @note This method requires tensor-product bases (gsTensorBasis).
     *       For non-tensor bases, use the standard assemble() method.
     */
    template <class E>
    void assembleSF(const Expr::BaseExpression<E>& expr)
    {
        // Check if this is matrix assembly
        GISMO_ASSERT(expr.Space == Expr::SpaceType::Both,
            "Sum factorization currently only supports bilinear forms (matrix assembly)");
        
        expr.parse(m_helper);
        
        // Get test and trial space data
        const E& derived_expr = static_cast<const E&>(expr);
        auto test_space_data = derived_expr.test().spaceData();
        auto trial_space_data = derived_expr.trial().spaceData();
        
        GISMO_ASSERT(test_space_data != nullptr && trial_space_data != nullptr,
            "Space data must be available for sum factorization");
        
        // Get the basis - must be tensor product
        const gsFunctionSet<T>* test_fs = test_space_data->fs;
        
        // Iterate over patches
        for (index_t patch = 0; patch < test_fs->nPieces(); ++patch)
        {
            const gsBasis<T>& test_basis = test_fs->basis(patch);
            
            // Get domain dimension
            const short_t d = test_basis.domainDim();
            
            // Try to cast to tensor basis
            // For now, we handle the common case of gsTensorBSplineBasis
            assembleSF_patch(expr, test_space_data, trial_space_data, patch, d);
        }
        
        m_fmatrix.toSparseMatrix_into(m_matrix);
    }

private:
    /**
     * @brief Sum factorization assembly for a single patch
     */
    template <class E>
    void assembleSF_patch(const Expr::BaseExpression<E>& expr,
                          Expr::FeSpaceData<T>* test_space_data,
                          Expr::FeSpaceData<T>* trial_space_data,
                          index_t patch,
                          short_t d)
    {
        const gsFunctionSet<T>* test_fs = test_space_data->fs;
        const gsFunctionSet<T>* trial_fs = trial_space_data->fs;
        const gsBasis<T>& test_basis = test_fs->basis(patch);
        const gsBasis<T>& trial_basis = trial_fs->basis(patch);
        
        // Get 1D quadrature rule
        gsGaussRule<T> quad1D;
        
        // Store 1D mass matrices for each direction
        std::vector<gsMatrix<T>> mass1D(d);
        
        // For each direction, compute the 1D "mass" contributions
        for (short_t dir = 0; dir < d; ++dir)
        {
            // Get the 1D component basis (assuming tensor structure)
            // This requires the basis to support component access
            const gsBasis<T>* test_comp = nullptr;
            const gsBasis<T>* trial_comp = nullptr;
            
            // Try to get component basis via dynamic cast to tensor basis
            if (auto* tb = dynamic_cast<const gsTensorBasis<2,T>*>(&test_basis))
            {
                test_comp = &tb->component(dir);
                trial_comp = &(dynamic_cast<const gsTensorBasis<2,T>*>(&trial_basis))->component(dir);
            }
            else if (auto* tb3 = dynamic_cast<const gsTensorBasis<3,T>*>(&test_basis))
            {
                test_comp = &tb3->component(dir);
                trial_comp = &(dynamic_cast<const gsTensorBasis<3,T>*>(&trial_basis))->component(dir);
            }
            else
            {
                GISMO_ERROR("Sum factorization requires tensor-product bases. "
                           "Use assemble() for non-tensor bases.");
            }
            
            // Setup 1D quadrature
            const index_t deg = math::max(test_comp->maxDegree(), trial_comp->maxDegree());
            gsVector<index_t> numNodes(1);
            numNodes[0] = static_cast<index_t>(m_options.getReal("quA") * deg + m_options.getInt("quB"));
            quad1D.setNodes(numNodes);
            
            // Get number of 1D basis functions
            index_t n_test = test_comp->size();
            index_t n_trial = trial_comp->size();
            
            // Initialize 1D mass matrix
            mass1D[dir].setZero(n_test, n_trial);
            
            // Get domain for this direction and iterate over elements
            auto domain1D = test_comp->domain();
            
            for (auto elemIt = domain1D->beginAll(); elemIt != domain1D->endAll(); ++elemIt)
            {
                T a = elemIt.lowerCorner().value();
                T b = elemIt.upperCorner().value();
                
                // Map quadrature to element [a,b]
                gsMatrix<T> nodes1D;
                gsVector<T> weights1D;
                quad1D.mapTo(a, b, nodes1D, weights1D);
                
                // Evaluate 1D basis functions on this element
                gsMatrix<T> test_vals, trial_vals;
                gsMatrix<index_t> test_actives, trial_actives;
                test_comp->active_into(nodes1D.col(0), test_actives);
                trial_comp->active_into(nodes1D.col(0), trial_actives);
                test_comp->eval_into(nodes1D, test_vals);
                trial_comp->eval_into(nodes1D, trial_vals);
                
                // Accumulate 1D mass matrix: M^dir_{ij} = sum_k w_k * B_i(x_k) * B_j(x_k)
                for (index_t k = 0; k < nodes1D.cols(); ++k)
                {
                    T w = weights1D[k];
                    for (index_t i = 0; i < test_actives.rows(); ++i)
                    {
                        for (index_t j = 0; j < trial_actives.rows(); ++j)
                        {
                            mass1D[dir](test_actives(i,0), trial_actives(j,0)) += w * test_vals(i, k) * trial_vals(j, k);
                        }
                    }
                }
            }
        }
        
        // Now compute the full matrix via Kronecker product
        // M = M^{d-1} ⊗ M^{d-2} ⊗ ... ⊗ M^0
        // We use the property that (A ⊗ B) vec(X) = vec(B X A^T)
        
        // For assembly, we need to map local to global indices
        const gsDofMapper& test_mapper = test_space_data->mapper;
        const gsDofMapper& trial_mapper = trial_space_data->mapper;
        
        // Get total sizes
        gsVector<index_t, -1> test_sizes(d), trial_sizes(d);
        index_t total_test = 1, total_trial = 1;
        for (short_t dir = 0; dir < d; ++dir)
        {
            test_sizes[dir] = mass1D[dir].rows();
            trial_sizes[dir] = mass1D[dir].cols();
            total_test *= test_sizes[dir];
            total_trial *= trial_sizes[dir];
        }
        
        // Compute full local matrix using Kronecker product structure
        // For efficiency, we iterate over tensor indices
        gsVector<index_t, -1> test_idx(d), trial_idx(d);
        test_idx.setZero();
        trial_idx.setZero();
        
        // Iterate over all (test, trial) pairs
        for (index_t i = 0; i < total_test; ++i)
        {
            // Convert flat index to tensor index for test
            index_t tmp = i;
            for (short_t dir = 0; dir < d; ++dir)
            {
                test_idx[dir] = tmp % test_sizes[dir];
                tmp /= test_sizes[dir];
            }
            
            // Get global test index
            index_t global_test = test_mapper.index(i, patch, 0);
            if (!test_mapper.is_free_index(global_test))
                continue;
            
            for (index_t j = 0; j < total_trial; ++j)
            {
                // Convert flat index to tensor index for trial
                tmp = j;
                for (short_t dir = 0; dir < d; ++dir)
                {
                    trial_idx[dir] = tmp % trial_sizes[dir];
                    tmp /= trial_sizes[dir];
                }
                
                // Get global trial index
                index_t global_trial = trial_mapper.index(j, patch, 0);
                if (!trial_mapper.is_free_index(global_trial))
                    continue;
                
                // Compute Kronecker product entry
                // M_{i,j} = prod_{dir} M^{dir}_{test_idx[dir], trial_idx[dir]}
                T val = mass1D[0](test_idx[0], trial_idx[0]);
                for (short_t dir = 1; dir < d; ++dir)
                {
                    val *= mass1D[dir](test_idx[dir], trial_idx[dir]);
                }
                
                if (math::abs(val) > 1e-15)
                {
                    m_fmatrix.coeffRef(global_test, global_trial) += val;
                }
            }
        }
    }

public:
    /// @brief Assemble multiple expressions at once (for block systems)
    template <class... Exprs>
    void assemble(const Exprs&... exprs)
    {
        // Use a dummy array to force expansion in order
        int dummy[] = {(assemble_single(exprs), 0)...};
        (void)dummy; // Suppress unused variable warning
    }

    // Getters for results
    const gsSparseMatrix<T>& matrix() const { return m_matrix; }
    const gsVector<T>& rhs() const { return m_rhs; }
    
    /// @brief Returns a block view of the system matrix
    /// Each block corresponds to a different test/trial space pair
    typename gsSparseMatrix<T>::BlockView matrixBlockView()
    {
        gsVector<index_t> rowSizes, colSizes;
        _blockDims(rowSizes, colSizes);
        return m_matrix.blockView(rowSizes, colSizes);
    }
    
    /// @brief Returns a const block view of the system matrix
    typename gsSparseMatrix<T>::constBlockView matrixBlockView() const
    {
        gsVector<index_t> rowSizes, colSizes;
        _blockDims(rowSizes, colSizes);
        return m_matrix.blockView(rowSizes, colSizes);
    }

private:
    /// @brief Helper to assemble a single expression (used by variadic assemble)
    template <class E>
    void assemble_single(const Expr::BaseExpression<E>& expr)
    {
        assemble(expr);
    }

public:

    // Allow evaluator to access the helper
    const ExpressionHelper<T>& helper() const { return m_helper; }
    
    // Options
    gsOptionList & options() {return m_options;}
    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
        opt.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
        return opt;
    }

    index_t numBlocks() const
    {
        index_t nb = 0;
        for (const auto* space_data : m_test_spaces)
            if (space_data)
                nb += space_data->dim;
        return nb;
    }

private:
    const gsDomain<T>& m_domain;
    ExpressionHelper<T> m_helper;

    std::list<Expr::FeSpaceData<T>> m_space_data_list;
    std::vector<Expr::FeSpaceData<T>*> m_test_spaces;
    std::vector<Expr::FeSpaceData<T>*> m_trial_spaces;

    gsFiberMatrix<T,ColMajor> m_fmatrix;
    gsSparseMatrix<T> m_matrix;
    gsVector<T> m_rhs;
    gsOptionList m_options;

    void clearMatrix(const bool& save_sparsety_pattern = true)
    {
        if (m_fmatrix.nonZeros() && save_sparsety_pattern)
        {
            m_fmatrix.assignZero();
        }
        else
        {
            m_fmatrix.resize(numTestDofs(), numDofs());

            if (0 == m_fmatrix.rows() || 0 == m_fmatrix.cols())
                gsWarn << " No internal DOFs, zero sized system.\n";
            // Not implementing reservation for now as domain().degree() is not straightforward
        }
    }

    // Helper method to reset dimensions of DoF mappers
    void resetDimensions()
    {
        // Set shifts for all spaces - shifts are used to compute global indices
        // for multi-unknown systems
        index_t current_shift = 0;
        
        // First, set shifts for trial spaces
        for (auto* space_data : m_trial_spaces)
        {
            if (space_data)
            {
                if (!space_data->mapper.isFinalized()) 
                    space_data->init();
                space_data->mapper.setShift(current_shift);
                current_shift += space_data->dim * space_data->mapper.freeSize();
            }
        }
        
        // Reset shift for test spaces if they don't share mappers with trial
        current_shift = 0;
        for (auto* space_data : m_test_spaces)
        {
            if (space_data)
            {
                if (!space_data->mapper.isFinalized()) 
                    space_data->init();
                space_data->mapper.setShift(current_shift);
                current_shift += space_data->dim * space_data->mapper.freeSize();
            }
        }
    }
    
    /// @brief Compute block dimensions for block matrix view
    void _blockDims(gsVector<index_t>& rowSizes, gsVector<index_t>& colSizes) const
    {
        // For multiple spaces, each space corresponds to a block
        rowSizes.resize(m_test_spaces.size());
        for (size_t r = 0; r < m_test_spaces.size(); ++r)
        {
            if (m_test_spaces[r])
                rowSizes[r] = m_test_spaces[r]->mapper.freeSize();
            else
                rowSizes[r] = 0;
        }
        
        colSizes.resize(m_trial_spaces.size());
        for (size_t c = 0; c < m_trial_spaces.size(); ++c)
        {
            if (m_trial_spaces[c])
                colSizes[c] = m_trial_spaces[c]->mapper.freeSize();
            else
                colSizes[c] = 0;
        }
    }
    
    // Helper to get total number of DoFs (sum of free DoFs from all trial spaces)
    index_t numDofs() const
    {
        index_t total_dofs = 0;
        for(const auto* space_data : m_trial_spaces) {
            if (space_data && space_data->mapper.isFinalized()) {
                total_dofs += space_data->mapper.freeSize();
            }
        }
        return total_dofs;
    }

    // Helper to get total number of Test DoFs (sum of free DoFs from all test spaces)
    index_t numTestDofs() const
    {
        index_t total_test_dofs = 0;
        for(const auto* space_data : m_test_spaces) {
            if (space_data && space_data->mapper.isFinalized()) {
                total_test_dofs += space_data->mapper.freeSize();
            }
        }
        return total_test_dofs;
    }

    // _assembler_op struct for evaluation and pushing contributions
    struct _assembler_op
    {
        ExpressionHelper<T>& m_helper;
        gsFiberMatrix<T,ColMajor> & m_fmatrix;
        gsVector<T>               & m_rhs;
        bool m_is_matrix_assembly;
        gsMatrix<T>                 localMat; // Local matrix/vector for an element

        _assembler_op(ExpressionHelper<T>& helper, gsFiberMatrix<T,ColMajor> & fmatrix, gsVector<T> & rhs, bool is_matrix)
            : m_helper(helper), m_fmatrix(fmatrix), m_rhs(rhs), m_is_matrix_assembly(is_matrix) {}

        // Push local contributions to global system
        void push(Expr::FeSpaceData<T>* row_space_data, Expr::FeSpaceData<T>* col_space_data, index_t patch_id)
        {
            // This is where the local contributions (localMat) get added to the global system (m_fmatrix/m_rhs)
            // using the dof mappers and actives from FeSpaceData.
            GISMO_ASSERT(row_space_data != nullptr, "Row space data is null!");
            
            const gsDofMapper& row_mapper = row_space_data->mapper;

            // Get active basis functions for the row space at this patch
            const gsMatrix<index_t>& row_actives = m_helper.funcData(row_space_data->fs).actives;

            if (m_is_matrix_assembly)
            {
                GISMO_ASSERT(col_space_data != nullptr, "Col space data is null for matrix assembly!");
                const gsDofMapper& col_mapper = col_space_data->mapper;
                
                // Get active basis functions for the col space at this patch
                const gsMatrix<index_t>& col_actives = m_helper.funcData(col_space_data->fs).actives;

                // localMat has dimensions (row_actives.rows(), col_actives.rows())
                GISMO_ASSERT(localMat.rows() == row_actives.rows(), "Local matrix row size mismatch");
                GISMO_ASSERT(localMat.cols() == col_actives.rows(), "Local matrix col size mismatch");

                for (index_t r = 0; r < row_actives.rows(); ++r)
                {
                    for (index_t c_comp = 0; c_comp < row_space_data->dim; ++c_comp)
                    {
                        index_t global_row_idx = row_mapper.index(row_actives(r,0), patch_id, c_comp);
                        if (row_mapper.is_free_index(global_row_idx))
                        {
                            for (index_t c = 0; c < col_actives.rows(); ++c)
                            {
                                for (index_t r_comp = 0; r_comp < col_space_data->dim; ++r_comp)
                                {
                                    index_t global_col_idx = col_mapper.index(col_actives(c,0), patch_id, r_comp);
                                    if (col_mapper.is_free_index(global_col_idx))
                                    {
                                        // Interior-interior coupling: add to matrix
                                        m_fmatrix.coeffRef(global_row_idx, global_col_idx) += localMat(r, c);
                                    }
                                    else
                                    {
                                        // Interior-boundary coupling: Dirichlet lift
                                        // Subtract K(row, boundary) * dirichlet_value from RHS
                                        index_t bnd_idx = col_mapper.global_to_bindex(global_col_idx);
                                        if (bnd_idx < col_space_data->fixedDofs.rows())
                                        {
                                            m_rhs(global_row_idx) -= localMat(r, c) * col_space_data->fixedDofs(bnd_idx, 0);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else // Vector assembly (RHS)
            {
                // localMat is a column vector with dimensions (row_actives.rows(), 1)
                GISMO_ASSERT(localMat.rows() == row_actives.rows() && localMat.cols() == 1, "localMat must be a column vector with one entry per active basis function");
                for (index_t r = 0; r < row_actives.rows(); ++r)
                {
                    for (index_t c_comp = 0; c_comp < row_space_data->dim; ++c_comp)
                    {
                        index_t global_row_idx = row_mapper.index(row_actives(r,0), patch_id, c_comp);
                        if (row_mapper.is_free_index(global_row_idx))
                        {
                            m_rhs(global_row_idx) += localMat(r, 0);
                        }
                    }
                }
            }
        }
    };
};

} // namespace gismo
