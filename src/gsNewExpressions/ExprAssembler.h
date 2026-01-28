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

    Expr::SpaceObject<T, Expr::SpaceType::Test, 1> getVectorTestSpace(const gsFunctionSet<T>& basis, index_t id = 0, std::string label="v")
    {
        index_t target_dim = basis.targetDim();
        GISMO_ASSERT(id < m_test_spaces.size(), "Test space ID exceeds pre-allocated size.");

        m_space_data_list.emplace_back(basis, target_dim, id); // Vector space
        Expr::FeSpaceData<T>* space_data = &m_space_data_list.back();
        space_data->init();
        
        std::array<size_t,1> sizes = {(size_t)target_dim};
        Expr::SpaceObject<T, Expr::SpaceType::Test, 1> expr(basis.domainDim(), sizes, id, label);
        expr.setSource(basis);
        m_helper.add(expr); // Register with helper for eval()
        expr.setSpaceData(*space_data);
        
        m_test_spaces[id] = space_data;
        return expr;
    }

    Expr::SpaceObject<T, Expr::SpaceType::Trial, 1> getVectorTrialSpace(const gsFunctionSet<T>& basis, index_t id = 0, std::string label="u")
    {
        index_t target_dim = basis.targetDim();
        GISMO_ASSERT(id < m_trial_spaces.size(), "Trial space ID exceeds pre-allocated size.");
        
        m_space_data_list.emplace_back(basis, target_dim, id); // Vector space
        Expr::FeSpaceData<T>* space_data = &m_space_data_list.back();
        space_data->init();

        std::array<size_t,1> sizes = {(size_t)target_dim};
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
        initMatrix();
        m_rhs.setZero(numTestDofs(), numRhs);
    }
    
    void initMatrix()
    {
        resetDimensions();
        clearMatrix(false);
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


    // Getters for results
    const gsSparseMatrix<T>& matrix() const { return m_matrix; }
    const gsVector<T>& rhs() const { return m_rhs; }

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
        index_t current_shift = 0;
        for (auto* space_data : m_trial_spaces)
        {
            if (space_data && !space_data->mapper.isFinalized()) 
            {
                space_data->init();
                space_data->mapper.setShift(current_shift);
                current_shift += space_data->dim * space_data->mapper.freeSize();
            }
        }
        for (auto* space_data : m_test_spaces) // Test spaces might share mappers with trial spaces
        {
            if (space_data && !space_data->mapper.isFinalized()) 
            { // Only init if not already done by trial
                space_data->init();
                space_data->mapper.setShift(current_shift);
                current_shift += space_data->dim * space_data->mapper.freeSize();
            }
        }
        // Finalize all mappers if needed
        for (auto* space_data : m_trial_spaces) 
        {
            if (space_data && !space_data->mapper.isFinalized()) 
                space_data->mapper.finalize();
        }
        for (auto* space_data : m_test_spaces) 
        {
            if (space_data && !space_data->mapper.isFinalized()) 
                space_data->mapper.finalize();
        }
    }
    
    // Helper to get total number of DoFs (sum of free DoFs from all trial spaces)
    index_t numDofs() const
    {
        index_t total_dofs = 0;
        for(const auto* space_data : m_trial_spaces) {
            if (space_data && space_data->mapper.isFinalized()) {
                total_dofs += space_data->mapper.firstIndex() + space_data->mapper.freeSize();
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
                total_test_dofs += space_data->mapper.firstIndex() + space_data->mapper.freeSize();
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
                                        m_fmatrix.coeffRef(global_row_idx, global_col_idx) += localMat(r, c);
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
