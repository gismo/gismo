/** @file ExpressionResult.h

    @brief Container for expression values with test/trial space cardinalities

    This file provides the ExpressionResult class, which handles the storage of evaluated
    expression results with proper support for test and trial space cardinalities.

    ## Motivation
    
    In finite element assembly, expressions can involve:
    - **Scalar expressions**: Independent of basis functions (e.g., constant, function values)
    - **Test space expressions**: Depend on test basis functions φ_i (i = 1..N)
    - **Trial space expressions**: Depend on trial basis functions ψ_j (j = 1..M)
    - **Bilinear expressions**: Depend on both test and trial basis functions (φ_i, ψ_j pairs)

    Each expression evaluation produces results for all relevant basis function combinations.
    ExpressionResult provides a unified container for these cases.

    ## Cardinality Mapping

    | Space Type       | Cardinality | Description                           | Example Expression |
    |------------------|-------------|---------------------------------------|-------------------|
    | SpaceType::None  | (1, 1)      | Single matrix (scalar expression)     | 2.0, sin(x)       |
    | SpaceType::Test  | (N, 1)      | N matrices (one per test function)    | φ, ∇φ             |
    | SpaceType::Trial | (1, M)      | M matrices (one per trial function)   | ψ, ∇ψ             |
    | SpaceType::Both  | (N, M)      | N×M matrices (one per basis pair)     | φ·ψ, ∇φ·∇ψ        |

    ## Usage Examples

    ```cpp
    // Scalar expression: value(0,0) contains the result at all quadrature points
    ExpressionResult<real_t> scalar_val = makeExpressionResult<real_t>(SpaceType::None);
    scalar_val(0, 0) = evaluateScalarExpression();

    // Test space: value(i,0) contains result for test basis function i
    ExpressionResult<real_t> test_val = makeExpressionResult<real_t>(SpaceType::Test, N);
    for (index_t i = 0; i < N; ++i)
        test_val(i, 0) = evaluateTestBasis(i);

    // Bilinear form: value(i,j) contains result for basis pair (φ_i, ψ_j)
    ExpressionResult<real_t> bilinear_val = makeExpressionResult<real_t>(SpaceType::Both, N, M);
    for (index_t i = 0; i < N; ++i)
        for (index_t j = 0; j < M; ++j)
            bilinear_val(i, j) = evaluateBilinearForm(i, j);
    ```

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <vector>

namespace gismo
{
namespace Expr
{

/**
 * @brief Container for expression evaluation results with test/trial space cardinalities
 * 
 * This class stores the evaluated values of an expression, accounting for the cardinality
 * of test and trial spaces:
 * - Scalar expressions (Space::None): cardinality (1,1) - single matrix
 * - Test space expressions (Space::Test): cardinality (N,1) - N matrices (one per test basis function)
 * - Trial space expressions (Space::Trial): cardinality (1,M) - M matrices (one per trial basis function)
 * - Bilinear expressions (Space::Both): cardinality (N,M) - N×M matrices (one per test-trial pair)
 * 
 * @tparam T Scalar type (e.g., real_t, double)
 */
template <typename T>
class ExpressionResult
{
public:
    typedef gsMatrix<T> MatrixType;

    /**
     * @brief Default constructor - creates empty value with cardinality (1,1)
     */
    ExpressionResult()
        : row_cardinality_(1), col_cardinality_(1)
    {
        data_.resize(1);
    }

    /**
     * @brief Constructor with specified cardinalities
     * @param row_card Row cardinality (number of test basis functions, or 1 for none)
     * @param col_card Column cardinality (number of trial basis functions, or 1 for none)
     */
    ExpressionResult(index_t row_card, index_t col_card)
        : row_cardinality_(row_card), col_cardinality_(col_card)
    {
        GISMO_ASSERT(row_card > 0 && col_card > 0, "Cardinalities must be positive");
        data_.resize(row_card * col_card);
    }

    /**
     * @brief Constructor with cardinalities and matrix dimensions
     * @param row_card Row cardinality
     * @param col_card Column cardinality
     * @param rows Number of rows in each matrix
     * @param cols Number of columns in each matrix
     */
    ExpressionResult(index_t row_card, index_t col_card, index_t rows, index_t cols)
        : row_cardinality_(row_card), col_cardinality_(col_card)
    {
        GISMO_ASSERT(row_card > 0 && col_card > 0, "Cardinalities must be positive");
        data_.resize(row_card * col_card);
        for (auto& mat : data_)
            mat.resize(rows, cols);
    }

    explicit ExpressionResult(const T & value)
        : row_cardinality_(1), col_cardinality_(1)
    {
        data_.resize(1);
        data_[0].resize(1, 1);
        data_[0](0, 0) = value;
    }

    template <int _Rows>
    explicit ExpressionResult(const gsVector<T, _Rows> & value)
        : row_cardinality_(1), col_cardinality_(1)
    {
        data_.resize(1);
        data_[0] = value;
    }

    template <int _Rows, int _Cols>
    explicit ExpressionResult(const gsMatrix<T, _Rows, _Cols> & value)
        : row_cardinality_(1), col_cardinality_(1)
    {
        data_.resize(1);
        data_[0] = value;
    }

    /**
     * @brief Access matrix for basis function pair (i,j)
     * @param i Row index (test basis function index, 0 for scalar/trial-only)
     * @param j Column index (trial basis function index, 0 for scalar/test-only)
     * @return Reference to the matrix for this basis pair
     */
    MatrixType& operator()(index_t i, index_t j)
    {
        GISMO_ASSERT(i >= 0 && i < row_cardinality_, "Row index out of bounds");
        GISMO_ASSERT(j >= 0 && j < col_cardinality_, "Column index out of bounds");
        return data_[i * col_cardinality_ + j];
    }

    /**
     * @brief Access matrix for basis function pair (i,j) - const version
     */
    const MatrixType& operator()(index_t i, index_t j) const
    {
        GISMO_ASSERT(i >= 0 && i < row_cardinality_, "Row index out of bounds");
        GISMO_ASSERT(j >= 0 && j < col_cardinality_, "Column index out of bounds");
        return data_[i * col_cardinality_ + j];
    }

    /**
     * @brief Access matrix for scalar expressions (cardinality 1,1)
     * @return Reference to the single matrix
     */
    MatrixType& operator()()
    {
        GISMO_ASSERT(row_cardinality_ == 1 && col_cardinality_ == 1, 
                    "Parameterless access only valid for scalar expressions");
        return data_[0];
    }

    /**
     * @brief Access matrix for scalar expressions - const version
     */
    const MatrixType& operator()() const
    {
        GISMO_ASSERT(row_cardinality_ == 1 && col_cardinality_ == 1, 
                    "Parameterless access only valid for scalar expressions");
        return data_[0];
    }

    /**
     * @brief Linear index access for advanced use
     * @param idx Linear index into the flattened data array
     */
    MatrixType& operator[](index_t idx)
    {
        GISMO_ASSERT(idx >= 0 && idx < static_cast<index_t>(data_.size()), "Index out of bounds");
        return data_[idx];
    }

    /**
     * @brief Linear index access - const version
     */
    const MatrixType& operator[](index_t idx) const
    {
        GISMO_ASSERT(idx >= 0 && idx < static_cast<index_t>(data_.size()), "Index out of bounds");
        return data_[idx];
    }

    /**
     * @brief Get row cardinality (number of test basis functions)
     */
    index_t rowCardinality() const { return row_cardinality_; }

    /**
     * @brief Get column cardinality (number of trial basis functions)
     */
    index_t colCardinality() const { return col_cardinality_; }

    /**
     * @brief Get total number of matrices stored
     */
    index_t size() const { return static_cast<index_t>(data_.size()); }

    /**
     * @brief Resize all matrices to new dimensions
     * @param rows New number of rows
     * @param cols New number of columns
     */
    void resizeMatrices(index_t rows, index_t cols)
    {
        for (auto& mat : data_)
            mat.resize(rows, cols);
    }

    /**
     * @brief Set all matrices to zero
     */
    void setZero()
    {
        for (auto& mat : data_)
            mat.setZero();
    }

    /**
     * @brief Set all matrices to a constant value
     */
    void setConstant(T value)
    {
        for (auto& mat : data_)
            mat.setConstant(value);
    }

    /**
     * @brief Resize the cardinality structure
     * @param row_card New row cardinality
     * @param col_card New column cardinality
     * @param preserve If true, preserve existing data where possible
     */
    void resize(index_t row_card, index_t col_card, bool preserve = false)
    {
        GISMO_ASSERT(row_card > 0 && col_card > 0, "Cardinalities must be positive");
        
        if (preserve && !data_.empty())
        {
            std::vector<MatrixType> new_data(row_card * col_card);
            index_t min_rows = std::min(row_cardinality_, row_card);
            index_t min_cols = std::min(col_cardinality_, col_card);
            
            for (index_t i = 0; i < min_rows; ++i)
                for (index_t j = 0; j < min_cols; ++j)
                    new_data[i * col_card + j] = data_[i * col_cardinality_ + j];
            
            data_ = std::move(new_data);
        }
        else
        {
            data_.resize(row_card * col_card);
        }
        
        row_cardinality_ = row_card;
        col_cardinality_ = col_card;
    }

    /**
     * @brief Flatten the data into a single matrix (for advanced use)
     * @return A single matrix containing all data concatenated
     */
    MatrixType flatten() const
    {
        // We assyme all matrices have the same size
        GISMO_ASSERT(!data_.empty(), "No data to flatten");
        index_t mat_rows = data_[0].rows();
        index_t mat_cols = data_[0].cols();
        index_t total_rows = mat_rows * row_cardinality_;
        index_t total_cols = mat_cols * col_cardinality_;
        MatrixType flat(total_rows, total_cols);
        for (index_t i = 0; i < row_cardinality_; ++i)
            for (index_t j = 0; j < col_cardinality_; ++j)
                flat.block(i * mat_rows, j * mat_cols, mat_rows, mat_cols) = data_[i * col_cardinality_ + j];

        return flat;
    }

    /**
     * @brief Get direct access to underlying data vector (advanced use)
     */
    std::vector<MatrixType>& data() { return data_; }
    const std::vector<MatrixType>& data() const { return data_; }

    /**
     * @brief Iterator support for range-based loops
     */
    typename std::vector<MatrixType>::iterator begin() { return data_.begin(); }
    typename std::vector<MatrixType>::iterator end() { return data_.end(); }
    typename std::vector<MatrixType>::const_iterator begin() const { return data_.begin(); }
    typename std::vector<MatrixType>::const_iterator end() const { return data_.end(); }

    std::ostream &print(std::ostream &os) const
    {
        os << "ExpressionResult with cardinality (" << row_cardinality_ << ", " << col_cardinality_ << "):\n";
        for (index_t i = 0; i < row_cardinality_; ++i)
        {
            for (index_t j = 0; j < col_cardinality_; ++j)
            {
                os << "Matrix (" << i << ", " << j << "):\n" << data_[i * col_cardinality_ + j] << "\n";
            }
        }
        return os;
    }

private:
    index_t row_cardinality_;    ///< Number of test basis functions (or 1 for scalar/trial-only)
    index_t col_cardinality_;    ///< Number of trial basis functions (or 1 for scalar/test-only)
    std::vector<MatrixType> data_; ///< Flattened array of matrices (row-major storage)
};

/// Helper class for scalar ExpressionResult. 
/// It can still have cardinality > 1, but each entry must be a scalar matrix (1x1).
template <typename T>
class ScalarExpressionResult : public ExpressionResult<T>
{
public:
    /// Constructor from scalar value (creates single-element result)
    explicit ScalarExpressionResult(T scalar_val)
        : ExpressionResult<T>(scalar_val)  // Use base class scalar constructor
    {
    }
    
    ScalarExpressionResult(const ExpressionResult<T>& expr_val)
        : ExpressionResult<T>(expr_val)
    {
        // Ensure all matrices are scalar (1x1)
        for (index_t i = 0; i < this->size(); ++i)
        {
            GISMO_ASSERT(this->operator[](i).rows() == 1 && this->operator[](i).cols() == 1,
                        "ScalarExpressionResult: All matrices must be scalar (1x1)");
        }
    }
};


template <typename T>
ExpressionResult<T> operator*(const T& scalar, const ExpressionResult<T>& expr)
{
    ExpressionResult<T> result = expr;
    for (index_t i = 0; i < result.size(); ++i)
        result[i] *= scalar;
    return result;
}

template <typename T>
ExpressionResult<T> operator*(const ExpressionResult<T>& expr, const T& scalar)
{
    ExpressionResult<T> result = expr;
    for (index_t i = 0; i < result.size(); ++i)
        result[i] *= scalar;
    return result;
}

template <typename T>
ExpressionResult<T> operator*(const ScalarExpressionResult<T>& lhs, const ExpressionResult<T>& rhs)
{
    ExpressionResult<T> result(math::max(lhs.rowCardinality(), rhs.rowCardinality()),
                              math::max(lhs.colCardinality(), rhs.colCardinality()));

    for (index_t i = 0; i < result.rowCardinality(); ++i)
        for (index_t j = 0; j < result.colCardinality(); ++j)
        {
            if (lhs.rowCardinality() == 1 && lhs.colCardinality() == 1)
                result(i, j) = lhs(0, 0)(0,0) * rhs(i, j);
            else if (rhs.rowCardinality() == 1 && rhs.colCardinality() == 1)
                result(i, j) = lhs(i, j)(0,0) * rhs(0, 0);
            else
                result(i, j) = lhs(i, j)(0,0) * rhs(i, j);
        }
    return result;
}

template <typename T>
ExpressionResult<T> operator*(const ExpressionResult<T>& lhs, const ScalarExpressionResult<T>& rhs)
{
    ExpressionResult<T> result(math::max(lhs.rowCardinality(), rhs.rowCardinality()),
                              math::max(lhs.colCardinality(), rhs.colCardinality()));
    for (index_t i = 0; i < result.rowCardinality(); ++i)
        for (index_t j = 0; j < result.colCardinality(); ++j)
        {
            if (lhs.rowCardinality() == 1 && lhs.colCardinality() == 1)
                result(i, j) = lhs(0, 0) * rhs(i, j)(0,0);
            else if (rhs.rowCardinality() == 1 && rhs.colCardinality() == 1)
                result(i, j) = lhs(i, j) * rhs(0, 0)(0,0);
            else
                result(i, j) = lhs(i, j) * rhs(i, j)(0,0);
        }
    return result;
}

template <typename T>
ExpressionResult<T> operator*(const ExpressionResult<T>& lhs, const ExpressionResult<T>& rhs)
{
    // This product is defined for:
    // - lhs.rowCardinality() == rhs.rowCardinality() && lhs.colCardinality() == rhs.colCardinality()
    // - (lhs.rowCardinality() == 1 && lhs.colCardinality() == 1) && rhs arbitrary
    // - lhs arbitrary && (rhs.rowCardinality() == 1 && rhs.colCardinality() == 1)
    // - lhs.rowCardinality() == rhs.rowCardinality() && (lhs.colCardinality() == 1 || rhs.colCardinality() == 1)
    // - lhs.colCardinality() == rhs.colCardinality() && (lhs.rowCardinality() == 1 || rhs.rowCardinality() == 1)
    // - lhs.rowCardinality() == 1 && rhs.colCardinality() == 1
    // - lhs.colCardinality() == 1 && rhs.rowCardinality() == 1

    ExpressionResult<T> result(math::max(lhs.rowCardinality(), rhs.rowCardinality()),
                              math::max(lhs.colCardinality(), rhs.colCardinality()));
    for (index_t i = 0; i < result.rowCardinality(); ++i)
        for (index_t j = 0; j < result.colCardinality(); ++j)
        {
            index_t lhs_i = (lhs.rowCardinality() == 1) ? 0 : i;
            index_t lhs_j = (lhs.colCardinality() == 1) ? 0 : j;
            index_t rhs_i = (rhs.rowCardinality() == 1) ? 0 : i;
            index_t rhs_j = (rhs.colCardinality() == 1) ? 0 : j;
            result(i, j) = lhs(lhs_i, lhs_j) * rhs(rhs_i, rhs_j);
        }
    return result;
}

template <typename T>
ExpressionResult<T> operator+(const ExpressionResult<T>& lhs, const ExpressionResult<T>& rhs)
{
    ExpressionResult<T> result(math::max(lhs.rowCardinality(), rhs.rowCardinality()),
                              math::max(lhs.colCardinality(), rhs.colCardinality()));
    for (index_t i = 0; i < result.rowCardinality(); ++i)
        for (index_t j = 0; j < result.colCardinality(); ++j)
        {
            if (lhs.rowCardinality() == 1 && lhs.colCardinality() == 1)
                result(i, j) = lhs(0, 0) + rhs(i, j);
            else if (rhs.rowCardinality() == 1 && rhs.colCardinality() == 1)
                result(i, j) = lhs(i, j) + rhs(0, 0);
            else
                result(i, j) = lhs(i, j) + rhs(i, j); 
        }
    return result;
}

template <typename T>
ExpressionResult<T> operator-(const ExpressionResult<T>& lhs, const ExpressionResult<T>& rhs)
{
    ExpressionResult<T> result(math::max(lhs.rowCardinality(), rhs.rowCardinality()),
                              math::max(lhs.colCardinality(), rhs.colCardinality()));
    for (index_t i = 0; i < result.rowCardinality(); ++i)
        for (index_t j = 0; j < result.colCardinality(); ++j)
        {
            if (lhs.rowCardinality() == 1 && lhs.colCardinality() == 1)
                result(i, j) = lhs(0, 0) - rhs(i, j);
            else if (rhs.rowCardinality() == 1 && rhs.colCardinality() == 1)
                result(i, j) = lhs(i, j) - rhs(0, 0);
            else
                result(i, j) = lhs(i, j) - rhs(i, j);
        }
    return result;
}

template <typename T>
ExpressionResult<T> operator/(const ExpressionResult<T>& expr, const T& scalar)
{
    ExpressionResult<T> result = expr;
    for (index_t i = 0; i < result.size(); ++i)
        result[i] /= scalar;
    return result;
}

template <typename T>  
ExpressionResult<T> operator/(const ExpressionResult<T>& expr, const ScalarExpressionResult<T>& scalar_expr)
{
    ExpressionResult<T> result(expr.rowCardinality(), expr.colCardinality());
    for (index_t i = 0; i < result.rowCardinality(); ++i)
        for (index_t j = 0; j < result.colCardinality(); ++j)
        {
            if (scalar_expr.rowCardinality() == 1 && scalar_expr.colCardinality() == 1)
                result(i, j) = expr(i, j) / scalar_expr(0, 0)(0,0);
            else
                result(i, j) = expr(i, j) / scalar_expr(i, j)(0,0);
        }
    return result;
}

/**
 * @brief Helper function to create ExpressionResult with appropriate cardinality based on space type
 * 
 * @tparam T Scalar type
 * @param space_type The space type (None, Test, Trial, Both)
 * @param test_card Number of test basis functions (ignored if not Test or Both)
 * @param trial_card Number of trial basis functions (ignored if not Trial or Both)
 * @return ExpressionResult with appropriate cardinality
 */
template <typename T>
ExpressionResult<T> makeExpressionResult(SpaceType space_type, index_t test_card = 1, index_t trial_card = 1)
{
    switch (space_type)
    {
        case SpaceType::None:
            return ExpressionResult<T>(1, 1);
        case SpaceType::Test:
            return ExpressionResult<T>(test_card, 1);
        case SpaceType::Trial:
            return ExpressionResult<T>(1, trial_card);
        case SpaceType::Both:
            return ExpressionResult<T>(test_card, trial_card);
        default:
            GISMO_ERROR("Unknown space type");
    }
}

/**
 * @brief Helper function to create ExpressionResult with matrix dimensions
 */
template <typename T>
ExpressionResult<T> makeExpressionResult(SpaceType space_type, index_t test_card, index_t trial_card, 
                                       index_t rows, index_t cols)
{
    switch (space_type)
    {
        case SpaceType::None:
            return ExpressionResult<T>(1, 1, rows, cols);
        case SpaceType::Test:
            return ExpressionResult<T>(test_card, 1, rows, cols);
        case SpaceType::Trial:
            return ExpressionResult<T>(1, trial_card, rows, cols);
        case SpaceType::Both:
            return ExpressionResult<T>(test_card, trial_card, rows, cols);
        default:
            GISMO_ERROR("Unknown space type");
    }
}

} // namespace Expr

template<class T>
std::ostream &operator<<(std::ostream &os, const Expr::ExpressionResult<T> &result)
{
    return result.print(os);
}

} // namespace gismo
