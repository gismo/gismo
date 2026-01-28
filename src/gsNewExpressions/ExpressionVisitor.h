/** @file ExpressionVisitor.h

    @brief Base class for user-defined expression visitors

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

namespace gismo
{
namespace Expr
{

/**
 * @brief Base class for user-defined expression visitors
 * 
 * ExpressionVisitor provides a framework for creating custom bilinear forms
 * that can be assembled into matrices or vectors. Users can define their own
 * expressions by inheriting from this class and implementing the pure virtual
 * methods.
 * 
 * ## Purpose
 * 
 * Visitors allow users to:
 * - Define custom weak forms without modifying the core expression library
 * - Encapsulate complex PDEs into reusable components
 * - Add problem-specific operators (e.g., custom stabilization terms)
 * - Create libraries of standard forms (mass, stiffness, etc.)
 * 
 * ## Architecture
 * 
 * The visitor pattern separates the expression structure from the operations
 * performed on it. Each visitor:
 * - Knows how to evaluate itself at quadrature points
 * - Manages its own local matrix contributions
 * - Reports which function spaces it depends on (test/trial)
 * 
 * ## Usage Pattern
 * 
 * 1. **Define Your Visitor:**
 *    ```cpp
 *    class MyVisitor : public ExpressionVisitor<T> {
 *        // Implement pure virtual methods
 *    };
 *    ```
 * 
 * 2. **Register with Assembler:**
 *    ```cpp
 *    MyVisitor visitor(parameters);
 *    assembler.assemble(visitor);
 *    ```
 * 
 * 3. **Assembler Calls Methods:**
 *    - `parse()` - register dependencies
 *    - `eval()` - compute local contributions
 *    - `test()/trial()` - determine matrix structure
 * 
 * ## Example: Mass Matrix Visitor
 * 
 * ```cpp
 * template<class T>
 * class MassMatrixVisitor : public ExpressionVisitor<T> 
 * {
 *     SpaceObject<T, SpaceType::Test, 0> m_v;
 *     SpaceObject<T, SpaceType::Trial, 0> m_u;
 *     
 * public:
 *     MassMatrixVisitor(const SpaceObject<T, SpaceType::Test, 0> & v,
 *                       const SpaceObject<T, SpaceType::Trial, 0> & u)
 *     : m_v(v), m_u(u) {}
 *     
 *     T eval(index_t k) const override {
 *         return m_v.eval(k)(0,0) * m_u.eval(k)(0,0);
 *     }
 *     
 *     void parse(ExpressionHelper<T> & helper) const override {
 *         m_v.parse(helper);
 *         m_u.parse(helper);
 *     }
 *     
 *     auto test() const override { return m_v; }
 *     auto trial() const override { return m_u; }
 * };
 * ```
 * 
 * ## Virtual Method Contract
 * 
 * ### Pure Virtual (MUST implement):
 * 
 * 1. **`eval(index_t k)`**
 *    - Evaluates the expression at quadrature point k
 *    - Returns scalar value for bilinear form
 *    - Called by assembler during integration
 * 
 * 2. **`parse(ExpressionHelper<T> &)`**
 *    - Registers all dependencies with helper
 *    - Sets required flags (derivatives, etc.)
 *    - Called once before assembly
 * 
 * 3. **`test()`**
 *    - Returns the test space object
 *    - Determines row structure of matrix
 *    - Must match eval() behavior
 * 
 * 4. **`trial()`**
 *    - Returns the trial space object
 *    - Determines column structure of matrix
 *    - Must match eval() behavior
 * 
 * ### Optional (default implementations provided):
 * 
 * 1. **`print(std::ostream &)`**
 *    - Human-readable description
 *    - Used for debugging and logging
 *    - Default: prints class name
 * 
 * 2. **`order()`**
 *    - Returns tensor order (default: 0 for scalars)
 *    - 0 = scalar, 1 = vector, 2 = matrix
 * 
 * 3. **`space()`**
 *    - Space type: None, Test, Trial, Both
 *    - Default: Both (for bilinear forms)
 * 
 * ## Thread Safety
 * 
 * Visitors should be **const-correct** and **thread-safe**:
 * - eval() is const - use mutable for temporaries
 * - Store parameters as const members
 * - Avoid shared mutable state
 * 
 * ## Performance Tips
 * 
 * 1. **Cache Expensive Operations:**
 *    ```cpp
 *    mutable gsMatrix<T> tmp; // Reusable temporary
 *    ```
 * 
 * 2. **Minimize parse() Work:**
 *    - Only register what you need
 *    - Don't recompute in parse()
 * 
 * 3. **Inline eval() if Simple:**
 *    - Helps compiler optimize
 *    - Consider marking as `inline` or defining in header
 * 
 * @tparam T Scalar type (typically real_t or double)
 */
template<class T>
class ExpressionVisitor : public BaseExpression<ExpressionVisitor<T>>
{
    using Base = BaseExpression<ExpressionVisitor<T>>;
    
public:
    typedef T Scalar;
    static constexpr size_t Order = 0; // Visitors produce scalars by default
    static constexpr SpaceType Space = SpaceType::Both; // Test and trial
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = false;

    /**
     * @brief Virtual destructor for proper cleanup
     */
    virtual ~ExpressionVisitor() = default;

    // ============================================================
    // PURE VIRTUAL METHODS - Must be implemented by derived classes
    // ============================================================

    /**
     * @brief Evaluate the expression at quadrature point k
     * 
     * This is the core computation method. It should return the scalar
     * value of the bilinear form at the given quadrature point.
     * 
     * **Threading:** This method must be thread-safe and const.
     * Use mutable members for temporary storage if needed.
     * 
     * @param k Quadrature point index
     * @return Scalar value of the expression
     * 
     * @example
     * ```cpp
     * T eval(index_t k) const override {
     *     auto v_val = m_v.eval(k);
     *     auto u_val = m_u.eval(k);
     *     return v_val(0,0) * u_val(0,0); // Mass matrix
     * }
     * ```
     */
    virtual ExpressionResult<T> eval(index_t k) const = 0;

    /**
     * @brief Register dependencies with the expression helper
     * 
     * This method is called once before assembly begins. It should:
     * - Call parse() on all sub-expressions
     * - Set any required flags (NEED_DERIV, etc.)
     * - Register function sources
     * 
     * **Note:** This is NOT called repeatedly - cache what you can.
     * 
     * @param helper The expression helper managing evaluation data
     * 
     * @example
     * ```cpp
     * void parse(ExpressionHelper<T> & helper) const override {
     *     m_v.parse(helper);
     *     m_u.parse(helper);
     *     // If you need gradients:
     *     // m_v.setDerivative(1);
     *     // m_u.setDerivative(1);
     * }
     * ```
     */
    virtual void parse(gismo::ExpressionHelper<T> & helper) const = 0;

    /**
     * @brief Return the test space object
     * 
     * The test space determines the row structure of the assembled matrix.
     * Must be implemented by derived classes.
     * 
     * @return Test space object (determines matrix rows)
     * 
     * @example
     * ```cpp
     * SpaceObject<T, SpaceType::Test, 0> test() const override {
     *     return m_v;  // Return your test space
     * }
     * ```
     */
    virtual const SpaceObject<T, SpaceType::Test, 0>& test() const = 0;

    /**
     * @brief Return the trial space object
     * 
     * The trial space determines the column structure of the assembled matrix.
     * Must be implemented by derived classes.
     * 
     * @return Trial space object (determines matrix columns)
     * 
     * @example
     * ```cpp
     * SpaceObject<T, SpaceType::Trial, 0> trial() const override {
     *     return m_u;  // Return your trial space
     * }
     * ```
     */
    virtual const SpaceObject<T, SpaceType::Trial, 0>& trial() const = 0;

    // ============================================================
    // VIRTUAL METHODS - Can be overridden, have default implementations
    // ============================================================

    /**
     * @brief Print a human-readable representation
     * 
     * Used for debugging and logging. Default prints the class name.
     * Override to provide more detailed information about parameters.
     * 
     * @param os Output stream
     * 
     * @example
     * ```cpp
     * void print(std::ostream & os) const override {
     *     os << "MassMatrix(v=" << m_v.label() 
     *        << ", u=" << m_u.label() << ")";
     * }
     * ```
     */
    virtual void print(std::ostream & os) const
    {
        os << "ExpressionVisitor";
    }

    /**
     * @brief Return the tensor order of the result
     * 
     * - 0: Scalar (default for bilinear forms)
     * - 1: Vector
     * - 2: Matrix
     * 
     * Most visitors return 0 (scalars). Override if your visitor
     * produces vector or matrix results.
     * 
     * @return Tensor order
     */
    virtual size_t order() const { return Order; }

    /**
     * @brief Return the space type
     * 
     * - SpaceType::None: No space dependence
     * - SpaceType::Test: Only test space
     * - SpaceType::Trial: Only trial space
     * - SpaceType::Both: Test and trial (default)
     * 
     * @return Space type
     */
    virtual SpaceType space() const { return Space; }

    // ============================================================
    // REQUIRED INTERFACE - Already implemented, do not override
    // ============================================================

    /**
     * @brief Return the sizes array for the result
     * 
     * For scalar visitors (Order=0), this returns an empty array.
     * Automatically handled by the expression system.
     */
    std::array<size_t, Order> sizes() const
    {
        return std::array<size_t, Order>();
    }

    /**
     * @brief Return the domain dimension
     * 
     * Default implementation queries test space.
     * Override if needed for different behavior.
     */
    virtual size_t domainDim() const
    {
        return test().domainDim();
    }
};

// ============================================================
// TRAIT SPECIALIZATION for ExpressionVisitor
// ============================================================

/**
 * @brief Specialization of ExpressionTraits for ExpressionVisitor<T>
 * 
 * Every user-defined visitor inheriting from ExpressionVisitor<T> 
 * will use this trait specialization.
 */
template <typename T>
struct ExpressionTraits<ExpressionVisitor<T>>
{
    typedef T Scalar;
    static constexpr size_t Order = 0;           // Visitors produce scalars
    static constexpr SpaceType Space = SpaceType::Both;  // Contains both test and trial
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = false;
};

// ============================================================
// HELPER MACROS for common visitor patterns
// ============================================================

/**
 * @brief Macro to define a simple bilinear form visitor
 * 
 * Usage:
 * ```cpp
 * DEFINE_BILINEAR_VISITOR(MassMatrix, v * u)
 * ```
 */
#define DEFINE_BILINEAR_VISITOR(ClassName, Implementation) \
template<class T> \
class ClassName : public ExpressionVisitor<T> { \
    /* Implementation here */ \
};

} // namespace Expr
} // namespace gismo
