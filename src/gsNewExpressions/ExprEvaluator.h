/** @file ExprEvaluator.h

    @brief Evaluates expressions from the gsNewExpressions module.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsNewExpressions/ExpressionHelper.h>
#include <gsNewExpressions/ExpressionResult.h>
#include <gsNewExpressions/NewExpressions.h>
#include <gsDomain/gsDomain.h> // For gsDomain
#include <gsAssembler/gsQuadrature.h> // For gsQuadrature

#include <vector>

namespace gismo
{
    template<class T>
    class ExprAssembler;

    template<class T>
    class ExprEvaluator
    {
    public:
        ExprEvaluator(ExpressionHelper<T>& helper, gsOptionList options = defaultOptions()) : m_helper(helper), m_options(options) {}
        
        // New constructor that takes an ExprAssembler
        ExprEvaluator(const ExprAssembler<T>& assembler, gsOptionList options = defaultOptions());

        static gsOptionList defaultOptions()
        {
            gsOptionList opt;
            opt.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
            opt.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
            return opt;
        }

        // template<class E>
        // gismo::Expr::ExpressionResult<T> eval(const Expr::BaseExpression<E>& expr,
        //                               const gsVector<T>& p)
        // {
        //     // single point evaluation
        //     m_helper.points() = gsMatrix<T>(p);
        //     expr.derived().parse(m_helper);
        //     m_helper.precompute();
        //     return expr.derived().eval(0);
        // }

        template<class E>
        std::vector<gismo::Expr::ExpressionResult<T>> eval(const Expr::BaseExpression<E>& expr,
                                                   const gsMatrix<T>& points)
        {
            m_helper.points() = points;
            expr.derived().parse(m_helper);
            m_helper.precompute();
            std::vector<gismo::Expr::ExpressionResult<T>> results;
            results.reserve(points.cols());
            for (index_t i = 0; i < points.cols(); ++i)
            {
                results.push_back(expr.derived().eval(i));
            }
            return results;
        }

        template<class E>
        gsMatrix<T> computeIntegral(const Expr::BaseExpression<E>& expr,
                                    const gsDomain<T>& domain) // Quadrature options
        {
            GISMO_ASSERT(expr.Space == Expr::SpaceType::None,
                         "computeIntegral: Only expressions with SpaceType::None are supported.");
            gsMatrix<T> integral_result;
            expr.derived().parse(m_helper);

            index_t rows = 1, cols = 1;
            if (expr.derived().order() == 1) {
                rows = expr.derived().sizes()[0];
            } else if (expr.derived().order() == 2) {
                rows = expr.derived().sizes()[0];
                cols = expr.derived().sizes()[1];
            }
            integral_result.setZero(rows, cols);

            #pragma omp parallel
            {
                gsMatrix<T> local_integral = gsMatrix<T>::Zero(rows, cols);
                typename gsQuadRule<T>::uPtr quad_rule;
                index_t current_patch = -1;

                #pragma omp for
                for (auto & elem : domain.allElements())
                {
                    if (current_patch != elem.patch())
                    {
                        current_patch = elem.patch();
                        // get Degree of the domain
                        quad_rule = gsQuadrature::getPtr(*domain.subdomain(current_patch), m_options);
                    }
                    
                    quad_rule->mapTo(elem.lowerCorner(), elem.upperCorner(), m_helper.points(), m_helper.weights());
                    
                    m_helper.precompute(current_patch);

                    for (index_t k = 0; k < m_helper.points().cols(); ++k)
                    {
                        gismo::Expr::ExpressionResult<T> eval_val = expr.derived().eval(k);
                        local_integral += eval_val(0,0) * m_helper.weights()[k];
                    }
                }
                
                #pragma omp critical
                integral_result += local_integral;
            }
            return integral_result;
        }


    private:
        ExpressionHelper<T>& m_helper;
        gsOptionList m_options;
    };


    template<class T>
    ExprEvaluator<T>::ExprEvaluator(const ExprAssembler<T>& assembler, gsOptionList options)
        : m_helper(const_cast<ExpressionHelper<T>&>(assembler.helper())), m_options(options) {}

} // namespace gismo
